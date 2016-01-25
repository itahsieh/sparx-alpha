#! /usr/bin/env python

# Validate treatment of line radiative transfer in SPARX by doing synthetic
# observations on a spherical gas cloud containing two-level o-H2O in LTE, which is
# optically thin.
#
# Since the cloud is optically thin, the brightness temperature Tb for a cloud of
# kinetic temperature Tk and optical depth tau can be estimated as
#                    Tb = Tk * (1 - exp(-tau)) ~ Tk * tau
#
# Note that since the apbove expression is based on the Rayleigh-Jeans approximation,
# the difference between such an estimate and that calculated from synthetic
# observations will become smaller and smaller as the observation wavelegth
# becomes longer and longer.
#
# The optical depth tau at the center of a cloud with diameter D, molecular upper level
# density n1, lower level density n2, and at frequency nu can be calculated as
#            tau = alpha_nu * D
#                = ((h * nu)/(4 * pi)) * (n1 * B12 - n2 * B21) * phi(nu) * D
#
# where h is Planck's constant, B12 and B21 are the Einstein B coefficients, and phi(nu) 
# the line profile function, usually assumed to be Gaussian
#        phi(nu) = (1 / (delta_nu * sqrt(pi))) * exp(-(nu - nu0)**2 / delta_nu**2))
#
# where delta_nu is the Gaussian line width, and nu0 the line central frequency.
#
# In LTE, the level densities n1 and n2 are populated according to the Boltzmann distribution
#        n2 / n1 = (g2 * exp(-h * E2 / (kB * Tk))) / (g1 * exp(-h * E1 / (kB * Tk)))
#                = (g2 / g1) * exp(-(E2 - E1) / (kB * Tk))
#                = (g2 / g1) * exp(-h * nu / (kB * Tk))
#
# where g1 and g2 are the statistical weights of levels 1 and 2, E1 and E2 are the energies
# of levels 1 and 2, and kB is Boltzmann's constant.
#
# This calculation should produce the following table:

import sparx
import sparx.physics as phys

def Estimate(molec):
	from math import exp
	# Constants
	C = phys.Const
	U = phys.Units
	pc = U.pc # m

	# Cloud properties
	Tk = 40 # [K] Cloud kinetic temperature
	n_H2 = 1e9 * 1e6 # [m^-3] Molecular gas number density
	D = 2 * 0.1 * pc # [m] Cloud diameter
	Xmol = 1e-10

	# Molecule
	mol = phys.Molecule(molec)
	Blu = mol.line_Blu[0] # Inu^-1s^-1
	Bul = mol.line_Bul[0] # Inu^-1s^-1
	nu0 = mol.line_freq[0] # Hz
	delta_nu = mol.get_thermal_fwidth(0, Tk)
	phi_nu = phys.gaussian_fprofile(nu0, nu0, delta_nu)

	# Fractional level population
	u = mol.get_boltzmann_ratio(0, Tk)
	n_u = u / (1 + u)
	n_l = 1.0 - n_u

	# Calculate tau at center of cloud
	alpha_nu = (C.h * nu0 / (4.0 * C.pi)) * n_H2 * Xmol * (n_l * Blu - n_u * Bul) * phi_nu
	tau = alpha_nu * D

	# Calculate brightness temperature at center of cloud
	Tb = Tk * (1.0 - exp(-tau))

	return Tb, tau


##
## Main
##
if __name__ == "__main__":
	# Some necessary imports
	from sparx import tasks, utils, inputs
	from glob import glob

	# Remove old files
	name = "uniform"
	if not utils.confirm_remove_files(glob("./%s.*"%name)):
		print "Aborted"
		exit()

	# Table header
	table = ['# %20s %20s %20s %20s %20s %20s %20s %20s %20s'%('Molec', 'Tau_S1D', 'Tau_S3D', 'Tau_Est', 'Tau_diff(%)', 'Tb_S1D(K)', 'Tb_S3D(K)', 'Tb_Est(K)', 'Tb_Diff(%)')]

	# Calculate image for different wavelengths
	for molec in 'o-h2o_2lev', 'co_2lev', 'hco+_2lev', 'sio_2lev', 'hnco_2lev':
		# Setup file names
		src1d = "%s.%s.src1d"%(name, molec)
		out1d = "%s.%s.xyv1d"%(name, molec)
		tau1d = "%s.%s.tau1d"%(name, molec)
		src3d = "%s.%s.src3d"%(name, molec)
		out3d = "%s.%s.xyv3d"%(name, molec)
		tau3d = "%s.%s.tau3d"%(name, molec)

		# Reset inputs (safer)
		inputs.reset_inputs()

		# Calculate analytical solution
		Tb_est, tau_est = Estimate(molec)

		# Generate model
		tasks.task_valline1d(out=src1d, molec=molec, xmol=1e-10)
		tasks.task_valline3d(out=src3d, molec=molec, xmol=1e-10)

		# Generate synthetic dust observations
		tasks.task_lineobs(source=src1d, out=out1d, tau=tau1d, line="0", cell="['0.5asec', '0.5asec']", unit='K', chan="[32, '0.2kms^-1']")
		tasks.task_lineobs(source=src3d, out=out3d, tau=tau3d, line="0", cell="['0.5asec', '0.5asec']", unit='K', chan="[32, '0.2kms^-1']")

		from subprocess import Popen, PIPE
		# Get brightness temperature
		p = Popen("imspec in=%s options=eformat,nohead region='rel,box(0,0,0,0)' plot=sum | head -n 18 | tail -n 1 | awk '{print $3}'"%out1d, shell=True, stdout=PIPE, stderr=PIPE)
		Tb_s1d = float(p.communicate()[0])
		p = Popen("imspec in=%s options=eformat,nohead region='rel,box(0,0,0,0)' plot=sum | head -n 18 | tail -n 1 | awk '{print $3}'"%out3d, shell=True, stdout=PIPE, stderr=PIPE)
		Tb_s3d = float(p.communicate()[0])

		# Get peak tau
		p = Popen("imspec in=%s options=eformat,nohead region='rel,box(0,0,0,0)' plot=sum | head -n 18 | tail -n 1 | awk '{print $3}'"%tau1d, shell=True, stdout=PIPE, stderr=PIPE)
		tau_s1d = float(p.communicate()[0])
		p = Popen("imspec in=%s options=eformat,nohead region='rel,box(0,0,0,0)' plot=sum | head -n 18 | tail -n 1 | awk '{print $3}'"%tau3d, shell=True, stdout=PIPE, stderr=PIPE)
		tau_s3d = float(p.communicate()[0])

		tau_diff = 100 * abs(tau_s3d - tau_s1d) / tau_s1d
		Tb_diff = 100 * abs(Tb_s3d - Tb_s1d) / Tb_s3d

		table += ['  %20s %20.5e %20.5e %20.5e %20.5f %20.5e %20.5e %20.5e %20.5f'%(molec, tau_s1d, tau_s3d, tau_est, tau_diff, Tb_s1d, Tb_s3d, Tb_est, Tb_diff)]

# Print out results
for i in table:
	print i




