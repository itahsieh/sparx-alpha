#! /usr/bin/env python

# Validate treatment of dust radiative transfer in SPARX by doing synthetic
# observations on a spherical gas cloud with given gas to dust ratio and
# dust opacity, and comparing the results with a direct calculation.
#
# Since dust is mostly optically thin at mm/sub-mm wavelengths, the brightness
# temperature Tb for a cloud of kinetic temperature Tk and optical depth tau
# can be estimated as
#                    Tb = Tk * (1 - exp(-tau)) ~ Tk * tau
#
# Note that since the apbove expression is based on the Rayleigh-Jeans approximation,
# the difference between such an estimate and that calculated from synthetic
# observations will become smaller and smaller as the observation wavelegth
# becomes longer and longer.
#
# The optical depth tau at the center of a cloud with diameter D, column density
# N_H2, mean molecular weight mu, gas to dust ratio X, and dust mass opacity kappa
# can be calculated as
#                    tau = kappa * dust_surface_density
#                        = kappa * ((N_H2 / X) * mu * amu)
#                        = kappa * ((n_H2 * D / X) * mu * amu)
#
# This calculation should produce the following table:
#           Lambda(mm)            Tau_SPARX              Tau_Est          Tb_SPARX(K)            Tb_Est(K)           Tb_Diff(%)
#                  0.1          2.86900e-04          2.86938e-04          1.28400e-02          2.86938e-02             55.25165
#                    1          2.86900e-04          2.86938e-04          2.66700e-02          2.86938e-02              7.05307
#                   10          2.86900e-04          2.86938e-04          2.84800e-02          2.86938e-02              0.74508
#                  100          2.86900e-04          2.86938e-04          2.86700e-02          2.86938e-02              0.08292
#                 1000          2.86900e-04          2.86938e-04          2.86900e-02          2.86938e-02              0.01322
#                10000          2.86900e-04          2.86938e-04          2.86900e-02          2.86938e-02              0.01322
#
# Note that the kappa used here is *independent* of wavelength!!!

import sparx.physics as phys

def Estimate():
	# Constants
	amu = phys.Units.amu # kg
	pc = phys.Units.pc # m

	# Cloud properties
	Tk = 100 # [K] Cloud kinetic temperature
	gas2dust = 100.0 # Gas to dust ratio
	kappa = 0.1 / gas2dust # [m^2kg^-1] Dust mass opacity
	n_H2 = 1e4 * 1e6 # [m^-3] Molecular gas number density
	D = 2 * 0.1 * pc # [m] Cloud diameter

	# Calculate tau at center of cloud
	tau = kappa * n_H2 * D * 2.8 * amu

	# Calculate brightness temperature at center of cloud
	Tb = Tk * tau

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

	# Generate 1D model
	src1d = "%s1d.src"%(name)
	src3d = "%s3d.src"%(name)
	tasks.task_valdust1d(out=src1d)
	tasks.task_valdust3d(out=src3d)

	# Table header
	table = ['# %20s %20s %20s %20s %20s %20s %20s %20s %20s'%('Lambda(mm)', 'Tau_S1D', 'Tau_S3D', 'Tau_Est', 'Tau_Diff(%)', 'Tb_S1D(K)', 'Tb_S3D(K)', 'Tb_Est(K)', 'Tb_Diff(%)')]

	# Calculate image for different wavelengths
	for wavlen in 0.1, 1, 10, 100, 1000, 10000:
		# Setup file names
		out1d = "%s1d.%7.2emm.xyv"%(name, wavlen)
		tau1d = "%s1d.%7.2emm.tau"%(name, wavlen)
		out3d = "%s3d.%7.2emm.xyv"%(name, wavlen)
		tau3d = "%s3d.%7.2emm.tau"%(name, wavlen)

		# Reset inputs (safer)
		inputs.reset_inputs()

		# Generate synthetic dust observations
		tasks.task_contobs(source=src1d, out=out1d, tau=tau1d, wavelen="%7.2emm"%wavlen, cell="['0.5asec', '0.5asec']", unit='K', chan="[32, '0.2kms^-1']")
		tasks.task_contobs(source=src3d, out=out3d, tau=tau3d, wavelen="%7.2emm"%wavlen, cell="['0.5asec', '0.5asec']", unit='K', chan="[32, '0.2kms^-1']")

		from subprocess import Popen, PIPE
		# Extract brightness temperature
		p = Popen("imspec in=%s options=eformat,nohead region='rel,box(0,0,0,0)' plot=sum | head -n 2 | tail -n 1 | awk '{print $3}'"%out1d, shell=True, stdout=PIPE, stderr=PIPE)
		Tb_s1d = float(p.communicate()[0])

		p = Popen("imspec in=%s options=eformat,nohead region='rel,box(0,0,0,0)' plot=sum | head -n 2 | tail -n 1 | awk '{print $3}'"%out3d, shell=True, stdout=PIPE, stderr=PIPE)
		Tb_s3d = float(p.communicate()[0])

		# Get peak tau
		p = Popen("imspec in=%s options=eformat,nohead region='rel,box(0,0,0,0)' plot=sum | head -n 2 | tail -n 1 | awk '{print $3}'"%tau1d, shell=True, stdout=PIPE, stderr=PIPE)
		tau_s1d = float(p.communicate()[0])

		p = Popen("imspec in=%s options=eformat,nohead region='rel,box(0,0,0,0)' plot=sum | head -n 2 | tail -n 1 | awk '{print $3}'"%tau3d, shell=True, stdout=PIPE, stderr=PIPE)
		tau_s3d = float(p.communicate()[0])

		Tb_est, tau_est = Estimate()
		Tb_diff = 100 * abs(Tb_s3d - Tb_s1d) / Tb_s1d
		tau_diff = 100 * abs(tau_s3d - tau_s1d) / tau_s1d

		table += ['  %20g %20.5e %20.5e %20.5e %20.5e %20.5e %20.5e %20.5e %20.5f'%(wavlen, tau_s1d, tau_s3d, tau_est, tau_diff, Tb_s1d, Tb_s3d, Tb_est, Tb_diff)]

# Print out results
for i in table:
	print i



