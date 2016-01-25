#! /usr/bin/env python

# Validate non-LTE excitation calculation by solving the water
# radiative transfer benchmark problem presented at the RT workshop
# at Leiden in 2004. This problem considers the excitation of
# 2-level ortho-water in an isothermal, uniform gas cloud. Two
# dynamical cases are considered -- a static cloud and a dynmaical
# cloud where the cloud is expanding with a linear velocity gradient
# of 100 km/s/pc.
#
# This program does the necessary work to reproduce all the results
# written in Dr. David Neufeld's document for benchmarking non-LTE
# line radiative transfer codes.
# 
# CAVEAT: Be ware that the line width used in Neufeld's documents
# is 1 km/s FWHM, and if a direct comparison is to be made, some
# adjustments must be applied first.

# Global filename prefix
NAME = "leiden"

# Some necessary imports
from os.path import exists
from math import sqrt, log
import sparx
from sparx import tasks, utils, inputs, physics, miriad, grid
Unit = physics.Units
Cnst = physics.Const
import pylab as pl
import numpy as np

# Setup plotting
import matplotlib.font_manager
LEGND_PROP = matplotlib.font_manager.FontProperties(size=8)

# Image type
IMG_TYP = "png"

class Cloud(object):
	"""
	# The uniform cloud model
	"""
	# Cloud parameters
	r = 0.1 * Unit.pc # [m] cloud radius
	n_H2 = 1e4 * 1e6 # [m^-3] cloud density
	T_k = 40.0 # [K] cloud kinetic temperature

	# Distance from source [m]
	dist = 2000.0 * Unit.pc

	# Calculate projected beam radius [m]
	fwhm = (10.0 * Unit.pc / dist) # rad
	sigma = fwhm / (2.0 * sqrt(2.0 * log(2.0)))
	R_beam = sigma * dist

	# Channel number and width
	chan_n = 100
	chan_width = 0.05 # km/s

	def __init__(self, molec):
		# The model molecule should contain only two states
		#self.mol = physics.Molecule("o-h2o_2lev")
		self.mol = physics.Molecule(molec)

		# Calculate transition frequency and wavelength
		self.freq = self.mol.line_freq[0] # Hz
		self.lmda = Cnst.c / self.freq # m
		return

################################################################################

class Static(Cloud):
	"""
	Analytic solution for the static cloud:
	Since there is no accurate analytic solution for the excitation of a
	static cloud with arbitrary molecular abundance, only two limiting cases
	are solved:
	  a) the optcially thick LTE limit
	  b) the optically thin case
	
	Formulating detailed balance in terms of photon escape probability,
	we may arrive at the following solution for the ratio of the upper
	and lower level populations
	
	n_u / n_l = u / (1 + s * beta)
	where u = (g_u / g_l) * exp(-(E_u - E_l) / (k * T_k))
	      s = n_cr / n_H2
	      beta = escape probability
	"""
	# Filename prefix
	name = NAME+"-static"

	def __init__(self, ndiv, Xmol_list, molec, tcmb="0K", s1d_path=None, s3d_path=None, r1d_path=None, fwhm=0.32e3):
		"""
		Initialize static problem: the only free parameter is
		the fwhm of line width.
		"""
		# Initialize parent class
		Cloud.__init__(self, molec)

		# Number of divisions along radial direction
		self.ndiv = ndiv

		# List of molecular abundances
		self.Xmol_list = Xmol_list

		# Background radiation
		self.tcmb = tcmb

		# Paths to calculation results
		self.s1d_path = s1d_path
		self.s3d_path = s3d_path
		self.r1d_path = r1d_path

		# Setup file names
		if s1d_path != None:
			name = s1d_path+"/"+"%s-s1d"%(self.name)
			self.s1d_srcs = [name+"-X=%.2e.src"%Xmol for Xmol in Xmol_list]
			self.s1d_pops = [name+"-X=%.2e.pop"%Xmol for Xmol in Xmol_list]
			self.s1d_imgs = [name+"-X=%.2e.img"%Xmol for Xmol in Xmol_list]
			self.s1d_cnvs = [name+"-X=%.2e.cnv"%Xmol for Xmol in Xmol_list]
			self.s1d_tmbs = [name+"-X=%.2e.tmb"%Xmol for Xmol in Xmol_list]
		if s3d_path != None:
			name = s3d_path+"/"+"%s-s3d"%(self.name)
			self.s3d_srcs = [name+"-X=%.2e.src"%Xmol for Xmol in Xmol_list]
			self.s3d_pops = [name+"-X=%.2e.pop"%Xmol for Xmol in Xmol_list]
			self.s3d_imgs = [name+"-X=%.2e.img"%Xmol for Xmol in Xmol_list]
			self.s3d_cnvs = [name+"-X=%.2e.cnv"%Xmol for Xmol in Xmol_list]
			self.s3d_tmbs = [name+"-X=%.2e.tmb"%Xmol for Xmol in Xmol_list]

		# Free parameters:
		# Convert line width from FWHM to sqrt(2) * sigma:
		# FWHM = 2 * sqrt(2 * log(2)) * sigma
		self.width = fwhm / (2.0 * sqrt(log(2.0))) # m/s
		self.nu_D = physics.Doppler_vel2frq(self.freq, self.width) - self.freq # Hz

		# Excitation parameters:
		# The 'u' parameter for excitation: this should be 0.512 for the H2O model
		self.exc_u = self.mol.get_boltzmann_ratio(0, self.T_k)

		# The 's' parameter for excitation:
		# This should be 1588 for the H2O model
		self.exc_s = self.mol.col[0].get_crit_dens(0, self.T_k) / self.n_H2

		# Level populations:
		# Upper level fractional density at optically thick condition (beta = 0)
		#    n_u / n_l = u / (1 + s * beta)
		# => n_u / (1 - n_u) = u / 1
		# => n_u = u / (u + 1)
		# This should be 0.339 for the H2O model
		self.n_u_thick = self.exc_u / (self.exc_u + 1.0)

		# Lower level fractional density at optically thin condition (beta = 1)
		#    n_u / n_l = u / (1 + s * beta)
		# => n_u / n_l = u / (1 + s)
		# This should be 3.22e-4 for the H2O model
		self.n_u_thin = 1.0 / (1.0 + (1.0 + self.exc_s) / self.exc_u)
		return

	def phi_nu(self, nu):
		"""
		Line profile function
		"""
		return phys.gaussian_fprofile(nu, self.freq, self.nu_D)

	def calc_luminosity(self, X_mol):
		"""
		Calculate theoretical luminosity assuming collisional de-excitation
		can be neglected.
		L = h * nu * q_12 * n_H2 * n_H2O * (4/3 * pi * r**3)
		"""
		return Cnst.h * self.freq * self.mol.col[0].get_up_rate(0, self.T_k) * self.n_H2 * (X_mol * self.n_H2) * (Cnst.pi * 4.0 / 3.0) * self.r**3.0

	def _pipeline_sparx(self, Xmol, src, pop, img, cnv, tmb, genimg=True):
		"""
		Pipeline for generating SPARX solutions
		"""
		# Calculate excitation
		if not exists(pop):
			#tasks.task_amc(source=src, molec='o-h2o_2lev', out=pop)
			utils.call("mpirun -np %d sparx --parallel run task_amc source='%s' molec='%s' out='%s' fixiter=10 lte=%s trace=%s nrays=%d snr=%.0f tolerance=%.2e maxiter=%d"%(NPROC, src, self.mol.name, pop, FROMLTE, TRACE, NRAYS, SNR, TOLERANCE, MAXITER))
			#utils.call("sparx run task_amc source='%s' molec='o-h2o_2lev' out='%s'"%(src, pop))

		if genimg:
			# Generate synthesized map
			if not exists(img):
				tasks.task_lineobs(
					source=pop,
					chan="[%d,'%gkms^-1']"%(self.chan_n, self.chan_width),
					dist="%gpc"%(self.dist / Unit.pc),
					line=0,
					out=img,
					npix="[256,256]",
					cell="['0.5asec','0.5asec']",
					unit="JY/PIXEL")

			# To avoid mysterious 'invalid argument' bug in miriad
			from time import sleep
			sleep(0.5)

			# Convolve with beam
			if not exists(cnv):
				miriad.convol(img, cnv, self.fwhm)

			# To avoid mysterious 'invalid argument' bug in miriad
			from time import sleep
			sleep(0.5)

			# Convert to brightness temperature
			if not exists(tmb):
				miriad.convert_flux_to_tb(cnv, tmb, self.freq, self.fwhm)
		return

	def _pipeline_s1d(self, iX, vgrad, genimg=True):
		"""
		SPARX-1D pipeline
		"""
		# Setup filenames
		src = self.s1d_srcs[iX]
		pop = self.s1d_pops[iX]
		img = self.s1d_imgs[iX]
		cnv = self.s1d_cnvs[iX]
		tmb = self.s1d_tmbs[iX]
		Xmol = self.Xmol_list[iX]

		# Reset inputs (safer)
		inputs.reset_inputs()

		# Generate model grid
		if not exists(src):
			if self.ndiv == None:
				tasks.task_leiden1d(out=src, xmol=Xmol, vgrad=vgrad, tk=opts.tk, tcmb=self.tcmb)
			else:
				tasks.task_leiden1d(out=src, xmol=Xmol, vgrad=vgrad, tk=opts.tk, tcmb=self.tcmb, ndiv=self.ndiv)

		self._pipeline_sparx(Xmol, src, pop, img, cnv, tmb, genimg)
		return

	def _pipeline_s3d(self, iX, vgrad, genimg=True):
		"""
		SPARX-3D pipeline
		"""
		# Setup filenames
		src = self.s3d_srcs[iX]
		pop = self.s3d_pops[iX]
		img = self.s3d_imgs[iX]
		cnv = self.s3d_cnvs[iX]
		tmb = self.s3d_tmbs[iX]
		Xmol = self.Xmol_list[iX]

		# Reset inputs (safer)
		inputs.reset_inputs()

		# Generate model grid
		if not exists(src):
			if self.ndiv == None:
				ndiv = 16
			else:
				ndiv = self.ndiv * 2
			tasks.task_leiden3d(out=src, xmol=Xmol, vgrad=vgrad, tk=opts.tk, tcmb=self.tcmb, ndiv=ndiv)

		self._pipeline_sparx(Xmol, src, pop, img, cnv, tmb, genimg)
		return

	def plot_figure1(self, iX_range=None):
		"""
		Figure 1: shows the upper level population calculated with SPARX as a
		function of radius, for different molecular abundances
		"""
		# Fall back to full range if iX_range not given
		if iX_range == None:
			iX_range = range(len(self.Xmol_list))

		# Clear current figure
		pl.cla()

		for i in iX_range:
			Xmol = self.Xmol_list[i]
			if self.s1d_path :
				# Get and plot S1D results
				h5f = grid.SPARXH5(self.s1d_pops[i])
				xlist = h5f.GetRadii()
				ylist = h5f.GetRadial("lev1")
				h5f.Close()
				pl.plot(xlist, ylist, "-o", label="X(H2O)=%.2e (S1D)"%Xmol)
			if self.s3d_path :
				# Get and plot S3D results
				h5f = grid.SPARXH5(self.s3d_pops[i])
				xlist = h5f.GetRadii()
				ylist = h5f.GetRadial("lev1")
				h5f.Close()
				pl.plot(xlist, ylist, "-^", label="X(H2O)=%.2e (S3D)"%Xmol)

		# Plot analytic solution
		pl.axhline(self.n_u_thick, color="r", ls=":", label="LTE limit")
		pl.axhline(self.n_u_thin, color="b", ls=":", label="Optically thin limit")

		# Setup plot
		pl.xscale("log")
		pl.yscale("log")
		pl.xlabel("Radius [pc]")
		pl.ylabel("Upper level fractional density")
		pl.legend(loc="best", prop=LEGND_PROP)
		pl.xlim((0, self.r / Unit.pc)) # pc

		# Save plot
		pl.savefig(self.name+"-fig1."+IMG_TYP)
		return

	def plot_figure2(self, iX_range=None):
		"""
		Figure 2: shows the upper level population at the center of the cloud
		calculated with SPARX as a function of abundance
		"""
		# Fall back to full range if iX_range not given
		if iX_range == None:
			iX_range = range(len(self.Xmol_list))

		# Clear figure
		pl.cla()

		# Allocate arrays for plotting
		Xmol_arr = np.array(self.Xmol_list)
		if self.s1d_path :
			n_arr_s1d = np.zeros(shape=(len(self.Xmol_list)))
		if self.s3d_path :
			n_arr_s3d = np.zeros(shape=(len(self.Xmol_list)))

		# Loop through abundance list and gather n_u at
		# central zone
		for i in iX_range:
			if self.s1d_path :
				# Get S1D results
				h5f = grid.SPARXH5(self.s1d_pops[i])
				levarr = h5f.GetRadial("lev1")
				n_arr_s1d[i] = levarr[0]
				h5f.Close()
			if self.s3d_path :
				# Get S3D results
				h5f = grid.SPARXH5(self.s3d_pops[i])
				levarr = h5f.GetRadial("lev1")
				n_arr_s3d[i] = levarr[0]
				h5f.Close()

		# Plot analytic solution
		pl.axhline(self.n_u_thick, color="r", ls=":", label="LTE limit")
		pl.axhline(self.n_u_thin, color="b", ls=":", label="Optically thin limit")

		# Plot SPARX solutions
		if self.s1d_path:
			pl.plot(Xmol_arr, n_arr_s1d, "c-o", label="S1D")

		if self.s3d_path:
			pl.plot(Xmol_arr, n_arr_s3d, "m-^", label="S3D")

		# Setup plot
		pl.xscale("log")
		pl.yscale("log")
		pl.xlabel("Molecular abundance")
		pl.ylabel("Upper level fractional density")
		pl.legend(loc="best", prop=LEGND_PROP)
		pl.xlim((Xmol_arr[0], Xmol_arr[-1]))

		# Save plot
		pl.savefig(self.name+"-fig2."+IMG_TYP)
		return

	def plot_figure3(self, iX_range=None):
		"""
		Figure 3: shows the emergent spectrum (T_A vs. v) for a Gaussian beam of
		projected size 10pc (FWHM)
		"""
		# Fall back to full range if iX_range not given
		if iX_range == None:
			iX_range = range(len(self.Xmol_list))

		# Clear figure
		pl.cla()

		# Loop through abundance list
		for i in iX_range:
			Xmol = self.Xmol_list[i]
			if self.s1d_path :
				# Get and plot S1D results
				tmb = miriad.MirXYV(self.s1d_tmbs[i])
				velo = tmb.v_list / 1e3 # [m/s] -> [km/s]
				spec = tmb.GetSpecOffASec(0, 0) # [K]
				pl.plot(velo, spec, ":", label="X(H2O)=%.2e (S1D)"%Xmol)
			if self.s3d_path :
				# Get and plot S3D results
				tmb = miriad.MirXYV(self.s3d_tmbs[i])
				velo = tmb.v_list / 1e3 # [m/s] -> [km/s]
				spec = tmb.GetSpecOffASec(0, 0) # [K]
				pl.plot(velo, spec, "-.", label="X(H2O)=%.2e (S3D)"%Xmol)

		# Setup plot
		pl.yscale("log")
		pl.xlabel("Projected velocity (km/s)")
		pl.ylabel("Antenna temperature (K)")
		pl.legend(loc="best", prop=LEGND_PROP)
		pl.xlim(-1.5, 1.5)
		#pl.ylim(1e-7, 1e-2)

		# Save plot
		pl.savefig(self.name+'-fig3.'+IMG_TYP)
		return

	def plot_figure4(self, iX_range=None):
		"""
		Figure 4: shows the total line luminosity as a function of the molecular
		abundance
		"""
		# Fall back to full range if iX_range not given
		if iX_range == None:
			iX_range = range(len(self.Xmol_list))

		# Clear figure
		pl.cla()

		# Setup data arrays
		L_list_theory = []
		L_list_s1d = []
		L_list_s3d = []

		# Loop through abundance list
		for i in iX_range:
			Xmol = self.Xmol_list[i]
			# Calculate theoretical luminosity limit
			L_list_theory += [self.calc_luminosity(Xmol) * 1e7] # [J s^-1] -> [erg s^-1]

			if self.s1d_path:
				cnv = miriad.MirXYV(self.s1d_cnvs[i])
				F_nu = cnv.GetSpecOffASec(0, 0) # [Jy/Beam]
				F_line = sum(F_nu) * 1e-26 # [J/m^2]
				L_line = F_line * 4.0 * Cnst.pi * self.R_beam**2 * 1e7 # [J s^-1] -> [erg s^-1]
				L_list_s1d += [L_line]
			if self.s3d_path:
				cnv = miriad.MirXYV(self.s3d_cnvs[i])
				F_nu = cnv.GetSpecOffASec(0, 0) # [Jy/Beam]
				F_line = sum(F_nu) * 1e-26 # [J/m^2]
				L_line = F_line * 4.0 * Cnst.pi * self.R_beam**2 * 1e7 # [J s^-1] -> [erg s^-1]
				L_list_s3d += [L_line]

		# Plot all solutions
		Xmol_list = [self.Xmol_list[i] for i in iX_range]
		pl.plot(Xmol_list, L_list_theory, ":", label="Theoretical limit")
		if len(L_list_s1d):
			pl.plot(Xmol_list, L_list_s1d, "o", label="S1D")
		if len(L_list_s3d):
			pl.plot(Xmol_list, L_list_s3d, "^", label="S3D")

		# Setup plot
		pl.xscale("log")
		pl.yscale("log")
		pl.xlabel("Molecular fractional abundance")
		pl.ylabel("Line luminosity (erg/s)")
		pl.legend(loc="best", prop=LEGND_PROP)
		pl.xlim(min(Xmol_list), max(Xmol_list))
		# pl.ylim(1e26, 1e32)

		# Save plot
		pl.savefig(self.name+'-fig4.'+IMG_TYP)
		return

	def run(self, exc_only=False, nofig=False, no_intermediate=False):
		"""
		Run the benchmark problems
		"""
		# Don't generate images if excitation only is requested
		genimg = not exc_only

		# Calculate static problem for all abundances for all pipelines
		for i in range(len(self.Xmol_list)):
			print "Running static problem for Xmol=%g..."%self.Xmol_list[i]
			# S1D problem
			if self.s1d_path != None:
				self._pipeline_s1d(i, '0kms^-1', genimg)
			# S3D problem
			if self.s3d_path != None:
				self._pipeline_s3d(i, '0kms^-1', genimg)

			# Plot intermediate results
			if not nofig and not no_intermediate:
				print "Generating intermediate figures..."
				self.plot_figure1(range(0,i+1))
				self.plot_figure2(range(0,i+1))

				if genimg:
					self.plot_figure3(range(0,i+1))
					self.plot_figure4(range(0,i+1))

		if not nofig:
			print "Generating final figures..."
			self.plot_figure1(range(0,i+1))
			self.plot_figure2(range(0,i+1))

			if genimg:
				self.plot_figure3(range(0,i+1))
				self.plot_figure4(range(0,i+1))
		print "Static problem completed"
		print
		return

################################################################################

class LVG(Static):
	"""
	Analytic solution for the expanding cloud:
	The expanding cloud problem is solved using Sobolev's method,
	which gives analytical solutions for the level populations.
	"""
	def __init__(self, ndiv, Xmol_list, molec, alpha=100.0e3 / Unit.pc, fwhm=0.32e3, **kwargs):
		self.alpha = alpha
		self.name = NAME+"-lvg"
		Static.__init__(self, ndiv, Xmol_list, molec, fwhm=fwhm, **kwargs)
		return

	def exc_t(self, Xmol):
		"""
		Calculate the `t' optical depth parameter for the LVG problem.
		alpha = v / r (velocity gradient)
		"""
		return (1.0 / (8.0 * Cnst.pi * self.alpha)) * (Cnst.c**3.0 / self.freq**3.0) * self.mol.line_Aul[0] * self.n_H2 * Xmol

	def n_u_hightau(self, x_mol):
		"""
		Calculate upper level fractional population for
		high-optical depth conditions
		"""
		from math import sqrt
		t = self.exc_t(x_mol)
		u = self.exc_u
		s = self.exc_s
		return (t * (3.0 * u + 1.0) + s - sqrt((1.0 - u)**2.0 * t**2.0 + 2.0 * t * s * (3.0 * u + 1.0) + s**2.0)) / (4.0 * t * (u + 1.0))

	def n_u_lowtau(self, x_mol):
		"""
		Calculate upper level fractional population for
		low-optical depth conditions
		"""
		from math import exp
		t = self.exc_t(x_mol)
		u = self.exc_u
		s = self.exc_s
		return u * t / (s * (1.0 - exp(-t)))

	def n_u_analytic(self, X_mol):
		if X_mol <= 5.591e-8:
			return self.n_u_lowtau(X_mol)
		else:
			return self.n_u_hightau(X_mol)

	def plot_figure2(self, iX_range=None):
		"""
		Figure 2: shows the upper level population at the center of the cloud
		calculated with SPARX as a function of abundance
		"""
		# Fall back to full range if iX_range not given
		if iX_range == None:
			iX_range = range(len(self.Xmol_list))

		# Clear figure
		pl.cla()

		# Allocate arrays for plotting
		Xmol_arr = np.array(self.Xmol_list)
		if self.s1d_path :
			n_arr_s1d = np.zeros(shape=(len(self.Xmol_list)))
		if self.s3d_path :
			n_arr_s3d = np.zeros(shape=(len(self.Xmol_list)))

		# Loop through abundance list and gather n_u at
		# central zone
		for i in iX_range:
			if self.s1d_path :
				# Get S1D results
				h5f = grid.SPARXH5(self.s1d_pops[i])
				levarr = h5f.GetRadial("lev1")
				n_arr_s1d[i] = levarr[0]
				h5f.Close()
			if self.s3d_path :
				# Get S3D results
				h5f = grid.SPARXH5(self.s3d_pops[i])
				levarr = h5f.GetRadial("lev1")
				n_arr_s3d[i] = levarr[0]
				h5f.Close()

		# Plot analytic solutions
		pl.axhline(self.n_u_thick, color="r", ls=":", label="LTE limit")
		pl.axhline(self.n_u_thin, color="b", ls=":", label="Optically thin limit")
		xarr = utils.generate_log_points(Xmol_arr[0], Xmol_arr[-1], 100)
		pl.plot(xarr, [self.n_u_analytic(Xmol) for Xmol in xarr], color="g", ls="-.", label="Sobolev approximation")

		# Plot SPARX solutions
		if self.s1d_path:
			pl.plot(Xmol_arr, n_arr_s1d, "c-o", label="S1D")
		if self.s3d_path:
			pl.plot(Xmol_arr, n_arr_s3d, "m-^", label="S3D")

		# Setup plot
		pl.xscale("log")
		pl.yscale("log")
		pl.xlabel("Molecular abundance")
		pl.ylabel("Upper level fractional density")
		pl.legend(loc="best", prop=LEGND_PROP)
		pl.xlim((Xmol_arr[0], Xmol_arr[-1]))

		# Save plot
		pl.savefig(self.name+"-fig2."+IMG_TYP)
		return

	def run(self, vgrad="100kms^-1", nofig=False, no_intermediate=False):
		"""
		Run the benchmark problems
		"""
		# Calculate static problem for all abundances for all pipelines
		for i in range(len(self.Xmol_list)):
			print "Running LVG problem for Xmol=%g..."%self.Xmol_list[i]
			# S1D problem
			if self.s1d_path != None:
				self._pipeline_s1d(i, vgrad, genimg=False)
			# S3D problem
			if self.s3d_path != None:
				self._pipeline_s3d(i, vgrad, genimg=False)

			# Plot intermediate results
			if not nofig and not no_intermediate:
				print "Generating intermediate figures..."
				self.plot_figure1(range(0,i+1))
				self.plot_figure2(range(0,i+1))

		if not nofig:
			print "Generating final figures..."
			self.plot_figure1(range(len(self.Xmol_list)))
			self.plot_figure2(range(len(self.Xmol_list)))
		print "LVG problem completed"
		print
		return

################################################################################

##
## Main
##
if __name__ == "__main__":
	##
	## Parse inputs
	##
	# Init option parser
	from optparse import OptionParser
	parser = OptionParser()
	 
	# Setup options
	parser.add_option("--ndiv", metavar="POSINT", dest="ndiv", nargs=1, default="8", help="Number of zones along radial direction")
	parser.add_option("--xmol", metavar="XLIST", dest="xmol", nargs=1, default="1e-10,1e-9,1e-8,1e-7,1e-6", help="List of abundances to calcualte")
	parser.add_option("--tk", metavar="TK", dest="tk", nargs=1, default="40K", help="Kinetic temperature")
	parser.add_option("--vgrad", metavar="VELGRAD", dest="vgrad", nargs=1, default="100kms^-1", help="Velocity gradient of LVG problem")
	parser.add_option("--tcmb", metavar="TCMB", dest="tcmb", nargs=1, default="0K", help="Brightness temperature of background radiation")
	parser.add_option("--snr", metavar="SNR", dest="snr", nargs=1, default="20", help="Final Monte Carlo S/N ratio")
	parser.add_option("--tolerance", metavar="TOLERANCE", dest="tolerance", nargs=1, default="1e-9", help="Convergence criterion for fixed rays stage")
	parser.add_option("--nrays", metavar="NRAYS", dest="nrays", nargs=1, default="1000", help="Number of initial rays")
	parser.add_option("--maxiter", metavar="MAXITER", dest="maxiter", nargs=1, default="1000", help="Maximum number of iterations for converging Jbar and n")
	parser.add_option("--molec", metavar="MOLEC", dest="molec", nargs=1, default="o-h2o_2lev", help="Molecule used in the benchmark")
	parser.add_option("--clear", dest="clear", action="store_true", default=False, help="Whether to remove old files")
	parser.add_option("--trace", dest="trace", action="store_true", default=False, help="Whether to trace convergence history")
	parser.add_option("--from-lte", dest="from_lte", action="store_true", default=False, help="Start convergence from LTE conditions")
	parser.add_option("--static-only", dest="static_only", action="store_true", default=False, help="Do static problem only")
	parser.add_option("--lvg-only", dest="lvg_only", action="store_true", default=False, help="Do LVG problem only")
	parser.add_option("--exc-only", dest="exc_only", action="store_true", default=False, help="Calculate excitation only")
	parser.add_option("--nofig", dest="nofig", action="store_true", default=False, help="Do not plot figures")
	parser.add_option("--no-intermediate", dest="no_intermediate", action="store_true", default=False, help="Do not plot intermediate figures")
	parser.add_option("--orig", dest="orig", action="store_true", default=False, help="Use original problem description")
	parser.add_option("--1d-only", dest="only1d", action="store_true", default=False, help="Do 1D problem only")
	parser.add_option("--np", metavar="NPROC", dest="nproc", nargs=1, default="1", help="Number of parallel processes")
	 
	# The actual parsing
	(opts, args) = parser.parse_args()

	# The list of molecular abundances to be considered
	#Xmol_LIST = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5][:]
	Xmol_LIST = [float(i) for i in opts.xmol.split(",")]
	for i in Xmol_LIST:
		if i <= 0:
			raise Exception, "xmol must be > 0"

	# SNR
	SNR = float(opts.snr)

	# NRAYS
	NRAYS = int(opts.nrays)

	# FROMLTE
	FROMLTE = opts.from_lte

	# TRACE
	TRACE = opts.trace

	# TOLERANCE
	TOLERANCE = float(opts.tolerance)

	# NPROC
	NPROC = int(opts.nproc)

	# MAXITER
	MAXITER = int(opts.maxiter)

	# NDIV
	if opts.orig:
		NDIV = None
	else:
		NDIV = int(opts.ndiv)

	# Clear old files?
	if opts.clear:
		from glob import glob
		utils.confirm_remove_files(glob("./%s*"%NAME))

	if opts.only1d:
		s3d_path = None
	else:
		s3d_path = "."

	##
	## Run validation
	##
	# Calculate static problem
	if not opts.lvg_only:
		static = Static(NDIV, Xmol_LIST, opts.molec, tcmb=opts.tcmb, s1d_path=".", s3d_path=s3d_path)
		static.run(opts.exc_only, opts.nofig, opts.no_intermediate)

	# Calculate LVG problem
	if not opts.static_only:
		lvg = LVG(NDIV, Xmol_LIST, opts.molec, tcmb=opts.tcmb, s1d_path=".", s3d_path=s3d_path)
		lvg.run(opts.vgrad, opts.nofig, opts.no_intermediate)







