"""
Various physical models
"""

from sparx import physics, utils
U = physics.Units
C = physics.Const
import numpy as np
from math import exp

class Shu1977:
	'''
	Shu (1997) self-similar collapse model
	Given isothermal sound speed and collapse age, calculate
	mass density, gas velocity and number density distribution
	of a collapsing singular isothermal sphere given by Shu (1977).
	'''

	def __init__(self, a_sound, t_age):
		'''Initialize physical parameters and create interpolation objects'''
		from numpy import array
		from scipy.interpolate import interpolate as interp
		from math import log10

		# Init sound speed and collapse age
		self.a_sound = float(a_sound)
		self.t_age = float(t_age)

		# Dimensionless parameters given in table 2 of the paper.
		# These describe the "expansion-wave collapse solution" where A = 2+
		# (see section II-c for details)
		table2_x = array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1])
		table2_alpha = array([71.5, 27.8, 16.4, 11.5, 8.76, 7.09, 5.95, 5.14, 4.52, 4.04, 3.66, 3.35, 3.08, 2.86, 2.67, 2.5, 2.35, 2.22, 2.1, 2])
		table2_negv = array([5.44, 3.47, 2.58, 2.05, 1.68, 1.4, 1.18, 1.01, 0.861, 0.735, 0.625, 0.528, 0.442, 0.363, 0.291, 0.225, 0.163, 0.106, 0.051, 0])

		# Set interpolation limits
		self.xmin = table2_x[0]
		self.xmax = table2_x[-1]

		# Init interpolation objects:
		#   alpha is best interpolated linearly in log space
		log10_x = [log10(x) for x in table2_x]
		log10_alpha = [log10(alpha) for alpha in table2_alpha]
		self.alpha_interp = interp.interp1d(log10_x, log10_alpha, kind='linear')

		#   neg_v is best interpolated cubically
		self.negv_interp = interp.interp1d(table2_x, table2_negv, kind='cubic')

	def u(self, r):
		'''Given radius r, calculate gas veocity u'''
		# x = r / (a * t)
		# u(r, t) = a * v(x)
		x = r / (self.a_sound * self.t_age)

		if x < self.xmin:
			# Solution for small radii (eq 17)
			v = -(2.0 * 0.975 / x)**0.5
		elif x > self.xmax:
			# Solution for radii outside of the expansion wave
			# (i.e. static zone, eq 13)
			v = 0
		else:
			# Interpolate from table 2
			v = -(self.negv_interp(x))
		return self.a_sound * v

	def rho(self, r):
		'''Given radius r, calculate mass density rho'''
		# x = r / (a * t)
		# rho(r, t) = alpha(x) / (4 * pi * G * t**2)
		x = r / (self.a_sound * self.t_age)

		if x < self.xmin:
			# Solution for small radii (eq 17)
			alpha = (0.975 / (2.0 * x**3.0))**0.5
		elif x > self.xmax:
			# Solution for radii outside of the expansion wave
			# (i.e. static zone, eq 13)
			alpha = 2.0 / (x**2.0)
		else:
			# Interpolate from table 2
			from math import log10
			alpha = 10.0**(self.alpha_interp(log10(x)))
		return alpha / (4.0 * CONST.pi * CONST.G * self.t_age**2.0)

	def n(self, r):
		return self.rho(r) / (2.8 * CONST.amu)


################################################################################

class Myers1996:
	"""
	The Myers (1996) two-layer model for spectra from infalling cores.
	This is only valid for the Shu (1977) model case where p=-3/2 and q=-1/2.
	"""
	def __init__(self, mol, iline, T_b, tau0, T_k, V_in, sigma_nt, trapping=False):
		from math import exp
		##
		## Free parameters set by user
		##
		self.mol = mol # Molecule
		self.iline = iline = int(iline) # Line index
		self.T_b = T_b = float(T_b) # Brightness temperature of cosmic background
		self.tau0 = tau0 = float(tau0) # Peak optical depth
		self.T_k = T_k = float(T_k) # Kinetic temperature (K)
		self.V_in = V_in = float(V_in) # Gas velocity at optically thin limit (m/s)
		self.sigma_nt = sigma_nt = float(sigma_nt) # Non-thermal line width (m/s)
		self.sigma = (sigma_nt**2.0 + (C.k * T_k / mol.mass))**0.5

		##
		## Calculate derived parameters (see sec. 2 of Myers 1996)
		##
		# <Q> = attenuation-weighted mean of Q
		# Velocities
		# V_in = 1/2 * V_max --> V_max = 2 * V_in
		self.V_max = V_max = 2.0 * V_in

		# <Vf> and <Vr> are mean gas velocities at the front and rear of the cloud
		self.Vf_mean = V_max * (1.0 - (1.0 + tau0) * exp(-tau0)) / (tau0 * (1.0 - exp(-tau0)))
		self.Vr_mean = V_max * (tau0 - (1.0 - exp(-tau0))) / (tau0 * (1.0 - exp(-tau0)))

		# Densities
		# n_cr = n_thin = 1/4 * n_max --> n_max = 4 * n_cr
		self.n_cr = mol.col[0].get_crit_dens(iline, T_k)
		self.n_max = n_max = 4.0 * self.n_cr

		# <nf> and <nr> are mean gas densities at the front and rear of the cloud
		self.nf_mean = nf_mean = n_max * (1.0 - exp(-tau0))**-1.0 * ((6.0 / tau0**3.0) - exp(-tau0) * (1.0 + 3.0 / tau0 + 6.0 / tau0**2.0 + 6.0 / tau0**3.0))
		self.nr_mean = nr_mean = n_max * (1.0 - exp(-tau0))**-1.0 * (1.0 - 3.0 / tau0 + 6.0 / tau0**2.0 - (6.0 / tau0**3.0) * (1.0 - exp(-tau0)))

		# To = h * nu / k
		self.To = To = C.h * mol.line_freq[iline] / C.k

		# Excitation temperatures at the front and rear of the cloud
		if trapping:
			beta = (1.0 - exp(-tau0)) / tau0
		else:
			beta = 1.0
		self.Tf_mean = T_k * (T_b + (4.0 * To / beta) * nf_mean / n_max) / (T_k + (4.0 * To / beta) * nf_mean / n_max)
		self.Tr_mean = T_k * (T_b + (4.0 * To / beta) * nr_mean / n_max) / (T_k + (4.0 * To / beta) * nr_mean / n_max)
		return

	def calc_tau_f(self, velo):
		"""
		Calculate the optical depth of the forward cloud
		"""
		from math import exp
		return self.tau0 * exp(-(velo - self.Vf_mean)**2.0 / (2.0 * self.sigma**2.0))

	def calc_tau_r(self, velo):
		"""
		Calculate the optical depth of the rear cloud
		"""
		from math import exp
		return self.tau0 * exp(-(velo + self.Vr_mean)**2.0 / (2.0 * self.sigma**2.0))

	def J(self, Tex):
		"""
		J(T) = (h*nu/k) / (exp(h*nu/(k*T))-1)
		"""
		from math import exp
		return self.To / (exp(self.To / Tex) - 1.0)

	def delta_TB(self, tau_f, tau_r):
		"""
		Brightness temperature of the spectral line, obtained from the equation
		of radiative transfer
		"""
		from math import exp
		return self.J(self.Tf_mean) * (1.0 - exp(-tau_f)) +\
		       self.J(self.Tr_mean) * (1.0 - exp(-tau_r)) * exp(-tau_f) -\
		       self.J(self.T_b) * (1.0 - exp(-tau_r - tau_f))

	def calc_spectrum(self, bandwidth, nchan):
		"""
		Given bandwidth and number of channels, calculate spectrum
		"""
		chan0 = nchan // 2
		delta_v = float(bandwidth) / nchan
		v_arr = np.array([(i - chan0) * delta_v for i in range(nchan)])
		return np.array([self.delta_TB(self.calc_tau_f(i), self.calc_tau_r(i)) for i in v_arr]), v_arr

################################################################################

class TwoLayer:
	"""
	The 'two-layer' model for infall spectra described in de Vries (2005)
	"""
	def __init__(self, tau0, sigma, Tf, Tr, vlsr, vin, mol, iline, Tb):
		# Free parameters
		self.tau0 = float(tau0) # Dimensionless
		self.sigma = float(sigma) # [m/s]
		self.Tf = float(Tf) # [K]
		self.Tr = float(Tr) # [K]
		self.vlsr = float(vlsr) # [m/s]
		self.vin = float(vin) # [m/s]

		# Molecule and background radiation
		self.mol = mol
		self.iline = int(iline)
		self.Tb = float(Tb)

		# T0 = h * nu / k
		self.T0 = C.h * mol.line_freq[iline] / C.k
		return

	def get_J(self, Tex):
		"""
		J(T) = (h*nu/k) / (exp(h*nu/(k*T))-1)
		"""
		try:
			x = abs(self.T0 / Tex)
			# exp overflows at around 700, this must be prevented
			if abs(Tex) > 1.0e-3 and abs(x) > 1.0e-4 and abs(x) <= 700.0:
				return self.T0 / (exp(self.T0 / Tex) - 1.0)
			else:
				return 0.0
		except:
			print "T0=", self.T0
			print "Tex=", Tex
			print "T0/Tex=", self.T0 / Tex
			raise
		

	def get_tauf(self, velo):
		return self.tau0 * exp(-((velo - self.vlsr - self.vin) / (2.0 * self.sigma))**2.0)

	def get_taur(self, velo):
		return self.tau0 * exp(-((velo - self.vlsr + self.vin) / (2.0 * self.sigma))**2.0)

	def get_deltaTB(self, velo):
		Jf = self.get_J(self.Tf)
		Jr = self.get_J(self.Tr)
		Jb = self.get_J(self.Tb)
		tauf = self.get_tauf(velo)
		taur = self.get_taur(velo)
		return Jf * (1.0 - exp(-tauf)) + Jr * (1.0 - exp(-taur)) * exp(-tauf) - Jb * (1.0 - exp(-taur - tauf))

################################################################################

class Hill(TwoLayer):
	"""
	The 'hill' model for infall spectra described in de Vries (2005)
	"""
	def __init__(self, tauC, sigma, To, Tp, vlsr, vC, mol, iline, Tb):
		# Free parameters
		self.tauC = float(tauC) # Dimensionless
		self.sigma = float(sigma) # [m/s]
		self.To = float(To) # [K]
		self.Tp = float(Tp) # [K]
		self.vlsr = float(vlsr) # [m/s]
		self.vC = float(vC) # [m/s]

		# Molecule and background radiation
		self.mol = mol
		self.iline = int(iline)
		self.Tb = float(Tb)

		# T0 = h * nu / k
		self.T0 = C.h * mol.line_freq[iline] / C.k
		return

	def get_tauf(self, velo):
		return self.tauC * exp(-((velo - self.vlsr - self.vC) / (2.0 * self.sigma))**2.0)

	def get_taur(self, velo):
		return self.tauC * exp(-((velo - self.vlsr + self.vC) / (2.0 * self.sigma))**2.0)

	def get_TB(self, J1, J2, tau0):
		"""
		Solution to the equation of radiative transfer assuming a linear J.
		Equations (6) and (9) are either WRONG or cannot be implemented numerically; below
		is the result of my own re-derivation.

		Since tau0 is at the denominator, a lower limit must be set so the result
		does not diverge for small tau0, otherwise poor fits can easily occur.
		"""
		if abs(tau0) > 1e-3:
			return J1 * (1.0 - exp(-tau0)) + ((J2 - J1) / tau0) * (1.0 - (tau0 + 1.0) * exp(-tau0))
		else:
			return 0.0

	def get_deltaTB(self, velo):
		Jp = self.get_J(self.Tp)
		Jo = self.get_J(self.To)
		Jb = self.get_J(self.Tb)
		tauf = self.get_tauf(velo)
		taur = self.get_taur(velo)
		TBf = self.get_TB(Jo, Jp, tauf)
		TBr = self.get_TB(Jp, Jo, taur)
		return TBf + TBr * exp(-tauf) + Jb * exp(-taur-tauf) - Jb

	def get_spectrum(self, vmin, vmax, n=100):
		velo_arr = utils.generate_linear_points(vmin, vmax, n)
		return np.array([self.get_deltaTB(velo) for velo in velo_arr])

	def get_spectrum_varr(self, velo_arr):
		return np.array([self.get_deltaTB(velo) for velo in velo_arr])

################################################################################

class Hill5(Hill):
	"""
	The 'hill' model for infall spectra described in de Vries (2005)
	This special case has 5 free parameters: tauC, vlsr, vC, sigma, Tp
	"""
	def __init__(self, tauC, vlsr, vC, sigma, Tp, Tb, mol, iline):
		To = Tb
		Hill.__init__(self, tauC, sigma, To, Tp, vlsr, vC, mol, iline, Tb)
		return

################################################################################

class Hill6b(Hill):
	"""
	The 'hill' model for infall spectra described in de Vries (2005)
	This special case has 6 free parameters: tauC, vlsr, vC, sigma, Tp, To
	"""
	def __init__(self, tauC, vlsr, vC, sigma, Tp, To, Tb, mol, iline):
		Hill.__init__(self, tauC, sigma, To, Tp, vlsr, vC, mol, iline, Tb)
		return





