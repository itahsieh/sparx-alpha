##
## List of available geometries
## DO NOT move or rename these symbols unless you know absolutely what you're doing!
##
SPH1D = "sph1d"
REC3D = "rec3d"
SPH3D = "sph3d"
CYL3D = "cyl3d"
GEOM_DICT = {SPH1D: 0, REC3D: 1, SPH3D: 2, CYL3D: 3}

# Some necessary imports
from math import sqrt

class Grid(object):
	'''
	Grid class used for describing source model with Python data structures.
	This is used mainly in grid generation routines defined by the user.
	'''
	def __init__(self, geom, shape, min, max, init_grid=True):
		import numpy as np
		assert len(shape) == 3
		assert len(min) == 3
		assert len(max) == 3
		assert geom in GEOM_DICT

		# Declare attributes
		self.geom = GEOM_DICT[geom]
		self.shape = []
		self.min = []
		self.max = []
		self.cen = []
		self.delta = []

		# Sanity checks
		if geom == SPH1D:
			assert shape[0] > 0 and shape[1] == 1 and shape[2] == 1
		#elif geom == REC3D:
		else:
			assert shape[0] * shape[1] * shape[2] > 0

		# Grid shape and dimensions
		for i in range(3):
			assert max[i] > min[i] and min[i] == 0
			self.shape += [int(shape[i])]
			self.min += [float(min[i])]
			self.max += [float(max[i])]
			self.cen += [self.min[i] + 0.5 * (self.max[i] - self.min[i])]
			self.delta += [float(max[i]) - float(min[i])]

		# Global physical parameters
		self.gas_to_dust = 0.0
		self.T_cmb = 0.0

		# Zone dimensions
		self.min_i = np.zeros(shape=self.shape)
		self.min_j = np.zeros(shape=self.shape)
		self.min_k = np.zeros(shape=self.shape)

		self.max_i = np.zeros(shape=self.shape)
		self.max_j = np.zeros(shape=self.shape)
		self.max_k = np.zeros(shape=self.shape)

		self.cen_i = np.zeros(shape=self.shape)
		self.cen_j = np.zeros(shape=self.shape)
		self.cen_k = np.zeros(shape=self.shape)

		if init_grid:
			# Init zone dimensions (linear by default)
			self.Grid_linear()

		# Gas density
		self.n_H2 = np.zeros(shape=self.shape)

		# Molecular and collisional partner abundances
		self.X_mol = np.zeros(shape=self.shape)
		self.X_pH2 = np.zeros(shape=self.shape)
		self.X_oH2 = np.zeros(shape=self.shape)
		self.X_e = np.zeros(shape=self.shape)
		self.X_H = np.zeros(shape=self.shape)
		self.X_He = np.zeros(shape=self.shape)

		# Kinetic temperature
		self.T_k = np.zeros(shape=self.shape)
		# Dust temperature
		self.T_d = np.zeros(shape=self.shape)
		# Free-free brightness temperature
		self.T_ff = np.zeros(shape=self.shape)
		# Continuum brightness temperature
		self.T_bb = np.zeros(shape=self.shape)

		# Turbulent velocity
		self.V_t = np.zeros(shape=self.shape)

		# Gas velocity
		self.V_i = np.zeros(shape=self.shape)
		self.V_j = np.zeros(shape=self.shape)
		self.V_k = np.zeros(shape=self.shape)

		# Dust opacity
		self.kapp_d = np.zeros(shape=self.shape, dtype='|S64') # NumPy string array
		# Free-free opacity
		self.kapp_ff = np.zeros(shape=self.shape, dtype='|S64') # NumPy string array

		# debug
		#print "Grid init complete"

	def SetZoneBoundary(self, pos, min, max):
		self.min_i[pos] = min[0]
		self.min_j[pos] = min[1]
		self.min_k[pos] = min[2]
		self.max_i[pos] = max[0]
		self.max_j[pos] = max[1]
		self.max_k[pos] = max[2]
		self.cen_i[pos] = min[0] + 0.5 * (max[0] - min[0])
		self.cen_j[pos] = min[1] + 0.5 * (max[1] - min[1])
		self.cen_k[pos] = min[2] + 0.5 * (max[2] - min[2])

	def SetZoneBoundary_sph1d(self, pos, min, max):
		from sparx.physics import Const as C
		i = int(pos)
		self.SetZoneBoundary((pos, 0, 0), (min, 0, 0), (max, C.pi, 2.0*C.pi));

	def Grid_linear(self):
		import numpy as np
		delta = [self.delta[i] / float(self.shape[i]) for i in range(3)]
		for pos in np.ndindex(self.shape[0], self.shape[1], self.shape[2]):
			min = [delta[i] * pos[i] for i in range(3)]
			max = [delta[i] * (pos[i] + 1) for i in range(3)]
			self.SetZoneBoundary(pos, min, max)

	def Grid_sph1d_linear(self):
		import numpy as np
		delta = (self.max[0] - self.min[0]) / float(self.shape[0])
		for pos in np.ndindex(self.shape[0], 1, 1):
			i = pos[0]
			self.SetZoneBoundary_sph1d(i, delta * i, delta * (i + 1))

	def Grid_sph1d_log10(self, min):
		import numpy as np
		from math import log10
		assert self.max[0] > min
		logmax = log10(self.max[0])
		logmin = log10(min)
		delta = (logmax - logmin) / float(self.shape[0] - 1)
		for pos in np.ndindex(self.shape[0], 1, 1):
			i = pos[0]
			if i == 0:
				min = 0
			else:
				min = 10.0**(logmin + delta * (i - 1))
			max = 10.0**(logmin + delta * i)
			self.SetZoneBoundary_sph1d(i, min, max)

	def GetZoneCen(self, pos):
		return (self.cen_i[pos], self.cen_j[pos], self.cen_k[pos])

################################################################################

# Grid_sph1d subclass
class Grid_sph1d(Grid):
	def __init__(self, nshell, r_max, **kwargs):
		from sparx.physics import Const as C
		super(Grid_sph1d, self).__init__("sph1d", (nshell, 1, 1), (0, 0, 0), (r_max, C.pi, 2.0*C.pi), **kwargs)

################################################################################

# Grid_rec3d subclass
class Grid_rec3d(Grid):
	def __init__(self, shape, dims, **kwargs):
		super(Grid_rec3d, self).__init__("rec3d", shape, (0, 0, 0), dims, **kwargs)

################################################################################

class SPARXH5(object):
	"""
	SPARXH5: interface to sparx model files
	"""
	def __init__(self, fname):
		# Remember filename
		self.fname = fname

		# Open file
		from tables import openFile
		self.h5f = openFile(fname)

		# Load parameters for root zone (although there is only one
		# row in the ZONE table, it still has to be loaded this way)
		#self.root = [i for i in self.h5f.root.ZONE][0] --> deprecated in newer versions?
		self.root = self.h5f.root.ZONE[0]
		self.shape = list(self.root['NAXES'])
		self.naxis = len(self.shape)
		assert self.naxis == 3

		# Get coordinate system
		self.coord = self.root['geom']
		assert self.coord in GEOM_DICT

		# Get total number of zones (rows)
		self.nzone = reduce(lambda x,y: x*y, self.shape)

		# Get dimensions
		self.X_max = [i for i in self.root['X_max']]
		self.X_min = [i for i in self.root['X_min']]
		self.X_cen = [i for i in self.root['X_cen']]
		self.X_delta = [(self.X_max[i] - self.X_min[i]) / self.shape[i] for i in range(self.naxis)]

		# Get central position and number of radial data points
		self.center = [i // 2 for i in self.shape]
		self.nradial = self.shape[0] // 2

		# Get grid table
		self.grid_table = getattr(self.h5f.root, "GRID")

		# Get pops table if available
		if hasattr(self.h5f.root, "POPS"):
			self.pops_table = getattr(self.h5f.root, "POPS")
		else:
			self.pops_table = None

		# Get tau table if available
		if hasattr(self.h5f.root, "TAU"):
			self.tau_table = getattr(self.h5f.root, "TAU")
		else:
			self.tau_table = None
		return

	def GetParmData(self, parm, index=0):
		# Identify correct table to load
		from re import match
		m_lev = match("lev([0-9]+)", parm)
		m_tau = match("tau([0-9]+)", parm)
		if m_lev is not None:
			table = self.pops_table
		elif m_tau is not None:
			table = self.tau_table
			parm = "line"+m_tau.group(1)
		else:
			table = self.grid_table

		# Load table data
		try:
			data = [i[parm][index] for i in table]
		except TypeError:
			data = [i[parm] for i in table]

		return data

	def GetSlice(self, slice_axis, slice_index, parm, index=0):
		'''Get slice at (slice_axis, slice_index)'''
		from numpy import ndarray
		# Get data cube
		data = self.GetParmData(parm, index)

		# Get shape of slice and init slice array
		shape = self.shape
		shape.pop(slice_axis)
		slice = ndarray(shape)

		# Iterate over slice and load values
		for i in ndindex(self.shape[0], self.shape[1]):
			index = [j for j in i]
			index.insert(slice_axis, slice_index)
			slice[i] = data[index[2] + naxes[2] * (index[1] + naxes[1] * index[0])]
		return slice

	def GetPixel(self, pixel_axis, pixel_pos, parm, index=0):
		'''Get pixel at (pixel_axis, pixel_pos)'''
		from numpy import ndarray
		# Get data table
		data = self.GetParmData(parm, index)

		# Get shape of slice and init slice array
		shape = self.shape
		npix = shape[slice_axis]
		pixel = ndarray([npix])

		# Iterate over slice and load values
		for i in range(npix):
			index = [pixel_pos[0], pixel_pos[1]]
			index.insert(pixel_axis, i)
			pixel[i] = data[index[2] + naxes[2] * (index[1] + naxes[1] * index[0])]
		return pixel

	def GetRadii(self):
		"""
		Ger array of radial coordinates in [pc]
		"""
		from numpy import array
		if self.coord == SPH1D:
			return array([i['X_cen'][0] for i in self.grid_table])

		elif self.coord == REC3D:
			delta = self.X_delta[0]
			return array([delta * (i + 0.5) for i in range(self.nradial)])

	def GetRadial(self, parm, index=0):
		"""
		Get array of radial physical parameters
		"""
		import numpy as np
		# Get data cube
		data = self.GetParmData(parm, index)
		if self.coord == SPH1D:
			return np.array(data)
			
		elif self.coord == REC3D:
			# Get zone positions
			Xcen = self.GetParmData("X_cen", 0)
			Ycen = self.GetParmData("X_cen", 1)
			Zcen = self.GetParmData("X_cen", 2)
			n_H2 = self.GetParmData("n_H2")

			# Generate list of radii
			rmax = (self.X_max[0] - self.X_min[0]) / 2.0
			delta =  rmax / self.nradial
			radii = np.array([i * delta for i in range(self.nradial)])

			# Init radial data array
			radial = np.zeros(self.nradial)
			navg = np.zeros(self.nradial)

			# Get averaged data
			from numpy import ndindex
			for i in ndindex(self.shape[0], self.shape[1], self.shape[2]):
				# Get location in table for current grid position
				arrpos = i[2] + self.shape[2] * (i[1] + self.shape[1] * i[0])

				# Calculate distance from center
				dx = Xcen[arrpos] - self.X_cen[0]
				dy = Ycen[arrpos] - self.X_cen[1]
				dz = Zcen[arrpos] - self.X_cen[2]
				r = sqrt(dx**2 + dy**2 + dz**2)
				rpos = np.searchsorted(radii, r) - 1

				# If within bounds, search for position in radial array
				if rpos < self.nradial and n_H2[arrpos] > 0:
					radial[rpos] += data[arrpos]
					navg[rpos] += 1

			# Calculate average
			for i in range(self.nradial):
				if navg[i] > 0:
					radial[i] /= navg[i]
			return radial

	def Close(self):
		self.h5f.close()
		return

	def __del__(self):
		self.Close()
		return

