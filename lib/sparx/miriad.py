##
## Interface to miriad data
##

import ctypes as ct
import sparx._sparx as _sparx
import numpy as np
from sparx import physics as phys
from sparx.utils import clear_path
Unit = phys.Units
Cnst = phys.Const

################################################################################

class MirXYV:
	"""
	Interface for Miriad image datasets (must be in XYV order)
	"""
	class axis:
		"""
		Axis parameters
		"""
		def __init__(self, n, cdelt, crpix, crval):
			self.n = int(n)
			self.cdelt = float(cdelt)
			self.crpix = float(crpix)
			self.crval = float(crval)

	def __init__(self, fname):
		from os.path import exists
		if not exists(fname):
			raise Exception, "File '%s' does not exist"%fname
		self.fname = fname
		image = _sparx.load_mir_xyv2(fname)
		self.x = x = self.axis(image["xcsize"], image["xcdelt"], image["xcrpix"], image["xcrval"])
		self.y = y = self.axis(image["ycsize"], image["ycdelt"], image["ycrpix"], image["ycrval"])
		self.v = v = self.axis(image["vcsize"], image["vcdelt"], image["vcrpix"], image["vcrval"])
		self.x_list = np.array([(i - x.crpix) * x.cdelt + x.crval for i in range(x.n)])
		self.y_list = np.array([(i - y.crpix) * y.cdelt + y.crval for i in range(y.n)])
		self.v_list = np.array([(i - v.crpix) * v.cdelt + v.crval for i in range(v.n)])
		self.cube = image["cube"]
		self.bunit = image["bunit"]
		return

	def pix_to_offasec(self, ix, iy):
		from sparx.physics import Units as U
		xasec = (self.x_list[ix] - self.x_list[self.x.n // 2]) / U.asec
		yasec = (self.y_list[iy] - self.y_list[self.y.n // 2]) / U.asec
		return xasec, yasec

	def offasec_to_pix(self, xasec, yasec):
		from sparx.physics import Units as U
		# Locate nearest offset pixel
		x_offset = self.x_list[self.x.n // 2] + xasec * U.asec
		y_offset = self.y_list[self.y.n // 2] + yasec * U.asec
		if self.x.cdelt >= 0:
			x_list = self.x_list
			x_pix = x_list.searchsorted(x_offset)
		else:
			x_list = self.x_list[::-1]
			x_pix = self.x.n - x_list.searchsorted(x_offset, side='right')
		if self.y.cdelt >= 0:
			y_list = self.y_list
			y_pix = y_list.searchsorted(y_offset)
		else:
			y_list = self.y_list[::-1]
			y_pix = self.y.n - y_list.searchsorted(y_offset, side='right')
		return x_pix, y_pix

	def GetPix(self, ix, iy, iv):
		return self.cube[ix, iy, iv]

	def GetSpectrum(self, x, y):
		import numpy as np
		spec = np.zeros(self.v.n)
		for i in range(self.v.n):
			spec[i] = self.cube[x, y, i]
		return spec

	def GetSpecOffASec(self, x_offasec, y_offasec):
		x_pix, y_pix = self.offasec_to_pix(x_offasec, y_offasec)
		return self.GetSpectrum(x_pix, y_pix)

################################################################################

class MirTask:
	"""
	General Miriad command class
	"""
	def __init__(self, name, shell=True, **popen_args):
		self.name = name
		self.shell = shell
		self.popen_args = popen_args
		return

	def _buildcommand(self, **task_args):
		"""
		Build command line
		"""
		# Replace 'inp' keys with 'in'
		if "inp" in task_args.keys():
			task_args["in"] = task_args["inp"]
			del task_args["inp"]

		# Build command line
		return self.name + " " + " ".join(["%s='%s'"%(key, task_args[key]) for key in task_args])

	def __call__(self, **task_args):
		"""
		Launch command
		"""
		command = self._buildcommand(**task_args)

		from subprocess import call
		ret = call(command, shell=self.shell, **self.popen_args)

		# Check for errors
		if ret != 0:
			raise Exception, "Command '%s' failed (returncode=%s)"%(command, ret)

		return ret


	def Popen(self, **task_args):
		"""
		Popen the Miriad command
		"""
		from subprocess import Popen
		command = self._buildcommand(**task_args)
		return Popen(command, shell=self.shell, **self.popen_args)

################################################################################

def convol(map, out, fwhm):
	"""
	Convolve with Gaussian beam
	"""
	MirTask("convol")(map=map, out=out, fwhm=fwhm / Unit.asec)
	return

################################################################################

def subcont(inp, out, cont_chan=0):
	"""
	Subtract continuum from line maps
	"""
	# Clear temporary space for the continuum channel
	cont = "/tmp/cont"
	clear_path(cont)

	# Do the subtraction
	MirTask("imsub")(inp=inp, out=cont, region='image(%d)'%(cont_chan+1))
	MirTask("maths")(exp='<%s>-<%s>'%(inp, cont), out=out, options='grow')
	return

################################################################################

def subcont_convol(inp, out, fwhm):
	"""
	Subtract continuum and do convolution in one call
	"""
	# Clear temporary space for the continuum-subtracted map
	sub = "/tmp/sub"
	clear_path(sub)

	# Subtract continuum
	subcont(inp, sub)

	# Convolve with beam
	convol(sub, out, fwhm)
	return

################################################################################

def convert_flux_to_tb(inp, out, freq, fwhm):
	"""
	Convert flux to brightness temperature
	"""
	mir_freq = freq / 1.0e9 # [Hz] -> [GHz]
	mir_fwhm = fwhm / Unit.asec # [radian] -> [arcsec]
	MirTask("maths")(out=out, exp="1.22e6*<%s>/((%g**2)*(%.2f**2))"%(inp, mir_fwhm, mir_freq))
	MirTask("puthd")(inp=out+"/bunit", value="K")
	return

################################################################################

def get_rms(inp, region):
	"""Get RMS noise from line free channels using the 'histo' task"""
        p = MirTask("histo").Popen(inp=inp, region=region)
        lines = p.stdout.readlines()
        #print "".join(lines)
        # Line 3 contains statistical parameters
        if VERSION == "4.1.3":
                cols = lines[2].strip().split()
        elif VERSION == "4.0.4":
                cols = lines[2].strip().split()
        rms = float(cols[3])
        return rms
















