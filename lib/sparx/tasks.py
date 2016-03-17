##
## This module handles most things related to sparx tasks
##
# SOP for adding new tasks:
# 1) Make a copy of 'class Task_Template' and the corresponding 'install_task()' call in this file.
# 2) Change the names 'Task_Template' and 'task_template' to 'Task_Something' and 'task_something'
#    where 'Something' and 'something' is the name of your task.
# 3) Add the proper keywords to self.keys.
# 4) Modify the 'self.configure' and 'self.main' methods to your liking.
#    self.configure() is run during self.__init__(). When the task is called,
#    the execution order is: slef._proc_inps() -> self.main() -> self.cfunc()
# 5) Point self.cfunc to a proper C function if needed.
# 6) All done!

# Some necessary imports
from math import sqrt
from sparx.utils import MESG as Mesg, MPI_SIZE, MPI_RANK
from sparx.inputs import Key, Type, INP_DICT
from sparx import _sparx
from sparx import physics as phys
import numpy as np

Unit = phys.Units
Cnst = phys.Const

##
## Global task dictionary
##
TASK_DICT = {}

##
## Task class
##
class Task(object):
	##
	## Object constructor
	##
	def __init__(self, name, desc="", expl="", cfunc=None, keys=[], parallel=False):
		# Name
		self.name = str(name)

		# Configuration overridden by child
		if hasattr(self, "configure"):
			self.configure()

		# Description
		if not hasattr(self, "desc"):
			self.desc = str(desc)

		# Explanation
		if not hasattr(self, "expl"):
			self.expl = str(expl)

		# C function to call
		if not hasattr(self, "cfunc"):
			self.cfunc = cfunc

		# Keywords for this task
		if not hasattr(self, "keys"):
			assert type(keys) is list
			self.keys = keys

		# Store keys in dictionary for easy lookup
		self.keydict = {}
		for key in self.keys:
			self.keydict[key.name] = key

		# Is this a parallelized task?
		if not hasattr(self, "parallel"):
			self.parallel = parallel

		# Generate docstring
		dochead = "Name:\n"+\
		"  %s\n"%self.name+\
		"\n"+\
		"Description:\n"+\
		"  %s\n"%self.desc

		dockey = "Keywords:\n"
		for key in self.keys:
			dockey += "%12s=%s\n"%(key.name, repr(key.typ))
			if key.deflt is None:
				dockey += " "*13+"(No default)\n"
			elif key.deflt is Type.Optional:
				dockey += " "*13+"(Optional)\n"
			else:
				dockey += " "*13+"(default: \"%s\")\n"%key.deflt
			dockey += " "*13+key.desc+"\n\n"

		doctail = \
		"Explanation:\n"+\
		"  %s"%self.expl
		self.__doc__ = "\n".join([dochead, dockey, doctail])

		return

	##
	## Inputs processor
	##
	def _proc_inps(self, **kwargs):
		# args is a dictionary of user-specified keyword=value entries,
		# so convert them to the proper values and store in the
		# sparx.inputs.INP_DICT dictionary for later use by the task
		for kwrd in kwargs:
			# Retrieve value
			valu = kwargs[kwrd]

			# Convert value to proper format if kwrd is a valid key
			if kwrd in self.keydict:
				key = self.keydict[kwrd]
				try:
					# Convert value string
					INP_DICT[kwrd] = key.typ(valu)
				except:
					# Something went wrong with the conversion
					raise
			else:
				raise Exception, "Key '%s' is not a valid keyword for task '%s'" % (kwrd, self.name)

		# Make sure input is valid by checking
		# whether all required keys have been set
		for key in self.keys:
			if key.name not in INP_DICT:
				if key.deflt is None:
					# This key is mandatory
					raise Exception, "Key '%s' must be given for task '%s'"%(key.name, self.name)
				elif key.deflt is Type.Optional:
					# Set key to None (VERY IMPORTANT!)
					INP_DICT[key.name] = None
				else:
					# This key has a default value
					INP_DICT[key.name] = key.typ(key.deflt)
					Mesg("Key '%s' not given, using default=\"%s\"" % (key.name, key.deflt))
		return

	##
	## What to do when the task is executed
	##
	def __call__(self, **kwargs):
		# If this is not a parallel task, do not allow non-root processes to execute
		# the task. THIS IS VERY IMPORTANT! Failing to check this may result in scary
		# race conditions!!!
		if not self.parallel and MPI_RANK != 0:
			return

		# Process user inputs
		self._proc_inps(**kwargs)

		# Run the main() hook function if defined
		if hasattr(self, "main"):
			self.main()

		# Run C function if it isn't None. This function must
		# be able to raise Python exceptions since no error
		# checking is done after execution
		if self.cfunc is not None:
			self.cfunc()
		return

	##
	## What to do when called from the command line
	##
	def run(self, args):
		# args is a list of key=val pairs, so break them up
		# and turn them into a dictionary
		kwargs = {}
		for arg in args:
			toks = arg.split("=")
			if len(toks) < 2:
				raise Exception, "'%s' must be in the 'keyword=value' format"%arg
			kwargs[toks[0]] = "=".join(toks[1:])
		self.__call__(**kwargs)
		return

##
## Function to install task in this module and add entry to
## TASK_DICT dictionary
##
def install_task(task):
	import sys
	mod = sys.modules[__name__]

	if not hasattr(mod, task.name):
		setattr(mod, task.name, task)
		TASK_DICT[task.name] = task
	else:
		raise Exception, "Task '%s' already exists" % task.name

##
## Task definitions
##

class Task_Template(Task):
	"""
	Task template
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "A brief description of what the task does"

		# Explanation
		self.expl = "Some documentation on the task"

		# Keys
		self.keys = [
			Key("pos_int", Type.PosInt, 10, "Keyword for testing purposes only"),
			Key("int", Type.Integer, -10, "Keyword for testing purposes only"),
			Key("angle", Type.Angle, '0.1asec', "Keyword for testing purposes only"),
			Key("velo", Type.Velo, '100kms^-1', "Keyword for testing purposes only"),
			Key("length", Type.Length, '10pc', "Keyword for testing purposes only"),
		]

		# C function to call
		self.cfunc = _sparx.task_template

	##
	## Task procedures
	##
	def main(self):
		return

install_task(Task_Template("task_template"))

################################################################################


class Task_input1D(Task):
	"""
	Spherical AGB envelope, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "1-D AGB envelope"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("in", Type.OldFile, None, "Input model file"),
			Key("Xmol", Type.Fraction, None, "Molecular abundance"),
			Key("Vt", Type.Velo, None, "turbulent speed"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		R_min=[]
		R_cen=[]
		R_max=[]
		n_H2=[]
		T_k=[]
		V_i=[]
		filename = INP_DICT["in"]
		finput = open(filename,'r')
		for line in finput.readline():
			columns = line.split()
			R_min = columns[0]
			R_cen = columns[1]
			R_max = columns[2]
			n_H2 = columns[3]
			T_k = columns[4]
			V_i = columns[5]
	
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		radius_out = INP_DICT["R_out"] # [m] outer radius
		radius_in = 0.0                # [m] inner radius
		r0 = INP_DICT["R_in"]          # [m] wind radius
		n = INP_DICT["ndiv"]
		dr =(radius_out-radius_in)/n   # [m]
		Vwind=INP_DICT["velocity"]     # [m/s]
		
		from math import pi
		from numpy import zeros
		manual_grid=zeros(n+1)
		for pos in range(0,n+1):
			manual_grid[pos] = (radius_in+pos*dr) / 3.086e16  # [pc]

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius_out / 3.086e16) # [pc]

		# Setup model parameters
		grid.gas_to_dust = 100.
		grid.T_cmb = 2.73 # K
		mH2 = 2.*1.660538921e-27 # kg
		X0=INP_DICT["abundance"]
		MassLoss = (INP_DICT["massloss"]/mH2)/31536000.  # (s^-1)
		T0 = INP_DICT["temperature"] 
		P_temp = INP_DICT["P_temp"]
		P_abund = INP_DICT["P_abund"]
		Vt = INP_DICT["Vt"]
			
		for pos in np.ndindex(n, 1, 1):
			# Manually regrid if requested
			id = pos[0]
			if id == 0:
				grid.SetZoneBoundary_sph1d(id, manual_grid[id], manual_grid[id+1])
				assert manual_grid[id] < manual_grid[id+1]
			# Set physical parameters
			
			R = grid.cen_i[pos] * 3.086e16 # [m]
			if ( R < r0 ):
				factor=1.0
			else:
				factor=R/r0
			grid.n_H2[pos] = MassLoss/(4.0*pi* Vwind *(r0*r0)) *factor**(-2) # m^-3
			grid.T_k[pos] = T0 * factor**(-P_temp) # K
			grid.X_mol[pos] = X0 * factor**(-P_abund) # Fraction
			grid.V_i[pos] = Vwind  # m/s
			grid.V_t[pos] = Vt # [m/s]
			grid.X_pH2[pos] = 0.25
			grid.X_oH2[pos] = 0.75
			grid.T_d[pos] = grid.T_k[pos]
			grid.kapp_d[pos] = "table,jena_thin_e5"
		return

install_task(Task_input1D("task_input1d"))


################################################################################
class Task_AFGL3068_Woods_model(Task):
	"""
	Spherical AGB envelope, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "AFGL3068 1-D model used in Woods et. al. (2003)"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("ndiv", Type.PosInt, 64, "Number of shells"),
			Key("abundance", Type.Fraction, '3e-5', "Molecular abundance at inner boundary"),
			Key("Vt", Type.Velo, '800ms^-1', "turbulent speed"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		radius_out = 0.2 # [pc] outer radius
		radius_in = 0.0                # [pc] inner radius
		r0 = 100./206264.806          # [pc] wind radius
		n = INP_DICT["ndiv"]
		dr =(radius_out-radius_in)/n   # [m]
		Vwind=14000.     # [m/s]
		
		from math import pi
		from numpy import zeros
		manual_grid=zeros(n+1)
		for pos in range(0,n+1):
			manual_grid[pos] = (radius_in+pos*dr)  # [pc]

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius_out) # [pc]

		# Setup model parameters
		grid.gas_to_dust = 100.
		grid.T_cmb = 2.73 # [K]
		Msun2MKS=1.98892e30 # [kg]
		mH2 = 2.*1.660538921e-27 # [kg]
		X0=INP_DICT["abundance"]
		MassLoss = (2e-5*Msun2MKS/mH2)/31536000.  # [s^-1]
		Vt = INP_DICT["Vt"]
			
		for pos in np.ndindex(n, 1, 1):
			# Manually regrid if requested
			id = pos[0]
			if id == 0:
				grid.SetZoneBoundary_sph1d(id, manual_grid[id], manual_grid[id+1])
				assert manual_grid[id] < manual_grid[id+1]
			# Set physical parameters
			
			R = grid.cen_i[pos] #[pc]
			grid.n_H2[pos] = MassLoss/(4.0*pi* Vwind ) * (R*3.08567758e16)**(-2) # [m^-3]
			grid.T_k[pos] = 300.*(R/r0)**(-1) + 3.434e-7*(R/r0)**3 + 5.4623 # [K]
			grid.X_mol[pos] = INP_DICT["abundance"] # Fraction
			grid.V_i[pos] = Vwind  # [m/s]
			grid.V_t[pos] = Vt # [m/s]
			grid.X_pH2[pos] = 0.25
			grid.X_oH2[pos] = 0.75
			grid.T_d[pos] = grid.T_k[pos]
			grid.kapp_d[pos] = "table,jena_thin_e5"
		return

install_task(Task_AFGL3068_Woods_model("task_afgl3068_woods_model"))

################################################################################


class Task_AGB1D(Task):
	"""
	Spherical AGB envelope, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "1-D AGB envelope"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("massloss", Type.Mass, '1e-5Msun', "mass loss rate per year"),
			Key("velocity", Type.Velo, '10kms^-1', "wind velocity"),
			Key("ndiv", Type.PosInt, 64, "Number of shells"),
			Key("abundance", Type.Fraction, '3e-5', "Molecular abundance at inner boundary"),
			Key("R_in", Type.Length, '1au', "radius of inner boundary"),
			Key("R_out", Type.Length, '0.1pc', "radius of outer boundary"),
			Key("temperature", Type.Temp, '1000K', "Temperature at inner boundary"),
			Key("P_temp", Type.Float, '0.0', "power law for temperature profile"),
			Key("P_abund", Type.Float, '0.0', "power law for abundance profile"),
			Key("Vt", Type.Velo, '800ms^-1', "turbulent speed"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		radius_out = INP_DICT["R_out"] # [m] outer radius
		radius_in = 0.0                # [m] inner radius
		r0 = INP_DICT["R_in"]          # [m] wind radius
		n = INP_DICT["ndiv"]
		dr =(radius_out-radius_in)/n   # [m]
		Vwind=INP_DICT["velocity"]     # [m/s]
		
		from math import pi
		from numpy import zeros
		manual_grid=zeros(n+1)
		for pos in range(0,n+1):
			manual_grid[pos] = (radius_in+pos*dr) / 3.086e16  # [pc]

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius_out / 3.086e16) # [pc]

		# Setup model parameters
		grid.gas_to_dust = 100.
		grid.T_cmb = 2.73 # K
		mH2 = 2.*1.660538921e-27 # kg
		X0=INP_DICT["abundance"]
		MassLoss = (INP_DICT["massloss"]/mH2)/31536000.  # (s^-1)
		T0 = INP_DICT["temperature"] 
		P_temp = INP_DICT["P_temp"]
		P_abund = INP_DICT["P_abund"]
		Vt = INP_DICT["Vt"]
			
		for pos in np.ndindex(n, 1, 1):
			# Manually regrid if requested
			id = pos[0]
			if id == 0:
				grid.SetZoneBoundary_sph1d(id, manual_grid[id], manual_grid[id+1])
				assert manual_grid[id] < manual_grid[id+1]
			# Set physical parameters
			
			R = grid.cen_i[pos] * 3.086e16 # [m]
			if ( R < r0 ):
				factor=1.0
			else:
				factor=R/r0
			grid.n_H2[pos] = MassLoss/(4.0*pi* Vwind *(r0*r0)) *factor**(-2) # m^-3
			grid.T_k[pos] = T0 * factor**(-P_temp) # K
			grid.X_mol[pos] = X0 * factor**(-P_abund) # Fraction
			grid.V_i[pos] = Vwind  # m/s
			grid.V_t[pos] = Vt # [m/s]
			grid.X_pH2[pos] = 0.25
			grid.X_oH2[pos] = 0.75
			grid.T_d[pos] = grid.T_k[pos]
			grid.kapp_d[pos] = "table,jena_thin_e5"
		return

install_task(Task_AGB1D("task_agb1d"))

################################################################################
class Task_AGB1D_Xdissociated(Task):
	"""
	Spherical AGB envelope, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "1-D AGB envelope"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("massloss", Type.Mass, '1e-5Msun', "mass loss rate per year"),
			Key("velocity", Type.Velo, '10kms^-1', "wind velocity"),
			Key("ndiv", Type.PosInt, 64, "Number of shells"),
			Key("abundance", Type.Fraction, '3e-5', "Molecular abundance at inner boundary"),
			Key("rp", Type.Length, '1.9e14cm', "radius of molecular photodissociation"),
			Key("R_in", Type.Length, '1au', "radius of inner boundary"),
			Key("R_out", Type.Length, '0.1pc', "radius of outer boundary"),
			Key("temperature", Type.Temp, '1000K', "Temperature at inner boundary"),
			Key("P_temp", Type.Float, '0.0', "power law for temperature profile"),
			Key("P_abund", Type.Float, '2.5', "power law for abundance profile"),
			Key("Vt", Type.Velo, '800ms^-1', "turbulent speed"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		radius_out = INP_DICT["R_out"] # [m] outer radius
		radius_in = 0.0                # [m] inner radius
		r0 = INP_DICT["R_in"]          # [m] wind radius
		rp = INP_DICT["rp"]            # [m] radius of molecular photodissociation
		n = INP_DICT["ndiv"]
		dr =(radius_out-radius_in)/n   # [m]
		Vwind=INP_DICT["velocity"]     # [m/s]
		
		from math import pi,exp,log
		from numpy import zeros
		manual_grid=zeros(n+1)
		for pos in range(0,n+1):
			manual_grid[pos] = (radius_in+pos*dr) / 3.086e16  # [pc]

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius_out / 3.086e16) # [pc]

		# Setup model parameters
		grid.gas_to_dust = 100.
		grid.T_cmb = 2.73 # K
		mH2 = 2.*1.660538921e-27 # kg
		X0=INP_DICT["abundance"]
		MassLoss = (INP_DICT["massloss"]/mH2)/31536000.  # (s^-1)
		T0 = INP_DICT["temperature"] 
		P_temp = INP_DICT["P_temp"]
		P_abund = INP_DICT["P_abund"]
		Vt = INP_DICT["Vt"]
			
		for pos in np.ndindex(n, 1, 1):
			# Manually regrid if requested
			id = pos[0]
			if id == 0:
				grid.SetZoneBoundary_sph1d(id, manual_grid[id], manual_grid[id+1])
				assert manual_grid[id] < manual_grid[id+1]
			# Set physical parameters
			
			R = grid.cen_i[pos] * 3.086e16 #[m]
			if ( R < r0 ):
				factor=1.0
			else:
				factor=R/r0
			grid.n_H2[pos] = MassLoss/(4.0*pi* Vwind *(r0*r0)) *factor**(-2) # m^-3
			grid.T_k[pos] = T0 * factor**(-P_temp) # K
			grid.X_mol[pos] = X0 * exp( log(2.0)*(R/rp)**(P_abund) ) # Fraction
			grid.V_i[pos] = Vwind  # m/s
			grid.V_t[pos] = Vt # [m/s]
			grid.X_pH2[pos] = 0.25
			grid.X_oH2[pos] = 0.75
			grid.T_d[pos] = grid.T_k[pos]
			grid.kapp_d[pos] = "table,jena_thin_e5"
		return

install_task(Task_AGB1D_Xdissociated("task_agb1d_Xdissociated"))



class Task_Empty(Task):
	"""
	An empty task mainly for generating valgrind suppressions
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "A task that does nothing"

		# Explanation
		self.expl = "Not much explanation here"

		# Keys
		self.keys = [
			Key("pos_int", Type.PosInt, 10, "Keyword for testing purposes only"),
			Key("int", Type.Integer, -10, "Keyword for testing purposes only"),
			Key("angle", Type.Angle, '0.1asec', "Keyword for testing purposes only"),
			Key("velo", Type.Velo, '100kms^-1', "Keyword for testing purposes only"),
			Key("length", Type.Length, '10pc', "Keyword for testing purposes only"),
		]

		# C function to call
		self.cfunc = None

install_task(Task_Empty("task_empty"))

################################################################################

class Task_ValDust1D(Task):
	"""
	Validate dust radiative transfer, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Generate sph1d grid for testing LTE dust radiative xfer"

		# Explanation
		self.expl = "Generate a uniform sph1d grid full of dust for verifying dust radiative transfer"

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Name of output file (HDF5 file)"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions
		n = 100
		length = 0.1 # pc

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, length)

		# Set model parameters
		grid.gas_to_dust = 100.0
		for i in range(n):
			pos = i,0,0
			grid.n_H2[pos] = 1e4 * 1e6 # m^-3
			#grid.X_mol[pos] = 1e-4
			grid.T_k[pos] = 100.0 # K
			#grid.T_ff[pos] = grid.T_k[pos]
			#grid.kapp_ff[pos] = Type.KappFLaw('100GHz', '1.0e-2cm^2g^-1', 0)
			grid.T_d[pos] = grid.T_k[pos]
			grid.kapp_d[pos] = Type.KappLLaw("['1mm', '0.1m^2kg^-1', 0]")
		return

install_task(Task_ValDust1D("task_valdust1d"))

################################################################################

class Task_ValDust3D(Task):
	"""
	Validate dust radiative transfer, 3D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Generate rec3d grid for testing LTE dust radiative xfer"

		# Explanation
		self.expl = "Generate a uniform rec3d grid full of dust for verifying dust radiative transfer"

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Name of output file (HDF5 file)"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_rec3d
		# Setup dimensions
		radius = 0.1 # [pc] cloud radius
		length = radius * 2 # [pc]
		n = 32

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_rec3d((n, n, n), (length, length, length))

		# Set model parameters
		grid.gas_to_dust = 100.0
		grid.T_cmb = 0
		for pos in np.ndindex(n, n, n):
			# Calculate radius at pos
			di = grid.cen_i[pos] - grid.cen[0]
			dj = grid.cen_j[pos] - grid.cen[1]
			dk = grid.cen_k[pos] - grid.cen[2]
			r = sqrt(di**2 + dj**2 + dk**2)

			# Set parameters only if within cloud radius
			if r <= radius:
				# Set physical parameters
				grid.n_H2[pos] = 1e4 * 1e6 # m^-3
				grid.T_k[pos] = 100.0 # K
				grid.T_d[pos] = grid.T_k[pos]
				grid.kapp_d[pos] = Type.KappLLaw("['1mm', '0.1m^2kg^-1', 0]")
		return

install_task(Task_ValDust3D("task_valdust3d"))

################################################################################

class Task_ValLine1D(Task):
	"""
	Validate line radiative transfer, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Generate sph1d grid for testing LTE molecular line radiative xfer"

		# Explanation
		self.expl = "Generate a uniform sph1d grid full of molecular gas for verifying line radiative transfer"

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Name of output file (HDF5 file)"),
			Key("molec", Type.Molec, None, "Molecule to generate LTE populations for"),
			Key("xmol", Type.Fraction, None, "Molecular abundance"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions
		n = 100
		length = 0.1 # pc

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, length)

		# Set model parameters
		grid.gas_to_dust = 0
		for i in range(n):
			pos = i,0,0
			grid.n_H2[pos] = 1e9 * 1e6 # m^-3
			grid.T_k[pos] = 40.0 # K
			grid.X_mol[pos] = INP_DICT["xmol"] # fraction
		return

install_task(Task_ValLine1D("task_valline1d"))

################################################################################

class Task_ValLine3D(Task):
	"""
	Validate line radiative transfer, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Generate sph1d grid for testing LTE molecular line radiative xfer"

		# Explanation
		self.expl = "Generate a uniform sph1d grid full of molecular gas for verifying line radiative transfer"

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Name of output file (HDF5 file)"),
			Key("molec", Type.Molec, None, "Molecule to generate LTE populations for"),
			Key("xmol", Type.Fraction, None, "Molecular abundance"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_rec3d
		# Setup dimensions
		radius = 0.1 # [pc] cloud radius
		length = radius * 2 # [pc]
		n = 32

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_rec3d((n, n, n), (length, length, length))

		# Set model parameters
		grid.gas_to_dust = 100.0
		grid.T_cmb = 0
		for pos in np.ndindex(n, n, n):
			# Calculate radius at pos
			di = grid.cen_i[pos] - grid.cen[0]
			dj = grid.cen_j[pos] - grid.cen[1]
			dk = grid.cen_k[pos] - grid.cen[2]
			r = sqrt(di**2 + dj**2 + dk**2)

			# Set parameters only if within cloud radius
			if r <= radius:
				# Set physical parameters
				grid.n_H2[pos] = 1e9 * 1e6 # m^-3
				grid.T_k[pos] = 40.0 # K
				grid.X_mol[pos] = INP_DICT["xmol"] # fraction
		return

install_task(Task_ValLine3D("task_valline3d"))

################################################################################
class Task_P2a(Task):
	"""
	Leiden benchmark problem, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "2002 Benchmark Problem 2a"

		# Explanation
		self.expl = "Some documentation on the task"

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file")
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		n = 50
		manual_grid = [3.242862E-03, 3.506789E-03, 3.792235E-03, 4.100885E-03, 4.434650E-03, 4.795573E-03, 
                 5.185890E-03, 5.607998E-03, 6.064426E-03, 6.558026E-03, 7.091778E-03, 7.668989E-03, 8.293159E-03, 
                 8.968143E-03, 9.698059E-03, 1.048741E-02, 1.134099E-02, 1.226402E-02, 1.326221E-02, 1.434164E-02, 
                 1.550890E-02, 1.677117E-02, 1.813621E-02, 1.961231E-02, 2.120858E-02, 2.293476E-02, 2.480144E-02, 
                 2.682004E-02, 2.900295E-02, 3.136352E-02, 3.391613E-02, 3.667660E-02, 3.966167E-02, 4.288978E-02, 
                 4.638072E-02, 5.015556E-02, 5.423794E-02, 5.865249E-02, 6.342613E-02, 6.858833E-02, 7.417085E-02, 
                 8.020773E-02, 8.673591E-02, 9.379525E-02, 1.014295E-01, 1.096850E-01, 1.186123E-01, 1.282662E-01, 
                 1.387056E-01, 1.499951E-01]
		manual_nh2 = [2.548344E+11, 2.356551E+11, 2.179171E+11, 2.015157E+11, 1.863490E+11, 1.723241E+11, 
                1.593541E+11, 1.473597E+11, 1.362689E+11, 1.260124E+11, 1.165283E+11, 1.077577E+11, 9.964752E+10, 
                9.214758E+10, 8.521217E+10, 7.879855E+10, 7.286777E+10, 6.738348E+10, 6.231184E+10, 5.762191E+10, 
                5.328508E+10, 4.927460E+10, 4.556591E+10, 4.213643E+10, 3.896503E+10, 3.603232E+10, 3.332036E+10, 
                3.081251E+10, 2.849341E+10, 2.634885E+10, 2.262470E+10, 1.881360E+10, 1.602180E+10, 1.387110E+10, 
                1.185840E+10, 1.010000E+10, 8.648850E+09, 7.395170E+09, 6.319050E+09, 5.402330E+09, 4.617710E+09, 
                3.946950E+09, 3.373380E+09, 2.882940E+09, 2.463690E+09, 2.105340E+09, 1.799070E+09, 1.537290E+09, 
                1.313610E+09, 1.122970E+09]
		manual_T = [1.890000E+01, 1.836352E+01, 1.784226E+01, 1.733580E+01, 1.684372E+01, 1.636560E+01, 
              1.590106E+01, 1.544970E+01, 1.501115E+01, 1.458505E+01, 1.417105E+01, 1.376880E+01, 1.339650E+01, 
              1.309480E+01, 1.284920E+01, 1.263430E+01, 1.244060E+01, 1.226860E+01, 1.211470E+01, 1.196620E+01, 
              1.182410E+01, 1.173050E+01, 1.164030E+01, 1.153480E+01, 1.140280E+01, 1.141230E+01, 1.129260E+01, 
              1.129980E+01, 1.119581E+01, 1.109880E+01, 1.110150E+01, 1.122620E+01, 1.128450E+01, 1.130700E+01, 
              1.135870E+01, 1.147290E+01, 1.160570E+01, 1.167560E+01, 1.180180E+01, 1.199910E+01, 1.213550E+01, 
              1.227160E+01, 1.244650E+01, 1.262950E+01, 1.278760E+01, 1.304070E+01, 1.316610E+01, 1.344790E+01, 
              1.361030E+01, 1.390000E+01]
		manual_Vr = 		[-7.659916E+02,-7.274061E+02,-6.896902E+02,-6.528437E+02,-6.168665E+02,-5.817588E+02,-5.475205E+02,
   -5.141516E+02,-4.816521E+02,-4.500222E+02,-4.192615E+02,-3.893703E+02,-3.603486E+02,-3.321962E+02,-3.049133E+02,
   -2.784998E+02,-2.529557E+02,-2.282811E+02,-2.044758E+02,-1.815400E+02,-1.594736E+02,-1.382767E+02,-1.179491E+02,
   -9.849095E+01,-7.990225E+01,-6.218296E+01,-4.533308E+01,-2.935264E+01,-1.424161E+01, 0.000000E+00, 0.000000E+00, 
   0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 
   0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 
   0.000000E+00, 0.000000E+00, 0.000000E+00]
		manual_Vt = [1.590000E+02, 1.583720E+02, 1.576660E+02, 1.568800E+02, 1.560160E+02, 1.550770E+02, 
               1.540750E+02, 1.530300E+02, 1.519720E+02, 1.509490E+02, 1.500280E+02, 1.492890E+02, 1.487420E+02, 
               1.483550E+02, 1.480620E+02, 1.477610E+02, 1.473860E+02, 1.468940E+02, 1.463170E+02, 1.459610E+02, 
               1.459890E+02, 1.459750E+02, 1.461140E+02, 1.452810E+02, 1.449930E+02, 1.449940E+02, 1.450000E+02, 
               1.450000E+02, 1.450000E+02, 1.450000E+02, 1.450020E+02, 1.449990E+02, 1.449930E+02, 1.450080E+02, 
               1.450260E+02, 1.449240E+02, 1.450270E+02, 1.458040E+02, 1.461190E+02, 1.458800E+02, 1.463480E+02, 
               1.471000E+02, 1.469980E+02, 1.468980E+02, 1.479300E+02, 1.479680E+02, 1.480280E+02, 1.483530E+02, 
               1.489640E+02, 1.500000E+02]
		

		# Generate grid and attach to input
		
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, manual_grid[n-1])

		# Setup model parameters
		grid.gas_to_dust = 100.
		grid.T_cmb = 2.728 # K
		for pos in np.ndindex(n, 1, 1):
			# Manually regrid if requested
			id = pos[0]

			if id == 0:
				grid.SetZoneBoundary_sph1d(id, 0., 0.5*(manual_grid[id]+manual_grid[id+1]) )
			elif id == n-1:
				grid.SetZoneBoundary_sph1d(id, 0.5*(manual_grid[id-1]+manual_grid[id]), manual_grid[id] + 0.5*(manual_grid[id]-manual_grid[id-1]) )
			else:
				grid.SetZoneBoundary_sph1d(id, 0.5*(manual_grid[id-1]+manual_grid[id]), 0.5*(manual_grid[id]+manual_grid[id+1]))
				assert manual_grid[id - 1] < manual_grid[id]

			# Set physical parameters
			grid.cen_i[pos] = manual_grid[id]
			grid.n_H2[pos] = manual_nh2[id] # m^-3
			grid.T_k[pos] = manual_T[id] # K
			grid.X_mol[pos] = 1e-9 # Fraction
			grid.V_i[pos] = manual_Vr[id] # m/s
			grid.V_t[pos] = manual_Vt[id] # [m/s]
                        grid.kapp_d[pos] = "table,jena_thin_e5" 
			grid.T_d[pos] = manual_T[id] # K
		return

install_task(Task_P2a("task_p2a"))

################################################################################
class Task_PetersAGB(Task):
	"""
	Leiden benchmark problem, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Peter's AGB model(2016 Jan)"

		# Explanation
		self.expl = "a spike-like 1D model"

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file")
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		n = 100
		GridRadius = [275000000000000.12,291763090668161.44,309548003913591.8,328417026661780.5,
		348436242643232.75,369675763834149.0,392209976004962.1,416117799236702.3,441482964317577.2,
		468394305987774.1,496946074059496.1,527238263501845.7,559376964646592.44,593474734741325.9,
		629650992151259.8,668032434590275.2,708753482845945.5,751956751552569.9,797793548660976.2,
		846424405354351.1,898019638265985.8,952759945967957.1,1010837041819787.4,1072454325393463.9,
		1137827594826309.0,1207185802596524.8,1280771857368316.0,1358843474714841.2,1441674079698420.2,
		1529553764469047.8,1622790304234946.2,1721710235163322.5,1826659997986392.0,1938007151317843.5,
		2056141658929065.5,2181477255493471.0,2314452895582084.5,2455534290985106.0,2605215541743515.0,
		2764020866602968.0,2932506438950441.5,3111262334663507.5,3300914598694051.0,3502127437624108.0,
		3715605545872647.5,3942096573700249.5,4182393745655088.0,4437338638630986.0,4707824129266264.0,
		4994797521006822.0,5299263861784053.0,5622289463928064.0,5965005638642491.0,6328612658121271.0,
		6714383959182190.0,7123670603140833.0,7557906007542857.0,8018610966327915.0,8507398976005250.0,
		9025981886496154.0,9576175896432578.0,1.0159907913909524e+16,1.0779222304967788e+16,
		1.1436288053442168e+16,1.2133406357249614e+16,1.2873018687723194e+16,1.3657715340214914e+16,
		1.4490244505916646e+16,1.5373521896667518e+16,1.6310640956459586e+16,1.7304883695401178e+16,
		1.8359732184084876e+16,1.9478880748611972e+16,2.0666248908987636e+16,2.192599510619536e+16,
		2.3262531266031144e+16,2.4680538250698012e+16,2.618498225228166e+16,2.7781132185514972e+16,
		2.9474578140751316e+16,3.1271250961766532e+16,3.3177443016962736e+16,3.519983023671205e+16,
		3.734549549402817e+16,3.962195341044183e+16,4.203717667396527e+16,4.459962396130762e+16,
		4.7318269562141624e+16,5.020263480916206e+16,5.326282142402255e+16,5.650954689592405e+16,
		5.995418201677188e+16,6.360879070434521e+16,6.74861722529596e+16,7.1599906159590264e+16,
		7.596439968244531e+16,8.059493829855208e+16,8.550773923707354e+16,9.072000827583883e+16,
		9.625000000000005e+16]
                for i in range(n):
                        GridRadius[i] /= 3.08567758e18
		density = [1.2764256315130103e-15,6.707988766139324e-16,4.692645667454499e-16,3.577309614525716e-16,2.84421281385618e-16,2.3186116511070704e-16,1.9217238391105812e-16,1.6116420902792403e-16,1.3635252487609542e-16,1.1614784053937355e-16,9.947379326219946e-17,8.556912080769166e-17,7.387670095869923e-17,6.397748623289185e-17,5.554917531359566e-17,4.833927507335331e-17,4.214692328806059e-17,3.6810254431867896e-17,3.2197391090083405e-17,2.819987630447864e-17,2.4727791459620792e-17,2.170606447799231e-17,1.9071635536054977e-17,1.6771251679578362e-17,1.4759730124859977e-17,1.2998575935015263e-17,1.1454871168616828e-17,1.0100374474200584e-17,8.910785591449585e-18,7.865140349376977e-18,6.945309860814119e-18,5.112986330461622e-20,4.518600409117765e-20,3.994713003533849e-20,3.532710366823901e-20,3.125077498727052e-20,2.765249317640863e-20,2.4474836483718636e-20,2.1667525063705502e-20,1.9186487901400747e-20,1.6993059986869342e-20,1.5053289979702155e-20,1.3337341903197279e-20,1.1818977098314568e-20,1.0475104872354933e-20,9.285392093233065e-21,8.231923482740108e-21,7.298905610827523e-21,6.472408634837501e-21,5.740140700305221e-21,5.09125065359201e-21,4.516155335390327e-21,4.0063882477294764e-21,3.554466831440934e-21,3.153775969183295e-21,2.798465651760191e-21,2.4833610214297322e-21,2.2038832424791197e-21,1.9559798525958787e-21,1.736063423541765e-21,1.540957510557839e-21,1.3678490003302894e-21,1.214246080209689e-21,1.0779411492109056e-21,9.56978076254279e-22,8.496232849425145e-22,7.543402084381674e-22,6.697667140172861e-22,5.946951457465603e-22,5.280546764164605e-22,4.688956971852112e-22,4.163760060498952e-22,3.6974858487152473e-22,3.283507797616772e-22,2.915947216440156e-22,2.589588431241944e-22,2.299803647773627e-22,2.0424863888654456e-22,1.8139925179484185e-22,1.6110879759096209e-22,1.4309024602600717e-22,1.2708883652778872e-22,1.1287843808560541e-22,1.0025832175154618e-22,8.905029865742015e-23,7.909618187742442e-23,7.025553526265568e-23,6.240367660986099e-20,5.542990627035478e-20,4.923593561381131e-20,4.3734492687326296e-20,3.884808499753359e-20,3.450790163260961e-20,3.065283896639541e-20,2.7228635979555116e-20,2.4187106819650593e-20,1.2891275776308163e-21,1.1451415139442073e-21,1.0172432258039388e-21,0.0]
                CGSDensity2MKSN_h2 = 1E6 * 1. / (1.36986301369863 * 2. * 1.6726219E-24)
                for i in range(n):
                        density[i] *= CGSDensity2MKSN_h2
		temperature = [894.4271909999156,868.35271199125,843.0383602041342,818.4619762929202,794.6020469046122,
   771.4376858467884,748.9486158045208,727.1151505902841,705.9181779113185,685.339142639362,
   665.360030568105,645.96335264415,627.1321296576734,608.8498773793879,591.1005921307924,
   573.8687367750828,557.1392271164561,540.8974186959051,525.1290939719432,509.8204498750394,
   494.9580857248685,480.52899149979856,466.52053644834876,452.92045803264796,439.7168511942142,
   426.89815793266155,414.4531571882083,402.3709550191329,390.6409750655784,379.2529492913563,
   368.19690899564824,357.463176086735,347.04235461011535,336.9253225235984,327.10322371217006,
   317.5674602356428,308.3096848023037,299.32179346197125,290.59591851206596,282.12442161048364,
   273.8998870892437,265.9151154630588,258.16311712714383,250.63710623874718,243.3304947770485,
   236.23688677622306,229.35007272662725,222.6640241391941,216.17288826830065,209.87098298846064,
   203.75279182038736,197.81295910204307,192.0462853004766,186.44772246031772,181.01236978497042,
   175.73546934661178,170.61240192126516,165.63868294527893,160.80995858969288,156.1220019490369,
   151.57070934123936,147.1520967153995,142.86229616428125,138.697552538474,134.6542201592621,
   130.72875962731425,126.91773472441481,123.21780940550792,119.62574487843813,116.13839676881504,
   112.75271236753589,109.46572795854121,106.27456622447966,103.176433727996,100.16861846645205,
   97.24848749792588,94.41348463642497,91.66112821428236,88.98900890978999,86.39478763815531,
   83.87619350394671,81.43102181322405,79.05713214362525,76.7524464707088,74.51494734892252,
   72.34267614559589,70.23373132642047,68.18626679090713,66.19849025637326,64.26866168903554,
   62.39509178084405,60.57614047071595,58.810215508881996,57.09577106308307,55.43130636539927,
   53.81536439852793,52.246530620359245,50.7234317257325,49.2447344442898,47.80914437337573]
		Vr = [3.1176914536240035,5.270368822273995,6.692991265396655,7.799848677540593,8.715359463301386,
               9.497822849454762,10.180422150479457,10.784322077659105,11.324080617258556,11.810266104988418,
               12.25086712345085,12.652114830531321,13.018992815996272,13.355568975686275,13.665220281783784,
               13.950790130643508,14.214701618033823,14.459041067952562,14.685620923469099,14.896027971480422,
               15.091660923258056,15.273760124138654,15.443431345297377,15.601665058849512,15.749352218783185,
               15.88729730532513,16.016229201910757,16.136810337818126,16.249644429815596,16.355283082181074,
               16.45423144888715,16.546953119560662,16.633874358468148,16.71538780072406,16.791855690346647,
               16.863612729369915,16.93096859498308,16.99421017188439,17.053603539159134,17.109395744609266,
               17.161816394258928,17.211079080494173,17.25738266877828,17.30091245996946,17.34184124284083,
               17.380330249371795,17.416530023673293,17.450581213968892,17.482615295833043,17.512755233849177,
               17.541116087963374,17.56780557004934,17.592924555546578,17.61656755446986,17.6388231455996,
               17.659774377239188,17.679499137555947,17.698070497200238,17.7155570266148,17.73202309019876,
               17.747529119272304,17.76213186559581,17.7758846370265,17.7888375167451,17.801037567350345,
               17.81252902099998,17.82335345667,17.833549965508528,17.84315530517518,17.85220404398002,
               17.86072869556702,17.868759844824844,17.876326265651365,17.883455031147598,17.89017161677065,
               17.896499996933446,17.902462735501043,17.908081070598726,17.913374994115557,17.91836332625835,
               17.923063785484555,17.92749305411868,17.931666839934625,17.935599933966106,17.939306264788797,
               17.94279894950059,17.946090341610674,17.94919207603358,17.952115111370947,17.95486976965131,
               17.957465773686835,17.959912282195265,17.962217922825516,17.96439082321632,17.966438640208764,
               17.968368587325774,17.97018746062438,17.971901663019615,17.97351722717283,17.97503983703116]

		

		# Generate grid and attach to input
		
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, GridRadius[n-1])

		
		
		
		# Setup model parameters
		grid.gas_to_dust = 100.
		grid.T_cmb = 2.728 # K
		for pos in np.ndindex(n, 1, 1):
			# Manually regrid if requested
			id = pos[0]

			if id == 0:
				grid.SetZoneBoundary_sph1d(id, 0., 0.5*(GridRadius[id]+GridRadius[id+1]) )
			elif id == n-1:
				grid.SetZoneBoundary_sph1d(id, 0.5*(GridRadius[id-1]+GridRadius[id]), GridRadius[id] + 0.5*(GridRadius[id]-GridRadius[id-1]) )
			else:
				grid.SetZoneBoundary_sph1d(id, 0.5*(GridRadius[id-1]+GridRadius[id]), 0.5*(GridRadius[id]+GridRadius[id+1]))
				assert GridRadius[id - 1] < GridRadius[id]

			# Set physical parameters
			grid.cen_i[pos] = GridRadius[id] #cm2pc
			grid.n_H2[pos] = density[id] # m^-3
			grid.T_k[pos] = temperature[id] # K
			grid.X_mol[pos] = 3e-4/30. # Fraction
			grid.V_i[pos] = Vr[id] * 1e3  # m/s
			grid.V_t[pos] = 1e3 # [m/s]
                        grid.kapp_d[pos] = "table,jena_thin_e5" 
			grid.T_d[pos] = temperature[id] # K
		return

install_task(Task_PetersAGB("task_petersagb"))

################################################################################

class Task_Benchmark_Ronny(Task):
	"""
	Leiden benchmark problem, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "2002 Benchmark Problem 2a"

		# Explanation
		self.expl = "Some documentation on the task"

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file")
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		n = 119
		manual_grid = [1.298701E-05, 1.300326E-05, 1.301951E-05, 1.306850E-05, 1.311749E-05, 1.319996E-05, 1.328243E-05, 1.339963E-05, 1.351684E-05, 1.367057E-05, 1.382431E-05, 1.401696E-05, 1.420962E-05, 1.444423E-05, 1.467885E-05, 1.495921E-05, 1.523956E-05, 1.557027E-05, 1.590097E-05, 1.628760E-05, 1.667422E-05, 1.712346E-05, 1.757269E-05, 1.809253E-05, 1.861237E-05, 1.921235E-05, 1.981234E-05, 2.050384E-05, 2.119535E-05, 2.199193E-05, 2.278852E-05, 2.370637E-05, 2.462421E-05, 2.568266E-05, 2.674111E-05, 2.796330E-05, 2.918550E-05, 3.059922E-05, 3.201294E-05, 3.365160E-05, 3.529025E-05, 3.719413E-05, 3.909801E-05, 4.131584E-05, 4.353367E-05, 4.612456E-05, 4.871545E-05, 5.175130E-05, 5.478716E-05, 5.835577E-05, 6.192437E-05, 6.613322E-05, 7.034207E-05, 7.532325E-05, 8.030442E-05, 8.622075E-05, 9.213708E-05, 9.919002E-05, 1.062430E-04, 1.146826E-04, 1.231223E-04, 1.332603E-04, 1.433982E-04, 1.556242E-04, 1.678502E-04, 1.826531E-04, 1.974560E-04, 2.154519E-04, 2.334478E-04, 2.554154E-04, 2.773830E-04, 3.043107E-04, 3.312384E-04, 3.643853E-04, 3.975321E-04, 4.385083E-04, 4.794844E-04, 5.303569E-04, 5.812293E-04, 6.446620E-04, 7.080946E-04, 7.875340E-04, 8.669733E-04, 9.668964E-04, 1.066819E-03, 1.193065E-03, 1.319310E-03, 1.479522E-03, 1.639734E-03, 1.843964E-03, 2.048193E-03, 2.309706E-03, 2.571218E-03, 2.907598E-03, 3.243977E-03, 3.678624E-03, 4.113271E-03, 4.677458E-03, 5.241646E-03, 5.977338E-03, 6.713031E-03, 7.676781E-03, 8.640530E-03, 9.908863E-03, 1.117720E-02, 1.285411E-02, 1.453102E-02, 1.675843E-02, 1.898585E-02, 2.195828E-02, 2.493071E-02, 2.891589E-02, 3.290107E-02, 3.826910E-02, 4.363712E-02, 5.090181E-02, 5.816650E-02, 6.804429E-02, 7.792207E-02]
		manual_nh2 = [6.335209E+11, 6.179611E+11, 6.029676E+11, 5.609226E+11, 5.231267E+11, 4.676037E+11, 4.204731E+11, 3.648513E+11, 3.195783E+11, 2.719305E+11, 2.341957E+11, 1.970475E+11, 1.680849E+11, 1.406842E+11, 1.194758E+11, 9.987248E+10, 8.472592E+10, 7.090516E+10, 6.020965E+10, 5.051243E+10, 4.298234E+10, 3.617070E+10, 3.085873E+10, 2.605237E+10, 2.228734E+10, 1.887457E+10, 1.618962E+10, 1.374936E+10, 1.182195E+10, 1.006467E+10, 8.671959E+09, 7.397944E+09, 6.385368E+09, 5.455996E+09, 4.715688E+09, 4.034019E+09, 3.490141E+09, 2.987822E+09, 2.586629E+09, 2.215057E+09, 1.918156E+09, 1.642483E+09, 1.422233E+09, 1.217274E+09, 1.053627E+09, 9.010477E+08, 7.793620E+08, 6.657217E+08, 5.752387E+08, 4.906267E+08, 4.233975E+08, 3.604665E+08, 3.105902E+08, 2.638689E+08, 2.269485E+08, 1.923485E+08, 1.650976E+08, 1.395553E+08, 1.195126E+08, 1.007286E+08, 8.604869E+07, 7.229574E+07, 6.159456E+07, 5.157524E+07, 4.381568E+07, 3.655663E+07, 3.096271E+07, 2.573512E+07, 2.172777E+07, 1.798754E+07, 1.513609E+07, 1.247856E+07, 1.046409E+07, 8.589661E+06, 7.177221E+06, 5.865296E+06, 4.882769E+06, 3.971911E+06, 3.294044E+06, 2.666909E+06, 2.203201E+06, 1.775124E+06, 1.460689E+06, 1.171072E+06, 9.597686E+05, 7.656016E+05, 6.249068E+05, 4.959360E+05, 4.031309E+05, 3.182717E+05, 2.576372E+05, 2.023358E+05, 1.631019E+05, 1.274117E+05, 1.022725E+05, 7.946443E+04, 6.351502E+04, 4.908319E+04, 3.906458E+04, 3.002368E+04, 2.379336E+04, 1.818631E+04, 1.435074E+04, 1.090830E+04, 8.570822E+03, 6.478669E+03, 5.068580E+03, 3.809944E+03, 2.967935E+03, 2.218428E+03, 1.720749E+03, 1.278962E+03, 9.877978E+02, 7.300429E+02, 5.614355E+02, 4.125840E+02, 3.159425E+02, 2.308582E+02, 1.760308E+02]
		manual_T = [2.000000E+03, 1.998250E+03, 1.996504E+03, 1.991262E+03, 1.986054E+03, 1.977360E+03, 1.968757E+03, 1.956687E+03, 1.944795E+03, 1.929460E+03, 1.914415E+03, 1.895958E+03, 1.877928E+03, 1.856523E+03, 1.835701E+03, 1.811551E+03, 1.788157E+03, 1.761486E+03, 1.735761E+03, 1.706815E+03, 1.679015E+03, 1.648058E+03, 1.618451E+03, 1.585758E+03, 1.554624E+03, 1.520478E+03, 1.488098E+03, 1.452786E+03, 1.419443E+03, 1.383254E+03, 1.349228E+03, 1.312445E+03, 1.278006E+03, 1.240905E+03, 1.206316E+03, 1.169162E+03, 1.134670E+03, 1.097714E+03, 1.063552E+03, 1.027029E+03, 9.934091E+02, 9.575344E+02, 9.246518E+02, 8.896206E+02, 8.576474E+02, 8.236333E+02, 7.927200E+02, 7.598740E+02, 7.301487E+02, 6.985988E+02, 6.701675E+02, 6.400188E+02, 6.129653E+02, 5.843004E+02, 5.586874E+02, 5.315673E+02, 5.074373E+02, 4.819022E+02, 4.592788E+02, 4.353496E+02, 4.142389E+02, 3.919187E+02, 3.723109E+02, 3.515867E+02, 3.334581E+02, 3.143023E+02, 2.976166E+02, 2.799892E+02, 2.646998E+02, 2.485499E+02, 2.346014E+02, 2.198692E+02, 2.071992E+02, 1.938179E+02, 1.823585E+02, 1.702558E+02, 1.599353E+02, 1.490352E+02, 1.397794E+02, 1.300032E+02, 1.217370E+02, 1.130051E+02, 1.056532E+02, 9.788599E+01, 9.137404E+01, 8.449319E+01, 7.874876E+01, 7.267775E+01, 6.763088E+01, 6.229596E+01, 5.787978E+01, 5.321044E+01, 4.936160E+01, 4.529106E+01, 4.195001E+01, 3.841552E+01, 3.552674E+01, 3.246980E+01, 2.998190E+01, 2.734834E+01, 2.521410E+01, 2.295414E+01, 2.113042E+01, 1.919861E+01, 1.764629E+01, 1.600137E+01, 1.468518E+01, 1.328994E+01, 1.217827E+01, 1.099937E+01, 1.006404E+01, 9.071747E+00, 8.287805E+00, 7.455775E+00, 6.801226E+00, 6.106231E+00, 5.561801E+00, 4.983476E+00, 4.532358E+00]
		manual_Vr = [1.500000E-01, 1.533928E-01, 1.568150E-01, 1.673079E-01, 1.780584E-01, 1.967196E-01, 2.160615E-01, 2.446634E-01, 2.745006E-01, 3.153838E-01, 3.581006E-01, 4.139920E-01, 4.722559E-01, 5.460549E-01, 6.225962E-01, 7.171462E-01, 8.145342E-01, 9.323964E-01, 1.052827E+00, 1.196075E+00, 1.341187E+00, 1.511231E+00, 1.681962E+00, 1.879425E+00, 2.075914E+00, 2.300555E+00, 2.522101E+00, 2.772795E+00, 3.017870E+00, 3.292640E+00, 3.558945E+00, 3.855047E+00, 4.139615E+00, 4.453656E+00, 4.752989E+00, 5.081078E+00, 5.391300E+00, 5.729218E+00, 6.046240E+00, 6.389612E+00, 6.709299E+00, 7.053765E+00, 7.372088E+00, 7.713455E+00, 8.026626E+00, 8.361005E+00, 8.665589E+00, 8.989498E+00, 9.282501E+00, 9.592945E+00, 9.871865E+00, 1.016638E+01, 1.042924E+01, 1.070593E+01, 1.095127E+01, 1.120878E+01, 1.143566E+01, 1.167316E+01, 1.188110E+01, 1.209825E+01, 1.228721E+01, 1.248408E+01, 1.265437E+01, 1.283143E+01, 1.298367E+01, 1.314167E+01, 1.327672E+01, 1.341664E+01, 1.353556E+01, 1.365855E+01, 1.376250E+01, 1.386985E+01, 1.396008E+01, 1.405313E+01, 1.413090E+01, 1.421102E+01, 1.427761E+01, 1.434613E+01, 1.440278E+01, 1.446100E+01, 1.450889E+01, 1.455805E+01, 1.459827E+01, 1.463954E+01, 1.467311E+01, 1.470753E+01, 1.473539E+01, 1.476393E+01, 1.478692E+01, 1.481044E+01, 1.482929E+01, 1.484857E+01, 1.486394E+01, 1.487964E+01, 1.489210E+01, 1.490483E+01, 1.491487E+01, 1.492513E+01, 1.493318E+01, 1.494139E+01, 1.494781E+01, 1.495436E+01, 1.495945E+01, 1.496463E+01, 1.496864E+01, 1.497273E+01, 1.497588E+01, 1.497908E+01, 1.498154E+01, 1.498404E+01, 1.498594E+01, 1.498788E+01, 1.498934E+01, 1.499084E+01, 1.499197E+01, 1.499311E+01, 1.499397E+01, 1.499485E+01, 1.499550E+01]

		# Generate grid and attach to input
		
		INP_DICT["pygrid"] = grid = Grid_sph1d(n+1, manual_grid[n-1] + 0.5*(manual_grid[n-1]-manual_grid[n-2]))

		# Setup model parameters
		grid.gas_to_dust = 100.0
		grid.T_cmb = 2.728 # K
		for pos in np.ndindex(n+1, 1, 1):
			# Manually regrid if requested
			id = pos[0]-1

			if id == -1:
				grid.SetZoneBoundary_sph1d(pos[0], 0., manual_grid[id+1] )
			elif id == 0:
				grid.SetZoneBoundary_sph1d(pos[0], manual_grid[id], 0.5*(manual_grid[id]+manual_grid[id+1]) )
			elif id == n-1:
				grid.SetZoneBoundary_sph1d(pos[0], 0.5*(manual_grid[id-1]+manual_grid[id]), manual_grid[id] + 0.5*(manual_grid[id]-manual_grid[id-1]) )
			else:
				grid.SetZoneBoundary_sph1d(pos[0], 0.5*(manual_grid[id-1]+manual_grid[id]), 0.5*(manual_grid[id]+manual_grid[id+1]))
				assert manual_grid[id - 1] < manual_grid[id]
				
				


			if id == -1:
				grid.cen_i[pos] = 0.5*manual_grid[id+1]
				grid.X_mol[pos] = 0.0 # Fraction# Set physical parameters
				grid.n_H2[pos] = manual_nh2[id+1]*1e6 # m^-3
				grid.T_k[pos] = manual_T[id+1] # K
				grid.V_i[pos] = 0. # m/s
				grid.kapp_d[pos] = Type.KappLLaw("['0.25mm','1e6cm^2g^-1',-1]")
			else:
				grid.cen_i[pos] = manual_grid[id]
				grid.n_H2[pos] = manual_nh2[id]*1e6 # m^-3
				grid.T_k[pos] = manual_T[id] # K
				grid.X_mol[pos] = 3e-4 # Fraction
				grid.V_i[pos] = manual_Vr[id]*1e3 # m/s
				grid.kapp_d[pos] = Type.KappLLaw("['0.25mm','10cm^2g^-1',-1]")
			grid.T_d[pos] = grid.T_k[pos]
			grid.V_t[pos] = 1000. # [m/s]Ft
 			grid.X_pH2[pos] = 0.25
 			grid.X_oH2[pos] = 0.75 
			
				
		return

install_task(Task_Benchmark_Ronny("task_benchmark_ronny"))

################################################################################

class Task_N1333(Task):
	"""
	N1333I4A modeling (1D)
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "N1333"

		# Explanation
		self.expl = "1D model for N1333"

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("abd1", Type.Fraction, None, "Inner molecular abundance"),
			Key("abd2", Type.Fraction, None, "Outer molecular abundance"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		radius_out = 0.035 # [pc] cloud radius
		radius_in = 0.0 # [pc] cloud radius
		n = 200
		# stretching grid
		stretch_ratior=1.02
		dr = (radius_out-radius_in)*(stretch_ratior-1.)/(stretch_ratior**(n)-1.)
		
		abd1 = INP_DICT["abd1"]
		abd2 = INP_DICT["abd2"]
		
		from numpy import zeros
		manual_grid=zeros(n+1)
		for pos in range(0,n+1):
			if (pos==0):
				manual_grid[pos] = radius_in
			else:
				manual_grid[pos] = manual_grid[pos-1]+dr
				dr=dr*stretch_ratior

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius_out)

		# Setup model parameters
		grid.gas_to_dust = 100
		grid.T_cmb = 2.73 # K
		
		
		for pos in np.ndindex(n, 1, 1):
			# Manually regrid if requested
			id = pos[0]

			grid.SetZoneBoundary_sph1d(id, manual_grid[id], manual_grid[id+1])
			assert manual_grid[id] < manual_grid[id+1]

			# Set physical parameters
			if(grid.cen_i[pos] < 23.9/2.06e5):			
				grid.n_H2[pos] = 5e15 # m^-3
				grid.kapp_d[pos] = Type.KappLLaw("['0.16mm','0.23cm^2g^-1',-2]")
			else:
				grid.n_H2[pos] = 6.3e12*(grid.cen_i[pos]/0.005)**-1.8 # m^-3
				grid.kapp_d[pos] = Type.KappLLaw("['0.16mm','0.23cm^2g^-1',-2]")
			grid.T_k[pos] = 71.*(grid.cen_i[pos]/0.00032467532)**-1 +60.7*(grid.cen_i[pos]/0.00032467532)**-0.4 # K
			if(grid.T_k[pos]>100.):
				grid.X_mol[pos] = abd1 # Fraction
				if(grid.T_k[pos]>250.):
					grid.T_k[pos] = 250.
			else:
				grid.X_mol[pos] = abd2 # Fraction
			grid.V_t[pos] = 400. # [m/s]
			grid.V_i[pos] = -1100.*(grid.cen_i[pos]/0.005)**-0.5 # m/s
			grid.T_d[pos] = grid.T_k[pos]
			
		return

install_task(Task_N1333("task_n1333"))

################################################################################
class Task_L1498(Task):
	"""
	L1498-1D
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Leiden H2O benchmark problem generator"

		# Explanation
		self.expl = "Some documentation on the task"

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		radius_out = 0.27 # [pc] cloud radius
		radius_in = 0.0 # [pc] cloud radius
		n = 200
		dr =(radius_out-radius_in)/n
		
		from numpy import zeros
		manual_grid=zeros(n+1)
		for pos in range(0,n+1):
			manual_grid[pos] = radius_in+pos*dr

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius_out)

		# Setup model parameters
		grid.gas_to_dust = 100
		grid.T_cmb = 2.73 # K
		xmol=INP_DICT["abundance"]
		
		
		for pos in np.ndindex(n, 1, 1):
			# Manually regrid if requested
			id = pos[0]
			if id == 0:
				grid.SetZoneBoundary_sph1d(id, manual_grid[id], manual_grid[id+1])
				assert manual_grid[id] < manual_grid[id+1]

			# Set physical parameters
			if(grid.cen_i[pos] < 4e17/3.08e18):			
				grid.n_H2[pos] = 0.94e11/(1+(grid.cen_i[pos]/(1.51875e17/3.08e18))**3.5) # m^-3
			else:
				grid.n_H2[pos] = 1e9 # m^-3
			grid.T_k[pos] = 10.0 # K
			if(grid.cen_i[pos]<(1.15e17/3.08e18)):
				grid.X_mol[pos] = 0.0 # Fraction
			elif(grid.cen_i[pos]<(4e17/3.08e18)):
				grid.X_mol[pos] = 3e-9 # Fraction
			else:
				grid.X_mol[pos] = 12e-9 # Fraction
			grid.V_t[pos] = 200 # [m/s]
			grid.V_i[pos] = 0.0 # m/s
		return

install_task(Task_L1498("task_l1498"))

################################################################################
class Task_Leiden1D(Task):
	"""
	Leiden benchmark problem, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Leiden H2O benchmark problem generator"

		# Explanation
		self.expl = "Some documentation on the task"

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("xmol", Type.Fraction, None, "Molecular abundance"),
			Key("vgrad", Type.Velo, None, "Radial velocity gradient"),
			Key("tk", Type.Temp, None, "Kinetic temperature"),
			Key("tcmb", Type.Temp, "0.0K", "Brightness temperature of background radiation field"),
			Key("ndiv", Type.PosInt, Type.Optional, "Number of shells (standard 200 shell grid used if not given)"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		radius = 0.1 # [pc] cloud radius
		n = INP_DICT["ndiv"]
		if n == None:
			n = 200
			manual_grid = [.100000000E-02, .626483225E-02, .112496854E-01, .159694484E-01, .204382186E-01, .246693436E-01, .286754610E-01, .324685367E-01, .360599000E-01, .394602777E-01, .426798264E-01, .457281623E-01, .486143904E-01, .513471314E-01, .539345477E-01, .563843674E-01, .587039079E-01, .609000973E-01, .629794953E-01, .649483128E-01, .668124302E-01, .685774156E-01, .702485405E-01, .718307966E-01, .733289096E-01, .747473544E-01, .760903674E-01, .773619603E-01, .785659310E-01, .797058755E-01, .807851988E-01, .818071247E-01, .827747054E-01, .836908311E-01, .845582380E-01, .853795169E-01, .861571210E-01, .868933728E-01, .875904713E-01, .882504988E-01, .888754266E-01, .894671213E-01, .900273502E-01, .905577866E-01, .910600148E-01, .915355350E-01, .919857675E-01, .924120570E-01, .928156768E-01, .931978325E-01, .935596654E-01, .939022564E-01, .942266287E-01, .945337512E-01, .948245412E-01, .950998672E-01, .953605517E-01, .956073731E-01, .958410688E-01, .960623368E-01, .962718380E-01, .964701980E-01, .966580095E-01, .968358333E-01, .970042006E-01, .971636142E-01, .973145504E-01, .974574599E-01, .975927697E-01, .977208838E-01, .978421848E-01, .979570352E-01, .980657780E-01, .981687379E-01, .982662225E-01, .983585230E-01, .984459150E-01, .985286595E-01, .986070038E-01, .986811818E-01, .987514151E-01, .988179134E-01, .988808755E-01, .989404892E-01, .989969328E-01, .990503747E-01, .991009746E-01, .991488837E-01, .991942450E-01, .992371940E-01, .992778591E-01, .993163616E-01, .993528166E-01, .993873329E-01, .994200137E-01, .994509566E-01, .994802539E-01, .995079932E-01, .995342574E-01, .995591249E-01, .995826699E-01, .996049629E-01, .996260703E-01, .996460553E-01, .996649774E-01, .996828933E-01, .996998565E-01, .997159175E-01, .997311245E-01, .997455227E-01, .997591553E-01, .997720629E-01, .997842841E-01, .997958554E-01, .998068113E-01, .998171846E-01, .998270063E-01, .998363056E-01, .998451104E-01, .998534470E-01, .998613403E-01, .998688138E-01, .998758898E-01, .998825896E-01, .998889375E-01, .998949433E-01, .999006298E-01, .999060139E-01, .999111116E-01, .999159383E-01, .999205083E-01, .999248353E-01, .999289321E-01, .999328111E-01, .999364838E-01, .999399612E-01, .999432537E-01, .999463711E-01, .999493227E-01, .999521173E-01, .999547633E-01, .999572686E-01, .999596407E-01, .999618867E-01, .999640132E-01, .999660266E-01, .999679329E-01, .999697379E-01, .999714469E-01, .999730650E-01, .999745970E-01, .999760476E-01, .999774211E-01, .999787215E-01, .999799527E-01, .999811185E-01, .999822223E-01, .999832674E-01, .999842569E-01, .999851937E-01, .999860808E-01, .999869207E-01, .999877159E-01, .999884689E-01, .999891818E-01, .999898568E-01, .999905415E-01, .999911429E-01, .999917114E-01, .999922498E-01, .999927596E-01, .999932425E-01, .999936997E-01, .999941328E-01, .999945429E-01, .999949314E-01, .999952993E-01, .999956477E-01, .999959777E-01, .999962903E-01, .999965863E-01, .999968667E-01, .999971323E-01, .999973838E-01, .999977677E-01, .999979928E-01, .999981857E-01, .999983665E-01, .999985387E-01, .999987031E-01, .999988601E-01, .999990103E-01, .999991538E-01, .999992911E-01, .999994225E-01, .999995483E-01, .999996687E-01, .999997841E-01, .999998947E-01, .100000000E+00]
		else:
			manual_grid = None

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius)

		# Setup model parameters
		grid.gas_to_dust = 0
		grid.T_cmb = INP_DICT["tcmb"] # K
		for pos in np.ndindex(n, 1, 1):
			# Manually regrid if requested
			id = pos[0]
			if manual_grid != None:
				if id == 0:
					grid.SetZoneBoundary_sph1d(id, 0, manual_grid[id])
				else:
					grid.SetZoneBoundary_sph1d(id, manual_grid[id - 1], manual_grid[id])
					assert manual_grid[id - 1] < manual_grid[id]

			# Set physical parameters
			grid.n_H2[pos] = 1e4 * 1e6 # m^-3
			grid.T_k[pos] = INP_DICT["tk"] # K
			grid.X_mol[pos] = INP_DICT["xmol"] # Fraction
			grid.V_i[pos] = INP_DICT["vgrad"] * grid.cen_i[pos] # m/s
		return

install_task(Task_Leiden1D("task_leiden1d"))

################################################################################

class Task_Leiden3D(Task_Leiden1D):
	"""
	Leiden benchmark problem, 3D version
	"""
	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_rec3d
		# Setup dimensions and gridding
		radius = 0.1 # [pc] cloud radius
		length = radius * 2 # [pc] length of each side of the box should be 2*radius
		n = INP_DICT["ndiv"] # number of divisions per axis

		# Set n=16 if not given by user
		if n == None:
			n = 16

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_rec3d((n, n, n), (length, length, length))

		# Setup model parameters
		grid.gas_to_dust = 0
		grid.T_cmb = INP_DICT["tcmb"]
		for pos in np.ndindex(n, n, n):
			# Calculate radius at pos
			di = grid.cen_i[pos] - grid.cen[0]
			dj = grid.cen_j[pos] - grid.cen[1]
			dk = grid.cen_k[pos] - grid.cen[2]
			r = sqrt(di**2 + dj**2 + dk**2)

			# Set parameters only if within cloud radius
			if r <= radius:
				# Set physical parameters
				grid.n_H2[pos] = 1e4 * 1e6 # m^-3
				grid.T_k[pos] = INP_DICT["tk"] # K
				grid.X_mol[pos] = INP_DICT["xmol"] # Fraction

				# Project radial velocity onto Cartesian coordinates
				vel = phys.Vec3_Normalize([di, dj, dk])
				vel = phys.Vec3_Scale(vel, INP_DICT["vgrad"] * r) # m/s
				grid.V_i[pos] = vel[0]
				grid.V_j[pos] = vel[1]
				grid.V_k[pos] = vel[2]
		return

install_task(Task_Leiden3D("task_leiden3d"))

################################################################################

class Task_W3OH1D(Task):
	"""
	1D toy model for molecular gas in the W3(OH) region
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "W3(OH) model generator"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output file name"),
			Key("xmol", Type.Fraction, None, "Molecular abundance"),
			Key("nshell", Type.PosInt, None, "Number of shells"),
			Key("radius", Type.Length, None, "Radius of model"),
			Key("n0", Type.NumDens, None, "Gas density"),
			Key("t0", Type.Temp, None, "Kinetic temperature"),
			Key("v0", Type.Velo, None, "Radial velocity gradient [velocity/pc]"),
			Key("vt", Type.Velo, None, "Turbulent line width"),
			#Key("tcen", Type.Temp, None, "Temperature of central source"),
			Key("rcen", Type.Length, None, "Radius of central cavity"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions
		n = INP_DICT["nshell"]
		length = INP_DICT["radius"] / Unit.pc # [pc]

		# Generate grid and attach to input
		grid = INP_DICT["pygrid"] = Grid_sph1d(n, length)

		# Radius of central "hole"
		rcen = INP_DICT["rcen"] / Unit.pc # [pc]

		# Regrid with inner radius arbitrarily set at 100 au (esc 09Dec30)
		# grid.Grid_sph1d_log10(100.0 * Unit.au / Unit.pc)
		grid.Grid_sph1d_log10(rcen)

		# Reference radius arbitrarily set at 500 au
		# DO NOT CHANGE unless for good reason
		r0 = 500 * Unit.au / Unit.pc # [pc]

		# Other model parameters
		grid.gas_to_dust = 100.0
		grid.T_cmb = 2.728 # K

		def T_k(r):
			return INP_DICT["t0"] * (r / r0)**-0.4 # [K]

		# Model description
		for pos in np.ndindex(n, 1, 1):
			# Get radius
			r = grid.cen_i[pos]

			# Set parameters
			grid.n_H2[pos] = INP_DICT["n0"] * (r / r0)**-1.5 # [m^-3]
			grid.X_mol[pos] = INP_DICT["xmol"] # Fraction
			grid.T_k[pos] =  T_k(r) # [K]
			grid.V_i[pos] = INP_DICT["v0"] * (r / r0)**-0.5 # [m/s]
			grid.V_t[pos] = INP_DICT["vt"] # [m/s]

			# Dust properties
			grid.T_d[pos] = grid.T_k[pos]
			grid.kapp_d[pos] = Type.KappLLaw("['0.25mm','10cm^2g^-1',-1]")

			if r <= rcen:
				# No molecules or dust present inside dust destruction radius,
				# which harbors an extremely optically thick continuum source
				grid.T_d[pos] = T_k(rcen) # K
				grid.kapp_d[pos] = Type.KappLLaw("['0.25mm','1e6cm^2g^-1',-1]")
				grid.X_mol[pos] = 0
		return

install_task(Task_W3OH1D("task_w3oh1d"))

################################################################################

class Task_ConstInfall1D(Task):
	"""
	Spherical caloud with constant infall, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Spherical caloud with constant infall, 1D version"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("xmol", Type.Fraction, None, "Molecular abundance"),
			Key("vin", Type.Velo, None, "Infall velocity"),
			Key("ndiv", Type.PosInt, 64, "Number of shells"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		# Model parameters (Tafalla et al. 2004)
		n0 = 6.0e4 * 1e6 # [cm^-3] -> [m^-3]
		r0 = 5e16 * 1e-2 / Unit.pc # [cm] -> [pc]
		radius = 2e17 * 1e-2 / Unit.pc # [cm] -> [pc]
		alpha = 4.0

		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		n = INP_DICT["ndiv"]

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius)

		# Setup model parameters
		grid.gas_to_dust = 0
		grid.T_cmb = 2.728 # K

		for pos in np.ndindex(n, 1, 1):
			r = grid.cen_i[pos]

			# Set physical parameters
			grid.n_H2[pos] = n0 / (1.0 + (r / r0)**alpha) # m^-3
			grid.T_k[pos] = 10.0 # K
			grid.X_mol[pos] = INP_DICT["xmol"] # Fraction
			grid.V_i[pos] = -INP_DICT["vin"] # m/s
			grid.V_t[pos] = 0.1e3 # [m/s]
		return

install_task(Task_ConstInfall1D("task_constinfall1d"))

################################################################################

class Task_ConstInfall3D(Task):
	"""
	Spherical caloud with constant infall, 3D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Spherical caloud with constant infall, 3D version"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("xmol", Type.Fraction, None, "Molecular abundance"),
			Key("vin", Type.Velo, None, "Infall velocity"),
			Key("ndiv", Type.PosInt, 64, "Number of shells"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		# Model parameters (Tafalla et al. 2004)
		n0 = 6.0e4 * 1e6 # [cm^-3] -> [m^-3]
		r0 = 5e16 * 1e-2 / Unit.pc # [cm] -> [pc]
		radius = 2e17 * 1e-2 / Unit.pc # [cm] -> [pc]
		length = radius * 2.0 # [pc]
		alpha = 4.0

		from sparx.grid import Grid_rec3d
		# Setup dimensions and gridding
		n = INP_DICT["ndiv"]

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_rec3d((n, n, n), (length, length, length))

		# Setup model parameters
		grid.gas_to_dust = 0
		grid.T_cmb = 2.728 # K

		for pos in np.ndindex(n, n, n):
			# Calculate radius at pos
			di = grid.cen_i[pos] - grid.cen[0]
			dj = grid.cen_j[pos] - grid.cen[1]
			dk = grid.cen_k[pos] - grid.cen[2]
			r = sqrt(di**2 + dj**2 + dk**2)

			# Set parameters only if within cloud radius
			if r <= radius:
				grid.n_H2[pos] = n0 / (1.0 + (r / r0)**alpha) # m^-3
				grid.T_k[pos] = 10.0 # K
				grid.X_mol[pos] = INP_DICT["xmol"] # Fraction
				grid.V_t[pos] = 0.1e3 # [m/s]

				# Project radial velocity onto Cartesian coordinates
				vel = phys.Vec3_Normalize([di, dj, dk])
				vel = phys.Vec3_Scale(vel, -INP_DICT["vin"]) # m/s
				grid.V_i[pos] = vel[0]
				grid.V_j[pos] = vel[1]
				grid.V_k[pos] = vel[2]
		return

install_task(Task_ConstInfall3D("task_constinfall3d"))
################################################################################



class Task_Hollow1D(Task):
	"""
	Spherical caloud with constant infall, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Spherical caloud with constant infall, 1D version"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("n_max", Type.NumDens, None, "Gas density"),
			Key("Vin", Type.Velo, None, "Infall velocity"),
			Key("ndiv", Type.PosInt, 64, "Number of shells"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		radius_out = 0.1 # [pc] cloud radius
		radius_in = 0.0 # [pc] cloud radius
		n = INP_DICT["ndiv"]
		dr =(radius_out-radius_in)/n
		v_max=2*INP_DICT["Vin"]
		from numpy import zeros
		manual_grid=zeros(n+1)
		for pos in range(0,n+1):
			manual_grid[pos] = radius_in+pos*dr

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius_out)

		# Setup model parameters
		grid.gas_to_dust = 0
		grid.T_cmb = 2.73 # K
		for pos in np.ndindex(n, 1, 1):
			# Manually regrid if requested
			id = pos[0]
			if id == 0:
				grid.SetZoneBoundary_sph1d(id, manual_grid[id], manual_grid[id+1])
				assert manual_grid[id] < manual_grid[id+1]

			# Set physical parameters
			grid.n_H2[pos] = 1e11 # m^-3
			grid.T_k[pos] = 10.0 # K
			if grid.cen_i[pos] <= 0.5*radius_out:
				grid.X_mol[pos] = 1e-9 # Fraction
			else:
				grid.X_mol[pos] = 1e-8 # Fraction
			grid.V_i[pos] = -v_max * (grid.cen_i[pos]/0.01)**-0.5 # m/s
			grid.V_t[pos] = 200 # [m/s]
		return

install_task(Task_Hollow1D("task_hollow1d"))

################################################################################


class Task_Hollow3D(Task):
	"""
	Spherical caloud with constant infall, 3D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Spherical caloud with constant infall, 3D version"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("n_max", Type.NumDens, None, "Gas density"),
			Key("Vin", Type.Velo, None, "Infall velocity"),
			Key("ndiv", Type.PosInt, 64, "Number of shells"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		# Model parameters (Tafalla et al. 2004)
		r0 = 0.01 # [cm] -> [pc]
		radius = 0.1  # [cm] -> [pc]
		length = radius * 2.0 # [pc]
		v_max=2*INP_DICT["Vin"]

		from sparx.grid import Grid_rec3d
		# Setup dimensions and gridding
		n = INP_DICT["ndiv"]

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_rec3d((n, n, n), (length, length, length))

		# Setup model parameters
		grid.gas_to_dust = 0
		grid.T_cmb = 2.73 # K

		for pos in np.ndindex(n, n, n):
			# Calculate radius at pos
			di = grid.cen_i[pos] - grid.cen[0]
			dj = grid.cen_j[pos] - grid.cen[1]
			dk = grid.cen_k[pos] - grid.cen[2]
			r = sqrt(di**2 + dj**2 + dk**2)

			# Set parameters only if within cloud radius


			grid.n_H2[pos] = 1e11 # m^-3
			grid.T_k[pos] = 10.0 # K
			if r <= 0.5*radius:
				grid.X_mol[pos] = 1e-9 # Fraction
			else:
				grid.X_mol[pos] = 1e-8 # Fraction
			grid.V_t[pos] = 200 # [m/s]
			# Project radial velocity onto Cartesian coordinates
			vel = phys.Vec3_Normalize([di, dj, dk])
			vel = phys.Vec3_Scale(vel, -v_max * (r/0.01)**-0.5) # m/s
			grid.V_i[pos] = vel[0]
			grid.V_j[pos] = vel[1]
			grid.V_k[pos] = vel[2]
		return

install_task(Task_Hollow3D("task_hollow3d"))


################################################################################


class Task_Shu1D(Task):
	"""
	Spherical cloud with constant infall, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Shu's 1-D spherical collapsing cloud"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("n_max", Type.NumDens, '1e12m^-3', "Gas density"),
			Key("Vin", Type.Velo, '0.1kms^-1', "Infall velocity"),
			Key("ndiv", Type.PosInt, 64, "Number of shells"),
			Key("abundance", Type.Fraction, None, "Molecular abundance")
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		radius_out = 0.1 # [pc] cloud radius
		radius_in = 0.0 # [pc] cloud radius
		n = INP_DICT["ndiv"]
		dr =(radius_out-radius_in)/n
		v_max=2*INP_DICT["Vin"]
		r0=0.01
		from numpy import zeros
		manual_grid=zeros(n+1)
		for pos in range(0,n+1):
			manual_grid[pos] = radius_in+pos*dr

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius_out)

		# Setup model parameters
		grid.gas_to_dust = 0
		grid.T_cmb = 2.73 # K
		xmol=INP_DICT["abundance"]
			
		for pos in np.ndindex(n, 1, 1):
			# Manually regrid if requested
			id = pos[0]
			if id == 0:
				grid.SetZoneBoundary_sph1d(id, manual_grid[id], manual_grid[id+1])
				assert manual_grid[id] < manual_grid[id+1]
			# Set physical parameters
			grid.n_H2[pos] = INP_DICT["n_max"]*(grid.cen_i[pos]/r0)**-1.5 # m^-3
			grid.T_k[pos] = 10.0 # K
			grid.X_mol[pos] = xmol # Fraction
			grid.V_t[pos] = 200. # [m/s]
			grid.V_i[pos] = -v_max * (grid.cen_i[pos]/r0)**-0.5 # m/s
		return

install_task(Task_Shu1D("task_shu1d"))

################################################################################

class Task_Shu3D(Task):
	"""
	Spherical caloud with constant infall, 3D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Spherical caloud with constant infall, 3D version"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("n_max", Type.NumDens, None, "Gas density"),
			Key("Vin", Type.Velo, None, "Infall velocity"),
			Key("ndiv", Type.PosInt, 64, "Number of shells"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		# Model parameters (Tafalla et al. 2004)
		r0 = 0.01 # [cm] -> [pc]
		radius = 0.1  # [cm] -> [pc]
		length = radius * 2.0 # [pc]
		v_max=2*INP_DICT["Vin"]

		from sparx.grid import Grid_rec3d
		# Setup dimensions and gridding
		n = INP_DICT["ndiv"]

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_rec3d((n, n, n), (length, length, length))

		# Setup model parameters
		grid.gas_to_dust = 0
		grid.T_cmb = 2.73 # K

		for pos in np.ndindex(n, n, n):
			# Calculate radius at pos
			di = grid.cen_i[pos] - grid.cen[0]
			dj = grid.cen_j[pos] - grid.cen[1]
			dk = grid.cen_k[pos] - grid.cen[2]
			r = sqrt(di**2 + dj**2 + dk**2)

			# Set parameters only if within cloud radius
			if r <= radius:
				grid.n_H2[pos] = INP_DICT["n_max"]*(r/r0)**-1.5 # m^-3
			else:
				grid.n_H2[pos] = 0
			grid.T_k[pos] = 10.0 # K
			grid.X_mol[pos] = 1e-9 # Fraction
			grid.V_t[pos] = 200 # [m/s]

			# Project radial velocity onto Cartesian coordinates
			vel = phys.Vec3_Normalize([di, dj, dk])
			vel = phys.Vec3_Scale(vel, -v_max * (r/0.01)**-0.5) # m/s
			grid.V_i[pos] = vel[0]
			grid.V_j[pos] = vel[1]
			grid.V_k[pos] = vel[2]
		return

install_task(Task_Shu3D("task_shu3d"))

################################################################################

class Task_Toy1D(Task):
	"""
	Toy, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Toy, 1D version"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("xmol", Type.Fraction, None, "Molecular abundance"),
			Key("ndiv", Type.PosInt, 64, "Number of shells"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		n0 = 1e6 * 1e6 # [cm^-3] -> [m^-3]
		radius = 0.1  # [cm] -> [pc]
 
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		n = INP_DICT["ndiv"]

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius)

		# Setup model parameters
		grid.gas_to_dust = 0
		grid.T_cmb = 2.728 # K	

		for pos in np.ndindex(n, 1, 1):
			r = grid.cen_i[pos]

			# Set physical parameters
			grid.n_H2[pos] = n0 # m^-3
			grid.T_k[pos] = 50.0 # K
			grid.X_mol[pos] = INP_DICT["xmol"] # Fraction
			grid.V_i[pos] = 0.0 # m/s
			grid.V_t[pos] = 400.0 # [m/s]
		return

install_task(Task_Toy1D("task_toy1d"))




################################################################################
class Task_Toy21D(Task):
	"""
	Toy, 1D version
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Toy, 1D version"

		# Explanation
		self.expl = ""

		# Keys
		self.keys = [
			Key("out", Type.NewFile, None, "Output model file"),
			Key("xmol", Type.Fraction, None, "Molecular abundance"),
			Key("ndiv", Type.PosInt, 64, "Number of shells"),
		]

		# C function to call
		self.cfunc = _sparx.task_pygrid

	##
	## Task procedures
	##
	def main(self):
		n0 = 1e6 * 1e6 # [cm^-3] -> [m^-3]
		radius = 0.1  # [cm] -> [pc]
 
		from sparx.grid import Grid_sph1d
		# Setup dimensions and gridding
		n = INP_DICT["ndiv"]

		# Generate grid and attach to input
		INP_DICT["pygrid"] = grid = Grid_sph1d(n, radius)

		# Setup model parameters
		grid.gas_to_dust = 0
		grid.T_cmb = 2.728 # K	

		for pos in np.ndindex(n, 1, 1):
			r = grid.cen_i[pos]

			# Set physical parameters
			grid.n_H2[pos] = n0 # m^-3
			grid.T_k[pos] = 50.0*(r/radius)**-0.5 # K
			grid.X_mol[pos] = INP_DICT["xmol"] # Fraction
			grid.V_i[pos] = 0.0 # m/s
			grid.V_t[pos] = 400.0 # [m/s]
		return

install_task(Task_Toy21D("task_toy21d"))




################################################################################
class Task_AMC(Task):
	"""
	Accelerated Monte Carlo method for solving detailed balance
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Non-LTE molecular excitation solver"

		# Explanation
		self.expl = "Calculate non-LTE molecular excitation by solving detailed balance with the Accelerated Monte Carlo method (Hogerheijde & van der Tak 2000)"

		# Keys
		self.keys = [
			Key("source", Type.OldFile, None, "Name of input source model (HDF5 file)"),
			Key("pops", Type.OldFile, Type.Optional, "Name of initial population file (HDF5 file)"),
			Key("out", Type.NewFile, None, "Name of output file (HDF5 file)"),
			Key("molec", Type.Molec, None, "Molecule to calculate"),
			Key("overlap", Type.Velo, '0kms^-1', "overlapping calculation (for hyperfine splitting)"),
			Key("lte", Type.Bool, "False", "Whether to start convergence from LTE conditions"),
			Key("trace", Type.Bool, "False", "Whether to trace convergence history"),
			Key("tolerance", Type.Fraction, "5e-3", "Convergence criterion for fixed rays stage"),
			Key("snr", Type.Float, "20", "Upper limit of Monte Carlo noise level"),
			Key("minpop", Type.Fraction, "1e-6", "Minimum pops to test for convergence"),
			Key("nrays", Type.PosInt, "1000", "Number of initial rays per zone"),
			Key("maxiter", Type.PosInt, "1000", "Maximum number of iterations for solving detailed balance"),
			# esc 09Sep29: It seems it would be best for raniter to be set to at least 5 to
			#              prevent false convergence (just an empirical guess)
			Key("fixiter", Type.PosInt, "5", "Minimum number of iterations for fixed rays stage"),
			Key("raniter", Type.PosInt, "5", "Minimum number of iterations for random rays stage"),
		]

		# C function to call
		self.cfunc = _sparx.task_amc
			##
	## Task procedures
	##
	def main(self):
		class amc:
			if INP_DICT["pops"]==None:
				popsold=0
			else:
				popsold=1
			
		        
		        overlap_vel = INP_DICT["overlap"]
			if ( overlap_vel == 0.0):
			        overlap_int=0
			else:
			        overlap_int=1
		INP_DICT["amc"] = amc
		return

install_task(Task_AMC("task_amc"))

##
## telsim related tasks:
## task_contobs and task_lineobs
##
# Keys both contobs and lineobs share:
telsim_keys = [
	Key("dist", Type.Length, "1kpc", "Distance to source"),
	Key("cell", Type.Custom([Type.Angle, Type.Angle]), "['1asec', '1asec']", "Angular size of each pixel"),
	Key("npix", Type.Custom([Type.PosInt, Type.PosInt]), "[128, 128]", "Image dimensions in number of pixels"),
	Key("chan", Type.Custom([Type.PosInt, Type.Velo]), "[64, '0.1kms^-1']", "Number of spectral channels and width of each channel (in velocity units)"),
	Key("unit", Type.Option(['JY/PIXEL', 'K']), "JY/PIXEL", "Image brightness unit"),
	Key("rotate", Type.Custom([Type.Angle, Type.Angle, Type.Angle]), "['0deg', '0deg', '0deg']", "Rotation of model about its x, y and z axes, in x-y-z order."),
	Key("tau", Type.NewFile, Type.Optional, "Name of output tau cube (Miriad image dataset)"),
	Key("subres", Type.Custom([[Type.Angle, Type.Angle, Type.Angle, Type.Angle, Type.PosInt]]), Type.Optional, "Boxed regions for sub-resolution averaging. Meaning of values are [[blc_x, blc_y, trc_x, trc_y, nsub], ...]")	
]

##
## task_contobs
##
class Task_ContObs(Task):
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Synthetic continuum map generator"

		# Explanation
		self.expl = "Generates synthetic continuum images"

		# Keys
		self.keys = [
			Key("source", Type.OldFile, None, "Input source model (HDF5 file), must contain dust information"),
			Key("out", Type.NewFile, None, "Name of output image (Miriad image dataset)"),
			Key("wavelen", Type.Length, None, "Wavelength of observation")
		]+telsim_keys

		# C function to call
		self.cfunc = _sparx.task_telsim

	##
	## Task procedures
	##
	def main(self):
		# Create 'obs' class used by C code (consider changing this
		# eventually)
		class obs:
                        # Enable continuum mode 
			cont = 1 
			wavelen = INP_DICT['wavelen']
			coldens=0
		INP_DICT["obs"] = obs
		return

install_task(Task_ContObs("task_contobs"))

##
## task_lineobs
##
class Task_LineObs(Task):
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Synthetic line map generator"

		# Explanation
		self.expl = "Generates synthetic line images"

		# Keys
		self.keys = [
			Key("source", Type.OldFile, None, "Input source model (HDF5 file), must contain dust information"),
			Key("out", Type.NewFile, None, "Name of output image (Miriad image dataset)"),
			Key("line", Type.Index, None, "Line index of observing line"),
			Key("overlap", Type.Velo, '0kms^-1', "line overlaping calculation"),
			Key("lte", Type.Bool, "False", "Init model to LTE pops of Molec"),
			Key("molec", Type.Molec, Type.Optional, "Molecule to calculate"),
			Key("excit", Type.Bool, "False", "Excitation temperature"),
			Key("zeeman", Type.Bool, "False", "Excitation temperature")
		]+telsim_keys

		# C function to call
		self.cfunc = _sparx.task_telsim

	##
	## Task procedures
	##
	def main(self):
		# Create 'obs' class used by C code (consider changing this
		# eventually)
		class obs:
			# Disable continuum and coldens mode -- line observations
			cont = 0 
			coldens=0
			line = INP_DICT["line"]
			overlap_vel = INP_DICT["overlap"]
			if ( overlap_vel == 0.0):
			        overlap_int=0
			else:
			        overlap_int=1

		INP_DICT["obs"] = obs
		return

install_task(Task_LineObs("task_lineobs"))


##
## task_coldens
##
class Task_ColDens(Task):
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Synthetic column density map generator"

		# Explanation
		self.expl = "Generates synthetic column density images"

		# Keys
		self.keys = [
			Key("source", Type.OldFile, None, "Input source model (HDF5 file), must contain dust information"),
			Key("out", Type.NewFile, None, "Name of output image (Miriad image dataset)")
		]+telsim_keys

		# C function to call
		self.cfunc = _sparx.task_telsim

	##
	## Task procedures
	##
	def main(self):
		# Create 'obs' class used by C code (consider changing this
		# eventually)
		class obs:
			# Enable coldens mode 
			cont = 0 
			coldens=1
		INP_DICT["obs"] = obs
		
		
		return

install_task(Task_ColDens("task_coldens"))



############################################
# Line Fitting Task (Dec 2013 added by I-Ta)
############################################
class Task_LineFitting(Task):
	"""
	Line Fitting module 
	"""
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Abundance fitting for a target spectrum from LTE/non-LTE line calculation"

		# Explanation
		self.expl = "Simulated Annealing approach to minimize CHI-square between an analytical model and desired spectrum"

		# Keys
		self.keys = [
			Key("target", Type.OldFile, None, "File to approach through the fitting"),
			Key("kmax", Type.PosInt, 100, "SA iteration number"),
			Key("X0", Type.Fraction, 1e-9, "The initial guess for molecular abundance"),
			Key("molecule", Type.Molec, None, "Molecule to fit"),
			Key("StopNR", Type.Fraction, 1e-2, "Stopping noise ratio"),
		]
	##
	## Task procedures
	##
	
	

	def main(self):
		global target,molecule
		
		target=INP_DICT["target"]
		kmax=INP_DICT["kmax"]
		X0=INP_DICT["X0"]
		molecule=INP_DICT["molecule"]
		StopNR=INP_DICT["StopNR"]
		
		def ProduceXmol(iteration,BestX,FurnaceTemperature):
			global kmax
			import random
			from math import tan,pi,exp
			
			RandomFloat = random.uniform(-1.0, 1.0) 
			RandomDistribution = FurnaceTemperature * tan( 0.5 * pi * RandomFloat )
			factor = exp( RandomDistribution )
			
			Xmol = BestX * factor
			if ( Xmol > 1.0):
				Xmol=ProduceXmol(iteration,BestX,FurnaceTemperature)
			return Xmol
		
		
		def SimulatedAnnealing(StopNR,kmax,X0):
			global target,molecule
			global targetFrequency,targetFlux
			from os import remove
			from os.path import exists
			
			#from sparx.plotter import GUIPlotter, mainloop
			#pltr = GUIPlotter()	
			
			if exists('fitting_log.dat'):
				remove('fitting_log.dat')
	
			# load target spectrum
			targetFrequency=[]
			targetFlux=[]
			targetfile = open(target, 'r')
			for line in targetfile:
				column = line.split()
				targetFrequency.append(float(column[1]))
				targetFlux.append(float(column[2]))
			targetfile.close()
			# calculate the signal intensity (to compare to the noise)
			signal=0.0
			for i in range(len(targetFrequency)):
				signal += targetFlux[i]
			signal /= len(targetFrequency)
			

			# iteration variables

			#IterArray=[]
			#ErrorArray=[]
			
			MinimumError=float("inf")
			BestX=X0
			
			
			for iteration in range(kmax):
				# initialize the abundance
				if (iteration == 0):
					Xmol=BestX
				else:
					# the Furnace temperature represen the range of the sampling abundance, and it's related to the noise to signal ratio
					FurnaceTemperature = NoiseRatio 
					Xmol=ProduceXmol(iteration,BestX,FurnaceTemperature)
					
				# call the physical simulation routine
				FSAE=ModelToSpectrumProcess(Xmol)
				  
				# update the best result if it is
				if (FSAE<MinimumError):
					BestX=Xmol
					MinimumError=FSAE
					NoiseRatio = MinimumError/signal
					# jump out/stop the iterative loop when achieve deserved StopNR
					if ( NoiseRatio < StopNR ):
						break
				# save the log file
				LOGfile = open('fitting_log.dat', 'a')
				print >>LOGfile,"%(0)d %(1)e %(2)e %(3)e %(4)e" %{ '0':iteration, '1':FSAE/signal, '2':NoiseRatio, '3':Xmol,'4':BestX}
				LOGfile.close()
				#IterArray.append(iteration)
				#ErrorArray.append(error)	
	
				#pltr.plot(
				#	IterArray,
				#	ErrorArray,
				#	name='ERROR',
				#	xlab='Iteration',
				#	ylab='deviation',
				#	logx=False,
				#	logy=True,
				#	beg=None,
				#	end=None,
				#	lsty="-",
				#	msty="o"
				#)
		
				#pltr.show()
				#mainloop()

		def ModelToSpectrumProcess(Xmol):
			global molecule
			global targetFrequency,targetFlux
			import subprocess
			from os import remove
			from os.path import exists
			import shutil
			
			fwhm=10
			freq=89.9
			
			if exists('model'):
				remove('model')
			if exists('pops'):
				remove('pops')
			if exists('mapJY/'):
				shutil.rmtree('mapJY/')
			if exists('subJY/'):
				shutil.rmtree('subJY/')
			if exists('contJY/'):
				shutil.rmtree('contJY/')
			if exists('conv/'):
				shutil.rmtree('conv/')
			if exists('convK/'):
				shutil.rmtree('convK/')
			
			# generate physical model
			TASK_DICT["task_shu1d"].run([
				"out=model",
				"abundance="+str(Xmol),
				"ndiv=16"
			])
			# compute non-LTE exitation
			#TASK_DICT["task_amc"].run([
			#	"source=model",
			#	"out=pops",
			#	"molec="+molecule,
			#	"lte=True",
			#	"tolerance=5e-3"
			#])		
			# compute line image
			TASK_DICT["task_lineobs"].run([
				#"source=pops",
				"source=model",
				"lte=True",
				"molec="+molecule,
				"out=mapJY",
				"line=0",
				"dist=100pc",			
				"subres=[['-20asec', '-20asec', '20asec', '20asec', 4]]",
				"chan=[100, '0.02kms^-1']"
			])
			# subtract continuum emission 
			subprocess.call(
				"contsub in=mapJY out=subJY cont=contJY contchan='(1,5),(96,100)' mode=mean",
				shell=True
			)
			#convolve image with beam size
			subprocess.call(
				"convol map=subJY fwhm="+str(fwhm)+" out=conv mode=mean",
				shell=True
			)
			# convert to the unit of Kelvin from Jansky
			subprocess.call(
				"maths exp=\"1.22e6*<conv>/(("+str(freq)+"**2)*("+str(fwhm)+"**2))\" out=convK",
				shell=True
			)
			# write the convolved spectrum into ASCII
			subprocess.call(
				"imspect in=convK log=temp region='rel,box(0,0,0,0)'",
				shell=True
			)
			subprocess.call(
				"sed '1,4d' temp > spectrum_model.dat",
				shell=True
			)
			# compute the error between the target spectrum
			modelFrequency=[]
			modelFlux=[]
			modelfile = open('spectrum_model.dat', 'r')
			for line in modelfile:
				column = line.split()
				modelFrequency.append(float(column[1]))
				modelFlux.append(float(column[2]))
			modelfile.close()
			error=0.0
			for i in range(len(modelFrequency)):
				if ( modelFrequency[i] == targetFrequency[i] ):
					diff = modelFlux[i]-targetFlux[i]
					error += diff * diff
			error /= len(modelFrequency)
			error = sqrt(error)
			
			
			
			return error
		
		SimulatedAnnealing(StopNR,kmax,X0)
		
		return

install_task(Task_LineFitting("task_linefitting"))
####################################
####################################



##
## Tasks for testing purposes only
##

# test_fft
install_task(Task("test_fft", "Test fast fourier transform functionality.", "", _sparx.test_fft))
install_task(Task("test_planck", "Test Planck function.", "", _sparx.test_planck))









