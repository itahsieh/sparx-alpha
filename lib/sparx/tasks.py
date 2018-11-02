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
			Key("overlap", Type.Velo, '0kms^-1', "overlapping calculation (for hyperfine splitting)"),
			Key("lte", Type.Bool, "True", "Whether to start convergence from LTE conditions"),
			Key("trace", Type.Integer, 0, "Write out the temporary result in every n stetp of iteration"),
			#Key("tolerance", Type.Fraction, "5e-3", "Convergence criterion for fixed rays stage"),
			Key("snr", Type.Float, "20", "Upper limit of Monte Carlo noise level"),
			Key("minpop", Type.Fraction, "1e-6", "Minimum pops to test for convergence"),
			Key("nrays", Type.PosInt, "1000", "Number of initial rays per zone"),
			Key("maxiter", Type.PosInt, "1000", "Maximum number of iterations for solving detailed balance"),
			# esc 09Sep29: It seems it would be best for raniter to be set to at least 5 to
			#              prevent false convergence (just an empirical guess)
			Key("fixiter", Type.PosInt, "5", "Minimum number of iterations for fixed rays stage"),
			Key("raniter", Type.PosInt, "5", "Minimum number of iterations for random rays stage"),
			Key("qmc", Type.Bool, "True", "Quasi-Monte-Carlo method"),
			Key("ali", Type.Bool, "False", "Lambda iteration only / Three-staged Monte-Carlo convergent automation"),
			Key("dat", Type.Bool, "False", "1-D level populations ascii file ouput"),
			Key("sor", Type.Float, "1.0", "successive and over-relaxation method"),
		]

		# C function to call
		self.cfunc = _sparx.task_amc
			##
	## Task procedures
	##
	def main(self):
                class obs: 
                        task = 'amc'
		return

install_task(Task_AMC("task_amc"))




################################################
#    The following tasks are preprocessors     #
################################################
# Keys all tasks share:


# postprocessing keys are for all post tasks, including telsim, contribution, and model2vtk
postprocess_keys = [
    Key("source", Type.OldFile, None, 
        "Input source model (HDF5 file), must contain dust information"
        )
]

# observer keys are for all post tasks except model2vtk
observer_keys = [
    Key("dist", Type.Length, "1kpc", 
        "Distance to source"
        ),
    Key("rotate", Type.Custom([Type.Angle, Type.Angle, Type.Angle]), 
        "['0deg', '0deg', '0deg']", 
        "Rotation of model about its x, y and z axes, in x-y-z order."
        )
]

# for lineobs, zeeman, contobs, coldens
telsim_keys = [
    Key("out", Type.NewFile, None, 
        "Name of output image (FITS image dataset)"
        ),
    Key("dist", Type.Length, "1kpc", 
        "Distance to source"
        ),
    Key("cell", Type.Custom([Type.Angle, Type.Angle]), "['1asec', '1asec']", 
        "Angular size of each pixel"
        ),
    Key("npix", Type.Custom([Type.PosInt, Type.PosInt]), "[128, 128]", 
        "Image dimensions in number of pixels"
        ),
    Key("rotate", Type.Custom([Type.Angle, Type.Angle, Type.Angle]), 
        "['0deg', '0deg', '0deg']", 
        "Rotation of model about its x, y and z axes, in x-y-z order."
        ),
    Key("subres", 
        Type.Custom([[Type.Angle, Type.Angle, Type.Angle, Type.Angle, Type.PosInt]]), 
        Type.Optional, 
        "Boxed regions for sub-resolution averaging. Meaning of values are [[blc_x, blc_y, trc_x, trc_y, nsub], ...]"
        )
]

# for lineobs, contobs, zeeman
radiation_keys = [
    Key("unit", Type.Option(['JY/PIXEL', 'K']), "JY/PIXEL", 
        "Image brightness unit"
        ),
    Key("tau", Type.NewFile, Type.Optional, 
        "Name of output tau cube (Miriad image dataset)"
        )
]
    
# for linectb, contctb, zeemanctb
contribution_keys = [
    Key("unit", Type.Option(['JY/PC', 'K/PC']), "K/PC", 
        "Intensity contribution per length"
        )
]

# for linectb, contctb, and vtk
vtk_keys = [
    Key("slice", Type.Bool, "False", 
        "Slice cut at X-Y plane, only used in sph1d model"
        ),
    Key("out", Type.NewFile, None, 
        "Name of output VTK file "
        )
]

# for lineobs, zeeman, and linectb
line_keys = [
    Key("line", Type.Index, None, 
            "Line index of observing line"
            ),
    Key("overlap", Type.Velo, '0kms^-1', 
            "line overlaping calculation"
            ),
    Key("chan", Type.Custom([Type.PosInt, Type.Velo]), "[64, '0.1kms^-1']", 
            "Number of spectral channels and width of each channel (in velocity units)"
            ),
    Key("lte", Type.Bool, "False", 
            "Init model to LTE pops of Molec"
            )
]

# for contobs, contctb
cont_keys = [
    Key("chan", Type.Custom([Type.PosInt, Type.Velo]), "[1, '0.1kms^-1']",
        "Number of spectral channels and width of each channel (in velocity units)"
        ),
    Key("wavelen", Type.Length, None, 
        "Wavelength of observation"
        )
]

# for coldens
coldens_keys = [
    Key("chan", Type.Custom([Type.PosInt, Type.Velo]), "[1, '0.1kms^-1']", 
        "Number of spectral channels and width of each channel (in velocity units)"
        ),
    Key("unit", Type.Option(['MKS', 'CGS']), "MKS", 
        "Image brightness unit"
        ),
    Key("tracer", Type.Bool, "False", 
        "Molecular tracer"
        )
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
		self.keys = postprocess_keys + observer_keys + telsim_keys + radiation_keys + cont_keys

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
			task = 'cont' 
			wavelen = INP_DICT['wavelen']
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
		self.keys = postprocess_keys + observer_keys + telsim_keys + radiation_keys + line_keys

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
			task='line'
			line = INP_DICT["line"]
			overlap_vel = INP_DICT["overlap"]

		INP_DICT["obs"] = obs
		return

install_task(Task_LineObs("task_lineobs"))



##
## task_zeeman
##
class Task_Zeeman(Task):
        ##
        ## Task configuration
        ##
        def configure(self):
                # Name: defined when task is registered

                # Description
                self.desc = "Zeeman effect : Stokes-V image"

                # Explanation
                self.expl = "Zeeman effect polarization to see B-strength along the line of sight"

                # Keys
                self.keys = postprocess_keys + observer_keys + telsim_keys + radiation_keys + line_keys

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
                        task='zeeman'
                        line = INP_DICT["line"]


                INP_DICT["obs"] = obs
                return

install_task(Task_Zeeman("task_zeeman"))


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
                self.keys = postprocess_keys + observer_keys + telsim_keys + coldens_keys

		# C function to call
		self.cfunc = _sparx.task_coldens

	##
	## Task procedures
	##
	def main(self):
		# Create 'obs' class used by C code (consider changing this
		# eventually)
		class obs:
			# Enable coldens mode 
			task='coldens'
		INP_DICT["obs"] = obs
		
		
		return

install_task(Task_ColDens("task_coldens"))



##
## task_linectb
##
class Task_LineCtb(Task):
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Compute the contribution function and excitation temperature of the molecular line and write out VTK file for visualization"

		# Explanation
		self.expl = "Visualizing the distribution of intensity-contribution could help to analyze the region where the molecule absorbs or emits radiation "

		# Keys
		self.keys = postprocess_keys + observer_keys + vtk_keys + contribution_keys + line_keys

		# C function to call
		self.cfunc = _sparx.task_visual

	##
	## Task procedures
	##
	def main(self):
		# Create 'obs' class used by C code (consider changing this
		# eventually)
		class obs:
			# Disable continuum and coldens mode -- line observations
			task='linectb'
			line = INP_DICT["line"]
			overlap_vel = INP_DICT["overlap"]

		INP_DICT["obs"] = obs
		return

install_task(Task_LineCtb("task_linectb"))


##
## task_contctb
##
class Task_ContCtb(Task):
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "Compute the contribution function of continnum emission and write out VTK file for visualization"

		# Explanation
		self.expl = "Visualizing the distribution of intensity-contribution could help to analyze the region where the dust/free-free absorbs or emits radiation "

		# Keys
		self.keys = postprocess_keys + observer_keys + vtk_keys + contribution_keys + cont_keys

		# C function to call
		self.cfunc = _sparx.task_visual

	##
	## Task procedures
	##
	def main(self):
		# Create 'obs' class used by C code (consider changing this
		# eventually)
		class obs:
			# Disable continuum and coldens mode -- line observations
			task='contctb'
			wavelen = INP_DICT['wavelen']

		INP_DICT["obs"] = obs
		return

install_task(Task_ContCtb("task_contctb"))

##
## task_model2vtk
##
class Task_Model2VTK(Task):
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "write out VTK file to visualize physical properties, including density, temperature, velocity field, B-field"

		# Explanation
		self.expl = "converting the model to VTK format"

		# Keys
		self.keys = postprocess_keys + vtk_keys

		# C function to call
		self.cfunc = _sparx.task_visual

	##
	## Task procedures
	##
	def main(self):
		# Create 'obs' class used by C code (consider changing this
		# eventually)
		class obs:
			# Disable continuum and coldens mode -- line observations
			task='model2vtk'

		INP_DICT["obs"] = obs
		return

install_task(Task_Model2VTK("task_model2vtk"))


##
## task_pops2ascii
##
class Task_Pops2ASCII(Task):
	##
	## Task configuration
	##
	def configure(self):
		# Name: defined when task is registered

		# Description
		self.desc = "write out level population to ASCII data file"

		# Explanation
		self.expl = "converting HDF5 population data to ASCII file"

		# Keys
		self.keys = postprocess_keys + [
                    Key("out", Type.NewFile, None, 
                        "Name of output file (.dat)"
                        )
                    ]

		# C function to call
		self.cfunc = _sparx.task_pops2ascii


install_task(Task_Pops2ASCII("task_pops2ascii"))

































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
		
		def ProduceXmol(BestX,FurnaceTemperature):
			import random
			from math import tan,pi,exp
			
			RandomFloat = random.uniform(-1.0, 1.0) 
			RandomDistribution = FurnaceTemperature * tan( 0.5 * pi * RandomFloat )
			factor = exp( RandomDistribution )
			
			Xmol = BestX * factor
			if ( Xmol > 1.0):
				Xmol=ProduceXmol(BestX,FurnaceTemperature)
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
					Xmol=ProduceXmol(BestX,FurnaceTemperature)
					
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









