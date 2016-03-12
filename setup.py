#! /bin/env python
# Python Distutils setup script for SPARX

# Version of SPARX
Version = '2.3.1'
dev = 1
# MIRIAD support option
MIRSUPPORT = 0  
# number of Thread using in per job
NumberOfThread = 24 

# Define macros
macros = [
        ('NTHREAD', str(NumberOfThread)),
        ('MIRSUPPORT',MIRSUPPORT),
        ('DEV_VERSION',dev)
]


# Test for MPI by checking whether mpicc can be called
from subprocess import call, Popen, PIPE
HAVE_MPI = (call("mpicc src/mpi-test.c -o/tmp/a.out", shell=True, stdout=PIPE, stderr=PIPE) == 0)
# Get svn revision and update VERSION
import time
p = Popen("svnversion", shell=True, stdout=PIPE)
REV = p.communicate()[0].strip()
fo = file("lib/sparx/VERSION", "w")
fo.write("%s (r%s, %s)"%(Version,REV, time.asctime()))
fo.close()

##
## Gather information for setting up the package
##                              
# Some necessary imports
import os, glob
from os.path import exists, realpath, expanduser

# Get Python paths
import sys
PYINC = sys.prefix+"/include"
PYLIB = sys.prefix+"/lib"
if not (exists(PYINC) and exists(PYLIB)):
	raise Exception, "Cannot locate Python include and lib directories"

# Get NumPy paths
import numpy
NPINC = numpy.get_include()

# Get Miriad paths
if MIRSUPPORT:
	MIRINC = os.getenv("MIRINC")
	MIRLIB = os.getenv("MIRLIB")
	if not (MIRINC and MIRLIB):
		raise Exception, "MIRIAD environment variables not present, cannot locate Miriad headers or libraries"
	MIRINC1 = realpath(MIRINC+"/../pgplot-miriad-remix")
	MIRINC2 = realpath(MIRINC+"/../miriad-c")
	if not (exists(MIRINC1) and exists(MIRINC2)):
		raise Exception, "MIRIAD include paths '%s' and '%s' not present, cannot continue"%(MIRINC1, MIRINC2)





# Check for additional search paths specified by user
USER_INCLUDE = []
USER_LIB = []
args = sys.argv[:]
mpi_libs = []
for arg in args:
	if arg.find('--with-include=') == 0:
		USER_INCLUDE += [expanduser(arg.split('=')[1])]
		sys.argv.remove(arg)
	elif arg.find('--with-lib=') == 0:
		USER_LIB += [expanduser(arg.split('=')[1])]
		sys.argv.remove(arg)
	elif arg.find('--no-mpi') == 0:
		HAVE_MPI = False
		sys.argv.remove(arg)
	elif arg.find('--lam') == 0:
		mpi_libs = ['lammpio', 'mpi', 'lam']
		sys.argv.remove(arg)
	elif arg.find('--mpich') == 0:
		mpi_libs = ['mpich', 'pgc', 'pgftnrtl', 'pgftnrtl', 'nspgc', 'pgc', 'rt']
		sys.argv.remove(arg)


if not HAVE_MPI:
	print\
'''



NO MPI SUPPORT!



'''
else:
	print\
'''



MPI support available



'''

# Compiler flags
compiler_flags = [
	'-std=c99',
	'-pedantic',
	'-fshort-enums',
	'-fno-common',
	'-Dinline=',
	'-g',
	'-rdynamic',
	'-O3',
	'-pthread',
#	'-Werror',
	'-Wall',
	'-W',
	'-Wmissing-prototypes',
	'-Wstrict-prototypes',
	'-Wpointer-arith',
	'-Wcast-qual',
	'-Wcast-align',
	'-Wwrite-strings',
	'-Wnested-externs',
]



# Header directories
header_dirs = [
	'src',
	PYINC,
	NPINC,
]+USER_INCLUDE

if MIRSUPPORT:
	header_dirs += [ MIRINC1, MIRINC2]


# Library directories
lib_dirs = [
	PYLIB,
]+USER_LIB

if MIRSUPPORT:
	lib_dirs += [
		MIRLIB
		]

# Libraries to link to
libs = [
	'X11',
	'm',
	'gsl',
	'gslcblas',
	'fftw3',
	'hdf5',
	'hdf5_hl',
	'cfitsio'
]

if MIRSUPPORT:
	libs += [
	'cpgplot',
	'pgplot',
	'mir',
	'mir_uvio',
	'mir_linpack'
	]

# Base source files
sources_base = [
	'src/data_structs.c',
	'src/debug.c',
	'src/error.c',
	'src/geometry.c',
	'src/kappa.c',
	'src/memory.c',
	'src/molec.c',
	'src/numerical.c',
	'src/physics.c',
	'src/python-wrappers.c',
	'src/zone.c',
	'src/zone-hdf5.c',
	'src/fits-and-miriad-wrappers.c',
]

if MIRSUPPORT:
	sources_base += [
		'src/cpgplot-wrappers.c'
		]

# Base dependencies
depends_base = [
	'src/data_structs.h',
	'src/debug.h',
	'src/error.h',
	'src/geometry.h',
	'src/kappa.h',
	'src/memory.h',
	'src/fits-and-miriad-wrappers.h',
	'src/molec.h',
	'src/numerical.h',
	'src/physics.h',
	'src/python-wrappers.h',
	'src/zone.h',
	'src/zone-hdf5.h',
]

if MIRSUPPORT:
	depends_base += [
		'src/cpgplot-wrappers.h'
		]


# SPARX sources files
sources_sparx = [
	'src/sparx-python.c',
	'src/sparx-test.c',
	'src/sparx-model.c',
	'src/sparx-physics.c',
	'src/sparx-inputs.c',
	'src/sparx-io.c',
	'src/sparx-utils.c',
]

# SPARX dependencies
depends_sparx = ['src/sparx.h']

##
## Distutils setup
##
from distutils.core import setup, Extension
from distutils import ccompiler, unixccompiler

# Things to include if MPI is available
if HAVE_MPI:
	macros += [('HAVE_MPI', None)]
	libs += mpi_libs
	#libs += ['lammpio', 'mpi', 'lam']
	#libs += ['mpich', 'pgc', 'pgftnrtl', 'pgftnrtl', 'nspgc', 'pgc', 'rt']
	#import os
	#os.environ['CC'] = "mpicc"

# Definition for the _sparx extension module

ext_sparx = Extension('sparx._sparx' if dev==0 else 'sparxdev._sparx',
	sources = sources_base+sources_sparx+[
		'src/sparx-pyext-_sparx.c',
		'src/sparx-task-amc.c',
		'src/sparx-task-telsim.c',
		'src/sparx-task-pygrid.c',
		'src/sparx-task-template.c',
	],
	depends = ['setup.py']+depends_base+depends_sparx,
	extra_compile_args = compiler_flags,
	define_macros = macros,
	include_dirs = header_dirs,
	library_dirs = lib_dirs,
	libraries = libs
)

# The main setup call
setup(
	name = 'sparx',
	version = Version,
	author = 'Eric Chung & I-Ta Hsieh',
	author_email = 'schung@asiaa.sinica.edu.tw / ita.hsieh@gmail.com',
	url = 'http://esclab.tw/wiki/index.php/Category:SPARX',
	description = 'SPARX Platform for Astrophysical Radiative Xfer',
	packages = ['sparx' if dev==0 else 'sparxdev'],
	package_dir = {'sparx' if dev==0 else 'sparxdev': "lib/sparx"},
	package_data = {'sparx' if dev==0 else 'sparxdev': [
		'data/molec/*.dat', # Molecular data files
		'data/opacity/*.tab', # Opacity data files
		'VERSION', # Program version
	]},
	ext_modules = [ext_sparx],
	scripts = [
		'bin/sparx', # Main sparx command line driver
		'bin/sparx-plot', # Model plotter
		'bin/sparx-plot.py', # Model plotter
		'bin/sparx-validate-dust.py', # Script for validating dust radiative transfer
		'bin/sparx-validate-line.py', # Script for validating line radiative transfer
		'bin/sparx-validate-leiden.py', # Script for validating with the Leiden 2004 benchmark problems
	],
)




