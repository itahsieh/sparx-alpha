#! /bin/env python

# Python Distutils setup script for SPARX

# VERSION_NUMBER of SPARX
VERSION_NUMBER = '3.1'

# number of Thread using in per job
ENABLE_MULTITREADING = 1
if ENABLE_MULTITREADING:
        import multiprocessing
        NumberOfThread = 2 * multiprocessing.cpu_count()
else:
        NumberOfThread = 1
print 'Number Of Thread =',NumberOfThread


# Some necessary imports
import os, glob
from os.path import exists, realpath, expanduser
# Test for MPI by checking whether mpicc can be called
from subprocess import call, Popen, PIPE

if not exists('unit_tests/tmp'):
    os.makedirs('unit_tests/tmp')
HAVE_MPI = (call("mpicc src/mpi-test.c -o unit_tests/tmp/a.out", shell=True, stdout=PIPE, stderr=PIPE) == 0)
# Get svn revision and update VERSION
import time
p = Popen("svnversion", shell=True, stdout=PIPE)
REV = p.communicate()[0].strip()
fo = file("lib/sparx/VERSION", "w")
fo.write("%s (r%s, %s)"%(VERSION_NUMBER,REV, time.asctime()))
fo.close()

##
## Gather information for setting up the package
##                              


# Get Python paths
import sys
PYINC = sys.prefix+"/include"
PYLIB = sys.prefix+"/lib"
if not (exists(PYINC) and exists(PYLIB)):
	raise Exception, "Cannot locate Python include and lib directories"

# Get NumPy paths
import numpy
NPINC = numpy.get_include()


# MIRIAD support option
MIRSUPPORT = 1 if os.environ.get('MIRLIB') is not None else 0

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
SPARXVERSION='sparx'
SPARX_VERSION='sparx'
USER_INCLUDE = []
USER_LIB = []
args = sys.argv[:]
mpi_libs = ['mpi']
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
	elif arg.find('--version') == 0:
                SPARXVERSION += '-'+expanduser(arg.split('=')[1])
                SPARX_VERSION += '_'+expanduser(arg.split('=')[1])
                sys.argv.remove(arg)

# Define macros
macros = [
        ('NTHREAD', str(NumberOfThread)),
        ('MIRSUPPORT', MIRSUPPORT),
        ('SPARXVERSION', '\"' + SPARXVERSION + '\"' ),
        ('SPARX_VERSION', '\"' + SPARX_VERSION + '\"' ),
]

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
#	'-pedantic',
	'-fshort-enums',
	'-fno-common',
	'-Dinline=',
	'-g',
	'-rdynamic',
	'-O3',
	'-pthread',
#	'-Werror',
	'-Wall',
#	'-W',
	'-Wmissing-prototypes',
	'-Wstrict-prototypes',
	'-Wpointer-arith',
	'-Wcast-qual',
	'-Wcast-align',
	'-Wwrite-strings',
	'-Wnested-externs',
	#'-finline-limit=600',
	#'-fwhole-program',
	'-ftree-vectorize',
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
	'src/vtk-wrapper.c',
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
	'src/vtk-wrapper.h',
]

if MIRSUPPORT:
	depends_base += [
		'src/cpgplot-wrappers.h'
		]


# SPARX sources files
sources_sparx = [
	'src/sparx-python.c',
	#'src/sparx-test.c',
	'src/sparx-model.c',
	'src/sparx-physics.c',
	'src/sparx-inputs.c',
	'src/sparx-io.c',
	'src/sparx-utils.c',
	'src/sparx-ImageTracing.c',
]

# SPARX dependencies
depends_sparx = ['src/sparx.h']

##
## Distutils setup
##
from distutils.core import setup, Extension
from distutils import ccompiler, unixccompiler

#os.environ["CXX"] = "icc"
#os.environ["CC"] = "icc"

# Things to include if MPI is available
if HAVE_MPI:
	macros += [('HAVE_MPI', None)]
	libs += mpi_libs
	#libs += ['lammpio', 'mpi', 'lam']
	#libs += ['mpich', 'pgc', 'pgftnrtl', 'pgftnrtl', 'nspgc', 'pgc', 'rt']
	#os.environ['CC'] = "mpicc"

# Definition for the _sparx extension module

ext_sparx = Extension( SPARX_VERSION + '._sparx'  ,
	sources = sources_base+sources_sparx+[
		'src/sparx-pyext-_sparx.c',
		'src/sparx-task-amc.c',
		'src/sparx-task-telsim.c',
		'src/sparx-task-coldens.c',
		'src/sparx-task-visual.c',
		'src/sparx-task-pops2ascii.c',
		#'src/sparx-task-pygrid.c',
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
	version = VERSION_NUMBER,
	author = 'Eric Chung & I-Ta Hsieh',
	author_email = 'schung@asiaa.sinica.edu.tw / ita.hsieh@gmail.com',
	#url = 'http://esclab.tw/wiki/index.php/Category:SPARX',
        url = 'https://github.com/itahsieh/sparx-alpha',
        description = 'SPARX Platform for Astrophysical Radiative Xfer',
	packages = [SPARX_VERSION],
	package_dir = { SPARX_VERSION : "lib/sparx"},
	package_data = { SPARX_VERSION : [
		'data/molec/*.dat', # Molecular data files
		'data/opacity/*.tab', # Opacity data files
		'VERSION', # Program version
	]},
	ext_modules = [ext_sparx],
	scripts = [
                'bin/presparx', # SPARX preprocessor
		'bin/sparx', # Main sparx command line driver
		'bin/sparx-plot', # Model plotter
		'bin/sparx-plot.py', # Model plotter
		'bin/sparx-validate-dust.py', # Script for validating dust radiative transfer
		'bin/sparx-validate-line.py', # Script for validating line radiative transfer
		'bin/sparx-validate-leiden.py', # Script for validating with the Leiden 2004 benchmark problems
	],
)




