##
## The main sparx module
##

# Some necessary imports
import os
from math import sqrt, exp



# Get root package directory
ROOT_DIR = os.path.realpath(os.path.dirname(__file__))

##
## List of available molecules
##
import re
MOLEC_DIR = ROOT_DIR+"/data/molec" # Used by C modules, DO NOT RENAME!
MOLEC_LIST = []
for i in os.listdir(MOLEC_DIR):
	match = re.match("(.*).dat", i)
	if match is not None:
		MOLEC_LIST.append(match.groups()[0])
MOLEC_LIST.sort()

##
## List of available opacities
##
import re
KAPPA_DIR = ROOT_DIR+"/data/opacity" # Used by C modules, DO NOT RENAME!
KAPPA_LIST = []
for i in os.listdir(KAPPA_DIR):
	match = re.match("(.*).tab", i)
	if match is not None:
		KAPPA_LIST.append(match.groups()[0])
KAPPA_LIST.sort()













