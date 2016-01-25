##
## Miscellaneous utilities
## All contents in this module should be independent from the sparx
## module so circular dependencies can be avoided
##

##
## Rank of process and size of process pool
##
## DO NOT delete or rename these two attributes unless you know
## absolutely what you're doing!
##
MPI_RANK = 0
MPI_SIZE = 1

# Some necessary imports
import numpy as np
from math import log10
import os

def confirm_remove_files(files):
	"""
	Convenience function for removing a bunch of files
	interactively
	"""
	from shutil import rmtree
	from os import remove
	# ans is whether the files have been removed
	ans = True
	if len(files) > 0:
		# Print list of files to be removed
		for i in files:
			print i

		# Prompt user for decision if files are to be removed
		while True:
			ans = raw_input("Remove the above files (y/N)? ")
			if ans == '' or ans == 'N' or ans == 'n':
				ans = False
				break
			elif ans == 'y' or ans == 'Y':
				ans = True
				break
			else:
				print "Please enter 'y' or 'n' to continue"

		# If ans is True, remove all files
		if ans:
			for i in files:
				try:
					rmtree(i)
				except:
					try:
						remove(i)
					except:
						raise
	return ans

################################################################################

class Message:
	"""
	Class for handling terminal messages for the user.
	This class is tied to the MPI_RANK attribute in this module,
	so be careful when moving it to somewhere else!
	"""
	from sys import stdout, stderr
	def __init__(self, prompt, outf=stdout, errf=stderr, debug=False):
		self.prompt = prompt
		self.outf = outf
		self.errf = errf
		self.debug = debug
		return

	def __call__(self, text, stream=None):
		text += "\n"
		if stream is None:
			stream = self.outf
		if self.debug:
			pos = self._GetSrcPos(currentframe().f_back)
			self.Raw(pos+":"+self.prompt+": "+text, stream)
		else:
			self.Raw(self.prompt+": "+text, stream)
		return

	def _GetSrcPos(self, frame):
		info = getframeinfo(frame)
		return "%s:%d" % (basename(info[0]), info[1])

	def _Raw(self, text, stream):
		if MPI_RANK == 0:
			print >>stream, text,
		return

	def SetPrompt(self, prompt):
		self.prompt = str(prompt)
		return

	def Raw(self, text, stream=None):
		if stream is None:
			stream = self.outf
		self._Raw(text, stream)
		return

	def Err(self, text):
		text += "\n"
		if self.debug:
			pos = self._GetSrcPos(currentframe().f_back)
			self.Raw(pos+":"+self.prompt+" error: "+text, self.errf)
		else:
			self.Raw(self.prompt+" error: "+text, self.errf)
		return

	def HR(self, char="="):
		self.Raw(char * 80 + "\n")
		return

	def Bonk(self, text, errno=1):
		self.Err(text)
		exit(errno)
		return

##
## Instantiate message object
##
MESG = Message("sparx", debug=False)

################################################################################

def call(prog, shell=True, **kwargs):
	"""
	Convenience function for calling shell commands
	"""
	from subprocess import call
	args = " ".join(["%s=\"%s\""%(key, str(kwargs[key])) for key in kwargs])
	ret = call(prog+" "+args, shell=shell)
	if ret != 0:
		raise Exception, "'%s' failed"%prog
	return

################################################################################

def generate_linear_points(min, max, n):
	"""
	Generate n linearly spaced data points between min and max
	"""
	delta = (max - min) / (n - 1)
	return np.array([min + delta * i for i in range(n)])

################################################################################

def generate_log_points(min, max, n):
	"""
	Generate n logarithmically spaced data points between min and max
	"""
	logmax = log10(max)
	logmin = log10(min)
	delta = (logmax - logmin) / (n - 1)
	return np.array([10.0**(logmin + delta * i) for i in range(n)])

################################################################################

def clear_path(path):
	"""
	Make sure path is available by checking for files/directories present
	and remove accordingly
	"""
	from shutil import rmtree
	from os.path import isdir, exists
	from os import remove
	if exists(path):
		if isdir(path):
			rmtree(path)
		else:
			remove(path)
	return

################################################################################

class ParmSpace:
	def __init__(self, parms, ondisk=True):
		"""
		'parms' is a list of parameter names
		"""
		# Keep track of sorted list of parameters
		self.parms = sorted(parms)

		self.ondisk = ondisk

		# Init new parmspace: this is a dictionary of dictionaries
		self.prm_dic = {}

		# This is for keeping track of parameters currently in the
		# parameter space (searching linear string lists is much faster
		# than searching multi-dimensional arrays
		self.prm_nlist = []
		self.prm_slist = []
		self.prm_ilist = []

		if self.ondisk:
			# Search in current directory for directories with names following
			# the pattern "prm-[0-9]+" and build parameter space on the fly.
			# These are directories that contain the "PARMS" text file and 
			# any other data defined by the user pipeline
			from os.path import exists, isdir
			import re
			for name in os.listdir(os.curdir):
				m = re.match(r"prm-([0-9]+)", name)
				if m and isdir(name):

					# Check whether PARMS file is present, it is
					# a fatal error if it doesn't exist since it should've
					# been saved by the insert() method
					fo = file(name+"/PARMS", "r")

					# Read parameters from PARMS: this is a list of parameters
					# inserted by the insert() method, with one key parameter per
					# line, and each line in the form of key=value
					string = fo.read()
					pdic = self.parse_parmstr(string)

					# Store parameter name, parmstr and id in lists
					self.prm_dic[name] = pdic
					self.prm_nlist += [name]
					self.prm_slist += [self.build_parmstr(**pdic)]
					self.prm_ilist += [int(m.groups()[0])]

		# Find gaps in the list of ids
		if len(self.prm_ilist) > 0:
			N = max(self.prm_ilist)+1
			dup_ilist = self.prm_ilist + range(N)
			self.prm_gaps = [i for i in dup_ilist if dup_ilist.count(i) == 1]
		else:
			self.prm_gaps = []
		return

	def build_parmstr(self, **parms):
		try:
			return "\n".join(["%10s=%s"%(key, parms[key]) for key in self.parms])
		except:
			print "blah", parms
			print "self.parms=", self.parms
			for key in self.parms:
				print key
				print parms[key]
			raise

	def parse_parmstr(self, parmstr):
		parms = {}
		lines = parmstr.split("\n")
		if len(lines) == 0:
			raise Exception, "Illegal parmstr '%s'"%parmstr
		for line in lines:
			pair = line.strip().split("=")
			if len(pair) >= 2:
				parms[pair[0]] = "=".join(pair[1:])
		return parms

	def insert(self, **parms):
		"""
		Insert parms combination into self.space and obtain unique id
		"""
		# Check whether input parms are valid
		for key in parms.keys():
			if key not in self.parms:
				raise Exception, "'%s' is not a valid parameter (%s)"%(key, ",".join(self.parms))

		# Check for missing parms
		missing_parms = []
		for key in self.parms:
			if key not in parms.keys():
				missing_parms += [key]
		if len(missing_parms) > 0:
			raise Exception, "Parameters '%s' missing"%(",".join(missing_parms))

		# Check for repeated parameter combinations
		parmstr = self.build_parmstr(**parms)
		if parmstr not in self.prm_slist:
			# Get unique id from next item in self.prm_gaps if there are gaps,
			# otherwise id is length of self.prm_ilist
			if len(self.prm_gaps) > 0:
				id = self.prm_gaps.pop(0)
			else:
				id = len(self.prm_ilist)
			name = "prm-%05d"%id
			self.prm_ilist += [id]
			self.prm_slist += [parmstr]

			assert name not in self.prm_dic # Just in case

			# Insert parameters
			self.prm_dic[name] = parms

			if self.ondisk:
				# Make directory
				os.mkdir(name)

				# Save parameter string in dir
				fo = file(name+"/PARMS", "w")
				print >>fo, parmstr
				fo.flush()
				os.fsync(fo.fileno())
				fo.close()
		else:
			# Repeated parameter, return its name
			name = self.prm_nlist[self.prm_slist.index(parmstr)]

			if self.ondisk:
				# Check whether parms file is present and
				# consistent with input
				fo = file(name+"/PARMS", "r")
				fparmstr = fo.read()
				fo.close()
				if not parmstr.strip() == fparmstr.strip():
					print "parmstr2=", parmstr
					print "fparmstr2=", fparmstr
					raise Exception, "Prams file in directory '%s' inconsistent with parms"%id

		# Return name
		return name

	def save(self, fname):
		"""
		Save parameter space to aother file
		"""
		# Open file for writing
		fo = file(fname, "w")

		# Write table header
		print >>fo, "\t".join(["name"] + self.parms)

		# Loop through names
		for name in sorted(self.prm_dic.keys()):
			print >>fo, "\t".join([name] + [self.prm_dic[name][key] for key in self.parms])

		# Close file
		os.fsync(fo.fileno())
		fo.close()
		return

	def filter(self, parm, min, max, convertor=None):
		"""
		Search for all models with parameter parm within the range [min, max].
		Returns a new ParmSpace instance with search results.
		"""
		ps = ParmSpace(self.parms, ondisk=False)

		# Convert min and max
		if convertor != None:
			min = convertor(min)
			max = convertor(max)

		# Loop through parm space and search for matching models
		for model in self.prm_dic:
			value = self.prm_dic[model][parm]
			if convertor != None:
				value = convertor(value)
			if value >= min and value <= max:
				ps.insert(**self.prm_dic[model])
		return ps

################################################################################

def calc_reduced_chisquare(theory, data, error):
	"""
	Calculate reduced Chi^2 of data
	"""
	if error <= 0.0:
		raise Exception, "Error must be > 0"
	n_t = len(theory)
	n_d = len(data)
	if n_t != n_d:
		raise Exception, "Data length inequal (%d vs %d)"%(n_t, n_d)
	chisquare = 0.0
	for i in range(n_d):
		chisquare += ((data[i] - theory[i]) / error)**2.0

	return chisquare / n_d

################################################################################

class RATRANDat:
	"""
	Interface to RATRAN data files (still under construction)
	"""
	def __init__(self, fname):
		# Open file
		fo = file(fname, "r")

		# Loop through lines and search for '@'
		for i in fo:
			line = i.strip()
			if line == '@':
				break

		# Load grid into table
		cols = [i.strip().split() for i in fo]
		# ignore last line which seems to be erroneous
		cols.pop(-1)
		r_list = [float(i[1]) for i in cols]
		return


