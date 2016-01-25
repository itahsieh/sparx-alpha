#! /usr/bin/env python

from tables import openFile
from numpy import ndarray, ndenumerate, concatenate, array, min
from math import sqrt, log10
import matplotlib as mpl
#LEGEND_PROP = mpl.font_manager.FontProperties(size=4)

################################################################################

def PylabShow():
	from pylab import show, savefig
	if opts.out:
		savefig(opts.out)
	else:
		show()

################################################################################

def slice_scalar(fname_list, tname, pname, slice_axis):
	# Open HDF5 file
	try:
		h5f = openFile(fname_list[0])
	except:
		raise

	# Get dimensions
	rows = [i for i in h5f.root.ZONE]
	shape = [i for i in rows[0]['NAXES']]
	naxes = [i for i in rows[0]['NAXES']]

	axis_titles = ["X", "Y", "Z"]
	slice_axis_title = axis_titles.pop(slice_axis)

	# Load data
	table = getattr(h5f.root, tname)
	data = [i[pname] for i in table]
	iter_max = shape.pop(slice_axis) # --> xy, xz or yz
	grid = ndarray(shape)

	if opts.index is not None:
		indices = [opts.index]
	else:
		indices = range(iter_max)

	for slice_index in indices:
		for i, val in ndenumerate(grid):
			index = [j for j in i]
			index.insert(slice_axis, slice_index)
			grid[i] = data[index[2] + naxes[2] * (index[1] + naxes[1] * index[0])]

		if opts.ascii:
			print "# %s, %s index=%d" % (pname, slice_axis_title, slice_index)
			for i in grid.transpose():
				for j in i:
					print "%10.4g" % j,
				print
			print
		else:
			from pylab import pcolor, colorbar, axis, title, xlabel, ylabel
			pcolor(grid.transpose())
			colorbar()
			axis([0, shape[0], 0, shape[1]])
			title("%s, %s index=%d" % (pname, slice_axis_title, slice_index))
			xlabel(axis_titles[0])
			ylabel(axis_titles[1])
			PylabShow()

	# Close HDF5 file
	h5f.close()

################################################################################

def radial_scalar(fname_list, tname, pname):
	for fname in fname_list:
		# Open HDF5 file
		try:
			h5f = openFile(fname)
		except:
			raise
		# Get dimensions
		rows = [i for i in h5f.root.ZONE]
		shape = [i for i in rows[0]['NAXES']]
		center = [i // 2 for i in rows[0]['NAXES']]
		X_max = [i for i in rows[0]['X_max']]
		X_min = [i for i in rows[0]['X_min']]
		X_delta = [(X_max[i] - X_min[i]) / shape[i] for i in range(len(shape))]

		# Load data
		table = getattr(h5f.root, tname)
		data = [i[pname] for i in table]
		grid = ndarray(shape)
		nprofile = min(shape) // 2
		profile = [0 for i in range(nprofile)]
		navg = [0 for i in range(nprofile)]
		radii = [(min(X_delta) * i) for i in range(nprofile)]

		for i, val in ndenumerate(grid):
			grid[i] = data[i[2] + shape[2] * (i[1] + shape[1] * i[0])]
			dx = i[0] - center[0]
			dy = i[1] - center[1]
			dz = i[2] - center[2]
			radius = int(sqrt(dx**2 + dy**2 + dz**2))

			if radius < len(profile):
				profile[radius] += grid[i]
				navg[radius] += 1

		for i in range(len(profile)):
			if navg[i] > 0:
				profile[i] /= navg[i]

		if opts.ascii:
			print "# ", fname
			for i in range(len(profile)):
				print "%10d %10g" % (i, profile[i])
			print
		else:
			if False:
				from pylab import plot as plot
				for i in range(len(profile)):
					profile[i] = log10(profile[i])
			else:
				if opts.logx and opts.logy:
					from pylab import loglog as plot
				elif opts.logx:
					from pylab import semilogx as plot
				elif opts.logy:
					from pylab import semilogy as plot
				else:
					from pylab import plot as plot
			plot(radii, profile, "o-", label=fname)
			import pylab as pl
			if opts.ymin:
				pl.ylim(ymin=float(opts.ymin))
			if opts.ymax:
				pl.ylim(ymax=float(opts.ymax))

		# Close HDF5 file
		h5f.close()


	if opts.ascii:
		pass
	else:
		from pylab import axis, title, xlabel, ylabel, legend
		xlabel("Radius (pc)")
		ylabel(pname)
		if not opts.nolegend:
			leg = legend(loc='best', prop=LEGND_PROP)
			for t in leg.get_texts():
				t.set_fontsize('small')
				t.set_linespacing(0)
		PylabShow()

################################################################################

def slice_velocity(fname_list, slice_axis):
	# Open HDF5 file
	try:
		h5f = openFile(fname_list[0])
	except:
		raise
	# Get dimensions
	rows = [i for i in h5f.root.ZONE]
	shape = [i for i in rows[0]['NAXES']]
	naxes = [i for i in rows[0]['NAXES']]

	axis_titles = ["X", "Y", "Z"]
	slice_axis_title = axis_titles.pop(slice_axis)

	axis_indices = [0, 1, 2]
	axis_indices.pop(slice_axis)

	# Load data
	table = getattr(h5f.root, "GRID")
	data = [i['V_cen'] for i in table]
	iter_max = shape.pop(slice_axis) # --> xy, xz or yz
	grid_i = ndarray(shape)
	grid_j = ndarray(shape)
	intensity = ndarray(shape)
	pos_x = ndarray(shape)
	pos_y = ndarray(shape)

	if opts.index is not None:
		indices = [opts.index]
	else:
		indices = range(iter_max)

	for slice_index in indices:
		for i, val in ndenumerate(grid_i):
			index = [j for j in i]
			index.insert(slice_axis, slice_index)
			row_id = index[2] + naxes[2] * (index[1] + naxes[1] * index[0])
			grid_i[i] = data[row_id][axis_indices[0]]
			grid_j[i] = data[row_id][axis_indices[1]]
			intensity[i] = sqrt(grid_i[i]**2 + grid_j[i]**2)
			pos_x[i] = i[0] + 0.5
			pos_y[i] = i[1] + 0.5

		from pylab import quiver, colorbar, axis, title, xlabel, ylabel, pcolor
		pcolor(intensity.transpose())
		colorbar()
		#quiver(pos_x.transpose(), pos_y.transpose(), grid_i.transpose(), grid_j.transpose(), scale=10000.0)
		quiver(pos_x.transpose(), pos_y.transpose(), grid_i.transpose(), grid_j.transpose(), scale=None)
		axis([0, shape[0], 0, shape[1]])
		title("Velocity, %s index=%d" % (slice_axis_title, slice_index))
		xlabel(axis_titles[0])
		ylabel(axis_titles[1])

		PylabShow()
	# Close HDF5 file
	h5f.close()

################################################################################

def radial_velocity(fname_list):
	for fname in fname_list:
		# Open HDF5 file
		try:
			h5f = openFile(fname)
		except:
			raise
		# Get dimensions
		rows = [i for i in h5f.root.ZONE]
		shape = [i for i in rows[0]['NAXES']]
		center = [i // 2 for i in rows[0]['NAXES']]
		X_max = [i for i in rows[0]['X_max']]
		X_min = [i for i in rows[0]['X_min']]
		X_delta = [(X_max[i] - X_min[i]) / shape[i] for i in range(len(shape))]

		# Load data
		table = getattr(h5f.root, "GRID")
		data = [i["V_cen"] for i in table]
		grid = ndarray(shape)
		nprofile = min(shape) // 2
		profile = [0 for i in range(nprofile)]
		navg = [0 for i in range(min(nprofile))]
		radii = [(min(X_delta) * i) for i in range(nprofile)]

		for i, val in ndenumerate(grid):
			vel = data[i[2] + shape[2] * (i[1] + shape[1] * i[0])]
			vel_mag = sqrt(vel[0]**2 + vel[1]**2 + vel[2]**2) / 1000.0
			dx = i[0] - center[0]
			dy = i[1] - center[1]
			dz = i[2] - center[2]
			radius = int(sqrt(dx**2 + dy**2 + dz**2))

			if radius < len(profile):
				profile[radius] += vel_mag
				navg[radius] += 1

		for i in range(len(profile)):
			if navg[i] > 0:
				profile[i] /= navg[i]

		if opts.ascii:
			print "# ", fname
			for i in range(len(profile)):
				print "%10d %10g" % (i, profile[i])
			print
		else:
			from pylab import semilogx as plot
			plot(radii, profile, label=fname)

		# Close HDF5 file
		h5f.close()


	if opts.ascii:
		pass
	else:
		from pylab import axis, title, xlabel, ylabel, legend
		xlabel("Radius (pc)")
		ylabel("Velocity (km/s)")
		if not opts.nolegend:
			legend(loc="best", prop=LEGND_PROP)
		PylabShow()



################################################################################

if __name__ == "__main__":
	# Configure options parser
	from optparse import OptionParser
	parser = OptionParser(usage="%prog <mode> <parm> <file>")
	parser.add_option("-o", metavar="File", dest="out", nargs=1, help="Output file name")
	parser.add_option("-x", metavar="Axis", dest="axis", nargs=1, help="Axis to slice along")
	parser.add_option("-i", metavar="Index", dest="index", nargs=1, help="Index along axis to slice at")
	parser.add_option("-a", dest="ascii", action="store_true", help="Print out ascii table instead of plotting")
	parser.add_option("--ymin", dest="ymin", action="store", help="Y axis lower limit")
	parser.add_option("--ymax", dest="ymax", action="store", help="Y axis upper limit")
	parser.add_option("--logx", dest="logx", action="store_true", help="log x axis")
	parser.add_option("--logy", dest="logy", action="store_true", help="log y axis")
	parser.add_option("--nolegend", dest="nolegend", action="store_true", help="Don't draw legend")
	(opts, args) = parser.parse_args()

	# Input parameters
	if opts.axis:
		slice_axis = int(opts.axis)
	else:
		slice_axis = 0

	if opts.index is not None:
		opts.index = int(opts.index)

	if len(args) < 3:
		parser.error("Not enough arguments")

	mode = args[0]
	parm = args[1]
	fname = args[2:]

	import matplotlib as mpl
	if opts.out:
		mpl.use('Agg')
	else:
		mpl.use('tkagg')

	if mode == "slice":
		if parm == "velo":
			slice_velocity(fname, slice_axis)
		elif parm == "pops":
			slice_scalar(fname, "POPS", "lev1", slice_axis)
		elif parm == "tau":
			slice_scalar(fname, "TAU", "line0", slice_axis)
		else:
			slice_scalar(fname, "GRID", parm, slice_axis)
	elif mode == "radial":
		if parm == "velo":
			radial_velocity(fname)
		elif parm == "pops":
			radial_scalar(fname, "POPS", "lev1")
		elif parm == "tau":
			radial_scalar(fname, "TAU", "line0")
		else:
			radial_scalar(fname, "GRID", parm)
	else:
		parser.error("`%s' is not a recongnized mode" % mode)




