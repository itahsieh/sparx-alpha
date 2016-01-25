#! /usr/bin/env python

from tables import openFile
from numpy import ndarray, ndenumerate, concatenate, array, min
from math import sqrt, log10

################################################################################

def GetRadialData(fname, tname, pname, pindex=0):
	'''
	Return arrays of values for plotting.
	'''
	# Open HDF5 file
	try:
		h5f = openFile(fname)
	except:
		raise

	# Get dimensions
	rows = [i for i in h5f.root.ZONE]
	shape = [i for i in rows[0]['NAXES']]
	assert shape[0] > 0 and shape[1] == 1 and shape[2] == 1
	nzone = shape[0]

	# Get table data
	table = getattr(h5f.root, tname)
	tdata = [i[pname] for i in table]

	# Init and load plot data
	plot_data = [0 for i in range(nzone)]
	try:
		for i in range(nzone):
			plot_data[i] = tdata[i][pindex]
	except TypeError:
		for i in range(nzone):
			plot_data[i] = tdata[i]

	# Close HDF5 file
	h5f.close()

	return plot_data

################################################################################

def PlotRadial(fname_list, parm):
	for fname in fname_list:
		# Get X axis (radii)
		x_axis = GetRadialData(fname, "GRID", "X_cen")

		# Get Y axis
		if parm[0:3] == "lev" and len(parm) > 3:
			y_axis = GetRadialData(fname, "POPS", parm)
		elif parm[0:3] == "tau" and len(parm) > 3:
			y_axis = GetRadialData(fname, "TAU", "line"+parm[3:])
		elif parm[0:5] == "ratio" and len(parm) > 5:
			iline = int(parm[5:])
			assert iline > 0
			lev1 = GetRadialData(fname, "POPS", "lev%s" % (iline - 1))
			lev2 = GetRadialData(fname, "POPS", "lev%s" % (iline))
			n = len(lev1)
			y_axis = [0 for i in range(n)]
			for i in range(n):
				if lev1[i] > 0:
					y_axis[i] = lev2[i] / lev1[i]
		else:
			y_axis = GetRadialData(fname, "GRID", parm)

		assert len(x_axis) == len(y_axis)

		# Get subset if requested
		if opts.beg:
			x_axis = x_axis[int(opts.beg):]
			y_axis = y_axis[int(opts.beg):]
		if opts.end:
			x_axis = x_axis[:int(opts.end)+1]
			y_axis = y_axis[:int(opts.end)+1]

		# Output ascii if requested
		if opts.ascii:
			print "# ", fname
			for i in range(len(y_axis)):
				print "%20g %20g" % (x_axis[i], y_axis[i])
			print
		# Otherwise plot with Matplotlib
		else:
			# Check data range
			if opts.logy:
				for i in range(len(y_axis)):
					if y_axis[i] <= 0:
						raise Exception, "log Y scale requires values to be > 0"
			if opts.logx:
				for i in range(len(y_axis)):
					if x_axis[i] <= 0:
						raise Exception, "log X scale requires values to be > 0"

			# Plot the data
			from pylab import show, axis, xlim, ylim
			if opts.logx and opts.logy:
				from pylab import loglog as plot
			elif opts.logx:
				from pylab import semilogx as plot
			elif opts.logy:
				from pylab import semilogy as plot
			else:
				from pylab import plot as plot

			plot(x_axis, y_axis, "o-", label=fname)

			if opts.tight:
				axis('tight')
				xmin, xmax = xlim()
				ymin, ymax = ylim()
				delta_x = xmax - xmin
				delta_y = ymax - ymin
				v = [
					xmin,
					xmax,
					ymin - 0.1 * delta_y,
					ymax + 0.1 * delta_y
				]
				axis(v)

			from pylab import axis, title, xlabel, ylabel, legend
			title(parm)
			xlabel("Radius [pc]")
			ylabel(parm)
			if not opts.nolegend:
				leg = legend(loc='best')
				for t in leg.get_texts():
					t.set_fontsize('small')
					t.set_linespacing(0)


	if not opts.ascii:
		from pylab import show, savefig
		if opts.out:
			savefig(opts.out)
		else:
			show()

################################################################################

if __name__ == "__main__":
	# Configure options parser
	from optparse import OptionParser
	parser = OptionParser(usage="%prog <mode> <parm> <file>")
	parser.add_option("-o", metavar="File", dest="out", nargs=1, help="Output file name")
	parser.add_option("-a", "--ascii", dest="ascii", action="store_true", help="Print out ascii table instead of plotting")
	parser.add_option("--logx", dest="logx", action="store_true", help="log x axis")
	parser.add_option("--logy", dest="logy", action="store_true", help="log y axis")
	parser.add_option("--beg", dest="beg", action="store", help="Beginning index of X axis to plot")
	parser.add_option("--end", dest="end", action="store", help="Ending index of X axis to plot")
	parser.add_option("--nolegend", dest="nolegend", action="store_true", help="Don't draw legend")
	parser.add_option("--tight", dest="tight", action="store_true", help="Fit plot scale to data")

	(opts, args) = parser.parse_args()

	if len(args) < 2:
		parser.error("Not enough arguments")

	import matplotlib as mpl
	if opts.out:
		mpl.use('Agg')
	else:
		mpl.use('tkagg')

	parm = args[0]
	fname = args[1:]

	PlotRadial(fname, parm)





