#! /bin/env python

import pysparx as ps

################################################################################
################################################################################

if __name__ == "__main__":
	# Configure options parser and parse options
	from optparse import OptionParser
	parser = OptionParser(usage="%prog <parm> <file>")
	parser.add_option("-o", metavar="File", dest="out", nargs=1, help="Output file name")
	parser.add_option("-a", dest="ascii", action="store_true", help="Print out ascii table instead of plotting")
	parser.add_option("--ymin", dest="ymin", action="store", help="Y axis lower limit")
	parser.add_option("--ymax", dest="ymax", action="store", help="Y axis upper limit")
	parser.add_option("--logx", dest="logx", action="store_true", help="log x axis")
	parser.add_option("--logy", dest="logy", action="store_true", help="log y axis")
	parser.add_option("--nolegend", dest="nolegend", action="store_true", help="Don't draw legend")
	parser.add_option("--tight", dest="tight", action="store_true", help="Fit plot scale to data")
	parser.add_option("--slice", metavar="Axis Index", dest="slice", nargs=2, help="Position to slice")
	(opts, args) = parser.parse_args()

	# Check number of arguments
	if len(args) < 2:
		parser.error("Not enough arguments")

	# Set plotting device
	import matplotlib as mpl
	if opts.out:
		mpl.use('agg')
	else:
		mpl.use('tkagg')

	# Set plot parameter
	parm = args[0]

	# Set number of subplots if plotting slices
	if opts.slice:
		nsub = int(len(args[1:])**0.5)
		slice_axis = int(opts.slice[0])
		slice_index = int(opts.slice[1])

	# Load files
	for i in range(len(args[1:])):
		fo = ps.SparxH5_3D(i)
		if opts.slice:
			slice = fo.GetSlice(slice_axis, slice_index, parm
			i_row = i // nsub
			i_col = i % nsub
	files = [ps.SparxH5_3D(i) for i in args[1:]]

	# Close files
	for i in files:
		i.Close()

