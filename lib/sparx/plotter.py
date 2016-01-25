def __init__():
	import sys
	mod = sys.modules[__name__]
	mod.mod = mod

	import matplotlib as mpl
	mpl.use('TkAgg')
	import Tkinter as Tk

	# Attach ROOT widget and mainloop method
	mod.Tk = Tk
	mod.ROOT = ROOT = Tk.Tk()
	mod.mainloop = Tk.mainloop

	# Window title
	ROOT.wm_title("PySPARX Plotter")

	# Preven user from deleting the plot window
	def callback(): return
	ROOT.bind("<Destroy>", callback)
	ROOT.protocol("WM_DELETE_WINDOW", callback)

	# The quit program button
	button = Tk.Button(ROOT, text='Quit Program', command=sys.exit)
	button.pack(side=Tk.BOTTOM)
	button.pack(fill=Tk.BOTH, expand=1)
	return

# Init module
__init__()

################################################################################
# GUIPlotter class
class GUIPlotter:
	def __init__(self, name=None, figsize=(6, 4), position=(0.1, 0.1, 0.8, 0.8)):
		import matplotlib.backends.backend_tkagg as mpTkAgg
		import matplotlib.figure as mpfig
		import sys

		# Window widget
		self.window = Tk.Toplevel()
		if name is not None:
			assert type(name) is str
			self.window.title(name)
		self.name = self.window.title()

		def callback(): self.window.withdraw()
		self.window.bind("<Destroy>", callback)
		self.window.protocol("WM_DELETE_WINDOW", callback)

		button = Tk.Button(ROOT, text='Show %s window'%self.window.title(), command=self.window.deiconify)
		button.pack(side=Tk.TOP)
		button.pack(fill=Tk.BOTH, expand=1)

		# Instantiate figure and plot
		self.f = mpfig.Figure(figsize=figsize, dpi=100)
		self.ax = self.f.add_subplot(111)
		self.ax.set_position(position)
 
		# Instantiate canvas
		self.canvas = canvas = mpTkAgg.FigureCanvasTkAgg(self.f, self.window)
 
		# Pack canvas into window
		canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
		canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
 
		# Instantiate and pack toolbar
		self.toolbar = toolbar = mpTkAgg.NavigationToolbar2TkAgg(canvas, self.window)
 
		# Instantiate and pack quit button
		#self.button = button = Tk.Button(self.window, text='Quit', command=sys.exit)
		#button.pack(side=Tk.BOTTOM)
 
		# Show canvas and toolbar
		toolbar.update()
		canvas.show()

		# Init nplots to zero
		self.nplots = 0

		# Init names
		self.names = []
		return

	def plot(self, x_axis, y_axis, name=None, xlab=None, ylab=None, logx=False, logy=False, beg=None, end=None, lsty=None, msty=None, color=None, xfontsiz=8, yfontsiz=8):
		# Lists of lines, markers and colors
		line_styles = ('-', ':', '-.', '--')
		marker_styles = ('^', 'd', 'h', 'o', 'p', 's', 'v', 'x')
		color_styles = ('b', 'g', 'r', 'c', 'm', 'k')

		# Set plotting function according to logx and logy
		if logx and logy:
			plot = self.ax.loglog
		elif logx:
			plot = self.ax.semilogx
		elif logy:
			plot = self.ax.semilogy
		else:
			plot = self.ax.plot

		# Get subset if requested
		if beg is not None:
			x_axis = x_axis[int(beg):]
			y_axis = y_axis[int(beg):]
		if end is not None:
			x_axis = x_axis[:int(end)+1]
			y_axis = y_axis[:int(end)+1]

		# Select line, marker and color
		if lsty is not None:
			ls = lsty
		else:
			ls = line_styles[self.nplots % len(line_styles)]

		if msty is not None:
			ms = msty
		else:
			ms = marker_styles[self.nplots % len(marker_styles)]

		if color is not None:
			clr = color
		else:
			clr = color_styles[self.nplots % len(color_styles)]

		# Do the plotting
		plot(x_axis, y_axis, clr+ls+ms, markersize=3, label=name)

		# Set labels
		self.ax.set_xlabel(xlab)
		self.ax.set_ylabel(ylab)

		self.ax.xaxis.get_label().set_fontsize(xfontsiz)
		self.ax.yaxis.get_label().set_fontsize(yfontsiz)

		for tick in self.ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(xfontsiz)

		for tick in self.ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(yfontsiz)

		# Add name to namelist
		self.names += [name]

		# Increment nplots
		self.nplots += 1
		return

	def legend(self, legpos='best', fontsiz=10):
		import matplotlib.font_manager as fm
		leg_prop = fm.FontProperties(size=fontsiz)
		#self.ax.legend(prop=leg_prop, labelspacing=0.5, loc=legpos, shadow=True)
		self.ax.legend(prop=leg_prop, loc=legpos, shadow=True)
		return

	def figlegend(self, legpos='upper right', fontsiz=10):
		import matplotlib.font_manager as fm
		leg_prop = fm.FontProperties(size=fontsiz)
		self.f.legend(self.ax.lines, self.names, prop=leg_prop, labelspacing=0.5, borderaxespad=1, shadow=True)
		return

	def save(self, fname=None):
		if fname is None:
			fname = self.name+".png"
		self.canvas.print_figure(fname)
		return

	def refresh(self):
		self.canvas.show()
		return

	def show(self):
		self.refresh()
		return














