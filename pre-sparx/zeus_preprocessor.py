#!/usr/bin/env python
#
##################################################################
# Reading in ZeusTW output and converting into SPARX HDF5 input. #
##################################################################
#

# Python library import 
import numpy as np
import tables as pt
import scipy.interpolate as scintp
import astropy.constants as const
import matplotlib.pyplot as pl
from argparse import ArgumentParser

# Unit conversion
au = const.au.to('cm').value
pc = const.pc.to('cm').value
mH = const.u.to('g').value

#
# Take care of parsed arguments:
#
# Prepare argument parser
parser = ArgumentParser(description='ZeusTW to SPARX import interface')
parser.add_argument('ztwfile', default='o_d__00100', type=str, \
                    help='ZeusTW output filename')
parser.add_argument('-l', dest='dir', default='../', type=str, \
                    help='directory for the variable and geometry files')
parser.add_argument('-o', dest='output', default='model_zeustw', type=str, \
                    help='SPARX input filename')
parser.add_argument('-xmol', dest='xmol', default='1.0d-4', type=float, \
                    help='fractional abundance of the molecule wrt molecular hydrogen')
parser.add_argument('-tk', dest='T_k', default='400.0', type=float, \
                    help='gas temperature for the run (K)')
parser.add_argument('-td', dest='T_d', default='400.0', type=float, \
                    help='dust temperature for the run (K)')
parser.add_argument('-vt', dest='V_t', default='0.0', type=float, \
                    help='turbulence velocity (m/s)')
parser.add_argument('-rmax', dest='Rmax', default='12500.0', type=float, \
                    help='maximum extension in the outflow axis (AU)')
parser.add_argument('-gtd', dest='gtd', default='100.0', type=float, \
                    help='gas-to-dust ratio for the run')
parser.add_argument('--dust', dest='dust', action='store_true', \
                    help='include dust in the calculation')
parser.add_argument('--cmb', dest='cmb', action='store_true', \
                    help='include CMB (2.73 K) in the calculation')
parser.add_argument('--debug', dest='debug', action='store_true', \
                    help='draw debugging plots')
parser.add_argument('--vtk', dest='writevtk', action='store_true', \
                    help='store VTK file for viewing') 
parser.add_argument('--mirror', dest='mirror', action='store_true', \
                    help='if only half domain is available, mirror the variables; or otherwise pad with zeros')
args = parser.parse_args()

# Some flow-control tags
debug = args.debug           # Whether to draw debugging plots
writevtk = args.writevtk     # Whether to write in VTK file for viewing

# Physical parameters setting for SPARX
X_mol = args.xmol            # Fractional abundance of the molecular species
V_t = args.V_t               # Turbulence velocity (m/s)
T_k = args.T_k               # Kinetic temperature (K)
if (args.cmb == True):
    T_cmb = 2.73             # CMB mean temperature (K) [set to zero for no CMB]
else:
    T_cmb = 0.00
if (args.dust == True):
    gas_to_dust = args.gtd   # Gas-to-dust ratio [set to zero for no dust emission]
    T_d = args.T_d           # Dust temperature (K) [set to zero for no dust emission]
else:
    gas_to_dust = 0.0
    T_d = 0.0
geom = "sph3d"               # Geometrical gridding used in SPARX
molec = ""                   # Name of the given molecular species [kept empty until AMC assigns]
root = "/"                   # Table root for the output HDF file

# Parameter settings for ZeusTW input
mirror = args.mirror         # Whether to mirror the hemispheric data or padded with zero
Rmax_AU = args.Rmax          # Maximum coverage of r of the ZeusTW results to be included in SPARX

#
# Class definition for SPARX input table, in the format of pyTables
#
class Particle(pt.IsDescription):
    """
      Standard SPARX table definitions:
        0. LEVEL:   refinement level for AMR grids [Int32]
        1. POS:     position counter for the grids [Int64]
        2. geom:    geometry characteristics -- sph1d, sph3d, rec3d [String Size of 6]
        3. X_max:   cell outer boundary -- (r,theta,phi) for sph; (x,y,z) for rec (pc) [Float64]
        4. X_min:   cell inner boundary (pc) [Float64]
        5. X_cen:   cell center (pc) [Float64]
        6. n_H2:    molecular hydrogen number density (m^-3) [Float64]
        7. T_k:     kinetic temperature for the cells (K) [Float64]
        8. X_mol:   fractional abundance of the molecular species (wrt H2) [Float64]
        9. X_pH2:   fraction of the para-H2 [Float64]
       10. X_oH2:   fraction of the ortho-H2 [Float64]
       11. X_e:     electron fraction [Float64]
       12. X_H:     atomic hydrogen fraction [Float64]
       13. X_He:    helium fraction [Float64]
       14. V_t:     turbulence velocity (m/s) [Float64]
       15. V_edge:  velocities at the surfaces of the cell (m/s) [Float32]
       16. V_cen:   velocity interpolated to the cell center (m/s) [Float32]
       17. B_cen:   magnetic field strength at the cell center [Float32]
       18. ds:      [Float32]
       19. NCHILDREN: number of child cells in the next level [Int64]
       20. NAXES:   number of axes in the next level [Int64]
       21. T_d:     precalculated dust temperature [Float64]
       22. kapp_d:  name of the dust opacity table [String size of 64]
       23. T_ff:    free-free emission temperature [Float64]
       24. kapp_ff: name of the free-free transition opacity table [String size of 64]
       25. T_bb:    black-body temperature [Float64]
      Physical properties: MKS unit; Length unit: pc.
    """
    LEVEL     = pt.Int32Col(pos=0)
    POS       = pt.Int64Col(pos=1)
    geom      = pt.StringCol(itemsize=6,pos=2)
    X_max     = pt.Float64Col(shape=3,pos=3)
    X_min     = pt.Float64Col(shape=3,pos=4)
    X_cen     = pt.Float64Col(shape=3,pos=5)
    n_H2      = pt.Float64Col(pos=6)
    T_k       = pt.Float64Col(pos=7)
    X_mol     = pt.Float64Col(pos=8)
    X_pH2     = pt.Float64Col(pos=9)
    X_oH2     = pt.Float64Col(pos=10)
    X_e       = pt.Float64Col(pos=11)
    X_H       = pt.Float64Col(pos=12)
    X_He      = pt.Float64Col(pos=13)
    V_t       = pt.Float64Col(pos=14)
    V_edge    = pt.FloatCol(shape=(6,3),pos=15)
    V_cen     = pt.FloatCol(shape=3,pos=16)
    B_cen     = pt.FloatCol(shape=3,pos=17)
    ds        = pt.FloatCol(pos=18)
    NCHILDREN = pt.Int64Col(pos=19)
    NAXES     = pt.Int64Col(shape=3,pos=20)
    T_d       = pt.Float64Col(pos=21)
    kapp_d    = pt.StringCol(itemsize=64,pos=22)
    T_ff      = pt.Float64Col(pos=23)
    kapp_ff   = pt.StringCol(itemsize=64,pos=24)
    T_bb      = pt.Float64Col(pos=25)

# 
# Function definition for generating both sides of equatorial plane if necessary
#
def quadrantset(naxes,R_ap,R_bp,T_ap,T_bp,density,Vr,Vth,Rmax=1e99,mirror=True):

    # Set ghost zone size for ZeusTW
    gh = 3

    # Determine whether the input data cover one or two quadrant(s)
    Theta_min = T_ap[gh] / (np.pi/2.0)
    Theta_max = T_ap[-gh] / (np.pi/2.0)

    # Padding, mirroring, or leaving as is
    if (np.abs(Theta_max - 1.0) <= 0.01):
        print "[ZeusTW2SPARX] Only one quadrant exists. Generate the other quardrant."
        if (mirror == True):
            print "[ZeusTW2SPARX] Mirror into two opposite quadrants."
	    quad = 1.0
        elif (mirror == False):
            print "[ZeusTW2SPARX] Pad zeros to the other quadrant."
	    quad = 0.0
        naxes_new = [ naxes[0], naxes[1]*2, naxes[2] ]
        R_bounds = R_ap[gh:-gh+1]
        T_bounds = np.concatenate((T_ap[gh:-gh+1], np.pi-T_ap[-gh-1:gh-1:-1]))
        R_cells = R_bp[gh:-gh]
        T_cells = np.concatenate((T_bp[gh:-gh], np.pi-T_bp[-gh-1:gh-1:-1]))
        Density_cells = np.concatenate( ( density[gh:-gh,gh:-gh], quad*density[-gh:gh:-1,gh:-gh]), axis=0)
        Vr_bounds = np.concatenate((Vr[gh:-gh+1,gh:-gh+1], Vr[-gh-1:gh-1:-1,gh:-gh+1]), axis=0)
        Vth_bounds = np.concatenate((Vth[gh:-gh+1,gh:-gh+1], -Vth[-gh-1:gh-1:-1,gh:-gh+1]), axis=0)
    elif (np.abs(Theta_max - 2.0) <= 0.01):
        print "[ZeusTW2SPARX] Two quadrants exist. No need for mirroring."
        naxes_new = naxes
        R_bounds = R_ap[gh:-gh+1]
        T_bounds = T_ap[gh:-gh+1]
        R_cells = R_bp[gh:-gh]
        T_cells = T_bp[gh:-gh]
        Density_cells = density[gh:-gh,gh:-gh]
        Vr_bounds = Vr[gh:-gh+1,gh:-gh+1]
        Vth_bounds = Vth[gh:-gh+1,gh:-gh+1]

    # Made-up Azimuthal properties for axisymmetric data set
    P_bounds = np.array([0.0, 2.0*np.pi])
    P_cells = np.array([np.pi])
    Vph_cells = np.zeros(Density_cells.shape)

    # Two-dimensional interpolation for the velocities onto cell center
    fspline_r = \
        scintp.RegularGridInterpolator((T_cells, R_bounds), Vr_bounds[1:,:], fill_value=0.0)
    fspline_th = \
        scintp.RegularGridInterpolator((T_bounds, R_cells), Vth_bounds[:,:-1], fill_value=0.0)
    T_map, R_map = np.meshgrid(T_cells, R_cells)
    Vr_cells = fspline_r((T_map.transpose(), R_map.transpose()))
    Vth_cells = fspline_th((T_map.transpose(), R_map.transpose()))

    # Determine the outermost index of the R-direction
    rind_max = np.argmax(np.where(R_cells < Rmax))
    naxes_new[0] = rind_max + 1
    R_bounds = R_bounds[:naxes_new[0]+1]
    R_cells = R_cells[:naxes_new[0]]
    Density_cells = Density_cells[:,:naxes_new[0]]
    Vr_bounds = Vr_bounds[:,:naxes_new[0]+1]
    Vr_cells = Vr_cells[:,:naxes_new[0]]
    Vth_bounds = Vth_bounds[:,:naxes_new[0]+1]
    Vth_cells = Vth_cells[:,:naxes_new[0]]
    Vph_cells = Vph_cells[:,:naxes_new[0]]

    # Return properties 
    return (naxes_new, R_bounds, R_cells, T_bounds, T_cells, P_bounds, P_cells, \
                       Density_cells, Vr_bounds, Vr_cells, Vth_bounds, Vth_cells, Vph_cells)

#
# Function definition for writing the "ZONE" properties in HDF table
#
def writezone(direc,lev,position,xmax,xmin,naxes):

    # Write a row in the grid table
    table = h5file.create_table(direc, 'ZONE', Particle, "Grid table")
    zone = table.row
    zone['LEVEL']  = lev
    zone['POS'] = position
    zone['geom'] = geom
    zone['X_max'] = [ xmax[0]/pc,xmax[1],xmax[2] ]
    zone['X_min'] = [ xmin[0]/pc,xmin[1],xmin[2] ]
    zone['X_cen'] = [ 0.5*(xmin[0]+xmax[0])/pc,0.5*(xmin[1]+xmax[1]),0.5*(xmin[2]+xmax[2]) ]
    zone['NCHILDREN'] = naxes[0]*naxes[1]*naxes[2]
    zone['NAXES'] = naxes
    zone.append()

    # Flush table and delete temporary variables
    table.flush()
    for p in range(25+1):
        exec 'del table.attrs.FIELD_%d_FILL' % p
    del table.attrs.NROWS

#
# Function definition for writing the "GRID" properties in HDF table
#
def writegrid(direc,lev,naxes,\
              R_bounds,R_cells,T_bounds,T_cells,P_bounds,P_cells,\
              Density_cells,Vr_cells,Vth_cells,Vph_cells,writevtk):

    # Initialize the grid table
    table = h5file.create_table(direc, 'GRID', Particle, "Grid table")
    cell = table.row
    Vx = np.zeros(naxes)
    Vy = np.zeros(naxes)
    Vz = np.zeros(naxes)
    
    # Start looping and putting values into rows
    nc = 0
    for i in range(naxes[0]):
        for j in range(naxes[1]):
            for k in range(naxes[2]):
                
                # Write a row of grid table
                cell['LEVEL']  = lev + 1
                cell['POS'] = naxes[1]*naxes[2]*i + naxes[2]*j + k
                cell['geom'] = geom
                cell['X_max'] = [ R_bounds[i+1]/pc, T_bounds[j+1], P_bounds[k+1] ]
                cell['X_min'] = [ R_bounds[i]/pc, T_bounds[j], P_bounds[k] ]
                cell['X_cen'] = [ R_cells[i]/pc, T_cells[j], P_cells[k] ]
                cell['NCHILDREN'] = 0
                cell['NAXES'] = np.zeros((1,3))
                if (Density_cells[j,i] <= 1e-30): # empty (toroid region or ignored half)
                    cell['n_H2'] = 0.0
                    cell['T_k'] = 0.0
                    cell['X_mol'] = 0.0
                else: 
                    cell['n_H2'] = Density_cells[j,i]/(2.0*mH)*1e6 # g/cm^3 -> m^-3
                    cell['T_k'] = T_k
                    cell['X_mol'] = X_mol
                cell['X_pH2'] = 0.0
                cell['X_oH2'] = 0.0
                cell['X_e']   = 0.0
                cell['X_H']     = 0.0
                cell['X_He']    = 0.0
                cell['V_t']     = V_t
                cell['V_edge']  = np.zeros((6,3))
                cell['V_cen']   = [ Vr_cells[j,i]*1e-2, Vth_cells[j,i]*1e-2, Vph_cells[j,i]*1e-2 ]
                cell['B_cen']   = np.zeros((1,3))
                cell['ds']      = 0.0
                cell['T_d']     = T_d
                cell['kapp_d']  = 'table,jena_bare_e6'
                cell['T_ff']    = 0.0
                cell['kapp_ff'] = ''
                cell['T_bb']    = 0.0

                cell.append()
                nc = nc + 1

    # Flush table and delete temporary variables
    table.flush()
    for p in range(25+1):
        exec 'del table.attrs.FIELD_%d_FILL' % p
    del table.attrs.NROWS

    return nc

################################## 
# 
# Main Routine Starts
#
##################################

#
# Reading ZeusTW output
# 
#  data directory
datadir = args.dir
datafile = args.ztwfile
# Load radius and get dimensions
tmpfile = datadir+'z_x1ap'
with open(tmpfile, 'rb') as f:
    ra = np.fromfile(f, dtype=np.float64)
tmpfile = datadir+'z_x1bp'
with open(tmpfile, 'rb') as f:
    rb = np.fromfile(f, dtype=np.float64)
nr = (ra.shape)[0]
# Load theta and get dimensions
tmpfile = datadir+'z_x2ap'
with open(tmpfile, 'rb') as f:
    thetaa = np.fromfile(f, dtype=np.float64)
tmpfile = datadir+'z_x2bp'
with open(tmpfile, 'rb') as f:
    thetab = np.fromfile(f, dtype=np.float64)
nt = (thetaa.shape)[0]
# Load density
tmpfile = datadir+datafile
with open(tmpfile, 'rb') as f:
    density = np.fromfile(f, dtype=np.float64)
density = np.reshape(density, (nt,nr))
# Load velocity
tmpfile = datadir+'o_v1_00100'
with open(tmpfile, 'rb') as f:
    Vr = np.fromfile(f, dtype=np.float64)
Vr = np.reshape(Vr, (nt,nr))
tmpfile = datadir+'o_v2_00100'
with open(tmpfile, 'rb') as f:
    Vt = np.fromfile(f, dtype=np.float64)
Vt = np.reshape(Vt, (nt,nr))

# Replicate ZeusTW output to form a full two-quadrant matrix if needed
naxes = [nr-6, nt-6, 1]
naxes_new, r_bounds, r_cells, theta_bounds, theta_cells, phi_bounds, phi_cells, \
density_cells, Vr_bounds, Vr_cells, Vth_bounds, Vth_cells, Vph_cells = \
    quadrantset(naxes,ra,rb,thetaa,thetab,density,Vr,Vt,Rmax=Rmax_AU*au,mirror=mirror)
print "[ZeusTW2SPARX] Include the inner R = %d AU into SPARX HDF5 table" % Rmax_AU

# Convert R-theta to X-Z 
Xab = np.outer(np.sin(theta_bounds), r_cells)
Zab = np.outer(np.cos(theta_bounds), r_cells)
Xba = np.outer(np.sin(theta_cells), r_bounds)
Zba = np.outer(np.cos(theta_cells), r_bounds)
Xbb = np.outer(np.sin(theta_cells), r_cells)
Zbb = np.outer(np.cos(theta_cells), r_cells)
                
# Debug plotting: Data read in
if debug == True:

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.ticker import MultipleLocator, LinearLocator
    print('[ZeusTW2SPARX][debug] Plotting input ZeusTW arrays...')

    f, (ax1, ax2, ax3) = pl.subplots(1, 3, sharex=True, sharey=True)

    im1 = ax1.contour(Xbb/au*0.01, Zbb/au*0.01, np.log10(density_cells), 50)
    ax1.set_title(r'$\log_{10}$(Density)')
    ax1.set_xlim((0,30))
    ax1.set_ylim((-120,120))
    div1 = make_axes_locatable(ax1)
    cax1 = div1.append_axes("right", size="20%", pad=0.05)
    cbar1 = pl.colorbar(im1, cax=cax1, ticks=LinearLocator(), format="%d")

    im2 = ax2.contour(Xba/au*0.01, Zba/au*0.01, Vr_bounds[1:,:]*1e-5, 50)
    ax2.set_title(r'$V_r$')
    ax2.set_xlim((0,30))
    ax2.set_ylim((-120,120))
    div2 = make_axes_locatable(ax2)
    cax2 = div2.append_axes("right", size="20%", pad=0.05)
    cbar2 = pl.colorbar(im2, cax=cax2, ticks=LinearLocator(), format="%.1f")

    im3 = ax3.contour(Xab/au*0.01, Zab/au*0.01, Vth_bounds[:,:-1]*1e-5, 50)
    ax3.set_title(r'$V_\theta$')
    ax3.set_xlim((0,30))
    ax3.set_ylim((-120,120))
    div3 = make_axes_locatable(ax3)
    cax3 = div3.append_axes("right", size="20%", pad=0.05)
    cbar3 = pl.colorbar(im3, cax=cax3, ticks=MultipleLocator(1.0), format="%.1f")

    ax1.set_ylabel(r'$z$ (100 AU)')
    ax2.set_xlabel(r'$\varpi$ (100 AU)')

   #pl.show()
    pl.savefig('ZeusTW_View_Vbounds.eps')
    ax1.set_xlim((0,6))
    ax1.set_ylim((-24,24))
    pl.savefig('ZeusTW_ZoomView_Vbounds.eps')

    f, (ax1, ax2, ax3) = pl.subplots(1, 3, sharex=True, sharey=True)

    im1 = ax1.contour(Xbb/au*0.01, Zbb/au*0.01, np.log10(density_cells), 50)
    ax1.set_title(r'$\log_{10}$(Density)')
    ax1.set_xlim((0,30))
    ax1.set_ylim((-120,120))
    div1 = make_axes_locatable(ax1)
    cax1 = div1.append_axes("right", size="20%", pad=0.05)
    cbar1 = pl.colorbar(im1, cax=cax1, ticks=LinearLocator(), format="%d")

    im2 = ax2.contour(Xbb/au*0.01, Zbb/au*0.01, Vr_cells*1e-5, 50)
    ax2.set_title(r'$V_r$ (cell ctr intp.)')
    ax2.set_xlim((0,30))
    ax2.set_ylim((-120,120))
    div2 = make_axes_locatable(ax2)
    cax2 = div2.append_axes("right", size="20%", pad=0.05)
    cbar2 = pl.colorbar(im2, cax=cax2, ticks=LinearLocator(), format="%.1f")

    im3 = ax3.contour(Xbb/au*0.01, Zbb/au*0.01, Vth_cells*1e-5, 50)
    ax3.set_title(r'$V_\theta$ (cell ctr intp.)')
    ax3.set_xlim((0,30))
    ax3.set_ylim((-120,120))
    div3 = make_axes_locatable(ax3)
    cax3 = div3.append_axes("right", size="20%", pad=0.05)
    cbar3 = pl.colorbar(im3, cax=cax3, ticks=MultipleLocator(1.0), format="%.1f")

    ax1.set_ylabel(r'$z$ (100 AU)')
    ax2.set_xlabel(r'$\varpi$ (100 AU)')

   #pl.show()
    pl.savefig('ZeusTW_View_Vcells.eps')
    ax1.set_xlim((0,6))
    ax1.set_ylim((-24,24))
    pl.savefig('ZeusTW_ZoomView_Vcells.eps')

# HDF table output
print "[ZeusTW2SPARX] Writing output HDF5 table"
outfilename = args.output
h5file = pt.open_file(outfilename, mode = "w", title = "ZeusTW file")
h5file.del_node_attr("/", "TITLE", name=None)
h5file.del_node_attr("/", "CLASS", name=None)
h5file.del_node_attr("/", "VERSION", name=None)
h5file.del_node_attr("/", "PYTABLES_FORMAT_VERSION", name=None)
h5file.set_node_attr("/", "molec", molec, name=None)
h5file.set_node_attr("/", "T_cmb", T_cmb, name=None)
h5file.set_node_attr("/", "gas_to_dust", gas_to_dust, name=None)
h5file.set_node_attr("/", "velfield", "grid ", name=None)	
writezone(root,-1,0,\
          [r_bounds[-1],theta_bounds[-1],phi_bounds[-1]],\
          [r_bounds[0],theta_bounds[0],phi_bounds[0]],naxes_new)
ncell = writegrid(root,-1,naxes_new,\
                  r_bounds,r_cells,theta_bounds,theta_cells,phi_bounds,phi_cells,\
                  density_cells,Vr_cells,Vth_cells,Vph_cells,writevtk)
h5file.close()

print '[ZeusTW2SPARX] Total cells = ', ncell
# if (writevtk):
#    print "Wrote out",vtkfilea,vtkfileb

