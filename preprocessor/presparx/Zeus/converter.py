import numpy as np
import scipy.interpolate as scintp
from sparx_tc.pre_unit import AU2cm



import astropy.constants as const
mH = const.u.to('g').value
pc = const.pc.to('cm').value


# directory for the variable and geometry files
zeus_data_dir='/tiara/home/ithsieh/zeustw2sparx/zeus_data'
# Time stamp of the files
time_stamp='00100'

def FetchZeusData(filename):
    with open(zeus_data_dir+'/'+filename, 'rb') as f:
        data = np.fromfile(f, dtype=np.float64)
    return data

# Load radius and get dimensions
ra = FetchZeusData('z_x1ap')
rb  = FetchZeusData('z_x1bp')
nr = (ra.shape)[0]

# Load theta and get dimensions
thetaa = FetchZeusData('z_x2ap')
thetab  = FetchZeusData('z_x2bp')
nt = (thetaa.shape)[0]

# Load density
density = FetchZeusData('o_d__'+time_stamp+'_vThr')
density = np.reshape(density, (nt,nr))
density = np.rollaxis(density,1,0)

# Load velocity
Vr = FetchZeusData('o_v1_'+time_stamp)
Vr = np.reshape(Vr, (nt,nr))
Vr = Vr.transpose()


Vt = FetchZeusData('o_v2_'+time_stamp)
Vt = np.reshape(Vt, (nt,nr))
Vt = Vt.transpose()

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
        Density_cells = np.concatenate( ( density[gh:-gh,gh:-gh], quad*density[-gh:gh:-1,gh:-gh]), axis=1)
        Vr_bounds = np.concatenate((Vr[gh:-gh+1,gh:-gh+1], Vr[gh:-gh+1,-gh-1:gh-1:-1]), axis=1)
        Vth_bounds = np.concatenate((Vth[gh:-gh+1,gh:-gh+1], -Vth[gh:-gh+1,-gh-1:gh-1:-1]), axis=1)
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
        scintp.RegularGridInterpolator(( R_bounds, T_cells), Vr_bounds[:,1:], fill_value=0.0)
    fspline_th = \
        scintp.RegularGridInterpolator((R_cells, T_bounds), Vth_bounds[:-1,:], fill_value=0.0)
    
    
    R_map, T_map = np.meshgrid(R_cells,T_cells)
    
    Vr_cells = fspline_r(( R_map.transpose(), T_map.transpose()))
    Vth_cells = fspline_th((R_map.transpose(), T_map.transpose()))
    

    # Determine the outermost index of the R-direction
    rind_max = np.argmax(np.where(R_cells < Rmax))
    naxes_new[0] = rind_max + 1
    R_bounds = R_bounds[:naxes_new[0]+1]
    R_cells = R_cells[:naxes_new[0]]
    Density_cells = Density_cells[:naxes_new[0],:]
    Density_cells = Density_cells.reshape(naxes_new)
    #Vr_bounds = Vr_bounds[:naxes_new[0]+1,:]
    Vr_cells = Vr_cells[:naxes_new[0],:]
    Vr_cells = Vr_cells.reshape(naxes_new)
    #Vth_bounds = Vth_bounds[:naxes_new[0]+1,:]
    Vth_cells = Vth_cells[:naxes_new[0],:]
    Vth_cells = Vth_cells.reshape(naxes_new)
    Vph_cells = Vph_cells[:naxes_new[0],:]
    Vph_cells = Vph_cells.reshape(naxes_new)

    # Return properties 
    return (naxes_new, \
            R_bounds/pc, \
            T_bounds, \
            P_bounds, \
            Density_cells/(2.0*mH)*1e6, \
            Vr_cells*1e-2, \
            Vth_cells*1e-2, \
            Vph_cells*1e-2 \
            )

# Replicate ZeusTW output to form a full two-quadrant matrix if needed
naxes = [nr-6, nt-6, 1]



# maximum extension in the outflow axis (AU)
Rmax_AU = 12500.0


### data to input to SPARX: ###
# GridType: the geometry-type of the data (sph3d/cyl3d/rec3d)
# naxes : resolution of 3-D input, [n1,n2,n3]
# x1    : bounding points of first  dimension, 
#         to be (n1+1)-elements array, unit in [parsec]
# x2    : bounding points of second dimension, 
#         to be (n2+1)-elements array, unit in [parsec]
# x3    : bounding points of third  dimension, 
#         to be (n3+1)-elements array, unit in [parsec]
# n_H2  : number density of H2 molecule, 
#         3-D numpy array in the size of "naxes", 
#         value in the unit [number of H2 per cubic meter]
# T_k   : Kinetic temperature of gas,
#         3-D numpy array in the size of "naxes", 
#         value in the unit [Kelvin]
# T_d   : dust temperature,
#         3-D numpy array in the size of "naxes", 
#         value in the unit [Kelvin]
# v1    : gas velocity on the first  dimension, unit in [m/s]
# v2    : gas velocity on the second dimension, unit in [m/s]
# v3    : gas velocity on the third  dimension, unit in [m/s]
# Vt    : turbulent velocity

# grid type
GridType = "SPH3D" 
naxes, x1, x2, x3, n_H2, v1, v2, v3 = \
quadrantset( naxes, ra, rb, thetaa, thetab, \
    density, Vr, Vt, Rmax = Rmax_AU*AU2cm, mirror=True )
T_k = 10.0 * np.ones(naxes)
Vt  = 500.0 * np.ones(naxes)
# dust attributes
# T_d   : dust temperature
# kapp_d: dust opacity profile
# dust_to_gas : dust-to-gas ratio
T_d = T_k
kapp_d = 'table,jena_thin_e5'
dust_to_gas = 0.01 * np.ones(naxes)


# model attribute: molecular name
molec='co@xpol'
X_mol  = 1e-3 * np.ones(naxes)
# model attribute: CMB temperature
T_cmb=2.73


print "[ZeusTW2SPARX] Include the inner R = %d AU into SPARX HDF5 table" % Rmax_AU


