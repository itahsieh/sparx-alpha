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

# POLARIZATION attributes
# b1    : B-field on the first  dimension, unit in [m/s]
# b2    : B-field on the second dimension, unit in [m/s]
# b3    : B-field on the third  dimension, unit in [m/s]

# DUST attributes:
# T_d   : dust temperature
# kapp_d: dust opacity profile
# dust_to_gas : dust-to-gas ratio

# MOLECULAR attributes: 
# molec : molecular name
# X_mol : molecular abundance

# BOUNDARY attribute: 
# T_cmb : CMB temperature
TIARA_CLUSTER = ['oc','tc','px','xl']

import numpy as np
import scipy.interpolate as scintp
from os.path import isfile
import zeus_parameter as ZeusPar

# the number of the Ngzost zone for ZeusTW
Ngz = 3

# Pre-check : grid type
if      ZeusPar.GridType == 'SPH':
    GridType = "SPH3D"
elif    ZeusPar.GridType == 'CYL':
    GridType = "CYL3D"
elif    ZeusPar.GridType == 'REC':
    GridType = "REC3D"
else:
    print 'GridType in zeus_parameter.py is not recognized. It should be SPH/CYL/REC .'
    exit(1)

# Pre-check : data type
if ZeusPar.DataType == 'BINARY':
    read_mode = 'rb'
elif ZeusPar.DataType == 'ASCII':
    read_mode = 'r'
else:
    print 'DataType in zeus_parameter.pyis not recognized. It should be BINARY/ASCII'
    exit(1)


# function to fetch Zeus data
def FetchZeusData(filename):
    filepath = ZeusPar.ZeusDataDir+'/'+filename
    if isfile(filepath):
        with open(filepath, read_mode) as f:
            data = np.fromfile(f, dtype=np.float64)
        return data
    else:
        print 'no %s data' % filename
        return None


# Load the griding 
# the bounding points on the first dimension
x1a = FetchZeusData('z_x1ap')
# the cell     points on the first dimension
x1b = FetchZeusData('z_x1bp')
# the bounding points on the first dimension
x2a = FetchZeusData('z_x2ap')
# the cell     points on the first dimension
x2b = FetchZeusData('z_x2bp')
# the bounding points on the first dimension
x3a = FetchZeusData('z_x3ap')
# the cell     points on the first dimension
x3b = FetchZeusData('z_x3bp')

# the resolution of the original data
n1 = x1a.shape[0]
n2 = x2a.shape[0]
n3 = x3a.shape[0] if x3a is not None else 1

# reverse the index order from (k,j,i) to (i,j,k)
def ReverseIndex(data):
    naxes = data.shape
    len_naxes = len(naxes)
    new_naxes = ()
    for i in range(len_naxes):
        new_naxes += (naxes[len_naxes-1-i],) 
    new_data = np.zeros(new_naxes)
    for i in range(new_naxes[0]):
        for j in range(new_naxes[1]):
            for k in range(new_naxes[2]):
                new_data[i,j,k] = data[k,j,i] 
    return new_data

# Function to fetch Zeus Physical data
def FetchZeusPhys(filename):
    data = FetchZeusData(filename)
    if data is None:
        return None
    data = np.reshape( data, (n3,n2,n1))
    return ReverseIndex(data)

def CheckAbdSetAtrr(attr):
    if hasattr(ZeusPar,attr):
        return getattr(ZeusPar,attr)
    else:
        print("{0} has no {1} data".format(ZeusPar.__name__,attr) )
        return None

def CheckAndSetArray(attr):
    if hasattr(ZeusPar,attr):
        return getattr(ZeusPar,attr) * np.ones(naxes)
    else:
        print("{0} has no {1} data".format(ZeusPar.__name__,attr) )
        return None
    
# Time stamp of the files
time_stamp = str('%.5d' % ZeusPar.TimeStamp)

# Load density
if ZeusPar.PostProcessedDensity:
    density = FetchZeusPhys('o_d__'+time_stamp+'_vThr')
else:
    density = FetchZeusPhys('o_d__'+time_stamp)

# Load temperature if desired
if hasattr( ZeusPar, 'ConstantTemperature'):
    if ZeusPar.ConstantTemperature > 0.0:
        temperature = ZeusPar.ConstantTemperature
    else:
        temperature = FetchZeusPhys('o_T__'+time_stamp)
else:
    temperature = FetchZeusPhys('o_T__'+time_stamp)

# Load velocity
V1 = FetchZeusPhys('o_v1_'+time_stamp)
V2 = FetchZeusPhys('o_v2_'+time_stamp)
V3 = FetchZeusPhys('o_v3_'+time_stamp)

# Load B-field
B1 = FetchZeusPhys('o_b1_'+time_stamp)
B2 = FetchZeusPhys('o_b2_'+time_stamp)
B3 = FetchZeusPhys('o_b3_'+time_stamp)

# 
# Function definition for generating both sides of equatorial plane if necessary
#
def Sph_MirrorNTrimR( 
    naxes, R_ap, R_bp, T_ap, T_bp, P_ap, P_bp, density, Vr, Vth, Vp,
    Rmax=1e99, T = temperature ):

    R_bounds = R_ap[Ngz:-Ngz+1]
    R_cells =  R_bp[Ngz:-Ngz]
    
    # Determine whether the input data cover one or two quadrant(s)
    Theta_max = T_ap[-Ngz] / (np.pi/2.0)
    # Padding, mirroring, or leaving as is
    
    if (np.abs(Theta_max - 1.0) <= 0.01):
        print "[ZeusTW2SPARX] Only one quadrant exists. Generate the other quardrant."
        print "[ZeusTW2SPARX] Mirror into two opposite quadrants."        
        naxes_new = [ naxes[0], naxes[1]*2, naxes[2] ]
        T_bounds = np.concatenate((T_ap[Ngz:-Ngz+1], np.pi-T_ap[-Ngz-1:Ngz-1:-1]))
        T_cells =  np.concatenate((T_bp[Ngz:-Ngz], np.pi-T_bp[-Ngz-1:Ngz-1:-1]))
        
        upper_hemisphere_index = slice(Ngz,-Ngz),slice(Ngz,-Ngz),   slice(0,naxes[2]) if naxes[2] == 1 else slice(Ngz,-Ngz)
        lower_hemisphere_index = slice(Ngz,-Ngz),slice(-Ngz,Ngz,-1),slice(0,naxes[2]) if naxes[2] == 1 else slice(Ngz,-Ngz)
        
        Density_cells = np.concatenate((density[upper_hemisphere_index],density[lower_hemisphere_index]), axis=1)
        if np.size(T) > 1:
            Temperature_cells = np.concatenate( ( T[upper_hemisphere_index], T[lower_hemisphere_index]), axis=1)
        
        upper_hemisphere_index = slice(Ngz,-Ngz+1),slice(Ngz,-Ngz+1),    slice(0,naxes[2]) if naxes[2] == 1 else slice(Ngz,-Ngz+1)
        lower_hemisphere_index = slice(Ngz,-Ngz+1),slice(-Ngz-1,Ngz-1,-1),slice(0,naxes[2]) if naxes[2] == 1 else slice(Ngz,-Ngz+1)
        
        Vr_bounds =     np.concatenate(( Vr[upper_hemisphere_index],   Vr[lower_hemisphere_index]), axis=1)
        Vth_bounds =    np.concatenate((Vth[upper_hemisphere_index], -Vth[lower_hemisphere_index]), axis=1)
        if Vp is not None:
            Vp_bounds = np.concatenate(( Vp[upper_hemisphere_index],   Vp[lower_hemisphere_index]), axis=1)

    elif (np.abs(Theta_max - 2.0) <= 0.01):
        print "[ZeusTW2SPARX] Two quadrants exist. No need for mirroring."
        naxes_new = naxes
        T_bounds = T_ap[Ngz:-Ngz+1]
        T_cells =  T_bp[Ngz:-Ngz]
        
        index = slice(Ngz,-Ngz), slice(Ngz,-Ngz), slice(0,naxes[2]) if naxes[2] == 1 else slice(Ngz,-Ngz)
        
        Density_cells = density[index]
        if np.size(T) > 1:
            Temperature_cells =   T[index] 
        
        index = slice(Ngz,-Ngz+1),slice(Ngz,-Ngz+1),    slice(0,naxes[2]) if naxes[2] == 1 else slice(Ngz,-Ngz+1)
        
        Vr_bounds =          Vr[index]
        Vth_bounds =        Vth[index]
        if Vp is not None:
            Vp_bounds =      Vp[index]
        

    # Made-up Azimuthal properties for axisymmetric data set
    if naxes_new[2] == 1:
        P_bounds = np.array([0.0, 2.0*np.pi])
        P_cells =  np.array([np.pi])
    else:
        P_bounds = P_ap[Ngz:-Ngz+1]
        P_cells =  P_bp[Ngz:-Ngz]
    
    R_map, T_map, P_map = np.meshgrid(R_cells,T_cells, P_cells, indexing='ij')
    mesh = (R_map, T_map, P_map)


    # Two-dimensional interpolation for the velocities onto cell center
    if Vp is None:
        fspline_r = \
            scintp.RegularGridInterpolator(( R_bounds, T_cells, P_cells), Vr_bounds[:,:-1,:], fill_value=0.0)
        Vr_cells =  fspline_r(mesh)
        fspline_th = \
            scintp.RegularGridInterpolator(( R_cells, T_bounds, P_cells), Vth_bounds[:-1,:,:], fill_value=0.0)
        Vth_cells = fspline_th(mesh)
    
        Vph_cells = np.zeros(naxes_new)
    else:
        fspline_r = \
            scintp.RegularGridInterpolator(( R_bounds, T_cells, P_cells), Vr_bounds[:,:-1,:-1], fill_value=0.0)
        Vr_cells =  fspline_r(mesh)
        fspline_th = \
            scintp.RegularGridInterpolator(( R_cells, T_bounds, P_cells), Vth_bounds[:-1,:,:-1], fill_value=0.0)
        Vth_cells = fspline_th(mesh)
        fspline_p = \
            scintp.RegularGridInterpolator(( R_cells, T_cells, P_bounds), Vp_bounds[:-1,:-1,:], fill_value=0.0)
        Vph_cells = fspline_p(mesh)

    # Determine the outermost index of the R-direction
    rind_max = np.argmax(np.where(R_cells < Rmax))
    naxes_new[0] = rind_max + 1
    R_bounds = R_bounds[:naxes_new[0]+1]
    R_cells =   R_cells[:naxes_new[0]]

    Density_cells = Density_cells[:naxes_new[0],:,:]
    Density_cells = Density_cells.reshape(naxes_new)
    Vr_cells =           Vr_cells[:naxes_new[0],:,:]
    Vr_cells =           Vr_cells.reshape(naxes_new)
    Vth_cells =         Vth_cells[:naxes_new[0],:,:]
    Vth_cells =         Vth_cells.reshape(naxes_new)
    Vph_cells =         Vph_cells[:naxes_new[0],:,:]
    Vph_cells =         Vph_cells.reshape(naxes_new)
    if np.size(T) > 1:
        Temperature_cells = Temperature_cells[:naxes_new[0],:]
        Temperature_cells = Temperature_cells.reshape(naxes_new)
    else:
        Temperature_cells = T * np.ones(naxes_new)
    
    
    
    import astropy.constants as const
    mH = const.u.to('g').value
    pc = const.pc.to('cm').value
    # Return properties 
    return (naxes_new, R_bounds/pc, T_bounds, P_bounds, \
            Density_cells/(2.0*mH)*1e6, Vr_cells*1e-2, Vth_cells*1e-2, Vph_cells*1e-2, Temperature_cells)

# Replicate ZeusTW output to form a full two-quadrant matrix if needed
naxes = [ n1 - 2*Ngz, n2 - 2*Ngz, n3 if n3 == 1 else n3 - 2*Ngz]  # sizes of the active zones from the ZeusTW data
print "[ZeusTW2SPARX] Total number of cells is %d x %d x %d = %d" % (naxes[0], naxes[1], naxes[2], np.prod(naxes))





if ZeusPar.GridType == 'SPH':
    import importlib
    import socket
    hostname = socket.gethostname()
    library_version = ''
    for cluster_name in TIARA_CLUSTER:
        if hostname[:2] == cluster_name:
            library_version = '_'+cluster_name
    pre_unit = "sparx"+library_version+".pre_unit"
    AU2cm = importlib.import_module(pre_unit).AU2cm

    naxes, x1, x2, x3, n_H2, v1, v2, v3, T_k = \
        Sph_MirrorNTrimR( naxes, x1a, x1b, x2a, x2b, x3a, x3b, density, V1, V2, V3, Rmax = CheckAbdSetAtrr('Rmax_AU')*AU2cm, T = temperature )
    print "[ZeusTW2SPARX] Include the inner R = %d AU into SPARX HDF5 table" % CheckAbdSetAtrr('Rmax_AU')
    
    







alpha = CheckAndSetArray('DustAlpha')
Vt = CheckAndSetArray('TurbulentVelocity')

T_d = T_k
kapp_d = CheckAbdSetAtrr('DustKappa')
dust_to_gas = CheckAndSetArray('DustToGas')


# model attribute: molecular name
molec = CheckAbdSetAtrr('MolecularSpecie')
X_mol = CheckAndSetArray('MolecularAbundance')
# model attribute: CMB temperature
T_cmb = CheckAbdSetAtrr('T_cmb')



