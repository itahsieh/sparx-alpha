# Grid Type : SPH1D / SPH3D / CYL3D /REC3D / REC1D
GridType = 'SPH1D'

# spacing type : uniform / stretch -> stretch_ratio_r
spacing = 'uniform'

# resolution of the radius
nr = 64

# inner radius (pc)
from model import R_star
Rin = R_star

# outer radius (pc)
from sparx_tc.pre_unit import m2pc
Rout = 3.0e8 * m2pc

