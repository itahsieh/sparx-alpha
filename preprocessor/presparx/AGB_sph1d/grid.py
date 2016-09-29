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
Rout = 0.1 
