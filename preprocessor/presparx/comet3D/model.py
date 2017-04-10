# Model Type : Function / Constant / TABLE / ZEUS 
ModelType = 'Function'

# Molecule
molec = 'hcn@hfs'

# inner Boundary condition
T_in = 190.

# CMB temperature (Kelvin, outer B.C.)
T_cmb = 2.73

# Outer Source 
from math import acos
nray = 1000
angle_of_view = 2. * acos(1. - 2.* 1./nray)
# Surface Temperature, theta, phi, angle of view
OuterSource = [[5500., 0., 0., angle_of_view]]


# Mass Loss Rate (Msun yr-1 to H2 number per second)
#from sparx_tc.pre_unit import Msun2kg, mH2, yr2sec
#massloss = (1e-5 * Msun2kg / mH2) / yr2sec 

# Radius of the star
from sparx_tc.pre_unit import m2pc
R_star = 2.0e4 * m2pc

# velocity of the wind
V_wind = 1.2 * 1e3 # (ms^-1)


from sympy import pi, exp
from sparx_tc.pre_unit import pc2m
# Gas Density (number/m^3)
def Density1D(r):
    factor = r / R_star
    n_H2 = 1.5e30 / ( 4.0 * pi * V_wind * (r*r*pc2m**2) ) * exp(-(r*pc2m - 2.0e4) / 8.0e7)
    return n_H2

class model:
    def __init__(self,r,theta,phi):
        # Gas Density (number/m^3)
        self.n_H2 = Density1D(r)
        
        self._Tgas1D(r)
        
        self.V_cen = [ V_wind, 0., 0.]
        # turbulent speed (m/s)
        self.Vt = 800.
        
        self._MolecAbd1D(r)
        
        # fraction of para-H2
        self.X_pH2 = 0.25
        
        # fraction of ortho-H2
        self.X_oH2 = 0.75
        
        # gas-to-dust ratio
        self.dust_to_gas = 0.01
        
        # dust temperature
        self.T_d = self.T_k
        
        # dust kappa
        self.kapp_d = 'table,jena_thin_e5'

    # Temperature on surface of the star (Kelvin)
    def _Tgas1D(self,r):
        #T_in = 0.
        # power law for temperature profile
        #P_temp = 0.0
        #factor = r / R_star
        self.T_k = T_in

    # Molecular Abundance (fraction)
    def _MolecAbd1D(self,r):
        self.X_mol = 1e-2 

