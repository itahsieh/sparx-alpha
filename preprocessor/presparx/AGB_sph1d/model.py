# Model Type : Function / Constant / TABLE / ZEUS 
ModelType = 'Function'

# Molecule
molec = '13co@xpol'

# inner Boundary condition
T_in = 1000.

# CMB temperature (Kelvin, outer B.C.)
T_cmb = 2.73


# Mass Loss Rate (Msun yr-1 to H2 number per second)
from pre_unit import Msun2kg, mH2, yr2sec
from math import pi
massloss = (1e-5 * Msun2kg / mH2) / yr2sec 

# Radius of the star
from pre_unit import AU2pc
R_star = 1. * AU2pc

# velocity of the wind
V_wind = 10. * 1e3 # (ms^-1)

# Gas Density (number/m^3)
def Density1D(r):
        factor = r / R_star
        n_H2 = massloss / ( 4.0 * pi * V_wind * ( R_star * R_star ) ) * factor**(-2)
        return n_H2

class model:
        def __init__(self,r):
                # Gas Density (number/m^3)
                self.n_H2 = Density1D(r)
                
                self._Tgas1D(r)
                
                self.Vr = V_wind
                # turbulent speed (m/s)
                self.Vt = 800.
                
                self._MolecAbd1D(r)
                
                # gas-to-dust ratio
                self.dust_to_gas = 0.01
                
                # dust temperature
                self.T_d = self.T_k
                
                # dust kappa
                self.kapp_d = 'table,jena_thin_e5'

        # Temperature on surface of the star (Kelvin)
        def _Tgas1D(self,r):
                T_in = 1000.
                # power law for temperature profile
                P_temp = 0.0
                factor = r / R_star
                self.T_k = T_in * factor**(-P_temp)

        # Molecular Abundance (fraction)
        def _MolecAbd1D(self,r):
                from math import exp,log
                # molecular abundance at inner boundary
                X0 = 3e-5
                # radius of molecular photodissociation
                from pre_unit import m2pc
                rp = 1.9e12 * m2pc
                # power law for abundance profile
                P_abund = 2.5
                self.X_mol = X0 * exp( log(2.0)*(r/rp)**(P_abund) )






