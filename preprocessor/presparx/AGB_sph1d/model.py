# Model Type : Function / Constant / TABLE / ZEUS 
ModelType = 'Function'

# Molecule
molec = '13co@xpol'

# inner Boundary condition
T_in = 1000.

# CMB temperature (Kelvin, outer B.C.)
T_cmb = 2.73


from sparx_xl.pre_unit import Msun2kg, mH2, yr2sec, AU2pc, m2pc, pc2m

# Mass Loss Rate (Msun yr-1 to H2 number per second)
massloss = (1e-5 * Msun2kg / mH2) / yr2sec 

# Radius of the star
R_star = 1. * AU2pc

# velocity of the wind
V_wind = 10. * 1e3 # (ms^-1)

# Gas Density (number/m^3)
def Density1D(r):
        factor = r / R_star
        from math import pi
        n_H2 = massloss / ( 4.0 * pi * V_wind * ( R_star * R_star * pc2m * pc2m) ) * factor**(-2)
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
                T_in = 1000.
                # power law for temperature profile
                P_temp = 0.0
                factor = r / R_star
                self.T_k = T_in * factor**(-P_temp)

        # Molecular Abundance (fraction)
        def _MolecAbd1D(self,r):
                # molecular abundance at inner boundary
                X0 = 3e-5
                self.X_mol = X0 






