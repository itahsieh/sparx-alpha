# Model Type : Function / Constant / TABLE / ZEUS 
ModelType = 'Function'

# Molecule
molec = 'hco+'

# CMB temperature (Kelvin, outer B.C.)
T_cmb = 2.73

# inner Boundary condition
T_in = 0.0

# reference radius
r_ref = 0.01

# reference H2 number density (m^-3)
n_H2_ref = 1e10

# reference velocity (m/s)
V_ref = -200.0

# Gas Density (number/m^3)
def Density1D(r):
        n_H2 = n_H2_ref * ( r / r_ref )**-2.
        return n_H2

class model:
        def __init__(self,r):
                self._Density1D(r)
                self._Tgas1D(r)
                self._Velo1D(r)
                self._Vt1D(r)
                self._MolecAbd1D(r)
                self._DustToGas1D(r)
                
                # Dust Temperature (Kelvin)
                self.T_d = self.T_k
                # dust kappa
                self.kapp_d = 'table,jena_thin_e5'
                
        # Gas Density (number/m^3)
        def _Density1D(self,r):
                self.n_H2 = Density1D(r)
        
        # Temperature (Kelvin)
        def _Tgas1D(self,r):
                self.T_k = 10.

        # Velocity (m/s)
        def _Velo1D(self,r):
                self.Vr = V_ref * ( r / r_ref )**-0.5
                
        # turbulent speed (m/s)
        def _Vt1D(self,r):
                self.Vt = 200.

        # Molecular Abundance (fraction)
        def _MolecAbd1D(self,r):
                self.X_mol = 1e-9

        # gas-to-dust ratio
        def _DustToGas1D(self,r):
                self.dust_to_gas = 0.01



