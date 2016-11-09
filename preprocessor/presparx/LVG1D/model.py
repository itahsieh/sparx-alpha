# Model Type : Function / Constant / TABLE / ZEUS 
ModelType = 'Function'

# Molecule
molec = 'hco+'

# CMB temperature (Kelvin, outer B.C.)
T_cmb = 2.73

# inner Boundary condition

# reference radius
r_ref = 0.01
# reference H2 number density (m^-3)
n_H2_ref = 1e10
# reference velocity (m/s)
V_ref = 200.0

# Gas Density (number/m^3)
def Density1D(r):
        n_H2 = n_H2_ref
        return n_H2

class model:
        def __init__(self,r):
                self._Density1D(r)
                self._Tgas1D(r)
                self._Velo1D(r)
                self._Vt1D(r)
                self._MolecAbd1D(r)

                
        # Gas Density (number/m^3)
        def _Density1D(self,r):
                self.n_H2 = Density1D(r)
        
        # Temperature (Kelvin)
        def _Tgas1D(self,r):
                self.T_k = 10.

        # Velocity (m/s)
        def _Velo1D(self,r):
                self.Vr = V_ref * ( r / r_ref )
                
        # turbulent speed (m/s)
        def _Vt1D(self,r):
                self.Vt = 100.

        # Molecular Abundance (fraction)
        def _MolecAbd1D(self,r):
                self.X_mol = 1e-9



