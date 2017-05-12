# Model Type : Function / Constant / TABLE / ZEUS 
ModelType = 'Function'

# Molecule
molec = 'p-h2co'
# CMB temperature (Kelvin, outer B.C.)
T_cmb = 2.73
# inner Boundary condition
T_in = 0.0
# reference radius
r_ref = 0.005
# reference H2 number density (m^-3)
n_H2_ref = 6.3e12
# reference velocity (m/s)
V_ref = -1100

abd1 = 0.6667e-8
abd2 = 0.3333e-9

# Gas Density (number/m^3)
def Density1D(r):
        if ( r < 23.9/2.06e5 ):
                n_H2 = 5e15
        else:
                n_H2 = n_H2_ref * ( r / r_ref )**-1.8
        return n_H2

class model:
        def __init__(self,r):
                # Gas Density (number/m^3)
                self.n_H2 = Density1D(r)
                # Velocity (m/s)
                self.Vr = V_ref * ( r / r_ref )**-0.5
                # Temperature (Kelvin)
                self.T_k = 71.*(r/0.00032467532)**-1 +60.7*(r/0.00032467532)**-0.4
                # turbulent speed (m/s)
                self.Vt = 400.
                # Molecular Abundance (fraction)
                if( self.T_k > 100. ):
                        self.X_mol = abd1 # Fraction
                        if( self.T_k > 250. ):
                                self.T_k = 250.
                else:
                        self.X_mol = abd2 # Fraction
                # dust-to-gas ratio
                self.dust_to_gas = 0.01
                # Dust Temperature (Kelvin)
                self.T_d = self.T_k
                # dust kappa
                self.kapp_d = 'powerlaw, 1.874e+12, 2.300e-02, +2.000e+00'

