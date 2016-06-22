# Model Type : Function / Constant / TABLE / ZEUS 
ModelType = 'Function'

# reference radius
r_ref = 0.01
# reference H2 number density (m^-3)
n_h2_ref = 1e12
# reference velocity (m/s)
V_ref = -200.0

# Gas Density (number/m^3)
def DensityFunc1D(r):
        n_h2 = n_h2_ref * ( r / r_ref )**-1.5
        return n_h2
        
# Temperature (Kelvin)
def TgasFunc1D(r):
        Tk = 10.
        return Tk

# Dust Temperature (Kelvin)
def TdustFunc1D(r):
        Td = TgasFunc1D(r)
        return Td

# Velocity (m/s)
def VeloFunc1D(r):
        Vr = V_ref * ( r / r_ref )**-0.5
        return Vr
        
# Molecular Abundance (fraction)
molec = 'hco+'
def MolecAbdFunc1D(r):
        X_mol = 1e-9
        return X_mol

# CMB temperature (Kelvin, outer B.C.)
Tcmb = 2.73

# gas-to-dust ratio
gas_to_dust = 100.

# dust kappa
kappa_d = 'table,jena_thin_e5'
#kappa_d = 'powerlaw, 1.874e+12, 2.300e-02,-2.0'

# turbulent speed (m/s)
Vt = 200.


