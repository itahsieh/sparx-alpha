# Model Type : Function / Constant / TABLE / ZEUS 
ModelType = 'Function'

# reference radius
r_ref = 0.01
# reference H2 number density (m^-3)
n_h2_ref = 1e10
# reference velocity (m/s)
V_ref = -200.0

# Gas Density (number/m^3)
def DensityFunc1D(r):
        n_h2 = n_h2_ref * ( r / r_ref )**-2.
        return n_h2
        
# Temperature (Kelvin)
def TgasFunc1D(r):
        Tk = 10.
        return Tk

# Velocity (m/s)
def VeloFunc1D(r):
        Vr = V_ref * ( r / r_ref )**-0.5
        return Vr
        
# turbulent speed (m/s)
def Vt(r):
        Vt = 200.
        return Vt


# Molecular data
# Molecule
molec = 'hco+'

# Molecular Abundance (fraction)
def MolecAbdFunc1D(r):
        X_mol = 1e-9
        return X_mol

# CMB temperature (Kelvin, outer B.C.)
T_cmb = 2.73

# inner Boundary condition
T_in = 0.0

# Enable / disable dust emission
dust = 1

# gas-to-dust ratio
def DustToGasFunc1D(r):
        dust_to_gas = 0.01
        return dust_to_gas

# Dust Temperature (Kelvin)
def TdustFunc1D(r):
        Td = TgasFunc1D(r)
        return Td

# dust kappa
def kappa_d_Func1D(r):
        kapp_d = 'table,jena_thin_e5'
        #kappa_d = 'powerlaw, 1.874e+12, 2.300e-02,-2.0'
        return kapp_d

