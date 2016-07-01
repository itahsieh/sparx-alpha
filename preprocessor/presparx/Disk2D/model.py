# Model Type : Function / Constant / TABLE / ZEUS 
ModelType = 'Function'

# enable/disable disk model
disk = 1

# enable/disable envelope model
env = 1

# Physical parameter
#parameters from Keto&Zhang
rho_e0 = 7.9e4 *1e6 # Envelope density at Rd (m^-3)
Rd = 6900./206260. # AU to pc
Ap = 5.1
Mt = 10.7 # stellar mass (Msun)
Rt = 26.*0.0046491/206260
p = -1
BT = 15.0
vk = 1.2 # Keplerian velocity (km/s)

from scipy.optimize import brentq

class Disk_Env:
        def __init__(self,r,theta):
                self.n_H2       = _DensityFuncSph2D(r,theta)
                self.T_k        = _TgasFuncSph2D(r,theta)
        
        def _CubicEq(xx):
                global pp,qq
                return xx*xx*xx+pp*xx+qq
        
        # Gas Density (number/m^3)
        def _DensityFuncSph2D(r,theta):
                if (env == 1):
                        pp=r/Rd-1.
                        qq=-cos(theta)*r/Rd
                        cos_theta0 = brentq(_CubicEq, -1.,1.)
                        
                        density_env = rho_e0 * ((r/Rd)**(-1.5))*((1+cos(theta)/cos_theta0)**(-0.5))*((1 + ((r/Rd)**(-1)) * (3*cos_theta0**2-1.0))**(-1))
                        mass_env += density_env*volume
                if (r<=Rd and disk==1):
                        rho_0 = rho_d0*(Rd/R)**2.25
                        H=H0*(R/Rt)**1.2
                        density_disc = rho_0 * exp(-(r*r-R*R)/(2.*H*H))
        
                n_h2 = density_disc + density_env
                return n_h2
        
        # Temperature (Kelvin)
        def TgasFuncSph2D(r,theta):
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

