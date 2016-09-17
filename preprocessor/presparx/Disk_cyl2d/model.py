####################################################
# The model of disk and envelope (Keto&Zhang 2010) #
####################################################

# Model Type : Function / Constant / TABLE / ZEUS 
ModelType = 'Function'

# Molecule
molec = 'hco+'

# Enable / disable dust emission
dust = 1

# CMB temperature (Kelvin, outer B.C.)
T_cmb = 2.73

# inner Boundary condition
T_in = 0.0

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

rho_d0 = Ap*rho_e0
H0 = 0.01*Rt
Tt = 30000.

# unit converter
km2m=1e3
pc2m=30.857e15
pc2km=30.857e12
# mean mollecular mass (kg)
mean_molecular_mass = 2*1.67*1.67e-27 # kg



from math import sin,cos,pi,exp,sqrt,isnan,acos
Mar = rho_e0*4.*pi*Rd*Rd*vk*(mean_molecular_mass*pc2m**2*km2m) # mass accretion rate (kg/s)

G = 4.302e-3 # gravitational constant (pc Msun^-1 (km/s)^2)
sigma = 5.67037321e-8 # Stefan-Boltzmann constant (W m^-2 K^-4)

pp = 0.0
qq = 0.0
def CubicEq(x):
        global pp,qq
        return x*x*x + pp*x + qq

from scipy.optimize import brentq

class model:        
        def __init__(self, Rc, z):
                self.Rc         = Rc
                self.z          = z
                self.r          = sqrt( Rc*Rc + z*z )
                self.sin_theta  = Rc / self.r
                self.cos_theta  = z / self.r
                self.theta      = acos( self.cos_theta )
                
                global pp,qq
                pp = self.r / Rd - 1.
                qq = -self.cos_theta * self.r / Rd
                self.cos_theta0 = brentq(CubicEq, -1.,1.)

                self._DensitySph2D
                self._TgasSph2D
                self._VeloSph2
                self._VtSph2D
                self._MolecAbdSph2D
                self._DustToGasSph2D
                self._TdustSph2D
                self._kappa_d_Sph2D
        

        
        # Gas Density (number/m^3)
        def _DensitySph2D(self):
                # envelope density
                if (env == 1):
                        self.n_H2_env = rho_e0 * ((self.r/Rd)**(-1.5)) * ((1. + self.cos_theta / self.cos_theta0)**(-0.5)) * ((1 + ((self.r/Rd)**(-1)) * (3 * self.cos_theta0**2 - 1.0))**(-1))
                else:
                        self.n_H2_env = 0.0
                
                # disk density
                if (self.r<=Rd and disk==1):
                        rho_0 = rho_d0*(Rd/self.Rc)**2.25
                        H=H0*(self.Rc/Rt)**1.2
                        self.n_H2_disc = rho_0 * exp(-(self.r*self.r-self.Rc*self.Rc)/(2.*H*H))
                else:
                        self.n_H2_disc = 0.0
                
                # total density
                self.n_H2 = self.n_H2_env + self.n_H2_disc

        
        # Temperature (Kelvin)
        def _TgasSph2D(self):
                if ( self.n_H2 != 0.0 ):
                        T_env = Tt*(Rt/(2.*self.r))**(2./(4+p))
                        T_disc = BT * ( (3.*G*Mt*Mar/(4.*pi * pc2km*pc2km * (self.Rc**3) * sigma)) * (1.-sqrt(Rt/self.Rc)) )**0.25
                        self.T_k = ( self.n_H2_disc*T_disc + self.n_H2_env*T_env ) / self.n_H2
                else:
                        self.T_k = 0.0

        # Velocity (m/s)
        def _VeloSph2D(self):
                # Keplerian velocity (km/s)
                Vkep = sqrt(G*Mt/self.r)
                # disk velocity (km/s)
                Vp_disc = sqrt(G*Mt/self.Rc)
                
                Vr_env = -Vkep * sqrt( 1. + self.cos_theta / self.cos_theta0 )
                
                if self.sin_theta == 0. or self.cos_theta0 == 0.:
                        print self.r, self.theta, self.sin_theta, self.cos_theta0
                        import sys
                        sys.exit(0)
                Vt_env = Vkep * ( (self.cos_theta0 - self.cos_theta / self.sin_theta ) * sqrt( 1. + self.cos_theta / self.cos_theta0 )
                Vp_env = Vkep * ( sqrt( 1. - self.cos_theta0 * self.cos_theta0) / self.sin_theta ) * sqrt( 1. + self.cos_theta / self.cos_theta0 )
                
                if self.n_H2 != 0.:
                        Vr = (self.n_H2_env * Vr_env) / self.n_H2
                        Vt = (self.n_H2_env * Vt_env) / self.n_H2
                        Vp = (self.n_H2_env * Vp_env + self.n_H2_disc * Vp_disc) / self.n_H2
                else:
                        Vr = 0.0
                        Vt = 0.0
                        Vp = 0.0
                
                self.V_cen = [km2m*Vr,km2m*Vt,km2m*Vp]

        # turbulent speed (m/s)
        def _VtSph2D(self):
                self.Vt = 200.

        # Molecular Abundance (fraction)
        def _MolecAbdSph2D(self):
                self.X_mol = 1e-9
3
        # gas-to-dust ratio
        def _DustToGasSph2D(self):
                self.dust_to_gas = 0.01

        # Dust Temperature (Kelvin)
        def _TdustSph2D(self):
                self.T_d = self.T_k

        # dust kappa
        def _kappa_d_Sph2D(self,r,theta):
                self.kapp_d = 'table,jena_thin_e5'
                #kappa_d = 'powerlaw, 1.874e+12, 2.300e-02,-2.0'

