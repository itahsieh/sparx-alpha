from numpy import zeros
from math import pi,cos,sqrt

from sympy import Symbol,integrate


pc2m = 30.857e15 # m per pc
volume_pc2m =  pc2m**3

# meam molecular mass, mass / mass_H2 ~ 1.67
MeanMolecularMass = 2. * 1.67 * 1.6726219e-27 # kg
kg2Msun = 1./ 1.98892E+30 # kg

class profile:
        def __init__(self, mesh, model):
                self.model = model
                self._mapping(mesh)
        
        def _mapping(self, mesh):
                ModelType = self.model.ModelType
                if ModelType == 'Function':
                        self._MappingFunction(mesh)
                elif ModelType == 'Constant':
                        pass
                elif ModelType == 'TABLE':
                        pass
                elif ModelType == 'ZEUS':
                        pass
                else:
                        raise RuntimeError('Model Type not defined : %s' % ModelType) 
        
        def _MappingFunction(self, mesh):
                gr = mesh.grid
                md = self.model
                
                self.molec = md.molec
                self.T_cmb = md.T_cmb
                self.T_in = md.T_in
                self.dust = md.dust
                
                
                # accumulated mass
                self.mass = 0.0
                # accumulated volume
                self.volume = 0.0
                # Max Velocity Dispersion to Vt
                self.MVD2Vt = 0.0 
                
                GridType = gr.GridType
                if GridType == 'SPH1D':
                        self._MappingFunction_sph1d(mesh)

                elif GridType == 'SPH2D':
                        self._MappingFunction_sph2d(mesh)
                
                elif GridType == 'SPH3D':
                        pass
                else:
                        raise RuntimeError('Grid Type not defined : %s' % GridType) 
                
                self.mass *= kg2Msun # Msun
                
                
        def _MappingFunction_sph1d(self,mesh):
                gr = mesh.grid
                md = self.model
                
                nr = gr.nr
                self.n_H2 = zeros(nr)
                self.T_k = zeros(nr)
                self.V_gas = zeros((nr,3))
                self.Vt = zeros(nr)
                for i in range(nr):
                        r       = mesh.R_c[i]
                        r_in    = mesh.R_p[i]
                        r_out   = mesh.R_p[i+1]
                        
                        phys = md.model(r)
                        self.n_H2[i]         = phys.n_H2
                        self.T_k[i]          = phys.T_k
                        self.V_gas[i]        = [phys.Vr, 0., 0.]
                        self.Vt[i]           = phys.Vt
                        
                        # volume
                        dVolume = 4. / 3. * pi * (r_out**3-r_in**3) # pc^3
                        self.volume += dVolume # pc^3
                        dVolume *= volume_pc2m # m^3
                        
                        # delta mass
                        dMass = self.n_H2[i] * dVolume * MeanMolecularMass # kg
                        # accumulated mass
                        self.mass += dMass #kg
                        
                        # max delta V (m/s)
                        if i == 0:
                                VeloDispersion = abs(self.V_gas[i+1,0] - self.V_gas[i,0])
                        elif i == nr-1:
                                VeloDispersion = abs(self.V_gas[i,0] - self.V_gas[i-1,0])
                        else:
                                VeloDispersion = max(
                                        abs(self.V_gas[i,0] - self.V_gas[i-1,0]),
                                        abs(self.V_gas[i+1,0] - self.V_gas[i,0])
                                                     )
                        # Velocity Dispersion to Vt
                        VeloDispersion2Vt = VeloDispersion / self.Vt[i]
                        # update Maximum VD2Vt
                        if VeloDispersion2Vt > self.MVD2Vt :
                                self.MVD2Vt = VeloDispersion2Vt
                                self.MVD2Vt_index = i
                
                if md.molec:
                        self.X_mol = zeros(nr)
                        for i in range(nr):
                                self.X_mol[i] = phys.X_mol
                                
                if md.dust:
                        self.T_d = zeros(nr)
                        self.dust_to_gas = zeros(nr)
                        self.kapp_d = []
                        for i in range(nr):
                                self.T_d[i]          = phys.T_d
                                self.dust_to_gas[i]  = phys.dust_to_gas
                                self.kapp_d.append(phys.kapp_d)

        def _MappingFunction_sph2d(self,mesh):
                gr = mesh.grid
                md = self.model
                
                nr = gr.nr
                nt = gr.nt
                self.n_H2 = zeros((nr,nt))
                self.T_k = zeros((nr,nt))
                self.V_gas = zeros((nr,nt,3))
                self.Vt = zeros((nr,nt))
                for i in range(nr):
                    for j in range(nt):
                        r       = mesh.R_c[i]
                        r_in    = mesh.R_p[i]
                        r_out   = mesh.R_p[i+1]
                        theta   = mesh.theta_c[j]
                        theta_in        = mesh.theta_p[j] 
                        theta_out       = mesh.theta_p[j+1]
                        
                        phys = md.model(r,theta)
                        self.n_H2[i,j]  = phys.n_H2
                        self.T_k[i,j]   = phys.T_k
                        self.V_gas[i,j] = phys.V_cen
                        self.Vt[i,j]    = phys.Vt
                        
                        # volume
                        dVolume = 2. / 3. * pi * (r_out**3-r_in**3) * (cos(theta_in)-cos(theta_out)) # pc^3
                        
                        self.volume += dVolume # pc^3
                        dVolume *= volume_pc2m # m^3
                        
                        # delta mass
                        dMass = self.n_H2[i,j] * dVolume * MeanMolecularMass # kg
                        # accumulated mass
                        self.mass += dMass #kg
                        
                        # max delta V (m/s)
                        # no need to concern the velocity deviation along theta & phi
                        # because sparx tracer would take care of it
                        if      i == 0  :
                                Vr_r = self.V_gas[i+1,j,0] - self.V_gas[i,j,0]
                        elif    i == nr-1:
                                Vr_r = self.V_gas[i,j,0] - self.V_gas[i-1,j,0]
                        else:
                                Vr_r = max( abs(self.V_gas[i+1,j,0] - self.V_gas[i,j,0]), abs(self.V_gas[i,j,0] - self.V_gas[i-1,j,0]) )
                                
                        if      j == 0:
                                Vt_t = self.V_gas[i,j+1,1] - self.V_gas[i,j,1]
                        elif    j == nt-1:
                                Vt_t = self.V_gas[i,j,1] - self.V_gas[i,j-1,1]
                        else:
                                Vt_t = max( abs(self.V_gas[i,j+1,1] - self.V_gas[i,j,1]), abs(self.V_gas[i,j,1] - self.V_gas[i,j-1,1]) )
                        
                        VeloDispersion = sqrt(Vr_r * Vr_r + Vt_t * Vt_t)

                        # Velocity Dispersion to Vt
                        VeloDispersion2Vt = VeloDispersion / self.Vt[i,j]
                        # update Maximum VD2Vt
                        if VeloDispersion2Vt > self.MVD2Vt :
                                self.MVD2Vt = VeloDispersion2Vt
                                self.MVD2Vt_index = [i,j]
                
                if md.molec:
                        self.X_mol = zeros((nr,nt))
                        for i in range(nr):
                            for j in range(nt):
                                self.X_mol[i,j] = phys.X_mol
                                
                if md.dust:
                        self.T_d = zeros((nr,nt))
                        self.dust_to_gas = zeros((nr,nt))
                        self.kapp_d = []
                        for i in range(nr):
                            for j in range(nt):
                                self.T_d[i,j]           = phys.T_d
                                self.dust_to_gas[i,j]   = phys.dust_to_gas
                                self.kapp_d.append(phys.kapp_d)

def calc_exact_mass(grid,model):
        GridType = grid.GridType
        if GridType =='SPH1D':
                r_min = grid.Rin
                r_max = grid.Rout
                r = Symbol('r')
                mass = integrate(model.Density1D(r) * 4.*pi*r**2,(r,r_min,r_max))
                mass *= volume_pc2m * MeanMolecularMass * kg2Msun
                return mass