from numpy import zeros
from math import pi

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

                elif GridType == SPH2D:
                        pass
                
                elif GridType == SPH3D:
                        pass

                
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
                        r = mesh.R_c[i]
                        phys = md.model(r)
                        self.n_H2[i]         = phys.n_H2
                        self.T_k[i]          = phys.T_k
                        self.V_gas[i]        = [phys.Vr, 0., 0.]
                        self.Vt[i]           = phys.Vt
                        
                        # volume
                        dVolume = 4. / 3. * pi * (mesh.R_p[i+1]**3-mesh.R_p[i]**3) # pc^3
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


def calc_exact_mass(grid,model):
        GridType = grid.GridType
        if GridType =='SPH1D':
                r_min = grid.Rin
                r_max = grid.Rout
                r = Symbol('r')
                mass = integrate(model.Density1D(r) * 4.*pi*r**2,(r,r_min,r_max))
                mass *= volume_pc2m * MeanMolecularMass *kg2Msun
                return mass