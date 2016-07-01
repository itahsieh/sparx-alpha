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
                self.Tcmb = md.Tcmb
                self.T_in = md.T_in
                
                # accumulated mass
                mass = 0.0
                # accumulated volume
                volume = 0.0
                # Max Velocity Dispersion to Vt
                MVD2Vt = 0.0 
                
                GridType = gr.GridType
                if GridType == 'SPH1D':
                        _MappingFunction_sph1d(mesh)

                elif GridType == SPH2D:
                        pass
                
                elif GridType == SPH3D:
                        pass
                
                self.n_H2 = n_H2
                self.T_k = T_k
                self.T_d = T_d
                self.V_gas = V_gas
                self.Vt = Vt

                if gr.molec:
                        self.X_mol = X_mol
                        
                if gr.dust:
                        self.T_d = T_d
                        self.dust_to_gas = dust_to_gas
                        self.kapp_d = kapp_d
                
                mass *= kg2Msun # Msun
                self.mass = mass
                self.volume = volume
                self.MVD2Vt = MVD2Vt
                self.MVD2Vt_index = MVD2Vt_index
                
                
        def _MappingFunction_sph1d(self,mesh):
                gr = mesh.grid
                md = self.model
                
                nr = gr.nr
                n_H2 = zeros(nr)
                T_k = zeros(nr)
                V_gas = zeros((nr,3))
                Vt = zeros(nr)
                for i in range(nr):
                        r = mesh.R_c[i]
                        n_H2[i]         = md.DensityFunc1D(r)
                        T_k[i]          = md.TgasFunc1D(r)
                        V_gas[i]        = [md.VeloFunc1D(r), 0., 0.]
                        Vt[i]           = md.Vt(r)
                        
                        # volume
                        dVolume = 4. / 3. * pi * (mesh.R_p[i+1]**3-mesh.R_p[i]**3) # pc^3
                        volume += dVolume # pc^3
                        dVolume *= volume_pc2m # m^3
                        
                        # delta mass
                        dMass = n_H2[i] * dVolume * MeanMolecularMass # kg
                        # accumulated mass
                        mass += dMass #kg
                        
                        # max delta V (m/s)
                        if i == 0:
                                VeloDispersion = abs(V_gas[i+1,0] - V_gas[i,0])
                        elif i == nr-1:
                                VeloDispersion = abs(V_gas[i,0] - V_gas[i-1,0])
                        else:
                                VeloDispersion = max(abs(V_gas[i,0] - V_gas[i-1,0]),abs(V_gas[i+1,0] - V_gas[i,0]))
                        # Velocity Dispersion to Vt
                        VeloDispersion2Vt = VeloDispersion / Vt[i]
                        # update Maximum VD2Vt
                        if VeloDispersion2Vt > MVD2Vt :
                                MVD2Vt = VeloDispersion2Vt
                                MVD2Vt_index = i
                
                if gr.molec:
                        X_mol = zeros(nr)
                        for i in range(nr):
                                r = mesh.R_c[i]
                                X_mol[i] = md.MolecAbdFunc1D(r)
                                
                if gr.dust:
                        T_d = zeros(nr)
                        dust_to_gas = zeros(nr)
                        kapp_d = zeros(nr)
                        for i in range(nr):
                                r = mesh.R_c[i]
                                T_d[i]          = md.TdustFunc1D(r)
                                dust_to_gas[i]  = md.DustToGasFunc1D(r)
                                kapp_d[i]       = md.kappa_d_Func1D(r)


def calc_exact_mass(grid,model):
        r_min = grid.Rin
        r_max = grid.Rout
        r = Symbol('r')
        mass = integrate(model.DensityFunc1D(r) * 4.*pi*r**2,(r,r_min,r_max))
        mass *= volume_pc2m * MeanMolecularMass *kg2Msun
        return mass