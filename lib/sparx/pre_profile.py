from numpy import zeros
from math import pi

pc2m = 30.857e15 # m per pc
volume_pc2m =  pc2m**3
mean_molecular_mass = 2. * 1.67 * 1.67e-27 # kg
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
                self.gas_to_dust = md.gas_to_dust 
                self.kappa_d = md.kappa_d
                
                mass = 0.0
                # Max Velocity Dispersion to Vt
                MVD2Vt = 0.0 
                
                GridType = gr.GridType
                if GridType == 'SPH1D':
                        nr = gr.nr
                        n_H2 = zeros(nr)
                        T_k = zeros(nr)
                        T_d = zeros(nr)
                        V_gas = zeros((nr,3))
                        X_mol = zeros(nr)
                        Vt = zeros(nr)
                        for i in range(nr):
                                r = mesh.R_c[i]
                                n_H2[i]         = md.DensityFunc1D(r)
                                T_k[i]          = md.TgasFunc1D(r)
                                T_d[i]          = md.TdustFunc1D(r)
                                V_gas[i]        = [md.VeloFunc1D(r), 0., 0.]
                                X_mol[i]        = md.MolecAbdFunc1D(r)
                                Vt[i]           = md.Vt(r)
                                
                                # volume
                                volume = 4. / 3. * pi * (mesh.R_p[i+1]**3-mesh.R_p[i]**3) # pc^3
                                volume *= volume_pc2m # m^3
                                # delta mass
                                dMass = n_H2[i] * volume * mean_molecular_mass # kg
                                dMass *= kg2Msun # Msun
                                # accumulated mass
                                mass += dMass
                                
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

                elif GridType == SPH2D:
                        pass
                elif GridType == SPH3D:
                        pass
                
                self.n_H2 = n_H2
                self.T_k = T_k
                self.T_d = T_d
                self.V_gas = V_gas
                self.X_mol = X_mol
                
                self.mass = mass
                self.MVD2Vt = MVD2Vt
                self.MVD2Vt_index = MVD2Vt_index
                