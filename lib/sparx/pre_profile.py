from numpy import zeros

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
                GridType = gr.GridType
                if GridType == 'SPH1D':
                        nr = gr.nr
                        n_H2 = zeros(nr)
                        T_k = zeros(nr)
                        T_d = zeros(nr)
                        V_gas = zeros((nr,3))
                        X_mol = zeros(nr)
                        for i in range(nr):
                                r = mesh.R_c[i]
                                n_H2[i]         = md.DensityFunc1D(r)
                                T_k[i]          = md.TgasFunc1D(r)
                                T_d[i]          = md.TdustFunc1D(r)
                                V_gas[i]        = [md.VeloFunc1D(r), 0., 0.]
                                X_mol[i]        = md.MolecAbdFunc1D(r)
                elif GridType == SPH3D:
                        pass
                
                self.n_H2 = n_H2
                self.T_k = T_k
                self.T_d = n_H2
                self.n_H2 = n_H2
                self.n_H2 = n_H2
                