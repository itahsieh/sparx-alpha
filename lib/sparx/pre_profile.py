from numpy import zeros
from math import pi,cos,sqrt
from pre_unit import *

class profile:
        def __init__(self, mesh, model):
                md = self.model = model

                self.molec = md.molec
                self.T_cmb = md.T_cmb
                self.T_in = md.T_in                
                
                # accumulated mass
                self.mass = 0.0
                # accumulated volume
                self.volume = 0.0
                # Max Velocity Dispersion to Vt
                self.MVD2Vt = 0.0 
                
                ModelType = self.model.ModelType
                if ModelType == 'Function':
                        self._MappingFunction(mesh)
                elif ModelType == 'Constant':
                        pass
                elif ModelType == 'user_defined':
                        self._MappingUserDefined(mesh)
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
                
                # accumulated mass
                self.mass = 0.0
                # accumulated volume
                self.volume = 0.0
                # Max Velocity Dispersion to Vt
                self.MVD2Vt = 0.0 
                # MVD2Vt index
                self.MVD2Vt_index = 0
                
                GridType = gr.GridType
                if GridType == 'SPH1D':
                        self._MappingFunction_sph1d(mesh)
                        self._Mass_VeloDisp_sph1d(mesh)
                elif GridType == 'SPH2D':
                        self._MappingFunction_sph2d(mesh)
                        self._Mass_VeloDisp_sph2d(mesh)
                elif GridType == 'SPH3D':
                        pass
                elif GridType == 'CYL2D':
                        self._MappingFunction_cyl2d(mesh)
                        self._Mass_VeloDisp_cyl2d(mesh)
                elif GridType == 'CYL3D':
                        pass
                else:
                        raise RuntimeError('Grid Type not defined : %s' % GridType) 
                
                self.mass *= kg2Msun # Msun
                
                
        def _MappingFunction_sph1d(self,mesh):
                md = self.model
                
                nr = mesh.grid.nr
                self.n_H2 = zeros(nr)
                self.T_k = zeros(nr)
                self.V_gas = zeros((nr,3))
                self.Vt = zeros(nr)
                for i in range(nr):
                        r    = mesh.R_c[i]
                        phys = md.model(r)
                        self.n_H2[i]         = phys.n_H2
                        self.T_k[i]          = phys.T_k
                        self.V_gas[i]        = [phys.Vr, 0., 0.]
                        self.Vt[i]           = phys.Vt
                
                if md.molec:
                        self.X_mol = zeros(nr)
                        for i in range(nr):
                                self.X_mol[i] = phys.X_mol
                
                if hasattr(phys, 'X_pH2'):
                        self.X_pH2 = zeros(nr)
                        for i in range(nr):
                                self.X_pH2[i] = phys.X_pH2
                if hasattr(phys, 'X_oH2'):
                        self.X_oH2 = zeros(nr)
                        for i in range(nr):
                                self.X_oH2[i] = phys.X_oH2
                                
                if hasattr(phys, 'T_d'):
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
                        theta   = mesh.theta_c[j]

                        phys = md.model(r,theta)
                        self.n_H2[i,j]  = phys.n_H2
                        self.T_k[i,j]   = phys.T_k
                        self.V_gas[i,j] = phys.V_cen
                        self.Vt[i,j]    = phys.Vt
                
                if md.molec:
                        self.X_mol = zeros((nr,nt))
                        for i in range(nr):
                            for j in range(nt):
                                self.X_mol[i,j] = phys.X_mol
                
                if hasattr(phys, 'X_pH2'):
                        self.X_pH2 = zeros((nr,nt))
                        for i in range(nr):
                            for j in range(nt):
                                self.X_pH2[i,j] = phys.X_pH2
                if hasattr(phys, 'X_oH2'):
                        self.X_oH2 = zeros((nr,nt))
                        for i in range(nr):
                            for j in range(nt):
                                self.X_oH2[i,j] = phys.X_oH2
                                
                if hasattr(phys, 'T_d'):
                        self.T_d = zeros((nr,nt))
                        self.dust_to_gas = zeros((nr,nt))
                        self.kapp_d = []
                        for i in range(nr):
                            for j in range(nt):
                                self.T_d[i,j]           = phys.T_d
                                self.dust_to_gas[i,j]   = phys.dust_to_gas
                                self.kapp_d.append(phys.kapp_d)
                
                if hasattr(phys, 'B_field'):
                        self.B_field = zeros((nr,nt,3))
                        self.alpha = zeros((nr,nt))
                        self.z     = zeros((nr,nt))
                        for i in range(nr):
                            for j in range(nt):
                                self.B_field[i,j] = phys.B_field
                                self.alpha[i,j]   = phys.alpha
                                self.z[i,j]       = phys.z

        def _MappingFunction_cyl2d(self,mesh):
                gr = mesh.grid
                md = self.model
                
                nrc = gr.nrc
                nz = gr.nz
                self.n_H2 = zeros((nrc,nz))
                self.T_k = zeros((nrc,nz))
                self.V_gas = zeros((nrc,nz,3))
                self.Vt = zeros((nrc,nz))
                for i in range(nrc):
                    for j in range(nz):
                        rc      = mesh.Rc_c[i]
                        z       = mesh.z_c[j]

                        phys = md.model(rc,z)
                        self.n_H2[i,j]  = phys.n_H2
                        self.T_k[i,j]   = phys.T_k
                        self.V_gas[i,j] = phys.V_cen
                        self.Vt[i,j]    = phys.Vt
                        
                
                if md.molec:
                        self.X_mol = zeros((nrc,nz))
                        for i in range(nrc):
                            for j in range(nz):
                                self.X_mol[i,j] = phys.X_mol

                if hasattr(phys, 'X_pH2'):
                        self.X_pH2 = zeros((nrc,nz))
                        for i in range(nrc):
                            for j in range(nz):
                                self.X_pH2[i,j] = phys.X_pH2
                if hasattr(phys, 'X_oH2'):
                        self.X_oH2 = zeros((nrc,nz))
                        for i in range(nrc):
                            for j in range(nz):
                                self.X_oH2[i,j] = phys.X_oH2
                                
                if hasattr(phys, 'T_d'):
                        self.T_d = zeros((nrc,nz))
                        self.dust_to_gas = zeros((nrc,nz))
                        self.kapp_d = []
                        for i in range(nrc):
                            for j in range(nz):
                                self.T_d[i,j]           = phys.T_d
                                self.dust_to_gas[i,j]   = phys.dust_to_gas
                                self.kapp_d.append(phys.kapp_d)
                
                if hasattr(phys, 'B_field'):
                        self.B_field = zeros((nrc,nz,3))
                        self.alpha = zeros((nrc,nz))
                        self.z     = zeros((nrc,nz))
                        for i in range(nrc):
                            for j in range(nz):
                                self.B_field[i,j] = phys.B_field
                                self.alpha[i,j]   = phys.alpha
                                self.z[i,j]       = phys.z



        def _MappingUserDefined(self, mesh):
                

                
                gr = mesh.grid
                GridType = gr.GridType
                if GridType == 'SPH1D':
                        self._MappingUserDefined_sph1d(mesh)
                        self._Mass_VeloDisp_sph1d(mesh)
                elif GridType == 'SPH2D':
                        pass
                elif GridType == 'SPH3D':
                        pass
                elif GridType == 'CYL2D':
                        pass
                else:
                        raise RuntimeError('Grid Type not defined : %s' % GridType) 
                
                self.mass *= kg2Msun # Msun

        def _MappingUserDefined_sph1d(self,mesh):
                phys = self.model
                nr = mesh.grid.nr
                self.n_H2 = zeros(nr)
                self.T_k = zeros(nr)
                self.V_gas = zeros((nr,3))
                self.Vt = zeros(nr)
                for i in range(nr):                        
                        self.n_H2[i]         = phys.n_H2[i]
                        self.T_k[i]          = phys.T_k[i]
                        self.V_gas[i]        = [ phys.Vr[i], 0., 0.]
                        self.Vt[i]           = phys.Vt[i]

                if phys.molec:
                        self.X_mol = phys.X_mol

                if hasattr(phys, 'T_d'):
                        self.T_d = phys.T_d
                        self.dust_to_gas = phys.dust_to_gas
                        self.kapp_d = []
                        for i in range(nr):
                                self.kapp_d.append(phys.kapp_d)


        def _Mass_VeloDisp_sph1d(self,mesh):                
                nr = mesh.grid.nr
                for i in range(nr):
                        r       = mesh.R_c[i]
                        r_in    = mesh.R_p[i]
                        r_out   = mesh.R_p[i+1]
                        
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

        def _Mass_VeloDisp_sph2d(self,mesh):     
                gr = mesh.grid                
                nr = gr.nr
                nt = gr.nt

                for i in range(nr):
                    for j in range(nt):
                        r_in    = mesh.R_p[i]
                        r_out   = mesh.R_p[i+1]
                        theta_in        = mesh.theta_p[j] 
                        theta_out       = mesh.theta_p[j+1]

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
                
                
        def _Mass_VeloDisp_cyl2d(self,mesh):
                gr = mesh.grid
                nrc = gr.nrc
                nz = gr.nz
                
                for i in range(nrc):
                    for j in range(nz):
                        rc_in   = mesh.Rc_p[i]
                        rc_out  = mesh.Rc_p[i+1]
                        z_min   = mesh.z_p[j] 
                        z_max   = mesh.z_p[j+1]

                        # volume
                        dVolume = 4. * pi * ( rc_out**2 - rc_in**2 ) * ( z_max -z_min ) # pc^3
                        
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
                                Vrc_rc = self.V_gas[i+1,j,0] - self.V_gas[i,j,0]
                        elif    i == nrc-1:
                                Vrc_rc = self.V_gas[i,j,0] - self.V_gas[i-1,j,0]
                        else:
                                Vrc_rc = max( abs(self.V_gas[i+1,j,0] - self.V_gas[i,j,0]), abs(self.V_gas[i,j,0] - self.V_gas[i-1,j,0]) )
                                
                        if      j == 0:
                                Vz_z = self.V_gas[i,j+1,1] - self.V_gas[i,j,2]
                        elif    j == nz-1:
                                Vz_z = self.V_gas[i,j,1] - self.V_gas[i,j-1,2]
                        else:
                                Vz_z = max( abs(self.V_gas[i,j+1,1] - self.V_gas[i,j,2]), abs(self.V_gas[i,j,1] - self.V_gas[i,j-1,1]) )
                        
                        VeloDispersion = sqrt(Vrc_rc * Vrc_rc + Vz_z * Vz_z)

                        # Velocity Dispersion to Vt
                        VeloDispersion2Vt = VeloDispersion / self.Vt[i,j]
                        # update Maximum VD2Vt
                        if VeloDispersion2Vt > self.MVD2Vt :
                                self.MVD2Vt = VeloDispersion2Vt
                                self.MVD2Vt_index = [i,j]