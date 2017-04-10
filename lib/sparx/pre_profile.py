from numpy import zeros
from math import pi,cos,sqrt
from pre_unit import *

class profile:
    def __init__(self, mesh, model):
        md = self.model = model

        self.molec = md.molec
        if hasattr(md, 'T_cmb'):
            self.T_cmb = md.T_cmb
        
        if hasattr(md, 'T_in'):
            self.T_in = md.T_in
        
        if hasattr(md, 'OuterSource'):
            self.OuterSource = md.OuterSource
        
        
        GridType = mesh.grid.GridType
        if GridType == 'SPH1D':
            phys = md.model(1.)
        elif GridType == 'SPH2D' or GridType == 'CYL2D':
            phys = md.model(1.,1.)
        elif GridType == 'SPH3D' or GridType == 'CYL3D' or GridType == 'REC3D':
            phys = md.model(1.,1.,1.)
        
        if hasattr(phys, 'B_cen'):
            self.B_field = 1
        
        if hasattr(phys, 'T_d'):
            self.dust = 1

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
                self._MappingFunction_sph3d(mesh)
                self._Mass_VeloDisp_sph3d(mesh)
        elif GridType == 'CYL2D':
                self._MappingFunction_cyl2d(mesh)
                self._Mass_VeloDisp_cyl2d(mesh)
        elif GridType == 'CYL3D':
                pass
        else:
                raise RuntimeError('Grid Type not defined : %s' % GridType) 
        
        self.mass *= kg2Msun # Msun
    
    def _InitPhys(self, n, phys):
        # n is a tuple
        self.n_H2 = zeros(n)
        self.T_k = zeros(n)
        self.Vt = zeros(n)
        self.V_gas = zeros(n + (3,))
        
        md = self.model
        if md.molec:
            self.X_mol = zeros(n)
            
        if hasattr(phys, 'X_pH2'):
            self.X_pH2 = zeros(n)
        
        if hasattr(phys, 'X_oH2'):
            self.X_oH2 = zeros(n)
        
        if hasattr(phys, 'X_e'):
            self.X_e = zeros(n)
        
        if hasattr(phys, 'X_H'):
            self.X_H = zeros(n)
        
        if hasattr(phys, 'X_He'):
            self.X_He = zeros(n)
        
        if hasattr(phys, 'T_d'):
            self.T_d = zeros(n)
            self.dust_to_gas = zeros(n)
            self.kapp_d = []
    
    def _MappingPhys(self, phys, molec, i):
        self.n_H2[i]         = phys.n_H2
        self.T_k[i]          = phys.T_k
        self.Vt[i]           = phys.Vt
        
        if molec:
            self.X_mol[i] = phys.X_mol
    
        if hasattr(phys, 'X_pH2'):
            self.X_pH2[i] = phys.X_pH2
        
        if hasattr(phys, 'X_oH2'):
            self.X_oH2[i] = phys.X_oH2
        
        if hasattr(phys, 'X_e'):
            self.X_e[i] = phys.X_e
        
        if hasattr(phys, 'X_H'):
            self.X_H[i] = phys.X_H    
        
        if hasattr(phys, 'X_He'):
            self.X_He[i] = phys.X_He    
        
        if hasattr(phys, 'T_d'):
            self.T_d[i]          = phys.T_d
            self.dust_to_gas[i]  = phys.dust_to_gas
            self.kapp_d.append(phys.kapp_d)

        
    def _MappingFunction_sph1d(self,mesh):
            md = self.model
            nr = mesh.grid.nr
            
            r    = mesh.R_c[0]
            phys = md.model(r)
            
            self._InitPhys( (nr,), phys )
                
            for i in range(nr):
                r    = mesh.R_c[i]
                phys = md.model(r)
                
                self._MappingPhys(phys, md.molec, i)
                
                self.V_gas[i]        = [phys.Vr, 0., 0.]

    def _MappingFunction_sph2d(self,mesh):
            gr = mesh.grid
            md = self.model
            
            nr = gr.nr
            nt = gr.nt
            
            r       = mesh.R_c[0]
            theta   = mesh.theta_c[0]
            phys    = md.model(r,theta)
            
            self._InitPhys( (nr,nt), phys )
            
            for i in range(nr):
                for j in range(nt):
                    r       = mesh.R_c[i]
                    theta   = mesh.theta_c[j]
                    phys = md.model(r,theta)
                    
                    self._MappingPhys(phys, md.molec, (i,j))

                    self.V_gas[i,j] = phys.V_cen

    def _MappingFunction_sph3d(self,mesh):
            gr = mesh.grid
            md = self.model
            
            nr = gr.nr
            nt = gr.nt
            np = gr.np
            
            r       = mesh.R_c[0]
            theta   = mesh.theta_c[0]
            phi     = mesh.phi_c[0]
            phys    = md.model(r,theta,phi)
            
            self._InitPhys( (nr,nt,np), phys )
            
            for i in range(nr):
                for j in range(nt):
                    for k in range(np):
                        r       = mesh.R_c[i]
                        theta   = mesh.theta_c[j]
                        phi     = mesh.phi_c[j]
                        phys    = md.model(r,theta,phi)
                        
                        self._MappingPhys(phys, md.molec, (i,j,k))

                        self.V_gas[i,j,k] = phys.V_cen
                        

    def _MappingFunction_cyl2d(self,mesh):
            gr = mesh.grid
            md = self.model
            
            nrc = gr.nrc
            nz = gr.nz
            
            rc      = mesh.Rc_c[0]
            z       = mesh.z_c[0]
            phys    = md.model(rc,z)
                    
            self._InitPhys( (nrc,nz), phys )
            
            for i in range(nrc):
                for j in range(nz):
                    rc      = mesh.Rc_c[i]
                    z       = mesh.z_c[j]
                    phys = md.model(rc,z)
                    
                    self._MappingPhys(phys, md.molec, (i,j))
                    
                    self.V_gas[i,j] = phys.V_cen


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
            
    def _Mass_VeloDisp_sph3d(self,mesh):     
            gr = mesh.grid                
            nr = gr.nr
            nt = gr.nt
            np = gr.np

            for i in range(nr):
                for j in range(nt):
                    for k in range(np):
                        r_in    = mesh.R_p[i]
                        r_out   = mesh.R_p[i+1]
                        theta_in        = mesh.theta_p[j] 
                        theta_out       = mesh.theta_p[j+1]
                        phi_in        = mesh.phi_p[k] 
                        phi_out       = mesh.phi_p[k+1]
                        
                        # volume
                        dVolume = (phi_out - phi_in) / 3. * (r_out**3-r_in**3) * (cos(theta_in)-cos(theta_out)) # pc^3
                        
                        self.volume += dVolume # pc^3
                        dVolume *= volume_pc2m # m^3
                        
                        # delta mass
                        dMass = self.n_H2[i,j,k] * dVolume * MeanMolecularMass # kg
                        # accumulated mass
                        self.mass += dMass #kg
                        
                        # max delta V (m/s)
                        # no need to concern the velocity deviation along theta & phi
                        # because sparx tracer would take care of it
                        if      i == 0  :
                                Vr_r = self.V_gas[i+1,j,k,0] - self.V_gas[i,j,k,0]
                        elif    i == nr-1:
                                Vr_r = self.V_gas[i,j,k,0] - self.V_gas[i-1,j,k,0]
                        else:
                                Vr_r = max( abs(self.V_gas[i+1,j,k,0] - self.V_gas[i,j,k,0]), abs(self.V_gas[i,j,k,0] - self.V_gas[i-1,j,k,0]) )
                                
                        if      j == 0:
                                Vt_t = self.V_gas[i,j+1,k,1] - self.V_gas[i,j,k,1]
                        elif    j == nt-1:
                                Vt_t = self.V_gas[i,j,k,1] - self.V_gas[i,j-1,k,1]
                        else:
                                Vt_t = max( abs(self.V_gas[i,j+1,k,1] - self.V_gas[i,j,k,1]), abs(self.V_gas[i,j,k,1] - self.V_gas[i,j-1,k,1]) )
                        
                        if      k == 0:
                                Vp_p = self.V_gas[i,j,k+1,2] - self.V_gas[i,j,k,2]
                        elif    k == np-1:
                                Vp_p = self.V_gas[i,j,k,2] - self.V_gas[i,j,k-1,2]
                        else:
                                Vp_p = max( abs(self.V_gas[i,j,k+1,2] - self.V_gas[i,j,k,2]), abs(self.V_gas[i,j,k,2] - self.V_gas[i,j,k-1,2]) )
                        
                        VeloDispersion = sqrt(Vr_r * Vr_r + Vt_t * Vt_t + Vp_p*Vp_p)

                        # Velocity Dispersion to Vt
                        VeloDispersion2Vt = VeloDispersion / self.Vt[i,j,k]
                        # update Maximum VD2Vt
                        if VeloDispersion2Vt > self.MVD2Vt :
                                self.MVD2Vt = VeloDispersion2Vt
                                self.MVD2Vt_index = [i,j,k]
 
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
