#import time
#start = time.time()

from numpy import zeros
from math import pi,ceil,floor
import sys      

#end = time.time()
#print(end - start)


class from_dataset:
    class grid(object):
        pass
    
    def __init__(self, converter):
        self.grid.naxes = converter.naxes
        self.grid.GridType = converter.GridType
        self.grid.x1 = converter.x1
        self.grid.x2 = converter.x2
        self.grid.x3 = converter.x3
        self._gen_mesh()
        
    def _gen_mesh(self):
        GridType = self.grid.GridType

        if GridType == 'SPH3D':
            self._gen_mesh_sph3d()
        elif GridType == 'REC3D':
            pass
        elif GridType == 'CYL3D':
            pass
        else:
            raise RuntimeError('Grid Type not defined : %s' % GridType)
            sys.exit(2)
            
    def _gen_mesh_sph3d(self):
        n = self.grid.naxes
        self.grid.Rin = self.grid.x1[0]
        self.grid.Rout = self.grid.x1[-1]
        self.grid.nr = n[0]
        self.grid.nt = n[1]
        self.grid.np = n[2]
        
        self.R_p        = self.grid.x1
        self.theta_p    = self.grid.x2
        self.phi_p      = self.grid.x3
        
        naxes = self.grid.naxes
        self.R_c        = zeros(naxes[0])
        self.theta_c    = zeros(naxes[1])
        self.phi_c      = zeros(naxes[2])
        
        for i in range(naxes[0]):
            self.R_c[i]     = 0.5 * ( self.R_p[i] + self.R_p[i+1] )
        for j in range(naxes[1]):
            self.theta_c[j] = 0.5 * ( self.theta_p[j] + self.theta_p[j+1] )
        for k in range(naxes[2]):
            self.phi_c[k]   = 0.5 * ( self.phi_p[k] + self.phi_p[k+1] )
        
        
class mesh:
    def __init__(self, grid):
        self.grid = grid
        self._gen_mesh()
            
    def _gen_mesh(self):
        GridType = self.grid.GridType

        if   GridType == 'SPH1D':
            self._gen_mesh_sph1d()
        elif GridType == 'SPH2D':
            self._gen_mesh_sph2d()
        elif GridType == 'SPH3D':
            self._gen_mesh_sph3d()
        elif GridType == 'REC3D':
            pass
        elif GridType == 'CYL2D':
            self._gen_mesh_cyl2d()
        elif GridType == 'CYL3D':
            pass
        else:
            raise RuntimeError('Grid Type not defined : %s' % GridType)
            sys.exit(2)

    def _gen_mesh_sph1d(self):
        self._sph_r_gridding()
        
        theta_p = [0., pi]
        theta_c = [0.5*pi]
        self.theta_p = theta_p
        self.theta_c = theta_c
                
        phi_p = [0., 2. * pi]
        phi_c = [pi]
        self.phi_p = phi_p
        self.phi_c = phi_c
            
    def _gen_mesh_sph2d(self):
        self._sph_r_gridding()
        self._sph_theta_gridding()
        
        phi_p = [0., 2. * pi]
        phi_c = [pi]
        self.phi_p = phi_p
        self.phi_c = phi_c
            
    def _gen_mesh_sph3d(self):
        self._sph_r_gridding()
        self._sph_theta_gridding()
        self._sph_phi_gridding()
            
       
    def _sph_r_gridding(self):
        gr = self.grid
        
        nr = gr.nr
        Rin = gr.Rin
        Rout = gr.Rout
        
        R_p = zeros(nr+1)
        R_c = zeros(nr)
        R_p[0] = Rin
        
        spacing = gr.spacing
        if ( spacing == 'uniform' ):
            dr = (Rout-Rin)/nr
            for i in range(1,nr+1):
                R_p[i] =  R_p[i-1] + dr
                R_c[i-1] = R_p[i-1] + 0.5 * dr
        elif ( spacing == 'stretch'):
            stretch_ratio_r = gr.stretch_ratio_r

            if stretch_ratio_r == 1.0:
                dr = (Rout-Rin)/nr
            else:
                dr = (Rout-Rin)*(stretch_ratio_r-1.)/(stretch_ratio_r**(nr)-1.)

            for i in range(1,nr+1):
                R_p[i] =  R_p[i-1] + dr
                R_c[i-1] = R_p[i-1] + 0.5 * dr
                dr *= stretch_ratio_r
        elif ( spacing == 'user_defined' ):
            for i in range(1,nr):
                R_p[i] = 0.5 * (gr.grid[i-1] + gr.grid[i])
            R_p[nr] = Rout
            for i in range(0,nr):    
                R_c[i] = gr.grid[i]
        else:
            raise RuntimeError('Spacing Type not defined : %s' % spacing)
            sys.exit(2)
       
        self.R_p = R_p
        self.R_c = R_c
        
        
    def _sph_theta_gridding(self):
        gr = self.grid
        
        nt = gr.nt
        
        theta_p = zeros(nt+1)
        theta_c = zeros(nt)
        theta_p[0] = 0.0
        
        spacing = gr.spacing
        if ( spacing == 'uniform' ):            
            dt = pi / nt
            for j in range(nt):
                theta_p[j+1] = theta_p[j] + dt
                theta_c[j] = theta_p[j] + 0.5 * dt
        elif ( spacing == 'stretch'):
            stretch_ratio_t = gr.stretch_ratio_t
            # resolution is even
            if stretch_ratio_t == 1.0:
                dt0 = pi / nt
                dt = dt0
            else:
                # resolution is even
                if nt % 2 == 0:        
                    dt0 = 0.5 * pi * (stretch_ratio_t - 1.) / (stretch_ratio_t**(nt/2) - 1.)
                    dt = dt0 * stretch_ratio_t**(nt/2-1)
                # resolution is odd
                else:
                    dt0 = pi * (stretch_ratio_t - 1.) / (2 * stretch_ratio_t**(nt/2) - stretch_ratio_t - 1.)
                    dt = dt0 * stretch_ratio_t**(nt/2)
            
            # north semi-sphere
            for j in range(nt/2):
                theta_p[j+1] = theta_p[j] + dt
                theta_c[j] = theta_p[j] + 0.5 * dt
                dt /= stretch_ratio_t
            # south semi-sphere
            dt = dt0
            for j in range(nt/2,nt):
                theta_p[j+1] = theta_p[j] + dt
                theta_c[j] = theta_p[j] + 0.5 * dt
                dt *= stretch_ratio_t
        else:
            raise RuntimeError('Spacing Type not defined : %s' % spacing)
            sys.exit(2)
            
        self.theta_p = theta_p
        self.theta_c = theta_c
        
    def _sph_phi_gridding(self):
        gr = self.grid
        
        np = gr.np
        
        phi_p = zeros(np+1)
        phi_c = zeros(np)
        phi_p[0] = 0.0
        
        spacing = gr.spacing
        if ( spacing == 'uniform' ):            
            dp = 2. * pi / np
            for k in range(np):
                phi_p[k+1] = phi_p[k] + dp
                phi_c[k] = phi_p[k] + 0.5 * dp
        
        elif ( spacing == 'stretch'):
            stretch_ratio_p = gr.stretch_ratio_p
            # resolution is even
            if stretch_ratio_p == 1.0:
                dp = 2. * pi / np
            else:
                # resolution is even
                if np % 2 == 0:        
                    dp = pi * (stretch_ratio_p - 1.) / (stretch_ratio_p**(np/2) - 1.)
                # resolution is odd
                else:
                    dp = 2. * pi * (stretch_ratio_p - 1.) / (2 * stretch_ratio_p**(np/2) - stretch_ratio_p - 1.)
            dp0 = dp
            
            # back semi-sphere
            for k in range(np/2):
                phi_p[k+1] = phi_p[k] + dp
                phi_c[k] = phi_p[k] + 0.5 * dp
                dt *= stretch_ratio_p
            # front semi-sphere
            dp = dp0
            for j in range(np/2,np):
                phi_p[k+1] = phi_p[k] + dp
                phi_c[k] = phi_p[k] + 0.5 * dp
                dp /= stretch_ratio_p
        else:
            raise RuntimeError('Spacing Type not defined : %s' % spacing)
            sys.exit(2)
            
        self.phi_p = phi_p
        self.phi_c = phi_c
            
    def _gen_mesh_cyl2d(self):
            gr = self.grid
            
            nrc = gr.nrc
            nz = gr.nz
            Rc_in = gr.Rc_in
            Rc_out = gr.Rc_out
            z_max = gr.z_max
            
            Rc_p = zeros(nrc+1)
            Rc_c = zeros(nrc)
            Rc_p[0] = Rc_in
            
            z_p = zeros(nz+1)
            z_c = zeros(nz)
            z_p[0] = - z_max
            

            spacing = gr.spacing
            if ( spacing == 'uniform' ):
                    drc = (Rc_out-Rc_in)/nrc
                    for i in range(1,nrc+1):
                            Rc_p[i] =  Rc_p[i-1] + drc
                            Rc_c[i-1] = Rc_p[i-1] + 0.5 * drc
                    
                    dz = pi / nt
                    for j in range(1,nz+1):
                            z_p[j] = z_p[j-1] + dz
                            z_c[j-1] = z_p[j-1] + 0.5 * dz

            elif ( spacing == 'stretch'):
                    stretch_ratio_rc = gr.stretch_ratio_rc
                    drc = (Rc_out-Rc_in)*(stretch_ratio_rc-1.)/(stretch_ratio_rc**(nrc)-1.)
                    for i in range(1,nrc+1):
                            Rc_p[i] =  Rc_p[i-1] + drc
                            Rc_c[i-1] = Rc_p[i-1] + 0.5 * drc
                            drc *= stretch_ratio_rc
                    
                    stretch_ratio_z = gr.stretch_ratio_z
                    # resolution is even
                    if nz % 2 == 0:        
                            dz0 = z_max * (stretch_ratio_z - 1.) / (stretch_ratio_z**(nz/2) - 1.)
                            dz = dz0 * stretch_ratio_z**(nz/2-1)
                    # resolution is odd
                    else:
                            dz0 = 2.0 * z_max * (stretch_ratio_z - 1.) / (2 * stretch_ratio_z**(nz/2) - stretch_ratio_z - 1.)
                            dz = dz0 * stretch_ratio_z**(nz/2)
                    
                    # north semi-sphere
                    for j in range(1,nz/2+1):
                            z_p[j] = z_p[j-1] + dz
                            z_c[j-1] = z_p[j-1] + 0.5 * dz
                            dz /= stretch_ratio_z
                    # south semi-sphere
                    dz = dz0
                    for j in range(nz/2+1,nz+1):
                            z_p[j] = z_p[j-1] + dz
                            z_c[j-1] = z_p[j-1] + 0.5 * dz
                            dz *= stretch_ratio_z
            
            else:
                    raise RuntimeError('Spacing Type not defined : %s' % spacing)
                    sys.exit(2)

            phi_p = [0., 2. * pi]
            phi_c = [pi]
            self.Rc_p = Rc_p
            self.Rc_c = Rc_c
            self.phi_p = phi_p
            self.phi_c = phi_c
            self.z_p = z_p
            self.z_c = z_c
