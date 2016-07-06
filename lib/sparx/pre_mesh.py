#import time
#start = time.time()

from numpy import zeros
from math import pi,ceil,floor
import sys      

#end = time.time()
#print(end - start)

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
                        pass
                elif GridType == 'REC3D':
                        pass
                elif GridType == 'CYL2D':
                        pass
                elif GridType == 'CYL3D':
                        pass
                else:
                        raise RuntimeError('Grid Type not defined : %s' % GridType)
                        sys.exit(2)

        def _gen_mesh_sph1d(self):
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
                        dr = (Rout-Rin)*(stretch_ratio_r-1.)/(stretch_ratio_r**(nr)-1.)
                        for i in range(1,nr+1):
                                R_p[i] =  R_p[i-1] + dr
                                R_c[i-1] = R_p[i-1] + 0.5 * dr
                                dr *= stretch_ratio_r
                else:
                        raise RuntimeError('Spacing Type not defined : %s' % spacing)
                        sys.exit(2)
                theta_p = [0., pi]
                theta_c = [0.5*pi]
                phi_p = [0., 2. * pi]
                phi_c = [pi]
                self.R_p = R_p
                self.R_c = R_c
                self.theta_p = theta_p
                self.theta_c = theta_c
                self.phi_p = phi_p
                self.phi_c = phi_c
                
        def _gen_mesh_sph2d(self):
                gr = self.grid
                
                nr = gr.nr
                nt = gr.nt
                Rin = gr.Rin
                Rout = gr.Rout
                
                R_p = zeros(nr+1)
                R_c = zeros(nr)
                R_p[0] = Rin
                
                theta_p = zeros(nt+1)
                theta_c = zeros(nt)
                theta_p[0] = 0.0
                

                spacing = gr.spacing
                if ( spacing == 'uniform' ):
                        dr = (Rout-Rin)/nr
                        for i in range(1,nr+1):
                                R_p[i] =  R_p[i-1] + dr
                                R_c[i-1] = R_p[i-1] + 0.5 * dr
                                
                        
                        dt = pi / nt
                        for j in range(1,nt+1):
                                theta_p[j] = theta_p[j-1] + dt
                                theta_c[j-1] = theta_p[j-1] + 0.5 * dt

                elif ( spacing == 'stretch'):
                        stretch_ratio_r = gr.stretch_ratio_r
                        dr = (Rout-Rin)*(stretch_ratio_r-1.)/(stretch_ratio_r**(nr)-1.)
                        for i in range(1,nr+1):
                                R_p[i] =  R_p[i-1] + dr
                                R_c[i-1] = R_p[i-1] + 0.5 * dr
                                dr *= stretch_ratio_r
                        
                        # resolution is even
                        if nt % 2 == 0:        
                                dt0 = 0.5 * pi * (stretch_ratio_t - 1.) / (stretch_ratio_t**(nt/2) - 1.)
                                dt = t0 * stretch_ratio_t**(nt/2-1)
                        # resolution is odd
                        else:
                                dt0 = pi * (stretch_ratio_t - 1.) / (2 * stretch_ratio_t**(nt/2) - stretch_ratio_t - 1.)
                                dt = dt0 * stretch_ratio_t**(nt/2)
                        
                        for j in range(1,nt/2+1):
                                theta_p[j] = theta_p[j-1] + dt
                                theta_c[j-1] = theta_p[j-1] + 0.5 * dt
                                dt /= stretch_ratio_t
                        dt = t0
                        for j in range(nt/2+1,nt):
                                theta_p[j] = theta_p[j-1] + dt
                                theta_c[j-1] = theta_p[j-1] + 0.5 * dt
                                dt *= stretch_ratio_t
                        theta_p[nt] = pi
                
                else:
                        raise RuntimeError('Spacing Type not defined : %s' % spacing)
                        sys.exit(2)

                phi_p = [0., 2. * pi]
                phi_c = [pi]
                self.R_p = R_p
                self.R_c = R_c
                self.theta_p = theta_p
                self.theta_c = theta_c
                self.phi_p = phi_p
                self.phi_c = phi_c
                
