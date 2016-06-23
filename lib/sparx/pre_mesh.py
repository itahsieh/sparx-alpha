#import time
#start = time.time()

from numpy import zeros
from math import pi
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
                        pass
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
