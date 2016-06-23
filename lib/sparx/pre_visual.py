import matplotlib.pyplot as plt
from pylab import savefig
from numpy import amax,amin,zeros
from math import pi,sin,cos


class plot:
        def __init__(self,mesh,phys):
                GridType = mesh.grid.GridType
                if   GridType == 'SPH1D':
                        self._plot_sph1d(mesh,phys)
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
                
        def _plot_sph1d(self,mesh,phys):     
                r = mesh.R_c
                
                # Density plot
                plt.subplot(221)
                plt.plot(r,phys.n_H2)
                plt.xscale('log')
                plt.yscale('log')
                plt.ylabel('H2 density (m^-3)')
                
                # Temperature plot
                plt.subplot(222)
                plt.plot(r,phys.T_k)
                plt.xscale('log')
                plt.yscale('linear')
                plt.ylabel('Kinetic Temperature (Kelvin)')
                
                # Velocity plot
                plt.subplot(223)        
                plt.plot(r,phys.V_gas[:,0])
                plt.xscale('log')
                plt.yscale('linear')
                plt.ylabel('Infalling Velocity (m/s)')
                
                # Abundance plot
                plt.subplot(224)
                plt.plot(r,phys.X_mol)
                X_mol_min = amin(phys.X_mol)
                X_mol_max = amax(phys.X_mol)
                if X_mol_min == X_mol_max:
                        plt.axis([0,0,0.1*amin(phys.X_mol),10.*amax(phys.X_mol)])
                        plt.autoscale(True,'x',False)
                plt.xscale('log')
                plt.yscale('log')
                plt.ylabel('Molecular Abundance (Fraction)')
                
                savefig('profile.png')


class vtk_output:
        def __init__(self, mesh, phys):
                GridType = mesh.grid.GridType
                if   GridType == 'SPH1D':
                        self._vtk_sph1d(mesh,phys)
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
                
        def _vtk_sph1d(self,mesh,phys):
                nr = mesh.grid.nr
                nt = 45
                np = 90
                
                # theta grid
                theta_p = zeros(nt+1)
                dtheta = pi / nt
                theta_p[0] = 0.0
                for j in range(1,nt+1):
                        theta_p[j] = theta_p[j-1] + dtheta
                # phi grid
                phi_p = zeros(np+1)
                dphi = 2.* pi / np
                for k in range(1,np+1):
                        phi_p[k] = phi_p[k-1] + dphi
                
                fvtk1=open('visual.vtk', mode = "w")
                print >>fvtk1,'# vtk DataFile Version 3.0'
                print >>fvtk1,'ENV_DISK'
                print >>fvtk1,'ASCII'
                print >>fvtk1,'DATASET STRUCTURED_GRID'
                print >>fvtk1,'DIMENSIONS %(0)d %(1)d %(2)d'%{'0':np+1,'1':nt+1,'2':nr+1}
                print >>fvtk1,'POINTS %(0)d float'%{'0':(nr+1)*(nt+1)*(np+1)}
                for i in range(nr+1):
                  for j in range(nt+1):
                    for k in range(np+1):
                        x = mesh.R_p[i] * sin(theta_p[j]) * cos(phi_p[k])
                        y = mesh.R_p[i] * sin(theta_p[j]) * sin(phi_p[k])
                        z = mesh.R_p[i] * cos(theta_p[j]) 
                        print >>fvtk1,'%(0)e %(1)e %(2)e'%{'0':x,'1':y,'2':z}
                print >>fvtk1,'CELL_DATA %(0)d'%{'0':nr * nt * np}
                print >>fvtk1,'SCALARS density float 1'
                print >>fvtk1,'LOOKUP_TABLE default'
                for i in range(nr):
                  for j in range(nt):
                    for k in range(np):
                        print >>fvtk1,'%(0)8.2e'%{'0':phys.n_H2[i]},
                    print >>fvtk1
                print >>fvtk1,'SCALARS temperature float 1'
                print >>fvtk1,'LOOKUP_TABLE default'
                for i in range(nr):
                  for j in range(nt):
                    for k in range(np):
                        print >>fvtk1,'%(0)7.1f'%{'0':phys.T_k[i]},
                    print >>fvtk1
                print >>fvtk1,'VECTORS velocity float'
                for i in range(nr):
                  for j in range(nt):
                    for k in range(np):
                        Vr = phys.V_gas[i][0]
                        Vx = Vr * sin(theta_p[j]) * cos(phi_p[k])
                        Vy = Vr * sin(theta_p[j]) * sin(phi_p[k])
                        Vz = Vr * cos(theta_p[j]) 
                        print >>fvtk1,'%(0)7.2e %(1)7.2e %(2)7.2e'%{'0':Vx,'1':Vy,'2':Vz}
                fvtk1.close()
                

                
                