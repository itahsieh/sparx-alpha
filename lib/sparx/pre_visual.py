import matplotlib
import matplotlib.pyplot as plt
from pylab import savefig
from numpy import amax,amin,zeros
import numpy as np
from math import pi,sin,cos


class plot:
        def __init__(self,mesh,phys):
                GridType = mesh.grid.GridType
                if   GridType == 'SPH1D':
                        self._plot_sph1d(mesh,phys)
                elif GridType == 'SPH2D':
                        self._plot_sph2d(mesh,phys)
                elif GridType == 'SPH3D':
                        pass
                elif GridType == 'REC3D':
                        pass
                elif GridType == 'CYL2D':
                        self._plot_cyl2d(mesh,phys)
                elif GridType == 'CYL3D':
                        pass
                
                filename='profile.png'
                savefig(filename)
                print filename,'generated'
                
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
                
                

        def _plot_sph2d(self,mesh,phys):     
                r = mesh.R_c
                theta = mesh.theta_c
                
                theta_grid, r_grid, = np.meshgrid(theta, r)
                
                x = r_grid * np.sin(theta_grid)
                y = r_grid * np.cos(theta_grid)
                
                # Density plot
                plt.subplot(231, aspect=1)
                plt.pcolormesh(x, y, phys.n_H2, norm=matplotlib.colors.LogNorm())
                plt.title('H2 density (m^-3)')
                plt.colorbar()
                
                # Temperature plot
                plt.subplot(232, aspect=1)
                plt.pcolormesh(x, y, phys.T_k, norm=matplotlib.colors.LogNorm())
                plt.title('Kinetic Temperature (Kelvin)')
                plt.colorbar()
                
                # Velocity plot
                plt.subplot(233, aspect=1)
                plt.pcolormesh(x, y, phys.X_mol, norm=matplotlib.colors.LogNorm(vmin=0.1*amin(phys.X_mol), vmax=10.*amax(phys.X_mol)))
                plt.title('Molecular Abundance (Fraction)')
                plt.colorbar()
                
                # Velocity plot
                V_gas_max = amax(phys.V_gas[:,:,:])
                V_gas_min = amin(phys.V_gas[:,:,:])
                absmax = 0.05 * max(abs(V_gas_max),abs(V_gas_min))
                
                plt.subplot(234, aspect=1)
                plt.pcolormesh(x, y, phys.V_gas[:,:,0], norm=matplotlib.colors.Normalize(vmin=-absmax, vmax=absmax))
                plt.title('V_r (m/s)')
                plt.colorbar()
                
                plt.subplot(235, aspect=1)
                plt.pcolormesh(x, y, phys.V_gas[:,:,1], norm=matplotlib.colors.Normalize(vmin=-absmax, vmax=absmax))
                plt.title('V_theta (m/s)')
                plt.colorbar()
                
                plt.subplot(236, aspect=1)
                plt.pcolormesh(x, y, phys.V_gas[:,:,2], norm=matplotlib.colors.Normalize(vmin=-absmax, vmax=absmax))
                plt.title('V_phi (m/s)')
                plt.colorbar()
                
                
        def _plot_cyl2d(self,mesh,phys):     
                rc = mesh.Rc_c
                z = mesh.z_c
                
                y, x, = np.meshgrid(z, rc)

                # Density plot
                plt.subplot(231, aspect=1)
                plt.pcolormesh(x, y, phys.n_H2, norm=matplotlib.colors.LogNorm())
                plt.title('H2 density (m^-3)')
                plt.colorbar()
                
                # Temperature plot
                plt.subplot(232, aspect=1)
                plt.pcolormesh(x, y, phys.T_k, norm=matplotlib.colors.LogNorm())
                plt.title('Kinetic Temperature (Kelvin)')
                plt.colorbar()
                
                # Velocity plot
                plt.subplot(233, aspect=1)
                plt.pcolormesh(x, y, phys.X_mol, norm=matplotlib.colors.LogNorm(vmin=0.1*amin(phys.X_mol), vmax=10.*amax(phys.X_mol)))
                plt.title('Molecular Abundance (Fraction)')
                plt.colorbar()
                
                # Velocity plot
                V_gas_max = amax(phys.V_gas[:,:,:])
                V_gas_min = amin(phys.V_gas[:,:,:])
                absmax = 0.05 * max(abs(V_gas_max),abs(V_gas_min))
                
                plt.subplot(234, aspect=1)
                plt.pcolormesh(x, y, phys.V_gas[:,:,0], norm=matplotlib.colors.Normalize(vmin=-absmax, vmax=absmax))
                plt.title('V_rc (m/s)')
                plt.colorbar()
                
                plt.subplot(235, aspect=1)
                plt.pcolormesh(x, y, phys.V_gas[:,:,1], norm=matplotlib.colors.Normalize(vmin=-absmax, vmax=absmax))
                plt.title('V_z (m/s)')
                plt.colorbar()
                
                plt.subplot(236, aspect=1)
                plt.pcolormesh(x, y, phys.V_gas[:,:,2], norm=matplotlib.colors.Normalize(vmin=-absmax, vmax=absmax))
                plt.title('V_phi (m/s)')
                plt.colorbar()
                
                
class vtk_output:
        def __init__(self, mesh, phys):
                self.filename='visual.vtk'
                
                GridType = mesh.grid.GridType
                if   GridType == 'SPH1D':
                        self._vtk_sph1d(mesh,phys)
                elif GridType == 'SPH2D':
                        self._vtk_sph2d(mesh,phys)
                elif GridType == 'SPH3D':
                        pass
                elif GridType == 'REC3D':
                        pass
                elif GridType == 'CYL2D':
                        self._vtk_cyl2d(mesh,phys)
                elif GridType == 'CYL3D':
                        pass
                
                print self.filename,'generated'
                
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
                
                fvtk1=open(self.filename, mode = "w")
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
                
        def _vtk_sph2d(self,mesh,phys):
                nr = mesh.grid.nr
                nt = mesh.grid.nt
                np = 90
                
                # phi grid
                phi_p = zeros(np+1)
                dphi = 2.* pi / np
                for k in range(1,np+1):
                        phi_p[k] = phi_p[k-1] + dphi
                
                fvtk1=open(self.filename, mode = "w")
                print >>fvtk1,'# vtk DataFile Version 3.0'
                print >>fvtk1,'ENV_DISK'
                print >>fvtk1,'ASCII'
                print >>fvtk1,'DATASET STRUCTURED_GRID'
                print >>fvtk1,'DIMENSIONS %(0)d %(1)d %(2)d'%{'0':np+1,'1':nt+1,'2':nr+1}
                print >>fvtk1,'POINTS %(0)d float'%{'0':(nr+1)*(nt+1)*(np+1)}
                for i in range(nr+1):
                  for j in range(nt+1):
                    for k in range(np+1):
                        r       = mesh.R_p[i]
                        theta   = mesh.theta_p[j]
                        phi     = phi_p[k]
                        x = r * sin(theta) * cos(phi)
                        y = r * sin(theta) * sin(phi)
                        z = r * cos(theta) 
                        print >>fvtk1,'%(0)e %(1)e %(2)e'%{'0':x,'1':y,'2':z}
                print >>fvtk1,'CELL_DATA %(0)d'%{'0':nr * nt * np}
                print >>fvtk1,'SCALARS density float 1'
                print >>fvtk1,'LOOKUP_TABLE default'
                for i in range(nr):
                  for j in range(nt):
                    for k in range(np):
                        print >>fvtk1,'%(0)8.2e'%{'0':phys.n_H2[i,j]},
                    print >>fvtk1
                print >>fvtk1,'SCALARS temperature float 1'
                print >>fvtk1,'LOOKUP_TABLE default'
                for i in range(nr):
                  for j in range(nt):
                    for k in range(np):
                        print >>fvtk1,'%(0)7.1f'%{'0':phys.T_k[i,j]},
                    print >>fvtk1
                print >>fvtk1,'VECTORS velocity float'
                for i in range(nr):
                  for j in range(nt):
                    for k in range(np):
                        Vr = phys.V_gas[i,j][0]
                        Vt = phys.V_gas[i,j][1]
                        Vp = phys.V_gas[i,j][2]
                        theta   = mesh.theta_p[j]
                        phi     = phi_p[k]
                        Vx = Vr * sin(theta) * cos(phi) + Vt * cos(theta) * cos(phi) - Vp * sin(phi)
                        Vy = Vr * sin(theta) * sin(phi) + Vt * cos(theta) * sin(phi) + Vp * cos(phi)
                        Vz = Vr * cos(theta)            - Vt * sin(theta)
                        print >>fvtk1,'%(0)7.2e %(1)7.2e %(2)7.2e'%{'0':Vx,'1':Vy,'2':Vz}
                fvtk1.close()
                
        def _vtk_cyl2d(self,mesh,phys):
                nrc = mesh.grid.nrc
                np = 90
                nz = mesh.grid.nz
                
                # phi grid
                phi_p = zeros(np+1)
                phi_c = zeros(np)
                phi_p[0] = 0.
                dphi = 2.* pi / np
                for j in range(1,np+1):
                        phi_p[j] = phi_p[j-1] + dphi
                        phi_c[j-1] = 0.5 * (phi_p[j-1] + phi_p[j])
                
                fvtk1=open(self.filename, mode = "w")
                print >>fvtk1,'# vtk DataFile Version 3.0'
                print >>fvtk1,'ENV_DISK'
                print >>fvtk1,'ASCII'
                print >>fvtk1,'DATASET STRUCTURED_GRID'
                print >>fvtk1,'DIMENSIONS %(0)d %(1)d %(2)d'%{'0':nz+1,'1':np+1,'2':nrc+1}
                print >>fvtk1,'POINTS %(0)d float'%{'0':(nrc+1)*(np+1)*(nz+1)}
                for i in range(nrc+1):
                  for j in range(np+1):
                    for k in range(nz+1):
                        rc  = mesh.Rc_p[i]
                        phi = phi_p[j]
                        z   = mesh.z_p[k]
                        
                        x = rc * cos(phi)
                        y = rc * sin(phi)
                        print >>fvtk1,'%(0)e %(1)e %(2)e'%{'0':x,'1':y,'2':z}
                print >>fvtk1,'CELL_DATA %(0)d'%{'0':nrc * np * nz}
                print >>fvtk1,'SCALARS density float 1'
                print >>fvtk1,'LOOKUP_TABLE default'
                for i in range(nrc):
                  for j in range(np):
                    for k in range(nz):
                        print >>fvtk1,'%(0)8.2e'%{'0':phys.n_H2[i,k]},
                    print >>fvtk1
                print >>fvtk1,'SCALARS temperature float 1'
                print >>fvtk1,'LOOKUP_TABLE default'
                for i in range(nrc):
                  for j in range(np):
                    for k in range(nz):
                        print >>fvtk1,'%(0)7.1f'%{'0':phys.T_k[i,k]},
                    print >>fvtk1
                print >>fvtk1,'VECTORS velocity float'
                for i in range(nrc):
                  for j in range(np):
                    for k in range(nz):
                        Vrc = phys.V_gas[i,k][0]
                        Vp  = phys.V_gas[i,k][1]
                        Vz  = phys.V_gas[i,k][2]
                        phi = phi_c[j]
                        Vx  = Vrc * cos(phi) - Vp * sin(phi)
                        Vy  = Vrc * sin(phi) + Vp * cos(phi)
                        print >>fvtk1,'%(0)7.2e %(1)7.2e %(2)7.2e'%{'0':Vx,'1':Vy,'2':Vz}
                fvtk1.close()