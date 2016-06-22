import matplotlib.pyplot as plt
from numpy import amax,amin

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
                
                plt.show()


class vtk_output:
        def __init__(self, mesh, phys, OutputFile):
                GridType = mesh.grid.GridType
                if   GridType == 'SPH1D':
                        self._vtk_sph1d(mesh,phys,OutputFile)
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
                
                