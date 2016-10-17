from sympy import Symbol,integrate
from math import pi
from pre_unit import *

class error:
        def __init__(self,mesh,phys):
                ModelType = phys.model.ModelType
                if ModelType == 'Function':
                        self.calc_exact_mass(mesh,phys)
                        print 'Analytical Mass : %f Msun' %self.exact_mass
                        MassError = 100. * ( phys.mass - self.exact_mass ) / self.exact_mass
                        print 'Mass Error of the Gridding : %f %%' %MassError
                else:
                        print 'Total mass = %f MSun' %phys.mass
                print 'Largest Velocity Dispersion to Turbulent Velocity : %f' %phys.MVD2Vt
                print 'Largest Dispersion occurs at spatial index = %d' %phys.MVD2Vt_index

        
        def calc_exact_mass(self,mesh,phys):
                gr = mesh.grid
                md = phys.model
                GridType = gr.GridType
                if GridType =='SPH1D':
                        r_min = gr.Rin
                        r_max = gr.Rout
                        r = Symbol('r')
                        self.exact_mass = integrate(md.Density1D(r) * 4.*pi*r**2,(r,r_min,r_max))
                        self.exact_mass *= volume_pc2m * MeanMolecularMass * kg2Msun
