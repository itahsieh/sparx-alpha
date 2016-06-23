from tables import * 


# Define a user record to characterize some kind of particles
class Particle(IsDescription):
        LEVEL=Int32Col(pos=0)
        POS=Int64Col(pos=1)
        geom=StringCol(itemsize=6,pos=2)
        X_max=Float64Col(shape=3,pos=3)
        X_min=Float64Col(shape=3,pos=4)
        X_cen=Float64Col(shape=3,pos=5)
        n_H2=Float64Col(pos=6)
        T_k=Float64Col(pos=7)
        X_mol=Float64Col(pos=8)
        X_pH2=Float64Col(pos=9)
        X_oH2=Float64Col(pos=10)
        X_e=Float64Col(pos=11)
        X_H=Float64Col(pos=12)
        X_He=Float64Col(pos=13)
        V_t=Float64Col(pos=14)
        V_edge=FloatCol(shape=(6,3),pos=15)
        V_cen=FloatCol(shape=3,pos=16)
        B_cen=FloatCol(shape=3,pos=17)
        ds=FloatCol(pos=18)
        NCHILDREN=Int64Col(pos=19)
        NAXES=Int64Col(shape=3,pos=20)
        T_d=Float64Col(pos=21)
        kapp_d=StringCol(itemsize=64,pos=22)
        T_ff=Float64Col(pos=23)
        kapp_ff=StringCol(itemsize=64,pos=24)
        T_bb=Float64Col(pos=25) 
        
def DelAttrs(table):
        del table.attrs.FIELD_0_FILL 
        del table.attrs.FIELD_1_FILL
        del table.attrs.FIELD_2_FILL
        del table.attrs.FIELD_3_FILL
        del table.attrs.FIELD_4_FILL
        del table.attrs.FIELD_5_FILL
        del table.attrs.FIELD_6_FILL
        del table.attrs.FIELD_7_FILL
        del table.attrs.FIELD_8_FILL
        del table.attrs.FIELD_9_FILL
        del table.attrs.FIELD_10_FILL
        del table.attrs.FIELD_11_FILL
        del table.attrs.FIELD_12_FILL
        del table.attrs.FIELD_13_FILL
        del table.attrs.FIELD_14_FILL
        del table.attrs.FIELD_15_FILL
        del table.attrs.FIELD_16_FILL
        del table.attrs.FIELD_17_FILL
        del table.attrs.FIELD_18_FILL
        del table.attrs.FIELD_19_FILL
        del table.attrs.FIELD_20_FILL
        del table.attrs.FIELD_21_FILL
        del table.attrs.FIELD_22_FILL
        del table.attrs.FIELD_23_FILL
        del table.attrs.FIELD_24_FILL
        del table.attrs.FIELD_25_FILL
        del table.attrs.NROWS

class export:
        def __init__(self,mesh,phys,FileName):
                h5file = openFile(FileName, mode = "w", title = "SPARX MODEL FILE")
                h5file.delNodeAttr("/", "TITLE", name=None)
                h5file.delNodeAttr("/", "CLASS", name=None)
                h5file.delNodeAttr("/", "VERSION", name=None)
                h5file.delNodeAttr("/", "PYTABLES_FORMAT_VERSION", name=None)
                h5file.setNodeAttr("/", "molec", phys.molec, name=None)
                h5file.setNodeAttr("/", "T_cmb", phys.Tcmb, name=None)
                h5file.setNodeAttr("/", "gas_to_dust", phys.gas_to_dust, name=None)
                h5file.setNodeAttr("/", "velfield", "grid ", name=None)
                
                GridType = mesh.grid.GridType
                if GridType == 'SPH1D':
                        self._export_sph1d(mesh,phys,h5file)
                else:
                        pass
                h5file.close()        
        
        def _export_sph1d(sellf,mesh,phys,h5file):
                nr = mesh.grid.nr
                nt = 1
                np = 1
                
                # Create ZONE table
                table = h5file.createTable("/", 'ZONE', Particle, "Grid table")
                particle = table.row
                particle['LEVEL']       = -1
                particle['POS']         = 0
                particle['geom']        = 'sph1d'
                particle['X_max']       = [ mesh.R_p[nr], mesh.theta_p[nt], mesh.phi_p[np] ]
                particle['X_min']       = [ mesh.R_p[0] , mesh.theta_p[0] , mesh.phi_p[0]  ]
                particle['X_cen']       = [ 0.5*( mesh.grid.Rin + mesh.grid.Rout ), mesh.theta_c[0], mesh.phi_c[0] ]
                particle['NCHILDREN']   = nr * nt * np
                particle['NAXES']       = [nr,nt,np]
                #Insert a new particle record
                particle.append()
                table.flush()
                DelAttrs(table)
                
                # Create GRID table
                table = h5file.createTable("/", 'GRID', Particle, "Grid table")
                particle = table.row
                for i in range(nr):
                        particle['LEVEL']       = 0
                        particle['POS']         = i
                        particle['geom']        = 'sph1d'
                        particle['X_max']       = [ mesh.R_p[i+1], mesh.theta_p[nt], mesh.phi_p[np] ]
                        particle['X_min']       = [ mesh.R_p[i]  , mesh.theta_p[0] , mesh.phi_p[0]  ]
                        particle['X_cen']       = [ mesh.R_c[i]  , mesh.theta_c[0] , mesh.phi_c[0]  ]
                        particle['n_H2']        = phys.n_H2[i]
                        particle['V_cen']       = phys.V_gas[i]
                        particle['T_k']         = phys.T_k[i]
                        particle['X_mol']       = phys.X_mol[i]
                        particle['V_t']         = phys.Vt[i]
                        particle['T_d']         = phys.T_d[i]
                        particle['kapp_d']      = phys.kappa_d
                        # Insert a new particle record
                        particle.append()
                table.flush()
                DelAttrs(table)
                
                 
               