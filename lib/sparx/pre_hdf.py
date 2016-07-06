from tables import * 


# Define a user record to characterize some kind of particles
class Particle(IsDescription):
        LEVEL   = Int32Col(pos=0)
        POS     = Int64Col(pos=1)
        X_max   = Float64Col(shape=3,pos=2)
        X_min   = Float64Col(shape=3,pos=3)
        X_cen   = Float64Col(shape=3,pos=4)
        n_H2    = Float64Col(pos=5)
        T_k     = Float64Col(pos=6)
        X_pH2   = Float64Col(pos=7)
        X_oH2   = Float64Col(pos=8)
        X_e     = Float64Col(pos=9)
        X_H     = Float64Col(pos=10)
        X_He    = Float64Col(pos=11)
        V_t     = Float64Col(pos=12)
        V_cen   = FloatCol(shape=3,pos=13)
        NCHILDREN=Int64Col(pos=14)
        NAXES   = Int64Col(shape=3,pos=15)

class Particle_molec(IsDescription):
        X_mol   = Float64Col(pos=0)

class Particle_dust(IsDescription):
        T_d     = Float64Col(pos=0)
        kapp_d  = StringCol(itemsize=64,pos=1)
        dust_to_gas = Float64Col(pos=2)
        
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
        del table.attrs.NROWS
        
def DelAttrs_molec(table):
        del table.attrs.FIELD_0_FILL 
        del table.attrs.NROWS
        
def DelAttrs_dust(table):
        del table.attrs.FIELD_0_FILL 
        del table.attrs.FIELD_1_FILL
        del table.attrs.FIELD_2_FILL
        del table.attrs.NROWS        

class export:
        def __init__(self,mesh,phys,FileName):
                h5file = openFile(FileName, mode = "w", title = "SPARX MODEL FILE")
                
                # delete attributes
                h5file.delNodeAttr("/", "TITLE", name=None)
                h5file.delNodeAttr("/", "CLASS", name=None)
                h5file.delNodeAttr("/", "VERSION", name=None)
                h5file.delNodeAttr("/", "PYTABLES_FORMAT_VERSION", name=None)
                
                # create attributes
                h5file.setNodeAttr("/", "format", "SPARX format v3", name=None)
                h5file.setNodeAttr("/", "molec", phys.molec, name=None)
                h5file.setNodeAttr("/", "pops", 0, name=None)
                h5file.setNodeAttr("/", "T_cmb", phys.T_cmb, name=None)
                h5file.setNodeAttr("/", "T_in", phys.T_in, name=None)
                h5file.setNodeAttr("/", "dust", phys.dust, name=None)
                
                GridType = mesh.grid.GridType
                if GridType == 'SPH1D':
                        h5file.setNodeAttr("/", "geom", "sph1d", name=None)
                        self._export_sph1d(mesh,phys,h5file)
                if GridType == 'SPH2D':
                        h5file.setNodeAttr("/", "geom", "sph3d", name=None)
                        self._export_sph2d(mesh,phys,h5file)
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
                        particle['X_max']       = [ mesh.R_p[i+1], mesh.theta_p[nt], mesh.phi_p[np] ]
                        particle['X_min']       = [ mesh.R_p[i]  , mesh.theta_p[0] , mesh.phi_p[0]  ]
                        particle['X_cen']       = [ mesh.R_c[i]  , mesh.theta_c[0] , mesh.phi_c[0]  ]
                        particle['n_H2']        = phys.n_H2[i]
                        particle['V_cen']       = phys.V_gas[i]
                        particle['T_k']         = phys.T_k[i]
                        particle['V_t']         = phys.Vt[i]
                        # Insert a new particle record
                        particle.append()
                table.flush()
                DelAttrs(table)
                
                if phys.molec:
                        # Create GRID table
                        table = h5file.createTable("/", 'molec', Particle_molec, "molecular table")
                        particle = table.row
                        for i in range(nr):
                                particle['X_mol']       = phys.X_mol[i]
                                particle.append()
                        table.flush()
                        DelAttrs_molec(table)
                
                if phys.dust:
                        # Create DUST table
                        table = h5file.createTable("/", 'DUST', Particle_dust, "dust table")
                        particle = table.row
                        for i in range(nr):
                                particle['T_d']         = phys.T_d[i]
                                particle['kapp_d']      = phys.kapp_d[i]
                                particle['dust_to_gas']      = phys.dust_to_gas[i]
                                particle.append()
                        table.flush()
                        DelAttrs_dust(table)
                

        def _export_sph2d(sellf,mesh,phys,h5file):
                pass
               