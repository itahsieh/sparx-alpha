from tables import * 


# Define a user record to characterize some kind of particles
class ParticleZone(IsDescription):
        LEVEL   = Int32Col(pos=0)
        X_max   = Float64Col(shape=3,pos=1)
        X_min   = Float64Col(shape=3,pos=2)
        X_cen   = Float64Col(shape=3,pos=3)
        NCHILDREN=Int64Col(pos=4)
        NAXES   = Int64Col(shape=3,pos=5)

def DelAttrsZone(table):
        del table.attrs.FIELD_0_FILL 
        del table.attrs.FIELD_1_FILL
        del table.attrs.FIELD_2_FILL
        del table.attrs.FIELD_3_FILL
        del table.attrs.FIELD_4_FILL
        del table.attrs.FIELD_5_FILL
        del table.attrs.NROWS

class ParticleGrid(IsDescription):
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

def DelAttrsGrid(table):
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
        
class Particle_molec(IsDescription):
        X_mol   = Float64Col(pos=0)

def DelAttrs_molec(table):
        del table.attrs.FIELD_0_FILL 
        del table.attrs.NROWS
        
class Particle_dust(IsDescription):
        T_d     = Float64Col(pos=0)
        kapp_d  = StringCol(itemsize=64,pos=1)
        dust_to_gas = Float64Col(pos=2)

def DelAttrs_dust(table):
        del table.attrs.FIELD_0_FILL 
        del table.attrs.FIELD_1_FILL
        del table.attrs.FIELD_2_FILL
        del table.attrs.NROWS    
        
class Particle_polariz(IsDescription):
        B_cen   = Float64Col(shape=3,pos=0)
        alpha   = Float64Col(pos=1)

def DelAttrs_polariz(table):
        del table.attrs.FIELD_0_FILL 
        del table.attrs.FIELD_1_FILL
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
                
                h5file.setNodeAttr("/", "molec", 
                                   phys.molec if hasattr(phys, 'molec') else 0,
                                   name=None)
                h5file.setNodeAttr("/", "pops", 0, name=None)
                h5file.setNodeAttr("/", "T_cmb", 
                                   phys.T_cmb if hasattr(phys, 'T_cmb') else 0,
                                   name=None)
                h5file.setNodeAttr("/", "T_in", 
                                   phys.T_in if hasattr(phys, 'T_in') else 0,
                                   name=None)
                h5file.setNodeAttr("/", "dust", 
                                   phys.dust if hasattr(phys, 'dust') else 0, 
                                   name=None)
                h5file.setNodeAttr("/", "polariz", 
                                   phys.polariz if hasattr(phys, 'polariz') else 0,
                                   name=None)
                h5file.setNodeAttr("/", "z", 
                                   phys.z if hasattr(phys, 'z') else 0.,
                                   name=None)
                
                GridType = mesh.grid.GridType
                if GridType == 'SPH1D':
                        h5file.setNodeAttr("/", "geom", "sph1d", name=None)
                        self._export_sph1d(mesh,phys,h5file)
                elif GridType == 'SPH2D':
                        h5file.setNodeAttr("/", "geom", "sph3d", name=None)
                        self._export_sph2d(mesh,phys,h5file)
                elif GridType == 'CYL2D':
                        h5file.setNodeAttr("/", "geom", "cyl3d", name=None)
                        self._export_cyl2d(mesh,phys,h5file)
                else:
                        pass
                
                h5file.close()        
        
        def _export_sph1d(self,mesh,phys,h5file):
                nr = mesh.grid.nr
                nt = 1
                np = 1
                
                # Create ZONE table
                table = h5file.createTable("/", 'ZONE', ParticleZone, "Grid table")
                particle = table.row
                particle['LEVEL']       = -1
                particle['X_max']       = [ mesh.R_p[nr], mesh.theta_p[nt], mesh.phi_p[np] ]
                particle['X_min']       = [ mesh.R_p[0] , mesh.theta_p[0] , mesh.phi_p[0]  ]
                particle['X_cen']       = [ 0.5*( mesh.grid.Rin + mesh.grid.Rout ), mesh.theta_c[0], mesh.phi_c[0] ]
                particle['NCHILDREN']   = nr * nt * np
                particle['NAXES']       = [nr,nt,np]
                #Insert a new particle record
                particle.append()
                table.flush()
                DelAttrsZone(table)
                
                # Create GRID table
                table = h5file.createTable("/", 'GRID', ParticleGrid, "Grid table")
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
                DelAttrsGrid(table)
                
                if phys.molec:
                        # Create MOLEC table
                        table = h5file.createTable("/", 'MOLEC', Particle_molec, "molecular table")
                        particle = table.row
                        for i in range(nr):
                                particle['X_mol']       = phys.X_mol[i]
                                particle.append()
                        table.flush()
                        DelAttrs_molec(table)
                
                if hasattr(phys, 'T_d'):
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
                

        def _export_sph2d(self,mesh,phys,h5file):
                nr = mesh.grid.nr
                nt = mesh.grid.nt
                np = 1
                
                # Create ZONE table
                table = h5file.createTable("/", 'ZONE', ParticleZone, "Grid table")
                particle = table.row
                particle['LEVEL']       = -1
                particle['X_max']       = [ mesh.R_p[nr], mesh.theta_p[nt], mesh.phi_p[np] ]
                particle['X_min']       = [ mesh.R_p[0] , mesh.theta_p[0] , mesh.phi_p[0]  ]
                particle['X_cen']       = [ 0.5*( mesh.grid.Rin + mesh.grid.Rout ), 0.5 * (mesh.theta_p[0] + mesh.theta_p[nt]), mesh.phi_c[0] ]
                particle['NCHILDREN']   = nr * nt * np
                particle['NAXES']       = [nr,nt,np]
                #Insert a new particle record
                particle.append()
                table.flush()
                DelAttrsZone(table)
                
                # Create GRID table
                table = h5file.createTable("/", 'GRID', ParticleGrid, "Grid table")
                particle = table.row
                for i in range(nr):
                    for j in range(nt):
                        particle['LEVEL']       = 0
                        particle['POS']         = i * nt + j
                        particle['X_max']       = [ mesh.R_p[i+1], mesh.theta_p[j+1], mesh.phi_p[np] ]
                        particle['X_min']       = [ mesh.R_p[i]  , mesh.theta_p[j] , mesh.phi_p[0]  ]
                        particle['X_cen']       = [ mesh.R_c[i]  , mesh.theta_c[j] , mesh.phi_c[0]  ]
                        particle['n_H2']        = phys.n_H2[i,j]
                        particle['V_cen']       = phys.V_gas[i,j]
                        particle['T_k']         = phys.T_k[i,j]
                        particle['V_t']         = phys.Vt[i,j]
                        # Insert a new particle record
                        particle.append()
                table.flush()
                DelAttrsGrid(table)
                
                if phys.molec:
                        # Create MOLEC table
                        table = h5file.createTable("/", 'MOLEC', Particle_molec, "molecular table")
                        particle = table.row
                        for i in range(nr):
                            for j in range(nt):
                                particle['X_mol']       = phys.X_mol[i,j]
                                particle.append()
                        table.flush()
                        DelAttrs_molec(table)
                
                if hasattr(phys, 'T_d'):
                        # Create DUST table
                        table = h5file.createTable("/", 'DUST', Particle_dust, "dust table")
                        particle = table.row
                        for i in range(nr):
                            for j in range(nt):
                                particle['T_d']         = phys.T_d[i,j]
                                particle['kapp_d']      = phys.kapp_d[i * nt + j]
                                particle['dust_to_gas'] = phys.dust_to_gas[i,j]
                                particle.append()
                        table.flush()
                        DelAttrs_dust(table)
                        
                if hasattr(phys, 'B_field'):
                        # Create POLARIZ table
                        table = h5file.createTable("/", 'POLARIZ', Particle_polariz, "polarization table")
                        particle = table.row
                        for i in range(nr):
                            for j in range(nt):
                                particle['B_cen']       = phys.B_field[i,j]
                                particle['alpha']       = phys.alpha[i,j]
                                particle.append()
                        table.flush()
                        DelAttrs_molec(table)



        def _export_cyl2d(self,mesh,phys,h5file):
                nrc = mesh.grid.nrc
                nz = mesh.grid.nz
                np = 1
                
                # Create ZONE table
                table = h5file.createTable("/", 'ZONE', ParticleZone, "Grid table")
                particle = table.row
                particle['LEVEL']       = -1
                particle['X_max']       = [ mesh.Rc_p[nrc], mesh.phi_p[np], mesh.z_p[nz] ]
                particle['X_min']       = [ mesh.Rc_p[0], mesh.phi_p[0], mesh.z_p[0] ]
                particle['X_cen']       = [ 0.5*( mesh.grid.Rc_in + mesh.grid.Rc_out ), mesh.phi_c[0], 0.5 * (mesh.z_p[0] + mesh.z_p[nz]) ]
                particle['NCHILDREN']   = nrc * np * nz
                particle['NAXES']       = [nrc, np, nz]
                #Insert a new particle record
                particle.append()
                table.flush()
                DelAttrsZone(table)
                
                # Create GRID table
                table = h5file.createTable("/", 'GRID', ParticleGrid, "Grid table")
                particle = table.row
                for i in range(nrc):
                    for j in range(nz):
                        particle['LEVEL']       = 0
                        particle['POS']         = i * nz + j
                        particle['X_max']       = [ mesh.Rc_p[i+1], mesh.phi_p[np], mesh.z_p[j+1] ]
                        particle['X_min']       = [ mesh.Rc_p[i]  , mesh.phi_p[0] , mesh.z_p[j]   ]
                        particle['X_cen']       = [ mesh.Rc_c[i]  , mesh.phi_c[0] , mesh.z_c[j]   ]
                        particle['n_H2']        = phys.n_H2[i,j]
                        particle['V_cen']       = phys.V_gas[i,j]
                        particle['T_k']         = phys.T_k[i,j]
                        particle['V_t']         = phys.Vt[i,j]
                        # Insert a new particle record
                        particle.append()
                table.flush()
                DelAttrsGrid(table)
                
                if phys.molec:
                        # Create MOLEC table
                        table = h5file.createTable("/", 'MOLEC', Particle_molec, "molecular table")
                        particle = table.row
                        for i in range(nrc):
                            for j in range(nz):
                                particle['X_mol']       = phys.X_mol[i,j]
                                particle.append()
                        table.flush()
                        DelAttrs_molec(table)
                
                if hasattr(phys, 'T_d'):
                        # Create DUST table
                        table = h5file.createTable("/", 'DUST', Particle_dust, "dust table")
                        particle = table.row
                        for i in range(nrc):
                            for j in range(nz):
                                particle['T_d']         = phys.T_d[i,j]
                                particle['kapp_d']      = phys.kapp_d[i * nz + j]
                                particle['dust_to_gas'] = phys.dust_to_gas[i,j]
                                particle.append()
                        table.flush()
                        DelAttrs_dust(table)
                        
                if hasattr(phys, 'B_field'):
                        # Create POLARIZ table
                        table = h5file.createTable("/", 'POLARIZ', Particle_polariz, "polarization table")
                        particle = table.row
                        for i in range(nrc):
                            for j in range(nz):
                                particle['B_cen']       = phys.B_field[i,j]
                                particle['alpha']       = phys.alpha[i,j]
                                particle.append()
                        table.flush()
                        DelAttrs_molec(table)