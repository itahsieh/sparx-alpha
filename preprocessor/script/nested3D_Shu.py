#!/usr/bin/env python
R_in=0.0
R_out=0.1
mi=8;mj=mi;mk=mi      # Dimension of the coarsest grid
Lx=2*R_out;Ly=Lx;Lz=Lx # Length of the Domain
Critical_density=37    # Density Criterion
Max_level=3            # Max hierachical level

import numpy as np
import array
from math import *
from tables import *
from tables.nodes import filenode
writegrid=1

# Physical parameter
n_max=1e12
vin=40
r0=0.01
T_k=10.
X_mol=1e-9
V_t=200.
nref=[2,2,2]
T_cmb=2.73
gas_to_dust=0.0
v_max=2*vin
molec=""
geom='rec3d'
root="/"


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
        ds=FloatCol(pos=17)
        NCHILDREN=Int64Col(pos=18)
        NAXES=Int64Col(shape=3,pos=19)
        T_d=Float64Col(pos=20)
        kapp_d=StringCol(itemsize=64,pos=21)
        T_ff=Float64Col(pos=22)
        kapp_ff=StringCol(itemsize=64,pos=23)
        T_bb=Float64Col(pos=24)


def main(pfile,direc,xaxis,yaxis,zaxis,cen_xx,cen_yy,cen_zz,density,Vxx,Vyy,Vzz,n1,n2,level,position,pf2,npart,cl):
	# Create ZONE table
	table = pfile.createTable(direc, 'ZONE', Particle, "Grid table")
	particle = table.row
	particle['LEVEL']  = level-1
	particle['POS'] = position
	particle['geom'] = geom
	particle['X_max'] =[ xaxis[n1[0]],yaxis[n1[1]],zaxis[n1[2]] ]
	particle['X_min'] =[ xaxis[0],yaxis[0],zaxis[0] ]
	particle['X_cen'] =[ 0.5*(particle['X_max'][0]+particle['X_min'][0]),0.5*(particle['X_max'][1]+particle['X_min'][1]),0.5*(particle['X_max'][2]+particle['X_min'][2])]
	particle['NCHILDREN'] =n1[0]*n1[1]*n1[2]
	particle['NAXES'] =n1
	#Insert a new particle record
	particle.append()
	table.flush()
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
	del table.attrs.NROWS


	# Create GRID table
	table = pfile.createTable(direc, 'GRID', Particle, "Grid table")
	particle = table.row
	rho2=np.zeros(n2)
        Vx2=np.zeros(n2)
        Vy2=np.zeros(n2)
        Vz2=np.zeros(n2)	
	for i in range(n1[0]):
            for j in range(n1[1]):
                for k in range(n1[2]):
                    # generate subgrid
                    dx2=(xaxis[i+1]-xaxis[i])/float(n2[0])
                    dy2=(yaxis[j+1]-yaxis[j])/float(n2[1])
                    dz2=(zaxis[k+1]-zaxis[k])/float(n2[2])
                    x2=np.zeros(n2[0]+1)
                    y2=np.zeros(n2[1]+1)
                    z2=np.zeros(n2[2]+1)
                    for tempi in range(n2[0]+1):
                            x2[tempi]=xaxis[i]+float(tempi)*dx2
                    for tempj in range(n2[1]+1):
                            y2[tempj]=yaxis[j]+float(tempj)*dy2
                    for tempk in range(n2[2]+1):
                            z2[tempk]=zaxis[k]+float(tempk)*dz2
                    cen_x2=np.zeros(n2[0])
                    cen_y2=np.zeros(n2[1])
                    cen_z2=np.zeros(n2[2])
                    for ii in range(n2[0]):
                            cen_x2[ii]=0.5*(x2[ii]+x2[ii+1])
                    for jj in range(n2[1]):
                            cen_y2[jj]=0.5*(y2[jj]+y2[jj+1])
                    for kk in range(n2[2]):
                            cen_z2[kk]=0.5*(z2[kk]+z2[kk+1])
                            
                    max_rho=0.0
                    for ii in range(n2[0]):
                        for jj in range(n2[1]):
                            for kk in range(n2[2]):
                                    # interpolate
                                    tempR = sqrt( cen_x2[ii]*cen_x2[ii] + cen_y2[jj]*cen_y2[jj] + cen_z2[kk]*cen_z2[kk] )
                                    if ( cen_z2[kk] == 0. ):
                                            tempTheta=0.5*pi
                                    else:
                                            tempTheta=atan( sqrt(cen_x2[ii]*cen_x2[ii]+cen_y2[jj]*cen_y2[jj]) / cen_z2[kk] )
                                    if (tempTheta<0):
                                            tempTheta=tempTheta+pi
                                    if (cen_x2[ii]==0.):
                                            if(cen_y2[jj]>0):
                                                    tempPhi=0.5*pi
                                            else:
                                                    tempPhi=-0.5*pi
                                    else:
                                            tempPhi=atan(cen_y2[jj]/cen_x2[ii])
                                    if (cen_x2[ii]<0): 
                                            tempPhi=tempPhi+pi
                                    elif (tempPhi<0):
                                            tempPhi=tempPhi+2*pi
                                    rho2[ii,jj,kk]=n_max*(tempR/r0)**-1.5
                                    tempVr=-v_max * (tempR/r0)**-0.5
                                    Vx2[ii,jj,kk]=sin(tempTheta)*cos(tempPhi)*tempVr
                                    Vy2[ii,jj,kk]=sin(tempTheta)*sin(tempPhi)*tempVr
                                    Vz2[ii,jj,kk]=cos(tempTheta)*tempVr
                                    max_rho=max(max_rho,rho2[ii,jj,kk])

                    # write a row of grid table
                    particle['LEVEL']  = level
                    particle['POS'] = n1[1]*n1[2]*i+n1[2]*j+k
                    particle['geom'] = geom
                    particle['X_max'] =[ xaxis[i+1],yaxis[j+1],zaxis[k+1]]
                    particle['X_min'] =[ xaxis[i],yaxis[j],zaxis[k]]
                    particle['X_cen'] =[ cen_xx[i],cen_yy[j],cen_zz[k]]
                    if ( ((log(max_rho)/log(2)-level)>Critical_density or abs(tempR-R_out)<0.5*sqrt(3.)*(xaxis[i+1]-xaxis[i])) and level<Max_level ):
                            particle['NCHILDREN'] =n2[0]*n2[1]*n2[2]
                            particle['NAXES'] =n2
                            gdir='grid'+'%(0)d'%{'0':n1[1]*n1[2]*i+n1[2]*j+k}
                            group = pfile.createGroup(direc,gdir,gdir)
                            if (direc=="/"):
                                    path=direc+gdir
                            else:
                                    path=direc+"/"+gdir
                            h5file.delNodeAttr(path, "TITLE", name=None)
                            h5file.delNodeAttr(path, "CLASS", name=None)
                            h5file.delNodeAttr(path, "VERSION", name=None)
                            h5file.setNodeAttr(path, "molec", molec, name=None)
                            h5file.setNodeAttr(path, "T_cmb", T_cmb, name=None)
                            h5file.setNodeAttr(path, "gas_to_dust", gas_to_dust, name=None)
                            h5file.setNodeAttr(path, "velfield", "grid ", name=None)
                            if (level+1>cl):
                                    cl=level+1
                                    print cl
                            (npart,cl)=main(pfile,path,x2,y2,z2,cen_x2,cen_y2,cen_z2,rho2,Vx2,Vy2,Vz2,n2,nref,level+1,particle['POS'],pf2,npart,cl)  
                    else:
                            if (tempR<=R_out):
                                    # H2 number density (m^-3)
                                    particle['n_H2'] = density[i,j,k]
                                    # velocity at cell center (m/s)
                                    particle['V_cen'] = [Vxx[i,j,k],Vyy[i,j,k],Vzz[i,j,k]]
                                    # magnetic field at cell center 
                                    # unit : Guass, used in Zeeman effect 
                                    # but the magnitude doesn't matter when calculating dust polarization
                                    particle['B_cen'] = [0.,0.,0.]
                                    # kinetic temperature (Kelvin)
                                    particle['T_k'] =T_k
                                    # fractional abundance (fraction related to H2 numebr)
                                    particle['X_mol'] =X_mol
                                    # turbulent speed (m/s)
                                    particle['V_t'] =V_t
                    # Insert a new particle record
                    particle.append()
                    
	table.flush()
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
	del table.attrs.NROWS
	# Write Grid for visualization
	if (writegrid):
        	fname='multiblock/post_'+str(npart)+'.vtr'
	        f=open(fname,'w')                
                # write in VTR format
                print >>f,'<?xml version="1.0"?>'
                print >>f,'<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">'
                print >>f,'  <RectilinearGrid WholeExtent="%(0)d %(1)d %(2)d %(3)d %(4)d %(5)d">'\
                          %{'0':0,'1':n1[0],'2':0,'3':n1[1],'4':0,'5':n1[2]}
                print >>f,'    <Piece Extent="%(0)d %(1)d %(2)d %(3)d %(4)d %(5)d">'\
                          %{'0':0,'1':n1[0],'2':0,'3':n1[1],'4':0,'5':n1[2]}
                print >>f,'      <Coordinates> '
                print >>f,'        <DataArray type="Float32" Name="X_COORDINATES" NumberOfComponents="1">'
                for i in range(n1[0]+1):
                        print >>f,'%(0)12.6e'%{'0':xaxis[i]},
                print >>f,'\n        </DataArray>'
                print >>f,'        <DataArray type="Float32" Name="Y_COORDINATES" NumberOfComponents="1">'
                for j in range(n1[1]+1):
                        print >>f,'%(0)12.6e'%{'0':yaxis[j]},
                print >>f,'\n        </DataArray>'
                print >>f,'        <DataArray type="Float32" Name="Z_COORDINATES" NumberOfComponents="1">'
                for k in range(n1[2]+1):
                        print >>f,'%(0)12.6e'%{'0':zaxis[k]},
                print >>f,'\n        </DataArray>'
                print >>f,'      </Coordinates>'
                print >>f,'      <CellData>'
                print >>f,'        <DataArray type="Float32" Name="density" NumberOfComponents="1" format="ascii">'
                for k in range(n1[2]):
                        for j in range(n1[1]):
                                for i in range(n1[0]):
                                        print >>f,density[i,j,k],
                print >>f,'\n        </DataArray>'
                print >>f,'      </CellData>'
                print >>f,'    </Piece>'
                print >>f,'  </RectilinearGrid>'
                print >>f,'</VTKFile>'
                f.close()
                print >>pf2,'    <DataSet group="%(0)d" dataset="0" file="%(1)s"/>'%{'0':npart,'1':fname}
        return (npart+1,cl)
 



# Create Cartesian Grid
naxe=[mi,mj,mk]
dx=Lx/float(naxe[0]);dy=Ly/float(naxe[1]);dz=Lz/float(naxe[2])
x=[];y=[];z=[]
for i in range(naxe[0]+1):
        x.append(-0.5*Lx+float(i)*dx)
for j in range(naxe[1]+1):        
        y.append(-0.5*Ly+float(j)*dy)
for k in range(naxe[2]+1):        
        z.append(-0.5*Lz+float(k)*dz)
cen_x=np.zeros(naxe[0])
cen_y=np.zeros(naxe[1])
cen_z=np.zeros(naxe[2])
for ii in range(naxe[0]):
        cen_x[ii]=0.5*(x[ii]+x[ii+1])
for jj in range(naxe[1]):
        cen_y[jj]=0.5*(y[jj]+y[jj+1])
for kk in range(naxe[2]):
        cen_z[kk]=0.5*(z[kk]+z[kk+1])
rho=np.zeros(naxe)
Vx=np.zeros(naxe)
Vy=np.zeros(naxe)
Vz=np.zeros(naxe)
for i in range(naxe[0]):
        for j in range(naxe[1]):
                for k in range(naxe[2]):
                        # interpolate
                        tempR=sqrt(cen_x[i]*cen_x[i]+cen_y[j]*cen_y[j]+cen_z[k]*cen_z[k])
                        if (cen_z[k]==0.):
                                tempTheta=0.5*pi
                        else:
                                tempTheta=atan( sqrt(cen_x[i]*cen_x[i]+cen_y[j]*cen_y[j]) / cen_z[k] )
                        if (tempTheta<0):
                                tempTheta=tempTheta+pi
                        if (cen_x[i]==0.):
                                if(cen_y[j]>0):
                                        tempPhi=0.5*pi
                                else:
                                        tempPhi=-0.5*pi
                        else:
                                tempPhi=atan(cen_y[j]/cen_x[i])
                        if (x[i]<0): 
                                tempPhi=tempPhi+pi
                        elif (tempPhi<0):
                                tempPhi=tempPhi+2*pi
                        if (tempR<=R_out):
                                rho[i,j,k]=n_max*(tempR/r0)**-1.5
                                tempVr=-v_max * (tempR/r0)**-0.5
                                Vx[i,j,k]=sin(tempTheta)*cos(tempPhi)*tempVr
                                Vy[i,j,k]=sin(tempTheta)*sin(tempPhi)*tempVr
                                Vz[i,j,k]=cos(tempTheta)*tempVr 
current_level=0

# Open a file in "w"rite mode
filename = "model"
h5file = openFile(filename, mode = "w", title = "Test file")
h5file.delNodeAttr("/", "TITLE", name=None)
h5file.delNodeAttr("/", "CLASS", name=None)
h5file.delNodeAttr("/", "VERSION", name=None)
h5file.delNodeAttr("/", "PYTABLES_FORMAT_VERSION", name=None)
h5file.setNodeAttr("/", "molec", molec, name=None)
h5file.setNodeAttr("/", "T_cmb", T_cmb, name=None)
h5file.setNodeAttr("/", "gas_to_dust", gas_to_dust, name=None)
h5file.setNodeAttr("/", "velfield", "grid ", name=None) 
if (writegrid):
        fmb=open('multiblock.pvd', mode = "w")
        print >>fmb,'<?xml version="1.0"?>'
        print >>fmb,'<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
        print >>fmb,'  <Collection>'
(part,current_level) = main(h5file,root,x,y,z,cen_x,cen_y,cen_z,rho,Vx,Vy,Vz,naxe,nref,0,0,fmb,0,current_level)
# Close (and flush) the file
h5file.close()
if (writegrid):
        print >>fmb,'  </Collection>'
        print >>fmb,'</VTKFile>'
        fmb.close()



