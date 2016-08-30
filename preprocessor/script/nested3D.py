#!/usr/bin/env python
mi=16;mj=mi;mk=mi       # Dimension of the coarsest grid
Lx=3e17;Ly=Lx;Lz=Lx # Length of the Domain
Critical_value=-119  # Density Criterion
Max_level=4            # Max hierachical level

import numpy as np
import array
from math import *
from tables import *
from tables.nodes import filenode
writegrid=1

# Physical parameter
T_k=10.
X_mol=1e-9
V_t=200.
nref=[2,2,2]
T_cmb=2.73
gas_to_dust=0.0
molec=""
geom='rec3d'
root="/"

# unit conversion
m_cm=1e-2
Nn_gcm=1e6*6.022e23/2.
pc_cm=1/3.08568025e18

# dimension of input data
ni=406
nj=186
# Create Cartesian Grid
dx=Lx/float(mi);dy=Ly/float(mj);dz=Lz/float(mk)
x=[];y=[];z=[]
for i in range(mi+1):
        x.append(-0.5*Lx+float(i)*dx)
for j in range(mj+1):        
        y.append(-0.5*Ly+float(j)*dy)
for k in range(mk+1):        
        z.append(-0.5*Lz+float(k)*dz)
x2 = array.array('d')
y2 = array.array('d')
z2 = array.array('d')

datadir='/home/vandine/work/GridConversion/'
# Load radius
tmpfile=datadir+'z_x1ap'
f=open( tmpfile,'rb')
ra = array.array('d')
ra.read(f,ni)
tmpfile=datadir+'z_x1bp'
f=open( tmpfile,'rb')
rb = array.array('d')
rb.read(f,ni)
# Load theta
tmpfile=datadir+'z_x2ap'
f=open( tmpfile,'rb')
thetaa = array.array('d')
thetaa.read(f,nj)
tmpfile=datadir+'z_x2bp'
f=open( tmpfile,'rb')
thetab = array.array('d')
thetab.read(f,nj)  
# Load density
tmpfile=datadir+'o_d__00100'
f=open( tmpfile,'rb')
density = array.array('d')
density.read(f,ni*nj)
density=np.reshape(density,(nj,ni))
# Load velocity
tmpfile=datadir+'o_v1_00100'
f=open( tmpfile,'rb')
Vr = array.array('d')
Vr.read(f,ni*nj)
Vr=np.reshape(Vr,(nj,ni))
tmpfile=datadir+'o_v2_00100'
f=open( tmpfile,'rb')
Vt = array.array('d')
Vt.read(f,ni*nj)
Vt=np.reshape(Vt,(nj,ni))
tmpfile=datadir+'o_v3_00100'
f=open( tmpfile,'rb')
Vp = array.array('d')
Vp.read(f,ni*nj)
Vp=np.reshape(Vp,(nj,ni))
f.close()

# geometry
R_in=ra[3]
R_out=ra[ni-3]

# convert R-theta to X-Z
Xb=np.zeros((nj,ni),np.float64)
Zb=np.zeros((nj,ni),np.float64)
for j in range(3,nj-2):
        for i in range(3,ni-3):                
                Xb[j,i]=rb[i]*sin(thetab[j])
                Zb[j,i]=rb[i]*cos(thetab[j])
                
# compute density gradient
grad_den=np.zeros((nj,ni),np.float64)
for j in xrange(3,nj-2):
        for i in xrange(3,ni-3):
                den_r=(density[j,i+1]-density[j,i-1])/(rb[i+1]-rb[i-1])
                den_theta=(density[j+1,i]-density[j-1,i])/(rb[i]*(thetab[j+1]-thetab[j-1]))
                grad_den[j,i]=sqrt(den_r*den_r+den_theta*den_theta)
                #print grad_den[j,i]

# write original 2D VTK file 
if (writegrid):
        fmb=open('2D.vtk', mode = "w")
        print >>fmb,'# vtk DataFile Version 3.0'
        print >>fmb,'2DHydro'
        print >>fmb,'ASCII'
        print >>fmb,'DATASET STRUCTURED_GRID'
        print >>fmb,'DIMENSIONS %(0)5d %(1)5d %(2)5d'%{'0':nj-5,'1':ni-5,'2':1}
        print >>fmb,'POINTS %(0)8d float'%{'0':(ni-5)*(nj-5)}
        for i in range(3,ni-2):
                for j in range(3,nj-2):                
                        print >>fmb,'%(0)11.4e %(1)1d %(2)11.4e'%{'0':ra[i]*sin(thetaa[j])*pc_cm,'1':0,'2':ra[i]*cos(thetaa[j])*pc_cm}
        print >>fmb,'CELL_DATA %(0)8d'%{'0':(ni-6)*(nj-6)}
        print >>fmb,'SCALARS density float 1'
        print >>fmb,'LOOKUP_TABLE default'
        for i in range(3,ni-3):
                for j in range(3,nj-3):
                        print >>fmb,'%(0)11.4e'%{'0':density[j,i]}
        fmb.close()




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


def main(pfile,direc,xaxis,yaxis,zaxis,n1,n2,level,position,pf2,npart,nzone):

	# Create ZONE table
	table = pfile.createTable(direc, 'ZONE', Particle, "Grid table")
	particle = table.row
	particle['LEVEL']  = level-1
	particle['POS'] = position
	particle['geom'] = geom
	particle['X_max'] =[ xaxis[n1[0]]*pc_cm,yaxis[n1[1]]*pc_cm,zaxis[n1[2]]*pc_cm ]
	particle['X_min'] =[ xaxis[0]*pc_cm,yaxis[0]*pc_cm,zaxis[0]*pc_cm ]
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

	cen_x=np.zeros(n1[0],np.float64)
        cen_y=np.zeros(n1[1],np.float64)
        cen_z=np.zeros(n1[2],np.float64)
        rho=np.zeros(n1,np.float64)
        Vx=np.zeros(n1,np.float64)
        Vy=np.zeros(n1,np.float64)
        Vz=np.zeros(n1,np.float64)
        for i in range(n1[0]):
                cen_x[i]=0.5*(xaxis[i]+xaxis[i+1])
        for j in range(n1[1]):
                cen_y[j]=0.5*(yaxis[j]+yaxis[j+1])
        for k in range(n1[2]):
                cen_z[k]=0.5*(zaxis[k]+zaxis[k+1])
	# Create GRID table
	table = pfile.createTable(direc, 'GRID', Particle, "Grid table")
	particle = table.row	
	for i in range(n1[0]):
	        for j in range(n1[1]):
	                for k in range(n1[2]):
	                        # write a row of grid table
		                particle['LEVEL']  = level
        	                particle['POS'] = n1[1]*n1[2]*i+n1[2]*j+k
        	                particle['geom'] = geom
        	                particle['X_max'] =[ xaxis[i+1]*pc_cm,yaxis[j+1]*pc_cm,zaxis[k+1]*pc_cm]
        	                particle['X_min'] =[ xaxis[i]*pc_cm,yaxis[j]*pc_cm,zaxis[k]*pc_cm]
        	                particle['X_cen'] =[ cen_x[i]*pc_cm,cen_y[j]*pc_cm,cen_z[k]*pc_cm]
	                        
	                        # project the cuboid zone to the R-theta plane
	                        abs_minx=min(abs(xaxis[i]),abs(xaxis[i+1]))
	                        abs_maxx=max(abs(xaxis[i]),abs(xaxis[i+1]))
	                        abs_miny=min(abs(yaxis[j]),abs(yaxis[j+1]))
	                        abs_maxy=max(abs(yaxis[j]),abs(yaxis[j+1]))
	                        minX=sqrt( abs_minx*abs_minx+abs_miny*abs_miny)
	                        maxX=sqrt( abs_maxx*abs_maxx+abs_maxy*abs_maxy)
	                        minZ=zaxis[k]
	                        maxZ=zaxis[k+1]
	                        abs_minz=min(abs(minZ),abs(maxZ))
	                        abs_maxz=max(abs(minZ),abs(maxZ))
	                        minR=sqrt( minX*minX+abs_minz*abs_minz)
	                        maxR=sqrt( maxX*maxX+abs_maxz*abs_maxz)
	                        if (maxZ>0.0):
	                                minTheta=atan(minX/maxZ)
	                        elif (maxZ<0.0):
	                                minTheta=atan(maxX/maxZ)+pi
	                        else:
	                                minTheta=0.5*pi
	                        if (minZ>0.0):
	                                maxTheta=atan(maxX/minZ)
	                        elif (minZ<0.0):
	                                maxTheta=atan(minX/minZ)+pi
	                        else:
	                                maxTheta=0.5*pi
	                        
	                        # narrow down the searching domain
	                        for tempi in xrange(3,ni-3):
	                                if (rb[tempi]>minR):
	                                        break
	                        i1=tempi
	                        for tempi in xrange(i1,ni-3):
	                                if (rb[tempi]>maxR):
	                                        break
	                        i2=tempi
	                        for tempj in xrange(3,nj-2):
	                                if (thetab[tempj]>minTheta):
	                                        break
	                        j1=tempj
	                        for tempj in xrange(j1,nj-2):
	                                if (thetab[tempj]>maxTheta):
	                                        break
	                        j2=tempj
	                        # search for maximum density gradient
	                        max_grad=1E-100
	                        for tempj in range(j1,j2):
                                        for tempi in range(i1,i2):
                                                if ( (minX<=Xb[tempj,tempi]<maxX) and (minZ<=Zb[tempj,tempi]<maxZ) ):
	                                                if ( max_grad<grad_den[tempj,tempi] ):
	                                                        max_grad=grad_den[tempj,tempi]
        	                
        	                # divide higher level        
        	                if ( ((log(max_grad)/log(2)-level)>Critical_value or minR<R_out<maxR) and level<Max_level ):
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
				        dx2=(xaxis[i+1]-xaxis[i])/float(n2[0])
				        dy2=(yaxis[j+1]-yaxis[j])/float(n2[1])
				        dz2=(zaxis[k+1]-zaxis[k])/float(n2[2])
				        x2=np.zeros(n2[0]+1,np.float64)
                                        y2=np.zeros(n2[1]+1,np.float64)
                                        z2=np.zeros(n2[2]+1,np.float64)
				        for tempi in range(n2[0]+1):
				                x2[tempi]=xaxis[i]+float(tempi)*dx2
				        for tempj in range(n2[1]+1):
				                y2[tempj]=yaxis[j]+float(tempj)*dy2
				        for tempk in range(n2[2]+1):
				                z2[tempk]=zaxis[k]+float(tempk)*dz2
				        #if (current_level<level+1):
				        #        current_level=level+1
				        #        print 'current level=',current_level
				                # recursive create the nested grid
			                (npart,nzone)=main(pfile,path,x2,y2,z2,n2,nref,level+1,particle['POS'],pf2,npart,nzone)
		                # the leaf zone
		                else:
		                        tempR=sqrt(cen_x[i]*cen_x[i]+cen_y[j]*cen_y[j]+cen_z[k]*cen_z[k])
			                if (tempR<=R_out):  # inside the boundary/ non-empty
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
                                                if (tempR<R_in):
                                                        print 'Mesh inside the boundary!'
                                                elif (tempR>R_out):
                                                        print 'Mesh outside the boundary!'
                                                # interpolation
                                                for tempai in xrange(3,ni-2):
                                                        if (ra[tempai]>tempR):
                                                                break
                                                for tempaj in xrange(3,nj-2):
                                                        if (thetaa[tempaj]>tempTheta):
                                                                break
                                                for tempbi in xrange(3,ni-3):
                                                        if (rb[tempbi]>tempR):
                                                                break
                                                for tempbj in xrange(3,nj-2):
                                                        if (thetab[tempbj]>tempTheta):
                                                                break              
                                                alpha=(tempR-rb[tempbi-1])/(rb[tempbi]-rb[tempbi-1])
                                                beta=(tempTheta-thetab[tempbj-1])/(thetab[tempbj]-thetab[tempbj-1])
                                                rho[i,j,k]=(1.-alpha)*(1.-beta)*density[tempbj-1,tempbi-1]+alpha*(1.-beta)*density[tempbj-1,tempbi]+\
                                                           beta*(1.-alpha)*density[tempbj,tempbi-1]+alpha*beta*density[tempbj,tempbi]
                                                tempVp=(1.-alpha)*(1.-beta)*Vp[tempbj-1,tempbi-1]+alpha*(1.-beta)*Vp[tempbj-1,tempbi]+\
                                                       beta*(1.-alpha)*Vp[tempbj,tempbi-1]+alpha*beta*Vp[tempbj,tempbi]
                                                alpha=(tempR-ra[tempai-1])/(ra[tempai]-ra[tempai-1])
                                                tempVr=(1.-alpha)*(1.-beta)*Vr[tempbj-1,tempai-1]+alpha*(1.-beta)*Vr[tempbj-1,tempai]+\
                                                       beta*(1.-alpha)*Vr[tempbj,tempai-1]+alpha*beta*Vr[tempbj,tempai]
                                                alpha=(tempR-rb[tempbi-1])/(rb[tempbi]-rb[tempbi-1])
                                                beta=(tempTheta-thetab[tempaj-1])/(thetab[tempaj]-thetab[tempaj-1])           
                                                tempVt=(1.-alpha)*(1.-beta)*Vt[tempaj-1,tempbi-1]+alpha*(1.-beta)*Vt[tempaj-1,tempbi]+\
                                                       beta*(1.-alpha)*Vt[tempaj,tempbi-1]+alpha*beta*Vt[tempaj,tempbi]
                                                Vx[i,j,k]=sin(tempTheta)*cos(tempPhi)*tempVr+cos(tempTheta)*cos(tempPhi)*tempVt-sin(tempPhi)*tempVp
                                                Vy[i,j,k]=sin(tempTheta)*sin(tempPhi)*tempVr+cos(tempTheta)*sin(tempPhi)*tempVt+cos(tempPhi)*tempVp
                                                Vz[i,j,k]=cos(tempTheta)*tempVr-sin(tempTheta)*tempVt
                                                # write out the non-empty-leaf zone
			                        particle['n_H2'] =rho[i,j,k]*Nn_gcm
			                        particle['V_cen'] =[Vx[i,j,k]*m_cm,Vy[i,j,k]*m_cm,Vz[i,j,k]*m_cm]      
        		                        particle['T_k'] =T_k
        		                        particle['X_mol'] =X_mol
        		                        particle['V_t'] =V_t
        		                        nzone=nzone+1
                                
                                # Insert a new particle record
        	                particle.append()
        	                if (level==0):
        	                        print n1[1]*n1[2]*i+n1[2]*j+k,'/',n1[0]*n1[1]*n1[2]
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
                        print >>f,'%(0)12.6e'%{'0':xaxis[i]*pc_cm},
                print >>f,'\n        </DataArray>'
                print >>f,'        <DataArray type="Float32" Name="Y_COORDINATES" NumberOfComponents="1">'
                for j in range(n1[1]+1):
                        print >>f,'%(0)12.6e'%{'0':yaxis[j]*pc_cm},
                print >>f,'\n        </DataArray>'
                print >>f,'        <DataArray type="Float32" Name="Z_COORDINATES" NumberOfComponents="1">'
                for k in range(n1[2]+1):
                        print >>f,'%(0)12.6e'%{'0':zaxis[k]*pc_cm},
                print >>f,'\n        </DataArray>'
                print >>f,'      </Coordinates>'
                print >>f,'    </Piece>'
                print >>f,'  </RectilinearGrid>'
                print >>f,'</VTKFile>'
                f.close()
                print >>fmb,'    <DataSet group="%(0)d" dataset="0" file="%(1)s"/>'%{'0':npart,'1':fname}
        return npart+1,nzone
 
# Timer
import time
tStart = time.time()
# Open a file in "w"rite mode
if (writegrid):
        fmb=open('multiblock.pvd', mode = "w")
        print >>fmb,'<?xml version="1.0"?>'
        print >>fmb,'<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
        print >>fmb,'  <Collection>'
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
naxe=[mi,mj,mk]
zone=0
(part,zone)=main(h5file,root,x,y,z,naxe,nref,0,0,fmb,0,zone)
# Close (and flush) the file
h5file.close()

if (writegrid):
        print >>fmb,'  </Collection>'
        print >>fmb,'</VTKFile>'
        fmb.close()
tEnd = time.time()

# print out meta information
total_time=tEnd-tStart
hh=int(total_time/60/60)
mm=int(total_time/60%60)
ss=int(total_time%60)
print 'Elapsing time = %(0)d h %(1)d m %(2)d s'%{'0':hh,'1':mm,'2':ss}
print 'max radius=',Lx*0.5*sqrt(3)
print 'original radius=',R_out
print zone,'zones'



