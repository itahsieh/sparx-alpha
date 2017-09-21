#include "sparx.h"
#include "task.h"
#include "vtk-wrapper.h"
#include "unit.h"

/* Global parameter struct */
static struct glb {
	DatINode *task;
        DatINode *unit;

	double ucon, overlap_vel;
	double dist, rotate[3], I_norm;
	SpModel model;
	size_t line;
	double lamb, freq;
        MirImg_Axis v;
        VtkData vtkdata;
        VtkFile vtkfile;

        int overlap, lte, tau, slice;
        int tracer;
} glb;


#define RELVEL(i,j)\
	glb.model.parms.mol->OL[NRAD*(i)+(j)]->RelativeVel
#define OVERLAP(i,j)\
	glb.model.parms.mol->OL[NRAD*(i)+(j)]->overlap

#define NRAD\
	glb.model.parms.mol->nrad
#define FREQ(i)\
	glb.model.parms.mol->rad[(i)]->freq


/* Subroutine prototypes */
static int InitModel(void);
static void *InitModelThread(void *tid_p);

static void InitRay(double *dx, double *dy, GeRay *ray);

static void Vtk_nested_hyosun(void);


static int CalcContrib();

static void *VtkContributionSph1dTread(void *tid_p);
static void *VtkContributionSph3dTread(void *tid_p);
static void *VtkContributionCyl3dTread(void *tid_p);

static void ContributionSubSamp(size_t idx, Zone *SampZone);

static void ContributionTracer( double * contrib, double * contrib_dust, double * tau, double * tau_dust, Zone * SampZone, GeVec3_d * SampCartPos);

static int HitVoxel( const GeRay * ray, const GeVox * SampVp, double * tSamp, size_t * side_Samp);
static int HitVoxel_sph3d( const GeRay *, const GeVox *, double * t, size_t * side);
static int HitVoxel_cyl3d( const GeRay *, const GeVox *, double * t, size_t * side);

static void CalcOpticalDepth( Zone *, const GeRay *, const double, double *, double *);
static void ContributionOfCell( Zone *, const GeRay *, const GeVec3_d *, double *, double *, double *, double *);

int ReachZone( const Zone * zp, const Zone * SampZone);
int ReachCell(const GeVec3_d * GeomPos, const GeVox * voxel, size_t side);




/*----------------------------------------------------------------------------*/

int SpTask_Visual(void)
{
	Mem_BZERO(&glb);
        int sts = 0;

// 1. GET PARAMETERS FROM PYTHON
        PyObject *o;
/*    1-1 those are the general parameters */
        /* out (mandatory) */
        if(!sts){
                sts = SpPy_GetInput_PyObj("out", &o);
                glb.vtkfile.FileName = Sp_PYSTR(o);
                SpPy_XDECREF(o);
        }


/*    1-2 get the task-based parameters */
	/* obs */
	if(!sts && !(sts = SpPy_GetInput_PyObj("obs", &o))) {
                
                /* task */
                PyObject *o_task = PyObject_GetAttrString(o, "task");
                glb.task = Dat_IList_NameLookup(TASKS, Sp_PYSTR(o_task));
                SpPy_XDECREF(o_task);
                

                switch (glb.task->idx){
                    // for task-lineobs
                    case TASK_LINECTB:
                    case TASK_ZEEMANCTB: 
                    { // get line transition
                        PyObject *o_line = PyObject_GetAttrString(o, "line");
                        glb.line = Sp_PYSIZE(o_line);
                        SpPy_XDECREF(o_line);
                    }
                        if(!sts) sts = SpPy_GetInput_bool("lte", &glb.lte);
                    
                    
                    { // get overlap velocity
                        PyObject *o3 = PyObject_GetAttrString(o, "overlap_vel");
                        glb.overlap_vel = Sp_PYDBL(o3);
                        SpPy_XDECREF(o3);
                        if (glb.overlap_vel == 0.0)
                            glb.overlap = 0;
                        else
                            glb.overlap = 1;
                    }
                    
                        /* chan */
                        if(!sts && !(sts = SpPy_GetInput_PyObj("chan", &o))) {
                            glb.v.n = Sp_PYSIZE(Sp_PYLST(o, 0));
                            glb.v.crpix = MirWr_CRPIX(glb.v.n);
                            glb.v.delt = Sp_PYDBL(Sp_PYLST(o, 1));
                        }
                    
                        break;
                    // for task-contobs 
                    case TASK_CONTCTB:
                    {    
                        PyObject *o_wavelen;
                        o_wavelen = PyObject_GetAttrString(o, "wavelen");
                        glb.lamb = Sp_PYDBL(o_wavelen);
                        glb.freq = PHYS_CONST_MKS_LIGHTC / glb.lamb;
                        SpPy_XDECREF(o_wavelen);
                    }
                        break;
                    case TASK_MODEL2VTK:
                        break;
                    default: 
                        /* Shouldn't reach here */
                        Deb_ASSERT(0);
                }
                SpPy_XDECREF(o);
                
                switch (glb.task->idx){
                      // for all radiation-related contrib tasks, except model2vtk
                      case TASK_LINECTB:
                      case TASK_ZEEMANCTB:
                      case TASK_CONTCTB:
                          /* radiation options */
                          /* unit */
                          if(!sts && !(sts = SpPy_GetInput_PyObj("unit", &o))) {
                              glb.unit = Dat_IList_NameLookup(UNITS, Sp_PYSTR(o));
                              Deb_ASSERT(glb.unit != NULL);
                              SpPy_XDECREF(o);
                          }
                          /* observer options  */
                          /* dist */
                          if(!sts) sts = SpPy_GetInput_dbl("dist", &glb.dist);
                          glb.dist /= Sp_LENFAC;
                          
                          /* rotate */
                          if(!sts && !(sts = SpPy_GetInput_PyObj("rotate", &o))) {
                              glb.rotate[0] = Sp_PYDBL(Sp_PYLST(o, 0));
                              glb.rotate[1] = Sp_PYDBL(Sp_PYLST(o, 1));
                              glb.rotate[2] = Sp_PYDBL(Sp_PYLST(o, 2));
                              SpPy_XDECREF(o);
                          }
                          break;
                }
                if(!sts) sts = SpPy_GetInput_bool("slice", &glb.slice);
	}

	
        TASK_TYPE task = glb.task->idx;
/*    1-3 read the source model */
        /* source */
        if(!sts) {
                int popsold = 0;
                /* only read populations data when non-LTE LINE/ZEEMAN mapping tasks*/ 
                if( task == TASK_LINE || task == TASK_ZEEMAN ){
                        if (!glb.lte)
                                popsold = 1;
                }
                sts = SpPy_GetInput_model("source","source", &glb.model, &popsold, task);
                
                
        }
        
        
/* 2. Initialize model */
        /* initialize sparx model */
        if(!sts) {
            sts = InitModel();
            // the dimension and grid
            Vtk_InitializeGrid( glb.v.n, 
                                glb.model.grid, 
                                &glb.vtkdata, 
                                glb.model.grid->voxel.geom, 
                                glb.slice
                              );
        }
/* 3. Computing the analyzed properties */
        if(!sts)
            if( task == TASK_LINECTB || task == TASK_CONTCTB)
                sts = CalcContrib();

/* 4. I/O : OUTPUT */
	if(!sts){
            // VTK visualization
            double scale_factor = 
                (task == TASK_MODEL2VTK) ? 
                1.0 : glb.I_norm/glb.ucon;
            
            Vtk_Output( &glb.vtkfile, 
                        &glb.vtkdata, 
                        &glb.model, 
                        glb.line, 
                        glb.v.n, 
                        task, 
                        scale_factor,
                        glb.unit
                      );
                
        }

#if 0	
#define MEAN_INT(ix, iy)\
	mean_int[(size_t)(iy) + glb.y.n * (size_t)(ix) ]
#define MEAN_VEL(ix, iy)\
	mean_vel[(size_t)(iy) + glb.y.n * (size_t)(ix) ]
#define MEAN_DEV(ix, iy)\
	mean_dev[(size_t)(iy) + glb.y.n * (size_t)(ix) ]
	
        // this visualization is ... not useful...
        // calculate three moments and output to VTK file
	if(glb.task->idx == TASK_LINE){
                double * mean_int = Mem_CALLOC(glb.x.n*glb.y.n, mean_int);
                double * mean_vel = Mem_CALLOC(glb.x.n*glb.y.n, mean_vel);
                double * mean_dev = Mem_CALLOC(glb.x.n*glb.y.n, mean_dev);


                
                for(size_t ix = 0; ix < glb.x.n; ix++) 
                        for(size_t iy = 0; iy < glb.y.n; iy++) {
                                MEAN_INT(ix,iy)=0.0;
                                MEAN_VEL(ix,iy)=0.0;
                                MEAN_DEV(ix,iy)=0.0;
                                for(size_t iv = 0; iv < glb.v.n; iv++) {
                                        double tempINT = (MirImg_PIXEL(*glb.image, iv, ix, iy)-glb.I_cmb)*glb.I_norm/glb.ucon;
                                        MEAN_INT(ix,iy) += tempINT;
                                        MEAN_VEL(ix,iy) += tempINT * ((double)iv - glb.v.crpix) * glb.v.delt;
                                }
                                MEAN_INT(ix,iy) /= (double)glb.v.n;
                                if (MEAN_INT(ix,iy)<1e-20){
                                        MEAN_VEL(ix,iy)=0.;
                                }
                                else{
                                        MEAN_VEL(ix,iy) /= MEAN_INT(ix,iy) * (double)glb.v.n;
                                }
                                for(size_t iv = 0; iv < glb.v.n; iv++) {
                                        double tempINT = (MirImg_PIXEL(*glb.image, iv, ix, iy)-glb.I_cmb)*glb.I_norm/glb.ucon;
                                        double dv = ((double)iv - glb.v.crpix) * glb.v.delt - MEAN_VEL(ix,iy);
                                        MEAN_DEV(ix,iy) += tempINT*dv*dv;
                                }
                                if (MEAN_INT(ix,iy)<1e-20){
                                        MEAN_DEV(ix,iy)=0.;
                                }
                                else{
                                        MEAN_DEV(ix,iy) /= MEAN_INT(ix,iy) * (double)glb.v.n;
                                        MEAN_DEV(ix,iy) = sqrt(MEAN_DEV(ix,iy));
                                }	
                        }

                char filename[32];
                sprintf(filename,"0thmov_%5s.vtk",glb.imgf->name);
                FILE *fp=fopen(filename,"w");
                fprintf(fp,"# vtk DataFile Version 3.0\n");
                fprintf(fp,"0thmov\n");
                fprintf(fp,"ASCII\n");
                fprintf(fp,"DATASET STRUCTURED_POINTS\n");
                fprintf(fp,"DIMENSIONS %zu %zu %d\n",glb.x.n,glb.y.n,2);
                fprintf(fp,"ORIGIN %f %f %d\n",-glb.x.crpix,-glb.y.crpix,1);
                fprintf(fp,"SPACING %d %d %d\n",1,1,1);
                fprintf(fp,"POINT_DATA %zu\n",glb.x.n * glb.y.n * 2);
                fprintf(fp,"SCALARS Intensity float 1\n");
                fprintf(fp,"LOOKUP_TABLE default\n");	
                for (size_t iz = 0; iz < 2; iz++)
                        for(size_t iy = 0; iy < glb.y.n; iy++) 
                                for(size_t ix = 0; ix < glb.x.n; ix++) 
                                        fprintf(fp,"%11.4e\n",MEAN_INT(ix,iy));
                fclose(fp);
                sprintf(filename,"1stmov_%5s.vtk",glb.imgf->name);
                fp=fopen(filename,"w");
                fprintf(fp,"# vtk DataFile Version 3.0\n");
                fprintf(fp,"1stmov\n");
                fprintf(fp,"ASCII\n");
                fprintf(fp,"DATASET STRUCTURED_POINTS\n");
                fprintf(fp,"DIMENSIONS %zu %zu %d\n",glb.x.n,glb.y.n,2);
                fprintf(fp,"ORIGIN %f %f %d\n",-glb.x.crpix,-glb.y.crpix,-1);
                fprintf(fp,"SPACING %d %d %d\n",1,1,1);
                fprintf(fp,"POINT_DATA %zu\n",glb.x.n * glb.y.n * 2);
                fprintf(fp,"SCALARS Velocity float 1\n");
                fprintf(fp,"LOOKUP_TABLE default\n");		
                for (size_t iz = 0; iz < 2; iz++)
                        for(size_t iy = 0; iy < glb.y.n; iy++) 
                                for(size_t ix = 0; ix < glb.x.n; ix++) 
                                        fprintf(fp,"%11.4e\n",MEAN_VEL(ix,iy));
                fclose(fp);
                sprintf(filename,"2ndmov_%5s.vtk",glb.imgf->name);
                fp=fopen(filename,"w");
                fprintf(fp,"# vtk DataFile Version 3.0\n");
                fprintf(fp,"2ndmov\n");
                fprintf(fp,"ASCII\n");
                fprintf(fp,"DATASET STRUCTURED_POINTS\n");
                fprintf(fp,"DIMENSIONS %zu %zu %d\n",glb.x.n,glb.y.n,2);
                fprintf(fp,"ORIGIN %f %f %d\n",-glb.x.crpix,-glb.y.crpix,0);
                fprintf(fp,"SPACING %d %d %d\n",1,1,1);
                fprintf(fp,"POINT_DATA %zu\n",glb.x.n * glb.y.n * 2);
                fprintf(fp,"SCALARS Deviation float 1\n");
                fprintf(fp,"LOOKUP_TABLE default\n");	
                for (size_t iz = 0; iz < 2; iz++)
                        for(size_t iy = 0; iy < glb.y.n; iy++) 
                                for(size_t ix = 0; ix < glb.x.n; ix++) 
                                        fprintf(fp,"%11.4e\n",MEAN_DEV(ix,iy));
                fclose(fp);
        }
#undefine MEAN_INT(ix, iy)
#undefine MEAN_VEL(ix, iy)
#undefine MEAN_DEV(ix, iy)
#endif

	

/* 5. Cleanup */
        SpModel_Cleanup( glb.model );
        Vtk_Mem_FREE( glb.model.grid->voxel.geom, &glb.vtkdata);

	return sts;
}

/*----------------------------------------------------------------------------*/

static int InitModel(void)
{
	Zone *root = glb.model.grid, *zp;
        SpPhysParm *parms = &glb.model.parms;
	
	int sts = 0;
        int task_id = glb.task->idx; 
        
        /* initialize line profile if LINE or ZEEMAN task */
        if( task_id == TASK_LINECTB || task_id == TASK_ZEEMANCTB ){
            
                glb.freq = parms->mol->rad[glb.line]->freq;
                glb.lamb = PHYS_CONST_MKS_LIGHTC / glb.freq;
                Deb_ASSERT(glb.line < parms->mol->nrad);
        }
        
        /* set the reference of the intensity: glb.ucon */
        if(glb.task->idx != TASK_MODEL2VTK){
            switch (glb.unit->idx){
                case UNIT_KPC:
                    glb.ucon = Phys_RayleighJeans(glb.freq, 1.0);
                    break;
                case UNIT_JYPC:
                    glb.ucon = (PHYS_UNIT_MKS_JY );
                    break;
                default:
                    Deb_ASSERT(0);
            }
            /* Sanity check */
            Deb_ASSERT((glb.ucon > 0) && (!Num_ISNAN(glb.ucon)) && (glb.ucon < HUGE_VAL));
        }
        
	/* initialization : construct overlapping table */
        if(glb.overlap){
            parms->mol->OL = Mem_CALLOC(NRAD*NRAD,parms->mol->OL);
            for(size_t i=0; i<NRAD; i++){
                for(size_t j=0; j<NRAD; j++){
                    parms->mol->OL[NRAD*i+j] = Mem_CALLOC(1, parms->mol->OL[NRAD*i+j]);
                    RELVEL(i,j) = ( 1e0 - FREQ(j)/FREQ(i) )*CONSTANTS_MKS_LIGHT_C;
                    if( fabs(RELVEL(i,j)) < glb.overlap_vel ){
                        OVERLAP(i,j) = 1;
                    }
                    else{
                        OVERLAP(i,j) = 0;
                    }
                }
            }
        }

        /* Set normalization intensity to 20K -- normalization prevents rounding
	   errors from creeping in when flux values are very small */
        if(task_id != TASK_MODEL2VTK){
		glb.I_norm = Phys_PlanckFunc(glb.freq, 10.0);
		Deb_ASSERT(glb.I_norm > 0); /* Just in case */

	}

	for(zp = Zone_GetMinLeaf(root); zp; zp = Zone_AscendTree(zp)) {
                SpPhys *pp;
		/* Pointer to physical parameters */
		pp = zp->data;

		if((pp->n_H2 > 1e-200) && !zp->children) {
			/* This is a non-empty leaf zone */
			pp->non_empty_leaf = 1;

			if(pp->X_mol > 0) {
				/* This zone contains tracer molecules */
				pp->has_tracer = 1;
          
			}
			else{
				pp->has_tracer = 0;
			}
		}
		else{
			pp->non_empty_leaf = 0;
			pp->has_tracer = 0;
		}
	}


	sts = SpUtil_Threads2(Sp_NTHREAD, InitModelThread);
	//SpUtil_Threads(InitModelThread);

	return sts;
}

/*----------------------------------------------------------------------------*/

static void *InitModelThread(void *tid_p)
{
    size_t tid = *((size_t *)tid_p);
    size_t zone_id;
    Zone *root = glb.model.grid, *zp;
    SpPhys *pp;
    int task_id = glb.task->idx; 
    
    if(task_id != TASK_MODEL2VTK){
        for(zp = Zone_GetMinLeaf(root), zone_id = 0; zp; zp = Zone_AscendTree(zp), zone_id++) {
            if(zone_id % Sp_NTHREAD == tid) {
                /* Check for thread termination */
                Sp_CHECKTERMTHREAD();
                
                /* Init zone parameters */
                pp = zp->data;
                
                if(task_id == TASK_CONTCTB) {
                    SpPhys_InitContWindows(pp, &glb.freq, (size_t)1);
                }
                else {
                    size_t nrad = glb.model.parms.mol->nrad;
                    double freq[nrad];
                    /* Set initial pops to either optically thin or LTE */
                    if(glb.lte){
                        // initialize level population for LTE condition
                        pp->mol = glb.model.parms.mol;
                        pp->pops_preserve = Mem_CALLOC(pp->mol->nlev, pp->pops_preserve);

                        // LTE : Boltzmann distribution
                        for(size_t j = 0; j < pp->mol->nlev; j++) 
                            pp->pops_preserve[j] = SpPhys_BoltzPops(pp->mol, j, pp->T_k);
                        
                        /* Allocate continuum emission/absorption */
                        for(size_t i = 0; i < nrad; i++) {
                            freq[i] = pp->mol->rad[i]->freq;
                        }
                        SpPhys_InitContWindows(pp, freq, nrad);
                    }
                        
                    pp->width = SpPhys_CalcLineWidth(pp);
                }
                
                /* Add dust emission/absorption if T_d > 0 */
                if(pp->T_d > 0) {
                    SpPhys_AddContinuum_d(pp, task_id == TASK_CONTCTB, pp->dust_to_gas);
                }
                /* Add free-free emission/absorption if T_ff > 0 */
                if(pp->X_e > 0) {
                    SpPhys_AddContinuum_ff(pp, task_id == TASK_CONTCTB);
                }
                
            }
        }
    }
    
    pthread_exit(NULL);
}

/*---------------------------------------------------------------------------- */
/* Initialize the position of the ray and its shooting direction */
static void InitRay(double *dx, double *dy, GeRay *ray)
{       
        Zone *root = glb.model.grid;

        /* Reset ray */
        Mem_BZERO(ray);
        

        double ModelSize = Zone_ZoneSize(root);
        double Model2DistanceRatio = ModelSize / glb.dist;
        
        if ( Model2DistanceRatio > 1.0 ){
            // distance is too close
            Sp_PRINT("The observing distance is too close!\n");
            Deb_ASSERT(0);
        }
        else{ 
            if ( Model2DistanceRatio > 1e-4 ){
                // stereopsis projection
                /* Init ray position to <dist, 0, 0> */
                GeRay_E(*ray, 0) = glb.dist;
                GeRay_E(*ray, 1) = 0.0;
                GeRay_E(*ray, 2) = 0.0;

                /* Set "reversed" direction of ray according to pixel position
                * !!! backward tracing !!! :
                *   theta = PI/2 - dy
                *   phi = PI - dx
                */
                
                double theta = 0.5 * M_PI - (*dy);
                double phi = M_PI - (*dx);
                
                /* Convert to Cartesian coordinates */
                GeRay_D(*ray, 0) = sin(theta) * cos(phi);
                GeRay_D(*ray, 1) = sin(theta) * sin(phi);
                GeRay_D(*ray, 2) = cos(theta);
            }
            /* parallel projection mode 
                if the ratio of domain size and observing distance is less than 10^-4,
                then switch to parallel projection mode to avoid precision lost
                particularly in cylindrical and spherical coordinate
                the ray shooting solves the second order equation with the coefficient 
                c = Ex^2 - r^2
            */
            else{
                GeRay_E(*ray, 0) = ModelSize;
                GeRay_E(*ray, 1) = (*dx) * glb.dist;
                GeRay_E(*ray, 2) = (*dy) * glb.dist;

                GeRay_D(*ray, 0) = -1.0;
                GeRay_D(*ray, 1) = 0.0;
                GeRay_D(*ray, 2) = 0.0;
            }
        }
        

        /* Rotate ray:
         * Since what we REALLY want to rotate is the model and that the rays
         * are pointed towards the model, rays should be rotated in the negative
         * direction, and in the opposite order of what we would've done to
         * rotate the model. */
        *ray = GeRay_Rotate(ray, 0, -glb.rotate[0]);
        *ray = GeRay_Rotate(ray, 1, -glb.rotate[1]);
        *ray = GeRay_Rotate(ray, 2, -glb.rotate[2]);

        /* Coordinate-dependent offset */
        switch(root->voxel.geom) {
                case GEOM_SPH1D:
                        break;
                        
                case GEOM_SPH3D:
                        break;

                case GEOM_REC3D:
                        ray->e = GeVec3_Add(&ray->e, &root->voxel.cen);
                        break;
                        
                case GEOM_CYL3D:
                        break;

                default: 
                        /* Shouldn't reach here */
                        Deb_ASSERT(0);
        }
        return;
}



/*----------------------------------------------------------------------------*/
/*   
write out VTK file of the excitation temperature for the nested Cartesian Cell in multiblock format 
*/
static void Vtk_nested_hyosun(void)
{
	Zone * root = glb.model.grid;
	static double k = PHYS_CONST_MKS_BOLTZK;
	
        FILE * fp=fopen("vis.pvd","w");
	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n");
	fprintf(fp,"  <Collection>\n");
        
	mkdir("multiblock", S_IRWXU | S_IRWXG | S_IRWXO);
	
	Zone * zp1 = root;
	for(size_t ix1=0; ix1 < zp1->naxes.x[0]; ix1++){
	 for(size_t iy1=0; iy1 < zp1->naxes.x[1]; iy1++){
	  for(size_t iz1=0; iz1 < zp1->naxes.x[2]; iz1++){
		size_t izone1 = iz1 + iy1*(zp1->naxes.x[2]) + ix1*(zp1->naxes.x[2])*(zp1->naxes.x[1]);
		Zone *zp2 = zp1->children[izone1];
		for(size_t ix2=0; ix2 < zp2->naxes.x[0]; ix2++){
		 for(size_t iy2=0; iy2 < zp2->naxes.x[1]; iy2++){
		  for(size_t iz2=0; iz2 < zp2->naxes.x[2]; iz2++){
			size_t izone2 = iz2 + iy2*(zp2->naxes.x[2]) + ix2*(zp2->naxes.x[2])*(zp2->naxes.x[1]);
			Zone *zp3 = zp2->children[izone2];
			char filename[22];
			sprintf(filename,"multiblock/block%2.2zu%2.2zu.vtr",izone1,izone2);
			fprintf(fp,"    <DataSet group=\"%zu\" dataset=\"0\" file=\"%22s\"/>\n",izone1*8+izone2,filename);
			
			FILE *fp2=fopen(filename,"w");
			fprintf(fp2,"<?xml version=\"1.0\"?>\n");
			fprintf(fp2,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
			fprintf(fp2,"<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">\n",0,64,0,64,0,32);
			fprintf(fp2,"    <Piece Extent=\"%d %d %d %d %d %d\">\n",0,64,0,64,0,32);
			fprintf(fp2,"      <Coordinates>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"X_COORDINATES\" NumberOfComponents=\"1\">\n");
			fprintf(fp2," %12.6e",zp3->voxel.min.x[0]);
			for(size_t ix=0; ix < zp3->naxes.x[0]; ix++){
				size_t izone = ix*(zp3->naxes.x[2])*(zp3->naxes.x[1]);
				Zone *zp = zp3->children[izone];
				fprintf(fp2," %12.6e",zp->voxel.max.x[0]);
			}
			fprintf(fp2,"\n        </DataArray>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"Y_COORDINATES\" NumberOfComponents=\"1\">\n");
			fprintf(fp2," %12.6e",zp3->voxel.min.x[1]);
			for(size_t iy=0; iy < zp3->naxes.x[1]; iy++){
				size_t izone = iy*(zp3->naxes.x[2]);
				Zone *zp = zp3->children[izone];
				fprintf(fp2," %12.6e",zp->voxel.max.x[1]);
			}
			fprintf(fp2,"\n        </DataArray>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"Z_COORDINATES\" NumberOfComponents=\"1\">\n");
			fprintf(fp2," %12.6e",zp3->voxel.min.x[2]);
			for(size_t iz=0; iz < zp3->naxes.x[2]; iz++){
				size_t izone = iz;
				Zone *zp = zp3->children[izone];
				fprintf(fp2," %12.6e",zp->voxel.max.x[2]);
			}
			fprintf(fp2,"\n        </DataArray>\n");
			fprintf(fp2,"      </Coordinates>\n");
			fprintf(fp2,"      <CellData>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"Tex*\" NumberOfComponents=\"1\">\n");
			
			for(size_t iz=0; iz < zp3->naxes.x[2]; iz++){
			 for(size_t iy=0; iy < zp3->naxes.x[1]; iy++){
			  for(size_t ix=0; ix < zp3->naxes.x[0]; ix++){
				size_t izone = iz + iy*(zp3->naxes.x[2]) + ix*(zp3->naxes.x[2])*(zp3->naxes.x[1]);
				Zone *zp = zp3->children[izone];
				SpPhys *pp = zp->data;
				if(pp->X_mol == 0.0){
					fprintf(fp2,"0.0\n");
				}
				else{
					MolTrRad *trans = pp->mol->rad[glb.line];
					size_t up = trans->up;
					size_t lo = trans->lo;
					double n_u = pp->pops_preserve[up];
					double n_l = pp->pops_preserve[lo];
					double E_u = pp->mol->lev[up]->E;
					double E_l = pp->mol->lev[lo]->E;
					double g_u = pp->mol->lev[up]->g;
					double g_l = pp->mol->lev[lo]->g;
					double Tex = (E_l-E_u)/(k*log((n_u*g_l)/(n_l*g_u))); 
					fprintf(fp2,"%12.6e\n",Tex/(pp->T_k));
				}
			  }
			 }
			}
			fprintf(fp2,"        </DataArray>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"Tk\" NumberOfComponents=\"1\">\n");
			
			for(size_t iz=0; iz < zp3->naxes.x[2]; iz++){
			 for(size_t iy=0; iy < zp3->naxes.x[1]; iy++){
			  for(size_t ix=0; ix < zp3->naxes.x[0]; ix++){
				size_t izone = iz + iy*(zp3->naxes.x[2]) + ix*(zp3->naxes.x[2])*(zp3->naxes.x[1]);
				Zone *zp = zp3->children[izone];
				SpPhys *pp = zp->data;
				fprintf(fp2,"%12.6e\n",pp->T_k);
			  }
			 }
			}
			fprintf(fp2,"        </DataArray>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"N_h2\" NumberOfComponents=\"1\">\n");
			
			for(size_t iz=0; iz < zp3->naxes.x[2]; iz++){
			 for(size_t iy=0; iy < zp3->naxes.x[1]; iy++){
			  for(size_t ix=0; ix < zp3->naxes.x[0]; ix++){
				size_t izone = iz + iy*(zp3->naxes.x[2]) + ix*(zp3->naxes.x[2])*(zp3->naxes.x[1]);
				Zone *zp = zp3->children[izone];
				SpPhys *pp = zp->data;
				fprintf(fp2,"%12.6e\n",pp->n_H2);
			  }
			 }
			}
			fprintf(fp2,"        </DataArray>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"contribution\" NumberOfComponents=\"1\">\n");
			
			for(size_t iz=0; iz < zp3->naxes.x[2]; iz++){
			 for(size_t iy=0; iy < zp3->naxes.x[1]; iy++){
			  for(size_t ix=0; ix < zp3->naxes.x[0]; ix++){
				size_t izone = iz + iy*(zp3->naxes.x[2]) + ix*(zp3->naxes.x[2])*(zp3->naxes.x[1]);
				Zone *zp = zp3->children[izone];
				SpPhys *pp = zp->data;
				if(pp->X_mol == 0.0){
					fprintf(fp2,"0.0\n");
				}
				else{
					MolTrRad *trans = pp->mol->rad[glb.line];
					size_t up = trans->up;
					size_t lo = trans->lo;
					double n_u = pp->pops_preserve[up];
					double n_l = pp->pops_preserve[lo];
					double E_u = pp->mol->lev[up]->E;
					double E_l = pp->mol->lev[lo]->E;
					double g_u = pp->mol->lev[up]->g;
					double g_l = pp->mol->lev[lo]->g;
					double Tex = (E_l-E_u)/(k*log((n_u*g_l)/(n_l*g_u))); 
					fprintf(fp2,"%12.6e\n", Tex * pp->n_H2 * pp->X_mol);
				}
			  }
			 }
			}
			fprintf(fp2,"        </DataArray>\n");
			fprintf(fp2,"      </CellData>\n");
			fprintf(fp2,"    </Piece>\n");
			fprintf(fp2,"  </RectilinearGrid>\n");
			fprintf(fp2,"</VTKFile>");
			fclose(fp2);
			
                  }
		 }
		}
	  }
	 }
	}
		
	fprintf(fp,"  </Collection>\n");
	fprintf(fp,"</VTKFile>");
	
	fclose(fp);
        
        return;
}



/*----------------------------------------------------------------------------*/

static int CalcContrib()
{
        GEOM_TYPE geom = glb.model.grid->voxel.geom;
        /* multithreading calculation of the contribution of the cells */
        switch (geom){
            case GEOM_SPH1D :
                return SpUtil_Threads2( Sp_NTHREAD, VtkContributionSph1dTread);
            case GEOM_SPH3D :
                return SpUtil_Threads2( Sp_NTHREAD, VtkContributionSph3dTread);
            case GEOM_CYL3D :
                return SpUtil_Threads2( Sp_NTHREAD, VtkContributionCyl3dTread);
            default:
                Deb_ASSERT(0);
        }
}

/*----------------------------------------------------------------------------*/

static void *VtkContributionSph1dTread(void *tid_p)
{
        size_t tid = *((size_t *)tid_p);
        Zone * root = glb.model.grid;
        VtkData *visual = &glb.vtkdata;
        
        // Dimension of the visualized resolution
        size_t nr = visual->sph3d->nr;
        size_t nt = visual->sph3d->nt;
        size_t np = visual->sph3d->np;
        
        // link to the global pointer
        double * theta  = visual->sph3d->theta;
        double * phi    = visual->sph3d->phi;

        size_t cell_id = 0;
        // calculate the contribution of the cells
        for( size_t i = 0; i < nr; i++){
            for( size_t j = 0; j < nt; j++){
                for( size_t k = 0; k < np; k++){
                    if ( cell_id % Sp_NTHREAD == tid ){
                        /* Check for thread termination */
                        Sp_CHECKTERMTHREAD();
                        // copy grid i to sampling zone
                        Zone SampZone = *root->children[i];
                        // the voxel pointer
                        GeVox *vp = &SampZone.voxel;
                        // change the GEOM of the voxel to SPH3D
                        vp->geom = GEOM_SPH3D;
                        
                        vp->min.x[1] = theta[j];
                        vp->max.x[1] = theta[j+1];
                        SampZone.index.x[1] = j;
                        vp->min.x[2] = phi[k];
                        vp->max.x[2] = phi[k+1];
                        SampZone.index.x[2] = k;
                        
                        size_t idx = ( i * nt + j ) * np + k ;
                        ContributionSubSamp( idx, &SampZone);
                        
                    }
                    cell_id += 1;
                }
            }
        }
        //pthread_exit(NULL);
        return NULL;
}

/*----------------------------------------------------------------------------*/

static void *VtkContributionSph3dTread(void *tid_p)
{
        size_t tid = *((size_t *)tid_p);
        Zone * root = glb.model.grid;
        VtkData *visual = &glb.vtkdata;
        
        // Dimension of the visualized resolution
        size_t nr = visual->sph3d->nr;
        size_t nt = visual->sph3d->nt;
        size_t np = visual->sph3d->np;
        
        // link to the global pointer
        double * theta  = visual->sph3d->theta;
        double * phi    = visual->sph3d->phi;
        
        size_t cell_id = 0;
        // calculate the contribution of the cells
        for( size_t i = 0; i < nr; i++){
          for( size_t j = 0; j < nt; j++){
            for( size_t k = 0; k < np; k++){
              if ( cell_id % Sp_NTHREAD == tid ){
                        /* Check for thread termination */
                        Sp_CHECKTERMTHREAD();
                        // izone : to coresponding zone index
                        size_t izone = ZoneIndex( GEOM_SPH3D, i, j, k, root);
                        // copy grid i_zone to sampling zone
                        Zone SampZone = *root->children[izone];
                        // the voxel pointer
                        GeVox *vp = &SampZone.voxel;
                        // change the GEOM of the voxel to SPH3D
                        vp->geom = GEOM_SPH3D;

                        if (root->naxes.x[1] == 1){
                                vp->min.x[1] = theta[j];
                                vp->max.x[1] = theta[j+1];
                                SampZone.index.x[1] = j;
                        }
                        if(root->naxes.x[2] == 1){
                                vp->min.x[2] = phi[k];
                                vp->max.x[2] = phi[k+1];
                                SampZone.index.x[2] = k;
                        }
                        
                        size_t idx = ( i * nt + j ) * np + k ;
                        ContributionSubSamp( idx, &SampZone);
                        
              }
              cell_id += 1;
            }
          }
        }
        //pthread_exit(NULL);
        return NULL;
}

/*----------------------------------------------------------------------------*/

static void *VtkContributionCyl3dTread(void *tid_p)
{
        size_t tid = *((size_t *)tid_p);
        Zone * root = glb.model.grid;
        VtkData *visual = &glb.vtkdata;
        
        // Dimension of the visualized resolution
        size_t nr = visual->cyl3d->nr;
        size_t np = visual->cyl3d->np;
        size_t nz = visual->cyl3d->nz;
        
        // link to the global pointer
        double * phi    = visual->cyl3d->phi;
        double * Z      = visual->cyl3d->Z;
        
        size_t cell_id = 0;
        // calculate the contribution of the cells
        for( size_t i = 0; i < nr; i++){
          for( size_t j = 0; j < np; j++){
            for( size_t k = 0; k < nz; k++){
              if ( cell_id % Sp_NTHREAD == tid ){
                        /* Check for thread termination */
                        Sp_CHECKTERMTHREAD();
                        // izone : to coresponding zone index
                        size_t izone = ZoneIndex( GEOM_CYL3D, i, j, k, root);
                        // copy grid i_zone to sampling zone
                        Zone SampZone = *root->children[izone];
                        // the voxel pointer
                        GeVox *vp = &SampZone.voxel;
                        // change the GEOM of the voxel to CYL3D
                        vp->geom = GEOM_CYL3D;

                        if (root->naxes.x[1] == 1){
                                vp->min.x[1] = phi[j];
                                vp->max.x[1] = phi[j+1];
                                SampZone.index.x[1] = j;
                        }
                        if(root->naxes.x[2] == 1){
                                vp->min.x[2] = Z[k];
                                vp->max.x[2] = Z[k+1];
                                SampZone.index.x[2] = k;
                        }
                        
                        size_t idx = ( i * np + j ) * nz + k ;
                        ContributionSubSamp( idx, &SampZone);
                        

              }
              cell_id += 1;
            }
          }
        }
    
        //pthread_exit(NULL);
        return NULL;
}

/*----------------------------------------------------------------------------*/

void ContributionSubSamp(size_t idx, Zone *SampZone){
        // initialize numer of sampling
        static int nSamp_initialized;
        static int nSamp1D;
        static double Devide_nSampCube;
        static double Devide_2nSamp1D;
        if (!nSamp_initialized) {
                nSamp1D = 1;
                Devide_nSampCube = 1. / ((double) (nSamp1D * nSamp1D * nSamp1D));
                Devide_2nSamp1D = 1. / ((double) (2 * nSamp1D));
                nSamp_initialized = 1;
        }
        
        // the contribution of the cell
        // sampling n * n * n points in a cell
        // extract the geometry of the voxel
        

        GeVox * vp = &SampZone->voxel;
        size_t nvelo = glb.v.n;
        VtkData *visual = &glb.vtkdata;
        
        double * contrib_dust = visual->contrib_dust + idx;
        double * tau_dust = visual->tau_dust + idx;
        
        double * contrib = visual->contrib[idx];
        double * tau   = visual->tau[idx];
        double * tau_dev = visual->tau_dev[idx];
        
        Mem_BZERO2( contrib_dust, 1);
        Mem_BZERO2( tau_dust, 1);

        Mem_BZERO2( contrib, nvelo);
        Mem_BZERO2( tau, nvelo);
        Mem_BZERO2( tau_dev, nvelo);
        
        for(int i = 0; i < nSamp1D; i++){
            for(int j = 0; j < nSamp1D; j++){
                for(int k = 0; k < nSamp1D; k++){
                    // sub-sampling position for the geom-coordinate
                    GeVec3_d SampGeomPos = GeSubSampPos(i, j, k, Devide_2nSamp1D, vp);
                    
                    // Sampling position in Cartesian (x,y,z) format
                    GeVec3_d SampPosXYZ;
                    switch (vp->geom){
                        case GEOM_SPH3D:
                                SampPosXYZ = GeVec3_Sph2Cart( &SampGeomPos);
                                break;
                        case GEOM_CYL3D:
                                SampPosXYZ = GeVec3_Cyl2Cart( &SampGeomPos);
                                break;
                        default:
                                Deb_ASSERT(0);
                            
                    }
                    // Call the contribution tracer
                    double * contrib_sub = Mem_CALLOC( nvelo, contrib_sub);
                    double * tau_sub = Mem_CALLOC( nvelo, tau_sub);
                    double contrib_dust_sub; 
                    double tau_dust_sub; 
                    
                    ContributionTracer( contrib_sub, &contrib_dust_sub, tau_sub, &tau_dust_sub, SampZone, &SampPosXYZ);

                    // statistics : sum up
                    for (size_t l = 0; l < nvelo; l++){
                            contrib[l] += contrib_sub[l];
                            tau[l] += tau_sub[l];
                            tau_dev[l] += tau_sub[l] * tau_sub[l];
                    }
                    *contrib_dust += contrib_dust_sub;
                    *tau_dust += tau_dust_sub;
                    
                    // free sub memory
                    free(contrib_sub);
                    free(tau_sub);
        }}}
        
        for (size_t l = 0; l < nvelo; l++){
                // average contrib, tau
                contrib[l] *= Devide_nSampCube;
                tau[l] *= Devide_nSampCube;
                /* evaluate tau_dev
                   tau_dev = ( 1/n * sum ( tau_i - tau_mean )^2 ) ^ 1/2
                           = ( 1/n * ( sum tau_i^2 - sum tau_mean^2 ) ) ^ 1/2
                           = ( 1/n * sum tau_i^2 - tau_mean^2 ) ^ 1/2         */
                tau_dev[l] *= Devide_nSampCube;
                tau_dev[l] -= tau[l] * tau[l];
                tau_dev[l] = sqrt( tau_dev[l] );
        }
        *contrib_dust *= Devide_nSampCube;
        *tau_dust *= Devide_nSampCube;

        return;
}

/*----------------------------------------------------------------------------*/


static void ContributionTracer( double *contrib, double *contrib_dust, double *tau_nu, double * tau_dust, Zone *SampZone, GeVec3_d *SampCartPos)
{
        GeVox * SampVp = &SampZone->voxel;
        GeRay ray;
        size_t side;
        Zone * root = glb.model.grid;
        double t;


#if 1
        GeVec3_d SampCartPosRotate = *SampCartPos;
        
        SampCartPosRotate = GeVec3_Rotate_x( &SampCartPosRotate, glb.rotate[0]);
        SampCartPosRotate = GeVec3_Rotate_y( &SampCartPosRotate, glb.rotate[1]);
        SampCartPosRotate = GeVec3_Rotate_z( &SampCartPosRotate, glb.rotate[2]);
        
        double dx = atan( SampCartPosRotate.x[1] / ( glb.dist - SampCartPosRotate.x[0] ) );
        double dy = atan( SampCartPosRotate.x[2] / ( glb.dist - SampCartPosRotate.x[0] ) );
        InitRay( &dx, &dy, &ray);
#else
        Mem_BZERO(&ray);
        GeRay_E(ray, 0) = glb.dist;
        GeRay_E(ray, 1) = 0.0;
        GeRay_E(ray, 2) = 0.0;
        
        ray = GeRay_Rotate(&ray, 0, -glb.rotate[0]);
        ray = GeRay_Rotate(&ray, 1, -glb.rotate[1]);
        ray = GeRay_Rotate(&ray, 2, -glb.rotate[2]);
        
        GeRay_D(ray, 0) = SampCartPos->x[0] - GeRay_E(ray, 0);
        GeRay_D(ray, 1) = SampCartPos->x[1] - GeRay_E(ray, 1);
        GeRay_D(ray, 2) = SampCartPos->x[2] - GeRay_E(ray, 2);
        ray.d = GeVec3_Normalize(&ray.d);

#endif
        /* Shoot ray at model and see what happens! */
        if(GeRay_IntersectVoxel(&ray, &root->voxel, &t, &side)) {
            /* Reset tau for all channels */
            Mem_BZERO2(tau_dust, 1);
            Mem_BZERO2(tau_nu, glb.v.n);
                
            
            /* Calculate intersection */
            ray = GeRay_Inc(&ray, t);

            /* Locate starting leaf zone according to intersection */
            Zone *zp = Zone_GetLeaf(root, side, &ray.e, &ray);
            
            #if 0
            if (
                (SampZone->index.x[0] == 48) && (SampZone->index.x[1] == 13) && (SampZone->index.x[2] == 17) ||
                (SampZone->index.x[0] == 47) && (SampZone->index.x[1] == 11) && (SampZone->index.x[2] == 73)
            )
            {
//                 printf("geom: %s\n",GEOM_TYPES[root->voxel.geom]);
                //printf("SampPos x y z : %E %E %E\n",SampCartPosRotate.x[0], SampCartPosRotate.x[1], SampCartPosRotate.x[2]);
                //printf("dist : %E\n", glb.dist);
                //printf("dx \t= %E, dy = %E\n",dx,dy);
                
                
                printf("ray.e \t= %E %E %E\n", ray.e.x[0], ray.e.x[1], ray.e.x[2]);
                printf("ray.d \t= %E %E %E\n", ray.d.x[0], ray.d.x[1], ray.d.x[2]);
                
//                 GeVec3_d tmp;
//                 tmp.x[0] = SampCartPosRotate.x[0] - glb.dist;
//                 tmp.x[1] = SampCartPosRotate.x[1];
//                 tmp.x[2] = SampCartPosRotate.x[2];
//                 tmp = GeVec3_Normalize(&tmp);
//                 printf("true direction \t= %E %E %E\n", tmp.x[0], tmp.x[1], tmp.x[2]);
                
                printf("SP index : %zu %zu %zu\n",SampZone->index.x[0],SampZone->index.x[1],SampZone->index.x[2]);
                
                printf("SampPos \t= %E %E %E\n", 
                SampCartPos->x[0], SampCartPos->x[1], SampCartPos->x[2]);
                printf("VoxelMin \t= %E %E %E\n", 
                SampVp->min.x[0], SampVp->min.x[1], SampVp->min.x[2]);
                printf("VoxelMax \t= %E %E %E\n", 
                SampVp->max.x[0], SampVp->max.x[1], SampVp->max.x[2]);
            }
            #endif
            
            int reached_sampling_zone = 0;
            
            /* Keep going until there's no next zone to traverse to */
            while(zp) {
                
                GeVec3_d RayGeomPos;
                
                #if 0
                if ((SampZone->index.x[0] == 47) && (SampZone->index.x[1] == 12) && (SampZone->index.x[2] == 73)){
                    printf("zp index : %zu %zu %zu\n",zp->index.x[0],zp->index.x[1],zp->index.x[2]);
                
                    // debugging
                    //printf("dx \t= %E, dy = %E\n",dx,dy);
                    printf("ray.e \t= %E %E %E\n", ray.e.x[0], ray.e.x[1], ray.e.x[2]);
                    //printf("ray.d \t= %E %E %E\n", ray.d.x[0], ray.d.x[1], ray.d.x[2]);
                    printf("SampPos \t= %E %E %E\n", 
                    SampCartPos->x[0], SampCartPos->x[1], SampCartPos->x[2]);
                    
                    printf("VoxelMin \t= %E %E %E\n", 
                    zp->voxel.min.x[0], zp->voxel.min.x[1], zp->voxel.min.x[2]);
                    
                    RayGeomPos = GeGeomPos(SampVp->geom, &ray.e);
                    printf("RayGeomPos \t= %E %E %E\n", 
                    RayGeomPos.x[0], RayGeomPos.x[1], RayGeomPos.x[2]);
                    printf("VoxelMax \t= %E %E %E\n", 
                    zp->voxel.max.x[0], zp->voxel.max.x[1], zp->voxel.max.x[2]);
                
                    if ( !ReachCell(&RayGeomPos, &zp->voxel, side) )
                        Deb_ASSERT(0);
                }
                #endif
                
                // see if the ray reach the target zone
                int reach_the_target_zone = ReachZone( zp, SampZone);
                
                
                // the ray reaches the taget zone
                if( reach_the_target_zone ){
                    
                    // the ordinates of the ray position
                    RayGeomPos = GeGeomPos(SampVp->geom, &ray.e);
                    
 

                    // see if the photon reach the sampling cell
                    int reached_cell = ReachCell(&RayGeomPos, SampVp, side);
                    
                    // if reached the sampling cell, escape the tau tracer.
                    if( reached_cell ){
                            reached_sampling_zone = 1;
                            break;
                    }
                    // the ray is not inside the sampling voxel but inside the target zone
                    else{
                            size_t side_Samp;
                            double tSamp;
                            int hit = HitVoxel( &ray, SampVp, &tSamp, &side_Samp);

                            GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &side);

                            if ( tSamp < t){
                                    CalcOpticalDepth( zp, &ray, tSamp, tau_nu, tau_dust);
                                    /* Calculate next position */
                                    ray = GeRay_Inc(&ray, tSamp);
                                    side = side_Samp;
                            }
                            else{
                                    // will not reach in this interval
                                    CalcOpticalDepth( zp, &ray, t, tau_nu, tau_dust);
                                    /* Calculate next position */
                                    ray = GeRay_Inc(&ray, t);
                                    /* Get next zone to traverse to */
                                    zp = Zone_GetNext(zp, &side, &ray);
                            }
                            #if 0
                            printf("hit = %d\n", hit);
                            printf("tSamp = %E, t = %E\n", tSamp, t);
                            #endif
                    }
                }
                // not inside the target voxel, keep tracing
                else{
                    /* Calculate path to next boundary */
                    GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &side);
                    
                    CalcOpticalDepth( zp, &ray, t, tau_nu, tau_dust);
                    
                    /* Calculate next position */
                    ray = GeRay_Inc(&ray, t);
                    /* Get next zone to traverse to */
                    zp = Zone_GetNext(zp, &side, &ray);
                }
            }
            if(!reached_sampling_zone){
                printf("didnot reach the sampling zone\n");
                printf("SampZone : %zu %zu %zu\n",SampZone->index.x[0],SampZone->index.x[1],SampZone->index.x[2]); 
                Deb_ASSERT(0);
            }
            
            ContributionOfCell( zp, &ray, SampCartPos, contrib, contrib_dust, tau_nu, tau_dust);
        }
        else{
            Deb_ASSERT(0);
        }
#if 0
        printf("contribution\n");
        for(size_t l = 0; l < glb.v.n;l++)
                printf("%E ",contrib[l]);
        printf("\n");
        printf("OK\n");exit(0);
#endif
        

        return;
}

/*----------------------------------------------------------------------------*/

static int HitVoxel( const GeRay * ray, const GeVox * voxel, double * t, size_t * side)
{
        switch (voxel->geom){
            case GEOM_SPH3D :
                return HitVoxel_sph3d( ray, voxel, t, side);
            case GEOM_CYL3D :
                return HitVoxel_cyl3d( ray, voxel, t, side);
            default:
                Deb_ASSERT(0);
        }
}

/*----------------------------------------------------------------------------*/

static int HitVoxel_sph3d( const GeRay *ray, const GeVox *voxel, double *t, size_t *side)
{
        static const double half_pi = 0.5 * 3.1415926535897932384626433832795;
        
        Deb_ASSERT( voxel->geom == GEOM_SPH3D );
        
        double R_in = voxel->min.x[0];
        double R_out = voxel->max.x[0];
        double theta_in = voxel->min.x[1];
        double theta_out = voxel->max.x[1];
        double phi_in = voxel->min.x[2];
        double phi_out = voxel->max.x[2];
        
        GeRay ray2;
        GeVec3_d RaySphPos;
        
        // t_Tl : distance to lower theta
        double t_Tl1, t_Tl2;
        GeRay_IntersectTheta(ray, theta_in, &t_Tl1, &t_Tl2);
        double t_Tl;
        if ( theta_in < half_pi && fabs(theta_in-half_pi) > 1e-10 ){
                t_Tl = Num_MAX(t_Tl1,t_Tl2);
                if(t_Tl <= 0.0) t_Tl = HUGE_VAL;
        }
        else{
                if ( t_Tl1 <= 0.0 ) t_Tl1 = HUGE_VAL;
                if ( t_Tl2 <= 0.0 ) t_Tl2 = HUGE_VAL;
                t_Tl = Num_MIN(t_Tl1,t_Tl2);
        }
        ray2 = GeRay_Inc(ray, t_Tl);
        RaySphPos = GeVec3_Cart2Sph(&ray2.e);
        #if 0
        printf("x y z \t= %E %E %E\n", ray2.e.x[0], ray2.e.x[1], ray2.e.x[2]);
        printf("t_Tl = %E, t_Tl1 = %E, t_Tl2 = %E\n", t_Tl, t_Tl1, t_Tl2);
        printf("R Theta Phi \t= %E %E %E\n", RaySphPos.x[0], RaySphPos.x[1], RaySphPos.x[2]);
        printf("%E %E %E %E %E %E\n", R_in, R_out, theta_in, theta_out, phi_in, phi_out);
        #endif
        if( ( RaySphPos.x[0] > R_out   ) ||
            ( RaySphPos.x[0] < R_in    ) ||
            ( RaySphPos.x[2] < phi_in  ) ||
            ( RaySphPos.x[2] > phi_out )    )
                t_Tl = HUGE_VAL;
        
        // t_Tu : distance to upper theta
        double t_Tu1, t_Tu2;
        GeRay_IntersectTheta(ray, theta_out, &t_Tu1, &t_Tu2);
        double t_Tu;
        if ( theta_out < half_pi || fabs(theta_in-half_pi) < 1e-10 ){
                if ( t_Tu1 <= 0.0 ) t_Tu1 = HUGE_VAL;
                if ( t_Tu2 <= 0.0 ) t_Tu2 = HUGE_VAL;
                t_Tu = Num_MIN(t_Tu1,t_Tu2);
        }
        else{
                t_Tu = Num_MAX(t_Tu1,t_Tu2);
                if ( t_Tu <= 0.0 ) t_Tu = HUGE_VAL;
        }
        ray2 = GeRay_Inc(ray, t_Tu);
        RaySphPos = GeVec3_Cart2Sph(&ray2.e);
        #if 0
        printf("x y z \t= %E %E %E\n", ray2.e.x[0], ray2.e.x[1], ray2.e.x[2]);
        printf("R Theta Phi \t= %E %E %E\n", RaySphPos.x[0], RaySphPos.x[1], RaySphPos.x[2]);
        printf("%E %E %E %E %E %E\n", R_in, R_out, theta_in, theta_out, phi_in, phi_out);
        #endif
        if( ( RaySphPos.x[0] > R_out   ) ||
            ( RaySphPos.x[0] < R_in    ) ||
            ( RaySphPos.x[2] < phi_in  ) ||
            ( RaySphPos.x[2] > phi_out )    )
                t_Tu = HUGE_VAL;

        // t_Pin : distance to inner phi
        double t_Pin;
        if( -sin(phi_in) * ray->d.x[0] + cos(phi_in) * ray->d.x[1] < 0.0 ){
                t_Pin = HUGE_VAL;
        }
        else{
                GeRay_IntersectPhi(ray, phi_in, &t_Pin);
                ray2 = GeRay_Inc(ray, t_Pin);
                RaySphPos = GeVec3_Cart2Sph(&ray2.e);
                if( ( RaySphPos.x[0] > R_out  ) ||
                    ( RaySphPos.x[0] < R_in   ) ||
                    ( RaySphPos.x[1] < theta_in ) ||
                    ( RaySphPos.x[1] > theta_out )   )
                        t_Pin = HUGE_VAL;
        }
        
        // distance to outer phi
        double t_Pout;
        if( sin(phi_out) * ray->d.x[0] - cos(phi_out) * ray->d.x[1] < 0.0 ){
                t_Pout=HUGE_VAL;
        }
        else{
                GeRay_IntersectPhi(ray, phi_out, &t_Pout);
                ray2 = GeRay_Inc(ray, t_Pout);
                RaySphPos = GeVec3_Cart2Sph(&ray2.e);
                if( ( RaySphPos.x[0] > R_out  ) ||
                    ( RaySphPos.x[0] < R_in   ) ||
                    ( RaySphPos.x[1] < theta_in ) ||
                    ( RaySphPos.x[1] > theta_out )   )
                        t_Pout = HUGE_VAL;
        }
        
        
        /* Init tmin */
        *t = HUGE_VAL;
        /* Find minimal intersection  */
        if( t_Tl < *t ){
                *t = t_Tl;
                *side = 2;
        }
        if( t_Tu < *t ){
                *t = t_Tu;
                *side = 3;
        }
        if( t_Pin < *t ){
                *t = t_Pin;
                *side = 4;
        }
        if( t_Pout < *t ){
                *t = t_Pout;
                *side = 5;
        }
#if 0
printf("side = %zu, t = %E\n", *side, *t);
printf("t_Tl = %E, t_Tu = %E, t_Pin = %E, t_Pout = %E\n", t_Tl, t_Tu, t_Pin, t_Pout);
#endif
        
        return ( *t == HUGE_VAL ) ? 0 : 1;
}

/*----------------------------------------------------------------------------*/

static int HitVoxel_cyl3d( const GeRay *ray, const GeVox *voxel, double *t, size_t *side)
{        
        Deb_ASSERT( voxel->geom == GEOM_CYL3D );
        
        double Rc_in = voxel->min.x[0];
        double Rc_out = voxel->max.x[0];
        double phi_in = voxel->min.x[1];
        double phi_out = voxel->max.x[1];
        double Z_in = voxel->min.x[2];
        double Z_out = voxel->max.x[2];
        
        GeRay ray2;
        GeVec3_d RayCylPos;

        // t_Pin : distance to inner phi
        double t_Pin;
        if( -sin(phi_in) * ray->d.x[0] + cos(phi_in) * ray->d.x[1] < 0.0 ){
                t_Pin = HUGE_VAL;
        }
        else{
                GeRay_IntersectPhi(ray, phi_in, &t_Pin);
                ray2 = GeRay_Inc(ray, t_Pin);
                RayCylPos = GeVec3_Cart2Cyl(&ray2.e);
                if( ( RayCylPos.x[0] > Rc_out  ) ||
                    ( RayCylPos.x[0] < Rc_in   ) ||
                    ( RayCylPos.x[2] < Z_in ) ||
                    ( RayCylPos.x[2] > Z_out )   )
                        t_Pin = HUGE_VAL;
        }
        
        // distance to outer phi
        double t_Pout;
        if( sin(phi_out) * ray->d.x[0] - cos(phi_out) * ray->d.x[1] < 0.0 ){
                t_Pout=HUGE_VAL;
        }
        else{
                GeRay_IntersectPhi(ray, phi_out, &t_Pout);
                ray2 = GeRay_Inc(ray, t_Pout);
                RayCylPos = GeVec3_Cart2Cyl(&ray2.e);
                if( ( RayCylPos.x[0] > Rc_out  ) ||
                    ( RayCylPos.x[0] < Rc_in   ) ||
                    ( RayCylPos.x[2] < Z_in ) ||
                    ( RayCylPos.x[2] > Z_out )   )
                        t_Pout = HUGE_VAL;
        }
        
        /* Init tmin */
        *t = HUGE_VAL;
        /* Find minimal intersection  */
        if( t_Pin < *t ){
                *t = t_Pin;
                *side = 2;
        }
        if( t_Pout < *t ){
                *t = t_Pout;
                *side = 3;
        }

        return ( *t == HUGE_VAL ) ? 0 : 1;
}

/*----------------------------------------------------------------------------*/

static void CalcOpticalDepth( Zone *zp, const GeRay *ray, const double t, double *tau_nu, double *tau_dust)
{       /* Pointer to physical parameters associated with this zone */
        SpPhys * pp = zp->data;
        
        /* Do radiative transfer only if gas is present in this zone */
        if(pp->non_empty_leaf) {
                /* Do calculations on all channels at this pixel. Try to minimize the 
                * amount of operations in this loop, since everything here is repeated 
                * for ALL channels, and can significantly increase computation time.
                */
                
                double k_dust = pp->cont[glb.line].k;
                double dtau_dust = k_dust * t * Sp_LENFAC;
                *tau_dust += dtau_dust;
                
                for(size_t iv = 0; iv < glb.v.n; iv++) {
                        /* Calculate velocity associated with this channel */
                        double dv = ((double)iv - glb.v.crpix) * glb.v.delt;
                        /* Reset emission and absorption coeffs */
                        double j_nu = 0.0;
                        double k_nu = 0.0;

                        if(pp->has_tracer) {
                                /* Calculate velocity line profile factor for this channel:
                                * This version averages over the line profile in steps
                                * of the local line width -- very time consuming! */
                                double vfac = SpPhys_GetVfac(ray, t, dv, zp, 0);
                                /* Calculate molecular line emission and absorption coefficients */
                                SpPhys_GetMoljk(pp, glb.line, vfac, &j_nu, &k_nu);
                                
                        }
                        
                        /* Add continuum absorption */
                        k_nu += pp->cont[glb.line].k;

                        /* Calculate source function and optical depth if
                        * absorption is NOT zero */
                        double dtau_nu = k_nu * t * Sp_LENFAC;

                        /* Accumulate total optical depth for this channel (must be done
                        * AFTER calculation of intensity!) */
                        tau_nu[iv] += dtau_nu;
                }
        }
        return;
}

/*----------------------------------------------------------------------------*/

static void ContributionOfCell( Zone *zp, const GeRay *ray, const GeVec3_d *SampCartPos, double * contrib, double * contrib_dust, double *tau_nu, double *tau_dust)
{
        /* Pointer to physical parameters associated with this zone */
        SpPhys *pp = zp->data;
        if(pp->non_empty_leaf) {
                /* Do calculations on all channels at this pixel. Try to minimize the 
                * amount of operations in this loop, since everything here is repeated 
                * for ALL channels, and can significantly increase computation time.
                */
                // the distance of the path to sampling
                double t = GeVec3_Mag2( &ray->e, SampCartPos);
                
                double j_dust = pp->cont[glb.line].j;
                double k_dust = pp->cont[glb.line].k;
                double S_dust = (fabs(k_dust) > 0.0) ?
                        j_dust / ( k_dust * glb.I_norm ) : 
                        0.;
                double dtau_dust = k_dust * t * Sp_LENFAC;
                *contrib_dust = 
                        S_dust * (1.0 - exp(-dtau_dust)) 
                        * exp(-(*tau_dust)) / t;
                *tau_dust += dtau_dust;
                
                for(size_t iv = 0; iv < glb.v.n; iv++) {
                        /* Calculate velocity associated with this channel */
                        double dv = ((double)iv - glb.v.crpix) * glb.v.delt;
                        /* Reset emission and absorption coeffs */
                        double j_nu = 0.0;
                        double k_nu = 0.0;

                        if(pp->has_tracer) {
                                /* Calculate velocity line profile factor for this channel:
                                * This version averages over the line profile in steps
                                * of the local line width -- very time consuming! */
                                double vfac = SpPhys_GetVfac(ray, t, dv, zp, 0);
                                /* Calculate molecular line emission and absorption coefficients */
                                SpPhys_GetMoljk(pp, glb.line, vfac, &j_nu, &k_nu);
                        }
                        
                        /* Add continuum emission/absorption */
                        j_nu += pp->cont[glb.line].j;
                        k_nu += pp->cont[glb.line].k;

                        /* Calculate source function and optical depth if
                        * absorption is NOT zero */
                        double dtau_nu = k_nu * t * Sp_LENFAC;
                        double S_nu = (fabs(k_nu) > 0.0) ?
                                j_nu / ( k_nu * glb.I_norm ) : 0.;

                        /* Calculate the contribution per length */
                        contrib[iv] = 
                                S_nu * (1.0 - exp(-dtau_nu))
                                * exp(-tau_nu[iv]) / t;
                        contrib[iv] -= *contrib_dust;
                                
                        /* Accumulate total optical depth for this channel (must be done
                        * AFTER calculation of intensity!) */
                        tau_nu[iv] += dtau_nu;
                }
        }
        else{
                Mem_BZERO2( contrib, glb.v.n);
        }
        
        return;
}

/*----------------------------------------------------------------------------*/

int ReachCell(const GeVec3_d * GeomPos, const GeVox * voxel, size_t side)
{
        static double threshold = 1E-6;
        static double threshold_0 = 1E-12;
        int reached_cell = 1;
        
        for (size_t i = 0; i < 3; i++){
            size_t side_in = 2 * i;
            size_t side_out = side_in + 1;
            if ( side_in == side ){
                if ( GeomPos->x[i] == 0.0 )
                    reached_cell *= 
                    ( fabs( voxel->min.x[i] - GeomPos->x[i] )  < threshold_0) ? 
                    1 : 0;
                else
                    reached_cell *= 
                    ( fabs( voxel->min.x[i] - GeomPos->x[i] ) / GeomPos->x[i] < threshold) ? 
                    1 : 0;
            }
            else{
                reached_cell *= 
                (voxel->min.x[i] <= GeomPos->x[i]) ? 1 : 0;
            }
            if ( side_out == side ){
                if ( GeomPos->x[i] == 0.0 )
                    reached_cell *= 
                    ( fabs( voxel->max.x[i] - GeomPos->x[i] ) < threshold_0 ) ? 
                    1 : 0;
                else
                    reached_cell *= 
                    ( fabs( voxel->max.x[i] - GeomPos->x[i] ) / GeomPos->x[i] < threshold ) ? 
                    1 : 0;
            }
            else{
                reached_cell *= 
                (voxel->max.x[i] >= GeomPos->x[i]) ? 1 : 0;
            }
        }
        
        return reached_cell;
}

/*----------------------------------------------------------------------------*/

int ReachZone( const Zone * zp, const Zone * SampZone)
{
        int reach_zone = 1;
        switch(zp->voxel.geom){
            case GEOM_SPH1D:
                reach_zone = 
                        ( zp->pos == SampZone->pos ) ? 1 : 0;
                break;
            case GEOM_SPH3D:
                for (int i =0; i < 3; i++){
                    if ( (i == 1) && (zp->parent->naxes.x[1] == 1) )
                        ;
                    else if ( (i == 2) && (zp->parent->naxes.x[2] == 1) )
                        ;
                    else
                        reach_zone *= 
                                (zp->index.x[i] == SampZone->index.x[i]) ? 1 : 0;
                }
                break;
            case GEOM_CYL3D:
                for (int i =0; i < 3; i++){
                    if ( (i == 1) && (zp->parent->naxes.x[1] == 1) )
                        ;
                    else
                        reach_zone *= 
                                (zp->index.x[i] == SampZone->index.x[i]) ? 1 : 0;
                }
                break;
            case GEOM_REC3D:
                for (int i =0; i < 3; i++)
                        reach_zone *= 
                                (zp->index.x[i] == SampZone->index.x[i]) ? 1 : 0;
                break;
            default:
                Deb_ASSERT(0);
        }
        
        return reach_zone;
}
