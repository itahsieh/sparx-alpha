#include "sparx.h"

/* Global parameter struct */
static struct glb {
	DatINode *task;
        int     overlap,
                lte,
                excit;
	DatINode *unit;
	double ucon, overlap_vel;
	MirImg_Axis x, y, v;
        MirFile *imgf, *StxQf, *StxUf, *tau_imgf;
        MirImg *image, *StokesQ, *StokesU, *StokesV, *sigma2, *tau_img;
	double dist, rotate[3], I_norm, I_cmb;
	SpModel model;
	size_t line;
	double lamb, freq;
	size_t nsubres;
	struct {
		double blc_x, blc_y, trc_x, trc_y;
		size_t nsub;
	} *subres;
        struct {
                struct {
                        size_t nr, nt, np;
                        double *radius, *theta, *phi;
                } *sph3d;
                struct {
                        size_t nr, np, nz;
                        double *Rc, *phi, *Z;
                } *cyl3d;
                double **contrib, **tau, **tau_dev;
        } *visual;
} glb;

enum {
	UNIT_K,
	UNIT_JYPX,
	UNIT_MKS,
        UNIT_CGS,
};

static DatINode UNITS[] = {
	{"K", UNIT_K},
	{"JY/PIXEL", UNIT_JYPX},
        {"MKS", UNIT_MKS},
        {"CGS", UNIT_MKS},
	{0, 0}
};


enum {
        TASK_LINE,
        TASK_ZEEMAN,
        TASK_CONT,
        TASK_COLDENS
};
static DatINode TASKS[] = {
        {"line", TASK_LINE},
        {"zeeman", TASK_ZEEMAN},
        {"cont", TASK_CONT},
        {"coldens", TASK_COLDENS},
        {0, 0}
};

#define MEAN_INT(ix, iy)\
	mean_int[(size_t)(iy) + glb.y.n * (size_t)(ix) ]
#define MEAN_VEL(ix, iy)\
	mean_vel[(size_t)(iy) + glb.y.n * (size_t)(ix) ]
#define MEAN_DEV(ix, iy)\
	mean_dev[(size_t)(iy) + glb.y.n * (size_t)(ix) ]

#define RELVEL(i,j)\
	(glb.model.parms.mol->OL[NRAD*i+j]->RelativeVel)
#define OVERLAP(i,j)\
	(glb.model.parms.mol->OL[NRAD*i+j]->overlap)

#define NRAD\
	(glb.model.parms.mol->nrad)
#define FREQ(i)\
	(glb.model.parms.mol->rad[i]->freq)


/* Subroutine prototypes */
static int InitModel(void);
static void *InitModelThread(void *tid_p);
static int CalcImage(void);
static void *CalcImageThreadLine(void *tid_p);
static void *CalcImageThreadZeeman(void *tid_p);
static void *CalcImageThreadCont(void *tid_p);
static void *CalcImageThreadColdens(void *tid_p);
static void RadiativeXferLine(double dx, double dy, double *I_nu, double *tau_nu, size_t tid);
static void RadiativeXferOverlap(double dx, double dy, double *I_nu, double *tau_nu, size_t tid);
static void RadiativeXferZeeman(double dx, double dy, double *V_nu, double *tau_nu, size_t tid);
static void RadiativeXferCont(double dx, double dy, double *I_nu, double *Q_nu, double *U_nu, double *sigma2, double *tau_nu);
static void ColumnDensityTracer(double dx, double dy, double *CD);
static void InitRay(double *dx, double *dy, GeRay *ray);
static void InitLOSCoord( double *dx, 
                          double *dy, 
                          GeRay *ray, 
                          GeVec3_d *z, 
                          GeVec3_d *n, 
                          GeVec3_d *e
                        );
static size_t Init_nsub(size_t ix, size_t iy);
static void InitSubPixel( double *dx, 
                          double *dy, 
                          size_t ix, 
                          size_t iy, 
                          size_t isub, 
                          size_t jsub, 
                          size_t nsub
                        );
static void visualization(void);


enum{
        VISUAL_CONTRIBUTION,
        VISUAL_TAU,
        VISUAL_TAU_DEV,
        VISUAL_EXITATION
};

DatINode VISUAL_TYPES[] = {
        {"contribution", VISUAL_CONTRIBUTION},
        {"optical_depth", VISUAL_TAU},
        {"optical_depth_deviation", VISUAL_TAU_DEV},
        {"exitation_temperature", VISUAL_EXITATION},
        {0, 0}
};

static int generic_vtk(int);

static int vtk_rec3d(void);
static int vtk_cyl3d(void);
static int vtk_sph3d(void);
static int vtk_sph1d(void);

static int CalcContrib(GEOM_TYPE geom_type);

static void *VtkContributionSph1dTread(void *tid_p);
static void *VtkContributionSph3dTread(void *tid_p);
static void *VtkContributionCyl3dTread(void *tid_p);

static void ContributionSubSamp(double *contrib, double *tau, double *tau_dev, Zone *SampZone, size_t nvelo, size_t tid);

static void ContributionTracer( double *, double *, Zone *, GeVec3_d *, size_t tid);

static int HitSph3dVoxel( const GeRay *, const GeVox *, double * t, size_t * side);
static int HitCyl3dVoxel( const GeRay *, const GeVox *, double * t, size_t * side);

static void CalcOpticalDepth( Zone *, const GeRay *, const double, double *, size_t tid);
static void ContributionOfCell( Zone *, const GeRay *, const GeVec3_d *, double *, double *, size_t tid);


/*----------------------------------------------------------------------------*/

int SpTask_Telsim(void)
{
	Mem_BZERO(&glb);
        int sts = 0;

// 1. GET PARAMETERS FROM PYTHON
        PyObject *o;
/*    1-1 those are the general parameters */
        /* source */
        if(!sts) {
                sts = SpPy_GetInput_model("source","source", &glb.model);
        }

        /* npix */
        if(!sts && !(sts = SpPy_GetInput_PyObj("npix", &o))) {
                glb.x.n = Sp_PYSIZE(Sp_PYLST(o, 0));
                glb.x.crpix = MirWr_CRPIX(glb.x.n);
                glb.y.n = Sp_PYSIZE(Sp_PYLST(o, 1));
                glb.y.crpix = MirWr_CRPIX(glb.y.n);
                SpPy_XDECREF(o);
        }
        /* cell */
        if(!sts && !(sts = SpPy_GetInput_PyObj("cell", &o))) {
                glb.x.delt = Sp_PYDBL(Sp_PYLST(o, 0));
                glb.y.delt = Sp_PYDBL(Sp_PYLST(o, 1));
                SpPy_XDECREF(o);
        }
        /* chan */
        if(!sts && !(sts = SpPy_GetInput_PyObj("chan", &o))) {
                glb.v.n = Sp_PYSIZE(Sp_PYLST(o, 0));
                glb.v.crpix = MirWr_CRPIX(glb.v.n);
                glb.v.delt = Sp_PYDBL(Sp_PYLST(o, 1));
                SpPy_XDECREF(o);
        }

        /* out (mandatory) */
        if(!sts){
                sts = SpPy_GetInput_mirxy_new("out", glb.x.n, glb.y.n, glb.v.n, &glb.imgf);
        }
        
        /* dist */
        if(!sts) sts = SpPy_GetInput_dbl("dist", &glb.dist);
          
        /* rotate */
        if(!sts && !(sts = SpPy_GetInput_PyObj("rotate", &o))) {
                glb.rotate[0] = Sp_PYDBL(Sp_PYLST(o, 0));
                glb.rotate[1] = Sp_PYDBL(Sp_PYLST(o, 1));
                glb.rotate[2] = Sp_PYDBL(Sp_PYLST(o, 2));
                SpPy_XDECREF(o);
        }
        /* subres */
        if(!sts && !(sts = SpPy_GetInput_PyObj("subres", &o))) {
                if(o != Py_None) {
                        /* Get number of boxes */
                        glb.nsubres = (size_t)PyList_Size(o);
                        glb.subres = Mem_CALLOC(glb.nsubres, glb.subres);
                        for(size_t i = 0; i < glb.nsubres; i++) {
                                PyObject *o_sub;
                                o_sub = PyList_GetItem(o, (Py_ssize_t)i);
                                glb.subres[i].blc_x = Sp_PYDBL(PyList_GetItem(o_sub, (Py_ssize_t)0));
                                glb.subres[i].blc_y = Sp_PYDBL(PyList_GetItem(o_sub, (Py_ssize_t)1));
                                glb.subres[i].trc_x = Sp_PYDBL(PyList_GetItem(o_sub, (Py_ssize_t)2));
                                glb.subres[i].trc_y = Sp_PYDBL(PyList_GetItem(o_sub, (Py_ssize_t)3));
                                glb.subres[i].nsub = Sp_PYSIZE(PyList_GetItem(o_sub, (Py_ssize_t)4));
                        }
                }
                SpPy_XDECREF(o);
        }

/*    1-2 get the task-based parameters */
	/* obs */
	if(!sts && !(sts = SpPy_GetInput_PyObj("obs", &o))) {
                PyObject *o_task;
                /* task */
                o_task = PyObject_GetAttrString(o, "task");
                glb.task = Dat_IList_NameLookup(TASKS, Sp_PYSTR(o_task));
                SpPy_XDECREF(o_task);	
                
                PyObject *o_line;
                PyObject *o1, *o2, *o3;
                switch (glb.task->idx){
                  // for task-contobs 
                  case TASK_CONT:
                          /* tau (optional) */
                          if(!sts && SpPy_CheckOptionalInput("tau")) {
                                  sts = SpPy_GetInput_mirxy_new("tau", glb.x.n, glb.y.n, glb.v.n, &glb.tau_imgf);
                          }
                          PyObject *o_wavelen;
                          o_wavelen = PyObject_GetAttrString(o, "wavelen");
                          glb.lamb = Sp_PYDBL(o_wavelen);
                          glb.freq = PHYS_CONST_MKS_LIGHTC / glb.lamb;
                          SpPy_XDECREF(o_wavelen);
                          break;
                  
                  // for task-lineobs
                  case TASK_ZEEMAN:
                          // get line transition
                          o_line = PyObject_GetAttrString(o, "line");
                          glb.line = Sp_PYSIZE(o_line);
                          SpPy_XDECREF(o_line);
                          
                          if(!sts) sts = SpPy_GetInput_bool("lte", &glb.lte);
                          if(!sts){
                                  if(glb.lte) sts = SpPy_GetInput_molec("molec", &glb.model.parms.mol);
                          }
                          Deb_ASSERT(glb.model.parms.mol != NULL);
                          
                          glb.freq = glb.model.parms.mol->rad[glb.line]->freq;
                          glb.lamb = PHYS_CONST_MKS_LIGHTC / glb.freq;
                          Deb_ASSERT(glb.line < glb.model.parms.mol->nrad);
                          
                          break;
                  // for task-lineobs
                  case TASK_LINE:
                          // get line transition
                          o1 = PyObject_GetAttrString(o, "line");
                          glb.line = Sp_PYSIZE(o1);
                          SpPy_XDECREF(o1);
                          
                          if(!sts) sts = SpPy_GetInput_bool("lte", &glb.lte);
                          if(!sts){
                                  if(glb.lte) sts = SpPy_GetInput_molec("molec", &glb.model.parms.mol);
                          }
                          Deb_ASSERT(glb.model.parms.mol != NULL);
                          
                          glb.freq = glb.model.parms.mol->rad[glb.line]->freq;
                          glb.lamb = PHYS_CONST_MKS_LIGHTC / glb.freq;
                          Deb_ASSERT(glb.line < glb.model.parms.mol->nrad);
                          
                          /* tau (optional) */
                          if(!sts && SpPy_CheckOptionalInput("tau")) {
                                  sts = SpPy_GetInput_mirxy_new("tau", glb.x.n, glb.y.n, glb.v.n, &glb.tau_imgf);
                          }

                          // get overlap switch
                          o2 = PyObject_GetAttrString(o, "overlap_int");
                          glb.overlap = Sp_PYINT(o2);
                          
                          // get overlap velocity
                          o3 = PyObject_GetAttrString(o, "overlap_vel");
                          glb.overlap_vel = Sp_PYDBL(o3);
                          SpPy_XDECREF(o2);
                          SpPy_XDECREF(o3);
                          
                          /* excitation visualization */
                          if(!sts) 
                                  sts = SpPy_GetInput_bool("excit", &glb.excit);
                          break;
                  default: 
                          /* Shouldn't reach here */
                          Deb_ASSERT(0);
                }
	}
	
	/* unit */
        if(!sts && !(sts = SpPy_GetInput_PyObj("unit", &o))) {
                glb.unit = Dat_IList_NameLookup(UNITS, Sp_PYSTR(o));
                Deb_ASSERT(glb.unit != NULL);
                if(glb.task->idx != TASK_COLDENS){
                        switch (glb.unit->idx){
                                case UNIT_K:
                                        glb.ucon = Phys_RayleighJeans(glb.freq, 1.0);
                                        break;
                                case UNIT_JYPX:
                                        glb.ucon = (PHYS_UNIT_MKS_JY / (glb.x.delt * glb.y.delt));
                                        break;
                                default:
                                        Deb_ASSERT(0);
                        }
                        /* Sanity check */
                        Deb_ASSERT((glb.ucon > 0) && (!Num_ISNAN(glb.ucon)) && (glb.ucon < HUGE_VAL));
                        SpPy_XDECREF(o);
                }
        }
	

/* 2. Initialize model */
	if(!sts) sts = InitModel();

/* 3. Synthesize image */
	if(!sts) {
		/* Allocate image */
		glb.image = MirImg_Alloc(glb.x, glb.y, glb.v);
		glb.image->restfreq = glb.freq;
		if(glb.task->idx == TASK_CONT){
			glb.StokesQ = MirImg_Alloc(glb.x, glb.y, glb.v);
			glb.StokesU = MirImg_Alloc(glb.x, glb.y, glb.v);
			glb.sigma2 = MirImg_Alloc(glb.x, glb.y, glb.v);
			glb.StokesQ->restfreq = glb.freq;
			glb.StokesU->restfreq = glb.freq;
			glb.sigma2->restfreq = glb.freq;
		}
		if(glb.tau_imgf)
			glb.tau_img = MirImg_Alloc(glb.x, glb.y, glb.v);
		/* Calculate image */
		sts = CalcImage();
	}

/* 4. I/O : OUTPUT */
	if(!sts){
                double scale_factor;
                int stokes;
                switch(glb.task->idx){
                  // output column density image
                  case TASK_COLDENS:
                          switch(glb.unit->idx){
                                  case UNIT_CGS:
                                          scale_factor = 1e-7;
                                          break;
                                  case UNIT_MKS:
                                          scale_factor = 1.0;
                                          break;
                          }
                          stokes = 0;
                          
                          char filename[32];
                          sprintf(filename,"%s.vtk",glb.imgf->name);
                          {
                            FILE *fp=fopen(filename,"w");
                            fprintf(fp,"# vtk DataFile Version 3.0\n");
                            fprintf(fp,"Column density\n");
                            fprintf(fp,"ASCII\n");
                            fprintf(fp,"DATASET STRUCTURED_POINTS\n");
                            fprintf(fp,"DIMENSIONS %zu %zu %d\n",glb.x.n,glb.y.n,1);
                            fprintf(fp,"ORIGIN %f %f %d\n",-glb.x.crpix,-glb.y.crpix,1);
                            fprintf(fp,"SPACING %d %d %d\n",1,1,1);
                            fprintf(fp,"POINT_DATA %zu\n",glb.x.n * glb.y.n );
                            fprintf(fp,"SCALARS CD float 1\n");
                            fprintf(fp,"LOOKUP_TABLE default\n");
                            size_t iv=0;
                            for(size_t iy = 0; iy < glb.y.n; iy++) 
                                    for(size_t ix = 0; ix < glb.x.n; ix++)
                                            fprintf(fp,"%11.4e\n",MirImg_PIXEL(*glb.image, iv, ix, iy));
                            
                            fclose(fp);
                          }
                          #if Sp_MIRSUPPORT
                          MirImg_WriteXY( glb.imgf, glb.image, glb.unit->name, scale_factor);
                          Sp_PRINT("Wrote Miriad image to `%s'\n", glb.imgf->name);
                          #endif
                          
                          break;
                  // output dust emission and its polarization image
                  case TASK_CONT:
                          scale_factor = glb.I_norm/glb.ucon;
                          stokes = 1;
                          #if Sp_MIRSUPPORT
                          MirImg_WriteXY(glb.imgf, glb.image, glb.unit->name, scale_factor);
                          Sp_PRINT("Wrote Miriad image to `%s'\n", glb.imgf->name);
                          
                          char SQFName[64],SUFName[64];
                          
                          sprintf( SQFName, "stokesQ_%s", glb.imgf->name);
                          glb.StxQf = MirXY_Open_new(SQFName, glb.x.n, glb.y.n, glb.v.n);
                          MirImg_WriteXY(glb.StxQf, glb.StokesQ, glb.unit->name, scale_factor);
                          Sp_PRINT("Wrote Miriad image to `%s'\n", glb.StxQf->name);
                          
                          sprintf( SUFName, "stokesU_%s", glb.imgf->name);
                          glb.StxUf = MirXY_Open_new(SUFName, glb.x.n, glb.y.n, glb.v.n);
                          MirImg_WriteXY(glb.StxUf, glb.StokesU, glb.unit->name, scale_factor);
                          Sp_PRINT("Wrote Miriad image to `%s'\n", glb.StxUf->name);
                          #endif
                          {
                            char filename[32];
                            sprintf(filename,"stokesIQU_%s.vtk",glb.imgf->name);
                            FILE *fp=fopen(filename,"w");
                            fprintf(fp,"# vtk DataFile Version 3.0\n");
                            fprintf(fp,"Stokes parameters\n");
                            fprintf(fp,"ASCII\n");
                            fprintf(fp,"DATASET STRUCTURED_POINTS\n");
                            fprintf(fp,"DIMENSIONS %zu %zu %d\n",glb.x.n,glb.y.n,1);
                            fprintf(fp,"ORIGIN %f %f %d\n",-glb.x.crpix,-glb.y.crpix,0);
                            fprintf(fp,"SPACING %11.4e %11.4e %d\n",glb.x.delt,glb.y.delt,1);
                            fprintf(fp,"POINT_DATA %zu\n",glb.x.n * glb.y.n);
                            fprintf(fp,"SCALARS Stokes_I float 1\n");
                            fprintf(fp,"LOOKUP_TABLE default\n");
                            for(size_t iy = 0; iy < glb.y.n; iy++) 
                                    for(size_t ix = 0; ix < glb.x.n; ix++) 
                                            fprintf(fp,"%11.4e ",MirImg_PIXEL(*glb.image, 0, ix, iy));
                            fprintf(fp,"\n");
                            fprintf(fp,"SCALARS Stokes_Q float 1\n");
                            fprintf(fp,"LOOKUP_TABLE default\n");
                            for(size_t iy = 0; iy < glb.y.n; iy++) 
                                    for(size_t ix = 0; ix < glb.x.n; ix++) 
                                            fprintf(fp,"%11.4e ",MirImg_PIXEL(*glb.StokesQ, 0, ix, iy));
                            fprintf(fp,"\n");
                            fprintf(fp,"SCALARS Stokes_U float 1\n");
                            fprintf(fp,"LOOKUP_TABLE default\n");
                            for(size_t iy = 0; iy < glb.y.n; iy++) 
                                    for(size_t ix = 0; ix < glb.x.n; ix++) 
                                            fprintf(fp,"%11.4e ",MirImg_PIXEL(*glb.StokesU, 0, ix, iy));
                            fclose(fp);
                            
                            
                            sprintf(filename,"stokesI_%s.dat",glb.imgf->name);
                            fp=fopen(filename,"w"); 
                            for(size_t iy = 0; iy < glb.y.n; iy++) {
                                    for(size_t ix = 0; ix < glb.x.n; ix++) 
                                            fprintf(fp,"%5zu %5zu %11.4e\n",ix,iy,MirImg_PIXEL(*glb.image, 0, ix, iy));
                                    fprintf(fp,"\n");
                            }
                            fclose(fp);
                            sprintf(filename,"stokesQ_%s.dat",glb.imgf->name);
                            fp=fopen(filename,"w"); 
                            for(size_t iy = 0; iy < glb.y.n; iy++) {
                                    for(size_t ix = 0; ix < glb.x.n; ix++) 
                                            fprintf(fp,"%5zu %5zu %11.4e\n",ix,iy,MirImg_PIXEL(*glb.StokesQ, 0, ix, iy));
                                    fprintf(fp,"\n");
                            }
                            fclose(fp);
                            sprintf(filename,"stokesU_%s.dat",glb.imgf->name);
                            fp=fopen(filename,"w"); 
                            for(size_t iy = 0; iy < glb.y.n; iy++) {
                                    for(size_t ix = 0; ix < glb.x.n; ix++) 
                                            fprintf(fp,"%5zu %5zu %11.4e\n",ix,iy,MirImg_PIXEL(*glb.StokesU, 0, ix, iy));
                                    fprintf(fp,"\n");
                            }
                            fclose(fp);
                            
                            sprintf(filename,"vector_%s.dat",glb.imgf->name);
                            fp=fopen(filename,"w");         
                            for(size_t iy = 0; iy < glb.y.n; iy++) {
                                    for(size_t ix = 0; ix < glb.x.n; ix++) {
                                            double StxI = MirImg_PIXEL(*glb.image, 0, ix, iy);
                                            double StxQ = MirImg_PIXEL(*glb.StokesQ, 0, ix, iy);
                                            double StxU = MirImg_PIXEL(*glb.StokesU, 0, ix, iy);
                                            double Sigma2 = MirImg_PIXEL(*glb.sigma2, 0, ix, iy);
                                            double pd = sqrt( StxQ*StxQ + StxU*StxU ) / (StxI - Sigma2);
                                            
                                            double xi=atan2(StxU,StxQ); xi = 0.5 * xi;
                                            double vecx=-pd*sin(xi);
                                            double vecy=pd*cos(xi);
                                            
                                            fprintf(fp,"%5zu %5zu %11.4e %11.4e %11.4e %11.4e\n",ix,iy,vecx,vecy,pd,xi);
                                    }
                                    fprintf(fp,"\n");
                            }
                            fclose(fp);
                          }
                          
                          break;
                  // output line emission or zeeman effect (stokes V) image
                  case TASK_LINE:
                          scale_factor = glb.I_norm/glb.ucon;
                          stokes = 0;
                          #if Sp_MIRSUPPORT
                          MirImg_WriteXY(glb.imgf, glb.image, glb.unit->name, glb.I_norm/glb.ucon);
                          Sp_PRINT("Wrote Miriad image to `%s'\n", glb.imgf->name);
                          #endif
                          /* write excitation visualization to VTK */
                          if(glb.excit && glb.task->idx == TASK_LINE)
                                visualization();
                          #if 1
                          sts = generic_vtk(glb.model.grid->voxel.geom);
                          #endif
                          break;
                  case TASK_ZEEMAN:
                          scale_factor = glb.I_norm/glb.ucon;
                          stokes = 0;
                          
                          #if Sp_MIRSUPPORT
                          MirImg_WriteXY(glb.imgf, glb.image, glb.unit->name, glb.I_norm/glb.ucon);
                          Sp_PRINT("Wrote Miriad image to `%s'\n", glb.imgf->name);
                          #endif
                          
                          break;
                  default:
                          Deb_ASSERT(0);
		}
		FITSoutput( glb.imgf, glb.image, glb.StokesQ, glb.StokesU, glb.unit->name, scale_factor, stokes);
                Sp_PRINT("Wrote FITS image to `%s'\n", glb.imgf->name);
			
		// output tau image
                if(glb.tau_imgf){
			FITSoutput( glb.tau_imgf, glb.tau_img, glb.StokesQ, glb.StokesU, "Optical depth", 1., 0);
                        
                        #if Sp_MIRSUPPORT
                        MirImg_WriteXY(glb.tau_imgf, glb.tau_img, "Optical depth", 1.0);
                        Sp_PRINT("Wrote Miriad image to `%s'\n", glb.tau_imgf->name);
                        MirXY_Close(glb.tau_imgf);
                        #endif
                }
                
        }

	



#if 0	
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
#endif

	

/* 5. Cleanup */
        #if Sp_MIRSUPPORT
        /* Miriad images must always be closed! */
        if(glb.imgf)
                MirXY_Close(glb.imgf);
        if(glb.StxQf)
                MirXY_Close(glb.StxQf);
        if(glb.StxUf)
                MirXY_Close(glb.StxUf);
        #endif
	if(glb.image)
		MirImg_Free(glb.image);
        if(glb.StokesQ)
                MirImg_Free(glb.StokesQ);
        if(glb.StokesU)
                MirImg_Free(glb.StokesU);
        if(glb.sigma2)
		MirImg_Free(glb.sigma2);
	if(glb.tau_img)
		MirImg_Free(glb.tau_img);
	if(glb.subres)
		free(glb.subres);
	
        SpModel_Cleanup(glb.model);

	return sts;
}

/*----------------------------------------------------------------------------*/

static int InitModel(void)
{
	Zone *root = glb.model.grid, *zp;
	SpPhys *pp;
	int sts = 0;

// 	FILE *fp;
// 	double radius;

	/* initialization : construct overlapping table */
        if(glb.overlap){
                glb.model.parms.mol->OL = Mem_CALLOC(NRAD*NRAD,glb.model.parms.mol->OL);
                for(size_t i=0; i<NRAD; i++){
                 for(size_t j=0; j<NRAD; j++){
                        glb.model.parms.mol->OL[NRAD*i+j] = Mem_CALLOC(1,glb.model.parms.mol->OL[NRAD*i+j]);
                        RELVEL(i,j) = ( 1e0-FREQ(j)/FREQ(i) )*CONSTANTS_MKS_LIGHT_C;
                        if( fabs(RELVEL(i,j)) < glb.overlap_vel ){
                                OVERLAP(i,j)=1;
                        }
                        else{
                                OVERLAP(i,j)=0;
                        }
                 }
                }
        }

        /* Set normalization intensity to 20K -- normalization prevents rounding
	   errors from creeping in when flux values are very small */
	if(glb.task->idx != TASK_COLDENS){
		glb.I_norm = Phys_PlanckFunc(glb.freq, 10.0);

		Deb_ASSERT(glb.I_norm > 0); /* Just in case */

		/* Calculate CMB intensity */
		if(glb.model.parms.T_cmb > 0) {
			glb.I_cmb = Phys_PlanckFunc(glb.freq, glb.model.parms.T_cmb);
			Deb_ASSERT(glb.I_cmb > 0); /* Just in case */

			/* Normalize CMB */
			glb.I_cmb /= glb.I_norm;
		}
	}

// 	fp=fopen("pops.dat","w");
	for(zp = Zone_GetMinLeaf(root); zp; zp = Zone_AscendTree(zp)) {
		/* Pointer to physical parameters */
		pp = zp->data;

		if((pp->n_H2 > 1e-200) && !zp->children) {
			/* This is a non-empty leaf zone */
			pp->non_empty_leaf = 1;

			if(pp->X_mol > 0) {
				/* This zone contains tracer molecules */
				pp->has_tracer = 1;
                                /*
 				radius=sqrt(zp->voxel.cen.x[0] * zp->voxel.cen.x[0] + zp->voxel.cen.x[1] * zp->voxel.cen.x[1] + zp->voxel.cen.x[2] * zp->voxel.cen.x[2]);
 				radius=zp->voxel.cen.x[0];
 				fprintf(fp,"%g %g %g %g %g %g %g %g %g %g %g\n",radius,pp->pops[0][0],pp->pops[0][1],pp->pops[0][2],pp->pops[0][3],pp->pops[0][4],pp->pops[0][5],pp->pops[0][6],pp->pops[0][7],pp->pops[0][8],pp->pops[0][9]);
                                */
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
// 	fclose(fp);

	sts = SpUtil_Threads2(Sp_NTHREAD, InitModelThread);
	//SpUtil_Threads(InitModelThread);

	return sts;
}

/*----------------------------------------------------------------------------*/

static void *InitModelThread(void *tid_p)
{
	size_t tid = *((size_t *)tid_p);
        size_t zone_id,j,k;
	Zone *root = glb.model.grid, *zp;
	SpPhys *pp;
        
        if(glb.task->idx != TASK_COLDENS){
          for(zp = Zone_GetMinLeaf(root), zone_id = 0; zp; zp = Zone_AscendTree(zp), zone_id++) {
            if(zone_id % Sp_NTHREAD == tid) {
                /* Check for thread termination */
                Sp_CHECKTERMTHREAD();

                /* Init zone parameters */
                pp = zp->data;

                if(glb.task->idx == TASK_CONT) {
                        SpPhys_InitContWindows(pp, &glb.freq, (size_t)1);
                }
                else {
                        size_t nrad = glb.model.parms.mol->nrad;
                        double freq[nrad];
                        /* Set initial pops to either optically thin or LTE */
                        if(glb.lte){
                                // initialize level population for LTE condition
                                pp->mol = glb.model.parms.mol;
                                for(k = 0; k < Sp_NTHREAD; k++) {
                                        pp->pops[k] = Mem_CALLOC(pp->mol->nlev, pp->pops[k]);
                                }
                                for(j = 0; j < pp->mol->nlev; j++) {
                                        pp->pops[0][j] = SpPhys_BoltzPops(pp->mol, j, pp->T_k);
                                        for(k = 1; k < Sp_NTHREAD; k++) {
                                                pp->pops[k][j] = pp->pops[0][j];
                                        }
                                }
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
                        SpPhys_AddContinuum_d(pp, glb.task->idx == TASK_CONT, glb.model.parms.gas_to_dust);
                }
                /* Add free-free emission/absorption if T_ff > 0 */
                if(pp->T_ff > 0) {
                        SpPhys_AddContinuum_ff(pp, glb.task->idx == TASK_CONT);
                }

                /* Set continuum flux */
                if(pp->T_bb > 0) {
                        //debug
                        SpPhys_SetContinuumIntens_bb(pp, glb.task->idx == TASK_CONT, pp->T_bb, glb.I_norm);
                        Deb_PRINT("T_bb=%g, F_nu=%g\n", pp->T_bb, pp->cont[0].I_bb);
                }
            }
          }
        }

	pthread_exit(NULL);
}

/*----------------------------------------------------------------------------*/

static int CalcImage(void)
{
	switch(glb.task->idx){
                case TASK_ZEEMAN:
                        return SpUtil_Threads2(Sp_NTHREAD, CalcImageThreadZeeman);
                case TASK_CONT:
                        return SpUtil_Threads2(Sp_NTHREAD, CalcImageThreadCont);
                case TASK_COLDENS:
                        return SpUtil_Threads2(Sp_NTHREAD, CalcImageThreadColdens);
                case TASK_LINE:
                        return SpUtil_Threads2(Sp_NTHREAD, CalcImageThreadLine);
                default:
                        /* Shouldn't reach here */
                        Deb_ASSERT(0);
        }
}

/*----------------------------------------------------------------------------*/

static void *CalcImageThreadLine(void *tid_p)
{
	size_t tid = *((size_t *)tid_p);

        /* pix_id is used for distributing work to threads */
        size_t pix_id = 0;
        for(size_t ix = 0; ix < glb.x.n; ix++) {
         for(size_t iy = 0; iy < glb.y.n; iy++) {
                if(pix_id % Sp_NTHREAD == tid) {
                        /* Check for thread termination */
                        Sp_CHECKTERMTHREAD();
                        /* check sub-sampling region */
                        size_t nsub = Init_nsub( ix, iy);

			/* I_nu is the brightness for all channels at pixel (ix, iy) */
			double *I_nu = Mem_CALLOC(glb.v.n, I_nu);
			/* tau_nu is the total optical depth for all channels at pixel (ix, iy) */
			double *tau_nu = Mem_CALLOC(glb.v.n, tau_nu);

			/* Loop through sub-resolution positions */
			for(size_t isub = 0; isub < nsub; isub++) {
			 for(size_t jsub = 0; jsub < nsub; jsub++) {
                                double dx, dy;
				/* Calculate sub-pixel position */
                                InitSubPixel( &dx, &dy, ix, iy, isub, jsub, nsub);
						
                                /* I_sub is the brightness for all channels at each sub-resolution */
				double *I_sub = Mem_CALLOC(glb.v.n, I_sub);
				double *tau_sub = Mem_CALLOC(glb.v.n, tau_sub);

				/* Calculate radiative transfer for this sub-los */
                                if(glb.overlap)
                                        RadiativeXferOverlap(dx, dy, I_sub, tau_sub, tid);
                                else
                                        RadiativeXferLine(dx, dy, I_sub, tau_sub, tid);

				/* Add I_sub to I_nu */
				for(size_t iv = 0; iv < glb.v.n; iv++) {
					I_nu[iv] += I_sub[iv];
					tau_nu[iv] += tau_sub[iv];
				}

				/* Cleanup */
				free(I_sub);
                                free(tau_sub);
			 }
			}
			
			/* Save averaged I_nu to map */
                        double DnsubSquare = 1. / (double)(nsub * nsub);
			for(size_t iv = 0; iv < glb.v.n; iv++) {
				MirImg_PIXEL(*glb.image, iv, ix, iy) = I_nu[iv] * DnsubSquare;
				if(glb.tau_img)
					MirImg_PIXEL(*glb.tau_img, iv, ix, iy) = tau_nu[iv] * DnsubSquare;
			}

                        /* Cleanup */
                        free(I_nu);
                        free(tau_nu);				
		}
                /* Increment pix_id */
                pix_id += 1;
	 }
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/


static void *CalcImageThreadZeeman(void *tid_p)
{
        size_t tid = *((size_t *)tid_p);
        /* pix_id is used for distributing work to threads */
        size_t pix_id = 0;
        for(size_t ix = 0; ix < glb.x.n; ix++) {
         for(size_t iy = 0; iy < glb.y.n; iy++) {
                if(pix_id % Sp_NTHREAD == tid) {
                        /* Check for thread termination */
                        Sp_CHECKTERMTHREAD();
                        /* check sub-sampling region */
                        size_t nsub = Init_nsub( ix, iy);
                        /* I_nu is the brightness for all channels at pixel (ix, iy) */
                        double *V_nu = Mem_CALLOC(glb.v.n, V_nu);
                        /* tau_nu is the total optical depth for all channels at pixel (ix, iy) */
                        double *tau_nu = Mem_CALLOC(glb.v.n, tau_nu);
                        /* Loop through sub-resolution positions */
                        for(size_t isub = 0; isub < nsub; isub++) {
                         for(size_t jsub = 0; jsub < nsub; jsub++) {
                                double dx, dy;
                                /* Calculate sub-pixel position */
                                InitSubPixel( &dx, &dy, ix, iy, isub, jsub, nsub);
                                /* V_sub is Stokes V for all channels at each sub-resolution */
                                double *V_sub = Mem_CALLOC(glb.v.n, V_sub);
                                double *tau_sub = Mem_CALLOC(glb.v.n, tau_sub);
                                /* Calculate radiative transfer for this sub-los (Zeeman) */
                                RadiativeXferZeeman(dx, dy, V_sub, tau_sub, tid);
                                /* Add I_sub to I_nu */
                                for(size_t iv = 0; iv < glb.v.n; iv++) {
                                        V_nu[iv] += V_sub[iv];
                                        tau_nu[iv] += tau_sub[iv];
                                }
                                /* Cleanup */
                                free(V_sub);
                                free(tau_sub);
                         }
                        }
                        /* Save averaged I_nu to map */
                        double DnsubSquare = 1. / (double)(nsub * nsub);
                        for(size_t iv = 0; iv < glb.v.n; iv++) {
                                MirImg_PIXEL(*glb.image, iv, ix, iy) = V_nu[iv] / (double)(nsub * nsub);
                                if(glb.tau_img)
                                        MirImg_PIXEL(*glb.tau_img, iv, ix, iy) = tau_nu[iv] * DnsubSquare;
                        }
                        /* Cleanup */
                        free(V_nu);
                        free(tau_nu);
                }
                /* Increment pix_id */
                pix_id += 1;
         }
        }

        return NULL;
}


/*----------------------------------------------------------------------------*/

static void *CalcImageThreadColdens(void *tid_p)
{
        /* Column density image tracer
           CD is the column density at pixel (ix, iy)
           glb.v.n = 1                                     */
        
        size_t tid = *((size_t *)tid_p);
        /* pix_id is used for distributing work to threads */
        size_t pix_id = 0;
        double *CD_sub = Mem_CALLOC( 1, CD_sub);
        double *CD = Mem_CALLOC( 1, CD);
        for(size_t ix = 0; ix < glb.x.n; ix++) {
         for(size_t iy = 0; iy < glb.y.n; iy++) {
                if(pix_id % Sp_NTHREAD == tid) {
                        /* Check for thread termination */
                        Sp_CHECKTERMTHREAD();
                        /* check sub-sampling region */
                        size_t nsub = Init_nsub( ix, iy);
                        Mem_BZERO2(CD,1);
                        /* Loop through sub-resolution positions */
                        for(size_t isub = 0; isub < nsub; isub++) {
                         for(size_t jsub = 0; jsub < nsub; jsub++) {
                                double dx, dy;
                                /* Calculate sub-pixel position */
                                InitSubPixel( &dx, &dy, ix, iy, isub, jsub, nsub);
                                /* Column density tracer */
                                Mem_BZERO2(CD_sub,1);
                                ColumnDensityTracer(dx,dy, CD_sub);
                                /* Add CD_sub to CD */
                                *CD += *CD_sub;
                         }
                        }
                        double DnsubSquare = 1. / (double)(nsub * nsub);
                        /* Save averaged I_nu to map */
                        MirImg_PIXEL(*glb.image, 0, ix, iy) = (*CD) * DnsubSquare;
                }
         /* Increment pix_id */
         pix_id += 1;
         }
        }

        return NULL;
}

/*----------------------------------------------------------------------------*/

static void *CalcImageThreadCont(void *tid_p)
{
        size_t tid = *((size_t *)tid_p);

         /* pix_id is used for distributing work to threads */
        size_t pix_id = 0;
        for(size_t ix = 0; ix < glb.x.n; ix++) {
         for(size_t iy = 0; iy < glb.y.n; iy++) {
                if(pix_id % Sp_NTHREAD == tid) {
                        /* Check for thread termination */
                        Sp_CHECKTERMTHREAD();
                        /* check sub-sampling region */
                        size_t nsub = Init_nsub( ix, iy);

                        /* I_nu is the brightness for all channels at pixel (ix, iy) */
                        double *I_nu = Mem_CALLOC(glb.v.n, I_nu);
                        double *Q_nu = Mem_CALLOC(glb.v.n, Q_nu);
                        double *U_nu = Mem_CALLOC(glb.v.n, U_nu);
                        double *sigma2 = Mem_CALLOC(glb.v.n, sigma2);
                        /* tau_nu is the total optical depth for all channels at pixel (ix, iy) */
                        double *tau_nu = Mem_CALLOC(glb.v.n, tau_nu);

                        /* Loop through sub-resolution positions */
                        for(size_t isub = 0; isub < nsub; isub++) {
                         for(size_t jsub = 0; jsub < nsub; jsub++) {
                                double dx, dy;
                                /* Calculate sub-pixel position */
                                InitSubPixel( &dx, &dy, ix, iy, isub, jsub, nsub);
                                /* I_sub is the brightness for all channels at each sub-resolution */
                                double *I_sub = Mem_CALLOC(glb.v.n, I_sub);
                                double *Q_sub = Mem_CALLOC(glb.v.n, Q_sub);
                                double *U_sub = Mem_CALLOC(glb.v.n, U_sub);
                                double *sigma2_sub = Mem_CALLOC(glb.v.n, sigma2_sub);
                                double *tau_sub = Mem_CALLOC(glb.v.n, tau_sub);
                                /* Calculate radiative transfer for this sub-los */
                                RadiativeXferCont(dx, dy, I_sub, Q_sub, U_sub, sigma2_sub, tau_sub);
                                /* Add I_sub to I_nu */
                                for(size_t iv = 0; iv < glb.v.n; iv++) {
                                        I_nu[iv] += I_sub[iv];
                                        Q_nu[iv] += Q_sub[iv];
                                        U_nu[iv] += U_sub[iv];
                                        sigma2[iv] += sigma2_sub[iv];
                                        tau_nu[iv] += tau_sub[iv];
                                }
                                /* Cleanup */
                                free(I_sub);
                                free(tau_sub);
                                free(Q_sub);
                                free(U_sub);
                                free(sigma2_sub);
                         }
                        }
                        /* Save averaged I_nu to map */
                        double DnsubSquare = 1. / (double)(nsub * nsub);
                        for(size_t iv = 0; iv < glb.v.n; iv++) {
                                static const double alpha=0.15; // polarized efficiency
                                MirImg_PIXEL(*glb.image, iv, ix, iy) = I_nu[iv] * DnsubSquare;
                                MirImg_PIXEL(*glb.StokesQ, iv, ix, iy) = alpha * Q_nu[iv] * DnsubSquare;
                                MirImg_PIXEL(*glb.StokesU, iv, ix, iy) = alpha * U_nu[iv] * DnsubSquare;
                                MirImg_PIXEL(*glb.sigma2, iv, ix, iy) = 0.5 * alpha * sigma2[iv] * DnsubSquare;
                                if(glb.tau_img)
                                        MirImg_PIXEL(*glb.tau_img, iv, ix, iy) = tau_nu[iv] * DnsubSquare;
                        }
                        /* Cleanup */
                        free(I_nu);
                        free(tau_nu);
                        free(Q_nu);
                        free(U_nu);
                        free(sigma2);
                                
                }
                /* Increment pix_id */
                pix_id += 1;
         }
        }

        return NULL;
}


/*----------------------------------------------------------------------------*/

static void RadiativeXferLine(double dx, double dy, double *I_nu, double *tau_nu, size_t tid)
{
	GeRay ray;
	size_t side;
	Zone *root = glb.model.grid;
        double t;

	InitRay( &dx, &dy, &ray);

	/* Reset tau for all channels */
	Mem_BZERO2(tau_nu, glb.v.n);
	/* Shoot ray at model and see what happens! */
	if(GeRay_IntersectVoxel(&ray, &root->voxel, &t, &side)) {
		/* Calculate intersection */
		ray = GeRay_Inc(&ray, t);
		/* Locate starting leaf zone according to intersection */
		Zone *zp = Zone_GetLeaf(root, side, &ray.e, &ray);
		/* Keep going until there's no next zone to traverse to */
		while(zp) {
			/* Calculate path to next boundary */
			GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &side);
                        #if 0
                        Deb_PRINT("traveling distance = %e\n",t);
                        Deb_PRINT("side = %zu\n",side);
                        #endif

			/* Pointer to physical parameters associated with this zone */
			SpPhys *pp = zp->data;
			/* Do radiative transfer only if gas is present in this zone */
			if(pp->non_empty_leaf) {
				/* Do calculations on all channels at this pixel. Try to minimize the 
				 * amount of operations in this loop, since everything here is repeated 
				 * for ALL channels, and can significantly increase computation time.
				 */
				for(size_t iv = 0; iv < glb.v.n; iv++) {
					/* Calculate velocity associated with this channel */
					double dv = ((double)iv - glb.v.crpix) * glb.v.delt;
					/* Reset emission and absorption coeffs */
					double j_nu = 0;
					double k_nu = 0;

					if(pp->has_tracer) {
						/* Calculate velocity line profile factor for this channel:
                                                 * This version averages over the line profile in steps
                                                 * of the local line width -- very time consuming! */
						double vfac = SpPhys_GetVfac(&ray, t, dv, zp, 0);
						/* Calculate molecular line emission and absorption coefficients */
						SpPhys_GetMoljk(tid, pp, glb.line, vfac, &j_nu, &k_nu);
						
					}
					
                                        /* Add continuum emission/absorption */
					j_nu += pp->cont[glb.line].j;
					k_nu += pp->cont[glb.line].k;

					/* Calculate source function and optical depth if
					 * absorption is NOT zero */
                                        double dtau_nu = k_nu * t * Sp_LENFAC;
                                        double S_nu = (fabs(k_nu) > 0.0) ?
                                                j_nu / ( k_nu * glb.I_norm ) : 0.;

					/* Calculate intensity contributed by this step */
					//debug
					//I_nu[iv] += S_nu * (1.0 - exp(-dtau_nu)) * exp(-tau_nu[iv]);
					I_nu[iv] += (S_nu * (1.0 - exp(-dtau_nu)) + pp->cont[glb.line].I_bb) * exp(-tau_nu[iv]);

					/* Accumulate total optical depth for this channel (must be done
					 * AFTER calculation of intensity!) */
					tau_nu[iv] += dtau_nu;
				}
			}
// 			Deb_PRINT("checkpoint: non_empty_leaf\n");
			/* Calculate next position */
			ray = GeRay_Inc(&ray, t);
// 			Deb_PRINT("checkpoint: GeRay_Inc\n");
			/* Get next zone to traverse to */
			zp = Zone_GetNext(zp, &side, &ray);
			#if 0
			if(zp)
				Deb_PRINT("zone index: %d %d %d\n",GeVec3_X(zp->index,0),GeVec3_X(zp->index,1),GeVec3_X(zp->index,2));
			else
				Deb_PRINT("shoot out!\n");
			#endif
			
			#if 0
			if(I_nu[0] != I_nu[0] && zp ) { // check if intensity become NaN
				Deb_PRINT("zone index: %d %d %d\n",GeVec3_X(zp->index,0),GeVec3_X(zp->index,1),GeVec3_X(zp->index,2));
				Deb_PRINT("dnesity=%g\n",pp->n_H2);
				if(pp->non_empty_leaf) Deb_PRINT("non empty leaf\n");
				
			}
			else if(zp == NULL)
				Deb_PRINT("shoot out!\n");
			#endif
		}
	}

	/* Add CMB to all channels -- this is done even if the ray misses the source */
	for(size_t iv = 0; iv < glb.v.n; iv++)
		I_nu[iv] += glb.I_cmb * exp(-tau_nu[iv]);

	return;
}

/*----------------------------------------------------------------------------*/

static void RadiativeXferOverlap(double dx, double dy, double *I_nu, double *tau_nu, size_t tid)
{
        GeRay ray;
        double t;
        size_t side;
        Zone *root = glb.model.grid;

        InitRay( &dx, &dy, &ray);

        /* Reset tau for all channels */
        Mem_BZERO2(tau_nu, glb.v.n);
        /* Shoot ray at model and see what happens! */
        if(GeRay_IntersectVoxel(&ray, &root->voxel, &t, &side)) {
                /* Calculate intersection */
                ray = GeRay_Inc(&ray, t);
                /* Locate starting leaf zone according to intersection */
                Zone *zp = Zone_GetLeaf(root, side, &ray.e, &ray);
                /* Keep going until there's no next zone to traverse to */
                while(zp) {
                        /* Calculate path to next boundary */
                        GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &side);
//                      Deb_PRINT("checkpoint: GeRay_TraverseVoxel\n");
                        /* Pointer to physical parameters associated with this zone */
                        SpPhys *pp = zp->data;
                        /* Do radiative transfer only if gas is present in this zone */
                        if(pp->non_empty_leaf) {
                                /* Do calculations on all channels at this pixel. Try to minimize the 
                                 * amount of operations in this loop, since everything here is repeated 
                                 * for ALL channels, and can significantly increase computation time.
                                 */
                                for(size_t iv = 0; iv < glb.v.n; iv++) {
                                        /* Calculate velocity associated with this channel */
                                        double dv = ((double)iv - glb.v.crpix) * glb.v.delt;
                                        /* Reset emission and absorption coeffs */
                                        double j_nu = 0;
                                        double k_nu = 0;

                                        if(pp->has_tracer) {
                                                /* Calculate velocity line profile factor for this channel:
                                                 * This version averages over the line profile in steps
                                                 * of the local line width -- very time consuming! */
                                                double vfac = SpPhys_GetVfac(&ray, t, dv, zp, 0);
                                                /* Calculate molecular line emission and absorption coefficients */                                             
                                                size_t i=glb.line;
                                                for(size_t j = 0; j < NRAD; j++) {
                                                        if(OVERLAP(i,j)){
                                                                double tempj_nu, tempk_nu;
                                                                if(i==j){
                                                                        /* Calculate molecular line emission and absorption coefficients */
                                                                        SpPhys_GetMoljk(tid, pp, j, vfac, &tempj_nu, &tempk_nu);
                                                                }
                                                                else{
                                                                        /* Calculate velocity line profile factor */
                                                                        double vfac2 = pp->has_tracer ? SpPhys_GetVfac(&ray, t, dv-RELVEL(i,j), zp, 0) : 0.0;
                                                                        /* Calculate molecular line emission and absorption coefficients */
                                                                        SpPhys_GetMoljk(tid, pp, j, vfac2, &tempj_nu, &tempk_nu);
                                                                }
                                                                j_nu += tempj_nu;
                                                                k_nu += tempk_nu;
                                                       }
                                                }
                                        }
                                        /* Add continuum emission/absorption */
                                        j_nu += pp->cont[glb.line].j;
                                        k_nu += pp->cont[glb.line].k;

                                        /* Calculate source function and optical depth if
                                         * absorption is NOT zero */
                                        double dtau_nu = k_nu * t * Sp_LENFAC;
                                        double S_nu = (fabs(k_nu) > 0.0) ?
                                                j_nu / ( k_nu * glb.I_norm ) : 0.;

                                        /* Calculate intensity contributed by this step */
                                        //debug
                                        I_nu[iv] += (S_nu * (1.0 - exp(-dtau_nu)) + pp->cont[glb.line].I_bb) * exp(-tau_nu[iv]);

                                        /* Accumulate total optical depth for this channel (must be done
                                         * AFTER calculation of intensity!) */
                                        tau_nu[iv] += dtau_nu;
                                }
                        }
                        /* Calculate next position */
                        ray = GeRay_Inc(&ray, t);
                        /* Get next zone to traverse to */
                        zp = Zone_GetNext(zp, &side, &ray);
                }
        }

        /* Add CMB to all channels -- this is done even if the ray misses the source */
        for(size_t iv = 0; iv < glb.v.n; iv++)
                I_nu[iv] += glb.I_cmb * exp(-tau_nu[iv]);

        return;
}

/*----------------------------------------------------------------------------*/

static void RadiativeXferZeeman(double dx, double dy, double *V_nu, double *tau_nu, size_t tid)
{
        GeRay ray;
        double t;
        Zone *root = glb.model.grid;
        InitRay( &dx, &dy, &ray);
        GeVec3_d z,n,e;
        InitLOSCoord( &dx, &dy, &ray, &z, &n, &e);

        /* Reset tau for all channels */
        Mem_BZERO2(tau_nu, glb.v.n);
        
        size_t side;
        /* Shoot ray at model and see what happens! */
        if(GeRay_IntersectVoxel(&ray, &root->voxel, &t, &side)) {
                /* Calculate intersection */
                ray = GeRay_Inc(&ray, t);
                /* Locate starting leaf zone according to intersection */
                Zone *zp = Zone_GetLeaf(root, side, &ray.e, &ray);
                /* Keep going until there's no next zone to traverse to */
                while(zp) {
                        /* Calculate path to next boundary */
                        GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &side);
//                      Deb_PRINT("checkpoint: GeRay_TraverseVoxel\n");
                        /* Pointer to physical parameters associated with this zone */
                        SpPhys *pp = zp->data;
                        /* Do radiative transfer only if gas is present in this zone */
                        if(pp->non_empty_leaf) {
                                /* Do calculations on all channels at this pixel. Try to minimize the 
                                 * amount of operations in this loop, since everything here is repeated 
                                 * for ALL channels, and can significantly increase computation time.
                                 */
    
                                 GeVec3_d B = SpPhys_GetBfac(&ray, t, zp, 0);
                                 double zproduct = GeVec3_DotProd(&z,&B);
                                 double B_Mag = GeVec3_Mag(&B);
                                 double costheta; 
                                 
                                 if (B_Mag == 0.)
                                         // preventing costheta overflow
                                         costheta = 0.;
                                 else
                                         costheta = zproduct / B_Mag;
                                 
                                 static const double g = 2.18/1.4;
                                 double dnu = 1.4e6 * g * B_Mag;
                                 double deltav = dnu * PHYS_CONST_MKS_LIGHTC / glb.freq;
                                 //printf("delta V = %e %e %e %e\n",deltav,dnu,PHYS_CONST_MKS_LIGHTC,glb.freq);

                                for(size_t iv = 0; iv < glb.v.n; iv++) {
                                        /* Calculate velocity associated with this channel */
                                        double dv = ((double)iv - glb.v.crpix) * glb.v.delt;

                                        /* Reset emission and absorption coeffs */
                                        double j_nu = 0;
                                        double k_nu = 0;

                                        if(pp->has_tracer) {
                                                double vfac = SpPhys_GetVfac(&ray, t, dv+deltav, zp, 0);
                                                double vfac2 = SpPhys_GetVfac(&ray, t, dv-deltav, zp, 0);
                                                vfac = vfac-vfac2;
                                                double tempj_nu, tempk_nu;
                                                SpPhys_GetMoljk(tid, pp, glb.line, vfac, &tempj_nu, &tempk_nu);
                                                j_nu = 0.5*tempj_nu;
                                                k_nu += tempk_nu;
                                        }

                                        /* Calculate source function and optical depth if
                                         * absorption is NOT zero */
                                        double dtau_nu = k_nu * t * Sp_LENFAC;
                                        double S_nu = (fabs(k_nu) > 0.0) ?
                                                j_nu / ( k_nu * glb.I_norm ) : 0.;

                                        /* Calculate intensity contributed by this step */
                                        //debug
                                        double temp = S_nu * (1.0 - exp(-dtau_nu)) * exp(-tau_nu[iv]);
                                        V_nu[iv] += temp*costheta;
        
                                        /* Accumulate total optical depth for this channel (must be done
                                         * AFTER calculation of intensity!) */
                                        tau_nu[iv] += dtau_nu;
                                }
                        }
                        /* Calculate next position */
                        ray = GeRay_Inc(&ray, t);
                        //Deb_PRINT("checkpoint: GeRay_Inc\n");
                        /* Get next zone to traverse to */
                        zp = Zone_GetNext(zp, &side, &ray);
                } 
        }

        return;
}

/*----------------------------------------------------------------------------*/

static void RadiativeXferCont(double dx, double dy, double *I_nu, double *Q_nu, double *U_nu, double *sigma2, double *tau_nu)
{
        GeRay ray;
        double t;
        Zone *root = glb.model.grid;
        GeVec3_d z,n,e;        
        
        InitRay( &dx, &dy, &ray);
        InitLOSCoord( &dx, &dy, &ray, &z, &n, &e);

        /* Reset tau for all channels */
        Mem_BZERO2(tau_nu, glb.v.n);

        size_t side;
        /* Shoot ray at model and see what happens! */
        if(GeRay_IntersectVoxel(&ray, &root->voxel, &t, &side)) {
                /* Calculate intersection */
                ray = GeRay_Inc(&ray, t);

                /* Locate starting leaf zone according to intersection */
                Zone *zp = Zone_GetLeaf(root, side, &ray.e, &ray);

                /* Keep going until there's no next zone to traverse to */
                while(zp) {
                        /* Calculate path to next boundary */
                        GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &side);
//                      Deb_PRINT("checkpoint: GeRay_TraverseVoxel\n");
                        /* Pointer to physical parameters associated with this zone */
                        SpPhys *pp = zp->data;
                        /* Do radiative transfer only if gas is present in this zone */
                        if(pp->non_empty_leaf) {
                                /* 
                                   Do calculations on all channels at this pixel. Try to minimize the
                                   amount of operations in this loop, since everything here is repeated
                                   for ALL channels, and can significantly increase computation time.
                                */
                                
                                /* 
                                   Reference : ARTIST(DustPol) -- http://arxiv.org/pdf/1204.6668.pdf
                                */
                                
                                GeVec3_d B = SpPhys_GetBfac(&ray, t, zp, 0);
                                double nproduct = GeVec3_DotProd(&n,&B);
                                double eproduct = GeVec3_DotProd(&e,&B);
                                double zproduct = GeVec3_DotProd(&z,&B);
                                
                                // psi is the angle between the projected B-field on p.o.s. and the north of the image 
                                double B_Mag = GeVec3_Mag(&B);
                                double psi = atan2( -eproduct, nproduct); 
                                // gamma ia the angle bettwen B-field an the plane of sky
                                double cosgammasquare = 1.0 - zproduct * zproduct / GeVec3_Mag(&B); 
 
                                for(size_t iv = 0; iv < glb.v.n; iv++) {
                                        /* Reset emission and absorption coeffs */
                                        double j_nu = 0;
                                        double k_nu = 0;

                                        /* Add continuum emission/absorption */
                                        j_nu += pp->cont[glb.line].j;
                                        k_nu += pp->cont[glb.line].k;

                                        /* Calculate source function and optical depth if
                                         * absorption is NOT zero */
                                        double dtau_nu = k_nu * t * Sp_LENFAC;
                                        double S_nu = (fabs(k_nu) > 0.0) ?
                                                j_nu / ( k_nu * glb.I_norm ) : 0.;

                                        /* Calculate intensity contributed by this step */
                                        //debug
                                        double contribution = (S_nu * (1.0 - exp(-dtau_nu)) + pp->cont[glb.line].I_bb) * exp(-tau_nu[iv]);
                                        
                                        I_nu[iv] += contribution;
                                        if (B_Mag == 0.){
                                                // do nothing
                                                // preventing undefined psi
                                        }
                                        else{
                                                Q_nu[iv] += contribution * cos(2.0 * psi) * cosgammasquare;
                                                U_nu[iv] += contribution * sin(2.0 * psi) * cosgammasquare;
                                        }
                                        static const double d23 = 2. / 3. ;
                                        sigma2[iv] += contribution * (cosgammasquare - d23);

                                        /* Accumulate total optical depth for this channel (must be done
                                         * AFTER calculation of intensity!) */
                                        tau_nu[iv] += dtau_nu;
                                }
                        }
//                      Deb_PRINT("checkpoint: non_empty_leaf\n");
                        /* Calculate next position */
                        ray = GeRay_Inc(&ray, t);
//                      Deb_PRINT("checkpoint: GeRay_Inc\n");
                        /* Get next zone to traverse to */
                        zp = Zone_GetNext(zp, &side, &ray);
                }
        }

        /* Add CMB to all channels -- this is done even if the ray misses the source */
        for(size_t iv = 0; iv < glb.v.n; iv++)
                I_nu[iv] += glb.I_cmb * exp(-tau_nu[iv]);

        return;
}

/*----------------------------------------------------------------------------*/
static void ColumnDensityTracer(double dx, double dy, double *CD)
{
	GeRay ray;
	double t;
	size_t side;
	Zone *root = glb.model.grid;

        InitRay( &dx, &dy, &ray);

//      GeVec3_d z;
// 	z.x[0]=-GeRay_D(ray, 0);
// 	z.x[1]=-GeRay_D(ray, 1);
// 	z.x[2]=-GeRay_D(ray, 2);
        
        Mem_BZERO(CD);
        
	/* Shoot ray at model and see what happens! */
	if(GeRay_IntersectVoxel(&ray, &root->voxel, &t, &side)) {

		/* Calculate intersection */
		ray = GeRay_Inc(&ray, t);

		/* Locate starting leaf zone according to intersection */
		Zone *zp = Zone_GetLeaf(root, side, &ray.e, &ray);

		/* Keep going until there's no next zone to traverse to */
		while(zp) {
			/* Calculate path to next boundary */
			GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &side);

			/* Pointer to physical parameters associated with this zone */
			SpPhys *pp = zp->data;

			/* Do radiative transfer only if gas is present in this zone */
			if(pp->non_empty_leaf) {
				/* Do calculations on all channels at this pixel. Try to minimize the 
				 * amount of operations in this loop, since everything here is repeated 
				 * for ALL channels, and can significantly increase computation time.
				 */ 
// 				GeVec3_d B = SpPhys_GetBfac(&ray, t, zp, 0);
// 				double zproduct = GeVec3_DotProd(&z,&B);
//				GeVec3_d V = SpPhys_GetVfac2(&ray, t, zp, 0);
// 				double zproduct = GeVec3_DotProd(&z,&V);
// 				*CD += pp->n_H2 * pp->X_mol * t * Sp_LENFAC * zproduct;
				*CD += pp->n_H2 * pp->X_mol * t * Sp_LENFAC;
			}
			/* Calculate next position */
			ray = GeRay_Inc(&ray, t);

			/* Get next zone to traverse to */
			zp = Zone_GetNext(zp, &side, &ray);
		}
	}

	return;
}


/*---------------------------------------------------------------------------- */
/* Initialize the position of the ray and its shooting direction */
static void InitRay(double *dx, double *dy, GeRay *ray)
{       
        Zone *root = glb.model.grid;

        /* Reset ray */
        Mem_BZERO(ray);

        /* Init ray position to <dist, 0, 0> */
        GeRay_E(*ray, 0) = glb.dist / Sp_LENFAC;
        GeRay_E(*ray, 1) = 0;
        GeRay_E(*ray, 2) = 0;

        /* Set direction of ray according to pixel position:
         *   theta = PI/2 + dy
         *   phi = -dx
         */
        double phi = -(*dx);
        double theta = 0.5 * M_PI + (*dy);

        /* Convert to Cartesian coordinates */
        GeRay_D(*ray, 0) = sin(theta) * cos(phi);
        GeRay_D(*ray, 1) = sin(theta) * sin(phi);
        GeRay_D(*ray, 2) = cos(theta);

        /* Rotate ray:
         * Since what we REALLY want to rotate is the model and that the rays
         * are pointed towards the model, rays should be rotated in the negative
         * direction, and in the opposite order of what we would've done to
         * rotate the model. */
        *ray = GeRay_Rotate(ray, 2, -glb.rotate[2]);
        *ray = GeRay_Rotate(ray, 1, -glb.rotate[1]);
        *ray = GeRay_Rotate(ray, 0, -glb.rotate[0]);

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
/* initialize line-of-sight coordinate */
static void InitLOSCoord( double *dx, double *dy, GeRay *ray, GeVec3_d *z, GeVec3_d *n, GeVec3_d *e)
{
        double phi = M_PI + *dx;   
        double theta = 0.5 * M_PI - *dy;
        
        /* line of sight coordinate */
        GeVec3_X(*z,0) = -GeRay_D(*ray, 0);
        GeVec3_X(*z,1) = -GeRay_D(*ray, 1);
        GeVec3_X(*z,2) = -GeRay_D(*ray, 2);
        
        GeVec3_X(*n,0) = -cos(theta)*cos(phi);
        GeVec3_X(*n,1) = -cos(theta)*sin(phi);
        GeVec3_X(*n,2) =  sin(theta);
        *n = GeVec3_Rotate_z(n, -glb.rotate[2]);
        *n = GeVec3_Rotate_y(n, -glb.rotate[1]);
        *n = GeVec3_Rotate_x(n, -glb.rotate[0]);
        
        GeVec3_X(*e,0) = sin(phi);
        GeVec3_X(*e,1) = -cos(phi);
        GeVec3_X(*e,2) = 0.0;
        *e = GeVec3_Rotate_z(e, -glb.rotate[2]);
        *e = GeVec3_Rotate_y(e, -glb.rotate[1]);
        *e = GeVec3_Rotate_x(e, -glb.rotate[0]);
        
        return;
}

/*----------------------------------------------------------------------------*/
/* determine nsub if the pixel is in the sub-sampling position */
static size_t Init_nsub(size_t ix, size_t iy)
{
        /* Determine angular offsets from the pointing center */
        double dx = ((int)ix - (int)glb.x.crpix) * glb.x.delt;
        double dy = ((int)iy - (int)glb.y.crpix) * glb.y.delt;

        /* nsub is by default 1 */
        size_t nsub = 1;

        /* Check if position is within any of the subres boxes */
        for(size_t ibox = 0; ibox < glb.nsubres; ibox++) {
                if((dx >= glb.subres[ibox].blc_x) && (dx <= glb.subres[ibox].trc_x) &&
                        (dy >= glb.subres[ibox].blc_y) && (dy <= glb.subres[ibox].trc_y)) {
                        /* Position within subres box, set nsub to subres[ibox].nsub */
                        nsub = glb.subres[ibox].nsub;
                        break;
                }
        }
        
        return nsub;
}

/*----------------------------------------------------------------------------*/
static void InitSubPixel( double *dx, double *dy, 
                          size_t ix, size_t iy, 
                          size_t isub,size_t jsub, 
                          size_t nsub)
{
        /* Determine sub-resolution angular offsets from the pointing center */
        double subix = (double)ix + (2.*(double)isub-(double)nsub+1.) / ((double)(2*nsub));
        double subiy = (double)iy + (2.*(double)jsub-(double)nsub+1.) / ((double)(2*nsub));
        *dx = (subix - (double)glb.x.crpix) * glb.x.delt;
        *dy = (subiy - (double)glb.y.crpix) * glb.y.delt;

        return;
}


/*----------------------------------------------------------------------------*/
/*   
write out VTK file of the excitation temperature for the nested Cartesian Cell in multiblock format 
*/
static void visualization(void)
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
					double n_u = pp->pops[0][up];
					double n_l = pp->pops[0][lo];
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
					double n_u = pp->pops[0][up];
					double n_l = pp->pops[0][lo];
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

static int generic_vtk(int geom)
{
        int sts = 0;
        
        glb.visual = Mem_CALLOC(1, glb.visual);
        switch (geom){
                case GEOM_SPH1D:
                        sts = vtk_sph1d();
                        break;
                case GEOM_SPH3D:
                        sts = vtk_sph3d();
                        break;
                case GEOM_REC3D:
                        sts = vtk_rec3d();
                        break;
                case GEOM_CYL3D:
                        sts = vtk_cyl3d();
                        break;
                default:
                        /* Should not happen */
                        Deb_ASSERT(0);
        }
        free(glb.visual);
        
        return sts;
}

/*----------------------------------------------------------------------------*/

static int vtk_sph1d(void)
{
        int sts = 0;
        Zone * root = glb.model.grid;
        size_t nvelo = glb.v.n;
        
        glb.visual->sph3d = Mem_CALLOC(1,glb.visual->sph3d);
        // Dimension of the visualized resolution
        size_t nr = glb.visual->sph3d->nr = root->nchildren;
        size_t nt = glb.visual->sph3d->nt = 45;
        size_t np = glb.visual->sph3d->np = 90;
        
        size_t nelement =  nr * np * nt;
        // declare the memory
        double * radius = Mem_CALLOC( nr+1, radius);
        double * theta = Mem_CALLOC( nt+1, theta);
        double * phi = Mem_CALLOC( np+1, phi);
        double ** contrib = Mem_CALLOC( nelement, contrib);
        double ** tau = Mem_CALLOC( nelement, tau);
        double ** tau_dev = Mem_CALLOC( nelement, tau_dev);
        for (size_t idx = 0; idx < nelement; idx++){
                contrib[idx] = Mem_CALLOC(nvelo, contrib[idx]);
                tau[idx]     = Mem_CALLOC(nvelo, tau[idx]);
                tau_dev[idx] = Mem_CALLOC(nvelo, tau_dev[idx]);
        }
        // link to the global pointer
        glb.visual->sph3d->radius = radius;
        glb.visual->sph3d->theta = theta;
        glb.visual->sph3d->phi = phi;
        glb.visual->contrib = contrib;
        glb.visual->tau = tau;
        glb.visual->tau_dev = tau_dev;
        
        // construct the expanding SPH3D mesh
        radius[0] = root->children[0]->voxel.min.x[0];
        for(size_t i = 1; i < nr+1; i++)
                radius[i] = root->children[i-1]->voxel.max.x[0];
        double delta_theta = M_PI / (double) nt;
        theta[0] = 0.;
        for ( size_t j = 1; j < nt+1; j++)
                theta[j] = theta[j-1] + delta_theta;
        double delta_phi = 2.0 * M_PI / (double) np;
        phi[0] = 0.;
        for (size_t k = 1; k < np+1; k ++)
                phi[k] = phi[k-1] + delta_phi;

        sts = CalcContrib(GEOM_SPH1D);     
        
        // open VTK file
        FILE *fp;
        char filename[32];
        sprintf(filename,"vis.vtk");
        fp=fopen(filename,"w");
        
        // write the header
        fprintf(fp,"# vtk DataFile Version 3.0\n");
        fprintf(fp,"%s\n", "POSTPROCESSING VISUALIZATION");
        fprintf(fp,"ASCII\n");
        
        // define the type of the gridding
        fprintf(fp,"DATASET STRUCTURED_GRID\n");
        fprintf(fp,"DIMENSIONS %zu %zu %zu\n", np+1, nt+1, nr+1);
        fprintf(fp,"POINTS %zu float\n", (nr+1) * (nt+1) * (np+1) );
        for( size_t i = 0; i < nr + 1; i++)
         for( size_t j = 0; j < nt + 1; j++)
          for( size_t k = 0; k < np + 1; k++){
                double x = radius[i] * sin(theta[j]) * cos(phi[k]);
                double y = radius[i] * sin(theta[j]) * sin(phi[k]);
                double z = radius[i] * cos(theta[j]);
                fprintf(fp,"%E %E %E\n", x, y, z);
          }
        
        // write the artributes
        fprintf(fp,"CELL_DATA %zu\n", nelement );
        
        // H2 number density
        fprintf(fp,"SCALARS H2_Number_Density float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++){
                Zone *zp = root->children[i];
                SpPhys *pp = zp->data;
                for( size_t j = 0; j < nt; j++)
                  for( size_t k = 0; k < np; k++){
                        fprintf(fp,"%E ", pp->n_H2);
                  }fprintf(fp,"\n");
        }
        
        // molecular number density
        fprintf(fp,"SCALARS Molecular_Number_Density float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++){
                Zone *zp = root->children[i];
                SpPhys *pp = zp->data;
                for( size_t j = 0; j < nt; j++)
                  for( size_t k = 0; k < np; k++){
                        fprintf(fp,"%E ", pp->n_H2 * pp->X_mol);
                  }fprintf(fp,"\n");
        }
        
        // kinetic temperature
        fprintf(fp,"SCALARS Temperature float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++){
                Zone *zp = root->children[i];
                SpPhys *pp = zp->data;
                for( size_t j = 0; j < nt; j++)
                  for( size_t k = 0; k < np; k++){
                        fprintf(fp,"%E ", pp->T_k);
                  }fprintf(fp,"\n");
        }
        
        // exitation temperature
        fprintf(fp,"SCALARS Tex float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++){
                Zone *zp = root->children[i];
                SpPhys *pp = zp->data;
                if(pp->X_mol == 0.0){
                        for( size_t j = 0; j < nt; j++)
                          for( size_t k = 0; k < np; k++)
                                fprintf(fp,"%E ", 0.0);
                }
                else{
                        MolTrRad *trans = pp->mol->rad[glb.line];
                        size_t up = trans->up;
                        size_t lo = trans->lo;
                        double n_u = pp->pops[0][up];
                        double n_l = pp->pops[0][lo];
                        double E_u = pp->mol->lev[up]->E;
                        double E_l = pp->mol->lev[lo]->E;
                        double g_u = pp->mol->lev[up]->g;
                        double g_l = pp->mol->lev[lo]->g;
                        double Tex = (E_l-E_u)
                                / ( PHYS_CONST_MKS_BOLTZK * log((n_u*g_l)/(n_l*g_u)) );
                        double Tex_ratio = Tex/(pp->T_k);
                        for( size_t j = 0; j < nt; j++)
                          for( size_t k = 0; k < np; k++){
                                fprintf(fp,"%E ", Tex_ratio);
                          }fprintf(fp,"\n");
                }
        }
        
        fprintf(fp,"FIELD CONTRIBUTION_TAU 3\n");
        // write the contribution of the cells
        fprintf(fp,"%s %zu %zu float\n", VISUAL_TYPES[0].name, nvelo, nelement);
        for (size_t idx = 0; idx < nelement; idx++){
                for ( size_t l = 0; l < nvelo; l++)
                        fprintf(fp,"%E ", contrib[idx][l]);
                fprintf(fp,"\n");
        }
        
        fprintf(fp,"%s %zu %zu float\n", VISUAL_TYPES[1].name, nvelo, nelement);
        for (size_t idx = 0; idx < nelement; idx++){
                for ( size_t l = 0; l < nvelo; l++)
                        fprintf(fp,"%E ", tau[idx][l]);
                fprintf(fp,"\n");
        }
        fprintf(fp,"%s %zu %zu float\n", VISUAL_TYPES[2].name, nvelo, nelement);
        for (size_t idx = 0; idx < nelement; idx++){
                for ( size_t l = 0; l < nvelo; l++)
                        fprintf(fp,"%E ", tau_dev[idx][l]);
                fprintf(fp,"\n");
        }
        fclose(fp);
        printf("wrote %s\n",filename);
        

        for (size_t idx = 0; idx < nelement; idx++){
                free(contrib[idx]);
                free(tau[idx]);
                free(tau_dev[idx]);
        }
        free(contrib);
        free(tau);
        free(tau_dev);
        free(radius);
        free(theta);
        free(phi);
        
        free(glb.visual->sph3d);

        return sts;
}

/*----------------------------------------------------------------------------*/

static int vtk_sph3d(void){
        int sts = 0;
        Zone * root = glb.model.grid;
        size_t nvelo = glb.v.n;
        
        glb.visual->sph3d = Mem_CALLOC(1,glb.visual->sph3d);
        // Dimension of the visualized resolution
        size_t nr = glb.visual->sph3d->nr = root->naxes.x[0];
        size_t nt = glb.visual->sph3d->nt = (root->naxes.x[1] == 1) ? 45 : root->naxes.x[1];
        size_t np = glb.visual->sph3d->np = (root->naxes.x[2] == 1) ? 90 : root->naxes.x[2];
        
        size_t nelement =  nr * np * nt;
        // declare the memory
        double * radius = Mem_CALLOC( nr+1, radius);
        double * theta = Mem_CALLOC( nt+1, theta);
        double * phi = Mem_CALLOC( np+1, phi);
        double ** contrib = Mem_CALLOC( nelement, contrib);
        double ** tau = Mem_CALLOC( nelement, tau);
        double ** tau_dev = Mem_CALLOC( nelement, tau_dev);
        for (size_t idx = 0; idx < nelement; idx++){
                contrib[idx] = Mem_CALLOC(nvelo, contrib[idx]);
                tau[idx]     = Mem_CALLOC(nvelo, tau[idx]);
                tau_dev[idx] = Mem_CALLOC(nvelo, tau_dev[idx]);
        }
        
        // link to the global pointer
        glb.visual->sph3d->radius = radius;
        glb.visual->sph3d->theta = theta;
        glb.visual->sph3d->phi = phi;
        glb.visual->contrib = contrib;
        glb.visual->tau = tau;
        glb.visual->tau_dev = tau_dev;
        
        // construct the expanding SPH3D mesh
        radius[0] = root->children[0]->voxel.min.x[0];
        for(size_t i = 1; i < nr+1; i++)
                radius[i] = root->children[(i-1)*root->naxes.x[1]*root->naxes.x[2]]->voxel.max.x[0];
        if (root->naxes.x[1] == 1){
                double delta_theta = M_PI / (double) nt;
                theta[0] = 0.;
                for ( size_t j = 1; j < nt+1; j++)
                        theta[j] = theta[j-1] + delta_theta;
        }
        else{
                theta[0] = root->children[0]->voxel.min.x[1];
                for ( size_t j = 1; j < nt+1; j++)
                        theta[j] = root->children[ (j-1)*root->naxes.x[2] ]->voxel.max.x[1];
        }
        if (root->naxes.x[2] == 1){
                double delta_phi = 2. * M_PI / (double) np;
                phi[0] = 0.;
                for (size_t k = 1; k < np+1; k ++)
                        phi[k] = phi[k-1] + delta_phi;
        }
        else{
                phi[0] = root->children[0]->voxel.min.x[2];
                for (size_t k = 1; k < np+1; k ++)
                        phi[k] = root->children[k-1]->voxel.max.x[2];
        }
        #if 0
        for(size_t i = 0; i < nr+1; i++) printf("%E ",radius[i]);printf("\n");
        for(size_t j = 0; j < nr+1; j++) printf("%E ",theta[j]);printf("\n");
        for(size_t k = 0; k < nr+1; k++) printf("%E ",phi[k]);printf("\n");
        printf("%zu %zu %zu\n",root->naxes.x[0],root->naxes.x[1],root->naxes.x[2]);
        printf("%zu %zu %zu\n",nr,nt,np);                                 
        #endif
         /* multithreading calculation of the contribution of the cells */
        sts = CalcContrib(GEOM_SPH3D);

        // open VTK file
        FILE *fp;
        char filename[32];
        sprintf(filename,"vis.vtk");
        fp=fopen(filename,"w");
        
        // write the header
        fprintf(fp,"# vtk DataFile Version 3.0\n");
        fprintf(fp,"%s\n", "POSTPROCESSING VISUALIZATION");
        fprintf(fp,"ASCII\n");
        
        // define the type of the gridding
        fprintf(fp,"DATASET STRUCTURED_GRID\n");
        fprintf(fp,"DIMENSIONS %zu %zu %zu\n", np+1, nt+1, nr+1);
        fprintf(fp,"POINTS %zu float\n", (nr+1) * (nt+1) * (np+1) );
        for( size_t i = 0; i < nr + 1; i++)
         for( size_t j = 0; j < nt + 1; j++)
          for( size_t k = 0; k < np + 1; k++){
                double x = radius[i] * sin(theta[j]) * cos(phi[k]);
                double y = radius[i] * sin(theta[j]) * sin(phi[k]);
                double z = radius[i] * cos(theta[j]);
                fprintf(fp,"%E %E %E\n", x, y, z);
          }
        // write the artributes
        fprintf(fp,"CELL_DATA %zu\n", nelement );
        
        // H2 number density
        fprintf(fp,"SCALARS H2_Number_Density float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++)
          for( size_t j = 0; j < nt; j++)
            for( size_t k = 0; k < np; k++){
                // izone : to coresponding zone index
                size_t izone = 
                  (root->naxes.x[1]==1) ?
                     (root->naxes.x[2]==1) ?
                             i :
                             i * root->naxes.x[2] + k :
                     (root->naxes.x[2]==1) ?
                             i * root->naxes.x[1] + j :
                             (i * root->naxes.x[1] + j) * root->naxes.x[2] + k;
                Zone *zp = root->children[izone];
                SpPhys *pp = zp->data;
                fprintf(fp,"%E ", pp->n_H2);
        }fprintf(fp,"\n");
        
        // molecular number density
        fprintf(fp,"SCALARS Molecular_Number_Density float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++)
          for( size_t j = 0; j < nt; j++)
            for( size_t k = 0; k < np; k++){
                // izone : to coresponding zone index
                size_t izone = 
                  (root->naxes.x[1]==1) ?
                     (root->naxes.x[2]==1) ?
                             i :
                             i * root->naxes.x[2] + k :
                     (root->naxes.x[2]==1) ?
                             i * root->naxes.x[1] + j :
                             (i * root->naxes.x[1] + j) * root->naxes.x[2] + k;
                Zone *zp = root->children[izone];
                SpPhys *pp = zp->data;
                fprintf(fp,"%E ", pp->n_H2 * pp->X_mol);
        }fprintf(fp,"\n");
        
        // kinetic temperature
        fprintf(fp,"SCALARS Temperature float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++)
          for( size_t j = 0; j < nt; j++)
            for( size_t k = 0; k < np; k++){
                // izone : to coresponding zone index
                size_t izone = 
                  (root->naxes.x[1]==1) ?
                     (root->naxes.x[2]==1) ?
                             i :
                             i * root->naxes.x[2] + k :
                     (root->naxes.x[2]==1) ?
                             i * root->naxes.x[1] + j :
                             (i * root->naxes.x[1] + j) * root->naxes.x[2] + k;
                Zone *zp = root->children[izone];
                SpPhys *pp = zp->data;
                fprintf(fp,"%E ", pp->T_k);
        }fprintf(fp,"\n");
        
        // exitation temperature
        fprintf(fp,"SCALARS Tex float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++)
          for( size_t j = 0; j < nt; j++)
            for( size_t k = 0; k < np; k++){
                // izone : to coresponding zone index
                size_t izone = 
                  (root->naxes.x[1]==1) ?
                     (root->naxes.x[2]==1) ?
                             i :
                             i * root->naxes.x[2] + k :
                     (root->naxes.x[2]==1) ?
                             i * root->naxes.x[1] + j :
                             (i * root->naxes.x[1] + j) * root->naxes.x[2] + k;
                Zone *zp = root->children[izone];
                SpPhys *pp = zp->data;
                if(pp->X_mol == 0.0){
                        fprintf(fp,"%E ", 0.0);
                }
                else{
                        MolTrRad *trans = pp->mol->rad[glb.line];
                        size_t up = trans->up;
                        size_t lo = trans->lo;
                        double n_u = pp->pops[0][up];
                        double n_l = pp->pops[0][lo];
                        double E_u = pp->mol->lev[up]->E;
                        double E_l = pp->mol->lev[lo]->E;
                        double g_u = pp->mol->lev[up]->g;
                        double g_l = pp->mol->lev[lo]->g;
                        double Tex = (E_l-E_u)
                                / ( PHYS_CONST_MKS_BOLTZK * log((n_u*g_l)/(n_l*g_u)) ); 
                        fprintf(fp,"%E ", Tex/(pp->T_k));
                }
        }fprintf(fp,"\n");
        
        fprintf(fp,"FIELD CONTRIBUTION_TAU 3\n");
        // write the contribution of the cells
        fprintf(fp,"%s %zu %zu float\n", VISUAL_TYPES[0].name, nvelo, nelement);
        for (size_t idx = 0; idx < nelement; idx++){
                for ( size_t l = 0; l < nvelo; l++)
                        fprintf(fp,"%E ", contrib[idx][l]);
                fprintf(fp,"\n");
        }
        
        fprintf(fp,"%s %zu %zu float\n", VISUAL_TYPES[1].name, nvelo, nelement);
        for (size_t idx = 0; idx < nelement; idx++){
                for ( size_t l = 0; l < nvelo; l++)
                        fprintf(fp,"%E ", tau[idx][l]);
                fprintf(fp,"\n");
        }
        fprintf(fp,"%s %zu %zu float\n", VISUAL_TYPES[2].name, nvelo, nelement);
        for (size_t idx = 0; idx < nelement; idx++){
                for ( size_t l = 0; l < nvelo; l++)
                        fprintf(fp,"%E ", tau_dev[idx][l]);
                fprintf(fp,"\n");
        }
        fclose(fp);
        printf("wrote %s\n",filename);  
        
        for (size_t idx = 0; idx < nelement; idx++){
                free(contrib[idx]);
                free(tau[idx]);
                free(tau_dev[idx]);
        }
        free(contrib);
        free(tau);
        free(tau_dev);
        free(radius);
        free(theta);
        free(phi);
        free(glb.visual->sph3d);
        return sts;
}

/*----------------------------------------------------------------------------*/

static int vtk_rec3d(void){
        int sts = 0;
        return sts;
}

/*----------------------------------------------------------------------------*/

static int vtk_cyl3d(void){
        int sts = 0;
        Zone * root = glb.model.grid;
        size_t nvelo = glb.v.n;
        
        glb.visual->cyl3d = Mem_CALLOC(1,glb.visual->cyl3d);
        // Dimension of the visualized resolution
        size_t nr = glb.visual->cyl3d->nr = root->naxes.x[0];
        size_t np = glb.visual->cyl3d->np = (root->naxes.x[1] == 1) ? 72 : root->naxes.x[1];
        size_t nz = glb.visual->cyl3d->nz = root->naxes.x[2];
        
        size_t nelement =  nr * np * nz;
        // declare the memory
        double * Rc = Mem_CALLOC( nr+1, Rc);
        double * phi = Mem_CALLOC( np+1, phi);
        double * Z = Mem_CALLOC( nz+1, Z);
        double ** contrib = Mem_CALLOC( nelement, contrib);
        double ** tau = Mem_CALLOC( nelement, tau);
        double ** tau_dev = Mem_CALLOC( nelement, tau_dev);
        for (size_t idx = 0; idx < nelement; idx++){
                contrib[idx] = Mem_CALLOC(nvelo, contrib[idx]);
                tau[idx]     = Mem_CALLOC(nvelo, tau[idx]);
                tau_dev[idx] = Mem_CALLOC(nvelo, tau_dev[idx]);
        }
        
        // link to the global pointer
        glb.visual->cyl3d->Rc   = Rc;
        glb.visual->cyl3d->phi  = phi;
        glb.visual->cyl3d->Z    = Z;
        glb.visual->contrib     = contrib;
        glb.visual->tau         = tau;
        glb.visual->tau_dev     = tau_dev;
        
        // construct the expanding CYL3D mesh
        // for Rc
        Rc[0] = root->children[0]->voxel.min.x[0];
        for(size_t i = 1; i < nr+1; i++)
                Rc[i] = root->children[ (i-1) * root->naxes.x[1] * root->naxes.x[2] ]->voxel.max.x[0];
        // for phi
        if (root->naxes.x[1] == 1){
                double delta_phi = 2.0 * M_PI / (double) np;
                phi[0] = 0.;
                for ( size_t j = 1; j < np+1; j++)
                        phi[j] = phi[j-1] + delta_phi;
        }
        else{
                phi[0] = root->children[0]->voxel.min.x[1];
                for ( size_t j = 1; j < np+1; j++)
                        phi[j] = root->children[ (j-1)*root->naxes.x[2] ]->voxel.max.x[1];
        }
        // for Z
        Z[0] = root->children[0]->voxel.min.x[2];
        for (size_t k = 1; k < nz+1; k ++)
                Z[k] = root->children[k-1]->voxel.max.x[2];

        #if 0
        for(size_t i = 0; i < nr+1; i++) printf("%E ",Rc[i]);printf("\n");
        for(size_t j = 0; j < np+1; j++) printf("%E ",phi[j]);printf("\n");
        for(size_t k = 0; k < nz+1; k++) printf("%E ",Z[k]);printf("\n");
        printf("%zu %zu %zu\n",root->naxes.x[0],root->naxes.x[1],root->naxes.x[2]);
        printf("%zu %zu %zu\n",nr,np,nz);                                 
        #endif

        sts = CalcContrib(GEOM_CYL3D);

        // open VTK file
        FILE *fp;
        char filename[32];
        sprintf(filename,"vis.vtk");
        fp=fopen(filename,"w");
        
        // write the header
        fprintf(fp,"# vtk DataFile Version 3.0\n");
        fprintf(fp,"%s\n", "POSTPROCESSING VISUALIZATION");
        fprintf(fp,"ASCII\n");
        
        // define the type of the gridding
        fprintf(fp,"DATASET STRUCTURED_GRID\n");
        fprintf(fp,"DIMENSIONS %zu %zu %zu\n", nz+1, np+1, nr+1);
        fprintf(fp,"POINTS %zu float\n", (nr+1) * (np+1) * (nz+1) );
        for( size_t i = 0; i < nr + 1; i++)
         for( size_t j = 0; j < np + 1; j++)
          for( size_t k = 0; k < nz + 1; k++){
                double x = Rc[i] * cos(phi[j]);
                double y = Rc[i] * sin(phi[j]);
                double z = Z[k];
                fprintf(fp,"%E %E %E\n", x, y, z);
          }
        // write the artributes
        fprintf(fp,"CELL_DATA %zu\n", nelement );
        
        // H2 number density
        fprintf(fp,"SCALARS H2_Number_Density float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++) 
         for( size_t j = 0; j < np; j++){ 
          for( size_t k = 0; k < nz; k++){
                // izone : to coresponding zone index
                size_t izone = (root->naxes.x[1]==1) ?
                        i * root->naxes.x[2] + k :
                        (i * root->naxes.x[1] + j) * root->naxes.x[2] + k;
                Zone *zp = root->children[izone];
                SpPhys *pp = zp->data;
                fprintf(fp,"%E ", pp->n_H2);
           }
           fprintf(fp,"\n");
         }
        // molecular number density
        fprintf(fp,"SCALARS Molecular_Number_Density float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++)
         for( size_t j = 0; j < np; j++)
          for( size_t k = 0; k < nz; k++){
                // izone : to coresponding zone index
                size_t izone = (root->naxes.x[1]==1) ?
                        i * root->naxes.x[2] + k :
                        (i * root->naxes.x[1] + j) * root->naxes.x[2] + k;
                Zone *zp = root->children[izone];
                SpPhys *pp = zp->data;
                fprintf(fp,"%E ", pp->n_H2 * pp->X_mol);
        }fprintf(fp,"\n");
        
        // kinetic temperature
        fprintf(fp,"SCALARS Temperature float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++)
         for( size_t j = 0; j < np; j++){
          for( size_t k = 0; k < nz; k++){
                // izone : to coresponding zone index
                size_t izone = (root->naxes.x[1]==1) ?
                        i * root->naxes.x[2] + k :
                        (i * root->naxes.x[1] + j) * root->naxes.x[2] + k;
                Zone *zp = root->children[izone];
                SpPhys *pp = zp->data;
                fprintf(fp,"%E ", pp->T_k);
           }
           fprintf(fp,"\n");
         }
        
        // exitation temperature
        fprintf(fp,"SCALARS Tex float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < nr; i++)
         for( size_t j = 0; j < np; j++){
          for( size_t k = 0; k < nz; k++){
                // izone : to coresponding zone index
                size_t izone = (root->naxes.x[1]==1) ?
                        i * root->naxes.x[2] + k :
                        (i * root->naxes.x[1] + j) * root->naxes.x[2] + k;
                Zone *zp = root->children[izone];
                SpPhys *pp = zp->data;
                if(pp->X_mol == 0.0){
                        fprintf(fp,"%E ", 0.0);
                }
                else{
                        MolTrRad *trans = pp->mol->rad[glb.line];
                        size_t up = trans->up;
                        size_t lo = trans->lo;
                        double n_u = pp->pops[0][up];
                        double n_l = pp->pops[0][lo];
                        double E_u = pp->mol->lev[up]->E;
                        double E_l = pp->mol->lev[lo]->E;
                        double g_u = pp->mol->lev[up]->g;
                        double g_l = pp->mol->lev[lo]->g;
                        double Tex = (E_l-E_u)
                                / ( PHYS_CONST_MKS_BOLTZK * log((n_u*g_l)/(n_l*g_u)) ); 
                        fprintf(fp,"%E ", Tex/(pp->T_k));
                }
           }
           fprintf(fp,"\n");
         }
        
        fprintf(fp,"FIELD CONTRIBUTION_TAU 3\n");
        // write the contribution of the cells
        fprintf(fp,"%s %zu %zu float\n", VISUAL_TYPES[0].name, nvelo, nelement);
        for (size_t idx = 0; idx < nelement; idx++){
                for ( size_t l = 0; l < nvelo; l++)
                        fprintf(fp,"%E ", contrib[idx][l]);
                fprintf(fp,"\n");
        }
        
        fprintf(fp,"%s %zu %zu float\n", VISUAL_TYPES[1].name, nvelo, nelement);
        for (size_t idx = 0; idx < nelement; idx++){
                for ( size_t l = 0; l < nvelo; l++)
                        fprintf(fp,"%E ", tau[idx][l]);
                fprintf(fp,"\n");
        }
        fprintf(fp,"%s %zu %zu float\n", VISUAL_TYPES[2].name, nvelo, nelement);
        for (size_t idx = 0; idx < nelement; idx++){
                for ( size_t l = 0; l < nvelo; l++)
                        fprintf(fp,"%E ", tau_dev[idx][l]);
                fprintf(fp,"\n");
        }        
        fclose(fp);
        printf("wrote %s\n",filename);  
        
        for (size_t idx = 0; idx < nelement; idx++){
                free(contrib[idx]);
                free(tau[idx]);
                free(tau_dev[idx]);
        }
        free(contrib);
        free(tau);
        free(tau_dev);
        free(Rc);
        free(phi);
        free(Z);
        free(glb.visual->cyl3d);

        return sts;
}

/*----------------------------------------------------------------------------*/

static int CalcContrib(GEOM_TYPE geom_type)
{
        /* multithreading calculation of the contribution of the cells */
        switch (geom_type){
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
        size_t nvelo = glb.v.n;
        
        // Dimension of the visualized resolution
        size_t nr = glb.visual->sph3d->nr;
        size_t nt = glb.visual->sph3d->nt;
        size_t np = glb.visual->sph3d->np;
        
        // link to the global pointer
        double * radius = glb.visual->sph3d->radius;
        double * theta = glb.visual->sph3d->theta;
        double * phi = glb.visual->sph3d->phi;
        double ** contrib = glb.visual->contrib;
        double ** tau = glb.visual->tau;
        double ** tau_dev = glb.visual->tau_dev;
        
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
                        ContributionSubSamp( contrib[idx], tau[idx], tau_dev[idx], &SampZone, nvelo, tid);
                        
                        #if 0
                        printf("zone pos = %zu %zu %zu\n",i,j,k);
                        //printf("%E %E %E\n",contrib[idx + 15], tau[idx + 15], tau_dev[idx + 15]);
                        if( j==3 && k==41){                                
                                //printf("OK\n");exit(0);
                        }
                        #endif
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
        size_t nvelo = glb.v.n;
        
        // Dimension of the visualized resolution
        size_t nr = glb.visual->sph3d->nr;
        size_t nt = glb.visual->sph3d->nt;
        size_t np = glb.visual->sph3d->np;
        
        // link to the global pointer
        double * radius = glb.visual->sph3d->radius;
        double * theta = glb.visual->sph3d->theta;
        double * phi = glb.visual->sph3d->phi;
        double ** contrib = glb.visual->contrib;
        double ** tau = glb.visual->tau;
        double ** tau_dev = glb.visual->tau_dev;
        
        size_t cell_id = 0;
        // calculate the contribution of the cells
        for( size_t i = 0; i < nr; i++){
          for( size_t j = 0; j < nt; j++){
            for( size_t k = 0; k < np; k++){
              if ( cell_id % Sp_NTHREAD == tid ){
                        /* Check for thread termination */
                        Sp_CHECKTERMTHREAD();
                        // izone : to coresponding zone index
                        size_t izone = 
                          (root->naxes.x[1]==1) ?
                             (root->naxes.x[2]==1) ?
                                     i :
                                     i * root->naxes.x[2] + k :
                             (root->naxes.x[2]==1) ?
                                     i * root->naxes.x[1] + j :
                                     (i * root->naxes.x[1] + j) * root->naxes.x[2] + k;
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
                        ContributionSubSamp( contrib[idx], tau[idx], tau_dev[idx], &SampZone, nvelo, tid);
                        
                        #if 0
                        printf("zone pos = %zu %zu %zu\n",i,j,k);
                        //printf("%E %E %E\n",contrib[idx + 15], tau[idx + 15], tau_dev[idx + 15]);
                        if( j==3 && k==41){                                
                                //printf("OK\n");exit(0);
                        }
                        #endif
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
        size_t nvelo = glb.v.n;
        
        // Dimension of the visualized resolution
        size_t nr = glb.visual->cyl3d->nr;
        size_t np = glb.visual->cyl3d->np;
        size_t nz = glb.visual->cyl3d->nz;
        
        // link to the global pointer
        double * Rc = glb.visual->cyl3d->Rc;
        double * phi = glb.visual->cyl3d->phi;
        double * Z = glb.visual->cyl3d->Z;
        double ** contrib = glb.visual->contrib;
        double ** tau = glb.visual->tau;
        double ** tau_dev = glb.visual->tau_dev;
        
        size_t cell_id = 0;
        // calculate the contribution of the cells
        for( size_t i = 0; i < nr; i++){
          for( size_t j = 0; j < np; j++){
            for( size_t k = 0; k < nz; k++){
              if ( cell_id % Sp_NTHREAD == tid ){
                        /* Check for thread termination */
                        Sp_CHECKTERMTHREAD();
                        // izone : to coresponding zone index
                        size_t izone = (root->naxes.x[1]==1) ?
                                     i * root->naxes.x[2] + k :
                                     (i * root->naxes.x[1] + j) * root->naxes.x[2] + k;
                        // copy grid i_zone to sampling zone
                        Zone SampZone = *root->children[izone];
                        // the voxel pointer
                        GeVox *vp = &SampZone.voxel;
                        // change the GEOM of the voxel to SPH3D
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
                        ContributionSubSamp( contrib[idx], tau[idx], tau_dev[idx], &SampZone, nvelo, tid);
                        
                        #if 0
                        printf("zone pos = %zu %zu %zu\n",i,j,k);
                        printf("%E %E %E\n",contrib[idx + 15], tau[idx + 15], tau_dev[idx + 15]);
                        if( j==0 && k==0){                                
                                //printf("OK\n");exit(0);
                        }
                        #endif
              }
              cell_id += 1;
            }
          }
        }
    
        //pthread_exit(NULL);
        return NULL;
}

/*----------------------------------------------------------------------------*/

void ContributionSubSamp(double *contrib, double *tau, double *tau_dev, Zone *SampZone, size_t nvelo, size_t tid){
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
        
        switch (vp->geom){
            case GEOM_SPH3D:{
                double R_in = vp->min.x[0];
                double R_out = vp->max.x[0];
                double theta_in = vp->min.x[1];
                double theta_out = vp->max.x[1];
                double phi_in = vp->min.x[2];
                double phi_out = vp->max.x[2];
                Mem_BZERO2(contrib, nvelo); 
                Mem_BZERO2(tau, nvelo);
                Mem_BZERO2(tau_dev, nvelo);
                for(int i = 0; i < nSamp1D; i++){
                        // radius position
                        double RFrac = (double) ( 2 * i + 1 ) * Devide_2nSamp1D;
                        #define ONE_THIRD (0.3333333333333333)
                        double R = ( R_in == 0. ) ? 
                                R_out * pow(RFrac, ONE_THIRD) : 
                                R_in * pow((1. + RFrac * (pow(R_out/R_in, 3.) - 1.)), ONE_THIRD);
                        #undef ONE_THIRD
                        for(int j = 0; j < nSamp1D; j++){
                                // theta position
                                double ThetaFrac = (double) ( 2 * j + 1 ) * Devide_2nSamp1D;
                                double theta = ( theta_in == 0. ) ?
                                        acos( 1. + ThetaFrac * ( cos(theta_out) - 1. ) ) :
                                        acos( ( cos(theta_in) - 1. ) * ( 1. + ThetaFrac * ( ( cos(theta_out) - 1. ) / ( cos(theta_in) - 1. ) - 1. ) ) + 1. );
                                for(int k = 0; k < nSamp1D; k++){
                                        // phi position
                                        double PhiFrac = (double) ( 2 * k + 1 ) * Devide_2nSamp1D;
                                        double phi = phi_in + PhiFrac * ( phi_out - phi_in );
                                        
                                        // Sampling position in Cartesian (x,y,z) format
                                        GeVec3_d SampPosXYZ = GeVec3_Sph2Cart( R, theta, phi);
                                        
                                        // Call the contribution tracer
                                        double * contrib_sub = Mem_CALLOC( nvelo, contrib_sub); 
                                        double * tau_sub = Mem_CALLOC( nvelo, tau_sub);

                                        ContributionTracer( contrib_sub, tau_sub, SampZone, &SampPosXYZ, tid);

                                        // statistics : sum up
                                        for (size_t l = 0; l < nvelo; l++){
                                                contrib[l] += contrib_sub[l];
                                                tau[l] += tau_sub[l];
                                                tau_dev[l] += tau_sub[l] * tau_sub[l];
                                        }
                                        
                                        free(contrib_sub);
                                        free(tau_sub);
                                }
                        }
                }
                break;
            }
            case GEOM_CYL3D:{
                double Rc_in = vp->min.x[0];
                double Rc_out = vp->max.x[0];
                double phi_in = vp->min.x[1];
                double phi_out = vp->max.x[1];
                double Z_in = vp->min.x[2];
                double Z_out = vp->max.x[2];
                Mem_BZERO2(contrib, nvelo); 
                Mem_BZERO2(tau, nvelo);
                Mem_BZERO2(tau_dev, nvelo);
                for(int i = 0; i < nSamp1D; i++){
                        // radius position
                        double RcFrac = (double) ( 2 * i + 1 ) * Devide_2nSamp1D;
                        double Rc = ( Rc_in == 0. ) ? 
                                Rc_out * pow(RcFrac, 0.5) : 
                                Rc_in * pow((1. + RcFrac * (pow(Rc_out/Rc_in, 2.0) - 1.)), 0.5);
                        for(int j = 0; j < nSamp1D; j++){
                                // phi position
                                double PhiFrac = (double) ( 2 * j + 1 ) * Devide_2nSamp1D;
                                double phi = phi_in + PhiFrac * ( phi_out - phi_in );
                                for(int k = 0; k < nSamp1D; k++){
                                        // Z position
                                        double ZFrac = (double) ( 2 * k + 1 ) * Devide_2nSamp1D;
                                        double Z = Z_in + ZFrac * ( Z_out - Z_in );
                                        
                                        // Sampling position in Cartesian (x,y,z) format
                                        GeVec3_d SampPosXYZ = GeVec3_Cyl2Cart( Rc, phi, Z);
                                        
                                        // Call the contribution tracer
                                        double * contrib_sub = Mem_CALLOC( nvelo, contrib_sub); 
                                        double * tau_sub = Mem_CALLOC( nvelo, tau_sub);

                                        ContributionTracer( contrib_sub, tau_sub, SampZone, &SampPosXYZ, tid);
                                        
                                        // statistics : sum up
                                        for (size_t l = 0; l < nvelo; l++){
                                                contrib[l] += contrib_sub[l];
                                                tau[l] += tau_sub[l];
                                                tau_dev[l] += tau_sub[l] * tau_sub[l];
                                        }
                                        
                                        free(contrib_sub);
                                        free(tau_sub);
                                }
                        }
                }
                break;    
            }    
        }
        
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

        return;
}

/*----------------------------------------------------------------------------*/


static void ContributionTracer( double *contrib, double *tau_nu, Zone *SampZone, GeVec3_d *SampCartPos, size_t tid)
{
        GeVox * SampVp = &SampZone->voxel;
        GeRay ray;
        size_t side;
        Zone * root = glb.model.grid;
        double t;
        int reached_sampling_zone = 0;
        static double threshold = 1E-6;

        double dx = atan( SampCartPos->x[1] / ( glb.dist / Sp_LENFAC - SampCartPos->x[0] ) );
        double dy = atan( SampCartPos->x[2] / ( glb.dist / Sp_LENFAC - SampCartPos->x[0] ) );
        InitRay( &dx, &dy, &ray);
        
        /* Reset tau for all channels */
        Mem_BZERO2(tau_nu, glb.v.n);
        /* Shoot ray at model and see what happens! */
        if(GeRay_IntersectVoxel(&ray, &root->voxel, &t, &side)) {
            /* Calculate intersection */
            ray = GeRay_Inc(&ray, t);
            /* Locate starting leaf zone according to intersection */
            Zone *zp = Zone_GetLeaf(root, side, &ray.e, &ray);

            /* Keep going until there's no next zone to traverse to */
            while(zp) {
                // see if the ray reach the target 1-D layer
                if( zp->pos == SampZone->pos ){
                        // the ordinates of the ray position
                        GeVec3_d RayCartPos = ray.e;
                        GeVec3_d RayGeomPos;
                        switch (SampVp->geom){
                            case GEOM_SPH3D:
                                RayGeomPos = GeVec3_Cart2Sph(&RayCartPos);
                                break;
                            case GEOM_CYL3D:
                                RayGeomPos = GeVec3_Cart2Cyl(&RayCartPos);
                                break;
                            default:
                                Deb_ASSERT(0);
                        }
                        #if 0
                        // debugging
                        //printf("dx \t= %E, dy = %E\n",dx,dy);
                        printf("ray.e \t= %E %E %E\n", ray.e.x[0], ray.e.x[1], ray.e.x[2]);
                        //printf("ray.d \t= %E %E %E\n", ray.d.x[0], ray.d.x[1], ray.d.x[2]);
                        printf("SampPos \t= %E %E %E\n", 
                        SampCartPos->x[0], SampCartPos->x[1], SampCartPos->x[2]);
                        printf("VoxelMin \t= %E %E %E\n", 
                        SampVp->min.x[0], SampVp->min.x[1], SampVp->min.x[2]);
                        printf("RayGeomPos \t= %E %E %E\n", 
                        RayGeomPos.x[0], RayGeomPos.x[1], RayGeomPos.x[2]);
                        printf("VoxelMax \t= %E %E %E\n", 
                        SampVp->max.x[0], SampVp->max.x[1], SampVp->max.x[2]);
                        #endif

                        // see if the photon reach the sampling cell
                        int reached_cell = 1;
                        for (size_t i = 0; i < 3; i++){
                                size_t side_in = 2 * i;
                                size_t side_out = side_in + 1;
                                if ( side_in == side ){
                                        reached_cell *= 
                                        ( fabs( SampVp->min.x[i] - RayGeomPos.x[i] ) / RayGeomPos.x[i] < threshold) ? 
                                        1 : 0;
                                        //printf("diff_in = %E\n",fabs( SampVp->min.x[i] - RayGeomPos.x[i] ));
                                        
                                }else{
                                        reached_cell *= 
                                        (SampVp->min.x[i] <= RayGeomPos.x[i]) ? 1 : 0;
                                }
                                if ( side_out == side ){
                                        reached_cell *= 
                                        ( fabs( SampVp->max.x[i] - RayGeomPos.x[i] ) / RayGeomPos.x[i] < threshold ) ? 
                                        1 : 0;
                                        //printf("diff_out = %E\n",fabs( SampVp->max.x[i] - RayGeomPos.x[i] ));
                                }
                                else{
                                        reached_cell *= 
                                        (SampVp->max.x[i] >= RayGeomPos.x[i]) ? 1 : 0;
                                }
                        }
                        
                        // if reached the sampling cell, escape the tau tracer.
                        if( reached_cell ){
                                reached_sampling_zone = 1;
                                break;
                        }
                        // the ray is not inside the sampling voxel but inside the target zone
                        else{
                                size_t side_Samp;
                                double tSamp;
                                
                                int hit;
                                switch (SampVp->geom){
                                    case GEOM_SPH3D :
                                        hit = HitSph3dVoxel( &ray, SampVp, &tSamp, &side_Samp);
                                        break;
                                    case GEOM_CYL3D :
                                        hit = HitCyl3dVoxel( &ray, SampVp, &tSamp, &side_Samp);
                                        break;
                                    default:
                                            Deb_ASSERT(0);
                                
                                }
                                GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &side);

                                if ( tSamp < t){
                                        CalcOpticalDepth( zp, &ray, tSamp, tau_nu, tid);
                                        /* Calculate next position */
                                        ray = GeRay_Inc(&ray, tSamp);
                                        side = side_Samp;
                                }
                                else{
                                        // will not reach in this interval
                                        CalcOpticalDepth( zp, &ray, t, tau_nu, tid);
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
                        
                        CalcOpticalDepth( zp, &ray, t, tau_nu, tid);
                        
                        /* Calculate next position */
                        ray = GeRay_Inc(&ray, t);
                        /* Get next zone to traverse to */
                        zp = Zone_GetNext(zp, &side, &ray);
                }
            } 
            ContributionOfCell( zp, &ray, SampCartPos, contrib, tau_nu, tid);
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
        Deb_ASSERT(reached_sampling_zone);

        return;
}



/*----------------------------------------------------------------------------*/

static void CalcOpticalDepth( Zone *zp, const GeRay *ray, const double t, double *tau_nu, size_t tid)
{       /* Pointer to physical parameters associated with this zone */
        SpPhys * pp = zp->data;
        
        /* Do radiative transfer only if gas is present in this zone */
        if(pp->non_empty_leaf) {
                /* Do calculations on all channels at this pixel. Try to minimize the 
                * amount of operations in this loop, since everything here is repeated 
                * for ALL channels, and can significantly increase computation time.
                */
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
                                SpPhys_GetMoljk(tid, pp, glb.line, vfac, &j_nu, &k_nu);
                                
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

static int HitSph3dVoxel( const GeRay *ray, const GeVox *voxel, double *t, size_t *side)
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

static int HitCyl3dVoxel( const GeRay *ray, const GeVox *voxel, double *t, size_t *side)
{
        static const double half_pi = 0.5 * 3.1415926535897932384626433832795;
        
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

static void ContributionOfCell( Zone *zp, const GeRay *ray, const GeVec3_d *SampCartPos, double * contrib, double *tau_nu, size_t tid)
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
                                SpPhys_GetMoljk(tid, pp, glb.line, vfac, &j_nu, &k_nu);
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
                        contrib[iv] = (S_nu * (1.0 - exp(-dtau_nu)) + pp->cont[glb.line].I_bb)
                                * exp(-tau_nu[iv]) / t;

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
