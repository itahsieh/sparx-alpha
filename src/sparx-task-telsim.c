#include "sparx.h"
#include "task.h"
#include "unit.h"

/* Miriad support */
#define Sp_MIRSUPPORT MIRSUPPORT

/* Global parameter struct */
static struct glb {
    DatINode *task;
    
    int overlap,lte;
    int tracer;
    
    DatINode *unit;
    double ucon, overlap_vel;
    
    double I_norm, I_cmb, I_in;
    SpModel model;
    
    size_t line;
    double lamb, freq;

    SpTelsim tel_parms;
    
    MirImg_Axis x, y, v;
    MirFile *imgf, *StxQf, *StxUf, *tau_imgf;
    MirImg *image, *StokesQ, *StokesU, *StokesV, *sigma2, *tau_img;
} glb;



#define RELVEL(i,j) glb.model.parms.mol->OL[NRAD*(i)+(j)]->RelativeVel
#define OVERLAP(i,j) glb.model.parms.mol->OL[NRAD*(i)+(j)]->overlap

#define NRAD glb.model.parms.mol->nrad
#define FREQ(i) glb.model.parms.mol->rad[(i)]->freq


/* Subroutine prototypes */
static int InitModel(void);
static void *InitModelThread(void *tid_p);
static int CalcImage(void);
static void *CalcImageThreadLine(void *tid_p);
static void *CalcImageThreadZeeman(void *tid_p);
static void *CalcImageThreadCont(void *tid_p);
static void RadiativeXferLine(double dx, double dy, double *I_nu, double *tau_nu, size_t tid);
static void RadiativeXferOverlap(double dx, double dy, double *I_nu, double *tau_nu, size_t tid);
static void RadiativeXferZeeman(double dx, double dy, double *V_nu, double *tau_nu, size_t tid);
static void RadiativeXferCont(double dx, double dy, double *I_nu, double *Q_nu, double *U_nu, double *sigma2, double *tau_nu);

static void IntensityBC( size_t side, double *I_nu, double *tau_nu, GeRay *ray);
static void InitLOSCoord( double *dx, double *dy, GeRay *ray, GeVec3_d *z, GeVec3_d *n, GeVec3_d *e);


/*----------------------------------------------------------------------------*/

int SpTask_Telsim(void)
{
    Mem_BZERO(&glb);
    int sts = 0;
    
    // 1. GET PARAMETERS FROM PYTHON
    PyObject *o;
    /*    1-1 those are the general parameters */
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
    if(!sts) sts = SpPy_GetInput_dbl("dist", &glb.tel_parms.dist);
    glb.tel_parms.dist /= Sp_LENFAC;
    
    /* rotate */
    if(!sts && !(sts = SpPy_GetInput_PyObj("rotate", &o))) {
        glb.tel_parms.rotate[0] = Sp_PYDBL(Sp_PYLST(o, 0));
        glb.tel_parms.rotate[1] = Sp_PYDBL(Sp_PYLST(o, 1));
        glb.tel_parms.rotate[2] = Sp_PYDBL(Sp_PYLST(o, 2));
        SpPy_XDECREF(o);
    }
    /* subres */
    if(!sts && !(sts = SpPy_GetInput_PyObj("subres", &o))) {
        if(o != Py_None) {
            /* Get number of boxes */
            glb.tel_parms.nsubres = (size_t)PyList_Size(o);
            glb.tel_parms.subres = Mem_CALLOC(glb.tel_parms.nsubres, glb.tel_parms.subres);
            for(size_t i = 0; i < glb.tel_parms.nsubres; i++) {
                PyObject *o_sub;
                o_sub = PyList_GetItem(o, (Py_ssize_t)i);
                glb.tel_parms.subres[i].blc_x = Sp_PYDBL(PyList_GetItem(o_sub, (Py_ssize_t)0));
                glb.tel_parms.subres[i].blc_y = Sp_PYDBL(PyList_GetItem(o_sub, (Py_ssize_t)1));
                glb.tel_parms.subres[i].trc_x = Sp_PYDBL(PyList_GetItem(o_sub, (Py_ssize_t)2));
                glb.tel_parms.subres[i].trc_y = Sp_PYDBL(PyList_GetItem(o_sub, (Py_ssize_t)3));
                glb.tel_parms.subres[i].nsub = Sp_PYSIZE(PyList_GetItem(o_sub, (Py_ssize_t)4));
            }
        }
        SpPy_XDECREF(o);
    }
    
    
    /*    1-2 get the task-based parameters */
    /* obs */
    if(!sts && !(sts = SpPy_GetInput_PyObj("obs", &o))) {
        
        /* task */
        PyObject *o_task;
        o_task = PyObject_GetAttrString(o, "task");
        glb.task = Dat_IList_NameLookup(TASKS, Sp_PYSTR(o_task));
        SpPy_XDECREF(o_task);	
        
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
            case TASK_LINE:
            case TASK_ZEEMAN:
                { // get line transition
                    PyObject *o_line;
                    o_line = PyObject_GetAttrString(o, "line");
                    glb.line = Sp_PYSIZE(o_line);
                    SpPy_XDECREF(o_line);
                }
                if(!sts) sts = SpPy_GetInput_bool("lte", &glb.lte);
                
                /* tau (optional) */
                if(!sts && SpPy_CheckOptionalInput("tau")) {
                    sts = SpPy_GetInput_mirxy_new("tau", glb.x.n, glb.y.n, glb.v.n, &glb.tau_imgf);
                }
                
                { // get overlap velocity
                    PyObject *o3;
                    o3 = PyObject_GetAttrString(o, "overlap_vel");
                    glb.overlap_vel = Sp_PYDBL(o3);
                    if (glb.overlap_vel == 0.0)
                        glb.overlap = 0;
                    else
                        glb.overlap = 1;
                    SpPy_XDECREF(o3);
                }
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
        SpPy_XDECREF(o);
        
    }
    
    /*    1-3 read the source model */
    /* source */
    if(!sts) {
        int task_id = glb.task->idx; 
        int popsold = 0;
        /* only read populations data when non-LTE LINE/ZEEMAN mapping tasks*/ 
        if( task_id == TASK_LINE || task_id == TASK_ZEEMAN ){
            if (!glb.lte)
                popsold = 1;
        }
        sts = SpPy_GetInput_model("source","source", &glb.model, &popsold, task_id);
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
        if(glb.tau_imgf){
            glb.tau_img = MirImg_Alloc(glb.x, glb.y, glb.v);
            glb.tau_img->restfreq = glb.freq;
        }
        /* Calculate image */
        sts = CalcImage();
    }
    
    /* 4. I/O : OUTPUT */
    if(!sts){
        double scale_factor = 1.0;
        int stokes;
        switch(glb.task->idx){
            // output dust emission and its polarization image
            case TASK_CONT:
                scale_factor = glb.I_norm/glb.ucon;
                stokes = 1;
                
                #if Sp_MIRSUPPORT
                MirImg_WriteXY(glb.imgf, glb.image, glb.unit->name, scale_factor);
                Sp_PRINT("Wrote Miriad image to `%s'\n", glb.imgf->name);
                #endif
                
                
                if (glb.model.parms.polariz){
                    
                    #if Sp_MIRSUPPORT
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
                            
                            double pd = sqrt( StxQ*StxQ + StxU*StxU ) / StxI;
                            
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
            case TASK_ZEEMAN:
                scale_factor = glb.I_norm/glb.ucon;
                stokes = 0;
                #if Sp_MIRSUPPORT
                MirImg_WriteXY(glb.imgf, glb.image, glb.unit->name, scale_factor);
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
            Sp_PRINT("Wrote FITS image to `%s'\n", glb.tau_imgf->name);
            
            #if Sp_MIRSUPPORT
            MirImg_WriteXY(glb.tau_imgf, glb.tau_img, "Optical depth", 1.0);
            Sp_PRINT("Wrote Miriad image to `%s'\n", glb.tau_imgf->name);
            MirXY_Close(glb.tau_imgf);
            #endif
        }
    }
    
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
    if(glb.tel_parms.subres)
        free(glb.tel_parms.subres);
    
    SpModel_Cleanup(glb.model);
    
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
    if( task_id == TASK_LINE || task_id == TASK_ZEEMAN ){
        Deb_ASSERT(glb.line < parms->mol->nrad);
        glb.freq = parms->mol->rad[glb.line]->freq;
        glb.lamb = PHYS_CONST_MKS_LIGHTC / glb.freq;
    }
    /* set the reference of the intensity: glb.ucon */
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
    /* initialization : construct overlapping table */
    if(glb.overlap){
        parms->mol->OL = Mem_CALLOC(NRAD*NRAD,parms->mol->OL);
        for(size_t i=0; i<NRAD; i++){
            for(size_t j=0; j<NRAD; j++){
                parms->mol->OL[NRAD*i+j] = Mem_CALLOC(1, parms->mol->OL[NRAD*i+j]);
                RELVEL(i,j) = ( 1e0-FREQ(j)/FREQ(i) )*CONSTANTS_MKS_LIGHT_C;
                if( fabs(RELVEL(i,j)) < glb.overlap_vel ){
                    OVERLAP(i,j)=1;
                    //if ( i != j && i == glb.line )
                    //printf("%zu %zu %e %e %e\n",i,j,FREQ(i),FREQ(j),RELVEL(i,j));
                }
                else{
                    OVERLAP(i,j)=0;
                }
            }
        }
    }
    
    /* Set normalization intensity to 20K -- normalization prevents rounding
     *	   errors from creeping in when flux values are very small */
    glb.I_norm = Phys_PlanckFunc(glb.freq, 10.0);
    Deb_ASSERT(glb.I_norm > 0); /* Just in case */
    
    /* Calculate CMB intensity */
    if(parms->T_cmb > 0) {
        glb.I_cmb = Phys_PlanckFunc(glb.freq, parms->T_cmb) / glb.I_norm;
        Deb_ASSERT(glb.I_cmb > 0); /* Just in case */
    }
    /* Calculate inner boundary intensity */
    if(parms->T_in > 0) {
        glb.I_in = Phys_PlanckFunc(glb.freq, parms->T_in) / glb.I_norm;
        Deb_ASSERT(glb.I_in > 0); /* Just in case */
    }
    /* Calculate Source intensity and the dim factor */
    int nSource = parms->Outer_Source;
    if(nSource){
        for (int source_id = 0; source_id < nSource; source_id++){
            SourceData *source = &parms->source[source_id];
            
            source->beta = source->radius / source->distance;
            
            source->intensity = Mem_CALLOC( 1, source->intensity);
            source->intensity[0] = Phys_PlanckFunc( glb.freq, source->temperature) / glb.I_norm;
            
            GeVec3_X(source->pt_sph,0) = source->distance;
            GeVec3_X(source->pt_sph,1) = source->theta;
            GeVec3_X(source->pt_sph,2) = source->phi;
            source->pt_cart = GeVec3_Sph2Cart(&source->pt_sph);
        }
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
    
    for(zp = Zone_GetMinLeaf(root), zone_id = 0; zp; zp = Zone_AscendTree(zp), zone_id++) {
        if(zone_id % Sp_NTHREAD == tid) {
            /* Check for thread termination */
            Sp_CHECKTERMTHREAD();
            
            /* Init zone parameters */
            pp = zp->data;
            
            if(task_id == TASK_CONT) {
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
                /*
                    *                    if (!(pp->width > 0.0)) {
                    *                        printf("width = %g pos=%zu temp=%g\n",pp->width,zp->pos,pp->T_k);
                    *                        exit(0);
            }
            */
            }
            
            /* Add dust emission/absorption if T_d > 0 */
            if(pp->T_d > 0) {
                SpPhys_AddContinuum_d(pp, task_id == TASK_CONT, pp->dust_to_gas);
            }
            /* Add free-free emission/absorption if T_ff > 0 */
            if(pp->X_e > 0) {
                SpPhys_AddContinuum_ff(pp, task_id == TASK_CONT);
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
                size_t nsub = SpImgTrac_Init_nsub( ix, iy, &glb.tel_parms, &glb.x, &glb.y);
                
                /* I_nu is the brightness for all channels at pixel (ix, iy) */
                double *I_nu = Mem_CALLOC(glb.v.n, I_nu);
                /* tau_nu is the total optical depth for all channels at pixel (ix, iy) */
                double *tau_nu = Mem_CALLOC(glb.v.n, tau_nu);
                
                /* Loop through sub-resolution positions */
                for(size_t isub = 0; isub < nsub; isub++) {
                    for(size_t jsub = 0; jsub < nsub; jsub++) {
                        double dx, dy;
                        /* Calculate sub-pixel position */
                        SpImgTrac_InitSubPixel( &dx, &dy, ix, iy, isub, jsub, nsub, &glb.x, &glb.y);
                        
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
                size_t nsub = SpImgTrac_Init_nsub( ix, iy, &glb.tel_parms, &glb.x, &glb.y);
                /* I_nu is the brightness for all channels at pixel (ix, iy) */
                double *V_nu = Mem_CALLOC(glb.v.n, V_nu);
                /* tau_nu is the total optical depth for all channels at pixel (ix, iy) */
                double *tau_nu = Mem_CALLOC(glb.v.n, tau_nu);
                /* Loop through sub-resolution positions */
                for(size_t isub = 0; isub < nsub; isub++) {
                    for(size_t jsub = 0; jsub < nsub; jsub++) {
                        double dx, dy;
                        /* Calculate sub-pixel position */
                        SpImgTrac_InitSubPixel( &dx, &dy, ix, iy, isub, jsub, nsub, &glb.x, &glb.y);
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

static void *CalcImageThreadCont(void *tid_p)
{
    size_t tid = *((size_t *)tid_p);
    size_t nvelo = glb.v.n;
    
    /* pix_id is used for distributing work to threads */
    size_t pix_id = 0;
    for(size_t ix = 0; ix < glb.x.n; ix++) {
        for(size_t iy = 0; iy < glb.y.n; iy++) {
            if(pix_id % Sp_NTHREAD == tid) {
                /* Check for thread termination */
                Sp_CHECKTERMTHREAD();
                /* check sub-sampling region */
                size_t nsub = SpImgTrac_Init_nsub( ix, iy, &glb.tel_parms, &glb.x, &glb.y);
                
                /* I_nu is the brightness for all channels at pixel (ix, iy) */
                double *I_nu = Mem_CALLOC(nvelo, I_nu);
                double *Q_nu = Mem_CALLOC(nvelo, Q_nu);
                double *U_nu = Mem_CALLOC(nvelo, U_nu);
                double *sigma2 = Mem_CALLOC(nvelo, sigma2);
                /* tau_nu is the total optical depth for all channels at pixel (ix, iy) */
                double *tau_nu = Mem_CALLOC(nvelo, tau_nu);
                
                /* Loop through sub-resolution positions */
                for(size_t isub = 0; isub < nsub; isub++) {
                    for(size_t jsub = 0; jsub < nsub; jsub++) {
                        double dx, dy;
                        /* Calculate sub-pixel position */
                        SpImgTrac_InitSubPixel( &dx, &dy, ix, iy, isub, jsub, nsub, &glb.x, &glb.y);
                        /* I_sub is the brightness for all channels at each sub-resolution */
                        double *I_sub = Mem_CALLOC(nvelo, I_sub);
                        double *Q_sub = Mem_CALLOC(nvelo, Q_sub);
                        double *U_sub = Mem_CALLOC(nvelo, U_sub);
                        double *sigma2_sub = Mem_CALLOC(nvelo, sigma2_sub);
                        double *tau_sub = Mem_CALLOC(nvelo, tau_sub);
                        /* Calculate radiative transfer for this sub-los */
                        RadiativeXferCont(dx, dy, I_sub, Q_sub, U_sub, sigma2_sub, tau_sub);
                        /* Add I_sub to I_nu */
                        for(size_t iv = 0; iv < nvelo; iv++) {
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
                for(size_t iv = 0; iv < nvelo; iv++) {
                    MirImg_PIXEL(*glb.image, iv, ix, iy) = ( I_nu[iv] - sigma2[iv] ) * DnsubSquare;
                    MirImg_PIXEL(*glb.StokesQ, iv, ix, iy) = Q_nu[iv] * DnsubSquare;
                    MirImg_PIXEL(*glb.StokesU, iv, ix, iy) = U_nu[iv] * DnsubSquare;
                    MirImg_PIXEL(*glb.sigma2, iv, ix, iy) = sigma2[iv] * DnsubSquare;
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
    
    SpImgTrac_InitRay(root, &dx, &dy, &ray, &glb.tel_parms);
    
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
                    
                    /* Calculate intensity contributed by this step */
                    //debug
                    //I_nu[iv] += S_nu * (1.0 - exp(-dtau_nu)) * exp(-tau_nu[iv]);
                    I_nu[iv] += S_nu * (1.0 - exp(-dtau_nu)) * exp(-tau_nu[iv]);
                    
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
    
    IntensityBC( side, I_nu, tau_nu, &ray);
    
    return;
}

/*----------------------------------------------------------------------------*/

static void RadiativeXferOverlap(double dx, double dy, double *I_nu, double *tau_nu, size_t tid)
{
    GeRay ray;
    double t;
    size_t side;
    Zone *root = glb.model.grid;
    
    SpImgTrac_InitRay(root, &dx, &dy, &ray, &glb.tel_parms);
    
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
            //              Deb_PRINT("checkpoint: GeRay_TraverseVoxel\n");
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
                        size_t i = glb.line;
                        for(size_t j = 0; j < NRAD; j++) {
                            if(OVERLAP(i,j)){
                                double tempj_nu, tempk_nu;
                                if(i==j){
                                    /* Calculate molecular line emission and absorption coefficients */
                                    SpPhys_GetMoljk(pp, j, vfac, &tempj_nu, &tempk_nu);
                                }
                                else{
                                    /* Calculate velocity line profile factor */
                                    double vfac2 = pp->has_tracer ? SpPhys_GetVfac(&ray, t, dv-RELVEL(i,j), zp, 0) : 0.0;
                                    /* Calculate molecular line emission and absorption coefficients */
                                    SpPhys_GetMoljk(pp, j, vfac2, &tempj_nu, &tempk_nu);
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
                    I_nu[iv] += S_nu * (1.0 - exp(-dtau_nu))  * exp(-tau_nu[iv]);
                    
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
    
    IntensityBC( side, I_nu, tau_nu, &ray);
    
    return;
}

/*----------------------------------------------------------------------------*/

static void RadiativeXferZeeman(double dx, double dy, double *V_nu, double *tau_nu, size_t tid)
{
    GeRay ray;
    double t;
    Zone *root = glb.model.grid;
    SpImgTrac_InitRay(root, &dx, &dy, &ray, &glb.tel_parms);
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
                
                /* Z is Zeeman splitting factor (Hz/G)
                 *                                   see the reference for CN 
                 *                                   http://adsabs.harvard.edu/abs/1996ApJ...456..217C */
                //static const double Z = 2.18 * 1e6;
                
                #define Z glb.model.parms.z
                double dnu = Z * B_Mag;
                #undef Z
                
                double deltav = dnu * PHYS_CONST_MKS_LIGHTC / glb.freq;
                
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
                        SpPhys_GetMoljk(pp, glb.line, vfac, &tempj_nu, &tempk_nu);
                        j_nu = 0.5 * tempj_nu;
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
                    V_nu[iv] += temp * costheta;
                    
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
    
    SpImgTrac_InitRay(root, &dx, &dy, &ray, &glb.tel_parms);
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
                 *                                   Do calculations on all channels at this pixel. Try to minimize the
                 *                                   amount of operations in this loop, since everything here is repeated
                 *                                   for ALL channels, and can significantly increase computation time.
                 */
                
                /* 
                 *                                   Reference : ARTIST(DustPol) -- http://arxiv.org/pdf/1204.6668.pdf
                 */
                
                GeVec3_d B = SpPhys_GetBfac(&ray, t, zp, 0);
                double nproduct = GeVec3_DotProd(&n,&B);
                double eproduct = GeVec3_DotProd(&e,&B);
                double zproduct = GeVec3_DotProd(&z,&B);
                
                // psi is the angle between the projected B-field on p.o.s. and the north of the image 
                double B_Mag = GeVec3_Mag(&B);
                double psi = atan2( -eproduct, nproduct); 
                // gamma is the angle bettwen B-field an the plane of sky
                double cosgammasquare = 1.0 - zproduct * zproduct / GeVec3_Mag(&B);
                
                double alpha = pp->alpha;
                
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
                    double contribution = S_nu * (1.0 - exp(-dtau_nu)) * exp(-tau_nu[iv]);
                    
                    I_nu[iv] += contribution;
                    if (B_Mag == 0.){
                        // do nothing
                        // preventing undefined psi
                    }
                    else{
                        Q_nu[iv] += alpha * contribution * cos(2.0 * psi) * cosgammasquare;
                        U_nu[iv] += alpha * contribution * sin(2.0 * psi) * cosgammasquare;
                        static const double d23 = 2. / 3. ;
                        sigma2[iv] += 0.5 * alpha * contribution * (cosgammasquare - d23);
                    }
                    
                    
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
    
    IntensityBC( side, I_nu, tau_nu, &ray);
    
    
    return;
}


/*---------------------------------------------------------------------------- */
static void IntensityBC( size_t side, double *I_nu, double *tau_nu, GeRay *ray){
    
    int geom = glb.model.grid->voxel.geom;
    
    #define ADD_BC( INTENSITY ) \
    for(size_t iv = 0; iv < glb.v.n; iv++)\
        I_nu[iv] += (INTENSITY) * exp(-tau_nu[iv]);
    
    if ( (geom == GEOM_SPH1D || geom == GEOM_SPH3D) && ( side == 0 ) ){
        /* Ray has been reached inner boundary, to give inner B.C. T_in */
        ADD_BC(glb.I_in)
    }
    else if ( (geom == GEOM_SPH1D || geom == GEOM_SPH3D) && (side >= 2)) {
        Deb_ASSERT(0);
    }
    /* Ray escaped cloud, add CMB to all lines */
    else{
        /* Add CMB to all channels -- this is done even if the ray misses the source */
        /* 
         *            calculate the difference of the angle (alpha) between the ray and the source
         *            cos alpha = unit vec d1  cross  unit vec d2                    (1)
         *            where unit ve =  ( sin theta cos phi, 
         *                            sin theta sin phi, 
         *                            cos theta           )                       (2)
         *            cos alpha =  sin theta1 cos phi1 sin theta2 cosphi2 
         *                    + sin theta1 sin phi1 sin theta2 sin phi2 
         *                    + cos theta1 ocs theta2
         *                    =  sin theta1 sin theta2 cos(phi1-phi2)
         *                    + cos theta1 cos theta2
         *                    =  (-cos(theta1+theta2)+cos(theta1-theta2))/2 cos(phi1-phi2)
         *                    + ( cos(theta1+theta2)+cos(theta1-theta2))/2        (3)
         */            
        SpPhysParm * parms = &glb.model.parms;
        SourceData *target = NULL; 
        for ( int source_id = 0; source_id < parms->Outer_Source; source_id++){
            SourceData *candidate = &parms->source[source_id];
            
            GeVec3_d source_vec = GeVec3_Sub( &ray->e, &candidate->pt_cart);
            GeVec3_d source_d = GeVec3_Normalize( &source_vec) ;
            double cos_alpha = GeVec3_DotProd( &ray->d, &source_d);
            double alpha = acos(cos_alpha);
            
            #if 0
            //printf("%e %e %e\n", candidate->pt_sph.x[0], candidate->pt_sph.x[1], candidate->pt_sph.x[2]);
            //printf("%e %e %e\n", candidate->pt_cart.x[0], candidate->pt_cart.x[1], candidate->pt_cart.x[2]);
            printf("%e %e\n", alpha, candidate->beta);
            exit(0);
            #endif         
            // the ray is inside the view of angle of the source
            // alpha < beta (angle of view of the radius of the source)
            // cos alpha > cos beta
            if ( alpha < candidate->beta ){
                // printf("%e %e\n",alpha,candidate->beta);
                if ( !target || (target && candidate->distance < target->distance) )
                    target = candidate;
            }
        }
        if (target){
            // printf("OK %e\n",target->intensity[0]);
            ADD_BC(target->intensity[0])
        }
        else{
            /* Add CMB to all channels -- this is done even if the ray misses the source */
            ADD_BC(glb.I_cmb)
        }
    }
    #undef ADD_BC
    
}


/*----------------------------------------------------------------------------*/
/* initialize line-of-sight coordinate */
static void InitLOSCoord( double *dx, double *dy, GeRay *ray, GeVec3_d *z, GeVec3_d *n, GeVec3_d *e)
{
    double phi = M_PI - *dx;   
    double theta = 0.5 * M_PI - *dy;
    
    SpTelsim * tel_parms = &glb.tel_parms;
    
    /* line of sight coordinate */
    GeVec3_X(*z,0) = GeRay_D(*ray, 0);
    GeVec3_X(*z,1) = GeRay_D(*ray, 1);
    GeVec3_X(*z,2) = GeRay_D(*ray, 2);
    
    GeVec3_X(*n,0) = -cos(theta)*cos(phi);
    GeVec3_X(*n,1) = -cos(theta)*sin(phi);
    GeVec3_X(*n,2) =  sin(theta);
    *n = GeVec3_Rotate_x(n, -tel_parms->rotate[0]);
    *n = GeVec3_Rotate_y(n, -tel_parms->rotate[1]);
    *n = GeVec3_Rotate_z(n, -tel_parms->rotate[2]);
    
    GeVec3_X(*e,0) = sin(phi);
    GeVec3_X(*e,1) = -cos(phi);
    GeVec3_X(*e,2) = 0.0;
    *e = GeVec3_Rotate_x(e, -tel_parms->rotate[0]);
    *e = GeVec3_Rotate_y(e, -tel_parms->rotate[1]);
    *e = GeVec3_Rotate_z(e, -tel_parms->rotate[2]);
    
    return;
}



