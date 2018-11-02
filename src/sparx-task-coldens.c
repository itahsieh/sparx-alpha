#include "sparx.h"
#include "task.h"
#include "unit.h"

/* Miriad support */
#define Sp_MIRSUPPORT MIRSUPPORT

/* Global parameter struct */
static struct glb {
    DatINode *task;
    int tracer;
    
    DatINode *unit;
    MirImg_Axis x, y, v;
    MirFile *imgf;
    MirImg *image;
    double I_norm, I_cmb, I_in;
    SpModel model;
    SpTelsim tel_parms;
} glb;

/* Subroutine prototypes */
static int InitModel(void);
static void *CalcImageThreadColdens(void *tid_p);
static void ColumnDensityTracer(double dx, double dy, double *CD);




/*----------------------------------------------------------------------------*/

int SpTask_ColDens(void)
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
            case TASK_COLDENS:
                if(!sts) sts = SpPy_GetInput_bool("tracer", &glb.tracer);
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
        sts = SpPy_GetInput_model("source","source", &glb.model, &popsold, task_id);
    }
    
    
    /* 2. Initialize model */
    if(!sts) sts = InitModel();
    
    /* 3. Synthesize image */
    if(!sts) {
        /* Allocate image */
        glb.image = MirImg_Alloc(glb.x, glb.y, glb.v);

        /* Calculate image */
        switch(glb.task->idx){
            case TASK_COLDENS:
                return SpUtil_Threads2(Sp_NTHREAD, CalcImageThreadColdens);
            default:
                /* Shouldn't reach here */
                Deb_ASSERT(0);
        }
    }
    
    /* 4. I/O : OUTPUT */
    if(!sts){
        double scale_factor = 1.0;
        int stokes;
        switch(glb.task->idx){
            // output column density image
            case TASK_COLDENS:
                switch(glb.unit->idx){
                    case UNIT_CGS:
                        scale_factor = 1e-4;
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
   
            default:
                Deb_ASSERT(0);
        }
        FITSoutput( glb.imgf, glb.image, NULL, NULL, glb.unit->name, scale_factor, stokes);
        Sp_PRINT("Wrote FITS image to `%s'\n", glb.imgf->name);
        
    }
    
    /* 5. Cleanup */
    #if Sp_MIRSUPPORT
    /* Miriad images must always be closed! */
    if(glb.imgf)
        MirXY_Close(glb.imgf);
    #endif
    if(glb.image)
        MirImg_Free(glb.image);

    if(glb.tel_parms.subres)
        free(glb.tel_parms.subres);
    
    SpModel_Cleanup(glb.model);
    
    return sts;
}

/*----------------------------------------------------------------------------*/

static int InitModel(void)
{
    Zone *root = glb.model.grid, *zp;
    
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

    return 0;
}






/*----------------------------------------------------------------------------*/

static void *CalcImageThreadColdens(void *tid_p)
{
    /* Column density image tracer
     *           CD is the column density at pixel (ix, iy)
     *           glb.v.n = 1                                     */
    
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
                size_t nsub = SpImgTrac_Init_nsub( ix, iy, &glb.tel_parms, &glb.x, &glb.y);
                Mem_BZERO2(CD,1);
                /* Loop through sub-resolution positions */
                for(size_t isub = 0; isub < nsub; isub++) {
                    for(size_t jsub = 0; jsub < nsub; jsub++) {
                        double dx, dy;
                        /* Calculate sub-pixel position */
                        SpImgTrac_InitSubPixel( &dx, &dy, ix, iy, isub, jsub, nsub, &glb.x, &glb.y);
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
static void ColumnDensityTracer(double dx, double dy, double *CD)
{
    GeRay ray;
    double t;
    size_t side;
    Zone *root = glb.model.grid;
    
    SpImgTrac_InitRay(root, &dx, &dy, &ray, &glb.tel_parms);
    
    //      GeVec3_d z;
    // 	z.x[0]=-GeRay_D(ray, 0);
    // 	z.x[1]=-GeRay_D(ray, 1);
    // 	z.x[2]=-GeRay_D(ray, 2);
    
    Mem_BZERO(CD);
    
    /* Shoot ray at model and see what happens! */
    if ( GeRay_IntersectVoxel(&ray, &root->voxel, &t, &side) ) {
        
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
                double column_density_per_cell = pp->n_H2 * t * Sp_LENFAC;
                if (glb.tracer)
                    column_density_per_cell *= pp->X_mol;
                *CD += column_density_per_cell;
            }
            /* Calculate next position */
            ray = GeRay_Inc(&ray, t);
            
            /* Get next zone to traverse to */
            zp = Zone_GetNext(zp, &side, &ray);
        }
    }
    
    return;
}





