#include "sparx.h"
#include "vtk-wrapper.h"

/* Global parameter struct */
static struct glb {
    DatINode *task;


    SpModel model;
    char *FileName;
} glb;


#define RELVEL(i,j)\
glb.model.parms.mol->OL[NRAD*(i)+(j)]->RelativeVel
#define OVERLAP(i,j)\
glb.model.parms.mol->OL[NRAD*(i)+(j)]->overlap

#define NRAD\
glb.model.parms.mol->nrad
#define FREQ(i)\
glb.model.parms.mol->rad[(i)]->freq


/*----------------------------------------------------------------------------*/

int SpTask_Pops2ASCII(void)
{
    Mem_BZERO(&glb);
    int sts = 0;
    
    // 1. GET PARAMETERS FROM PYTHON
    PyObject *o;
    /*    1-1 those are the general parameters */
    /* out (mandatory) */
    if(!sts){
        sts = SpPy_GetInput_PyObj("out", &o);
        glb.FileName = Sp_PYSTR(o);
        SpPy_XDECREF(o);
    }
    

    /*    1-2 read the source model */
    /* source */
    if(!sts) {
        int popsold = 1;
        sts = SpPy_GetInput_model("source","source", &glb.model, &popsold, TASK_POPS2ASCII);
    }

    /* 2. I/O : OUTPUT */
    if(!sts){
        Zone *root = glb.model.grid, *zp;
        
        FILE * fp = fopen( glb.FileName, "w");
        for(zp = Zone_GetMinLeaf(root); zp; zp = Zone_AscendTree(zp)) {
            SpPhys *pp = zp->data;
            if((pp->n_H2 > 1e-200) && !zp->children) {
                /* This is a non-empty leaf zone */
                if(pp->X_mol > 0) {
                    /* This zone contains tracer molecules */
                    GeVec3_d SphPos = GeVec3_Geom2Sph( glb.model.parms.geom, &zp->voxel.cen);
                    
                    fprintf(fp,"%g %g %g", SphPos.x[0], SphPos.x[1], SphPos.x[2]);
                    for (size_t l=0; l < glb.model.parms.mol->nlev; l++)
                        fprintf(fp," %g", pp->pops_preserve[l]);
                    fprintf(fp," \n");
                }
            }
        }
        fclose(fp);
        
    }
    
/* 3. Cleanup */

SpModel_Cleanup(glb.model);

return sts;
}







