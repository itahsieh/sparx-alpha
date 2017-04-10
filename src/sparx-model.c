#include "sparx.h"

/*----------------------------------------------------------------------------*/

void SpModel_PrintModel(SpModel model)
/* Display model contents on stdout */
{
    Sp_PRINTF("T_cmb: %g K\n", model.parms.T_cmb);
    
    if(model.parms.mol)
        Mol_Fprintf(stdout, model.parms.mol);
    
    if(model.grid)
        SpZone_FPRINTF(stdout, model.grid);
    
    return;
}

/*----------------------------------------------------------------------------*/

void SpModel_Cleanup(SpModel model)
{
    if(model.grid)
        SpZone_FREE(model.grid);
    
    /* Molecule must be freed AFTER grid, since it
     * requires information on the molecule to free
     * collisional rates and things */
    if(model.parms.mol)
        Mol_Free(model.parms.mol);
    
    if(model.parms.source)
        free(model.parms.source);
    
    return;
}

/*----------------------------------------------------------------------------*/

void SpModel_InitGrid(SpModel *model, int geom, GeVec3_d min, GeVec3_d max, GeVec3_s ndiv)
/* Allocate root zone and set geometry/dimensions of bounding volume */
{
    model->grid = SpZone_ALLOC(&model->parms);
    model->grid->voxel = GeVox_Init(geom,
                                    GeVec3_X(min, 0),
                                    GeVec3_X(min, 1),
                                    GeVec3_X(min, 2),
                                    GeVec3_X(max, 0),
                                    GeVec3_X(max, 1),
                                    GeVec3_X(max, 2));
    SpZone_GROW(model->grid, ndiv, &model->parms);
    
    return;
}

/*----------------------------------------------------------------------------*/

void SpModel_GenModel_UniSphere(SpModel *model, double n_H2, double T_k, double X_mol)
/* Generate a uniform, static, spherical model */
{
    Zone *grid = model->grid, *zp;
    SpPhys *pp;
    GeVec3_d delta = grid->voxel.delta;
    double dx, dy, dz, radius;
    double dim_min = Num_MIN(GeVec3_X(delta, 0), Num_MIN(GeVec3_X(delta, 1), GeVec3_X(delta, 2)));
    
    /* Fill in values */
    for(zp = Zone_GetMinLeaf(grid); zp; zp = Zone_AscendTree(zp)) {
        dx = zp->voxel.cen.x[0] - grid->voxel.cen.x[0];
        dy = zp->voxel.cen.x[1] - grid->voxel.cen.x[1];
        dz = zp->voxel.cen.x[2] - grid->voxel.cen.x[2];
        
        radius = sqrt(dx * dx + dy * dy + dz * dz);
        
        if(radius <= (dim_min / 2.0)) {
            pp = zp->data;
            pp->n_H2 = n_H2;
            pp->T_k = T_k;
            pp->X_mol = X_mol;
            
            pp->v_cen.x[0] = 0;
            pp->v_cen.x[1] = 0;
            pp->v_cen.x[2] = 0;
        }
    }

}

