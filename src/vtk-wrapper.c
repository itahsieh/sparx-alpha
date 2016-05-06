#include "vtk-wrapper.h"
#include "geometry.h"
#include "zone.h"
#include "memory.h"
#include "debug.h"
#include "sparx.h"

#include <stdio.h>



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


/*----------------------------------------------------------------------------*/

void Vtk_Mem_CALL(GEOM_TYPE geom, VtkData * visual, size_t nvelo)
{
        
        
        size_t nelement;
        
        switch(geom){
                case GEOM_SPH1D:
                case GEOM_SPH3D:
                {        
                        // Dimension of the visualized resolution
                        size_t nr = visual->sph3d->nr;
                        size_t nt = visual->sph3d->nt;
                        size_t np = visual->sph3d->np;
                        nelement =  nr * np * nt;
                        
                        // declare the memory
                        double * radius = Mem_CALLOC( nr+1, radius);
                        double * theta = Mem_CALLOC( nt+1, theta);
                        double * phi = Mem_CALLOC( np+1, phi);
 
                        // link to the global pointer
                        visual->sph3d->radius = radius;
                        visual->sph3d->theta = theta;
                        visual->sph3d->phi = phi;
                }
                        break;
                case GEOM_CYL3D:
                {
                        // Dimension of the visualized resolution
                        size_t nr = visual->cyl3d->nr;
                        size_t np = visual->cyl3d->np;
                        size_t nz = visual->cyl3d->nz;
                        nelement =  nr * np * nz;
                        
                        // declare the memory
                        double * Rc = Mem_CALLOC( nr+1, Rc);
                        double * phi = Mem_CALLOC( np+1, phi);
                        double * Z = Mem_CALLOC( nz+1, Z);

                        
                        // link to the global pointer
                        visual->cyl3d->Rc   = Rc;
                        visual->cyl3d->phi  = phi;
                        visual->cyl3d->Z    = Z;
                }
                        break;
                default:
                        Deb_ASSERT(0);
        }
        
        double ** contrib       = Mem_CALLOC( nelement, contrib);
        double * contrib_dust   = Mem_CALLOC( nelement, contrib_dust);
        double ** tau           = Mem_CALLOC( nelement, tau);
        double ** tau_dev       = Mem_CALLOC( nelement, tau_dev);
        
        for (size_t idx = 0; idx < nelement; idx++){
                contrib[idx] = Mem_CALLOC(nvelo, contrib[idx]);
                tau[idx]     = Mem_CALLOC(nvelo, tau[idx]);
                tau_dev[idx] = Mem_CALLOC(nvelo, tau_dev[idx]);
        }
        visual->contrib_dust = contrib_dust;
        visual->contrib = contrib;
        visual->tau = tau;
        visual->tau_dev = tau_dev;

        return;
}

/*----------------------------------------------------------------------------*/

void Vtk_Mem_FREE(GEOM_TYPE geom, VtkData * visual)
{
        size_t nelement;
        
        switch(geom){
                case GEOM_SPH1D:
                case GEOM_SPH3D:
                {
                        size_t nr = visual->sph3d->nr;
                        size_t nt = visual->sph3d->nt;
                        size_t np = visual->sph3d->np;
                        nelement =  nr * nt * np;

                        free(visual->sph3d->radius);
                        free(visual->sph3d->theta);
                        free(visual->sph3d->phi);
                        free(visual->sph3d);
                }
                        break;
                case GEOM_CYL3D:
                {
                        size_t nr = visual->cyl3d->nr;
                        size_t np = visual->cyl3d->np;
                        size_t nz = visual->cyl3d->nz;
                        nelement =  nr * np * nz;

                        free(visual->cyl3d->Rc);
                        free(visual->cyl3d->phi);
                        free(visual->cyl3d->Z);
                        free(visual->cyl3d);
                }
                        break;
                default:
                        Deb_ASSERT(0);
        }
        
        double **contrib        = visual->contrib;
        double *contrib_dust    = visual->contrib_dust;
        double **tau            = visual->tau;
        double **tau_dev        = visual->tau_dev;
        
        for (size_t idx = 0; idx < nelement; idx++){
                free(contrib[idx]);
                free(    tau[idx]);
                free(tau_dev[idx]);
        }
        free(contrib);
        free(contrib_dust);
        free(tau);
        free(tau_dev);
        
        return;
}

/*----------------------------------------------------------------------------*/

GeVec3_d Vtk_Index2GeomPos(size_t i, size_t j, size_t k, GEOM_TYPE geom, VtkData * visual)
{
        GeVec3_d GeomPos;
        
        switch (geom){
            case GEOM_SPH1D:
                {
                // link to the global pointer
                double * radius = visual->sph3d->radius;
                double * theta = visual->sph3d->theta;
                double * phi = visual->sph3d->phi;

                GeomPos.x[0] = radius[i];
                GeomPos.x[1] = theta[j];
                GeomPos.x[2] = phi[k];
                }
                break;
            case GEOM_SPH3D:
                {
                // link to the global pointer
                double * radius = visual->sph3d->radius;
                double * theta = visual->sph3d->theta;
                double * phi = visual->sph3d->phi;
                
                GeomPos.x[0] = radius[i];
                GeomPos.x[1] = theta[j];
                GeomPos.x[2] = phi[k];
                }
                break;
            case GEOM_REC3D:
                {
                // link to the global pointer
                double * x = visual->rec3d->x;
                double * y = visual->rec3d->y;
                double * z = visual->rec3d->z;
                
                GeomPos.x[0] = x[i];
                GeomPos.x[1] = y[j];
                GeomPos.x[2] = z[k];
                }
                break;
            case GEOM_CYL3D:
                {
                // link to the global pointer
                double * Rc     = visual->cyl3d->Rc;
                double * phi    = visual->cyl3d->phi;
                double * Z      = visual->cyl3d->Z;
                
                GeomPos.x[0] = Rc[i];
                GeomPos.x[1] = phi[j];
                GeomPos.x[2] = Z[k];
                }
                break;
            default:
                /* Should not happen */
                Deb_ASSERT(0);
        }
        
        return GeomPos;
}

/*----------------------------------------------------------------------------*/

GeVec3_d Vtk_Geom2CartPos( GEOM_TYPE geom, GeVec3_d * GeomPos)
{
        GeVec3_d CartPos;
        
        switch (geom){
            case GEOM_SPH1D:
                {
                double radius   = GeomPos->x[0];
                double theta    = GeomPos->x[1];
                double phi      = GeomPos->x[2];

                CartPos.x[0] = radius * sin(theta) * cos(phi);
                CartPos.x[1] = radius * sin(theta) * sin(phi);
                CartPos.x[2] = radius * cos(theta);
                }
                break;
            case GEOM_SPH3D:
                {
                double radius   = GeomPos->x[0];
                double theta    = GeomPos->x[1];
                double phi      = GeomPos->x[2];

                CartPos.x[0] = radius * sin(theta) * cos(phi);
                CartPos.x[1] = radius * sin(theta) * sin(phi);
                CartPos.x[2] = radius * cos(theta);
                }
                break;
            case GEOM_REC3D:
                CartPos = *GeomPos;
                break;
            case GEOM_CYL3D:
                {
                // link to the global pointer
                double Rc  = GeomPos->x[0];
                double phi = GeomPos->x[1];
                double Z   = GeomPos->x[2];
                
                CartPos.x[0] = Rc * cos(phi);
                CartPos.x[1] = Rc * sin(phi);
                CartPos.x[2] = Z;
                }
                break;
            default:
                /* Should not happen */
                Deb_ASSERT(0);
        }
        
        return CartPos;
}

/*----------------------------------------------------------------------------*/

void Vtk_Output(size_t n1, size_t n2, size_t n3, VtkData * visual, Zone * root, size_t line, size_t nvelo)
{
        GEOM_TYPE geom = root->voxel.geom;
        size_t nelement = n1 * n2 * n3;
        size_t npoint = (n3+1) * (n2+1) * (n1+1);
        
        FILE *fp;
        char filename[32];
        
        // open VTK file
        sprintf(filename,"vis.vtk");
        fp=fopen(filename,"w");
        
        // write the header
        fprintf(fp,"# vtk DataFile Version 3.0\n");
        fprintf(fp,"%s\n", "POSTPROCESSING VISUALIZATION");
        fprintf(fp,"ASCII\n");
        
        // define the type of the gridding
        fprintf(fp,"DATASET STRUCTURED_GRID\n"); 
        fprintf(fp,"DIMENSIONS %zu %zu %zu\n", n3+1, n2+1, n1+1);
        fprintf(fp,"POINTS %zu float\n", npoint );
        for( size_t i = 0; i < n1 + 1; i++)
          for( size_t j = 0; j < n2 + 1; j++)
            for( size_t k = 0; k < n3 + 1; k++){
                GeVec3_d GeomPos = Vtk_Index2GeomPos(i, j, k, geom, visual);
                GeVec3_d CartPos = Vtk_Geom2CartPos(geom, &GeomPos);
                fprintf(fp,"%E %E %E\n", CartPos.x[0],  CartPos.x[1],  CartPos.x[2]);
            }
        
        // write the artributes
        fprintf(fp,"CELL_DATA %zu\n", nelement );
        // H2 number density
        fprintf(fp,"SCALARS H2_Number_Density float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < n1; i++) 
         for( size_t j = 0; j < n2; j++){ 
          for( size_t k = 0; k < n3; k++){
                size_t izone = ZoneIndex( geom, i, j, k, root);
                SpPhys *pp = root->children[izone]->data;
                fprintf(fp,"%E ", pp->n_H2);
           }
           fprintf(fp,"\n");
         }
        // molecular number density
        fprintf(fp,"SCALARS Molecular_Number_Density float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < n1; i++)
         for( size_t j = 0; j < n2; j++)
          for( size_t k = 0; k < n3; k++){
                size_t izone = ZoneIndex( geom, i, j, k, root);
                SpPhys *pp = root->children[izone]->data;
                fprintf(fp,"%E ", pp->n_H2 * pp->X_mol);
        }fprintf(fp,"\n");
        
        // kinetic temperature
        fprintf(fp,"SCALARS Temperature float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < n1; i++)
         for( size_t j = 0; j < n2; j++){
          for( size_t k = 0; k < n3; k++){
                size_t izone = ZoneIndex( geom, i, j, k, root);
                SpPhys *pp = root->children[izone]->data;
                fprintf(fp,"%E ", pp->T_k);
           }
           fprintf(fp,"\n");
         }
        
        // exitation temperature
        fprintf(fp,"SCALARS Tex float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for( size_t i = 0; i < n1; i++)
         for( size_t j = 0; j < n2; j++){
          for( size_t k = 0; k < n3; k++){
                size_t izone = ZoneIndex( geom, i, j, k, root);
                SpPhys *pp = root->children[izone]->data;
                if(pp->X_mol == 0.0){
                        fprintf(fp,"%E ", 0.0);
                }
                else{
                        MolTrRad *trans = pp->mol->rad[line];
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
        
        // the global pointer
        double ** contrib       = visual->contrib;
        double * contrib_dust   = visual->contrib_dust;
        double ** tau           = visual->tau;
        double ** tau_dev       = visual->tau_dev;
        
        // write the dust contribution of the cells
        fprintf(fp,"SCALARS DUST_CONTRIBUTION float 1\n");
        fprintf(fp,"LOOKUP_TABLE default\n");
        for (size_t idx = 0; idx < nelement; idx++)
                fprintf(fp,"%E ", contrib_dust[idx]);
        fprintf(fp,"\n");
        
        
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
        
        
        // for seperate contribution channel as a file
        for( size_t l = 0; l < nvelo; l++){
                // open VTK file
                sprintf(filename,"contribution_%04zu.vtk", l);
                fp=fopen(filename,"w");
                
                // write the header
                fprintf(fp,"# vtk DataFile Version 3.0\n");
                fprintf(fp,"%s\n", "POSTPROCESSING VISUALIZATION");
                fprintf(fp,"ASCII\n");
                
                // define the type of the gridding
                fprintf(fp,"DATASET STRUCTURED_GRID\n"); 
                fprintf(fp,"DIMENSIONS %zu %zu %zu\n", n3+1, n2+1, n1+1);
                fprintf(fp,"POINTS %zu float\n", npoint );
                for( size_t i = 0; i < n1 + 1; i++)
                  for( size_t j = 0; j < n2 + 1; j++)
                    for( size_t k = 0; k < n3 + 1; k++){
                        GeVec3_d GeomPos = Vtk_Index2GeomPos(i, j, k, geom, visual);
                        GeVec3_d CartPos = Vtk_Geom2CartPos(geom, &GeomPos);
                        fprintf(fp,"%E %E %E\n", CartPos.x[0],  CartPos.x[1],  CartPos.x[2]);
                }
                
                // write the artributes
                fprintf(fp,"CELL_DATA %zu\n", nelement );
                // write the contribution of the cells
                fprintf(fp,"SCALARS DUST_CONTRIBUTION float 1\n");
                fprintf(fp,"LOOKUP_TABLE default\n");
                for (size_t idx = 0; idx < nelement; idx++)
                        fprintf(fp,"%E ", contrib[idx][l]);
                fprintf(fp,"\n");
                
                fclose(fp);
                printf("wrote %s\n",filename);
        }
        
        return;
}
