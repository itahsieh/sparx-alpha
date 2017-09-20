#include "vtk-wrapper.h"
#include "geometry.h"
#include "zone.h"
#include "memory.h"
#include "debug.h"
#include "sparx.h"

#include <stdio.h>




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
        
        double * contrib_dust   = Mem_CALLOC( nelement, contrib_dust);
        double * tau_dust       = Mem_CALLOC( nelement, tau_dust);
        
        double ** contrib       = Mem_CALLOC( nelement, contrib);
        double ** tau           = Mem_CALLOC( nelement, tau);
        double ** tau_dev       = Mem_CALLOC( nelement, tau_dev);
        
        for (size_t idx = 0; idx < nelement; idx++){
                contrib[idx] = Mem_CALLOC(nvelo, contrib[idx]);
                tau[idx]     = Mem_CALLOC(nvelo, tau[idx]);
                tau_dev[idx] = Mem_CALLOC(nvelo, tau_dev[idx]);
        }
        visual->contrib_dust = contrib_dust;
        visual->tau_dust = tau_dust;
        
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

void Vtk_InitializeGrid(size_t nvelo, Zone * root, VtkData *visual, GEOM_TYPE geom)
{
        size_t n1, n2, n3;
        switch (geom){
            case GEOM_SPH1D:
                visual->sph3d = Mem_CALLOC( 1, visual->sph3d);
                    
                // Dimension of the visualized resolution
                n1 = visual->sph3d->nr = root->nchildren;
                n2 = visual->sph3d->nt = 45;
                n3 = visual->sph3d->np = 90;
                
                // initialize memory
                Vtk_Mem_CALL(geom, visual, nvelo);
                {
                // link to the global pointer
                double * radius = visual->sph3d->radius;
                double * theta = visual->sph3d->theta;
                double * phi = visual->sph3d->phi;
                
                // construct the expanding SPH3D mesh
                // radius
                radius[0] = root->children[0]->voxel.min.x[0];
                for(size_t i = 1; i < n1 + 1; i++)
                        radius[i] = root->children[i-1]->voxel.max.x[0];
                // theta
                double delta_theta = M_PI / (double) n2;
                theta[0] = 0.;
                for ( size_t j = 1; j < n2 + 1; j++)
                        theta[j] = theta[j-1] + delta_theta;
                //phi
                double delta_phi = 2.0 * M_PI / (double) n3;
                phi[0] = 0.;
                for (size_t k = 1; k < n3 + 1; k ++)
                        phi[k] = phi[k-1] + delta_phi;

                }
                break;
            case GEOM_SPH3D:
                visual->sph3d = Mem_CALLOC( 1, visual->sph3d);
                
                // Dimension of the visualized resolution
                n1 = visual->sph3d->nr = root->naxes.x[0];
                n2 = visual->sph3d->nt = (root->naxes.x[1] == 1) ? 
                        45 : root->naxes.x[1];
                n3 = visual->sph3d->np = (root->naxes.x[2] == 1) ? 
                        90 : root->naxes.x[2];
                        
                // initialize memory
                Vtk_Mem_CALL(geom, visual, nvelo);
                {
                // link to the global pointer
                double * radius = visual->sph3d->radius;
                double * theta = visual->sph3d->theta;
                double * phi = visual->sph3d->phi;
                
                // construct the expanding SPH3D mesh
                // radius
                radius[0] = root->children[0]->voxel.min.x[0];
                for(size_t i = 1; i < n1 + 1; i++)
                        radius[i] = root->children[ (i-1) * root->naxes.x[1] * root->naxes.x[2] ]->voxel.max.x[0];
                // theta
                if (root->naxes.x[1] == 1){
                        double delta_theta = M_PI / (double) n2;
                        theta[0] = 0.;
                        for ( size_t j = 1; j < n2 + 1; j++)
                                theta[j] = theta[j-1] + delta_theta;
                }
                else{
                        theta[0] = root->children[0]->voxel.min.x[1];
                        for ( size_t j = 1; j < n2 + 1; j++)
                                theta[j] = root->children[ (j-1)*root->naxes.x[2] ]->voxel.max.x[1];
                }
                //phi
                if (root->naxes.x[2] == 1){
                        double delta_phi = 2. * M_PI / (double) n3;
                        phi[0] = 0.;
                        for (size_t k = 1; k < n3 + 1; k ++)
                                phi[k] = phi[k-1] + delta_phi;
                }
                else{
                        phi[0] = root->children[0]->voxel.min.x[2];
                        for (size_t k = 1; k < n3 + 1; k ++)
                                phi[k] = root->children[k-1]->voxel.max.x[2];
                }

                }
                break;
            case GEOM_REC3D:
                visual->rec3d = Mem_CALLOC( 1, visual->rec3d);
                
                // Dimension of the visualized resolution
                n1 = visual->rec3d->nx = root->naxes.x[0];
                n2 = visual->rec3d->ny = root->naxes.x[1];
                n3 = visual->rec3d->nz = root->naxes.x[2];
                
                // initialize memory
                Vtk_Mem_CALL(geom, visual, nvelo);
                {
                // link to the global pointer
                double * x = visual->rec3d->x;
                double * y = visual->rec3d->y;
                double * z = visual->rec3d->z;
                
                // construct the expanding CYL3D mesh
                // for x
                x[0] = root->children[0]->voxel.min.x[0];
                for(size_t i = 1; i < n1 + 1; i++)
                        x[i] = root->children[ (i-1)*root->naxes.x[1]*root->naxes.x[2] ]->voxel.max.x[0];
                // for y
                y[0] = root->children[0]->voxel.min.x[0];
                for(size_t j = 1; j < n2 + 1; j++)
                        y[j] = root->children[ (j-1)*root->naxes.x[2] ]->voxel.max.x[0];
                // for z
                z[0] = root->children[0]->voxel.min.x[2];
                for (size_t k = 1; k < n3 + 1; k ++)
                        z[k] = root->children[k-1]->voxel.max.x[2];
                
                }
                break;
            case GEOM_CYL3D:
                visual->cyl3d = Mem_CALLOC( 1, visual->cyl3d);
                // Dimension of the visualized resolution
                n1 = visual->cyl3d->nr = root->naxes.x[0];
                n2 = visual->cyl3d->np = (root->naxes.x[1] == 1) ? 
                        72 : root->naxes.x[1];
                n3 = visual->cyl3d->nz = root->naxes.x[2];

                // initialize memory
                Vtk_Mem_CALL(geom, visual, nvelo);
                {
                // link to the global pointer
                double * Rc     = visual->cyl3d->Rc;
                double * phi    = visual->cyl3d->phi;
                double * Z      = visual->cyl3d->Z;
                
                // construct the expanding CYL3D mesh
                // for Rc
                Rc[0] = root->children[0]->voxel.min.x[0];
                for(size_t i = 1; i < n1 + 1; i++)
                        Rc[i] = root->children[ (i-1) * root->naxes.x[1] * root->naxes.x[2] ]->voxel.max.x[0];
                // for phi
                if (root->naxes.x[1] == 1){
                        double delta_phi = 2.0 * M_PI / (double) n2;
                        phi[0] = 0.;
                        for ( size_t j = 1; j < n2 + 1; j++)
                                phi[j] = phi[j-1] + delta_phi;
                }
                else{
                        phi[0] = root->children[0]->voxel.min.x[1];
                        for ( size_t j = 1; j < n2 + 1; j++)
                                phi[j] = root->children[ (j-1)*root->naxes.x[2] ]->voxel.max.x[1];
                }
                // for Z
                Z[0] = root->children[0]->voxel.min.x[2];
                for (size_t k = 1; k < n3 + 1; k ++)
                        Z[k] = root->children[k-1]->voxel.max.x[2];
                
                }
                break;
            default:
                /* Should not happen */
                Deb_ASSERT(0);
        }
}

/*----------------------------------------------------------------------------*/

void Vtk_Output(VtkFile *vtkfile, VtkData * visual, SpModel *model, size_t line, size_t nvelo, TASK_TYPE task)
{
        // Dimension of the visualized resolution
    Zone * root = model->grid;
    GEOM_TYPE geom = root->voxel.geom;
    
        size_t n1, n2, n3;
        switch (geom){            
            case GEOM_SPH1D:
            case GEOM_SPH3D:
                n1 = visual->sph3d->nr;
                n2 = visual->sph3d->nt;
                n3 = visual->sph3d->np;
                break;
            case GEOM_CYL3D:
                n1 = visual->cyl3d->nr;
                n2 = visual->cyl3d->np;
                n3 = visual->cyl3d->nz;
                break;
            case GEOM_REC3D:
                n1 = visual->rec3d->nx;
                n2 = visual->rec3d->ny;
                n3 = visual->rec3d->nz;
                break;
            default:
                Deb_ASSERT(0);
        }
        
        size_t nelement = n1 * n2 * n3;
        size_t npoint = (n3+1) * (n2+1) * (n1+1);
        
        
        #define WRITE_HEADER_AND_GRID() \
        fprintf(fp,"# vtk DataFile Version 3.0\n"); \
        fprintf(fp,"%s\n", "POSTPROCESSING VISUALIZATION"); \
        fprintf(fp,"ASCII\n"); \
        fprintf(fp,"DATASET STRUCTURED_GRID\n"); \
        fprintf(fp,"DIMENSIONS %zu %zu %zu\n", n3+1, n2+1, n1+1); \
        fprintf(fp,"POINTS %zu float\n", npoint ); \
        for( size_t i = 0; i < n1 + 1; i++) \
            for( size_t j = 0; j < n2 + 1; j++) \
                for( size_t k = 0; k < n3 + 1; k++){ \
                    GeVec3_d GeomPos = Vtk_Index2GeomPos(i, j, k, geom, visual); \
                    GeVec3_d CartPos = Vtk_Geom2CartPos(geom, &GeomPos); \
                    fprintf(fp,"%E %E %E\n", CartPos.x[0], CartPos.x[1], CartPos.x[2]); \
                }

        FILE *fp = vtkfile->fp;

        // open VTK file
        char filename[64];
        sprintf(filename, "%s.vtk", vtkfile->FileName);
        fp = fopen( filename, "w");
        
        WRITE_HEADER_AND_GRID()
        
        // write the artributes
        fprintf(fp,"CELL_DATA %zu\n", nelement );
        

        #define WRITE_SCALAR_MODEL_DATA( SCALAR_NAME, MODEL_DATA ) \
        fprintf(fp,"SCALARS %s float 1\n", SCALAR_NAME); \
        fprintf(fp,"LOOKUP_TABLE default\n"); \
        for( size_t i = 0; i < n1; i++)  \
            for( size_t j = 0; j < n2; j++){  \
                for( size_t k = 0; k < n3; k++){ \
                    size_t izone = ZoneIndex( geom, i, j, k, root); \
                    SpPhys *pp = root->children[izone]->data; \
                    fprintf(fp,"%E ", MODEL_DATA); \
                } \
                fprintf(fp,"\n"); \
            } 
        
        WRITE_SCALAR_MODEL_DATA( "H2_Number_Density", pp->n_H2 ) 
        WRITE_SCALAR_MODEL_DATA( "Molecular_Number_Density", pp->n_H2 * pp->X_mol ) 
        WRITE_SCALAR_MODEL_DATA( "Temperature", pp->T_k ) 

        #undef WRITE_SCALAR_MODEL_DATA( SCALAR_NAME, MODEL_DATA )
        // write the velocity field
        
        #define WRITE_VECTOR_MODEL_DATA( VECTOR_NAME, VECTOR_FUNCTION) \
        fprintf(fp,"VECTORS %s float\n", VECTOR_NAME); \
        for( size_t i = 0; i < n1; i++) for( size_t j = 0; j < n2; j++) for( size_t k = 0; k < n3; k++){ \
            GeVec3_d GeomPos_min = Vtk_Index2GeomPos(i, j, k, geom, visual); \
            GeVec3_d GeomPos_max = Vtk_Index2GeomPos(i+1, j+1, k+1, geom, visual); \
            GeVec3_d GeomPos; \
            GeomPos.x[0] = 0.5 * ( GeomPos_min.x[0] + GeomPos_max.x[0] ); \
            GeomPos.x[1] = 0.5 * ( GeomPos_min.x[1] + GeomPos_max.x[1] ); \
            GeomPos.x[2] = 0.5 * ( GeomPos_min.x[2] + GeomPos_max.x[2] ); \
            GeVec3_d CartPos = Vtk_Geom2CartPos(geom, &GeomPos); \
            size_t izone = ZoneIndex( geom, i, j, k, root); \
            Zone * zp = root->children[izone]; \
            GeVec3_d VecCart = VECTOR_FUNCTION(&CartPos, zp); \
            fprintf(fp,"%E %E %E\n", VecCart.x[0], VecCart.x[1], VecCart.x[2]); \
        }
        WRITE_VECTOR_MODEL_DATA("velocity", SpPhys_GetVgas)
        if( model->parms.polariz ){
            WRITE_VECTOR_MODEL_DATA("b-field", SpPhys_GetBgas)
        }
        #undef WRITE_VECTOR_MODEL_DATA

        
        if(task == TASK_LINECTB || task == TASK_ZEEMANCTB){
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
                          double n_u = pp->pops_preserve[up];
                          double n_l = pp->pops_preserve[lo];
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
        }
        if(task == TASK_LINECTB || task == TASK_CONTCTB || task == TASK_ZEEMANCTB){
            #define WRITE_SCALAR_VISUAL_DATA( SCALAR_NAME, visual_data) \
            fprintf(fp,"SCALARS %s float 1\n", SCALAR_NAME); \
            fprintf(fp,"LOOKUP_TABLE default\n"); \
            double * visual_data   = visual->visual_data; \
            for (size_t idx = 0; idx < nelement; idx++) \
                fprintf(fp,"%E ", visual_data[idx]); \
            fprintf(fp,"\n"); 
            WRITE_SCALAR_VISUAL_DATA( "DUST_CONTRIBUTION", contrib_dust )
            WRITE_SCALAR_VISUAL_DATA( "DUST_TAU", tau_dust )
            #undef WRITE_SCALAR_VISUAL_DATA( SCALAR_NAME, visual_data)
        }
        

        fclose(fp);
        printf("wrote %s\n",filename);
        
        
        
        
        
        if(task == TASK_LINECTB || task == TASK_ZEEMANCTB){
            // for seperate contribution channel as a file
            for( size_t l = 0; l < nvelo; l++){
                // open VTK file
                sprintf( filename, "%s_%04zu.vtk",vtkfile->FileName, l);
                fp=fopen(filename,"w");
                
                WRITE_HEADER_AND_GRID()

                // write the artributes
                fprintf(fp,"CELL_DATA %zu\n", nelement );
                
                
                #define WRITE_SCALAR_VISUAL_DATA( SCALAR_NAME, visual_data) \
                fprintf(fp,"SCALARS %s float 1\n", SCALAR_NAME); \
                fprintf(fp,"LOOKUP_TABLE default\n"); \
                double ** visual_data   = visual->visual_data; \
                for (size_t idx = 0; idx < nelement; idx++) \
                    fprintf(fp,"%E ", visual_data[idx][l]); \
                    fprintf(fp,"\n"); 
                WRITE_SCALAR_VISUAL_DATA( "CONTRIBUTION", contrib)
                WRITE_SCALAR_VISUAL_DATA( "TAU", tau)
                #undef WRITE_SCALAR_VISUAL_DATA( SCALAR_NAME, visual_data)

                fclose(fp);
                printf("wrote %s\n",filename);
            }
        }
        
        #undef WRITE_HEADER_AND_GRID()
 
        return;
}

void FreeVtkFile(VtkFile * vtkfile){
    free(vtkfile->FileName);
    free(vtkfile->fp);
}




