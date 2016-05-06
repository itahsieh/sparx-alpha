#include "vtk-wrapper.h"
#include "geometry.h"
#include "memory.h"
#include "debug.h"
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

