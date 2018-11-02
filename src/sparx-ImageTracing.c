#include "sparx.h"

/*---------------------------------------------------------------------------- */
/* Initialize the position of the ray and its shooting direction */
void SpImgTrac_InitRay(Zone *root, double *dx, double *dy, GeRay *ray, SpTelsim *tel_parms)
{           
    /* Reset ray */
    Mem_BZERO(ray);
    
    
    double ModelSize = Zone_ZoneSize(root);
    double Model2DistanceRatio = ModelSize / tel_parms->dist;
    
    if ( Model2DistanceRatio > 1.0 ){
        // distance is too close
        Sp_PRINT("The observing distance is too close!\n");
        Deb_ASSERT(0);
    }
    else{ 
        if ( Model2DistanceRatio > 1e-4 ){
            // stereopsis projection
            /* Init ray position to <dist, 0, 0> */
            GeRay_E(*ray, 0) = tel_parms->dist;
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
         *                if the ratio of domain size and observing distance is less than 10^-4,
         *                then switch to parallel projection mode to avoid precision lost
         *                particularly in cylindrical and spherical coordinate
         *                the ray shooting solves the second order equation with the coefficient 
         *                c = Ex^2 - r^2
         */
        else{
            GeRay_E(*ray, 0) = ModelSize;
            GeRay_E(*ray, 1) = (*dx) * tel_parms->dist;
            GeRay_E(*ray, 2) = (*dy) * tel_parms->dist;
            
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
    *ray = GeRay_Rotate(ray, 0, -tel_parms->rotate[0]);
    *ray = GeRay_Rotate(ray, 1, -tel_parms->rotate[1]);
    *ray = GeRay_Rotate(ray, 2, -tel_parms->rotate[2]);
    
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
/* determine nsub if the pixel is in the sub-sampling position */
size_t SpImgTrac_Init_nsub(size_t ix, size_t iy, SpTelsim *tel_parms, MirImg_Axis *x, MirImg_Axis *y)
{
    /* Determine angular offsets from the pointing center */
    double dx = ((int)ix - (int)x->crpix) * x->delt;
    double dy = ((int)iy - (int)y->crpix) * y->delt;
    
    /* nsub is by default 1 */
    size_t nsub = 1;
    
    /* Check if position is within any of the subres boxes */
    for(size_t ibox = 0; ibox < tel_parms->nsubres; ibox++) {
        if((dx >= tel_parms->subres[ibox].blc_x) && (dx <= tel_parms->subres[ibox].trc_x) &&
            (dy >= tel_parms->subres[ibox].blc_y) && (dy <= tel_parms->subres[ibox].trc_y)) {
            /* Position within subres box, set nsub to subres[ibox].nsub */
            nsub = tel_parms->subres[ibox].nsub;
            break;
        }
    }
    
    return nsub;
}

/*----------------------------------------------------------------------------*/
void SpImgTrac_InitSubPixel( double *dx, double *dy, 
                          size_t ix, size_t iy, 
                          size_t isub,size_t jsub, 
                          size_t nsub,
                          MirImg_Axis *x, MirImg_Axis *y
                                )
{
    /* Determine sub-resolution angular offsets from the pointing center */
    double subix = (double)ix + (2.*(double)isub-(double)nsub+1.) / ((double)(2*nsub));
    double subiy = (double)iy + (2.*(double)jsub-(double)nsub+1.) / ((double)(2*nsub));
    *dx = (subix - (double)x->crpix) * x->delt;
    *dy = (subiy - (double)y->crpix) * y->delt;
    
    return;
}
