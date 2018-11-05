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



/*---------------------------------------------------------------------------- */
void SpImgtrac_IntensityBC( size_t side, double *I_nu, double *tau_nu, GeRay *ray,
                            int geom, size_t vn, double I_in, double I_cmb, SpPhysParm * parms){
        
    #define ADD_BC( INTENSITY ) \
    for(size_t iv = 0; iv < vn; iv++)\
        I_nu[iv] += (INTENSITY) * exp(-tau_nu[iv]);
    
    if ( (geom == GEOM_SPH1D || geom == GEOM_SPH3D) && ( side == 0 ) ){
        /* Ray has been reached inner boundary, to give inner B.C. T_in */
        ADD_BC(I_in)
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
            ADD_BC(I_cmb)
        }
    }
    #undef ADD_BC
    
}


/*----------------------------------------------------------------------------*/
/* initialize line-of-sight coordinate */
void SpImgtrac_InitLOSCoord( double *dx, double *dy, GeRay *ray, GeVec3_d *z, GeVec3_d *n, GeVec3_d *e, SpTelsim * tel_parms)
{
    double phi = M_PI - *dx;   
    double theta = 0.5 * M_PI - *dy;
    
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

