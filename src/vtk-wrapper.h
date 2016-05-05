#include "geometry.h"
#include <stdio.h>


typedef struct VtkData{
        struct {
                size_t nr, nt, np;
                double *radius, *theta, *phi;
        } *sph3d;
        struct {
                size_t nr, np, nz;
                double *Rc, *phi, *Z;
        } *cyl3d;
        struct {
                size_t nx, ny, nz;
                double *x, *y, *z;
        } *rec3d;
        double 
                **contrib, 
                *contrib_dust, 
                **tau, 
                **tau_dev;
} VtkData;


void Mem_CALL_VISUAL(GEOM_TYPE geom, VtkData * visual, size_t nvelo);
void Mem_FREE_VISUAL(GEOM_TYPE geom, VtkData * visual);