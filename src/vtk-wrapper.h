#ifndef __VTK_WRAPPER_H__
#define __VTK_WRAPPER_H__


#include "geometry.h"
#include "zone.h"

#include <math.h>
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




void Vtk_Mem_CALL(GEOM_TYPE geom, VtkData * visual, size_t nvelo);
void Vtk_Mem_FREE(GEOM_TYPE geom, VtkData * visual);

GeVec3_d Vtk_Index2GeomPos(size_t i, size_t j, size_t k, GEOM_TYPE geom, VtkData * visual);
GeVec3_d Vtk_Geom2CartPos( GEOM_TYPE geom, GeVec3_d * GeomPos);

void Vtk_Output(size_t n1, size_t n2, size_t n3, VtkData * visual, Zone * root, size_t line, size_t nvelo);

#endif