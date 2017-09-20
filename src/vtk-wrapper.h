#ifndef __VTK_WRAPPER_H__
#define __VTK_WRAPPER_H__

#include "task.h"
#include "geometry.h"
#include "zone.h"


typedef enum{
    VTK_GRID,
    VTK_DESITY,
    VTK_TEMPERATURE,
    VTK_ABUNDANCE,
    VTK_V_FIELD,
    VTK_B_FIELD,
    VTK_MC_NOISE,
    VTK_LI_NOISE,
    VTK_EXITATION,
    VTK_CONTRIBUTION,
    N_VTK_OUTPUT_TYPE
} VTK_OUTPUT_TYPE;

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
                **tau, 
                **tau_dev,
                *contrib_dust,
                *tau_dust,
                *tau_dust_dev;
} VtkData;

typedef struct VtkFile{
    char *FileName;
    FILE *fp
} VtkFile;

void Vtk_Mem_CALL(GEOM_TYPE geom, VtkData * visual, size_t nvelo);
void Vtk_Mem_FREE(GEOM_TYPE geom, VtkData * visual);

GeVec3_d Vtk_Index2GeomPos(size_t i, size_t j, size_t k, GEOM_TYPE geom, VtkData * visual);
GeVec3_d Vtk_Geom2CartPos( GEOM_TYPE geom, GeVec3_d * GeomPos);

void Vtk_InitializeGrid(size_t nvelo, Zone * root, VtkData *visual, GEOM_TYPE geom);

void Vtk_Output(VtkFile *vtkfile, VtkData * visual, Zone * root, size_t line, size_t nvelo, TASK_TYPE task);

void FreeVtkFile(VtkFile * vtkfile);

#endif
