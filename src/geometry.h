#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__
#define Sp_MIRSUPPORT MIRSUPPORT

/* The geometry.c/.h interface provides geometrical objects,
 * basic vector operations, and ray tracing operations. */

#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include "data_structs.h"

typedef enum {
	GEOM_SPH1D,
	GEOM_SPH3D,
	GEOM_REC3D,
	GEOM_CYL3D
} GEOM_TYPE;

extern DatINode GEOM_TYPES[];

typedef struct GeVec3_i {
	int x[3];
} GeVec3_i;

typedef struct GeVec3_s {
	size_t x[3];
} GeVec3_s;

typedef struct GeVec3_d {
	double x[3];
} GeVec3_d;

typedef struct GeMat3_d {
	double m[3 * 3];
} GeMat3_d;

typedef struct GeVox {
	/* Remember to update GeVox_INIT and GeVox_Get_subvox() if the
	 * GeVox struct changes!!! */
	int geom;
	GeVec3_d delta, min, max, cen;
} GeVox;

/*
 * Graphical representation of a voxel:
 *
            Z
	    |
	    |
	    |     /
	    4----------7
	   /|   /     /|
	  / |  /     / |
	 /  | /     /  |
	5----------6   |  
--------|---0------|---3---------Y
	|  /|      |  /
	| / |      | /
	|/  |      |/
	1----------2
       /    |
      /     |
     /
    X

*/

typedef struct GeCam {
	GeVec3_d min, max, phi;
} GeCam;

#define GeCam_INIT(phix, phiy, phiz) {\
	GeVec3_INIT(0, 0, 0),\
	GeVec3_INIT(0, 0, 0),\
	GeVec3_INIT(phix, phiy, phiz)\
}

typedef struct GeRay {
	GeVec3_d e, d, tMax, tDelta;
	GeVec3_s pos;
	double t;
} GeRay;

#define RayAW_INIT(ex, ey, ez, dx, dy, dz) {\
	GeVec3_INIT((ex), (ey), (ez)),\
	GeVec3_INIT((dx), (dy), (dz)),\
	GeVec3_INIT(0, 0, 0),\
	0\
}

#define GeMat3_INIT(ii, ij, ik, ji, jj, jk, ki, kj, kk)\
	{{(ii), (ij), (ik), (ji), (jj), (jk), (ki), (kj), (kk)}}

#define GeMat3_X(mat, i, j)\
	(mat).m[(j)+3*(i)]

/* The GeVec3_* init macro */
#define GeVec3_INIT(x1, x2, x3)\
	{{(x1), (x2), (x3)}}

#define GeVec3_PRINT(fp, vec)\
	{fprintf((fp), "<%g, %g, %g>", (double)(vec).x[0], (double)(vec).x[1], (double)(vec).x[2]);}

#define GeVec3_X(vec, axis)\
	((vec).x[(axis)])

#define GeRay_INIT(e1, e2, e3, d1, d2, d3) {\
	GeVec3_INIT((e1), (e2), (e3)), /* origin */\
	GeVec3_INIT((d1), (d2), (d3)), /* direction */\
	GeVec3_INIT(0.0, 0.0, 0.0), /* tMax */\
	GeVec3_INIT(0.0, 0.0, 0.0), /* tDelta */\
	GeVec3_INIT(0, 0, 0), /* pos */\
	0.0 /* t */\
}

#define GeRay_D(ray, i)\
	(GeVec3_X((ray).d, (i)))

#define GeRay_E(ray, i)\
	(GeVec3_X((ray).e, (i)))

/* The voxel init macro */
#define GeVox_INIT(geom, xmin, ymin, zmin, xmax, ymax, zmax)\
	{\
		(geom),\
		GeVec3_INIT(\
			(xmax) - (xmin), /* dx */ \
			(ymax) - (ymin), /* dy */ \
			(zmax) - (zmin)  /* dz */ \
		),\
		GeVec3_INIT((xmin), (ymin), (zmin)), /* Min position */ \
		GeVec3_INIT((xmax), (ymax), (zmax)), /* Max position */ \
		GeVec3_INIT(\
			(xmin) + ((xmax) - (xmin)) * 0.5, /* Center */ \
			(ymin) + ((ymax) - (ymin)) * 0.5,\
			(zmin) + ((zmax) - (zmin)) * 0.5\
		)\
	}

#define GeVox_FPRINTF(fp, voxel)\
	 fprintf(fp, "min: "); GeVec3_PRINT(fp, voxel.min); fprintf(fp, "\n");\
	 fprintf(fp, "max: "); GeVec3_PRINT(fp, voxel.max); fprintf(fp, "\n");\
	 fprintf(fp, "cen: "); GeVec3_PRINT(fp, voxel.cen); fprintf(fp, "\n");}

#define Geom_CodeToName(geom)\
	(Dat_IList_IdxLookup(GEOM_TYPES, (geom))->name)
GeVec3_d GeVec3_d_Init(double x, double y, double z);
GeVec3_s GeVec3_s_Init(size_t x, size_t y, size_t z);
double GeVec3_Mag(const GeVec3_d *a);
double GeVec3_Mag2(const GeVec3_d *a, const GeVec3_d *b);
GeVec3_d GeVec3_Add(const GeVec3_d *a, const GeVec3_d *b);
GeVec3_d GeVec3_Sub(const GeVec3_d *a, const GeVec3_d *b);
GeVec3_d GeVec3_Normalize(const GeVec3_d *a);
GeVec3_d GeVec3_Scale(const GeVec3_d *a, double fac);
GeVec3_d GeVec3_InterpLinear(const GeVec3_d *xa, const GeVec3_d *xb, const GeVec3_d *a, 
	const GeVec3_d *b, const GeVec3_d *pos);
double GeVec3_DotProd(const GeVec3_d *a, const GeVec3_d *b);
GeVec3_d GeVec3_CrossProd(const GeVec3_d *a, const GeVec3_d *b);
GeVec3_d GeVec3_MatOp(const GeVec3_d *vec, const GeMat3_d *mat);
GeVec3_d GeVec3_Rotate_x(const GeVec3_d *vec, double phi);
GeVec3_d GeVec3_Rotate_y(const GeVec3_d *vec, double phi);
GeVec3_d GeVec3_Rotate_z(const GeVec3_d *vec, double phi);
GeVec3_d GeVec3_Rotate(const GeVec3_d *vec, const GeCam *cam);

GeVec3_d GeVec3_Sph2Cart( double R, double theta, double phi);
GeVec3_d GeVec3_Cyl2Cart( double Rc, double phi, double Z);

GeVec3_d GeVec3_Cart2Sph( const GeVec3_d *);
GeVec3_d GeVec3_Cart2Cyl( const GeVec3_d *);

int point_in_voxel2(const GeVec3_d *pt, const GeVox *voxel, size_t axis);
int point_in_voxel_cyl3d(const GeVec3_d *pt, const GeVox *voxel, size_t axis);

size_t Ge_PosToIelem(size_t i, size_t j, size_t k, const GeVec3_s *naxes);
size_t Ge_IndexToIelem(const GeVec3_s *idx, const GeVec3_s *naxes);

GeVec3_s Ge_IelemToIndex(size_t ielem, const GeVec3_s *naxes);
GeVox GeVox_GetSubVox(const GeVox *voxel, const GeVec3_s *idx, const GeVec3_s *naxes);

GeCam GeCam_Init(double phix, double phiy, double phiz);

void GeRay_AWInit(GeRay *ray, const GeVox *voxel);
void GeRay_AWTraverse(GeRay *ray, double *tmin, size_t *plane);
GeVec3_d GeRay_AWPos(const GeRay *ray);

GeRay GeRay_Init(double ex, double ey, double ez, double dx, double dy, double dz);
double GeRay_IntersectPlane(const GeRay *ray, const GeVec3_d *n, const GeVec3_d *q);
void GeRay_IntersectSphere(const GeRay *ray, double r, double *t1, double *t2);

void GeRay_IntersectSphere3d(const GeRay *ray, double r, double *t1, double *t2);
void GeRay_IntersectTheta(const GeRay *ray, double theta, double *t1, double *t2);
void GeRay_IntersectPhi(const GeRay *ray, double phi, double *t);

void GeRay_IntersectRc(const GeRay *ray, double Rc, double *t1, double *t2);
void GeRay_IntersectHz(const GeRay *ray, double Hz, double *t);

GeRay GeRay_Inc(const GeRay *ray, double t);
GeRay GeRay_Rotate(const GeRay *ray, int axis, double phi);
int GeRay_IntersectVoxel(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side);
int GeRay_IntersectVoxel_sph1d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side);
int GeRay_IntersectVoxel_sph3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side);
int GeRay_IntersectVoxel_rec3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side);
int GeRay_IntersectVoxel_cyl3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side);

void GeRay_TraverseVoxel(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side);
void GeRay_TraverseVoxel_sph1d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side);
void GeRay_TraverseVoxel_sph3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side);
void GeRay_TraverseVoxel_rec3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side);
void GeRay_TraverseVoxel_cyl3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side);

GeRay GeRay_Rand(gsl_rng *rng, const GeVox *voxel);
GeRay GeRay_Rand_sph1d(gsl_rng *rng, const GeVox *voxel);
GeRay GeRay_Rand_sph3d(gsl_rng *rng, const GeVox *voxel);
GeRay GeRay_Rand_rec3d(gsl_rng *rng, const GeVox *voxel);
GeRay GeRay_Rand_cyl3d(gsl_rng *rng, const GeVox *voxel);

#if Sp_MIRSUPPORT
void GeVec3_Cpgpt1(const GeVec3_d *vec, int sty, const GeCam *cam);
void GeVec3_Cpgline2(const GeVec3_d *v1, const GeVec3_d *v2, const GeCam *cam);
void GeVec3_Cpgarro2(const GeVec3_d *v1, const GeVec3_d *v2, const GeCam *cam);
void GeVox_Cpgenv(const GeVox *voxel);
void GeVox_Cpgplot(const GeVox *voxel, const GeCam *cam);
void GeRay_Cpgarro(const GeRay *ray, const GeCam *cam);
#endif

size_t GeVox_VertIndex2Pos(size_t i, size_t j, size_t k);
GeVox GeVox_Init(int cosys, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax);
GeVox GeVox_Init2(int cosys, GeVec3_d min, GeVec3_d max);

#endif



