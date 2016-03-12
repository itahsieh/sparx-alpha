#ifndef __MIRIAD_WRAPPERS_H__
#define __MIRIAD_WRAPPERS_H__

#include "fitsio.h"


enum {
	MirWr_NEW,
	MirWr_OLD
};

typedef struct {
	char type;
	const char *name;
	void *value;
	int n;
} MirVar;


typedef struct MirImg_Axis {
	size_t n;
	double delt, crpix, crval;
} MirImg_Axis;

typedef struct MirImg {
	MirImg_Axis x, y, v;
	double *cube, restfreq;
} MirImg;

typedef struct MirFile {
        char *name;
        int tno;
} MirFile;

#define MirImg_PIXEL(img, iv, ix, iy)\
	(img).cube[(size_t)(iy) + (img).y.n * ((size_t)(ix) + (img).x.n * (size_t)(iv))]

#define MirImg_CHAN(img, iv)\
	&MirImg_PIXEL((img), (iv), 0, 0)

/* Calculate reference pixel position
 * Note: CRPIX is assumed to be a floating point number! */
#define MirWr_CRPIX(n)\
	((n) == 1 ? 0.0 : (double)(n) / 2.0)

#define MirWr_IMGAXIS(n, delt, crpix, crval) {\
		(n), /* size_t, number of pixels */\
		(delt), /* double, pixel size -- values must be physically reasonable, or\
		        Miriad will have trouble displaying the image. */\
		(crpix), /* double, position of reference pixel */\
		(crval) /* double, value of reference pixel */\
	}

	
MirVar *MirVar_ListLookup(MirVar vars[], const char *name);
void MirVar_ListRead(MirVar vars[], MirFile *fp);
void MirVar_ListWrite(MirVar vars[], MirFile *fp);
void MirVar_SetWidth(MirVar vars[], const char *name, int n);
void MirVar_SetValue_a(MirVar vars[], const char *name, const char *value);
void MirVar_SetValue_i(MirVar vars[], const char *name, int value);
void MirVar_SetValue_r(MirVar vars[], const char *name, double value);
void MirVar_SetValue_d(MirVar vars[], const char *name, double value);
void MirVar_SetValue(MirVar vars[], const char *name, const void *value);
void MirVar_ListCleanup(MirVar vars[]);
double *MirImg_GetChannel(MirImg img, size_t i);
MirImg *MirImg_Alloc(MirImg_Axis x, MirImg_Axis y, MirImg_Axis v);
void MirImg_Free(void *img_p);
MirFile *MirUV_Open(const char *name, const char *mode);
void MirUV_Close(void *fp);
MirFile *MirXY_Open_new(const char *name, size_t nx, size_t ny, size_t nv);
MirFile *MirXY_Open_old(const char *name, size_t *nx, size_t *ny, size_t *nv);
void MirXY_Close(void *fp);
void MirWr_SetBugLabel(const char *name);
double *MirImg_LoadXYV(const char *fname, MirImg_Axis *x, MirImg_Axis *y, MirImg_Axis *v, char **bunit);
void MirImg_ReadXY(MirFile *fp, MirImg *image, char **bunit);
void MirImg_WriteXY(MirFile *fp, const MirImg *image, const char *bunit, double bfac);
void MirXY_WriteImgCube(int tno, const MirImg *image, char *bunit, double bfac);
void MirWr_WriteImgCube(int tno, const MirImg_Axis *x, const MirImg_Axis *y,
	const MirImg_Axis *v, const char *bunit, double bfac, double *cube);
void MirImg_UVResamp(MirImg *image, MirFile *uvin, MirFile *uvout);
void FITSoutput( MirFile *fp, MirImg *image, MirImg *StokesQ, MirImg *StokesU, const char *bunit, double scale, int Stokes);
#endif
