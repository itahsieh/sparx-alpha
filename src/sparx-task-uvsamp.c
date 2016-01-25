#include "sparx.h"

/* Global parameter struct */
static struct glb {
	MirImg_Axis x, y, v;
	MirImg *img;
	char *bunit;
	MirFile *xyin, *uvin, *uvout;
} glb;


#if Sp_MIRSUPPORT
int SpTask_UVSamp(void)
{
	int status = 0;

	/* Inputs --------------------------------------------------------------------*/
	Mem_BZERO(&glb);
	glb.xyin = SpInp_GetKey_mirxy_old("xyin", &glb.x.n, &glb.y.n, &glb.v.n);
	glb.uvin = SpInp_GetKey_miruv_old("uvin");
	glb.uvout = SpInp_GetKey_miruv_new("uvout");
	/*----------------------------------------------------------------------------*/

	/* Calculate uv-dataset sampled from image -----------------------------------*/
	/* Read image into memory */
	glb.img = MirImg_Alloc(glb.x, glb.y, glb.v);
	MirImg_ReadXY(glb.xyin, glb.img, &glb.bunit);

	/* Sanity check */
	if(glb.img->restfreq <= 0) {
		status = Err_SETSTRING("Rest frequency in image `%s' must be > 0", glb.xyin->name);
	}
	else if(glb.img->v.delt <= 0) {
		status = Err_SETSTRING("Channel width in image `%s' must be > 0", glb.xyin->name);
	}
	else {
		/* Do UV-resampling */
		MirImg_UVResamp(glb.img, glb.uvin, glb.uvout);
	}

	if(!status) {
		Sp_PRINT("Wrote UV-sampled visibilities to `%s'\n", glb.uvout->name);
	}
	/*----------------------------------------------------------------------------*/

	/* Cleanup -------------------------------------------------------------------*/
	MirXY_Close(glb.xyin);
	MirUV_Close(glb.uvin);
	MirUV_Close(glb.uvout);
	/*----------------------------------------------------------------------------*/

	return status;
}
#endif

















