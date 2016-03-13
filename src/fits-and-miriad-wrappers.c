/*
 * Wrapper routines for accessing Miriad datasets
 * All functions must begin with the prefix `MirWr'
 */
#define Sp_MIRSUPPORT MIRSUPPORT

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>


#if Sp_MIRSUPPORT
#include <maxdimc.h> /* The Miriad header for defining special dimensions */
#include <miriad.h>
#endif

#include "physics.h"
#include "numerical.h"
#include "memory.h"
#include "fits-and-miriad-wrappers.h"
#include "debug.h"

/*----------------------------------------------------------------------------*/

MirVar *MirVar_ListLookup(MirVar vars[], const char *name)
{
	int i;

	for(i = 0; vars[i].name; i++) {
		if(!strcmp(vars[i].name, name))
			return &vars[i];
	}

	return NULL;
}



/*----------------------------------------------------------------------------*/



void MirVar_SetWidth(MirVar vars[], const char *name, int n)
{
	MirVar *vp;

	vp = MirVar_ListLookup(vars, name);
	vp->n = n;

	return;
}

/*----------------------------------------------------------------------------*/

void MirVar_SetValue_a(MirVar vars[], const char *name, const char *value)
{
	MirVar_SetValue(vars, name, value);

	return;
}

/*----------------------------------------------------------------------------*/

void MirVar_SetValue_i(MirVar vars[], const char *name, int value)
{
	MirVar_SetValue(vars, name, &value);

	return;
}

/*----------------------------------------------------------------------------*/

void MirVar_SetValue_r(MirVar vars[], const char *name, double value)
{
	float fdum = (float)value;

	MirVar_SetValue(vars, name, &fdum);

	return;
}

/*----------------------------------------------------------------------------*/

void MirVar_SetValue_d(MirVar vars[], const char *name, double value)
{
	MirVar_SetValue(vars, name, &value);

	return;
}

/*----------------------------------------------------------------------------*/

void MirVar_SetValue(MirVar vars[], const char *name, const void *value)
{
	int i;
	MirVar *var = NULL;

	/*locate variable*/
	var = MirVar_ListLookup(vars, name);
	Deb_ASSERT(var != NULL);

	/*free original values*/
	if(var->value) {
		free(var->value);
		var->value = NULL;
	}

	/*write new values*/
	switch(var->type) {
			case 'a':
				var->value = Mem_STRDUP((const char *)value);
				break;

			case 'j':
				var->value = Mem_CALLOC2(var->n, sizeof(int));
				for(i = 0; i < var->n; i++) {
					*((int *)var->value+i) = *((const int *)value+i);
				}
				break;

			case 'i':
				var->value = Mem_CALLOC2(var->n,sizeof(int));
				for(i = 0; i < var->n; i++) {
					*((int *)var->value+i) = *((const int *)value+i);
				}
				break;

			case 'r':
				var->value = Mem_CALLOC2(var->n,sizeof(float));
				for(i = 0; i < var->n; i++) {
					*((float *)var->value+i) = *((const float *)value+i);
				}
				break;

			case 'd':
				var->value = Mem_CALLOC2(var->n,sizeof(double));
				for(i = 0; i < var->n; i++) {
					*((double *)var->value+i) = *((const double *)value+i);
				}
				break;

			case 'c':
				var->value = Mem_CALLOC2(var->n*2,sizeof(float));
				for(i = 0; i < var->n; i++) {
					*((float *)var->value+i*2) = *((const float *)value+i*2);
					*((float *)var->value+i*2+1) = *((const float *)value+i*2+1);
				}
				break;

		default: /* Shouldn't happen */
			Deb_ASSERT(0);
	}

	return;
}

/*----------------------------------------------------------------------------*/

void MirVar_ListCleanup(MirVar vars[])
{
	MirVar *var = NULL;

	for(var = &vars[0]; var->name; var++) {
		if(var->value) {
			free(var->value);
			var->value = NULL;
		}
	}

	return;
}

/*----------------------------------------------------------------------------*/

MirImg *MirImg_Alloc(MirImg_Axis x, MirImg_Axis y, MirImg_Axis v)
{
	MirImg *img = Mem_CALLOC(1, img);

	/* Some error checking */
	Deb_ASSERT(x.n > 0);
	Deb_ASSERT(y.n > 0);
	Deb_ASSERT(v.n > 0);

	img->x = x;
	img->y = y;
	img->v = v;
	/* First index of image is V!!! */
	img->cube = Mem_CALLOC(v.n * x.n * y.n, img->cube);

	return img;
}

/*----------------------------------------------------------------------------*/

void MirImg_Free(void *img_p)
{
	MirImg *img = img_p;

	if(img->cube)
		free(img->cube);

	free(img);
}


/*----------------------------------------------------------------------------*/

MirFile *MirXY_Open_new(const char *name, size_t nx, size_t ny, size_t nv)
{
	int nsize[3] = {(int)nx, (int)ny, (int)nv}, tno;
	MirFile *fp;
        #if Sp_MIRSUPPORT
	xyopen_c(&tno, name, "new", 3, nsize);
        #endif
	fp = Mem_CALLOC(1, fp);
	fp->name = Mem_STRDUP(name);
	fp->tno = tno;

	return fp;
}

#if Sp_MIRSUPPORT
/*----------------------------------------------------------------------------*/

void MirVar_ListRead(MirVar vars[], MirFile *fp)
{
	size_t i;
	char sbuff[BUFSIZ] = "";

	for(i = 0; vars[i].name; i++) {
		switch(vars[i].type) {
			case 'a':
				uvgetvra_c(fp->tno, vars[i].name, sbuff, BUFSIZ);
				vars[i].value = Mem_STRDUP(sbuff);
				break;

			case 'j':
				vars[i].value = Mem_CALLOC2(vars[i].n, sizeof(int));
				uvgetvrj_c(fp->tno, vars[i].name, vars[i].value, vars[i].n);
				break;

			case 'i':
				vars[i].value = Mem_CALLOC2(vars[i].n, sizeof(int));
				uvgetvri_c(fp->tno, vars[i].name, vars[i].value, vars[i].n);
				break;

			case 'r':
				vars[i].value = Mem_CALLOC2(vars[i].n, sizeof(float));
				uvgetvrr_c(fp->tno, vars[i].name, vars[i].value, vars[i].n);
				break;

			case 'd':
				vars[i].value = Mem_CALLOC2(vars[i].n, sizeof(double));
				uvgetvrd_c(fp->tno, vars[i].name, vars[i].value, vars[i].n);
				break;

			case 'c':
				vars[i].value = Mem_CALLOC2(vars[i].n*2, sizeof(float));
				uvgetvrc_c(fp->tno, vars[i].name, vars[i].value, vars[i].n);
				return;

			default: /* Should not happen */
				Deb_ASSERT(0);
				return;
		}
	}

	return;
}

/*----------------------------------------------------------------------------*/
void MirVar_ListWrite(MirVar vars[], MirFile *fp)
{
	size_t i;

	for(i = 0; vars[i].name; i++) {
		switch(vars[i].type) {
			case 'a':
				//uvputvra_c(fp->tno, vars[i].name, vars[i].value);
				uvputvr_c(fp->tno, H_BYTE, vars[i].name, vars[i].value, (int)strlen(vars[i].value));
				break;

			case 'j':
				uvputvrj_c(fp->tno,vars[i].name,vars[i].value,vars[i].n);
				break;

			case 'i':
				uvputvri_c(fp->tno,vars[i].name,vars[i].value,vars[i].n);
				break;

			case 'r':
				uvputvrr_c(fp->tno,vars[i].name,vars[i].value,vars[i].n);
				break;

			case 'd':
				uvputvrd_c(fp->tno,vars[i].name,vars[i].value,vars[i].n);
				break;

			case 'c':
				uvputvrc_c(fp->tno,vars[i].name,vars[i].value,vars[i].n);
				return;

			default: /* Shouldn't happen */
				Deb_ASSERT(0);
		}
	}

	return;
}

/*----------------------------------------------------------------------------*/

MirFile *MirUV_Open(const char *name, const char *mode)
{
	int tno;
	MirFile *fp;

	uvopen_c(&tno, name, mode);
	fp = Mem_CALLOC(1, fp);
	fp->name = Mem_STRDUP(name);
	fp->tno = tno;

	return fp;
}

/*----------------------------------------------------------------------------*/

void MirUV_Close(void *fp)
{
	MirFile *f = fp;

	Deb_ASSERT(f->tno > 0);

	free(f->name);
	uvclose_c(f->tno);
	free(fp);

	return;
}



/*----------------------------------------------------------------------------*/

MirFile *MirXY_Open_old(const char *name, size_t *nx, size_t *ny, size_t *nv)
/* In "old" mode of xyopen_c, the cube dimensions are stored in the nsize
 * parameter upon success */
{
	int nsize[3] = {0, 0, 0}, tno;
	MirFile *fp;

	xyopen_c(&tno, name, "old", 3, nsize);
	
	*nx = (size_t)nsize[0];
	*ny = (size_t)nsize[1];
	*nv = (size_t)nsize[2];

	fp = Mem_CALLOC(1, fp);
	fp->name = Mem_STRDUP(name);
	fp->tno = tno;

	return fp;
}

/*----------------------------------------------------------------------------*/

void MirXY_Close(void *fp)
{
	MirFile *f = fp;

	Deb_ASSERT(f->tno > 0);

	free(f->name);
	xyclose_c(f->tno);
	free(fp);

	return;
}

/*----------------------------------------------------------------------------*/

void MirWr_SetBugLabel(const char *name)
{
	buglabel_c(name);

	return;
}

/*----------------------------------------------------------------------------*/

double *MirImg_LoadXYV(const char *fname, MirImg_Axis *x, MirImg_Axis *y, MirImg_Axis *v, char **bunit)
{
	MirFile *xy;
	MirImg *img;
	size_t i, j, k;
	double *cube = 0;

	Mem_BZERO(x);
	Mem_BZERO(y);
	Mem_BZERO(v);

	xy = MirXY_Open_old(fname, &x->n, &y->n, &v->n);
	img = MirImg_Alloc(*x, *y, *v);
	MirImg_ReadXY(xy, img, bunit);

	*x = img->x;
	*y = img->y;
	*v = img->v;
	cube = Mem_CALLOC(x->n * y->n * v->n, cube);

	for(k = 0; k < img->v.n; k++) {
		for(j = 0; j < img->y.n; j++) {
			for(i = 0; i < img->x.n; i++) {
				cube[k + v->n * (j + y->n * i)] = MirImg_PIXEL(*img, k, i, j);
			}
		}
	}

	MirImg_Free(img);
	MirXY_Close(xy);

	return cube;
}

/*----------------------------------------------------------------------------*/

void MirImg_ReadXY(MirFile *fp, MirImg *image, char **bunit)
{
	char sval[BUFSIZ] = "";
	int i, j, k;
	int ind;
	float *row = NULL;

	/************************************************************************/
	/*                           axis 1--RA                                 */
	/************************************************************************/
	rdhda_c(fp->tno,"ctype1", sval, "", BUFSIZ);
	rdhdd_c(fp->tno,"crpix1", &image->x.crpix, 0.0);
	Deb_ASSERT(image->x.crpix > 0);
	image->x.crpix -= 1; /* Indices are offset-1 in Miriad */
	rdhdd_c(fp->tno,"crval1", &image->x.crval, 0.0);
	rdhdd_c(fp->tno,"cdelt1", &image->x.delt, 0.0);

	/************************************************************************/
	/*                          axis 2--DEC                                 */
	/************************************************************************/
	rdhda_c(fp->tno,"ctype2", sval, "", BUFSIZ);
	rdhdd_c(fp->tno,"crpix2", &image->y.crpix, 0.0);
	Deb_ASSERT(image->y.crpix > 0);
	image->y.crpix -= 1; /* Indices are offset-1 in Miriad */
	rdhdd_c(fp->tno,"crval2", &image->y.crval, 0.0);
	rdhdd_c(fp->tno,"cdelt2", &image->y.delt, 0.0);

	/************************************************************************/
	/*                          axis 3--VEL                                 */
	/************************************************************************/
	rdhda_c(fp->tno,"ctype3", sval, "", BUFSIZ);
	rdhdd_c(fp->tno,"crpix3", &image->v.crpix, 0.0);
	Deb_ASSERT(image->v.crpix > 0);
	image->v.crpix -= 1; /* Indices are offset-1 in Miriad */
	rdhdd_c(fp->tno,"crval3", &image->v.crval, 0.0);
	rdhdd_c(fp->tno,"cdelt3", &image->v.delt, 0.0);
	image->v.crval *= 1000.0; /* km/s -> m/s */
	image->v.delt *= 1000.0; /* km/s -> m/s */

	/* bunit */
	strcpy(sval,"");
	rdhda_c(fp->tno,"bunit",sval,"",BUFSIZ);
	if(bunit && strcmp(sval,"")) {
		*bunit = Mem_STRDUP(sval);
	}

	/* Rest frequency */
	rdhdd_c(fp->tno, "restfreq", &image->restfreq, 0.0);
	image->restfreq *= 1.0e9; /* GHz -> Hz */

	Deb_ASSERT(image->cube != NULL);

	/* Loop through channels */
	for(k = 0; k < (int)image->v.n; k++) {
		ind = k + 1;
		xysetpl_c(fp->tno, 1, &ind);

		/* Loop through rows */
		for(j = 0; j < (int)image->y.n; j++) {
			/* Allocate row */
			row = Mem_CALLOC(image->x.n, row);

			/* Read row from file */
			xyread_c(fp->tno, j + 1, row);

			/* Load row into cube */
			for(i = 0; i < (int)image->x.n; i++) {
				MirImg_PIXEL(*image, k, i, j) = row[i];
			}

			/* Free row */
			free(row);
		}
	}

	return;
}

/*----------------------------------------------------------------------------*/

void MirImg_WriteXY(MirFile *fp, const MirImg *image, const char *bunit, double bfac)
/* Write cube of data (with dimensions x->n * y->n * v->n) to a
 * Miriad image dataset. All error handling is done by Miriad routines, which
 * are all fatal -- so I guess I don't need to bother with error handling here.
 *
 * Values of ImgAxis->delt must be physically reasonable, or Miriad will have
 * trouble displaying the image. */
{
	size_t i, j, k;
	int idx;
	float *row;
        //printf("xcrp=%d xcrv=%e x\n");
	/************************************************************************/
	/*                           Axis 1--RA                                 */
	/************************************************************************/
	wrhda_c(fp->tno, "ctype1", "RA---SIN");
	wrhdd_c(fp->tno, "crpix1", image->x.crpix + 1.0); /* Indices are offset-1 in Miriad */
	wrhdd_c(fp->tno, "crval1", image->x.crval);
	wrhdd_c(fp->tno, "cdelt1", image->x.delt);

	/************************************************************************/
	/*                           Axis 2--DEC                                */
	/************************************************************************/
	wrhda_c(fp->tno, "ctype2", "DEC--SIN");
	wrhdd_c(fp->tno, "crpix2", image->y.crpix + 1.0); /* Indices are offset-1 in Miriad */
	wrhdd_c(fp->tno, "crval2", image->y.crval);
	wrhdd_c(fp->tno, "cdelt2", image->y.delt);

	/************************************************************************/
	/*                           Axis 3--VEL                                */
	/************************************************************************/
	wrhda_c(fp->tno, "ctype3", "VELO-LSR");
	wrhdd_c(fp->tno, "crpix3", image->v.crpix + 1.0); /* Indices are offset-1 in Miriad */
	wrhdd_c(fp->tno, "crval3", image->v.crval);
	wrhdd_c(fp->tno, "cdelt3", image->v.delt / 1000.0); /* Convert to km/s */

	/* Brightness unit */
	wrhda_c(fp->tno, "bunit", bunit);

	/* Rest frequency */
	wrhdd_c(fp->tno, "restfreq", image->restfreq / 1.0e9); /* Hz -> GHz */

	/* Loop through all channels and write cube data */
	for(k = 0; k < image->v.n; k++) {
		/* Select appropriate plane (indices start with 1 as always) */
		idx = (int)k + 1;
		xysetpl_c(fp->tno, 1, &idx);

		/* Loop through rows */
		for(j = 0; j < image->y.n; j++) {
			/* Allocate row buffer */
			row = Mem_CALLOC(image->x.n, row);

			/* Load cube data into row */
			for(i = 0; i < image->x.n; i++)
				row[i] = (float)(MirImg_PIXEL(*image, k, i, j) * bfac);

			/* Write row of data to file */
			idx = (int)j + 1;
			xywrite_c(fp->tno, idx, row);

			/* Free row buffer */
			free(row);
		}
	}

	return;
}

/*----------------------------------------------------------------------------*/

void MirImg_UVResamp(MirImg *image, MirFile *uvin, MirFile *uvout)
/* this program simulates the effect of interferometric observation
   by doing an FFT on the 'image' channel maps, and producing
   nvis visbilities according to the Miriad dataset (tno) specified by
   the user. */
{
	size_t i, nx = image->x.n, ny = image->y.n, nchan = image->v.n;
	int *flags, /* nwide, */ nspect, nants, nread = 0;
	double *preamble, *chan, *chan_real, *chan_imag, *uarr, *varr,
	       umax, vmax, uu, vv, restfreq, dfreq, sfreq;
	MirImg *fft_real, *fft_imag;
	float *data;
	MirVar uv_vars[] = {
		/* Variables that are commented out, are those that are not
		 * absolutely necessary */
		//{'a', "observer", NULL, 0},
		{'a', "telescop", NULL, 0},
		{'a', "version", NULL, 0},
		//{'d', "latitud", NULL, 1},
		{'i', "nants", NULL, 1},
		//{'r', "pbfwhm", NULL, 1},
		{'a', "source", NULL, 0},
		{'d', "ra", NULL, 1},
		{'d', "dec", NULL, 1},
		{'r', "epoch", NULL, 1},
		{'d', "obsra", NULL, 1},
		{'d', "obsdec", NULL, 1},
		{'a', "veltype", NULL, 1},
		{'r', "veldop", NULL, 1},
		{'r', "vsource", NULL, 1},
		{'d', "lo1", NULL, 1},
		//{'i', "nwide", NULL, 1},
		//{'r', "wfreq", NULL, 0},
		//{'r', "wwidth", NULL, 0},
		//{'i', "numchan", NULL, 1},
		{'i', "nspect", NULL, 1},
		{'d', "restfreq", NULL, 0},
		{'d', "sfreq", NULL, 0},
		{'d', "sdf", NULL, 0},
		{'i', "nschan", NULL, 0},
		{'i', "ischan", NULL, 0},
		{'i', "npol", NULL, 1},
		{'i', "cormode", NULL, 1},
		{'r', "corfin", NULL, 4},
		//{'r', "corbw", NULL, 2},
		{'d', "antpos", NULL, 0},
		{0, 0, 0, 0}
	},
	integ_vars[] = {
		{'r', "inttime", NULL, 1},
		{'d', "ut", NULL, 1},
		{'d', "lst", NULL, 1},
		{'r', "wsystemp", NULL, 0},
		{'r', "systemp", NULL, 0},
		{'i', "pol", NULL, 1},
		{'r', "jyperk", NULL, 1},
		{0, 0, 0, 0}
	};

	/* Sanity checks */
	Deb_ASSERT(image->restfreq > 0);
	Deb_ASSERT(image->v.delt > 0);

	/* Allocate arrays */
	flags = Mem_CALLOC(MAXCHAN, flags);
	preamble = Mem_CALLOC(MAXCHAN, preamble);
	data = Mem_CALLOC(MAXCHAN, data);

	/* Acquire info from input dataset */
	/* Read first record (required by Miriad UV format) */
	uvread_c(uvin->tno, preamble, data, flags, MAXCHAN, &nread);

	/************************************************************************/
	/*                           Setup Header                               */
	/************************************************************************/
	/* Set data width of UV variables (must be set before reading
	 * from file or Miriad will complain about it */
	uvgetvri_c(uvin->tno, "nspect", &nspect, 1);
	uvgetvri_c(uvin->tno, "nants", &nants, 1);
	MirVar_SetWidth(uv_vars, "restfreq", nspect);
	MirVar_SetWidth(uv_vars, "sfreq", nspect);
	MirVar_SetWidth(uv_vars, "sdf", nspect);
	MirVar_SetWidth(uv_vars, "nschan", nspect);
	MirVar_SetWidth(uv_vars, "ischan", nspect);
	MirVar_SetWidth(uv_vars, "antpos", nants * 3);

	/* Wide-band variables -- currently defunct */
	//uvgetvri_c(uvin->tno, "nwide", &nwide, 1);
	//MirVar_SetWidth(uv_vars, "wfreq", nwide);
	//MirVar_SetWidth(uv_vars, "wwidth", nwide);

	/* Set data width of integration variables */
	MirVar_SetWidth(integ_vars, "wsystemp",  nants * (2 + nspect));
	MirVar_SetWidth(integ_vars, "systemp", nants * nspect);

	/* Read UV variables from input UV dataset */
	MirVar_ListRead(uv_vars, uvin);

	/* One image contains only one spectral window */
	nspect = 1;

	/* Widths must be reset if nspect has changed */
	MirVar_SetWidth(uv_vars, "restfreq", nspect);
	MirVar_SetWidth(uv_vars, "sfreq", nspect);
	MirVar_SetWidth(uv_vars, "sdf", nspect);
	MirVar_SetWidth(uv_vars, "nschan", nspect);
	MirVar_SetWidth(uv_vars, "ischan", nspect);

	/* Calculate frequency spacing from image velocity difference */
	restfreq = image->restfreq;

	/* Channel width is minus Doppler-shifted velocity width */
	dfreq = -(Phys_DoppShiftFreq(restfreq, image->v.delt) - restfreq);
	
	/* Channel 1 frequency is restfreq + bw/2:
	 * high-freq --> blueshift --> negative velocity
	 * low-freq  --> redshift  --> positive velocity
	 */
	sfreq = restfreq + 0.5 * fabs((int)nchan * dfreq);

	/* Modify UV variables */
	//MirVar_SetValue_i(uv_vars, "numchan", (int)nchan);
	//MirVar_SetValue_r(uv_vars, "veldop", 0.0); /* Observatory velocity */
	MirVar_SetValue_i(uv_vars, "nspect", nspect); /* Number of spectral windos */
	MirVar_SetValue_i(uv_vars, "nschan", (int)nchan); /* Number of spectral channels */
	MirVar_SetValue_i(uv_vars, "ischan", 1); /* Starting channel */
	MirVar_SetValue_d(uv_vars, "restfreq", restfreq / 1.0e9); /* Window rest frequency: Hz -> GHz */
	MirVar_SetValue_d(uv_vars, "sfreq", sfreq / 1.0e9); /* Frequency of starting channel: Hz -> GHz */
	MirVar_SetValue_d(uv_vars, "sdf", dfreq / 1.0e9); /* Channel width: Hz -> GHz */

	/* Write UV variables to output UV dataset */
	MirVar_ListWrite(uv_vars, uvout);

	/* `npol' must be present in the header, or the file can't be opened
	 * by Miriad. We can handle only one polarization for now */
	wrhdi_c(uvout->tno, "npol", 1);

	/* Read integration variables from input dataset */
	MirVar_ListRead(integ_vars, uvin);

	/************************************************************************/
	/*                              FFT                                     */
	/************************************************************************/
	/* Datasets all set, now FFT all channel maps to prepare
	 * for interpolation in uv-plane */
	fft_real = MirImg_Alloc(image->x, image->y, image->v);
	fft_imag = MirImg_Alloc(image->x, image->y, image->v);
	for(i = 0; i < nchan; i++) {
		chan = MirImg_CHAN(*image, i);
		chan_real = MirImg_CHAN(*fft_real, i);
		chan_imag = MirImg_CHAN(*fft_imag, i);
		NumFFT_Xform2d(chan, nx, ny, chan_real, chan_imag);
	}

	/************************************************************************/
	/*                   Read and interpolate UV data                       */
	/************************************************************************/
	/* Setup uv-plane axes */
	uarr = Mem_CALLOC(nx, uarr);
	varr = Mem_CALLOC(ny, varr);
	Deb_ASSERT((image->x.delt > 0) && (image->y.delt > 0));

	/* U_max and V_max are reciprocal of pixel size:
	 * radians -> n_lambda */
	umax = 1.0 / image->x.delt;
	vmax = 1.0 / image->y.delt;
	for(i = 0; i < nx; i++) {
		uarr[i] = -umax + (umax * 2.0 / (double)nx) * (double)i;
	}
	for(i = 0; i < ny; i++) {
		varr[i] = -vmax + (vmax * 2.0 / (double)ny) * (double)i;
	}

	/* Loop through visibilities */
	while(nread > 0) {
		/* Unpack the preamble: preamble[0] and preamble[1] are uv coordinates
		 * in nanoseconds; multiply be restfreq to convert to n_lambda */
		uu = preamble[0] * 1.0E-9 * restfreq;
		vv = preamble[1] * 1.0E-9 * restfreq;

		/* interpolate visibilities */
		for(i = 0; i < nchan; i++) {
			/* Pointer to real and imag fft maps */
			chan_real = MirImg_CHAN(*fft_real, i);
			chan_imag = MirImg_CHAN(*fft_imag, i);

			/* Interpolate ONLY if u and v within range */
			if(fabs(uu) <= umax && fabs(vv) <= vmax) {
				/* Data is in complex format -- i*2=real, i*2+1=imag */
				data[i*2] = (float)Num_InterpPoly2d(uu, vv, uarr, varr, chan_real, nx, ny, (size_t)3, (size_t)3);
				data[i*2+1] = (float)Num_InterpPoly2d(uu, vv, uarr, varr, chan_imag, nx, ny, (size_t)3, (size_t)3);
			}
			else {
				/* Otherwise flag as bad data -- set flag=0 */
				data[i*2] = 0.0;
				data[i*2+1] = 0.0;
				flags[i] = 0;
			}
		}

		/* Write integration variables to output dataset */
		MirVar_ListWrite(integ_vars, uvout);

		/* Write UV data to output dataset */
		uvwrite_c(uvout->tno, preamble, data, flags, (int)nchan);

		/* Read next record from input dataset */
		uvread_c(uvin->tno, preamble, data, flags, MAXCHAN, &nread);

		/* Read next set of integration variables from input dataset */
		MirVar_ListRead(integ_vars, uvin);
	}

	/* Cleanup */
	MirImg_Free(fft_real);
	MirImg_Free(fft_imag);
	free(flags);
	free(preamble);
	free(data);
	MirVar_ListCleanup(uv_vars);

	return;
}

#endif

void FITSoutput( MirFile *fp, MirImg *image, MirImg *StokesQ, MirImg *StokesU, const char *bunit, double scale, int Stokes)
{
        #define M_PI 3.14159265358979323846264338327950288419716939937510
        /* output FITS file  */
        int status = 0;
        fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
        
        char FileName[32];
        sprintf(FileName,"%s.fits", fp->name);
        
        /* createing new FITS file  */
        // if the file exist, remove it.
        if(access( FileName, F_OK ) != -1){ 
                fits_delete_file( fptr, &status);
                status = 0;
        }
        fits_create_file(&fptr, FileName, &status);   /* create new file */
        Deb_ASSERT(status == 0);

        long naxis;
        if (Stokes) naxis = 4;
        else naxis =3;
        
        long naxes[naxis];
        
        naxes[0] = (long) image->x.n;
        naxes[1] = (long) image->y.n;
        naxes[2] = (long) image->v.n;
        if(Stokes) naxes[3] = (long) 3;
        long nelements = 1;
        for (int i = 0; i < naxis; i++)
                nelements *= naxes[i];

        char 
        ctype1[32],
        ctype2[32],
        ctype3[32],
        ctype4[32],
        cellscal[32];
        sprintf(ctype1,"RA---SIN");
        sprintf(ctype2,"DEC--SIN");
        sprintf(ctype3,"VELO-LSR");
        sprintf(ctype4,"STOKES");
        sprintf(cellscal,"1/F     ");
        
        double 
        crpix1 = image->x.crpix + 1. ,
        crpix2 = image->y.crpix + 1. ,
        crpix3 = image->v.crpix + 1. , 
        crpix4 = 1.,
        cdelt1 = 180. / M_PI * image->x.delt,
        cdelt2 = 180. / M_PI * image->y.delt,
        cdelt3 = 0.001 * image->v.delt,
        cdelt4 = 1.,
        crval = 0.0,
        crval4 = 1.,
        beam = 0.0,
        bscale = 1.,
        bzero = 0.,
        restfreq = image->restfreq;

        double *array = Mem_CALLOC(nelements, array);;
        
        fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
        Deb_ASSERT(status == 0);
            /* Write a keyword; must pass the ADDRESS of the value */
         
        fits_write_key(fptr,TSTRING,"CTYPE1",&ctype1,"Type of first axis",&status);
        fits_write_key(fptr,TDOUBLE,"CRPIX1",&crpix1,"Reference of first axis",&status);
        fits_write_key(fptr,TDOUBLE,"CDELT1",&cdelt1,"Increment value of first axis",&status);
        fits_write_key(fptr,TDOUBLE,"CRVAL1",&crval,"Offset of first axis ",&status);
        fits_write_key(fptr,TSTRING,"CTYPE2",&ctype2,"Type of second axis",&status);
        fits_write_key(fptr,TDOUBLE,"CRPIX2",&crpix2,"Reference of second axis",&status);
        fits_write_key(fptr,TDOUBLE,"CDELT2",&cdelt2,"Increment value of second axis",&status);
        fits_write_key(fptr,TDOUBLE,"CRVAL2",&crval,"Offset of second axis ",&status);
        fits_write_key(fptr,TSTRING,"CTYPE3",&ctype3,"Type of third axis",&status);
        fits_write_key(fptr,TDOUBLE,"CRPIX3",&crpix3,"Reference of third axis",&status);
        fits_write_key(fptr,TDOUBLE,"CDELT3",&cdelt3,"Increment value of third axis",&status);
        fits_write_key(fptr,TDOUBLE,"CRVAL3",&crval,"Offset of third axis ",&status);
        if(Stokes){
                fits_write_key(fptr,TSTRING,"CTYPE4",&ctype4,"Type of third axis",&status);
                fits_write_key(fptr,TINT,"CRPIX4",&crpix4,"Reference of third axis",&status);
                fits_write_key(fptr,TINT,"CDELT4",&cdelt4,"Increment value of third axis",&status);
                fits_write_key(fptr,TINT,"CRVAL4",&crval4,"Offset of third axis ",&status);
        }
        fits_write_key(fptr,TDOUBLE,"BMAJ",&beam,"Major beam axis",&status);
        fits_write_key(fptr,TDOUBLE,"BMIN",&beam,"Minor beam axis ",&status);
        fits_write_key(fptr,TDOUBLE,"BPA",&beam,"PA",&status);
        fits_write_key(fptr,TDOUBLE,"BSCALE",&bscale,"",&status);
        fits_write_key(fptr,TDOUBLE,"BZERO",&bzero,"",&status);
        fits_write_key(fptr,TDOUBLE,"RESTFREQ",&restfreq,"",&status);
        fits_write_key(fptr,TSTRING,"BUNIT",bunit,"",&status);
        fits_write_key(fptr,TSTRING,"CELLSCAL",&cellscal,"",&status);
        Deb_ASSERT(status == 0);
        
        
        if(Stokes) {
                for (int iv = 0; iv < naxes[2]; iv++) 
                 for (int iy = 0; iy < naxes[1]; iy++)      
                  for (int ix = 0; ix < naxes[0]; ix++) {
                        int idx = ( ( 0 * naxes[2] + iv) * naxes[1] + iy) * naxes[0] + ix;
                        array[idx] = scale * MirImg_PIXEL( *image, iv, ix, iy);
                }
                for (int iv = 0; iv < naxes[2]; iv++) 
                 for (int iy = 0; iy < naxes[1]; iy++) 
                  for (int ix = 0; ix < naxes[0]; ix++) {
                        
                        int idx = ( ( 1 * naxes[2] + iv) * naxes[1] + iy) * naxes[0] + ix;
                        array[idx] = scale * MirImg_PIXEL( *StokesQ, iv, ix, iy);
                }
                for (int iv = 0; iv < naxes[2]; iv++) 
                 for (int iy = 0; iy < naxes[1]; iy++) 
                  for (int ix = 0; ix < naxes[0]; ix++) {
                        int idx = ( ( 2 * naxes[2] + iv) * naxes[1] + iy) * naxes[0] + ix;
                        array[idx] = scale * MirImg_PIXEL( *StokesU, iv, ix, iy);
                }
        }
        else{
                for (int iv = 0; iv < naxes[2]; iv++)
                 for (int iy = 0; iy < naxes[1]; iy++)
                  for (int ix = 0; ix < naxes[0]; ix++) {
                        int idx = (iv * naxes[1] + iy) * naxes[0] + ix;
                        array[idx] = scale * MirImg_PIXEL(*image, iv, ix, iy);
                }
        }
        
        Deb_ASSERT(status == 0);
        long *fpixel;
        if (Stokes){
                fpixel = Mem_CALLOC( 4, fpixel);
                fpixel[0] = fpixel[1] = fpixel[2] = fpixel[3] = 1;
        }
        else{
                fpixel = Mem_CALLOC( 3, fpixel);
                fpixel[0] = fpixel[1] = fpixel[2] = 1;
        }
        
        fits_write_pix(fptr, TDOUBLE, fpixel, nelements, array, &status);
        free(fpixel);
        fits_close_file(fptr, &status);            /* close the file */
        Deb_ASSERT(status == 0);
        
        free(array);

        return;
}




