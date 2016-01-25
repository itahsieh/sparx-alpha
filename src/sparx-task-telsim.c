#include "sparx.h"

/* Global parameter struct */
static struct glb {
	int test, cont,overlap,lte,coldens,excit,zeeman;
	DatINode *unit;
	double ucon, overlap_vel;
	MirImg_Axis x, y, v;
	MirImg *xyv_img, *stokesv, *stokesq, *stokesu, *sigma2, *tau_img;
	double dist, rotate[3], I_norm, I_cmb;
	char *imgname;
	MirFile *xyv_imgf, *tau_imgf, *stokesvf, *stokesqf, *stokesuf;
	SpModel model;
	size_t line;
	double lamb, freq;
	size_t nsubres;
	struct {
		double blc_x, blc_y, trc_x, trc_y;
		size_t nsub;
	} *subres;
} glb;

enum {
	UNIT_K,
	UNIT_JYPX
};

static DatINode UNITS[] = {
	{"K", UNIT_K},
	{"JY/PIXEL", UNIT_JYPX},
	{0, 0}
};

#define MEAN_INT(ix, iy)\
	mean_int[(size_t)(iy) + glb.y.n * (size_t)(ix) ]
#define MEAN_VEL(ix, iy)\
	mean_vel[(size_t)(iy) + glb.y.n * (size_t)(ix) ]
#define MEAN_DEV(ix, iy)\
	mean_dev[(size_t)(iy) + glb.y.n * (size_t)(ix) ]

#define RELVEL(i,j)\
	(glb.model.parms.mol->OL[NRAD*i+j]->RelativeVel)
#define OVERLAP(i,j)\
	(glb.model.parms.mol->OL[NRAD*i+j]->overlap)

#define NRAD\
	(glb.model.parms.mol->nrad)
#define FREQ(i)\
	(glb.model.parms.mol->rad[i]->freq)


/* Subroutine prototypes */
static int InitModel(void);
static void *InitModelThread(void *tid_p);
static int CalcImage(void);
static void *CalcImageThread(void *tid_p);
static void RadiativeXfer(double dx, double dy, double *I_nu, double *V_nu, double *Q_nu, double *U_nu, double *sigma2, double *tau_nu);
static void ColumnDensityTracer(double dx, double dy, double *I_nu);
static void call_fits_output( const char *FileName, const char *bunit, const double scale, const int Stokes);
static void visualization(void);
/*----------------------------------------------------------------------------*/

int SpTask_Telsim(void)
{
	size_t i,j;
	int sts = 0;
	PyObject *o, *o1, *o2, *o3;
	

	Mem_BZERO(&glb);
	/* npix */
	if(!sts && !(sts = SpPy_GetInput_PyObj("npix", &o))) {
		glb.x.n = Sp_PYSIZE(Sp_PYLST(o, 0));
		glb.x.crpix = MirWr_CRPIX(glb.x.n);
		glb.y.n = Sp_PYSIZE(Sp_PYLST(o, 1));
		glb.y.crpix = MirWr_CRPIX(glb.y.n);
		SpPy_XDECREF(o);
	}
	/* cell */
	if(!sts && !(sts = SpPy_GetInput_PyObj("cell", &o))) {
		glb.x.delt = Sp_PYDBL(Sp_PYLST(o, 0));
		glb.y.delt = Sp_PYDBL(Sp_PYLST(o, 1));
		SpPy_XDECREF(o);
	}
	/* chan */
	if(!sts && !(sts = SpPy_GetInput_PyObj("chan", &o))) {
		glb.v.n = Sp_PYSIZE(Sp_PYLST(o, 0));
		glb.v.crpix = MirWr_CRPIX(glb.v.n);
		glb.v.delt = Sp_PYDBL(Sp_PYLST(o, 1));
		SpPy_XDECREF(o);
	}
	/* rotate */
	if(!sts && !(sts = SpPy_GetInput_PyObj("rotate", &o))) {
		glb.rotate[0] = Sp_PYDBL(Sp_PYLST(o, 0));
		glb.rotate[1] = Sp_PYDBL(Sp_PYLST(o, 1));
		glb.rotate[2] = Sp_PYDBL(Sp_PYLST(o, 2));
		SpPy_XDECREF(o);
	}
	/* subres */
	if(!sts && !(sts = SpPy_GetInput_PyObj("subres", &o))) {
		if(o != Py_None) {
			/* Get number of boxes */
			glb.nsubres = (size_t)PyList_Size(o);
			glb.subres = Mem_CALLOC(glb.nsubres, glb.subres);
			for(i = 0; i < glb.nsubres; i++) {
				o1 = PyList_GetItem(o, (Py_ssize_t)i);
				glb.subres[i].blc_x = Sp_PYDBL(PyList_GetItem(o1, (Py_ssize_t)0));
				glb.subres[i].blc_y = Sp_PYDBL(PyList_GetItem(o1, (Py_ssize_t)1));
				glb.subres[i].trc_x = Sp_PYDBL(PyList_GetItem(o1, (Py_ssize_t)2));
				glb.subres[i].trc_y = Sp_PYDBL(PyList_GetItem(o1, (Py_ssize_t)3));
				glb.subres[i].nsub = Sp_PYSIZE(PyList_GetItem(o1, (Py_ssize_t)4));
			}
		}
		SpPy_XDECREF(o);
	}

	/* source */
	if(!sts) {
		sts = SpPy_GetInput_model("source","source", &glb.model);
	}

	/* obs */
	if(!sts && !(sts = SpPy_GetInput_PyObj("obs", &o))) {
		o1 = PyObject_GetAttrString(o, "cont");
		o2 = PyObject_GetAttrString(o, "data");

		/* Temporary measure */
		Deb_ASSERT(o1 != NULL);
		//Deb_ASSERT(o2 != NULL);
		//Deb_ASSERT(o3 != NULL);
		glb.cont = Sp_PYINT(o1);
		if(glb.cont) {
			glb.lamb = Sp_PYDBL(o2);
			glb.freq = PHYS_CONST_MKS_LIGHTC / glb.lamb;
		}
		else {
			o3 = PyObject_GetAttrString(o, "coldens");
			glb.coldens = Sp_PYINT(o3);

			if(!glb.coldens){

				if(!sts) sts = SpPy_GetInput_bool("lte", &glb.lte);
				if(!sts){
					if(glb.lte) sts = SpPy_GetInput_molec("molec", &glb.model.parms.mol);
				}
				if(!sts) sts = SpPy_GetInput_bool("zeeman", &glb.zeeman);
							
				Deb_ASSERT(glb.model.parms.mol != NULL);

				glb.line = Sp_PYSIZE(o2);
				glb.freq = glb.model.parms.mol->rad[glb.line]->freq;
				glb.lamb = PHYS_CONST_MKS_LIGHTC / glb.freq;
				Deb_ASSERT(glb.line < glb.model.parms.mol->nrad);
				o2 = PyObject_GetAttrString(o, "overlap_int");
				o3 = PyObject_GetAttrString(o, "overlap_vel");
				glb.overlap = Sp_PYINT(o2);
				glb.overlap_vel = Sp_PYDBL(o3);
				
				/* construct overlapping table */
				if(glb.overlap){
					glb.model.parms.mol->OL = Mem_CALLOC(NRAD*NRAD,glb.model.parms.mol->OL);
					for(i=0; i<NRAD; i++){
						for(j=0; j<NRAD; j++){
							glb.model.parms.mol->OL[NRAD*i+j] = Mem_CALLOC(1,glb.model.parms.mol->OL[NRAD*i+j]);
							RELVEL(i,j) = ( 1e0-FREQ(j)/FREQ(i) )*CONSTANTS_MKS_LIGHT_C;
							if( fabs(RELVEL(i,j)) < glb.overlap_vel ){
								OVERLAP(i,j)=1;
							}
							else{
								OVERLAP(i,j)=0;
							}
						}
					}
				}
				
			}
			
		}
		
		SpPy_XDECREF(o1);
		SpPy_XDECREF(o2);
		SpPy_XDECREF(o3);
		SpPy_XDECREF(o);
	}
	if(!sts && !(sts = SpPy_GetInput_PyObj("unit", &o))) {
		/* unit */
		glb.unit = Dat_IList_NameLookup(UNITS, Sp_PYSTR(o));
		Deb_ASSERT(glb.unit != NULL);
		if(glb.unit->idx == UNIT_K)
			glb.ucon = Phys_RayleighJeans(glb.freq, 1.0);
		else if(glb.unit->idx == UNIT_JYPX)
			glb.ucon = (PHYS_UNIT_MKS_JY / (glb.x.delt * glb.y.delt));
		else
			Deb_ASSERT(0);
		/* Sanity check */
		if(!glb.coldens){
			Deb_ASSERT((glb.ucon > 0) && (!Num_ISNAN(glb.ucon)) && (glb.ucon < HUGE_VAL));
			SpPy_XDECREF(o);
		}

	}
	/* dist */
	if(!sts) sts = SpPy_GetInput_dbl("dist", &glb.dist);
	
	/* out (mandatory) */
	if(!sts){
		if(glb.coldens) glb.v.n=1;
		
		sts = SpPy_GetInput_mirxy_new("out", glb.x.n, glb.y.n, glb.v.n, &glb.xyv_imgf);
		
		const char *tempchar1[64],*tempchar2[64];
// 		if(glb.zeeman){
// 			sprintf(tempchar1,"stokesv_%s",glb.xyv_imgf->name);
// 			glb.stokesvf = MirXY_Open_new(tempchar1, glb.x.n, glb.y.n, glb.v.n);
// 		}
// 		else{
#if Sp_MIRSUPPORT
			if(glb.cont){
				sprintf(tempchar1,"stokesq_%s",glb.xyv_imgf->name);
				sprintf(tempchar2,"stokesu_%s",glb.xyv_imgf->name);
				glb.stokesqf = MirXY_Open_new(tempchar1, glb.x.n, glb.y.n, glb.v.n);
				glb.stokesuf = MirXY_Open_new(tempchar2, glb.x.n, glb.y.n, glb.v.n);
			}
#endif
// 		}
	}
	/* tau (optional) */
	if(!sts && SpPy_CheckOptionalInput("tau")) {
		sts = SpPy_GetInput_mirxy_new("tau", glb.x.n, glb.y.n, glb.v.n, &glb.tau_imgf);
	}
	/* excitation visualization */
	if(!sts && !glb.coldens && !glb.cont) sts = SpPy_GetInput_bool("excit", &glb.excit);

	/*
	 * Initialize model
	 */
	if(!sts) sts = InitModel();

        
	/*
	 * Synthesize image
	 */
	if(!sts) {
		/* Allocate image */	
// 		if(glb.zeeman){
// 			glb.stokesv = MirImg_Alloc(glb.x, glb.y, glb.v);
// 			glb.stokesv->restfreq = glb.freq;
// 		}
// 		else{
			glb.xyv_img = MirImg_Alloc(glb.x, glb.y, glb.v);
			glb.xyv_img->restfreq = glb.freq;
			if(glb.cont){
				glb.stokesq = MirImg_Alloc(glb.x, glb.y, glb.v);
				glb.stokesu = MirImg_Alloc(glb.x, glb.y, glb.v);
				glb.sigma2 = MirImg_Alloc(glb.x, glb.y, glb.v);
				glb.stokesq->restfreq = glb.freq;
				glb.stokesu->restfreq = glb.freq;
				glb.sigma2->restfreq = glb.freq;
			}
// 		}
		if(glb.tau_imgf)
			glb.tau_img = MirImg_Alloc(glb.x, glb.y, glb.v);
		/* Calculate image */

		sts = CalcImage();

	}

#if Sp_MIRSUPPORT
	/* writing out MIRIAD file */
	if(!sts) {
		/* Denormalize and convert image to proper units, then write cube to
		   Miriad image dataset */
		if(glb.coldens) {
			MirImg_WriteXY(glb.xyv_imgf, glb.xyv_img, glb.unit->name, 1e-7);
		}
		else {
			MirImg_WriteXY(glb.xyv_imgf, glb.xyv_img, glb.unit->name, glb.I_norm/glb.ucon);
			if(glb.cont){
				MirImg_WriteXY(glb.stokesqf, glb.stokesq, glb.unit->name, glb.I_norm/glb.ucon);
				MirImg_WriteXY(glb.stokesuf, glb.stokesu, glb.unit->name, glb.I_norm/glb.ucon);
			}
			Sp_PRINT("Wrote Miriad image to `%s'\n", glb.xyv_imgf->name);
		}
		if(glb.tau_imgf)
			MirImg_WriteXY(glb.tau_imgf, glb.tau_img, "Optical depth", 1.0);
	}
#endif

	char filename[32];
	
	
	if(!sts){
		
		sprintf(filename,"%s.fits", glb.xyv_imgf->name);

		if(glb.coldens) {
			// column density image
			call_fits_output( &filename, glb.unit->name, 1e-7, 0);
		}
		else{
			call_fits_output( &filename, glb.unit->name, glb.I_norm/glb.ucon, 0);
			if(glb.cont)
				call_fits_output( &filename, glb.unit->name, glb.I_norm/glb.ucon, 1);
		}
		Sp_PRINT("Wrote FITS image to `%s'\n", filename);		
		if(glb.tau_imgf){
			sprintf(filename,"%s.fits", glb.tau_imgf->name);
			call_fits_output( &filename, "Optical depth", 1, 0);
		}
	}

	/* write excitation visualization to VTK */
	if(glb.excit && !glb.cont)
		visualization();
	
	FILE *fp;
	size_t ix,iy,iz,iv;
	/* column density map to VTK */
	if(glb.coldens && !glb.cont){
		char filename[32];
		sprintf(filename,"%s.vtk",glb.xyv_imgf->name);
		fp=fopen(filename,"w");
	
		fprintf(fp,"# vtk DataFile Version 3.0\n");
		fprintf(fp,"Column density\n");
		fprintf(fp,"ASCII\n");
		fprintf(fp,"DATASET STRUCTURED_POINTS\n");
		fprintf(fp,"DIMENSIONS %zu %zu %d\n",glb.x.n,glb.y.n,1);
		fprintf(fp,"ORIGIN %f %f %d\n",-glb.x.crpix,-glb.y.crpix,1);
		fprintf(fp,"SPACING %d %d %d\n",1,1,1);
		fprintf(fp,"POINT_DATA %zu\n",glb.x.n * glb.y.n );
		fprintf(fp,"SCALARS CD float 1\n");
		fprintf(fp,"LOOKUP_TABLE default\n");	
		
		iv=0;
		for(iy = 0; iy < glb.y.n; iy++) {
			for(ix = 0; ix < glb.x.n; ix++) {
				fprintf(fp,"%11.4e\n",MirImg_PIXEL(*glb.xyv_img, iv, ix, iy));
			}	
		}
	
		fclose(fp);
		
	}

	
	// calculate three moments and output to VTK file
	/*
	if(!glb.cont){
	
	double *mean_int,*mean_vel,*mean_dev,dv,tempINT;
	mean_int = Mem_CALLOC(glb.x.n*glb.y.n, mean_int);
	mean_vel = Mem_CALLOC(glb.x.n*glb.y.n, mean_vel);
	mean_dev = Mem_CALLOC(glb.x.n*glb.y.n, mean_dev);
	
	
	for(ix = 0; ix < glb.x.n; ix++) {
		for(iy = 0; iy < glb.y.n; iy++) {
			MEAN_INT(ix,iy)=0.0;
			MEAN_VEL(ix,iy)=0.0;
			MEAN_DEV(ix,iy)=0.0;
			for(iv = 0; iv < glb.v.n; iv++) {
				tempINT = (MirImg_PIXEL(*glb.xyv_img, iv, ix, iy)-glb.I_cmb)*glb.I_norm/glb.ucon;
				MEAN_INT(ix,iy) += tempINT;
				MEAN_VEL(ix,iy) += tempINT * ((double)iv - glb.v.crpix) * glb.v.delt;
			}
			MEAN_INT(ix,iy) /= (double)glb.v.n;
			if (MEAN_INT(ix,iy)<1e-20){
				MEAN_VEL(ix,iy)=0.;
			}
			else{
				MEAN_VEL(ix,iy) /= MEAN_INT(ix,iy) * (double)glb.v.n;
			}
			for(iv = 0; iv < glb.v.n; iv++) {
				tempINT = (MirImg_PIXEL(*glb.xyv_img, iv, ix, iy)-glb.I_cmb)*glb.I_norm/glb.ucon;
			 	dv = ((double)iv - glb.v.crpix) * glb.v.delt - MEAN_VEL(ix,iy);
				MEAN_DEV(ix,iy) += tempINT*dv*dv;
			}
			if (MEAN_INT(ix,iy)<1e-20){
				MEAN_DEV(ix,iy)=0.;
			}
			else{
				MEAN_DEV(ix,iy) /= MEAN_INT(ix,iy) * (double)glb.v.n;
				MEAN_DEV(ix,iy) = sqrt(MEAN_DEV(ix,iy));
			}	
		}
	}
	
	char filename[32];
	
	sprintf(filename,"0thmov_%5s.vtk",glb.xyv_imgf->name);
	fp=fopen(filename,"w");
	fprintf(fp,"# vtk DataFile Version 3.0\n");
	fprintf(fp,"0thmov\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET STRUCTURED_POINTS\n");
	fprintf(fp,"DIMENSIONS %zu %zu %d\n",glb.x.n,glb.y.n,2);
	fprintf(fp,"ORIGIN %f %f %d\n",-glb.x.crpix,-glb.y.crpix,1);
	fprintf(fp,"SPACING %d %d %d\n",1,1,1);
	fprintf(fp,"POINT_DATA %zu\n",glb.x.n * glb.y.n * 2);
	fprintf(fp,"SCALARS Intensity float 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");	
	for (iz = 0; iz < 2; iz++){
		for(iy = 0; iy < glb.y.n; iy++) {
			for(ix = 0; ix < glb.x.n; ix++) {
				fprintf(fp,"%11.4e\n",MEAN_INT(ix,iy));
			}	
		}
	}
	fclose(fp);
	sprintf(filename,"1stmov_%5s.vtk",glb.xyv_imgf->name);
	fp=fopen(filename,"w");
	fprintf(fp,"# vtk DataFile Version 3.0\n");
	fprintf(fp,"1stmov\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET STRUCTURED_POINTS\n");
	fprintf(fp,"DIMENSIONS %zu %zu %d\n",glb.x.n,glb.y.n,2);
	fprintf(fp,"ORIGIN %f %f %d\n",-glb.x.crpix,-glb.y.crpix,-1);
	fprintf(fp,"SPACING %d %d %d\n",1,1,1);
	fprintf(fp,"POINT_DATA %zu\n",glb.x.n * glb.y.n * 2);
	fprintf(fp,"SCALARS Velocity float 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");		
	for (iz = 0; iz < 2; iz++){
		for(iy = 0; iy < glb.y.n; iy++) {
			for(ix = 0; ix < glb.x.n; ix++) {
				fprintf(fp,"%11.4e\n",MEAN_VEL(ix,iy));
			}	
		}
	}
	fclose(fp);
	sprintf(filename,"2ndmov_%5s.vtk",glb.xyv_imgf->name);
	fp=fopen(filename,"w");
	fprintf(fp,"# vtk DataFile Version 3.0\n");
	fprintf(fp,"2ndmov\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET STRUCTURED_POINTS\n");
	fprintf(fp,"DIMENSIONS %zu %zu %d\n",glb.x.n,glb.y.n,2);
	fprintf(fp,"ORIGIN %f %f %d\n",-glb.x.crpix,-glb.y.crpix,0);
	fprintf(fp,"SPACING %d %d %d\n",1,1,1);
	fprintf(fp,"POINT_DATA %zu\n",glb.x.n * glb.y.n * 2);
	fprintf(fp,"SCALARS Deviation float 1\n");
	fprintf(fp,"LOOKUP_TABLE default\n");	
	for (iz = 0; iz < 2; iz++){
		for(iy = 0; iy < glb.y.n; iy++) {
			for(ix = 0; ix < glb.x.n; ix++) {
				fprintf(fp,"%11.4e\n",MEAN_DEV(ix,iy));
			}	
		}
	}
	fclose(fp);
	
	}
	*/
	if(glb.cont){
		char filename[32];
		sprintf(filename,"stokesIQU_%s.vtk",glb.xyv_imgf->name);
		fp=fopen(filename,"w");
		fprintf(fp,"# vtk DataFile Version 3.0\n");
		fprintf(fp,"Stokes parameters\n");
		fprintf(fp,"ASCII\n");
		fprintf(fp,"DATASET STRUCTURED_POINTS\n");
		fprintf(fp,"DIMENSIONS %zu %zu %d\n",glb.x.n,glb.y.n,1);
		fprintf(fp,"ORIGIN %f %f %d\n",-glb.x.crpix,-glb.y.crpix,0);
		fprintf(fp,"SPACING %11.4e %11.4e %d\n",glb.x.delt,glb.y.delt,1);
		fprintf(fp,"POINT_DATA %zu\n",glb.x.n * glb.y.n);
		fprintf(fp,"SCALARS Stokes_I float 1\n");
		fprintf(fp,"LOOKUP_TABLE default\n");
		for(iy = 0; iy < glb.y.n; iy++) {
			for(ix = 0; ix < glb.x.n; ix++) {
				fprintf(fp,"%11.4e ",MirImg_PIXEL(*glb.xyv_img, 0, ix, iy));
			}
		}
		fprintf(fp,"\n");
		fprintf(fp,"SCALARS Stokes_Q float 1\n");
		fprintf(fp,"LOOKUP_TABLE default\n");
		for(iy = 0; iy < glb.y.n; iy++) {
			for(ix = 0; ix < glb.x.n; ix++) {
				fprintf(fp,"%11.4e ",MirImg_PIXEL(*glb.stokesq, 0, ix, iy));
			}
		}
		fprintf(fp,"\n");
		fprintf(fp,"SCALARS Stokes_U float 1\n");
		fprintf(fp,"LOOKUP_TABLE default\n");
		for(iy = 0; iy < glb.y.n; iy++) {
			for(ix = 0; ix < glb.x.n; ix++) {
				fprintf(fp,"%11.4e ",MirImg_PIXEL(*glb.stokesu, 0, ix, iy));
			}
		}
		fclose(fp);

		/*
		sprintf(filename,"stokesI_%s.dat",glb.xyv_imgf->name);
		fp=fopen(filename,"w");	
		for(iy = 0; iy < glb.y.n; iy++) {
			for(ix = 0; ix < glb.x.n; ix++) {
				fprintf(fp,"%5zu %5zu %11.4e\n",ix,iy,MirImg_PIXEL(*glb.xyv_img, 0, ix, iy));
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		sprintf(filename,"stokesQ_%s.dat",glb.xyv_imgf->name);
		fp=fopen(filename,"w");	
		for(iy = 0; iy < glb.y.n; iy++) {
			for(ix = 0; ix < glb.x.n; ix++) {
				fprintf(fp,"%5zu %5zu %11.4e\n",ix,iy,MirImg_PIXEL(*glb.stokesq, 0, ix, iy));
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		sprintf(filename,"stokesU_%s.dat",glb.xyv_imgf->name);
		fp=fopen(filename,"w");	
		for(iy = 0; iy < glb.y.n; iy++) {
			for(ix = 0; ix < glb.x.n; ix++) {
				fprintf(fp,"%5zu %5zu %11.4e\n",ix,iy,MirImg_PIXEL(*glb.stokesu, 0, ix, iy));
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
		*/
		sprintf(filename,"vector_%s.dat",glb.xyv_imgf->name);
		fp=fopen(filename,"w");		
		for(iy = 0; iy < glb.y.n; iy++) {
			for(ix = 0; ix < glb.x.n; ix++) {
				double pd,xi,vecx,vecy;
				pd=sqrt( MirImg_PIXEL(*glb.stokesq, 0, ix, iy)*MirImg_PIXEL(*glb.stokesq, 0, ix, iy)+ \
				         MirImg_PIXEL(*glb.stokesu, 0, ix, iy)*MirImg_PIXEL(*glb.stokesu, 0, ix, iy) )\
				   /(MirImg_PIXEL(*glb.xyv_img, 0, ix, iy) - MirImg_PIXEL(*glb.sigma2, 0, ix, iy));
				
// 				if (MirImg_PIXEL(*glb.stokesq, 0, ix, iy)==0.0){
//                   			if(MirImg_PIXEL(*glb.stokesu, 0, ix, iy)>0) 
//                   				xi=0.5*M_PI;
//                                         else xi=-0.5*M_PI;
// 				}
//                                else{
// 					xi=atan(MirImg_PIXEL(*glb.stokesu, 0, ix, iy)/MirImg_PIXEL(*glb.stokesq, 0, ix, iy));
// 				} 
// 				if (MirImg_PIXEL(*glb.stokesq, 0, ix, iy)<0){
// 					xi=xi+M_PI;
// 				}
// 				else if (xi<0){ 
// 					xi=xi+2*M_PI;
// 				}
				   
				xi=atan2(MirImg_PIXEL(*glb.stokesu, 0, ix, iy),MirImg_PIXEL(*glb.stokesq, 0, ix, iy));
				
				xi = 0.5 * xi;
				vecx=-pd*sin(xi);
				vecy=pd*cos(xi);
// 				vecx=pd*cos(xi);
// 				vecy=pd*sin(xi);
				
				fprintf(fp,"%5zu %5zu %11.4e %11.4e %11.4e %11.4e\n",ix,iy,vecx,vecy,pd,xi);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
	}
	



	/*
	 * Cleanup
	 */
	if(glb.xyv_img)
		MirImg_Free(glb.xyv_img);
 
	if(glb.stokesv)
		MirImg_Free(glb.stokesv);
	if(glb.stokesu)
		MirImg_Free(glb.stokesu);
        if(glb.stokesq)
		MirImg_Free(glb.stokesq);
        if(glb.sigma2)
		MirImg_Free(glb.sigma2);
	
	
	if(glb.tau_img)
		MirImg_Free(glb.tau_img);

	if(glb.imgname)
		free(glb.imgname);

	if(glb.subres)
		free(glb.subres);

	SpModel_Cleanup(glb.model);

	
#if Sp_MIRSUPPORT
	/* Miriad images must always be closed! */
	if(glb.xyv_imgf)
		MirXY_Close(glb.xyv_imgf);
	
	if(glb.stokesvf)
		MirXY_Close(glb.stokesvf);
	if(glb.stokesuf)
		MirXY_Close(glb.stokesuf);
	if(glb.stokesqf)
		MirXY_Close(glb.stokesqf);

	if(glb.tau_imgf)
		MirXY_Close(glb.tau_imgf);
#endif
	return sts;
}

/*----------------------------------------------------------------------------*/

static int InitModel(void)
{
	Zone *root = glb.model.grid, *zp;
	SpPhys *pp;
	int sts = 0;

// 	FILE *fp;
// 	double radius;

	/* Set normalization intensity to 20K -- normalization prevents rounding
	   errors from creeping in when flux values are very small */
	if(!glb.coldens){
		glb.I_norm = Phys_PlanckFunc(glb.freq, 10.0);

		Deb_ASSERT(glb.I_norm > 0); /* Just in case */

		/* Calculate CMB intensity */
		if(glb.model.parms.T_cmb > 0) {
			glb.I_cmb = Phys_PlanckFunc(glb.freq, glb.model.parms.T_cmb);
			Deb_ASSERT(glb.I_cmb > 0); /* Just in case */

			/* Normalize CMB */
			glb.I_cmb /= glb.I_norm;
		}
	}

// 	fp=fopen("pops.dat","w");

	for(zp = Zone_GetMinLeaf(root); zp; zp = Zone_AscendTree(zp)) {
		/* Pointer to physical parameters */
		pp = zp->data;

		if((pp->n_H2 > 1e-200) && !zp->children) {
			/* This is a non-empty leaf zone */
			pp->non_empty_leaf = 1;

			if(pp->X_mol > 0) {
				/* This zone contains tracer molecules */
				pp->has_tracer = 1;
/*
 				radius=sqrt(zp->voxel.cen.x[0] * zp->voxel.cen.x[0] + zp->voxel.cen.x[1] * zp->voxel.cen.x[1] + zp->voxel.cen.x[2] * zp->voxel.cen.x[2]);
 				radius=zp->voxel.cen.x[0];
 				fprintf(fp,"%g %g %g %g %g %g %g %g %g %g %g\n",radius,pp->pops[0][0],pp->pops[0][1],pp->pops[0][2],pp->pops[0][3],pp->pops[0][4],pp->pops[0][5],pp->pops[0][6],pp->pops[0][7],pp->pops[0][8],pp->pops[0][9]);
*/
				
			}
			else{
				pp->has_tracer = 0;
			}
		}
		else{
			pp->non_empty_leaf = 0;
			pp->has_tracer = 0;
		}
	}

// 	fclose(fp);

	sts = SpUtil_Threads2(Sp_NTHREAD, InitModelThread);
	//SpUtil_Threads(InitModelThread);

	return sts;
}

/*----------------------------------------------------------------------------*/

static void *InitModelThread(void *tid_p)
{
	size_t tid = *((size_t *)tid_p), zone_id,j,k;
	Zone *root = glb.model.grid, *zp;
	SpPhys *pp;
        
        size_t nrad = glb.model.parms.mol->nrad;
        double freq[nrad];

        

	for(zp = Zone_GetMinLeaf(root), zone_id = 0; zp; zp = Zone_AscendTree(zp), zone_id++) {
		if(zone_id % Sp_NTHREAD == tid) {

			/* Check for thread termination */
			Sp_CHECKTERMTHREAD();

			/* Init zone parameters */
			pp = zp->data;

			if(!glb.coldens){
				/* Set initial pops to either optically thin or LTE */
				if(glb.lte) {

					pp->mol = glb.model.parms.mol;
					for(k = 0; k < Sp_NTHREAD; k++) {
//printf("telsim.c:688 checkpoint zone_id=%d pp=%d pops=%d nlev=%d size=%zd\n",zone_id,pp,pp->pops[k],pp->mol->nlev,sizeof(*(pp->pops[k])) );
						pp->pops[k] = Mem_CALLOC(pp->mol->nlev, pp->pops[k]);
					}


					for(j = 0; j < pp->mol->nlev; j++) {
						pp->pops[0][j] = SpPhys_BoltzPops(pp->mol, j, pp->T_k);
						for(k = 1; k < Sp_NTHREAD; k++) {
							pp->pops[k][j] = pp->pops[0][j];
						}
					}
					
				}

				if(glb.cont) {
// 				if(glb.cont || glb.lte) {
					SpPhys_InitContWindows(pp, &glb.freq, (size_t)1);
				}
				else {

// 				if(!glb.cont) {
 					if(glb.lte){
                                                /* Allocate continuum emission/absorption */
                                                for(size_t i = 0; i < nrad; i++) {
                                                        freq[i] = pp->mol->rad[i]->freq;
                                                }
                                                SpPhys_InitContWindows(pp, freq, nrad);
                                        }
                                        pp->width = SpPhys_CalcLineWidth(pp);
				}

				/* Add dust emission/absorption if T_d > 0 */
				if(pp->T_d > 0) {
					SpPhys_AddContinuum_d(pp, glb.cont, glb.model.parms.gas_to_dust);
				}
				/* Add free-free emission/absorption if T_ff > 0 */

				if(pp->T_ff > 0) {
					SpPhys_AddContinuum_ff(pp, glb.cont);
				}

				/* Set continuum flux */
				if(pp->T_bb > 0) {
					//debug
					SpPhys_SetContinuumIntens_bb(pp, glb.cont, pp->T_bb, glb.I_norm);
					Deb_PRINT("T_bb=%g, F_nu=%g\n", pp->T_bb, pp->cont[0].I_bb);
				}

			}
		}
	}

	pthread_exit(NULL);
}

/*----------------------------------------------------------------------------*/

static int CalcImage(void)
{
	return SpUtil_Threads2(Sp_NTHREAD, CalcImageThread);
}

/*----------------------------------------------------------------------------*/

static void *CalcImageThread(void *tid_p)
{
	size_t tid = *((size_t *)tid_p), pix_id, ix, iy, iv, nsub, isub, jsub, ibox;
	double dx, dy, subix, subiy,
		*I_sub, *I_nu, *tau_sub, *tau_nu,
		*V_sub, *Q_sub, *U_sub, *V_nu, *Q_nu, *U_nu, 
		*sigma2, *sigma2_sub,
		alpha;
	
	alpha=0.15; // polarized efficiency

	/* pix_id is used for distributing work to threads */
	pix_id = 0;
	

	
	for(ix = 0; ix < glb.x.n; ix++) {
		for(iy = 0; iy < glb.y.n; iy++) {
			if(pix_id % Sp_NTHREAD == tid) {

				/* Check for thread termination */
				Sp_CHECKTERMTHREAD();

				/* Determine angular offsets from the pointing center */
				dx = ((int)ix - (int)glb.x.crpix) * glb.x.delt;
				dy = ((int)iy - (int)glb.y.crpix) * glb.y.delt;

				/* nsub is by default 1 */
				nsub = 1;

				/* Check if position is within any of the subres boxes */
				for(ibox = 0; ibox < glb.nsubres; ibox++) {
					if((dx >= glb.subres[ibox].blc_x) && (dx <= glb.subres[ibox].trc_x) &&
					   (dy >= glb.subres[ibox].blc_y) && (dy <= glb.subres[ibox].trc_y)) {
						/* Position within subres box, set nsub to subres[ibox].nsub */
						nsub = glb.subres[ibox].nsub;
						break;
					}
				}

				/* I_nu is the brightness for all channels at pixel (ix, iy) */
				I_nu = Mem_CALLOC(glb.v.n, I_nu);
				if(glb.zeeman) {
					V_nu = Mem_CALLOC(glb.v.n, V_nu);
				}
				else if(glb.cont) {
					Q_nu = Mem_CALLOC(glb.v.n, Q_nu);
					U_nu = Mem_CALLOC(glb.v.n, U_nu);
					sigma2 = Mem_CALLOC(glb.v.n, sigma2);
				}

				/* tau_nu is the total optical depth for all channels at pixel (ix, iy) */
				tau_nu = Mem_CALLOC(glb.v.n, tau_nu);

				/* Loop through sub-resolution positions */
				for(isub = 0; isub < nsub; isub++) {
					for(jsub = 0; jsub < nsub; jsub++) {
						
						/* Determine sub-resolution angular offsets from the pointing center */
// 						subix = (double)ix + ((double)isub) / (double)nsub;
// 						subiy = (double)iy + ((double)jsub) / (double)nsub;
						subix = (double)ix + (2.*(double)isub-(double)nsub+1.) / ((double)(2*nsub));
 						subiy = (double)iy + (2.*(double)jsub-(double)nsub+1.) / ((double)(2*nsub));
						dx = (subix - (double)glb.x.crpix) * glb.x.delt;
						dy = (subiy - (double)glb.y.crpix) * glb.y.delt;
						/* I_sub is the brightness for all channels at each sub-resolution */
						I_sub = Mem_CALLOC(glb.v.n, I_sub);
						tau_sub = Mem_CALLOC(glb.v.n, tau_sub);
						if(glb.zeeman) {
							V_sub = Mem_CALLOC(glb.v.n, V_sub);
						}
						else if(glb.cont) {
							Q_sub = Mem_CALLOC(glb.v.n, Q_sub);
							U_sub = Mem_CALLOC(glb.v.n, U_sub);
							sigma2_sub = Mem_CALLOC(glb.v.n, sigma2_sub);
						}

						/* Column density tracer */
						if(glb.coldens){
							ColumnDensityTracer(dx,dy,I_sub);
						}
						/* intensity tracer */
						else{
							/* Calculate radiative transfer for this sub-los */
							RadiativeXfer(dx, dy, I_sub, V_sub, Q_sub, U_sub, sigma2_sub, tau_sub);
							
						}

						/* Add I_sub to I_nu */
						for(iv = 0; iv < glb.v.n; iv++) {
							

							if (glb.zeeman){
								V_nu[iv] += V_sub[iv];
							}
							else{
								I_nu[iv] += I_sub[iv];
								if(glb.cont){
									Q_nu[iv] += Q_sub[iv];
									U_nu[iv] += U_sub[iv];
									sigma2[iv] += sigma2_sub[iv];
								}
							}
							
							tau_nu[iv] += tau_sub[iv];
						}

						/* Cleanup */
						free(I_sub);
                                                free(tau_sub);
						if (glb.zeeman)
                                                        free(V_sub);
						else if(glb.cont) {
                                                        free(Q_sub);
                                                        free(U_sub);
                                                        free(sigma2_sub);
                                                }

					}
				}
			
				/* Save averaged I_nu to map */
				for(iv = 0; iv < glb.v.n; iv++) {
					
					if (glb.zeeman){
// 						MirImg_PIXEL(*glb.stokesv, iv, ix, iy) = V_nu[iv] / (double)(nsub * nsub);
						MirImg_PIXEL(*glb.xyv_img, iv, ix, iy) = V_nu[iv] / (double)(nsub * nsub);
					}
					else{
						MirImg_PIXEL(*glb.xyv_img, iv, ix, iy) = I_nu[iv] / (double)(nsub * nsub);
						if(glb.cont){
							MirImg_PIXEL(*glb.stokesq, iv, ix, iy) = alpha*Q_nu[iv] / (double)(nsub * nsub);
							MirImg_PIXEL(*glb.stokesu, iv, ix, iy) = alpha*U_nu[iv] / (double)(nsub * nsub);
							MirImg_PIXEL(*glb.sigma2, iv, ix, iy) = 0.5*alpha*sigma2[iv] / (double)(nsub * nsub);
						}
					}

					if(glb.tau_img)
						MirImg_PIXEL(*glb.tau_img, iv, ix, iy) = tau_nu[iv] / (double)(nsub * nsub);
				}


				/* Cleanup */
				free(I_nu);
                                free(tau_nu);
				if (glb.zeeman)
                                        free(V_nu);
				else if(glb.cont) {
                                        free(Q_nu);
                                        free(U_nu);
                                        free(sigma2);
                                }
				
			}
			/* Increment pix_id */
			pix_id += 1;
		}
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/

static void RadiativeXfer(double dx, double dy, double *I_nu, double *V_nu,double *Q_nu, double *U_nu, double *sigma2, double *tau_nu)
{
	GeRay ray;
	double dv, theta, phi, t, vfac, vfac2, j_nu, k_nu, S_nu, dtau_nu, tempj_nu, tempk_nu;
	SpPhys *pp;
	size_t iv, side,i,j;
	Zone *zp, *root = glb.model.grid;
	
	GeVec3_d z,n,e;
	double nproduct,zproduct,eproduct,psi,cosgammasquare,costheta,temp;
	double Vr,Vt,Vp,thetac,phic;
	GeVec3_d B;
	double d23 = 2./3. ;
	
	double dnu, deltav;
	double g = 2.18/1.4;

	/* Reset ray */
	Mem_BZERO(&ray);

	/* Init ray position to <dist, 0, 0> */
	GeRay_E(ray, 0) = glb.dist / Sp_LENFAC;
	GeRay_E(ray, 1) = 0;
	GeRay_E(ray, 2) = 0;

	/* Set direction of ray according to pixel position:
	 *   theta = PI/2 + dy
	 *   phi = -dx
	 */
	phi = M_PI + dx;   //
	theta = (M_PI / 2.0) - dy;

	/*  the coordinate
	 *                          north   eastern
	 *                             ^   A 
	 *                             |  / 
	 *                             | /
	 *                             |/
	 * tracing direction D <-------|-------> light traversing direction z
	 */


	/* Convert to Cartesian coordinates */
	GeRay_D(ray, 0) = sin(theta) * cos(phi);
	GeRay_D(ray, 1) = sin(theta) * sin(phi);
	GeRay_D(ray, 2) = cos(theta);
	
	/* line of sight coordinate */
	if (glb.cont || glb.zeeman){
		
		z.x[0]=-GeRay_D(ray, 0);
		z.x[1]=-GeRay_D(ray, 1);
		z.x[2]=-GeRay_D(ray, 2);
		n.x[0]=-cos(theta)*cos(phi);
		n.x[1]=-cos(theta)*sin(phi);
		n.x[2]= sin(theta);
		e.x[0]= sin(phi);
		e.x[1]= -cos(phi);
		e.x[2]= 0.0;
		
		z = GeVec3_Rotate_z(&z, -glb.rotate[2]);
		n = GeVec3_Rotate_z(&n, -glb.rotate[2]);
		e = GeVec3_Rotate_z(&e, -glb.rotate[2]);
		z = GeVec3_Rotate_y(&z, -glb.rotate[1]);
		n = GeVec3_Rotate_y(&n, -glb.rotate[1]);
		e = GeVec3_Rotate_y(&e, -glb.rotate[1]);
		z = GeVec3_Rotate_x(&z, -glb.rotate[0]);
		n = GeVec3_Rotate_x(&n, -glb.rotate[0]);
		e = GeVec3_Rotate_x(&e, -glb.rotate[0]);
	
	}

	/* Rotate ray:
	 * Since what we REALLY want to rotate is the model and that the rays
	 * are pointed towards the model, rays should be rotated in the negative
	 * direction, and in the opposite order of what we would've done to
	 * rotate the model. */
	ray = GeRay_Rotate(&ray, 2, -glb.rotate[2]);
	ray = GeRay_Rotate(&ray, 1, -glb.rotate[1]);
	ray = GeRay_Rotate(&ray, 0, -glb.rotate[0]);

	/* Coordinate-dependent offset */
	switch(root->voxel.geom) {
		case GEOM_SPH1D:
			break;
			
		case GEOM_SPH3D:
			break;

		case GEOM_REC3D:
			ray.e = GeVec3_Add(&ray.e, &root->voxel.cen);
			break;
			
		case GEOM_CYL3D:
			break;

		default: /* Shouldn't reach here */
			Deb_ASSERT(0);
	}

	/* Reset tau for all channels */
	Mem_BZERO2(tau_nu, glb.v.n);

	/* Shoot ray at model and see what happens! */
	if(GeRay_IntersectVoxel(&ray, &root->voxel, &t, &side)) {
		/* Calculate intersection */
		ray = GeRay_Inc(&ray, t);

		/* Locate starting leaf zone according to intersection */
		zp = Zone_GetLeaf(root, side, &ray.e, &ray);

		/* Keep going until there's no next zone to traverse to */
		while(zp) {
			/* Calculate path to next boundary */
			GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &side);
// 			Deb_PRINT("checkpoint: GeRay_TraverseVoxel\n");

			/* Pointer to physical parameters associated with this zone */
			pp = zp->data;

			/* Do radiative transfer only if gas is present in this zone */
			if(pp->non_empty_leaf) {
				/* Do calculations on all channels at this pixel. Try to minimize the 
				 * amount of operations in this loop, since everything here is repeated 
				 * for ALL channels, and can significantly increase computation time.
				 */
				if(glb.zeeman){
					B = SpPhys_GetBfac(&ray, t, zp, 0);
					zproduct = GeVec3_DotProd(&z,&B);
					costheta = zproduct/GeVec3_Mag(&B);
					dnu = 1.4e6 * g * GeVec3_Mag(&B);
					deltav = dnu * PHYS_CONST_MKS_LIGHTC / glb.freq;
					//printf("delta V = %e %e %e %e\n",deltav,dnu,PHYS_CONST_MKS_LIGHTC,glb.freq);
				}
				else if (glb.cont){
					B = SpPhys_GetBfac(&ray, t, zp, 0);
					
					nproduct = GeVec3_DotProd(&n,&B);
					eproduct = GeVec3_DotProd(&e,&B);
					zproduct = GeVec3_DotProd(&z,&B);
					
// 					if (nproduct==0.0){
// 						if(-eproduct>0) psi=0.5*M_PI;
// 						else psi=-0.5*M_PI;
// 					}
// 					else
// 						psi=atan(-eproduct/nproduct);
// 					if (nproduct<0)
// 						psi=psi+M_PI;
// 					else if (psi<0)
// 						psi=psi+2*M_PI;

					psi=atan2(-eproduct,nproduct); // psi is the angle between the projected B-field on p.o.s. and the north of the image
					
					cosgammasquare=1.0-zproduct*zproduct/( GeVec3_X(B, 0)*GeVec3_X(B, 0) + GeVec3_X(B, 1)*GeVec3_X(B, 1) + GeVec3_X(B, 2)*GeVec3_X(B, 2) ); // gamma ia the angle bettwen B-field an the plane of sky
				}
				

				for(iv = 0; iv < glb.v.n; iv++) {
					/* Calculate velocity associated with this channel */
					dv = ((double)iv - glb.v.crpix) * glb.v.delt;

					/* Reset emission and absorption coeffs */
					j_nu = 0;
					k_nu = 0;

					if(!glb.cont && pp->has_tracer) {
						if(glb.zeeman){
							
// 							vfac = SpPhys_GetVfac(&ray, t, dv+deltav, zp, 0); 
// 							SpPhys_GetMoljk((size_t)0, pp, j, vfac, &tempj_nu, &tempk_nu);
// 							
// 							j_nu += 0.5*tempj_nu;
// 							k_nu += 0.5*tempk_nu;
// 							
// 							vfac = SpPhys_GetVfac(&ray, t, dv-deltav, zp, 0);  
// 							SpPhys_GetMoljk((size_t)0, pp, j, vfac, &tempj_nu, &tempk_nu);
// 							
// 							j_nu -= 0.5*tempj_nu;
// 							k_nu += 0.5*tempk_nu;


							vfac = SpPhys_GetVfac(&ray, t, dv+deltav, zp, 0);
							vfac2 = SpPhys_GetVfac(&ray, t, dv-deltav, zp, 0);
							vfac = vfac-vfac2;
							SpPhys_GetMoljk((size_t)0, pp, j, vfac, &tempj_nu, &tempk_nu);

							j_nu = 0.5*tempj_nu;
							k_nu += tempk_nu;
							//if(deltav>1.0)printf("dv= %e vfac = %e j_nu = %e\n",dv,vfac,j_nu);
							
						
						}
						else{ // discard zeeman splitting
							/* Calculate velocity line profile factor for this channel:
							 * This version averages over the line profile in steps
							 * of the local line width -- very time consuming! */
							vfac = SpPhys_GetVfac(&ray, t, dv, zp, 0);

							/* Calculate molecular line emission and absorption coefficients */
							//SpPhys_GetMoljk((size_t)0, pp, glb.line, vfac, &j_nu, &k_nu);
						
							if(glb.overlap){
								i=glb.line;
								for(j = 0; j < NRAD; j++) {
									if(OVERLAP(i,j)){
										if(i==j){
											/* Calculate molecular line emission and absorption coefficients */
											SpPhys_GetMoljk((size_t)0, pp, j, vfac, &tempj_nu, &tempk_nu);
										}
										else{
											/* Calculate velocity line profile factor */
											vfac2 = pp->has_tracer ? SpPhys_GetVfac(&ray, t, dv-RELVEL(i,j), zp, 0) : 0.0;
											/* Calculate molecular line emission and absorption coefficients */
											SpPhys_GetMoljk((size_t)0, pp, j, vfac2, &tempj_nu, &tempk_nu);
										}
										j_nu += tempj_nu;
										k_nu += tempk_nu;
									}
								}
							}
							else{
								/* Calculate molecular line emission and absorption coefficients */
								SpPhys_GetMoljk((size_t)0, pp, glb.line, vfac, &j_nu, &k_nu);
							}
						}
						
						
					}
					
                                        /* Add continuum emission/absorption */
					if(!glb.zeeman) j_nu += pp->cont[glb.line].j;
					k_nu += pp->cont[glb.line].k;

					/* Calculate source function and optical depth if
					 * absorption is NOT zero */
					if(fabs(k_nu) > 0.0) {
						//S_nu = j_nu / ( k_nu * glb.I_norm );
						S_nu = j_nu / ( k_nu * glb.I_norm );
						dtau_nu = k_nu * t * Sp_LENFAC;
					}
					else {
						S_nu = dtau_nu = 0.0;
					}


					/* Calculate intensity contributed by this step */
					//debug
					temp = S_nu * (1.0 - exp(-dtau_nu)) * exp(-tau_nu[iv]);
					
					//I_nu[iv] += (S_nu * (1.0 - exp(-dtau_nu)) + pp->cont[glb.line].I_bb) * exp(-tau_nu[iv]);
					
					if(glb.zeeman){
						V_nu[iv] += temp*costheta;
					}
					else{
						I_nu[iv] += temp;
						if (glb.cont){
							Q_nu[iv] += temp*cos(2.0 * psi) * cosgammasquare;
							U_nu[iv] += temp*sin(2.0 * psi) * cosgammasquare;
							sigma2[iv] += temp * (cosgammasquare - d23);
						}
					}
					/* Accumulate total optical depth for this channel (must be done
					 * AFTER calculation of intensity!) */
					tau_nu[iv] += dtau_nu;
				}
			}
// 			Deb_PRINT("checkpoint: non_empty_leaf\n");
			/* Calculate next position */
			ray = GeRay_Inc(&ray, t);
// 			Deb_PRINT("checkpoint: GeRay_Inc\n");
			/* Get next zone to traverse to */
			zp = Zone_GetNext(zp, &side, &ray);
			#if 0
			if(zp){
				Deb_PRINT("zone index: %d %d %d\n",GeVec3_X(zp->index,0),GeVec3_X(zp->index,1),GeVec3_X(zp->index,2));
				Deb_PRINT("traveling distance = %e\n",t);
			}
			else
				Deb_PRINT("shoot out!\n");
			#endif
			
			#if 0
			if(I_nu[0] != I_nu[0] && zp ) { // check if intensity become NaN
				Deb_PRINT("zone index: %d %d %d\n",GeVec3_X(zp->index,0),GeVec3_X(zp->index,1),GeVec3_X(zp->index,2));
				Deb_PRINT("dnesity=%g\n",pp->n_H2);
				if(pp->non_empty_leaf) Deb_PRINT("non empty leaf\n");
				
			}
			else if(zp == NULL)
				Deb_PRINT("shoot out!\n");
			#endif
		}
	}

	/* Add CMB to all channels -- this is done even if the ray misses the source */
	

	for(iv = 0; iv < glb.v.n; iv++)
		I_nu[iv] += glb.I_cmb * exp(-tau_nu[iv]);

	return;
}





static void ColumnDensityTracer(double dx, double dy, double *I_nu)
{
	GeRay ray;
	double dv, theta, phi, t, vfac;
	SpPhys *pp;
	size_t side;
	Zone *zp, *root = glb.model.grid;
	
	GeVec3_d z;
// 	GeVec3_d B;
// 	GeVec3_d V;
	double zproduct;

	/* Reset ray */
	Mem_BZERO(&ray);

	/* Init ray position to <dist, 0, 0> */
	GeRay_E(ray, 0) = glb.dist / Sp_LENFAC;
	GeRay_E(ray, 1) = 0;
	GeRay_E(ray, 2) = 0;

	/* Set direction of ray according to pixel position:
	 *   theta = PI/2 + dy
	 *   phi = -dx
	 */
	phi = -dx;
	theta = (M_PI / 2.0) + dy;

	/* Convert to Cartesian coordinates */
	GeRay_D(ray, 0) = sin(theta) * cos(phi);
	GeRay_D(ray, 1) = sin(theta) * sin(phi);
	GeRay_D(ray, 2) = cos(theta);
	
	
			
	/* Rotate ray:
	 * Since what we REALLY want to rotate is the model and that the rays
	 * are pointed towards the model, rays should be rotated in the negative
	 * direction, and in the opposite order of what we would've done to
	 * rotate the model. */
	ray = GeRay_Rotate(&ray, 2, -glb.rotate[2]);
	ray = GeRay_Rotate(&ray, 1, -glb.rotate[1]);
	ray = GeRay_Rotate(&ray, 0, -glb.rotate[0]);

	/* Coordinate-dependent offset */
	switch(root->voxel.geom) {
		case GEOM_SPH1D:
			break;
			
		case GEOM_SPH3D:
			break;

		case GEOM_REC3D:
			ray.e = GeVec3_Add(&ray.e, &root->voxel.cen);
			break;
			
		case GEOM_CYL3D:
			break;

		default: /* Shouldn't reach here */
			Deb_ASSERT(0);
	}


	z.x[0]=-GeRay_D(ray, 0);
	z.x[1]=-GeRay_D(ray, 1);
	z.x[2]=-GeRay_D(ray, 2);


	

	/* Shoot ray at model and see what happens! */
	if(GeRay_IntersectVoxel(&ray, &root->voxel, &t, &side)) {

		/* Calculate intersection */
		ray = GeRay_Inc(&ray, t);

		/* Locate starting leaf zone according to intersection */
		zp = Zone_GetLeaf(root, side, &ray.e, &ray);

		/* Keep going until there's no next zone to traverse to */
		while(zp) {
			/* Calculate path to next boundary */
			GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &side);


			/* Pointer to physical parameters associated with this zone */
			pp = zp->data;

			/* Do radiative transfer only if gas is present in this zone */
			if(pp->non_empty_leaf) {
				/* Do calculations on all channels at this pixel. Try to minimize the 
				 * amount of operations in this loop, since everything here is repeated 
				 * for ALL channels, and can significantly increase computation time.
				 */ 
				
				
// 				B = SpPhys_GetBfac(&ray, t, zp, 0);
// 				zproduct = GeVec3_DotProd(&z,&B);
//				V = SpPhys_GetVfac2(&ray, t, zp, 0);
// 				zproduct = GeVec3_DotProd(&z,&V);
				
				
                      
// 				I_nu[0] += pp->n_H2 * pp->X_mol * t * Sp_LENFAC * zproduct;
				I_nu[0] += pp->n_H2 * pp->X_mol * t * Sp_LENFAC;
			}

			/* Calculate next position */
			ray = GeRay_Inc(&ray, t);

			/* Get next zone to traverse to */
			zp = Zone_GetNext(zp, &side, &ray);


		}
	}

	return;
}


static void visualization(void)
{
	Zone *zp, *zp1, *zp2, *zp3, *root = glb.model.grid;
	FILE *fp, *fp2;
	size_t ix,iy,iz, izone;
	size_t ix1,iy1,iz1, izone1;
	size_t ix2,iy2,iz2, izone2;
	char filename[22];
	double Tex;
	SpPhys *pp;
	static const double c = PHYS_CONST_MKS_LIGHTC,
		h = PHYS_CONST_MKS_PLANCKH,
		k = PHYS_CONST_MKS_BOLTZK;
	MolTrRad *trans;
	double n_u, n_l, g_u, g_l, E_u, E_l;
	size_t up, lo;
	
	
	fp=fopen("vis.pvd","w");
	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n");
	fprintf(fp,"  <Collection>\n");
	mkdir("multiblock", S_IRWXU | S_IRWXG | S_IRWXO);
	
	zp1 = root;
	for(ix1=0; ix1 < zp1->naxes.x[0]; ix1++){
	for(iy1=0; iy1 < zp1->naxes.x[1]; iy1++){
	for(iz1=0; iz1 < zp1->naxes.x[2]; iz1++){
		izone1 = iz1 + iy1*(zp1->naxes.x[2]) + ix1*(zp1->naxes.x[2])*(zp1->naxes.x[1]);
		zp2 = zp1->children[izone1];
		for(ix2=0; ix2 < zp2->naxes.x[0]; ix2++){
		for(iy2=0; iy2 < zp2->naxes.x[1]; iy2++){
		for(iz2=0; iz2 < zp2->naxes.x[2]; iz2++){
			izone2 = iz2 + iy2*(zp2->naxes.x[2]) + ix2*(zp2->naxes.x[2])*(zp2->naxes.x[1]);
			zp3 = zp2->children[izone2];
			
			sprintf(filename,"multiblock/block%2.2zu%2.2zu.vtr",izone1,izone2);
			fprintf(fp,"    <DataSet group=\"%zu\" dataset=\"0\" file=\"%22s\"/>\n",izone1*8+izone2,filename);
			
			fp2=fopen(filename,"w");
			fprintf(fp2,"<?xml version=\"1.0\"?>\n");
			fprintf(fp2,"<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
			fprintf(fp2,"<RectilinearGrid WholeExtent=\"%d %d %d %d %d %d\">\n",0,64,0,64,0,32);
			fprintf(fp2,"    <Piece Extent=\"%d %d %d %d %d %d\">\n",0,64,0,64,0,32);
			fprintf(fp2,"      <Coordinates>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"X_COORDINATES\" NumberOfComponents=\"1\">\n");
			fprintf(fp2," %12.6e",zp3->voxel.min.x[0]);
			for(ix=0; ix < zp3->naxes.x[0]; ix++){
				izone = ix*(zp3->naxes.x[2])*(zp3->naxes.x[1]);
				zp = zp3->children[izone];
				fprintf(fp2," %12.6e",zp->voxel.max.x[0]);
			}
			fprintf(fp2,"\n        </DataArray>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"Y_COORDINATES\" NumberOfComponents=\"1\">\n");
			fprintf(fp2," %12.6e",zp3->voxel.min.x[1]);
			for(iy=0; iy < zp3->naxes.x[1]; iy++){
				izone = iy*(zp3->naxes.x[2]);
				zp = zp3->children[izone];
				fprintf(fp2," %12.6e",zp->voxel.max.x[1]);
			}
			fprintf(fp2,"\n        </DataArray>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"Z_COORDINATES\" NumberOfComponents=\"1\">\n");
			fprintf(fp2," %12.6e",zp3->voxel.min.x[2]);
			for(iz=0; iz < zp3->naxes.x[2]; iz++){
				izone = iz;
				zp = zp3->children[izone];
				fprintf(fp2," %12.6e",zp->voxel.max.x[2]);
			}
			fprintf(fp2,"\n        </DataArray>\n");
			fprintf(fp2,"      </Coordinates>\n");
			fprintf(fp2,"      <CellData>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"Tex*\" NumberOfComponents=\"1\">\n");
			
			for(iz=0; iz < zp3->naxes.x[2]; iz++){
			for(iy=0; iy < zp3->naxes.x[1]; iy++){
			for(ix=0; ix < zp3->naxes.x[0]; ix++){
				izone = iz + iy*(zp3->naxes.x[2]) + ix*(zp3->naxes.x[2])*(zp3->naxes.x[1]);
				zp = zp3->children[izone];
				pp = zp->data;
				if(pp->X_mol == 0.0){
					fprintf(fp2,"0.0\n");
				}
				else{
					trans = pp->mol->rad[glb.line];
					up = trans->up;
					lo = trans->lo;
					n_u = pp->pops[0][up];
					n_l = pp->pops[0][lo];
					E_u = pp->mol->lev[up]->E;
					E_l = pp->mol->lev[lo]->E;
					g_u = pp->mol->lev[up]->g;
					g_l = pp->mol->lev[lo]->g;
					Tex = (E_l-E_u)/(k*log((n_u*g_l)/(n_l*g_u))); 
					fprintf(fp2,"%12.6e\n",Tex/(pp->T_k));
				}
			}
			}
			}
			fprintf(fp2,"        </DataArray>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"Tk\" NumberOfComponents=\"1\">\n");
			
			for(iz=0; iz < zp3->naxes.x[2]; iz++){
			for(iy=0; iy < zp3->naxes.x[1]; iy++){
			for(ix=0; ix < zp3->naxes.x[0]; ix++){
				izone = iz + iy*(zp3->naxes.x[2]) + ix*(zp3->naxes.x[2])*(zp3->naxes.x[1]);
				zp = zp3->children[izone];
				pp = zp->data;
				fprintf(fp2,"%12.6e\n",pp->T_k);
			}
			}
			}
			fprintf(fp2,"        </DataArray>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"N_h2\" NumberOfComponents=\"1\">\n");
			
			for(iz=0; iz < zp3->naxes.x[2]; iz++){
			for(iy=0; iy < zp3->naxes.x[1]; iy++){
			for(ix=0; ix < zp3->naxes.x[0]; ix++){
				izone = iz + iy*(zp3->naxes.x[2]) + ix*(zp3->naxes.x[2])*(zp3->naxes.x[1]);
				zp = zp3->children[izone];
				pp = zp->data;
				fprintf(fp2,"%12.6e\n",pp->n_H2);
			}
			}
			}
			fprintf(fp2,"        </DataArray>\n");
			fprintf(fp2,"        <DataArray type=\"Float32\" Name=\"contribution\" NumberOfComponents=\"1\">\n");
			
			for(iz=0; iz < zp3->naxes.x[2]; iz++){
			for(iy=0; iy < zp3->naxes.x[1]; iy++){
			for(ix=0; ix < zp3->naxes.x[0]; ix++){
				izone = iz + iy*(zp3->naxes.x[2]) + ix*(zp3->naxes.x[2])*(zp3->naxes.x[1]);
				zp = zp3->children[izone];
				pp = zp->data;
				if(pp->X_mol == 0.0){
					fprintf(fp2,"0.0\n");
				}
				else{
					trans = pp->mol->rad[glb.line];
					up = trans->up;
					lo = trans->lo;
					n_u = pp->pops[0][up];
					n_l = pp->pops[0][lo];
					E_u = pp->mol->lev[up]->E;
					E_l = pp->mol->lev[lo]->E;
					g_u = pp->mol->lev[up]->g;
					g_l = pp->mol->lev[lo]->g;
					Tex = (E_l-E_u)/(k*log((n_u*g_l)/(n_l*g_u))); 
					fprintf(fp2,"%12.6e\n", Tex * pp->n_H2 * pp->X_mol);
				}
			}
			}
			}
			fprintf(fp2,"        </DataArray>\n");
			fprintf(fp2,"      </CellData>\n");
			fprintf(fp2,"    </Piece>\n");
			fprintf(fp2,"  </RectilinearGrid>\n");
			fprintf(fp2,"</VTKFile>");
			fclose(fp2);
			
		}
		}
		}
	}
	}
	}
		
		fprintf(fp,"  </Collection>\n");
		fprintf(fp,"</VTKFile>");
		
		
		fclose(fp);
}

static void call_fits_output( const char *FileName, const char *bunit, const double scale, const int Stokes){
	size_t ix,iy,iz,iv;

	/* output FITS file  */
	int status = 0;
	fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
	
	/* createing new FITS file  */
 	fits_create_file(&fptr, FileName, &status);   /* create new file */
	Deb_ASSERT(status == 0);
        
	long naxis;
	if (Stokes) naxis = 4;
	else naxis =3;
	
	long naxes[naxis];
	long  fpixel[3];
	naxes[0] = (long) glb.x.n;
	naxes[1] = (long) glb.y.n;
	naxes[2] = (long) glb.v.n;
	if(Stokes) naxes[3] = (long) 3;
	long nelements = 1;
	for (int i = 0; i < naxis; i++)
		nelements *= naxes[i];
	
	fpixel[0] = fpixel[1] = fpixel[2] = 1;
	
	double 
	crpix1 = (double) naxes[0]/2. + 0.5,
	crpix2 = (double) naxes[1]/2. + 0.5,
	crpix3 = (double) naxes[2]/2. + 0.5,
	crpix4 = 1.,
	cdelt1 = 180. / M_PI * glb.x.delt,
	cdelt2 = 180. / M_PI * glb.y.delt,
	cdelt3 = 180. / M_PI * glb.v.delt,
	cdelt4 = 1.,
	crval = 0.0,
	crval4 = 1.,
	beam = 0.0,
	bscale = 1.,
	bzero = 0.,
	restfreq = glb.freq;
	
	double *array = malloc(nelements*sizeof(*array));
	
	fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
	Deb_ASSERT(status == 0);
	    /* Write a keyword; must pass the ADDRESS of the value */
	 
	fits_write_key(fptr,TSTRING,"CTYPE1","RA---SIN","Type of first axis",&status);
	fits_write_key(fptr,TDOUBLE,"CRPIX1",&crpix1,"Reference of first axis",&status);
	fits_write_key(fptr,TDOUBLE,"CDELT1",&cdelt1,"Increment value of first axis",&status);
	fits_write_key(fptr,TDOUBLE,"CRVAL1",&crval,"Offset of first axis ",&status);
	fits_write_key(fptr,TSTRING,"CTYPE2","DEC--SIN","Type of second axis",&status);
	fits_write_key(fptr,TDOUBLE,"CRPIX2",&crpix2,"Reference of second axis",&status);
	fits_write_key(fptr,TDOUBLE,"CDELT2",&cdelt2,"Increment value of second axis",&status);
	fits_write_key(fptr,TDOUBLE,"CRVAL2",&crval,"Offset of second axis ",&status);
	fits_write_key(fptr,TSTRING,"CTYPE3","VELO-LSR","Type of third axis",&status);
	fits_write_key(fptr,TDOUBLE,"CRPIX3",&crpix3,"Reference of third axis",&status);
	fits_write_key(fptr,TDOUBLE,"CDELT3",&cdelt3,"Increment value of third axis",&status);
	fits_write_key(fptr,TDOUBLE,"CRVAL3",&crval,"Offset of third axis ",&status);
	if(Stokes){
		fits_write_key(fptr,TSTRING,"CTYPE4","STOKES","Type of third axis",&status);
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
	fits_write_key(fptr,TSTRING,"CELLSCAL","1/F     ","",&status);
        Deb_ASSERT(status == 0);
	
	
	if(Stokes) {
		for (iv = 0; iv < naxes[2]; iv++) for (iy = 0; iy < naxes[1]; iy++) for (ix = 0; ix < naxes[0]; ix++) {
			int idx = ( ( 0 * naxes[2] + iv) * naxes[1] + iy) * naxes[0] + ix;
			array[idx] = scale * MirImg_PIXEL(*glb.xyv_img, iv, ix, iy);
		}
		for (iv = 0; iv < naxes[2]; iv++) for (iy = 0; iy < naxes[1]; iy++) for (ix = 0; ix < naxes[0]; ix++) {
			
			int idx = ( ( 1 * naxes[2] + iv) * naxes[1] + iy) * naxes[0] + ix;
			array[idx] = scale * MirImg_PIXEL(*glb.stokesq, iv, ix, iy);
		}
		for (iv = 0; iv < naxes[2]; iv++) for (iy = 0; iy < naxes[1]; iy++) for (ix = 0; ix < naxes[0]; ix++) {
			int idx = ( ( 2 * naxes[2] + iv) * naxes[1] + iy) * naxes[0] + ix;
			array[idx] = scale * MirImg_PIXEL(*glb.stokesu, iv, ix, iy);
		}
	}
	else{
		for (iv = 0; iv < naxes[2]; iv++)
		  for (iy = 0; iy < naxes[1]; iy++)
		    for (ix = 0; ix < naxes[0]; ix++) {
			int idx = (iv * naxes[1] + iy) * naxes[0] + ix;
			array[idx] = scale * MirImg_PIXEL(*glb.xyv_img, iv, ix, iy);
		}
	}
	
	Deb_ASSERT(status == 0);
	fits_write_pix(fptr, TDOUBLE, fpixel, nelements, array, &status);
	fits_close_file(fptr, &status);            /* close the file */
	Deb_ASSERT(status == 0);
	
	free(array);
	
	
	
}










