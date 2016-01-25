#include "sparx.h"

/* Global parameter struct */
static struct glb {
	size_t nzone;
	SpFile *outf;
	PyObject *r, *nh2, *xmol, *x_pH2, *x_oH2, *x_e, *x_H, *x_He, *tk, *vt, *vr;
	SpModel model;
	double rmax;
	Kappa *dust;
} glb;

int SpTask_AsciiGrid(void)
{
	int status = 0;
	size_t i;
	GeVec3_d min, max;
	GeVec3_s ndiv;
	Zone *zone;
	SpPhys *pp;
	double ra, rb;

	/* Inputs --------------------------------------------------------------------*/
	Mem_BZERO(&glb);
	glb.r = SpInp_GetKey_obj("r");
	glb.nh2 = SpInp_GetKey_obj("nh2");
	glb.xmol = SpInp_GetKey_obj("xmol");
	glb.x_pH2 = SpInp_GetKey_obj("x_pH2");
	glb.x_oH2 = SpInp_GetKey_obj("x_oH2");
	glb.x_e = SpInp_GetKey_obj("x_e");
	glb.x_H = SpInp_GetKey_obj("x_H");
	glb.x_He = SpInp_GetKey_obj("x_He");
	glb.tk = SpInp_GetKey_obj("tk");
	glb.vt = SpInp_GetKey_obj("vt");
	glb.vr = SpInp_GetKey_obj("vr");
	glb.rmax = SpInp_GetKey_dbl("rmax");
	glb.nzone = SpInp_GetKey_size_t("nzone");
	glb.model.parms.T_cmb = SpInp_GetKey_dbl("tcmb");
	glb.model.parms.gas_to_dust = SpInp_GetKey_dbl("gas2dust");

	/* outf */
	if(!status && !(glb.outf = SpInp_GetKey_spfile("out", Sp_NEW)))
		status = 1;
	/*----------------------------------------------------------------------------*/

	/* Generate grid -------------------------------------------------------------*/
	if(!status) {
		/* Generate grid */
		min = GeVec3_d_Init(0.0, 0.0, 0.0);
		max = GeVec3_d_Init(glb.rmax, 0.0, 0.0);
		ndiv = GeVec3_s_Init(glb.nzone, (size_t)1, (size_t)1);
		SpModel_InitGrid(&glb.model, GEOM_SPH1D, min, max, ndiv);

		/* Fill in values */
		for(i = 0; i < glb.nzone; i++) {
			zone = Zone_CHILD2(glb.model.grid, i, 0, 0);
			pp = zone->data;

			/* voxel */
			ra = (i == 0 ? 0.0 : Sp_PYDBL(Sp_PYLST(glb.r, i - 1)));
			rb = Sp_PYDBL(Sp_PYLST(glb.r, i));

			if(!(rb > ra)) {
				Deb_PRINT("rb=%.5e, ra=%.5e\n", rb, ra);
				assert(0);
			}

			zone->voxel = GeVox_Init(GEOM_SPH1D, ra, 0.0, 0.0, rb, PI, TWOPI);
			pp->n_H2 = Sp_PYDBL(Sp_PYLST(glb.nh2, i)); /* m^-3 */
			pp->X_mol = Sp_PYDBL(Sp_PYLST(glb.xmol, i)); /* fraction */
			pp->X_pH2 = Sp_PYDBL(Sp_PYLST(glb.x_pH2, i)); /* fraction */
			pp->X_oH2 = Sp_PYDBL(Sp_PYLST(glb.x_oH2, i)); /* fraction */
			pp->X_e = Sp_PYDBL(Sp_PYLST(glb.x_e, i)); /* fraction */
			pp->X_H = Sp_PYDBL(Sp_PYLST(glb.x_H, i)); /* fraction */
			pp->X_He = Sp_PYDBL(Sp_PYLST(glb.x_He, i)); /* fraction */
			pp->T_k = Sp_PYDBL(Sp_PYLST(glb.tk, i)); /* K */
			pp->V_t = Sp_PYDBL(Sp_PYLST(glb.vt, i)); /* m/s */
			GeVec3_X(pp->v_cen, 0) = Sp_PYDBL(Sp_PYLST(glb.vr, i)); /* m/s */
		}

		/* Write model to file */
		Sp_PRINT("Wrote model to `%s'\n", glb.outf->name);
		status = SpIO_FwriteModel(glb.outf, glb.model);
	}
	/*----------------------------------------------------------------------------*/

	/* Cleanup -------------------------------------------------------------------*/
	Py_DECREF(glb.r);
	Py_DECREF(glb.nh2);
	Py_DECREF(glb.xmol);
	Py_DECREF(glb.x_pH2);
	Py_DECREF(glb.x_oH2);
	Py_DECREF(glb.x_e);
	Py_DECREF(glb.x_H);
	Py_DECREF(glb.x_He);
	Py_DECREF(glb.tk);
	Py_DECREF(glb.vt);
	Py_DECREF(glb.vr);
	SpModel_Cleanup(glb.model);
	SpIO_CloseFile(glb.outf);
	/*----------------------------------------------------------------------------*/

	return status;
}

















