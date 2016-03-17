#include "sparx.h"

/* Global parameter struct */
static struct glb {
	size_t nzone;
	SpFile *outf;
	PyObject *pygrid;
	SpModel model;
	Kappa *dust;
} glb;

int SpTask_PyGrid(void)
{
	int sts = 0, geom=0;
	double gas_to_dust=0., T_cmb=0.;
	size_t i, j, k, l, m;
	GeVec3_d min, max;
	GeVec3_s ndiv;
	PyObject *obj=0, *n_H2=0, *X_mol=0, *X_pH2=0, *X_oH2=0, *X_e=0, *X_H=0, *X_He=0, *T_k=0, *T_d=0, *T_ff=0, *T_bb=0, *V_t=0, *V_i=0, *V_j=0, *V_k=0,
		*max_i=0, *max_j=0, *max_k=0, *min_i=0, *min_j=0, *min_k=0, *kapp_d=0, *kapp_ff=0;
	Zone *zone;
	SpPhys *pp;

	/*
	 * Inputs
	 */
	Mem_BZERO(&glb);

	/* Get molecule (optional) */
	if(!sts && SpPy_CheckOptionalInput("molec")) {
		sts = SpPy_GetInput_molec("molec", &glb.model.parms.mol);
	}

	/* Retrieve pygrid object */
	if(!sts) sts = SpPy_GetInput_PyObj("pygrid", &glb.pygrid);

	/* Get pointers to pygrid contents */
	if(!sts) {
		/* geom */
		obj = PyObject_GetAttrString(glb.pygrid, "geom");
		geom = Sp_PYINT(obj);
		Py_DECREF(obj);

		/* gas_to_dust */
		obj = PyObject_GetAttrString(glb.pygrid, "gas_to_dust");
		gas_to_dust = Sp_PYDBL(obj);
		Py_DECREF(obj);

		/* T_cmb */
		obj = PyObject_GetAttrString(glb.pygrid, "T_cmb");
		T_cmb = Sp_PYDBL(obj);
		Py_DECREF(obj);

		/* shape */
		obj = PyObject_GetAttrString(glb.pygrid, "shape");
		ndiv = GeVec3_s_Init(
			(size_t)Sp_PYINT(Sp_PYLST(obj, 0)),
			(size_t)Sp_PYINT(Sp_PYLST(obj, 1)),
			(size_t)Sp_PYINT(Sp_PYLST(obj, 2))
		);
		Py_DECREF(obj);

		/* min */
		obj = PyObject_GetAttrString(glb.pygrid, "min");
		min = GeVec3_d_Init(
			Sp_PYDBL(Sp_PYLST(obj, 0)),
			Sp_PYDBL(Sp_PYLST(obj, 1)),
			Sp_PYDBL(Sp_PYLST(obj, 2))
		);
		Py_DECREF(obj);

		/* max */
		obj = PyObject_GetAttrString(glb.pygrid, "max");
		max = GeVec3_d_Init(
			Sp_PYDBL(Sp_PYLST(obj, 0)),
			Sp_PYDBL(Sp_PYLST(obj, 1)),
			Sp_PYDBL(Sp_PYLST(obj, 2))
		);
		Py_DECREF(obj);

		min_i = PyObject_GetAttrString(glb.pygrid, "min_i");
		min_j = PyObject_GetAttrString(glb.pygrid, "min_j");
		min_k = PyObject_GetAttrString(glb.pygrid, "min_k");

		max_i = PyObject_GetAttrString(glb.pygrid, "max_i");
		max_j = PyObject_GetAttrString(glb.pygrid, "max_j");
		max_k = PyObject_GetAttrString(glb.pygrid, "max_k");

		n_H2 = PyObject_GetAttrString(glb.pygrid, "n_H2");
		X_mol = PyObject_GetAttrString(glb.pygrid, "X_mol");
		X_pH2 = PyObject_GetAttrString(glb.pygrid, "X_pH2");
		X_oH2 = PyObject_GetAttrString(glb.pygrid, "X_oH2");
		X_e = PyObject_GetAttrString(glb.pygrid, "X_e");
		X_H = PyObject_GetAttrString(glb.pygrid, "X_H");
		X_He = PyObject_GetAttrString(glb.pygrid, "X_He");
		T_k = PyObject_GetAttrString(glb.pygrid, "T_k");
		T_d = PyObject_GetAttrString(glb.pygrid, "T_d");
		T_ff = PyObject_GetAttrString(glb.pygrid, "T_ff");
		T_bb = PyObject_GetAttrString(glb.pygrid, "T_bb");
		V_t = PyObject_GetAttrString(glb.pygrid, "V_t");
		V_i = PyObject_GetAttrString(glb.pygrid, "V_i");
		V_j = PyObject_GetAttrString(glb.pygrid, "V_j");
		V_k = PyObject_GetAttrString(glb.pygrid, "V_k");

		kapp_d = PyObject_GetAttrString(glb.pygrid, "kapp_d");
		kapp_ff = PyObject_GetAttrString(glb.pygrid, "kapp_ff");
	}

	/* outf */
	if(!sts) sts = SpPy_GetInput_spfile("out", &glb.outf, Sp_NEW);

	/*
	 * Generate grid
	 */
	if(!sts) {
		glb.model.parms.gas_to_dust = gas_to_dust;
		glb.model.parms.T_cmb = T_cmb;
		SpModel_InitGrid(&glb.model, geom, min, max, ndiv);

		for(i = 0; i < GeVec3_X(ndiv, 0); i++) {
			for(j = 0; j < GeVec3_X(ndiv, 1); j++) {
				for(k = 0; k < GeVec3_X(ndiv, 2); k++) {
					zone = Zone_CHILD2(glb.model.grid, i, j, k);
					pp = zone->data;

					#define ARRAY(obj, i, j, k)\
						(*((double *)PyWrArray_GetPtr3((obj), (i), (j), (k))))
					#define ARRAY_STR(obj, i, j, k)\
						((char *)PyWrArray_GetPtr3((obj), (i), (j), (k)))


					zone->voxel = GeVox_Init(
						geom,
						ARRAY(min_i, i, j, k),
						ARRAY(min_j, i, j, k),
						ARRAY(min_k, i, j, k),
						ARRAY(max_i, i, j, k),
						ARRAY(max_j, i, j, k),
						ARRAY(max_k, i, j, k)
					);

					pp->n_H2 = ARRAY(n_H2, i, j, k);
					pp->X_mol = ARRAY(X_mol, i, j, k);
					pp->X_pH2 = ARRAY(X_pH2, i, j, k);
					pp->X_oH2 = ARRAY(X_oH2, i, j, k);
					pp->X_e = ARRAY(X_e, i, j, k);
					pp->X_H = ARRAY(X_H, i, j, k);
					pp->X_He = ARRAY(X_He, i, j, k);
					pp->T_k = ARRAY(T_k, i, j, k);
					pp->T_d = ARRAY(T_d, i, j, k);
					pp->T_ff = ARRAY(T_ff, i, j, k);
					pp->T_bb = ARRAY(T_bb, i, j, k);
					pp->V_t = ARRAY(V_t, i, j, k);

					GeVec3_X(pp->v_cen, 0) = ARRAY(V_i, i, j, k); /* m/s */
					GeVec3_X(pp->v_cen, 1) = ARRAY(V_j, i, j, k); /* m/s */
					GeVec3_X(pp->v_cen, 2) = ARRAY(V_k, i, j, k); /* m/s */

					strncpy(pp->kapp_d, ARRAY_STR(kapp_d, i, j, k), ZoneH5_KAPPLEN);
					strncpy(pp->kapp_ff, ARRAY_STR(kapp_ff, i, j, k), ZoneH5_KAPPLEN);

					#undef ARRAY
					#undef ARRAY_STR

					/* Set initial pops to either optically thin or LTE */
					if(pp->mol) {
						for(l = 0; l < pp->mol->nlev; l++) {
							for(m = 0; m < Sp_NTHREAD; m++) {
								pp->pops[m][l] = SpPhys_BoltzPops(pp->mol, l, pp->T_k);
							}
						}
					}
				}
			}
		}

		/* Write model to file */
		Sp_PRINT("Wrote model to `%s'\n", glb.outf->name);
		sts = SpIO_FwriteModel(glb.outf, glb.model);
	}

	/*
	 * Cleanup
	 */
	Py_XDECREF(glb.pygrid);
	Py_XDECREF(min_i);
	Py_XDECREF(min_j);
	Py_XDECREF(min_k);
	Py_XDECREF(max_i);
	Py_XDECREF(max_j);
	Py_XDECREF(max_k);
	Py_XDECREF(n_H2);
	Py_XDECREF(X_mol);
	Py_XDECREF(X_pH2);
	Py_XDECREF(X_oH2);
	Py_XDECREF(X_e);
	Py_XDECREF(X_H);
	Py_XDECREF(X_He);
	Py_XDECREF(T_k);
	Py_XDECREF(V_t);
	Py_XDECREF(V_i);
	Py_XDECREF(V_j);
	Py_XDECREF(V_k);
	SpModel_Cleanup(glb.model);
	SpIO_CloseFile(glb.outf);

	return sts;
}

















