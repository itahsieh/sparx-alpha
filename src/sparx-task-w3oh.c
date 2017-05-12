#include "sparx.h"

/* Global parameter struct */
static struct glb {
	int geom;
	size_t ndiv;
	double dims, kcen, tcen, rcen, rref, nref, nidx, tref, tidx, vref, vidx, xmol;
	SpModel model;
	SpFile *out_fp;
	char *fname;
	Kappa *kap;
} glb;

/* Keywords */
static SpKey
K_OUTF = Sp_KEY("out", "NewFile", 0, "Name of output file"),
K_GEOM = Sp_KEY("geom", "Geometry", "'sph1d'", "Geometry (coordinate system) of model"),
K_NDIV = Sp_KEY("ndiv", "Size_t", "32", "Number of divisions along each axis"),
K_DIMS = Sp_KEY("dims", "Length", 0, "Dimension of each axis"),
K_KCEN = Sp_KEY("kapcen", "Opacity", "'0.01cm^2g^-1'", "Opacity of central continuum source at 1.1mm"),
K_TCEN = Sp_KEY("tcen", "Temperature", "'5000K'", "Brightness temperature of central continuum source"),
K_RCEN = Sp_KEY("rcen", "Length", "'0.01pc'", "Radius of central continuum source"),
K_RREF = Sp_KEY("r0", "Length", 0, "Reference radius"),
K_NREF = Sp_KEY("n0", "NumDensity", 0, "Molecular gas number density at r0"),
K_NIDX = Sp_KEY("na", "Float", 0, "Molecular gas power index"),
K_TREF = Sp_KEY("t0", "Temperature", 0, "Kinetic temperature at r0"),
K_TIDX = Sp_KEY("ta", "Float", 0, "Kinetic temperature power index"),
K_VREF = Sp_KEY("v0", "Velocity", 0, "Radial velocity gradient (velocity/pc)"),
K_VIDX = Sp_KEY("va", "Float", 0, "Radial velocity gradient (velocity/pc)"),
K_XMOL = Sp_KEY("xmol", "Fraction", 0, "Molecular fractional abundance"),
K_TCMB = Sp_KEY("tcmb", "Temperature", "Optional", "Cosmic microwave background"),
K_GASD = Sp_KEY("gas2dust", "Float", "100.0", "Gas-to-dust ratio"),
K_DUST = Sp_KEY("dust", "Kappa", "KapPLaw('300GHz', '1cm^2g^-1', 2.0)", "Dust opacity model"),
K_LTEM = Sp_KEY("lte", "Molec", "Optional", "Init with given molecule in lte");

static SpKey *keys[] = {
	&K_OUTF,
	&K_NDIV,
	&K_DIMS,
	&K_KCEN,
	&K_TCEN,
	&K_RCEN,
	&K_RREF,
	&K_GEOM,
	&K_NREF,
	&K_NIDX,
	&K_TREF,
	&K_TIDX,
	&K_VREF,
	&K_VIDX,
	&K_XMOL,
	&K_TCMB,
	&K_GASD,
	&K_DUST,
	&K_LTEM,
	0
};

static int TaskMain(void);
static int ProcInps(void);
static int GenModel(void);
static void *GenModelThread(void *arg);
static void GenZone(Zone *zp);
static void GenZone_sph1d(Zone *zp);
static void Cleanup(void);

/* Task definition */
SpTask SpTask_w3oh = Sp_TASK("w3oh", "W3(OH) model generator", TaskMain, keys);

/*----------------------------------------------------------------------------*/

static int TaskMain(void)
{
	int status = 0;

	Mem_BZERO(&glb);

	status = ProcInps();

	/* Check for unused keys */
	if(!status)
		status = SpInp_CheckKeys();

	/* Generate model */
	if(!status)
		status = GenModel();

	/* Cleanup */
	Cleanup();

	return status;
}

/*----------------------------------------------------------------------------*/

static int ProcInps(void)
{
	int status = 0;
	PyObject *o;

	/* K_GEOM: geometry */
	if(!status && !(o = SpInp_TASKGETKEY(K_GEOM)))
		status = 1;
	if(!status)
		glb.geom = Sp_PYINT(o);
	SpPy_XDECREF(o);

	/* K_OUTF: out_fp */
	if(!status && !(o = SpInp_TASKGETKEY(K_OUTF)))
		status = 1;
	if(!status && !(glb.out_fp = SpIO_OpenFile(Sp_PYSTR(o), Sp_NEW)))
		status = 1;
	SpPy_XDECREF(o);

	/* K_NDIV: ndiv */
	if(!status && !(o = SpInp_TASKGETKEY(K_NDIV)))
		status = 1;
	if(!status)
		glb.ndiv = Sp_PYSIZE(o);
	SpPy_XDECREF(o);

	/* K_DIMS: dims */
	if(!status && !(o = SpInp_TASKGETKEY(K_DIMS)))
		status = 1;
	if(!status)
		glb.dims = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/* K_KCEN: kcen */
	if(!status && !(o = SpInp_TASKGETKEY(K_KCEN)))
		status = 1;
	if(!status)
		glb.kcen = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/* K_TCEN: tcen */
	if(!status && !(o = SpInp_TASKGETKEY(K_TCEN)))
		status = 1;
	if(!status)
		glb.tcen = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/* K_RCEN: rcen */
	if(!status && !(o = SpInp_TASKGETKEY(K_RCEN)))
		status = 1;
	if(!status)
		glb.rcen = Sp_PYDBL(o) / PHYS_UNIT_MKS_PC;
	SpPy_XDECREF(o);

	/* K_RREF: rref */
	if(!status && !(o = SpInp_TASKGETKEY(K_RREF)))
		status = 1;
	if(!status)
		glb.rref = Sp_PYDBL(o) / PHYS_UNIT_MKS_PC;
	SpPy_XDECREF(o);

	/* K_NREF: nref */
	if(!status && !(o = SpInp_TASKGETKEY(K_NREF)))
		status = 1;
	if(!status)
		glb.nref = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/* K_NIDX: nidx */
	if(!status && !(o = SpInp_TASKGETKEY(K_NIDX)))
		status = 1;
	if(!status)
		glb.nidx = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/* K_TREF: tref */
	if(!status && !(o = SpInp_TASKGETKEY(K_TREF)))
		status = 1;
	if(!status)
		glb.tref = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/* K_TIDX: tidx */
	if(!status && !(o = SpInp_TASKGETKEY(K_TIDX)))
		status = 1;
	if(!status)
		glb.tidx = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/* K_VREF: vref */
	if(!status && !(o = SpInp_TASKGETKEY(K_VREF)))
		status = 1;
	if(!status)
		glb.vref = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/* K_VIDX: vidx */
	if(!status && !(o = SpInp_TASKGETKEY(K_VIDX)))
		status = 1;
	if(!status)
		glb.vidx = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/* K_XMOL: xmol */
	if(!status && !(o = SpInp_TASKGETKEY(K_XMOL)))
		status = 1;
	if(!status)
		glb.xmol = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/* K_LTEM: model.parms.mol */
	if(!status && !(o = SpInp_TASKGETKEY(K_LTEM)))
		status = 1;
	if(!status && o != Py_None) {
		glb.model.parms.mol = SpIO_FreadMolec(Sp_PYSTR(o));
		if(!glb.model.parms.mol)
			status = 1;
	}
	SpPy_XDECREF(o);

	/* K_TCMB: model.parms.T_cmb */
	if(!status && !(o = SpInp_TASKGETKEY(K_TCMB)))
		status = 1;
	if(!status && o != Py_None) {
		glb.model.parms.T_cmb = Sp_PYDBL(o);
	}
	SpPy_XDECREF(o);

	/* K_GASD: model.parms.gas_to_dust */
	if(!status && !(o = SpInp_TASKGETKEY(K_GASD)))
		status = 1;
	if(!status && o != Py_None) {
		glb.model.parms.gas_to_dust = Sp_PYDBL(o);
	}
	SpPy_XDECREF(o);

	/* K_DUST: glb.kap */
	glb.kap = SpInp_GetKey_kappa(K_DUST.name);

	return status;
}

/*----------------------------------------------------------------------------*/

static void Cleanup(void)
{
	if(glb.fname)
		free(glb.fname);

	SpModel_Cleanup(glb.model);

	if(glb.out_fp)
		SpIO_CloseFile(glb.out_fp);

	return;
}

/*----------------------------------------------------------------------------*/

#define NI (glb.ndiv)
#define NJ (glb.ndiv)
#define NK (glb.ndiv)
#define PC (PHYS_UNIT_MKS_PC)

static int GenModel(void)
{
	int status = 0;
	GeVec3_d min, max;
	GeVec3_s ndiv;

	switch(glb.geom) {
		case GEOM_SPH1D:
			min = GeVec3_d_Init(0.0, 0.0, 0.0);
			max = GeVec3_d_Init(glb.dims / PC, PI, TWOPI);
			ndiv = GeVec3_s_Init(glb.ndiv, (size_t)1, (size_t)1);
			break;

		case GEOM_REC3D:
			min = GeVec3_d_Init(0.0, 0.0, 0.0);
			max = GeVec3_d_Init(glb.dims/PC, glb.dims/PC, glb.dims/PC);
			ndiv = GeVec3_s_Init(glb.ndiv, glb.ndiv, glb.ndiv);
			break;

		default: /* Shouldn't happen */
			assert(0);
	}

	/* Allocate grid */
	SpModel_InitGrid(&glb.model, glb.geom, min, max, ndiv);

	/* Fill in grid parameters */
	SpUtil_Threads(GenModelThread);

	/* Write model to file */
	status = SpIO_FwriteModel(glb.out_fp, glb.model);

	if(!status)
		Sp_PRINT("Wrote source model to `%s'\n", glb.out_fp->name);

	return status;
}

/*----------------------------------------------------------------------------*/

static void *GenModelThread(void *tid_p)
{
	size_t tid = *((size_t *)tid_p), zone_id;
	Zone *zp, *root = glb.model.grid;

	/* Setup model */
	for(zp = Zone_GetMinLeaf(root), zone_id = 0; zp; zp = Zone_AscendTree(zp), zone_id++) {
		/* Skip when zone_id % Sp_NTHREAD != tid */
		if(zone_id % Sp_NTHREAD != tid)
			continue;

		GenZone(zp);
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/

static void GenZone(Zone *zp)
{
	switch(zp->voxel.geom) {
		case GEOM_SPH1D:
			GenZone_sph1d(zp);
			break;

		default: /* This should never happen */
			assert(0);
	}

	return;
}

/*----------------------------------------------------------------------------*/

static void GenZone_sph1d(Zone *zp)
{
	size_t i;
	SpPhys *pp;

	pp = zp->data;

	/* n_H2 */
	pp->n_H2 = glb.nref * pow(GeVec3_X(zp->voxel.cen, 0) / glb.rref, glb.nidx); /* m^-3 */

	/* T_k */
	pp->T_k = glb.tref * pow(GeVec3_X(zp->voxel.cen, 0) / glb.rref, glb.tidx); /* K */

	/* V_cen */
	GeVec3_X(pp->v_cen, 0) = glb.vref * pow(GeVec3_X(zp->voxel.cen, 0) / glb.rref, glb.vidx); /* m/s */

	/* X_mol */
	pp->X_mol = glb.xmol; /* fraction */

	strncpy(pp->kapp_d, glb.kap->name, ZoneH5_KAPPLEN);

	/* Init molecular level populations if requested */
	if(glb.model.parms.mol) {
		for(i = 0; i < pp->mol->nlev; i++) {
			pp->pops_preserve[i] = SpPhys_BoltzPops(glb.model.parms.mol, i, pp->T_k);
		}
	}

	/* Add central continuum source */
	if(GeVec3_X(zp->voxel.cen, 0) <= glb.rcen) {
		pp->T_ff = glb.tcen;
		snprintf(pp->kapp_ff, ZoneH5_KAPPLEN, "powerlaw,%10.3e,%10.3e,%10.3e", 
			PHYS_CONST_MKS_LIGHTC / 1.1e-3,
			glb.kcen,
			2.0);
	}


	return;
}










