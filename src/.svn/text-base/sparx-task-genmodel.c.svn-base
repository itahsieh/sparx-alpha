#include "sparx.h"

/* Global parameter struct */
static struct glb {
	int geom;
	size_t ndiv;
	double dims, ngas, tkin, xmol;
	SpModel model;
	SpFile *out_fp;
	char *fname;
} glb;

/* Keywords */
static SpKey
K_OUTF = Sp_KEY("out", "NewFile", 0, "Name of output file"),
K_NDIV = Sp_KEY("ndiv", "Size_t", "32", "Number of divisions along each axis"),
K_DIMS = Sp_KEY("dims", "Length", 0, "Dimension of each axis"),
K_LTEM = Sp_KEY("lte", "Molec", "Optional", "Init with given molecule in lte");

static SpKey *keys[] = {
	&K_OUTF,
	&K_NDIV,
	&K_DIMS,
	&K_LTEM,
	0
};

static int TaskMain(void);
static int ProcInps(void);
static int GenModel(void);
static void *GenModelThread(void *arg);
static double n_H2_Profile(const GeVec3_d *cen, const GeVec3_d *pos);
static double T_k_Profile(const GeVec3_d *cen, const GeVec3_d *pos);
static double X_mol_Profile(const GeVec3_d *cen, const GeVec3_d *pos);
static GeVec3_d v_Profile(const GeVec3_d *cen, const GeVec3_d *pos);
static void Cleanup(void);

/* Task definition */
SpTask SpTask_genmodel = Sp_TASK("genmodel", "Generic model generator", TaskMain, keys);

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

	/* K_LTEM: model.parms.mol */
	if(!status && !(o = SpInp_TASKGETKEY(K_LTEM)))
		status = 1;
	if(!status && o != Py_None) {
		glb.model.parms.mol = SpIO_FreadMolec(Sp_PYSTR(o));
		if(!glb.model.parms.mol)
			status = 1;
	}
	SpPy_XDECREF(o);

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
#define LENFAC (Sp_LENFAC)

static int GenModel(void)
{
	int status = 0;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, glb.dims / LENFAC, glb.dims / LENFAC, glb.dims / LENFAC);
	GeVec3_s naxes = GeVec3_INIT(NI, NJ, NK);
	Zone *root;

	/* Allocate grid and write to file */
	glb.model.grid = root = SpZone_ALLOC(&glb.model.parms);
	root->voxel = voxel;
	SpZone_GROW(root, naxes, &glb.model.parms);

	/* Multi-threading! */
	SpUtil_Threads(GenModelThread);

	/* Write model to file */
	status = SpIO_FwriteModel(glb.out_fp, glb.model);

	if(!status)
		Sp_PRINT("Wrote source model to `%s'\n", glb.out_fp->name);

	return status;
}

/*----------------------------------------------------------------------------*/

static void *GenModelThread(void *arg)
{
	size_t i, id = *((size_t *)arg), zone_id;
	Zone *zp, *root = glb.model.grid;
	SpPhys *pp;

	/* Setup cloud model */
	for(zp = Zone_GetMinLeaf(root), zone_id = 0; zp; zp = Zone_AscendTree(zp), zone_id++) {
		/* Skip when zone_id % Sp_NTHREAD != id */
		if(zone_id % Sp_NTHREAD != id)
			continue;

		pp = zp->data;

		pp->n_H2 = n_H2_Profile(&root->voxel.cen, &zp->voxel.cen);
		pp->T_k = T_k_Profile(&root->voxel.cen, &zp->voxel.cen);
		pp->X_mol = X_mol_Profile(&root->voxel.cen, &zp->voxel.cen);

		if(glb.model.parms.mol) {
			for(i = 0; i < glb.model.parms.mol->nlev; i++) {
				pp->pops[0][i] = SpPhys_BoltzPops(pp->mol, i, pp->T_k);
			}
		}

		pp->v_cen = v_Profile(&root->voxel.cen, &zp->voxel.cen);
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/

static double n_H2_Profile(const GeVec3_d *cen, const GeVec3_d *pos)
{
	double value = 0.0, radius;

	radius = GeVec3_Mag2(cen, pos);

	if(radius < 0.1)
		value = 1e4 * PHYS_UNIT_MKS_PERCC * (0.1 - radius) / 0.1;

	return value;
}

/*----------------------------------------------------------------------------*/

static double T_k_Profile(const GeVec3_d *cen, const GeVec3_d *pos)
{
	double value = 0.0, radius;

	radius = GeVec3_Mag2(cen, pos);

	if(radius < 0.1)
		value = 40.0 * (0.1 - radius) / 0.1;

	return value;
}

/*----------------------------------------------------------------------------*/

static double X_mol_Profile(const GeVec3_d *cen, const GeVec3_d *pos)
{
	double
		value = 0.0,
		radius,
		x = GeVec3_X(*pos, 0),
		y = GeVec3_X(*pos, 1),
		z = GeVec3_X(*pos, 2),
		x0 = GeVec3_X(*cen, 0),
		y0 = GeVec3_X(*cen, 1),
		z0 = GeVec3_X(*cen, 2);

	radius = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0));

	if(radius < 0.1)
		value = 1.0e-9;

	return value;
}

/*----------------------------------------------------------------------------*/

static GeVec3_d v_Profile(const GeVec3_d *cen, const GeVec3_d *pos)
{
	GeVec3_d vel = GeVec3_INIT(0, 0, 0);
	double
		radius,
		x = GeVec3_X(*pos, 0),
		y = GeVec3_X(*pos, 1),
		z = GeVec3_X(*pos, 2),
		x0 = GeVec3_X(*cen, 0),
		y0 = GeVec3_X(*cen, 1),
		z0 = GeVec3_X(*cen, 2);

	radius = GeVec3_Mag2(cen, pos);

	if(radius <= 0.1 && 1) {
		if(radius > 0) {
			GeVec3_X(vel, 0) = 10.0 * (x - x0) / radius; /* km/s */
			GeVec3_X(vel, 1) = 10.0 * (y - y0) / radius; /* km/s */
			GeVec3_X(vel, 2) = 10.0 * (z - z0) / radius; /* km/s */
		}
		else {
			GeVec3_X(vel, 0) = 0.0;
			GeVec3_X(vel, 1) = 0.0;
			GeVec3_X(vel, 2) = 0.0;
		}
	}

	return vel;
}









