#include "sparx.h"

/* Global parameter struct */
static struct glb {
	int geom;
	size_t ndiv;
	double dims, ngas, tkin, xmol, vrad;
	Kappa *dust;
	SpModel model;
	SpFile *outf;
} glb;

static int GenModel(void);
static void *GenModelThread(void *arg);

/*----------------------------------------------------------------------------*/

int SpTask_Uniform(void)
{
	int status = 0;

	Mem_BZERO(&glb);

	glb.geom = SpInp_GetKey_int("geom");
	glb.dims = SpInp_GetKey_dbl("dims");
	glb.ndiv = SpInp_GetKey_size_t("ndiv");
	glb.ngas = SpInp_GetKey_dbl("nh2");
	glb.tkin = SpInp_GetKey_dbl("tk");
	glb.xmol = SpInp_GetKey_dbl("xmol");
	glb.vrad = SpInp_GetKey_dbl("vrad");
	glb.dust = SpInp_GetKey_kappa("dust");
	glb.model.parms.gas_to_dust = SpInp_GetKey_dbl("gas2dust");
	glb.outf = SpInp_GetKey_spfile("out", Sp_NEW);

	/* Generate model */
	if(!status)
		status = GenModel();

	/* Cleanup */
	SpModel_Cleanup(glb.model);

	if(glb.outf)
		SpIO_CloseFile(glb.outf);


	return status;
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

	/* Models of different coordinate systems differ only in their
	 * dimensions */
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
	status = SpIO_FwriteModel(glb.outf, glb.model);

	if(!status)
		Sp_PRINT("Wrote source model to `%s'\n", glb.outf->name);

	return status;
}

/*----------------------------------------------------------------------------*/

static void *GenModelThread(void *tid_p)
{
	size_t i, tid = *((size_t *)tid_p), zone_id;
	Zone *zp, *root = glb.model.grid;
	double radius, velo;
	SpPhys *pp;
	GeVec3_d pos;

	/* Setup model */
	for(zp = Zone_GetMinLeaf(root), zone_id = 0; zp; zp = Zone_AscendTree(zp), zone_id++) {
		/* Skip when zone_id % Sp_NTHREAD != tid */
		if(zone_id % Sp_NTHREAD != tid)
			continue;

		if(zp->children)
			continue;

		pp = NULL;

		switch(glb.geom) {
			case GEOM_SPH1D:
				pp = zp->data;
				GeVec3_X(pp->v_cen, 0) = glb.vrad * GeVec3_X(zp->voxel.cen, 0);
				break;

			case GEOM_REC3D:
				pos = GeVec3_Sub(&root->voxel.cen, &zp->voxel.cen);
				radius = GeVec3_Mag(&pos);
				if(radius <= 0.5 * glb.dims / PHYS_UNIT_MKS_PC) {
					pp = zp->data;

					velo = radius * glb.vrad;
					pos = GeVec3_Normalize(&pos);
					pp->v_cen = GeVec3_Scale(&pos, velo);
				}
				break;

			default:
				assert(0);
		}

		if(pp) {
			pp->T_k = glb.tkin; /* K */
			pp->n_H2 = glb.ngas; /* m^-3 */
			pp->X_mol = glb.xmol; /* fraction */
			if(glb.dust) {
				strncpy(pp->kapp_d, glb.dust->name, ZoneH5_KAPPLEN);
			}

			/* Init molecular level populations if requested */
			if(glb.model.parms.mol) {
				for(i = 0; i < pp->mol->nlev; i++) {
					pp->pops_preserve[i] = SpPhys_BoltzPops(glb.model.parms.mol, i, pp->T_k);
				}
			}
		}
	}

	return NULL;
}












