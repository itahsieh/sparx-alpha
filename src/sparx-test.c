#include <gsl/gsl_rng.h>

#include "sparx.h"

#if Sp_MIRSUPPORT
#include "miriad.h"
#endif

/* Prototypes of test tasks should no longer be kept here
   but in the sparx.h header, since they should be called from PyMain() */

void SpTest_MirImg_LoadXYV(const char *fname);
static int SpTest_Kappa(void);
static int SpTest_RayTracing3(void);
static int SpTest_RayTracing2(void);
static int SpTest_Math(void);
static int SpTest_RayAW(void);
static int SpTest_Nothing(void);
static int SpTest_Model(void);
static int SpTest_HDF5(void);
static int SpTest_RandRay(void);
static int SpTest_Phys(void);
static int SpTest_Zone(void);
static int SpTest_MirImg(void);
static int SpTest_Pkeys(void);
static int SpTest_Llst(void);
static int SpTest_RayTracing(void);
static int SpTest_Molec(void);
static int SpTest_VecInterp(void);

/* Task definitions */
SpTask
	SpTask_t_kappa = Sp_TASK("t_kappa", "Test opacities loading and interpolation", SpTest_Kappa, 0),
#if Sp_MIRSUPPORT
	SpTask_t_raytracing3 = Sp_TASK("t_raytracing3", "Test ray tracing in various geometries (zone by zone)", SpTest_RayTracing3, 0),
	SpTask_t_raytracing2 = Sp_TASK("t_raytracing2", "Test ray tracing in various geometries", SpTest_RayTracing2, 0),
	SpTask_t_math = Sp_TASK("t_math", "Test QR decomposition", SpTest_Math, 0),
	SpTask_t_rayaw = Sp_TASK("t_rayaw", "Test ray tracing with Amanatides-Woo algorithm", SpTest_RayAW, 0),
	SpTask_t_randray = Sp_TASK("t_randray", "Test physics.", SpTest_RandRay, 0),
	SpTask_t_mirimg = Sp_TASK("t_mirimg", "A demonstration of writing Miriad images to disk.", SpTest_MirImg, 0),
	SpTask_t_raytracing = Sp_TASK("t_raytracing", "For debugging purposes", SpTest_RayTracing, 0),
#endif
	
	
	SpTask_t_nothing = Sp_TASK("t_nothing", "This is for generating valgrind suppressions.", SpTest_Nothing, 0),
	SpTask_t_model = Sp_TASK("t_model", "Test model allocation and i/o.", SpTest_Model, 0),
	SpTask_t_hdf5 = Sp_TASK("t_hdf5", "Test HDF5 access.", SpTest_HDF5, 0),	
	SpTask_t_phys = Sp_TASK("t_phys", "Test physics.", SpTest_Phys, 0),
	SpTask_t_zone = Sp_TASK("t_zone", "Test zone manipulation.", SpTest_Zone, 0),
	SpTask_t_pkeys = Sp_TASK("t_pkeys", "For debugging purposes", SpTest_Pkeys, 0),
	SpTask_t_llst = Sp_TASK("t_llst", "For debugging purposes", SpTest_Llst, 0),
	SpTask_t_molec = Sp_TASK("t_molec", "Test molecule access", SpTest_Molec, 0),
	SpTask_t_vecinterp = Sp_TASK("t_vecinterp", "Test vector interpolation", SpTest_VecInterp, 0);



/*----------------------------------------------------------------------------*/

int SpTest_Key(void)
{
	int boolean;
	boolean = SpInp_GetKey_TF("test");

	printf(boolean ? "True\n" : "False\n");

	return 0;
}

/*----------------------------------------------------------------------------*/

int SpTest_Gaussian(void)
{
	size_t i, n = 1001;
	double v_min = -100.0, v_max = 100.0, v_delta, vel, width = 50.0;

	v_delta = (v_max - v_min) / (double)(n - 1);

	for(i = 0; i < 1001; i++) {
		vel = v_min + v_delta * (double)i;
		printf("%20g %20g\n", vel, Num_GaussNormal(vel, width));
	}

	return 0;
}


/*----------------------------------------------------------------------------*/

int SpTest_Interp2d(void)
{
	#define idim 6
	#define jdim 6
	#define idim2 50
	#define jdim2 50
	#define coord {0, 10, 20, 30, 40, 50}

	int i, j;
	int icen = idim / 2, jcen = jdim / 2;
	double *array2d = Mem_CALLOC(idim * jdim, array2d),
		x1a[idim] = coord,
		x2a[jdim] = coord;
	double width = 1.0, interp_value, x, y;
	#define ARR(i, j)\
		array2d[j + jdim * i]
	FILE *fp = fopen("SpTest_Interp2d.txt", "w");

	/* Setup a 2D gaussian profile exp(-(x^2/width^2 + y^2/width^2)) */
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim ; j++) {
			ARR(i, j) = exp(-(pow((i - icen) / width, 2.0) + pow((j - jcen) / width, 2.0)));
		}
	}

	fprintf(fp, "# Before:\n");
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim ; j++) {
			fprintf(fp, j == 0 ? "%10.3e" : " %10.3e", ARR(i, j));
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "# After:\n");
	for(i = 0; i < idim2; i++) {
		for(j = 0; j < jdim2; j++) {
			x = x1a[0] + (double)i * ((x1a[idim-1] - x1a[0]) / (double)idim2);
			y = x2a[0] + (double)j * ((x2a[jdim-1] - x2a[0]) / (double)jdim2);
			interp_value = Num_InterpPoly2d(x, y, x1a, x2a, array2d, (size_t)idim, (size_t)jdim, (size_t)4, (size_t)4);
			//interp_value = num_interp_poly_2d(x, y, x1a, x2a, array2d, idim, jdim, 4, 4);
			fprintf(fp, j == 0 ? "%10.3e" : " %10.3e", interp_value);
		}
		fprintf(fp, "\n");
	}

	#undef idim
	#undef jdim

	return 0;
}

/*----------------------------------------------------------------------------*/

int SpTest_FFT(void)
{
	int idim = 128, jdim = 128, i, j;
	int icen = idim / 2, jcen = jdim / 2;
	double *array2d = Mem_CALLOC(idim * jdim, array2d),
		*real = Mem_CALLOC(idim*jdim, real),
		*imag = Mem_CALLOC(idim*jdim, imag);
	double width = 1.0;
	#define ARR(i, j)\
		array2d[j + jdim * i]
	#define REAL(i, j)\
		real[j + jdim * i]
	#define IMAG(i, j)\
		imag[j + jdim * i]
	FILE *fp = fopen("SpTest_FFT.txt", "w");

	/* Setup a 2D gaussian profile exp(-(x^2/width^2 + y^2/width^2)) */
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim ; j++) {
			ARR(i, j) = exp(-(pow((i - icen) / width, 2.0) + pow((j - jcen) / width, 2.0)));
		}
	}

	fprintf(fp, "# Before:\n");
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim ; j++) {
			fprintf(fp, j == 0 ? "%10.3e" : " %10.3e", ARR(i, j));
		}
		fprintf(fp, "\n");
	}

	NumFFT_Xform2d(array2d, (size_t)idim, (size_t)jdim, real, imag);

	fprintf(fp, "\n\n# After (real):\n");
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim ; j++) {
			fprintf(fp, j == 0 ? "%10.3e" : " %10.3e", REAL(i, j));
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n\n# After (imag):\n");
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim ; j++) {
			fprintf(fp, j == 0 ? "%10.3e" : " %10.3e", IMAG(i, j));
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	free(array2d);
	free(real);
	free(imag);

	Sp_PRINT("Wrote data to 'SpTest_FFT.txt'\n");
	Sp_PRINT("You may now plot the data, e.g. with gnuplot,\n");
	Sp_PRINT("gnuplot> set pm3d map\n");
	Sp_PRINT("gnuplot> splot 'SpTest_FFT.txt' i 0 matrix w pm3d # Input Gaussian\n");
	Sp_PRINT("gnuplot> splot 'SpTest_FFT.txt' i 1 matrix w pm3d # Output (real)\n");
	Sp_PRINT("gnuplot> splot 'SpTest_FFT.txt' i 2 matrix w pm3d # Output (imag)\n");

	return 0;
}

/*----------------------------------------------------------------------------*/

int SpTest_Test(void)
{
	int status = 0;
	PyObject *o;

	/* xmol */
	if(!status && !(o = SpInp_GETKEYSTR("xmol")))
		status = 1;

	if(!status)
		Deb_PRINT("xmol=%g!!!\n", Sp_PYDBL(o));

	if(!status) {
		Deb_PRINT("Hooray!\n");
	}

	return status;
}

/*----------------------------------------------------------------------------*/

static int SpTest_Kappa(void)
{
	Kappa *kap;
	
	size_t i, n = 1000;
	double c = PHYS_CONST_MKS_LIGHTC;
	double nu_min, nu_max, nu_del, nu;
	FILE *fp = fopen("kappa.dat", "w");

	nu_min = c / 2.0e-6;
	nu_max = c / 1.0e-3;
	nu_del = (nu_max - nu_min) / (double)n;

	//kap = Kap_New_Table(Sp_ROOT"/data/opacity/jena_bare_e5.tab");
	kap = SpIO_FreadKappa("jena_bare_e5");

	fprintf(fp, "# Table:\n");
	for(i = 0; i < kap->nrows; i++) {
		fprintf(fp, "%20g %20g\n", kap->freq[i], kap->kappa[i]);
	}

	fprintf(fp, "\n\n# Interpolated:\n");

	for(i = 0; i < n; i++) {
		nu = nu_min + nu_del * (double)i;
		fprintf(fp, "%20g %20g\n", nu, Kap_FromFreq(kap, nu));
	}

	Kap_Free(kap);

	return 0;
}

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

static int SpTest_Math(void)
{
	double A[9] = {1, 3, 5, 2, 5, 1, 2, 3, 8}, b[3] = {10, 8, 3}, x[3];
	size_t M = 3, N = 2;

	/* Test quadratic equation solver */
	double x1, x2;
	/* Distinct real roots -2 and 1 */
	Num_QuadraticRoots(1.0, 1.0, -2.0, &x1, &x2);
	printf("/* Distinct real roots -2 and 1 */ x1=%g, x2=%g\n", x1, x2);

	/* Double real roots 5 */
	Num_QuadraticRoots(1.0, -10.0, 25.0, &x1, &x2);
	printf("/* Double real roots 5 */ x1=%g, x2=%g\n", x1, x2);

	/* Complex conjugate roots +/-2i */
	Num_QuadraticRoots(1.0, 0.0, 4.0, &x1, &x2);
	printf("/* Complex conjugate roots +/-2i */ x1=%g, x2=%g\n", x1, x2);

	Num_QRDecompSolve(A, M, N, b, x);

	/* Normalized position vector */
	GeVec3_d a = GeVec3_INIT(-1, 2, 3);
	Deb_PRINT("Before normalization:"); GeVec3_PRINT(stdout, a); printf("\n");
	a = GeVec3_Normalize(&a);
	Deb_PRINT("After normalization:"); GeVec3_PRINT(stdout, a); printf("\n");
	Deb_PRINT("mag=%g\n", GeVec3_Mag(&a));
	a = GeVec3_Scale(&a, 3.0);
	Deb_PRINT("After scaling by 3:"); GeVec3_PRINT(stdout, a); printf("\n");
	Deb_PRINT("mag=%g\n", GeVec3_Mag(&a));

	return 0;
}

/*----------------------------------------------------------------------------*/

static int SpTest_Nothing(void)
/* This task is just for generating valgrind supressions */
{
	return 0;
}

/*----------------------------------------------------------------------------*/

#define NDIV ((size_t)2)

static int SpTest_Model(void)
{
	size_t i;
	Zone *zp;
	SpPhys *pp;
	int status = 0;
	SpModel model;
	const char *fname = "test.h5";
	SpFile *sfp = 0;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, 0.1, 0.1, 0.1);
	GeVec3_s naxes = GeVec3_INIT(NDIV, NDIV, NDIV);

	/* Create and write model */
	Mem_BZERO(&model);

	if(1) {
		model.parms.mol = SpIO_FreadMolec("o-h2o_2lev");
		if(!model.parms.mol)
			status = 1;
	}

	model.grid = SpZone_ALLOC(&model.parms);
	model.grid->voxel = voxel;
	SpZone_GROW(model.grid, naxes, &model.parms);

	if(model.parms.mol) {
		for(zp = Zone_GetMinLeaf(model.grid); zp; zp = Zone_AscendTree(zp)) {
			pp = zp->data;
			pp->T_k = 40.0;
			for(i = 0; i < model.parms.mol->nlev; i++) {
				pp->pops[0][i] = SpPhys_BoltzPops(model.parms.mol, i, pp->T_k);
			}
		}
	}

	if(!status && !(sfp = SpIO_OpenFile(fname, Sp_TRUNC)))
		status = 1;

	if(!status)
		status = SpIO_FwriteModel(sfp, model);

	if(sfp)
		SpIO_CloseFile(sfp);

	SpModel_Cleanup(model);

	/* Open and read model */
	Mem_BZERO(&model);

	sfp = 0;

	if(!status && !(sfp = SpIO_OpenFile(fname, Sp_OLD)))
		status = 1;

	if(!status)
		status = SpIO_FreadModel(sfp, sfp, &model);

	if(sfp)
		SpIO_CloseFile(sfp);

	if(!status)
		SpModel_PrintModel(model);

	SpModel_Cleanup(model);

	return status;
}

/*----------------------------------------------------------------------------*/

#define FNAME "test.h5"

static int SpTest_HDF5(void)
{
	Zone *root;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, 0.1, 0.1, 0.1);
	GeVec3_s naxes = GeVec3_INIT(NDIV, NDIV, NDIV);
	hid_t file_id;

	root = SpZone_ALLOC(0);
	root->voxel = voxel;
	SpZone_GROW(root, naxes, 0);

	/* Open HDF5 file */
	file_id = H5Fcreate(FNAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* Write grid to file */
	SpIO_H5WriteGrid(file_id, root);

	/* Cleanup */
	SpZone_FREE(root);

	/* Close the file */
	H5Fclose(file_id);

	/* Reopen file */
	file_id = H5Fopen(FNAME, H5F_ACC_RDONLY, H5P_DEFAULT);

	root = 0;

	/* Load grid data */
	SpIO_H5ReadGrid(file_id, file_id, &root, 0);

	/* Close file */
	H5Fclose(file_id);
	
	/* Print grid */
	SpZone_FPRINTF(stdout, root);

	/* Free grid */
	SpZone_FREE(root);

	return 0;
}

#undef NDIV
#undef FNAME


/*----------------------------------------------------------------------------*/

static int SpTest_Phys(void)
{
	int status = 0;
	size_t i;
	double T_k;
	const char *molec = "o-h2o_2lev";
	Molec *mol;

	/* 2-level molecules are SIGNIFICANTLY faster when calculating LTE pops! */
	mol = SpIO_FreadMolec(molec);

	Deb_ASSERT(mol != 0);

	for(i = 0; i < 10; i++) {
		T_k = 100.0 * ((int)i + 1) / 10.0;
		printf("Botlzmann ratio for %s 1-0 @ %g K=%g\n", molec, T_k, SpPhys_BoltzRatio(mol, 1, 0, T_k));
	}

	Mol_Free(mol);

	return status;
}

/*----------------------------------------------------------------------------*/

#define NDIV 2
#define NI NDIV
#define NJ NDIV
#define NK NDIV

static int SpTest_Zone(void)
{
	int status = 0, non_empty = 0;
	Zone *root, *zp;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, 0.1, 0.1, 0.1);
	GeVec3_s naxes = GeVec3_INIT(NDIV, NDIV, NDIV);
	SpFile *sfp;
	SpModel model;
	int test = 0;
	size_t i;
	double radius, ipos, jpos, kpos;
	SpPhys *pp;
	SpPhysParm prm;
	size_t num = 0;
	time_t tm0, tm1;

	Mem_BZERO(&model);

	/* 2-level molecules are SIGNIFICANTLY faster when calculating LTE pops! */
	model.parms.mol = SpIO_FreadMolec("o-h2o_2lev");

	Deb_ASSERT(model.parms.mol != 0);

	prm = model.parms;

	/* Allocate grid and write to file */
	if(!status) {
		/* Allocate model grid, assume dimensions of root node is 1pc x 1pc x 1pc */
		model.grid = root = Zone_Alloc(0, 0, SpPhys_Alloc, &prm);
		root->voxel = voxel;
		Zone_GrowChildren(root, naxes, SpPhys_Alloc, &prm);
		Zone_GrowChildren(Zone_CHILD2(root, 0, 0, 0), naxes, SpPhys_Alloc, &prm);

		Deb_PRINT("Setting up model...\n");
		time(&tm0);
		num = 0;
		for(zp = Zone_GetMinLeaf(root), num = 0; zp; zp = Zone_AscendTree(zp)) {
			/* Calculate positions */
			ipos = (int)GeVec3_X(zp->index, 0) - (int)NI / 2 + 0.5;
			jpos = (int)GeVec3_X(zp->index, 1) - (int)NJ / 2 + 0.5;
			kpos = (int)GeVec3_X(zp->index, 2) - (int)NK / 2 + 0.5;

			switch(test) {
				case 0:
					/* Cubic model */
					non_empty = 1;
					break;

				case 1:
					/* Plane parallel model */
					non_empty = fabs(ipos) <= 2;
					break;


				case 2:
					/* Spherical model */
					radius = sqrt(ipos * ipos + jpos * jpos + kpos * kpos);
					non_empty = radius <= NI / 2;
					break;

				case 3:
					/* Disk model */
					radius = sqrt(ipos * ipos + jpos * jpos);
					non_empty = (radius <= NI / 2 && fabs(kpos) <= 2);
					break;

				default:
					/* No such case! */
					Deb_ASSERT(0);
			}

			if(non_empty) {
				pp = zp->data;
				pp->T_k = 40.0; /* K */
				pp->n_H2 = 1.0e4 * PHYS_UNIT_MKS_PERCC; /* m^-3 */
				pp->X_mol = 1.0e-8;
				pp->width = 192.0;

				for(i = 0; i < prm.mol->nlev; i++) {
					pp->pops[0][i] = SpPhys_BoltzPops(prm.mol, i, pp->T_k);
				}

				if(0) {
					GeVec3_X(pp->v_cen, 0) = -1000;
					GeVec3_X(pp->v_cen, 1) = 0;
					GeVec3_X(pp->v_cen, 2) = 0;

					for(i = 0; i < 8; i++) {
						GeVec3_X(pp->v_edge[i], 0) = -1000;
						GeVec3_X(pp->v_edge[i], 1) = 0;
						GeVec3_X(pp->v_edge[i], 2) = 0;
					}
				}
			}
			num++;
		}
		time(&tm1);
		Deb_PRINT("Done, num=%d, delta_t=%d\n", num, tm1 - tm0);
	}

	sfp = SpIO_OpenFile("zone_test.sparx", Sp_NEW);
	if(sfp) {
		SpIO_FwriteModel(sfp, model);
		SpIO_CloseFile(sfp);
	}
	else
		status = 1;

	SpModel_Cleanup(model);
	Mem_BZERO(&model);

	/* Load grid from file */
	if(!status) {
		sfp = SpIO_OpenFile("zone_test.sparx", Sp_OLD);
	}

	if(sfp) {
		SpIO_FreadModel(sfp, sfp, &model);
		SpIO_CloseFile(sfp);

		SpZone_FPRINTF(stdout, model.grid);

		if(model.parms.mol)
			Mol_Fprintf(stdout, model.parms.mol);
	}

	SpModel_Cleanup(model);
	Mem_BZERO(&model);

	return status;
}

#undef NDIV
#undef NI
#undef NJ
#undef NK

/*----------------------------------------------------------------------------*/

#define LEN 100.0

static int SpTest_VecInterp(void)
{
	size_t i, j, NI = 20, NJ = 20;
	float
		*mapx = Mem_CALLOC(NI * NJ, mapx),
		*mapy = Mem_CALLOC(NI * NJ, mapy),
		tr[6] = {-0.5, 1, 0, -0.5, 0, 1};
	GeVec3_d
		a = GeVec3_INIT(5.0, 4.0, 0.0), 
		b = GeVec3_INIT(-4.0, -1.0, 0.0),
		xa = GeVec3_INIT(0.0, 0.0, 0.0),
		xb = GeVec3_INIT(LEN, LEN, LEN),
		pos;
	GeVec3_d *map = Mem_CALLOC(NI * NJ, map);

	tr[1] = (float)(LEN / ((double)NI - 1));
	tr[5] = (float)(LEN / ((double)NJ - 1));

	#define MAP(i, j)\
		map[(j) + NI * (i)]

	#define MAPX(i, j)\
		mapx[(j) + NI * (i)]

	#define MAPY(i, j)\
		mapy[(j) + NI * (i)]

	for(i = 0; i < NI; i++) {
		for(j = 0; j < NJ; j++) {
			GeVec3_X(pos, 0) = GeVec3_X(xa, 0) + (double)i * LEN / ((double)NI - 1);
			GeVec3_X(pos, 1) = GeVec3_X(xa, 1) + (double)j * LEN / ((double)NJ - 1);
			GeVec3_X(pos, 2) = 0.0;

			MAP(i, j) = GeVec3_InterpLinear(&xa, &xb, &a, &b, &pos);
			MAPX(i, j) = (float)GeVec3_X(MAP(i, j), 0);
			MAPY(i, j) = (float)GeVec3_X(MAP(i, j), 1);
		}
	}
#if Sp_MIRSUPPORT
	cpgopen("/xs");
	Cpg_Env(-5.0, LEN + 10.0, -5.0, LEN + 10.0, 1, 0);
	Cpg_Sch(0.6);
	Cpg_Vect(mapx, mapy, NI, NJ, 1, NI, 1, NJ, 0.0, 0, tr, HUGE_VAL);
	cpgclos();
#endif
	free(map);
	free(mapx);
	free(mapy);

	#undef MAP
	#undef MAPX
	#undef MAPY

	return 0;
}

#undef LEN

/*----------------------------------------------------------------------------*/

static int SpTest_Molec(void)
{
	int status = 0;
	char *mname = 0;
	Molec *mol = 0;
	FILE *fp = 0;
	PyObject *o;

	/* K_MOLE: model.parms.mol */
	if(!status && !(o = SpInp_GETKEYSTR("molec")))
		status = 1;
	if(!status) {
		mol = SpIO_FreadMolec(Sp_PYSTR(o));
		if(!mol)
			status = 1;
	}
	SpPy_XDECREF(o);

	if(fp)
		fclose(fp);

	if(1 && !status)
		Mol_Fprintf(stdout, mol);
	else {
		if(!status) {
			fp = fopen("mol_test.mol", "wb");
			Mol_FwriteBinary(mol, fp);
			fclose(fp);
		}

		if(mol) {
			Mol_Free(mol);
			mol = 0;
		}

		if(!status) {
			fp = fopen("mol_test.mol", "rb");
			mol = Mol_FreadBinary(fp);
			fclose(fp);
			Mol_Fprintf(stdout, mol);
		}
	}

	if(mol)
		Mol_Free(mol);

	if(mname)
		free(mname);

	return status;
}



/*----------------------------------------------------------------------------*/

static int SpTest_Pkeys(void)
{
	SpInp_PrintKeys();

	return 0;
}

/*----------------------------------------------------------------------------*/

static int SpTest_Llst(void)
{
	size_t i;
	const char *strings[] = {
		"one",
		"two",
		"three",
		"four",
		"five",
		0
	};
	LNode *list = 0, *lp;
	char *sp;

	for(i = 0; strings[i]; i++) {
		printf("string[%lu]='%s'\n", (unsigned long)i, strings[i]);
		sp = Mem_STRDUP(strings[i]);

		if(1)
			list = Dat_Llst_Push(list, sp, free);
	}

	printf("iterating through list:\n");
	for(lp = list; lp; lp = lp->next) {
		printf("<%p>: %s->%s\n", (void *)lp, lp->name, (char *)lp->value);
	}

	Err_SETSTRING("This is a test error condition");
	Err_SETSTRING("This is yet another test error condition");

	if(Err_Occurred())
		Err_Fprintf(stderr, Sp_parm.debug);

	Err_Clear();

	if(list)
		Dat_Llst_Free(list);

	return 0;
}

/*----------------------------------------------------------------------------*/

#define NDIV ((size_t)2)
#if Sp_MIRSUPPORT
static int SpTest_RayTracing3(void)
{
	SpModel model;
	#if 1
		GeVec3_d min = GeVec3_INIT(0.0, 0.0, 0.0),
			 max = GeVec3_INIT(0.1, PI, TWOPI);
		GeVec3_s ndiv = GeVec3_INIT(10, 1, 1);
		int geom = GEOM_SPH1D;
	#else
		GeVec3_d min = GeVec3_INIT(0.0, 0.0, 0.0),
			 max = GeVec3_INIT(0.1, 0.1, 0.1),
		GeVec3_s ndiv = GeVec3_INIT(10, 10, 10);
		int geom = GEOM_REC3D;
	#endif
	size_t i, side;
	double r, theta, phi, t;
	GeRay ray1, ray2;
	Zone *root, *zp, *zone;
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	GeCam cam = GeCam_INIT(0, 0.0 * PI, 0);

	Mem_BZERO(&model);
	SpModel_InitGrid(&model, geom, min, max, ndiv);

	root = model.grid;

	cpgopen("/xs");
	GeVox_Cpgenv(&model.grid->voxel);
	Zone_Cpgplot(model.grid, 0);
	cpgupdt();

	#define NRAY 100

	for(zone = Zone_GetMinLeaf(model.grid); zone; zone = Zone_AscendTree(zone)) {
		if(zone->children)
			continue;

		for(i = 0; i < NRAY; i++) {
			if(1) {
				ray1 = GeRay_Rand(rng, &zone->voxel);
			}
			else {
				/* Reset ray */
				Mem_BZERO(&ray1);

				/* This samples a random number uniformly in the
				 * interval [0, 1) */
				#define RAND()\
					gsl_rng_uniform(rng)
				/* This samples a random number uniformly in the
				 * interval (0, 1) */
				#define PRAND()\
					gsl_rng_uniform_pos(rng)

				/* Init ray origin to center of voxel */
				ray1.e = zone->voxel.cen;

				r = GeRay_E(ray1, 0);
				theta = PRAND() * PI;
				phi = RAND() * TWOPI;
				//phi = 0.5 * PI;

				#if 0
				Deb_PRINT("r=%g\n", GeRay_E(ray1, 0));
				Deb_PRINT("theta=%g\n", GeRay_E(ray1, 1));
				Deb_PRINT("phi=%g\n", GeRay_E(ray1, 2));
				#endif

				/* Convert to rectangular coordinates */
				GeRay_E(ray1, 0) = r * sin(theta) * cos(phi);
				GeRay_E(ray1, 1) = r * sin(theta) * sin(phi);
				GeRay_E(ray1, 2) = r * cos(theta);

				#if 0
				Deb_PRINT("x=%g\n", GeRay_E(ray1, 0));
				Deb_PRINT("y=%g\n", GeRay_E(ray1, 1));
				Deb_PRINT("z=%g\n", GeRay_E(ray1, 2));
				#endif

				GeVec3_Cpgpt1(&ray1.e, 0, &cam);

				theta = PRAND() * PI;
				phi = RAND() * TWOPI;
				//phi = 0.5 * PI;

				/* Convert to rectangular coordinates */
				GeRay_D(ray1, 0) = sin(theta) * cos(phi);
				GeRay_D(ray1, 1) = sin(theta) * sin(phi);
				GeRay_D(ray1, 2) = cos(theta);
			}

			if(1) {
				zp = zone;

				while(zp) {
					if(0) {
						Deb_PRINT("Traversing zone: [level %d] ", zp->level);
						GeVec3_PRINT(stdout, zp->index);
						printf("\n");
					}

					/* Calculate path to next boundary */
					GeRay_TraverseVoxel(&ray1, &zp->voxel, &t, &side);

					//Deb_PRINT("t=%g, side=%d\n", t, (int)side);

					/* Calculate next ray */
					ray2 = GeRay_Inc(&ray1, t);

					if(1) GeVec3_Cpgarro2(&ray1.e, &ray2.e, &cam);
					if(0) { cpgupdt(); Deb_PAUSE(); }

					/* Get next zone to traverse to */
					zp = Zone_GetNext(zp, &side, &ray2);

					ray1 = ray2;
				}

				if(0) {
					fprintf(stderr, "GAR GAR GAR!!! NO INTERSECTION!!!\n");
					exit(1);
				}

				if(0) {
					cpgupdt();
					Deb_PAUSE();
				}

				
			}
		}
	}

	cpgclos();

	#undef NRAY

	return 0;
}

/*----------------------------------------------------------------------------*/

static int SpTest_RayTracing2(void)
{
	SpModel model;
	#if 1
		GeVec3_d min = GeVec3_INIT(0.0, 0.0, 0.0),
			 max = GeVec3_INIT(0.1, PI, TWOPI),
			 offset = GeVec3_INIT(0.0, 0.0, 0.0);
		GeVec3_s ndiv = GeVec3_INIT(10, 1, 1);
		int geom = GEOM_SPH1D;
	#else
		GeVec3_d min = GeVec3_INIT(0.0, 0.0, 0.0),
			 max = GeVec3_INIT(0.1, 0.1, 0.1),
			 offset = GeVec3_INIT(0.05, 0.05, 0.05);
		GeVec3_s ndiv = GeVec3_INIT(10, 10, 10);
		int geom = GEOM_REC3D;
	#endif
	size_t i, side;
	double theta, phi, t;
	GeRay ray0 = GeRay_INIT(0, 0, 1, 0, 0, -1), ray1, ray2;
	Zone *root, *zp;
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

	Mem_BZERO(&model);
	SpModel_InitGrid(&model, geom, min, max, ndiv);

	root = model.grid;

	cpgopen("/xs");
	GeVox_Cpgenv(&model.grid->voxel);
	Zone_Cpgplot(model.grid, 0);
	cpgupdt();

	#define NRAY 1000

	for(i = 0; i < NRAY; i++) {
		/* This samples a random number uniformly in the
		 * interval [0, 1) */
		#define RAND()\
			gsl_rng_uniform(rng)
		/* This samples a random number uniformly in the
		 * interval (0, 1) */
		#define PRAND()\
			gsl_rng_uniform_pos(rng)

		theta = PRAND() * PI;
		phi = RAND() * TWOPI;

		/* debug */
		if(0) Deb_PRINT("i=%d\n", i);

		/* For every single ray, there *must* be two intersections with the voxel */
		/* Reset ray */
		ray1 = ray0;

		/* Rotate theta rad about y axis */
		ray1 = GeRay_Rotate(&ray1, 0, theta);

		/* Rotate phi rad about z axis */
		ray1 = GeRay_Rotate(&ray1, 2, phi);

		/* Offset ray to direct at center of model */
		ray1.e = GeVec3_Add(&ray1.e, &offset);

		if(0) { GeVec3_Cpgpt1(&ray1.e, 0, NULL); continue; }

		/* Test for intersection with root zone */
		if(GeRay_IntersectVoxel(&ray1, &root->voxel, &t, &side)) {
			/* GeRay hits root zone, locate starting leaf zone for traversal. */
			ray2 = GeRay_Inc(&ray1, t);
			GeVec3_Cpgarro2(&ray1.e, &ray2.e, 0);

			zp = Zone_GetLeaf(root, side, &ray2.e, &ray2);

			if(0) {
				Deb_PRINT("starting from leaf: [level %d] ", zp->level);
				GeVec3_PRINT(stdout, zp->index);
				printf("\n");
			}

			if(0) { cpgupdt(); Deb_PAUSE(); }

			while(zp) {
				if(0) {
					Deb_PRINT("Traversing zone: [level %d] ", zp->level);
					GeVec3_PRINT(stdout, zp->index);
					printf("\n");
				}

				/* Put ray1 at starting position */
				ray1 = ray2;

				/* Calculate path to next boundary */
				GeRay_TraverseVoxel(&ray1, &zp->voxel, &t, &side);

				//Deb_PRINT("t=%g\n", t);

				/* Calculate next ray */
				ray2 = GeRay_Inc(&ray1, t);

				if(1) GeVec3_Cpgarro2(&ray1.e, &ray2.e, 0);
				if(0) { cpgupdt(); Deb_PAUSE(); }

				/* Get next zone to traverse to */
				zp = Zone_GetNext(zp, &side, &ray2);
			}
		}
		else if(1) {
			fprintf(stderr, "GAR GAR GAR!!! NO INTERSECTION!!!\n");
			exit(1);
		}

		if(0) {
			cpgupdt();
			Deb_PAUSE();
		}
	}
	cpgclos();

	#undef NRAY

	return 0;
}

/*----------------------------------------------------------------------------*/
static int SpTest_RandRay(void)
{
	int status = 0;
	SpModel model;
	PyObject *o;
	GeCam cam;
	size_t i, j, nzone = 0, plane;
	GeRay ray;
	GeVec3_d tail;
	gsl_rng *rng;
	double theta, phi, t;
	Zone **zones = 0, *zp;
	SpPhys *pp;

	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, (unsigned long)time(NULL));

	switch(4) {
		case 0: /* No rotation */
			cam = GeCam_Init(0.0, 0.0, 0.0);
			break;

		case 1: /* Rotate 90 deg about Y axis */
			cam = GeCam_Init(0.0, PI / 2.0, 0.0);
			break;

		case 2: /* Rotate 90 deg about Z axis */
			cam = GeCam_Init(0.0, 0.0, PI/2.0);
			break;

		case 3: /* Rotate 45 deg about Y axis, then 45 deg aobut Z axis */
			cam = GeCam_Init(0.0, PI/4.0, PI/4.0);
			break;

		case 4: /* Rotate 45 deg about Y axis, then 45 deg aobut Z axis */
			cam = GeCam_Init(PI / 4.0, PI/4.0, PI/4.0);
			break;
	}

	Mem_BZERO(&model);

	/* Load source model */
	if(!status && !(o = SpInp_GETKEYSTR("source")))
		status = 1;
	if(!status)
		status = SpIO_OpenModel(Sp_PYSTR(o), Sp_PYSTR(o), &model);
	SpPy_XDECREF(o);

	if(!status) {
		for(zp = Zone_GetMinLeaf(model.grid); zp; zp = Zone_AscendTree(zp)) {
			/* Pointer to physical parameters */
			pp = zp->data;

			/* Collect non-empty zones for easy looping and multi-threading */
			if(pp->n_H2 * pp->X_mol > DBL_EPSILON && !zp->children) {
				nzone += 1;
				zones = Mem_REALLOC(nzone, zones);
				zones[nzone - 1] = zp;
			}
		}
	}

	if(!status) {
		/* Print model */
		SpModel_PrintModel(model);

		/* Open PGPLOT device */
		cpgopen("/xs");

		/* Set plot window */
		GeVox_Cpgenv(&model.grid->voxel);
		//Cpg_Env(-1.0, 1.0, -1.0, 1.0, 1, 0);

		/* Plot voxel */
		Zone_Cpgplot(model.grid, &cam);
		//GeVox_Cpgplot(&model.grid->voxel, &cam);

		#define NRAY 5

		for(i = 0; i < nzone; i++) {
			/* Loop through all rays */
			for(j = 0; j < NRAY; j++) {
				zp = zones[i];

				/* Init ray origin to center of voxel */
				ray.e = zp->voxel.cen;

				/* Set random ray direction: first obtrain direction
				 * in spherical coordinates then convert to rectangular coordinates */
				theta = gsl_rng_uniform_pos(rng) * PHYS_CONST_PI;
				phi = gsl_rng_uniform_pos(rng) * PHYS_CONST_TWOPI;

				GeVec3_X(ray.d, 0) = sin(theta) * cos(phi);
				GeVec3_X(ray.d, 1) = sin(theta) * sin(phi);
				GeVec3_X(ray.d, 2) = cos(theta);

				while(zp) {
					tail = ray.e;

					/* Calculate path to next boundary */
					GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &plane);
					ray = GeRay_Inc(&ray, t);

					GeVec3_Cpgarro2(&tail, &ray.e, &cam);

					zp = Zone_GetNext_rec3d(zp, plane, &ray.e);
				} //while(zp)
				cpgupdt();
				Deb_PAUSE();
			} //for(NRAY)
		} //for(nzone)

		#undef NRAY
	} //if(!status)

	/* Cleanup */
	/* Close PGPLOT device */
	cpgclos();
	SpModel_Cleanup(model);
	gsl_rng_free(rng);

	return status;
}
/*----------------------------------------------------------------------------*/

#define NX 128
#define NY 128
#define NV 64

#define CRV (NV / 2)

static int SpTest_MirImg(void)
{
	double deg = PI / 3.14159,
	       min = deg / 60.0,
	       sec = min / 60.0;
	MirImg_Axis x = MirWr_IMGAXIS(NX, sec, 64.0, 0.0),
		    y = MirWr_IMGAXIS(NY, sec, 64.0, 0.0),
		    v = MirWr_IMGAXIS(NV, 10.0, CRV, 0.0);
	MirImg *img;
	double *cube = Mem_CALLOC(NX * NY * NV, cube);
	double xx, yy, zz, width = 30;
	size_t i, j, k;
	MirFile *fp;

	img = MirImg_Alloc(x, y, v);

	for(i = 0; i < NX; i++) {
		for(j = 0; j < NY; j++) {
			for(k = 0; k < NV; k++) {
				xx = (double)i - NX / 2.0;
				yy = (double)j - NY / 2.0;
				zz = (double)k - CRV;
				MirImg_PIXEL(*img, k, i, j) = exp(-1.0 * (xx * xx + yy * yy + zz * zz) / (width * width));
			}
		}
	}

	fp = MirXY_Open_new("test.xy", x.n, y.n, v.n);
	MirImg_WriteXY(fp, img, "K", 1.0);
	MirXY_Close(fp);

	free(cube);

	return 0;
}

#undef NX
#undef NY
#undef NV

/*----------------------------------------------------------------------------*/
static int SpTest_RayAW(void)
{
	Zone *root = 0, *zp;
	GeRay ray, tmp_ray;
	GeCam cam;
	double t, theta, phi;
	size_t plane, i, iter, axis;
	int inc;
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlux);
	GeVec3_s pos;
	GeVec3_d pt;

	gsl_rng_set(rng, (unsigned long)time(NULL));

	/* Setup grid */
	root = SpZone_ALLOC(0);
	root->voxel = GeVox_Init(GEOM_REC3D, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1);
	SpZone_GROW(root, GeVec3_s_Init(NDIV, NDIV, NDIV), 0);

	//cam = GeCam_Init(PI/4.0, PI/4.0, PI/4.0);
	//cam = GeCam_Init(0.25 * PI, 0.25 * PI, 0.25 * PI);
	cam = GeCam_Init(0.0 * PI, 0.0 * PI, 0.0 * PI);

	//print grid
	//SpZone_FPRINTF(stdout, root);

	/* Setup plot */
	cpgopen("/xs");
	//GeVox_Cpgenv(&root->voxel);
	Cpg_Env(-0.2, 0.3, -0.2, 0.3, 1, 0);
	Zone_Cpgplot(root, &cam);

	for(iter = 0; iter < 10; iter++) {
		/* Start ray 300 pc away */
		ray = GeRay_Init(0.2, 0.0, 0.0, -1.0, 0.0, 0.0);

		/* This samples a random number uniformly in the
		 * interval [0, 1) */
		#define RAND()\
			gsl_rng_uniform(rng)
		/* This samples a random number uniformly in the
		 * interval (0, 1) */
		#define PRAND()\
			gsl_rng_uniform_pos(rng)

		theta = asin(2.0 * PRAND() - 1.0) + 0.5 * PI;
		phi = (RAND() * TWOPI);

		if(1) {
			/* Rotate phi rad about z axis */
			ray = GeRay_Rotate(&ray, 2, phi);

			/* Rotate theta rad about y axis */
			ray = GeRay_Rotate(&ray, 1, theta);

			/* Translate to grid-centric coordinates */
			ray.e = GeVec3_Add(&ray.e, &root->voxel.cen);
		}

		/* Calculate first intersection */
		if(GeRay_IntersectVoxel(&ray, &root->voxel, &t, &plane)) {
			printf("================================================\n");
			/* Draw first arrow */
			tmp_ray = GeRay_Inc(&ray, t);
			GeVec3_Cpgarro2(&ray.e, &tmp_ray.e, &cam);
			cpgupdt();
			Deb_PAUSE();

			ray = tmp_ray;

			/* Locate zone of entry */
			zp = Zone_GetLeaf_rec3d(root, plane, &ray.e);

			/* Init ray for AW traversal */
			GeRay_AWInit(&ray, &zp->voxel);
			pos = zp->index;

			printf("entering from zone <%lu,%lu,%lu>\n", (unsigned long)pos.x[0], (unsigned long)pos.x[1], (unsigned long)pos.x[2]);

			while(zp) {
				GeRay_AWTraverse(&ray, &t, &plane);

				printf("Now in zone <%lu,%lu,%lu>\n", (unsigned long)zp->index.x[0], (unsigned long)zp->index.x[1], (unsigned long)zp->index.x[2]);

				if(1) {
					printf("t=%g\n", t);
					for(i = 0; i < 3; i++) {
						printf("tMax[%lu]=%g\n", (unsigned long)i, ray.tMax.x[i]);
						printf("tDelta[%lu]=%g\n", (unsigned long)i, ray.tDelta.x[i]);
					}
					printf("\n");
				}

				tmp_ray = GeRay_Inc(&ray, ray.t);

				GeVec3_Cpgarro2(&ray.e, &tmp_ray.e, &cam);

				if(1) {
					zp = Zone_GetNext(zp, &plane, &ray);
				}
				else if(1) {
					pt = GeRay_AWPos(&ray);
					zp = Zone_GetNext_rec3d(zp, plane, &pt);
				}
				else {
					/* Locate next zone */
					axis = plane / 2;
					inc = (plane % 2 == 0) ? -1 : 1;
					//printf("axis=%d, inc=%d\n", axis, inc);
					//
					printf("moving %d along axis %lu\n", inc, (unsigned long)axis);

					if(inc > 0) {
						if(GeVec3_X(pos, axis) == (GeVec3_X(root->naxes, axis) - 1))
							zp = NULL;
						else {
							GeVec3_X(pos, axis) += (size_t)inc;
							zp = Zone_CHILD(root, pos);
						}
					}
					else {
						if(GeVec3_X(pos, axis) == 0)
							zp = NULL;
						else {
							GeVec3_X(pos, axis) += (size_t)inc;
							zp = Zone_CHILD(root, pos);
						}
					}
				}

				if(!zp)
					printf("exiting from zone <%lu,%lu,%lu>\n", (unsigned long)pos.x[0], (unsigned long)pos.x[1], (unsigned long)pos.x[2]);

				cpgupdt();
				Deb_PAUSE();
			}
			printf("================================================\n");
		}
	}

	cpgclos();

	return 0;
}

/*----------------------------------------------------------------------------*/

static int SpTest_RayTracing(void)
{
	Zone *root = 0, *zp;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, 10, 10, 10);
	GeVec3_s naxes = GeVec3_INIT(NDIV, NDIV, NDIV);
	GeVec3_d offset = GeVec3_INIT(5, 5, 5);
	GeCam cam;
	GeRay ray0 = GeRay_INIT(0, 0, 15, 0, 0, -1), ray1, ray2;
	double t = 0, theta, phi;
	int i;
	size_t plane;
	GeVec3_s cidx = GeVec3_INIT(1, 0, 0);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlux);

	/* Init rng */
	gsl_rng_set(rng, (unsigned long)time(NULL));

	/* Allocate root zone */
	root = Zone_Alloc(0, NULL, 0, 0);
	root->voxel = voxel;

	/* Grow L0 grid */
	if(1)
		Zone_GrowChildren(root, naxes, 0, 0);

	/* Grow L1 grid */
	if(1)
		Zone_GrowChildren(Zone_CHILD(root, cidx), naxes, 0, 0);

	/* Print to stdout */
	if(0)
		Zone_Fprintf(stdout, root, 0);

	/* Setup plot device on X-server */
	cpgopen("/xs");

	if(0)
		GeVox_Cpgenv(&root->voxel);
	else
		Cpg_Env(-15.0, 25.0, -15.0, 25.0, 1, 0);

	/* Setup camera */
	if(1)
		cam = GeCam_Init(0.0, PI/4.0, PI/4.0);
	else
		cam = GeCam_Init(0.0, 0.0, 0.0);

	/* Plot grid */
	if(1)
		Zone_Cpgplot(root, &cam);

	#define NRAY 1000

	for(i = 0; i < NRAY; i++) {
		/* This samples a random number uniformly in the
		 * interval [0, 1) */
		#define RAND()\
			gsl_rng_uniform(rng)
		/* This samples a random number uniformly in the
		 * interval (0, 1) */
		#define PRAND()\
			gsl_rng_uniform_pos(rng)

		theta = PRAND() * PI;
		phi = RAND() * TWOPI;

		/* debug */
		if(0) Deb_PRINT("i=%d\n", i);

		/* For every single ray, there *must* be two intersections with the voxel */
		/* Reset ray */
		ray1 = ray0;

		/* Rotate theta rad about y axis */
		ray1 = GeRay_Rotate(&ray1, 0, theta);

		/* Rotate phi rad about z axis */
		ray1 = GeRay_Rotate(&ray1, 2, phi);

		/* Translate to grid-centric coordinates */
		ray1.e = GeVec3_Add(&ray1.e, &offset);

		if(0) { GeVec3_Cpgpt1(&ray1.e, 0, &cam); continue; }

		/* Test for intersection with root zone */
		if(GeRay_IntersectVoxel(&ray1, &root->voxel, &t, &plane)) {
			/* GeRay hits root zone, locate starting leaf zone for traversal. */
			ray2 = GeRay_Inc(&ray1, t);
			zp = Zone_GetLeaf_rec3d(root, plane, &ray2.e);

			if(0) {
				Deb_PRINT("starting from leaf: [level %d] ", zp->level);
				GeVec3_PRINT(stdout, zp->index);
				printf("\n");
			}

			if(1) GeVec3_Cpgarro2(&ray1.e, &ray2.e, &cam);
			if(0) { cpgupdt(); Deb_PAUSE(); }

			while(zp) {
				if(0) {
					Deb_PRINT("Traversing zone: [level %d] ", zp->level);
					GeVec3_PRINT(stdout, zp->index);
					printf("\n");
				}

				/* Put ray1 at starting position */
				ray1 = ray2;

				/* Calculate path to next boundary */
				GeRay_TraverseVoxel(&ray1, &zp->voxel, &t, &plane);

				/* Calculate next ray */
				ray2 = GeRay_Inc(&ray1, t);

				if(1) GeVec3_Cpgarro2(&ray1.e, &ray2.e, &cam);
				if(0) { cpgupdt(); Deb_PAUSE(); }

				/* Get next zone to traverse to */
				zp = Zone_GetNext_rec3d(zp, plane, &ray2.e);
			}
		}
		else if(1) {
			fprintf(stderr, "GAR GAR GAR!!! NO INTERSECTION!!!\n");
			exit(1);
		}

		if(0) {
			cpgupdt();
			Deb_PAUSE();
		}
	}

	/* Close device */
	cpgclos();

	Zone_Free(root, 0);


	return 0;
}

/*----------------------------------------------------------------------------*/



int SpTest_UVResamp(void)
{
	MirFile *uvin, *uvout;
	MirImg_Axis x, y, v;
	MirImg *img;
	int i, j, k, icen, jcen;
	double width = 5.0, *chan;
	FILE *fp = fopen("SpTest_UVResamp.txt", "w");

	x.n = 16;
	icen = (int)x.n / 2 - 1;
	x.delt = PHYS_UNIT_ASEC;
	y.n = 16;
	jcen = (int)y.n / 2 - 1;
	y.delt = PHYS_UNIT_ASEC;
	v.n = 4;
	v.delt = 100.0;

	img = MirImg_Alloc(x, y, v);
	img->restfreq = 300.0e9; /* Hz */

	for(i = 0; i < (int)img->v.n; i++) {
		for(j = 0; j < (int)img->x.n; j++) {
			for(k = 0; k < (int)img->y.n; k++) {
				MirImg_PIXEL(*img, i, j, k) = exp(-(pow((j - icen) / width, 2.0) + pow((k - jcen) / width, 2.0)));
			}
		}
	}

	for(i = 0; i < (int)img->v.n; i++) {
		chan = MirImg_CHAN(*img, i);
		fprintf(fp, "# Channel %d\n", i);
		for(j = 0; j < (int)img->x.n; j++) {
			for(k = 0; k < (int)img->y.n; k++) {
				fprintf(fp, j == 0 ? "%10.3e" : " %10.3e", chan[(size_t)k + img->y.n * (size_t)j]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n");
	}

	uvin = MirUV_Open(SpInp_GetKey_str("uvin"), "old");
	uvout = MirUV_Open(SpInp_GetKey_str("uvout"), "new");

	MirImg_UVResamp(img, uvin, uvout);

	MirUV_Close(uvin);
	MirUV_Close(uvout);

	return 0;
}
/*----------------------------------------------------------------------------*/
int SpTest_XYReadWrite(void)
{
	MirFile *xy;
	MirImg_Axis x, y, v;
	MirImg *img;
	int i, j, k, icen, jcen;
	double width = 5.0, *chan;
	char *bunit;
	FILE *fp;
	
	Mem_BZERO(&x);
	Mem_BZERO(&y);
	Mem_BZERO(&v);
	x.n = 16;
	icen = (int)x.n / 2 - 1;
	x.delt = PHYS_UNIT_ASEC;
	y.n = 16;
	jcen = (int)y.n / 2 - 1;
	y.delt = PHYS_UNIT_ASEC;
	v.n = 4;
	v.delt = 100.0;
	img = MirImg_Alloc(x, y, v);

	for(i = 0; i < (int)img->v.n; i++) {
		for(j = 0; j < (int)img->x.n; j++) {
			for(k = 0; k < (int)img->y.n; k++) {
				MirImg_PIXEL(*img, i, j, k) = exp(-(pow((j - icen) / width, 2.0) + pow((k - jcen) / width, 2.0)));
			}
		}
	}

	fp = fopen("SpTest_XYReadWrite_original.txt", "w");
	fprintf(fp, "################################################################################\n");
	for(i = 0; i < (int)img->v.n; i++) {
		chan = MirImg_CHAN(*img, i);
		fprintf(fp, "# Channel %d\n", i);
		for(j = 0; j < (int)img->x.n; j++) {
			for(k = 0; k < (int)img->y.n; k++) {
				fprintf(fp, k == 0 ? "%10.3e" : " %10.3e", chan[(size_t)k + img->y.n * (size_t)j]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n");
	}

	fclose(fp);

	/* Write image to Miriad image */
	xy = MirXY_Open_new("SpTest_XYReadWrite.xy", img->x.n, img->y.n, img->v.n);
	MirImg_WriteXY(xy, img, "JY/PIXEL", 1.0);
	MirXY_Close(xy);
	MirImg_Free(img);

	/* Read from Miriad image */
	Mem_BZERO(&x);
	Mem_BZERO(&y);
	Mem_BZERO(&v);
	xy = MirXY_Open_old("SpTest_XYReadWrite.xy", &x.n, &y.n, &v.n);
	img = MirImg_Alloc(x, y, v);
	MirImg_ReadXY(xy, img, &bunit);

	fp = fopen("SpTest_XYReadWrite_fromfile.txt", "w");
	fprintf(fp, "################################################################################\n");
	for(i = 0; i < (int)img->v.n; i++) {
		chan = MirImg_CHAN(*img, i);
		fprintf(fp, "# Channel %d\n", i);
		for(j = 0; j < (int)img->x.n; j++) {
			for(k = 0; k < (int)img->y.n; k++) {
				fprintf(fp, k == 0 ? "%10.3e" : " %10.3e", chan[(size_t)k + img->y.n * (size_t)j]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n");
	}

	return 0;
}
/*----------------------------------------------------------------------------*/

void SpTest_MirImg_LoadXYV(const char *fname)
{
	MirImg_Axis x, y, v;
	double *cube = 0;
	char *bunit;
	size_t i, j, k;

	cube = MirImg_LoadXYV(fname, &x, &y, &v, &bunit);

	for(i = 0; i < x.n; i++) {
		for(j = 0; j < y.n; j++) {
			for(k = 0; k < v.n; k++) {
				printf("cube(%3lu,%3lu,%3lu)=%g\n", (unsigned long)i, (unsigned long)j, (unsigned long)k, cube[k + v.n * (j * y.n * i)]);
			}
		}
	}

	return;
}
/*----------------------------------------------------------------------------*/
#endif
