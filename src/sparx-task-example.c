#include <gsl/gsl_interp.h>
#include "sparx.h"

/* Global parameter struct */
static struct glb {
	int geom, pops;
	double xmol;
	SpModel model;
	SpFile *outf;
} glb;

static void Cleanup(void);
static int GenModel(void);

#define NZONE 16
#define NLEV 21
#define RMAX 1.196800E+15 /* m */
#define TCMB 2.728 /* k */
#define GAS2DUST 100.0

/*----------------------------------------------------------------------------*/

static void Cleanup(void)
{
	SpModel_Cleanup(glb.model);

	if(glb.outf)
		SpIO_CloseFile(glb.outf);

	return;
}

/*----------------------------------------------------------------------------*/

int SpTask_Example(void)
{
	int status = 0;

	/* Reset parameters */
	Mem_BZERO(&glb);

	/* Process inputs */
	glb.xmol = SpInp_GetKey_dbl("xmol");
	glb.geom = SpInp_GetKey_int("geom");
	glb.model.parms.T_cmb = SpInp_GetKey_dbl("tcmb");
	glb.model.parms.gas_to_dust = SpInp_GetKey_dbl("gas2dust");
	glb.pops = SpInp_GetKey_int("pops");

	if(!status && !(glb.outf = SpInp_GetKey_spfile("out", Sp_NEW)))
		status = 1;
			
	/* Generate model */
	if(!status)
		status = GenModel();

	/* Write model to file */
	if(!status)
		status = SpIO_FwriteModel(glb.outf, glb.model);

	if(!status)
		Sp_PRINT("Wrote example model to `%s'\n", glb.outf->name);

	Cleanup();

	return status;
}

/*----------------------------------------------------------------------------*/

static double
	/* In meters */
	ra[NZONE] = { 0.000000E+00, 1.027200E+13, 1.410600E+13, 1.937200E+13, 2.660300E+13, 3.653300E+13, 5.017000E+13, 6.889700E+13, 9.461500E+13, 1.299300E+14, 1.784300E+14, 2.450400E+14, 3.365000E+14, 4.621100E+14, 6.346100E+14, 8.714900E+14 },
	rb[NZONE] = { 1.027200E+13, 1.410600E+13, 1.937200E+13, 2.660300E+13, 3.653300E+13, 5.017000E+13, 6.889700E+13, 9.461500E+13, 1.299300E+14, 1.784300E+14, 2.450400E+14, 3.365000E+14, 4.621100E+14, 6.346100E+14, 8.714900E+14, 1.196800E+15 },
	/* In cm^-3 */
	nh[NZONE] = { 0.000000E+00, 1.177500E+07, 7.316900E+06, 4.546600E+06, 2.825200E+06, 1.825000E+06, 1.184600E+06, 7.467600E+05, 4.979800E+05, 3.365300E+05, 2.298300E+05, 1.608800E+05, 1.144600E+05, 8.281800E+04, 6.013800E+04, 3.211900E+04 },
	/* In K */
	tk[NZONE] = { 1.155700E+02, 8.179100E+01, 7.204400E+01, 6.346000E+01, 5.589800E+01, 4.923700E+01, 4.337000E+01, 3.820200E+01, 3.365000E+01, 2.964000E+01, 2.610800E+01, 2.299700E+01, 2.025700E+01, 1.784300E+01, 1.571700E+01, 1.384400E+01 },
	/* In cm^-3 */
	nm[NZONE] = { 0.000000E+00, 2.355000E-02, 1.463400E-02, 9.093300E-03, 5.650500E-03, 3.649900E-03, 2.369200E-03, 1.493500E-03, 9.959700E-04, 6.730600E-04, 4.596700E-04, 3.217700E-04, 2.289200E-04, 1.656400E-04, 1.202800E-04, 6.423900E-05 },
	/* In km s^-1 */
	vr[NZONE] = { -4.069800E+00, -2.641800E+00, -2.254300E+00, -1.923700E+00, -1.641600E+00, -1.211100E+00, -9.889100E-01, -7.874800E-01, -6.246000E-01, -4.847600E-01, -3.615400E-01, -2.545500E-01, -1.616200E-01, -7.861900E-02, -1.309900E-03, 0.000000E+00 },
	db[NZONE] = { 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01, 1.201200E-01 },
	/* In K */
	//td[NZONE] = { 1.155700E+02, 8.179100E+01, 7.204400E+01, 6.346000E+01, 5.589800E+01, 4.923700E+01, 4.337000E+01, 3.820200E+01, 3.365000E+01, 2.964000E+01, 2.610800E+01, 2.299700E+01, 2.025700E+01, 1.784300E+01, 1.571700E+01, 1.384400E+01 },
	/* Fraction */
	lp[NZONE][NLEV] = {
		{ 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00 },
		{ 2.885190E-02, 8.330146E-02, 1.271923E-01, 1.540454E-01, 1.607279E-01, 1.482958E-01, 1.208282E-01, 8.537473E-02, 5.064383E-02, 2.466711E-02, 1.019875E-02, 3.836174E-03, 1.360526E-03, 4.624391E-04, 1.485214E-04, 4.603741E-05, 1.366917E-05, 3.764840E-06, 1.060774E-06, 2.834690E-07, 4.463728E-08 },
		{ 3.290926E-02, 9.486093E-02, 1.428834E-01, 1.686648E-01, 1.697912E-01, 1.492601E-01, 1.132515E-01, 7.155140E-02, 3.584022E-02, 1.406519E-02, 4.737401E-03, 1.518959E-03, 4.705438E-04, 1.407043E-04, 3.976218E-05, 1.089223E-05, 2.879561E-06, 6.949589E-07, 1.742112E-07, 4.124663E-08, 5.630955E-09 },
		{ 3.804131E-02, 1.098033E-01, 1.631627E-01, 1.869208E-01, 1.792136E-01, 1.458854E-01, 9.783016E-02, 5.088956E-02, 1.969518E-02, 6.092346E-03, 1.769308E-03, 5.071978E-04, 1.400705E-04, 3.695553E-05, 9.187804E-06, 2.216759E-06, 5.184459E-07, 1.085263E-07, 2.371974E-08, 4.881794E-09, 5.618941E-10 },
		{ 4.425079E-02, 1.279979E-01, 1.862953E-01, 2.048445E-01, 1.838207E-01, 1.345220E-01, 7.583151E-02, 3.041626E-02, 8.904131E-03, 2.308918E-03, 6.041860E-04, 1.547374E-04, 3.784986E-05, 8.755418E-06, 1.902354E-06, 3.997424E-07, 8.165781E-08, 1.456430E-08, 2.707685E-09, 4.739667E-10, 4.501453E-11 },
		{ 5.107038E-02, 1.465504E-01, 2.067664E-01, 2.177695E-01, 1.829225E-01, 1.193223E-01, 5.468307E-02, 1.614115E-02, 3.663087E-03, 8.501670E-04, 2.026896E-04, 4.597093E-05, 9.872497E-06, 1.983079E-06, 3.720104E-07, 6.698178E-08, 1.173908E-08, 1.747915E-09, 2.697062E-10, 3.916525E-11, 3.005922E-12 },
		{ 5.940654E-02, 1.685350E-01, 2.288502E-01, 2.279816E-01, 1.750468E-01, 9.735379E-02, 3.384167E-02, 7.209124E-03, 1.397256E-03, 2.989947E-04, 6.367935E-05, 1.255557E-05, 2.330440E-06, 3.991936E-07, 6.359669E-08, 9.640065E-09, 1.418699E-09, 1.730030E-10, 2.158575E-11, 2.515968E-12, 1.529783E-13 },
		{ 6.930146E-02, 1.922629E-01, 2.501625E-01, 2.341505E-01, 1.609660E-01, 7.233311E-02, 1.750262E-02, 2.716315E-03, 4.878731E-04, 9.536260E-05, 1.777593E-05, 3.003222E-06, 4.738069E-07, 6.799994E-08, 9.014166E-09, 1.124034E-09, 1.352747E-10, 1.315563E-11, 1.287148E-12, 1.168339E-13, 5.437550E-15 },
		{ 8.072132E-02, 2.178408E-01, 2.704286E-01, 2.344427E-01, 1.396827E-01, 4.785842E-02, 7.820213E-03, 9.957106E-04, 1.738598E-04, 3.010193E-05, 4.809461E-06, 6.851042E-07, 9.006307E-08, 1.060595E-08, 1.142552E-09, 1.140391E-10, 1.087407E-11, 8.204722E-13, 6.100052E-14, 4.166880E-15, 1.441671E-16 },
		{ 9.289141E-02, 2.425415E-01, 2.878713E-01, 2.297622E-01, 1.153640E-01, 2.815173E-02, 2.998563E-03, 3.498574E-04, 5.924063E-05, 8.843666E-06, 1.187830E-06, 1.399021E-07, 1.497156E-08, 1.411062E-09, 1.201924E-10, 9.295740E-12, 6.782266E-13, 3.845489E-14, 2.094131E-15, 1.027743E-16, 2.101771E-18 },
		{ 1.061413E-01, 2.679053E-01, 3.019953E-01, 2.181149E-01, 8.938619E-02, 1.512923E-02, 1.176638E-03, 1.290256E-04, 1.943202E-05, 2.424144E-06, 2.659313E-07, 2.519938E-08, 2.137434E-09, 1.572083E-10, 1.027780E-11, 5.909443E-13, 3.179438E-14, 1.312100E-15, 4.892530E-17, 6.586714E-19, 1.000000E-30 },
		{ 1.209381E-01, 2.961111E-01, 3.142674E-01, 1.992051E-01, 6.225366E-02, 6.757148E-03, 4.157652E-04, 4.504037E-05, 5.910077E-06, 6.028210E-07, 5.259417E-08, 3.896125E-09, 2.542783E-10, 1.412250E-11, 6.836791E-13, 2.808762E-14, 1.065678E-15, 2.939588E-17, 1.000000E-30, 1.000000E-30, 1.000000E-30 },
		{ 1.375849E-01, 3.255934E-01, 3.207242E-01, 1.735486E-01, 3.956781E-02, 2.813999E-03, 1.501997E-04, 1.510369E-05, 1.658168E-06, 1.345771E-07, 9.056223E-09, 5.076060E-10, 2.464401E-11, 9.941304E-13, 3.417066E-14, 9.579559E-16, 2.367980E-17, 1.000000E-30, 1.000000E-30, 1.000000E-30, 1.000000E-30 },
		{ 1.596636E-01, 3.633808E-01, 3.197070E-01, 1.362249E-01, 2.006407E-02, 9.059228E-04, 4.863944E-05, 4.550882E-06, 4.111055E-07, 2.541727E-08, 1.284767E-09, 5.215416E-11, 1.830106E-12, 5.079548E-14, 1.177247E-15, 2.155812E-17, 4.085761E-19, 5.648188E-20, 4.393560E-20, 3.728379E-20, 3.194698E-20 },
		{ 1.922841E-01, 4.100132E-01, 3.003015E-01, 8.952025E-02, 7.603424E-03, 2.604980E-04, 1.569869E-05, 1.234651E-06, 8.786356E-08, 3.970082E-09, 1.439578E-10, 4.069789E-12, 9.930135E-14, 1.797040E-15, 2.666614E-17, 3.013204E-19, 5.819199E-21, 2.504336E-21, 2.089612E-21, 1.774715E-21, 1.520623E-21 },
		{ 2.628169E-01, 4.724182E-01, 2.274972E-01, 3.565728E-02, 1.555402E-03, 5.147937E-05, 3.242127E-06, 2.047948E-07, 1.108228E-08, 3.458759E-10, 8.409739E-12, 1.605258E-13, 2.620520E-15, 2.904156E-17, 1.000000E-30, 1.000000E-30, 1.000000E-30, 1.000000E-30, 1.000000E-30, 1.000000E-30, 1.000000E-30 },
	};

/*----------------------------------------------------------------------------*/

static int GenModel(void)
{
	int status = 0;
	GeVec3_d min, max, pos;
	GeVec3_s ndiv;
	size_t i, j;
	Zone *grid, *zone;
	SpPhys *pp;
	double
		radius,
		lograd,
		logdum,
		percc = PHYS_UNIT_MKS_PERCC,
		km = PHYS_UNIT_MKS_KM,
		pc = PHYS_UNIT_MKS_PC;
	static double logrb[NZONE], lognh[NZONE], logtk[NZONE], lognm[NZONE], loglp[NLEV][NZONE];

	/* Load hco+ molecule if pops requested */
	if(glb.pops) {
		glb.model.parms.mol = SpIO_FreadMolec("hco+");
		assert(glb.model.parms.mol != NULL);
	}

	switch(glb.geom) {
		case GEOM_SPH1D:
			/* GEOM_SPH1D */
			min = GeVec3_d_Init(0.0, 0.0, 0.0);
			max = GeVec3_d_Init(RMAX/pc, 0.0, 0.0);
			ndiv = GeVec3_s_Init((size_t)NZONE, (size_t)1, (size_t)1);

			/* Allocate grid */
			SpModel_InitGrid(&glb.model, GEOM_SPH1D, min, max, ndiv);
			grid = glb.model.grid;

			for(i = 0; i < NZONE; i++) {
				zone = Zone_CHILD2(grid, i, 0, 0);
				pp = zone->data;

				/* Reset voxel */
				zone->voxel = GeVox_Init(GEOM_SPH1D, ra[i]/pc, 0.0, 0.0, rb[i]/pc, PI, TWOPI);

				/* Set values */
				pp->n_H2 = nh[i] * percc; /* cm^-3 -> m^-3 */
				pp->X_mol = (nh[i] > 0 ? nm[i] / nh[i] : 0.0); /* fraction */
				pp->T_k = tk[i]; /* K */
				pp->V_t = db[i] * km; /* km/s -> m/s */
				GeVec3_X(pp->v_cen, 0) = vr[i] * km; /* km/s -> m/s */

				if(glb.pops) {
					for(j = 0; j < NLEV; j++) {
						pp->pops_preserve[j] = lp[i][j]; /* fraction */
					}
				}
			}
			break;

		case GEOM_REC3D:
			#define NGRID 1
			/* GEOM_REC3D */
			min = GeVec3_d_Init(0.0, 0.0, 0.0);
			max = GeVec3_d_Init(2.0 * RMAX / pc, 2.0 * RMAX / pc, 2.0 * RMAX / pc);
			ndiv = GeVec3_s_Init((size_t)(NGRID*NZONE), (size_t)(NGRID*NZONE), (size_t)(NGRID*NZONE));

			#define LOGZERO -100

			/* Init log values for interpolation */
			for(i = 0; i < NZONE; i++) {
				logrb[i] = rb[i] > 0 ? log10(rb[i]) : LOGZERO;
				lognh[i] = nh[i] > 0 ? log10(nh[i]) : LOGZERO;
				logtk[i] = tk[i] > 0 ? log10(tk[i]) : LOGZERO;
				lognm[i] = nm[i] > 0 ? log10(nm[i]) : LOGZERO;

				for(j = 0; j < NLEV; j++) {
					loglp[j][i] = lp[i][j] > 0 ? log10(lp[i][j]) : LOGZERO;
				}
			}

			/* Allocate grid */
			SpModel_InitGrid(&glb.model, GEOM_REC3D, min, max, ndiv);
			grid = glb.model.grid;

			for(zone = Zone_GetMinLeaf(grid); zone; zone = Zone_AscendTree(zone)) {
				/* Calculate radius from center */
				radius = GeVec3_Mag2(&zone->voxel.cen, &grid->voxel.cen) * pc;

				if(radius <= RMAX) {
					/* Search for position in example grid */
					i = gsl_interp_bsearch(rb, radius, (size_t)0, (size_t)NZONE);

					if(i < 10 && 0)
						Deb_PRINT("radius=%10.3e, ra=%10.3e, rb=%10.3e, i=%2d %s\n", radius / pc, ra[i] / pc, rb[i] / pc, (int)i, i <= 5 ? "<--" : "");

					pp = zone->data;

					if(radius > 0) {
						lograd = log10(radius);

						logdum = Num_InterpPoly(lograd, logrb, lognh, (size_t)NZONE, (size_t)3);
						pp->n_H2 = pow(10.0, logdum) * percc; /* cm^-3 -> m^-3 */

						logdum = Num_InterpPoly(lograd, logrb, lognm, (size_t)NZONE, (size_t)3);
						logdum = pow(10.0, logdum) * percc;
						pp->X_mol = (nh[i] > 0 ? logdum / nh[i] : 0.0); /* fraction */

						logdum = Num_InterpPoly(lograd, logrb, logtk, (size_t)NZONE, (size_t)3);
						pp->T_k = pow(10.0, logdum); /* K */

						pp->V_t = db[i] * km; /* km/s -> m/s */

						/* Calculate position vector of voxel center */
						pos = GeVec3_Sub(&grid->voxel.cen, &zone->voxel.cen);
						pos = GeVec3_Normalize(&pos);
						pp->v_cen = GeVec3_Scale(&pos, vr[i] * km); /* km/s -> m/s */

						if(glb.pops) {
							for(j = 0; j < NLEV; j++) {
								logdum = Num_InterpPoly(lograd, logrb, loglp[j], (size_t)NZONE, (size_t)3);
								pp->pops_preserve[j] = pow(10.0, logdum); /* fraction */
							}
						}
					}
				}
			}
			#undef NGRID
			break;

		default:
			assert(0);
	}

	return status;
}

/*----------------------------------------------------------------------------*/









