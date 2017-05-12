#include "sparx.h"
#include "miriad.h"
#include <gsl/gsl_rng.h>

static int SpVali_WaterSphere1D(void);
static int SpVali_RayPath(void);
static int SpVali_LVGSphere(void);
static int SpVali_StaticSphere(void);
static int SpVali_StaticCube(void);
static int SpVali_PlaneParallel(void);
static int SpVali_ModelGen(void);
static int SpVali_KepRotate(void);
static int SpVali_PlaneParallelVelo(void);

/* Task definitions */
SpTask
	SpTask_v_watersphere1d = Sp_TASK("v_watersphere1d", "Program for generating the water benchmark problem (1D static)", SpVali_WaterSphere1D, 0),
	SpTask_v_raypath = Sp_TASK("v_raypath", "Program for validating ray path lengths for random rays", SpVali_RayPath, 0),
	SpTask_v_lvgsphere = Sp_TASK("v_lvgsphere", "Program for validating spherical, LVG non-lte excitation", SpVali_LVGSphere, 0),
	SpTask_v_statsphere = Sp_TASK("v_statsphere", "Program for validating spherical, static non-lte excitation", SpVali_StaticSphere, 0),
	SpTask_v_statcube = Sp_TASK("v_statcube", "Program for validating cubical, static non-lte excitation", SpVali_StaticCube, 0),
	SpTask_v_planepar = Sp_TASK("v_planepar", "Program for validating plane-parallel, static RT in 2-lev H2O", SpVali_PlaneParallel, 0),
	SpTask_v_genmodel = Sp_TASK("v_genmodel", "Free-style function for validating model generation in general", SpVali_ModelGen, 0),
	SpTask_v_keprotate = Sp_TASK("t_keprotate", "Generate keplerian disk model for RT validation", SpVali_KepRotate, 0),
	SpTask_v_planeparallelvelo = Sp_TASK("t_planeparallelvelo", "Generate plane-parallel model for RT validation (with systemic velocity)", SpVali_PlaneParallelVelo, 0);

#define NDIV 16
#define DIMS 0.2

/*----------------------------------------------------------------------------*/

static int SpVali_WaterSphere1D(void)
/* The 2-level o-H2O benchmark problem from the 2004 Leiden RT workshop */
{
	#define NZONE 200
	int status = 0;
	PyObject *o;
	double radii[NZONE] = {.100000000E-02, .626483225E-02, .112496854E-01, .159694484E-01, .204382186E-01, .246693436E-01, .286754610E-01, .324685367E-01, .360599000E-01, .394602777E-01, .426798264E-01, .457281623E-01, .486143904E-01, .513471314E-01, .539345477E-01, .563843674E-01, .587039079E-01, .609000973E-01, .629794953E-01, .649483128E-01, .668124302E-01, .685774156E-01, .702485405E-01, .718307966E-01, .733289096E-01, .747473544E-01, .760903674E-01, .773619603E-01, .785659310E-01, .797058755E-01, .807851988E-01, .818071247E-01, .827747054E-01, .836908311E-01, .845582380E-01, .853795169E-01, .861571210E-01, .868933728E-01, .875904713E-01, .882504988E-01, .888754266E-01, .894671213E-01, .900273502E-01, .905577866E-01, .910600148E-01, .915355350E-01, .919857675E-01, .924120570E-01, .928156768E-01, .931978325E-01, .935596654E-01, .939022564E-01, .942266287E-01, .945337512E-01, .948245412E-01, .950998672E-01, .953605517E-01, .956073731E-01, .958410688E-01, .960623368E-01, .962718380E-01, .964701980E-01, .966580095E-01, .968358333E-01, .970042006E-01, .971636142E-01, .973145504E-01, .974574599E-01, .975927697E-01, .977208838E-01, .978421848E-01, .979570352E-01, .980657780E-01, .981687379E-01, .982662225E-01, .983585230E-01, .984459150E-01, .985286595E-01, .986070038E-01, .986811818E-01, .987514151E-01, .988179134E-01, .988808755E-01, .989404892E-01, .989969328E-01, .990503747E-01, .991009746E-01, .991488837E-01, .991942450E-01, .992371940E-01, .992778591E-01, .993163616E-01, .993528166E-01, .993873329E-01, .994200137E-01, .994509566E-01, .994802539E-01, .995079932E-01, .995342574E-01, .995591249E-01, .995826699E-01, .996049629E-01, .996260703E-01, .996460553E-01, .996649774E-01, .996828933E-01, .996998565E-01, .997159175E-01, .997311245E-01, .997455227E-01, .997591553E-01, .997720629E-01, .997842841E-01, .997958554E-01, .998068113E-01, .998171846E-01, .998270063E-01, .998363056E-01, .998451104E-01, .998534470E-01, .998613403E-01, .998688138E-01, .998758898E-01, .998825896E-01, .998889375E-01, .998949433E-01, .999006298E-01, .999060139E-01, .999111116E-01, .999159383E-01, .999205083E-01, .999248353E-01, .999289321E-01, .999328111E-01, .999364838E-01, .999399612E-01, .999432537E-01, .999463711E-01, .999493227E-01, .999521173E-01, .999547633E-01, .999572686E-01, .999596407E-01, .999618867E-01, .999640132E-01, .999660266E-01, .999679329E-01, .999697379E-01, .999714469E-01, .999730650E-01, .999745970E-01, .999760476E-01, .999774211E-01, .999787215E-01, .999799527E-01, .999811185E-01, .999822223E-01, .999832674E-01, .999842569E-01, .999851937E-01, .999860808E-01, .999869207E-01, .999877159E-01, .999884689E-01, .999891818E-01, .999898568E-01, .999905415E-01, .999911429E-01, .999917114E-01, .999922498E-01, .999927596E-01, .999932425E-01, .999936997E-01, .999941328E-01, .999945429E-01, .999949314E-01, .999952993E-01, .999956477E-01, .999959777E-01, .999962903E-01, .999965863E-01, .999968667E-01, .999971323E-01, .999973838E-01, .999977677E-01, .999979928E-01, .999981857E-01, .999983665E-01, .999985387E-01, .999987031E-01, .999988601E-01, .999990103E-01, .999991538E-01, .999992911E-01, .999994225E-01, .999995483E-01, .999996687E-01, .999997841E-01, .999998947E-01, .100000000E+00};
	GeVec3_d min, max;
	GeVec3_s ndiv;
	size_t i;
	Zone *grid, *zone;
	SpPhys *pp;
	double
		vgrad =0,
		xmol = 0,
		percc = PHYS_UNIT_MKS_PERCC;
	SpFile *out_fp = 0;
	SpModel model;

	/* Reset model */
	Mem_BZERO(&model);

	/*-------------------------------------------------------*/

	/* out_fp */
	if(!status && !(o = SpInp_GETKEYSTR("out")))
		status = 1;
	if(!status && !(out_fp = SpIO_OpenFile(Sp_PYSTR(o), Sp_NEW)))
		status = 1;
	SpPy_XDECREF(o);

	/* xmol */
	if(!status && !(o = SpInp_GETKEYSTR("xmol")))
		status = 1;
	if(!status)
		xmol = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/* vgrad */
	if(!status && !(o = SpInp_GETKEYSTR("vgrad")))
		status = 1;
	if(!status)
		vgrad = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/*-------------------------------------------------------*/

	if(!status) {
		/* GEOM_SPH1D */
		min = GeVec3_d_Init(0.0, 0.0, 0.0);
		max = GeVec3_d_Init(1.0, 0.0, 0.0);
		ndiv = GeVec3_s_Init((size_t)NZONE, (size_t)1, (size_t)1);

		/* Allocate grid */
		SpModel_InitGrid(&model, GEOM_SPH1D, min, max, ndiv);
		grid = model.grid;

		for(i = 0; i < NZONE; i++) {
			zone = Zone_CHILD2(grid, i, 0, 0);
			pp = zone->data;

			/* Re-init voxel */
			if(i == 0)
				zone->voxel = GeVox_Init(GEOM_SPH1D, 0.0, 0.0, 0.0, radii[i], PI, TWOPI);
			else
				zone->voxel = GeVox_Init(GEOM_SPH1D, radii[i - 1], 0.0, 0.0, radii[i], PI, TWOPI);
			pp->n_H2 =  1.0e4 * percc; /* cm^-3 -> m^-3 */
			pp->X_mol = xmol; /* fraction */
			pp->T_k = 40.0; /* K */
			GeVec3_X(pp->v_cen, 0) = vgrad * 1000.0 * radii[i]; /* km/s -> m/s */
		}
	
		/* Write model to file */
		status = SpIO_FwriteModel(out_fp, model);
	}

	SpModel_Cleanup(model);

	return status;
}

/*----------------------------------------------------------------------------*/

static void *SpVali_RayPath_Thread(void *tid_p);

static SpModel v_raypath_model;

static int SpVali_RayPath(void)
/* This program is used to validate the total path length traced by a ray by
 * calculating the average path length traced by random rays in a spherical model,
 * which should approach the sphere radius for large numbers of rays */
{
	Zone *zone;
	SpPhys *zone_pp;
	#if 0
		int geom = GEOM_REC3D;
		GeVec3_d min = GeVec3_INIT(0.0, 0.0, 0.0),
		       max = GeVec3_INIT(2, 2, 2);
		size_t n_div = 32;
		GeVec3_s ndiv = GeVec3_INIT(n_div, n_div, n_div);
	#else
		int geom = GEOM_SPH1D;
		GeVec3_d min = GeVec3_INIT(0.0, 0.0, 0.0),
		       max = GeVec3_INIT(0.1, PI, TWOPI);
		size_t n_div = 32;
		GeVec3_s ndiv = GeVec3_INIT(n_div, 1, 1);
	#endif
	SpFile *sfp;

	/* Construct spherical model */
	Mem_BZERO(&v_raypath_model);
	SpModel_InitGrid(&v_raypath_model, geom, min, max, ndiv);
	SpModel_GenModel_UniSphere(&v_raypath_model, 1.0e10, 40.0, 1.0e-8);

	/* Propagate random rays out of the center of the model */
	for(zone = Zone_GetMinLeaf(v_raypath_model.grid); zone; zone = Zone_AscendTree(zone)) {
		zone_pp = zone->data;
		zone_pp->non_empty_leaf = (zone_pp->n_H2 * zone_pp->X_mol > 0.0) && !zone->children;
	}

	SpUtil_Threads(SpVali_RayPath_Thread);

	#ifdef HAVE_MPI
	size_t izone, zone_rank;
	MPI_Status mpi_status;

	for(izone = 0, zone = Zone_GetMinLeaf(v_raypath_model.grid); zone; izone++, zone = Zone_AscendTree(zone)) {
		zone_pp = zone->data;
		zone_rank = (izone % Sp_MPISIZE);

		if(Sp_MPIRANK != 0) {
			/* This is a slave process; send to master */
			if(zone_rank == Sp_MPIRANK) {
				MPI_Send(&zone_pp->ds, 1, MPI_DOUBLE, 0, Sp_MPITAG, MPI_COMM_WORLD);
			}
		}
		else {
			/* This is the master process; receive from slave */
			if(zone_rank != 0) {
				MPI_Recv(&zone_pp->ds, 1, MPI_DOUBLE, (int)zone_rank, Sp_MPITAG, MPI_COMM_WORLD, &mpi_status);
			}
		}
	}
	#endif

	if(Sp_MPIRANK == 0) {
		sfp = SpIO_OpenFile("v_raypath.h5", Sp_TRUNC);
		Deb_ASSERT(sfp != NULL);
		SpIO_FwriteModel(sfp, v_raypath_model);
		SpIO_CloseFile(sfp);
	}

	/* Cleanup */
	SpModel_Cleanup(v_raypath_model);

	return 0;
}

static void *SpVali_RayPath_Thread(void *tid_p)
{
	size_t tid = *((size_t *)tid_p);
	Zone *zp, *zone;
	SpPhys *pp, *zone_pp;
	size_t i, plane, izone;
	GeRay ray;
	double r, theta, phi, t;
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_ranlux);

	gsl_rng_set(rng, (unsigned long)time(NULL));

	for(izone = 0, zone = Zone_GetMinLeaf(v_raypath_model.grid); zone; izone++, zone = Zone_AscendTree(zone)) {
		zone_pp = zone->data;

		/* Skip empty and non-leaf zones */
		if(!zone_pp->non_empty_leaf)
			continue;

		if((izone % Sp_MPISIZE) != Sp_MPIRANK)
			continue;
			
		if(((izone / Sp_MPISIZE) % Sp_NTHREAD) != tid)
			continue;

		#define NRAY 5000

		/* This samples a random number uniformly in the
		 * interval [0, 1) */
		#define RAND()\
			gsl_rng_uniform(rng)
		/* This samples a random number uniformly in the
		 * interval (0, 1) */
		#define PRAND()\
			gsl_rng_uniform_pos(rng)

		for(i = 0; i < NRAY; i++) {
			if(1) {
				ray = GeRay_Rand(rng, &zone->voxel);
			}
			else {
				Mem_BZERO(&ray);

				/* Init random ray origin */
				ray.e = zone->voxel.cen;
				r = GeRay_E(ray, 0);
				theta = PRAND() * PI;
				phi = RAND() * TWOPI;

				/* Convert to rectangular coordinates */
				GeRay_E(ray, 0) = r * sin(theta) * cos(phi);
				GeRay_E(ray, 1) = r * sin(theta) * sin(phi);
				GeRay_E(ray, 2) = r * cos(theta);

				/* Set random ray direction: first obtain direction in spherical
				 * coordinates then convert to rectangular coordinates */

				/* theta must be sampled with PRAND() or there may be a
				 * bias towards the theta=0 direction */
				theta = asin(2.0 * RAND() - 1.0) + 0.5 * PI;

				/* phi must be sampled with RAND() to avoid insufficient
				 * samples in the phi=0 direction */
				phi = PRAND() * TWOPI;

				/* Convert to rectangular coordinates */
				GeRay_D(ray, 0) = sin(theta) * cos(phi);
				GeRay_D(ray, 1) = sin(theta) * sin(phi);
				GeRay_D(ray, 2) = cos(theta);
			}

			zp = zone;

			/* Calculate radiative transfer along this direction */
			/* Propagate the ray through the cloud until we've
			 * reached the edge */
			while(zp) {
				/* Calculate path to next boundary */
				GeRay_TraverseVoxel(&ray, &zp->voxel, &t, &plane);

				/* Pointer to physical parameters associated with this zone */
				pp = zp->data;

				/* Do radiative transfer only if there are molecules in
				 * this zone */
				if(pp->non_empty_leaf) {
					zone_pp->ds += t;
				}

				/* Propagate ray to next position */
				ray = GeRay_Inc(&ray, t);

				/* Get next zone to traverse to */
				zp = Zone_GetNext(zp, plane, &ray);
			}
		}
		zone_pp->ds /= (double)NRAY;
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/

static int SpVali_LVGSphere(void)
{
	size_t ndiv = 0;
	int status = 0;
	PyObject *o = 0;
	SpModel model;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, DIMS, DIMS, DIMS);
	GeVec3_s naxes = GeVec3_INIT(NDIV, NDIV, NDIV);
	Zone *grid, *zp;
	SpPhys *pp;
	SpFile *out_fp = 0;
	double dx, dy, dz, xmol = 0, radius, vel, vgrad = 0;
	char *fname;

	/*-------------------------------------------------------*/

	if(!status && !(o = SpInp_GETKEYSTR("ndiv")))
		status = 1;
	if(!status)
		ndiv = Sp_PYSIZE(o);
	SpPy_XDECREF(o);

	/*-------------------------------------------------------*/

	if(!status && !(o = SpInp_GETKEYSTR("vgrad")))
		status = 1;
	if(!status)
		vgrad = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/*-------------------------------------------------------*/

	if(!status && !(o = SpInp_GETKEYSTR("xmol")))
		status = 1;
	if(!status)
		xmol = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/*-------------------------------------------------------*/

	fname = Mem_Sprintf("v_lvgsphere_xmol_%.0e_ndiv_%02lu_vgrad_%.2e.h5", xmol, ndiv, vgrad);

	if(!status && !(out_fp = SpIO_OpenFile(fname, Sp_NEW)))
		status = 1;

	free(fname);

	/*-------------------------------------------------------*/

	if(status)
		return 1;

	Sp_PRINT("ndiv=%lu\n", ndiv);
	Sp_PRINT("xmol=%g\n", xmol);
	Sp_PRINT("vgrad=%g km/s/pc\n", vgrad);

	GeVec3_X(naxes, 0) = ndiv;
	GeVec3_X(naxes, 1) = ndiv;
	GeVec3_X(naxes, 2) = ndiv;

	/* Reset model */
	Mem_BZERO(&model);

	/* Allocate grid */
	model.grid = grid = SpZone_ALLOC(&model.parms);
	grid->voxel = voxel;
	SpZone_GROW(grid, naxes, &model.parms);

	/* Fill in values */
	for(zp = Zone_GetMinLeaf(grid); zp; zp = Zone_AscendTree(zp)) {
		dx = zp->voxel.cen.x[0] - grid->voxel.cen.x[0];
		dy = zp->voxel.cen.x[1] - grid->voxel.cen.x[1];
		dz = zp->voxel.cen.x[2] - grid->voxel.cen.x[2];

		radius = sqrt(dx * dx + dy * dy + dz * dz);
		vel = vgrad * PHYS_UNIT_MKS_KM;

		//if(1) {
		if(radius <= (DIMS / 2.0)) {
			pp = zp->data;
			pp->n_H2 = 1.0e4 * PHYS_UNIT_MKS_PERCC;
			pp->T_k = 40.0;
			pp->X_mol = xmol;

			pp->v_cen.x[0] = vel * dx;
			pp->v_cen.x[1] = vel * dy;
			pp->v_cen.x[2] = vel * dz;
		}
	}

	/* Write model to file */
	status = SpIO_FwriteModel(out_fp, model);

	return status;
}

/*----------------------------------------------------------------------------*/

static int SpVali_StaticSphere(void)
{
	size_t ndiv = 0;
	int status = 0;
	PyObject *o = 0;
	SpModel model;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, DIMS, DIMS, DIMS);
	GeVec3_s naxes = GeVec3_INIT(NDIV, NDIV, NDIV);
	Zone *grid, *zp;
	SpPhys *pp;
	SpFile *out_fp = 0;
	double dx, dy, dz, xmol = 0, radius;
	char *fname;

	/*-------------------------------------------------------*/

	if(!status && !(o = SpInp_GETKEYSTR("ndiv")))
		status = 1;
	if(!status)
		ndiv = Sp_PYSIZE(o);
	SpPy_XDECREF(o);

	/*-------------------------------------------------------*/

	if(!status && !(o = SpInp_GETKEYSTR("xmol")))
		status = 1;
	if(!status)
		xmol = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/*-------------------------------------------------------*/

	fname = Mem_Sprintf("v_statsphere_xmol_%.0e_ndiv_%lu.h5", xmol, ndiv);

	if(!status && !(out_fp = SpIO_OpenFile(fname, Sp_NEW)))
		status = 1;

	free(fname);

	/*-------------------------------------------------------*/

	if(status)
		return 1;

	Sp_PRINT("ndiv=%lu\n", ndiv);
	Sp_PRINT("xmol=%g\n", xmol);

	GeVec3_X(naxes, 0) = ndiv;
	GeVec3_X(naxes, 1) = ndiv;
	GeVec3_X(naxes, 2) = ndiv;

	/* Reset model */
	Mem_BZERO(&model);

	/* Allocate grid */
	model.grid = grid = SpZone_ALLOC(&model.parms);
	grid->voxel = voxel;
	SpZone_GROW(grid, naxes, &model.parms);

	/* Fill in values */
	for(zp = Zone_GetMinLeaf(grid); zp; zp = Zone_AscendTree(zp)) {
		dx = zp->voxel.cen.x[0] - grid->voxel.cen.x[0];
		dy = zp->voxel.cen.x[1] - grid->voxel.cen.x[1];
		dz = zp->voxel.cen.x[2] - grid->voxel.cen.x[2];

		radius = sqrt(dx * dx + dy * dy + dz * dz);

		#if 0
		if((fabs(dx) <= (DIMS / 4.0)) &&
		   (fabs(dy) <= (DIMS / 4.0)) &&
		   (fabs(dz) <= (DIMS / 4.0))) {
		#endif
		#if 0
		if(1) {
		#endif
		#if 1
		if(radius <= (DIMS / 2.0)) {
		#endif
			pp = zp->data;
			pp->n_H2 = 1.0e4 * PHYS_UNIT_MKS_PERCC;
			pp->T_k = 40.0;
			pp->X_mol = xmol;

			pp->v_cen.x[0] = 0;
			pp->v_cen.x[1] = 0;
			pp->v_cen.x[2] = 0;
		}
	}

	/* Write model to file */
	status = SpIO_FwriteModel(out_fp, model);

	return status;
}

/*----------------------------------------------------------------------------*/

static int SpVali_StaticCube(void)
{
	size_t ndiv = 0;
	int status = 0;
	PyObject *o = 0;
	SpModel model;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, DIMS, DIMS, DIMS);
	GeVec3_s naxes = GeVec3_INIT(NDIV, NDIV, NDIV);
	Zone *grid, *zp;
	SpPhys *pp;
	SpFile *out_fp = 0;
	double dx, dy, dz, xmol = 0, radius;
	char *fname;

	/*-------------------------------------------------------*/

	if(!status && !(o = SpInp_GETKEYSTR("ndiv")))
		status = 1;
	if(!status)
		ndiv = Sp_PYSIZE(o);
	SpPy_XDECREF(o);

	/*-------------------------------------------------------*/

	if(!status && !(o = SpInp_GETKEYSTR("xmol")))
		status = 1;
	if(!status)
		xmol = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/*-------------------------------------------------------*/

	fname = Mem_Sprintf("v_statcube_xmol_%.0e_ndiv_%lu.h5", xmol, ndiv);

	if(!status && !(out_fp = SpIO_OpenFile(fname, Sp_NEW)))
		status = 1;

	free(fname);

	/*-------------------------------------------------------*/

	if(status)
		return 1;

	Sp_PRINT("ndiv=%lu\n", ndiv);
	Sp_PRINT("xmol=%g\n", xmol);

	GeVec3_X(naxes, 0) = ndiv;
	GeVec3_X(naxes, 1) = ndiv;
	GeVec3_X(naxes, 2) = ndiv;

	/* Reset model */
	Mem_BZERO(&model);

	/* Allocate grid */
	model.grid = grid = SpZone_ALLOC(&model.parms);
	grid->voxel = voxel;
	SpZone_GROW(grid, naxes, &model.parms);

	/* Fill in values */
	for(zp = Zone_GetMinLeaf(grid); zp; zp = Zone_AscendTree(zp)) {
		dx = zp->voxel.cen.x[0] - grid->voxel.cen.x[0];
		dy = zp->voxel.cen.x[1] - grid->voxel.cen.x[1];
		dz = zp->voxel.cen.x[2] - grid->voxel.cen.x[2];

		radius = sqrt(dx * dx + dy * dy + dz * dz);

		pp = zp->data;
		//debug
		pp->n_H2 = 1.0e4 * PHYS_UNIT_MKS_PERCC;
		pp->T_k = 40.0;
		pp->X_mol = xmol;

		pp->v_cen.x[0] = 0;
		pp->v_cen.x[1] = 0;
		pp->v_cen.x[2] = 0;
	}

	/* Write model to file */
	status = SpIO_FwriteModel(out_fp, model);

	return status;
}

/*----------------------------------------------------------------------------*/

static int SpVali_PlaneParallel(void)
{
	size_t i;
	int status = 0;
	PyObject *o = 0;
	SpModel model;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, DIMS, DIMS, DIMS);
	GeVec3_s naxes = GeVec3_INIT(NDIV, NDIV, NDIV);
	Zone *grid, *zp;
	SpPhys *pp;
	SpFile *out_fp = 0;
	double dx, dy, dz, xmol = 0;
	char *fname;

	/*-------------------------------------------------------*/

	if(!status && !(o = SpInp_GETKEYSTR("xmol")))
		status = 1;
	if(!status)
		xmol = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	/*-------------------------------------------------------*/

	fname = Mem_Sprintf("v_planepar_xmol_%.0e.h5", xmol);

	if(!status && !(out_fp = SpIO_OpenFile(fname, Sp_NEW)))
		status = 1;

	free(fname);

	/*-------------------------------------------------------*/

	if(status)
		return 1;

	/* Reset model */
	Mem_BZERO(&model);

	/* Load molecule */
	model.parms.mol = SpIO_FreadMolec("o-h2o_2lev");

	/* Allocate grid */
	model.grid = grid = SpZone_ALLOC(&model.parms);
	grid->voxel = voxel;
	SpZone_GROW(grid, naxes, &model.parms);

	/* Fill in values */
	for(zp = Zone_GetMinLeaf(grid); zp; zp = Zone_AscendTree(zp)) {
		dx = zp->voxel.cen.x[0] - grid->voxel.cen.x[0];
		dy = zp->voxel.cen.x[1] - grid->voxel.cen.x[1];
		dz = zp->voxel.cen.x[2] - grid->voxel.cen.x[2];

		if(dx <= 0) {
			pp = zp->data;
			pp->n_H2 = 1.0e4 * PHYS_UNIT_MKS_PERCC;
			pp->T_k = 40.0;
			pp->X_mol = xmol;

			pp->v_cen.x[0] = 0;
			pp->v_cen.x[1] = 0;
			pp->v_cen.x[2] = 0;

			for(i = 0; i < pp->mol->nlev; i++) {
				pp->pops_preserve[i] = SpPhys_BoltzPops(pp->mol, i, pp->T_k);
			}
		}
	}

	/* Write model to file */
	status = SpIO_FwriteModel(out_fp, model);

	return status;
}

/*----------------------------------------------------------------------------*/

static int SpVali_ModelGen(void)
/* Free-style function for validating model generation in general */
{
	size_t i;
	int status = 0;
	PyObject *o = 0;
	SpModel model;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, DIMS, DIMS, DIMS);
	GeVec3_s naxes = GeVec3_INIT(NDIV, NDIV, NDIV);
	Zone *grid, *zp;
	SpPhys *pp;
	SpFile *out_fp = 0;
	double dx, dy, dz, vel, radius;

	if(!status && !(out_fp = SpIO_OpenFile("v_genmodel.h5", Sp_NEW)))
		status = 1;
	SpPy_XDECREF(o);

	/* Reset model */
	Mem_BZERO(&model);

	/* Load molecule */
	model.parms.mol = SpIO_FreadMolec("o-h2o_2lev");

	/* Allocate grid */
	model.grid = grid = SpZone_ALLOC(&model.parms);
	grid->voxel = voxel;
	SpZone_GROW(grid, naxes, &model.parms);

	/* Fill in values */
	for(zp = Zone_GetMinLeaf(grid); zp; zp = Zone_AscendTree(zp)) {
		dx = zp->voxel.cen.x[0] - grid->voxel.cen.x[0];
		dy = zp->voxel.cen.x[1] - grid->voxel.cen.x[1];
		dz = zp->voxel.cen.x[2] - grid->voxel.cen.x[2];

		if((sqrt(dx * dx + dy * dy) <= (DIMS / 2.0)) && (fabs(dz) < DIMS / 8.0)) {
			pp = zp->data;
			pp->n_H2 = 1.0e4 * PHYS_UNIT_MKS_PERCC;
			pp->T_k = 40.0;
			pp->X_mol = 1.0e-8;

			radius = sqrt(dx * dx + dy * dy);
			vel = sqrt(PHYS_CONST_MKS_GRAVG * 5.0 * PHYS_CONST_MKS_MSUN / (radius * Sp_LENFAC));
			pp->v_cen.x[0] = vel * dx / radius;
			pp->v_cen.x[1] = vel * dy / radius;
			pp->v_cen.x[2] = 0;

			pp->v_cen = GeVec3_Rotate_z(&pp->v_cen, PI / 2.0);

			for(i = 0; i < pp->mol->nlev; i++) {
				pp->pops_preserve[i] = SpPhys_BoltzPops(pp->mol, i, pp->T_k);
			}
		}
	}

	/* Write model to file */
	status = SpIO_FwriteModel(out_fp, model);

	return status;
}

#undef NDIV
#undef DIMS

/*----------------------------------------------------------------------------*/

#define NDIV 32
#define DIMS 0.0002

static int SpVali_KepRotate(void)
/* Generate a plane-parallel model with dimensions 0.1x0.1x0.1 pc */
{
	size_t i;
	int status = 0, non_empty;
	PyObject *o = 0;
	SpModel model;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, DIMS, DIMS, DIMS);
	GeVec3_s naxes = GeVec3_INIT(NDIV, NDIV, NDIV);
	Zone *grid, *zp;
	SpPhys *pp;
	double xmol = 0, dx, dy, dz, radius, vel;
	SpFile *out_fp = 0;

	/*-------------------------------------------------------*/
	if(!status && !(o = SpInp_GETKEYSTR("out")))
		status = 1;
	if(!status && !(out_fp = SpIO_OpenFile(Sp_PYSTR(o), Sp_NEW)))
		status = 1;
	SpPy_XDECREF(o);

	/*-------------------------------------------------------*/

	if(!status && !(o = SpInp_GETKEYSTR("xmol")))
		status = 1;
	if(!status)
		xmol = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	if(status)
		return 1;
	/*-------------------------------------------------------*/

	/* Reset model */
	Mem_BZERO(&model);

	/* Load molecule */
	model.parms.mol = SpIO_FreadMolec("o-h2o_2lev");

	/* Allocate grid */
	model.grid = grid = SpZone_ALLOC(&model.parms);
	grid->voxel = voxel;
	SpZone_GROW(grid, naxes, &model.parms);

	/* Fill in values */
	for(zp = Zone_GetMinLeaf(grid); zp; zp = Zone_AscendTree(zp)) {
		non_empty = 0;

		dx = zp->voxel.cen.x[0] - grid->voxel.cen.x[0];
		dy = zp->voxel.cen.x[1] - grid->voxel.cen.x[1];
		dz = zp->voxel.cen.x[2] - grid->voxel.cen.x[2];
		radius = sqrt(dx * dx + dy * dy);

		/* Radius is 0.05 pc */
		if((radius <= DIMS / 2.0) && (fabs(dz) <= DIMS / 8.0)) {
			pp = zp->data;
			pp->n_H2 = 1.0e4 * PHYS_UNIT_MKS_PERCC;
			pp->T_k = 40.0;
			pp->X_mol = xmol;

			if(radius > 0) {
				vel = 1000000.0 * PHYS_CONST_MKS_GRAVG * PHYS_CONST_MKS_MSUN / pow(radius * Sp_LENFAC, 2.0);
				pp->v_cen.x[0] = vel * dx / radius;
				pp->v_cen.x[1] = vel * dy / radius;
				pp->v_cen.x[2] = 0;

				pp->v_cen = GeVec3_Rotate_z(&pp->v_cen, PI / 2.0);
			}

			for(i = 0; i < pp->mol->nlev; i++) {
				pp->pops_preserve[i] = SpPhys_BoltzPops(pp->mol, i, pp->T_k);
			}
		}
	}

	/* Write model to file */
	status = SpIO_FwriteModel(out_fp, model);

	return status;
}

/*----------------------------------------------------------------------------*/

static int SpVali_PlaneParallelVelo(void)
/* Generate a plane-parallel model with dimensions 0.1x0.1x0.1 pc */
{
	size_t i;
	int status = 0;
	PyObject *o = 0;
	SpModel model;
	GeVox voxel = GeVox_INIT(GEOM_REC3D, 0, 0, 0, DIMS, DIMS, DIMS);
	GeVec3_s naxes = GeVec3_INIT(NDIV, NDIV, NDIV);
	Zone *grid, *zp;
	SpPhys *pp;
	double xmol = 0;
	SpFile *out_fp = 0;

	/*-------------------------------------------------------*/
	if(!status && !(o = SpInp_GETKEYSTR("out")))
		status = 1;
	if(!status && !(out_fp = SpIO_OpenFile(Sp_PYSTR(o), Sp_NEW)))
		status = 1;
	SpPy_XDECREF(o);

	/*-------------------------------------------------------*/

	if(!status && !(o = SpInp_GETKEYSTR("xmol")))
		status = 1;
	if(!status)
		xmol = Sp_PYDBL(o);
	SpPy_XDECREF(o);

	if(status)
		return 1;
	/*-------------------------------------------------------*/

	/* Reset model */
	Mem_BZERO(&model);

	/* Load molecule */
	model.parms.mol = SpIO_FreadMolec("o-h2o_2lev");

	/* Allocate grid */
	model.grid = grid = SpZone_ALLOC(&model.parms);
	grid->voxel = voxel;
	SpZone_GROW(grid, naxes, &model.parms);

	/* Fill in values */
	for(zp = Zone_GetMinLeaf(grid); zp; zp = Zone_AscendTree(zp)) {
		pp = zp->data;
		pp->n_H2 = 1.0e4 * PHYS_UNIT_MKS_PERCC;
		pp->T_k = 40.0;
		pp->X_mol = xmol;

		pp->v_cen.x[0] = 1000.0;
		pp->v_cen.x[1] = 0;
		pp->v_cen.x[2] = 0;

		for(i = 0; i < pp->mol->nlev; i++) {
			pp->pops_preserve[i] = SpPhys_BoltzPops(pp->mol, i, pp->T_k);
		}
	}

	/* Write model to file */
	status = SpIO_FwriteModel(out_fp, model);

	return status;
}

#undef DIMS
#undef NDIV

