#include "sparx.h"
#include <gsl/gsl_rng.h>

#if 0 //esc 09Dec15: failed attempt at dynamic load balancing

enum {
	STAGE_FIX,
	STAGE_RAN,
	STAGE_N
};

/* Global parameter struct */
static struct glb {
	SpFile *outf;
	SpModel model;
	size_t nzone, nray, nconv, maxi, rani, nray_tot;
	Zone **zones;
	size_t *zone_rank, *zone_tid, izone;
	gsl_rng *rng[Sp_NTHREAD];
	unsigned long seed;
	double tolerance, minpop, snr, *I_norm, *I_cmb, max_diff;
	int stage, fully_random, lte;

	pthread_mutex_t exc_mutex, izone_mutex;
} glb;

/* Subroutine prototypes */
int SpTask_Amc(void);
static int InitModel(void);
static void *InitModelThread(void *tid_p);
static int CalcExc(void);
static void *CalcExcThread(void *tid_p);
static void SyncProcs(void);
static void SyncPops(size_t izone);
static void CalcRays(size_t tid, Zone *zone, double *ds0, double *vfac0,
	double *intensity, double *tau);
static void RadiativeXfer(size_t tid, Zone *zone, GeRay *ray, double vel, double *ds0,
	double *vfac0, double *intensity, double *tau);
static void RadiativeXfer2(size_t tid, Zone *zone, GeRay *ray, double vel, double *ds0,
	double *vfac0, double *intensity, double *tau);
static void CalcDetailedBalance(size_t tid, SpPhys *pp, const double *ds0,
	const double *vfac0, const double *intensity, const double *tau);
static void CalcJbar(size_t tid, SpPhys *pp, const double *ds0, const double *vfac0,
	const double *intensity, const double *tau, double *J_bar);
static double CalcDiff(const double *hist, size_t *max_diff_lev);
static double CalcDiff2(const double *hist, size_t *max_diff_lev);
static void Cleanup(void);

/*----------------------------------------------------------------------------*/

/* Some useful definitions */
#define NLEV\
	(glb.model.parms.mol->nlev)
#define NRAD\
	(glb.model.parms.mol->nrad)
#define RAD(i)\
	(glb.model.parms.mol->rad[i])
#define FREQ(i)\
	(glb.model.parms.mol->rad[i]->freq)
#define NCOL\
	(glb.model.parms.mol->ncol)
#define COLL(i)\
	(glb.model.parms.mol->col[i])
#define NCTR(i)\
	(COLL(i)->ntr)
#define CTR(i, j)\
	(COLL(i)->tr[j])
#define CSPECIES(i)\
	(glb.model.parms.mol->col[i]->species)

/*----------------------------------------------------------------------------*/

int SpTask_Amc2(void)
{
	size_t i;
	int status = 0;

	/* Reset parms */
	Mem_BZERO(&glb);

	/* Read keywords */
	glb.nray = SpInp_GetKey_size_t("nrays");
	glb.maxi = SpInp_GetKey_size_t("maxiter");
	glb.rani = SpInp_GetKey_size_t("raniter") - 1;
	glb.tolerance = SpInp_GetKey_dbl("fixlev");
	glb.snr = SpInp_GetKey_dbl("ranlev");
	glb.minpop = SpInp_GetKey_dbl("minpop");
	glb.lte = SpInp_GetKey_TF("lte");
	/* Init RNG seed */
	glb.seed = (unsigned long)time(NULL);

	/* Init RNG */
	for(i = 0; i < Sp_NTHREAD; i++) {
		glb.rng[i] = gsl_rng_alloc(gsl_rng_ranlux);
		gsl_rng_set(glb.rng[i], glb.seed);
	}

	#ifdef HAVE_MPI
	//Sp_PRINT("RANK[%d]: SIZE=%d\n", Sp_MPIRANK, Sp_MPISIZE);
	//Sp_PRINT("-->RANK[%d]: SIZE=%d\n", Sp_parm.mpi_rank, Sp_parm.mpi_size);
	#endif


	/* Load source model */
	if(!status)
		status = SpInp_GetKey_model("source", &glb.model);

	/* Source model must not already contain pops */
	assert(glb.model.parms.mol == NULL);

	/* Load molecule */
	if(!status && !(glb.model.parms.mol = SpInp_GetKey_molec("molec")))
		status = 1;

	/* Open output file handle -- only the master process can write files! */
	if(Sp_MPIRANK == 0) {
		if(!status && !(glb.outf = SpInp_GetKey_spfile("out", Sp_NEW)))
			status = 1;
	}

	/* Initialize model */
	if(!status)
		status = InitModel();

	/* Calculate excitation */
	if(!status)
		status = CalcExc();

	/* Write model to file */
	if(Sp_MPIRANK == 0) {
		if(!status)
			status = SpIO_FwriteModel(glb.outf, glb.model);

		if(!status)
			Sp_PRINT("Wrote output model to `%s'\n", glb.outf->name);
	}

	/* Cleanup */
	Cleanup();

	return status;
}

/*----------------------------------------------------------------------------*/

static void Cleanup(void)
{
	size_t i;

	if(Sp_MPIRANK == 0 && glb.outf)
		SpIO_CloseFile(glb.outf);
	
	SpModel_Cleanup(glb.model);

	if(glb.zones)
		free(glb.zones);

	if(glb.zone_rank)
		free(glb.zone_rank);

	if(glb.zone_tid)
		free(glb.zone_tid);

	for(i = 0; i < Sp_NTHREAD; i++) {
		if(glb.rng[i])
			gsl_rng_free(glb.rng[i]);
	}

	if(glb.I_norm)
		free(glb.I_norm);

	if(glb.I_cmb)
		free(glb.I_cmb);

	return;
}

/*----------------------------------------------------------------------------*/

static int InitModel(void)
{
	Zone *root = glb.model.grid, *zp;
	SpPhys *pp;
	size_t i, nn, nm;

	/* Set normalization intensity to 20K -- normalization prevents rounding
	   errors from creeping in when flux values are very small */
	glb.I_norm = Mem_CALLOC(NRAD, glb.I_norm);
	for(i = 0; i < NRAD; i++) {
		glb.I_norm[i] = Phys_PlanckFunc(FREQ(i), 20.0);
		assert(glb.I_norm[i] > 0); /* Just in case */
	}

	/* Set I_cmb */
	glb.I_cmb = Mem_CALLOC(NRAD, glb.I_cmb);
	if(glb.model.parms.T_cmb > 0) {
		for(i = 0; i < NRAD; i++) {
			glb.I_cmb[i] = Phys_PlanckFunc(FREQ(i), glb.model.parms.T_cmb);
			assert(glb.I_cmb[i] > 0); /* Just in case */
		}
	}

	for(zp = Zone_GetMinLeaf(root); zp; zp = Zone_AscendTree(zp)) {
		/* Pointer to physical parameters */
		pp = zp->data;

		/* Init molecular data */
		SpPhys_InitMol(pp, glb.model.parms.mol);

		//debug
		//if(((pp->n_H2 * pp->X_mol) > 0) && !zp->children) {
		if((pp->n_H2 > 0) && !zp->children) {
			/* This is a non-empty leaf zone */
			pp->non_empty_leaf = 1;

			if(pp->X_mol > 0) {
				/* This zone requires excitation calculations */
				pp->has_tracer = 1;
				glb.nzone += 1;

				/* Collect in zones array */
				glb.zones = Mem_REALLOC(glb.nzone, glb.zones);
				glb.zones[glb.nzone - 1] = zp;
			}
		}
	}

	glb.zone_rank = Mem_CALLOC(glb.nzone, glb.zone_rank);
	glb.zone_tid = Mem_CALLOC(glb.nzone, glb.zone_tid);
	nn = Num_MAX(1, glb.nzone / Sp_MPISIZE);
	nm = Num_MAX(1, nn / Sp_NTHREAD);
	for(i = 0; i < glb.nzone; i++) {
		/* zone_rank and zone_tid are set so that the zones are
		 * distributed as e.g.
		 * NZONE = 8
		 * NPROC = 2
		 * NTHREAD = 2
		 * NN = NZONE / NPROC = 4
		 * NM = NN / NTHREAD = 2
		 * I		0	1	2	3	4	5	6	7
		 * (I/NN)%NPROC	0	0	0	0	1	1	1	1
		 * (I/NM)%NTHRD	0	0	1	1	0	0	1	1
		 */
		/* Calculate corresponding zone_rank */
		glb.zone_rank[i] = (i / nn) % Sp_MPISIZE;

		/* Calculate corresponding zone_tid */
		glb.zone_tid[i] = (i / nm) % Sp_NTHREAD;

		//debug
		/* OR
		 * I		0	1	2	3	4	5	6	7
		 * (I/NN)%NPROC	0	0	0	0	1	1	1	1
		 * (I/NPROC)%NTHRD	0	0	1	1	0	0	1	1
		 */
		//glb.zone_rank[i] = (i / nn) % Sp_MPISIZE;
		//glb.zone_tid[i] = (i / Sp_MPISIZE) % Sp_NTHREAD;

		//debug
		/* OR
		 * I			0	1	2	3	4	5	6	7
		 * i%NPROC		0	1	0	1	0	1	0	1
		 * I%NTHRD	0	1	0	1	0	1	0	1
		 */
		//glb.zone_rank[i] = i % Sp_MPISIZE;
		//glb.zone_tid[i] = i % Sp_NTHREAD;
	}

	SpUtil_Threads(InitModelThread);

	return 0;
}

/*----------------------------------------------------------------------------*/

static void *InitModelThread(void *tid_p)
{
	size_t tid = *((size_t *)tid_p), zone_id, j, k;
	Zone *zp, *root = glb.model.grid;
	SpPhys *pp;

	//for(zone_id = 0; zone_id < glb.nzone; zone_id++) {
	for(zp = Zone_GetMinLeaf(root), zone_id = 0; zp; zp = Zone_AscendTree(zp), zone_id++) {
		/* Skip if zone id and thread id don't match */
		if(zone_id % Sp_NTHREAD != tid)
			continue;

		/* Pointer to physical parameters */
		//zp = glb.zones[zone_id];
		pp = zp->data;

		/* Init RT parameters only if there's gas here */
		if(pp->non_empty_leaf) {
			/* Init excitation parameters only if tracers are present */
			if(pp->has_tracer) {
				/* Set starting number of rays */
				pp->nray = glb.nray;

				/* Calculate thermal line width */
				pp->sigma = SpPhys_CalcLineWidth(pp);

				/* Interpolate downward collisional coeffs and infer upward coeffs */
				SpPhys_InitCollRates(pp);

				/* Set initial pops to either optically thin or LTE */
				for(j = 0; j < pp->mol->nlev; j++) {
					for(k = 0; k < Sp_NTHREAD; k++) {
						if(glb.lte) {
							pp->pops[k][j] = SpPhys_BoltzPops(pp->mol, j, pp->T_k);
						}
						else {
							pp->pops[k][j] = (j == 0) ? 1.0 : 0.0;
						}
					}
				}
			}

			/* Add dust emission/absorption if T_d > 0 */
			if(pp->T_d > 0) {
				SpPhys_AddContinuum_d(pp, 0, glb.model.parms.gas_to_dust);
			}

			/* Add free-free emission/absorption if T_ff > 0 */
			if(pp->T_ff > 0) {
				SpPhys_AddContinuum_ff(pp, 0);
			}

		}
	}

	return NULL;
}

/*----------------------------------------------------------------------------*/

#define MAXRAYS ((size_t)1e8)
#define MINPOPS (glb.minpop)
#define TOLERANCE (glb.tolerance)
#define MCNOISE (1.0 / glb.snr)
#define NHIST ((size_t)3)
#define HIST(i, j)\
	hist[(j) + NLEV * (i)]
#define INTENSITY(i, j)\
	intensity[(j) + NRAD * (i)]
#define TAU(i, j)\
	tau[(j) + NRAD * (i)]

/*----------------------------------------------------------------------------*/

static int CalcExc(void)
{
	int status = 0;
	size_t iter;
	time_t t_start, t_iter;

	/* Make sure exc_mutex is initialized */
	pthread_mutex_init(&glb.exc_mutex, NULL);
	pthread_mutex_init(&glb.izone_mutex, NULL);

	/* Start timer */
	time(&t_start);

	Sp_PRINT("Model geometry is `%s'\n", Geom_CodeToName(glb.model.grid->voxel.geom));
	Sp_PRINT("Solving excitation for %s\n", glb.model.parms.mol->chemname);
	Sp_PRINT("Total %d levels, %d lines\n", glb.model.parms.mol->nlev, glb.model.parms.mol->nrad);
	Sp_PRINT("Beginning convergence from %s conditions\n", glb.lte ? "LTE" : "optically thin");

	for(glb.stage = 0; glb.stage < STAGE_N; glb.stage++) {
		SpPy_CHECKEXC(break);

		/* Sleep for 1 second just in case things went too fast */
		glb.fully_random = (glb.stage == STAGE_FIX) ? 0 : 1;

		/* Get new seed for each stage */
		glb.seed = (unsigned long)time(NULL);

		/* Reset number of converged zones */
		glb.nconv = 0;
		
		/* Some pretty output for the user's convenience */
		if(!glb.fully_random) {
			Sp_PRINT("\n");
			Sp_PRINT("Iterating for convergence with FIXED set of random rays, seed=%lu\n", glb.seed);
			Sp_PRINT("%6s|%15s|%10s|%10s|%10s|%9s\n",
				"Iter.", "Converged/Total", "Prcntg.", "Max diff.", "Goal", "Elapsed");
			Sp_PRINT("------|---------------|----------|----------|----------|---------\n");
		}
		else {
			Sp_PRINT("\n");
			Sp_PRINT("Iterating for convergence with FULLY RANDOM rays, initial seed=%lu\n", glb.seed);
			Sp_PRINT("%6s|%15s|%10s|%10s|%10s|%9s|%20s\n",
				"Iter.", "Converged/Total", "Prcntg.", "Max diff.", "Goal", "Elapsed", "Status");
			Sp_PRINT("------|---------------|----------|----------|----------|---------|--------------------\n");
		}

		for(iter = 0; glb.nconv < glb.nzone; iter++) {
			SpPy_CHECKEXC(break);

			/* Reset global parameters */
			glb.nconv = 0;
			glb.max_diff = 0;
			glb.nray_tot = 0;

			/* Calculate excitation */
			glb.izone = 0;
			SpUtil_Threads(CalcExcThread);

			/* Sync pops from all processes and threads */
			SyncProcs();

			/* Get time for this iteration */
			time(&t_iter);

			/* Pretty output for the user */
			Sp_PRINT("%6g|%7g/%-7g|%9.2f%|%10.4e|%10.4e|%9s",
				(double)(iter + 1),
				(double)glb.nconv,
				(double)glb.nzone,
				100.0 * (double)glb.nconv / (double)glb.nzone,
				glb.max_diff,
				!glb.fully_random ? TOLERANCE : MCNOISE,
				Phys_SecsToHMS_Str((int)difftime(t_iter, t_start)));
				
			if(!glb.fully_random) {
				Sp_PRINTF("\n");
			}
			else {
				if((glb.nconv == glb.nzone)) {
					if(iter < glb.rani) {
						Sp_PRINTF("|      %3d more iters\n", glb.rani - iter);
						glb.nconv = 0;
					}
					else {
						Sp_PRINTF("|%20s\n", "Converged!");
					}
				}
				else {
					Sp_PRINTF("|->%13g rays\n", (double)glb.nray_tot);
				}
			}
		}
	}

	/* Check for Python exceptions */
	if(!status && PyErr_Occurred()) {
		status = 1;
	}

	return status;
}

/*----------------------------------------------------------------------------*/

static void *CalcExcThread(void *tid_p)
{
	size_t tid = *((size_t *)tid_p);
	Zone *zp;
	SpPhys *pp;
	size_t i, ihist, izone;
	double *ds0, *vfac0, *intensity, *tau, *hist, diff, vel;

	hist = Mem_CALLOC(NHIST * NLEV, hist);

	while(glb.izone < glb.nzone) {
		/* Check for Python exceptions */
		SpPy_CHECKEXC(break);

		/* Lock izone mutex */
		pthread_mutex_lock(&glb.izone_mutex);
		izone = glb.izone;

		/* Remember which thread is responsible for this zone */
		glb.zone_tid[izone] = tid;

		/* Update glb.izone */
		glb.izone += 1;

		/* Unlock izone mutex */
		pthread_mutex_unlock(&glb.izone_mutex);

		/* Skip zones that don't belong to this rank */
		if(glb.zone_rank[izone] != Sp_MPIRANK) {
			continue;
		}

		/* Zone related pointers */
		zp = glb.zones[izone];
		pp = zp->data;

		/* Buffers for calculating excitation */
		ds0 = Mem_CALLOC(pp->nray, ds0);
		vfac0 = Mem_CALLOC(pp->nray, vfac0);
		intensity = Mem_CALLOC(pp->nray * NRAD, intensity);
		tau = Mem_CALLOC(pp->nray * NRAD, tau);

		/* Calculate NHIST times for statistics */
		for(ihist = 0; ihist < NHIST; ihist++) {
			if(!glb.fully_random)
				gsl_rng_set(glb.rng[tid], glb.seed);

			if(1) {
				CalcRays(tid, zp, ds0, vfac0, intensity, tau);
			}
			else {
				#define XI()\
					gsl_rng_uniform(glb.rng[tid])
				for(i = 0; i < pp->nray; i++) {
					ds0[i] = XI() * 0.1 * Sp_LENFAC;
					vel = (XI() - 0.5) * 4.3 * pp->sigma;
					vfac0[i] = Num_GaussNormal(vel, pp->sigma);
				}
				#undef XI
			}
			CalcDetailedBalance(tid, pp, ds0, vfac0, intensity, tau);

			for(i = 0; i < NLEV; i++) {
				HIST(ihist, i) = pp->pops[tid][i];
			}
		}
		diff = CalcDiff(hist, NULL);

		/* Lock mutex for global parameters */
		pthread_mutex_lock(&glb.exc_mutex);

		/* Update max_diff */
		if(diff > glb.max_diff)
			glb.max_diff = diff;

		/* Determine convergence and act accordingly */
		if(!glb.fully_random) {
			if(diff <= TOLERANCE) {
				glb.nconv += 1;
			}
		}
		else {
			if(diff <= MCNOISE) {
				glb.nconv += 1;
			}
			else {
				pp->nray *= 4;
				if(pp->nray > MAXRAYS) {
					pp->nray = MAXRAYS;
				}
			}
		}
		/* Update total number of rays */
		glb.nray_tot += pp->nray;

		/* Unlock mutex */
		pthread_mutex_unlock(&glb.exc_mutex);

		free(ds0);
		free(vfac0);
		free(intensity);
		free(tau);
	}

	free(hist);

	//return NULL;
	pthread_exit(NULL);
}

/*----------------------------------------------------------------------------*/

static void SyncProcs(void)
{
	size_t i;

	/* Sync pops from all processes and threads */
	for(i = 0; i < glb.nzone; i++)
		SyncPops(i);

	#ifdef HAVE_MPI
	int glb_nconv, nconv_buff, glb_nray_tot, nray_tot_buff;
	double max_diff_buff[Sp_MPISIZE];

	/* Sync total number of rays */
	glb_nray_tot = (int)glb.nray_tot;
	MPI_Reduce(&glb_nray_tot, &nray_tot_buff, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	glb_nray_tot = (int)nray_tot_buff;
	MPI_Bcast(&glb_nray_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
	glb.nray_tot = (size_t)glb_nray_tot;

	/* Sync total number of converged zones */
	glb_nconv = (int)glb.nconv;
	MPI_Reduce(&glb_nconv, &nconv_buff, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	glb_nconv = nconv_buff;
	MPI_Bcast(&glb_nconv, 1, MPI_INT, 0, MPI_COMM_WORLD);
	glb.nconv = (size_t)glb_nconv;

	/* Gather max_diff from all processes */
	MPI_Gather(&glb.max_diff, 1, MPI_DOUBLE, max_diff_buff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(Sp_MPIRANK == 0) {
		/* Sort max_diff_buff to get largest max_diff */
		Num_Qsort_d(max_diff_buff, Sp_MPISIZE);
		glb.max_diff = max_diff_buff[Sp_MPISIZE - 1];
	}

	/* Sync max_diff with master */
	MPI_Bcast(&glb.max_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	return;
}

/*----------------------------------------------------------------------------*/

static void SyncPops(size_t izone)
{
	size_t zone_tid = glb.zone_tid[izone];
	size_t i;
	Zone *zp;
	SpPhys *pp;

	zp = glb.zones[izone];
	pp = zp->data;

	/* First copy pops from all other threads to thread 0 */
	if(zone_tid != 0) {
		Mem_MEMCPY(pp->pops[0], pp->pops[zone_tid], NLEV);
	}

	#ifdef HAVE_MPI
	/* Then sync results from all processes if using MPI */

	size_t zone_rank = glb.zone_rank[izone];
	MPI_Status mpi_status;

	/* Gather all calculated results to master process */
	if(Sp_MPIRANK != 0) {
		/* This is a slave process; send to master */
		if(zone_rank == Sp_MPIRANK) {
			MPI_Send(pp->pops[0], (int)NLEV, MPI_DOUBLE, 0, Sp_MPITAG, MPI_COMM_WORLD);
			MPI_Send(pp->tau, (int)NRAD, MPI_DOUBLE, 0, Sp_MPITAG, MPI_COMM_WORLD);
			MPI_Send(&pp->ds, 1, MPI_DOUBLE, 0, Sp_MPITAG, MPI_COMM_WORLD);
		}
	}
	else {
		/* This is the master process; receive from a slave */
		if(zone_rank != 0) {
			MPI_Recv(pp->pops[0], (int)NLEV, MPI_DOUBLE, (int)zone_rank, Sp_MPITAG, MPI_COMM_WORLD, &mpi_status);
			MPI_Recv(pp->tau, (int)NRAD, MPI_DOUBLE, (int)zone_rank, Sp_MPITAG, MPI_COMM_WORLD, &mpi_status);
			MPI_Recv(&pp->ds, 1, MPI_DOUBLE, (int)zone_rank, Sp_MPITAG, MPI_COMM_WORLD, &mpi_status);
		}
	}

	/* Sync all processes with master */
	MPI_Bcast(pp->pops[0], (int)NLEV, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(pp->tau, (int)NRAD, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	#endif

	/* Copy thread 0 pops to all other threads */
	for(i = 1; i < Sp_NTHREAD; i++) {
		Mem_MEMCPY(pp->pops[i], pp->pops[0], NLEV);
	}

	return;
}

/*----------------------------------------------------------------------------*/

static void CalcRays(size_t tid, Zone *zone, double *ds0, double *vfac0,
	double *intensity, double *tau)
/* Collect `external' contribution to local mean radiation field (J_bar)
 * by shooting NRAY rays in random directions and calculating the
 * corresponding intensity.
 */
{
	size_t i;
	SpPhys *pp = zone->data;
	double vel;
	GeRay ray;
	GeVec3_d v_gas;

	assert(pp->nray > 0); /* Just in case */

	#if debug_ray
	cpgopen("/xs");
	Vox_Cpgenv(&zone->root->voxel);
	Zone_Cpgplot(zone->root, &cam);
	cpgupdt();
	#endif

	/* Reset pp->ds */
	pp->ds = 0;

	/* Reset tau */
	Mem_BZERO2(tau, pp->nray * NRAD);

	for(i = 0; i < pp->nray; i++) {
		/* Set random ray origin and direction */
		ray = GeRay_Rand(glb.rng[tid], &zone->voxel);

		/* This samples a random number uniformly in the
		 * interval [0, 1) */
		#define RAND()\
			gsl_rng_uniform(glb.rng[tid])

		/* This samples a random number uniformly in the
		 * interval (0, 1) */
		#define PRAND()\
			gsl_rng_uniform_pos(glb.rng[tid])

		/* Set random velocity within local linewidth: PRAND() is used
		 * so that the sampling is uniform within the line width */
		v_gas = SpPhys_GetVgas(&ray.e, zone);
		vel = (PRAND() - 0.5) * 4.3 * pp->sigma + GeVec3_DotProd(&v_gas, &ray.d);

		/* Calculate radiative transfer along this direction */
		//debug
		if(0) {
			RadiativeXfer(tid, zone, &ray, vel, &ds0[i], &vfac0[i], &INTENSITY(i, 0), &TAU(i, 0));
		}
		else {
			RadiativeXfer2(tid, zone, &ray, vel, &ds0[i], &vfac0[i], &INTENSITY(i, 0), &TAU(i, 0));
		}
	}

	#if debug_ray
	cpgclos();
	Deb_PAUSE();
	#endif

	/* Calculate average path length */
	pp->ds /= (double)pp->nray;

	#undef RAND
	#undef PRAND

	#if debug_ray
	cpgclos();
	Deb_PAUSE();
	#endif

	return;
}

/*----------------------------------------------------------------------------*/

static void RadiativeXfer(size_t tid, Zone *zone, GeRay *ray, double vel, double *ds0,
	double *vfac0, double *intensity, double *tau)
/* Given a previously initialized ray, calculate intensity for all lines
 * along the ray from the ray origin to the edge of the cloud */
{
	size_t i;
	Zone *zp = zone;
	SpPhys *pp, *zone_pp;
	double t, ds, vfac, j_nu, k_nu, S_nu, dtau_nu;
	size_t plane;

	/* Reset intensity, ds0 and vfac0 */
	Mem_BZERO2(intensity, NRAD);
	Mem_BZERO2(tau, NRAD);
	*ds0 = 0;
	*vfac0 = 0;

	#if debug_ray //debug
	GeRay tmp_ray;
	size_t iter = 0;
	printf("%12s %12s %12s %12s %12s %12s %12s %12s\n", "iter", "ds", "vfac", "n_H2", "X_mol", "sigma", "dtau", "tau_nu");
	#endif

	/* Reset counter for average path length */
	zone_pp = zone->data;

	/* Propagate the ray through the cloud until we've
	 * reached the edge */
	while(zp) {
		/* Calculate path to next boundary */
		GeRay_TraverseVoxel(ray, &zp->voxel, &t, &plane);

		/* Pointer to physical parameters associated with this zone */
		pp = zp->data;

		//debug
		dtau_nu = 0;

		/* Do calculations on non-empty leaf zones only */
		if(pp->non_empty_leaf) {
			/* Calculate propagation path length and velocity factor
			 * for line profile */
			ds = t * Sp_LENFAC;

			/* Increment average path length counter */
			zone_pp->ds += t;

			/* Calculate velocity line profile factor */
			vfac = pp->has_tracer ? SpPhys_GetVfac(ray, t, vel, zp, 0) : 0;

			if(zp == zone) {
				/* Save ds, vfac for later use and do nothing else */
				*ds0 = ds;
				*vfac0 = vfac;
			}
			else {
				/* Calculate radiative contribution from neighboring
				 * zones */
				for(i = 0; i < NRAD; i++) {
					/* Calculate molecular line emission and absorption coefficients */
					if(pp->has_tracer) {
						SpPhys_GetMoljk(tid, pp, i, vfac, &j_nu, &k_nu);
					}
					else {
						j_nu = 0;
						k_nu = 0;
					}

					/* Add continuum emission/absorption */
					j_nu += pp->cont[i].j;
					k_nu += pp->cont[i].k;

					/* Calculate source function and optical depth if
					 * absorption is NOT zero */
					if(fabs(k_nu) > 0) {
						S_nu = j_nu / k_nu / glb.I_norm[i];
						dtau_nu = k_nu * ds;
						SpPhys_LIMITTAU(dtau_nu);

						/* Calculate intensity contributed by this step */
						intensity[i] += S_nu * (1.0 - exp(-dtau_nu)) * exp(-tau[i]);

						/* Accumulate total optical depth for this line (must be done
						 * AFTER calculation of intensity!) */
						tau[i] += dtau_nu;

						SpPhys_LIMITTAU(tau[i]);
					}
				}
			}
			#if debug_ray
			printf("%12lu %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n", (unsigned long)iter, t, vfac, pp->n_H2, pp->X_mol, pp->sigma, dtau_nu, tau[0]);
			#endif
		}

		#if debug_ray //debug
		tmp_ray = GeRay_Inc(ray, t);
		GeVec3_Cpgarro2(&ray->e, &tmp_ray.e, &cam);
		if(1) {
			//printf("plane=%lu\n", (unsigned long)plane);
			cpgupdt();
			Deb_PAUSE();
		}
		#endif

		/* Propagate ray to next position */
		*ray = GeRay_Inc(ray, t);

		/* Get next zone to traverse to */
		//zp = Zone_Getnext_Rec3d(zp, plane, &ray->e);
		zp = Zone_GetNext(zp, plane, ray);

		#if debug_ray //debug
		iter++;
		#endif
	}

	/* Ray escaped cloud, add CMB to all lines */
	for(i = 0; i < NRAD; i++) {
		intensity[i] += glb.I_cmb[i] * exp(-tau[i]);
	}

	#if debug_ray //debug
	printf("ds0=%12.4e, vfac0=%12.4e\n", (*ds0) / Sp_LENFAC, *vfac0);
	//Deb_PAUSE();
	#endif

	return;
}

/*----------------------------------------------------------------------------*/

static void RadiativeXfer2(size_t tid, Zone *zone, GeRay *ray, double vel, double *ds0,
	double *vfac0, double *intensity, double *tau)
/* Given a previously initialized ray, calculate intensity for all lines
 * along the ray from the ray origin to the edge of the cloud */
{
	size_t i, firststep = 1;
	Zone *zp = zone;
	SpPhys *pp, *zone_pp;
	double t, ds, vfac, j_nu, k_nu, S_nu, dtau_nu;
	size_t plane;

	/* Reset intensity, ds0 and vfac0 */
	Mem_BZERO2(intensity, NRAD);
	Mem_BZERO2(tau, NRAD);
	*ds0 = 0;
	*vfac0 = 0;

	#if debug_ray //debug
	GeRay tmp_ray;
	size_t iter = 0;
	printf("%12s %12s %12s %12s %12s %12s %12s %12s\n", "iter", "ds", "vfac", "n_H2", "X_mol", "sigma", "dtau", "tau_nu");
	#endif

	/* Reset counter for average path length */
	zone_pp = zone->data;

	/* Propagate the ray through the cloud until we've
	 * reached the edge */
	while(zp) {
		/* Calculate path to next boundary */
		GeRay_TraverseVoxel(ray, &zp->voxel, &t, &plane);

		/* Pointer to physical parameters associated with this zone */
		pp = zp->data;

		//debug
		dtau_nu = 0;

		/* Do calculations on non-empty leaf zones only */
		if(pp->non_empty_leaf) {
			/* Calculate propagation path length and velocity factor
			 * for line profile */
			ds = t * Sp_LENFAC;

			/* Increment average path length counter */
			zone_pp->ds += t;

			if(pp->has_tracer) {
				/* Calculate velocity line profile factor */
				vfac = SpPhys_GetVfac(ray, t, vel, zp, 0);

				/* Calculate radiative contribution from neighboring
				 * zones */
				for(i = 0; i < NRAD; i++) {
					/* Calculate molecular line emission and absorption coefficients */
					SpPhys_GetMoljk(tid, pp, i, vfac, &j_nu, &k_nu);

					/* Add continuum emission/absorption */
					j_nu += pp->cont[i].j;
					k_nu += pp->cont[i].k;

					/* Calculate source function and optical depth if
					 * absorption is NOT zero */
					if(fabs(k_nu) > 0) {
						S_nu = j_nu / k_nu / glb.I_norm[i];
					}
					else {
						S_nu = 0;
					}

					dtau_nu = k_nu * ds;
					SpPhys_LIMITTAU(dtau_nu);

					if(!firststep) {
						/* Calculate intensity contributed by this step */
						intensity[i] += S_nu * (1.0 - exp(-dtau_nu)) * exp(-tau[i]);

						/* Accumulate total optical depth for this line (must be done
						 * AFTER calculation of intensity!) */
						tau[i] += dtau_nu;
						SpPhys_LIMITTAU(tau[i]);
					}
				}

				if(firststep) {
					/* Save ds, vfac for later use and do nothing else */
					*ds0 = ds;
					*vfac0 = vfac;
					firststep = 0;
				}
			}
			#if debug_ray
			printf("%12lu %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n", (unsigned long)iter, t, vfac, pp->n_H2, pp->X_mol, pp->sigma, dtau_nu, tau[0]);
			#endif
		}

		#if debug_ray //debug
		tmp_ray = GeRay_Inc(ray, t);
		GeVec3_Cpgarro2(&ray->e, &tmp_ray.e, &cam);
		if(1) {
			//printf("plane=%lu\n", (unsigned long)plane);
			cpgupdt();
			Deb_PAUSE();
		}
		#endif

		/* Propagate ray to next position */
		*ray = GeRay_Inc(ray, t);

		/* Get next zone to traverse to */
		//zp = Zone_Getnext_Rec3d(zp, plane, &ray->e);
		zp = Zone_GetNext(zp, plane, ray);

		#if debug_ray //debug
		iter++;
		#endif
	}

	/* Ray escaped cloud, add CMB to all lines */
	for(i = 0; i < NRAD; i++) {
		intensity[i] += glb.I_cmb[i] * exp(-tau[i]);
	}

	#if debug_ray //debug
	printf("ds0=%12.4e, vfac0=%12.4e\n", (*ds0) / Sp_LENFAC, *vfac0);
	//Deb_PAUSE();
	#endif

	return;
}

/*----------------------------------------------------------------------------*/

static void CalcJbar(size_t tid, SpPhys *pp, const double *ds0, const double *vfac0,
	const double *intensity, const double *tau, double *J_bar)
/* Calculate J_bar, the mean radiation field intensity, by averaging the
 * radiation field over direction and velocity */
{
	size_t i, j;
	double j_nu, k_nu, S_nu, dtau_nu, vfac0_sum;

	/* Reset Jbar (very important!) */
	Mem_BZERO2(J_bar, NRAD);

	/* Reset tau */
	Mem_BZERO2(pp->tau, NRAD);

	/* Reset vfac0_sum */
	vfac0_sum = 0;

	/* Loop through all rays and sum intensities for averaging */
	for(i = 0; i < pp->nray; i++) {
		/* Accumulate vfac_sum */
		vfac0_sum += vfac0[i];

		/* Loop through lines */
		for(j = 0; j < NRAD; j++) {
			/* Calculate local emission and absorption */
			SpPhys_GetMoljk(tid, pp, j, vfac0[i], &j_nu, &k_nu);

			if(fabs(k_nu) > 0) {
				dtau_nu = k_nu * ds0[i];
				S_nu = j_nu / k_nu / glb.I_norm[j];

				SpPhys_LIMITTAU(dtau_nu);

				/* Accumulate intensities weighted over vfac in J_bar: the average
				 * over velocity is extremely important -- failing to do so would
				 * result in a dependence of the excitation on size of grid! */
				J_bar[j] += vfac0[i] * (INTENSITY(i, j) * exp(-dtau_nu) + S_nu * (1.0 - exp(-dtau_nu)));

				/* Store total tau in zone for bookkeeping */
				pp->tau[j] += vfac0[i] * (TAU(i, j) + dtau_nu);
			}
		}
	}

	assert(vfac0_sum >= 0); /* Just in case */

	if(vfac0_sum > 0) {
		for(i = 0; i < NRAD; i++) {
			/* Denormalized and average J_bar */
			J_bar[i] = J_bar[i] * glb.I_norm[i] / vfac0_sum;

			/* Calculate averaged tau for this zone */
			pp->tau[i] /= vfac0_sum;

			#ifdef DEBUG
			J_bar[i] *= 1.5;
			#endif
		}
	}

	return;
}

/*----------------------------------------------------------------------------*/

static void CalcDetailedBalance(size_t tid, SpPhys *pp, const double *ds0,
	const double *vfac0, const double *intensity, const double *tau)
/* Calculate and invert detailed balance equations to solve for level
 * populations, and save results in pp */
{
	size_t iter, ihist, i, j, k, up, lo, max_diff_lev = 0;
	double *rmat, *J_bar, *rhs, *hist, diff = 0;

	/* Allocate J_bar array */
	J_bar = Mem_CALLOC(NRAD, J_bar);

	/* Allocate rates matrix: (NLEV + 1) x NLEV,
	 * where the additional row is for the constraint
	 * that all levels must sum to unity */
	rmat = Mem_CALLOC((NLEV + 1) * NLEV, rmat);

	#define RMAT(i, j)\
		rmat[(j) + NLEV * (i)]
	#define CMAT(i, j)\
		(pp->cmat[(j) + NLEV * (i)])

	/* RHS of rate equation */
	rhs = Mem_CALLOC(NLEV + 1, rhs);
	rhs[NLEV] = 1.0;

	/* Allocate hist array */
	hist = Mem_CALLOC(NHIST * NLEV, hist);

	for(iter = 0; iter < glb.maxi; iter++) {
		for(ihist = 0; ihist < NHIST; ihist++) {
			/* Calculate J_bar, the mean radiation field intensity */
			CalcJbar(tid, pp, ds0, vfac0, intensity, tau, J_bar);

			/* Reset rates matrix */
			Mem_BZERO2(rmat, (NLEV + 1) * NLEV);

			/* Add radiative terms to rates matrix */
			for(i = 0; i < NRAD; i++) {
				up = RAD(i)->up;
				lo = RAD(i)->lo;

				/* Diagonal terms are transitions `out of' row state */
				RMAT(up, up) -= (RAD(i)->A_ul + J_bar[i] * RAD(i)->B_ul);
				RMAT(lo, lo) -= (J_bar[i] * RAD(i)->B_lu);

				/* Off-diagonal terms are transitions `into' row state */
				RMAT(up, lo) += (J_bar[i] * RAD(i)->B_lu);
				RMAT(lo, up) += (RAD(i)->A_ul + J_bar[i] * RAD(i)->B_ul);
			}

			/* Add collisional terms to rates matrix */
			for(i = 0; i < NLEV; i++) {
				for(j = 0; j < NLEV; j++) {
					if(i == j) {
					/* Diagonal terms are minus the sums of the rates from all
					 * collisional transitions `out of' state i */
						for(k = 0; k < NLEV; k++)
							RMAT(i, j) -= CMAT(i, k);
					}
					else {
					/* Off-diagonal terms are sums of rates from state j
					 * `into' state i */
						RMAT(i, j) += CMAT(j, i);
					}
				}
				/* Last row is the constraint that all level densities
				 * must sum to unity */
				RMAT(NLEV, i) = 1.0;
			}

			/* Invert matrix with QR decomposition and solve for
			 * level populations */
			Num_QRDecompSolve(rmat, NLEV + 1, NLEV, rhs, pp->pops[tid]);

			/* Zero out negative/uncertain pops */
			for(i = 0; i < NLEV; i++) {
				if(pp->pops[tid][i] < 0)
					pp->pops[tid][i] = 0.0;
			}

			for(i = 0; i < NLEV; i++) {
				HIST(ihist, i) = pp->pops[tid][i];
			}
		}
		/* Calculate relative difference */
		//diff = CalcDiff(hist, &max_diff_lev);
		diff = CalcDiff2(hist, &max_diff_lev);

		/* Stop iterating if diff is less than TOLERANCE */
		if(diff <= TOLERANCE)
			break;
	}

	/* Warn for non-convergence */
	if(diff > TOLERANCE) {
		Sp_PWARN("non-convergence detected at level %lu in zone <%lu,%lu,%lu> (diff=%.3e)\n",
			max_diff_lev,
			GeVec3_X(pp->zp->index, 0),
			GeVec3_X(pp->zp->index, 1),
			GeVec3_X(pp->zp->index, 2),
			diff);
	}

	#undef RMAT
	#undef CMAT

	/* Cleanup */
	free(J_bar);
	free(rmat);
	free(rhs);
	free(hist);

	return;
}

/*----------------------------------------------------------------------------*/

static double CalcDiff(const double *hist, size_t *max_diff_lev)
{
	size_t i, j;
	double mean, *diffs, max_diff;

	max_diff = 0.0;
	diffs = Mem_CALLOC(NHIST, diffs);

	/* Loop through all levels */
	for(i = 0; i < NLEV; i++) {
		/* Calculate mean for this level */
		mean = 0;
		for(j = 0; j < NHIST; j++) {
			mean += HIST(j, i);
		}
		mean /= (double)NHIST;

		/* Calculate relative difference from mean if
		 * mean >= minpop */
		if(mean >= glb.minpop) {
			for(j = 0; j < NHIST; j++)
				diffs[j] = (HIST(j, i) - mean) / mean;

			/* Find max diff */
			Num_Qsort_d(diffs, NHIST);

			if(diffs[NHIST - 1] > max_diff) {
				max_diff = diffs[NHIST - 1];
				if(max_diff_lev)
					*max_diff_lev = i;
			}
		}
	}

	free(diffs);

	return max_diff;
}

/*----------------------------------------------------------------------------*/

static double CalcDiff2(const double *hist, size_t *max_diff_lev)
{
	size_t i, j;
	double *diffs, *pops, max_diff;

	#define NDIFF (NHIST - 1)
	max_diff = 0.0;
	diffs = Mem_CALLOC(NDIFF, diffs);
	pops = Mem_CALLOC(NHIST, pops);

	/* Loop through all levels */
	for(i = 0; i < NLEV; i++) {
		/* Load pops into array for sorting */
		for(j = 0; j < NHIST; j++) {
			pops[j] = HIST(j, i);
		}
		/* Get smallest pops */
		Num_Qsort_d(pops, NHIST);

		if(pops[0] >= glb.minpop) {
			/* Calculate difference of hist[0] between all
			   other hist values */
			for(j = 0; j < NDIFF; j++) {
				diffs[j] = fabs((HIST(j + 1, i) - HIST(0, i)) / HIST(0, i));
			}
			/* Get smallest diff */
			Num_Qsort_d(diffs, NDIFF);

			if(diffs[NDIFF - 1] > max_diff) {
				max_diff = diffs[NDIFF - 1];
				if(max_diff_lev)
					*max_diff_lev = i;
			}
		}
	}

	#undef NDIFF

	free(diffs);
	free(pops);

	return max_diff;
}

/*----------------------------------------------------------------------------*/

#if 0
static double CalcDiff(const double *hist)
{
	size_t i, j;
	double mean, *diffs, *max_diffs, max_diff;

	diffs = Mem_CALLOC(NHIST, diffs);
	max_diffs = Mem_CALLOC(NLEV, max_diffs);

	/* Loop through all levels */
	for(i = 0; i < NLEV; i++) {
		/* Calculate mean for this level */
		mean = 0;
		for(j = 0; j < NHIST; j++) {
			mean += HIST(j, i);
		}
		mean /= (double)NHIST;

		/* Calculate relative difference from mean if
		 * mean >= minpop */
		if(mean >= glb.minpop) {
			for(j = 0; j < NHIST; j++)
				diffs[j] = (HIST(j, i) - mean) / mean;

			/* Find max diff */
			Num_Qsort_d(diffs, NHIST);
			max_diffs[i] = diffs[NHIST - 1];
		}
		else
			max_diffs[i] = 0.0;
	}
	Num_Qsort_d(max_diffs, NLEV);
	max_diff = max_diffs[NLEV - 1];

	free(diffs);
	free(max_diffs);

	return max_diff;
}
#endif


#endif


