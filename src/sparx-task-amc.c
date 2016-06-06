#include "sparx.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>

#include <time.h>

enum {
	STAGE_FIX,
	STAGE_RAN,
	STAGE_FIN,
	STAGE_N
};

#define QRAN_DIM 6

/* Global parameter struct */
static struct glb {
	SpFile *outf;
	SpModel model,pops;
	size_t nzone, nray, nconv, maxi, fixi, rani, nray_tot;
	Zone **zones;
	size_t *zone_rank, *zone_tid;
	gsl_rng *rng[Sp_NTHREAD];
        gsl_qrng *qrng[Sp_NTHREAD];
	unsigned long seed;
	double tolerance, minpop, snr, *I_norm, *I_cmb, max_diff, overlap_vel, sor;
	int stage, fully_random, lte, overlap, trace, popsold, qmc, ali, dat;

	pthread_mutex_t exc_mutex;
} glb;

/* Subroutine prototypes */
int SpTask_Amc(void);
static int InitModel(void);
static void *InitModelThread(void *tid_p);
static int CalcExc(void);
static void *CalcExcThread(void *tid_p);
static void SyncProcs(void);
static void SyncPops(size_t izone);
static void CalcRays_RNG(
        size_t tid, 
        Zone *zone, 
        double *ds0, 
        double *vfac0,
	double *intensity, 
        double *tau);
static void CalcRays_QRNG(
        size_t tid, 
        Zone *zone, 
        double *ds0, 
        double *vfac0,
        double *intensity, 
        double *tau);
static void RadiativeXfer(size_t tid, Zone *zone, GeRay *ray, double vel, double *ds0,
	double *vfac0, double *intensity, double *tau);
static void CalcDetailedBalance(size_t tid, SpPhys *pp, const double *ds0,
	const double *vfac0, const double *intensity, const double *tau);
static void CalcJbar(size_t tid, SpPhys *pp, const double *ds0, const double *vfac0,
	const double *intensity, const double *tau, double *J_bar);
static double CalcDiff(const double *hist, size_t *max_diff_lev,size_t izone,size_t tid);
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

#define RELVEL(i,j)\
	(glb.model.parms.mol->OL[NRAD*i+j]->RelativeVel)
#define OVERLAP(i,j)\
	(glb.model.parms.mol->OL[NRAD*i+j]->overlap)
	
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

int SpTask_Amc(void)
{
	int sts = 0;

	/* Reset parms */
	Mem_BZERO(&glb);

	#ifdef HAVE_MPI
	/* Sync all processes first */
	if(Sp_MPISIZE > 1)
		MPI_Barrier(MPI_COMM_WORLD);
	#endif
	/*
	 * Process inputs
	 */
// 1. get parameters from python interface
	if(!sts) sts = SpPy_GetInput_sizt("nrays", &glb.nray);
	if(!sts) sts = SpPy_GetInput_sizt("maxiter", &glb.maxi);
	if(!sts) sts = SpPy_GetInput_sizt("fixiter", &glb.fixi);
	if(!sts) sts = SpPy_GetInput_sizt("raniter", &glb.rani);
	/* Deduct raniter and fixiter by 1 (remove at some point?) */
	if(!sts) glb.fixi -= 1;
	if(!sts) glb.rani -= 1;
	if(!sts) sts = SpPy_GetInput_bool("lte", &glb.lte);
	if(!sts) sts = SpPy_GetInput_bool("trace", &glb.trace);
	//if(!sts) sts = SpPy_GetInput_dbl("tolerance", &glb.tolerance);
        if(!sts) sts = SpPy_GetInput_bool("qmc", &glb.qmc);
        if(!sts) sts = SpPy_GetInput_bool("ali", &glb.ali);
        if(!sts) sts = SpPy_GetInput_bool("dat", &glb.dat);
	if(!sts) sts = SpPy_GetInput_dbl("snr", &glb.snr);
        glb.tolerance = 0.1/glb.snr;
	if(!sts) sts = SpPy_GetInput_dbl("minpop", &glb.minpop);
        if(!sts) sts = SpPy_GetInput_dbl("sor", &glb.sor);
        
        PyObject *o;
	sts = SpPy_GetInput_PyObj("amc", &o);
	PyObject *popsobj = PyObject_GetAttrString(o, "popsold");
	glb.popsold = Sp_PYINT(popsobj);
	PyObject *overlap_obj = PyObject_GetAttrString(o, "overlap");
	glb.overlap_vel = Sp_PYINT(overlap_obj);
        glb.overlap = (glb.overlap_vel == 0.0) ? 0 : 1;

	
	if( glb.popsold ){
		if(!sts) sts = SpPy_GetInput_model("source","pops", &glb.model);
	}
	else{
		if(!sts) sts = SpPy_GetInput_model("source","source", &glb.model);
	}
	/* Source model must not already contain pops */  // hidden by I-Ta 2012.10.26


	if(!sts && glb.model.parms.mol != NULL && glb.popsold==0) {
		PyWrErr_SetString(PyExc_Exception, "Model must not already contain pops");
		sts = 1;
	}

	
	/* Load molecule */
	if( glb.popsold==0 ){	
		if(!sts && glb.model.parms.mol == NULL) sts = SpPy_GetInput_molec("molec", &glb.model.parms.mol);
		//if(!sts && glb.model.parms.mol != NULL) sts = SpPy_GetInput_molec_hyper("hyperfine", &glb.model.parms.mol_hyper);
	}
	/* Open output file handle -- only the master process can write files! */
	if(Sp_MPIRANK == 0 && !sts) sts = SpPy_GetInput_spfile("out", &glb.outf, Sp_NEW);
	/*
	 * Initializtion
	 */
	/* Init RNG seed */
	glb.seed = (unsigned long)time(NULL);
	
	/* Init RNG */
	for(size_t i = 0; i < Sp_NTHREAD; i++) {
		glb.rng[i] = gsl_rng_alloc(gsl_rng_ranlux);
		gsl_rng_set(glb.rng[i], glb.seed);
	}
	
#if 0
	int n=1024;
        int dim=6;
	double v[dim];
        
        char filename[32];
        sprintf(filename,"qmc.dat");
        FILE *fp=fopen(filename,"w+");
        gsl_qrng * g = gsl_qrng_alloc(gsl_qrng_niederreiter_2, dim);
        gsl_qrng_get(g, v);
        for (int i = 0; i < n; i++){
                gsl_qrng_get(g, v);
                for (int j = 0; j < dim; j++)
                        fprintf(fp, "%f ", v[j]);
                fprintf(fp,"\n");
        }
        gsl_qrng_free(g);
        fclose(fp);
        
	exit(0);
#endif	
	

	/* Init model */
	if(!sts)
		sts = InitModel();
	
	/* Calculate excitation */       
	if(!sts)
		sts = CalcExc();
        
	/* Write model to file */
	if(Sp_MPIRANK == 0) {
		if(!sts)
			sts = SpIO_FwriteModel(glb.outf, glb.model);
		if(!sts)
			Sp_PRINT("Wrote output model to `%s'\n", glb.outf->name);
	}

	/* Cleanup */
	Cleanup();

	return sts;
}

/*----------------------------------------------------------------------------*/

static void Cleanup(void)
{
	if(Sp_MPIRANK == 0 && glb.outf)
		SpIO_CloseFile(glb.outf);
	
	SpModel_Cleanup(glb.model);

	if(glb.zones)
		free(glb.zones);

	if(glb.zone_rank)
		free(glb.zone_rank);

	if(glb.zone_tid)
		free(glb.zone_tid);

	for(size_t i = 0; i < Sp_NTHREAD; i++) {
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

	/* Set normalization intensity to 20K -- normalization prevents rounding
	   errors from creeping in when flux values are very small */
	glb.I_norm = Mem_CALLOC(NRAD, glb.I_norm);
	for(size_t i = 0; i < NRAD; i++) {
		glb.I_norm[i] = Phys_PlanckFunc(FREQ(i), 20.0);
		Deb_ASSERT(glb.I_norm[i] > 0); /* Just in case */
	}

	/* Set I_cmb -- DO NOT FORGET TO NORMALIZE! */
	glb.I_cmb = Mem_CALLOC(NRAD, glb.I_cmb);
	if(glb.model.parms.T_cmb > 0) {
		for(size_t i = 0; i < NRAD; i++) {
			glb.I_cmb[i] = Phys_PlanckFunc(FREQ(i), glb.model.parms.T_cmb) / glb.I_norm[i];
			Deb_ASSERT(glb.I_cmb[i] > 0); /* Just in case */
		}
	}
	


	/* construct overlapping table */
	if(glb.overlap){
		glb.model.parms.mol->OL = Mem_CALLOC(NRAD*NRAD,glb.model.parms.mol->OL);
		for(size_t i=0; i<NRAD; i++){
			for(size_t j=0; j<NRAD; j++){
				glb.model.parms.mol->OL[NRAD*i+j] = Mem_CALLOC(1,glb.model.parms.mol->OL[NRAD*i+j]);
				RELVEL(i,j) = ( 1e0-FREQ(j)/FREQ(i) )*CONSTANTS_MKS_LIGHT_C;
				if( fabs(RELVEL(i,j)) < glb.overlap_vel ){
					OVERLAP(i,j)=1;
				}
				else{
					OVERLAP(i,j)=0;
				}
// 				printf("%2d",OVERLAP(i,j));
			}
// 			printf("\n");
		}
	}

	
	for(zp = Zone_GetMinLeaf(root); zp; zp = Zone_AscendTree(zp)) {
		/* Pointer to physical parameters */
		SpPhys *pp = zp->data;
		
		/* Init molecular data */
		SpPhys_InitMol(pp, glb.model.parms.mol,glb.popsold);

		if((pp->n_H2 > 1e-200) && !zp->children) {
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
			else{
				pp->has_tracer = 0;	
			}
		}
		else{
			pp->non_empty_leaf = 0;
			pp->has_tracer = 0;
		}
	}

	glb.zone_rank = Mem_CALLOC(glb.nzone, glb.zone_rank);
	glb.zone_tid = Mem_CALLOC(glb.nzone, glb.zone_tid);
	size_t nn = Num_MAX(1, glb.nzone / Sp_MPISIZE);
	size_t nm = Num_MAX(1, nn / Sp_NTHREAD);
	for(size_t i = 0; i < glb.nzone; i++) {
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
	

	SpUtil_Threads2(Sp_NTHREAD, InitModelThread);

	return 0;
}

/*----------------------------------------------------------------------------*/

static void *InitModelThread(void *tid_p)
{
	size_t tid = *((size_t *)tid_p);
	Zone *root = glb.model.grid;	

	//for(zone_id = 0; zone_id < glb.nzone; zone_id++) {
        Zone *zp;
        size_t zone_id;
	for(zp = Zone_GetMinLeaf(root), zone_id = 0; zp; zp = Zone_AscendTree(zp), zone_id++) {
		/* Skip if zone id and thread id don't match */
		if(zone_id % Sp_NTHREAD != tid)
			continue;

		/* Pointer to physical parameters */
		//zp = glb.zones[zone_id];
		SpPhys *pp = zp->data;//printf("OK \n");exit(0);

		/* Init RT parameters only if there's gas here */
		if(pp->non_empty_leaf) {
			/* Init excitation parameters only if tracers are present */
			if(pp->has_tracer) {
				/* Set starting number of rays */
				pp->nray = glb.nray;
				pp->diff = 0.0;

				/* Calculate thermal line width */
				pp->width = SpPhys_CalcLineWidth(pp);

				/* Interpolate downward collisional coeffs and infer upward coeffs */
				SpPhys_InitCollRates(pp);

				/* Set initial pops to either optically thin or LTE */
				if (glb.popsold == 0){   // only if not using old pops (added by I-Ta 2012.10.26)
					for(size_t j = 0; j < pp->mol->nlev; j++) {
						for(size_t k = 0; k < Sp_NTHREAD; k++) {
							if(glb.lte) {
								pp->pops[k][j] = SpPhys_BoltzPops(pp->mol, j, pp->T_k);
							}
							else {
								pp->pops[k][j] = (j == 0) ? 1.0 : 0.0;
							}
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
	Zone *zone = glb.model.grid;
	/* set the minimal converged cell number for the fully-random stage */
	size_t ran_min_nconv = (size_t)(0.995*(double)glb.nzone+0.5);
	
	/* Make sure exc_mutex is initialized */
	pthread_mutex_init(&glb.exc_mutex, NULL);

	/* Start timer */
        time_t t_start;
	time(&t_start);

	Sp_PRINT("Model geometry is `%s'\n", Geom_CodeToName(glb.model.grid->voxel.geom));
	Sp_PRINT("Solving excitation for %s\n", glb.model.parms.mol->chemname);
	Sp_PRINT("Total %d levels, %d lines\n", glb.model.parms.mol->nlev, glb.model.parms.mol->nrad);
	Sp_PRINT("Beginning convergence from %s conditions\n", glb.lte ? "LTE" : "GROUND STATE");

        int sts = 0;
	for(glb.stage = 0;
            !sts && (glb.stage < (glb.ali ? 1 : STAGE_N)); 
            glb.stage++) {
		
		/* Sleep for 1 second just in case things went too fast */
		glb.fully_random = (glb.stage == STAGE_RAN) ? 1 : 0;

                if(glb.qmc){
                        if(glb.fully_random){
                                // reset QMC random seed
                                for (size_t i =0; i < Sp_NTHREAD; i++){
                                        glb.qrng[i] = gsl_qrng_alloc(gsl_qrng_niederreiter_2, QRAN_DIM);
                                        // skip the initial zero-array 
                                        double QRN[QRAN_DIM];
                                        gsl_qrng_get(glb.qrng[i], QRN);
                                }
                        }
                }
                else{
                        /* Get new seed for each stage */
                        glb.seed = (unsigned long)time(NULL);
                }
                
		/* Reset number of converged zones */
		glb.nconv = 0;
		
		/* Some pretty output for the user's convenience */
		if(!glb.fully_random) {
			Sp_PRINT("\n");
                        if(glb.qmc)
                                Sp_PRINT("Iterating for convergence with FIXED set of random rays, Quasi Random Ray\n");
                        else
                                Sp_PRINT("Iterating for convergence with FIXED set of random rays, seed=%lu\n", glb.seed);
		}
		else {
			Sp_PRINT("\n");
                        if(glb.qmc)
                                Sp_PRINT("Iterating for convergence with FULLY RANDOM rays,  Quasi Random Ray\n");
                        else
                                Sp_PRINT("Iterating for convergence with FULLY RANDOM rays,  initial seed=%lu\n", glb.seed);
		}

		Sp_PRINT("%6s|%15s|%10s|%10s|%10s|%9s|%20s\n", "Iter.", "Converged/Total", "Prcntg.", "Max diff.", "Goal", "Elapsed", "Status");
		Sp_PRINT("------|---------------|----------|----------|----------|---------|--------------------\n");
		

		//for(iter = 0; !sts && (glb.nconv < glb.nzone ); iter++) {
		for (size_t iter = 0; 
                     !sts && (glb.nconv < ( glb.fully_random ? ran_min_nconv : glb.nzone ) );
                     iter++){
			/* Check for thread termination */
			if(SpUtil_TermThread()) {
				sts = 1;
				break;
			}

			/* Reset global parameters */
			glb.nconv = 0;
			glb.max_diff = 0;
			glb.nray_tot = 0;

			/* Calculate excitation */
			SpUtil_Threads2(Sp_NTHREAD, CalcExcThread);

			/* Sync pops from all processes and threads */
			SyncProcs();

			/* Get time for this iteration */
                        time_t t_iter;
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
				
			Sp_PRINTF("|->%13g rays\n", (double)glb.nray_tot);
                        if (glb.nconv >= ((glb.fully_random) ? ran_min_nconv:glb.nzone)) {
				if(iter < (glb.fully_random ? glb.rani : glb.fixi)) {
					Sp_PRINTF("%5d more iters\n", (glb.fully_random ? glb.rani : glb.fixi) - iter);
					glb.nconv = 0;
				}
				else {
					Sp_PRINTF("%16s\n", "Converged!");
				}
			}


			#ifdef HAVE_MPI
			if(Sp_MPISIZE > 1)
				MPI_Barrier(MPI_COMM_WORLD);
			#endif

			/* Write out model if requested by user to trace convergence history */
			if(glb.trace && Sp_MPIRANK == 0) {
                                SpFile *tracefp = NULL;
				if(!sts) {
					char *stmp = Mem_Sprintf("%s.stage%d-%05d", glb.outf->name, glb.stage, iter);
					sts = SpIO_OpenFile2(stmp, Sp_NEW, &tracefp);
					free(stmp);
				}
				if(!sts) sts = SpIO_FwriteModel(tracefp, glb.model);
				if(!sts) SpIO_CloseFile(tracefp);
			}
			
			/* Redistribute parallelization by weighting nray (added by I-Ta)*/
			if(glb.stage == STAGE_RAN){
				size_t temp_nray_tot = 0;
				for(size_t izone = 0; izone < glb.nzone; izone++) {
					Zone *zp = glb.zones[izone];
					SpPhys *pp = zp->data;
					temp_nray_tot+=pp->nray;
				}
				if(temp_nray_tot!=glb.nray_tot){
					if(Sp_MPIRANK == 0){
						printf("total number of ray isn't consistent: %13g\n", (double)temp_nray_tot);
					}
					glb.nray_tot=temp_nray_tot;
				}
				temp_nray_tot=0;
				size_t ithread=0;
				double average_nray = (double)glb.nray_tot/(double)(Sp_MPISIZE*Sp_NTHREAD);
				double critical_nray = average_nray;		
				for(size_t izone = 0; izone < glb.nzone; izone++) {
					Zone *zp = glb.zones[izone];
					SpPhys *pp = zp->data;
					size_t temp_old = temp_nray_tot;
					temp_nray_tot=temp_nray_tot+pp->nray;
					if ( (double)temp_nray_tot > critical_nray ){	
						if( (double)temp_nray_tot-critical_nray > critical_nray-(double)temp_old ){
							ithread=ithread+1;
							glb.zone_tid[izone] = ithread % Sp_NTHREAD;
							glb.zone_rank[izone] = (ithread/Sp_NTHREAD) % Sp_MPISIZE;
						}
						else{
							glb.zone_tid[izone] = ithread % Sp_NTHREAD;
							glb.zone_rank[izone] = (ithread/Sp_NTHREAD) % Sp_MPISIZE;
							ithread=ithread+1;
						}
						critical_nray=critical_nray+average_nray;
					}
					else{
						glb.zone_tid[izone] = ithread % Sp_NTHREAD;
						glb.zone_rank[izone] = (ithread/Sp_NTHREAD) % Sp_MPISIZE;
					}
				}
			}
#if 0
			/* Write out convergent information for visualization */
			if(glb.stage==1 && Sp_MPIRANK == 0){
                                char filename[13];
				sprintf(filename,"SPARX%4.4zu.vtk",iter);
				FILE *fp = fopen(filename,"w");
				fprintf(fp,"# vtk DataFile Version 3.0\n");
				fprintf(fp,"SPARX\n");
				fprintf(fp,"ASCII\n");
				fprintf(fp,"DATASET STRUCTURED_GRID\n");
				fprintf(fp,"DIMENSIONS %zu %zu %zu \n",glb.model.grid->naxes.x[0],glb.model.grid->naxes.x[1],glb.model.grid->naxes.x[2]);
				//fprintf(fp,"POINTS %zu float\n",glb.nzone);
				fprintf(fp,"POINTS %zu float\n",zone->nchildren);
				//for(size_t izone = 0; izone < glb.nzone; izone++){					
				for(size_t izone=0; izone < zone->nchildren; izone++ ){
					//zp = glb.zones[izone];
					Zone *zp = zone->children[izone];
					fprintf(fp,"%11.4e %11.4e %11.4e\n",zp->voxel.cen.x[0],zp->voxel.cen.x[1],zp->voxel.cen.x[2]);
				}
				//fprintf(fp,"POINT_DATA %zu\n",glb.nzone);
				fprintf(fp,"POINT_DATA %zu\n",zone->nchildren);
				fprintf(fp,"SCALARS Intensity_1-0 float 1 \n");
				fprintf(fp,"LOOKUP_TABLE default\n");
				//for(size_t izone = 0; izone < glb.nzone; izone++){
				for(size_t izone=0; izone < zone->nchildren; izone++ ){
					//zp = glb.zones[izone];
					Zone *zp = zone->children[izone];
					SpPhys *pp = zp->data;
					if(pp->non_empty_leaf){fprintf(fp,"%11.4e\n",pp->J_bar[0]);}
				 	else{fprintf(fp,"%11.4e\n",glb.I_cmb[0]*glb.I_norm[0]);}
				}
				fprintf(fp,"SCALARS Intensity_4-3 float 1 \n");
				fprintf(fp,"LOOKUP_TABLE default\n");
				for(size_t izone=0; izone < zone->nchildren; izone++ ){
					Zone *zp = zone->children[izone];
					SpPhys *pp = zp->data;
					if(pp->non_empty_leaf){fprintf(fp,"%11.4e\n",pp->J_bar[3]);}
				 	else{fprintf(fp,"%11.4e\n",glb.I_cmb[3]*glb.I_norm[3]);}
				}
				fprintf(fp,"SCALARS Population float 9 \n");
				fprintf(fp,"LOOKUP_TABLE default\n");
				for(size_t izone=0; izone < zone->nchildren; izone++ ){
					Zone *zp = zone->children[izone];
					SpPhys *pp = zp->data;
					if(pp->non_empty_leaf){
						fprintf(fp,"%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",\
						pp->pops[0][0],pp->pops[0][1],pp->pops[0][2],pp->pops[0][3],pp->pops[0][4],pp->pops[0][5],\
						pp->pops[0][6],pp->pops[0][7],pp->pops[0][8]);
					}
				 	else{fprintf(fp,"0 0 0 0 0 0 0 0 0\n");}
				}
				
				fprintf(fp,"SCALARS noise float 1\n");
				fprintf(fp,"LOOKUP_TABLE default\n");
				for(size_t izone = 0; izone < glb.nzone; izone++){
					Zone *zp = glb.zones[izone];
					SpPhys *pp = zp->data;
				 	fprintf(fp,"%11.4e\n",pp->diff);
				}
				fprintf(fp,"SCALARS nray float 1\n");
				fprintf(fp,"LOOKUP_TABLE default\n");
				for(size_t izone = 0; izone < glb.nzone; izone++){
					Zone *zp = glb.zones[izone];
					SpPhys *pp = zp->data;
				 	fprintf(fp,"%8zu\n",pp->nray);
				}
				fprintf(fp,"SCALARS converged float 1\n");
				fprintf(fp,"LOOKUP_TABLE default\n");
				for(size_t izone = 0; izone < glb.nzone; izone++){
					Zone *zp = glb.zones[izone];
					SpPhys *pp = zp->data;
				 	fprintf(fp,"%1zu\n",pp->converged );
				}
				
				fclose(fp);
			}
#endif
#if 0
			if(glb.stage==1 && Sp_MPIRANK == 0){
                                char filename[13];
				sprintf(filename,"nray%4.4zu.dat",iter);
				FILE *fp = fopen(filename,"w");
				for(size_t izone = 0; izone < glb.nzone; izone++){					
					Zone *zp = glb.zones[izone];
					SpPhys *pp = zp->data;
					fprintf(fp,"%12.5e %12.5e %12.5e %zu\n",zp->voxel.cen.x[0], \
									      zp->voxel.cen.x[1],\
					                                 zp->voxel.cen.x[2],pp->nray);
				}
				fclose(fp);
			}
#endif
			#ifdef HAVE_MPI
			if(Sp_MPISIZE > 1)
				MPI_Barrier(MPI_COMM_WORLD);
			#endif
		}
		
		if(glb.qmc && glb.fully_random)
                        for (size_t i =0; i < Sp_NTHREAD; i++)
                                gsl_qrng_free(glb.qrng[i]);
	}
#if 0	
	if(Sp_MPIRANK == 0){  
                char filename[13];
	  	sprintf(filename,"SPARX_final.vtk");
		FILE *fp=fopen(filename,"w");
		fprintf(fp,"# vtk DataFile Version 3.0\n");
		fprintf(fp,"MODEL\n");
		fprintf(fp,"ASCII\n");
		fprintf(fp,"DATASET STRUCTURED_GRID\n");
		fprintf(fp,"DIMENSIONS %zu %zu %zu \n",glb.model.grid->naxes.x[0],glb.model.grid->naxes.x[1],glb.model.grid->naxes.x[2]);
		//fprintf(fp,"POINTS %zu float\n",glb.nzone);
		fprintf(fp,"POINTS %zu float\n",zone->nchildren);
		//for(size_t izone = 0; izone < glb.nzone; izone++){
		for(size_t izone=0; izone < zone->nchildren; izone++ ){
			//Zone *zp = glb.zones[izone];
			Zone *zp = zone->children[izone];
			fprintf(fp,"%11.4e %11.4e %11.4e\n",zp->voxel.cen.x[0],zp->voxel.cen.x[1],zp->voxel.cen.x[2]);
		}
		//fprintf(fp,"POINT_DATA %zu\n",glb.nzone);
		fprintf(fp,"POINT_DATA %zu\n",zone->nchildren);
		fprintf(fp,"SCALARS DENSITY float 1\n");
		fprintf(fp,"LOOKUP_TABLE default\n");
		//for(size_t izone = 0; izone < glb.nzone; izone++){
		for(size_t izone=0; izone < zone->nchildren; izone++ ){
			//Zone *zp = glb.zones[izone];
			Zone *zp = zone->children[izone];
			SpPhys *pp = zp->data;
		 	fprintf(fp,"%11.4e\n",pp->n_H2);
		}
		fprintf(fp,"VECTORS VELOCITY float\n");
		//for(size_t izone = 0; izone < glb.nzone; izone++){
		for(size_t izone=0; izone < zone->nchildren; izone++ ){
			//Zone *zp = glb.zones[izone];
			Zone *zp = zone->children[izone];
			SpPhys *pp = zp->data;
		 	fprintf(fp,"%11.4e %11.4e %11.4e\n",pp->v_cen.x[0],pp->v_cen.x[1],pp->v_cen.x[2]);
		}
		fclose(fp);
	}
#endif
        if(glb.dat){
            if(Sp_MPIRANK == 0){
                char filename[32];
                sprintf(filename,"%s.dat",glb.outf->name);
                FILE *fp=fopen(filename,"w");
                for(size_t izone = 0; izone < glb.nzone; izone++){
                        Zone *zp = glb.zones[izone];
                        SpPhys *pp = zp->data;
                        double tempR;
                        switch (zp->voxel.geom){
                                case GEOM_REC3D:{
                                        double x = zp->voxel.cen.x[0]-zone->voxel.cen.x[0];
                                        double y = zp->voxel.cen.x[1]-zone->voxel.cen.x[1];
                                        double z = zp->voxel.cen.x[2]-zone->voxel.cen.x[2];
                                        tempR = sqrt( x*x + y*y+ z*z );
                                }
                                        break;
                                case GEOM_SPH1D:
                                        tempR = zp->voxel.cen.x[0];
                                        break;
                                default:
                                        Deb_ASSERT(0);
                        
                        }
                        fprintf(fp,"%11.4e ", tempR);
                        for (size_t j = 0; j < pp->mol->nlev; j++ )
                                fprintf(fp,"%11.4e ", pp->pops[0][j]);
                        fprintf(fp,"\n");
                }
                fclose(fp);
            }
        }
	/* Check for Python exceptions */
	if(!sts && PyErr_Occurred()) {
		sts = 1;
	}

	return sts;
}

/*----------------------------------------------------------------------------*/

static void *CalcExcThread(void *tid_p)
{
	size_t tid = *((size_t *)tid_p);
	double *hist = Mem_CALLOC(NHIST * NLEV, hist);
	double *popsold = Mem_CALLOC(NLEV, popsold);
        
#define TIMER 1
#if  TIMER        
        // Timer
        float Tmc_thread = 0.0;
        float Tdb_thread = 0.0;
        float Tall_thread = 0.0;
#endif
	for(size_t izone = 0; izone < glb.nzone; izone++) {
		/* Skip zones that don't belong to this rank/thread */
		if((glb.zone_tid[izone] != tid) || (glb.zone_rank[izone] != Sp_MPIRANK))
			continue;

		/* Check for thread termination */
		if(SpUtil_TermThread()) break;

		/* Zone related pointers */
		Zone *zp = glb.zones[izone];
		SpPhys *pp = zp->data;
		
		// restore
		for(size_t i = 0; i < NLEV; i++) {
                        popsold[i]=pp->pops[tid][i];
		}

		/* Buffers for calculating excitation */
		double *ds0 = Mem_CALLOC(pp->nray, ds0);
		double *vfac0 = Mem_CALLOC(pp->nray, vfac0);
		double *intensity = Mem_CALLOC(pp->nray * NRAD, intensity);
		double *tau = Mem_CALLOC(pp->nray * NRAD, tau);
                
                if(!glb.fully_random){
                        if(glb.qmc){
                                // reset QMC random seed
                                glb.qrng[tid] = gsl_qrng_alloc(gsl_qrng_niederreiter_2, QRAN_DIM);
                                // skip the initial zero-array 
                                double QRN[QRAN_DIM];
                                gsl_qrng_get(glb.qrng[tid], QRN);
                        }
                        else{
                        // reset random seed while in the fixed-ray stage
                                gsl_rng_set(glb.rng[tid], glb.seed);
                        }
                }
                
		/* Calculate NHIST times for statistics */
		for(size_t ihist = 0; ihist < (glb.stage == STAGE_RAN? NHIST:1); ihist++) {
#if  TIMER
                        clock_t start = clock();
#endif
                        #if 1
                        if(glb.qmc)
                                CalcRays_QRNG(tid, zp, ds0, vfac0, intensity, tau);
                        else
                                CalcRays_RNG(tid, zp, ds0, vfac0, intensity, tau);
                        #else 
			#define XI()\
				gsl_rng_uniform(glb.rng[tid])
			for(size_t i = 0; i < pp->nray; i++) {
				ds0[i] = XI() * 0.1 * Sp_LENFAC;
				double vel = (XI() - 0.5) * 4.3 * pp->width;
				vfac0[i] = Num_GaussNormal(vel, pp->width);
			}
			#undef XI
                        #endif
#if  TIMER
                        float Tmc = (float)(clock() - start) / (float)CLOCKS_PER_SEC;
#endif
			CalcDetailedBalance(tid, pp, ds0, vfac0, intensity, tau);
#if  TIMER
                        float Tall = (float)(clock() - start) / (float)CLOCKS_PER_SEC;
                        float Tdb =  Tall - Tmc;
                        
                        Tmc_thread += Tmc;
                        Tdb_thread += Tdb;
                        Tall_thread += Tall;
#endif
			if(glb.stage == STAGE_RAN){
				for(size_t i = 0; i < NLEV; i++) {
					HIST(ihist, i) = pp->pops[tid][i];
				}
			}
		}

                if(glb.qmc && !glb.fully_random)
                        gsl_qrng_free(glb.qrng[tid]);
		
		/* ================ */
		/* modified by I-Ta */
                double diff;
		if (glb.stage == STAGE_RAN){  
			// Stage 2, calculate difference by inner loop: Monte Carlo noise
			diff = CalcDiff(hist, NULL,izone,tid);
			pp->diff=diff;
		}
		else{
			// Stage 1, calculate difference by outer loop: iterative noise
			diff = 0.0;
			for(size_t i = 0; i < NLEV; i++) {
			     	if(pp->pops[tid][i]>glb.minpop){
					double temp = fabs( pp->pops[tid][i] - popsold[i] ) / pp->pops[tid][i];
					diff = ( diff>temp? diff:temp );
				}
			} 
		}
		/* ================ */
		
		
		
		
		/* Lock mutex for global parameters */
		pthread_mutex_lock(&glb.exc_mutex);

		/* Update max_diff */
		if( diff > glb.max_diff )
			glb.max_diff = diff;
		
		/* Determine convergence and act accordingly */
		/* The nray-doubling critiria modified by I-Ta */
		if( glb.fully_random ){
			if( diff < MCNOISE )
				glb.nconv += 1;
			if( diff > 0.5 * MCNOISE )
				pp->nray *= 4;
		}
		else{
			if( diff < TOLERANCE )
				glb.nconv += 1;
		}
		
		
		/* The old critiria
		if(diff <= (glb.fully_random ? MCNOISE : TOLERANCE)) {
			glb.nconv += 1;
		}
		else if(!glb.fully_random) {
			//debug
			//pp->nray += (size_t)(0.01 * glb.nray);
		}
		else {
			pp->nray *= 4;
		}
		*/
		

		/* Make sure nray does not exceed upper limit */
		if(pp->nray > MAXRAYS) {
			pp->nray = MAXRAYS;
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
#if TIMER	
	printf("Tid = %zu MC : %f %% , Detailed Balance : %f %% \n", tid, 100*Tmc_thread/Tall_thread, 100*Tdb_thread/Tall_thread);
#endif
#undef TIMER	

	free(hist);
	free(popsold);
	pthread_exit(NULL);
}

/*----------------------------------------------------------------------------*/

static void SyncProcs(void)
{
	/* Sync pops from all processes and threads */
	for(size_t i = 0; i < glb.nzone; i++)
		SyncPops(i);

	#ifdef HAVE_MPI
	if(Sp_MPISIZE > 1) {
		/* Sync total number of rays */
		int glb_nray_tot = (int)glb.nray_tot;
                int nray_tot_buff;
		MPI_Reduce(&glb_nray_tot, &nray_tot_buff, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		glb_nray_tot = (int)nray_tot_buff;
		MPI_Bcast(&glb_nray_tot, 1, MPI_INT, 0, MPI_COMM_WORLD);
		glb.nray_tot = (size_t)glb_nray_tot;

		/* Sync total number of converged zones */
		int glb_nconv = (int)glb.nconv;
                int nconv_buff;
		MPI_Reduce(&glb_nconv, &nconv_buff, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		glb_nconv = nconv_buff;
		MPI_Bcast(&glb_nconv, 1, MPI_INT, 0, MPI_COMM_WORLD);
		glb.nconv = (size_t)glb_nconv;

		/* Gather max_diff from all processes */
                double max_diff_buff[Sp_MPISIZE];
		MPI_Gather(&glb.max_diff, 1, MPI_DOUBLE, max_diff_buff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		if(Sp_MPIRANK == 0) {
			/* Sort max_diff_buff to get largest max_diff */
			Num_Qsort_d(max_diff_buff, Sp_MPISIZE);
			glb.max_diff = max_diff_buff[Sp_MPISIZE - 1];
		}

		/* Sync max_diff with master */
		MPI_Bcast(&glb.max_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	#endif

	return;
}

/*----------------------------------------------------------------------------*/

static void SyncPops(size_t izone)
{
	size_t zone_tid = glb.zone_tid[izone];
	Zone *zp = glb.zones[izone];
	SpPhys *pp = zp->data;

	/* First copy pops from all other threads to thread 0 */
	if(zone_tid != 0) {
		Mem_MEMCPY(pp->pops[0], pp->pops[zone_tid], NLEV);
	}

	#ifdef HAVE_MPI
	/* Then sync results from all processes if using MPI */
	size_t zone_rank = glb.zone_rank[izone];
	MPI_Status mpi_status;

	if(Sp_MPISIZE > 1) {
		/* Gather all calculated results to master process */
		if(Sp_MPIRANK != 0) {
			/* This is a slave process; send to master */
			if(zone_rank == Sp_MPIRANK) {
				MPI_Send(pp->pops[0], (int)NLEV, MPI_DOUBLE, 0, Sp_MPITAG, MPI_COMM_WORLD);
				MPI_Send(pp->tau, (int)NRAD, MPI_DOUBLE, 0, Sp_MPITAG, MPI_COMM_WORLD);
				MPI_Send(&pp->ds, 1, MPI_DOUBLE, 0, Sp_MPITAG, MPI_COMM_WORLD);
				
				//MPI_Send(pp->J_bar, (int)NRAD, MPI_DOUBLE, 0, Sp_MPITAG, MPI_COMM_WORLD);
				MPI_Send(&pp->nray, sizeof(size_t), MPI_CHAR, 0, Sp_MPITAG, MPI_COMM_WORLD);
				//MPI_Send(&pp->diff, 1, MPI_DOUBLE, 0, Sp_MPITAG, MPI_COMM_WORLD);
			}
		}
		else {
			/* This is the master process; receive from a slave */
			if(zone_rank != 0) {
				MPI_Recv(pp->pops[0], (int)NLEV, MPI_DOUBLE, (int)zone_rank, Sp_MPITAG, MPI_COMM_WORLD, &mpi_status);
				MPI_Recv(pp->tau, (int)NRAD, MPI_DOUBLE, (int)zone_rank, Sp_MPITAG, MPI_COMM_WORLD, &mpi_status);
				MPI_Recv(&pp->ds, 1, MPI_DOUBLE, (int)zone_rank, Sp_MPITAG, MPI_COMM_WORLD, &mpi_status);
				
				//MPI_Recv(pp->J_bar, (int)NRAD, MPI_DOUBLE, (int)zone_rank, Sp_MPITAG, MPI_COMM_WORLD, &mpi_status);
				MPI_Recv(&pp->nray, sizeof(size_t), MPI_CHAR, (int)zone_rank, Sp_MPITAG, MPI_COMM_WORLD, &mpi_status);
				//MPI_Recv(&pp->diff, 1, MPI_DOUBLE, (int)zone_rank, Sp_MPITAG, MPI_COMM_WORLD, &mpi_status);
			}
		}

		/* Sync all processes with master */
		MPI_Bcast(pp->pops[0], (int)NLEV, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(pp->tau, (int)NRAD, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		//MPI_Bcast(pp->J_bar, (int)NRAD, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&pp->nray, sizeof(size_t), MPI_CHAR, 0, MPI_COMM_WORLD);
		//MPI_Bcast(&pp->diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	#endif

	/* Copy thread 0 pops to all other threads */
	for(size_t i = 1; i < Sp_NTHREAD; i++) {
		Mem_MEMCPY(pp->pops[i], pp->pops[0], NLEV);
	}

	return;
}

/*----------------------------------------------------------------------------*/

static void CalcRays_RNG(size_t tid, Zone *zone, double *ds0, double *vfac0,
	double *intensity, double *tau)
/* Collect `external' contribution to local mean radiation field (J_bar)
 * by shooting NRAY rays in random directions and calculating the
 * corresponding intensity.
 */
{
	SpPhys *pp = zone->data;
	Deb_ASSERT(pp->nray > 0); /* Just in case */

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

	for(size_t i = 0; i < pp->nray; i++) {
		/* Set random ray origin and direction */
                GeRay ray = GeRay_Rand(glb.rng[tid], &zone->voxel);

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
		GeVec3_d v_gas = SpPhys_GetVgas(&ray.e, zone);
                double vel = ( PRAND() - 0.5) * 4.3 * pp->width 
                        + GeVec3_DotProd(&v_gas, &ray.d);
                #undef RAND
                #undef PRAND
                        
		/* Calculate radiative transfer along this direction */
		RadiativeXfer(tid, zone, &ray, vel, &ds0[i], &vfac0[i], &INTENSITY(i, 0), &TAU(i, 0));
	}

	#if debug_ray
	cpgclos();
	Deb_PAUSE();
	#endif

	/* Calculate average path length */
	pp->ds /= (double)pp->nray;

	return;
}

/*----------------------------------------------------------------------------*/

static void CalcRays_QRNG(size_t tid, Zone *zone, double *ds0, double *vfac0,
        double *intensity, double *tau)
/* Collect `external' contribution to local mean radiation field (J_bar)
 * by shooting NRAY rays in random directions and calculating the
 * corresponding intensity.
 */
{
        SpPhys *pp = zone->data;
        Deb_ASSERT(pp->nray > 0); /* Just in case */

        /* Reset pp->ds */
        pp->ds = 0;

        /* Reset tau */
        Mem_BZERO2(tau, pp->nray * NRAD);

        for(size_t i = 0; i < pp->nray; i++) {
                /* Set random ray origin and direction */
                double QRanNumber[QRAN_DIM];
                gsl_qrng_get(glb.qrng[tid], QRanNumber);
                GeRay ray = GeRay_QRand(QRanNumber, &zone->voxel);

                GeVec3_d v_gas = SpPhys_GetVgas(&ray.e, zone);
                double vel = ( QRanNumber[5] - 0.5) * 4.3 * pp->width 
                        + GeVec3_DotProd(&v_gas, &ray.d);

                /* Calculate radiative transfer along this direction */
                RadiativeXfer(tid, zone, &ray, vel, &ds0[i], &vfac0[i], &INTENSITY(i, 0), &TAU(i, 0));
        }

        /* Calculate average path length */
        pp->ds /= (double)pp->nray;

        return;
}

/*----------------------------------------------------------------------------*/

static void RadiativeXfer(size_t tid, Zone *zone, GeRay *ray, double vel, double *ds0,
	double *vfac0, double *intensity, double *tau)
/* Given a previously initialized ray, calculate intensity for all lines
 * along the ray from the ray origin to the edge of the cloud */
{
	size_t firststep = 1;
	Zone *zp = zone;
	size_t plane=6; // the ray is shooting outward instead of coming from somewhere initially

	#if debug_ray //debug
	size_t iter = 0;
	printf("%12s %12s %12s %12s %12s %12s %12s %12s\n", "iter", "ds", "vfac", "n_H2", "X_mol", "width", "dtau", "tau_nu");
	#endif

	/* Allocate array for dtau */
	double *dtau = Mem_CALLOC(NRAD, dtau);

	/* Reset intensity, tau, ds0 and vfac0 */
	Mem_BZERO2(intensity, NRAD);
	Mem_BZERO2(tau, NRAD);
	*ds0 = 0;
	*vfac0 = 0;

	/* Reset counter for average path length */
	SpPhys *zone_pp = zone->data;

	/* Propagate the ray through the cloud until we've
	 * reached the edge */
	while(zp) {
		/* Check for thread termination */
		if(SpUtil_TermThread()) break;

		/* Reset t */
		double t = 0;

		/* Reset dtau */
		Mem_BZERO2(dtau, NRAD);

		/* Calculate path to next boundary */
		GeRay_TraverseVoxel(ray, &zp->voxel, &t, &plane);

		/* Pointer to physical parameters associated with this zone */
		SpPhys *pp = zp->data;

		/* Do calculations on non-empty leaf zones only */
		if(pp->non_empty_leaf) {
			/* Calculate propagation path length and velocity factor
			 * for line profile */
			double ds = t * Sp_LENFAC;

			/* Increment average path length counter */
			zone_pp->ds += t;

			/* Calculate velocity line profile factor */
			double vfac = pp->has_tracer ? SpPhys_GetVfac(ray, t, vel, zp, 0) : 0.0;

			/* Calculate radiative contribution from neighboring
			 * zones */
			for(size_t i = 0; i < NRAD; i++) {
				double j_nu = 0.;
                                double k_nu = 0.;
                                if(pp->has_tracer) {
				 if(glb.overlap){
				  for(size_t j = 0; j < NRAD; j++) {
				   if(OVERLAP(i,j)==1){
                                        double tempj_nu, tempk_nu;
					if(i==j){
                                                /* Calculate molecular line emission and absorption coefficients */
                                                SpPhys_GetMoljk(tid, pp, j, vfac, &tempj_nu, &tempk_nu);
                                        }
					else{
						/* Calculate velocity line profile factor */
						double vfac2 = pp->has_tracer ? SpPhys_GetVfac(ray, t, vel-RELVEL(i,j), zp, 0) : 0.0;
						/* Calculate molecular line emission and absorption coefficients */
						SpPhys_GetMoljk(tid, pp, j, vfac2, &tempj_nu, &tempk_nu);
					}
					j_nu += tempj_nu;
					k_nu += tempk_nu;
				   }
				  }
				 }
                                 else{
                                        /* Calculate molecular line emission and absorption coefficients */
					SpPhys_GetMoljk(tid, pp, i, vfac, &j_nu, &k_nu);
				 }
				}
				else {
					j_nu = k_nu = 0.;
				}

				/* Add continuum emission/absorption */
				j_nu += pp->cont[i].j;
				k_nu += pp->cont[i].k;

				/* Calculate source function and optical depth if
				 * absorption is NOT zero */
				double S_nu = fabs(k_nu) > 0.0 ? j_nu / k_nu / glb.I_norm[i] : 0.0;

				dtau[i] = k_nu * ds;
				SpPhys_LIMITTAU(dtau[i]);

				if(!firststep) {
					/* Calculate intensity contributed by this step */
					intensity[i] += S_nu * (1.0 - exp(-dtau[i])) * exp(-tau[i]);

					/* Accumulate total optical depth for this line (must be done
					 * AFTER calculation of intensity!) */
					tau[i] += dtau[i];
					SpPhys_LIMITTAU(tau[i]);
				}
			}

			/* If this is the first step, save ds, vfac for later use and do
			   nothing else */
			if(firststep) {
				*ds0 = ds;
				*vfac0 = vfac;
				firststep = 0;
			}
			#if debug_ray
			printf("%12lu %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n", (unsigned long)iter, t, vfac, pp->n_H2, pp->X_mol, pp->width, dtau_nu, tau[0]);
			#endif
		}

		#if debug_ray //debug
		GeRay tmp_ray = GeRay_Inc(ray, t);
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
		zp = Zone_GetNext(zp, &plane, ray);

		#if debug_ray //debug
		iter++;
		#endif
	}

	/* Ray escaped cloud, add CMB to all lines */
	for(size_t i = 0; i < NRAD; i++) {
		intensity[i] += glb.I_cmb[i] * exp(-tau[i]);
	}

	#if debug_ray //debug
	printf("ds0=%12.4e, vfac0=%12.4e\n", (*ds0) / Sp_LENFAC, *vfac0);
	//Deb_PAUSE();
	#endif

	free(dtau);

	return;
}

/*----------------------------------------------------------------------------*/

static void CalcJbar(size_t tid, SpPhys *pp, const double *ds0, const double *vfac0,
	const double *intensity, const double *tau, double *J_bar)
/* Calculate J_bar, the mean radiation field intensity, by averaging the
 * radiation field over direction and velocity */
{
	/* Reset Jbar (very important!) */
	//Mem_BZERO2(pp->J_bar, NRAD);
	Mem_BZERO2(J_bar, NRAD);

	/* Reset tau */
	Mem_BZERO2(pp->tau, NRAD);

	/* Reset vfac0_sum */
	double vfac0_sum = 0.;

	/* Loop through all rays and sum intensities for averaging */
	for(size_t i = 0; i < pp->nray; i++) {
		/* Loop through lines */
		for(size_t j = 0; j < NRAD; j++) {
			/* Calculate local emission and absorption */
                        double j_nu, k_nu;
			SpPhys_GetMoljk(tid, pp, j, vfac0[i], &j_nu, &k_nu);

			/* Add continuum emission/absorption */
			j_nu += pp->cont[j].j;
			k_nu += pp->cont[j].k;

			if(fabs(k_nu) > 0) {
				double dtau_nu = k_nu * ds0[i];
				double S_nu = j_nu / k_nu / glb.I_norm[j];

				SpPhys_LIMITTAU(dtau_nu);

				/* Accumulate intensities weighted over vfac in J_bar: the average
				 * over velocity is extremely important -- failure to do so would
				 * result in a dependence of the excitation on size of grid! */
				//pp->J_bar[j] += vfac0[i] * (INTENSITY(i, j) * exp(-dtau_nu) + S_nu * (1.0 - exp(-dtau_nu)));
				J_bar[j] += vfac0[i] * (INTENSITY(i, j) * exp(-dtau_nu) + S_nu * (1.0 - exp(-dtau_nu)));
				
				/* Store total tau in zone for bookkeeping */
				pp->tau[j] += vfac0[i] * (TAU(i, j) + dtau_nu);
			}
		}

		/* Accumulate vfac0_sum */
		vfac0_sum += vfac0[i];
	}

	Deb_ASSERT(vfac0_sum >= 0); /* Just in case */

	if(vfac0_sum > 0) {
		for(size_t i = 0; i < NRAD; i++) {
			/* Denormalized and average J_bar */
			//pp->J_bar[i] = pp->J_bar[i] * glb.I_norm[i] / vfac0_sum;
			J_bar[i] = J_bar[i] * glb.I_norm[i] / vfac0_sum;

			/* Calculate averaged tau for this zone */
			pp->tau[i] /= vfac0_sum;
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
	size_t  max_diff_lev = 0;
	double diff = 0;
	const double MAXDIFF = glb.minpop;//TOLERANCE * 0.1;

	/* Allocate J_bar array (no need now) */
	double *J_bar = Mem_CALLOC(NRAD, J_bar);

	/* Allocate rates matrix: (NLEV + 1) x NLEV,
	 * where the additional row is for the constraint
	 * that all levels must sum to unity */
	double *rmat = Mem_CALLOC((NLEV + 1) * NLEV, rmat);

	#define RMAT(i, j)\
		rmat[(j) + NLEV * (i)]
	#define CMAT(i, j)\
		(pp->cmat[(j) + NLEV * (i)])

	/* RHS of rate equation */
	double *rhs = Mem_CALLOC(NLEV + 1, rhs);
	rhs[NLEV] = 1.0;

	/* Allocate hist array */
	double *hist = Mem_CALLOC(NHIST * NLEV, hist);

	for(size_t iter = 0; iter < glb.maxi; iter++) {
		for(size_t ihist = 0; ihist < NHIST; ihist++) {
			/* Calculate J_bar, the mean radiation field intensity */
			CalcJbar(tid, pp, ds0, vfac0, intensity, tau, J_bar);

			/* Reset rates matrix */
			Mem_BZERO2(rmat, (NLEV + 1) * NLEV);

			/* Add radiative terms to rates matrix */
			for(size_t i = 0; i < NRAD; i++) {
				size_t up = RAD(i)->up;
				size_t lo = RAD(i)->lo;

				/* Diagonal terms are transitions `out of' row state */
				RMAT(up, up) -= (RAD(i)->A_ul + J_bar[i] * RAD(i)->B_ul);
				RMAT(lo, lo) -= (J_bar[i] * RAD(i)->B_lu);

				/* Off-diagonal terms are transitions `into' row state */
				RMAT(up, lo) += (J_bar[i] * RAD(i)->B_lu);
				RMAT(lo, up) += (RAD(i)->A_ul + J_bar[i] * RAD(i)->B_ul);
			}

			/* Add collisional terms to rates matrix */
			for(size_t i = 0; i < NLEV; i++) {
				for(size_t j = 0; j < NLEV; j++) {
					if(i == j) {
					/* Diagonal terms are minus the sums of the rates from all
					 * collisional transitions `out of' state i */
						for(size_t k = 0; k < NLEV; k++)
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
			for(size_t i = 0; i < NLEV; i++) {
				if(pp->pops[tid][i] < 0)
					pp->pops[tid][i] = 0.0;
			}

			for(size_t i = 0; i < NLEV; i++) {
				HIST(ihist, i) = pp->pops[tid][i];
			}
		}
		/* Calculate relative difference */
		//diff = CalcDiff(hist, &max_diff_lev);
		diff = CalcDiff2(hist, &max_diff_lev);

		/* Stop iterating if diff is less than MAXDIFF */
		if(diff <= MAXDIFF)
			break;
	}

	/* Warn for non-convergence */
	if(diff > MAXDIFF) {
		Sp_PWARN("non-convergence detected at level %lu in zone <%lu,%lu,%lu> (diff=%.3e MAXDIFF=%.3e %.3e %.3e )\n",
			max_diff_lev,
			GeVec3_X(pp->zp->index, 0),
			GeVec3_X(pp->zp->index, 1),
			GeVec3_X(pp->zp->index, 2),
			diff,MAXDIFF,TOLERANCE,glb.tolerance);
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

static double CalcDiff(const double *hist, size_t *max_diff_lev,size_t izone,size_t tid)
{
	double max_diff = 0.0;
	double *diffs = Mem_CALLOC(NHIST, diffs);
	
	Zone *zp = glb.zones[izone];
	SpPhys *pp = zp->data;
	
	/* Loop through all levels */
	for(size_t i = 0; i < NLEV; i++) {
		/* Calculate mean for this level */
		double mean = 0.;
		for(size_t j = 0; j < NHIST; j++) {
			mean += HIST(j, i);
		}
		mean /= (double)NHIST;
		
		/* ============= */
		/* added by I-Ta */
		pp->pops[tid][i]=mean;
		/* ============= */

		/* Calculate relative difference from mean if
		 * mean >= minpop */
		if(mean >= glb.minpop) {
			for(size_t j = 0; j < NHIST; j++)
				diffs[j] = fabs(HIST(j, i) - mean) / mean;

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
	#define NDIFF (NHIST - 1)
	double max_diff = 0.0;
	double *diffs = Mem_CALLOC(NDIFF, diffs);
	double *pops = Mem_CALLOC(NHIST, pops);

	/* Loop through all levels */
	for(size_t i = 0; i < NLEV; i++) {
		/* Load pops into array for sorting */
		for(size_t j = 0; j < NHIST; j++) {
			pops[j] = HIST(j, i);
		}
		/* Get smallest pops */
		Num_Qsort_d(pops, NHIST);

		if(pops[0] >= glb.minpop) {
			/* Calculate difference of hist[0] between all
			   other hist values */
			for(size_t j = 0; j < NDIFF; j++) {
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





