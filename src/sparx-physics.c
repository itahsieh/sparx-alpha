#include "sparx.h"
#include <gsl/gsl_interp.h>

/*----------------------------------------------------------------------------*/

void *SpPhys_Alloc(const Zone *zp, const void *parms_p)
{
	const SpPhysParm *parms = parms_p;
	SpPhys *pp = Mem_CALLOC(1, pp);

	/* Back-reference to zone */
	pp->zp = zp;

	/* Init according to parms */
	if(parms) {
		/* Assign molecule and allocate lines */
		if(parms->mol)
			SpPhys_InitMol(pp, parms->mol,0);

		/* Set velocity field */
		if(parms->velfield) {
			pp->velfield = parms->velfield;
		}
		else {
			pp->velfield = NULL;
		}
	}

	return pp;
}

/*----------------------------------------------------------------------------*/

void SpPhys_Free(void *ptr)
{
	SpPhys *pp = ptr;

	for(size_t i = 0; i < Sp_NTHREAD; i++) {
		if(pp->pops[i])
			free(pp->pops[i]);
	}

	if(pp->cmat)
		free(pp->cmat);

	if(pp->tau)
		free(pp->tau);

	if(pp->cont)
		free(pp->cont);

	free(pp);

	return;
}

/*----------------------------------------------------------------------------*/

void SpPhys_InitMol(SpPhys *pp, const Molec *mol,int popsold)
{
	Deb_ASSERT(mol != NULL);
	if(!popsold){
		Deb_ASSERT(pp->mol == NULL);printf("tau=%zu \n",pp->tau);
		Deb_ASSERT(pp->tau == NULL);
	}

	if(!popsold){
		for(size_t i = 0; i < Sp_NTHREAD; i++) {
			Deb_ASSERT(pp->pops[i] == NULL);
			pp->pops[i] = Mem_CALLOC(mol->nlev, pp->pops[i]);
		}
	}
	/* Allocate tau for book keeping */
	pp->tau = Mem_CALLOC(mol->nrad, pp->tau);
	//pp->J_bar = Mem_CALLOC(mol->nrad, pp->J_bar);

	/* Allocate continuum emission/absorption */
	double *freq = Mem_CALLOC(mol->nrad, freq);
	for(size_t i = 0; i < mol->nrad; i++) {
		freq[i] = mol->rad[i]->freq;
	}
	SpPhys_InitContWindows(pp, freq, mol->nrad);

	free(freq);

	return;
}

/*----------------------------------------------------------------------------*/

void SpPhys_InitContWindows(SpPhys *pp, double freq[], size_t nfreq)
{
	Deb_ASSERT(nfreq > 0);

	if(pp->cont)
		free(pp->cont);

	/* Allocate continuum emission/absorption */
	pp->ncont = nfreq;
	pp->cont = Mem_CALLOC(nfreq, pp->cont);

	for(size_t i = 0; i < nfreq; i++) {
		Deb_ASSERT(freq[i] > 0);
		pp->cont[i].freq = freq[i];
		pp->cont[i].lambda = PHYS_CONST_MKS_LIGHTC / freq[i];
	}

	return;
}

/*----------------------------------------------------------------------------*/

void SpPhys_InitCollRates(SpPhys *pp)
/* Allocate and initialize collisional rates:
 * Interpolate downward rates and infer upward rates from Boltzmann relation */
{
	Deb_ASSERT(pp->mol != NULL);
	Deb_ASSERT(pp->cmat == NULL);

	#define NLEV (pp->mol->nlev)
	#define COL(i) (pp->mol->col[(i)])
	#define TMP(i, j) (COL(i)->tmp[(j)])
	#define TR(i, j) (COL(i)->tr[(j)])
	#define CMAT(i, j) (pp->cmat[(j) + NLEV * (i)])

	/* Allocate collisional rates matrix */
	pp->cmat = Mem_CALLOC(NLEV * NLEV, pp->cmat);

	/* Fill in downward rates: collisional rates for each transition are the
	 * sum of collisional rates from ALL collisional partners */
	for(size_t i = 0; i < pp->mol->ncol; i++) {
		/* Locate nearest temperature available for this species */
		size_t itmp = gsl_interp_bsearch(COL(i)->tmp, pp->T_k, (size_t)0, COL(i)->ntmp);

		/* Loop through all collisional transitions and calculate
		 * collisional rats */
		for(size_t j = 0; j < COL(i)->ntr; j++) {
                        double K_ul;
			/* Interpolate downward rate coeffs */
			if(itmp == COL(i)->ntmp - 1) {
				/* T_k greater than available tempratures, coerce to
				 * upper end of K_ul */
				K_ul = TR(i, j)->K_ul[COL(i)->ntmp - 1];
			}
			else if(pp->T_k < TMP(i, 0)) {
				/* T_k less than available temperatures, coerce to
				 * lower end of K_ul */
				K_ul = TR(i, j)->K_ul[0];
			}
			else {
				/* T_k within range, linearly interpolate K_ul */
				K_ul = Num_InterpLinear(pp->T_k, TMP(i, itmp), TMP(i, itmp + 1), TR(i, j)->K_ul[itmp], TR(i, j)->K_ul[itmp + 1]);
			}

			/* Collisional rate is density of collisional partner multiplied
			 * by downard rate */
			CMAT(TR(i, j)->up, TR(i, j)->lo) += SpPhys_GetCollDens(pp, COL(i)->species) * K_ul;
		}
	}

	/* Calculate upward rates from downward rates using the Boltzmann
	 * relation (cf. Rohlfs & Wilson p. 309)
	 * 	C_lu / C_ul = N_u / N_l = (g_u / g_l) * exp(-(E_u - E_l) / (k * T)) = BoltzRatio
	 * -->	C_lu = C_ul * BoltzRatio
	 */
	for(size_t i = 0; i < NLEV; i++) {
		for(size_t j = i + 1; j < NLEV; j++) {
			CMAT(i, j) = CMAT(j, i) * SpPhys_BoltzRatio(pp->mol, j, i, pp->T_k);
		}
                
                /* Diagonal terms are minus the sums of the rates from all
                 * collisional transitions `out of' state i */
                for (size_t j = 0; j < NLEV; j++){
                        if( j != i )
                                CMAT(i,i) -= CMAT(i,j);
                }
               
	}

        /* matrix transpose 
         * Off-diagonal terms are sums of rates from state j
         * `into' state i */
        for(size_t i = 0; i < NLEV; i++) {
                for (size_t j = 0; j < NLEV; j++){
                        if( j != i ){
                                double temp = CMAT(i,j);
                                CMAT(i,j) = CMAT(j,i);
                                CMAT(j,i) = temp;
                        }
                }
        }

	#undef NLEV
	#undef COL
	#undef TMP
	#undef TR
	#undef CMAT

	return;
}

/*----------------------------------------------------------------------------*/

void SpPhys_SetContinuumIntens_bb(SpPhys *pp, int cont, double T_bb, double I_norm)
/* Add continuum abosrption and emission associated with T_bb, kap and rho to total
   continuum absorption and emission for all line or continuum windows. */
{
	size_t nrad;

	if(cont) {
		Deb_ASSERT(pp->cont != NULL);
		nrad = pp->ncont;
	}
	else {
		Deb_ASSERT(pp->mol != NULL);
		nrad = pp->mol->nrad;
	}

	for(size_t i = 0; i < nrad; i++) {
                #define FREQ(i)\
                        (cont ? pp->cont[i].freq : pp->mol->rad[i]->freq)
		pp->cont[i].I_bb = Phys_PlanckFunc(FREQ(i), T_bb) / I_norm;
	}

	#undef FREQ

	return;
}

/*----------------------------------------------------------------------------*/

void SpPhys_AddContinuum(SpPhys *pp, int cont, double T_bb, const Kappa *kap, double rho)
/* Add continuum abosrption and emission associated with T_bb, kap and rho to total
   continuum absorption and emission for all line or continuum windows. */
{
	size_t nrad;

	Deb_ASSERT(kap != NULL);

	if(cont) {
		Deb_ASSERT(pp->cont != NULL);
		nrad = pp->ncont;
	}
	else {
		Deb_ASSERT(pp->mol != NULL);
		nrad = pp->mol->nrad;
	}

	#define FREQ(i)\
		(cont ? pp->cont[i].freq : pp->mol->rad[i]->freq)

	for(size_t i = 0; i < nrad; i++) {
		/* k = kappa_dust * rho_dust */
		double k_nu = Kap_FromFreq(kap, FREQ(i)) * rho;

		/* j = B_nu * k */
		double j_nu = Phys_PlanckFunc(FREQ(i), T_bb) * k_nu;

		/* Accumulate j and k */
		pp->cont[i].j += j_nu;
		pp->cont[i].k += k_nu;
	}

	
	#undef FREQ

	return;
}

/*----------------------------------------------------------------------------*/

void SpPhys_AddContinuum_d(SpPhys *pp, int cont, double gas_to_dust)
{
	/* Dust mass density is
		rho_dust = (n_H2 * mu_H2 * amu) / gas_to_dust
	   where
	   	n_H2 = H2 number density
		gas_to_dust = gas-to-dust ratio
		mu_H2 = mean molecular weight per H2 --
		        assuming M_H : M_He : M_Z = 0.71 : 0.27 : 0.02,
		        this would be ~ 2.8
		amu = atomic mass unit
	*/
	Deb_ASSERT(gas_to_dust > 0);
	double rho = (pp->n_H2 * 2.8 * PHYS_UNIT_MKS_AMU) / gas_to_dust;

	Kappa *kap = SpIO_LoadKappa(pp->kapp_d);

        SpPhys_AddContinuum(pp, cont, pp->T_d, kap, rho);

        Kap_Free(kap);

	return;
}

/*----------------------------------------------------------------------------*/

void SpPhys_AddContinuum_ff(SpPhys *pp, int cont)
{
	/* Total mass density is
		rho_dust = (n_H2 * mu_H2 * amu)
	   where
	   	n_H2 = H2 number density
		mu_H2 = mean molecular weight per H2 --
		        assuming M_H : M_He : M_Z = 0.71 : 0.27 : 0.02,
		        this would be ~ 2.8
		amu = atomic mass unit
	*/
	double rho = (pp->n_H2 * 2.8 * PHYS_UNIT_MKS_AMU);
	Kappa *kap = SpIO_LoadKappa(pp->kapp_ff);
	SpPhys_AddContinuum(pp, cont, pp->T_ff, kap, rho);
	Kap_Free(kap);

	return;
}

/*----------------------------------------------------------------------------*/

#if 0
void SpPhys_AddCont(SpPhys *pp, int cont)
{
	if(cont)
		Deb_ASSERT(pp->cont != NULL);
	else
		Deb_ASSERT(pp->mol != NULL);

	Deb_ASSERT(pp->bbsrc.nu0 > 0);

	size_t i, nrad;
	double rho, j_nu, k_nu;

	if(cont)
		nrad = pp->ncont;
	else
		nrad = pp->mol->nrad;

	/* Total mass density is
		rho_dust = (n_H2 * mu_H2 * amu)
	   where
	   	n_H2 = H2 number density
		mu_H2 = mean molecular weight per H2 --
		        assuming M_H : M_He : M_Z = 0.71 : 0.27 : 0.02,
		        this would be ~ 2.8
		amu = atomic mass unit
	*/
	rho = (pp->n_H2 * 2.8 * PHYS_UNIT_MKS_AMU);

	#define FREQ(i)\
		(cont ? pp->cont[i].freq : pp->mol->rad[i]->freq)

	for(i = 0; i < nrad; i++) {
		/* k = kappa * rho */
		k_nu = pp->bbsrc.kappa0 * pow(FREQ(i) / pp->bbsrc.nu0, pp->bbsrc.beta) * rho;

		/* j = B_nu * k */
		j_nu = Phys_PlanckFunc(FREQ(i), pp->bbsrc.T_bb) * k_nu;

		/* Accumulate j and k */
		pp->cont[i].j += j_nu;
		pp->cont[i].k += k_nu;
	}

	#undef FREQ

	return;
}

/*----------------------------------------------------------------------------*/

void SpPhys_AddDust(SpPhys *pp, int cont, const Kappa *kap, double gas_to_dust)
{
	Deb_ASSERT(kap != NULL);
	if(cont)
		Deb_ASSERT(pp->cont != NULL);
	else
		Deb_ASSERT(pp->mol != NULL);

	size_t i, nrad;
	double rho_dust, j_nu, k_nu;

	if(cont)
		nrad = pp->ncont;
	else
		nrad = pp->mol->nrad;

	/* Dust mass density is
		rho_dust = (n_H2 * mu_H2 * amu) / gas_to_dust
	   where
	   	n_H2 = H2 number density
		gas_to_dust = gas-to-dust ratio
		mu_H2 = mean molecular weight per H2 --
		        assuming M_H : M_He : M_Z = 0.71 : 0.27 : 0.02,
		        this would be ~ 2.8
		amu = atomic mass unit
	*/
	Deb_ASSERT(gas_to_dust > 0);
	rho_dust = (pp->n_H2 * 2.8 * PHYS_UNIT_MKS_AMU) / gas_to_dust;

	#define FREQ(i)\
		(cont ? pp->cont[i].freq : pp->mol->rad[i]->freq)

	for(i = 0; i < nrad; i++) {
		/* k = kappa_dust * rho_dust */
		k_nu = Kap_FromFreq(kap, FREQ(i)) * rho_dust;

		/* j = B_nu * k */
		j_nu = Phys_PlanckFunc(FREQ(i), pp->T_k) * k_nu;

		/* Accumulate j and k */
		pp->cont[i].j += j_nu;
		pp->cont[i].k += k_nu;
	}

	#undef FREQ

	return;
}
#endif

/*----------------------------------------------------------------------------*/

#if 0
void SpPhys_AddKappa(SpPhys *pp)
/* Once SpPhys.cont has been initialized, add continuum emission/absorption
   specified by SpPhys.kappa to all frequencies in SpPhys.cont */
{
	size_t i, nrad;
	double rho_tot, j_nu, k_nu;

	Deb_ASSERT(pp->kappa != NULL);

	if(cont)
		nrad = pp->ncont;
	else
		nrad = pp->mol->nrad;

	/* Total mass density is
		rho_dust = (n_H2 * mu_H2 * amu)
	   where
	   	n_H2 = H2 number density
		mu_H2 = mean molecular weight per H2 --
		        assuming M_H : M_He : M_Z = 0.71 : 0.27 : 0.02,
		        this would be ~ 2.8
		amu = atomic mass unit
	*/
	rho_dust = (pp->n_H2 * 2.8 * PHYS_UNIT_MKS_AMU);

	#define FREQ(i)\
		(cont ? pp->cont[i].freq : pp->mol->rad[i]->freq)

	for(i = 0; i < nrad; i++) {
		/* k = kappa_dust * rho_dust */
		k_nu = Kap_FromFreq(ppkap, FREQ(i)) * rho_dust;

		/* j = B_nu * k */
		j_nu = Phys_PlanckFunc(FREQ(i), pp->T_k) * k_nu;

		/* Accumulate j and k */
		pp->cont[i].j += j_nu;
		pp->cont[i].k += k_nu;
	}

	#undef FREQ
}
#endif

/*----------------------------------------------------------------------------*/

double SpPhys_GetCollDens(const SpPhys *pp, int species)
{
	switch(species) {
		case 1: /* H2 */
			return pp->n_H2;
		case 2: /* para-H2 */
			return pp->n_H2 * pp->X_pH2;
		case 3: /* ortho-H2 */
			return pp->n_H2 * pp->X_oH2;
		case 4: /* electrons */
			return pp->n_H2 * pp->X_e;
		case 5: /* H */
			return pp->n_H2 * pp->X_H;
		case 6: /* He */
			return pp->n_H2 * pp->X_He;

		default: /* Illegal species code */
			Deb_ASSERT(0);
	}

	/* This should never happen */
	return 0.0;
}

/*----------------------------------------------------------------------------*/

void SpPhys_Fprintf(SpPhys *pp, FILE *fp)
{
	fprintf(fp, "n_H2=%g, T_k=%g, X_mol=%g, width=%g, dust=`%s'\n",
		pp->n_H2, pp->T_k, pp->X_mol, pp->width, strlen(pp->kapp_d) > 0 ? pp->kapp_d : "None");

	if(pp->mol) {
		Deb_ASSERT(pp->pops[0] != 0);
		fprintf(fp, " %5s|%20s\n", "Level", "Fractional density");
		fprintf(fp, " %5s|%20s\n", "-----", "--------------------");
		for(size_t i = 0; i < pp->mol->nlev; i++) {
			fprintf(fp, " %5lu|%20g\n", (unsigned long)i, pp->pops[0][i]);
		}
	}

	return;
}

/*----------------------------------------------------------------------------*/

size_t SpPhys_Fwrite(SpPhys *pp, FILE *fp)
{
	size_t nbytes = 0;

	nbytes += Mem_FWRITE(&pp->n_H2, 1, fp);
	nbytes += Mem_FWRITE(&pp->T_k, 1, fp);
	nbytes += Mem_FWRITE(&pp->X_mol, 1, fp);
	nbytes += Mem_FWRITE(&pp->width, 1, fp);

	if(pp->mol) {
		Deb_ASSERT(pp->pops[0] != 0);
		nbytes += Mem_FWRITE(pp->pops[0], pp->mol->nlev, fp);
	}

	return nbytes;
}

/*----------------------------------------------------------------------------*/

size_t SpPhys_Fread(SpPhys *pp, FILE *fp)
{
	size_t nbytes = 0;

	nbytes += Mem_FREAD(&pp->n_H2, 1, fp);
	nbytes += Mem_FREAD(&pp->T_k, 1, fp);
	nbytes += Mem_FREAD(&pp->X_mol, 1, fp);
	nbytes += Mem_FREAD(&pp->width, 1, fp);

	if(pp->mol) {
		Deb_ASSERT(pp->pops[0] != 0);
		nbytes += Mem_FREAD(pp->pops[0], pp->mol->nlev, fp);
	}

	return nbytes;
}

/*----------------------------------------------------------------------------*/

double SpPhys_Zfunc(const Molec *mol, double T_k)
/* Calculate Z(T), the partition function */
{
	size_t i;
	double Z = 0, k = PHYS_CONST_MKS_BOLTZK;

	for(i = 0; i < mol->nlev; i++) {
		Z += mol->lev[i]->g * exp(-mol->lev[i]->E / (k * T_k));
	}

	return Z;
}

/*----------------------------------------------------------------------------*/

void SpPhys_ProcLamda(Molec *mol)
/* Process a molecule loaded from a LAMDA comaptible data file:
 * Convert everything to SI units */
{
	static double
		m_u = PHYS_UNIT_MKS_AMU,
		cm = PHYS_UNIT_MKS_CM,
		cc = PHYS_UNIT_MKS_CC,
		h = PHYS_CONST_MKS_PLANCKH,
		c = PHYS_CONST_MKS_LIGHTC;
	size_t i, j, k, up, lo;
	double g_u, g_l, E_u, E_l, nu;
	MolTrRad *rad;

	/* Convert molecular weight to kg */
	mol->weight *= m_u;

	/* Convert energy of each level from cm^-1 to J,
	 * according to E = h * c / lambda
	 *
	 * where E = energy
	 *       h = Planck's constant
	 *       c = speed of light in vacuum
	 *       lambda = transition wavelength
	 */
	for(i = 0; i < mol->nlev; i++)
		mol->lev[i]->E = h * c * (mol->lev[i]->E * (1.0 / cm));

	/* Calculate line parameters */
	for(i = 0; i < mol->nrad; i++) {
		rad = mol->rad[i];
		up = rad->up;
		lo = rad->lo;
		g_u = mol->lev[up]->g;
		g_l = mol->lev[lo]->g;
		E_u = mol->lev[up]->E;
		E_l = mol->lev[lo]->E;

		/* Recalculate frequency based on level energies:
		 * 	nu = (E_u - E_l) / h
		 *
		 * where nu = frequency
		 *       E_u = uppler level energy
		 *       E_l = lower level energy
		 *       h = Planck's constant
		 */
		nu = rad->freq = (E_u - E_l) / h;

		/* Einstein B coefficient for stimulated emission:
		 * 	B_ul = A_ul * c^2 / (2.0 * h * nu^3)
		 *
		 * where A_ul = Einstein A coefficient
		 *       c = speed of light in vacuum
		 *       h = Planck's constant
		 *       nu = transition frequency
		 */
		rad->B_ul = rad->A_ul * c * c / (2.0 * h * pow(nu, 3.0));

		/* Einstein B coefficient for absorption:
		 *	B_lu * g_l = B_ul * g_u
		 *
		 * where g_u = upper level statistical weight
		 *       g_l = lower level statistical weight
		 */
		rad->B_lu = (g_u / g_l) * rad->B_ul;
	}

	/* Convert collisional rate coefficients from Gaussian units
	 * to SI units */
	for(i = 0; i < mol->ncol; i++) {
		for(j = 0; j < mol->col[i]->ntr; j++) {
			for(k = 0; k < mol->col[i]->ntmp; k++) {
				mol->col[i]->tr[j]->K_ul[k] *= cc; /* [cm^3 s^-1] -> [m^3 s^-1] */
			}
		}
	}

	return;
}

/*----------------------------------------------------------------------------*/

double SpPhys_BoltzPops(const Molec *mol, size_t lev, double T_k)
/* Given kinetic temperature, calculate thermal equilibrium fractional
 * density for a particular level (i.e. level population) */
{
	static double
		k = PHYS_CONST_MKS_BOLTZK; /* Boltzmann's constant */
	double
		g = mol->lev[lev]->g, /* Statistical weight */
		E = mol->lev[lev]->E, /* Level energy */
		Z = SpPhys_Zfunc(mol, T_k); /* Partition function */

	/* The fractional density of lev is
	 *          g * e^(-E/kT) / Z
	 */
	return g * exp(-E / (k * T_k)) / Z;
}

/*----------------------------------------------------------------------------*/

void SpPhys_GetMoljk(size_t tid, const SpPhys *pp, size_t tr, double vfac, double *j_nu, double *k_nu)
/* Calculate mlecular emission and absorption coefficients (cf. Rybicki &
 * Lightman 1985).
 *
 * Parameters:
 *   tr -- index of transition
 *   vfac -- exponential term of line profile
 */
{
	/* Constants */
	static const double
		pi = PHYS_CONST_PI,
		c = PHYS_CONST_MKS_LIGHTC,
		KONST = PHYS_CONST_MKS_PLANCKH / (4.0 * PHYS_CONST_PI);
	const MolTrRad *trans = pp->mol->rad[tr];
	double
		nu = trans->freq,
		n_u = pp->pops[tid][trans->up],
		n_l = pp->pops[tid][trans->lo],
		factor;
		/* Factor is angle-averaged photon energy multiplied by
		 * the line profile function:
		 * 	(h * nu / (4 * pi)) * phi
		 *
		 * where h = Planck's constant
		 *       nu = line center frequency
		 *       phi = line profile function
		 */
		factor = KONST * pp->n_H2 * pp->X_mol * nu * (c / (pp->width * nu * sqrt(pi))) * vfac;

	#if 0 //debug
	printf("factor=%10.4e, n_H2=%10.4e, X_mol=%10.4e\n", factor, pp->n_H2, pp->X_mol);
	#endif

	/* Emission coefficient */
	*j_nu = factor * (n_u * trans->A_ul);

	/* Absorption coefficient */
	*k_nu = factor * (n_l * trans->B_lu - n_u * trans->B_ul);

	#if 0 //debug
	printf("factor=%10.4e, j_nu=%10.4e, k_nu=%10.4e, vfac=%10.4e\n", factor, *j_nu, *k_nu, vfac);
	#endif
	return;
}

/*----------------------------------------------------------------------------*/

GeVec3_d SpPhys_GetVgas(const GeVec3_d *pos, const Zone *zone)
/* Retriev gas velocity at pos */
{
	#define USE_LVG 0
	#define USE_CONST 0
        
	GeVec3_d v_gas;

        #if USE_CONST
	double velo = -0.2e3; // [km/s]
	#endif

	/* Get zone physics */
	SpPhys *pp = zone->data;

	switch(zone->voxel.geom) {
		case GEOM_SPH1D:
			#if USE_LVG
			double radius = GeVec3_Mag(pos);
			v_gas = GeVec3_Normalize(pos);
			v_gas = GeVec3_Scale(&v_gas, 100.0e3 * radius); /* 100 km/s/pc velocity gradient */
			#elif USE_CONST
			v_gas = GeVec3_Normalize(pos);
			v_gas = GeVec3_Scale(&v_gas, velo); /* 100 km/s/pc velocity gradient */
			#else
			/* Project radial velocity onto position vector */
			v_gas = GeVec3_Normalize(pos);
			v_gas = GeVec3_Scale(&v_gas, GeVec3_X(pp->v_cen, 0));
			#endif
			break;
			
		case GEOM_SPH3D:{
			/* Project radial velocity onto position vector */
			/* get normalized vr, vt vp */
                        GeVec3_d vr = GeVec3_Normalize(pos);
                        
                        GeVec3_d vp;
                        GeVec3_d CylPos;
                        CylPos = GeVec3_Cart2Cyl(pos);
			if( CylPos.x[0] == 0.0 ){
				GeVec3_X(vp,0) = 0.0;
				GeVec3_X(vp,1) = 1.0;
				GeVec3_X(vp,2) = 0.0;
			}
			else{
				GeVec3_X(vp,0) = -pos->x[1];
				GeVec3_X(vp,1) = pos->x[0];
				GeVec3_X(vp,2) = 0.0;	
				vp = GeVec3_Normalize(&vp);
			}

			GeVec3_d vt = GeVec3_CrossProd(&vp,&vr);
                        
			/* Project directional velocity onto the normalized vector */
			vr = GeVec3_Scale(&vr, GeVec3_X(pp->v_cen, 0));
			vt = GeVec3_Scale(&vt, GeVec3_X(pp->v_cen, 1));
			vp = GeVec3_Scale(&vp, GeVec3_X(pp->v_cen, 2));
                        
			/* add vr, vt, vp to get v_gas */
			v_gas = GeVec3_Add(&vr,&vt);
			v_gas = GeVec3_Add(&v_gas,&vp);
                
			break;
                }        
		case GEOM_REC3D:
			#if USE_LVG
			/* Analytic expression for LVG expanding cloud */
			v_gas = GeVec3_Sub(&zone->root->voxel.cen, pos);
			double radius = GeVec3_Mag(&v_gas); /* [pc] */
			v_gas = GeVec3_Normalize(&v_gas); /* Direction vector */
			v_gas = GeVec3_Scale(&v_gas, 100.0e3 * radius); /* 100 km/s/pc velocity gradient */
			#elif USE_CONST
			v_gas = GeVec3_Sub(&zone->root->voxel.cen, pos);
			v_gas = GeVec3_Normalize(&v_gas); /* Direction vector */
			v_gas = GeVec3_Scale(&v_gas, velo);
			#else
			/* Obtain gas velocity from model grid */
			v_gas = pp->v_cen;
			#endif
			break;

		case GEOM_CYL3D:{
			GeVec3_d vRc = *pos;
			vRc.x[2] = 0.;
			vRc = GeVec3_Normalize(&vRc);

                        GeVec3_d vp;
                        GeVec3_d CylPos;
                        CylPos = GeVec3_Cart2Cyl( pos);
			if( CylPos.x[0] == 0.0 ){
				GeVec3_X(vp,0) = 0.0;
				GeVec3_X(vp,1) = 1.0;
				GeVec3_X(vp,2) = 0.0;
			}
			else{
				GeVec3_X(vp,0) = -pos->x[1];
				GeVec3_X(vp,1) = pos->x[0];
				GeVec3_X(vp,2) = 0.0;	
				vp = GeVec3_Normalize(&vp);
			}

			GeVec3_d vz;
			GeVec3_X(vz,0) = 0.0;
			GeVec3_X(vz,1) = 0.0;
			GeVec3_X(vz,2) = 1.0;
			
			/* Project directional velocity onto the normalized vector */
			vRc = GeVec3_Scale(&vRc, GeVec3_X(pp->v_cen, 0));
			vp = GeVec3_Scale(&vp, GeVec3_X(pp->v_cen, 1));
			vz = GeVec3_Scale(&vz, GeVec3_X(pp->v_cen, 2));
			/* add vr, vt, vp to get v_gas */
			v_gas = GeVec3_Add(&vRc,&vz);
			v_gas = GeVec3_Add(&v_gas,&vp);

			break;
                }
		default: /* Shouldn't reach here */
			Deb_ASSERT(0);
	}

	return v_gas;
}

/*----------------------------------------------------------------------------*/

GeVec3_d SpPhys_GetBgas(const GeVec3_d *pos, const Zone *zone)
/* Retriev gas velocity at pos */
{
	#define USE_LVG 0
	#define USE_CONST 0

	GeVec3_d b_gas;
        GeVec3_d br, bp, bt;
        GeVec3_d bRc, bz;

	#if USE_CONST
	double velo = -0.2e3; // [km/s]
	#endif

	/* Get zone physics */
	SpPhys *pp = zone->data;

	switch(zone->voxel.geom) {
		case GEOM_SPH1D:
			#if USE_LVG
			double radius = GeVec3_Mag(pos);
			b_gas = GeVec3_Normalize(pos);
			b_gas = GeVec3_Scale(&b_gas, 100.0e3 * radius); /* 100 km/s/pc velocity gradient */
			#elif USE_CONST
			b_gas = GeVec3_Normalize(pos);
			b_gas = GeVec3_Scale(&b_gas, velo); /* 100 km/s/pc velocity gradient */
			#else
			/* Project radial velocity onto position vector */
			b_gas = GeVec3_Normalize(pos);
			b_gas = GeVec3_Scale(&b_gas, GeVec3_X(pp->b_cen, 0));
			#endif
			break;
			
		case GEOM_SPH3D:
			
			#if USE_LVG
			/* Analytic expression for LVG expanding cloud */
			b_gas = GeVec3_Sub(&zone->root->voxel.cen, pos);
			double radius = GeVec3_Mag(&b_gas); /* [pc] */
			b_gas = GeVec3_Normalize(&b_gas); /* Direction vector */
			b_gas = GeVec3_Scale(&b_gas, 100.0e3 * radius); /* 100 km/s/pc velocity gradient */
			#elif USE_CONST
			GeVec3_X(b_gas,0) = 0.0;
			GeVec3_X(b_gas,1) = 0.0;
			GeVec3_X(b_gas,2) = 1.0;
			#else

			/* Project radial velocity onto position vector */
			/* get normalized vr, vt vp */
			br = GeVec3_Normalize(pos);
			if( fabs(pos->x[0])<1e-6 && fabs(pos->x[1])<1e-6 ){
				GeVec3_X(bp,0) = 0.0;
				GeVec3_X(bp,1) = 1.0;
				GeVec3_X(bp,2) = 0.0;
			}
			else{
				GeVec3_X(bp,0) = -pos->x[1];
				GeVec3_X(bp,1) = pos->x[0];
				GeVec3_X(bp,2) = 0.0;	
				bp = GeVec3_Normalize(&bp);
			}
			bt = GeVec3_CrossProd(&bp,&br);
			/* Project directional velocity onto the normalized vector */
			br = GeVec3_Scale(&br, GeVec3_X(pp->b_cen, 0));
			bt = GeVec3_Scale(&bt, GeVec3_X(pp->b_cen, 1));
			bp = GeVec3_Scale(&bp, GeVec3_X(pp->b_cen, 2));
			/* add vr, vt, vp to get v_gas */
			b_gas = GeVec3_Add(&br,&bt);
			b_gas = GeVec3_Add(&b_gas,&bp);
// 			if ( GeVec3_X(b_gas,2) < 0. ) printf("b_z < 0\n");
			#endif

			break;

		case GEOM_REC3D:
			#if USE_LVG
			/* Analytic expression for LVG expanding cloud */
			b_gas = GeVec3_Sub(&zone->root->voxel.cen, pos);
			double radius = GeVec3_Mag(&b_gas); /* [pc] */
			b_gas = GeVec3_Normalize(&b_gas); /* Direction vector */
			b_gas = GeVec3_Scale(&b_gas, 100.0e3 * radius); /* 100 km/s/pc velocity gradient */
			#elif USE_CONST
			b_gas = GeVec3_Sub(&zone->root->voxel.cen, pos);
			b_gas = GeVec3_Normalize(&b_gas); /* Direction vector */
			b_gas = GeVec3_Scale(&b_gas, velo);
			#else
			/* Obtain gas velocity from model grid */
			b_gas = pp->b_cen;
			#endif
			break;

		case GEOM_CYL3D:
			
			/* Project radial velocity onto position vector */
			/* get normalized vr, vt vp */
			bRc = *pos;
			bRc.x[2] = 0.;
			bRc = GeVec3_Normalize(&bRc);

			if( fabs(pos->x[0])<1e-6 && fabs(pos->x[1])<1e-6 ){
				GeVec3_X(bp,0) = 0.0;
				GeVec3_X(bp,1) = 1.0;
				GeVec3_X(bp,2) = 0.0;
			}
			else{
				GeVec3_X(bp,0) = -pos->x[1];
				GeVec3_X(bp,1) = pos->x[0];
				GeVec3_X(bp,2) = 0.0;	
				bp = GeVec3_Normalize(&bp);
			}
			
			GeVec3_X(bz,0) = 0.0;
			GeVec3_X(bz,1) = 0.0;
			GeVec3_X(bz,2) = 1.0;
			
			/* Project directional velocity onto the normalized vector */
			bRc = GeVec3_Scale(&bRc, GeVec3_X(pp->b_cen, 0));
			bp = GeVec3_Scale(&bp, GeVec3_X(pp->b_cen, 1));
			bz = GeVec3_Scale(&bz, GeVec3_X(pp->b_cen, 2));
			
			/* add vr, vt, vp to get v_gas */
			b_gas = GeVec3_Add(&bRc,&bz);
			b_gas = GeVec3_Add(&b_gas,&bp);

			break;
		
		default: /* Shouldn't reach here */
			Deb_ASSERT(0);
			
			
	}

	return b_gas;
}

/*----------------------------------------------------------------------------*/

GeVec3_d SpPhys_GetVfunc(const GeRay *ray, double dt, const Zone *zone)
/* Calculate projected velocity of ray in zone at offset dt. v_func is v_los
 * minus projected gas velocity -- this offsets v_los so that minus velocity
 * is blueshift and plus velocity is redshift. */
{
	GeRay offset_ray = GeRay_Inc(ray, dt);
	GeVec3_d v_gas = SpPhys_GetVgas(&offset_ray.e, zone);

	//return v_los - GeVec3_DotProd(&v_gas, &ray->d);
	return v_gas;
}



/*----------------------------------------------------------------------------*/

GeVec3_d SpPhys_GetBfunc(const GeRay *ray, double dt, const Zone *zone)
/* Calculate projected velocity of ray in zone at offset dt. v_func is v_los
 * minus projected gas velocity -- this offsets v_los so that minus velocity
 * is blueshift and plus velocity is redshift. */
{
	GeRay offset_ray = GeRay_Inc(ray, dt);
	GeVec3_d b_gas = SpPhys_GetBgas(&offset_ray.e, zone);

	//return v_los - GeVec3_DotProd(&v_gas, &ray->d);
	return b_gas;
}

/*----------------------------------------------------------------------------*/

double SpPhys_GetVfac(const GeRay *ray, double dt, double v_los, const Zone *zone, int debug)
{
	SpPhys *pp = zone->data;
	double vfac = 0;

	//debug
	//return Num_GaussNormal(SpPhys_GetVfunc(ray, dt, v_los, zone), pp->width);

	/* Number of steps is
	 * 	delta_v / width
	 * 
	 * where delta_v = |v_1 - v_1|
	 *       width = local gaussian line width
	 */
	GeVec3_d v_0 = SpPhys_GetVfunc(ray, 0.0, zone);
	GeVec3_d v_1 = SpPhys_GetVfunc(ray, dt, zone);

 	//Deb_PRINT("dt=%e, v_0=%e %e %e, v_1=%e %e %e\n", dt, v_0.x[0], v_0.x[1], v_0.x[2], v_1.x[0], v_1.x[1], v_1.x[2]);

	size_t n_step = Num_MAX((size_t)(GeVec3_Mag2(&v_1, &v_0) / pp->width), 1);

 	//Deb_PRINT("n_step=%d\n", (int)n_step);
        //Deb_PRINT("n_step=%d\n", (size_t)(GeVec3_Mag2(&v_1, &v_0) / pp->width));

	Deb_ASSERT(n_step > 0); /* Just in case */

	for(size_t i = 0; i < n_step; i++) {
		/* Check for velocity differences within each step
		 * in case there are large variations */
		double s_0 = dt * (double)i / (double)n_step;
		double s_1 = dt * (double)(i + 1) / (double)n_step;
		v_0 = SpPhys_GetVfunc(ray, s_0, zone);
		v_1 = SpPhys_GetVfunc(ray, s_1, zone);
		size_t n_avg = Num_MAX( (size_t)(GeVec3_Mag2(&v_1, &v_0) / pp->width) , 1);

		Deb_ASSERT(n_avg > 0); /* Just in case */

		//Deb_PRINT("  n_avg=%d\n", (int)n_avg);

		/* Average line profile over n_avg */
		for(size_t j = 0; j < n_avg; j++) {
			double s = s_0 + (s_1 - s_0) * ((double)j + 0.5) / (double)n_avg;
			GeVec3_d v = SpPhys_GetVfunc(ray, s, zone);
			vfac += Num_GaussNormal(v_los - GeVec3_DotProd(&v, &ray->d), pp->width) / (double)n_avg;
		}
	}
	//debug
	#if 0
		Deb_PRINT("zone=<%lu, %lu, %lu>, r=%g pc, dir=<%g, %g, %g>, v=%g km/s\n",
			(unsigned long)GeVec3_X(zone->index, 0), (unsigned long)GeVec3_X(zone->index, 1), (unsigned long)GeVec3_X(zone->index, 2),
			GeVec3_Mag2(&zone->parent->voxel.cen, &zone->voxel.cen),
			GeRay_D(*ray, 0), GeRay_D(*ray, 1), GeRay_D(*ray, 2),
			GeVec3_Mag(&pp->v_cen) / 1000.0);
	#endif

	//Deb_PRINT("getvfac done\n");

	return vfac / (double)n_step;
}

/*----------------------------------------------------------------------------*/

GeVec3_d SpPhys_GetVfac2(const GeRay *ray, double dt, const Zone *zone, int debug)
{
	//debug
        //SpPhys *pp = zone->data;
	//return Num_GaussNormal(SpPhys_GetVfunc(ray, dt, v_los, zone), pp->width);

	/* Number of steps is
	 * 	delta_v / width
	 *
	 * where delta_v = |v_1 - v_1|
	 *       width = local gaussian line width
	 */
	GeVec3_d v_0 = SpPhys_GetVfunc(ray, 0.0, zone);
	GeVec3_d v_1 = SpPhys_GetVfunc(ray, dt, zone);

// 	Deb_PRINT("dt=%g, v_0=%g, v_1=%g\n", dt, v_0, v_1);

	size_t n_step;
        double V_Mag = GeVec3_Mag(&v_0);
        if (V_Mag == 0.)
                n_step = 1;
        else
                n_step = Num_MAX((size_t)(GeVec3_Mag2(&v_1, &v_0) / (V_Mag * 0.01)), 1);

// 	Deb_PRINT("n_step=%d\n", (int)n_step);

	Deb_ASSERT(n_step > 0); /* Just in case */
	//Deb_PRINT("n_step=%d\n", (int)n_step);
        size_t n_avg;
        GeVec3_d v_field = GeVec3_INIT(0.,  0.,  0.);
	for(size_t i = 0; i < n_step; i++) {
		/* Check for velocity differences within each step
		 * in case there are large variations */
		double s_0 = dt * (double)i / (double)n_step;
		double s_1 = dt * (double)(i + 1) / (double)n_step;
		v_0 = SpPhys_GetVfunc(ray, s_0, zone);
		v_1 = SpPhys_GetVfunc(ray, s_1, zone);
                
                V_Mag = GeVec3_Mag(&v_0);
                if (V_Mag == 0.)
                        n_avg = 1;
                else
                        n_avg = Num_MAX((size_t)(GeVec3_Mag2(&v_1, &v_0) / (V_Mag*0.01)), 1);

		Deb_ASSERT(n_avg > 0); /* Just in case */

		//Deb_PRINT("  n_avg=%d\n", (int)n_step);

		/* Average line profile over n_avg */
		for(size_t j = 0; j < n_avg; j++) {
			double s = s_0 + (s_1 - s_0) * ((double)j + 0.5) / (double)n_avg;
			GeVec3_d v = SpPhys_GetBfunc(ray, s, zone);
			GeVec3_X(v_field,0) += GeVec3_X(v,0);
			GeVec3_X(v_field,1) += GeVec3_X(v,1);
			GeVec3_X(v_field,2) += GeVec3_X(v,2);
		}
	}
	v_field = GeVec3_Scale( &v_field, 1. / ( n_avg * n_step ));
	//debug
	#if 0
		printf("zone=<%lu, %lu, %lu>, r=%g pc, dir=<%g, %g, %g>, v=%g km/s\n",
			(unsigned long)GeVec3_X(zone->index, 0), (unsigned long)GeVec3_X(zone->index, 1), (unsigned long)GeVec3_X(zone->index, 2),
			GeVec3_Mag2(&zone->parent->voxel.cen, &zone->voxel.cen),
			GeRay_D(*ray, 0), GeRay_D(*ray, 1), GeRay_D(*ray, 2),
			GeVec3_Mag(&pp->v_cen) / 1000.0);
	#endif

	//Deb_PRINT("getvfac done\n");

	return v_field;
}

/*----------------------------------------------------------------------------*/

GeVec3_d SpPhys_GetBfac(const GeRay *ray, double dt, const Zone *zone, int debug)
{
	//debug
        //SpPhys *pp = zone->data;
	//return Num_GaussNormal(SpPhys_GetVfunc(ray, dt, v_los, zone), pp->width);

	/* Number of steps is
	 * 	delta_v / width
	 * 
	 * where delta_v = |v_1 - v_1|
	 *       width = local gaussian line width
	 */
	GeVec3_d b_0 = SpPhys_GetBfunc(ray, 0.0, zone);
	GeVec3_d b_1 = SpPhys_GetBfunc(ray, dt, zone);

// 	Deb_PRINT("dt=%g, v_0=%g, v_1=%g\n", dt, v_0, v_1);
        
        double B_Mag = GeVec3_Mag(&b_0);
        size_t n_step;
        if ( B_Mag == 0. )
                n_step = 1;
        else
                n_step = Num_MAX((size_t)(GeVec3_Mag2(&b_1, &b_0) / ( B_Mag * 0.01 )), 1);

// 	Deb_PRINT("n_step=%d\n", (int)n_step);

	Deb_ASSERT(n_step > 0); /* Just in case */
	//Deb_PRINT("n_step=%d\n", (int)n_step);
        size_t n_avg;
        GeVec3_d b_field = GeVec3_INIT(0.,  0.,  0.);
	for(size_t i = 0; i < n_step; i++) {
		/* Check for velocity differences within each step
		 * in case there are large variations */
		double s_0 = dt * (double)i / (double)n_step;
		double s_1 = dt * (double)(i + 1) / (double)n_step;
		b_0 = SpPhys_GetBfunc(ray, s_0, zone);
		b_1 = SpPhys_GetBfunc(ray, s_1, zone);
                B_Mag = GeVec3_Mag(&b_0);
                if ( B_Mag == 0. )
                        n_avg = 1;
                else
                        n_avg = Num_MAX((size_t)(GeVec3_Mag2(&b_1, &b_0) / (B_Mag * 0.01)), 1);
		Deb_ASSERT(n_avg > 0); /* Just in case */
		//Deb_PRINT("  n_avg=%d\n", (int)n_step);
		/* Average line profile over n_avg */
		for(size_t j = 0; j < n_avg; j++) {
			double s = s_0 + (s_1 - s_0) * ((double)j + 0.5) / (double)n_avg;
			GeVec3_d b = SpPhys_GetBfunc(ray, s, zone);
			GeVec3_X(b_field,0) += GeVec3_X(b,0);
			GeVec3_X(b_field,1) += GeVec3_X(b,1);
			GeVec3_X(b_field,2) += GeVec3_X(b,2);
		}
	}
	b_field = GeVec3_Scale(&b_field, 1./(n_avg*n_step));
	//debug
	#if 0
		printf("zone=<%lu, %lu, %lu>, r=%g pc, dir=<%g, %g, %g>, v=%g km/s\n",
			(unsigned long)GeVec3_X(zone->index, 0), (unsigned long)GeVec3_X(zone->index, 1), (unsigned long)GeVec3_X(zone->index, 2),
			GeVec3_Mag2(&zone->parent->voxel.cen, &zone->voxel.cen),
			GeRay_D(*ray, 0), GeRay_D(*ray, 1), GeRay_D(*ray, 2),
			GeVec3_Mag(&pp->v_cen) / 1000.0);
	#endif

	//Deb_PRINT("getvfac done\n");

	return b_field;
}

/*----------------------------------------------------------------------------*/

double _SpPhys_BoltzRatio(const Molec *mol, size_t up, size_t lo, double T_k)
/* Calculate the Boltzmann ratio = (g_u / g_l) * exp(-(E_u - E_l) / (k * T_k))
 */
{
	static double k = PHYS_CONST_MKS_BOLTZK;
	double
		g_u = mol->lev[up]->g,
		g_l = mol->lev[lo]->g,
		E_u = mol->lev[up]->E,
		E_l = mol->lev[lo]->E;

	Deb_ASSERT((up > lo) && (T_k > 0));

	return (g_u / g_l) * exp(-(E_u - E_l) / (k * T_k));
}

/*----------------------------------------------------------------------------*/

double SpPhys_CalcLineWidth(const SpPhys *pp)
/* Calculate total line width accounting for both thermal and
 * turbulent motion, where V_t is the RMS turbulent velocity. */
{
	double V_th = Phys_ThermLineWidth(pp->T_k, pp->mol->weight);
	double V_t = pp->V_t;

	return sqrt(V_th * V_th + V_t * V_t);
}

/*----------------------------------------------------------------------------*/

#define STACK(i, j)\
	stack[(j) + nlev * (i)]

void SpPhys_PushPopsStack(double *stack, size_t nstack, double *pops, size_t nlev)
{
	/* Stack:
	 *
	 *	1*	2	3
	 *-----------------------------
	 *	11	21	31
	 *	12	22	32
	 *	13	23	33
	 *	14	24	34
	 */

	for(size_t i = 0; i < nlev; i++) {
		for(size_t j = 0; j < nstack; j++) {
			size_t dest_id = (nstack - 1) - j;
			if(dest_id == 0) {
				STACK(dest_id, i) = pops[i];
			}
			else {
				STACK(dest_id, i) = STACK(dest_id - 1, i);
			}
		}
	}

	return;
}

/*----------------------------------------------------------------------------*/

double SpPhys_CalcPopsDiff(const double *stack, size_t nstack, size_t nlev, double minpop)
{
	double *pops = Mem_CALLOC(nstack, pops);
	double *diffs = Mem_CALLOC(nstack, diffs);
	double *max_diffs = Mem_CALLOC(nlev, max_diffs);

	for(size_t i = 0; i < nlev; i++) {
		double mean = 0;
		for(size_t j = 0; j < nstack; j++) {
			/* Load pops for sorting */
			pops[j] = STACK(j, i);

			/* Accumulate mean for later use */
			mean += pops[j];
		}
		Num_Qsort_d(pops, nstack);

		/* Calculate variance if smallest pops is
		 * greater than or euqal to minpop */
		if(pops[0] >= minpop) {
			/* Calculate mean */
			mean /= (double)nstack;

			/* Variance is relative difference from mean */
			for(size_t j = 0; j < nstack; j++) {
				diffs[j] = fabs(STACK(j, i) - mean) / mean;
			}
			Num_Qsort_d(diffs, nstack);

			/* Maximum variance for this level */
			max_diffs[i] = diffs[nstack - 1];
		}
	}
	/* Find Maximum variance of entire stack */
	Num_Qsort_d(max_diffs, nlev);
	double max_diff = max_diffs[nlev - 1];
	
	free(pops);
	free(diffs);
	free(max_diffs);

	return max_diff;
}

#undef STACK











