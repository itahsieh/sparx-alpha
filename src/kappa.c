#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "error.h"
#include "memory.h"
#include "numerical.h"
#include "kappa.h"
#include "debug.h"

#include "physics.h"
#define Kap_LIGHTC_MKS PHYS_CONST_MKS_LIGHTC

static double Kap_FromFreq_table(const Kappa *kap, double freq);
static double Kap_FromFreq_plaw(const Kappa *kap, double freq);

/*----------------------------------------------------------------------------*/

Kappa *Kap_New_Powerlaw(double freq_0, double kappa_0, double index)
{

	Kappa *kap = Mem_CALLOC(1, kap);

	
	Deb_ASSERT(freq_0 > 0);

	kap->type = Kap_POWERLAW;

	kap->name = Mem_Sprintf("powerlaw,%10.3e,%10.3e,%10.3e", freq_0, kappa_0, index);

	kap->freq_0 = freq_0;
	kap->kappa_0 = kappa_0;
	kap->index = index;

	return kap;
}

/*----------------------------------------------------------------------------*/

Kappa *Kap_New_Table(const char *name, const char *fname, FILE *fp)
{
	Kappa *kap = Mem_CALLOC(1, kap);
	char buffer[BUFSIZ];
	double c = Kap_LIGHTC_MKS, /* Speed of light in vacuum in MKS units */
	       um = 1.0e-6; /* 1 micron in meters */
	size_t row_id, ncol;

	/* Set opacity type */
	kap->type = Kap_TABLE;
	kap->name = Mem_Sprintf("table,%s", name);
	kap->fname = Mem_STRDUP(fname);
	kap->nrows = 0;

	while(fgets(buffer, BUFSIZ, fp)) {
		kap->nrows += 1;
		kap->freq = Mem_REALLOC(kap->nrows, kap->freq);
		kap->lambda = Mem_REALLOC(kap->nrows, kap->lambda);
		kap->kappa = Mem_REALLOC(kap->nrows, kap->kappa);

		kap->loglam = Mem_REALLOC(kap->nrows, kap->loglam);
		kap->logkap = Mem_REALLOC(kap->nrows, kap->logkap);

		row_id = kap->nrows - 1;
		ncol = (size_t)sscanf(buffer, "%lf %lf", &kap->lambda[row_id], &kap->kappa[row_id]);
		Deb_ASSERT(ncol == 2);
		Deb_ASSERT(kap->lambda[row_id] > 0);

		kap->lambda[row_id] *= um; /* um -> m */
		kap->freq[row_id] = c / kap->lambda[row_id]; /* nu = c / lambda */
		kap->kappa[row_id] *= 0.1; /* cm^2 g^-1 -> m^2 kg^-1 */

		/* Calculate log10 values (for interpolation) */
		Deb_ASSERT(kap->lambda[row_id] > 0);
		Deb_ASSERT(kap->freq[row_id] > 0);
		Deb_ASSERT(kap->kappa[row_id] > 0);

		kap->loglam[row_id] = log10(kap->lambda[row_id]);
		kap->logkap[row_id] = log10(kap->kappa[row_id]);
	};

	return kap;
}

/*----------------------------------------------------------------------------*/

void Kap_Free(void *ptr)
{
	Kappa *kap = ptr;

	#define FREE(p)\
		if(kap->p) { free(kap->p); }

	if(kap) {
		FREE(name);
		FREE(fname);
		FREE(freq);
		FREE(lambda);
		FREE(kappa);
		FREE(loglam);
		FREE(logkap);
		free(kap);
	}

	return;
}

/*----------------------------------------------------------------------------*/

double Kap_FromFreq(const Kappa *kap, double freq)
{       
	switch(kap->type) {
		case Kap_POWERLAW:
			return Kap_FromFreq_plaw(kap, freq);

		case Kap_TABLE:
			return Kap_FromFreq_table(kap, freq);

		default:
			Deb_ASSERT(0);
	}

	return HUGE_VAL;
}

/*----------------------------------------------------------------------------*/

static double Kap_FromFreq_plaw(const Kappa *kap, double freq)
{
	Deb_ASSERT(kap->type == Kap_POWERLAW);
	Deb_ASSERT((kap->freq_0 > 0) && (freq > 0));

	return kap->kappa_0 * pow(freq / kap->freq_0, kap->index);
}

/*----------------------------------------------------------------------------*/

static double Kap_FromFreq_table(const Kappa *kap, double freq)
{
	Deb_ASSERT(kap->type == Kap_TABLE);
	Deb_ASSERT(freq > 0);

	double loglam = log10(Kap_LIGHTC_MKS / freq), logkap;
        
        static int over_max = 0;
        static int under_min = 0;

        
        int m = kap->nrows;
        double min_loglam = kap->loglam[0];
        double max_loglam = kap->loglam[m-1];
        Deb_ASSERT( min_loglam < max_loglam);
        if (loglam > max_loglam ){
                if (!over_max){
                        printf("Warning : Wavelength exceeds the table of kappa, lambda adopts the upper limit.\n");
                        //printf("log(lambda) = %f, max log(lambda) = %f\n", loglam, max_loglam);
                }
                over_max = 1;
                loglam = max_loglam;
        }
        else if(loglam < min_loglam){
                if (!under_min){
                        printf("Warning : Wavelength falls short of the table, lambda adopts the lower limit.\n");
                        //printf("log(lambda) = %f, min log(lambda) = %f\n", loglam, min_loglam);
                }
                under_min = 1;
                loglam = min_loglam;
        }

	logkap = Num_InterpPoly(loglam, kap->loglam, kap->logkap, kap->nrows, (size_t)3);

	return pow(10.0, logkap);
}

/*----------------------------------------------------------------------------*/







