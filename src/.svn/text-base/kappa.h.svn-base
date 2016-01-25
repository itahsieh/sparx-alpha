#ifndef __KAPPA_H__
#define __KAPPA_H__

#include <stdlib.h>

/* This API defines an opacity object which can calculate the opacity
 * at a certain frequency/wavelength based either on an analytic description
 * or interpolation from an external table */

enum {
	Kap_POWERLAW,
	Kap_TABLE
};

typedef struct Kappa {
	int type;
	char *name;

	/* Power law parameters */
	double freq_0, kappa_0, index;

	/* Table parameters */
	char *fname;
	size_t nrows;
	double *freq; /* Frequency column */
	double *lambda, *loglam; /* Wavelength column */
	double *kappa, *logkap; /* Opacities column */
} Kappa;

Kappa *Kap_New_Powerlaw(double freq_0, double kappa_0, double index);
Kappa *Kap_New_Table(const char *name, const char *fname, FILE *fp);
void Kap_Free(void *ptr);

double Kap_FromFreq(const Kappa *kap, double freq);

#endif














