#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include "physics.h"
#include "debug.h"

/*----------------------------------------------------------------------------*/

double Phys_DoppShiftFreq(double freq, double vel)
/* Calculate the Doppler-shifted frequency of an object moving at velocity
 * vel and source frequency freq. Vel MUST NOT BE TOO CLOSE TO C! */
{
	double c = PHYS_CONST_MKS_LIGHTC;

	Deb_ASSERT(vel < c);

	return sqrt((1.0 + vel / c) / (1.0 - vel / c)) * freq;
}

/*----------------------------------------------------------------------------*/

double Phys_BrightTemp(double nu, double I_nu)
/* Calculate birghtness temperature according to the Rayleigh-Jeans
 * law for thermal radiation
 * 	T_b = I_nu * c^2 / (2 * nu^2 * k) */
{
	static double
		c = PHYS_CONST_MKS_LIGHTC,
		k = PHYS_CONST_MKS_BOLTZK;

	return (c * c) * I_nu / (2.0 * (nu * nu) * k);
}

/*----------------------------------------------------------------------------*/

double Phys_RayleighJeans(double nu, double T_k)
{
	static double
		c = PHYS_CONST_MKS_LIGHTC,
		k = PHYS_CONST_MKS_BOLTZK;

	return 2.0 * (nu * nu) * k * T_k / (c * c);
}

/*----------------------------------------------------------------------------*/

#if 1
double Phys_PlanckFunc(double nu, double T_k)
/* Return the intensity of black-body radiation according to the Planck
   function:
	B_nu(T) = (2 * h * nu^3 / c^2) / (exp(h * nu / (k * T_k)) - 1)
   This version should prevent floating-point overflows by taking the Wien
   approximation explicitly at high energy.
*/
{
	static const double
		h = PHYS_CONST_MKS_PLANCKH,
		c = PHYS_CONST_MKS_LIGHTC,
		k = PHYS_CONST_MKS_BOLTZK;

	if(T_k <= DBL_EPSILON) {
		return 0.0;
	}
	/* Explicitely take Wien approximation for h*nu >> k*T: */
	else if(h * nu >= 100.0 * k * T_k) {
		return 2.0 * h * (nu * nu * nu / (c * c)) * exp(-h * nu / (k * T_k));
	}
	else {
		return (2.0 * h * nu * nu * nu / (c * c)) / (exp(h * nu / k / T_k) - 1.0);
	}
}

#else

double Phys_PlanckFunc(double nu, double T_k)
/* Return the intensity of black-body radiation according to the Planck
   function:

	B_nu(T) = (2 * h * nu^3 / c^2) / (exp(h * nu / (k * T_k)) - 1)

*/
{
	static const double
		h = PHYS_CONST_MKS_PLANCKH,
		c = PHYS_CONST_MKS_LIGHTC,
		k = PHYS_CONST_MKS_BOLTZK;

	return (2.0 * h * (nu * nu * nu) / (c * c)) / (exp(h * nu / k / T_k) - 1.0);
}
#endif

/*----------------------------------------------------------------------------*/

double Phys_ThermLineWidth(double T_k, double m_a)
/* Calculate thermal line width of molecule with molecular mass m_a
 * at kinetic temperature T_k. This is the gaussian sigma multiplied
 * by sqrt(2), which can be converted to FWHM by multiplying 2 * sqrt(log(2)).
 * See e.g. Rybicki & Lightman (1979) p. 288 for details.
 */
{
	static const double k = PHYS_CONST_MKS_BOLTZK;

	return sqrt(2.0 * k * T_k / m_a);
}

/*----------------------------------------------------------------------------*/

void Phys_SecsToHMS(int t, int *hr, int *min, int *sec)
{
	*hr = (int)floor(t / 3600.0);
	*min = (int)floor((t % 3600) / 60.0);
	*sec = (t % 3600) % 60;

	return;
}

/*----------------------------------------------------------------------------*/

char *Phys_SecsToHMS_Str(int t)
{
	static char string[BUFSIZ] = "";
	int hr, min, sec;

	Phys_SecsToHMS(t, &hr, &min, &sec);
	snprintf(string, (size_t)BUFSIZ, "%03d:%02d:%02d", hr, min, sec);

	return string;
}







