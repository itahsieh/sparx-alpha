#ifndef __PHYSICS_H__
#define __PHYSICS_H__

#include <time.h>
#include "constants.h"

/* constants in MKS units */
#define PHYS_CONST_MKS_LIGHTC	    CONSTANTS_MKS_LIGHT_C /* speed of light in vacuum [m s^-1] */
#define PHYS_CONST_MKS_LIGHTCSQR    CONSTANTS_MKS_LIGHT_CSQR /* speed of light in vacuum [m s^-1] */
#define PHYS_CONST_MKS_PLANCKH	    CONSTANTS_MKS_PLANCK_H /* Planck's constant [J s] */
#define PHYS_CONST_MKS_PLANCKHBAR   CONSTANTS_MKS_PLANCK_HBAR /* Planck's constant divided by 2*PI [J s] */
#define PHYS_CONST_MKS_BOLTZK	    CONSTANTS_MKS_BOLTZ_K /* Boltzmann's constant [J K^-1] */
#define PHYS_CONST_MKS_GRAVG	    CONSTANTS_MKS_GRAV_G /* gravitational constant [N m^2/kg^2] */
#define PHYS_CONST_MKS_STEFBOLTZ    CONSTANTS_MKS_STEFBOLTZ_SIGMA /* Stefan-Boltzmann constant [W m^-2 K^4] */
#define PHYS_CONST_MKS_RSUN	    CONSTANTS_MKS_R_SUN /* solar radius [m] (see Brown et al. (1998), ApJ 500 L195) */
#define PHYS_CONST_MKS_MSUN	    CONSTANTS_MKS_M_SUN /* solar mass [kg] */
#define PHYS_CONST_MKS_LSUN	    CONSTANTS_MKS_L_SUN /* solar luminosity [W] */
#define PHYS_CONST_MKS_WHYDROGEN    CONSTANTS_MKS_W_HYDROGEN
#define PHYS_CONST_MKS_ECHARGE	    CONSTANTS_MKS_CHARGE_E /* elementary charge */

/* constants in CGS units */
#define PHYS_CONST_CGS_GRAVG	    6.67e-8 /* gravitational constant [ gm^-1 cm^3 s^-2] */
#define PHYS_CONST_CGS_LIGHTC	    3.00e10 /* speed of light [ cm s^-1 ] */

/* constants in units common to both MKS and CGS */
#define PHYS_CONST_PI \
        3.14159265358979323846264338327950288419716939937510
#define PHYS_CONST_TWOPI	    6.28318530
#define PHYS_CONST_2S2L2	    2.35482 /* Factor to convert from Gaussian 1 sigma width to FWHM */
#define PHYS_CONST_TCMB		    2.725 /* CMB temperature [K] (see Mather et al. (1999), ApJ 512 511) */
#define PHYS_CONST_GAS2DUST	    100.0 /* galactic gas:dust ratio (see Devereux et al. (1990) ApJ 359 42) */

/*
 * various unit conversion factors
 */
/* angles -> [rad] */
#define PHYS_UNIT_DEG		    1.74532925000000E-02
#define PHYS_UNIT_AMIN		    2.90888208333333E-04
#define PHYS_UNIT_ASEC		    4.84813680555556E-06

/* frequency -> [Hz] */
#define PHYS_UNIT_GHZ		    1.0E9
#define PHYS_UNIT_MHZ		    1.0E6

/* time -> [s] */
#define PHYS_UNIT_YEAR		    31557600 /* number of seconds in a Julian year */
#define PHYS_UNIT_DAY		    86400
#define PHYS_UNIT_HOUR		    3600
#define PHYS_UNIT_MINUTE	    60
#define PHYS_UNIT_NS		    1.0E-9

/* temperature -> [K] */
#define PHYS_UNIT_K		    1.0

/* units in MKS */
/* flux density */
#define PHYS_UNIT_MKS_JY	    1.0E-26 /* [W m^-2 Hz^-1] */

/* Mass */
#define PHYS_UNIT_MKS_AMU	    1.66053873E-27 /* atomic mass unit [kg] */
#define PHYS_UNIT_MKS_KG	    1.0

/* lengths -> [m] */
#define PHYS_UNIT_MKS_AU	    1.49598E11 /* astronomical unit [m] */
#define PHYS_UNIT_MKS_LY	    9.4605284E15 /* lightyear [m] */
#define PHYS_UNIT_MKS_PC	    3.08568025E16 /* parsec [m] */
#define PHYS_UNIT_MKS_KM	    1.0E3
#define PHYS_UNIT_MKS_M		    1.0
#define PHYS_UNIT_MKS_CM	    1.0E-2
#define PHYS_UNIT_MKS_MM	    1.0E-3
#define PHYS_UNIT_MKS_UM	    1.0E-6
#define PHYS_UNIT_MKS_CC	    1.0E-6
#define PHYS_UNIT_MKS_PERCC	    1.0E6
#define PHYS_UNIT_MKS_MM	    1.0E-3
#define PHYS_UNIT_MKS_UM	    1.0E-6
#define PHYS_UNIT_MKS_NM	    1.0E-9

/* units in CGS */
/* lengths -> [cm] */
#define PHYS_UNIT_CGS_PC	    3.08568025E18 /* parsec [cm] */
#define PHYS_UNIT_CGS_LY	    9.4605284E17 /* lightyear [cm] */
#define PHYS_UNIT_CGS_AU	    1.49598E13 /* astronomical unit [cm] */
#define PHYS_UNIT_CGS_KM	    1.0E5
#define PHYS_UNIT_CGS_MM	    1.0E-1
#define PHYS_UNIT_CGS_UM	    1.0E-4

double Phys_DoppShiftFreq(double nu, double vel);
double Phys_BrightTemp(double nu, double I_nu);
double Phys_RayleighJeans(double nu, double T_k);
double Phys_PlanckFunc(double nu, double T_k);
double Phys_ThermLineWidth(double T_k, double m_a);
void Phys_SecsToHMS(int t, int *hr, int *min, int *sec);
char *Phys_SecsToHMS_Str(int t);

#endif
