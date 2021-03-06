#ifndef __NUMERICAL_H__
#define __NUMERICAL_H__

#include <gsl/gsl_rng.h>

/* The numerical.c/.h interface should provide basic mathematical
 * constants and operations, such as PI, TWOPI, searching and sorting,
 * ... etc.
 */

#define Num_PI 3.141592653589793238462643383279502884197
#define Num_TWOPI  6.283185307179586
#define Num_HALFPI 1.5707963267948966

#define Num_MAX(a, b)\
	((a) > (b) ? (a) : (b))

#define Num_MIN(a, b)\
	((a) < (b) ? (a) : (b))

#define Num_SIGN(a)\
	((a) < 0 ? -1 : 1)

#define Num_ISNAN(a)\
	isnan(a) /* May not be portable */

#define Num_CUBE(a)\
	((a) * (a) * (a))

void *Num_Bsearch_d(double key, double *base, size_t n);
size_t Num_GetInterpIdx(double x, const double *xa, size_t n, size_t m);
double Num_InterpLinear(double x, double x_0, double x_1, double y_0, double y_1);
double Num_InterpPoly(double x, const double *xa, const double *ya, size_t n, size_t m);
double Num_InterpPoly2d(double x1, double x2, const double *x1a, const double *x2a,
	const double *ya, size_t n1, size_t n2, size_t m1, size_t m2);

void Num_QRDecompSolve(double *A, size_t M, size_t N, const double *b, double *x);
void Num_SVDecompSolve(double *A, size_t M, size_t N, const double *b, double *x);
void Num_LUDecompSolve(double *A, size_t N, const double *b, double *x);

void Num_EigenSolver(double *A, size_t N, double *eigen_value);

void Num_Qsort_d(double array[], size_t n);
void Num_QuadraticRoots(double a, double b, double c, double *x1, double *x2);
void NumFFT_Xform2d(double *arr, size_t idim, size_t jdim, double *real, double *imag);
void NumFFT_Swap1d(double *arr, size_t n);
void NumFFT_Swap2d(double *arr, size_t idim, size_t jdim);
double Num_GaussNormal(double x, double width);
double Num_GaussNormal2(double x, double width);
void Num_RanDir3D(gsl_rng *rng, double *cost, double *sint, double *phi);
void Num_QRanDir3D(const double *QRN, double *cost, double *sint, double *phi);

#endif
