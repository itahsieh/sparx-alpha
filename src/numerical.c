#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_eigen.h>
#include <fftw3.h>

#include "debug.h"

#include "memory.h"
#include "numerical.h"

static int compar_d(const void *a, const void *b);

/*----------------------------------------------------------------------------*/

static int compar_d(const void *a, const void *b)
{
	double aa = *((const double *)a), bb = *((const double *)b);

	if(aa > bb)
		return 1;
	else if(fabs(aa - bb) <= DBL_EPSILON)
		return 0;
	else
		return -1;
}

/*----------------------------------------------------------------------------*/

void *Num_Bsearch_d(double key, double *base, size_t n)
{
	return bsearch(&key, base, n, sizeof(*base), compar_d);
}

/*----------------------------------------------------------------------------*/

size_t Num_GetInterpIdx(double x, const double *xa, size_t n, size_t m)
{
	size_t idx;

	/* Search for position of x within xa */
	idx = gsl_interp_bsearch(xa, x, (size_t)0, n - 1);

	/* Try to place idx at the beginning of the nearest
	 * block of m elements; if not possible then is either
	 * at xa[0] or xa[n - m].
	 *
	 * Special care must be taken since ((size_t)0 - 1)
	 * could overflow and turn into garbage.
	 */
	if(idx > (m - 1) / 2) {
		idx = Num_MIN(idx - (m - 1) / 2, n - m);
	}
	else {
		idx = 0;
	}

	return idx;
}

/*----------------------------------------------------------------------------*/

double Num_InterpLinear(double x, double x_0, double x_1, double y_0, double y_1)
/* linearly interpolate between two points:
 * 	y = y_0 + m * dx
 *
 *   where m = (y_1 - y_0) / (x_1 - x_0)
 */
{
	if(fabs(x_1 - x_0) > 0) {
		return y_0 + ((y_1 - y_0) / (x_1 - x_0)) * (x - x_0);
	}
	else {
		return y_0;
	}
}

/*----------------------------------------------------------------------------*/

double Num_InterpPoly(double x, const double *xa, const double *ya, size_t n, size_t m)
/* Do polynomial interpolation of arrays xa[] and ya[] of size n and order m */
{
	size_t idx;
	double result = 0.0;
	gsl_interp_accel *acc;
	gsl_spline *spline;

	idx = Num_GetInterpIdx(x, xa, n, m);

	/*init interpolation objects*/
	acc = gsl_interp_accel_alloc(); //init accelerator
	spline = gsl_spline_alloc(gsl_interp_polynomial, m); //init spline object
	gsl_spline_init(spline, &xa[idx], &ya[idx], m);

	/*execute interpolation*/
	result = gsl_spline_eval(spline, x, acc);

	/*cleanup*/
	gsl_spline_free(spline);
	gsl_interp_accel_free(acc);

	return result;
}

/*----------------------------------------------------------------------------*/

#define YA(a,b) (ya[(b)+n2*(a)])
#define YSUB(a,b) (ysub[(b)+m2*(a)])

double Num_InterpPoly2d(double x1, double x2, const double *x1a, const double *x2a,
	const double *ya, size_t n1, size_t n2, size_t m1, size_t m2)
/*do an order m1 x m2 2d polynomial interpolation of 2d array ya[] at (x1,x2)
 in coordinates given by x1a[] & x2a[]*/
{
	size_t idx1, idx2, i, j;
	double *x1sub = NULL, *x2sub = NULL;
	double y, *ysub = NULL, *ytmp = NULL;

	/* Allocate sub-block and coordinates */
	ysub = Mem_CALLOC(m1*m2, ysub);
	x1sub = Mem_CALLOC(m1, x1sub);
	x2sub = Mem_CALLOC(m2, x2sub);

	/* Allocate temporary storage for interp along columns */
	ytmp = Mem_CALLOC(m1, ytmp);

	/* Locate position of m1 x m2 sub-block of ya[] */
	idx1 = Num_GetInterpIdx(x1, x1a, n1, m1);
	idx2 = Num_GetInterpIdx(x2, x2a, n2, m2);

	/* Load sub-block */
	for(i = 0; i < m1; i++) {
		for(j = 0; j < m2; j++) {
			YSUB(i,j) = YA(idx1+i, idx2+j);
		}
	}

	/* Load coordinates */
	for(i = 0; i < m1; i++)
		x1sub[i] = x1a[idx1+i];

	for(i = 0; i < m2; i++)
		x2sub[i] = x2a[idx2+i];

	/* Step 1: interpolate along each column */
	for(i = 0; i < m1; i++) {
		ytmp[i] = Num_InterpPoly(x2, x2sub, &YSUB(i, 0), m2, m2);
	}

	/* Step 2: interpolate along row */
	y = Num_InterpPoly(x1, x1sub, ytmp, m1, m1);

	/* Cleanup */
	free(x1sub);
	free(x2sub);
	free(ysub);
	free(ytmp);

	return y;
}

#undef YA
#undef YSUB

/*----------------------------------------------------------------------------*/

void Num_QRDecompSolve(double *A, size_t M, size_t N, const double *b, double *x)
/* Solve the linear system A x = b by doing QR decomposition and
 * least-squares solution, where A is a general M by N matrix.
 *
 * Note on gsl_linalg_QR_* routines:
 * in gsl_linalg_QR_decomp(), dim(tau) must == \min(M, N) and
 * gsl_linalg_QR_lssolve(), dim(x) must == N, while dim(b) & dim(residual)
 * must == M
 */
{
	gsl_matrix_view
		QR = gsl_matrix_view_array(A, M, N);
	gsl_vector_const_view
		bb = gsl_vector_const_view_array(b, M);
	gsl_vector
		*tau = gsl_vector_alloc(Num_MIN(M, N)),
		*xx = gsl_vector_alloc(N),
		*residual = gsl_vector_alloc(M);

	gsl_linalg_QR_decomp(&QR.matrix, tau);
	gsl_linalg_QR_lssolve(&QR.matrix, tau, &bb.vector, xx, residual);

	for(size_t i = 0; i < N; i++)
		x[i] = gsl_vector_get(xx, i);

	gsl_vector_free(tau);
	gsl_vector_free(xx);
	gsl_vector_free(residual);

	return;
}

/*----------------------------------------------------------------------------*/

void Num_SVDecompSolve(double *A, size_t M, size_t N, const double *b, double *x)
/* Solve the linear system A * x = b using the Singular Value Decomposition
 * method. This version uses one-sided Jacobi orthogonalization, which requires
 * M >= N.
 */
{
	Deb_ASSERT(M >= N);

	size_t i;
	gsl_matrix_view
		U = gsl_matrix_view_array(A, M, N);
	gsl_vector_const_view
		bb = gsl_vector_const_view_array(b, M);
	gsl_matrix
		*V = gsl_matrix_alloc(N, N);
	gsl_vector
		*S = gsl_vector_alloc(N),
		*xx = gsl_vector_alloc(N);

	gsl_linalg_SV_decomp_jacobi(&U.matrix, V, S);
	gsl_linalg_SV_solve(&U.matrix, V, S, &bb.vector, xx);

	for(i = 0; i < N; i++)
		x[i] = gsl_vector_get(xx, i);

	return;
}

/*----------------------------------------------------------------------------*/

void Num_LUDecompSolve(double *A, size_t N, const double *b, double *x)
/* Solve the linear system A x = b by doing QR decomposition and
 * least-squares solution, where A is a general M by N matrix.
 *
 * Note on gsl_linalg_QR_* routines:
 * in gsl_linalg_QR_decomp(), dim(tau) must == \min(M, N) and
 * gsl_linalg_QR_lssolve(), dim(x) must == N, while dim(b) & dim(residual)
 * must == M
 */
{
        gsl_matrix_view
                LU = gsl_matrix_view_array(A, N, N);
        gsl_vector_const_view
                bb = gsl_vector_const_view_array(b, N);
        gsl_vector
                *xx = gsl_vector_alloc(N);
        gsl_permutation 
                *p = gsl_permutation_alloc (N);
        int s;

        gsl_linalg_LU_decomp(&LU.matrix, p, &s);
        
        gsl_linalg_LU_solve(&LU.matrix, p, &bb.vector, xx);

        for(size_t i = 0; i < N; i++)
                x[i] = gsl_vector_get(xx, i);

        gsl_vector_free(xx);
        gsl_permutation_free (p);

        return;
}

/*----------------------------------------------------------------------------*/

void Num_EigenSolver(double *A, size_t N, double *eigen_value)
{
        gsl_vector *eval = gsl_vector_alloc (N);
        gsl_matrix *evec = gsl_matrix_alloc (N, N);
        
        gsl_matrix_view 
                m = gsl_matrix_view_array(A, N, N);
        gsl_eigen_symm_workspace 
                * w = gsl_eigen_symmv_alloc(N);
                
        gsl_eigen_symmv (&m.matrix, eval, evec, w);

        gsl_eigen_symmv_free (w);

        gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

        for (int i = 0; i < N; i++){
                double eval_i = gsl_vector_get (eval, i);
                eigen_value[i] = eval_i;
                gsl_vector_view 
                        evec_i = gsl_matrix_column (evec, i);

                printf ("eigenvalue = %g\n", eval_i);
                //printf ("eigenvector = \n");
                //gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
        }

        gsl_vector_free (eval);
        gsl_matrix_free (evec);
        
        return;
}

/*----------------------------------------------------------------------------*/

void Num_Qsort_d(double array[], size_t n)
{
	qsort(array, n, sizeof(double), (int (*)(const void *, const void *))compar_d);

	return;
}

/*----------------------------------------------------------------------------*/

void Num_QuadraticRoots(double a, double b, double c, double *x1, double *x2)
/* 
Solve the quadratic equation a * x^2 + b * x + c = 0 and store the roots in
x1 and x2.

Note: In Numerical Recipes in Pascal (Press, etal. , Cambridge University
Press, 1989) - section 5.5, there is a better way to compute the roots of
the quadratic equation.

The normal way:  x = [-b (+/-) (b2 - 4ac).5]/2a 
If either a or c small then we get -b (+/-) b' with b' ~ b and this may
give a numerical error. A better way is: q = -0.5 [b +sgn(b) (b2 - 4ac) .5]
then  x1 = q/a;  x2 = c/q;

### Evaluates to nan for coplex roots! ###
*/
{
	double q;

	q = -0.5 * (b + Num_SIGN(b) * sqrt(b * b - 4.0 * a * c));

	*x1 = q / a;
	*x2 = c / q;

	return;
}

/*----------------------------------------------------------------------------*/

#define IN(a,b)\
	(fftw_in[(b)+jdim*(a)])
#define OUT(a,b)\
	(fftw_out[(b)+jdim*(a)])
#define ARR(a,b)\
	(arr[(b)+jdim*(a)])
#define REAL(a,b)\
	(real[(b)+jdim*(a)])
#define IMAG(a,b)\
	(imag[(b)+jdim*(a)])
#define DUP1(a,b)\
	(dup1[(b)+jdim*(a)])
#define DUP2(a,b)\
	(dup2[(b)+jdim*(a)])

void NumFFT_Xform2d(double *arr, size_t idim, size_t jdim, double *real, double *imag)
{
	size_t i, j;
	double *dup1, *dup2;
	fftw_plan p;
	fftw_complex *fftw_in, *fftw_out;

	/* Sanity check */
	Deb_ASSERT(idim * jdim > 0);

	dup1 = Mem_CALLOC(idim*jdim, dup1);
	dup2 = Mem_CALLOC(idim*jdim, dup2);
	fftw_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*idim*jdim);
	fftw_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*idim*jdim);
	Deb_ASSERT((fftw_in != NULL) && (fftw_out != NULL));

	/* FFTW doesn't seem to use size_t for array dimensions -- consider change to
	 * another library? GSL perhaps? */
	p = fftw_plan_dft_2d((int)idim, (int)jdim, fftw_in, fftw_out, FFTW_FORWARD, FFTW_MEASURE);

	/* Load arr into dup1 */
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim; j++) {
			DUP1(i, j) = ARR(i, j);
		}
	}

	/* FFTSwap dup1 */
	NumFFT_Swap2d(&DUP1(0, 0), idim, jdim);

	/* Load dup1 into real part of input and reset output */
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim; j++) {
			IN(i,j)[0] = DUP1(i,j);
			IN(i,j)[1] = 0.0;
			OUT(i,j)[0] = 0.0;
			OUT(i,j)[1] = 0.0;
		}
	}

	/*execute FFT*/
	fftw_execute(p);

	/* Load output: real->dup1 imag->dup2 */
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim; j++) {
			DUP1(i,j) = OUT(i,j)[0];
			DUP2(i,j) = OUT(i,j)[1];
		}
	}

	/* FFTSwap dup1 and dup2 */
	NumFFT_Swap2d(&DUP1(0, 0), idim, jdim);
	NumFFT_Swap2d(&DUP2(0, 0), idim, jdim);

	/* Load results to user array: dup1->real, dup2->imag */
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim; j++) {
			REAL(i,j) = DUP1(i,j);
			IMAG(i,j) = DUP2(i,j);
		}
	}

	/* Cleanup */
	free(dup1);
	free(dup2);
	fftw_free(fftw_in);
	fftw_free(fftw_out);
	fftw_destroy_plan(p);

	return;
}

#undef IN
#undef OUT
#undef ARR
#undef DUP

/*----------------------------------------------------------------------------*/

void NumFFT_Swap1d(double *arr, size_t n)
{
	size_t i;
	double *dup;

	dup = Mem_CALLOC(n, dup);

	for(i = 0; i < n; i++) {
		if(i < n/2)
			dup[i] = arr[i+n/2];
		else
			dup[i] = arr[i-n/2];
	}

	for(i = 0; i < n; i++) {
		arr[i] = dup[i];
	}

	free(dup);

	return;
}

/*----------------------------------------------------------------------------*/

#define ARR(a,b) (arr[(b)+jdim*(a)])
#define DUP(a,b) (dup[(b)+jdim*(a)])

void NumFFT_Swap2d(double *arr, size_t idim, size_t jdim)
{
	size_t i, j;
	size_t icen = 0, jcen = 0;
	double *dup = NULL;

	dup = Mem_CALLOC(idim * jdim, dup);
	icen = idim/2;
	jcen = jdim/2;

	/* Copy and swap array to dup */
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim; j++) {
			#if 1
			if(i >= icen && j >= jcen) { //load Q1 of arr to Q3 of dup
				DUP(i-icen, j-jcen) = ARR(i,j);
			}
			else if(i < icen && j >= jcen) { //load Q2 of arr to Q4 of dup
				DUP(i+icen, j-jcen) = ARR(i,j);
			}
			else if(i < icen && j < jcen) { //load Q3 of arr to Q1 of dup
				DUP(i+icen, j+jcen) = ARR(i,j);
			}
			else if(i >= icen && j < jcen) { //load Q4 of arr to Q2 of dup
				DUP(i-icen, j+jcen) = ARR(i,j);
			}
			#endif

			#if 0 //debug
			if(i < idim/2) {
			if(j < jdim/2)
			  DUP(i,j) = ARR(i+idim/2,j+jdim/2); //i < idim/2 && j < jdim/2
			else
			  DUP(i,j) = ARR(i+idim/2,j-jdim/2); //i < idim/2 && j >= jdim/2
			}
			else {
			if(j < jdim/2)
			  DUP(i,j) = ARR(i-idim/2,j+jdim/2); //i >= idim/2 && j < jdim/2
			else
			  DUP(i,j) = ARR(i-idim/2,j-jdim/2); //i >= idim/2 && j >= jdim/2
			}
			#endif
		}
	}

	/* Copy swapped array back to input */
	for(i = 0; i < idim; i++) {
		for(j = 0; j < jdim; j++) {
			ARR(i,j) = DUP(i,j);
		}
	}

	/* Cleanup */
	free(dup);

	return;
}

#undef ARR
#undef DUP

/*----------------------------------------------------------------------------*/

double Num_GaussNormal(double x, double width)
/* Return Gauss normal at x, with width=sqrt(2)*sigma of the Gaussian.
 *
 * Since exp(x) *can overflow* for x = +\-(large number), indexing
 * a static array of precalculated gaussian normals not only saves time
 * but also prevents overflow of floating point numbers.
 *
 * MAXGAU = total number of data points -- must be odd numbered so that
 * the central data point is zero.
 * MAXWIDTH = maximum number of widths to trace Gaussian.
 * fac = number of data points per width.
 */
#define MAXGAU 801
#define MAXWIDTH 4
{
	static const double fac = (MAXGAU - 1.0) / MAXWIDTH;
	static int init = 0;
	static double gaunormal[MAXGAU];
	size_t i, igau;

	double xp, alpha, beta;

	/* If init==0, init gaunormal */
	if(!init) {
		for(i = 0; i < MAXGAU; i++) {
			gaunormal[i] = exp(-(double)(i * i) / (fac * fac));
		}
		init = 1;
	}
	if (!(width > 0.0)) printf("width=%g\n",width);
	/* Find position of x in Gaussian normal */
	Deb_ASSERT(width > 0.0); /* Just in case */
	//igau = (size_t)round(fac * fabs(x) / width);
        xp = fac * fabs(x) / width;
	igau = (size_t)xp;
	
	beta =  xp-(double)igau;
	alpha = 1.0-beta;
	/* If igau is >= MAXGAU, it is greater than 4*width, so return 0.
	 * Otherwise return Gaussian normal at igau. */
	if(igau >= MAXGAU-1)
	//if(xp >= MAXWIDTH)
		return 0.0;
	else
		//return gaunormal[igau];
		return alpha*gaunormal[igau] + beta*gaunormal[igau+1];
		//return;
}
#undef MAXGAU
#undef MAXWIDTH

/*----------------------------------------------------------------------------*/

double Num_GaussNormal2(double x, double width)
/* Return Gauss normal at x with width=2*sigma of Gaussian.
 * Note: exp() *may overflow* for large arguments!!! */
{
	return exp(-(x * x) / (width * width));
}

/*----------------------------------------------------------------------------*/

void Num_RanDir3D(gsl_rng *rng, double *cost, double *sint, double *phi)
/* Generate a random direction in 3D space */
{
	/* This samples a random number uniformly in the
	 * interval [0, 1) */
	#define RAND()\
		gsl_rng_uniform(rng)
	/* This samples a random number uniformly in the
	 * interval (0, 1) */
	#define PRAND()\
		gsl_rng_uniform_pos(rng)

	#if 1
	*cost = 2.0 * RAND() - 1.0;
	*sint = sqrt(1.0 - (*cost) * (*cost));
	*phi = PRAND() * Num_TWOPI;
	#else
		#if 0
		*cost = cos(0);
		*sint = sin(0);
		*phi = RAND() * Num_TWOPI;
		#endif

		#if 0
		*cost = 2.0 * RAND() - 1.0;
		*sint = sqrt(1.0 - (*cost) * (*cost));
		*phi = 0;
		#endif

		#if 0
		*cost = cos(0);
		*sint = sin(0);
		*phi = 0;
		#endif
	#endif

	#undef RAND
	#undef PRAND

	return;
}


/*----------------------------------------------------------------------------*/

void Num_QRanDir3D(const double *QRN, double *cost, double *sint, double *phi)
/* Generate a random direction in 3D space */
{
        *cost = 2.0 * QRN[3] - 1.0;
        *sint = sqrt(1.0 - (*cost) * (*cost));
        *phi = QRN[4] * Num_TWOPI;

        return;
}




























