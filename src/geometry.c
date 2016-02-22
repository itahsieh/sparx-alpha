#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <gsl/gsl_rng.h>

#include "memory.h"
#include "debug.h"

#if Sp_MIRSUPPORT
#include "cpgplot-wrappers.h"
#endif

#include "numerical.h"
#include "geometry.h"

DatINode GEOM_TYPES[] = {
	{"sph1d", GEOM_SPH1D},
	{"sph3d", GEOM_SPH3D},
	{"rec3d", GEOM_REC3D},
	{"cyl3d", GEOM_CYL3D},
	{0, 0}
};

DatINode *GEOM_TYPESP = GEOM_TYPES;
#if Sp_MIRSUPPORT
static void GeVox_Cpgplot_sph1d(const GeVox *voxel);
static void GeVox_Cpgplot_rec(const GeVox *voxel, const GeCam *cam);
#endif
/*----------------------------------------------------------------------------*/

size_t Ge_PosToIelem(size_t i, size_t j, size_t k, const GeVec3_s *naxes)
{
 	size_t ielem = 0;

	ielem = (k) + GeVec3_X(*naxes, 2) * ((j) + GeVec3_X(*naxes, 1) * (i));

	return ielem;
}

/*----------------------------------------------------------------------------*/

GeVec3_s Ge_IelemToIndex(size_t ielem, const GeVec3_s *naxes)
{
	size_t i;
	GeVec3_s idx;

	for(i = 0; i < 3; i++) {
		GeVec3_X(idx, 2 - i) = ielem % GeVec3_X(*naxes, 2 - i);
		ielem = ielem / GeVec3_X(*naxes, 2 - i);
	}

	return idx;
}

/*----------------------------------------------------------------------------*/

size_t Ge_IndexToIelem(const GeVec3_s *idx, const GeVec3_s *naxes)
{
	size_t i, ielem = 0;

	for(i = 0; i < 3; i++) {
		ielem = i == 0 ? GeVec3_X(*idx, i) : ielem * GeVec3_X(*naxes, i) + GeVec3_X(*idx, i);
	}

	return ielem;
}

/*----------------------------------------------------------------------------*/

GeVox GeVox_GetSubVox(const GeVox *voxel, const GeVec3_s *idx, const GeVec3_s *naxes)
{
	size_t i;
	double delta;
	GeVox subvox;

	for(i = 0; i < 3; i++) {
		/* Remember to update these if the GeVox struct changes!!! */
		subvox.geom = voxel->geom;
		delta = subvox.delta.x[i] = voxel->delta.x[i] / (double)naxes->x[i];
		subvox.min.x[i] = voxel->min.x[i] + (double)idx->x[i] * delta;
		subvox.max.x[i] = voxel->min.x[i] + (double)(idx->x[i] + 1) * delta;
		subvox.cen.x[i] = subvox.min.x[i] + 0.5 * delta; /* This is only valid if the gridding is uniform */
	}

	return subvox;
}

/*----------------------------------------------------------------------------*/

GeCam GeCam_Init(double phix, double phiy, double phiz)
/* Set plot bounds according to voxel and rotation for camera */
{
	GeCam cam = GeCam_INIT(phix, phiy, phiz);

	return cam;
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_Rotate(const GeVec3_d *vec, const GeCam *cam)
{
	GeVec3_d result;

	result = *vec;
	result = GeVec3_Rotate_x(&result, GeVec3_X(cam->phi, 0));
	result = GeVec3_Rotate_y(&result, GeVec3_X(cam->phi, 1));
	result = GeVec3_Rotate_z(&result, GeVec3_X(cam->phi, 2));

	return result;
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_d_Init(double x, double y, double z)
{
	GeVec3_d vec = GeVec3_INIT(x, y, z);

	return vec;
}

/*----------------------------------------------------------------------------*/

GeVec3_s GeVec3_s_Init(size_t x, size_t y, size_t z)
{
	GeVec3_s vec = GeVec3_INIT(x, y, z);

	return vec;
}

/*----------------------------------------------------------------------------*/

double GeVec3_Mag(const GeVec3_d *a)
/* Calculate the magnitude of <a> = |<a>| */
{
	size_t i;
	double mag = 0;

	for(i = 0; i < 3; i++) {
		mag += pow(GeVec3_X(*a, i), 2.0);
	}

	return sqrt(mag);
}

/*----------------------------------------------------------------------------*/

double GeVec3_Mag2(const GeVec3_d *a, const GeVec3_d *b)
/* Calculate the magnitude of <b> - <a> = |<b> - <a>| */
{
	GeVec3_d c;

	c = GeVec3_Sub(a, b);

	return GeVec3_Mag(&c);
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_Add(const GeVec3_d *a, const GeVec3_d *b)
/* Calculate the result of <c> = <a> + <b> */
{
	size_t i;
	GeVec3_d c = GeVec3_INIT(0,0,0);

	for(i = 0; i < 3; i++)
		c.x[i] = a->x[i] + b->x[i];
	
	return c;
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_Sub(const GeVec3_d *a, const GeVec3_d *b)
/* Calculate the result of <c> = <b> - <a> */
{
	size_t i;
	GeVec3_d c;

	for(i = 0; i < 3; i++)
		c.x[i] = b->x[i] - a->x[i];
	
	return c;
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_Normalize(const GeVec3_d *a)
/* Normalize <a> */
{
	GeVec3_d b;
	double mag = GeVec3_Mag(a);
	size_t i;

	for(i = 0; i < 3; i++) {
		GeVec3_X(b, i) = GeVec3_X(*a, i) / mag;
	}

	return b;
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_Scale(const GeVec3_d *a, double fac)
/* Scale <a> by fac, i.e. <b> = fac * <a> */
{
	size_t i;
	GeVec3_d b;

	Mem_BZERO(&b);

	for(i = 0; i < 3; i++) {
		GeVec3_X(b, i) = fac * GeVec3_X(*a, i);
	}

	return b;
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_InterpLinear(const GeVec3_d *xa, const GeVec3_d *xb, const GeVec3_d *a, 
	const GeVec3_d *b, const GeVec3_d *pos)
/* Interpolate linearly between <a> and <b>, i.e. infer value at pos on each
 * axis according to
 * 	result_i = a_i + m_i * (pos_i - xa_i)
 *
 * 	and
 *
 * 	m_i = (b_i - a_i) / (xb_i - xa_i) -- if |xb_i - xa_i| > 0
 * 	      = 0                           -- otherwise
 *
 * 	where origin_i = ith component of origin
 * 	      b_i and a_i = ith component of <a> and <b>
 * 	      xa_i and xb_i = ith component of positions of <a> and <b>
 * 	      pos_i = ith component of interpolation position
 */
{
	size_t i;
	double m, delta_x;
	//const double *xav = xa->x, *xbv = xb->x, *av = a->x, *bv = b->x, *posv = pos->x;
	GeVec3_d result;

	for(i = 0; i < 3; i++) {
		/* Calculate delta_x */
		delta_x = (GeVec3_X(*xb, i) - GeVec3_X(*xa, i));
		//delta_x = xbv[i] - xav[i];

		/* Calculate slope */
		m = fabs(delta_x) > 0 ? (GeVec3_X(*b, i) - GeVec3_X(*a, i)) / delta_x : 0;
		//m = fabs(delta_x) > 0 ? (bv[i] - av[i]) / delta_x : 0;

		/* Calculate linearly interpolated value */
		GeVec3_X(result, i) = GeVec3_X(*a, i) + m * (GeVec3_X(*pos, i) - GeVec3_X(*xa, i));
		//GeVec3_X(result, i) = av[i] + m * (posv[i] - xav[i]);
	}

	return result;
}

/*----------------------------------------------------------------------------*/

double GeVec3_DotProd(const GeVec3_d *a, const GeVec3_d *b)
/* Calculate prod = <a> dot <b> */
{
	size_t i;
	double prod = 0;

	for(i = 0; i < 3; i++)
		prod += a->x[i] * b->x[i];

	return prod;
}
/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_CrossProd(const GeVec3_d *a, const GeVec3_d *b)
/* Calculate prod = <a> cross <b> */
{
	GeVec3_d prod;

	GeVec3_X(prod,0) = a->x[1] * b->x[2] - a->x[2] * b->x[1];
	GeVec3_X(prod,1) = a->x[2] * b->x[0] - a->x[0] * b->x[2];
	GeVec3_X(prod,2) = a->x[0] * b->x[1] - a->x[1] * b->x[0];
	
	return prod;
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_MatOp(const GeVec3_d *vec, const GeMat3_d *mat)
/* Calculate <result> = [mat]<vec> */
{
	size_t i, j;
	GeVec3_d result = GeVec3_INIT(0, 0, 0);

	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			result.x[i] += GeMat3_X(*mat, i, j) * vec->x[j];
		}
	}

	return result;
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_Rotate_x(const GeVec3_d *vec, double phi)
/*
Rotate vector phi radians about x axis
Rotation matrix about x:
	1  	       0          0
	0 	cos(phi)  -sin(phi)
	0 	sin(phi)   cos(phi)
*/
{
	GeMat3_d matrix;

	GeMat3_X(matrix, 0, 0) = 1;
	GeMat3_X(matrix, 0, 1) = 0;
	GeMat3_X(matrix, 0, 2) = 0;

	GeMat3_X(matrix, 1, 0) = 0;
	GeMat3_X(matrix, 1, 1) = cos(phi);
	GeMat3_X(matrix, 1, 2) = -sin(phi);

	GeMat3_X(matrix, 2, 0) = 0;
	GeMat3_X(matrix, 2, 1) = sin(phi);
	GeMat3_X(matrix, 2, 2) = cos(phi);
	
	return GeVec3_MatOp(vec, &matrix);
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_Rotate_y(const GeVec3_d *vec, double phi)
/*
Rotate vector phi radians about y axis
Rotation matrix about y:
	cos(phi)  	0     sin(phi)
               0   	1 	     0
       -sin(phi) 	0     cos(phi)
*/
{
	GeMat3_d matrix;

	GeMat3_X(matrix, 0, 0) = cos(phi);
	GeMat3_X(matrix, 0, 1) = 0;
	GeMat3_X(matrix, 0, 2) = sin(phi);

	GeMat3_X(matrix, 1, 0) = 0;
	GeMat3_X(matrix, 1, 1) = 1;
	GeMat3_X(matrix, 1, 2) = 0;

	GeMat3_X(matrix, 2, 0) = -sin(phi);
	GeMat3_X(matrix, 2, 1) = 0;
	GeMat3_X(matrix, 2, 2) = cos(phi);
	
	return GeVec3_MatOp(vec, &matrix);
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeVec3_Rotate_z(const GeVec3_d *vec, double phi)
/*
Rotate vector phi radians about z axis
Rotation matrix about z:
	cos(phi)  -sin(phi)     0
        sin(phi)   cos(phi)     0
               0          0     1
*/
{
	GeMat3_d matrix;

	GeMat3_X(matrix, 0, 0) = cos(phi);
	GeMat3_X(matrix, 0, 1) = -sin(phi);
	GeMat3_X(matrix, 0, 2) = 0;

	GeMat3_X(matrix, 1, 0) = sin(phi);
	GeMat3_X(matrix, 1, 1) = cos(phi);
	GeMat3_X(matrix, 1, 2) = 0;

	GeMat3_X(matrix, 2, 0) = 0;
	GeMat3_X(matrix, 2, 1) = 0;
	GeMat3_X(matrix, 2, 2) = 1;
	
	return GeVec3_MatOp(vec, &matrix);
}

/*----------------------------------------------------------------------------*/

GeRay GeRay_Init(double ex, double ey, double ez, double dx, double dy, double dz)
{
	GeRay ray = GeRay_INIT(ex, ey, ez, dx, dy, dz);

	return ray;
}

/*----------------------------------------------------------------------------*/

double GeRay_IntersectPlane(const GeRay *ray, const GeVec3_d *n, const GeVec3_d *q)
/*
Calculate the intersection between ray and plane according to

       t = <n> dot (<q> - <e>) / <n> dot <d>

	where <n> is the normal vector,
	      <q> a vertex of the plane,
	      <e> the origin of the ray,
	  and <d> the direction of the ray
*/
{
	GeVec3_d q_m_e = GeVec3_Sub(&ray->e, q);

	return GeVec3_DotProd(n, &q_m_e) / GeVec3_DotProd(n, &ray->d);
}

/*----------------------------------------------------------------------------*/

void GeRay_IntersectSphere(const GeRay *ray, double r, double *t1, double *t2)
/*
Solve for the intersection between a ray
	<x> = <e> + <d> * t
and a sphere with radius R in Cartesian coordinates (where <e> is OUTSIDE
the sphere). The 3D equatin of a sphere is
	|<x>|^2 = R^2
resulting in
	|<e>|^2 + 2 * <d> dot <e> * t + |<d>|^2 * t^2 = R^2
and the intersection between a line and a sphere can be found by solving the
quadratic equaiton
	a * t^2 + b * t + c = 0
where
	a = |<d>|^2
	b = 2 * <e> dot <d>
	c = |<e>|^2 - R^2
*/
{
	double a, b, c;

	a = GeVec3_DotProd(&ray->d, &ray->d); /* Is <d> guaranteed to be unit vector? (questioned by I-Ta) */
	b = 2.0 * GeVec3_DotProd(&ray->e, &ray->d);
	c = GeVec3_DotProd(&ray->e, &ray->e) - r * r;

	Num_QuadraticRoots(a, b, c, t1, t2);

	#if 0
	Deb_PRINT("solve: r=%g, t1=%g, t2=%g\n", r, *t1, *t2);
	Deb_PAUSE();
	#endif

	if(Num_ISNAN(*t1))
		*t1 = HUGE_VAL;

	if(Num_ISNAN(*t2))
		*t2 = HUGE_VAL;

	return;
}

/*----------------------------------------------------------------------------*/
void GeRay_IntersectSphere3d(const GeRay *ray, double r, double *t1, double *t2)
/*
Solve for the intersection between a ray
	<x> = <e> + <d> * t
and a sphere with radius R in Cartesian coordinates (where <e> is OUTSIDE
the sphere). The 3D equatin of a sphere is
	|<x>|^2 = R^2
resulting in
	|<e>|^2 + 2 * <d> dot <e> * t + |<d>|^2 * t^2 = R^2
and the intersection between a line and a sphere can be found by solving the
quadratic equaiton
	a * t^2 + b * t + c = 0
where
	a = |<d>|^2
	b = 2 * <e> dot <d>
	c = |<e>|^2 - R^2
*/
{
	double a, b, c;

	a = GeVec3_DotProd(&ray->d, &ray->d); /* Is <d> guaranteed to be unit vector? (questioned by I-Ta) */
	b = 2.0 * GeVec3_DotProd(&ray->e, &ray->d);
	c = GeVec3_DotProd(&ray->e, &ray->e) - r * r;

	Num_QuadraticRoots(a, b, c, t1, t2);

	#if 0
	Deb_PRINT("solve: r=%g, t1=%g, t2=%g\n", r, *t1, *t2);
	Deb_PAUSE();
	#endif
	if(Num_ISNAN(*t1))
		*t1 = HUGE_VAL;
	
	if(Num_ISNAN(*t2))
		*t2 = HUGE_VAL;

	return;
}

/*----------------------------------------------------------------------------*/
void GeRay_IntersectTheta(const GeRay *ray, double theta, double *t1, double *t2)
/*
Solve for the intersection between a ray
	<x> = <e> + <d> * t
and a cone with the angle Theta in Cartesian coordinates (where <e> is OUTSIDE
the cone). The 3D equatin of a cone is
	cos(theta) = <x>_z / |<x>|
resulting in
	cos^2(theta)*[ |<e>|^2 + * <d> dot <e> * t + |<d>|^2 * t^2 ]
= e_z + d_z * t
&&
	cos(theta)*(e_z+d_z*t)>0
and the intersection between a line and a cone can be found by solving the
quadratic equaiton
	a * t^2 + b * t + c = 0
where
	a = cos^2(theta)*|<d>|^2
	b = cos^2(theta)* 2 * <e> dot <d> - d_z
	c = cos^2(theta)*|<e>|^2 - e_z
*/
{
	double a, b, c, cos2theta, costheta;
	
	costheta = cos(theta);
	cos2theta = costheta * costheta;

	a = cos2theta - GeRay_D(*ray, 2) * GeRay_D(*ray, 2); 
	b = cos2theta * 2.0 * GeVec3_DotProd(&ray->e,&ray->d) - 2.0 * GeRay_D(*ray, 2) * GeRay_E(*ray, 2);
	c = cos2theta * GeVec3_DotProd(&ray->e, &ray->e) - GeRay_E(*ray, 2) * GeRay_E(*ray, 2);

	

	
	if( cos2theta < 1e-24){
		*t1 = -0.5*b/a;
		*t2 = *t1;
	}
	else{
		Num_QuadraticRoots(a, b, c, t1, t2);
		if(Num_ISNAN(*t1) || Num_ISNAN(*t2)){
			*t1 = HUGE_VAL;
			*t2 = HUGE_VAL;
		}
		/* check cos(theta)*(e_z+d_z*t)>0 */
		else{
			if( costheta * (GeRay_E(*ray, 2) + GeRay_D(*ray, 2) * (*t1)) < 0 )
				*t1 = HUGE_VAL;
			if( costheta * (GeRay_E(*ray, 2) + GeRay_D(*ray, 2) * (*t2)) < 0 )
				*t2 = HUGE_VAL;
		}
	}
	
	#if 0
	Deb_PRINT("solve: theta=%g, t1=%g, t2=%g\n", theta, *t1, *t2);
	Deb_PAUSE();
	#endif
	
	return;
}

/*----------------------------------------------------------------------------*/
void GeRay_IntersectPhi(const GeRay *ray, double phi, double *t)
/*
Solve for the intersection between a ray
	<x> = <e> + <d> * t
and a cone with the angle Theta in Cartesian coordinates (where <e> is OUTSIDE
the cone). The 3D equatin of a cone is
	tan(phi) = <x>_y / <x>_x
resulting in
	tan(phi) = (<e>_y+<d>_y*t) / (<e>_x+<d>_x*t)
&&
	<e>_y+<d>_y*t >= 0, 0=<phi=<pi
	<e>_y+<d>_y*t <  0, pi<phi<2pi
and the intersection between a line and a cone can be found by solving the
linear equaiton
	t = (tan(phi)*<e>_x-<e>_y)/(<d>_y-tan(phi)*<d>_x)
*/
{
	double tanphi, a, b;
        static const double pi = 4.* atan(1.);
	
	tanphi=tan(phi);
	a=( tanphi*GeRay_E(*ray, 0)-GeRay_E(*ray, 1) );
	b=( GeRay_D(*ray, 1)-tanphi*GeRay_D(*ray, 0) );
	if(fabs(b)<1e-8)
		*t=HUGE_VAL;
	else
		*t = a/b;
	
	/* check consistency of phi and <x>_y */
	if( *t <= 0.0 )
		*t = HUGE_VAL;
	else{ 
		double RayPos_y = GeRay_E(*ray, 1) + GeRay_D(*ray, 1) * (*t);
		if( phi <= pi ){
			if ( RayPos_y < 0. )
				/* phi <= phi but <x>_y < 0 => inconsistency */
				*t = HUGE_VAL;
		}
		else{
			if ( RayPos_y >= 0. ) {
				/* phi > phi but <x>_y >= 0 => inconsistency */
				*t = HUGE_VAL;
			}
		}
	}


	return;
}
/*----------------------------------------------------------------------------*/
void GeRay_IntersectRc(const GeRay *ray, double Rc, double *t1, double *t2)
/*
Solve for the intersection between a ray
	<x> = <e> + <d> * t
and a cylinder with radius Rc in Cartesian coordinates (where <e> is OUTSIDE
the sphere). The 3D equatin of a sphere is
	<x>_x^2 + <x>_y^2 = Rc^2
resulting in
	<e>_x^2 + 2 * <d>_x * <e>_x * t + <d>_x^2 * t^2 
      + <e>_y^2 + 2 * <d>_y * <e>_y * t + <d>_y^2 * t^2 = Rc^2
and the intersection between a line and a sphere can be found by solving the
quadratic equaiton
	a * t^2 + b * t + c = 0
where
	a = <d>_x^2 + <d>_y^2
	b = 2 * <e>_x * <d>_x + 2 * <e>_y * <d>_y
	c = <e>_x^2 + <e>_y^2 - R_c^2
*/
{
	double a, b, c;

	a = ray->d.x[0] * ray->d.x[0] + ray->d.x[1] * ray->d.x[1]; 
	b = 2.0 * ray->d.x[0] * ray->e.x[0] + 2.0 * ray->d.x[1] * ray->e.x[1];
	c = ray->e.x[0] * ray->e.x[0] + ray->e.x[1] * ray->e.x[1] - Rc * Rc;

	Num_QuadraticRoots(a, b, c, t1, t2);

	#if 0
	Deb_PRINT("solve: r=%g, t1=%g, t2=%g\n", r, *t1, *t2);
	Deb_PAUSE();
	#endif
	if(Num_ISNAN(*t1))
		*t1 = HUGE_VAL;
	
	if(Num_ISNAN(*t2))
		*t2 = HUGE_VAL;

	return;
}
/*----------------------------------------------------------------------------*/
void GeRay_IntersectHz(const GeRay *ray, double Hz, double *t)
/*
Calculate the intersection between ray and plane according to

       t = (H_z - e_z) / d_z

	where H_z is the vertical height,
	      e_z is the z-component of the origin of the ray,
	  and d_z is the z-component of the direction of the ray
*/
{
	*t = ( Hz - ray->e.x[2]) / ray->d.x[2];
	return ;
}
/*----------------------------------------------------------------------------*/
GeRay GeRay_Inc(const GeRay *ray, double t)
{
	size_t i;
	GeRay newray = *ray;

	for(i = 0; i < 3; i++) {
		newray.e.x[i] += newray.d.x[i] * t;
	}

	return newray;
}

/*----------------------------------------------------------------------------*/

GeRay GeRay_Rotate(const GeRay *ray, int axis, double phi)
/* Rotate ray about x axis */
{
	GeRay ray2 = GeRay_INIT(0, 0, 0, 0, 0, 0);

	switch(axis) {
		case 0:
			ray2.e = GeVec3_Rotate_x(&ray->e, phi);
			ray2.d = GeVec3_Rotate_x(&ray->d, phi);
			break;

		case 1:
			ray2.e = GeVec3_Rotate_y(&ray->e, phi);
			ray2.d = GeVec3_Rotate_y(&ray->d, phi);
			break;

		case 2:
			ray2.e = GeVec3_Rotate_z(&ray->e, phi);
			ray2.d = GeVec3_Rotate_z(&ray->d, phi);
			break;

		default:
			/* Not a valid axis */
			Deb_ASSERT(0);
	}

	return ray2;
}

/*----------------------------------------------------------------------------*/

int GeRay_IntersectVoxel(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side)
/* Find the intersection of a ray OUTSIDE voxel in rectangular coordinates. */
{
	int hit = 0;
	switch(voxel->geom) {
		case GEOM_SPH1D:
			hit = GeRay_IntersectVoxel_sph1d(ray, voxel, tmin, side);
			break;
			
		case GEOM_SPH3D:
			hit = GeRay_IntersectVoxel_sph3d(ray, voxel, tmin, side);
			break;

		case GEOM_REC3D:
			hit = GeRay_IntersectVoxel_rec3d(ray, voxel, tmin, side);
			break;
			
		case GEOM_CYL3D:
			hit = GeRay_IntersectVoxel_cyl3d(ray, voxel, tmin, side);
			break;

		default: /* Shouldn't happen */
			Deb_ASSERT(0);
	}

	return hit;
}

/*----------------------------------------------------------------------------*/

int GeRay_IntersectVoxel_sph1d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side)
/*
Solve for intersection between ray and voxel in which ray is OUTSIDE of voxel.

A 1D spherical voxel is described by two concentric spheres, thus for a ray
to intersect the voxel, it merely needs to hit the outer sphere.

side=0 --> inner sphere
side=1 --> outer sphere

If there is a hit, side would always be 1 for this routine.
*/
{
	double r = GeVec3_X(voxel->max, 0), t1, t2;

	GeRay_IntersectSphere(ray, r, &t1, &t2);
	
	if(t1 < t2)
		*tmin = t1;
	else
		*tmin = t2;

	*side = 1;

	return (*tmin < HUGE_VAL ? 1 : 0);
}

/*----------------------------------------------------------------------------*/

int GeRay_IntersectVoxel_sph3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side)
/* Find the intersection of a ray OUTSIDE voxel in rectangular coordinates. */
{

	double r = GeVec3_X(voxel->max, 0);
	double t1, t2;

	GeRay_IntersectSphere3d(ray, r, &t1, &t2);
	
	if( t1 < t2 )
		*tmin = t1;
	else
		*tmin = t2;
	

	*side = 1;

	return (*tmin < HUGE_VAL ? 1 : 0);
}

/*----------------------------------------------------------------------------*/

int GeRay_IntersectVoxel_rec3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side)
/* Find the intersection of a ray OUTSIDE voxel in rectangular coordinates. */
{
	double t;
	size_t i;
	int within_box;
	static GeVec3_d normal[6] = {
		GeVec3_INIT(-1,  0,  0),
		GeVec3_INIT( 1,  0,  0),
		GeVec3_INIT( 0, -1,  0),
		GeVec3_INIT( 0,  1,  0),
		GeVec3_INIT( 0,  0, -1),
		GeVec3_INIT( 0,  0,  1)
	};
	const GeVec3_d *q = NULL;
	GeRay tstray = GeRay_INIT(0, 0, 0, 0, 0, 0);

	/* Init tmin */
	*tmin = HUGE_VAL;

	/* debug */
	for(i = 0; i < 6; i++) {
		/* Calculate intersections with side */
		q = i % 2 == 0 ? &voxel->min : &voxel->max;
		t = GeRay_IntersectPlane(ray, &normal[i], q);

		/* Check if intersection is inside box */
		tstray = GeRay_Inc(ray, t);
		within_box = point_in_voxel2(&tstray.e, voxel, i / 2);

		/* Find tmin if intersection is within box*/
		if(within_box && (t < *tmin)) {
			*tmin = t;
			*side = i;
		}
	}
	
	return (*tmin < HUGE_VAL ? 1 : 0);
}

/*----------------------------------------------------------------------------*/
int GeRay_IntersectVoxel_cyl3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side)
/* Find the intersection of a ray OUTSIDE voxel in rectangular coordinates. */
{
	double t[6];
	size_t i;
	int within_box;
	GeRay tstray = GeRay_INIT(0, 0, 0, 0, 0, 0);

	/* Init tmin */
	*tmin = HUGE_VAL;

	/* debug */

	/* Calculate intersections with side */
	GeRay_IntersectRc(ray, voxel->max.x[0], &t[0], &t[1]);
	GeRay_IntersectPhi(ray, voxel->min.x[1], &t[2]);
	GeRay_IntersectPhi(ray, voxel->max.x[1], &t[3]);
	GeRay_IntersectHz(ray, voxel->min.x[2], &t[4]);
	GeRay_IntersectHz(ray, voxel->max.x[2], &t[5]);

	for(i = 0; i < 6; i++) {
		/* Check if intersection is inside box */
		tstray = GeRay_Inc(ray, t[i]);
		within_box = point_in_voxel_cyl3d(&tstray.e, voxel, i / 2);

		/* Find tmin if intersection is within box*/
		if(within_box && (t[i] < *tmin)) {
			*tmin = t[i];
			*side = (i == 0) ? 1:i;
		}
	}
	
	return (*tmin < HUGE_VAL ? 1 : 0);
}

/*----------------------------------------------------------------------------*/
int point_in_voxel2(const GeVec3_d *pt, const GeVox *voxel, size_t axis)
/* Check if coordinates of pt NOT on axis are within the limits of the voxel */
{
	size_t i;
	int within_box = 1;

	for(i = 0; i < 3; i++) {
		if(i != axis) {
			if((GeVec3_X(*pt, i) < GeVec3_X(voxel->min, i)) || (GeVec3_X(*pt, i) > GeVec3_X(voxel->max, i))) {
				within_box = 0;
				break;
			}
		}
	}

	return within_box;
}

/*----------------------------------------------------------------------------*/
int point_in_voxel_cyl3d(const GeVec3_d *pt, const GeVox *voxel, size_t axis)
/* Check if coordinates of pt NOT on axis are within the limits of the voxel */
{
	int within_box = 1;
	
	double Rc, phi, Hz;
	static double pi = 4. * atan(1.); 
	Rc = sqrt( GeVec3_X(*pt, 0) * GeVec3_X(*pt, 0) + GeVec3_X(*pt, 1) * GeVec3_X(*pt, 1) );
	if( GeVec3_X(*pt, 1) >= 0. ) phi = acos(GeVec3_X(*pt, 0)/Rc);
	else phi = 2. * pi - acos(GeVec3_X(*pt, 0)/Rc);
	Hz = GeVec3_X(*pt, 2); 
	
	switch(axis) {
		case 0:
			/* Pt is the intersection on cylinder with the radius Rc 
			   It should be guaranteed that 
			   Phi_min <= phi <= Phi_max &&
			   Hz_min <= Hz <= Hz_max  */
			
			if(
				phi < GeVec3_X(voxel->min, 1) ||
				phi > GeVec3_X(voxel->max, 1) ||
				Hz < GeVec3_X(voxel->min, 2) ||
				Hz > GeVec3_X(voxel->max, 2)
			){
				within_box = 0;
				break;
			}
		case 1:
			/* Pt is the intersection on horizontal plane on z=Hz 
			   It should be guaranteed that 
			   Rc_min <= Rc <= Rc_max &&
			   Hz_min <= Hz <= Hz_max*/
			if( 
				Rc < GeVec3_X(voxel->min, 0) ||
				Rc > GeVec3_X(voxel->max, 0) ||
				Hz < GeVec3_X(voxel->min, 2) ||
				Hz > GeVec3_X(voxel->max, 2)
			){
				within_box = 0;
				break;
			}
		case 2:
			/* Pt is the intersection on horizontal plane on z=Hz 
			   It should be guaranteed that 
			   Rc_min <= Rc <= Rc_max && 
			   Phi_min <= phi <= Phi_max  */
			if( 
				Rc < GeVec3_X(voxel->min, 0) ||
				Rc > GeVec3_X(voxel->max, 0) ||
				phi < GeVec3_X(voxel->min, 1) ||
				phi > GeVec3_X(voxel->max, 1)
				
			){
				within_box = 0;
				break;
			}
	}
	return within_box;
}

/*----------------------------------------------------------------------------*/

void GeRay_TraverseVoxel(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side)
{
	switch(voxel->geom) {
		case GEOM_SPH1D:
			GeRay_TraverseVoxel_sph1d(ray, voxel, tmin, side);
			break;
        
                case GEOM_SPH3D:
			GeRay_TraverseVoxel_sph3d(ray, voxel, tmin, side);
			break;
        
		case GEOM_REC3D:
			GeRay_TraverseVoxel_rec3d(ray, voxel, tmin, side);
			break;

		case GEOM_CYL3D:
			GeRay_TraverseVoxel_cyl3d(ray, voxel, tmin, side);
			break;
			
		default: /* Shouldn't happen */
			Deb_ASSERT(0);
	}

	return;
}

/*----------------------------------------------------------------------------*/

void GeRay_TraverseVoxel_sph1d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side)
/*
Traverse a spherical shell voxel by calculating intersections with its inner
and outer spheres. Since the ray must originate from INSIDE the sphere1d voxel,
the larger distance of the outer sphere should be used, while the smaller of the
innter sphere should be used.

Side 0: inner sphere
Side 1: outer sphere
*/
{
	double t1, t2, t_in = HUGE_VAL, t_out = HUGE_VAL,
	       r_in = GeVec3_X(voxel->min, 0), r_out = GeVec3_X(voxel->max, 0);

	/* Find intersections with outer sphere */
	GeRay_IntersectSphere(ray, r_out, &t1, &t2);

	//Deb_PRINT("outer: t1=%g, t2=%g\n", t1, t2);

	t_out = (t1 > t2 ? t1 : t2);

	if(t_out <= 0)
		t_out = HUGE_VAL;

	/* Find intersection with inner sphere if r_in > 0 */
	if(r_in > 0) {
		GeRay_IntersectSphere(ray, r_in, &t1, &t2);

		//Deb_PRINT("inner: t1=%g, t2=%g\n", t1, t2);

		t_in = (t1 < t2 ? t1 : t2);

		if(t_in <= 0)
			t_in = HUGE_VAL;
	}

	Deb_ASSERT((t_in < HUGE_VAL) || (t_out < HUGE_VAL));

	/* Find final intersection */
	if(t_in < t_out) {
		*tmin = t_in;
		*side = 0;
	}
	else {
		*tmin = t_out;
		*side = 1;
	}

	return;
}

/*----------------------------------------------------------------------------*/
void GeRay_TraverseVoxel_sph3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side)
/*
Traverse a spherical shell voxel by calculating intersections with its inner
and outer spheres. Since the ray must originate from INSIDE the sphere1d voxel,
the larger distance of the outer sphere should be used, while the smaller of the
innter sphere should be used.

Side 0&1: inner/outer sphere
Side 2&3: lower/upper theta
Side 4&5: lower/upper phi
*/
{
	double t0, t1, t2, t3, t4, t5, 
	       r_in = GeVec3_X(voxel->min, 0), r_out = GeVec3_X(voxel->max, 0),
	       theta_in = GeVec3_X(voxel->min, 1), theta_out = GeVec3_X(voxel->max, 1),
	       phi_in = GeVec3_X(voxel->min, 2), phi_out = GeVec3_X(voxel->max, 2);
	size_t oldside=*side;
	static const double half_pi = 0.5*3.1415926535897932384626433832795;
	/* Find intersections with outer sphere */
		
	
	if( oldside == 0 )
		t0=HUGE_VAL;
	else{
		double t00, t01;
                GeRay_IntersectSphere3d(ray, r_in, &t00, &t01);
		t0 = Num_MIN(t00,t01);
		//if( t0 <= 0.0 ) t0=HUGE_VAL;
		if( t0 < 0.0 ) t0=HUGE_VAL;
	}
	
	double t10, t11;
	GeRay_IntersectSphere3d(ray, r_out, &t10, &t11);
	if( oldside == 1 )
		t1 = Num_MAX(t10,t11);
	else{
		if(t10<=0.0) t10=HUGE_VAL;
		if(t11<=0.0) t11=HUGE_VAL;
		t1 = Num_MIN(t10,t11);
	}
	
	double t20, t21;
	if( theta_in < half_pi || fabs(theta_in-half_pi)<1e-16 ){ 
		//printf("upper semi-sphere %g %g\n",fabs(half_pi-theta_out),fabs(half_pi-theta_in));
		if( oldside == 2 )
			t2 = HUGE_VAL;
		else{
			GeRay_IntersectTheta(ray, theta_in, &t20, &t21);
			t2 = Num_MIN(t20,t21);
			if(t2<=0.0) t2=HUGE_VAL;
		}
	}
	else{
		GeRay_IntersectTheta(ray, theta_in, &t20, &t21);
		if( oldside == 2 ){
			t2 = Num_MAX(t20,t21);
			if(t2<=0.0) t2=HUGE_VAL;
		}
		else{
			if(t20<=0.0) t20=HUGE_VAL;
			if(t21<=0.0) t21=HUGE_VAL;
			t2 = Num_MIN(t20,t21);
		}
	}
	
	double t30, t31;
	if(theta_out < half_pi && fabs(theta_in-half_pi)>1e-16 ){
                
		GeRay_IntersectTheta(ray, theta_out, &t30, &t31);
		if( oldside == 3 ){
			t3 = Num_MAX(t30,t31);
			if(t3<=0.0) t3=HUGE_VAL;
		}
		else{
			if(t30<=0.0) t30=HUGE_VAL;
			if(t31<=0.0) t31=HUGE_VAL;
			t3 = Num_MIN(t30,t31);
		}
	}
	else{
		if( oldside == 3 )
			t3 = HUGE_VAL;
		else{
			GeRay_IntersectTheta(ray, theta_out, &t30, &t31);
			t3 = Num_MIN(t30,t31);
			if(t3<=0.0) t3=HUGE_VAL;
		}
	}
	
	if( oldside == 4 )
		t4=HUGE_VAL;
	else if( sin(phi_in)*GeRay_D(*ray, 0) - cos(phi_in)*GeRay_D(*ray, 1) < 0.0 )
		t4=HUGE_VAL;
	else
		GeRay_IntersectPhi(ray, phi_in, &t4);

	if( oldside == 5 )
		t5 = HUGE_VAL;
	else if( -sin(phi_out)*GeRay_D(*ray, 0) + cos(phi_out)*GeRay_D(*ray, 1) < 0.0 )
		t5=HUGE_VAL;
	else
		GeRay_IntersectPhi(ray, phi_out, &t5);
	

	/* Init tmin */
	*tmin = HUGE_VAL;
	/* Find minimal intersection  */
	if( t0 < *tmin ){
		*tmin=t0;
		*side=0;
	}
	if( t1 < *tmin ){
		*tmin=t1;
		*side=1;
	}
	if( t2 < *tmin ){
		*tmin=t2;
		*side=2;
	}
	if( t3 < *tmin ){
		*tmin=t3;
		*side=3;
	}
	if( t4 < *tmin ){
		*tmin=t4;
		*side=4;
	}
	if( t5 < *tmin ){
		*tmin=t5;
		*side=5;
	}
	
	#if 0
	Deb_PRINT("solve: r=%g\n", r_out );
	Deb_PRINT("solve: t0=%g, t1=%g, t2=%g, t3=%g, t4=%g, t5=%g, tmin=%g\n",  t0, t1, t2, t3, t4, t5, *tmin);
	Deb_PRINT("solve: oldside=%d side=%d\n", oldside,*side);
	Deb_PRINT("solve: theta_in=%g t20=%g t21=%g\n", theta_in,t20,t21);
	Deb_PRINT("solve: theta_out=%g t30=%g t31=%g\n", theta_out,t30,t31);
	Deb_PRINT("solve: Ex=%g Ey=%g Ez=%g\n", GeRay_E(*ray, 0),GeRay_E(*ray, 1),GeRay_E(*ray, 2));
	Deb_PRINT("solve: Dx=%g Dy=%g Dz=%g\n", GeRay_D(*ray, 0),GeRay_D(*ray, 1),GeRay_D(*ray, 2));
	
	R=sqrt(GeRay_E(*ray, 0)*GeRay_E(*ray, 0)+GeRay_E(*ray, 1)*GeRay_E(*ray, 1)+GeRay_E(*ray, 2)*GeRay_E(*ray, 2));
	theta=acos( GeRay_E(*ray, 2) / R );
	
	if(GeRay_E(*ray, 1) >= 0.0)
		phi = acos(GeRay_E(*ray, 0)/sqrt(GeRay_E(*ray, 0)*GeRay_E(*ray, 0)+GeRay_E(*ray, 1)*GeRay_E(*ray, 1)));
	else 
		phi = 2.*pi-acos(GeRay_E(*ray, 0)/sqrt(GeRay_E(*ray, 0)*GeRay_E(*ray, 0)+GeRay_E(*ray, 1)*GeRay_E(*ray, 1)));
	
	if( oldside!=0 ){
		if(R<r_in){
			Deb_PRINT("solve: R_in is out of range. R=%g R_in=%g \n",R,r_in);
			Deb_ASSERT(0);
		}
	}
	else{
		if(fabs(R-r_in)>=1e-8){
			Deb_PRINT("solve: R_in is not on the side. R=%g R_in=%g\n",R,r_in);
			Deb_ASSERT(0);
		}
	}
	
	if( oldside!=1){
		if(R>r_out){
			Deb_PRINT("solve: R_out is not compatible. R=%g R_out=%g\n",R,r_out);
			Deb_ASSERT(0);			
		}
	}
	else{
		if(fabs(R-r_out)>=1e-8){
			Deb_PRINT("solve: R_out is not on the side. R=%g R_out=%g\n",R,r_out);
			Deb_ASSERT(0);
		}
	}

	if( oldside!=2){
		if(theta<theta_in){
			Deb_PRINT("solve: theta_in is not compatible. theta=%g theta_in=%g\n",theta,theta_in);
			Deb_ASSERT(0);			
		}
	}
	else{
		if(fabs(theta-theta_in)>=1e-8){
			Deb_PRINT("solve: theta_in is not on the side. theta=%g theta_in=%g\n",theta,theta_in);
			Deb_ASSERT(0);
		}
	}
		
	if( oldside!=3){
		if(theta>theta_out){
			Deb_PRINT("solve: theta_out is not compatible. theta=%g theta_out=%g\n",theta,theta_out);
			Deb_ASSERT(0);			
		}
	}
	else{
		if(fabs(theta-theta_out)>=1e-6){
			Deb_PRINT("solve: theta_out is not on the side. theta=%g theta_out=%g\n",theta,theta_out);
			Deb_ASSERT(0);
		}
	}
	
	if( oldside!=4){
		if(phi<phi_in){
			Deb_PRINT("solve: phi_in is not compatible. phi=%g phi_in=%g\n",phi,phi_in);
			Deb_ASSERT(0);			
		}
	}
	else{
		if( fabs(phi-phi_in)>=1e-8 && fabs(phi-phi_out)>=1e-8 ){
			Deb_PRINT("solve: phi_in is not on the side. phi=%g phi_in=%g\n",phi,phi_in);
			Deb_ASSERT(0);
		}
	}
	
	if( oldside!=5){
		if(phi>phi_out){
			Deb_PRINT("solve: phi_out is not compatible. phi=%g phi_out=%g\n",phi,phi_out);
			Deb_ASSERT(0);			
		}
	}
	else{
		if( fabs(phi-phi_out)>=1e-8 && fabs(phi-phi_in)>=1e-8 ){
			Deb_PRINT("solve: phi_out is not on the side. phi=%g phi_out=%g\n",phi,phi_out);
			Deb_ASSERT(0);
		}
	}
	Deb_PRINT("solve: \n");
	//Deb_PAUSE();
	#endif
	Deb_ASSERT( *tmin < HUGE_VAL);



	return;
}

/*----------------------------------------------------------------------------*/
void GeRay_TraverseVoxel_rec3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side)
/*
Traverse a rectangular voxel using a variant of the Amanatides-Woo algorithm. An
intersection is ALWAYS guranteed, since we're traveling out of a bounded volume.
	* sides 0 & 1: x axis
	* sides 2 & 3: y axis
	* sides 4 & 5: z axis
*/
{
	static const GeVec3_d normal[6] = {
		GeVec3_INIT(-1,  0,  0),
		GeVec3_INIT( 1,  0,  0),
		GeVec3_INIT( 0, -1,  0),
		GeVec3_INIT( 0,  1,  0),
		GeVec3_INIT( 0,  0, -1),
		GeVec3_INIT( 0,  0,  1)
	};
	const GeVec3_d *q = NULL;

	/* Init tmin */
	*tmin = HUGE_VAL;

	for(size_t i = 0; i < 3; i++) {
		/*
		Calculate intersections with planes on all 3 axes
	          i == 0: x axis
		  i == 1: y axis
		  i == 2: z axis
		*/
		if(GeRay_D(*ray, i) == 0.0) {
			/* GeRay does not intersect either plane on this axis */
			continue;
		}
		else {
                        size_t iside;
			if(GeRay_D(*ray, i) > 0) {
				iside = i * 2 + 1;
				q = &voxel->max;
			}
			else {
				iside = i * 2;
				q = &voxel->min;
			}
			double t = GeRay_IntersectPlane(ray, &normal[iside], q);

			/* Find tmin */
			if(t < *tmin) {
				*tmin = t;
				*side = iside;
			}
		}

		#if 0 //debug
		debug_printf("D="); GeVec3_PRINT(stdout, ray->d); printf("\n");
		debug_printf("iside=%d, t=%g\n", iside, t);
		#endif
	}

	return;
}

/*----------------------------------------------------------------------------*/
void GeRay_TraverseVoxel_cyl3d(const GeRay *ray, const GeVox *voxel, double *tmin, size_t *side)
/*
Traverse a cylindrical coordinate voxel by calculating intersections with its corresponding cylinder and vertical height. Since the ray must originate from INSIDE the CYL3D voxel,
the larger distance of the outer cylinder should be used, while the smaller of the
innter cylinder should be used.

Side 0&1: inner/outer cylinder
Side 2&3: minimum/maximum phi
Side 4&5: lower/upper Height
*/
{
	double t[4],
	       Rc_in = GeVec3_X(voxel->min, 0), 
	       Rc_out = GeVec3_X(voxel->max, 0),
	       Phi_min = GeVec3_X(voxel->min, 1),
	       Phi_max = GeVec3_X(voxel->max, 1),
	       Hz_l = GeVec3_X(voxel->min, 2), 
	       Hz_u = GeVec3_X(voxel->max, 2);
	size_t oldside=*side;
		
	
	/* if the ray is shooting from inside Rc */
	if( oldside == 0 )
		t[0]=HUGE_VAL;
	else{
                double t00, t01;
		GeRay_IntersectRc(ray, Rc_in, &t00, &t01);
		t[0] = Num_MIN(t00,t01);
		if( t[0] < 0.0 ) t[0]=HUGE_VAL; /* just in case */
	}
	
	/* if the ray is shooting from outside Rc */
	double t10, t11;
        GeRay_IntersectRc(ray, Rc_out, &t10, &t11);
	if( oldside == 1 )
		t[1] = Num_MAX(t10,t11);
	else{
		if(t10<=0.0) t10=HUGE_VAL;
		if(t11<=0.0) t11=HUGE_VAL;
		t[1] = Num_MIN(t10,t11);
	}
	
	/* if the ray is shooting from the side of phi_min  */
	if( oldside == 2 )
		t[2] = HUGE_VAL;
	else{
		GeRay_IntersectPhi(ray, Phi_min, &t[2]);
		if( t[2] < 0.0 ) t[2] = HUGE_VAL;
	}

	/* if the ray is shooting from the side of phi_max  */
	if( oldside == 3 )
		t[3] = HUGE_VAL;
	else{
		GeRay_IntersectPhi(ray, Phi_max, &t[3]);
		if( t[3] < 0.0 ) t[3] = HUGE_VAL;
	}
	
	/* if the ray is shooting from the side of Hz_l  */
	if( oldside == 4 )
		t[4] = HUGE_VAL;
	else{
		GeRay_IntersectHz(ray, Hz_l, &t[4]);
		if( t[4] < 0.0 ) t[4] = HUGE_VAL;
	}

	/* if the ray is shooting from the side of Hz_u  */
	if( oldside == 5 )
		t[5] = HUGE_VAL;
	else{
		GeRay_IntersectHz(ray, Hz_u, &t[5]);
		if( t[5] < 0.0 ) t[5] = HUGE_VAL;
	}
	

	/* Init tmin */
	*tmin = HUGE_VAL;
	/* Find minimal intersection  */
	for (int i = 0; i < 6; i++){
		if( t[i] < *tmin ){
			*tmin=t[i];
			*side=i;
		}
	}
	
	Deb_ASSERT( *tmin < HUGE_VAL);



	return;
}

/*----------------------------------------------------------------------------*/
/* This samples a random number uniformly in the
 * interval [0, 1) */
#define RAND()\
	gsl_rng_uniform(rng)
/* This samples a random number uniformly in the
 * interval (0, 1) */
#define PRAND()\
	gsl_rng_uniform_pos(rng)

GeRay GeRay_Rand(gsl_rng *rng, const GeVox *voxel)
{
	GeRay ray;

	switch(voxel->geom) {
		case GEOM_SPH1D:
			ray = GeRay_Rand_sph1d(rng, voxel);
			break;

		case GEOM_SPH3D:
			ray = GeRay_Rand_sph3d(rng, voxel);
			break;
						
		case GEOM_REC3D:
			ray = GeRay_Rand_rec3d(rng, voxel);
			break;
			
		case GEOM_CYL3D:
			ray = GeRay_Rand_cyl3d(rng, voxel);
			break;

		default: /* Shouldn't reach here */
			Deb_PRINT("Uh oh, voxel->geom holds an unidentified geometry code '%d'\n", voxel->geom);
			Deb_ASSERT(0);
	}

	return ray;
}

/*----------------------------------------------------------------------------*/

GeRay GeRay_Rand_sph1d(gsl_rng *rng, const GeVox *voxel)
{
	
	double r, phi, cost, sint,
		r_in = GeVec3_X(voxel->min, 0),
		r_out = GeVec3_X(voxel->max, 0);

	#define ONE_THIRD (0.3333333333333333)

	/* Calculate random ray origin in spherical coordinates */
	if(r_in > 0) {
		r = r_in * pow((1.0 + PRAND() * (pow(r_out/r_in, 3.0) - 1.0)), ONE_THIRD);
	}
	else {
		r = r_out * pow(PRAND(), ONE_THIRD);
	}
	
	#undef ONE_THIRD

	/* Reset ray */
        GeRay ray;
	Mem_BZERO(&ray);

	/* Since this is a 1D problem, starting every ray from the
	   +Z axis is good enough */
	GeRay_E(ray, 2) = r;

	/* Generate random 3D direction */
	Num_RanDir3D(rng, &cost, &sint, &phi);

	/* Convert to rectangular coordinates */
	GeRay_D(ray, 0) = sint * cos(phi);
	GeRay_D(ray, 1) = sint * sin(phi);
	GeRay_D(ray, 2) = cost;

	return ray;
}
/*----------------------------------------------------------------------------*/

GeRay GeRay_Rand_sph3d(gsl_rng *rng, const GeVox *voxel)
/* Generate a randomly directed ray starting from a random position within
   the voxel */
{
	
	double Er, Et, Ep;
	double r_in = GeVec3_X(voxel->min, 0),
		r_out = GeVec3_X(voxel->max, 0),
		theta_in = GeVec3_X(voxel->min, 1),
		theta_out = GeVec3_X(voxel->max, 1),
		phi_in = GeVec3_X(voxel->min, 2),
		phi_out = GeVec3_X(voxel->max, 2);
	double phi, cost, sint;
	#define ONE_THIRD (0.3333333333333333)
	#define USE_RNG 1

	/* Reset ray */
        GeRay ray;
	Mem_BZERO(&ray);
	
	

	/* Set random ray origin in rectangular coordinates */
	#if USE_RNG
	if( r_in != 0.0 )
		Er = r_in * pow((1.0 + PRAND() * (pow(r_out/r_in, 3.0) - 1.0)), ONE_THIRD);
	else
		Er = r_out * pow(PRAND(), ONE_THIRD);
	
	if( theta_in != 0.0)
		Et = acos( (cos(theta_in)-1.0) * (1.0+PRAND()*((cos(theta_out)-1.0)/(cos(theta_in)-1.0)-1.0)) + 1.0 );
	else
		Et = acos( 1.0+PRAND()*(cos(theta_out)-1.0) );
	Ep = phi_in + PRAND()*(phi_out-phi_in);
	#else
	Er = GeVec3_X(voxel->cen, 0);
	Et = GeVec3_X(voxel->cen, 1);
	Ep = GeVec3_X(voxel->cen, 2);
	#endif
	
	#if 0 //debug
	Deb_PRINT("Er=%g r_in=%g r_out=%g\n",Er,r_in,r_out);
	Deb_PRINT("Et=%g theta_in=%g theta_out=%g\n",Et,theta_in,theta_out);
	Deb_PRINT("Ep=%g phi_in=%g phi_out=%g\n",Ep,phi_in,phi_out);
	#endif
	GeRay_E(ray, 0) = Er * sin(Et) * cos(Ep);
	GeRay_E(ray, 1) = Er * sin(Et) * sin(Ep);
	GeRay_E(ray, 2) = Er * cos(Et);


	/* Generate random 3D direction */
	Num_RanDir3D(rng, &cost, &sint, &phi);

	/* Convert to rectangular coordinates */
	GeRay_D(ray, 0) = sint * cos(phi);
	GeRay_D(ray, 1) = sint * sin(phi);
	GeRay_D(ray, 2) = cost;

	return ray;
}


/*----------------------------------------------------------------------------*/

GeRay GeRay_Rand_rec3d(gsl_rng *rng, const GeVox *voxel)
/* Generate a randomly directed ray starting from a random position within
   the voxel */
{
	GeRay ray;
	double phi, cost, sint;
	size_t i;

	/* Reset ray */
	Mem_BZERO(&ray);

	#if 1
	/* Set random ray origin in rectangular coordinates */
	for(i = 0; i < 3; i++) {
		GeRay_E(ray, i) = GeVec3_X(voxel->cen, i) + (RAND() - 0.5) * GeVec3_X(voxel->delta, i);
	}
	#else
	for(i = 0; i < 3; i++) {
		GeRay_E(ray, i) = GeVec3_X(voxel->min, i) + RAND() * GeVec3_X(voxel->delta, i);
	}
	#endif

	/* Generate random 3D direction */
	Num_RanDir3D(rng, &cost, &sint, &phi);

	/* Convert to rectangular coordinates */
	GeRay_D(ray, 0) = sint * cos(phi);
	GeRay_D(ray, 1) = sint * sin(phi);
	GeRay_D(ray, 2) = cost;

	return ray;
}
/*----------------------------------------------------------------------------*/

GeRay GeRay_Rand_cyl3d(gsl_rng *rng, const GeVox *voxel)
/* Generate a randomly directed ray starting from a random position within
   the voxel */
{
	
	double ERc, Ez, Ep;
	double 
		Rc_in = GeVec3_X(voxel->min, 0),
		Rc_out = GeVec3_X(voxel->max, 0),
		phi_in = GeVec3_X(voxel->min, 1),
		phi_out = GeVec3_X(voxel->max, 1),
		Hz_l = GeVec3_X(voxel->min, 2),
		Hz_u = GeVec3_X(voxel->max, 2);
		
	double phi, cost, sint;
	#define USE_RNG 1

	/* Reset ray */
        GeRay ray;
	Mem_BZERO(&ray);
	
	

	/* Set random ray origin in rectangular coordinates */
	#if USE_RNG
	if( Rc_in != 0.0 )
		ERc = Rc_in * pow((1.0 + PRAND() * (pow(Rc_out/Rc_in, 2.0) - 1.0)), 0.5);
	else
		ERc = Rc_out * pow( PRAND() , 0.5);
	Ep = phi_in + PRAND()*(phi_out-phi_in);
	Ez = Hz_l + PRAND()*(Hz_u-Hz_l);
	#else
	ERc = GeVec3_X(voxel->cen, 0);
	Ep = GeVec3_X(voxel->cen, 1);
	Ez = GeVec3_X(voxel->cen, 2);
	#endif
	
	#if 0 //debug
	Deb_PRINT("Er=%g r_in=%g r_out=%g\n",Er,r_in,r_out);
	Deb_PRINT("Et=%g theta_in=%g theta_out=%g\n",Et,theta_in,theta_out);
	Deb_PRINT("Ep=%g phi_in=%g phi_out=%g\n",Ep,phi_in,phi_out);
	#endif
	GeRay_E(ray, 0) = ERc * cos(Ep);
	GeRay_E(ray, 1) = ERc * sin(Ep);
	GeRay_E(ray, 2) = Ez;


	/* Generate random 3D direction */
	Num_RanDir3D(rng, &cost, &sint, &phi);

	/* Convert to rectangular coordinates */
	GeRay_D(ray, 0) = sint * cos(phi);
	GeRay_D(ray, 1) = sint * sin(phi);
	GeRay_D(ray, 2) = cost;

	return ray;
}
/*----------------------------------------------------------------------------*/
#if Sp_MIRSUPPORT
void GeVox_Cpgplot(const GeVox *voxel, const GeCam *cam)
{
	int idum;

	switch(voxel->geom) {
		case GEOM_SPH1D:
			cpgqfs(&idum);
			cpgsfs(2);
			GeVox_Cpgplot_sph1d(voxel);
			cpgsfs(idum);
			break;

		case GEOM_REC3D:
			GeVox_Cpgplot_rec(voxel, cam);
			break;

		default: /* Shouldn't happen */
			Deb_ASSERT(0);
	}

	return;
}

/*----------------------------------------------------------------------------*/

static void GeVox_Cpgplot_sph1d(const GeVox *voxel)
/* Draw a circle to represent projection/cross section of sphere */
{
	Cpg_Circ(0.0, 0.0, GeVec3_X(voxel->max, 0));

	return;
}

/*----------------------------------------------------------------------------*/

static void GeVox_Cpgplot_rec(const GeVox *voxel, const GeCam *cam)
/*
            Z
	    |
	    |
	    |     /
	    4----------7
	   /|   /     /|
	  / |  /     / |
	 /  | /     /  |
	5----------6   |  
--------|---0------|---3---------Y
	|  /|      |  /
	| / |      | /
	|/  |      |/
	1----------2
       /    |
      /     |
     /
    X

  point 0: (x_min, y_min, z_min)
  point 6: (x_max, y_max, z_max)

*/
{
	size_t i;
	GeVec3_d min = voxel->min, max = voxel->max;
	GeVec3_d verts[8];

	verts[0] = min;
	verts[1] = min; verts[1].x[1] = GeVec3_X(max, 1);
	verts[2] = max; verts[2].x[2] = GeVec3_X(min, 2);
	verts[3] = min; verts[3].x[0] = GeVec3_X(max, 0);

	verts[4] = min; verts[4].x[2] = GeVec3_X(max, 2);
	verts[5] = max; verts[5].x[0] = GeVec3_X(min, 0);
	verts[6] = max;
	verts[7] = max; verts[7].x[1] = GeVec3_X(min, 1);

	for(i = 0; i < 4; i++) {
		GeVec3_Cpgline2(&verts[i], &verts[(i+1)%4], cam); /* lower plane */
		GeVec3_Cpgline2(&verts[4+i], &verts[4+(i+1)%4], cam); /* upper plane */
		GeVec3_Cpgline2(&verts[i], &verts[i+4], cam); /* sides */
	}

	return;
}

/*----------------------------------------------------------------------------*/

void GeVec3_Cpgpt1(const GeVec3_d *vec, int sty, const GeCam *cam)
/* By default project on y-z plane */
{
	GeVec3_d v;

	if(cam) {
		v = GeVec3_Rotate(vec, cam);
	}
	else {
		v = *vec;
	}

	Cpg_Pt1(GeVec3_X(v, 1), GeVec3_X(v, 2), sty);

	return;
}

/*----------------------------------------------------------------------------*/

void GeVec3_Cpgline2(const GeVec3_d *v1, const GeVec3_d *v2, const GeCam *cam)
/* By default project on y-z plane */
{
	float xbuf[2], zbuf[2];
	GeVec3_d v11, v22;
	
	/* Optionally rotate coordinates */
	if(cam) {
		v11 = GeVec3_Rotate(v1, cam);
		v22 = GeVec3_Rotate(v2, cam);
	}
	else {
		v11 = *v1;
		v22 = *v2;
	}

	xbuf[0] = (float)v11.x[1];
	xbuf[1] = (float)v22.x[1];
	zbuf[0] = (float)v11.x[2];
	zbuf[1] = (float)v22.x[2];

	cpgline(2, xbuf, zbuf);

	return;
}

/*----------------------------------------------------------------------------*/

void GeVec3_Cpgarro2(const GeVec3_d *v1, const GeVec3_d *v2, const GeCam *cam)
/* By default project on y-z plane */
{
	GeVec3_d v11, v22;

	if(cam) {
		v11 = GeVec3_Rotate(v1, cam);
		v22 = GeVec3_Rotate(v2, cam);
	}
	else {
		v11 = *v1;
		v22 = *v2;
	}

	Cpg_Arro(v11.x[1], v11.x[2], v22.x[1], v22.x[2]);

	return;
}

/*----------------------------------------------------------------------------*/

void GeRay_Cpgarro(const GeRay *ray, const GeCam *cam)
/* Plot ray as an arrow according to origin and direction vector */
{
	size_t i;
	GeVec3_d head = GeVec3_INIT(0, 0, 0);

	for(i = 0; i < 3; i++)
		head.x[i] = ray->e.x[i] + ray->d.x[i];

	GeVec3_Cpgarro2(&ray->e, &head, cam);

	return;
}

/*----------------------------------------------------------------------------*/

void GeVox_Cpgenv(const GeVox *voxel)
/* Set size of plotting window according to voxel size
 * default is to plot on y-z plane!!! */
{
	const double
		*min = voxel->min.x,
	      	*max = voxel->max.x,
	      	dx = max[1]-min[1],
	      	dz = max[2]-min[2];

	switch(voxel->geom) {
		case GEOM_SPH1D:
			Cpg_Env(-max[0]*1.1, max[0]*1.1, -max[0]*1.1, max[0]*1.1, 1, 0);
			break;

		case GEOM_REC3D:
			Cpg_Env(min[1]-dx*1.5, max[1]+dx*0.5, min[2]-dz*1.5, max[2]+dz*0.5, 1, 0);
			break;

		default: /* Shouldn't happen */
			Deb_ASSERT(0);
	}

	return;
}
#endif


/*----------------------------------------------------------------------------*/
size_t GeVox_VertIndex2Pos(size_t i, size_t j, size_t k)
/* 
            Z
	    |
	    |
	    |     /
	    4----------7
	   /|   /     /|
	  / |  /     / |
	 /  | /     /  |
	5----------6   |  
--------|---0------|---3---------Y
	|  /|      |  /
	| / |      | /
	|/  |      |/
	1----------2
       /    |
      /     |
     /
    X
 *
 * k | 0 1 0 1 0 1 0 1
 * j | 0 0 1 1 0 0 1 1
 * i | 0 0 0 0 1 1 1 1
 */
{
	static size_t pos[8] = {
		0, /* (0, 0, 0) */
		4, /* (0, 0, 1) */
		3, /* (0, 1, 0) */
		7, /* (0, 1, 1) */
		1, /* (1, 0, 0) */
		5, /* (1, 0, 1) */
		2, /* (1, 1, 0) */
		6, /* (1, 1, 1) */
	};

	return pos[k + 2 * (j + 2 * i)];
}


/*----------------------------------------------------------------------------*/

GeVox GeVox_Init(int geom, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
{
	GeVox vox = GeVox_INIT(geom, xmin, ymin, zmin, xmax, ymax, zmax);

	return vox;
}

/*----------------------------------------------------------------------------*/

GeVox GeVox_Init2(int geom, GeVec3_d min, GeVec3_d max)
{
	GeVox vox = GeVox_INIT(
		geom,
		GeVec3_X(min, 0),
		GeVec3_X(min, 1),
		GeVec3_X(min, 2),
		GeVec3_X(max, 0),
		GeVec3_X(max, 1),
		GeVec3_X(max, 2)
	);

	return vox;
}

/*----------------------------------------------------------------------------*/

void GeRay_AWInit(GeRay *ray, const GeVox *voxel)
{
	double t;
	size_t i, iplane;
	static const GeVec3_d normal[6] = {
		GeVec3_INIT(-1,  0,  0),
		GeVec3_INIT( 1,  0,  0),
		GeVec3_INIT( 0, -1,  0),
		GeVec3_INIT( 0,  1,  0),
		GeVec3_INIT( 0,  0, -1),
		GeVec3_INIT( 0,  0,  1)
	};
	const GeVec3_d *q = NULL;

	for(i = 0; i < 3; i++) {
		/*
		Calculate intersections with planes on all 3 axes
	          i == 0: x axis
		  i == 1: y axis
		  i == 2: z axis
		*/
		if(GeRay_D(*ray, i) == 0.0) {
			/* GeRay does not intersect either plane on this axis */
			GeVec3_X(ray->tMax, i) = HUGE_VAL;
			GeVec3_X(ray->tDelta, i) = 0;
			continue;
		}
		else {
			if(GeRay_D(*ray, i) > 0) {
				iplane = i * 2 + 1;
				q = &voxel->max;
			}
			else {
				iplane = i * 2;
				q = &voxel->min;
			}
			t = GeRay_IntersectPlane(ray, &normal[iplane], q);

			/* Set tMax for this axis */
			GeVec3_X(ray->tMax, i) = t;

			/* Calculate tDelta */
			GeVec3_X(ray->tDelta, i) = GeVec3_X(voxel->delta, i) / fabs(GeVec3_X(ray->d, i));
		}
	}

	return;
}

/*----------------------------------------------------------------------------*/

void GeRay_AWTraverse(GeRay *ray, double *dt, size_t *plane)
{
	double tmin;

	#define tMaxX (GeVec3_X(ray->tMax, 0))
	#define tMaxY (GeVec3_X(ray->tMax, 1))
	#define tMaxZ (GeVec3_X(ray->tMax, 2))
	#define tDeltaX (GeVec3_X(ray->tDelta, 0))
	#define tDeltaY (GeVec3_X(ray->tDelta, 1))
	#define tDeltaZ (GeVec3_X(ray->tDelta, 2))

	if(tMaxX < tMaxY) {
		if(tMaxX < tMaxZ) {
			if(GeRay_D(*ray, 0) > 0) {
				*plane = 1;
			}
			else {
				*plane = 0;
			}
			tmin = tMaxX;
			tMaxX += tDeltaX;
		}
		else {
			if(GeRay_D(*ray, 2) > 0) {
				*plane = 5;
			}
			else {
				*plane = 4;
			}
			tmin = tMaxZ;
			tMaxZ += tDeltaZ;
		}
	}
	else {
		if(tMaxY < tMaxZ) {
			if(GeRay_D(*ray, 1) > 0) {
				*plane = 3;
			}
			else {
				*plane = 2;
			}
			tmin = tMaxY;
			tMaxY += tDeltaY;
		}
		else {
			if(GeRay_D(*ray, 2) > 0) {
				*plane = 5;
			}
			else {
				*plane = 4;
			}
			tmin = tMaxZ;
			tMaxZ += tDeltaZ;
		}
	}

	*dt = tmin - ray->t;
	ray->t = tmin;

	return;
}

/*----------------------------------------------------------------------------*/

GeVec3_d GeRay_AWPos(const GeRay *ray)
/* Calculate current position based on t and ray origin */
{
	GeRay tmp_ray;

	tmp_ray = GeRay_Inc(ray, ray->t);

	return tmp_ray.e;
}





