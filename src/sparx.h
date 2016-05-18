#ifndef __SPARX_H__
#define __SPARX_H__

#define STRINGIZE_(x) #x
#define STRINGIZE(x) STRINGIZE_(x)

/* Root of the sparx directory, holds Python source files
 * needed at runtime */
#ifdef ROOT
	#define Sp_ROOT STRINGIZE(ROOT)
#else
	#define Sp_ROOT "./"
#endif

/* SPARX VERSION */
#define Sp_SPARXVERSION SPARXVERSION
#define Sp_SPARX_VERSION SPARX_VERSION

/* Miriad support */
#define Sp_MIRSUPPORT MIRSUPPORT

/* Maximum number threads for multi-threading */
#ifdef NTHREAD
#define Sp_NTHREAD ((size_t)NTHREAD)
#else
#define Sp_NTHREAD ((size_t)4)
#endif

/* Normalization scale for geometrical objects */
#define Sp_LENFAC PHYS_UNIT_MKS_PC

/* Input handling in sparx relies heavily on embedded Python. Check out
 * The Python/C API (http://www.python.org/doc/2.5/api/api.html) for
 * more details. The Python.h header must be placed at the very beginning
 * of every file that needs to access the API. */
#include <Python.h>

/* MPI header file */
#ifdef HAVE_MPI
#include <mpi.h>
#endif

/* HDF5! */
#include <hdf5.h>
#include <hdf5_hl.h>

/* Standard C library headers */
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <time.h>
#include <string.h>

/* Other in-house routines */
#include "memory.h"
/* Error handling is done through an in-house api. See the source files
 * for more details */
#include "error.h"
#include "fits-and-miriad-wrappers.h"

#if Sp_MIRSUPPORT
#include "cpgplot-wrappers.h"
#endif

#include "python-wrappers.h"
#include "geometry.h"
#include "zone.h"
#include "zone-hdf5.h"
#include "debug.h"
#include "numerical.h"
#include "physics.h"
#include "molec.h"
#include "kappa.h"
#include "vtk-wrapper.h"


/* Globally useful constants */
#define PI Num_PI
#define TWOPI Num_TWOPI
#define T_CMB 2.728
#define GAS_TO_DUST 100.0

/******************************************************************************
 * Global resources
 ******************************************************************************/

/* Definition of the sparx global parameter structure */
typedef struct SpParm {
	char prog[BUFSIZ];
	int debug, verbose, deprecated, user; /* Flags */
	const char *task; /* Task name */
	FILE *out_fp, *err_fp; /* FILE pointer for outputs */
	int mpi_size, mpi_rank; /* MPI stuff */
	char *molec_path, *kappa_path;
} SpParm;

#define Sp_MPIRANK ((size_t)Sp_parm.mpi_rank)
#define Sp_MPISIZE ((size_t)Sp_parm.mpi_size)
#define Sp_MPITAG 0



/******************************************************************************
 * Inputs
 ******************************************************************************/
/* The SpKey structure is meant for defining inputs to a task. */
typedef struct  SpKey {
	const char *name, *format, *deflt, *doc;
} SpKey;

#define Sp_KEY(name, format, deflt, doc)\
	{(name), (format), (deflt), (doc)}

typedef struct SpTask {
	const char *name, *doc;
	int (*func)(void);
	SpKey **keys;
} SpTask;

#define Sp_TASK(name, doc, func, keys)\
	{(name), (doc), (func), (keys)}

void SpInp_PrintKeys(void);

PyObject *SpInp_GetKey(const char *key, const char *file, int line, const char *func);
#define SpInp_GETKEY(key)\
	SpInp_GetKey((key).name, __FILE__, __LINE__, __FUNCTION__)

#define SpInp_GETKEYSTR(key)\
	SpInp_GetKey((key), __FILE__, __LINE__, __FUNCTION__)

PyObject *SpInp_GetTaskKey(const char *task, const char *key, const char *file, int line, const char *func);
#define SpInp_GETTASKKEY(task, key)\
	SpInp_GetTaskKey((task), (key).name, __FILE__, __LINE__, __FUNCTION__)

#define SpInp_TASKGETKEY(key)\
	SpInp_GetTaskKey(Sp_parm.task, (key).name, __FILE__, __LINE__, __FUNCTION__)

#define SpInp_GETTASKKEYSTR(task, key)\
	SpInp_GetTaskKey((task), (key), __FILE__, __LINE__, __FUNCTION__)

int SpPy_Initialize(void);
int SpInp_InitTasks(SpTask *tasks[]);

int SpInp_CheckKeys(void);

#define Sp_PYLST(o, i)\
	PyList_GetItem((o), (Py_ssize_t)(i))
#define Sp_PYDBL(o)\
	PyFloat_AsDouble(o)
#define Sp_PYINT(o)\
	(int)PyInt_AsLong(o)
#define Sp_PYSIZE(o)\
	(size_t)PyInt_AsLong(o)
#define Sp_PYSTR(o)\
	PyString_AsString(o)

/* sparx-argp routines */
void SpArg_Parse(int argc, char *argv[], SpParm *Sp_parm);

/* sparx-physics routines */
typedef struct SpPhys {
	int non_empty_leaf, has_tracer;
	PyObject *velfield;
	const Molec *mol;
	size_t nray; /* For use in amc only */
	size_t ncont; /* For use in telsim only */
	double
		n_H2, /* Molecular gas density (m^-3) */
		T_k, /* Kinetic temperature (K) */
		T_d, /* Dust temperature (K) */
		T_ff, /* Brightness temperature of free-free emission (K) */
		T_bb,
		X_mol, /* Molecular fractional abundance */
		X_pH2,
		X_oH2,
		X_e,
		X_H,
		X_He,
		V_t, /* RMS turbulent velocity */
		width, /* Local velocity width (=sqrt(2)*sigma) (m/s) */
		*pops[Sp_NTHREAD], /* Fractional density of each level */
		
		//*popsold, /* old level population */
		//*J_bar, /* mean intensity */
		
		*cmat, /* nlev x nlev matrix of collisional rates */
		*tau, /* nrad array of average tau */
		ds; /* Path length averaged over all directions */
	struct {
		double I_bb; /* Extra continuum flux */
		double freq, lambda;
		double j, k; /* Sum of all continuum emission and absorption (e.g. dust, fre-free... etc. */
	} *cont;
	GeVec3_d
		v_cen, /* Velocity at voxel center (km/s) */
		b_cen, /* B-field at voxel center (km/s) */
		v_edge[6]; /* Velocity at cell edges (km/s) */
	const Zone *zp; /* Zone containing this set of parameters */
	/* Strings describing continuum opacities: meant to be accessed by SpIO_LoadKappa(),
	   must be either 'powerlaw,%10.3e,%10.3e,%10.3e' or 'table,<filename>' */
	char kapp_d[ZoneH5_KAPPLEN]; /* Dust opacity */ 
	char kapp_ff[ZoneH5_KAPPLEN]; /* Free-free opacity */

	double diff;

} SpPhys;

typedef struct SpPhysParm {
	PyObject *velfield;
	Molec *mol,*mol_hyper;
	double T_cmb, gas_to_dust;

} SpPhysParm;

typedef struct SpModel {
	SpPhysParm parms;
	Zone *grid;

} SpModel;

void *SpPhys_Alloc(const Zone *zp, const void *parms_p);
void SpPhys_Free(void *ptr);
void SpPhys_Fprintf(SpPhys *pp, FILE *fp);
size_t SpPhys_Fwrite(SpPhys *pp, FILE *fp);
size_t SpPhys_Fread(SpPhys *pp, FILE *fp);
void SpPhys_InitMol(SpPhys *pp, const Molec *mol,int popsold);
void SpPhys_InitContWindows(SpPhys *pp, double freq[], size_t nfreq);
void SpPhys_AddDust(SpPhys *pp, int cont, const Kappa *kap, double gas_to_dust);
#define SpPhys_AddDust_mol(pp, kap, gas_to_dust)\
	SpPhys_AddDust((pp), 0, (kap), (gas_to_dust))
#define SpPhys_AddDust_cont(pp, kap, gas_to_dust)\
	SpPhys_AddDust((pp), 1, (kap), (gas_to_dust))
void SpPhys_AddCont(SpPhys *pp, int cont);
#define SpPhys_AddCont_mol(pp)\
	SpPhys_AddCont((pp), 0)
#define SpPhys_AddCont_cont(pp)\
	SpPhys_AddCont((pp), 1)

/* New 08 Aug 2009 */
void SpPhys_AddContinuum(SpPhys *pp, int cont, double T_bb, const Kappa *kap, double rho);
void SpPhys_AddContinuum_d(SpPhys *pp, int cont, double gas_to_dust);
void SpPhys_AddContinuum_ff(SpPhys *pp, int cont);
void SpPhys_SetContinuumIntens_bb(SpPhys *pp, int cont, double T_bb, double I_norm);

void SpPhys_InitCollRates(SpPhys *pp);
double SpPhys_GetCollDens(const SpPhys *pp, int species);
void SpPhys_ProcLamda(Molec *mol);
double SpPhys_Zfunc(const Molec *mol, double T_k);
double SpPhys_BoltzPops(const Molec *mol, size_t lev, double T_k);
void SpPhys_GetMoljk(size_t tid, const SpPhys *pp, size_t tr, double vfac, double *j_nu, double *k_nu);
GeVec3_d SpPhys_GetVgas(const GeVec3_d *pos, const Zone *zone);
GeVec3_d SpPhys_GetBgas(const GeVec3_d *pos, const Zone *zone);
GeVec3_d SpPhys_GetVfunc(const GeRay *ray, double dt, const Zone *zone);
GeVec3_d SpPhys_GetBfunc(const GeRay *ray, double dt, const Zone *zone);

GeVec3_d SpPhys_GetBfac(const GeRay *ray, double dt, const Zone *zone, int debug);
double SpPhys_GetVfac(const GeRay *ray, double dt, double v_los, const Zone *zone, int debug);
GeVec3_d SpPhys_GetVfac2(const GeRay *ray, double dt, const Zone *zone, int debug);
double _SpPhys_BoltzRatio(const Molec *mol, size_t up, size_t lo, double T_k);
#define SpPhys_BoltzRatio(mol, up, lo, T_k)\
	_SpPhys_BoltzRatio((mol), (size_t)(up), (size_t)(lo), (T_k))
void SpPhys_PushPopsStack(double *stack, size_t nstack, double *pops, size_t nlev);
double SpPhys_CalcPopsDiff(const double *stack, size_t nstack, size_t nlev, double minpop);
double SpPhys_CalcLineWidth(const SpPhys *pp);
#define SpPhys_LIMITTAU(tau)\
	{if((tau) < -30.0) {(tau) = -30.0;}}

/* sparx-model routines */
void SpModel_PrintModel(SpModel model);
void SpModel_Cleanup(SpModel model);
void SpModel_InitGrid(SpModel *model, int geom, GeVec3_d min, GeVec3_d max, GeVec3_s ndiv);
void SpModel_GenModel_UniSphere(SpModel *model, double n_H2, double T_k, double X_mol);

/* sparx-util routines */
void SpUtil_Threads(void *(*ThreadFunc)(void *));
int SpUtil_Threads2(size_t nthread, void *(*ThreadFunc)(void *));
int SpUtil_TermThread(void);

#define Sp_CHECKTERMTHREAD()\
	{if(SpUtil_TermThread()) break;}

/* sparx-zone routines */
#define SpZone_ALLOC(parms_ptr)\
	Zone_Alloc(0, 0, SpPhys_Alloc, (parms_ptr))

#define SpZone_FREE(zone)\
	Zone_Free((zone), SpPhys_Free);

#define SpZone_FPRINTF(fp, zone)\
	Zone_Fprintf((fp), (zone), (void (*)(void *, FILE *fp))SpPhys_Fprintf);

#define SpZone_GROW(zone, naxes, parms_ptr)\
	Zone_GrowChildren((zone), (naxes), SpPhys_Alloc, (parms_ptr))

#define SpZone_FWRITE(zone, fp)\
	Zone_Fwrite((zone), (size_t (*)(void *, FILE *))SpPhys_Fwrite, (fp))

#define SpZone_FREAD(parms_ptr, fp)\
	Zone_Fread(SpPhys_Alloc, (parms_ptr), (size_t (*)(void *, FILE *))SpPhys_Fread, (fp))

/* sparx-io routines */
#define Sp_PRINT(...)\
	SpIO_Print(__FILE__, __LINE__, __FUNCTION__, __VA_ARGS__)
#define Sp_PRINTF(...)\
	SpIO_Printf(__FILE__, __LINE__, __FUNCTION__, __VA_ARGS__)
#define Sp_PERR(...)\
	SpIO_Perror(__FILE__, __LINE__, __FUNCTION__, __VA_ARGS__)
#define Sp_PWARN(...)\
	SpIO_Pwarn(__FILE__, __LINE__, __FUNCTION__, __VA_ARGS__)

enum {
	Sp_NEW,
	Sp_OLD,
	Sp_TRUNC
};



typedef struct SpFile {
	char *name;
	FILE *fp;
	hid_t h5f_id;
} SpFile;

SpFile *SpIO_OpenFile(const char *fname, int mode); /* mode=Sp_NEW or Sp_OLD */
int SpIO_OpenFile2(const char *fname, int mode, SpFile **fp);
void SpIO_CloseFile(SpFile *sfp);
int SpIO_OpenModel(const char *fname, const char *popsfname, SpModel *model);
int SpIO_FwriteModel(SpFile *sfp, SpModel model);
int SpIO_FreadModel(const SpFile *sfp, const SpFile *popsfp, SpModel *model);
Molec *SpIO_FreadMolec(const char *molname);
Molec *SpIO_FreadMolec_hyper(const char *molname);
Kappa *SpIO_FreadKappa(const char *name);
Kappa *SpIO_LoadKappa(const char *string);
void SpIO_SetTaskName(const char *str);
void SpIO_Print(const char *file, int line, const char *func, const char *format, ...);
void SpIO_Printf(const char *file, int line, const char *func, const char *format, ...);
void SpIO_Pwarn(const char *file, int line, const char *func, const char *format, ...);
void SpIO_Perror(const char *file, int line, const char *func, const char *format, ...);
ZoneH5_Record SpIO_ZoneToH5Record(const Zone *zone);
void SpIO_ZoneFromH5Record(Zone *zone, ZoneH5_Record record);
int SpIO_H5WriteGrid(hid_t h5f_id, const Zone *zone);
int SpIO_H5ReadGrid(hid_t h5f_id, hid_t popsh5f_id, Zone **zone, SpPhysParm *parms);
int SpIO_H5WritePops(hid_t h5f_id, const Zone *zone);
int SpIO_H5ReadPops(hid_t h5f_id, Zone *zone);
int SpIO_H5WriteTau(hid_t h5f_id, const Zone *zone);
int SpIO_H5ReadTau(hid_t h5f_id, Zone *zone);
int SpIO_H5GetAttribute_string(hid_t h5f_id, const char *obj_name, const char *attr_name, char **attribute);

/* Global variables -- better keep these at the bottom! */
extern struct SpParm Sp_parm;
extern SpTask *Sp_tasks[];

/* Tasks */
int SpTask_Amc(void);
int SpTask_AsciiGrid(void);
int SpTask_Example(void);
int SpTask_Powerlaw(void);
int SpTask_PyGrid(void);
int SpTask_Telsim(void);
int SpTask_Uniform(void);
int SpTask_UVSamp(void);
int SpTask_Template(void);

/* Testing tasks */
int SpTest_Key(void);
int SpTest_Gaussian(void);
int SpTest_XYReadWrite(void);
int SpTest_UVResamp(void);
int SpTest_Interp2d(void);
int SpTest_FFT(void);
int SpTest_Test(void);

/* Keyword inputs */
Kappa *SpInp_GetKey_kappa(const char *name);
const char *SpInp_GetKey_str(const char *name);
MirFile *SpInp_GetKey_miruv(const char *name, const char *mode);
#define SpInp_GetKey_miruv_new(name)\
	SpInp_GetKey_miruv((name), "new");
#define SpInp_GetKey_miruv_old(name)\
	SpInp_GetKey_miruv((name), "old");
MirFile *SpInp_GetKey_mirxy_new(const char *name, size_t nx, size_t ny, size_t nv);
MirFile *SpInp_GetKey_mirxy_old(const char *name, size_t *nx, size_t *ny, size_t *nv);
PyObject *SpInp_GetKey_obj(const char *name);
Molec *SpInp_GetKey_molec(const char *name);
int SpInp_GetKey_model(const char *name, SpModel *model);
SpFile *SpInp_GetKey_spfile(const char *name, int mode);
int SpInp_GetKey_TF(const char *name);
int SpInp_GetKey_int(const char *name);
size_t SpInp_GetKey_size_t(const char *name);
double SpInp_GetKey_dbl(const char *name);

/******************************************************************************
 * Embedded Python
 ******************************************************************************/
PyObject *SpPy_GetMain(const char *symb, const char *file, int line, const char *func);

#define SpPy_GETMAIN(symb)\
	SpPy_GetMain((symb), __FILE__, __LINE__, __FUNCTION__)

#define SpPy_XDECREF(o)\
	{if(o) {Py_XDECREF(o); (o) = NULL;}}

#define SpPy_CHECKEXC(action)\
	{if(PyErr_Occurred()) action;}

void SpPy_CallVoidFunc(const char *func);
int SpPy_GetInput_PyObj(const char *name, PyObject **obj);
int SpPy_CheckOptionalInput(const char *name);
int SpPy_GetInput_int(const char *name, int *value);
int SpPy_GetInput_sizt(const char *name, size_t *value);
int SpPy_GetInput_dbl(const char *name, double *value);
int SpPy_GetInput_bool(const char *name, int *value);
int SpPy_GetInput_model(const char *name, const char *popsname, SpModel *model);
int SpPy_GetInput_molec(const char *name, Molec **molec);
int SpPy_GetInput_molec_hyper(const char *name, Molec **molec);
int SpPy_GetInput_spfile(const char *name, SpFile **fp, int mode);
int SpPy_GetInput_mirxy_new(const char *name, size_t nx, size_t ny, size_t nv, MirFile **fp);

#endif









