#include <Python.h>
#include <numpy/arrayobject.h>
#include "sparx.h"


/* Global parameters */
SpParm Sp_parm;

#define USEUP_SELF_ARGS()\
	{self = NULL; args = NULL;}

/* Function prototypes */
//PyMODINIT_FUNC PY_MOD_INIT_FUNC(void);
PyMODINIT_FUNC init_sparx (void);
static PyObject *get_prog(PyObject *self, PyObject *args);
static PyObject *set_prog(PyObject *self, PyObject *args);
static PyObject *get_mpi_info(PyObject *self, PyObject *args);
static PyObject *task_template(PyObject *self, PyObject *args);
static PyObject *task_amc(PyObject *self, PyObject *args);

static PyObject *task_telsim(PyObject *self, PyObject *args);
static PyObject *task_contobs(PyObject *self, PyObject *args);
static PyObject *task_coldens(PyObject *self, PyObject *args);
static PyObject *task_visual(PyObject *self, PyObject *args);
static PyObject *task_pops2ascii(PyObject *self, PyObject *args);

static PyObject *task_pygrid(PyObject *self, PyObject *args);
static PyObject *test_fft(PyObject *self, PyObject *args);
#if Sp_MIRSUPPORT
static PyObject *load_mir_xyv(PyObject *self, PyObject *args);
static PyObject *load_mir_xyv2(PyObject *self, PyObject *args);
#endif
static PyObject *test_array(PyObject *self, PyObject *args);
static PyObject *test_planck(PyObject *self, PyObject *args);
#ifdef HAVE_MPI
static PyObject *init_mpi(PyObject *self, PyObject *args);
static PyObject *finalize_mpi(PyObject *self, PyObject *args);
#endif

/* Python module method table */
static PyMethodDef _SPARXMethods[] = {
	{"get_prog", get_prog, METH_VARARGS, "Get program name."},
	{"set_prog", set_prog, METH_VARARGS, "Set program name."},
#ifdef HAVE_MPI
	{"init_mpi", init_mpi, METH_VARARGS, "Initialize MPI environment."},
	{"finalize_mpi", finalize_mpi, METH_VARARGS, "Finalize MPI environment."},
#endif
	{"get_mpi_info", get_mpi_info, METH_VARARGS, "Get MPI rank and size."},
	{"task_amc", task_amc, METH_VARARGS, "Accelerated Monte Carlo solver for non-LTE molecular excitation."},
	{"task_telsim", task_telsim, METH_VARARGS, "Observation synthesizer."},
        {"task_contobs", task_contobs, METH_VARARGS, "Observation synthesizer."},
    	{"task_coldens", task_coldens, METH_VARARGS, "Column density tracer."},
        {"task_visual", task_visual, METH_VARARGS, "Visualization postprocessing."},
        {"task_pops2ascii", task_pops2ascii, METH_VARARGS, "Level populations writer."},
	{"task_pygrid", task_pygrid, METH_VARARGS, "Python interface for generate SPARX models."},
	{"task_template", task_template, METH_VARARGS, "Template for C tasks."},
	{"test_fft", test_fft, METH_VARARGS, "Test Fast Fourier Transform."},
	{"test_array", test_array, METH_VARARGS, "Test NumPy array creation."},
	{"test_planck", test_planck, METH_VARARGS, "Test the Planck function."},
#if Sp_MIRSUPPORT
	{"load_mir_xyv", load_mir_xyv, METH_VARARGS, "Load a Miriad XYV datacube."},
	{"load_mir_xyv2", load_mir_xyv2, METH_VARARGS, "Load a Miriad XYV datacube."},
#endif
	{NULL, NULL, 0, NULL} /* Sentinel */
};
 
/*----------------------------------------------------------------------------*/

PyMODINIT_FUNC init_sparx (void)
/* Module init function */
{
	PyObject *o, *mod, *dic, *sparx;

	/* Reset global parameters */
	Mem_BZERO(&Sp_parm);

	/* Set default values */
        strncpy(Sp_parm.prog, Sp_SPARX_VERSION, BUFSIZ);
	Sp_parm.debug = 0;
	Sp_parm.verbose = 1;
	Sp_parm.out_fp = stdout;
	Sp_parm.err_fp = stderr;

	/* Default settings for single process runs */
	Sp_parm.mpi_rank = 0;
	Sp_parm.mpi_size = 1;

	/* Initialize module */
        char static_library[7] = "._sparx";
        char *static_library_path = 
          malloc(1 + strlen(Sp_SPARX_VERSION) + strlen(static_library) );
        strcpy(static_library_path, Sp_SPARX_VERSION);
        strcat(static_library_path, static_library);
        mod = Py_InitModule( static_library_path, _SPARXMethods);
        free(static_library_path);
	/* Necessary for NumPy */
	import_array();

	/* Add flag indicating whether the module has been compiled
	 * with MPI support */
	dic = PyModule_GetDict(mod);
	#ifdef HAVE_MPI
	o = PyBool_FromLong(1);
	#else
	o = PyBool_FromLong(0);
	#endif
	PyDict_SetItemString(dic, "HAVE_MPI", o);
	Py_DECREF(o);

	/* Import sparx module */
        sparx = PyImport_ImportModule(Sp_SPARX_VERSION);

	/* Import MOLEC_DIR from module sparx */
	o = PyObject_GetAttrString(sparx, "MOLEC_DIR");
	Deb_ASSERT(o != NULL);
	Sp_parm.molec_path = Mem_STRDUP(Sp_PYSTR(o));
	Py_DECREF(o);

	/* Import KAPPA_DIR from module sparx */
	o = PyObject_GetAttrString(sparx, "KAPPA_DIR");
	Deb_ASSERT(o != NULL);
	Sp_parm.kappa_path = Mem_STRDUP(Sp_PYSTR(o));
	Py_DECREF(o);

	/* Cleanup */
	Py_DECREF(sparx);

	return;
}

/*----------------------------------------------------------------------------*/

#ifdef HAVE_MPI
static PyObject *init_mpi(PyObject *self, PyObject *args)
/* Init MPI environment if not already done. If init_complete is true, do not
 * do any initialization and simply return rank and size.
 */
{
	int sts = 0;
	static int init_complete = 0;

	USEUP_SELF_ARGS();

	if(!init_complete) {
		sts = MPI_Init(NULL, NULL); /* Init MPI */
		if(!sts) sts = MPI_Comm_rank(MPI_COMM_WORLD, &Sp_parm.mpi_rank); /* Get rank of process */
		if(!sts) sts = MPI_Comm_size(MPI_COMM_WORLD, &Sp_parm.mpi_size); /* Get size of process pool */
		if(!sts) sts = MPI_Barrier(MPI_COMM_WORLD); /* Sync all processes */
		if(!sts) init_complete = 1;
	}

	/* Return rank and size if successful */
	if(!sts) return Py_BuildValue("(i,i)", Sp_parm.mpi_rank, Sp_parm.mpi_size);
	else return NULL;
}
#endif

/*----------------------------------------------------------------------------*/

#ifdef HAVE_MPI
static PyObject *finalize_mpi(PyObject *self, PyObject *args)
/* Finalize MPI environment */
{
	USEUP_SELF_ARGS();

	MPI_Finalize();

	Py_RETURN_NONE;
}
#endif

/*----------------------------------------------------------------------------*/

static PyObject *get_mpi_info(PyObject *self, PyObject *args)
{
	USEUP_SELF_ARGS();

	return Py_BuildValue("(i,i)", Sp_parm.mpi_rank, Sp_parm.mpi_size);
}

/*----------------------------------------------------------------------------*/

static PyObject *task_template(PyObject *self, PyObject *args)
{
	int status = 0;

	USEUP_SELF_ARGS();

	status = SpTask_Template();

	if(!status) Py_RETURN_NONE;
	else return NULL;
}

/*----------------------------------------------------------------------------*/

static PyObject *get_prog(PyObject *self, PyObject *args)
{
	USEUP_SELF_ARGS();

	return Py_BuildValue("s", Sp_parm.prog);
}

/*----------------------------------------------------------------------------*/

static PyObject *set_prog(PyObject *self, PyObject *args)
{
	const char *prog;

	if (!PyArg_ParseTuple(args, "s", &prog))
		return NULL;

	strncpy(Sp_parm.prog, prog, BUFSIZ);

	USEUP_SELF_ARGS();

	Py_RETURN_NONE;
}

/**
 ** Tasks
 **/
static PyObject *task_amc(PyObject *self, PyObject *args)
{
	int status = 0;

	USEUP_SELF_ARGS();

	status = SpTask_Amc();

	if(!status) Py_RETURN_NONE;
	else return NULL;
}

/*----------------------------------------------------------------------------*/

static PyObject *task_telsim(PyObject *self, PyObject *args)
{
	int status = 0;

	USEUP_SELF_ARGS();

	status = SpTask_Telsim();

	if(!status) Py_RETURN_NONE;
	else return NULL;
}

/*----------------------------------------------------------------------------*/

static PyObject *task_contobs(PyObject *self, PyObject *args)
{
    int status = 0;
    
    USEUP_SELF_ARGS();
    
    status = SpTask_ContObs();
    
    if(!status) Py_RETURN_NONE;
    else return NULL;
}
/*----------------------------------------------------------------------------*/

static PyObject *task_coldens(PyObject *self, PyObject *args)
{
	int status = 0;

	USEUP_SELF_ARGS();

	status = SpTask_ColDens();

	if(!status) Py_RETURN_NONE;
	else return NULL;
}
/*----------------------------------------------------------------------------*/

static PyObject *task_visual(PyObject *self, PyObject *args)
{
	int status = 0;

	USEUP_SELF_ARGS();

	status = SpTask_Visual();

	if(!status) Py_RETURN_NONE;
	else return NULL;
}
/*----------------------------------------------------------------------------*/

static PyObject *task_pops2ascii(PyObject *self, PyObject *args)
{
    int status = 0;
    
    USEUP_SELF_ARGS();
    
    status = SpTask_Pops2ASCII();
    
    if(!status) Py_RETURN_NONE;
    else return NULL;
}

/*----------------------------------------------------------------------------*/

static PyObject *task_pygrid(PyObject *self, PyObject *args)
{
	int status = 0;

	USEUP_SELF_ARGS();

	//status = SpTask_PyGrid();

	if(!status) Py_RETURN_NONE;
	else return NULL;
}

/**
 ** Testing
 **/
static PyObject *test_fft(PyObject *self, PyObject *args)
{
	//SpTest_FFT();

	USEUP_SELF_ARGS();

	return Py_BuildValue("i", 0);
}

/*----------------------------------------------------------------------------*/

static PyObject *test_planck(PyObject *self, PyObject *args)
{
	size_t i;
	const size_t n = 100;
	double T_k, nu_min, nu_max, nu_delt, *nu_arr = 0, *planck_arr = 0;

	nu_arr = Mem_CALLOC(n, nu_arr);
	planck_arr = Mem_CALLOC(n, planck_arr);

	/* Get arguments */
	PyArg_ParseTuple(args, "ddd", &T_k, &nu_min, &nu_max);

	T_k = 8e-3; /* [K] */
	nu_min = 1e3; /* [Hz] */
	nu_max = 1e12; /* [Hz] */

	/* Logarithmically space frequency */
	nu_delt = (log10(nu_max) - log10(nu_min)) / (double)n;

	for(i = 0; i < n; i++) {
		nu_arr[i] = pow(10.0, log10(nu_min) + (double)i * nu_delt);
		planck_arr[i] = Phys_PlanckFunc(nu_arr[i], T_k);

		printf("%20.5e %20.5e\n", nu_arr[i], planck_arr[i]);
	}

	free(nu_arr);
	free(planck_arr);

	USEUP_SELF_ARGS();

	Py_RETURN_NONE;
}


#if Sp_MIRSUPPORT
/**
 ** Utilities
 **/
static PyObject *load_mir_xyv(PyObject *self, PyObject *args)
/* Load a Miriad XYV image data cube and return a dictionary containing
 * all headers and data */
{
	int sts = 0;
	const char *fname;
	char *bunit = 0;
	MirImg_Axis x, y, v;
	double *cube = 0;
	PyObject *np_cube = 0, *ret = 0;
	size_t i, j, k;
	npy_intp dims[3];

	#define ARRAY(obj, i, j, k)\
		(*((double *)PyWrArray_GetPtr3((obj), (i), (j), (k))))

	/* The only argument is the filename */
	if (!PyArg_ParseTuple(args, "s", &fname))
		sts = 1;

	if(!sts) {
		/* Load Miriad image into C data structure */
		cube = MirImg_LoadXYV(fname, &x, &y, &v, &bunit);

                /* Create NumPy cubic array */
                dims[0] = x.n;
                dims[1] = y.n;
                dims[2] = v.n;
                np_cube = PyArray_ZEROS(3, dims, NPY_FLOAT, 0);
                
		/* Load image data into NumPy array */
		for(i = 0; i < x.n; i++) {
			for(j = 0; j < y.n; j++) {
				for(k = 0; k < v.n; k++) {
					ARRAY(np_cube, i, j, k) = cube[k + v.n * (j + y.n * i)];
				}
			}
		}

		/* Create dictionary for returning data */
		if(!(ret = PyDict_New()))
			sts = 1;

		free(cube);
	}

	/* Insert np_cube into the dictionary: the DECREF afterwords is VERY IMPORTANT! 
	   Neglecting to DECREF any objects inserted into a dictionary can cause serious
	   memory leaks! */
	if(!sts) sts = PyDict_SetItemString(ret, "cube", np_cube); Py_DECREF(np_cube);

	/* Insert non-Python values into dictionary */
	if(!sts) sts = PyDict_SetItemString(ret, "bunit", Py_BuildValue("s", bunit));
	if(!sts) sts = PyDict_SetItemString(ret, "xcrpix", Py_BuildValue("d", x.crpix));
	if(!sts) sts = PyDict_SetItemString(ret, "xcrval", Py_BuildValue("d", x.crval));
	if(!sts) sts = PyDict_SetItemString(ret, "xcdelt", Py_BuildValue("d", x.delt));
	if(!sts) sts = PyDict_SetItemString(ret, "xcsize", Py_BuildValue("k", (unsigned long)x.n));
	if(!sts) sts = PyDict_SetItemString(ret, "ycrpix", Py_BuildValue("d", y.crpix));
	if(!sts) sts = PyDict_SetItemString(ret, "ycrval", Py_BuildValue("d", y.crval));
	if(!sts) sts = PyDict_SetItemString(ret, "ycdelt", Py_BuildValue("d", y.delt));
	if(!sts) sts = PyDict_SetItemString(ret, "ycsize", Py_BuildValue("k", (unsigned long)y.n));
	if(!sts) sts = PyDict_SetItemString(ret, "vcrpix", Py_BuildValue("d", v.crpix));
	if(!sts) sts = PyDict_SetItemString(ret, "vcrval", Py_BuildValue("d", v.crval));
	if(!sts) sts = PyDict_SetItemString(ret, "vcdelt", Py_BuildValue("d", v.delt));
	if(!sts) sts = PyDict_SetItemString(ret, "vcsize", Py_BuildValue("k", (unsigned long)v.n));

	/* Cleanup */
	if(bunit) free(bunit);

	USEUP_SELF_ARGS();

	return ret;
}

/*----------------------------------------------------------------------------*/

static PyObject *load_mir_xyv2(PyObject *self, PyObject *args)
/* Load a Miriad XYV image data cube and return a dictionary containing
 * all headers and data */
{
	int sts = 0;
	const char *fname;
	char *bunit = 0;
	MirImg_Axis x, y, v;
	double *cube = 0;
	PyObject *np_cube = 0, *numpy = 0, *numpy_zeros = 0, *numpy_float = 0, *ret = 0;
	char format[BUFSIZ];
	size_t i, j, k;

	#define ARRAY(obj, i, j, k)\
		(*((double *)PyWrArray_GetPtr3((obj), (i), (j), (k))))

	/* The only argument is the filename */
	if (!PyArg_ParseTuple(args, "s", &fname))
		sts = 1;

	if(!sts && !(numpy = PyImport_ImportModule("numpy")))
		sts = 1;

	if(!sts && !(numpy_zeros = PyObject_GetAttrString(numpy, "zeros")))
		sts = 1;

	if(!sts && !(numpy_float = PyObject_GetAttrString(numpy, "float")))
		sts = 1;

	if(!sts) {
		/* Load Miriad image into C data structure */
		cube = MirImg_LoadXYV(fname, &x, &y, &v, &bunit);

		/* Create NumPy cubic array */
		snprintf(format, BUFSIZ, "(k,k,k),O");
		np_cube = PyObject_CallFunction(numpy_zeros, format, (unsigned long)x.n, (unsigned long)y.n, (unsigned long)v.n, numpy_float);

		/* Load image data into NumPy array */
		for(i = 0; i < x.n; i++) {
			for(j = 0; j < y.n; j++) {
				for(k = 0; k < v.n; k++) {
					ARRAY(np_cube, i, j, k) = cube[k + v.n * (j + y.n * i)];
				}
			}
		}

		/* Create dictionary for returning data */
		if(!(ret = PyDict_New()))
			sts = 1;
	}

	/* Insert np_cube into the dictionary: the DECREF afterwords is VERY IMPORTANT! 
	   Neglecting to DECREF any objects inserted into a dictionary can cause serious
	   memory leaks! */
	if(!sts) sts = PyDict_SetItemString(ret, "cube", np_cube); Py_DECREF(np_cube);

	/* Insert non-Python values into dictionary */
	if(!sts) sts = PyDict_SetItemString(ret, "bunit", Py_BuildValue("s", bunit));
	if(!sts) sts = PyDict_SetItemString(ret, "xcrpix", Py_BuildValue("d", x.crpix));
	if(!sts) sts = PyDict_SetItemString(ret, "xcrval", Py_BuildValue("d", x.crval));
	if(!sts) sts = PyDict_SetItemString(ret, "xcdelt", Py_BuildValue("d", x.delt));
	if(!sts) sts = PyDict_SetItemString(ret, "xcsize", Py_BuildValue("k", (unsigned long)x.n));
	if(!sts) sts = PyDict_SetItemString(ret, "ycrpix", Py_BuildValue("d", y.crpix));
	if(!sts) sts = PyDict_SetItemString(ret, "ycrval", Py_BuildValue("d", y.crval));
	if(!sts) sts = PyDict_SetItemString(ret, "ycdelt", Py_BuildValue("d", y.delt));
	if(!sts) sts = PyDict_SetItemString(ret, "ycsize", Py_BuildValue("k", (unsigned long)y.n));
	if(!sts) sts = PyDict_SetItemString(ret, "vcrpix", Py_BuildValue("d", v.crpix));
	if(!sts) sts = PyDict_SetItemString(ret, "vcrval", Py_BuildValue("d", v.crval));
	if(!sts) sts = PyDict_SetItemString(ret, "vcdelt", Py_BuildValue("d", v.delt));
	if(!sts) sts = PyDict_SetItemString(ret, "vcsize", Py_BuildValue("k", (unsigned long)v.n));

	/* Cleanup */
	if(cube) free(cube);
	if(bunit) free(bunit);
	Py_XDECREF(numpy_float);
	Py_XDECREF(numpy_zeros);
	Py_XDECREF(numpy);

	USEUP_SELF_ARGS();

	return ret;
}
#endif
/*----------------------------------------------------------------------------*/

static PyObject *test_array(PyObject *self, PyObject *args)
{
	PyObject *a;
	npy_intp dims[2];

	USEUP_SELF_ARGS();

	dims[0] = 500;
	dims[1] = 500;
	a = PyArray_ZEROS(2, dims, NPY_FLOAT, 0);

	return a;
}
















