#include "sparx.h"

/*----------------------------------------------------------------------------*/

Kappa *SpInp_GetKey_kappa(const char *name)
{
	PyObject *Inputs = 0, *o = 0, *o1;
	Kappa *kap = NULL;
	char *sp;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);

	if(o != Py_None) {
		o1 = PyObject_GetAttrString(o, "data");
		sp = Sp_PYSTR(Sp_PYLST(o1, 0));
		if(!strcmp(sp, "plaw"))
			kap = Kap_New_Powerlaw(
				Sp_PYDBL(Sp_PYLST(o1, 1)),
				Sp_PYDBL(Sp_PYLST(o1, 2)),
				Sp_PYDBL(Sp_PYLST(o1, 3))
			);
		else if(!strcmp(sp, "table"))
			kap = SpIO_FreadKappa(Sp_PYSTR(Sp_PYLST(o1, 1)));
		else /* ID string returned by python not recognized */
			Deb_ASSERT(0);

		Py_DECREF(o1);
		Py_DECREF(o);
		Deb_ASSERT(kap != NULL);
	}
        
	Py_DECREF(Inputs);

	return kap;
}

/*----------------------------------------------------------------------------*/

const char *SpInp_GetKey_str(const char *name)
{
	PyObject *Inputs = 0, *o = 0;
	static char buffer[BUFSIZ] = "";
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);

	Mem_BZERO2(buffer, BUFSIZ);
	strncpy(buffer, Sp_PYSTR(o), (size_t)BUFSIZ);

	Py_DECREF(Inputs);
	Py_DECREF(o);

	return buffer;
}

/*----------------------------------------------------------------------------*/

MirFile *SpInp_GetKey_miruv(const char *name, const char *mode)
/* Retrieve a keyword from user input */
{
	PyObject *Inputs = 0, *o = 0;
	MirFile *fp = 0;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);
#if Sp_MIRSUPPORT
	fp = MirUV_Open(Sp_PYSTR(o), mode);
#endif
	Py_DECREF(Inputs);
	Py_DECREF(o);

	return fp;
}

/*----------------------------------------------------------------------------*/

MirFile *SpInp_GetKey_mirxy_new(const char *name, size_t nx, size_t ny, size_t nv)
/* Retrieve a keyword from user input */
{
	PyObject *Inputs = 0, *o = 0;
	MirFile *fp = 0;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);
#if Sp_MIRSUPPORT
	if(o != Py_None) {
		fp = MirXY_Open_new(Sp_PYSTR(o), nx, ny, nv);
		Py_DECREF(o);
	}
#endif
	Py_DECREF(Inputs);

	return fp;
}

/*----------------------------------------------------------------------------*/

MirFile *SpInp_GetKey_mirxy_old(const char *name, size_t *nx, size_t *ny, size_t *nv)
/* Retrieve a keyword from user input */
{
	PyObject *Inputs = 0, *o = 0;
	MirFile *fp = 0;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);
#if Sp_MIRSUPPORT
	fp = MirXY_Open_old(Sp_PYSTR(o), nx, ny, nv);
#endif
	Py_DECREF(Inputs);
	Py_DECREF(o);

	return fp;
}

/*----------------------------------------------------------------------------*/

PyObject *SpInp_GetKey_obj(const char *name)
/* Retrieve a keyword from user input */
{
	PyObject *Inputs = 0, *o = 0;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	if(o == NULL) {
		Sp_PERR("cannot retrieve objct '%s' from Inputs\n", name);
	}
	Deb_ASSERT(o != NULL);

	Py_DECREF(Inputs);

	return o;
}

/*----------------------------------------------------------------------------*/

Molec *SpInp_GetKey_molec(const char *name)
/* Retrieve a keyword from user input */
{
	PyObject *Inputs = 0, *o = 0;
	Molec *mol = NULL;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);

	mol = SpIO_FreadMolec(Sp_PYSTR(o));

	Py_DECREF(Inputs);
	Py_DECREF(o);

	return mol;
}

/*----------------------------------------------------------------------------*/

int SpInp_GetKey_model(const char *name, SpModel *model)
/* Retrieve a keyword from user input */
{
	int status = 0;
	PyObject *Inputs = 0, *o = 0;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);

	status = SpIO_OpenModel(Sp_PYSTR(o), Sp_PYSTR(o), model);

	Py_DECREF(Inputs);
	Py_DECREF(o);

	return status;
}

/*----------------------------------------------------------------------------*/

SpFile *SpInp_GetKey_spfile(const char *name, int mode)
/* Retrieve a keyword from user input */
{
	PyObject *Inputs = 0, *o = 0;
	SpFile *fp;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);

	fp = SpIO_OpenFile(Sp_PYSTR(o), mode);

	Py_DECREF(Inputs);
	Py_DECREF(o);

	return fp;
}

/*----------------------------------------------------------------------------*/

size_t SpInp_GetKey_size_t(const char *name)
/* Retrieve a keyword from user input */
{
	size_t size;
	PyObject *Inputs = 0, *o = 0;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);

	/* Convert to integer */
	if(o != Py_None)
		size = (size_t)PyInt_AsLong(o);
	else
		size = 0;

	Py_DECREF(Inputs);
	Py_DECREF(o);

	return size;
}

/*----------------------------------------------------------------------------*/

int SpInp_GetKey_TF(const char *name)
/* Retrieve a keyword from user input */
{
	int intgr = 0;
	PyObject *Inputs = 0, *o = 0;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);

	/* Convert to integer */
	if(o != Py_None) {
		intgr = (int)PyInt_AsLong(o);
	}

	Py_DECREF(Inputs);
	Py_DECREF(o);

	return intgr;
}

/*----------------------------------------------------------------------------*/

int SpInp_GetKey_int(const char *name)
/* Retrieve a keyword from user input */
{
	int intgr;
	PyObject *Inputs = 0, *o = 0;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);

	/* Convert to integer */
	if(o != Py_None)
		intgr = (int)PyInt_AsLong(o);
	else
		intgr = 0;

	Py_DECREF(Inputs);
	Py_DECREF(o);

	return intgr;
}

/*----------------------------------------------------------------------------*/

double SpInp_GetKey_dbl(const char *name)
/* Retrieve a keyword from user input */
{
	double dbl;
	PyObject *Inputs = 0, *o = 0;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", __FILE__, __LINE__, name);
	Deb_ASSERT(Inputs != NULL);

	o = PyObject_GetAttrString(Inputs, name);
	Deb_ASSERT(o != NULL);

	/* Convert to floating point */
	if(o != Py_None)
		dbl = PyFloat_AsDouble(o);
	else
		dbl = 0;

	Py_DECREF(Inputs);
	Py_DECREF(o);

	return dbl;
}

/*----------------------------------------------------------------------------*/

PyObject *SpInp_GetKey(const char *key, const char *file, int line, const char *func)
/* Retrieve a keyword from user input by doing the following:
 * 1. Check for existence of key.name in Inp_Keys
 * 2. If key exists, check against key.format using Sp_CheckFormat
 * This returns a NEW reference to PyObject *.
 */
{
	int status = 0;
	PyObject *Inputs = 0, *op = 0;
	
	/* Get class Inputs */
	Inputs = SpPy_GetMain("Inputs", file, line, func);

	if(!Inputs)
		Err_SETSTRING("Internal Python error: class Inputs not found in __main__!");
	else {
		/* Check if key has been set in Inputs */
		if(!status && PyObject_HasAttrString(Inputs, key))
			op = PyObject_GetAttrString(Inputs, key);
		else
			Err_SetString(file, line, func, "Error retrieving keyword `%s'.", key);

		/* Remember to DECREF Inputs */
		Py_DECREF(Inputs);
	}

	return op;
}

/*----------------------------------------------------------------------------*/

PyObject *SpInp_GetTaskKey(const char *task, const char *key, const char *file, int line, const char *func)
/* Wrapper for the Python function SpInp_GetTaskKey() */
{
	PyObject *fo, *ko;
	char format[BUFSIZ] = "ss";

	/* Load Python function */
	fo = SpPy_GetMain("PySpInp_GetTaskKey", file, line, func);

	Deb_ASSERT(fo != NULL);

	/* Get and process input key */
	ko = PyObject_CallFunction(fo, format, task, key);

	if(PyErr_Occurred()) {
		Py_XDECREF(ko);
		ko = NULL;
		Err_SetString(file, line, func, "Error processing keyword `%s'", key);
	}

	Py_DECREF(fo);

	return ko;
}

/*----------------------------------------------------------------------------*/

void SpInp_PrintKeys(void)
{
	PyObject *Inputs,
		 *dir = 0, /* new ref */
		 *o; /* Generic pointer */
	Py_ssize_t i, n;

	/* Get class Inputs */
	Inputs = SpPy_GETMAIN("Inputs");

	if(Inputs) {
		/* Get list of attributes in Inputs */
		dir = PyObject_Dir(Inputs);
		n = PyList_Size(dir);

		Sp_PRINT("Variables in Inputs:\n");
		for(i = 0; i < n; i++) {
			o = PyList_GetItem(dir, i);
			Sp_PRINT("  ");
			PyObject_Print(o, Sp_parm.out_fp, 0);
			Sp_PRINT("\n");
		}
	}
	else {
		Sp_PRINT("Error retrieving __main__.Inputs\n");
	}

	Py_XDECREF(dir);
	Py_XDECREF(Inputs);

	return;
}

/*----------------------------------------------------------------------------*/

int SpInp_CheckKeys(void)
/* Check for redundant input */
{
	PyObject *op;

	op = SpPy_GETMAIN("SpInp_CheckKeys");

	Deb_ASSERT(op != NULL);

	PyObject_CallFunction(op, 0);
	Py_DECREF(op);

	if(PyErr_Occurred())
		return Err_SETSTRING("Internal Python error");

	return 0;
}

/*----------------------------------------------------------------------------*/

int SpInp_InitTasks(SpTask *tasks[])
/* Build dictionary of tasks in embedded Python environment. */
{
	size_t i, j;
	SpTask *tp;
	SpKey *kp;

	PyWrRun_SimpleString("Sp_tasks = {}");
	for(i = 0; tasks[i]; i++) {
		tp = tasks[i];
		PyWrRun_SimpleString("Sp_tasks['%s'] = {}", tp->name);
		PyWrRun_SimpleString("Sp_tasks['%s']['doc'] = \"%s\"", tp->name, tp->doc);
		if(tp->keys) {
			PyWrRun_SimpleString("Sp_tasks['%s']['keys'] = {}", tp->name);
			for(j = 0; tp->keys[j]; j++) {
				kp = tp->keys[j];
				PyWrRun_SimpleString("Sp_tasks['%s']['keys']['%s'] = {}", tp->name, kp->name);
				PyWrRun_SimpleString("Sp_tasks['%s']['keys']['%s']['format'] = \"%s\"", tp->name, kp->name, kp->format);
				if(kp->deflt) {
					PyWrRun_SimpleString("Sp_tasks['%s']['keys']['%s']['default'] = %s", tp->name, kp->name, kp->deflt);
				}
				else
					PyWrRun_SimpleString("Sp_tasks['%s']['keys']['%s']['default'] = None", tp->name, kp->name);
				PyWrRun_SimpleString("Sp_tasks['%s']['keys']['%s']['doc'] = \"%s\"", tp->name, kp->name, kp->doc);
			}
		}
		else
			PyWrRun_SimpleString("Sp_tasks['%s']['keys'] = None", tp->name);
	}

	return 0;
}

/*----------------------------------------------------------------------------*/














