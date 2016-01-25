/*
 * Wrapper routines for the Python/C API
 * All wrapper routines should have prefix `PyWr'.
 */

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "memory.h"
#include "python-wrappers.h"
#include "numpy/arrayobject.h"

/*----------------------------------------------------------------------------*/

int PyWrRun_SimpleString(const char *format, ...)
{
	int status = 0;
	char *string;
	va_list ap;

	va_start(ap, format);
	string = Mem_VSprintf(format, ap);
	va_end(ap);

	status = PyRun_SimpleString(string);

	free(string);

	return status;
}

/*----------------------------------------------------------------------------*/

int PyWrRun_SimpleFile(const char *fname)
{
	FILE *fp;

	fp = fopen(fname, "r");
	if(!fp) {
		return 1;
	}

	return PyRun_SimpleFile(fp, fname);
}

/*----------------------------------------------------------------------------*/

void* PyWrArray_GetPtr3(PyObject* obj, size_t i, size_t j, size_t k)
{
	return PyArray_GETPTR3(obj, (npy_intp)i, (npy_intp)j, (npy_intp)k);
}

/*----------------------------------------------------------------------------*/

void PyWrErr_SetString(PyObject *type, const char *format, ...)
{
	char buffer[BUFSIZ];
	va_list ap;

	va_start(ap, format);
	vsnprintf(buffer, BUFSIZ, format, ap);
	va_end(ap);

	PyErr_SetString(type, buffer);

	return;
}


