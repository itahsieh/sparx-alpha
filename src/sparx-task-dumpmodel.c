#include "sparx.h"

/* Global parameter struct */
static struct glb {
	SpModel model;
} glb;

/* Keywords */
static SpKey
K_SRCF = Sp_KEY("source", "OldFile", 0, "Input source model");

static SpKey *keys[] = {
	&K_SRCF,
	0
};

/* Subroutine prototypes */
static int TaskMain(void);
static int ProcInps(void);

/* Task definition */
SpTask SpTask_dumpmodel = Sp_TASK("dumpmodel", "Print out model", TaskMain, keys);

/*----------------------------------------------------------------------------*/

static int TaskMain(void)
{
	int status = 0;

	/* Reset parms */
	Mem_BZERO(&glb);

	/* Process inputs */
	status = ProcInps();

	if(!status)
		SpModel_PrintModel(glb.model);

	SpModel_Cleanup(glb.model);

	return status;
}

/*----------------------------------------------------------------------------*/

static int ProcInps(void)
{
	int status = 0;
	PyObject *o;

	/* K_SRCF: source model */
	if(!status && !(o = SpInp_TASKGETKEY(K_SRCF)))
		status = 1;
	if(!status)
		status = SpIO_OpenModel(Sp_PYSTR(o), &glb.model);
	SpPy_XDECREF(o);

	return status;
}







