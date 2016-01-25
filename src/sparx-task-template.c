#include "sparx.h"

/*
 * This task is meant as a template for C tasks, and
 * for testing all the constructs required for a C task
 * to work.
 *
 * SOP for adding a C task:
 *   1. Write an SpTask_* function which handles the main execution flow
 *   2. Add function prototype to src/sparx.h
 *   3. Add Python/C function to src/sparx-pyext-_sparx.c
 *   4. Add function to the module methods table (_SPARXMethods[])
 *   5. Add install_task entry to lib/sparx/tasks.py
 *
 * Inputs from the user are stored in the sparx.inputs.INPUTS dictionary,
 */

/* Global parameter struct */
static struct glb {
	size_t sizt;
	int intgr;
	double dbl;
} glb;

static void *ExecThread(void *tid_p);

/*----------------------------------------------------------------------------*/

int SpTask_Template(void)
{
	int sts = 0;
	#ifdef HAVE_MPI
	//debug
	printf("proc %d/%d: check!\n", (int)Sp_MPIRANK, (int)Sp_MPISIZE);
	#endif

	Mem_BZERO(&glb);

	#define TEST_SIZT(name)\
		if(!(sts = SpPy_GetInput_sizt(name, &glb.sizt))) {\
			printf("%s=%g\n", name, (double)glb.sizt);\
		}

	#define TEST_INT(name)\
		if(!(sts = SpPy_GetInput_int(name, &glb.intgr))) {\
			printf("%s=%g\n", name, (double)glb.intgr);\
		}

	#define TEST_DBL(name)\
		if(!(sts = SpPy_GetInput_dbl(name, &glb.dbl))) {\
			printf("%s=%g\n", name, (double)glb.dbl);\
		}

	TEST_INT("int");
	TEST_SIZT("pos_int");
	TEST_DBL("angle");
	TEST_DBL("velo");
	TEST_DBL("length");
	PyErr_Clear();
	sts = 0;

	/* This is for checking whether the KeyboardInterrupt exception can
	 * be caught and propagated to the top-level interpreter */
	printf("Send ctrl-c to continue...\n");
	while(1) {
		if(PyErr_CheckSignals()) {
			printf("Signal caught in PyErr_CheckSignals()!\n");
			break;
		}
		sleep(1);
	}
	PyErr_Clear();
	sts = 0;

	/* This is to test whether threads can be safely terminated by ctrl-c
	 * and propagated to the top-level interpreter.
	 * Threads in Python are a bit tricky to use, since CPython implements
	 * the "GIL" or "Global Interpreter Lock" which restricts execution
	 * of the Python VM to a single instruction at a time. Hence if cpu-intensive
	 * operations are to be done within a C extension module, the GIL needs to be
	 * released (using Py_BEGIN_ALLOW_THREADS) so that the top-level interpreter
	 * can tend to other tasks while the (possibly threaded) C code operates.
	 * Once the C code has completed its task, reacquire the GIL
	 * (using Py_END_ALLOW_THREADS) so the rest of the operation can be managed
	 * by the interpreter.
	 *
	 * Caveat: the top-level interpreter is no longer in charge of the
	 * activities of code following the release of the GIL, so signals
	 * and such must be managed explicitly!
	 */
	if(!sts) sts = SpUtil_Threads2(2, ExecThread);

	return sts;
}

/*----------------------------------------------------------------------------*/

static void *ExecThread(void *tid_p)
{
	int i;
	size_t tid;
	
	tid = *((size_t *)tid_p);

	for(i = 0; i < 5; i++) {
		printf("Thread %d: count=%d\n", (int)tid, i);
		if(SpUtil_TermThread()) {
			printf("Thread %d terminated by SpUtil_termthread\n", (int)tid);
			break;
		}

		sleep(1);
	}

	pthread_exit(NULL);
}
















