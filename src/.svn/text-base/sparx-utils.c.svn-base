#include "sparx.h"
#include <signal.h>

static int SpUtil_termthread = 0; /* Flag for terminating threads */

static void *SpUtil_SigMonThread(void);

/*----------------------------------------------------------------------------*/

void SpUtil_Threads(void *(*ThreadFunc)(void *))
/* Distribute ThreadFunc to Sp_NTHREAD threads, passing the id of each thread
 * as the argument */
{
	pthread_t threads[Sp_NTHREAD];
	pthread_attr_t attr;
	int status;
	size_t i, thread_id[Sp_NTHREAD];

	/* Init thread attribute */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	/* Create and execute threads */
	for(i = 0; i < Sp_NTHREAD; i++) {
		thread_id[i] = i;
		status = pthread_create(&threads[i], &attr, ThreadFunc, &thread_id[i]);
		Deb_ASSERT(status == 0); /* Temporary measure */
	}

	/* Join threads */
	for(i = 0; i < Sp_NTHREAD; i++) {
		status = pthread_join(threads[i], NULL);
		Deb_ASSERT(status == 0); /* Temporary measure */
	}

	return;
}

/*----------------------------------------------------------------------------*/

int SpUtil_Threads2(size_t nthread, void *(*ThreadFunc)(void *))
/* Distribute ThreadFunc to nthread threads, passing the id of each thread
 * as the argument */
{
	pthread_t threads[nthread], mon_thread;
	pthread_attr_t attr;
	int sts;
	size_t i, thread_id[nthread];
	sigset_t mask;

	/* Release the GIL */
	Py_BEGIN_ALLOW_THREADS

	/* Catch only SIGINT and SIGKILL */
	sigemptyset(&mask);
	sigaddset(&mask, SIGINT);
	pthread_sigmask(SIG_BLOCK, &mask, NULL);

	/* Reset termthread flag */
	SpUtil_termthread = 0;

	/* Init thread attribute to joinable */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	/* Create monitor thread */
	sts = pthread_create(&mon_thread, &attr, (void *(*)(void *))SpUtil_SigMonThread, &SpUtil_termthread);
	Deb_ASSERT(sts == 0);

	/* Create and execute threads */
	for(i = 0; i < nthread; i++) {
		thread_id[i] = i;
		sts = pthread_create(&threads[i], &attr, ThreadFunc, &thread_id[i]);
		Deb_ASSERT(sts == 0); /* Temporary measure */
	}

	/* Join threads */
	for(i = 0; i < nthread; i++) {
		sts = pthread_join(threads[i], NULL);
		Deb_ASSERT(sts == 0); /* Temporary measure */
	}

	/* All threads joined, terminate monitor thread */
	pthread_kill(mon_thread, SIGTERM);

	/* Unblock signals for main thread */
	pthread_sigmask(SIG_UNBLOCK, &mask, NULL);

	/* Destroy attribute */
	pthread_attr_destroy(&attr);

	/* Reacquire the GIL */
	Py_END_ALLOW_THREADS

	/* The only thing that could go wrong here should
	 * be interruption by the user */
	if(PyErr_CheckSignals())
		sts = 1;

	return sts;
}

/*----------------------------------------------------------------------------*/

static void *SpUtil_SigMonThread(void)
{
	sigset_t mask;
	struct timespec timeout;
	int sig;

	/* Catch SIGTERM and SIGKILL */
	sigemptyset(&mask);
	sigaddset(&mask, SIGINT);
	sigaddset(&mask, SIGTERM);

	/* Set timeout to 2.5ms */
	Mem_BZERO(&timeout);
	//timeout.tv_sec = 1;
	timeout.tv_nsec = (time_t)2.5e8;

	while(1) {
		sig = sigtimedwait(&mask, NULL, &timeout);
		if(sig == SIGINT) {
			printf("Signal %d caught in '%s'\n", sig, __FUNCTION__);
			PyErr_SetInterrupt();
			SpUtil_termthread = 1;
			break;
		}
		else if(sig == SIGTERM) {
			//debug
			//printf("SIGTERM caught, terminating '%s' thread       \n", __FUNCTION__);
			break;
		}
	}

	pthread_exit(NULL);
}

/*----------------------------------------------------------------------------*/

int SpUtil_TermThread(void)
{
	return SpUtil_termthread;
}

