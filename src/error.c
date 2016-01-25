#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>

#include "memory.h"
#include "error.h"

LNode *Err_stack = 0;

/*----------------------------------------------------------------------------*/

Error *Err_Alloc(const char *file, int line, const char *func, const char *mesg)
{
	Error *err = Mem_CALLOC(1, err);

	err->file = Mem_STRDUP(file);
	err->line = line;
	err->func = Mem_STRDUP(func);
	err->mesg = Mem_STRDUP(mesg);

	return err;
}

/*----------------------------------------------------------------------------*/

void Err_Free(void *ptr)
{
	Error *err = ptr;

	free(err->file);
	free(err->func);
	free(err->mesg);
	free(err);

	return;
}

/*----------------------------------------------------------------------------*/

int Err_SetString(const char *file, int line, const char *func, const char *format, ...)
{
	char string[BUFSIZ] = "";
	Error *err;
	va_list ap;

	/* Format message string */
	va_start(ap, format);
	vsnprintf(string, (size_t)BUFSIZ, format, ap);
	va_end(ap);

	/* Allocate error structure */
	err = Err_Alloc(file, line, func, string);

	/* Push onto error stack */
	Err_stack = Dat_Llst_Push(Err_stack, err, Err_Free); 

	/* Return a positive integer, which can be used by the caller
	 * to indicate an error condition */
	return 1;
}

/*----------------------------------------------------------------------------*/

int Err_Occurred(void)
{
	return Err_stack ? 1 : 0;
}

/*----------------------------------------------------------------------------*/

void Err_Fprintf(FILE *fp, int debug)
{
	LNode *end, *np;
	Error *err;

	assert(Err_stack != NULL);

	/* Seek to end of list */
	for(end = Err_stack; end; end = end->next)
		if(!end->next) break;

	/* Print out error messages from bottom of stack */
	fprintf(fp, "Error traceback (most recent last):\n");
	for(np = end; np; np = np->prev) {
		err = np->value;

		if(debug) {
			fprintf(fp, "  %s:%d: In function `%s':\n", err->file, err->line, err->func);
			fprintf(fp, "%s\n", err->mesg);
		}
		else {
			fprintf(fp, "  %s\n", err->mesg);
		}
	}

	return;
}

/*----------------------------------------------------------------------------*/

void Err_Clear(void)
{
	if(Err_stack) {
		Dat_Llst_Free(Err_stack);
		Err_stack = NULL;
	}

	return;
}










