#ifndef __ERROR_H__
#define __ERROR_H__

#include <stdio.h>
#include "data_structs.h"

typedef struct Error {
	char *file, *func, *mesg;
	int line;
} Error;

/* The error stack is just a linked list of messages */
extern LNode *Err_Stack;

Error *Err_Alloc(const char *file, int line, const char *func, const char *mesg);
void Err_Free(void *ptr);
int Err_SetString(const char *file, int line, const char *func, const char *format, ...);
#define Err_SETSTRING(...) \
	Err_SetString(__FILE__, __LINE__, __FUNCTION__, __VA_ARGS__)
int Err_Occurred(void);
void Err_Fprintf(FILE *fp, int debug);
void Err_Clear(void);


#endif
