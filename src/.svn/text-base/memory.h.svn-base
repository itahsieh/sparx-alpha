#ifndef __MEM_H__
#define __MEM_H__

#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#define Mem_BZERO(var)\
	memset((var), 0, sizeof(*(var)))

#define Mem_BZERO2(var, n)\
	memset((var), 0, sizeof(*(var)) * (n))

#define Mem_MEMSET(var, val, n)\
	memset((var), (val), sizeof(*(var)) * (n))

#define Mem_RESET(var, n)\
	Mem_MEMSET((var), 0, (n))

#define Mem_MEMCPY(a, b, n)\
	memcpy((a), (b), sizeof(*(a)) * (size_t)(n))

#define Mem_MALLOC(n)\
	Mem_Calloc(n, __FILE__, __LINE__)

#define Mem_CALLOC(a, b) \
	Mem_Calloc((size_t)(a) * sizeof(*(b)), __FILE__, __LINE__)

#define Mem_CALLOC2(a, b) \
	Mem_Calloc((size_t)(a) * (b), __FILE__, __LINE__)

#define Mem_REALLOC(a, b) \
	Mem_Realloc((b), (a) * sizeof(*(b)), __FILE__, __LINE__)

#define Mem_STRDUP(a) \
	Mem_Strdup((a), __FILE__, __LINE__)

#define Mem_FWRITE(ptr, len, fp)\
	fwrite((ptr), (size_t)len, sizeof(*(ptr)), fp)

#define Mem_FREAD(ptr, len, fp)\
	fread((ptr), (size_t)len, sizeof(*(ptr)), fp)

void *Mem_Calloc(size_t n, const char *file, int line);
void *Mem_Realloc(void *ptr, size_t n, const char *file, int line);
void *Mem_Strdup(const char *string, const char *file, int line);
char *Mem_VSprintf(const char *format, va_list ap);
char *Mem_Sprintf(const char *format, ...);
char *MemStr_Stripwhite(char *string);
char **MemStr_Split(const char *string, const char *delimiter, size_t *n);
void MemStr_FreeList(char **list, size_t n);
size_t Mem_FwriteStr(const char *str, FILE *fp);
size_t Mem_FreadStr(char *str, size_t len, FILE *fp);
size_t Mem_FwriteSize_t(size_t siz, FILE *fp);
size_t Mem_FreadSize_t(size_t *siz, FILE *fp);
size_t Mem_FwriteInt(int integer, FILE *fp);
size_t Mem_FreadInt(int *integer, FILE *fp);
char *Mem_FgetStr(FILE *fp);
void Mem_FputStr(const char *str, FILE *fp);

#endif
