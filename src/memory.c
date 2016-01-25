#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "memory.h"

/*----------------------------------------------------------------------------*/

void *Mem_Calloc(size_t n, const char *file, int line)
{
	void *ptr = malloc(n);

	if(!ptr) {
		fprintf(stderr, "Out of memory at %s:%d\n", file, line);
		exit(1);
	}
	else
		memset(ptr, 0, n);
// printf("Mem_Calloc ptr=%d\n",ptr);
	return ptr;
}

/*----------------------------------------------------------------------------*/

void *Mem_Realloc(void *ptr, size_t n, const char *file, int line)
{
  void *tmpptr = 0;

  if(!(tmpptr = realloc(ptr, n))) {
    fprintf(stderr, "Out of memory at %s:%d\n", file, line);
    exit(1);
  }
//   printf("memory.c:34 checkpoint ptr=%d\n",ptr);
  return tmpptr;
}

/*----------------------------------------------------------------------------*/

void *Mem_Strdup(const char *string, const char *file, int line)
{
	size_t n = strlen(string);
	char *ptr = malloc(n + 1 * sizeof(char));

	strncpy(ptr, string, n);
	ptr[n] = '\0';

	if(!ptr) {
		fprintf(stderr, "Out of memory at %s:%d\n", file, line);
		exit(1);
	}
// 	printf("Mem_Strdup ptr=%d\n",ptr);
	return ptr;
}

/*----------------------------------------------------------------------------*/

char *Mem_VSprintf(const char *format, va_list ap)
{
	size_t len;
// 	printf("memory.c:61 checkpoint\n");;
	char *string = Mem_CALLOC(BUFSIZ, string);
// 	printf("memory.c:63 checkpoint string=%d\n",string);

	len = (size_t)vsnprintf(string, (size_t)BUFSIZ, format, ap);

	if(len >= BUFSIZ) {
		string = Mem_REALLOC(len + 1, string);
		vsnprintf(string, len + 1, format, ap);
	}

	return string;
}

/*----------------------------------------------------------------------------*/

char *Mem_Sprintf(const char *format, ...)
{
	char *string;
	va_list ap;

	va_start(ap, format);
	string = Mem_VSprintf(format, ap);
	va_end(ap);

	return string;
}

/*----------------------------------------------------------------------------*/

char *MemStr_Stripwhite(char *string)
{
	char *s,*t;

	for (s = string; isspace (*s); s++)
		;

	if (*s == 0)
		return s;

	t = s + strlen(s) - 1;

	while (t > s && isspace (*t))
		t--;

	*++t = '\0';

	return s;
}

/*----------------------------------------------------------------------------*/

char **MemStr_Split(const char *string, const char *delimiter, size_t *n)
{
	char *buffer = NULL, *tok = NULL, **list = NULL;

	buffer = Mem_STRDUP(string);
	*n = 0;

	tok = strtok(buffer, delimiter);
	while(tok) {
		(*n)++;
		list = Mem_REALLOC((*n), list);
		list[(*n) - 1] = Mem_STRDUP(MemStr_Stripwhite(tok));
		tok = strtok(NULL, delimiter);
	}

	free(buffer);

	return list;
}

/*----------------------------------------------------------------------------*/

void MemStr_FreeList(char **list, size_t n)
{
	size_t i;

	for(i = 0; i < n; i++) {
		free(list[i]);
	}
	free(list);

	return;
}

/*----------------------------------------------------------------------------*/

size_t Mem_FwriteStr(const char *str, FILE *fp)
{
	size_t len = strlen(str), nbytes;

	nbytes = fwrite(str, len, sizeof(*str), fp);

	return nbytes;
}

/*----------------------------------------------------------------------------*/

size_t Mem_FreadStr(char *str, size_t len, FILE *fp)
{
	size_t nbytes;

	nbytes = fread(str, len, sizeof(*str), fp);

	return nbytes;
}

/*----------------------------------------------------------------------------*/

size_t Mem_FwriteInt(int integer, FILE *fp)
{
	size_t nbytes;

	nbytes = fwrite(&integer, (size_t)1, sizeof(integer), fp);

	return nbytes;
}

/*----------------------------------------------------------------------------*/

size_t Mem_FreadInt(int *integer, FILE *fp)
{
	size_t nbytes;

	nbytes = fread(integer, (size_t)1, sizeof(*integer), fp);

	return nbytes;
}

/*----------------------------------------------------------------------------*/

size_t Mem_FwriteSize_t(size_t siz, FILE *fp)
{
	size_t nbytes;

	nbytes = fwrite(&siz, (size_t)1, sizeof(siz), fp);

	return nbytes;
}

/*----------------------------------------------------------------------------*/

size_t Mem_FreadSize_t(size_t *siz, FILE *fp)
{
	size_t nbytes;

	nbytes = fread(siz, (size_t)1, sizeof(*siz), fp);

	return nbytes;
}

/*----------------------------------------------------------------------------*/

void Mem_FputStr(const char *str, FILE *fp)
{
	Mem_FwriteSize_t(strlen(str), fp);
	Mem_FwriteStr(str, fp);

	return;
}

/*----------------------------------------------------------------------------*/

char *Mem_FgetStr(FILE *fp)
{
	size_t len;
	char *buf;

	Mem_FreadSize_t(&len, fp);
	buf = Mem_CALLOC(len + 1, buf);
	Mem_FreadStr(buf, len, fp);

	return buf;
}











