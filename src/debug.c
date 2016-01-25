/*
 * Copyright (C) 2007 Eric Chung
 * e-mail: schung@asiaa.sinica.edu.tw
 *
 * some useful tool for debugging
 * 
 * History:
 * esc 23Jul07 genesis
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>
#include <debug.h>

/*----------------------------------------------------------------------------*/

void Deb_Pausef(const char *file, int line)
{
  char buffer[BUFSIZ] = "";

  fprintf(stderr, "%s:%d:press ENTER to continue...\n", file, line);
  fgets(buffer, BUFSIZ, stdin);

  return;
}

#if 0 /* esc 25Dec08: where is usleep() defined? */
/*----------------------------------------------------------------------------*/

void Deb_Sleep(double secs)
{
  usleep((int)(secs * 1.0e6));

  return;
}
#endif

/*----------------------------------------------------------------------------*/

void Deb_Fprintf(FILE *fp, const char *file, int line, const char *format, ...)
{
  va_list ap;

  fprintf(fp, "%s:%d: ", file, line);
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  
  return;
}
