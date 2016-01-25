#ifndef __DEBUG_H__
#define __DEBUG_H__

/* void Deb_Sleep(double secs); esc 25Dec08 */
void Deb_Fprintf(FILE *fp, const char *file, int line, const char *format, ...);
void Deb_Pausef(const char *file, int line);

#define Deb_PRINT(...) \
	Deb_Fprintf(stderr, __FILE__, __LINE__, __VA_ARGS__)

#define Deb_PAUSE() \
	Deb_Pausef(__FILE__, __LINE__)

#endif

/* Apparently assert() from assert.h does not work in Python extension
 * modules
 */
#define Deb_ASSERT(expr)\
	{if(!(expr)) {Deb_PRINT("Assertion failed\n"); exit(1);}}

