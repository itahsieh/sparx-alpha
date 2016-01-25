#ifndef __CPGPLOT_WRAPPERS_H__
#define __CPGPLOT_WRAPPERS_H__

#include <stdlib.h>
#include "cpgplot.h"

void Cpg_Circ(double xcent, double ycent, double radius);
void Cpg_Pt1(double xpt, double ypt, int symbol);
void Cpg_Arro(double x1, double y1, double x2, double y2);
void Cpg_Env(double xmin, double xmax, double ymin, double ymax, int just, int axis);
void Cpg_Swin(double x1, double x2, double y1, double y2);
void Cpg_Box(const char *xopt, double xtick, int nxsub, const char *yopt, double ytick, int nysub);
void Cpg_Reset(void);
void _Cpg_Vect(const float *a, const float *b, size_t idim, size_t jdim, size_t i1, 
	size_t i2, size_t j1, size_t j2, double c, int nc, const float *tr, double blank);
#define Cpg_Vect(a, b, idim, jdim, i1, i2, j1, j2, c, nc, tr, blank)\
	_Cpg_Vect((a), (b), (size_t)(idim), (size_t)(jdim), (size_t)i1, (size_t)i2, (size_t)j1, (size_t)j2, (c), (nc), (tr), (blank))
void Cpg_Sch(double ch);



#endif
