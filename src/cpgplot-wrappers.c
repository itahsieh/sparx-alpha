/* Wrappers and conveneice functions for cpgplot routines */

#include <assert.h>
#include "cpgplot-wrappers.h"

void Cpg_Circ(double xcent, double ycent, double radius)
/* void cpgcirc(float xcent, float ycent, float radius); */
{
	cpgcirc((float)xcent, (float)ycent, (float)radius);

	return;
}

/*----------------------------------------------------------------------------*/

void Cpg_Pt1(double xpt, double ypt, int symbol)
/* void cpgpt1(float xpt, float ypt, int symbol); */
{
	cpgpt1((float)xpt, (float)ypt, symbol);

	return;
}

/*----------------------------------------------------------------------------*/

void Cpg_Env(double xmin, double xmax, double ymin, double ymax, int just, int axis)
/* void cpgenv(float xmin, float xmax, float ymin, float ymax, int just, int axis); */
{
	cpgenv((float)xmin, (float)xmax, (float)ymin, (float)ymax, just, axis);

	return;
}

/*----------------------------------------------------------------------------*/

void Cpg_Arro(double x1, double y1, double x2, double y2)
/* void cpgarro(float x1, float y1, float x2, float y2); */
{
	cpgarro((float)x1, (float)y1, (float)x2, (float)y2);

	return;
}

/*----------------------------------------------------------------------------*/

void Cpg_Swin(double x1, double x2, double y1, double y2)
/* void cpgswin(float x1, float x2, float y1, float y2); */
{
	cpgswin((float)x1, (float)x2, (float)y1, (float)y2);

	return;
}

/*----------------------------------------------------------------------------*/

void Cpg_Box(const char *xopt, double xtick, int nxsub, const char *yopt, double ytick, int nysub)
/* void cpgbox(const char *xopt, float xtick, int nxsub, \
 const char *yopt, float ytick, int nysub);

Options (for parameters XOPT and YOPT):
 A : draw Axis (X axis is horizontal line Y=0, Y axis is vertical
     line X=0).
 B : draw bottom (X) or left (Y) edge of frame.
 C : draw top (X) or right (Y) edge of frame.
 G : draw Grid of vertical (X) or horizontal (Y) lines.
 I : Invert the tick marks; ie draw them outside the viewport
     instead of inside.
 L : label axis Logarithmically (see below).
 N : write Numeric labels in the conventional location below the
     viewport (X) or to the left of the viewport (Y).
 P : extend ("Project") major tick marks outside the box (ignored if
     option I is specified).
 M : write numeric labels in the unconventional location above the
     viewport (X) or to the right of the viewport (Y).
 T : draw major Tick marks at the major coordinate interval.
 S : draw minor tick marks (Subticks).
 V : orient numeric labels Vertically. This is only applicable to Y.
     The default is to write Y-labels parallel to the axis.
 1 : force decimal labelling, instead of automatic choice (see PGNUMB).
 2 : force exponential labelling, instead of automatic.

When PGENV calls PGBOX, it sets both XOPT and YOPT according to the value of its parameter
AXIS: -1: 'BC', 0: 'BCNST', 1: 'ABCNST', 2: 'ABCGNST'.

*/
{
	cpgbox(xopt, (float)xtick, nxsub, yopt, (float)ytick, nysub);

	return;
}

/*----------------------------------------------------------------------------*/

void Cpg_Reset(void)
{
	/* Erase scene */
	cpgeras();

	/* Redraw bounding box */
	Cpg_Box("BCNST", 0.0, 0, "BCNST", 0.0, 0);

	return;
}

/*----------------------------------------------------------------------------*/

void _Cpg_Vect(const float *a, const float *b, size_t idim, size_t jdim, size_t i1, 
	size_t i2, size_t j1, size_t j2, double c, int nc, const float *tr, double blank)
{
	cpgvect(a, b, idim, jdim, i1, i2, j1, j2, c, nc, tr, blank);

	return;
}

/*----------------------------------------------------------------------------*/

void Cpg_Sch(double ch)
{
	cpgsch(ch);

	return;
}










