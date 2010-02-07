/* vsurf.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"


/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
/*  .                                                             . */
/*  .                  copyright (c) 1998 by UCAR                 . */
/*  .                                                             . */
/*  .       University Corporation for Atmospheric Research       . */
/*  .                                                             . */
/*  .                      all rights reserved                    . */
/*  .                                                             . */
/*  .                                                             . */
/*  .                         SPHEREPACK                          . */
/*  .                                                             . */
/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */


/* Subroutine */ int vsurf_(real *xeye, real *yeye, real *zeye, integer *ntri,
	 real *x1, real *y1, real *z1, real *x2, real *y2, real *z2, real *x3,
	 real *y3, real *z3, integer *itype, real *work, integer *iwork)
{
    extern /* Subroutine */ int vsurf1_(real *, real *, real *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, integer *, real *, real *, real *, real *, real *, real *,
	     real *, real *, real *, real *, real *, real *, real *, real *, 
	    integer *, integer *, integer *, integer *);


/*    subroutine vsurf is like subroutine hidel except the triangles */
/*    are categorized. vsurf is also like solid except triangles rather */
/*    than lines are covered. */

/*     written by paul n. swarztrauber, national center for atmospheric */
/*     research, p.o. box 3000, boulder, colorado, 80307 */

/*    this program plots visible lines for the surface defined */
/*    by the input 3-d triangles with corners at (x1,y1,z1), (x2,y2,z2) */
/*    and (x3,y3,z3). the sides of these these triangles may or */
/*    may not be plotted depending on itype. if itype is 1 then the */
/*    side between points (x1,y1,z1) and (x2,y2,z2) is plotted if it */
/*    is visible. if itype is 2 then the side between (x2,y2,z2) */
/*    and (x3,y3,z3) is plotted. if itype is 3 then the visible portion */
/*    of the side between (x3,y3,z3) and (x1,y1,z1) is plotted. */
/*    any combination is possible by specifying itype to be one */
/*    of the following values: 0,1,2,3,12,13,23,123. */

/*    the length of real    array  work must be at least 19*ntri */

/*    the length of integer array iwork must be at least 19*ntri */


/*    the vertices of the triangles are renumbered by vsurf so that */
/*    their projections are orientated counterclockwise. the user need */
/*    only be aware that the vertices may be renumbered by vsurf. */


    /* Parameter adjustments */
    --iwork;
    --work;
    --itype;
    --z3;
    --y3;
    --x3;
    --z2;
    --y2;
    --x2;
    --z1;
    --y1;
    --x1;

    /* Function Body */
    vsurf1_(xeye, yeye, zeye, ntri, &x1[1], &y1[1], &z1[1], &x2[1], &y2[1], &
	    z2[1], &x3[1], &y3[1], &z3[1], &itype[1], &work[1], &work[*ntri + 
	    1], &work[(*ntri << 1) + 1], &work[*ntri * 3 + 1], &work[(*ntri <<
	     2) + 1], &work[*ntri * 5 + 1], &work[*ntri * 6 + 1], &work[*ntri 
	    * 7 + 1], &work[(*ntri << 3) + 1], &work[*ntri * 9 + 1], &work[*
	    ntri * 10 + 1], &work[*ntri * 11 + 1], &work[*ntri * 12 + 1], &
	    work[*ntri * 13 + 1], &iwork[*ntri * 14 + 1], &iwork[*ntri * 6 + 
	    1], &iwork[*ntri * 15 + 1], &iwork[*ntri * 17 + 1]);
    return 0;
} /* vsurf_ */

