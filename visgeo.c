/* visgeo.f -- translated by f2c (version 20061008).
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

/* Table of constant values */

static integer c__2 = 2;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b32 = 0.;
static doublereal c_b33 = 1.;


/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
/*  .                                                             . */
/*  .                  copyright (c) 1998 by UCAR                 . */
/*  .                                                             . */
/*  .       University Corporation for Atmospheric Research       . */
/*  .                                                             . */
/*  .                      all rights reserved                    . */
/*  .                                                             . */
/*  .                                                             . */
/*  .                         SPHEREPACK                       . */
/*  .                                                             . */
/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* ... visgeo.f */

/*     contains documentation and code for subroutine visgeo */

/* Subroutine */ int visgeo_(integer *m, integer *idp, integer *jdp, doublereal *x, 
	doublereal *y, doublereal *z__, doublereal *h__, doublereal *eyer, doublereal *eyelat, doublereal *eyelon,
	 doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer 
	*ierror)
{
    /* System generated locals */
    integer h_dim1, h_dim2, h_offset, x_dim1, x_dim2, x_offset, y_dim1, 
	    y_dim2, y_offset, z_dim1, z_dim2, z_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, 
	    i12, i13, i14, lg, lt, mmsq;
    extern /* Subroutine */ int visgeo1_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *);


/*     subroutine visgeo will display a function on the sphere */
/*     as a solid. ie. as a "lumpy" sphere. visgeo calls subroutine */
/*     vsurf to produce the visible surface rendering. X, Y, and Z */
/*     are the points on an icosahedral geodesic computed by */
/*     subroutine geopts available in spherepack. */

/*     requires routines visgeo1 ctos stoc vsurf vsurf1 */
/*                       prjct box */

/*     visgeo uses the ncar graphics package. */
/*     compile with: ncargf77 (all programs above) */

/*     execute with:  a.out */

/*     on screen display with:  ctrans -d x11 gmeta */

/*     print with:  ctrans -d ps.color gmeta > gmeta.ps */
/*                  lpr -P(your printer) gmeta.ps */


/*     input parameters */

/*     m        the number of points on one edge of the icosahedron */

/*     idp,jdp  the first and second dimensions of the three */
/*              dimensional arrays x, y, z, and h. */

/*     x,y,z    the coordinates of the geodesic points on */
/*              the unit sphere computed by subroutine geopts. */
/*              the indices are defined on the unfolded */
/*              icosahedron as follows for the case m=3 */

/*                north pole */

/*                 (5,1)          0      l */
/*        i     (4,1) (5,2)              a    (repeated for */
/*           (3,1) (4,2) (5,3)  theta1   t    k=2,3,4,5 in */
/*        (2,1) (3,2) (4,3)              i        --> */
/*     (1,1) (2,2) (3,3)        theta2   t    the longitudinal */
/*        (1,2) (2,3)                    u    direction) */
/*           (1,3)                pi     d */
/*      j                                e */
/*         south pole */

/*            total number of vertices is  10*(m-1)**2+2 */
/*            total number of triangles is 20*(m-1)**2 */

/*     h      a three dimensional array that contains the discrete */
/*            function to be displayed. h(i,j,k) is the distance from */
/*            the center of the sphere to the "lumpy" surface at the */
/*             point [x(i,j,k),y(i,j,k),z(i,j,k)] on the unit sphere. */

/*     eyer   the distance from the center of the sphere to the eye. */

/*     eyelat the colatitudinal coordinate of the eye (in degrees). */

/*     eyelon the longitudinal  coordinate of the eye (in degrees). */

/*     idp    the first dimension of the array h as it appears in */
/*            the program that calls visgeo */

/*     jdp    the second dimension of the array h as it appears in */
/*            the program that calls visgeo */

/*     work   a doublereal work array */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls visgeo. lwork must be at least */
/*                       480*(m-1)**2. */

/*     iwork  an integer work array */

/*     liwork the dimension of the array iwork as it appears in the */
/*            program that calls visgeo. liwork must be at least */
/*                       140*(m-1)**2. */

/*     input parameter */

/*     ierror = 0    no error */
/*            = 1    h(i,j,k) is less than zero for some i,j,k. */
/*            = 2    eyer is less than h(i,j,k) for some i,k,k. */
/*            = 3    lwork  is less than 480*(m-1)**2 */
/*            = 4    liwork is less than 140*(m-1)**2 */

    /* Parameter adjustments */
    h_dim1 = *idp;
    h_dim2 = *jdp;
    h_offset = 1 + h_dim1 * (1 + h_dim2);
    h__ -= h_offset;
    z_dim1 = *idp;
    z_dim2 = *jdp;
    z_offset = 1 + z_dim1 * (1 + z_dim2);
    z__ -= z_offset;
    y_dim1 = *idp;
    y_dim2 = *jdp;
    y_offset = 1 + y_dim1 * (1 + y_dim2);
    y -= y_offset;
    x_dim1 = *idp;
    x_dim2 = *jdp;
    x_offset = 1 + x_dim1 * (1 + x_dim2);
    x -= x_offset;
    --work;
    --iwork;

    /* Function Body */
/* Computing 2nd power */
    i__1 = *m - 1;
    mmsq = i__1 * i__1;
    *ierror = 3;
    if (*lwork < mmsq * 480) {
	return 0;
    }
    *ierror = 4;
    if (*liwork < mmsq * 140) {
	return 0;
    }
    for (k = 1; k <= 5; ++k) {
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m + *m - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		if (h__[i__ + (j + k * h_dim2) * h_dim1] >= 0.) {
		    goto L15;
		}
		*ierror = 1;
		return 0;
L15:
		if (*eyer > h__[i__ + (j + k * h_dim2) * h_dim1]) {
		    goto L10;
		}
		*ierror = 2;
		return 0;
L10:
		;
	    }
	}
    }
    *ierror = 0;
/* Computing 2nd power */
    i__2 = *m - 1;
    lt = i__2 * i__2 * 20;
    lg = *m * 5 * (*m + *m - 1);
    i1 = 1;
    i2 = i1 + lt;
    i3 = i2 + lt;
    i4 = i3 + lt;
    i5 = i4 + lt;
    i6 = i5 + lt;
    i7 = i6 + lt;
    i8 = i7 + lt;
    i9 = i8 + lt;
    i10 = i9 + lt;
    i11 = i10 + lt;
    i12 = i11;
    i13 = i12 + lg;
    i14 = i13 + lg;
    visgeo1_(m, idp, jdp, &h__[h_offset], eyer, eyelat, eyelon, &x[x_offset], 
	    &y[y_offset], &z__[z_offset], &work[i1], &work[i2], &work[i3], &
	    work[i4], &work[i5], &work[i6], &work[i7], &work[i8], &work[i9], &
	    iwork[1], &work[i11], &work[i12], &work[i13], &work[i14], &iwork[
	    lt + 1]);
    return 0;
} /* visgeo_ */

/* Subroutine */ int visgeo1_(integer *m, integer *idp, integer *jdp, doublereal *
	h__, doublereal *eyer, doublereal *eyelat, doublereal *eyelon, doublereal *xi, doublereal *yi, doublereal 
	*zi, doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, doublereal *z2, doublereal 
	*x3, doublereal *y3, doublereal *z3, integer *itype, doublereal *work, doublereal *x, doublereal *y,
	 doublereal *z__, integer *iwork)
{
    /* System generated locals */
    integer h_dim1, h_dim2, h_offset, xi_dim1, xi_dim2, xi_offset, yi_dim1, 
	    yi_dim2, yi_offset, zi_dim1, zi_dim2, zi_offset, x_dim1, x_dim2, 
	    x_offset, y_dim1, y_dim2, y_offset, z_dim1, z_dim2, z_offset, 
	    i__1, i__2;

    /* Builtin functions */
    double atan(doublereal), sin(doublereal), cos(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal pi, rad, dtr;
    extern /* Subroutine */ int ctos_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), stoc_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal xeye, yeye;
    static integer ntri;
    static doublereal zeye, theta;
    extern /* Subroutine */ int vsurf_(doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal elambda;


/*     the * above refers to 20*(m-1)**2 locations which is the */
/*     number of triangles */

    /* Parameter adjustments */
    z_dim1 = *m + *m - 1;
    z_dim2 = *m;
    z_offset = 1 + z_dim1 * (1 + z_dim2);
    z__ -= z_offset;
    y_dim1 = *m + *m - 1;
    y_dim2 = *m;
    y_offset = 1 + y_dim1 * (1 + y_dim2);
    y -= y_offset;
    x_dim1 = *m + *m - 1;
    x_dim2 = *m;
    x_offset = 1 + x_dim1 * (1 + x_dim2);
    x -= x_offset;
    zi_dim1 = *idp;
    zi_dim2 = *jdp;
    zi_offset = 1 + zi_dim1 * (1 + zi_dim2);
    zi -= zi_offset;
    yi_dim1 = *idp;
    yi_dim2 = *jdp;
    yi_offset = 1 + yi_dim1 * (1 + yi_dim2);
    yi -= yi_offset;
    xi_dim1 = *idp;
    xi_dim2 = *jdp;
    xi_offset = 1 + xi_dim1 * (1 + xi_dim2);
    xi -= xi_offset;
    h_dim1 = *idp;
    h_dim2 = *jdp;
    h_offset = 1 + h_dim1 * (1 + h_dim2);
    h__ -= h_offset;
    --x1;
    --y1;
    --z1;
    --x2;
    --y2;
    --z2;
    --x3;
    --y3;
    --z3;
    --itype;
    --work;
    --iwork;

    /* Function Body */
    for (k = 1; k <= 5; ++k) {
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m + *m - 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ctos_(&xi[i__ + (j + k * xi_dim2) * xi_dim1], &yi[i__ + (j + 
			k * yi_dim2) * yi_dim1], &zi[i__ + (j + k * zi_dim2) *
			 zi_dim1], &rad, &theta, &elambda);
		stoc_(&h__[i__ + (j + k * h_dim2) * h_dim1], &theta, &elambda,
			 &x[i__ + (j + k * x_dim2) * x_dim1], &y[i__ + (j + k 
			* y_dim2) * y_dim1], &z__[i__ + (j + k * z_dim2) * 
			z_dim1]);
/* L10: */
	    }
	}
    }
    ntri = 0;
    for (k = 1; k <= 5; ++k) {
	i__2 = *m - 1;
	for (j = 1; j <= i__2; ++j) {
	    i__1 = *m + *m - 2;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		++ntri;
		x1[ntri] = x[i__ + (j + k * x_dim2) * x_dim1];
		y1[ntri] = y[i__ + (j + k * y_dim2) * y_dim1];
		z1[ntri] = z__[i__ + (j + k * z_dim2) * z_dim1];
		x2[ntri] = x[i__ + 1 + (j + 1 + k * x_dim2) * x_dim1];
		y2[ntri] = y[i__ + 1 + (j + 1 + k * y_dim2) * y_dim1];
		z2[ntri] = z__[i__ + 1 + (j + 1 + k * z_dim2) * z_dim1];
		x3[ntri] = x[i__ + 1 + (j + k * x_dim2) * x_dim1];
		y3[ntri] = y[i__ + 1 + (j + k * y_dim2) * y_dim1];
		z3[ntri] = z__[i__ + 1 + (j + k * z_dim2) * z_dim1];
		itype[ntri] = 13;
		++ntri;
		x1[ntri] = x[i__ + (j + k * x_dim2) * x_dim1];
		y1[ntri] = y[i__ + (j + k * y_dim2) * y_dim1];
		z1[ntri] = z__[i__ + (j + k * z_dim2) * z_dim1];
		x2[ntri] = x[i__ + 1 + (j + 1 + k * x_dim2) * x_dim1];
		y2[ntri] = y[i__ + 1 + (j + 1 + k * y_dim2) * y_dim1];
		z2[ntri] = z__[i__ + 1 + (j + 1 + k * z_dim2) * z_dim1];
		x3[ntri] = x[i__ + (j + 1 + k * x_dim2) * x_dim1];
		y3[ntri] = y[i__ + (j + 1 + k * y_dim2) * y_dim1];
		z3[ntri] = z__[i__ + (j + 1 + k * z_dim2) * z_dim1];
		itype[ntri] = 3;
/* L20: */
	    }
	}
    }
/*      write(6,22) ntri */
/* L22: */
/*      write(6,23) (x1(l2),y1(l2),z1(l2),x2(l2),y2(l2),z2(l2), */
/*     1             x3(l2),y3(l2),z3(l2),l2=1,ntri) */
/* 23   format(9f10.7) */

    pi = atan(1.f) * 4.;
    dtr = pi / 180.;
    xeye = *eyer * sin(dtr * *eyelat);
    yeye = xeye * sin(dtr * *eyelon);
    xeye *= cos(dtr * *eyelon);
    zeye = *eyer * cos(dtr * *eyelat);
    vsurf_(&xeye, &yeye, &zeye, &ntri, &x1[1], &y1[1], &z1[1], &x2[1], &y2[1],
	     &z2[1], &x3[1], &y3[1], &z3[1], &itype[1], &work[1], &iwork[1]);
    return 0;
} /* visgeo1_ */

/* Subroutine */ int ctos_(doublereal *x, doublereal *y, doublereal *z__, doublereal *r__, doublereal *
	theta, doublereal *phi)
{
    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), atan2(doublereal, doublereal);

    /* Local variables */
    static doublereal r1;

    r1 = *x * *x + *y * *y;
    if (r1 != 0.) {
	goto L10;
    }
    *phi = 0.;
    *theta = 0.;
    if (*z__ < 0.) {
	*theta = atan(1.f) * 4.;
    }
    return 0;
L10:
    *r__ = sqrt(r1 + *z__ * *z__);
    r1 = sqrt(r1);
    *phi = atan2(*y, *x);
    *theta = atan2(r1, *z__);
    return 0;
} /* ctos_ */

/* Subroutine */ int stoc_(doublereal *r__, doublereal *theta, doublereal *phi, doublereal *x, doublereal *
	y, doublereal *z__)
{
    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static doublereal st;

    st = sin(*theta);
    *x = *r__ * st * cos(*phi);
    *y = *r__ * st * sin(*phi);
    *z__ = *r__ * cos(*theta);
    return 0;
} /* stoc_ */

/* Subroutine */ int vsurf_(doublereal *xeye, doublereal *yeye, doublereal *zeye, integer *ntri,
	 doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, doublereal *z2, doublereal *x3,
	 doublereal *y3, doublereal *z3, integer *itype, doublereal *work, integer *iwork)
{
    extern /* Subroutine */ int vsurf1_(doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, 
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

/*    the length of doublereal    array  work must be at least 14*ntri */

/*    the length of integer array iwork must be at least  6*ntri */


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
	    work[*ntri * 13 + 1], &iwork[1], &iwork[*ntri + 1], &iwork[(*ntri 
	    << 1) + 1], &iwork[(*ntri << 2) + 1]);
    return 0;
} /* vsurf_ */

/* Subroutine */ int vsurf1_(doublereal *xeye, doublereal *yeye, doublereal *zeye, integer *
	ntri, doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, doublereal *z2, 
	doublereal *x3, doublereal *y3, doublereal *z3, integer *itype, doublereal *px1, doublereal *py1, 
	doublereal *px2, doublereal *py2, doublereal *px3, doublereal *py3, doublereal *vx1, doublereal *vy1, 
	doublereal *vx2, doublereal *vy2, doublereal *vx3, doublereal *vy3, doublereal *tl, doublereal *tr, 
	integer *kh, integer *next, integer *istart, integer *ifinal)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal r__1, r__2;

    /* Builtin functions */
    double log(doublereal);
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static doublereal a, b, c__, d__;
    static integer i__, k, l;
    static doublereal x, y, z__;
    static integer i2, j2, j1, k1, k2;
    static doublereal x4, y4, z4, x5, y5, z5, c14, c25, c17, c27;
    static integer id, jd, if__, jf, kd, lf;
    static doublereal c37;
    static integer il;
    static doublereal c36;
    static integer kb;
    static doublereal xa, hr;
    static integer ir;
    static doublereal ya, xb;
    static integer ks, ls, ns;
    static doublereal x54, y54, yb;
    static integer id1, id2, i1f, j1f, id3;
    static doublereal l2e;
    static doublereal le2;
    static integer i2m, j2m, ip2[11], ir1, ir2, i1s, j1s;
    static doublereal tl1, tl2, px4, py4, px5, py5;
    static integer kdf, ijd;
    static doublereal den, apl, bpl;
    static integer ird[11], isd;
    static doublereal hgr;
    static integer icv, nct[11];
    static doublereal hdy;
    static integer ncv[11];
    static doublereal hdx, dmx, dmy;
    static integer ifx;
    extern /* Subroutine */ int set_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *), box_(integer *, integer *, 
	    integer *, integer *, integer *);
    static integer kcv, ltp[3];
    static doublereal dxt;
    static integer ist, ixh, isx, ixs, ixf, kds, ldo, ith;
    static doublereal xpl, ypl, zpl, xpr, ypr, zpr, cpl, dpl, til, tir, tim, tlh, 
	    trh, dum1, dum2;
    static integer irp1;
    static doublereal px1h, py1h, vx1t, vy1t, vz1t, vx2t, vy2t, vz2t, vx3t, vy3t, 
	    vz3t;
    extern /* Subroutine */ int line_(doublereal *, doublereal *, doublereal *, doublereal *);
    static integer nseg, irdp;
    static doublereal xmid, ymid;
    static integer last[11];
    static doublereal pmin, pmax, tmin;
    static integer maxs;
    static doublereal tmax, xmin;
    static integer isxm;
    static doublereal xmax, ymin, ymax;
    static integer ityp, irmp1;
    extern /* Subroutine */ int frame_(void);
    static doublereal cprod, thold;
    static integer nsegp, irmax;
    static doublereal fntri;
    static integer isize;
    extern /* Subroutine */ int prjct_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal x1hold, y1hold, z1hold;




/*     compute projections of 3-d points */

    /* Parameter adjustments */
    --ifinal;
    --istart;
    --next;
    --kh;
    --tr;
    --tl;
    --vy3;
    --vx3;
    --vy2;
    --vx2;
    --vy1;
    --vx1;
    --py3;
    --px3;
    --py2;
    --px2;
    --py1;
    --px1;
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
    le2 = .6931471805599453094172321;
    l2e = 1. / le2;
    fntri = (doublereal) (*ntri);
    irmax = l2e * .5 * log(fntri);
    irmax = min(irmax,10);
    irmp1 = irmax + 1;
    for (icv = 1; icv <= 11; ++icv) {
	ncv[icv - 1] = 0;
/* L4: */
    }
    nct[0] = 0;
    ip2[0] = 1;
    ird[0] = 0;
    isize = 4;
    i__1 = irmp1;
    for (irp1 = 2; irp1 <= i__1; ++irp1) {
	ir = irp1 - 1;
	nct[irp1 - 1] = 0;
	ip2[irp1 - 1] = pow_ii(&c__2, &ir);
	ird[irp1 - 1] = ird[ir - 1] + isize;
/* Computing 2nd power */
	i__2 = ip2[irp1 - 1] + 1;
	isize = i__2 * i__2;
/* L7: */
    }
    isxm = ird[irmp1 - 1] + isize + 1;
    i__1 = isxm;
    for (isx = 1; isx <= i__1; ++isx) {
	istart[isx] = 0;
	ifinal[isx] = 0;
/* L8: */
    }
    i__1 = *ntri;
    for (i__ = 1; i__ <= i__1; ++i__) {
	next[i__] = 0;
/* L6: */
    }
    prjct_(&c__0, xeye, yeye, zeye, &x, &y, &z__, &dum1, &dum2);
/*      write(6,127) ntri */
/* L127: */
    i__1 = *ntri;
    for (k = 1; k <= i__1; ++k) {
	prjct_(&c__1, xeye, yeye, zeye, &x1[k], &y1[k], &z1[k], &px1[k], &py1[
		k]);
	prjct_(&c__1, xeye, yeye, zeye, &x2[k], &y2[k], &z2[k], &px2[k], &py2[
		k]);
	prjct_(&c__1, xeye, yeye, zeye, &x3[k], &y3[k], &z3[k], &px3[k], &py3[
		k]);
	if (k < 3) {
/*          write(6,333) xeye,yeye,zeye,x1(k),y1(k),z1(k),px1(k),py1(k) */
/* L333: */
	}
/* L86: */
    }

/*     orientate triangles counter clockwise */

    i__1 = *ntri;
    for (k = 1; k <= i__1; ++k) {
	cprod = (px2[k] - px1[k]) * (py3[k] - py1[k]) - (py2[k] - py1[k]) * (
		px3[k] - px1[k]);
/*      if(cprod.eq.0.) write(6,79) k,px1(k),px2(k),px3(k), */
/*     -                              py1(k),py2(k),py3(k) */
/* L79: */
	if (cprod >= 0.) {
	    goto L70;
	}
	px1h = px1[k];
	py1h = py1[k];
	px1[k] = px2[k];
	py1[k] = py2[k];
	px2[k] = px1h;
	py2[k] = py1h;
	x1hold = x1[k];
	y1hold = y1[k];
	z1hold = z1[k];
	x1[k] = x2[k];
	y1[k] = y2[k];
	z1[k] = z2[k];
	x2[k] = x1hold;
	y2[k] = y1hold;
	z2[k] = z1hold;
	ityp = itype[k];
	if (ityp == 2) {
	    itype[k] = 3;
	}
	if (ityp == 3) {
	    itype[k] = 2;
	}
	if (ityp == 12) {
	    itype[k] = 13;
	}
	if (ityp == 13) {
	    itype[k] = 12;
	}
L70:
	;
    }

/*     set screen limits */

    pmax = px1[1];
    pmin = px1[1];
    i__1 = *ntri;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	r__1 = pmin, r__2 = px1[k], r__1 = min(r__1,r__2), r__2 = py1[k], 
		r__1 = min(r__1,r__2), r__2 = px2[k], r__1 = min(r__1,r__2), 
		r__2 = py2[k], r__1 = min(r__1,r__2), r__2 = px3[k], r__1 = 
		min(r__1,r__2), r__2 = py3[k];
	pmin = dmin(r__1,r__2);
/* Computing MAX */
	r__1 = pmax, r__2 = px1[k], r__1 = max(r__1,r__2), r__2 = py1[k], 
		r__1 = max(r__1,r__2), r__2 = px2[k], r__1 = max(r__1,r__2), 
		r__2 = py2[k], r__1 = max(r__1,r__2), r__2 = px3[k], r__1 = 
		max(r__1,r__2), r__2 = py3[k];
	pmax = dmax(r__1,r__2);
/* L87: */
    }
    pmin *= 1.1;
    pmax *= 1.1;
    set_(&c_b32, &c_b33, &c_b32, &c_b33, &pmin, &pmax, &pmin, &pmax, &c__1);
/* Computing MIN */
    r__1 = min(px1[1],px2[1]);
    xmin = dmin(r__1,px3[1]);
/* Computing MAX */
    r__1 = max(px1[1],px2[1]);
    xmax = dmax(r__1,px3[1]);
/* Computing MIN */
    r__1 = min(py1[1],py2[1]);
    ymin = dmin(r__1,py3[1]);
/* Computing MAX */
    r__1 = max(py1[1],py2[1]);
    ymax = dmax(r__1,py3[1]);
    i__1 = *ntri;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MIN */
	r__1 = xmin, r__2 = px1[i__], r__1 = min(r__1,r__2), r__2 = px2[i__], 
		r__1 = min(r__1,r__2), r__2 = px3[i__];
	xmin = dmin(r__1,r__2);
/* Computing MAX */
	r__1 = xmax, r__2 = px1[i__], r__1 = max(r__1,r__2), r__2 = px2[i__], 
		r__1 = max(r__1,r__2), r__2 = px3[i__];
	xmax = dmax(r__1,r__2);
/* Computing MIN */
	r__1 = ymin, r__2 = py1[i__], r__1 = min(r__1,r__2), r__2 = py2[i__], 
		r__1 = min(r__1,r__2), r__2 = py3[i__];
	ymin = dmin(r__1,r__2);
/* Computing MAX */
	r__1 = ymax, r__2 = py1[i__], r__1 = max(r__1,r__2), r__2 = py2[i__], 
		r__1 = max(r__1,r__2), r__2 = py3[i__];
	ymax = dmax(r__1,r__2);
/* L1: */
    }
    dmx = xmax - xmin;
    dmy = ymax - ymin;
    if (dmx > dmy) {
	goto L2;
    }
    c__ = ymin;
    d__ = ymax;
    xmid = (xmin + xmax) * .5;
    hdy = dmy * .5;
    a = xmid - hdy;
    b = xmid + hdy;
    goto L3;
L2:
    a = xmin;
    b = xmax;
    ymid = (ymin + ymax) * .5;
    hdx = dmx * .5;
    c__ = ymid - hdx;
    d__ = ymid + hdx;
L3:
    hgr = b - a;

/*     categorize triangles */

    i__1 = *ntri;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	r__1 = px1[i__], r__2 = px2[i__], r__1 = min(r__1,r__2), r__2 = px3[
		i__];
	xmin = dmin(r__1,r__2);
/* Computing MAX */
	r__1 = px1[i__], r__2 = px2[i__], r__1 = max(r__1,r__2), r__2 = px3[
		i__];
	xmax = dmax(r__1,r__2);
/* Computing MIN */
	r__1 = py1[i__], r__2 = py2[i__], r__1 = min(r__1,r__2), r__2 = py3[
		i__];
	ymin = dmin(r__1,r__2);
/* Computing MAX */
	r__1 = py1[i__], r__2 = py2[i__], r__1 = max(r__1,r__2), r__2 = py3[
		i__];
	ymax = dmax(r__1,r__2);
/* Computing MAX */
	r__1 = xmax - xmin, r__2 = ymax - ymin;
	dxt = dmax(r__1,r__2);
	if (dxt > 0.) {
	    goto L10;
	}
	ir = irmax;
	goto L20;
L10:
	ir = l2e * log(hgr / dxt);
	ir = min(ir,irmax);
L20:
	irp1 = ir + 1;
	++nct[irp1 - 1];
	hr = hgr / ip2[irp1 - 1];
	xmid = (xmin + xmax) * .5;
	id = (xmid - a) / hr + 1.5;
	ymid = (ymin + ymax) * .5;
	jd = (ymid - c__) / hr + 1.5;
	ijd = ip2[irp1 - 1] + 1;
	isx = id + (jd - 1) * ijd + ird[irp1 - 1];
	ifx = ifinal[isx];
	if (ifx > 0) {
	    goto L50;
	}
	istart[isx] = i__;
	goto L60;
L50:
	next[ifx] = i__;
L60:
	ifinal[isx] = i__;
/* L100: */
    }
/*      write(6,106) tcat,(irp1,nct(irp1),irp1=1,irmp1) */
/* L106: */

/*     sort triangles into boxes */

    l = 0;
    i__1 = irmp1;
    for (irp1 = 1; irp1 <= i__1; ++irp1) {
	if (nct[irp1 - 1] == 0) {
	    goto L30;
	}
	ist = ird[irp1 - 1] + 1;
	isd = ip2[irp1 - 1] + 1;
	box_(&isd, &istart[ist], &next[1], &l, &ifinal[1]);
	last[irp1 - 1] = l + 1;
L30:
	;
    }
    i__1 = irmp1;
    for (irp1 = 1; irp1 <= i__1; ++irp1) {
/* Computing 2nd power */
	i__2 = ip2[irp1 - 1] + 1;
	il = ird[irp1 - 1] + i__2 * i__2 + 1;
	if (istart[il] == 0) {
	    istart[il] = last[irp1 - 1];
	}
/* L35: */
    }
/*      write(6,31) tsort,l,ntri */
/* L31: */
    i__1 = *ntri;
    for (k = 1; k <= i__1; ++k) {
	vx1[k] = px2[k] - px1[k];
	vy1[k] = py2[k] - py1[k];
	vx2[k] = px3[k] - px2[k];
	vy2[k] = py3[k] - py2[k];
	vx3[k] = px1[k] - px3[k];
	vy3[k] = py1[k] - py3[k];
/* L90: */
    }
    tl1 = 0.;
    tl2 = 0.;
    maxs = 0;
    i__1 = irmp1;
    for (ir2 = 1; ir2 <= i__1; ++ir2) {
	if (nct[ir2 - 1] == 0) {
	    goto L500;
	}
	ist = ird[ir2 - 1];
	isd = ip2[ir2 - 1] + 1;
	i__2 = isd;
	for (j2 = 1; j2 <= i__2; ++j2) {
	    i__3 = isd;
	    for (i2 = 1; i2 <= i__3; ++i2) {
		++ist;
		ls = istart[ist];
		lf = istart[ist + 1] - 1;
		if (lf < ls) {
		    goto L480;
		}

/*     define coverings */

		kcv = 0;
		i2m = i2 - 1;
		j2m = j2 - 1;
		i__4 = irmp1;
		for (ir1 = 1; ir1 <= i__4; ++ir1) {
		    if (nct[ir1 - 1] == 0) {
			goto L300;
		    }
		    if (ir1 >= ir2) {
			goto L260;
		    }
		    i__5 = ir2 - ir1;
		    irdp = pow_ii(&c__2, &i__5);
		    i1s = (i2m - 1) / irdp;
		    i1f = (i2m + 1) / irdp;
		    if__ = i2m + 1 - i1f * irdp;
		    if (if__ > 0) {
			++i1f;
		    }
		    j1s = (j2m - 1) / irdp;
		    j1f = (j2m + 1) / irdp;
		    jf = j2m + 1 - j1f * irdp;
		    if (jf > 0) {
			++j1f;
		    }
		    goto L270;
L260:
		    i__5 = ir1 - ir2;
		    irdp = pow_ii(&c__2, &i__5);
		    i1s = irdp * (i2m - 1);
		    i1f = irdp * (i2m + 1);
		    j1s = irdp * (j2m - 1);
		    j1f = irdp * (j2m + 1);
L270:
		    ijd = ip2[ir1 - 1] + 1;
/* Computing MAX */
		    i__5 = i1s + 1;
		    i1s = max(i__5,1);
/* Computing MIN */
		    i__5 = i1f + 1;
		    i1f = min(i__5,ijd);
/* Computing MAX */
		    i__5 = j1s + 1;
		    j1s = max(i__5,1);
/* Computing MIN */
		    i__5 = j1f + 1;
		    j1f = min(i__5,ijd);
		    ixh = (j1s - 2) * ijd + ird[ir1 - 1];
		    ixs = i1s + ixh;
		    ixf = i1f + ixh;
		    i__5 = j1f;
		    for (j1 = j1s; j1 <= i__5; ++j1) {
			ixs += ijd;
			kds = istart[ixs];
			ixf += ijd;
			kdf = istart[ixf + 1] - 1;
			if (kdf < kds) {
			    goto L290;
			}
			i__6 = kdf;
			for (kd = kds; kd <= i__6; ++kd) {
			    ++kcv;
			    kh[kcv] = ifinal[kd];
/* L280: */
			}
L290:
			;
		    }
L300:
		    ;
		}
		for (icv = 1; icv <= 10; ++icv) {
		    if (kcv <= ncv[icv - 1]) {
			goto L310;
		    }
		    ncv[icv - 1] = kcv;
		    goto L320;
L310:
		    ;
		}


L320:
		i__4 = lf;
		for (ldo = ls; ldo <= i__4; ++ldo) {
		    l = ifinal[ldo];
		    ith = itype[l];
		    if (ith == 0) {
			goto L470;
		    }
		    ltp[0] = 0;
		    ltp[1] = 0;
		    ltp[2] = 0;
		    id1 = ith / 100;
		    ith -= id1 * 100;
		    id2 = ith / 10;
		    id3 = ith - id2 * 10;
		    if (id1 != 0) {
			ltp[id1 - 1] = 1;
		    }
		    if (id2 != 0) {
			ltp[id2 - 1] = 1;
		    }
		    if (id3 != 0) {
			ltp[id3 - 1] = 1;
		    }
/*     if((ith.eq.123) .or. (ith.eq.12) .or.(ith.eq.13)) ltp(1) = 1 */
/*     if((ith.eq.123) .or. (ith.eq.23) .or.(ith.eq.12)) ltp(2) = 1 */
/*     if((ith.eq.123) .or. (ith.eq.13) .or.(ith.eq.23)) ltp(3) = 1 */
		    for (ns = 1; ns <= 3; ++ns) {
			switch (ns) {
			    case 1:  goto L101;
			    case 2:  goto L102;
			    case 3:  goto L103;
			}
L101:
			if (ltp[ns - 1] == 0) {
			    goto L460;
			}
			px4 = px1[l];
			py4 = py1[l];
			px5 = px2[l];
			py5 = py2[l];
			x4 = x1[l];
			y4 = y1[l];
			z4 = z1[l];
			x5 = x2[l];
			y5 = y2[l];
			z5 = z2[l];
			goto L105;
L102:
			if (ltp[ns - 1] == 0) {
			    goto L460;
			}
			px4 = px2[l];
			py4 = py2[l];
			px5 = px3[l];
			py5 = py3[l];
			x4 = x2[l];
			y4 = y2[l];
			z4 = z2[l];
			x5 = x3[l];
			y5 = y3[l];
			z5 = z3[l];
			goto L105;
L103:
			if (ltp[ns - 1] == 0) {
			    goto L460;
			}
			px4 = px1[l];
			py4 = py1[l];
			px5 = px3[l];
			py5 = py3[l];
			x4 = x1[l];
			y4 = y1[l];
			z4 = z1[l];
			x5 = x3[l];
			y5 = y3[l];
			z5 = z3[l];
L105:
			x54 = px5 - px4;
			y54 = py5 - py4;
			nseg = 0;
			i__5 = kcv;
			for (kd = 1; kd <= i__5; ++kd) {
			    k = kh[kd];
			    c17 = vx1[k] * y54 - vy1[k] * x54;
			    c27 = vx2[k] * y54 - vy2[k] * x54;
			    c37 = vx3[k] * y54 - vy3[k] * x54;
			    c14 = vy1[k] * (px4 - px1[k]) - vx1[k] * (py4 - 
				    py1[k]);
			    c25 = vy2[k] * (px4 - px2[k]) - vx2[k] * (py4 - 
				    py2[k]);
			    c36 = vy3[k] * (px4 - px3[k]) - vx3[k] * (py4 - 
				    py3[k]);
			    tmin = 0.;
			    tmax = 1.;
			    if (c17 < 0.) {
				goto L151;
			    } else if (c17 == 0) {
				goto L152;
			    } else {
				goto L153;
			    }
L151:
/* Computing MIN */
			    r__1 = c14 / c17;
			    tmax = dmin(r__1,tmax);
			    goto L154;
L152:
			    if (c14 >= 0.) {
				goto L440;
			    } else {
				goto L154;
			    }
L153:
/* Computing MAX */
			    r__1 = c14 / c17;
			    tmin = dmax(r__1,tmin);
L154:
			    if (c27 < 0.) {
				goto L155;
			    } else if (c27 == 0) {
				goto L156;
			    } else {
				goto L157;
			    }
L155:
/* Computing MIN */
			    r__1 = c25 / c27;
			    tmax = dmin(r__1,tmax);
			    goto L158;
L156:
			    if (c25 >= 0.) {
				goto L440;
			    } else {
				goto L158;
			    }
L157:
/* Computing MAX */
			    r__1 = c25 / c27;
			    tmin = dmax(r__1,tmin);
L158:
			    if (c37 < 0.) {
				goto L159;
			    } else if (c37 == 0) {
				goto L160;
			    } else {
				goto L161;
			    }
L159:
/* Computing MIN */
			    r__1 = c36 / c37;
			    tmax = dmin(r__1,tmax);
			    goto L162;
L160:
			    if (c36 >= 0.) {
				goto L440;
			    } else {
				goto L162;
			    }
L161:
/* Computing MAX */
			    r__1 = c36 / c37;
			    tmin = dmax(r__1,tmin);
L162:
			    if (tmax - tmin < 1e-5f) {
				goto L440;
			    }
			    xpl = x4 + tmin * (x5 - x4);
			    ypl = y4 + tmin * (y5 - y4);
			    zpl = z4 + tmin * (z5 - z4);
			    xpr = x4 + tmax * (x5 - x4);
			    ypr = y4 + tmax * (y5 - y4);
			    zpr = z4 + tmax * (z5 - z4);

/*     the projections of line and plane intersect */
/*     now determine if plane covers line */

			    vx1t = x2[k] - x1[k];
			    vy1t = y2[k] - y1[k];
			    vz1t = z2[k] - z1[k];
			    vx2t = x3[k] - x1[k];
			    vy2t = y3[k] - y1[k];
			    vz2t = z3[k] - z1[k];
			    apl = vy1t * vz2t - vy2t * vz1t;
			    bpl = vx2t * vz1t - vx1t * vz2t;
			    cpl = vx1t * vy2t - vx2t * vy1t;
			    dpl = apl * x1[k] + bpl * y1[k] + cpl * z1[k];
			    vx3t = xpl - *xeye;
			    vy3t = ypl - *yeye;
			    vz3t = zpl - *zeye;
			    den = apl * vx3t + bpl * vy3t + cpl * vz3t;
			    til = 0.;
			    if (den == 0.) {
				goto L410;
			    }
			    til = (dpl - apl * *xeye - bpl * *yeye - cpl * *
				    zeye) / den;
L410:
			    vx3t = xpr - *xeye;
			    vy3t = ypr - *yeye;
			    vz3t = zpr - *zeye;
			    den = apl * vx3t + bpl * vy3t + cpl * vz3t;
			    tir = 0.;
			    if (den == 0.) {
				goto L412;
			    }
			    tir = (dpl - apl * *xeye - bpl * *yeye - cpl * *
				    zeye) / den;
L412:
			    if (til >= .99999 && tir >= .99999) {
				goto L440;
			    }
			    if (til < 1. && tir < 1.) {
				goto L164;
			    }
			    vx3t = xpr - xpl;
			    vy3t = ypr - ypl;
			    vz3t = zpr - zpl;
			    den = apl * vx3t + bpl * vy3t + cpl * vz3t;
			    tim = 0.;
			    if (den == 0.) {
				goto L414;
			    }
			    tim = (dpl - apl * xpl - bpl * ypl - cpl * zpl) / 
				    den;
L414:
			    thold = tmin + tim * (tmax - tmin);
			    if (til >= 1.) {
				goto L163;
			    }
			    tmax = thold;
			    goto L164;
L163:
			    tmin = thold;
L164:
			    ++nseg;
			    tl[nseg] = tmin;
			    tr[nseg] = tmax;
L440:
			    ;
			}
			maxs = max(maxs,nseg);
			if ((i__5 = nseg - 1) < 0) {
			    goto L171;
			} else if (i__5 == 0) {
			    goto L180;
			} else {
			    goto L172;
			}
L171:
			line_(&px4, &py4, &px5, &py5);
			goto L460;

/*     order the segments according to left end point tl(k) */

L172:
			i__5 = nseg;
			for (k = 2; k <= i__5; ++k) {
			    i__6 = nseg;
			    for (i__ = k; i__ <= i__6; ++i__) {
				if (tl[k - 1] <= tl[i__]) {
				    goto L173;
				}
				tlh = tl[k - 1];
				trh = tr[k - 1];
				tl[k - 1] = tl[i__];
				tr[k - 1] = tr[i__];
				tl[i__] = tlh;
				tr[i__] = trh;
L173:
				;
			    }
			}

/*     eliminate segment overlap */

			k1 = 1;
			k2 = 1;
L174:
			++k2;
			if (k2 > nseg) {
			    goto L176;
			}
			if (tr[k1] < tl[k2]) {
			    goto L175;
			}
/* Computing MAX */
			r__1 = tr[k1], r__2 = tr[k2];
			tr[k1] = dmax(r__1,r__2);
			goto L174;
L175:
			++k1;
			tl[k1] = tl[k2];
			tr[k1] = tr[k2];
			goto L174;
L176:
			nseg = k1;

/*     plot all segments of the line */

L180:
			i__6 = nseg;
			for (ks = 1; ks <= i__6; ++ks) {
			    kb = nseg - ks + 1;
			    tl[kb + 1] = tr[kb];
			    tr[kb] = tl[kb];
/* L181: */
			}
			tl[1] = 0.;
			tr[nseg + 1] = 1.;
			nsegp = nseg + 1;
			i__6 = nsegp;
			for (k = 1; k <= i__6; ++k) {
			    if ((r__1 = tr[k] - tl[k], dabs(r__1)) < 1e-6f) {
				goto L450;
			    }
			    xa = px4 + tl[k] * (px5 - px4);
			    ya = py4 + tl[k] * (py5 - py4);
			    xb = px4 + tr[k] * (px5 - px4);
			    yb = py4 + tr[k] * (py5 - py4);
			    line_(&xa, &ya, &xb, &yb);
L450:
			    ;
			}
L460:
			;
		    }
L470:
		    ;
		}
L480:
		;
	    }
/* L490: */
	}
L500:
	;
    }
/*      write(6,903) tl1,tl2 */
/* L903: */
/*      write(6,904) maxs */
/* L904: */
/*      write(6,250) (ncv(icv),icv=1,10) */
/* L250: */
    frame_();
    return 0;
} /* vsurf1_ */

/* Subroutine */ int prjct_(integer *init, doublereal *xeye, doublereal *yeye, doublereal *zeye,
	 doublereal *x, doublereal *y, doublereal *z__, doublereal *px, doublereal *py)
{
    /* System generated locals */
    doublereal r__1, r__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d1, d2, x1, y1, z1, cx1, cy1, cx2, cy2, cz2, cx3, cy3, cz3, 
	    rads1, rads2, ratio;


/*     subroutine prjct projects the point x,y,z onto a plane through */
/*     the origin that is perpendicular to a line between the origin */
/*     and the eye. the projection is along the line between the eye */
/*     and the point x,y,z. px and py are the coordinates of the */
/*     projection in the plane. */
/*     (version 2 , 12-10-82) */

    if (*init != 0) {
	goto L1;
    }
/* Computing 2nd power */
    r__1 = *xeye;
/* Computing 2nd power */
    r__2 = *yeye;
    rads1 = r__1 * r__1 + r__2 * r__2;
/* Computing 2nd power */
    r__1 = *zeye;
    rads2 = rads1 + r__1 * r__1;
    d1 = sqrt(rads1);
    d2 = sqrt(rads2);
    if (d1 != 0.) {
	goto L2;
    }
    cx1 = 1.;
    cy1 = 0.;
    cx2 = 0.;
    cy2 = 1.;
    cz2 = 0.;
    cx3 = 0.;
    cy3 = 0.;
    cz3 = 1.;
    return 0;
L2:
    cx1 = -(*yeye) / d1;
    cy1 = *xeye / d1;
    cx2 = -(*xeye) * *zeye / (d1 * d2);
    cy2 = -(*yeye) * *zeye / (d1 * d2);
    cz2 = d1 / d2;
    cx3 = *xeye / d2;
    cy3 = *yeye / d2;
    cz3 = *zeye / d2;
    return 0;
L1:
    x1 = cx1 * *x + cy1 * *y;
    y1 = cx2 * *x + cy2 * *y + cz2 * *z__;
    z1 = cx3 * *x + cy3 * *y + cz3 * *z__;
    ratio = d2 / (d2 - z1);
    *px = ratio * x1;
    *py = ratio * y1;
    return 0;
} /* prjct_ */

/* Subroutine */ int box_(integer *isd, integer *istart, integer *next, 
	integer *l, integer *list)
{
    /* System generated locals */
    integer istart_dim1, istart_offset, i__1, i__2;

    /* Local variables */
    static integer id, jd, idx;

    /* Parameter adjustments */
    istart_dim1 = *isd;
    istart_offset = 1 + istart_dim1;
    istart -= istart_offset;
    --next;
    --list;

    /* Function Body */
    i__1 = *isd;
    for (jd = 1; jd <= i__1; ++jd) {
	i__2 = *isd;
	for (id = 1; id <= i__2; ++id) {
	    idx = istart[id + jd * istart_dim1];
	    istart[id + jd * istart_dim1] = *l + 1;
	    if (idx == 0) {
		goto L10;
	    }
L20:
	    ++(*l);
	    list[*l] = idx;
	    if (next[idx] == 0) {
		goto L10;
	    }
	    idx = next[idx];
	    goto L20;
L10:
	    ;
	}
/* L30: */
    }
    return 0;
} /* box_ */

