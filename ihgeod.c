/* ihgeod.f -- translated by f2c (version 20061008).
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

static doublereal c_b2 = 1.;
static doublereal c_b8 = 0.;


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


/* Subroutine */ int ihgeod_(integer *m, integer *idp, integer *jdp, 
	doublereal *x, doublereal *y, doublereal *z__)
{
    /* System generated locals */
    integer x_dim1, x_dim2, x_offset, y_dim1, y_dim2, y_offset, z_dim1, 
	    z_dim2, z_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal), acos(doublereal);

    /* Local variables */
    integer i__, j, k;
    doublereal x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6,
	     y6, z6, pi, xs, ys, zs, rad, phi, dxi, dyi, dzi, dxj, dyj, dzj, 
	    beta, dphi;
    extern /* Subroutine */ int stoc_(doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *), ctos_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    doublereal hdphi, theta, tdphi, theta1, theta2;


/*     m         is the number of points on the edge of a */
/*               single geodesic triangle */

/*     x,y,z     the coordinates of the geodesic points on */
/*               the sphere are x(i,j,k), y(i,j,k), z(i,j,k) */
/*               where i=1,...,m+m-1; j=1,...,m; and k=1,...,5. */
/*               the indices are defined on the unfolded */
/*               icosahedron as follows for the case m=3 */

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

/*                total number of points is 10*(m-1)**2+2 */
/*                total number of triangles is 20*(m-1)**2 */
/*                total number of edges is 30*(m-1)**2 */

    /* Parameter adjustments */
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

    /* Function Body */
    pi = atan(1.) * 4.;
    dphi = pi * .4;
    beta = cos(dphi);
    theta1 = acos(beta / (1. - beta));
    theta2 = pi - theta1;
    hdphi = dphi / 2.;
    tdphi = hdphi * 3.;
    for (k = 1; k <= 5; ++k) {
	phi = (k - 1) * dphi;
	stoc_(&c_b2, &theta2, &phi, &x1, &y1, &z1);
	d__1 = phi + hdphi;
	stoc_(&c_b2, &pi, &d__1, &x2, &y2, &z2);
	d__1 = phi + dphi;
	stoc_(&c_b2, &theta2, &d__1, &x3, &y3, &z3);
	dxi = (x2 - x1) / (*m - 1);
	dyi = (y2 - y1) / (*m - 1);
	dzi = (z2 - z1) / (*m - 1);
	dxj = (x3 - x2) / (*m - 1);
	dyj = (y3 - y2) / (*m - 1);
	dzj = (z3 - z2) / (*m - 1);
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xs = x1 + (i__ - 1) * dxi;
	    ys = y1 + (i__ - 1) * dyi;
	    zs = z1 + (i__ - 1) * dzi;
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		x[j + (i__ + k * x_dim2) * x_dim1] = xs + (j - 1) * dxj;
		y[j + (i__ + k * y_dim2) * y_dim1] = ys + (j - 1) * dyj;
		z__[j + (i__ + k * z_dim2) * z_dim1] = zs + (j - 1) * dzj;
	    }
	}
	d__1 = phi + hdphi;
	stoc_(&c_b2, &theta1, &d__1, &x4, &y4, &z4);
	dxi = (x3 - x4) / (*m - 1);
	dyi = (y3 - y4) / (*m - 1);
	dzi = (z3 - z4) / (*m - 1);
	dxj = (x4 - x1) / (*m - 1);
	dyj = (y4 - y1) / (*m - 1);
	dzj = (z4 - z1) / (*m - 1);
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    xs = x1 + (j - 1) * dxj;
	    ys = y1 + (j - 1) * dyj;
	    zs = z1 + (j - 1) * dzj;
	    i__2 = j;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[j + (i__ + k * x_dim2) * x_dim1] = xs + (i__ - 1) * dxi;
		y[j + (i__ + k * y_dim2) * y_dim1] = ys + (i__ - 1) * dyi;
		z__[j + (i__ + k * z_dim2) * z_dim1] = zs + (i__ - 1) * dzi;
	    }
	}
	d__1 = phi + tdphi;
	stoc_(&c_b2, &theta1, &d__1, &x5, &y5, &z5);
	dxj = (x5 - x3) / (*m - 1);
	dyj = (y5 - y3) / (*m - 1);
	dzj = (z5 - z3) / (*m - 1);
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    xs = x4 + (i__ - 1) * dxi;
	    ys = y4 + (i__ - 1) * dyi;
	    zs = z4 + (i__ - 1) * dzi;
	    i__2 = i__;
	    for (j = 1; j <= i__2; ++j) {
		x[j + *m - 1 + (i__ + k * x_dim2) * x_dim1] = xs + (j - 1) * 
			dxj;
		y[j + *m - 1 + (i__ + k * y_dim2) * y_dim1] = ys + (j - 1) * 
			dyj;
		z__[j + *m - 1 + (i__ + k * z_dim2) * z_dim1] = zs + (j - 1) *
			 dzj;
	    }
	}
	d__1 = phi + dphi;
	stoc_(&c_b2, &c_b8, &d__1, &x6, &y6, &z6);
	dxi = (x5 - x6) / (*m - 1);
	dyi = (y5 - y6) / (*m - 1);
	dzi = (z5 - z6) / (*m - 1);
	dxj = (x6 - x4) / (*m - 1);
	dyj = (y6 - y4) / (*m - 1);
	dzj = (z6 - z4) / (*m - 1);
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    xs = x4 + (j - 1) * dxj;
	    ys = y4 + (j - 1) * dyj;
	    zs = z4 + (j - 1) * dzj;
	    i__2 = j;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		x[j + *m - 1 + (i__ + k * x_dim2) * x_dim1] = xs + (i__ - 1) *
			 dxi;
		y[j + *m - 1 + (i__ + k * y_dim2) * y_dim1] = ys + (i__ - 1) *
			 dyi;
		z__[j + *m - 1 + (i__ + k * z_dim2) * z_dim1] = zs + (i__ - 1)
			 * dzi;
	    }
	}
    }
    for (k = 1; k <= 5; ++k) {
	i__1 = *m + *m - 1;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ctos_(&x[j + (i__ + k * x_dim2) * x_dim1], &y[j + (i__ + k * 
			y_dim2) * y_dim1], &z__[j + (i__ + k * z_dim2) * 
			z_dim1], &rad, &theta, &phi);
		stoc_(&c_b2, &theta, &phi, &x[j + (i__ + k * x_dim2) * x_dim1]
			, &y[j + (i__ + k * y_dim2) * y_dim1], &z__[j + (i__ 
			+ k * z_dim2) * z_dim1]);
	    }
	}
    }
    return 0;
} /* ihgeod_ */

/* Subroutine */ int ctos_(doublereal *x, doublereal *y, doublereal *z__, 
	doublereal *r__, doublereal *theta, doublereal *phi)
{
    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), atan2(doublereal, doublereal);

    /* Local variables */
    doublereal r1;

    r1 = *x * *x + *y * *y;
    if (r1 != 0.) {
	goto L10;
    }
    *phi = 0.;
    *theta = 0.;
    if (*z__ < 0.) {
	*theta = atan(1.) * 4.;
    }
    return 0;
L10:
    *r__ = sqrt(r1 + *z__ * *z__);
    r1 = sqrt(r1);
    *phi = atan2(*y, *x);
    *theta = atan2(r1, *z__);
    return 0;
} /* ctos_ */

/* Subroutine */ int stoc_(doublereal *r__, doublereal *theta, doublereal *
	phi, doublereal *x, doublereal *y, doublereal *z__)
{
    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    doublereal st;

    st = sin(*theta);
    *x = *r__ * st * cos(*phi);
    *y = *r__ * st * sin(*phi);
    *z__ = *r__ * cos(*theta);
    return 0;
} /* stoc_ */

