/* visequ.f -- translated by f2c (version 20061008).
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

static integer c__16 = 16;
static integer c__0 = 0;
static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b70 = 0.;
static doublereal c_b71 = 1.;


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

/* ... file visgau.f */

/*     contains documentation and code for subroutine visgau */

/* Subroutine */ int visequ_(integer *nlat, integer *nlon, doublereal *h__, 
	integer *len, doublereal *eyer, doublereal *eyelat, doublereal *
	eyelon, doublereal *wk, integer *lwk, integer *iwk, integer *liwk, 
	integer *ierror)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2;

    /* Builtin functions */
    double atan(doublereal), sin(doublereal), cos(doublereal);

    /* Local variables */
    static integer i__, j, m, n;
    static doublereal pi;
    static integer mn, nx, ny, nz, nw1, nw2, nx1, ny1, nz1, nx2, ny2, nz2, 
	    nx3, ny3, nz3;
    static doublereal dtr;
    static integer nxp, nyp;
    extern /* Subroutine */ int diag_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static integer itri;
    extern /* Subroutine */ int sptc_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *);
    static doublereal xeye, yeye;
    static integer ntri;
    static doublereal zeye;
    static integer nwrk, nmst, nmfac, nclat, nslat;
    extern /* Subroutine */ int vsurf_(doublereal *, doublereal *, doublereal 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    static integer niflag;
    extern /* Subroutine */ int triang_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *), stride_(integer *, integer *, integer *, integer *), 
	    interp_(doublereal *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, integer *), projct_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer nitype;


/*     subroutine visequ will display a function on the sphere */
/*     as a solid. ie. as a "lumpy" sphere. visequ calls subroutine */
/*     vsurf to produce the visible surface rendering. */

/*     requires routines visequ1 interp sptc diag stride triang vsurf */
/*                       vsurf1 prjct box icvmg projct */

/*     visgeo uses the ncar graphics package. */
/*     compile with: ncargf77 (all programs above) */

/*     execute with:  a.out */

/*     on screen display with:  ctrans -d x11 gmeta */

/*     print with:  ctrans -d ps.color gmeta > gmeta.ps */
/*                  lpr -p(your printer) gmeta.ps */

/*     input parameters */

/*     nlat   the number of colatitudes on the full sphere including the */
/*            poles. for example, nlat = 37 for a five degree grid. */
/*            nlat determines the grid increment in colatitude as */
/*            pi/(nlat-1).  if nlat is odd the equator is located at */
/*            grid point i=(nlat+1)/2. if nlat is even the equator is */
/*            located half way between points i=nlat/2 and i=nlat/2+1. */
/*            nlat must be at least 3. note: on the half sphere, the */
/*            number of grid points in the colatitudinal direction is */
/*            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd. */

/*     nlon   the number of distinct londitude points.  nlon determines */
/*            the grid increment in longitude as 2*pi/nlon. for example */
/*            nlon = 72 for a five degree grid. nlon must be greater */
/*            than or equal to 4. the efficiency of the computation is */
/*            improved when nlon is a product of small prime numbers. */


/*     h      a two dimensional array that contains the discrete */
/*            function to be displayed. h(i,j) is the distance from the */
/*            center of the sphere to the surface at the colatitude */
/*            point theta(i) = (i-1)*pi/(nlat-1) and longitude point */
/*            phi(j) = (j-1)*2*pi/nlon. */

/*     len    the first dimension of the array h as it appears in the */
/*            program that calls visequ. */

/*     eyer   the distance from the center of the sphere to the eye. */

/*     eyelat the colatitudinal coordinate of the eye (in degrees). */

/*     eyelon the longitudinal  coordinate of the eye (in degrees). */

/*     wk     a real work array */

/*     lwk    the dimension of the array wk as it appears in the */
/*            program that calls visequ. lwk must be at least */
/*                       46*nlat*(nlon+1). */

/*     iwk    an integer work array */

/*     liwk   the dimension of the array iwk as it appears in the */
/*            program that calls visequ. liwk must be at least */
/*                       14*nlat*(nlon+1). */

/*     ierror = 0    no error */
/*            = 1    the eye is positioned inside the sphere */
/*            = 2    lwk  is less than 46*nlat*(nlon+1) */
/*            = 3    liwk is less than 14*nlat*(nlon+1) */


    /* Parameter adjustments */
    h_dim1 = *len;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --wk;
    --iwk;

    /* Function Body */
    n = *nlat;
    m = *nlon + 1;
    mn = m * n;
    *ierror = 2;
    if (*lwk < mn * 46) {
	return 0;
    }
    *ierror = 3;
    if (*liwk < mn * 14) {
	return 0;
    }
    *ierror = 1;
    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (*eyer <= h__[i__ + j * h_dim1]) {
		return 0;
	    }
/* L10: */
	}
    }
    *ierror = 0;
    pi = atan(1.) * 4.;
    dtr = pi / 180.;

/*     ****     set up pointers to sub work arrays in wk and iwk */

    ntri = mn + mn;
    nw1 = 1;
    nw2 = nw1 + mn;
    nclat = 1;
    nslat = nclat + n;
    nxp = 1;
    nyp = nxp + mn;
    nx1 = 1;
    ny1 = nx1 + ntri;
    nz1 = ny1 + ntri;
    nx2 = nz1 + ntri;
    ny2 = nx2 + ntri;
    nz2 = ny2 + ntri;
    nx3 = nz2 + ntri;
    ny3 = nx3 + ntri;
    nz3 = ny3 + ntri;
    nx = nz3 + ntri;
    ny = nx + mn;
    nz = ny + mn;
    nwrk = nx;
    nitype = 1;
    niflag = ntri + 1;
    nmst = niflag + mn;
    nmfac = nmst + n;
/*     total iwk is 7*ntri */
/*     total wk is 58*nlat*(nlon+1) */
/*     ****     mid-cell interpolation, calculation of polar values */
    interp_(&h__[h_offset], len, &m, &n, &wk[nw1], &wk[nw2], &iwk[niflag]);
/*     ****     transform grid points to cartesian coordinates */
    sptc_(&h__[h_offset], len, &m, &n, &wk[nclat], &wk[nslat], &wk[nx], &wk[
	    ny], &wk[nz]);
/*     ****     transform eye position to cartesian coordinates */
    xeye = *eyer * sin(dtr * *eyelat);
    yeye = xeye * sin(dtr * *eyelon);
    xeye *= cos(dtr * *eyelon);
    zeye = *eyer * cos(dtr * *eyelat);
/*     ****     project grid points */
    projct_(&m, &n, &xeye, &yeye, &zeye, &wk[nx], &wk[ny], &wk[nz], &wk[nxp], 
	    &wk[nyp]);
/*     ****     check for visibility of cell boundaries */
    diag_(&m, &n, &wk[nxp], &wk[nyp], &iwk[niflag]);
/*     ****     compute longitude stride as a function of latitude */
    stride_(&m, &n, &iwk[nmst], &iwk[nmfac]);
/*     ****     perform triangulation */
    triang_(&m, &n, &wk[nx], &wk[ny], &wk[nz], &itri, &wk[nx1], &wk[ny1], &wk[
	    nz1], &wk[nx2], &wk[ny2], &wk[nz2], &wk[nx3], &wk[ny3], &wk[nz3], 
	    &iwk[nitype], &iwk[niflag], &iwk[nmst]);
/*     ****     call surface plotting routine */
    vsurf_(&xeye, &yeye, &zeye, &itri, &wk[nx1], &wk[ny1], &wk[nz1], &wk[nx2],
	     &wk[ny2], &wk[nz2], &wk[nx3], &wk[ny3], &wk[nz3], &iwk[nitype], &
	    wk[nwrk], &iwk[niflag]);
    return 0;
} /* visequ_ */

/* Subroutine */ int interp_(doublereal *h__, integer *len, integer *m, 
	integer *n, doublereal *w1, doublereal *w2, integer *iflag)
{
    /* Initialized data */

    static doublereal sten[16]	/* was [4][4] */ = { .015625,-.078125,
	    -.078125,.015625,-.078125,.390625,.390625,-.078125,-.078125,
	    .390625,.390625,-.078125,.015625,-.078125,-.078125,.015625 };

    /* System generated locals */
    integer h_dim1, h_offset, w1_dim1, w1_offset, w2_dim1, w2_offset, 
	    iflag_dim1, iflag_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k, l, n1, n2, mm1;
    extern integer icvmg_(integer *, integer *, doublereal *);

/*     ****     interpolates to mid points of grid cells using second */
/*     ****     order formula */
    /* Parameter adjustments */
    h_dim1 = *len;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    iflag_dim1 = *n;
    iflag_offset = 1 + iflag_dim1;
    iflag -= iflag_offset;
    w2_dim1 = *n + 2;
    w2_offset = 1 + w2_dim1;
    w2 -= w2_offset;
    w1_dim1 = *n;
    w1_offset = 1 + w1_dim1;
    w1 -= w1_offset;

    /* Function Body */
/*     ****     copy h to w2 */
    mm1 = *m - 1;
    i__1 = mm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    w2[j + (i__ + 1) * w2_dim1] = h__[j + i__ * h_dim1];
/* L1: */
	}
    }
/*     ****     add periodic points */
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
	w2[j + w2_dim1] = w2[j + *m * w2_dim1];
	w2[j + (*m + 1) * w2_dim1] = w2[j + (w2_dim1 << 1)];
	w2[j + (*m + 2) * w2_dim1] = w2[j + w2_dim1 * 3];
/* L2: */
    }
    n1 = 2;
    n2 = *n - 2;
/*     ****     perform interpolation */
/*     ****     set w1 to zero */
    i__2 = *m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    w1[j + i__ * w1_dim1] = 0.;
/* L7: */
	}
    }
/*     ****     interpolate */
    for (k = 1; k <= 4; ++k) {
	for (l = 1; l <= 4; ++l) {
	    i__1 = *m - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = n2;
		for (j = n1; j <= i__2; ++j) {
		    w1[j + i__ * w1_dim1] += w2[j + l - 2 + (i__ + k - 1) * 
			    w2_dim1] * sten[k + (l << 2) - 5];
/* L8: */
		}
	    }
	}
    }
/*     ****     set up iflag array */
/*     ****     iflag(j,i)=0  if diagonal is (j,i) to (j+1,i+1) */
/*     ****     iflag(j,i)=16 if diagonal is (j+1,i), (j,i+1) */
    i__2 = *m - 1;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = n2;
	for (j = n1; j <= i__1; ++j) {
	    d__3 = (d__1 = (w2[j + (i__ + 1) * w2_dim1] + w2[j + 1 + (i__ + 2)
		     * w2_dim1]) * .5 - w1[j + i__ * w1_dim1], abs(d__1)) - (
		    d__2 = (w2[j + (i__ + 2) * w2_dim1] + w2[j + 1 + (i__ + 1)
		     * w2_dim1]) * .5 - w1[j + i__ * w1_dim1], abs(d__2));
	    iflag[j + i__ * iflag_dim1] = icvmg_(&c__16, &c__0, &d__3);
/* L9: */
	}
    }
    return 0;
} /* interp_ */

/* Subroutine */ int sptc_(doublereal *r__, integer *len, integer *m, integer 
	*n, doublereal *clat, doublereal *slat, doublereal *x, doublereal *y, 
	doublereal *z__)
{
    /* System generated locals */
    integer r_dim1, r_offset, x_dim1, x_offset, y_dim1, y_offset, z_dim1, 
	    z_offset, i__1, i__2;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal dp, dt, pi, clon, slon;

/*     ****     transforms from spherical to cartesian coordinates */
    /* Parameter adjustments */
    r_dim1 = *len;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    z_dim1 = *n;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    y_dim1 = *n;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --slat;
    --clat;

    /* Function Body */
    pi = atan(1.) * 4.;
    dt = pi / (*n - 1);
    dp = (pi + pi) / (*m - 1);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	clat[j] = cos((j - 1) * dt);
	slat[j] = sin((j - 1) * dt);
/* L10: */
    }
    i__1 = *m - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	clon = cos((i__ - 1) * dp);
	slon = sin((i__ - 1) * dp);
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    x[j + i__ * x_dim1] = r__[j + i__ * r_dim1] * slat[j];
	    y[j + i__ * y_dim1] = x[j + i__ * x_dim1] * slon;
	    x[j + i__ * x_dim1] *= clon;
	    z__[j + i__ * z_dim1] = r__[j + i__ * r_dim1] * clat[j];
/* L20: */
	}
    }
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
	x[j + *m * x_dim1] = x[j + x_dim1];
	y[j + *m * y_dim1] = y[j + y_dim1];
	z__[j + *m * z_dim1] = z__[j + z_dim1];
/* L30: */
    }
    return 0;
} /* sptc_ */

/* Subroutine */ int diag_(integer *m, integer *n, doublereal *xp, doublereal 
	*yp, integer *iflag)
{
    /* System generated locals */
    integer xp_dim1, xp_offset, yp_dim1, yp_offset, iflag_dim1, iflag_offset, 
	    i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j;


/*     ****     label visibility of cell sides */

/*     north side corresponds to j */
/*     south side corresponds to j+1 */
/*     west  side corresponds to i */
/*     east  side corresponds to i+1 */

/*     let iflag = b4 b3 b2 b1 b0 (in binary) then b0 through b3 are */
/*     either o or 1 depending on whether the east, south, north */
/*     or west side is either invisible or visible, respectively. */

/*     b4 is o if the diagonal is from (i,j) to (i+1,j+1) and 1 if */
/*     the diagonal is from (i,j+1) to (i+1,j). */

/*     ****     arithmetic statement function */
    /* Parameter adjustments */
    iflag_dim1 = *n;
    iflag_offset = 1 + iflag_dim1;
    iflag -= iflag_offset;
    yp_dim1 = *n;
    yp_offset = 1 + yp_dim1;
    yp -= yp_offset;
    xp_dim1 = *n;
    xp_offset = 1 + xp_dim1;
    xp -= xp_offset;

    /* Function Body */
    i__1 = *n - 2;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *m - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (iflag[j + i__ * iflag_dim1] >= 16) {
		goto L20;
	    }
	    i__3 = j + 1;
	    i__4 = i__ + 1;
	    i__5 = j + 1;
	    if ((xp[i__3 + i__4 * xp_dim1] - xp[i__5 + i__ * xp_dim1]) * (yp[
		    j + i__ * yp_dim1] - yp[i__5 + i__ * yp_dim1]) - (xp[j + 
		    i__ * xp_dim1] - xp[i__5 + i__ * xp_dim1]) * (yp[i__3 + 
		    i__4 * yp_dim1] - yp[i__5 + i__ * yp_dim1]) <= 0.) {
		goto L10;
	    }
/*     west and south are visible */
	    iflag[j + i__ * iflag_dim1] += 10;
L10:
	    i__3 = i__ + 1;
	    i__4 = j + 1;
	    i__5 = i__ + 1;
	    if ((xp[j + i__ * xp_dim1] - xp[j + i__3 * xp_dim1]) * (yp[i__4 + 
		    i__5 * yp_dim1] - yp[j + i__3 * yp_dim1]) - (xp[i__4 + 
		    i__5 * xp_dim1] - xp[j + i__3 * xp_dim1]) * (yp[j + i__ * 
		    yp_dim1] - yp[j + i__3 * yp_dim1]) <= 0.) {
		goto L100;
	    }
/*     east and north are visible */
	    iflag[j + i__ * iflag_dim1] += 5;
	    goto L100;
L20:
	    i__3 = j + 1;
	    i__4 = i__ + 1;
	    if ((xp[i__3 + i__ * xp_dim1] - xp[j + i__ * xp_dim1]) * (yp[j + 
		    i__4 * yp_dim1] - yp[j + i__ * yp_dim1]) - (xp[j + i__4 * 
		    xp_dim1] - xp[j + i__ * xp_dim1]) * (yp[i__3 + i__ * 
		    yp_dim1] - yp[j + i__ * yp_dim1]) <= 0.) {
		goto L30;
	    }
/*     west and north are visible */
	    iflag[j + i__ * iflag_dim1] += 12;
L30:
	    i__3 = i__ + 1;
	    i__4 = j + 1;
	    i__5 = i__ + 1;
	    i__6 = j + 1;
	    if ((xp[j + i__3 * xp_dim1] - xp[i__4 + i__5 * xp_dim1]) * (yp[
		    i__6 + i__ * yp_dim1] - yp[i__4 + i__5 * yp_dim1]) - (xp[
		    i__6 + i__ * xp_dim1] - xp[i__4 + i__5 * xp_dim1]) * (yp[
		    j + i__3 * yp_dim1] - yp[i__4 + i__5 * yp_dim1]) <= 0.) {
		goto L100;
	    }
/*     east and south are visible */
	    iflag[j + i__ * iflag_dim1] += 3;
L100:
	    ;
	}
    }

/*     classify the poles */

    i__2 = *m - 1;
    for (i__ = 1; i__ <= i__2; ++i__) {
	iflag[i__ * iflag_dim1 + 1] = 0;
	i__1 = i__ + 1;
	if ((xp[c__2 + i__1 * xp_dim1] - xp[c__2 + i__ * xp_dim1]) * (yp[c__1 
		+ i__ * yp_dim1] - yp[c__2 + i__ * yp_dim1]) - (xp[c__1 + i__ 
		* xp_dim1] - xp[c__2 + i__ * xp_dim1]) * (yp[c__2 + i__1 * 
		yp_dim1] - yp[c__2 + i__ * yp_dim1]) > 0.) {
	    iflag[i__ * iflag_dim1 + 1] = 15;
	}
	iflag[*n - 1 + i__ * iflag_dim1] = 0;
	i__1 = *n - 1;
	i__3 = *n - 1;
	i__4 = i__ + 1;
	if ((xp[*n + i__ * xp_dim1] - xp[i__1 + i__ * xp_dim1]) * (yp[i__3 + 
		i__4 * yp_dim1] - yp[i__1 + i__ * yp_dim1]) - (xp[i__3 + i__4 
		* xp_dim1] - xp[i__1 + i__ * xp_dim1]) * (yp[*n + i__ * 
		yp_dim1] - yp[i__1 + i__ * yp_dim1]) > 0.) {
	    iflag[*n - 1 + i__ * iflag_dim1] = 31;
	}
/* L200: */
    }
    i__2 = *n - 1;
    for (j = 1; j <= i__2; ++j) {
	iflag[j + *m * iflag_dim1] = iflag[j + iflag_dim1];
/* L250: */
    }
    return 0;
} /* diag_ */

/* Subroutine */ int stride_(integer *m, integer *n, integer *mst, integer *
	mfac)
{
    /* Initialized data */

    static integer mtryh[3] = { 2,3,5 };

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double atan(doublereal), sin(doublereal);

    /* Local variables */
    static integer j, jf, nf;
    static doublereal dt;
    static integer ml;
    static doublereal pi;
    static integer mq, mr;
    static doublereal st;
    static integer mf1, mf2, ns2, jdo;
    static doublereal tphi;
    static integer mtry;
    static doublereal theta;

    /* Parameter adjustments */
    --mst;
    --mfac;

    /* Function Body */

/*     find prime factors of m-1 */

    ml = *m - 1;
    nf = 0;
    j = 0;
L101:
    ++j;
    if (j - 3 <= 0) {
	goto L102;
    } else {
	goto L103;
    }
L102:
    mtry = mtryh[j - 1];
    goto L104;
L103:
    mtry += 2;
L104:
    mq = ml / mtry;
    mr = ml - mtry * mq;
    if (mr != 0) {
	goto L101;
    } else {
	goto L105;
    }
L105:
    ++nf;
    mfac[nf] = mtry;
    ml = mq;
    if (ml != 1) {
	goto L104;
    }
    if (mfac[nf] > 2) {
	goto L106;
    }
    --nf;
    mfac[nf] = 4;
L106:
    tphi = .707 / (doublereal) (*m - 1);
    ns2 = *n / 2;
    mf1 = mfac[nf];
    mst[1] = (*m - 1) / mf1;
    pi = atan(1.) * 4.;
    dt = pi / (doublereal) (*n - 1);
    jf = nf - 1;
    i__1 = ns2;
    for (jdo = 2; jdo <= i__1; ++jdo) {
	j = jdo;
	theta = (j - 1) * dt;
	st = sin(theta);
	mf2 = mf1 * mfac[jf];
	if ((d__1 = st / mf1 - tphi, abs(d__1)) > (d__2 = st / mf2 - tphi, 
		abs(d__2))) {
	    goto L115;
	}
	mst[j] = mst[j - 1];
	goto L110;
L115:
	mst[j] = (*m - 1) / mf2;
	mf1 = mf2;
	--jf;
	if (jf == 0) {
	    goto L120;
	}
L110:
	;
    }
L120:
    i__1 = ns2;
    for (jdo = j; jdo <= i__1; ++jdo) {
	mst[jdo] = 1;
/* L125: */
    }
    i__1 = ns2;
    for (jdo = 1; jdo <= i__1; ++jdo) {
	mst[*n - jdo] = mst[jdo];
/* L130: */
    }
/*      write (6,135) (mst(j),j=1,n) */
/* L135: */

    return 0;
} /* stride_ */

/* Subroutine */ int triang_(integer *m, integer *n, doublereal *x, 
	doublereal *y, doublereal *z__, integer *itri, doublereal *x1, 
	doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, 
	doublereal *z2, doublereal *x3, doublereal *y3, doublereal *z3, 
	integer *ityp, integer *iflag, integer *mst)
{
    /* Initialized data */

    static integer icl[8] = { 0,1,2,12,3,13,23,123 };

    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, z_dim1, z_offset, iflag_dim1, 
	    iflag_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, n1, n2, ityph;

/*     ****     performs triangulation */
    /* Parameter adjustments */
    --mst;
    iflag_dim1 = *n;
    iflag_offset = 1 + iflag_dim1;
    iflag -= iflag_offset;
    z_dim1 = *n;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    y_dim1 = *n;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --x1;
    --y1;
    --z1;
    --x2;
    --y2;
    --z2;
    --x3;
    --y3;
    --z3;
    --ityp;

    /* Function Body */
    *itri = 0;
    n1 = 2;
    n2 = *n - 2;
    i__1 = n2;
    for (j = n1; j <= i__1; ++j) {
	i__2 = *m - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (iflag[j + i__ * iflag_dim1] >= 16) {
		goto L50;
	    }
	    if (iflag[j + i__ * iflag_dim1] % 16 < 8) {
		goto L70;
	    }
	    ++(*itri);
	    x1[*itri] = x[j + i__ * x_dim1];
	    y1[*itri] = y[j + i__ * y_dim1];
	    z1[*itri] = z__[j + i__ * z_dim1];
	    x2[*itri] = x[j + 1 + i__ * x_dim1];
	    y2[*itri] = y[j + 1 + i__ * y_dim1];
	    z2[*itri] = z__[j + 1 + i__ * z_dim1];
	    x3[*itri] = x[j + 1 + (i__ + 1) * x_dim1];
	    y3[*itri] = y[j + 1 + (i__ + 1) * y_dim1];
	    z3[*itri] = z__[j + 1 + (i__ + 1) * z_dim1];
	    ityph = 3;
	    if ((i__ - 1) % mst[j] == 0) {
		goto L60;
	    }
	    if (iflag[j + (i__ - 1) * iflag_dim1] % 2 == 0) {
		goto L60;
	    }
	    --ityph;
L60:
	    if (iflag[j + i__ * iflag_dim1] % 2 == 0) {
		ityph += 4;
	    }
	    ityp[*itri] = icl[ityph];
L70:
	    if (iflag[j + i__ * iflag_dim1] % 2 == 0) {
		goto L100;
	    }
	    ++(*itri);
	    x1[*itri] = x[j + i__ * x_dim1];
	    y1[*itri] = y[j + i__ * y_dim1];
	    z1[*itri] = z__[j + i__ * z_dim1];
	    x2[*itri] = x[j + 1 + (i__ + 1) * x_dim1];
	    y2[*itri] = y[j + 1 + (i__ + 1) * y_dim1];
	    z2[*itri] = z__[j + 1 + (i__ + 1) * z_dim1];
	    x3[*itri] = x[j + (i__ + 1) * x_dim1];
	    y3[*itri] = y[j + (i__ + 1) * y_dim1];
	    z3[*itri] = z__[j + (i__ + 1) * z_dim1];
	    ityph = 0;
	    if (iflag[j + i__ * iflag_dim1] % 16 < 8) {
		++ityph;
	    }
	    if (iflag[j + (i__ + 1) * iflag_dim1] % 16 < 8) {
		ityph += 2;
	    }
	    if (iflag[j - 1 + i__ * iflag_dim1] % 4 < 2) {
		ityph += 4;
	    }
	    ityp[*itri] = icl[ityph];
	    goto L100;
L50:
	    if (iflag[j + i__ * iflag_dim1] % 16 < 8) {
		goto L20;
	    }
	    ++(*itri);
	    x1[*itri] = x[j + i__ * x_dim1];
	    y1[*itri] = y[j + i__ * y_dim1];
	    z1[*itri] = z__[j + i__ * z_dim1];
	    x2[*itri] = x[j + 1 + i__ * x_dim1];
	    y2[*itri] = y[j + 1 + i__ * y_dim1];
	    z2[*itri] = z__[j + 1 + i__ * z_dim1];
	    x3[*itri] = x[j + (i__ + 1) * x_dim1];
	    y3[*itri] = y[j + (i__ + 1) * y_dim1];
	    z3[*itri] = z__[j + (i__ + 1) * z_dim1];
	    ityph = 1;
	    if ((i__ - 1) % mst[j] == 0) {
		goto L10;
	    }
	    if (iflag[j + (i__ - 1) * iflag_dim1] % 2 == 0) {
		goto L10;
	    }
	    ityph = 0;
L10:
	    if (iflag[j + i__ * iflag_dim1] % 2 == 0) {
		ityph += 2;
	    }
	    if (iflag[j - 1 + i__ * iflag_dim1] % 4 < 2) {
		ityph += 4;
	    }
	    ityp[*itri] = icl[ityph];
L20:
	    if (iflag[j + i__ * iflag_dim1] % 2 == 0) {
		goto L100;
	    }
	    ++(*itri);
	    x1[*itri] = x[j + 1 + i__ * x_dim1];
	    y1[*itri] = y[j + 1 + i__ * y_dim1];
	    z1[*itri] = z__[j + 1 + i__ * z_dim1];
	    x2[*itri] = x[j + 1 + (i__ + 1) * x_dim1];
	    y2[*itri] = y[j + 1 + (i__ + 1) * y_dim1];
	    z2[*itri] = z__[j + 1 + (i__ + 1) * z_dim1];
	    x3[*itri] = x[j + (i__ + 1) * x_dim1];
	    y3[*itri] = y[j + (i__ + 1) * y_dim1];
	    z3[*itri] = z__[j + (i__ + 1) * z_dim1];
	    ityph = 1;
	    if (iflag[j + (i__ + 1) * iflag_dim1] % 16 < 8) {
		ityph += 2;
	    }
	    if (iflag[j + i__ * iflag_dim1] % 16 < 8) {
		ityph += 4;
	    }
	    ityp[*itri] = icl[ityph];
L100:
	    ;
	}
    }

/*     ****     triangles around north and south poles */

    i__2 = *m - 1;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (iflag[i__ * iflag_dim1 + 1] % 16 < 8) {
	    goto L250;
	}
	++(*itri);
	x1[*itri] = x[i__ * x_dim1 + 1];
	y1[*itri] = y[i__ * y_dim1 + 1];
	z1[*itri] = z__[i__ * z_dim1 + 1];
	x2[*itri] = x[i__ * x_dim1 + 2];
	y2[*itri] = y[i__ * y_dim1 + 2];
	z2[*itri] = z__[i__ * z_dim1 + 2];
	x3[*itri] = x[(i__ + 1) * x_dim1 + 2];
	y3[*itri] = y[(i__ + 1) * y_dim1 + 2];
	z3[*itri] = z__[(i__ + 1) * z_dim1 + 2];
	ityph = 3;
	if ((i__ - 1) % mst[1] == 0) {
	    goto L260;
	}
	if (iflag[(i__ - 1) * iflag_dim1 + 1] % 2 == 0) {
	    goto L260;
	}
	--ityph;
L260:
	if (iflag[(i__ + 1) * iflag_dim1 + 1] % 16 < 8) {
	    ityph += 4;
	}
	ityp[*itri] = icl[ityph];
L250:
	if (iflag[*n - 1 + i__ * iflag_dim1] % 16 < 8) {
	    goto L200;
	}
	++(*itri);
	x1[*itri] = x[*n - 1 + i__ * x_dim1];
	y1[*itri] = y[*n - 1 + i__ * y_dim1];
	z1[*itri] = z__[*n - 1 + i__ * z_dim1];
	x2[*itri] = x[*n + i__ * x_dim1];
	y2[*itri] = y[*n + i__ * y_dim1];
	z2[*itri] = z__[*n + i__ * z_dim1];
	x3[*itri] = x[*n - 1 + (i__ + 1) * x_dim1];
	y3[*itri] = y[*n - 1 + (i__ + 1) * y_dim1];
	z3[*itri] = z__[*n - 1 + (i__ + 1) * z_dim1];
	ityph = 1;
	if ((i__ - 1) % mst[*n - 1] == 0) {
	    goto L210;
	}
	if (iflag[*n - 1 + (i__ - 1) * iflag_dim1] % 2 == 0) {
	    goto L210;
	}
	ityph = 0;
L210:
	if (iflag[*n - 1 + (i__ + 1) * iflag_dim1] % 16 < 8) {
	    ityph += 2;
	}
	if (iflag[*n - 2 + i__ * iflag_dim1] % 4 < 2) {
	    ityph += 4;
	}
	ityp[*itri] = icl[ityph];
L200:
	;
    }
    return 0;
} /* triang_ */

/* Subroutine */ int vsurf_(doublereal *xeye, doublereal *yeye, doublereal *
	zeye, integer *ntri, doublereal *x1, doublereal *y1, doublereal *z1, 
	doublereal *x2, doublereal *y2, doublereal *z2, doublereal *x3, 
	doublereal *y3, doublereal *z3, integer *itype, doublereal *work, 
	integer *iwork)
{
    extern /* Subroutine */ int vsurf1_(doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *, integer *);


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

/*    the length of real    array  work must be at least 14*ntri */

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

/* Subroutine */ int vsurf1_(doublereal *xeye, doublereal *yeye, doublereal *
	zeye, integer *ntri, doublereal *x1, doublereal *y1, doublereal *z1, 
	doublereal *x2, doublereal *y2, doublereal *z2, doublereal *x3, 
	doublereal *y3, doublereal *z3, integer *itype, doublereal *px1, 
	doublereal *py1, doublereal *px2, doublereal *py2, doublereal *px3, 
	doublereal *py3, doublereal *vx1, doublereal *vy1, doublereal *vx2, 
	doublereal *vy2, doublereal *vx3, doublereal *vy3, doublereal *tl, 
	doublereal *tr, integer *kh, integer *next, integer *istart, integer *
	ifinal)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;

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
    static doublereal l2e, le2;
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
    extern /* Subroutine */ int set_(doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), box_(integer *, integer *, integer *, 
	    integer *, integer *);
    static integer kcv, ltp[3];
    static doublereal dxt;
    static integer ist, ixh, isx, ixs, ixf, kds, ldo, ith;
    static doublereal xpl, ypl, zpl, xpr, ypr, zpr, cpl, dpl, til, tir, tim, 
	    tlh, trh, dum1, dum2;
    static integer irp1;
    static doublereal px1h, py1h, vx1t, vy1t, vz1t, vx2t, vy2t, vz2t, vx3t, 
	    vy3t, vz3t;
    extern /* Subroutine */ int line_(doublereal *, doublereal *, doublereal *
	    , doublereal *);
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
    extern /* Subroutine */ int prjct_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
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
    irmax = (integer) (l2e * .5 * log(fntri));
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
	d__1 = pmin, d__2 = px1[k], d__1 = min(d__1,d__2), d__2 = py1[k], 
		d__1 = min(d__1,d__2), d__2 = px2[k], d__1 = min(d__1,d__2), 
		d__2 = py2[k], d__1 = min(d__1,d__2), d__2 = px3[k], d__1 = 
		min(d__1,d__2), d__2 = py3[k];
	pmin = min(d__1,d__2);
/* Computing MAX */
	d__1 = pmax, d__2 = px1[k], d__1 = max(d__1,d__2), d__2 = py1[k], 
		d__1 = max(d__1,d__2), d__2 = px2[k], d__1 = max(d__1,d__2), 
		d__2 = py2[k], d__1 = max(d__1,d__2), d__2 = px3[k], d__1 = 
		max(d__1,d__2), d__2 = py3[k];
	pmax = max(d__1,d__2);
/* L87: */
    }
    pmin *= 1.1;
    pmax *= 1.1;
    set_(&c_b70, &c_b71, &c_b70, &c_b71, &pmin, &pmax, &pmin, &pmax, &c__1);
/* Computing MIN */
    d__1 = min(px1[1],px2[1]);
    xmin = min(d__1,px3[1]);
/* Computing MAX */
    d__1 = max(px1[1],px2[1]);
    xmax = max(d__1,px3[1]);
/* Computing MIN */
    d__1 = min(py1[1],py2[1]);
    ymin = min(d__1,py3[1]);
/* Computing MAX */
    d__1 = max(py1[1],py2[1]);
    ymax = max(d__1,py3[1]);
    i__1 = *ntri;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MIN */
	d__1 = xmin, d__2 = px1[i__], d__1 = min(d__1,d__2), d__2 = px2[i__], 
		d__1 = min(d__1,d__2), d__2 = px3[i__];
	xmin = min(d__1,d__2);
/* Computing MAX */
	d__1 = xmax, d__2 = px1[i__], d__1 = max(d__1,d__2), d__2 = px2[i__], 
		d__1 = max(d__1,d__2), d__2 = px3[i__];
	xmax = max(d__1,d__2);
/* Computing MIN */
	d__1 = ymin, d__2 = py1[i__], d__1 = min(d__1,d__2), d__2 = py2[i__], 
		d__1 = min(d__1,d__2), d__2 = py3[i__];
	ymin = min(d__1,d__2);
/* Computing MAX */
	d__1 = ymax, d__2 = py1[i__], d__1 = max(d__1,d__2), d__2 = py2[i__], 
		d__1 = max(d__1,d__2), d__2 = py3[i__];
	ymax = max(d__1,d__2);
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
	d__1 = px1[i__], d__2 = px2[i__], d__1 = min(d__1,d__2), d__2 = px3[
		i__];
	xmin = min(d__1,d__2);
/* Computing MAX */
	d__1 = px1[i__], d__2 = px2[i__], d__1 = max(d__1,d__2), d__2 = px3[
		i__];
	xmax = max(d__1,d__2);
/* Computing MIN */
	d__1 = py1[i__], d__2 = py2[i__], d__1 = min(d__1,d__2), d__2 = py3[
		i__];
	ymin = min(d__1,d__2);
/* Computing MAX */
	d__1 = py1[i__], d__2 = py2[i__], d__1 = max(d__1,d__2), d__2 = py3[
		i__];
	ymax = max(d__1,d__2);
/* Computing MAX */
	d__1 = xmax - xmin, d__2 = ymax - ymin;
	dxt = max(d__1,d__2);
	if (dxt > 0.) {
	    goto L10;
	}
	ir = irmax;
	goto L20;
L10:
	ir = (integer) (l2e * log(hgr / dxt));
	ir = min(ir,irmax);
L20:
	irp1 = ir + 1;
	++nct[irp1 - 1];
	hr = hgr / ip2[irp1 - 1];
	xmid = (xmin + xmax) * .5;
	id = (integer) ((xmid - a) / hr + 1.5);
	ymid = (ymin + ymax) * .5;
	jd = (integer) ((ymid - c__) / hr + 1.5);
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
			    d__1 = c14 / c17;
			    tmax = min(d__1,tmax);
			    goto L154;
L152:
			    if (c14 >= 0.) {
				goto L440;
			    } else {
				goto L154;
			    }
L153:
/* Computing MAX */
			    d__1 = c14 / c17;
			    tmin = max(d__1,tmin);
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
			    d__1 = c25 / c27;
			    tmax = min(d__1,tmax);
			    goto L158;
L156:
			    if (c25 >= 0.) {
				goto L440;
			    } else {
				goto L158;
			    }
L157:
/* Computing MAX */
			    d__1 = c25 / c27;
			    tmin = max(d__1,tmin);
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
			    d__1 = c36 / c37;
			    tmax = min(d__1,tmax);
			    goto L162;
L160:
			    if (c36 >= 0.) {
				goto L440;
			    } else {
				goto L162;
			    }
L161:
/* Computing MAX */
			    d__1 = c36 / c37;
			    tmin = max(d__1,tmin);
L162:
			    if (tmax - tmin < 1e-5) {
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
			d__1 = tr[k1], d__2 = tr[k2];
			tr[k1] = max(d__1,d__2);
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
			    if ((d__1 = tr[k] - tl[k], abs(d__1)) < 1e-6) {
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

/* Subroutine */ int prjct_(integer *init, doublereal *xeye, doublereal *yeye,
	 doublereal *zeye, doublereal *x, doublereal *y, doublereal *z__, 
	doublereal *px, doublereal *py)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d1, d2, x1, y1, z1, cx1, cy1, cx2, cy2, cz2, cx3, cy3, 
	    cz3, rads1, rads2, ratio;


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
    d__1 = *xeye;
/* Computing 2nd power */
    d__2 = *yeye;
    rads1 = d__1 * d__1 + d__2 * d__2;
/* Computing 2nd power */
    d__1 = *zeye;
    rads2 = rads1 + d__1 * d__1;
    d1 = sqrt(rads1);
    d2 = sqrt(rads2);
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

integer icvmg_(integer *i1, integer *i2, doublereal *r__)
{
    /* System generated locals */
    integer ret_val;


/*     returns i1 if i3.ge.0 and returns i2 if i3.lt.0 . */

    ret_val = *i1;
    if (*r__ < 0.) {
	ret_val = *i2;
    }
    return ret_val;
} /* icvmg_ */

/* Subroutine */ int projct_(integer *m, integer *n, doublereal *xeye, 
	doublereal *yeye, doublereal *zeye, doublereal *x, doublereal *y, 
	doublereal *z__, doublereal *px, doublereal *py)
{
    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, z_dim1, z_offset, px_dim1, 
	    px_offset, py_dim1, py_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal rdum1, rdum2, rdum3, rdum4, rdum5;
    extern /* Subroutine */ int prjct_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*     ****     projects point (x,y,z) onto plane thru origin and perp */
/*     ****     to line joining origin and eye */
    /* Parameter adjustments */
    py_dim1 = *n;
    py_offset = 1 + py_dim1;
    py -= py_offset;
    px_dim1 = *n;
    px_offset = 1 + px_dim1;
    px -= px_offset;
    z_dim1 = *n;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    y_dim1 = *n;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    x_dim1 = *n;
    x_offset = 1 + x_dim1;
    x -= x_offset;

    /* Function Body */
    prjct_(&c__0, xeye, yeye, zeye, &rdum1, &rdum2, &rdum3, &rdum4, &rdum5);
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    prjct_(&c__1, xeye, yeye, zeye, &x[j + i__ * x_dim1], &y[j + i__ *
		     y_dim1], &z__[j + i__ * z_dim1], &px[j + i__ * px_dim1], 
		    &py[j + i__ * py_dim1]);
/* L100: */
	}
    }
    return 0;
} /* projct_ */

