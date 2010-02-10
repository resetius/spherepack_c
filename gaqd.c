/* gaqd.f -- translated by f2c (version 20061008).
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


/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
/*  .                                                             . */
/*  .                  copyright (c) 1998 by UCAR                 . */
/*  .                                                             . */
/*  .       University Corporation for Atmospheric Research       . */
/*  .                                                             . */
/*  .                      all rights reserved                    . */
/*  .                                                             . */
/*  .                                                             . */
/*  .                          SPHEREPACK                          . */
/*  .                                                             . */
/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/*                             August 2003 */


/*     This version of gaqd implements the method presented in: */
/*     P. N. swarztrauber, Computing the points and weights for */
/*     Gauss-Legendre quadrature, SIAM J. Sci. Comput., */
/*     24(2002) pp. 945-954. */

/*     It the version that is new to spherepack 3.1 */
/*     The w and lwork arrays are dummy and included only to */
/*     permit a simple pluggable exchange with the */
/*     old gaqd in spherepack 3.0. */



/* Subroutine */ int gaqd_(integer *nlat, doublereal *theta, doublereal *wts, 
	doublereal *w, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double acos(doublereal), sqrt(doublereal), atan(doublereal), cos(
	    doublereal), sin(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal x, pb, pi, cz;
    static integer it, ns2;
    static doublereal dpb;
    static integer idx;
    static doublereal eps;
    static integer nix;
    static doublereal sum, pis2;
    extern /* Subroutine */ int cpdp_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal dcor, cmax, sgnd;
    extern /* Subroutine */ int tpdp_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal zero;
    static integer nhalf, mnlat;
    static doublereal zhold;
    extern doublereal dzeps_(doublereal *);
    static doublereal zlast, zprev, dthalf, dtheta;


/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
/*  .                                                             . */
/*  .                  copyright (c) 2001 by ucar                 . */
/*  .                                                             . */
/*  .       university corporation for atmospheric research       . */
/*  .                                                             . */
/*  .                      all rights reserved                    . */
/*  .                                                             . */
/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/*                             February 2002 */

/*     gauss points and weights are computed using the fourier-newton */
/*     described in "on computing the points and weights for */
/*     gauss-legendre quadrature", paul n. swarztrauber, siam journal */
/*     on scientific computing that has been accepted for publication. */
/*     This routine is faster and more accurate than older program */
/*     with the same name. */

/*     subroutine gaqd computes the nlat gaussian colatitudes and weights */
/*     in double precision. the colatitudes are in radians and lie in the */
/*     in the interval (0,pi). */

/*     input parameters */

/*     nlat    the number of gaussian colatitudes in the interval (0,pi) */
/*             (between the two poles).  nlat must be greater than zero. */

/*     w       unused double precision variable that permits a simple */
/*             exchange with the old routine with the same name */
/*             in spherepack. */

/*     lwork   unused variable that permits a simple exchange with the */
/*             old routine with the same name in spherepack. */

/*     output parameters */

/*     theta   a double precision array with length nlat */
/*             containing the gaussian colatitudes in */
/*             increasing radians on the interval (0,pi). */

/*     wts     a double precision array with lenght nlat */
/*             containing the gaussian weights. */

/*     ierror = 0 no errors */
/*            = 1 if nlat.le.0 */

/*  ***************************************************************** */


/*     check work space length */

    /* Parameter adjustments */
    --wts;
    --theta;

    /* Function Body */
    *ierror = 1;
    if (*nlat <= 0) {
	return 0;
    }
    *ierror = 0;

/*     compute weights and points analytically when nlat=1,2 */

    if (*nlat == 1) {
	theta[1] = acos(0.);
	wts[1] = 2.;
	return 0;
    }
    if (*nlat == 2) {
	x = sqrt(.33333333333333331);
	theta[1] = acos(x);
	theta[2] = acos(-x);
	wts[1] = 1.;
	wts[2] = 1.;
	return 0;
    }
    eps = sqrt(dzeps_(&c_b2));
    eps *= sqrt(eps);
    pis2 = atan(1.) * 2.;
    pi = pis2 + pis2;
    mnlat = *nlat % 2;
    ns2 = *nlat / 2;
    nhalf = (*nlat + 1) / 2;
    idx = ns2 + 2;

    cpdp_(nlat, &cz, &theta[ns2 + 1], &wts[ns2 + 1]);

    dtheta = pis2 / nhalf;
    dthalf = dtheta / 2.;
    cmax = dtheta * .2;

/*     estimate first point next to theta = pi/2 */

    if (mnlat != 0) {
	zero = pis2 - dtheta;
	zprev = pis2;
	nix = nhalf - 1;
    } else {
	zero = pis2 - dthalf;
	nix = nhalf;
    }
L9:
    it = 0;
L10:
    ++it;
    zlast = zero;

/*     newton iterations */

    tpdp_(nlat, &zero, &cz, &theta[ns2 + 1], &wts[ns2 + 1], &pb, &dpb);
    dcor = pb / dpb;
    sgnd = 1.;
    if (dcor != 0.) {
	sgnd = dcor / abs(dcor);
    }
/* Computing MIN */
    d__1 = abs(dcor);
    dcor = sgnd * min(d__1,cmax);
    zero -= dcor;
    if ((d__1 = zero - zlast, abs(d__1)) > eps * abs(zero)) {
	goto L10;
    }
    theta[nix] = zero;
    zhold = zero;
/*      wts(nix) = (nlat+nlat+1)/(dpb*dpb) */

/*     yakimiw's formula permits using old pb and dpb */

/* Computing 2nd power */
    d__1 = dpb + pb * cos(zlast) / sin(zlast);
    wts[nix] = (*nlat + *nlat + 1) / (d__1 * d__1);
    --nix;
    if (nix == 0) {
	goto L30;
    }
    if (nix == nhalf - 1) {
	zero = zero * 3. - pi;
    }
    if (nix < nhalf - 1) {
	zero = zero + zero - zprev;
    }
    zprev = zhold;
    goto L9;

/*     extend points and weights via symmetries */

L30:
    if (mnlat != 0) {
	theta[nhalf] = pis2;
	tpdp_(nlat, &pis2, &cz, &theta[ns2 + 1], &wts[ns2 + 1], &pb, &dpb);
	wts[nhalf] = (*nlat + *nlat + 1) / (dpb * dpb);
    }
    i__1 = ns2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wts[*nlat - i__ + 1] = wts[i__];
	theta[*nlat - i__ + 1] = pi - theta[i__];
    }
    sum = 0.;
    i__1 = *nlat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum += wts[i__];
    }
    i__1 = *nlat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wts[i__] = wts[i__] * 2. / sum;
    }
    return 0;
} /* gaqd_ */

/* Subroutine */ int cpdp_(integer *n, doublereal *cz, doublereal *cp, 
	doublereal *dcp)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal t1, t2, t3, t4;
    static integer ncp;


/*     computes the fourier coefficients of the legendre */
/*     polynomial p_n^0 and its derivative. */
/*     n is the degree and n/2 or (n+1)/2 */
/*     coefficients are returned in cp depending on whether */
/*     n is even or odd. The same number of coefficients */
/*     are returned in dcp. For n even the constant */
/*     coefficient is returned in cz. */

    /* Parameter adjustments */
    --dcp;
    --cp;

    /* Function Body */
    ncp = (*n + 1) / 2;
    t1 = -1.;
    t2 = *n + 1.;
    t3 = 0.;
    t4 = *n + *n + 1.;
    if (*n % 2 == 0) {
	cp[ncp] = 1.;
	for (j = ncp; j >= 2; --j) {
	    t1 += 2.;
	    t2 += -1.;
	    t3 += 1.;
	    t4 += -2.;
	    cp[j - 1] = t1 * t2 / (t3 * t4) * cp[j];
	}
	t1 += 2.;
	t2 += -1.;
	t3 += 1.;
	t4 += -2.;
	*cz = t1 * t2 / (t3 * t4) * cp[1];
	i__1 = ncp;
	for (j = 1; j <= i__1; ++j) {
	    dcp[j] = (j + j) * cp[j];
	}
    } else {
	cp[ncp] = 1.;
	for (j = ncp - 1; j >= 1; --j) {
	    t1 += 2.;
	    t2 += -1.;
	    t3 += 1.;
	    t4 += -2.;
	    cp[j] = t1 * t2 / (t3 * t4) * cp[j + 1];
	}
	i__1 = ncp;
	for (j = 1; j <= i__1; ++j) {
	    dcp[j] = (j + j - 1) * cp[j];
	}
    }
    return 0;
} /* cpdp_ */

/* Subroutine */ int tpdp_(integer *n, doublereal *theta, doublereal *cz, 
	doublereal *cp, doublereal *dcp, doublereal *pb, doublereal *dpb)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal fn, chh, cdt;
    static integer kdo;
    static doublereal cth, sdt, sth;


/*     computes pn(theta) and its derivative dpb(theta) with */
/*     respect to theta */


    /* Parameter adjustments */
    --dcp;
    --cp;

    /* Function Body */
    fn = (doublereal) (*n);
    cdt = cos(*theta + *theta);
    sdt = sin(*theta + *theta);
    if (*n % 2 == 0) {

/*     n even */

	kdo = *n / 2;
	*pb = *cz * .5;
	*dpb = 0.;
	if (*n > 0) {
	    cth = cdt;
	    sth = sdt;
	    i__1 = kdo;
	    for (k = 1; k <= i__1; ++k) {
/*      pb = pb+cp(k)*cos(2*k*theta) */
		*pb += cp[k] * cth;
/*      dpb = dpb-(k+k)*cp(k)*sin(2*k*theta) */
		*dpb -= dcp[k] * sth;
		chh = cdt * cth - sdt * sth;
		sth = sdt * cth + cdt * sth;
		cth = chh;
/* L170: */
	    }
	}
    } else {

/*     n odd */

	kdo = (*n + 1) / 2;
	*pb = 0.;
	*dpb = 0.;
	cth = cos(*theta);
	sth = sin(*theta);
	i__1 = kdo;
	for (k = 1; k <= i__1; ++k) {
/*      pb = pb+cp(k)*cos((2*k-1)*theta) */
	    *pb += cp[k] * cth;
/*      dpb = dpb-(k+k-1)*cp(k)*sin((2*k-1)*theta) */
	    *dpb -= dcp[k] * sth;
	    chh = cdt * cth - sdt * sth;
	    sth = sdt * cth + cdt * sth;
	    cth = chh;
/* L190: */
	}
    }
    return 0;
} /* tpdp_ */

doublereal dzeps_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal a, b, c__, eps;


/*     estimate unit roundoff in quantities of size x. */


/*     this program should function properly on all systems */
/*     satisfying the following two assumptions, */
/*        1.  the base used in representing floating point */
/*            numbers is not a power of three. */
/*        2.  the quantity  a  in statement 10 is represented to */
/*            the accuracy used in floating point variables */
/*            that are stored in memory. */
/*     the statement number 10 and the go to 10 are intended to */
/*     force optimizing compilers to generate code satisfying */
/*     assumption 2. */
/*     under these assumptions, it should be true that, */
/*            a  is not exactly equal to four-thirds, */
/*            b  has a zero for its last bit or digit, */
/*            c  is not exactly equal to one, */
/*            eps  measures the separation of 1.0 from */
/*                 the next larger floating point number. */
/*     the developers of eispack would appreciate being informed */
/*     about any systems where these assumptions do not hold. */

/*     this version dated 4/6/83. */

    a = 1.3333333333333333;
L10:
    b = a - 1.;
    c__ = b + b + b;
    eps = (d__1 = c__ - 1., abs(d__1));
    if (eps == 0.) {
	goto L10;
    }
    ret_val = eps * abs(*x);
    return ret_val;
} /* dzeps_ */

