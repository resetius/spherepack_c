/* lfpt.f -- translated by f2c (version 20061008).
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

/* subroutine lfpt (n,m,theta,cp,pb) */

/* dimension of */
/* arguments */
/*                        cp((n/2)+1) */

/* purpose                routine lfpt uses coefficients computed by */
/*                        routine alfk to compute the single precision */
/*                        normalized associated legendre function pbar(n, */
/*                        m,theta) at colatitude theta. */

/* usage                  call lfpt(n,m,theta,cp,pb) */

/* arguments */

/* on input               n */
/*                          nonnegative integer specifying the degree of */
/*                          pbar(n,m,theta) */
/*                        m */
/*                          is the order of pbar(n,m,theta). m can be */
/*                          any integer however pbar(n,m,theta) = 0 */
/*                          if abs(m) is greater than n and */
/*                          pbar(n,m,theta) = (-1)**m*pbar(n,-m,theta) */
/*                          for negative m. */

/*                        theta */
/*                          single precision colatitude in radians */

/*                        cp */
/*                          single precision array of length (n/2)+1 */
/*                          containing coefficients computed by routine */
/*                          alfk */

/* on output              pb */
/*                          single precision variable containing */
/*                          pbar(n,m,theta) */

/* special conditions     calls to routine lfpt must be preceded by an */
/*                        appropriate call to routine alfk. */

/* precision              single */

/* algorithm              the trigonometric series formula used by */
/*                        routine lfpt to calculate pbar(n,m,th) at */
/*                        colatitude th depends on m and n as follows: */

/*                           1) for n even and m even, the formula is */
/*                              .5*cp(1) plus the sum from k=1 to k=n/2 */
/*                              of cp(k)*cos(2*k*th) */
/*                           2) for n even and m odd. the formula is */
/*                              the sum from k=1 to k=n/2 of */
/*                              cp(k)*sin(2*k*th) */
/*                           3) for n odd and m even, the formula is */
/*                              the sum from k=1 to k=(n+1)/2 of */
/*                              cp(k)*cos((2*k-1)*th) */
/*                           4) for n odd and m odd, the formula is */
/*                              the sum from k=1 to k=(n+1)/2 of */
/*                              cp(k)*sin((2*k-1)*th) */

/* accuracy               comparison between routines lfpt and double */
/*                        precision dlfpt on the cray1 indicates greater */
/*                        accuracy for greater values on input parameter */
/*                        n.  agreement to 13 places was obtained for */
/*                        n=10 and to 12 places for n=100. */

/* timing                 time per call to routine lfpt is dependent on */
/*                        the input parameter n. */

/* Subroutine */ int lfpt_(integer *n, integer *m, real *theta, real *cp, 
	real *pb)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k, ma;
    static real ct, st;
    static integer kp1, np1;
    static real cdt;
    static integer kdo;
    static real cth, sdt, sum;
    static integer mmod, nmod;


    /* Parameter adjustments */
    --cp;

    /* Function Body */
    *pb = 0.f;
    ma = abs(*m);
    if (ma > *n) {
	return 0;
    }
    if (*n <= 0) {
	goto L10;
    } else {
	goto L30;
    }
L10:
    if (ma <= 0) {
	goto L20;
    } else {
	goto L30;
    }
L20:
    *pb = sqrt(.5f);
    goto L140;
L30:
    np1 = *n + 1;
    nmod = *n % 2;
    mmod = ma % 2;
    if (nmod <= 0) {
	goto L40;
    } else {
	goto L90;
    }
L40:
    if (mmod <= 0) {
	goto L50;
    } else {
	goto L70;
    }
L50:
    kdo = *n / 2 + 1;
    cdt = cos(*theta + *theta);
    sdt = sin(*theta + *theta);
    ct = 1.f;
    st = 0.f;
    sum = cp[1] * .5f;
    i__1 = kdo;
    for (kp1 = 2; kp1 <= i__1; ++kp1) {
	cth = cdt * ct - sdt * st;
	st = sdt * ct + cdt * st;
	ct = cth;
	sum += cp[kp1] * ct;
/* L60: */
    }
    *pb = sum;
    goto L140;
L70:
    kdo = *n / 2;
    cdt = cos(*theta + *theta);
    sdt = sin(*theta + *theta);
    ct = 1.f;
    st = 0.f;
    sum = 0.f;
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
	cth = cdt * ct - sdt * st;
	st = sdt * ct + cdt * st;
	ct = cth;
	sum += cp[k] * st;
/* L80: */
    }
    *pb = sum;
    goto L140;
L90:
    kdo = (*n + 1) / 2;
    if (mmod <= 0) {
	goto L100;
    } else {
	goto L120;
    }
L100:
    cdt = cos(*theta + *theta);
    sdt = sin(*theta + *theta);
    ct = cos(*theta);
    st = -sin(*theta);
    sum = 0.f;
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
	cth = cdt * ct - sdt * st;
	st = sdt * ct + cdt * st;
	ct = cth;
	sum += cp[k] * ct;
/* L110: */
    }
    *pb = sum;
    goto L140;
L120:
    cdt = cos(*theta + *theta);
    sdt = sin(*theta + *theta);
    ct = cos(*theta);
    st = -sin(*theta);
    sum = 0.f;
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
	cth = cdt * ct - sdt * st;
	st = sdt * ct + cdt * st;
	ct = cth;
	sum += cp[k] * st;
/* L130: */
    }
    *pb = sum;
L140:
    return 0;
} /* lfpt_ */

