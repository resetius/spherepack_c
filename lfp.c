/* lfp.f -- translated by f2c (version 20061008).
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

/* subroutine lfp (init,n,m,l,cp,pb,w) */

/* dimension of           cp((n/2)+1), pb(l), w(5*l+41) */
/* arguments */

/* purpose                routine lfp uses coefficients computed by */
/*                        routine alfk to calculate the single precision */
/*                        normalized associated legendre function pbar(n, */
/*                        m,theta) at colatitudes theta=(i-1)*pi/(l-1), */
/*                        i=1,...,l. subroutine lfp evaluates pbar */
/*                        using one of the following trigonometric */
/*                        expansions */

/*                        1) for n even and m even, pbar(m,n,theta) = */
/*                           .5*cp(1) plus the sum from k=1 to k=n/2 */
/*                           of cp(k)*cos(2*k*th) */

/*                        2) for n even and m odd, pbar(m,n,theta) = */
/*                           the sum from k=1 to k=n/2 of */
/*                           cp(k)*sin(2*k*th) */

/*                        3) for n odd and m even, pbar(m,n,theta) = */
/*                           the sum from k=1 to k=(n+1)/2 of */
/*                           cp(k)*cos((2*k-1)*th) */

/*                        4) for n odd and m odd,  pbar(m,n,theta) = */
/*                           the sum from k=1 to k=(n+1)/2 of */
/*                           cp(k)*sin((2*k-1)*th) */


/* usage                  call lfp(init,n,m,l,cp,pb,w) */

/* arguments */

/* on input               init */
/*                          = 0 initialization only */
/*                          = 1 compute pbar(n,m,theta) */

/*                          lfp call with init = 0 initializes array w; */
/*                          no values of pbar(n,m,theta) are computed. */
/*                          init=0 should be used on the first call, or */
/*                          if l or w values differ from those in the */
/*                          previous call. */

/*                        n */
/*                          nonnegative integer, less than l, specifying */
/*                          the degree of pbar(n,m,theta) */

/*                        m */
/*                          is the order of pbar(n,m,theta). m can be */
/*                          any integer however pbar(n,m,theta) = 0 */
/*                          if abs(m) is greater than n and */
/*                          pbar(n,m,theta) = (-1)**m*pbar(n,-m,theta) */
/*                          for negative m. */

/*                        l */
/*                          number of colatitudes theta=(i-1)*pi/(l-1) */
/*                          for i=1,...,l where l is greater than 1. */
/*                          l must be an odd integer. */

/*                        cp */
/*                          single precision array of length (n/2)+1 */
/*                          containing coefficients computed by routine */
/*                          alfk */

/*                        w */
/*                          a single precision work array with at */
/*                          least 5*l+41 locations */

/* on output              pb */
/*                          single precision array of length l containing */
/*                          pbar(n,m,theta), theta=(i-1)*pi/(l-1) for i=1 */
/*                          ,...,l. */

/*                        w */
/*                          a single precision array containing values */
/*                          which must not be destroyed if the next call */
/*                          will have the same value of input parameter n */

/* special conditions     calls to routine lfp must be preceded by an */
/*                        appropriate call to routine alfk. */

/* precision              single */

/* algorithm              the trigonometric series formula used by */
/*                        routine lfp to calculate pbar(n,m,theta) for */
/*                        theta=(i-1)*pi/(l-1), i=1,...,n, depends on */
/*                        m and n as follows: */

/*                           1) for n even and m even, the formula is */
/*                              .5*cp(1) plus the sum from k=1 to k=n/2 */
/*                              of cp(k)*cos(2*k*theta) */
/*                           2) for n even and m odd. the formula is */
/*                              the sum from k=1 to k=n/2 of */
/*                              cp(k)*sin(2*k*theta) */
/*                           3) for n odd and m even, the formula is */
/*                              the sum from k=1 to k=(n+1)/2 of */
/*                              cp(k)*cos((2*k-1)*theta) */
/*                           4) for n odd and m odd, the formula is */
/*                              the sum from k=1 to k=(n+1)/2 of */
/*                              cp(k)*sin((2*k-1)*theta) */

/* accuracy               comparison between routines lfp and double */
/*                        precision dlfp on the cray1 indicates greater */
/*                        accuracy for smaller values of input parameter */
/*                        n.  agreement to 12 places was obtained for */
/*                        n=10 and to 11 places for n=100. */

/* timing                 time per call to routine lfp is dependent on */
/*                        the input parameters l and n. */

/* Subroutine */ int lfp_(integer *init, integer *n, integer *m, integer *l, 
	doublereal *cp, doublereal *pb, doublereal *w)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, ma, iw1, iw2;
    extern /* Subroutine */ int lfp1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *);


    /* Parameter adjustments */
    --w;
    --pb;
    --cp;

    /* Function Body */
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pb[i__] = 0.;
/* L10: */
    }
    ma = abs(*m);
    if (ma > *n) {
	return 0;
    }
    iw1 = *l + *l + 12;
    iw2 = iw1 + (*l + 1) * 3 / 2 + 15;
    lfp1_(init, n, &ma, l, &cp[1], &pb[1], &w[1], &w[iw1], &w[iw2]);
    return 0;
} /* lfp_ */

/* Subroutine */ int lfp1_(integer *init, integer *n, integer *m, integer *l, 
	doublereal *cp, doublereal *p, doublereal *wsave1, doublereal *wsave2,
	 doublereal *wsave3)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal), atan(doublereal);

    /* Local variables */
    integer i__;
    static integer lc;
    doublereal dt, pi;
    static integer lq, ls;
    integer lm1, np1, ls2, kdp, lmi, mmod, nmod;
    extern /* Subroutine */ int cost_(integer *, doublereal *, doublereal *), 
	    sint_(integer *, doublereal *, doublereal *), cosqb_(integer *, 
	    doublereal *, doublereal *), sinqb_(integer *, doublereal *, 
	    doublereal *), cosqi_(integer *, doublereal *), costi_(integer *, 
	    doublereal *), sinti_(integer *, doublereal *);
    doublereal ssqrt2;

    /* Parameter adjustments */
    --wsave3;
    --wsave2;
    --wsave1;
    --p;
    --cp;

    /* Function Body */
    if (*init != 0) {
	goto L41;
    }
    lc = (*l + 1) / 2;
    ls = lc - 2;
    lq = lc - 1;
    sinti_(&ls, &wsave1[1]);
    costi_(&lc, &wsave2[1]);
    cosqi_(&lq, &wsave3[1]);
    return 0;
L41:
    if (*n <= 0) {
	goto L10;
    } else {
	goto L40;
    }
L10:
    if (*m <= 0) {
	goto L20;
    } else {
	goto L40;
    }
L20:
    ssqrt2 = 1. / sqrt(2.);
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = ssqrt2;
/* L30: */
    }
    return 0;
L40:
    ls2 = (*l + 1) / 2;
    lm1 = *l - 1;
    np1 = *n + 1;
    pi = atan(1.) * 4.;
    dt = pi / lm1;
    nmod = *n % 2;
    mmod = *m % 2;
    if (nmod <= 0) {
	goto L50;
    } else {
	goto L120;
    }
L50:
    if (mmod <= 0) {
	goto L60;
    } else {
	goto L90;
    }
L60:
    kdp = *n / 2 + 1;
    i__1 = kdp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = cp[i__] * .5;
/* L70: */
    }
    p[lc] += p[lc];
    cost_(&lc, &p[1], &wsave2[1]);
    i__1 = lc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lmi = *l - i__;
	p[lmi + 1] = p[i__];
/* L80: */
    }
    goto L190;
L90:
    kdp = *n / 2;
    i__1 = kdp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__ + 1] = cp[i__] * .5;
/* L100: */
    }
    p[ls + 2] = 0.;
    sint_(&ls, &p[2], &wsave1[1]);
    i__1 = ls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lmi = *l - i__;
	p[lmi] = -p[i__ + 1];
/* L110: */
    }
    p[*l] = 0.;
    goto L190;
L120:
    kdp = (*n + 1) / 2;
    if (mmod <= 0) {
	goto L140;
    } else {
	goto L160;
    }
L140:
    i__1 = kdp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = cp[i__] * .25;
/* L130: */
    }
    cosqb_(&lq, &p[1], &wsave3[1]);
    i__1 = lq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lmi = *l - i__;
	p[lmi + 1] = -p[i__];
/* L150: */
    }
    goto L190;
L160:
    i__1 = kdp;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__ + 1] = cp[i__] * .25;
/* L180: */
    }
    sinqb_(&lq, &p[2], &wsave3[1]);
    i__1 = lq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	lmi = *l - i__;
	p[lmi] = p[i__ + 1];
/* L170: */
    }
    p[*l] = 0.;
L190:
    return 0;
} /* lfp1_ */

