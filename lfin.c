/* lfin.f -- translated by f2c (version 20061008).
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

static integer c__0 = 0;
static integer c__1 = 1;

/* subroutine lfin (init,theta,l,m,nm,pb,id,wlfin) */

/* dimension of           theta(l),  pb(id,nm+1),  wlfin(4*l*(nm+1)) */
/* arguments */

/* purpose                given m and l, routine lfin calculates */
/*                        the normalized associated legendre functions */
/*                        pbar(n,m,theta) for n=m,...,nm and theta(i) */
/*                        for i=1,...,l where */

/*                        pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m) */
/*                        /(2*factorial(n+m)))*sin(theta)**m/(2**n* */
/*                        factorial(n)) times the (n+m)th derivative of */
/*                        (x**2-1)**n with respect to x=cos(theta) */

/* usage                  call lfin (init,theta,l,m,nm,pb,id,wlfin) */

/* arguments */
/* on input               init */
/*                        = 0 */
/*                            initialization only - using parameters */
/*                            l, nm and the array theta, subroutine lfin */
/*                            initializes the array wlfin for subsequent */
/*                            use in the computation of the associated */
/*                            legendre functions pb. initialization does */
/*                            not have to be repeated unless l, nm or */
/*                            the array theta are changed. */
/*                        = 1 */
/*                            subroutine lfin uses the array wlfin that */
/*                            was computed with init = 0 to compute pb */

/*                        theta */
/*                          an array that contains the colatitudes */
/*                          at which the associated legendre functions */
/*                          will be computed. the colatitudes must be */
/*                          specified in radians. */

/*                        l */
/*                          the length of the theta array. lfin is */
/*                          vectorized with vector length l. */

/*                        m */
/*                          nonnegative integer, less than nm, specifying */
/*                          degree of pbar(n,m,theta). subroutine lfin */
/*                          must be called starting with n=0. n must be */
/*                          incremented by one in subsequent calls and */
/*                          must not exceed nm. */

/*                        nm */
/*                          the maximum value of n and m */

/*                        id */
/*                          the first dimension of the two dimensional */
/*                          array pb as it appears in the program that */
/*                          calls lfin. (see output parameter pb) */

/*                        wlfin */
/*                          an array with length 4*l*(nm+1) which */
/*                          must be initialized by calling lfin */
/*                          with init=0 (see parameter init)  it */
/*                          must not be altered between calls to */
/*                          lfin. */


/* on output              pb */
/*                          a two dimensional array with first */
/*                          dimension id in the program that calls */
/*                          lfin. the second dimension of pb must */
/*                          be at least nm+1. starting with m=0 */
/*                          lfin is called repeatedly with m being */
/*                          increased by one between calls. on each */
/*                          call, subroutine lfin computes pb(i,n+1) */
/*                          = pbar(m,n,theta(i)) for n=m,...,nm and */
/*                          i=1,...l. */

/*                        wlfin */
/*                          array containing values which must not */
/*                          be altered unless l, nm or the array theta */
/*                          are changed in which case lfin must be */
/*                          called with init=0 to reinitialize the */
/*                          wlfin array. */

/* special conditions     m must be increased by one between calls */
/*                        of lfin in which m is not zero. */

/* precision              single */

/* algorithm              routine lfin calculates pbar(n,m,theta) using */
/*                        a four term recurrence relation. (unpublished */
/*                        notes by paul n. swarztrauber) */

/* Subroutine */ int lfin_(integer *init, real *theta, integer *l, integer *m,
	 integer *nm, real *pb, integer *id, real *wlfin)
{
    static integer iw1, iw2, iw3, lnx;
    extern /* Subroutine */ int lfin1_(integer *, real *, integer *, integer *
	    , integer *, integer *, real *, real *, real *, real *, real *, 
	    real *);


/*     total length of wlfin is 4*l*(nm+1) */

    /* Parameter adjustments */
    --wlfin;
    --pb;

    /* Function Body */
    lnx = *l * (*nm + 1);
    iw1 = lnx + 1;
    iw2 = iw1 + lnx;
    iw3 = iw2 + lnx;
    lfin1_(init, theta, l, m, nm, id, &pb[1], &wlfin[1], &wlfin[iw1], &wlfin[
	    iw2], &wlfin[iw3], &wlfin[iw2]);
    return 0;
} /* lfin_ */

/* Subroutine */ int lfin1_(integer *init, real *theta, integer *l, integer *
	m, integer *nm, integer *id, real *p3, real *phz, real *ph1, real *p1,
	 real *p2, real *cp)
{
    /* System generated locals */
    integer p1_dim1, p1_offset, p2_dim1, p2_offset, p3_dim1, p3_offset, 
	    phz_dim1, phz_offset, ph1_dim1, ph1_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, n;
    static real cc, dd, ee, cn, fm, fn;
    static integer nh;
    static real tm, tn;
    static integer mp1, np1, mp3, nmp1;
    extern /* Subroutine */ int alfk_(integer *, integer *, real *);
    static real fnmm, fnpm;
    extern /* Subroutine */ int lfpt_(integer *, integer *, real *, real *, 
	    real *);
    static real temp, ssqrt2;

    /* Parameter adjustments */
    --theta;
    p2_dim1 = *l;
    p2_offset = 1 + p2_dim1;
    p2 -= p2_offset;
    p1_dim1 = *l;
    p1_offset = 1 + p1_dim1;
    p1 -= p1_offset;
    ph1_dim1 = *l;
    ph1_offset = 1 + ph1_dim1;
    ph1 -= ph1_offset;
    phz_dim1 = *l;
    phz_offset = 1 + phz_dim1;
    phz -= phz_offset;
    p3_dim1 = *id;
    p3_offset = 1 + p3_dim1;
    p3 -= p3_offset;
    --cp;

    /* Function Body */
    nmp1 = *nm + 1;
    if (*init != 0) {
	goto L5;
    }
    ssqrt2 = 1.f / sqrt(2.f);
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	phz[i__ + phz_dim1] = ssqrt2;
/* L10: */
    }
    i__1 = nmp1;
    for (np1 = 2; np1 <= i__1; ++np1) {
	nh = np1 - 1;
	alfk_(&nh, &c__0, &cp[1]);
	i__2 = *l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    lfpt_(&nh, &c__0, &theta[i__], &cp[1], &phz[i__ + np1 * phz_dim1])
		    ;
/* L16: */
	}
	alfk_(&nh, &c__1, &cp[1]);
	i__2 = *l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    lfpt_(&nh, &c__1, &theta[i__], &cp[1], &ph1[i__ + np1 * ph1_dim1])
		    ;
/* L17: */
	}
/* L15: */
    }
    return 0;
L5:
    mp1 = *m + 1;
    fm = (real) (*m);
    tm = fm + fm;
    if ((i__1 = *m - 1) < 0) {
	goto L25;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L35;
    }
L25:
    i__1 = nmp1;
    for (np1 = 1; np1 <= i__1; ++np1) {
	i__2 = *l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    p3[i__ + np1 * p3_dim1] = phz[i__ + np1 * phz_dim1];
	    p1[i__ + np1 * p1_dim1] = phz[i__ + np1 * phz_dim1];
/* L45: */
	}
    }
    return 0;
L30:
    i__2 = nmp1;
    for (np1 = 2; np1 <= i__2; ++np1) {
	i__1 = *l;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    p3[i__ + np1 * p3_dim1] = ph1[i__ + np1 * ph1_dim1];
	    p2[i__ + np1 * p2_dim1] = ph1[i__ + np1 * ph1_dim1];
/* L50: */
	}
    }
    return 0;
L35:
    temp = tm * (tm - 1.f);
    cc = sqrt((tm + 1.f) * (tm - 2.f) / temp);
    ee = sqrt(2.f / temp);
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p3[i__ + (*m + 1) * p3_dim1] = cc * p1[i__ + (*m - 1) * p1_dim1] - ee 
		* p1[i__ + (*m + 1) * p1_dim1];
/* L85: */
    }
    if (*m == *nm) {
	return 0;
    }
    temp = tm * (tm + 1.f);
    cc = sqrt((tm + 3.f) * (tm - 2.f) / temp);
    ee = sqrt(6.f / temp);
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p3[i__ + (*m + 2) * p3_dim1] = cc * p1[i__ + *m * p1_dim1] - ee * p1[
		i__ + (*m + 2) * p1_dim1];
/* L70: */
    }
    mp3 = *m + 3;
    if (nmp1 < mp3) {
	goto L80;
    }
    i__1 = nmp1;
    for (np1 = mp3; np1 <= i__1; ++np1) {
	n = np1 - 1;
	fn = (real) n;
	tn = fn + fn;
	cn = (tn + 1.f) / (tn - 3.f);
	fnpm = fn + fm;
	fnmm = fn - fm;
	temp = fnpm * (fnpm - 1.f);
	cc = sqrt(cn * (fnpm - 3.f) * (fnpm - 2.f) / temp);
	dd = sqrt(cn * fnmm * (fnmm - 1.f) / temp);
	ee = sqrt((fnmm + 1.f) * (fnmm + 2.f) / temp);
	i__2 = *l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    p3[i__ + np1 * p3_dim1] = cc * p1[i__ + (np1 - 2) * p1_dim1] + dd 
		    * p3[i__ + (np1 - 2) * p3_dim1] - ee * p1[i__ + np1 * 
		    p1_dim1];
/* L75: */
	}
    }
L80:
    i__2 = nmp1;
    for (np1 = *m; np1 <= i__2; ++np1) {
	i__1 = *l;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    p1[i__ + np1 * p1_dim1] = p2[i__ + np1 * p2_dim1];
	    p2[i__ + np1 * p2_dim1] = p3[i__ + np1 * p3_dim1];
/* L90: */
	}
    }
    return 0;
} /* lfin1_ */

