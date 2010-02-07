/* alf.f -- translated by f2c (version 20061008).
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

static real c_b10 = 2.f;
static integer c__0 = 0;
static integer c__1 = 1;


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



/*     file alf.f contains subroutines alfk,lfim,lfim1,lfin,lfin1,lfpt */
/*     for computing normalized associated legendre polynomials */

/* subroutine alfk (n,m,cp) */

/* dimension of           real cp(n/2 + 1) */
/* arguments */

/* purpose                routine alfk computes single precision fourier */
/*                        coefficients in the trigonometric series */
/*                        representation of the normalized associated */
/*                        legendre function pbar(n,m,theta) for use by */
/*                        routines lfp and lfpt in calculating single */
/*                        precision pbar(n,m,theta). */

/*                        first define the normalized associated */
/*                        legendre functions */

/*                        pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m) */
/*                        /(2*factorial(n+m)))*sin(theta)**m/(2**n* */
/*                        factorial(n)) times the (n+m)th derivative of */
/*                        (x**2-1)**n with respect to x=cos(theta) */

/*                        where theta is colatitude. */

/*                        then subroutine alfk computes the coefficients */
/*                        cp(k) in the following trigonometric */
/*                        expansion of pbar(m,n,theta). */

/*                        1) for n even and m even, pbar(m,n,theta) = */
/*                           .5*cp(1) plus the sum from k=1 to k=n/2 */
/*                           of cp(k+1)*cos(2*k*th) */

/*                        2) for n even and m odd, pbar(m,n,theta) = */
/*                           the sum from k=1 to k=n/2 of */
/*                           cp(k)*sin(2*k*th) */

/*                        3) for n odd and m even, pbar(m,n,theta) = */
/*                           the sum from k=1 to k=(n+1)/2 of */
/*                           cp(k)*cos((2*k-1)*th) */

/*                        4) for n odd and m odd,  pbar(m,n,theta) = */
/*                           the sum from k=1 to k=(n+1)/2 of */
/*                           cp(k)*sin((2*k-1)*th) */


/* usage                  call alfk(n,m,cp) */

/* arguments */

/* on input               n */
/*                          nonnegative integer specifying the degree of */
/*                          pbar(n,m,theta) */

/*                        m */
/*                          is the order of pbar(n,m,theta). m can be */
/*                          any integer however cp is computed such that */
/*                          pbar(n,m,theta) = 0 if abs(m) is greater */
/*                          than n and pbar(n,m,theta) = (-1)**m* */
/*                          pbar(n,-m,theta) for negative m. */

/* on output              cp */
/*                          single precision array of length (n/2)+1 */
/*                          which contains the fourier coefficients in */
/*                          the trigonometric series representation of */
/*                          pbar(n,m,theta) */


/* special conditions     none */

/* precision              single */

/* algorithm              the highest order coefficient is determined in */
/*                        closed form and the remainig coefficients are */
/*                        determined as the solution of a backward */
/*                        recurrence relation. */

/* accuracy               comparison between routines alfk and double */
/*                        precision dalfk on the cray1 indicates */
/*                        greater accuracy for smaller values */
/*                        of input parameter n.  agreement to 14 */
/*                        places was obtained for n=10 and to 13 */
/*                        places for n=100. */

/* Subroutine */ int alfk_(integer *n, integer *m, real *cp)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_ri(real *, integer *);

    /* Local variables */
    static integer i__, l;
    static real a1, b1, c1, t1, t2;
    static integer ma;
    static real fk, cp2, pm1;
    static integer nex;
    static real fden, fnmh, fnum, fnnp1;
    static integer nmms2;
    static real fnmsq;


    /* Parameter adjustments */
    --cp;

    /* Function Body */
    cp[1] = 0.f;
    ma = abs(*m);
    if (ma > *n) {
	return 0;
    }
    if ((i__1 = *n - 1) < 0) {
	goto L2;
    } else if (i__1 == 0) {
	goto L3;
    } else {
	goto L5;
    }
L2:
    cp[1] = sqrt(2.f);
    return 0;
L3:
    if (ma != 0) {
	goto L4;
    }
    cp[1] = sqrt(1.5f);
    return 0;
L4:
    cp[1] = sqrt(.75f);
    if (*m == -1) {
	cp[1] = -cp[1];
    }
    return 0;
L5:
    if ((*n + ma) % 2 != 0) {
	goto L10;
    }
    nmms2 = (*n - ma) / 2;
    fnum = (real) (*n + ma + 1);
    fnmh = (real) (*n - ma + 1);
    pm1 = 1.f;
    goto L15;
L10:
    nmms2 = (*n - ma - 1) / 2;
    fnum = (real) (*n + ma + 2);
    fnmh = (real) (*n - ma + 2);
    pm1 = -1.f;
L15:
    t1 = 9.5367431640625e-7f;
    nex = 20;
    fden = 2.f;
    if (nmms2 < 1) {
	goto L20;
    }
    i__1 = nmms2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t1 = fnum * t1 / fden;
	if (t1 > 1048576.f) {
	    t1 /= 1099511627776.f;
	    nex += 40;
	}
	fnum += 2.f;
	fden += 2.f;
/* L18: */
    }
L20:
    i__1 = *n - 1 - nex;
    t1 /= pow_ri(&c_b10, &i__1);
    if (ma / 2 % 2 != 0) {
	t1 = -t1;
    }
    t2 = 1.f;
    if (ma == 0) {
	goto L26;
    }
    i__1 = ma;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t2 = fnmh * t2 / (fnmh + pm1);
	fnmh += 2.f;
/* L25: */
    }
L26:
    cp2 = t1 * sqrt((*n + .5f) * t2);
    fnnp1 = (real) (*n * (*n + 1));
    fnmsq = fnnp1 - ma * 2.f * ma;
    l = (*n + 1) / 2;
    if (*n % 2 == 0 && ma % 2 == 0) {
	++l;
    }
    cp[l] = cp2;
    if (*m >= 0) {
	goto L29;
    }
    if (ma % 2 != 0) {
	cp[l] = -cp[l];
    }
L29:
    if (l <= 1) {
	return 0;
    }
    fk = (real) (*n);
    a1 = (fk - 2.f) * (fk - 1.f) - fnnp1;
    b1 = (fk * fk - fnmsq) * 2.f;
    cp[l - 1] = b1 * cp[l] / a1;
L30:
    --l;
    if (l <= 1) {
	return 0;
    }
    fk += -2.f;
    a1 = (fk - 2.f) * (fk - 1.f) - fnnp1;
    b1 = (fk * fk - fnmsq) * -2.f;
    c1 = (fk + 1.f) * (fk + 2.f) - fnnp1;
    cp[l - 1] = -(b1 * cp[l] + c1 * cp[l + 1]) / a1;
    goto L30;
} /* alfk_ */

/* subroutine lfim (init,theta,l,n,nm,pb,id,wlfim) */

/* dimension of           theta(l),  pb(id,nm+1),  wlfim(4*l*(nm+1)) */
/* arguments */

/* purpose                given n and l, routine lfim calculates */
/*                        the normalized associated legendre functions */
/*                        pbar(n,m,theta) for m=0,...,n and theta(i) */
/*                        for i=1,...,l where */

/*                        pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m) */
/*                        /(2*factorial(n+m)))*sin(theta)**m/(2**n* */
/*                        factorial(n)) times the (n+m)th derivative of */
/*                        (x**2-1)**n with respect to x=cos(theta) */

/* usage                  call lfim (init,theta,l,n,nm,pb,id,wlfim) */

/* arguments */
/* on input               init */
/*                        = 0 */
/*                            initialization only - using parameters */
/*                            l, nm and array theta, subroutine lfim */
/*                            initializes array wlfim for subsequent */
/*                            use in the computation of the associated */
/*                            legendre functions pb. initialization */
/*                            does not have to be repeated unless */
/*                            l, nm, or array theta are changed. */
/*                        = 1 */
/*                            subroutine lfim uses the array wlfim that */
/*                            was computed with init = 0 to compute pb. */

/*                        theta */
/*                          an array that contains the colatitudes */
/*                          at which the associated legendre functions */
/*                          will be computed. the colatitudes must be */
/*                          specified in radians. */

/*                        l */
/*                          the length of the theta array. lfim is */
/*                          vectorized with vector length l. */

/*                        n */
/*                          nonnegative integer, less than nm, specifying */
/*                          degree of pbar(n,m,theta). subroutine lfim */
/*                          must be called starting with n=0. n must be */
/*                          incremented by one in subsequent calls and */
/*                          must not exceed nm. */

/*                        nm */
/*                          the maximum value of n and m */

/*                        id */
/*                          the first dimension of the two dimensional */
/*                          array pb as it appears in the program that */
/*                          calls lfim. (see output parameter pb) */

/*                        wlfim */
/*                          an array with length 4*l*(nm+1) which */
/*                          must be initialized by calling lfim */
/*                          with init=0 (see parameter init)  it */
/*                          must not be altered between calls to */
/*                          lfim. */


/* on output              pb */
/*                          a two dimensional array with first */
/*                          dimension id in the program that calls */
/*                          lfim. the second dimension of pb must */
/*                          be at least nm+1. starting with n=0 */
/*                          lfim is called repeatedly with n being */
/*                          increased by one between calls. on each */
/*                          call, subroutine lfim computes */
/*                          = pbar(m,n,theta(i)) for m=0,...,n and */
/*                          i=1,...l. */

/*                        wlfim */
/*                          array containing values which must not */
/*                          be altered unless l, nm or the array theta */
/*                          are changed in which case lfim must be */
/*                          called with init=0 to reinitialize the */
/*                          wlfim array. */

/* special conditions     n must be increased by one between calls */
/*                        of lfim in which n is not zero. */

/* precision              single */


/* algorithm              routine lfim calculates pbar(n,m,theta) using */
/*                        a four term recurrence relation. (unpublished */
/*                        notes by paul n. swarztrauber) */

/* Subroutine */ int lfim_(integer *init, real *theta, integer *l, integer *n,
	 integer *nm, real *pb, integer *id, real *wlfim)
{
    static integer iw1, iw2, iw3, lnx;
    extern /* Subroutine */ int lfim1_(integer *, real *, integer *, integer *
	    , integer *, integer *, real *, real *, real *, real *, real *, 
	    real *);


/*     total length of wlfim is 4*l*(nm+1) */

    /* Parameter adjustments */
    --wlfim;
    --pb;

    /* Function Body */
    lnx = *l * (*nm + 1);
    iw1 = lnx + 1;
    iw2 = iw1 + lnx;
    iw3 = iw2 + lnx;
    lfim1_(init, theta, l, n, nm, id, &pb[1], &wlfim[1], &wlfim[iw1], &wlfim[
	    iw2], &wlfim[iw3], &wlfim[iw2]);
    return 0;
} /* lfim_ */

/* Subroutine */ int lfim1_(integer *init, real *theta, integer *l, integer *
	n, integer *nm, integer *id, real *p3, real *phz, real *ph1, real *p1,
	 real *p2, real *cp)
{
    /* System generated locals */
    integer p1_dim1, p1_offset, p2_dim1, p2_offset, p3_dim1, p3_offset, 
	    phz_dim1, phz_offset, ph1_dim1, ph1_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, m;
    static real cc, dd, ee, cn, fn, fm;
    static integer nh;
    static real tn;
    static integer nm1, mp1, np1, nmp1;
    extern /* Subroutine */ int alfk_(integer *, integer *, real *);
    static real fnmm, fnpm;
    extern /* Subroutine */ int lfpt_(integer *, integer *, real *, real *, 
	    real *);
    static real temp, sq1s6, sq5s6, ssqrt2;

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
    if (*n > 2) {
	goto L60;
    }
    if ((i__1 = *n - 1) < 0) {
	goto L25;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L35;
    }
L25:
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p3[i__ + p3_dim1] = phz[i__ + phz_dim1];
/* L45: */
    }
    return 0;
L30:
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p3[i__ + p3_dim1] = phz[i__ + (phz_dim1 << 1)];
	p3[i__ + (p3_dim1 << 1)] = ph1[i__ + (ph1_dim1 << 1)];
/* L50: */
    }
    return 0;
L35:
    sq5s6 = sqrt(.83333333333333337f);
    sq1s6 = sqrt(.16666666666666666f);
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p3[i__ + p3_dim1] = phz[i__ + phz_dim1 * 3];
	p3[i__ + (p3_dim1 << 1)] = ph1[i__ + ph1_dim1 * 3];
	p3[i__ + p3_dim1 * 3] = sq5s6 * phz[i__ + phz_dim1] - sq1s6 * p3[i__ 
		+ p3_dim1];
	p1[i__ + p1_dim1] = phz[i__ + (phz_dim1 << 1)];
	p1[i__ + (p1_dim1 << 1)] = ph1[i__ + (ph1_dim1 << 1)];
	p2[i__ + p2_dim1] = phz[i__ + phz_dim1 * 3];
	p2[i__ + (p2_dim1 << 1)] = ph1[i__ + ph1_dim1 * 3];
	p2[i__ + p2_dim1 * 3] = p3[i__ + p3_dim1 * 3];
/* L55: */
    }
    return 0;
L60:
    nm1 = *n - 1;
    np1 = *n + 1;
    fn = (real) (*n);
    tn = fn + fn;
    cn = (tn + 1.f) / (tn - 3.f);
    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p3[i__ + p3_dim1] = phz[i__ + np1 * phz_dim1];
	p3[i__ + (p3_dim1 << 1)] = ph1[i__ + np1 * ph1_dim1];
/* L65: */
    }
    if (nm1 < 3) {
	goto L71;
    }
    i__1 = nm1;
    for (mp1 = 3; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	fm = (real) m;
	fnpm = fn + fm;
	fnmm = fn - fm;
	temp = fnpm * (fnpm - 1.f);
	cc = sqrt(cn * (fnpm - 3.f) * (fnpm - 2.f) / temp);
	dd = sqrt(cn * fnmm * (fnmm - 1.f) / temp);
	ee = sqrt((fnmm + 1.f) * (fnmm + 2.f) / temp);
	i__2 = *l;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    p3[i__ + mp1 * p3_dim1] = cc * p1[i__ + (mp1 - 2) * p1_dim1] + dd 
		    * p1[i__ + mp1 * p1_dim1] - ee * p3[i__ + (mp1 - 2) * 
		    p3_dim1];
/* L70: */
	}
    }
L71:
    fnpm = fn + fn - 1.f;
    temp = fnpm * (fnpm - 1.f);
    cc = sqrt(cn * (fnpm - 3.f) * (fnpm - 2.f) / temp);
    ee = sqrt(6.f / temp);
    i__2 = *l;
    for (i__ = 1; i__ <= i__2; ++i__) {
	p3[i__ + *n * p3_dim1] = cc * p1[i__ + (*n - 2) * p1_dim1] - ee * p3[
		i__ + (*n - 2) * p3_dim1];
/* L75: */
    }
    fnpm = fn + fn;
    temp = fnpm * (fnpm - 1.f);
    cc = sqrt(cn * (fnpm - 3.f) * (fnpm - 2.f) / temp);
    ee = sqrt(2.f / temp);
    i__2 = *l;
    for (i__ = 1; i__ <= i__2; ++i__) {
	p3[i__ + (*n + 1) * p3_dim1] = cc * p1[i__ + (*n - 1) * p1_dim1] - ee 
		* p3[i__ + (*n - 1) * p3_dim1];
/* L80: */
    }
    i__2 = np1;
    for (mp1 = 1; mp1 <= i__2; ++mp1) {
	i__1 = *l;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    p1[i__ + mp1 * p1_dim1] = p2[i__ + mp1 * p2_dim1];
	    p2[i__ + mp1 * p2_dim1] = p3[i__ + mp1 * p3_dim1];
/* L90: */
	}
    }
    return 0;
} /* lfim1_ */

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

