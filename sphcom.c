/* sphcom.f -- translated by f2c (version 20061008).
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

static doublereal c_b10 = 2.;
static integer c__0 = 0;


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


/* ... file sphcom.f */

/*     this file must be loaded with all driver level files */
/*     in spherepack3.0.  it includes undocumented subroutines */
/*     called by some or all of the drivers */

/* Subroutine */ int dnlfk_(integer *m, integer *n, doublereal *cp)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__, l;
    static doublereal a1, b1, c1, t1, t2;
    static integer ma;
    static doublereal fk, cp2, pm1;
    static integer nex;
    static doublereal fden, fnmh, fnum, fnnp1;
    static integer nmms2;
    static doublereal fnmsq;


/*     cp requires n/2+1 double precision locations */


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
    cp[1] = sqrt(2.);
    return 0;
L3:
    if (ma != 0) {
	goto L4;
    }
    cp[1] = sqrt(1.5);
    return 0;
L4:
    cp[1] = sqrt(.75);
    if (*m == -1) {
	cp[1] = -cp[1];
    }
    return 0;
L5:
    if ((*n + ma) % 2 != 0) {
	goto L10;
    }
    nmms2 = (*n - ma) / 2;
    fnum = (doublereal) (*n + ma + 1);
    fnmh = (doublereal) (*n - ma + 1);
    pm1 = 1.;
    goto L15;
L10:
    nmms2 = (*n - ma - 1) / 2;
    fnum = (doublereal) (*n + ma + 2);
    fnmh = (doublereal) (*n - ma + 2);
    pm1 = -1.;
/*      t1 = 1. */
/*      t1 = 2.d0**(n-1) */
/*      t1 = 1.d0/t1 */
L15:
    t1 = 9.5367431640625e-7;
    nex = 20;
    fden = 2.;
    if (nmms2 < 1) {
	goto L20;
    }
    i__1 = nmms2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t1 = fnum * t1 / fden;
	if (t1 > 1048576.) {
	    t1 /= 1099511627776.;
	    nex += 40;
	}
	fnum += 2.f;
	fden += 2.f;
/* L18: */
    }
L20:
    i__1 = *n - 1 - nex;
    t1 /= pow_di(&c_b10, &i__1);
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
    cp2 = t1 * sqrt((*n + .5) * t2);
    fnnp1 = (doublereal) (*n * (*n + 1));
    fnmsq = fnnp1 - ma * 2. * ma;
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
    fk = (doublereal) (*n);
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
} /* dnlfk_ */

/* Subroutine */ int dnlft_(integer *m, integer *n, doublereal *theta, 
	doublereal *cp, doublereal *pb)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal chh, cdt;
    static integer kdo;
    static doublereal cth, sdt, sth;
    static integer mmod, nmod;

    /* Parameter adjustments */
    --cp;

    /* Function Body */
    cdt = cos(*theta + *theta);
    sdt = sin(*theta + *theta);
    nmod = *n % 2;
    mmod = *m % 2;
    if (nmod <= 0) {
	goto L1;
    } else {
	goto L2;
    }
L1:
    if (mmod <= 0) {
	goto L3;
    } else {
	goto L4;
    }

/*     n even, m even */

L3:
    kdo = *n / 2;
    *pb = cp[1] * .5f;
    if (*n == 0) {
	return 0;
    }
    cth = cdt;
    sth = sdt;
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k+1)*dcos(2*k*theta) */
	*pb += cp[k + 1] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L170: */
    }
    return 0;

/*     n even, m odd */

L4:
    kdo = *n / 2;
    *pb = 0.f;
    cth = cdt;
    sth = sdt;
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k)*dsin(2*k*theta) */
	*pb += cp[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L180: */
    }
    return 0;
L2:
    if (mmod <= 0) {
	goto L13;
    } else {
	goto L14;
    }

/*     n odd, m even */

L13:
    kdo = (*n + 1) / 2;
    *pb = 0.f;
    cth = cos(*theta);
    sth = sin(*theta);
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k)*dcos((2*k-1)*theta) */
	*pb += cp[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L190: */
    }
    return 0;

/*     n odd, m odd */

L14:
    kdo = (*n + 1) / 2;
    *pb = 0.f;
    cth = cos(*theta);
    sth = sin(*theta);
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k)*dsin((2*k-1)*theta) */
	*pb += cp[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L200: */
    }
    return 0;
} /* dnlft_ */

/* Subroutine */ int dnlftd_(integer *m, integer *n, doublereal *theta, 
	doublereal *cp, doublereal *pb)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal chh, cdt;
    static integer kdo;
    static doublereal cth, sdt, sth;
    static integer mmod, nmod;


/*     computes the derivative of pmn(theta) with respect to theta */

    /* Parameter adjustments */
    --cp;

    /* Function Body */
    cdt = cos(*theta + *theta);
    sdt = sin(*theta + *theta);
    nmod = *n % 2;
    mmod = abs(*m) % 2;
    if (nmod <= 0) {
	goto L1;
    } else {
	goto L2;
    }
L1:
    if (mmod <= 0) {
	goto L3;
    } else {
	goto L4;
    }

/*     n even, m even */

L3:
    kdo = *n / 2;
    *pb = 0.;
    if (*n == 0) {
	return 0;
    }
    cth = cdt;
    sth = sdt;
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k+1)*dcos(2*k*theta) */
	*pb -= k * 2. * cp[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L170: */
    }
    return 0;

/*     n even, m odd */

L4:
    kdo = *n / 2;
    *pb = 0.f;
    cth = cdt;
    sth = sdt;
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k)*dsin(2*k*theta) */
	*pb += k * 2. * cp[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L180: */
    }
    return 0;
L2:
    if (mmod <= 0) {
	goto L13;
    } else {
	goto L14;
    }

/*     n odd, m even */

L13:
    kdo = (*n + 1) / 2;
    *pb = 0.f;
    cth = cos(*theta);
    sth = sin(*theta);
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k)*dcos((2*k-1)*theta) */
	*pb -= (k * 2. - 1) * cp[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L190: */
    }
    return 0;

/*     n odd, m odd */

L14:
    kdo = (*n + 1) / 2;
    *pb = 0.f;
    cth = cos(*theta);
    sth = sin(*theta);
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k)*dsin((2*k-1)*theta) */
	*pb += (k * 2. - 1) * cp[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L200: */
    }
    return 0;
} /* dnlftd_ */

/* Subroutine */ int legin_(integer *mode, integer *l, integer *nlat, integer 
	*m, real *w, real *pmn, integer *km)
{
    static integer i1, i2, i3, i4, i5, late;
    extern /* Subroutine */ int legin1_(integer *, integer *, integer *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *, integer *);

/*     this subroutine computes legendre polynomials for n=m,...,l-1 */
/*     and  i=1,...,late (late=((nlat+mod(nlat,2))/2)gaussian grid */
/*     in pmn(n+1,i,km) using swarztrauber's recursion formula. */
/*     the vector w contains quantities precomputed in shigc. */
/*     legin must be called in the order m=0,1,...,l-1 */
/*     (e.g., if m=10 is sought it must be preceded by calls with */
/*     m=0,1,2,...,9 in that order) */
/*     set size of pole to equator gaussian grid */
    /* Parameter adjustments */
    --pmn;
    --w;

    /* Function Body */
    late = (*nlat + *nlat % 2) / 2;
/*     partition w (set pointers for p0n,p1n,abel,bbel,cbel,pmn) */
    i1 = *nlat + 1;
    i2 = i1 + *nlat * late;
    i3 = i2 + *nlat * late;
    i4 = i3 + ((*nlat << 1) - *l) * (*l - 1) / 2;
    i5 = i4 + ((*nlat << 1) - *l) * (*l - 1) / 2;
    legin1_(mode, l, nlat, &late, m, &w[i1], &w[i2], &w[i3], &w[i4], &w[i5], &
	    pmn[1], km);
    return 0;
} /* legin_ */

/* Subroutine */ int legin1_(integer *mode, integer *l, integer *nlat, 
	integer *late, integer *m, real *p0n, real *p1n, real *abel, real *
	bbel, real *cbel, real *pmn, integer *km)
{
    /* Initialized data */

    static integer km0 = 1;
    static integer km1 = 2;
    static integer km2 = 3;

    /* System generated locals */
    integer p0n_dim1, p0n_offset, p1n_dim1, p1n_offset, pmn_dim1, pmn_dim2, 
	    pmn_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, n, ms, np1, imn, kmt, ninc;

    /* Parameter adjustments */
    pmn_dim1 = *nlat;
    pmn_dim2 = *late;
    pmn_offset = 1 + pmn_dim1 * (1 + pmn_dim2);
    pmn -= pmn_offset;
    p1n_dim1 = *nlat;
    p1n_offset = 1 + p1n_dim1;
    p1n -= p1n_offset;
    p0n_dim1 = *nlat;
    p0n_offset = 1 + p0n_dim1;
    p0n -= p0n_offset;
    --abel;
    --bbel;
    --cbel;

    /* Function Body */
/*     define index function used in storing triangular */
/*     arrays for recursion coefficients (functions of (m,n)) */
/*     for 2.le.m.le.n-1 and 2.le.n.le.l-1 */
/*     for l.le.n.le.nlat and 2.le.m.le.l */
/*     set do loop indices for full or half sphere */
    ms = *m + 1;
    ninc = 1;
    if (*mode == 1) {
/*     only compute pmn for n-m odd */
	ms = *m + 2;
	ninc = 2;
    } else if (*mode == 2) {
/*     only compute pmn for n-m even */
	ms = *m + 1;
	ninc = 2;
    }
    if (*m > 1) {
	i__1 = *nlat;
	i__2 = ninc;
	for (np1 = ms; i__2 < 0 ? np1 >= i__1 : np1 <= i__1; np1 += i__2) {
	    n = np1 - 1;
	    imn = (n - 1) * (n - 2) / 2 + *m - 1;
	    if (n >= *l) {
		imn = *l * (*l - 1) / 2 + (n - *l - 1) * (*l - 1) + *m - 1;
	    }
	    i__3 = *late;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		pmn[np1 + (i__ + km0 * pmn_dim2) * pmn_dim1] = abel[imn] * 
			pmn[n - 1 + (i__ + km2 * pmn_dim2) * pmn_dim1] + bbel[
			imn] * pmn[n - 1 + (i__ + km0 * pmn_dim2) * pmn_dim1] 
			- cbel[imn] * pmn[np1 + (i__ + km2 * pmn_dim2) * 
			pmn_dim1];
/* L100: */
	    }
	}
    } else if (*m == 0) {
	i__3 = *nlat;
	i__2 = ninc;
	for (np1 = ms; i__2 < 0 ? np1 >= i__3 : np1 <= i__3; np1 += i__2) {
	    i__1 = *late;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		pmn[np1 + (i__ + km0 * pmn_dim2) * pmn_dim1] = p0n[np1 + i__ *
			 p0n_dim1];
/* L101: */
	    }
	}
    } else if (*m == 1) {
	i__1 = *nlat;
	i__2 = ninc;
	for (np1 = ms; i__2 < 0 ? np1 >= i__1 : np1 <= i__1; np1 += i__2) {
	    i__3 = *late;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		pmn[np1 + (i__ + km0 * pmn_dim2) * pmn_dim1] = p1n[np1 + i__ *
			 p1n_dim1];
/* L102: */
	    }
	}
    }
/*     permute column indices */
/*     km0,km1,km2 store m,m-1,m-2 columns */
    kmt = km0;
    km0 = km2;
    km2 = km1;
    km1 = kmt;
/*     set current m index in output param km */
    *km = kmt;
    return 0;
} /* legin1_ */

/* Subroutine */ int zfin_(integer *isym, integer *nlat, integer *nlon, 
	integer *m, real *z__, integer *i3, real *wzfin)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, iw3, iw4, lim, labc, imid, mmax;
    extern /* Subroutine */ int zfin1_(integer *, integer *, integer *, real *
	    , integer *, integer *, real *, real *, real *, real *, real *);

    /* Parameter adjustments */
    --wzfin;
    --z__;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    lim = *nlat * imid;
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mmax = min(i__1,i__2);
    labc = (mmax - 2) * (*nlat + *nlat - mmax - 1) / 2;
    iw1 = lim + 1;
    iw2 = iw1 + lim;
    iw3 = iw2 + labc;
    iw4 = iw3 + labc;

/*     the length of wzfin is 2*lim+3*labc */

    zfin1_(isym, nlat, m, &z__[1], &imid, i3, &wzfin[1], &wzfin[iw1], &wzfin[
	    iw2], &wzfin[iw3], &wzfin[iw4]);
    return 0;
} /* zfin_ */

/* Subroutine */ int zfin1_(integer *isym, integer *nlat, integer *m, real *
	z__, integer *imid, integer *i3, real *zz, real *z1, real *a, real *b,
	 real *c__)
{
    /* System generated locals */
    integer z_dim1, z_dim2, z_offset, zz_dim1, zz_offset, z1_dim1, z1_offset, 
	    i__1, i__2, i__3;

    /* Local variables */
    static integer i__, i1, i2, ns, np1, nstp, ihold, nstrt;

    /* Parameter adjustments */
    z1_dim1 = *imid;
    z1_offset = 1 + z1_dim1;
    z1 -= z1_offset;
    zz_dim1 = *imid;
    zz_offset = 1 + zz_dim1;
    zz -= zz_offset;
    z_dim1 = *imid;
    z_dim2 = *nlat;
    z_offset = 1 + z_dim1 * (1 + z_dim2);
    z__ -= z_offset;
    --a;
    --b;
    --c__;

    /* Function Body */
    ihold = i1;
    i1 = i2;
    i2 = *i3;
    *i3 = ihold;
    if ((i__1 = *m - 1) < 0) {
	goto L25;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L35;
    }
L25:
    i1 = 1;
    i2 = 2;
    *i3 = 3;
    i__1 = *nlat;
    for (np1 = 1; np1 <= i__1; ++np1) {
	i__2 = *imid;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__ + (np1 + *i3 * z_dim2) * z_dim1] = zz[i__ + np1 * zz_dim1]
		    ;
/* L45: */
	}
    }
    return 0;
L30:
    i__2 = *nlat;
    for (np1 = 2; np1 <= i__2; ++np1) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    z__[i__ + (np1 + *i3 * z_dim2) * z_dim1] = z1[i__ + np1 * z1_dim1]
		    ;
/* L50: */
	}
    }
    return 0;
L35:
    ns = (*m - 2) * (*nlat + *nlat - *m - 1) / 2 + 1;
    if (*isym == 1) {
	goto L36;
    }
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__ + (*m + 1 + *i3 * z_dim2) * z_dim1] = a[ns] * z__[i__ + (*m - 
		1 + i1 * z_dim2) * z_dim1] - c__[ns] * z__[i__ + (*m + 1 + i1 
		* z_dim2) * z_dim1];
/* L85: */
    }
L36:
    if (*m == *nlat - 1) {
	return 0;
    }
    if (*isym == 2) {
	goto L71;
    }
    ++ns;
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__ + (*m + 2 + *i3 * z_dim2) * z_dim1] = a[ns] * z__[i__ + (*m + 
		i1 * z_dim2) * z_dim1] - c__[ns] * z__[i__ + (*m + 2 + i1 * 
		z_dim2) * z_dim1];
/* L70: */
    }
L71:
    nstrt = *m + 3;
    if (*isym == 1) {
	nstrt = *m + 4;
    }
    if (nstrt > *nlat) {
	goto L80;
    }
    nstp = 2;
    if (*isym == 0) {
	nstp = 1;
    }
    i__1 = *nlat;
    i__2 = nstp;
    for (np1 = nstrt; i__2 < 0 ? np1 >= i__1 : np1 <= i__1; np1 += i__2) {
	ns += nstp;
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    z__[i__ + (np1 + *i3 * z_dim2) * z_dim1] = a[ns] * z__[i__ + (np1 
		    - 2 + i1 * z_dim2) * z_dim1] + b[ns] * z__[i__ + (np1 - 2 
		    + *i3 * z_dim2) * z_dim1] - c__[ns] * z__[i__ + (np1 + i1 
		    * z_dim2) * z_dim1];
/* L75: */
	}
    }
L80:
    return 0;
} /* zfin1_ */

/* Subroutine */ int zfinit_(integer *nlat, integer *nlon, real *wzfin, 
	doublereal *dwork)
{
    static integer iw1, imid;
    extern /* Subroutine */ int zfini1_(integer *, integer *, integer *, real 
	    *, real *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --dwork;
    --wzfin;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     the length of wzfin is 3*((l-3)*l+2)/2 + 2*l*imid */
/*     the length of dwork is nlat+2 */

    zfini1_(nlat, nlon, &imid, &wzfin[1], &wzfin[iw1], &dwork[1], &dwork[*
	    nlat / 2 + 1]);
    return 0;
} /* zfinit_ */

/* Subroutine */ int zfini1_(integer *nlat, integer *nlon, integer *imid, 
	real *z__, real *abc, doublereal *cz, doublereal *work)
{
    /* System generated locals */
    integer z_dim1, z_dim2, z_offset, i__1, i__2;

    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static integer i__, m, n;
    static doublereal dt, pi, th, zh;
    static integer mp1, np1;
    extern /* Subroutine */ int rabcp_(integer *, integer *, real *), dnzfk_(
	    integer *, integer *, integer *, doublereal *, doublereal *), 
	    dnzft_(integer *, integer *, integer *, doublereal *, doublereal *
	    , doublereal *);


/*     abc must have 3*((mmax-2)*(nlat+nlat-mmax-1))/2 locations */
/*     where mmax = min0(nlat,nlon/2+1) */
/*     cz and work must each have nlat+1 locations */

    /* Parameter adjustments */
    z_dim1 = *imid;
    z_dim2 = *nlat;
    z_offset = 1 + z_dim1 * (1 + z_dim2);
    z__ -= z_offset;
    --abc;
    --cz;
    --work;

    /* Function Body */
    pi = atan(1.) * 4.f;
    dt = pi / (*nlat - 1);
    for (mp1 = 1; mp1 <= 2; ++mp1) {
	m = mp1 - 1;
	i__1 = *nlat;
	for (np1 = mp1; np1 <= i__1; ++np1) {
	    n = np1 - 1;
	    dnzfk_(nlat, &m, &n, &cz[1], &work[1]);
	    i__2 = *imid;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		th = (i__ - 1) * dt;
		dnzft_(nlat, &m, &n, &th, &cz[1], &zh);
		z__[i__ + (np1 + mp1 * z_dim2) * z_dim1] = zh;
/* L165: */
	    }
	    z__[(np1 + mp1 * z_dim2) * z_dim1 + 1] *= .5f;
/* L160: */
	}
    }
    rabcp_(nlat, nlon, &abc[1]);
    return 0;
} /* zfini1_ */

/* Subroutine */ int dnzfk_(integer *nlat, integer *m, integer *n, doublereal 
	*cz, doublereal *work)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k;
    static doublereal t1, t2;
    static integer lc;
    static doublereal sc1;
    static integer kp1, kdo, idx;
    static doublereal sum;
    static integer mmod, nmod;
    extern /* Subroutine */ int dnlfk_(integer *, integer *, doublereal *);


/*     dnzfk computes the coefficients in the trigonometric */
/*     expansion of the z functions that are used in spherical */
/*     harmonic analysis. */


/*     cz and work must both have nlat/2+1 locations */

    /* Parameter adjustments */
    --work;
    --cz;

    /* Function Body */
    lc = (*nlat + 1) / 2;
    sc1 = 2. / (real) (*nlat - 1);
    dnlfk_(m, n, &work[1]);
    nmod = *n % 2;
    mmod = *m % 2;
    if (nmod <= 0) {
	goto L1;
    } else {
	goto L2;
    }
L1:
    if (mmod <= 0) {
	goto L3;
    } else {
	goto L4;
    }

/*     n even, m even */

L3:
    kdo = *n / 2 + 1;
    i__1 = lc;
    for (idx = 1; idx <= i__1; ++idx) {
	i__ = idx + idx - 2;
	sum = work[1] / (1. - i__ * i__);
	if (kdo < 2) {
	    goto L29;
	}
	i__2 = kdo;
	for (kp1 = 2; kp1 <= i__2; ++kp1) {
	    k = kp1 - 1;
/* Computing 2nd power */
	    i__3 = k + k + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - i__;
	    t2 = 1. - i__3 * i__3;
/* L8: */
	    sum += work[kp1] * (t1 + t2) / (t1 * t2);
/* L6: */
	}
L29:
	cz[idx] = sc1 * sum;
/* L5: */
    }
    return 0;

/*     n even, m odd */

L4:
    kdo = *n / 2;
    i__1 = lc;
    for (idx = 1; idx <= i__1; ++idx) {
	i__ = idx + idx - 2;
	sum = 0.f;
	i__2 = kdo;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    i__3 = k + k + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - i__;
	    t2 = 1. - i__3 * i__3;
/* L12: */
	    sum += work[k] * (t1 - t2) / (t1 * t2);
/* L101: */
	}
	cz[idx] = sc1 * sum;
/* L9: */
    }
    return 0;
L2:
    if (mmod <= 0) {
	goto L13;
    } else {
	goto L14;
    }

/*     n odd, m even */

L13:
    kdo = (*n + 1) / 2;
    i__1 = lc;
    for (idx = 1; idx <= i__1; ++idx) {
	i__ = idx + idx - 1;
	sum = 0.f;
	i__2 = kdo;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    i__3 = k + k - 1 + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - 1 - i__;
	    t2 = 1. - i__3 * i__3;
/* L18: */
	    sum += work[k] * (t1 + t2) / (t1 * t2);
/* L16: */
	}
	cz[idx] = sc1 * sum;
/* L15: */
    }
    return 0;

/*     n odd, m odd */

L14:
    kdo = (*n + 1) / 2;
    i__1 = lc;
    for (idx = 1; idx <= i__1; ++idx) {
	i__ = idx + idx - 3;
	sum = 0.f;
	i__2 = kdo;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    i__3 = k + k - 1 + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - 1 - i__;
	    t2 = 1. - i__3 * i__3;
/* L22: */
	    sum += work[k] * (t1 - t2) / (t1 * t2);
/* L20: */
	}
	cz[idx] = sc1 * sum;
/* L19: */
    }
    return 0;
} /* dnzfk_ */

/* Subroutine */ int dnzft_(integer *nlat, integer *m, integer *n, doublereal 
	*th, doublereal *cz, doublereal *zh)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k, lc, lq, ls;
    static doublereal chh, cdt, cth, sdt, sth;
    static integer lmod, mmod, nmod;

    /* Parameter adjustments */
    --cz;

    /* Function Body */
    *zh = 0.f;
    cdt = cos(*th + *th);
    sdt = sin(*th + *th);
    lmod = *nlat % 2;
    mmod = *m % 2;
    nmod = *n % 2;
    if (lmod <= 0) {
	goto L20;
    } else {
	goto L10;
    }
L10:
    lc = (*nlat + 1) / 2;
    lq = lc - 1;
    ls = lc - 2;
    if (nmod <= 0) {
	goto L1;
    } else {
	goto L2;
    }
L1:
    if (mmod <= 0) {
	goto L3;
    } else {
	goto L4;
    }

/*     nlat odd n even m even */

L3:
    *zh = (cz[1] + cz[lc] * cos((lq << 1) * *th)) * .5f;
    cth = cdt;
    sth = sdt;
    i__1 = lq;
    for (k = 2; k <= i__1; ++k) {
/*     zh = zh+cz(k)*dcos(2*(k-1)*th) */
	*zh += cz[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L201: */
    }
    return 0;

/*     nlat odd n even m odd */

L4:
    cth = cdt;
    sth = sdt;
    i__1 = ls;
    for (k = 1; k <= i__1; ++k) {
/*     zh = zh+cz(k+1)*dsin(2*k*th) */
	*zh += cz[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L202: */
    }
    return 0;

/*     nlat odd n odd, m even */

L2:
    if (mmod <= 0) {
	goto L5;
    } else {
	goto L6;
    }
L5:
    cth = cos(*th);
    sth = sin(*th);
    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
/*     zh = zh+cz(k)*dcos((2*k-1)*th) */
	*zh += cz[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L203: */
    }
    return 0;

/*     nlat odd n odd m odd */

L6:
    cth = cos(*th);
    sth = sin(*th);
    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
/*     zh = zh+cz(k+1)*dsin((2*k-1)*th) */
	*zh += cz[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L204: */
    }
    return 0;
L20:
    lc = *nlat / 2;
    lq = lc - 1;
    if (nmod <= 0) {
	goto L30;
    } else {
	goto L80;
    }
L30:
    if (mmod <= 0) {
	goto L40;
    } else {
	goto L60;
    }

/*     nlat even n even m even */

L40:
    *zh = cz[1] * .5f;
    cth = cdt;
    sth = sdt;
    i__1 = lc;
    for (k = 2; k <= i__1; ++k) {
/*     zh = zh+cz(k)*dcos(2*(k-1)*th) */
	*zh += cz[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L50: */
    }
    return 0;

/*     nlat even n even m odd */

L60:
    cth = cdt;
    sth = sdt;
    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
/*     zh = zh+cz(k+1)*dsin(2*k*th) */
	*zh += cz[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L70: */
    }
    return 0;

/*     nlat even n odd m even */

L80:
    if (mmod <= 0) {
	goto L90;
    } else {
	goto L110;
    }
L90:
    *zh = cz[lc] * .5f * cos((*nlat - 1) * *th);
    cth = cos(*th);
    sth = sin(*th);
    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
/*     zh = zh+cz(k)*dcos((2*k-1)*th) */
	*zh += cz[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L100: */
    }
    return 0;

/*     nlat even n odd m odd */

L110:
    cth = cos(*th);
    sth = sin(*th);
    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
/*     zh = zh+cz(k+1)*dsin((2*k-1)*th) */
	*zh += cz[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L120: */
    }
    return 0;
} /* dnzft_ */

/* Subroutine */ int alin_(integer *isym, integer *nlat, integer *nlon, 
	integer *m, real *p, integer *i3, real *walin)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, iw3, iw4, lim, labc, imid, mmax;
    extern /* Subroutine */ int alin1_(integer *, integer *, integer *, real *
	    , integer *, integer *, real *, real *, real *, real *, real *);

    /* Parameter adjustments */
    --walin;
    --p;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    lim = *nlat * imid;
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mmax = min(i__1,i__2);
    labc = (mmax - 2) * (*nlat + *nlat - mmax - 1) / 2;
    iw1 = lim + 1;
    iw2 = iw1 + lim;
    iw3 = iw2 + labc;
    iw4 = iw3 + labc;

/*     the length of walin is ((5*l-7)*l+6)/2 */

    alin1_(isym, nlat, m, &p[1], &imid, i3, &walin[1], &walin[iw1], &walin[
	    iw2], &walin[iw3], &walin[iw4]);
    return 0;
} /* alin_ */

/* Subroutine */ int alin1_(integer *isym, integer *nlat, integer *m, real *p,
	 integer *imid, integer *i3, real *pz, real *p1, real *a, real *b, 
	real *c__)
{
    /* System generated locals */
    integer p_dim1, p_dim2, p_offset, pz_dim1, pz_offset, p1_dim1, p1_offset, 
	    i__1, i__2, i__3;

    /* Local variables */
    static integer i__, i1, i2, ns, np1, nstp, ihold, nstrt;

    /* Parameter adjustments */
    p1_dim1 = *imid;
    p1_offset = 1 + p1_dim1;
    p1 -= p1_offset;
    pz_dim1 = *imid;
    pz_offset = 1 + pz_dim1;
    pz -= pz_offset;
    p_dim1 = *imid;
    p_dim2 = *nlat;
    p_offset = 1 + p_dim1 * (1 + p_dim2);
    p -= p_offset;
    --a;
    --b;
    --c__;

    /* Function Body */
    ihold = i1;
    i1 = i2;
    i2 = *i3;
    *i3 = ihold;
    if ((i__1 = *m - 1) < 0) {
	goto L25;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L35;
    }
L25:
    i1 = 1;
    i2 = 2;
    *i3 = 3;
    i__1 = *nlat;
    for (np1 = 1; np1 <= i__1; ++np1) {
	i__2 = *imid;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    p[i__ + (np1 + *i3 * p_dim2) * p_dim1] = pz[i__ + np1 * pz_dim1];
/* L45: */
	}
    }
    return 0;
L30:
    i__2 = *nlat;
    for (np1 = 2; np1 <= i__2; ++np1) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    p[i__ + (np1 + *i3 * p_dim2) * p_dim1] = p1[i__ + np1 * p1_dim1];
/* L50: */
	}
    }
    return 0;
L35:
    ns = (*m - 2) * (*nlat + *nlat - *m - 1) / 2 + 1;
    if (*isym == 1) {
	goto L36;
    }
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__ + (*m + 1 + *i3 * p_dim2) * p_dim1] = a[ns] * p[i__ + (*m - 1 + 
		i1 * p_dim2) * p_dim1] - c__[ns] * p[i__ + (*m + 1 + i1 * 
		p_dim2) * p_dim1];
/* L85: */
    }
L36:
    if (*m == *nlat - 1) {
	return 0;
    }
    if (*isym == 2) {
	goto L71;
    }
    ++ns;
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__ + (*m + 2 + *i3 * p_dim2) * p_dim1] = a[ns] * p[i__ + (*m + i1 *
		 p_dim2) * p_dim1] - c__[ns] * p[i__ + (*m + 2 + i1 * p_dim2) 
		* p_dim1];
/* L70: */
    }
L71:
    nstrt = *m + 3;
    if (*isym == 1) {
	nstrt = *m + 4;
    }
    if (nstrt > *nlat) {
	goto L80;
    }
    nstp = 2;
    if (*isym == 0) {
	nstp = 1;
    }
    i__1 = *nlat;
    i__2 = nstp;
    for (np1 = nstrt; i__2 < 0 ? np1 >= i__1 : np1 <= i__1; np1 += i__2) {
	ns += nstp;
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    p[i__ + (np1 + *i3 * p_dim2) * p_dim1] = a[ns] * p[i__ + (np1 - 2 
		    + i1 * p_dim2) * p_dim1] + b[ns] * p[i__ + (np1 - 2 + *i3 
		    * p_dim2) * p_dim1] - c__[ns] * p[i__ + (np1 + i1 * 
		    p_dim2) * p_dim1];
/* L75: */
	}
    }
L80:
    return 0;
} /* alin1_ */

/* Subroutine */ int alinit_(integer *nlat, integer *nlon, real *walin, 
	doublereal *dwork)
{
    static integer iw1, imid;
    extern /* Subroutine */ int alini1_(integer *, integer *, integer *, real 
	    *, real *, doublereal *);

    /* Parameter adjustments */
    --dwork;
    --walin;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     the length of walin is 3*((l-3)*l+2)/2 + 2*l*imid */
/*     the length of work is nlat+1 */

    alini1_(nlat, nlon, &imid, &walin[1], &walin[iw1], &dwork[1]);
    return 0;
} /* alinit_ */

/* Subroutine */ int alini1_(integer *nlat, integer *nlon, integer *imid, 
	real *p, real *abc, doublereal *cp)
{
    /* System generated locals */
    integer p_dim1, p_dim2, p_offset, i__1, i__2;

    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static integer i__, m, n;
    static doublereal dt, pi, ph, th;
    static integer mp1, np1;
    extern /* Subroutine */ int rabcp_(integer *, integer *, real *), dnlfk_(
	    integer *, integer *, doublereal *), dnlft_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);

    /* Parameter adjustments */
    p_dim1 = *imid;
    p_dim2 = *nlat;
    p_offset = 1 + p_dim1 * (1 + p_dim2);
    p -= p_offset;
    --abc;
    --cp;

    /* Function Body */
    pi = atan(1.) * 4.f;
    dt = pi / (*nlat - 1);
    for (mp1 = 1; mp1 <= 2; ++mp1) {
	m = mp1 - 1;
	i__1 = *nlat;
	for (np1 = mp1; np1 <= i__1; ++np1) {
	    n = np1 - 1;
	    dnlfk_(&m, &n, &cp[1]);
	    i__2 = *imid;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		th = (i__ - 1) * dt;
		dnlft_(&m, &n, &th, &cp[1], &ph);
		p[i__ + (np1 + mp1 * p_dim2) * p_dim1] = ph;
/* L160: */
	    }
	}
    }
    rabcp_(nlat, nlon, &abc[1]);
    return 0;
} /* alini1_ */

/* Subroutine */ int rabcp_(integer *nlat, integer *nlon, real *abc)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, labc, mmax;
    extern /* Subroutine */ int rabcp1_(integer *, integer *, real *, real *, 
	    real *);


/*     subroutine rabcp computes the coefficients in the recurrence */
/*     relation for the associated legendre fuctions. array abc */
/*     must have 3*((mmax-2)*(nlat+nlat-mmax-1))/2 locations. */

    /* Parameter adjustments */
    --abc;

    /* Function Body */
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mmax = min(i__1,i__2);
    labc = (mmax - 2) * (*nlat + *nlat - mmax - 1) / 2;
    iw1 = labc + 1;
    iw2 = iw1 + labc;
    rabcp1_(nlat, nlon, &abc[1], &abc[iw1], &abc[iw2]);
    return 0;
} /* rabcp_ */

/* Subroutine */ int rabcp1_(integer *nlat, integer *nlon, real *a, real *b, 
	real *c__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer m, n;
    static real cn, fm, fn;
    static integer ns;
    static real tm, tn;
    static integer mp1, np1, mp3;
    static real fnmm, fnpm;
    static integer mmax;
    static real temp;


/*     coefficients a, b, and c for computing pbar(m,n,theta) are */
/*     stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1 */

    /* Parameter adjustments */
    --c__;
    --b;
    --a;

    /* Function Body */
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mmax = min(i__1,i__2);
    i__1 = mmax;
    for (mp1 = 3; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	ns = (m - 2) * (*nlat + *nlat - m - 1) / 2 + 1;
	fm = (real) m;
	tm = fm + fm;
	temp = tm * (tm - 1.f);
	a[ns] = sqrt((tm + 1.f) * (tm - 2.f) / temp);
	c__[ns] = sqrt(2.f / temp);
	if (m == *nlat - 1) {
	    goto L215;
	}
	++ns;
	temp = tm * (tm + 1.f);
	a[ns] = sqrt((tm + 3.f) * (tm - 2.f) / temp);
	c__[ns] = sqrt(6.f / temp);
	mp3 = m + 3;
	if (mp3 > *nlat) {
	    goto L215;
	}
	i__2 = *nlat;
	for (np1 = mp3; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    ++ns;
	    fn = (real) n;
	    tn = fn + fn;
	    cn = (tn + 1.f) / (tn - 3.f);
	    fnpm = fn + fm;
	    fnmm = fn - fm;
	    temp = fnpm * (fnpm - 1.f);
	    a[ns] = sqrt(cn * (fnpm - 3.f) * (fnpm - 2.f) / temp);
	    b[ns] = sqrt(cn * fnmm * (fnmm - 1.f) / temp);
	    c__[ns] = sqrt((fnmm + 1.f) * (fnmm + 2.f) / temp);
/* L210: */
	}
L215:
	;
    }
    return 0;
} /* rabcp1_ */

/* Subroutine */ int sea1_(integer *nlat, integer *nlon, integer *imid, real *
	z__, integer *idz, real *zin, real *wzfin, doublereal *dwork)
{
    /* System generated locals */
    integer z_dim1, z_offset, zin_dim1, zin_dim2, zin_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, m, i3, mn, mp1, np1, mmax;
    extern /* Subroutine */ int zfin_(integer *, integer *, integer *, 
	    integer *, real *, integer *, real *), zfinit_(integer *, integer 
	    *, real *, doublereal *);

    /* Parameter adjustments */
    zin_dim1 = *imid;
    zin_dim2 = *nlat;
    zin_offset = 1 + zin_dim1 * (1 + zin_dim2);
    zin -= zin_offset;
    z_dim1 = *idz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --wzfin;
    --dwork;

    /* Function Body */
    zfinit_(nlat, nlon, &wzfin[1], &dwork[1]);
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mmax = min(i__1,i__2);
    i__1 = mmax;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	zfin_(&c__0, nlat, nlon, &m, &zin[zin_offset], &i3, &wzfin[1]);
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    mn = m * (*nlat - 1) - m * (m - 1) / 2 + np1;
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		z__[mn + i__ * z_dim1] = zin[i__ + (np1 + i3 * zin_dim2) * 
			zin_dim1];
/* L33: */
	    }
	}
    }
    return 0;
} /* sea1_ */

/* Subroutine */ int ses1_(integer *nlat, integer *nlon, integer *imid, real *
	p, real *pin, real *walin, doublereal *dwork)
{
    /* System generated locals */
    integer p_dim1, p_offset, pin_dim1, pin_dim2, pin_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, m, i3, mn, mp1, np1;
    extern /* Subroutine */ int alin_(integer *, integer *, integer *, 
	    integer *, real *, integer *, real *);
    static integer mmax;
    extern /* Subroutine */ int alinit_(integer *, integer *, real *, 
	    doublereal *);

    /* Parameter adjustments */
    pin_dim1 = *imid;
    pin_dim2 = *nlat;
    pin_offset = 1 + pin_dim1 * (1 + pin_dim2);
    pin -= pin_offset;
    p_dim1 = *imid;
    p_offset = 1 + p_dim1;
    p -= p_offset;
    --walin;
    --dwork;

    /* Function Body */
    alinit_(nlat, nlon, &walin[1], &dwork[1]);
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mmax = min(i__1,i__2);
    i__1 = mmax;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	alin_(&c__0, nlat, nlon, &m, &pin[pin_offset], &i3, &walin[1]);
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    mn = m * (*nlat - 1) - m * (m - 1) / 2 + np1;
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		p[i__ + mn * p_dim1] = pin[i__ + (np1 + i3 * pin_dim2) * 
			pin_dim1];
/* L10: */
	    }
	}
    }
    return 0;
} /* ses1_ */

/* Subroutine */ int zvinit_(integer *nlat, integer *nlon, real *wzvin, 
	doublereal *dwork)
{
    static integer iw1, imid;
    extern /* Subroutine */ int zvini1_(integer *, integer *, integer *, real 
	    *, real *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --dwork;
    --wzvin;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     the length of wzvin is */
/*         2*nlat*imid +3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 */
/*     the length of dwork is nlat+2 */

    zvini1_(nlat, nlon, &imid, &wzvin[1], &wzvin[iw1], &dwork[1], &dwork[*
	    nlat / 2 + 2]);
    return 0;
} /* zvinit_ */

/* Subroutine */ int zvini1_(integer *nlat, integer *nlon, integer *imid, 
	real *zv, real *abc, doublereal *czv, doublereal *work)
{
    /* System generated locals */
    integer zv_dim1, zv_dim2, zv_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static integer i__, m, n;
    static doublereal dt, pi, th;
    static integer mp1, np1, mdo;
    static doublereal zvh;
    extern /* Subroutine */ int dzvk_(integer *, integer *, integer *, 
	    doublereal *, doublereal *), dzvt_(integer *, integer *, integer *
	    , doublereal *, doublereal *, doublereal *), rabcv_(integer *, 
	    integer *, real *);


/*     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 */
/*     locations where mmax = min0(nlat,(nlon+1)/2) */
/*     czv and work must each have nlat/2+1  locations */

    /* Parameter adjustments */
    zv_dim1 = *imid;
    zv_dim2 = *nlat;
    zv_offset = 1 + zv_dim1 * (1 + zv_dim2);
    zv -= zv_offset;
    --abc;
    --czv;
    --work;

    /* Function Body */
    pi = atan(1.) * 4.f;
    dt = pi / (*nlat - 1);
/* Computing MIN */
    i__1 = min(2,*nlat), i__2 = (*nlon + 1) / 2;
    mdo = min(i__1,i__2);
    i__1 = mdo;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    dzvk_(nlat, &m, &n, &czv[1], &work[1]);
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		th = (i__ - 1) * dt;
		dzvt_(nlat, &m, &n, &th, &czv[1], &zvh);
		zv[i__ + (np1 + mp1 * zv_dim2) * zv_dim1] = zvh;
/* L165: */
	    }
	    zv[(np1 + mp1 * zv_dim2) * zv_dim1 + 1] *= .5f;
/* L160: */
	}
    }
    rabcv_(nlat, nlon, &abc[1]);
    return 0;
} /* zvini1_ */

/* Subroutine */ int zwinit_(integer *nlat, integer *nlon, real *wzwin, 
	doublereal *dwork)
{
    static integer iw1, imid;
    extern /* Subroutine */ int zwini1_(integer *, integer *, integer *, real 
	    *, real *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --dwork;
    --wzwin;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     the length of wzvin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2 */
/*     the length of dwork is nlat+2 */

    zwini1_(nlat, nlon, &imid, &wzwin[1], &wzwin[iw1], &dwork[1], &dwork[*
	    nlat / 2 + 2]);
    return 0;
} /* zwinit_ */

/* Subroutine */ int zwini1_(integer *nlat, integer *nlon, integer *imid, 
	real *zw, real *abc, doublereal *czw, doublereal *work)
{
    /* System generated locals */
    integer zw_dim1, zw_dim2, zw_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static integer i__, m, n;
    static doublereal dt, pi, th;
    static integer mp1, np1, mdo;
    static doublereal zwh;
    extern /* Subroutine */ int dzwk_(integer *, integer *, integer *, 
	    doublereal *, doublereal *), dzwt_(integer *, integer *, integer *
	    , doublereal *, doublereal *, doublereal *), rabcw_(integer *, 
	    integer *, real *);


/*     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 */
/*     locations where mmax = min0(nlat,(nlon+1)/2) */
/*     czw and work must each have nlat+1 locations */

    /* Parameter adjustments */
    zw_dim1 = *imid;
    zw_dim2 = *nlat;
    zw_offset = 1 + zw_dim1 * (1 + zw_dim2);
    zw -= zw_offset;
    --abc;
    --czw;
    --work;

    /* Function Body */
    pi = atan(1.) * 4.f;
    dt = pi / (*nlat - 1);
/* Computing MIN */
    i__1 = min(3,*nlat), i__2 = (*nlon + 1) / 2;
    mdo = min(i__1,i__2);
    if (mdo < 2) {
	return 0;
    }
    i__1 = mdo;
    for (mp1 = 2; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    dzwk_(nlat, &m, &n, &czw[1], &work[1]);
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		th = (i__ - 1) * dt;
		dzwt_(nlat, &m, &n, &th, &czw[1], &zwh);
		zw[i__ + (np1 + m * zw_dim2) * zw_dim1] = zwh;
/* L165: */
	    }
	    zw[(np1 + m * zw_dim2) * zw_dim1 + 1] *= .5f;
/* L160: */
	}
    }
    rabcw_(nlat, nlon, &abc[1]);
    return 0;
} /* zwini1_ */

/* Subroutine */ int zvin_(integer *ityp, integer *nlat, integer *nlon, 
	integer *m, real *zv, integer *i3, real *wzvin)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, iw3, iw4, lim, labc, imid, mmax;
    extern /* Subroutine */ int zvin1_(integer *, integer *, integer *, real *
	    , integer *, integer *, real *, real *, real *, real *, real *);

    /* Parameter adjustments */
    --wzvin;
    --zv;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    lim = *nlat * imid;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) / 2;
    iw1 = lim + 1;
    iw2 = iw1 + lim;
    iw3 = iw2 + labc;
    iw4 = iw3 + labc;

/*     the length of wzvin is 2*lim+3*labc */

    zvin1_(ityp, nlat, m, &zv[1], &imid, i3, &wzvin[1], &wzvin[iw1], &wzvin[
	    iw2], &wzvin[iw3], &wzvin[iw4]);
    return 0;
} /* zvin_ */

/* Subroutine */ int zvin1_(integer *ityp, integer *nlat, integer *m, real *
	zv, integer *imid, integer *i3, real *zvz, real *zv1, real *a, real *
	b, real *c__)
{
    /* System generated locals */
    integer zv_dim1, zv_dim2, zv_offset, zvz_dim1, zvz_offset, zv1_dim1, 
	    zv1_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, i1, i2, ns, np1, nstp, ihold, nstrt;

    /* Parameter adjustments */
    zv1_dim1 = *imid;
    zv1_offset = 1 + zv1_dim1;
    zv1 -= zv1_offset;
    zvz_dim1 = *imid;
    zvz_offset = 1 + zvz_dim1;
    zvz -= zvz_offset;
    zv_dim1 = *imid;
    zv_dim2 = *nlat;
    zv_offset = 1 + zv_dim1 * (1 + zv_dim2);
    zv -= zv_offset;
    --a;
    --b;
    --c__;

    /* Function Body */
    ihold = i1;
    i1 = i2;
    i2 = *i3;
    *i3 = ihold;
    if ((i__1 = *m - 1) < 0) {
	goto L25;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L35;
    }
L25:
    i1 = 1;
    i2 = 2;
    *i3 = 3;
    i__1 = *nlat;
    for (np1 = 1; np1 <= i__1; ++np1) {
	i__2 = *imid;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    zv[i__ + (np1 + *i3 * zv_dim2) * zv_dim1] = zvz[i__ + np1 * 
		    zvz_dim1];
/* L45: */
	}
    }
    return 0;
L30:
    i__2 = *nlat;
    for (np1 = 2; np1 <= i__2; ++np1) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    zv[i__ + (np1 + *i3 * zv_dim2) * zv_dim1] = zv1[i__ + np1 * 
		    zv1_dim1];
/* L50: */
	}
    }
    return 0;
L35:
    ns = (*m - 2) * (*nlat + *nlat - *m - 1) / 2 + 1;
    if (*ityp == 1) {
	goto L36;
    }
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zv[i__ + (*m + 1 + *i3 * zv_dim2) * zv_dim1] = a[ns] * zv[i__ + (*m - 
		1 + i1 * zv_dim2) * zv_dim1] - c__[ns] * zv[i__ + (*m + 1 + 
		i1 * zv_dim2) * zv_dim1];
/* L85: */
    }
L36:
    if (*m == *nlat - 1) {
	return 0;
    }
    if (*ityp == 2) {
	goto L71;
    }
    ++ns;
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zv[i__ + (*m + 2 + *i3 * zv_dim2) * zv_dim1] = a[ns] * zv[i__ + (*m + 
		i1 * zv_dim2) * zv_dim1] - c__[ns] * zv[i__ + (*m + 2 + i1 * 
		zv_dim2) * zv_dim1];
/* L70: */
    }
L71:
    nstrt = *m + 3;
    if (*ityp == 1) {
	nstrt = *m + 4;
    }
    if (nstrt > *nlat) {
	goto L80;
    }
    nstp = 2;
    if (*ityp == 0) {
	nstp = 1;
    }
    i__1 = *nlat;
    i__2 = nstp;
    for (np1 = nstrt; i__2 < 0 ? np1 >= i__1 : np1 <= i__1; np1 += i__2) {
	ns += nstp;
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    zv[i__ + (np1 + *i3 * zv_dim2) * zv_dim1] = a[ns] * zv[i__ + (np1 
		    - 2 + i1 * zv_dim2) * zv_dim1] + b[ns] * zv[i__ + (np1 - 
		    2 + *i3 * zv_dim2) * zv_dim1] - c__[ns] * zv[i__ + (np1 + 
		    i1 * zv_dim2) * zv_dim1];
/* L75: */
	}
    }
L80:
    return 0;
} /* zvin1_ */

/* Subroutine */ int zwin_(integer *ityp, integer *nlat, integer *nlon, 
	integer *m, real *zw, integer *i3, real *wzwin)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, iw3, iw4, lim, labc, imid, mmax;
    extern /* Subroutine */ int zwin1_(integer *, integer *, integer *, real *
	    , integer *, integer *, real *, real *, real *, real *, real *);

    /* Parameter adjustments */
    --wzwin;
    --zw;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    lim = *nlat * imid;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) / 2;
    iw1 = lim + 1;
    iw2 = iw1 + lim;
    iw3 = iw2 + labc;
    iw4 = iw3 + labc;

/*     the length of wzwin is 2*lim+3*labc */

    zwin1_(ityp, nlat, m, &zw[1], &imid, i3, &wzwin[1], &wzwin[iw1], &wzwin[
	    iw2], &wzwin[iw3], &wzwin[iw4]);
    return 0;
} /* zwin_ */

/* Subroutine */ int zwin1_(integer *ityp, integer *nlat, integer *m, real *
	zw, integer *imid, integer *i3, real *zw1, real *zw2, real *a, real *
	b, real *c__)
{
    /* System generated locals */
    integer zw_dim1, zw_dim2, zw_offset, zw1_dim1, zw1_offset, zw2_dim1, 
	    zw2_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, i1, i2, ns, np1, nstp, ihold, nstrt;

    /* Parameter adjustments */
    zw2_dim1 = *imid;
    zw2_offset = 1 + zw2_dim1;
    zw2 -= zw2_offset;
    zw1_dim1 = *imid;
    zw1_offset = 1 + zw1_dim1;
    zw1 -= zw1_offset;
    zw_dim1 = *imid;
    zw_dim2 = *nlat;
    zw_offset = 1 + zw_dim1 * (1 + zw_dim2);
    zw -= zw_offset;
    --a;
    --b;
    --c__;

    /* Function Body */
    ihold = i1;
    i1 = i2;
    i2 = *i3;
    *i3 = ihold;
    if ((i__1 = *m - 2) < 0) {
	goto L25;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L35;
    }
L25:
    i1 = 1;
    i2 = 2;
    *i3 = 3;
    i__1 = *nlat;
    for (np1 = 2; np1 <= i__1; ++np1) {
	i__2 = *imid;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    zw[i__ + (np1 + *i3 * zw_dim2) * zw_dim1] = zw1[i__ + np1 * 
		    zw1_dim1];
/* L45: */
	}
    }
    return 0;
L30:
    i__2 = *nlat;
    for (np1 = 3; np1 <= i__2; ++np1) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    zw[i__ + (np1 + *i3 * zw_dim2) * zw_dim1] = zw2[i__ + np1 * 
		    zw2_dim1];
/* L50: */
	}
    }
    return 0;
L35:
    ns = (*m - 2) * (*nlat + *nlat - *m - 1) / 2 + 1;
    if (*ityp == 1) {
	goto L36;
    }
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zw[i__ + (*m + 1 + *i3 * zw_dim2) * zw_dim1] = a[ns] * zw[i__ + (*m - 
		1 + i1 * zw_dim2) * zw_dim1] - c__[ns] * zw[i__ + (*m + 1 + 
		i1 * zw_dim2) * zw_dim1];
/* L85: */
    }
L36:
    if (*m == *nlat - 1) {
	return 0;
    }
    if (*ityp == 2) {
	goto L71;
    }
    ++ns;
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	zw[i__ + (*m + 2 + *i3 * zw_dim2) * zw_dim1] = a[ns] * zw[i__ + (*m + 
		i1 * zw_dim2) * zw_dim1] - c__[ns] * zw[i__ + (*m + 2 + i1 * 
		zw_dim2) * zw_dim1];
/* L70: */
    }
L71:
    nstrt = *m + 3;
    if (*ityp == 1) {
	nstrt = *m + 4;
    }
    if (nstrt > *nlat) {
	goto L80;
    }
    nstp = 2;
    if (*ityp == 0) {
	nstp = 1;
    }
    i__1 = *nlat;
    i__2 = nstp;
    for (np1 = nstrt; i__2 < 0 ? np1 >= i__1 : np1 <= i__1; np1 += i__2) {
	ns += nstp;
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    zw[i__ + (np1 + *i3 * zw_dim2) * zw_dim1] = a[ns] * zw[i__ + (np1 
		    - 2 + i1 * zw_dim2) * zw_dim1] + b[ns] * zw[i__ + (np1 - 
		    2 + *i3 * zw_dim2) * zw_dim1] - c__[ns] * zw[i__ + (np1 + 
		    i1 * zw_dim2) * zw_dim1];
/* L75: */
	}
    }
L80:
    return 0;
} /* zwin1_ */

/* Subroutine */ int vbinit_(integer *nlat, integer *nlon, real *wvbin, 
	doublereal *dwork)
{
    static integer iw1, imid;
    extern /* Subroutine */ int vbini1_(integer *, integer *, integer *, real 
	    *, real *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --dwork;
    --wvbin;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2 */
/*     the length of dwork is nlat+2 */

    vbini1_(nlat, nlon, &imid, &wvbin[1], &wvbin[iw1], &dwork[1], &dwork[*
	    nlat / 2 + 2]);
    return 0;
} /* vbinit_ */

/* Subroutine */ int vbini1_(integer *nlat, integer *nlon, integer *imid, 
	real *vb, real *abc, doublereal *cvb, doublereal *work)
{
    /* System generated locals */
    integer vb_dim1, vb_dim2, vb_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static integer i__, m, n;
    static doublereal pi, dt, th;
    static integer mp1, np1;
    static doublereal vbh;
    static integer mdo;
    extern /* Subroutine */ int dvbk_(integer *, integer *, doublereal *, 
	    doublereal *), dvbt_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), rabcv_(integer *, integer *, real *);


/*     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 */
/*     locations where mmax = min0(nlat,(nlon+1)/2) */
/*     cvb and work must each have nlat+1 locations */

    /* Parameter adjustments */
    vb_dim1 = *imid;
    vb_dim2 = *nlat;
    vb_offset = 1 + vb_dim1 * (1 + vb_dim2);
    vb -= vb_offset;
    --abc;
    --cvb;
    --work;

    /* Function Body */
    pi = atan(1.) * 4.f;
    dt = pi / (*nlat - 1);
/* Computing MIN */
    i__1 = min(2,*nlat), i__2 = (*nlon + 1) / 2;
    mdo = min(i__1,i__2);
    i__1 = mdo;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    dvbk_(&m, &n, &cvb[1], &work[1]);
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		th = (i__ - 1) * dt;
		dvbt_(&m, &n, &th, &cvb[1], &vbh);
		vb[i__ + (np1 + mp1 * vb_dim2) * vb_dim1] = vbh;
/* L165: */
	    }
/* L160: */
	}
    }
    rabcv_(nlat, nlon, &abc[1]);
    return 0;
} /* vbini1_ */

/* Subroutine */ int wbinit_(integer *nlat, integer *nlon, real *wwbin, 
	doublereal *dwork)
{
    static integer iw1, imid;
    extern /* Subroutine */ int wbini1_(integer *, integer *, integer *, real 
	    *, real *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --dwork;
    --wwbin;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2 */
/*     the length of dwork is nlat+2 */

    wbini1_(nlat, nlon, &imid, &wwbin[1], &wwbin[iw1], &dwork[1], &dwork[*
	    nlat / 2 + 2]);
    return 0;
} /* wbinit_ */

/* Subroutine */ int wbini1_(integer *nlat, integer *nlon, integer *imid, 
	real *wb, real *abc, doublereal *cwb, doublereal *work)
{
    /* System generated locals */
    integer wb_dim1, wb_dim2, wb_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static integer i__, m, n;
    static doublereal dt, pi, th;
    static integer mp1, np1, mdo;
    static doublereal wbh;
    extern /* Subroutine */ int dwbk_(integer *, integer *, doublereal *, 
	    doublereal *), dwbt_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), rabcw_(integer *, integer *, real *);


/*     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 */
/*     locations where mmax = min0(nlat,(nlon+1)/2) */
/*     cwb and work must each have nlat/2+1 locations */

    /* Parameter adjustments */
    wb_dim1 = *imid;
    wb_dim2 = *nlat;
    wb_offset = 1 + wb_dim1 * (1 + wb_dim2);
    wb -= wb_offset;
    --abc;
    --cwb;
    --work;

    /* Function Body */
    pi = atan(1.) * 4.f;
    dt = pi / (*nlat - 1);
/* Computing MIN */
    i__1 = min(3,*nlat), i__2 = (*nlon + 1) / 2;
    mdo = min(i__1,i__2);
    if (mdo < 2) {
	return 0;
    }
    i__1 = mdo;
    for (mp1 = 2; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    dwbk_(&m, &n, &cwb[1], &work[1]);
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		th = (i__ - 1) * dt;
		dwbt_(&m, &n, &th, &cwb[1], &wbh);
		wb[i__ + (np1 + m * wb_dim2) * wb_dim1] = wbh;
/* L165: */
	    }
/* L160: */
	}
    }
    rabcw_(nlat, nlon, &abc[1]);
    return 0;
} /* wbini1_ */

/* Subroutine */ int vbin_(integer *ityp, integer *nlat, integer *nlon, 
	integer *m, real *vb, integer *i3, real *wvbin)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, iw3, iw4, lim, labc, imid, mmax;
    extern /* Subroutine */ int vbin1_(integer *, integer *, integer *, real *
	    , integer *, integer *, real *, real *, real *, real *, real *);

    /* Parameter adjustments */
    --wvbin;
    --vb;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    lim = *nlat * imid;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) / 2;
    iw1 = lim + 1;
    iw2 = iw1 + lim;
    iw3 = iw2 + labc;
    iw4 = iw3 + labc;

/*     the length of wvbin is 2*lim+3*labc */

    vbin1_(ityp, nlat, m, &vb[1], &imid, i3, &wvbin[1], &wvbin[iw1], &wvbin[
	    iw2], &wvbin[iw3], &wvbin[iw4]);
    return 0;
} /* vbin_ */

/* Subroutine */ int vbin1_(integer *ityp, integer *nlat, integer *m, real *
	vb, integer *imid, integer *i3, real *vbz, real *vb1, real *a, real *
	b, real *c__)
{
    /* System generated locals */
    integer vb_dim1, vb_dim2, vb_offset, vbz_dim1, vbz_offset, vb1_dim1, 
	    vb1_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, i1, i2, ns, np1, nstp, ihold, nstrt;

    /* Parameter adjustments */
    vb1_dim1 = *imid;
    vb1_offset = 1 + vb1_dim1;
    vb1 -= vb1_offset;
    vbz_dim1 = *imid;
    vbz_offset = 1 + vbz_dim1;
    vbz -= vbz_offset;
    vb_dim1 = *imid;
    vb_dim2 = *nlat;
    vb_offset = 1 + vb_dim1 * (1 + vb_dim2);
    vb -= vb_offset;
    --a;
    --b;
    --c__;

    /* Function Body */
    ihold = i1;
    i1 = i2;
    i2 = *i3;
    *i3 = ihold;
    if ((i__1 = *m - 1) < 0) {
	goto L25;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L35;
    }
L25:
    i1 = 1;
    i2 = 2;
    *i3 = 3;
    i__1 = *nlat;
    for (np1 = 1; np1 <= i__1; ++np1) {
	i__2 = *imid;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    vb[i__ + (np1 + *i3 * vb_dim2) * vb_dim1] = vbz[i__ + np1 * 
		    vbz_dim1];
/* L45: */
	}
    }
    return 0;
L30:
    i__2 = *nlat;
    for (np1 = 2; np1 <= i__2; ++np1) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    vb[i__ + (np1 + *i3 * vb_dim2) * vb_dim1] = vb1[i__ + np1 * 
		    vb1_dim1];
/* L50: */
	}
    }
    return 0;
L35:
    ns = (*m - 2) * (*nlat + *nlat - *m - 1) / 2 + 1;
    if (*ityp == 1) {
	goto L36;
    }
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vb[i__ + (*m + 1 + *i3 * vb_dim2) * vb_dim1] = a[ns] * vb[i__ + (*m - 
		1 + i1 * vb_dim2) * vb_dim1] - c__[ns] * vb[i__ + (*m + 1 + 
		i1 * vb_dim2) * vb_dim1];
/* L85: */
    }
L36:
    if (*m == *nlat - 1) {
	return 0;
    }
    if (*ityp == 2) {
	goto L71;
    }
    ++ns;
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	vb[i__ + (*m + 2 + *i3 * vb_dim2) * vb_dim1] = a[ns] * vb[i__ + (*m + 
		i1 * vb_dim2) * vb_dim1] - c__[ns] * vb[i__ + (*m + 2 + i1 * 
		vb_dim2) * vb_dim1];
/* L70: */
    }
L71:
    nstrt = *m + 3;
    if (*ityp == 1) {
	nstrt = *m + 4;
    }
    if (nstrt > *nlat) {
	goto L80;
    }
    nstp = 2;
    if (*ityp == 0) {
	nstp = 1;
    }
    i__1 = *nlat;
    i__2 = nstp;
    for (np1 = nstrt; i__2 < 0 ? np1 >= i__1 : np1 <= i__1; np1 += i__2) {
	ns += nstp;
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    vb[i__ + (np1 + *i3 * vb_dim2) * vb_dim1] = a[ns] * vb[i__ + (np1 
		    - 2 + i1 * vb_dim2) * vb_dim1] + b[ns] * vb[i__ + (np1 - 
		    2 + *i3 * vb_dim2) * vb_dim1] - c__[ns] * vb[i__ + (np1 + 
		    i1 * vb_dim2) * vb_dim1];
/* L75: */
	}
    }
L80:
    return 0;
} /* vbin1_ */

/* Subroutine */ int wbin_(integer *ityp, integer *nlat, integer *nlon, 
	integer *m, real *wb, integer *i3, real *wwbin)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, iw3, iw4, lim, labc, imid, mmax;
    extern /* Subroutine */ int wbin1_(integer *, integer *, integer *, real *
	    , integer *, integer *, real *, real *, real *, real *, real *);

    /* Parameter adjustments */
    --wwbin;
    --wb;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    lim = *nlat * imid;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) / 2;
    iw1 = lim + 1;
    iw2 = iw1 + lim;
    iw3 = iw2 + labc;
    iw4 = iw3 + labc;

/*     the length of wwbin is 2*lim+3*labc */

    wbin1_(ityp, nlat, m, &wb[1], &imid, i3, &wwbin[1], &wwbin[iw1], &wwbin[
	    iw2], &wwbin[iw3], &wwbin[iw4]);
    return 0;
} /* wbin_ */

/* Subroutine */ int wbin1_(integer *ityp, integer *nlat, integer *m, real *
	wb, integer *imid, integer *i3, real *wb1, real *wb2, real *a, real *
	b, real *c__)
{
    /* System generated locals */
    integer wb_dim1, wb_dim2, wb_offset, wb1_dim1, wb1_offset, wb2_dim1, 
	    wb2_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, i1, i2, ns, np1, nstp, ihold, nstrt;

    /* Parameter adjustments */
    wb2_dim1 = *imid;
    wb2_offset = 1 + wb2_dim1;
    wb2 -= wb2_offset;
    wb1_dim1 = *imid;
    wb1_offset = 1 + wb1_dim1;
    wb1 -= wb1_offset;
    wb_dim1 = *imid;
    wb_dim2 = *nlat;
    wb_offset = 1 + wb_dim1 * (1 + wb_dim2);
    wb -= wb_offset;
    --a;
    --b;
    --c__;

    /* Function Body */
    ihold = i1;
    i1 = i2;
    i2 = *i3;
    *i3 = ihold;
    if ((i__1 = *m - 2) < 0) {
	goto L25;
    } else if (i__1 == 0) {
	goto L30;
    } else {
	goto L35;
    }
L25:
    i1 = 1;
    i2 = 2;
    *i3 = 3;
    i__1 = *nlat;
    for (np1 = 2; np1 <= i__1; ++np1) {
	i__2 = *imid;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    wb[i__ + (np1 + *i3 * wb_dim2) * wb_dim1] = wb1[i__ + np1 * 
		    wb1_dim1];
/* L45: */
	}
    }
    return 0;
L30:
    i__2 = *nlat;
    for (np1 = 3; np1 <= i__2; ++np1) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    wb[i__ + (np1 + *i3 * wb_dim2) * wb_dim1] = wb2[i__ + np1 * 
		    wb2_dim1];
/* L50: */
	}
    }
    return 0;
L35:
    ns = (*m - 2) * (*nlat + *nlat - *m - 1) / 2 + 1;
    if (*ityp == 1) {
	goto L36;
    }
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wb[i__ + (*m + 1 + *i3 * wb_dim2) * wb_dim1] = a[ns] * wb[i__ + (*m - 
		1 + i1 * wb_dim2) * wb_dim1] - c__[ns] * wb[i__ + (*m + 1 + 
		i1 * wb_dim2) * wb_dim1];
/* L85: */
    }
L36:
    if (*m == *nlat - 1) {
	return 0;
    }
    if (*ityp == 2) {
	goto L71;
    }
    ++ns;
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wb[i__ + (*m + 2 + *i3 * wb_dim2) * wb_dim1] = a[ns] * wb[i__ + (*m + 
		i1 * wb_dim2) * wb_dim1] - c__[ns] * wb[i__ + (*m + 2 + i1 * 
		wb_dim2) * wb_dim1];
/* L70: */
    }
L71:
    nstrt = *m + 3;
    if (*ityp == 1) {
	nstrt = *m + 4;
    }
    if (nstrt > *nlat) {
	goto L80;
    }
    nstp = 2;
    if (*ityp == 0) {
	nstp = 1;
    }
    i__1 = *nlat;
    i__2 = nstp;
    for (np1 = nstrt; i__2 < 0 ? np1 >= i__1 : np1 <= i__1; np1 += i__2) {
	ns += nstp;
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    wb[i__ + (np1 + *i3 * wb_dim2) * wb_dim1] = a[ns] * wb[i__ + (np1 
		    - 2 + i1 * wb_dim2) * wb_dim1] + b[ns] * wb[i__ + (np1 - 
		    2 + *i3 * wb_dim2) * wb_dim1] - c__[ns] * wb[i__ + (np1 + 
		    i1 * wb_dim2) * wb_dim1];
/* L75: */
	}
    }
L80:
    return 0;
} /* wbin1_ */

/* Subroutine */ int dzvk_(integer *nlat, integer *m, integer *n, doublereal *
	czv, doublereal *work)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k;
    static doublereal t1, t2;
    static integer id, lc;
    static doublereal sc1;
    static integer kdo;
    static doublereal sum;
    extern /* Subroutine */ int dvbk_(integer *, integer *, doublereal *, 
	    doublereal *);
    static integer mmod, nmod;


/*     subroutine dzvk computes the coefficients in the trigonometric */
/*     expansion of the quadrature function zvbar(n,m,theta) */

/*     input parameters */

/*     nlat      the number of colatitudes including the poles. */

/*     n      the degree (subscript) of wbarv(n,m,theta) */

/*     m      the order (superscript) of wbarv(n,m,theta) */

/*     work   a work array with at least nlat/2+1 locations */

/*     output parameter */

/*     czv     the fourier coefficients of zvbar(n,m,theta). */

    /* Parameter adjustments */
    --work;
    --czv;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    lc = (*nlat + 1) / 2;
    sc1 = 2. / (real) (*nlat - 1);
    dvbk_(m, n, &work[1], &czv[1]);
    nmod = *n % 2;
    mmod = *m % 2;
    if (nmod != 0) {
	goto L1;
    }
    if (mmod != 0) {
	goto L2;
    }

/*     n even, m even */

    kdo = *n / 2;
    i__1 = lc;
    for (id = 1; id <= i__1; ++id) {
	i__ = id + id - 2;
	sum = 0.f;
	i__2 = kdo;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    i__3 = k + k + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - i__;
	    t2 = 1. - i__3 * i__3;
	    sum += work[k] * (t1 - t2) / (t1 * t2);
/* L10: */
	}
	czv[id] = sc1 * sum;
/* L9: */
    }
    return 0;

/*     n even, m odd */

L2:
    kdo = *n / 2;
    i__1 = lc;
    for (id = 1; id <= i__1; ++id) {
	i__ = id + id - 2;
	sum = 0.f;
	i__2 = kdo;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    i__3 = k + k + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - i__;
	    t2 = 1. - i__3 * i__3;
	    sum += work[k] * (t1 + t2) / (t1 * t2);
/* L6: */
	}
	czv[id] = sc1 * sum;
/* L5: */
    }
    return 0;
L1:
    if (mmod != 0) {
	goto L3;
    }

/*     n odd, m even */

    kdo = (*n + 1) / 2;
    i__1 = lc;
    for (id = 1; id <= i__1; ++id) {
	i__ = id + id - 3;
	sum = 0.f;
	i__2 = kdo;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    i__3 = k + k - 1 + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - 1 - i__;
	    t2 = 1. - i__3 * i__3;
	    sum += work[k] * (t1 - t2) / (t1 * t2);
/* L20: */
	}
	czv[id] = sc1 * sum;
/* L19: */
    }
    return 0;

/*     n odd, m odd */

L3:
    kdo = (*n + 1) / 2;
    i__1 = lc;
    for (id = 1; id <= i__1; ++id) {
	i__ = id + id - 1;
	sum = 0.f;
	i__2 = kdo;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    i__3 = k + k - 1 + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - 1 - i__;
	    t2 = 1. - i__3 * i__3;
	    sum += work[k] * (t1 + t2) / (t1 * t2);
/* L16: */
	}
	czv[id] = sc1 * sum;
/* L15: */
    }
    return 0;
} /* dzvk_ */

/* Subroutine */ int dzvt_(integer *nlat, integer *m, integer *n, doublereal *
	th, doublereal *czv, doublereal *zvh)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k, lc, lq, ls;
    static doublereal chh, cdt, cth, sdt, sth;
    static integer lmod, mmod, nmod;


/*     subroutine dzvt tabulates the function zvbar(n,m,theta) */
/*     at theta = th in double precision */

/*     input parameters */

/*     nlat      the number of colatitudes including the poles. */

/*     n      the degree (subscript) of zvbar(n,m,theta) */

/*     m      the order (superscript) of zvbar(n,m,theta) */

/*     czv     the fourier coefficients of zvbar(n,m,theta) */
/*             as computed by subroutine zwk. */

/*     output parameter */

/*     zvh     zvbar(m,n,theta) evaluated at theta = th */

    /* Parameter adjustments */
    --czv;

    /* Function Body */
    *zvh = 0.f;
    if (*n <= 0) {
	return 0;
    }
    lc = (*nlat + 1) / 2;
    lq = lc - 1;
    ls = lc - 2;
    cth = cos(*th);
    sth = sin(*th);
    cdt = cth * cth - sth * sth;
    sdt = sth * 2.f * cth;
    lmod = *nlat % 2;
    mmod = *m % 2;
    nmod = *n % 2;
    if (lmod == 0) {
	goto L50;
    }
    if (nmod != 0) {
	goto L1;
    }
    cth = cdt;
    sth = sdt;
    if (mmod != 0) {
	goto L2;
    }

/*     nlat odd  n even  m even */

    i__1 = ls;
    for (k = 1; k <= i__1; ++k) {
	*zvh += czv[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L10: */
    }
    return 0;

/*     nlat odd  n even  m odd */

L2:
    *zvh = czv[1] * .5f;
    i__1 = lq;
    for (k = 2; k <= i__1; ++k) {
	*zvh += czv[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L20: */
    }
    *zvh += czv[lc] * .5f * cos((*nlat - 1) * *th);
    return 0;
L1:
    if (mmod != 0) {
	goto L3;
    }

/*     nlat odd  n odd  m even */

    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
	*zvh += czv[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L30: */
    }
    return 0;

/*     nlat odd  n odd  m odd */

L3:
    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
	*zvh += czv[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L40: */
    }
    return 0;
L50:
    if (nmod != 0) {
	goto L51;
    }
    cth = cdt;
    sth = sdt;
    if (mmod != 0) {
	goto L52;
    }

/*     nlat even  n even  m even */

    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
	*zvh += czv[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L55: */
    }
    return 0;

/*     nlat even  n even  m odd */

L52:
    *zvh = czv[1] * .5f;
    i__1 = lc;
    for (k = 2; k <= i__1; ++k) {
	*zvh += czv[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L57: */
    }
    return 0;
L51:
    if (mmod != 0) {
	goto L53;
    }

/*     nlat even  n odd  m even */

    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
	*zvh += czv[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L58: */
    }
    return 0;

/*     nlat even  n odd  m odd */

L53:
    *zvh = czv[lc] * .5f * cos((*nlat - 1) * *th);
    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
	*zvh += czv[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L60: */
    }
    return 0;
} /* dzvt_ */

/* Subroutine */ int dzwk_(integer *nlat, integer *m, integer *n, doublereal *
	czw, doublereal *work)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k;
    static doublereal t1, t2;
    static integer id, lc;
    static doublereal sc1;
    static integer kp1, kdo;
    static doublereal sum;
    extern /* Subroutine */ int dwbk_(integer *, integer *, doublereal *, 
	    doublereal *);
    static integer mmod, nmod;


/*     subroutine dzwk computes the coefficients in the trigonometric */
/*     expansion of the quadrature function zwbar(n,m,theta) */

/*     input parameters */

/*     nlat      the number of colatitudes including the poles. */

/*     n      the degree (subscript) of zwbar(n,m,theta) */

/*     m      the order (superscript) of zwbar(n,m,theta) */

/*     work   a work array with at least nlat/2+1 locations */

/*     output parameter */

/*     czw     the fourier coefficients of zwbar(n,m,theta). */

    /* Parameter adjustments */
    --work;
    --czw;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    lc = (*nlat + 1) / 2;
    sc1 = 2. / (real) (*nlat - 1);
    dwbk_(m, n, &work[1], &czw[1]);
    nmod = *n % 2;
    mmod = *m % 2;
    if (nmod != 0) {
	goto L1;
    }
    if (mmod != 0) {
	goto L2;
    }

/*     n even, m even */

    kdo = *n / 2;
    i__1 = lc;
    for (id = 1; id <= i__1; ++id) {
	i__ = id + id - 3;
	sum = 0.f;
	i__2 = kdo;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    i__3 = k + k - 1 + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - 1 - i__;
	    t2 = 1. - i__3 * i__3;
	    sum += work[k] * (t1 - t2) / (t1 * t2);
/* L20: */
	}
	czw[id] = sc1 * sum;
/* L19: */
    }
    return 0;

/*     n even, m odd */

L2:
    kdo = *n / 2;
    i__1 = lc;
    for (id = 1; id <= i__1; ++id) {
	i__ = id + id - 1;
	sum = 0.f;
	i__2 = kdo;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    i__3 = k + k - 1 + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - 1 - i__;
	    t2 = 1. - i__3 * i__3;
	    sum += work[k] * (t1 + t2) / (t1 * t2);
/* L16: */
	}
	czw[id] = sc1 * sum;
/* L15: */
    }
    return 0;
L1:
    if (mmod != 0) {
	goto L3;
    }

/*     n odd, m even */

    kdo = (*n - 1) / 2;
    i__1 = lc;
    for (id = 1; id <= i__1; ++id) {
	i__ = id + id - 2;
	sum = 0.f;
	i__2 = kdo;
	for (k = 1; k <= i__2; ++k) {
/* Computing 2nd power */
	    i__3 = k + k + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - i__;
	    t2 = 1. - i__3 * i__3;
	    sum += work[k] * (t1 - t2) / (t1 * t2);
/* L10: */
	}
	czw[id] = sc1 * sum;
/* L9: */
    }
    return 0;

/*     n odd, m odd */

L3:
    kdo = (*n + 1) / 2;
    i__1 = lc;
    for (id = 1; id <= i__1; ++id) {
	i__ = id + id - 2;
	sum = work[1] / (1. - i__ * i__);
	if (kdo < 2) {
	    goto L29;
	}
	i__2 = kdo;
	for (kp1 = 2; kp1 <= i__2; ++kp1) {
	    k = kp1 - 1;
/* Computing 2nd power */
	    i__3 = k + k + i__;
	    t1 = 1. - i__3 * i__3;
/* Computing 2nd power */
	    i__3 = k + k - i__;
	    t2 = 1. - i__3 * i__3;
	    sum += work[kp1] * (t1 + t2) / (t1 * t2);
/* L6: */
	}
L29:
	czw[id] = sc1 * sum;
/* L5: */
    }
    return 0;
} /* dzwk_ */

/* Subroutine */ int dzwt_(integer *nlat, integer *m, integer *n, doublereal *
	th, doublereal *czw, doublereal *zwh)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k, lc, lq, ls;
    static doublereal chh, cdt, cth, sdt, sth;
    static integer lmod, mmod, nmod;


/*     subroutine dzwt tabulates the function zwbar(n,m,theta) */
/*     at theta = th in double precision */

/*     input parameters */

/*     nlat      the number of colatitudes including the poles. */
/*            nlat must be an odd integer */

/*     n      the degree (subscript) of zwbar(n,m,theta) */

/*     m      the order (superscript) of zwbar(n,m,theta) */

/*     czw     the fourier coefficients of zwbar(n,m,theta) */
/*             as computed by subroutine zwk. */

/*     output parameter */

/*     zwh     zwbar(m,n,theta) evaluated at theta = th */

    /* Parameter adjustments */
    --czw;

    /* Function Body */
    *zwh = 0.f;
    if (*n <= 0) {
	return 0;
    }
    lc = (*nlat + 1) / 2;
    lq = lc - 1;
    ls = lc - 2;
    cth = cos(*th);
    sth = sin(*th);
    cdt = cth * cth - sth * sth;
    sdt = sth * 2.f * cth;
    lmod = *nlat % 2;
    mmod = *m % 2;
    nmod = *n % 2;
    if (lmod == 0) {
	goto L50;
    }
    if (nmod != 0) {
	goto L1;
    }
    if (mmod != 0) {
	goto L2;
    }

/*     nlat odd  n even  m even */

    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
	*zwh += czw[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L30: */
    }
    return 0;

/*     nlat odd  n even  m odd */

L2:
    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
	*zwh += czw[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L40: */
    }
    return 0;
L1:
    cth = cdt;
    sth = sdt;
    if (mmod != 0) {
	goto L3;
    }

/*     nlat odd  n odd  m even */

    i__1 = ls;
    for (k = 1; k <= i__1; ++k) {
	*zwh += czw[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L10: */
    }
    return 0;

/*     nlat odd  n odd  m odd */

L3:
    *zwh = czw[1] * .5f;
    i__1 = lq;
    for (k = 2; k <= i__1; ++k) {
	*zwh += czw[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L20: */
    }
    *zwh += czw[lc] * .5f * cos((*nlat - 1) * *th);
    return 0;
L50:
    if (nmod != 0) {
	goto L51;
    }
    if (mmod != 0) {
	goto L52;
    }

/*     nlat even  n even  m even */

    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
	*zwh += czw[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L55: */
    }
    return 0;

/*     nlat even  n even  m odd */

L52:
    *zwh = czw[lc] * .5f * cos((*nlat - 1) * *th);
    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
	*zwh += czw[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L60: */
    }
    return 0;
L51:
    cth = cdt;
    sth = sdt;
    if (mmod != 0) {
	goto L53;
    }

/*     nlat even  n odd  m even */

    i__1 = lq;
    for (k = 1; k <= i__1; ++k) {
	*zwh += czw[k + 1] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L65: */
    }
    return 0;

/*     nlat even  n odd  m odd */

L53:
    *zwh = czw[1] * .5f;
    i__1 = lc;
    for (k = 2; k <= i__1; ++k) {
	*zwh += czw[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L70: */
    }
    return 0;
} /* dzwt_ */

/* Subroutine */ int dvbk_(integer *m, integer *n, doublereal *cv, doublereal 
	*work)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer l;
    static doublereal cf, fk, fn;
    static integer ncv, modm, modn;
    static real srnp1;
    extern /* Subroutine */ int dnlfk_(integer *, integer *, doublereal *);

    /* Parameter adjustments */
    --work;
    --cv;

    /* Function Body */
    cv[1] = 0.f;
    if (*n <= 0) {
	return 0;
    }
    fn = (doublereal) (*n);
    srnp1 = sqrt(fn * (fn + 1.f));
    cf = *m * 2.f / srnp1;
    modn = *n % 2;
    modm = *m % 2;
    dnlfk_(m, n, &work[1]);
    if (modn != 0) {
	goto L70;
    }
    ncv = *n / 2;
    if (ncv == 0) {
	return 0;
    }
    fk = 0.f;
    if (modm != 0) {
	goto L60;
    }

/*     n even m even */

    i__1 = ncv;
    for (l = 1; l <= i__1; ++l) {
	fk += 2.f;
	cv[l] = -fk * work[l + 1] / srnp1;
/* L55: */
    }
    return 0;

/*     n even m odd */

L60:
    i__1 = ncv;
    for (l = 1; l <= i__1; ++l) {
	fk += 2.f;
	cv[l] = fk * work[l] / srnp1;
/* L65: */
    }
    return 0;
L70:
    ncv = (*n + 1) / 2;
    fk = -1.f;
    if (modm != 0) {
	goto L80;
    }

/*     n odd m even */

    i__1 = ncv;
    for (l = 1; l <= i__1; ++l) {
	fk += 2.f;
	cv[l] = -fk * work[l] / srnp1;
/* L75: */
    }
    return 0;

/*     n odd m odd */

L80:
    i__1 = ncv;
    for (l = 1; l <= i__1; ++l) {
	fk += 2.f;
	cv[l] = fk * work[l] / srnp1;
/* L85: */
    }
    return 0;
} /* dvbk_ */

/* Subroutine */ int dwbk_(integer *m, integer *n, doublereal *cw, doublereal 
	*work)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer l;
    static doublereal cf, fn;
    static integer modm, modn;
    static doublereal srnp1;
    extern /* Subroutine */ int dnlfk_(integer *, integer *, doublereal *);

    /* Parameter adjustments */
    --work;
    --cw;

    /* Function Body */
    cw[1] = 0.f;
    if (*n <= 0 || *m <= 0) {
	return 0;
    }
    fn = (doublereal) (*n);
    srnp1 = sqrt(fn * (fn + 1.f));
    cf = *m * 2.f / srnp1;
    modn = *n % 2;
    modm = *m % 2;
    dnlfk_(m, n, &work[1]);
    if (*m == 0) {
	goto L50;
    }
    if (modn != 0) {
	goto L30;
    }
    l = *n / 2;
    if (l == 0) {
	goto L50;
    }
    if (modm != 0) {
	goto L20;
    }

/*     n even m even */

    cw[l] = -cf * work[l + 1];
L10:
    --l;
    if (l <= 0) {
	goto L50;
    }
    cw[l] = cw[l + 1] - cf * work[l + 1];
    goto L10;

/*     n even m odd */

L20:
    cw[l] = cf * work[l];
L25:
    --l;
    if (l <= 0) {
	goto L50;
    }
    cw[l] = cw[l + 1] + cf * work[l];
    goto L25;
L30:
    if (modm != 0) {
	goto L40;
    }
    l = (*n - 1) / 2;
    if (l == 0) {
	goto L50;
    }

/*     n odd m even */

    cw[l] = -cf * work[l + 1];
L35:
    --l;
    if (l <= 0) {
	goto L50;
    }
    cw[l] = cw[l + 1] - cf * work[l + 1];
    goto L35;

/*     n odd m odd */

L40:
    l = (*n + 1) / 2;
    cw[l] = cf * work[l];
L45:
    --l;
    if (l <= 0) {
	goto L50;
    }
    cw[l] = cw[l + 1] + cf * work[l];
    goto L45;
L50:
    return 0;
} /* dwbk_ */

/* Subroutine */ int dvbt_(integer *m, integer *n, doublereal *theta, 
	doublereal *cv, doublereal *vh)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal chh, cdt, cth;
    static integer ncv;
    static doublereal sdt, sth;
    static integer mmod, nmod;

    /* Parameter adjustments */
    --cv;

    /* Function Body */
    *vh = 0.f;
    if (*n == 0) {
	return 0;
    }
    cth = cos(*theta);
    sth = sin(*theta);
    cdt = cth * cth - sth * sth;
    sdt = sth * 2.f * cth;
    mmod = *m % 2;
    nmod = *n % 2;
    if (nmod != 0) {
	goto L1;
    }
    cth = cdt;
    sth = sdt;
    if (mmod != 0) {
	goto L2;
    }

/*     n even  m even */

    ncv = *n / 2;
    i__1 = ncv;
    for (k = 1; k <= i__1; ++k) {
	*vh += cv[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L10: */
    }
    return 0;

/*     n even  m odd */

L2:
    ncv = *n / 2;
    i__1 = ncv;
    for (k = 1; k <= i__1; ++k) {
	*vh += cv[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L15: */
    }
    return 0;
L1:
    if (mmod != 0) {
	goto L3;
    }

/*     n odd m even */

    ncv = (*n + 1) / 2;
    i__1 = ncv;
    for (k = 1; k <= i__1; ++k) {
	*vh += cv[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L20: */
    }
    return 0;

/* case m odd and n odd */

L3:
    ncv = (*n + 1) / 2;
    i__1 = ncv;
    for (k = 1; k <= i__1; ++k) {
	*vh += cv[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L25: */
    }
    return 0;
} /* dvbt_ */

/* Subroutine */ int dwbt_(integer *m, integer *n, doublereal *theta, 
	doublereal *cw, doublereal *wh)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal chh, cdt, cth;
    static integer ncw;
    static doublereal sdt, sth;
    static integer mmod, nmod;

    /* Parameter adjustments */
    --cw;

    /* Function Body */
    *wh = 0.f;
    if (*n <= 0 || *m <= 0) {
	return 0;
    }
    cth = cos(*theta);
    sth = sin(*theta);
    cdt = cth * cth - sth * sth;
    sdt = sth * 2.f * cth;
    mmod = *m % 2;
    nmod = *n % 2;
    if (nmod != 0) {
	goto L1;
    }
    if (mmod != 0) {
	goto L2;
    }

/*     n even  m even */

    ncw = *n / 2;
    i__1 = ncw;
    for (k = 1; k <= i__1; ++k) {
	*wh += cw[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L10: */
    }
    return 0;

/*     n even  m odd */

L2:
    ncw = *n / 2;
    i__1 = ncw;
    for (k = 1; k <= i__1; ++k) {
	*wh += cw[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L8: */
    }
    return 0;
L1:
    cth = cdt;
    sth = sdt;
    if (mmod != 0) {
	goto L3;
    }

/*     n odd m even */

    ncw = (*n - 1) / 2;
    i__1 = ncw;
    for (k = 1; k <= i__1; ++k) {
	*wh += cw[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L20: */
    }
    return 0;

/* case m odd and n odd */

L3:
    ncw = (*n + 1) / 2;
    *wh = cw[1] * .5f;
    if (ncw < 2) {
	return 0;
    }
    i__1 = ncw;
    for (k = 2; k <= i__1; ++k) {
	*wh += cw[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L25: */
    }
    return 0;
} /* dwbt_ */

/* Subroutine */ int rabcv_(integer *nlat, integer *nlon, real *abc)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, labc, mmax;
    extern /* Subroutine */ int rabcv1_(integer *, integer *, real *, real *, 
	    real *);


/*     subroutine rabcp computes the coefficients in the recurrence */
/*     relation for the functions vbar(m,n,theta). array abc */
/*     must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 locations. */

    /* Parameter adjustments */
    --abc;

    /* Function Body */
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) / 2;
    iw1 = labc + 1;
    iw2 = iw1 + labc;
    rabcv1_(nlat, nlon, &abc[1], &abc[iw1], &abc[iw2]);
    return 0;
} /* rabcv_ */

/* Subroutine */ int rabcv1_(integer *nlat, integer *nlon, real *a, real *b, 
	real *c__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer m, n;
    static real cn, fm, fn;
    static integer ns;
    static real tm, tn;
    static integer mp1, np1, mp3;
    static real tpn, fnmm, fnpm;
    static integer mmax;
    static real temp;


/*     coefficients a, b, and c for computing vbar(m,n,theta) are */
/*     stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1 */

    /* Parameter adjustments */
    --c__;
    --b;
    --a;

    /* Function Body */
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
    if (mmax < 3) {
	return 0;
    }
    i__1 = mmax;
    for (mp1 = 3; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	ns = (m - 2) * (*nlat + *nlat - m - 1) / 2 + 1;
	fm = (real) m;
	tm = fm + fm;
	temp = tm * (tm - 1.f);
	tpn = (fm - 2.f) * (fm - 1.f) / (fm * (fm + 1.f));
	a[ns] = sqrt(tpn * (tm + 1.f) * (tm - 2.f) / temp);
	c__[ns] = sqrt(2.f / temp);
	if (m == *nlat - 1) {
	    goto L215;
	}
	++ns;
	temp = tm * (tm + 1.f);
	tpn = (fm - 1.f) * fm / ((fm + 1.f) * (fm + 2.f));
	a[ns] = sqrt(tpn * (tm + 3.f) * (tm - 2.f) / temp);
	c__[ns] = sqrt(6.f / temp);
	mp3 = m + 3;
	if (mp3 > *nlat) {
	    goto L215;
	}
	i__2 = *nlat;
	for (np1 = mp3; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    ++ns;
	    fn = (real) n;
	    tn = fn + fn;
	    cn = (tn + 1.f) / (tn - 3.f);
	    tpn = (fn - 2.f) * (fn - 1.f) / (fn * (fn + 1.f));
	    fnpm = fn + fm;
	    fnmm = fn - fm;
	    temp = fnpm * (fnpm - 1.f);
	    a[ns] = sqrt(tpn * cn * (fnpm - 3.f) * (fnpm - 2.f) / temp);
	    b[ns] = sqrt(tpn * cn * fnmm * (fnmm - 1.f) / temp);
	    c__[ns] = sqrt((fnmm + 1.f) * (fnmm + 2.f) / temp);
/* L210: */
	}
L215:
	;
    }
    return 0;
} /* rabcv1_ */

/* Subroutine */ int rabcw_(integer *nlat, integer *nlon, real *abc)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, labc, mmax;
    extern /* Subroutine */ int rabcw1_(integer *, integer *, real *, real *, 
	    real *);


/*     subroutine rabcw computes the coefficients in the recurrence */
/*     relation for the functions wbar(m,n,theta). array abc */
/*     must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 locations. */

    /* Parameter adjustments */
    --abc;

    /* Function Body */
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) / 2;
    iw1 = labc + 1;
    iw2 = iw1 + labc;
    rabcw1_(nlat, nlon, &abc[1], &abc[iw1], &abc[iw2]);
    return 0;
} /* rabcw_ */

/* Subroutine */ int rabcw1_(integer *nlat, integer *nlon, real *a, real *b, 
	real *c__)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer m, n;
    static real cn, fm, fn;
    static integer ns;
    static real tm, tn;
    static integer mp1, np1, mp3;
    static real tph, tpn, fnmm, fnpm;
    static integer mmax;
    static real temp;


/*     coefficients a, b, and c for computing wbar(m,n,theta) are */
/*     stored in location ((m-2)*(nlat+nlat-m-1))/2+n+1 */

    /* Parameter adjustments */
    --c__;
    --b;
    --a;

    /* Function Body */
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
    if (mmax < 4) {
	return 0;
    }
    i__1 = mmax;
    for (mp1 = 4; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	ns = (m - 2) * (*nlat + *nlat - m - 1) / 2 + 1;
	fm = (real) m;
	tm = fm + fm;
	temp = tm * (tm - 1.f);
	tpn = (fm - 2.f) * (fm - 1.f) / (fm * (fm + 1.f));
	tph = fm / (fm - 2.f);
	a[ns] = tph * sqrt(tpn * (tm + 1.f) * (tm - 2.f) / temp);
	c__[ns] = tph * sqrt(2.f / temp);
	if (m == *nlat - 1) {
	    goto L215;
	}
	++ns;
	temp = tm * (tm + 1.f);
	tpn = (fm - 1.f) * fm / ((fm + 1.f) * (fm + 2.f));
	tph = fm / (fm - 2.f);
	a[ns] = tph * sqrt(tpn * (tm + 3.f) * (tm - 2.f) / temp);
	c__[ns] = tph * sqrt(6.f / temp);
	mp3 = m + 3;
	if (mp3 > *nlat) {
	    goto L215;
	}
	i__2 = *nlat;
	for (np1 = mp3; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    ++ns;
	    fn = (real) n;
	    tn = fn + fn;
	    cn = (tn + 1.f) / (tn - 3.f);
	    fnpm = fn + fm;
	    fnmm = fn - fm;
	    temp = fnpm * (fnpm - 1.f);
	    tpn = (fn - 2.f) * (fn - 1.f) / (fn * (fn + 1.f));
	    tph = fm / (fm - 2.f);
	    a[ns] = tph * sqrt(tpn * cn * (fnpm - 3.f) * (fnpm - 2.f) / temp);
	    b[ns] = sqrt(tpn * cn * fnmm * (fnmm - 1.f) / temp);
	    c__[ns] = tph * sqrt((fnmm + 1.f) * (fnmm + 2.f) / temp);
/* L210: */
	}
L215:
	;
    }
    return 0;
} /* rabcw1_ */

/* Subroutine */ int vtinit_(integer *nlat, integer *nlon, real *wvbin, 
	doublereal *dwork)
{
    static integer iw1, imid;
    extern /* Subroutine */ int vtini1_(integer *, integer *, integer *, real 
	    *, real *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --dwork;
    --wvbin;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2 */
/*     the length of dwork is nlat+2 */

    vtini1_(nlat, nlon, &imid, &wvbin[1], &wvbin[iw1], &dwork[1], &dwork[*
	    nlat / 2 + 2]);
    return 0;
} /* vtinit_ */

/* Subroutine */ int vtini1_(integer *nlat, integer *nlon, integer *imid, 
	real *vb, real *abc, doublereal *cvb, doublereal *work)
{
    /* System generated locals */
    integer vb_dim1, vb_dim2, vb_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static integer i__, m, n;
    static doublereal pi, dt, th;
    static integer mp1, np1;
    static doublereal vbh;
    static integer mdo;
    extern /* Subroutine */ int dvtk_(integer *, integer *, doublereal *, 
	    doublereal *), dvtt_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), rabcv_(integer *, integer *, real *);


/*     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 */
/*     locations where mmax = min0(nlat,(nlon+1)/2) */
/*     cvb and work must each have nlat/2+1 locations */

    /* Parameter adjustments */
    vb_dim1 = *imid;
    vb_dim2 = *nlat;
    vb_offset = 1 + vb_dim1 * (1 + vb_dim2);
    vb -= vb_offset;
    --abc;
    --cvb;
    --work;

    /* Function Body */
    pi = atan(1.) * 4.f;
    dt = pi / (*nlat - 1);
/* Computing MIN */
    i__1 = min(2,*nlat), i__2 = (*nlon + 1) / 2;
    mdo = min(i__1,i__2);
    i__1 = mdo;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    dvtk_(&m, &n, &cvb[1], &work[1]);
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		th = (i__ - 1) * dt;
		dvtt_(&m, &n, &th, &cvb[1], &vbh);
		vb[i__ + (np1 + mp1 * vb_dim2) * vb_dim1] = vbh;
/* L165: */
	    }
/* L160: */
	}
    }
    rabcv_(nlat, nlon, &abc[1]);
    return 0;
} /* vtini1_ */

/* Subroutine */ int wtinit_(integer *nlat, integer *nlon, real *wwbin, 
	doublereal *dwork)
{
    static integer iw1, imid;
    extern /* Subroutine */ int wtini1_(integer *, integer *, integer *, real 
	    *, real *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --dwork;
    --wwbin;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2 */
/*     the length of dwork is nlat+2 */

    wtini1_(nlat, nlon, &imid, &wwbin[1], &wwbin[iw1], &dwork[1], &dwork[*
	    nlat / 2 + 2]);
    return 0;
} /* wtinit_ */

/* Subroutine */ int wtini1_(integer *nlat, integer *nlon, integer *imid, 
	real *wb, real *abc, doublereal *cwb, doublereal *work)
{
    /* System generated locals */
    integer wb_dim1, wb_dim2, wb_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static integer i__, m, n;
    static doublereal dt, pi, th;
    static integer mp1, np1, mdo;
    static doublereal wbh;
    extern /* Subroutine */ int dwtk_(integer *, integer *, doublereal *, 
	    doublereal *), dwtt_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), rabcw_(integer *, integer *, real *);


/*     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 */
/*     locations where mmax = min0(nlat,(nlon+1)/2) */
/*     cwb and work must each have nlat/2+1 locations */

    /* Parameter adjustments */
    wb_dim1 = *imid;
    wb_dim2 = *nlat;
    wb_offset = 1 + wb_dim1 * (1 + wb_dim2);
    wb -= wb_offset;
    --abc;
    --cwb;
    --work;

    /* Function Body */
    pi = atan(1.) * 4.f;
    dt = pi / (*nlat - 1);
/* Computing MIN */
    i__1 = min(3,*nlat), i__2 = (*nlon + 1) / 2;
    mdo = min(i__1,i__2);
    if (mdo < 2) {
	return 0;
    }
    i__1 = mdo;
    for (mp1 = 2; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    dwtk_(&m, &n, &cwb[1], &work[1]);
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		th = (i__ - 1) * dt;
		dwtt_(&m, &n, &th, &cwb[1], &wbh);
		wb[i__ + (np1 + m * wb_dim2) * wb_dim1] = wbh;
/* L165: */
	    }
/* L160: */
	}
    }
    rabcw_(nlat, nlon, &abc[1]);
    return 0;
} /* wtini1_ */

/* Subroutine */ int vtgint_(integer *nlat, integer *nlon, doublereal *theta, 
	real *wvbin, doublereal *work)
{
    static integer iw1, imid;
    extern /* Subroutine */ int vtgit1_(integer *, integer *, integer *, 
	    doublereal *, real *, real *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --work;
    --wvbin;
    --theta;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     theta is a double precision array with (nlat+1)/2 locations */
/*     nlat is the maximum value of n+1 */
/*     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2 */
/*     the length of work is nlat+2 */

    vtgit1_(nlat, nlon, &imid, &theta[1], &wvbin[1], &wvbin[iw1], &work[1], &
	    work[*nlat / 2 + 2]);
    return 0;
} /* vtgint_ */

/* Subroutine */ int vtgit1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *theta, real *vb, real *abc, doublereal *cvb, doublereal *
	work)
{
    /* System generated locals */
    integer vb_dim1, vb_dim2, vb_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, m, n, mp1, np1;
    static doublereal vbh;
    static integer mdo;
    extern /* Subroutine */ int dvtk_(integer *, integer *, doublereal *, 
	    doublereal *), dvtt_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), rabcv_(integer *, integer *, real *);


/*     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 */
/*     locations where mmax = min0(nlat,(nlon+1)/2) */
/*     cvb and work must each have nlat/2+1   locations */

    /* Parameter adjustments */
    vb_dim1 = *imid;
    vb_dim2 = *nlat;
    vb_offset = 1 + vb_dim1 * (1 + vb_dim2);
    vb -= vb_offset;
    --theta;
    --abc;
    --cvb;
    --work;

    /* Function Body */
/* Computing MIN */
    i__1 = min(2,*nlat), i__2 = (*nlon + 1) / 2;
    mdo = min(i__1,i__2);
    i__1 = mdo;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    dvtk_(&m, &n, &cvb[1], &work[1]);
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		dvtt_(&m, &n, &theta[i__], &cvb[1], &vbh);
		vb[i__ + (np1 + mp1 * vb_dim2) * vb_dim1] = vbh;
/* L165: */
	    }
/* L160: */
	}
    }
    rabcv_(nlat, nlon, &abc[1]);
    return 0;
} /* vtgit1_ */

/* Subroutine */ int wtgint_(integer *nlat, integer *nlon, doublereal *theta, 
	real *wwbin, doublereal *work)
{
    static integer iw1, imid;
    extern /* Subroutine */ int wtgit1_(integer *, integer *, integer *, 
	    doublereal *, real *, real *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --work;
    --wwbin;
    --theta;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     theta is a double precision array with (nlat+1)/2 locations */
/*     nlat is the maximum value of n+1 */
/*     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2 */
/*     the length of work is nlat+2 */

    wtgit1_(nlat, nlon, &imid, &theta[1], &wwbin[1], &wwbin[iw1], &work[1], &
	    work[*nlat / 2 + 2]);
    return 0;
} /* wtgint_ */

/* Subroutine */ int wtgit1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *theta, real *wb, real *abc, doublereal *cwb, doublereal *
	work)
{
    /* System generated locals */
    integer wb_dim1, wb_dim2, wb_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, m, n, mp1, np1, mdo;
    static doublereal wbh;
    extern /* Subroutine */ int dwtk_(integer *, integer *, doublereal *, 
	    doublereal *), dwtt_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), rabcw_(integer *, integer *, real *);


/*     abc must have 3*((nlat-3)*nlat+2)/2 locations */
/*     cwb and work must each have nlat/2+1 locations */

    /* Parameter adjustments */
    wb_dim1 = *imid;
    wb_dim2 = *nlat;
    wb_offset = 1 + wb_dim1 * (1 + wb_dim2);
    wb -= wb_offset;
    --theta;
    --abc;
    --cwb;
    --work;

    /* Function Body */
/* Computing MIN */
    i__1 = min(3,*nlat), i__2 = (*nlon + 1) / 2;
    mdo = min(i__1,i__2);
    if (mdo < 2) {
	return 0;
    }
    i__1 = mdo;
    for (mp1 = 2; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    dwtk_(&m, &n, &cwb[1], &work[1]);
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		dwtt_(&m, &n, &theta[i__], &cwb[1], &wbh);
		wb[i__ + (np1 + m * wb_dim2) * wb_dim1] = wbh;
/* L165: */
	    }
/* L160: */
	}
    }
    rabcw_(nlat, nlon, &abc[1]);
    return 0;
} /* wtgit1_ */

/* Subroutine */ int dvtk_(integer *m, integer *n, doublereal *cv, doublereal 
	*work)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer l;
    static doublereal cf, fk, fn;
    static integer ncv, modm, modn;
    static doublereal srnp1;
    extern /* Subroutine */ int dnlfk_(integer *, integer *, doublereal *);

    /* Parameter adjustments */
    --work;
    --cv;

    /* Function Body */
    cv[1] = 0.f;
    if (*n <= 0) {
	return 0;
    }
    fn = (doublereal) (*n);
    srnp1 = sqrt(fn * (fn + 1.f));
    cf = *m * 2.f / srnp1;
    modn = *n % 2;
    modm = *m % 2;
    dnlfk_(m, n, &work[1]);
    if (modn != 0) {
	goto L70;
    }
    ncv = *n / 2;
    if (ncv == 0) {
	return 0;
    }
    fk = 0.f;
    if (modm != 0) {
	goto L60;
    }

/*     n even m even */

    i__1 = ncv;
    for (l = 1; l <= i__1; ++l) {
	fk += 2.f;
	cv[l] = -fk * fk * work[l + 1] / srnp1;
/* L55: */
    }
    return 0;

/*     n even m odd */

L60:
    i__1 = ncv;
    for (l = 1; l <= i__1; ++l) {
	fk += 2.f;
	cv[l] = -fk * fk * work[l] / srnp1;
/* L65: */
    }
    return 0;
L70:
    ncv = (*n + 1) / 2;
    fk = -1.f;
    if (modm != 0) {
	goto L80;
    }

/*     n odd m even */

    i__1 = ncv;
    for (l = 1; l <= i__1; ++l) {
	fk += 2.f;
	cv[l] = -fk * fk * work[l] / srnp1;
/* L75: */
    }
    return 0;

/*     n odd m odd */

L80:
    i__1 = ncv;
    for (l = 1; l <= i__1; ++l) {
	fk += 2.f;
	cv[l] = -fk * fk * work[l] / srnp1;
/* L85: */
    }
    return 0;
} /* dvtk_ */

/* Subroutine */ int dwtk_(integer *m, integer *n, doublereal *cw, doublereal 
	*work)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer l;
    static doublereal cf, fn;
    static integer modm, modn;
    static doublereal srnp1;
    extern /* Subroutine */ int dnlfk_(integer *, integer *, doublereal *);

    /* Parameter adjustments */
    --work;
    --cw;

    /* Function Body */
    cw[1] = 0.f;
    if (*n <= 0 || *m <= 0) {
	return 0;
    }
    fn = (doublereal) (*n);
    srnp1 = sqrt(fn * (fn + 1.f));
    cf = *m * 2.f / srnp1;
    modn = *n % 2;
    modm = *m % 2;
    dnlfk_(m, n, &work[1]);
    if (*m == 0) {
	goto L50;
    }
    if (modn != 0) {
	goto L30;
    }
    l = *n / 2;
    if (l == 0) {
	goto L50;
    }
    if (modm != 0) {
	goto L20;
    }

/*     n even m even */

    cw[l] = -cf * work[l + 1];
L10:
    --l;
    if (l <= 0) {
	goto L50;
    }
    cw[l] = cw[l + 1] - cf * work[l + 1];
    cw[l + 1] = (l + l + 1) * cw[l + 1];
    goto L10;

/*     n even m odd */

L20:
    cw[l] = cf * work[l];
L25:
    --l;
    if (l < 0) {
	goto L50;
    } else if (l == 0) {
	goto L27;
    } else {
	goto L26;
    }
L26:
    cw[l] = cw[l + 1] + cf * work[l];
L27:
    cw[l + 1] = -(l + l + 1) * cw[l + 1];
    goto L25;
L30:
    if (modm != 0) {
	goto L40;
    }
    l = (*n - 1) / 2;
    if (l == 0) {
	goto L50;
    }

/*     n odd m even */

    cw[l] = -cf * work[l + 1];
L35:
    --l;
    if (l < 0) {
	goto L50;
    } else if (l == 0) {
	goto L37;
    } else {
	goto L36;
    }
L36:
    cw[l] = cw[l + 1] - cf * work[l + 1];
L37:
    cw[l + 1] = (l + l + 2) * cw[l + 1];
    goto L35;

/*     n odd m odd */

L40:
    l = (*n + 1) / 2;
    cw[l] = cf * work[l];
L45:
    --l;
    if (l < 0) {
	goto L50;
    } else if (l == 0) {
	goto L47;
    } else {
	goto L46;
    }
L46:
    cw[l] = cw[l + 1] + cf * work[l];
L47:
    cw[l + 1] = -(l + l) * cw[l + 1];
    goto L45;
L50:
    return 0;
} /* dwtk_ */

/* Subroutine */ int dvtt_(integer *m, integer *n, doublereal *theta, 
	doublereal *cv, doublereal *vh)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal chh, cdt, cth;
    static integer ncv;
    static doublereal sdt, sth;
    static integer mmod, nmod;

    /* Parameter adjustments */
    --cv;

    /* Function Body */
    *vh = 0.f;
    if (*n == 0) {
	return 0;
    }
    cth = cos(*theta);
    sth = sin(*theta);
    cdt = cth * cth - sth * sth;
    sdt = sth * 2.f * cth;
    mmod = *m % 2;
    nmod = *n % 2;
    if (nmod != 0) {
	goto L1;
    }
    cth = cdt;
    sth = sdt;
    if (mmod != 0) {
	goto L2;
    }

/*     n even  m even */

    ncv = *n / 2;
    i__1 = ncv;
    for (k = 1; k <= i__1; ++k) {
	*vh += cv[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L10: */
    }
    return 0;

/*     n even  m odd */

L2:
    ncv = *n / 2;
    i__1 = ncv;
    for (k = 1; k <= i__1; ++k) {
	*vh += cv[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L15: */
    }
    return 0;
L1:
    if (mmod != 0) {
	goto L3;
    }

/*     n odd m even */

    ncv = (*n + 1) / 2;
    i__1 = ncv;
    for (k = 1; k <= i__1; ++k) {
	*vh += cv[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L20: */
    }
    return 0;

/* case m odd and n odd */

L3:
    ncv = (*n + 1) / 2;
    i__1 = ncv;
    for (k = 1; k <= i__1; ++k) {
	*vh += cv[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L25: */
    }
    return 0;
} /* dvtt_ */

/* Subroutine */ int dwtt_(integer *m, integer *n, doublereal *theta, 
	doublereal *cw, doublereal *wh)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal chh, cdt, cth;
    static integer ncw;
    static doublereal sdt, sth;
    static integer mmod, nmod;

    /* Parameter adjustments */
    --cw;

    /* Function Body */
    *wh = 0.f;
    if (*n <= 0 || *m <= 0) {
	return 0;
    }
    cth = cos(*theta);
    sth = sin(*theta);
    cdt = cth * cth - sth * sth;
    sdt = sth * 2.f * cth;
    mmod = *m % 2;
    nmod = *n % 2;
    if (nmod != 0) {
	goto L1;
    }
    if (mmod != 0) {
	goto L2;
    }

/*     n even  m even */

    ncw = *n / 2;
    i__1 = ncw;
    for (k = 1; k <= i__1; ++k) {
	*wh += cw[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L10: */
    }
    return 0;

/*     n even  m odd */

L2:
    ncw = *n / 2;
    i__1 = ncw;
    for (k = 1; k <= i__1; ++k) {
	*wh += cw[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L8: */
    }
    return 0;
L1:
    cth = cdt;
    sth = sdt;
    if (mmod != 0) {
	goto L3;
    }

/*     n odd m even */

    ncw = (*n - 1) / 2;
    i__1 = ncw;
    for (k = 1; k <= i__1; ++k) {
	*wh += cw[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L20: */
    }
    return 0;

/* case m odd and n odd */

L3:
    ncw = (*n + 1) / 2;
    *wh = 0.f;
    if (ncw < 2) {
	return 0;
    }
    i__1 = ncw;
    for (k = 2; k <= i__1; ++k) {
	*wh += cw[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L25: */
    }
    return 0;
} /* dwtt_ */

/* Subroutine */ int vbgint_(integer *nlat, integer *nlon, doublereal *theta, 
	real *wvbin, doublereal *work)
{
    static integer iw1, imid;
    extern /* Subroutine */ int vbgit1_(integer *, integer *, integer *, 
	    doublereal *, real *, real *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --work;
    --wvbin;
    --theta;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     theta is a double precision array with (nlat+1)/2 locations */
/*     nlat is the maximum value of n+1 */
/*     the length of wvbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2 */
/*     the length of work is nlat+2 */

    vbgit1_(nlat, nlon, &imid, &theta[1], &wvbin[1], &wvbin[iw1], &work[1], &
	    work[*nlat / 2 + 2]);
    return 0;
} /* vbgint_ */

/* Subroutine */ int vbgit1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *theta, real *vb, real *abc, doublereal *cvb, doublereal *
	work)
{
    /* System generated locals */
    integer vb_dim1, vb_dim2, vb_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, m, n, mp1, np1;
    static doublereal vbh;
    static integer mdo;
    extern /* Subroutine */ int dvbk_(integer *, integer *, doublereal *, 
	    doublereal *), dvbt_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), rabcv_(integer *, integer *, real *);


/*     abc must have 3*(max0(mmax-2,0)*(nlat+nlat-mmax-1))/2 */
/*     locations where mmax = min0(nlat,(nlon+1)/2) */
/*     cvb and work must each have nlat/2+1 locations */

    /* Parameter adjustments */
    vb_dim1 = *imid;
    vb_dim2 = *nlat;
    vb_offset = 1 + vb_dim1 * (1 + vb_dim2);
    vb -= vb_offset;
    --theta;
    --abc;
    --cvb;
    --work;

    /* Function Body */
/* Computing MIN */
    i__1 = min(2,*nlat), i__2 = (*nlon + 1) / 2;
    mdo = min(i__1,i__2);
    i__1 = mdo;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    dvbk_(&m, &n, &cvb[1], &work[1]);
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		dvbt_(&m, &n, &theta[i__], &cvb[1], &vbh);
		vb[i__ + (np1 + mp1 * vb_dim2) * vb_dim1] = vbh;
/* L165: */
	    }
/* L160: */
	}
    }
    rabcv_(nlat, nlon, &abc[1]);
    return 0;
} /* vbgit1_ */

/* Subroutine */ int wbgint_(integer *nlat, integer *nlon, doublereal *theta, 
	real *wwbin, doublereal *work)
{
    static integer iw1, imid;
    extern /* Subroutine */ int wbgit1_(integer *, integer *, integer *, 
	    doublereal *, real *, real *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --work;
    --wwbin;
    --theta;

    /* Function Body */
    imid = (*nlat + 1) / 2;
    iw1 = (*nlat << 1) * imid + 1;

/*     theta is a double precision array with (nlat+1)/2 locations */
/*     nlat is the maximum value of n+1 */
/*     the length of wwbin is 2*nlat*imid+3*((nlat-3)*nlat+2)/2 */
/*     the length of work is nlat+2 */

    wbgit1_(nlat, nlon, &imid, &theta[1], &wwbin[1], &wwbin[iw1], &work[1], &
	    work[*nlat / 2 + 2]);
    return 0;
} /* wbgint_ */

/* Subroutine */ int wbgit1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *theta, real *wb, real *abc, doublereal *cwb, doublereal *
	work)
{
    /* System generated locals */
    integer wb_dim1, wb_dim2, wb_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, m, n, mp1, np1, mdo;
    static doublereal wbh;
    extern /* Subroutine */ int dwbk_(integer *, integer *, doublereal *, 
	    doublereal *), dwbt_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), rabcw_(integer *, integer *, real *);


/*     abc must have 3*((nlat-3)*nlat+2)/2 locations */
/*     cwb and work must each have nlat/2+1 locations */

    /* Parameter adjustments */
    wb_dim1 = *imid;
    wb_dim2 = *nlat;
    wb_offset = 1 + wb_dim1 * (1 + wb_dim2);
    wb -= wb_offset;
    --theta;
    --abc;
    --cwb;
    --work;

    /* Function Body */
/* Computing MIN */
    i__1 = min(3,*nlat), i__2 = (*nlon + 1) / 2;
    mdo = min(i__1,i__2);
    if (mdo < 2) {
	return 0;
    }
    i__1 = mdo;
    for (mp1 = 2; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    n = np1 - 1;
	    dwbk_(&m, &n, &cwb[1], &work[1]);
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		dwbt_(&m, &n, &theta[i__], &cwb[1], &wbh);
		wb[i__ + (np1 + m * wb_dim2) * wb_dim1] = wbh;
/* L165: */
	    }
/* L160: */
	}
    }
    rabcw_(nlat, nlon, &abc[1]);
    return 0;
} /* wbgit1_ */

