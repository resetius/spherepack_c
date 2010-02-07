/* shigs.f -- translated by f2c (version 20061008).
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



/* ... file shigs.f */

/*     this file contains code and documentation for subroutine shigs */

/* ... files which must be loaded with shigs.f */

/*     sphcom.f, hrfft.f, gaqd.f */

/*     3/6/98 */

/* *** shigs is functionally the same as shagsi or shsgsi.  It */
/*     it included in spherepack3.0 because legacy codes, using */
/*     the older version of spherepack, call shigs to initialize */
/*     the saved work space wshigs for either shags or shsgs */
/*     Its arguments are identical to those of shagsi or shsgsi. */

/* **************************************************************** */

/*     subroutine shigs(nlat,nlon,wshigs,lshigs,work,lwork,dwork,ldwork, */
/*    +                 ierror) */

/*     subroutine shigs initializes the array wshigs which can then */
/*     be used repeatedly by subroutines shags,shsgs. it precomputes */
/*     and stores in wshigs quantities such as gaussian weights, */
/*     legendre polynomial coefficients, and fft trigonometric tables. */

/*     input parameters */

/*     nlat   the number of points in the gaussian colatitude grid on the */
/*            full sphere. these lie in the interval (0,pi) and are compu */
/*            in radians in theta(1),...,theta(nlat) by subroutine gaqd. */
/*            if nlat is odd the equator will be included as the grid poi */
/*            theta((nlat+1)/2).  if nlat is even the equator will be */
/*            excluded as a grid point and will lie half way between */
/*            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3. */
/*            note: on the half sphere, the number of grid points in the */
/*            colatitudinal direction is nlat/2 if nlat is even or */
/*            (nlat+1)/2 if nlat is odd. */

/*     nlon   the number of distinct londitude points.  nlon determines */
/*            the grid increment in longitude as 2*pi/nlon. for example */
/*            nlon = 72 for a five degree grid. nlon must be greater */
/*            than or equal to 4. the efficiency of the computation is */
/*            improved when nlon is a product of small prime numbers. */

/*     wshigs an array which must be initialized by subroutine shigs . */
/*            once initialized, wshigs can be used repeatedly by shigs */
/*            as long as nlat and nlon remain unchanged.  wshigs must */
/*            not be altered between calls of shigs. */

/*     lshigs the dimension of the array wshigs as it appears in the */
/*            program that calls shigs. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshigs must be at least */

/*            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15 */

/*     work   a real work space which need not be saved */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls shigs. lwork must be at least */
/*            4*nlat*(nlat+2)+2 in the routine calling shigs */

/*     dwork   a double precision work array that does not have to be saved. */

/*     ldwork  the length of dwork in the calling routine.  ldwork must */
/*             be at least nlat*(nlat+4) */

/*     output parameter */

/*     wshags an array which must be initialized before calling shags or */
/*            once initialized, wshags can be used repeatedly by shags or */
/*            as long as nlat and nlon remain unchanged.  wshags must not */
/*            altered between calls of shasc. */

/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of lshags */
/*            = 4  error in the specification of lwork */
/*            = 5  error in the specification of ldwork */
/*            = 6  failure in gaqd to compute gaussian points */
/*                 (due to failure in eigenvalue routine) */


/* **************************************************************** */

/* Subroutine */ int shigs_(integer *nlat, integer *nlon, real *wshigs, 
	integer *lshigs, real *work, integer *lwork, doublereal *dwork, 
	integer *ldwork, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer l, l1, l2, lp, late, ipmnf;
    extern /* Subroutine */ int shigsp_(integer *, integer *, real *, integer 
	    *, doublereal *, integer *, integer *), shigss1_(integer *, 
	    integer *, integer *, real *, real *, real *);


/*     this subroutine must be called before calling shags or shsgs with */
/*     fixed nlat,nlon. it precomputes the gaussian weights, points */
/*     and all necessary legendre polys and stores them in wshigs. */
/*     these quantities must be preserved when calling shsgs or shags */
/*     repeatedly with fixed nlat,nlon. */

    /* Parameter adjustments */
    --wshigs;
    --work;
    --dwork;

    /* Function Body */
    *ierror = 1;
    if (*nlat < 3) {
	return 0;
    }
    *ierror = 2;
    if (*nlon < 4) {
	return 0;
    }
/*     set triangular truncation limit for spherical harmonic basis */
/* Computing MIN */
    i__1 = (*nlon + 2) / 2;
    l = min(i__1,*nlat);
/*     set equator or nearest point (if excluded) pointer */
    late = (*nlat + 1) / 2;
    l1 = l;
    l2 = late;
/*     check permanent work space length */
    *ierror = 3;
    lp = *nlat * ((l1 + l2) * 3 - 2) + (l1 - 1) * (l2 * ((*nlat << 1) - l1) - 
	    l1 * 3) / 2 + *nlon + 15;
    if (*lshigs < lp) {
	return 0;
    }
    *ierror = 4;
/*     check temporary work space */
    if (*lwork < (*nlat << 2) * (*nlat + 2) + 2) {
	return 0;
    }
/*     check temp double precision space */
    *ierror = 5;
    if (*ldwork < *nlat * (*nlat + 4)) {
	return 0;
    }
    *ierror = 0;
/*     set preliminary quantites needed to compute and store legendre polys */
    shigsp_(nlat, nlon, &wshigs[1], lshigs, &dwork[1], ldwork, ierror);
    if (*ierror != 0) {
	return 0;
    }
/*     set legendre poly pointer in wshigs */
    ipmnf = *nlat + (*nlat << 1) * late + (l * (l - 1) / 2 + (*nlat - l) * (l 
	    - 1)) * 3 + *nlon + 16;
    shigss1_(nlat, &l, &late, &wshigs[1], &work[1], &wshigs[ipmnf]);
    return 0;
} /* shigs_ */

/* Subroutine */ int shigss1_(integer *nlat, integer *l, integer *late, real *
	w, real *pmn, real *pmnf)
{
    /* System generated locals */
    integer pmn_dim1, pmn_dim2, pmn_offset, pmnf_dim1, pmnf_offset, i__1, 
	    i__2, i__3;

    /* Local variables */
    static integer i__, j, k, m, km, mn, mp1, np1, mml1, mode;
    extern /* Subroutine */ int legin_(integer *, integer *, integer *, 
	    integer *, real *, real *, integer *);

/*     compute and store legendre polys for i=1,...,late,m=0,...,l-1 */
/*     and n=m,...,l-1 */
    /* Parameter adjustments */
    pmnf_dim1 = *late;
    pmnf_offset = 1 + pmnf_dim1;
    pmnf -= pmnf_offset;
    pmn_dim1 = *nlat;
    pmn_dim2 = *late;
    pmn_offset = 1 + pmn_dim1 * (1 + pmn_dim2);
    pmn -= pmn_offset;
    --w;

    /* Function Body */
    i__1 = *nlat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *late;
	for (j = 1; j <= i__2; ++j) {
	    for (k = 1; k <= 3; ++k) {
		pmn[i__ + (j + k * pmn_dim2) * pmn_dim1] = 0.f;
	    }
	}
    }
    i__1 = *l;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	mml1 = m * ((*nlat << 1) - m - 1) / 2;
/*     compute pmn for n=m,...,nlat-1 and i=1,...,(l+1)/2 */
	mode = 0;
	legin_(&mode, l, nlat, &m, &w[1], &pmn[pmn_offset], &km);
/*     store above in pmnf */
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    mn = mml1 + np1;
	    i__3 = *late;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		pmnf[i__ + mn * pmnf_dim1] = pmn[np1 + (i__ + km * pmn_dim2) *
			 pmn_dim1];
/* L102: */
	    }
/* L101: */
	}
/* L100: */
    }
    return 0;
} /* shigss1_ */

/* Subroutine */ int shigsp_(integer *nlat, integer *nlon, real *wshigs, 
	integer *lshigs, doublereal *dwork, integer *ldwork, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer l, i1, i2, i3, l1, l2, i4, i5, i6, i7, iw, late, idth, 
	    idwts;
    extern /* Subroutine */ int shigsp1_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     doublereal *, doublereal *, doublereal *, integer *);

    /* Parameter adjustments */
    --wshigs;
    --dwork;

    /* Function Body */
    *ierror = 1;
    if (*nlat < 3) {
	return 0;
    }
    *ierror = 2;
    if (*nlon < 4) {
	return 0;
    }
/*     set triangular truncation limit for spherical harmonic basis */
/* Computing MIN */
    i__1 = (*nlon + 2) / 2;
    l = min(i__1,*nlat);
/*     set equator or nearest point (if excluded) pointer */
    late = (*nlat + *nlat % 2) / 2;
    l1 = l;
    l2 = late;
    *ierror = 3;
/*     check permanent work space length */
    if (*lshigs < *nlat * ((l2 << 1) + l1 * 3 - 2) + l1 * 3 * (1 - l1) / 2 + *
	    nlon + 15) {
	return 0;
    }
    *ierror = 4;
/*     if (lwork.lt.4*nlat*(nlat+2)+2) return */
    if (*ldwork < *nlat * (*nlat + 4)) {
	return 0;
    }
    *ierror = 0;
/*     set pointers */
    i1 = 1;
    i2 = i1 + *nlat;
    i3 = i2 + *nlat * late;
    i4 = i3 + *nlat * late;
    i5 = i4 + l * (l - 1) / 2 + (*nlat - l) * (l - 1);
    i6 = i5 + l * (l - 1) / 2 + (*nlat - l) * (l - 1);
    i7 = i6 + l * (l - 1) / 2 + (*nlat - l) * (l - 1);
/*     set indices in temp work for double precision gaussian wts and pts */
    idth = 1;
/*     idwts = idth+2*nlat */
/*     iw = idwts+2*nlat */
    idwts = idth + *nlat;
    iw = idwts + *nlat;
    shigsp1_(nlat, nlon, &l, &late, &wshigs[i1], &wshigs[i2], &wshigs[i3], &
	    wshigs[i4], &wshigs[i5], &wshigs[i6], &wshigs[i7], &dwork[idth], &
	    dwork[idwts], &dwork[iw], ierror);
    if (*ierror != 0) {
	*ierror = 5;
    }
    return 0;
} /* shigsp_ */

/* Subroutine */ int shigsp1_(integer *nlat, integer *nlon, integer *l, 
	integer *late, real *wts, real *p0n, real *p1n, real *abel, real *
	bbel, real *cbel, real *wfft, doublereal *dtheta, doublereal *dwts, 
	doublereal *work, integer *ier)
{
    /* System generated locals */
    integer p0n_dim1, p0n_offset, p1n_dim1, p1n_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, m, n;
    static doublereal pb;
    static integer lw, np1, imn;
    extern /* Subroutine */ int gaqd_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    static integer mlim;
    extern /* Subroutine */ int dnlfk_(integer *, integer *, doublereal *), 
	    dnlft_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *), hrffti_(integer *, real *);

    /* Parameter adjustments */
    --dwts;
    --dtheta;
    --wts;
    p1n_dim1 = *nlat;
    p1n_offset = 1 + p1n_dim1;
    p1n -= p1n_offset;
    p0n_dim1 = *nlat;
    p0n_offset = 1 + p0n_dim1;
    p0n -= p0n_offset;
    --abel;
    --bbel;
    --cbel;
    --wfft;
    --work;

    /* Function Body */
    hrffti_(nlon, &wfft[1]);
/*     compute double precision gaussian points and weights */
/*     lw = 4*nlat*(nlat+2) */
    lw = *nlat * (*nlat + 2);
    gaqd_(nlat, &dtheta[1], &dwts[1], &work[1], &lw, ier);
    if (*ier != 0) {
	return 0;
    }
/*     store gaussian weights single precision to save computation */
/*     in inner loops in analysis */
    i__1 = *nlat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wts[i__] = dwts[i__];
/* L100: */
    }
/*     initialize p0n,p1n using double precision dnlfk,dnlft */
    i__1 = *nlat;
    for (np1 = 1; np1 <= i__1; ++np1) {
	i__2 = *late;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    p0n[np1 + i__ * p0n_dim1] = 0.f;
	    p1n[np1 + i__ * p1n_dim1] = 0.f;
/* L101: */
	}
    }
/*     compute m=n=0 legendre polynomials for all theta(i) */
    np1 = 1;
    n = 0;
    m = 0;
    dnlfk_(&m, &n, &work[1]);
    i__2 = *late;
    for (i__ = 1; i__ <= i__2; ++i__) {
	dnlft_(&m, &n, &dtheta[i__], &work[1], &pb);
	p0n[i__ * p0n_dim1 + 1] = pb;
/* L103: */
    }
/*     compute p0n,p1n for all theta(i) when n.gt.0 */
    i__2 = *nlat;
    for (np1 = 2; np1 <= i__2; ++np1) {
	n = np1 - 1;
	m = 0;
	dnlfk_(&m, &n, &work[1]);
	i__1 = *late;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dnlft_(&m, &n, &dtheta[i__], &work[1], &pb);
	    p0n[np1 + i__ * p0n_dim1] = pb;
/* L105: */
	}
/*     compute m=1 legendre polynomials for all n and theta(i) */
	m = 1;
	dnlfk_(&m, &n, &work[1]);
	i__1 = *late;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dnlft_(&m, &n, &dtheta[i__], &work[1], &pb);
	    p1n[np1 + i__ * p1n_dim1] = pb;
/* L106: */
	}
/* L104: */
    }

/*     compute and store swarztrauber recursion coefficients */
/*     for 2.le.m.le.n and 2.le.n.le.nlat in abel,bbel,cbel */
    i__2 = *nlat;
    for (n = 2; n <= i__2; ++n) {
	mlim = min(n,*l);
	i__1 = mlim;
	for (m = 2; m <= i__1; ++m) {
	    imn = (n - 1) * (n - 2) / 2 + m - 1;
	    if (n >= *l) {
		imn = *l * (*l - 1) / 2 + (n - *l - 1) * (*l - 1) + m - 1;
	    }
	    abel[imn] = sqrt((real) (((n << 1) + 1) * (m + n - 2) * (m + n - 
		    3)) / (real) (((n << 1) - 3) * (m + n - 1) * (m + n)));
	    bbel[imn] = sqrt((real) (((n << 1) + 1) * (n - m - 1) * (n - m)) /
		     (real) (((n << 1) - 3) * (m + n - 1) * (m + n)));
	    cbel[imn] = sqrt((real) ((n - m + 1) * (n - m + 2)) / (real) ((n 
		    + m - 1) * (n + m)));
/* L107: */
	}
    }
    return 0;
} /* shigsp1_ */

