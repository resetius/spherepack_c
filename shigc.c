/* shigc.f -- translated by f2c (version 20061008).
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



/* ... file shigc.f */

/*     this file contains code and documentation for subroutine shigc */

/* ... files which must be loaded with shigc.f */

/*     sphcom.f, hrfft.f, gaqd.f */

/*     3/6/98 */

/* *** shigc is functionally the same as shagci or shsgci.  It */
/*     it included in spherepack3.0 because legacy codes, using */
/*     the older version of spherepack call shigc to initialize */
/*     the saved work space wshigc, for either shagc or shsgc */

/*     subroutine shigc(nlat,nlon,wshigc,lshigc,dwork,ldwork,ierror) */

/*     subroutine shigc initializes the array wshigc which can then */
/*     be used repeatedly by subroutines shsgc or shagc. it precomputes */
/*     and stores in wshigc quantities such as gaussian weights, */
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

/*     wshigc an array which must be initialized by subroutine shigc. */
/*            once initialized, wshigc can be used repeatedly by shsgc */
/*            or shagc as long as nlat and nlon remain unchanged.  wshigc */
/*            must not be altered between calls of shsgc or shagc. */

/*     lshigc the dimension of the array wshigc as it appears in the */
/*            program that calls shsgc or shagc. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshigc must be at least */

/*                  nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15 */

/*     dwork  a double precision work array that does not have to be saved. */

/*     ldwork the dimension of the array dwork as it appears in the */
/*            program that calls shigc. ldwork must be at least */

/*               nlat*(nlat+4) */

/*     output parameter */

/*     wshigc an array which must be initialized before calling shsgc or shagc. */
/*            once initialized, wshigc can be used repeatedly by shsgc or shagc */
/*            as long as nlat and nlon remain unchanged.  wshigc must not */
/*            altered between calls of shsgc or shagc */

/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of lshigc */
/*            = 4  error in the specification of ldwork */
/*            = 5  failure in gaqd to compute gaussian points */
/*                 (due to failure in eigenvalue routine) */


/* **************************************************************** */
/* Subroutine */ int shigc_(integer *nlat, integer *nlon, doublereal *wshigc, 
	integer *lshigc, doublereal *dwork, integer *ldwork, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer l, i1, i2, i3, l1, l2, i4, i5, i6, i7, iw, late, idth, idwts;
    extern /* Subroutine */ int shigc1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);

/*     this subroutine must be called before calling shsgc/shagc with */
/*     fixed nlat,nlon. it precomputes quantites such as the gaussian */
/*     points and weights, m=0,m=1 legendre polynomials, recursion */
/*     recursion coefficients. */
    /* Parameter adjustments */
    --wshigc;
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
    if (*lshigc < *nlat * ((l2 << 1) + l1 * 3 - 2) + l1 * 3 * (1 - l1) / 2 + *
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
    shigc1_(nlat, nlon, &l, &late, &wshigc[i1], &wshigc[i2], &wshigc[i3], &
	    wshigc[i4], &wshigc[i5], &wshigc[i6], &wshigc[i7], &dwork[idth], &
	    dwork[idwts], &dwork[iw], ierror);
    if (*ierror != 0) {
	*ierror = 5;
    }
    return 0;
} /* shigc_ */

/* Subroutine */ int shigc1_(integer *nlat, integer *nlon, integer *l, 
	integer *late, doublereal *wts, doublereal *p0n, doublereal *p1n, 
	doublereal *abel, doublereal *bbel, doublereal *cbel, doublereal *
	wfft, doublereal *dtheta, doublereal *dwts, doublereal *work, integer 
	*ier)
{
    /* System generated locals */
    integer p0n_dim1, p0n_offset, p1n_dim1, p1n_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, m, n;
    doublereal pb;
    integer lw, np1, imn;
    extern /* Subroutine */ int gaqd_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    integer mlim;
    extern /* Subroutine */ int dnlfk_(integer *, integer *, doublereal *), 
	    dnlft_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *), hrffti_(integer *, doublereal *);

/*     compute the nlat  gaussian points and weights, the */
/*     m=0,1 legendre polys for gaussian points and all n, */
/*     and the legendre recursion coefficients */
/*     define index function used in storing */
/*     arrays for recursion coefficients (functions of (m,n)) */
/*     the index function indx(m,n) is defined so that */
/*     the pairs (m,n) map to [1,2,...,indx(l-1,l-1)] with no */
/*     "holes" as m varies from 2 to n and n varies from 2 to l-1. */
/*     (m=0,1 are set from p0n,p1n for all n) */
/*     define for 2.le.n.le.l-1 */
/*     define index function for l.le.n.le.nlat */
/*     preset quantites for fourier transform */
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
/*     lw = 4*nlat*(nlat+1)+2 */
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
	    p0n[np1 + i__ * p0n_dim1] = 0.;
	    p1n[np1 + i__ * p1n_dim1] = 0.;
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
	    abel[imn] = sqrt((doublereal) (((n << 1) + 1) * (m + n - 2) * (m 
		    + n - 3)) / (doublereal) (((n << 1) - 3) * (m + n - 1) * (
		    m + n)));
	    bbel[imn] = sqrt((doublereal) (((n << 1) + 1) * (n - m - 1) * (n 
		    - m)) / (doublereal) (((n << 1) - 3) * (m + n - 1) * (m + 
		    n)));
	    cbel[imn] = sqrt((doublereal) ((n - m + 1) * (n - m + 2)) / (
		    doublereal) ((n + m - 1) * (n + m)));
/* L107: */
	}
    }
    return 0;
} /* shigc1_ */

