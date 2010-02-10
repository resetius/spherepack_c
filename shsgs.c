/* shsgs.f -- translated by f2c (version 20061008).
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



/* ... file shsgs.f */

/*     this file contains code and documentation for subroutines */
/*     shsgs and shsgsi */

/* ... files which must be loaded with shsgs.f */

/*     sphcom.f, hrfft.f, gaqd.f */

/*     subroutine shsgs(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab, */
/*    1                    wshsgs,lshsgs,work,lwork,ierror) */

/*     subroutine shsgs performs the spherical harmonic synthesis */
/*     on the arrays a and b and stores the result in the array g. */
/*     the synthesis is performed on an equally spaced longitude grid */
/*     and a gaussian colatitude grid.  the associated legendre functions */
/*     are stored rather than recomputed as they are in subroutine */
/*     shsgc.  the synthesis is described below at output parameter */
/*     g. */


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

/*     isym   = 0  no symmetries exist about the equator. the synthesis */
/*                 is performed on the entire sphere.  i.e. on the */
/*                 array g(i,j) for i=1,...,nlat and j=1,...,nlon. */
/*                 (see description of g below) */

/*            = 1  g is antisymmetric about the equator. the synthesis */
/*                 is performed on the northern hemisphere only.  i.e. */
/*                 if nlat is odd the synthesis is performed on the */
/*                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even the synthesis is performed on the */
/*                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */


/*            = 2  g is symmetric about the equator. the synthesis is */
/*                 performed on the northern hemisphere only.  i.e. */
/*                 if nlat is odd the synthesis is performed on the */
/*                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even the synthesis is performed on the */
/*                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*     nt     the number of syntheses.  in the program that calls shsgs, */
/*            the arrays g,a and b can be three dimensional in which */
/*            case multiple synthesis will be performed.  the third */
/*            index is the synthesis index which assumes the values */
/*            k=1,...,nt.  for a single synthesis set nt=1. the */
/*            discription of the remaining parameters is simplified */
/*            by assuming that nt=1 or that the arrays g,a and b */
/*            have only two dimensions. */

/*     idg    the first dimension of the array g as it appears in the */
/*            program that calls shagc. if isym equals zero then idg */
/*            must be at least nlat.  if isym is nonzero then idg must */
/*            be at least nlat/2 if nlat is even or at least (nlat+1)/2 */
/*            if nlat is odd. */

/*     jdg    the second dimension of the array g as it appears in the */
/*            program that calls shagc. jdg must be at least nlon. */

/*     a,b    two or three dimensional arrays (see the input parameter */
/*            nt) that contain the coefficients in the spherical harmonic */
/*            expansion of g(i,j) given below at the definition of the */
/*            output parameter g.  a(m,n) and b(m,n) are defined for */
/*            indices m=1,...,mmax and n=m,...,nlat where mmax is the */
/*            maximum (plus one) longitudinal wave number given by */
/*            mmax = min0(nlat,(nlon+2)/2) if nlon is even or */
/*            mmax = min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     mdab   the first dimension of the arrays a and b as it appears */
/*            in the program that calls shsgs. mdab must be at least */
/*            min0((nlon+2)/2,nlat) if nlon is even or at least */
/*            min0((nlon+1)/2,nlat) if nlon is odd. */

/*     ndab   the second dimension of the arrays a and b as it appears */
/*            in the program that calls shsgs. ndab must be at least nlat */

/*     wshsgs an array which must be initialized by subroutine shsgsi. */
/*            once initialized, wshsgs can be used repeatedly by shsgs */
/*            as long as nlat and nlon remain unchanged.  wshsgs must */
/*            not be altered between calls of shsgs. */

/*     lshsgs the dimension of the array wshsgs as it appears in the */
/*            program that calls shsgs. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshsgs must be at least */

/*            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15 */


/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls shsgs. define */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */


/*            if isym is zero then lwork must be at least */

/*                  nlat*nlon*(nt+1) */

/*            if isym is nonzero then lwork must be at least */

/*                  l2*nlon*(nt+1) */


/*     ************************************************************** */

/*     output parameters */

/*     g      a two or three dimensional array (see input parameter nt) */
/*            that contains the discrete function which is synthesized. */
/*            g(i,j) contains the value of the function at the gaussian */
/*            colatitude point theta(i) and longitude point */
/*            phi(j) = (j-1)*2*pi/nlon. the index ranges are defined */
/*            above at the input parameter isym.  for isym=0, g(i,j) */
/*            is given by the the equations listed below.  symmetric */
/*            versions are used when isym is greater than zero. */

/*     the normalized associated legendre functions are given by */

/*     pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m))) */
/*                       *sin(theta)**m/(2**n*factorial(n)) times the */
/*                       (n+m)th derivative of (x**2-1)**n with respect */
/*                       to x=cos(theta) */

/*     define the maximum (plus one) longitudinal wave number */
/*     as   mmax = min0(nlat,(nlon+2)/2) if nlon is even or */
/*          mmax = min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     then g(i,j) = the sum from n=0 to n=nlat-1 of */

/*                   .5*pbar(0,n,theta(i))*a(1,n+1) */

/*              plus the sum from m=1 to m=mmax-1 of */

/*                   the sum from n=m to n=nlat-1 of */

/*              pbar(m,n,theta(i))*(a(m+1,n+1)*cos(m*phi(j)) */
/*                                    -b(m+1,n+1)*sin(m*phi(j))) */


/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of isym */
/*            = 4  error in the specification of nt */
/*            = 5  error in the specification of idg */
/*            = 6  error in the specification of jdg */
/*            = 7  error in the specification of mdab */
/*            = 8  error in the specification of ndab */
/*            = 9  error in the specification of lshsgs */
/*            = 10 error in the specification of lwork */


/* **************************************************************** */

/*     subroutine shsgsi(nlat,nlon,wshsgs,lshsgs,work,lwork,dwork,ldwork, */
/*    +                  ierror) */

/*     subroutine shsgsi initializes the array wshsgs which can then */
/*     be used repeatedly by subroutines shsgs. it precomputes */
/*     and stores in wshsgs quantities such as gaussian weights, */
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

/*     wshsgs an array which must be initialized by subroutine shsgsi. */
/*            once initialized, wshsgs can be used repeatedly by shsgs */
/*            as long as nlat and nlon remain unchanged.  wshsgs must */
/*            not be altered between calls of shsgs. */

/*     lshsgs the dimension of the array wshsgs as it appears in the */
/*            program that calls shsgs. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshsgs must be at least */

/*            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15 */

/*     work   a doublereal work space which need not be saved */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls shsgsi. lwork must be at least */
/*            4*nlat*(nlat+2)+2 in the routine calling shsgsi */

/*     dwork   a double precision work array that does not have to be saved. */

/*     ldwork  the length of dwork in the calling routine.  ldwork must */
/*             be at least nlat*(nlat+4) */

/*     output parameter */

/*     wshsgs an array which must be initialized before calling shsgs or */
/*            once initialized, wshsgs can be used repeatedly by shsgs or */
/*            as long as nlat and nlon remain unchanged.  wshsgs must not */
/*            altered between calls of shsgs. */

/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of lshsgs */
/*            = 4  error in the specification of lwork */
/*            = 5  error in the specification of ldwork */
/*            = 5  failure in gaqd to compute gaussian points */
/*                 (due to failure in eigenvalue routine) */


/* Subroutine */ int shsgs_(integer *nlat, integer *nlon, integer *mode, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a, doublereal *b, 
	integer *mdab, integer *ndab, doublereal *wshsgs, integer *lshsgs, doublereal *
	work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer g_dim1, g_dim2, g_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, i__1;

    /* Local variables */
    static integer l, l1, l2, lp, iw, lat, late, ifft, ipmn;
    extern /* Subroutine */ int shsgs1_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *);

/*     check input parameters */
    /* Parameter adjustments */
    g_dim1 = *idg;
    g_dim2 = *jdg;
    g_offset = 1 + g_dim1 * (1 + g_dim2);
    g -= g_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    --wshsgs;
    --work;

    /* Function Body */
    *ierror = 1;
    if (*nlat < 3) {
	return 0;
    }
    *ierror = 2;
    if (*nlon < 4) {
	return 0;
    }
    *ierror = 3;
    if (*mode < 0 || *mode > 2) {
	return 0;
    }
    *ierror = 4;
    if (*nt < 1) {
	return 0;
    }
/*     set limit on m subscript */
/* Computing MIN */
    i__1 = (*nlon + 2) / 2;
    l = min(i__1,*nlat);
/*     set gaussian point nearest equator pointer */
    late = (*nlat + *nlat % 2) / 2;
/*     set number of grid points for analysis/synthesis */
    lat = *nlat;
    if (*mode != 0) {
	lat = late;
    }
    *ierror = 5;
    if (*idg < lat) {
	return 0;
    }
    *ierror = 6;
    if (*jdg < *nlon) {
	return 0;
    }
    *ierror = 7;
    if (*mdab < l) {
	return 0;
    }
    *ierror = 8;
    if (*ndab < *nlat) {
	return 0;
    }
    l1 = l;
    l2 = late;
    *ierror = 9;
/*     check permanent work space length */
    lp = *nlat * ((l1 + l2) * 3 - 2) + (l1 - 1) * (l2 * ((*nlat << 1) - l1) - 
	    l1 * 3) / 2 + *nlon + 15;
    if (*lshsgs < lp) {
	return 0;
    }
/*     check temporary work space length */
    *ierror = 10;
    if (*mode == 0 && *lwork < *nlat * *nlon * (*nt + 1)) {
	return 0;
    }
    if (*mode != 0 && *lwork < l2 * *nlon * (*nt + 1)) {
	return 0;
    }
    *ierror = 0;
/*     starting address for fft values and legendre polys in wshsgs */
    ifft = *nlat + (*nlat << 1) * late + (l * (l - 1) / 2 + (*nlat - l) * (l 
	    - 1)) * 3 + 1;
    ipmn = ifft + *nlon + 15;
/*     set pointer for internal storage of g */
    iw = lat * *nlon * *nt + 1;
    shsgs1_(nlat, nlon, &l, &lat, mode, &g[g_offset], idg, jdg, nt, &a[
	    a_offset], &b[b_offset], mdab, ndab, &wshsgs[ifft], &wshsgs[ipmn],
	     &late, &work[1], &work[iw]);
    return 0;
} /* shsgs_ */

/* Subroutine */ int shsgs1_(integer *nlat, integer *nlon, integer *l, 
	integer *lat, integer *mode, doublereal *gs, integer *idg, integer *jdg, 
	integer *nt, doublereal *a, doublereal *b, integer *mdab, integer *ndab, doublereal *
	wfft, doublereal *pmn, integer *late, doublereal *g, doublereal *work)
{
    /* System generated locals */
    integer gs_dim1, gs_dim2, gs_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, pmn_dim1, pmn_offset, g_dim1, g_dim2, g_offset, 
	    i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal t1, t2, t3, t4;
    static integer mn, is, ms, ns, lm1, nl2, lp1, mp1, np1, mp2, meo, mml1;
    extern /* Subroutine */ int hrfftb_(integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *);

/*     reconstruct fourier coefficients in g on gaussian grid */
/*     using coefficients in a,b */
/*     initialize to zero */
    /* Parameter adjustments */
    g_dim1 = *lat;
    g_dim2 = *nlon;
    g_offset = 1 + g_dim1 * (1 + g_dim2);
    g -= g_offset;
    gs_dim1 = *idg;
    gs_dim2 = *jdg;
    gs_offset = 1 + gs_dim1 * (1 + gs_dim2);
    gs -= gs_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    --wfft;
    pmn_dim1 = *late;
    pmn_offset = 1 + pmn_dim1;
    pmn -= pmn_offset;
    --work;

    /* Function Body */
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *lat;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		g[i__ + (j + k * g_dim2) * g_dim1] = 0.;
/* L100: */
	    }
	}
    }
    lm1 = *l;
    if (*nlon == *l + *l - 2) {
	lm1 = *l - 1;
    }
    if (*mode == 0) {
/*     set first column in g */
	m = 0;
	mml1 = m * ((*nlat << 1) - m - 1) / 2;
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
/*     n even */
	    i__2 = *nlat;
	    for (np1 = 1; np1 <= i__2; np1 += 2) {
		mn = mml1 + np1;
		i__1 = *late;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    g[i__ + (k * g_dim2 + 1) * g_dim1] += a[(np1 + k * a_dim2)
			     * a_dim1 + 1] * pmn[i__ + mn * pmn_dim1];
/* L102: */
		}
	    }
/*     n odd */
	    nl2 = *nlat / 2;
	    i__1 = *nlat;
	    for (np1 = 2; np1 <= i__1; np1 += 2) {
		mn = mml1 + np1;
		i__2 = nl2;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    is = *nlat - i__ + 1;
		    g[is + (k * g_dim2 + 1) * g_dim1] += a[(np1 + k * a_dim2) 
			    * a_dim1 + 1] * pmn[i__ + mn * pmn_dim1];
/* L103: */
		}
	    }
/* L101: */
	}
/*     restore m=0 coefficients from odd/even */
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__2 = nl2;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		is = *nlat - i__ + 1;
		t1 = g[i__ + (k * g_dim2 + 1) * g_dim1];
		t3 = g[is + (k * g_dim2 + 1) * g_dim1];
		g[i__ + (k * g_dim2 + 1) * g_dim1] = t1 + t3;
		g[is + (k * g_dim2 + 1) * g_dim1] = t1 - t3;
/* L112: */
	    }
	}
/*     sweep interior columns of g */
	i__2 = lm1;
	for (mp1 = 2; mp1 <= i__2; ++mp1) {
	    m = mp1 - 1;
	    mml1 = m * ((*nlat << 1) - m - 1) / 2;
	    mp2 = m + 2;
	    i__3 = *nt;
	    for (k = 1; k <= i__3; ++k) {
/*     for n-m even store (g(i,p,k)+g(nlat-i+1,p,k))/2 in g(i,p,k) p=2*m,2*m+1 */
/*     for i=1,...,late */
		i__1 = *nlat;
		for (np1 = mp1; np1 <= i__1; np1 += 2) {
		    mn = mml1 + np1;
		    i__4 = *late;
		    for (i__ = 1; i__ <= i__4; ++i__) {
			g[i__ + ((m << 1) + k * g_dim2) * g_dim1] += a[mp1 + (
				np1 + k * a_dim2) * a_dim1] * pmn[i__ + mn * 
				pmn_dim1];
			g[i__ + ((m << 1) + 1 + k * g_dim2) * g_dim1] += b[
				mp1 + (np1 + k * b_dim2) * b_dim1] * pmn[i__ 
				+ mn * pmn_dim1];
/* L107: */
		    }
/* L106: */
		}
/*     for n-m odd store g(i,p,k)-g(nlat-i+1,p,k) in g(nlat-i+1,p,k) */
/*     for i=1,...,nlat/2 (p=2*m,p=2*m+1) */
		i__1 = *nlat;
		for (np1 = mp2; np1 <= i__1; np1 += 2) {
		    mn = mml1 + np1;
		    i__4 = nl2;
		    for (i__ = 1; i__ <= i__4; ++i__) {
			is = *nlat - i__ + 1;
			g[is + ((m << 1) + k * g_dim2) * g_dim1] += a[mp1 + (
				np1 + k * a_dim2) * a_dim1] * pmn[i__ + mn * 
				pmn_dim1];
			g[is + ((m << 1) + 1 + k * g_dim2) * g_dim1] += b[mp1 
				+ (np1 + k * b_dim2) * b_dim1] * pmn[i__ + mn 
				* pmn_dim1];
/* L109: */
		    }
/* L108: */
		}
/*     now set fourier coefficients using even-odd reduction above */
		i__1 = nl2;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    is = *nlat - i__ + 1;
		    t1 = g[i__ + ((m << 1) + k * g_dim2) * g_dim1];
		    t2 = g[i__ + ((m << 1) + 1 + k * g_dim2) * g_dim1];
		    t3 = g[is + ((m << 1) + k * g_dim2) * g_dim1];
		    t4 = g[is + ((m << 1) + 1 + k * g_dim2) * g_dim1];
		    g[i__ + ((m << 1) + k * g_dim2) * g_dim1] = t1 + t3;
		    g[i__ + ((m << 1) + 1 + k * g_dim2) * g_dim1] = t2 + t4;
		    g[is + ((m << 1) + k * g_dim2) * g_dim1] = t1 - t3;
		    g[is + ((m << 1) + 1 + k * g_dim2) * g_dim1] = t2 - t4;
/* L110: */
		}
/* L105: */
	    }
/* L104: */
	}
/*     set last column (using a only) if necessary */
	if (*nlon == *l + *l - 2) {
	    m = *l - 1;
	    mml1 = m * ((*nlat << 1) - m - 1) / 2;
	    i__2 = *nt;
	    for (k = 1; k <= i__2; ++k) {
/*     n-m even */
		i__3 = *nlat;
		for (np1 = *l; np1 <= i__3; np1 += 2) {
		    mn = mml1 + np1;
		    i__1 = *late;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			g[i__ + (*nlon + k * g_dim2) * g_dim1] += a[*l + (np1 
				+ k * a_dim2) * a_dim1] * 2. * pmn[i__ + mn *
				 pmn_dim1];
/* L131: */
		    }
		}
		lp1 = *l + 1;
/*     n-m odd */
		i__1 = *nlat;
		for (np1 = lp1; np1 <= i__1; np1 += 2) {
		    mn = mml1 + np1;
		    i__3 = nl2;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			is = *nlat - i__ + 1;
			g[is + (*nlon + k * g_dim2) * g_dim1] += a[*l + (np1 
				+ k * a_dim2) * a_dim1] * 2. * pmn[i__ + mn *
				 pmn_dim1];
/* L132: */
		    }
		}
		i__3 = nl2;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    is = *nlat - i__ + 1;
		    t1 = g[i__ + (*nlon + k * g_dim2) * g_dim1];
		    t3 = g[is + (*nlon + k * g_dim2) * g_dim1];
		    g[i__ + (*nlon + k * g_dim2) * g_dim1] = t1 + t3;
		    g[is + (*nlon + k * g_dim2) * g_dim1] = t1 - t3;
/* L133: */
		}
/* L111: */
	    }
	}
    } else {
/*     half sphere (mode.ne.0) */
/*     set first column in g */
	m = 0;
	mml1 = m * ((*nlat << 1) - m - 1) / 2;
	meo = 1;
	if (*mode == 1) {
	    meo = 2;
	}
	ms = m + meo;
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *nlat;
	    for (np1 = ms; np1 <= i__3; np1 += 2) {
		mn = mml1 + np1;
		i__1 = *late;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    g[i__ + (k * g_dim2 + 1) * g_dim1] += a[(np1 + k * a_dim2)
			     * a_dim1 + 1] * pmn[i__ + mn * pmn_dim1];
/* L113: */
		}
	    }
	}
/*     sweep interior columns of g */
	i__1 = lm1;
	for (mp1 = 2; mp1 <= i__1; ++mp1) {
	    m = mp1 - 1;
	    mml1 = m * ((*nlat << 1) - m - 1) / 2;
	    ms = m + meo;
	    i__3 = *nt;
	    for (k = 1; k <= i__3; ++k) {
		i__2 = *nlat;
		for (np1 = ms; np1 <= i__2; np1 += 2) {
		    mn = mml1 + np1;
		    i__4 = *late;
		    for (i__ = 1; i__ <= i__4; ++i__) {
			g[i__ + ((m << 1) + k * g_dim2) * g_dim1] += a[mp1 + (
				np1 + k * a_dim2) * a_dim1] * pmn[i__ + mn * 
				pmn_dim1];
			g[i__ + ((m << 1) + 1 + k * g_dim2) * g_dim1] += b[
				mp1 + (np1 + k * b_dim2) * b_dim1] * pmn[i__ 
				+ mn * pmn_dim1];
/* L115: */
		    }
		}
	    }
/* L114: */
	}
	if (*nlon == *l + *l - 2) {
/*     set last column */
	    m = *l - 1;
	    mml1 = m * ((*nlat << 1) - m - 1) / 2;
	    ns = *l;
	    if (*mode == 1) {
		ns = *l + 1;
	    }
	    i__1 = *nt;
	    for (k = 1; k <= i__1; ++k) {
		i__4 = *nlat;
		for (np1 = ns; np1 <= i__4; np1 += 2) {
		    mn = mml1 + np1;
		    i__2 = *late;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			g[i__ + (*nlon + k * g_dim2) * g_dim1] += a[*l + (np1 
				+ k * a_dim2) * a_dim1] * 2. * pmn[i__ + mn *
				 pmn_dim1];
/* L116: */
		    }
		}
	    }
	}
    }
/*     do inverse fourier transform */
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	hrfftb_(lat, nlon, &g[(k * g_dim2 + 1) * g_dim1 + 1], lat, &wfft[1], &
		work[1]);
/* L120: */
    }
/*     scale output in gs */
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__4 = *nlon;
	for (j = 1; j <= i__4; ++j) {
	    i__1 = *lat;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		gs[i__ + (j + k * gs_dim2) * gs_dim1] = g[i__ + (j + k * 
			g_dim2) * g_dim1] * .5;
/* L122: */
	    }
	}
    }
    return 0;
} /* shsgs1_ */

/* Subroutine */ int shsgsi_(integer *nlat, integer *nlon, doublereal *wshsgs, 
	integer *lshsgs, doublereal *work, integer *lwork, doublereal *dwork, 
	integer *ldwork, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer l, l1, l2, lp, ldw, late, ipmnf;
    extern /* Subroutine */ int shsgsp_(integer *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, integer *), shsgss1_(integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *);


/*     this subroutine must be called before calling shags or shsgs with */
/*     fixed nlat,nlon. it precomputes the gaussian weights, points */
/*     and all necessary legendre polys and stores them in wshsgs. */
/*     these quantities must be preserved when calling shsgs */
/*     repeatedly with fixed nlat,nlon. */

    /* Parameter adjustments */
    --wshsgs;
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
    if (*lshsgs < lp) {
	return 0;
    }
    *ierror = 4;
/*     check temporary work space */
    if (*lwork < (*nlat << 2) * (*nlat + 2) + 2) {
	return 0;
    }
    *ierror = 5;
    if (*ldwork < *nlat * (*nlat + 4)) {
	return 0;
    }
    *ierror = 0;
/*     set preliminary quantites needed to compute and store legendre polys */
    ldw = *nlat * (*nlat + 4);
    shsgsp_(nlat, nlon, &wshsgs[1], lshsgs, &dwork[1], ldwork, ierror);
    if (*ierror != 0) {
	return 0;
    }
/*     set legendre poly pointer in wshsgs */
    ipmnf = *nlat + (*nlat << 1) * late + (l * (l - 1) / 2 + (*nlat - l) * (l 
	    - 1)) * 3 + *nlon + 16;
    shsgss1_(nlat, &l, &late, &wshsgs[1], &work[1], &wshsgs[ipmnf]);
    return 0;
} /* shsgsi_ */

/* Subroutine */ int shsgss1_(integer *nlat, integer *l, integer *late, doublereal *
	w, doublereal *pmn, doublereal *pmnf)
{
    /* System generated locals */
    integer pmn_dim1, pmn_dim2, pmn_offset, pmnf_dim1, pmnf_offset, i__1, 
	    i__2, i__3;

    /* Local variables */
    static integer i__, j, k, m, km, mn, mp1, np1, mml1, mode;
    extern /* Subroutine */ int legin_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *);

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
		pmn[i__ + (j + k * pmn_dim2) * pmn_dim1] = 0.;
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
} /* shsgss1_ */

/* Subroutine */ int shsgsp_(integer *nlat, integer *nlon, doublereal *wshsgs, 
	integer *lshsgs, doublereal *dwork, integer *ldwork, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer l, i1, i2, i3, l1, l2, i4, i5, i6, i7, iw, late, idth, 
	    idwts;
    extern /* Subroutine */ int shsgsp1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *);

    /* Parameter adjustments */
    --wshsgs;
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
    if (*lshsgs < *nlat * ((l2 << 1) + l1 * 3 - 2) + l1 * 3 * (1 - l1) / 2 + *
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
    shsgsp1_(nlat, nlon, &l, &late, &wshsgs[i1], &wshsgs[i2], &wshsgs[i3], &
	    wshsgs[i4], &wshsgs[i5], &wshsgs[i6], &wshsgs[i7], &dwork[idth], &
	    dwork[idwts], &dwork[iw], ierror);
    if (*ierror != 0) {
	*ierror = 6;
    }
    return 0;
} /* shsgsp_ */

/* Subroutine */ int shsgsp1_(integer *nlat, integer *nlon, integer *l, 
	integer *late, doublereal *wts, doublereal *p0n, doublereal *p1n, doublereal *abel, doublereal *
	bbel, doublereal *cbel, doublereal *wfft, doublereal *dtheta, doublereal *dwts, 
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
	    doublereal *), hrffti_(integer *, doublereal *);

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
	    abel[imn] = sqrt((doublereal) (((n << 1) + 1) * (m + n - 2) * (m + n - 
		    3)) / (doublereal) (((n << 1) - 3) * (m + n - 1) * (m + n)));
	    bbel[imn] = sqrt((doublereal) (((n << 1) + 1) * (n - m - 1) * (n - m)) /
		     (doublereal) (((n << 1) - 3) * (m + n - 1) * (m + n)));
	    cbel[imn] = sqrt((doublereal) ((n - m + 1) * (n - m + 2)) / (doublereal) ((n 
		    + m - 1) * (n + m)));
/* L107: */
	}
    }
    return 0;
} /* shsgsp1_ */

