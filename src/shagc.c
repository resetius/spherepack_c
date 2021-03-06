/* shagc.f -- translated by f2c (version 20061008).
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



/* ... file shagc.f */

/*     this file contains code and documentation for subroutines */
/*     shagc and shagci */

/* ... files which must be loaded with shagc.f */

/*     sphcom.f, hrfft.f, gaqd.f */


/*     subroutine shagc(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab, */
/*    +                 wshagc,lshagc,work,lwork,ierror) */

/*     subroutine shagc performs the spherical harmonic analysis */
/*     on the array g and stores the result in the arrays a and b. */
/*     the analysis is performed on a gaussian grid in colatitude */
/*     and an equally spaced grid in longitude.  the associated */
/*     legendre functions are recomputed rather than stored as they */
/*     are in subroutine shags.  the analysis is described below */
/*     at output parameters a,b. */

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

/*     isym   = 0  no symmetries exist about the equator. the analysis */
/*                 is performed on the entire sphere.  i.e. on the */
/*                 array g(i,j) for i=1,...,nlat and j=1,...,nlon. */
/*                 (see description of g below) */

/*            = 1  g is antisymmetric about the equator. the analysis */
/*                 is performed on the northern hemisphere only.  i.e. */
/*                 if nlat is odd the analysis is performed on the */
/*                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even the analysis is performed on the */
/*                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */


/*            = 2  g is symmetric about the equator. the analysis is */
/*                 performed on the northern hemisphere only.  i.e. */
/*                 if nlat is odd the analysis is performed on the */
/*                 array g(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even the analysis is performed on the */
/*                 array g(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*     nt     the number of analyses.  in the program that calls shagc, */
/*            the arrays g,a and b can be three dimensional in which */
/*            case multiple analyses will be performed.  the third */
/*            index is the analysis index which assumes the values */
/*            k=1,...,nt.  for a single analysis set nt=1. the */
/*            discription of the remaining parameters is simplified */
/*            by assuming that nt=1 or that the arrays g,a and b */
/*            have only two dimensions. */

/*     g      a two or three dimensional array (see input parameter */
/*            nt) that contains the discrete function to be analyzed. */
/*            g(i,j) contains the value of the function at the gaussian */
/*            point theta(i) and longitude point phi(j) = (j-1)*2*pi/nlon */
/*            the index ranges are defined above at the input parameter */
/*            isym. */

/*     idg    the first dimension of the array g as it appears in the */
/*            program that calls shagc. if isym equals zero then idg */
/*            must be at least nlat.  if isym is nonzero then idg must */
/*            be at least nlat/2 if nlat is even or at least (nlat+1)/2 */
/*            if nlat is odd. */

/*     jdg    the second dimension of the array g as it appears in the */
/*            program that calls shagc. jdg must be at least nlon. */

/*     mdab   the first dimension of the arrays a and b as it appears */
/*            in the program that calls shagc. mdab must be at least */
/*            min0((nlon+2)/2,nlat) if nlon is even or at least */
/*            min0((nlon+1)/2,nlat) if nlon is odd */

/*     ndab   the second dimension of the arrays a and b as it appears */
/*            in the program that calls shaec. ndab must be at least nlat */

/*     wshagc an array which must be initialized by subroutine shagci. */
/*            once initialized, wshagc can be used repeatedly by shagc. */
/*            as long as nlat and nlon remain unchanged.  wshagc must */
/*            not be altered between calls of shagc. */

/*     lshagc the dimension of the array wshagc as it appears in the */
/*            program that calls shagc. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshagc must be at least */

/*                  nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls shagc. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            if isym is zero then lwork must be at least */

/*                      nlat*(nlon*nt+max0(3*l2,nlon)) */

/*            if isym is not zero then lwork must be at least */

/*                      l2*(nlon*nt+max0(3*nlat,nlon)) */

/*     ************************************************************** */

/*     output parameters */

/*     a,b    both a,b are two or three dimensional arrays (see input */
/*            parameter nt) that contain the spherical harmonic */
/*            coefficients in the representation of g(i,j) given in the */
/*            discription of subroutine shagc. for isym=0, a(m,n) and */
/*            b(m,n) are given by the equations listed below. symmetric */
/*            versions are used when isym is greater than zero. */

/*     definitions */

/*     1. the normalized associated legendre functions */

/*     pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m))) */
/*                       *sin(theta)**m/(2**n*factorial(n)) times the */
/*                       (n+m)th derivative of (x**2-1)**n with respect */
/*                       to x=cos(theta). */

/*     2. the fourier transform of g(i,j). */

/*     c(m,i)          = 2/nlon times the sum from j=1 to j=nlon of */
/*                       g(i,j)*cos((m-1)*(j-1)*2*pi/nlon) */
/*                       (the first and last terms in this sum */
/*                       are divided by 2) */

/*     s(m,i)          = 2/nlon times the sum from j=2 to j=nlon of */
/*                       g(i,j)*sin((m-1)*(j-1)*2*pi/nlon) */


/*     3. the gaussian points and weights on the sphere */
/*        (computed by subroutine gaqd). */

/*        theta(1),...,theta(nlat) (gaussian pts in radians) */
/*        wts(1),...,wts(nlat) (corresponding gaussian weights) */

/*     4. the maximum (plus one) longitudinal wave number */

/*            mmax = min0(nlat,(nlon+2)/2) if nlon is even or */
/*            mmax = min0(nlat,(nlon+1)/2) if nlon is odd. */


/*     then for m=0,...,mmax-1 and n=m,...,nlat-1 the arrays a,b */
/*     are given by */

/*     a(m+1,n+1)     =  the sum from i=1 to i=nlat of */
/*                       c(m+1,i)*wts(i)*pbar(m,n,theta(i)) */

/*     b(m+1,n+1)      = the sum from i=1 to nlat of */
/*                       s(m+1,i)*wts(i)*pbar(m,n,theta(i)) */

/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of isym */
/*            = 4  error in the specification of nt */
/*            = 5  error in the specification of idg */
/*            = 6  error in the specification of jdg */
/*            = 7  error in the specification of mdab */
/*            = 8  error in the specification of ndab */
/*            = 9  error in the specification of lshagc */
/*            = 10 error in the specification of lwork */


/* **************************************************************** */

/*     subroutine shagci(nlat,nlon,wshagc,lshagc,dwork,ldwork,ierror) */

/*     subroutine shagci initializes the array wshagc which can then */
/*     be used repeatedly by subroutines shagc. it precomputes */
/*     and stores in wshagc quantities such as gaussian weights, */
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

/*     wshagc an array which must be initialized by subroutine shagci. */
/*            once initialized, wshagc can be used repeatedly by shagc */
/*            as long as nlat and nlon remain unchanged.  wshagc must */
/*            not be altered between calls of shagc. */

/*     lshagc the dimension of the array wshagc as it appears in the */
/*            program that calls shagc. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshagc must be at least */

/*                  nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15 */

/*     dwork   a double precision work array that does not have to be saved. */

/*     ldwork  the dimension of the array dwork as it appears in the */
/*            program that calls shagci. ldwork must be at least */

/*                nlat*(nlat+4) */

/*     output parameter */

/*     wshagc an array which must be initialized before calling shagc or */
/*            once initialized, wshagc can be used repeatedly by shagc or */
/*            as long as nlat and nlon remain unchanged.  wshagc must not */
/*            altered between calls of shagc. */

/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of lshagc */
/*            = 4  error in the specification of ldwork */
/*            = 5  failure in gaqd to compute gaussian points */
/*                 (due to failure in eigenvalue routine) */


/* **************************************************************** */
/* Subroutine */ int shagc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a,
	 doublereal *b, integer *mdab, integer *ndab, doublereal *wshagc, 
	integer *lshagc, doublereal *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer g_dim1, g_dim2, g_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, i__1;

    /* Local variables */
    integer l, l1, l2, lat, late, ifft, ipmn, iwts;
    extern /* Subroutine */ int shagc1_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    );

/*     subroutine shagc performs the spherical harmonic analysis on */
/*     a gaussian grid on the array(s) in g and returns the coefficients */
/*     in array(s) a,b. the necessary legendre polynomials are computed */
/*     as needed in this version. */

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
    --wshagc;
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
    if (*isym < 0 || *isym > 2) {
	return 0;
    }
    *ierror = 4;
    if (*nt < 1) {
	return 0;
    }
/*     set upper limit on m for spherical harmonic basis */
/* Computing MIN */
    i__1 = (*nlon + 2) / 2;
    l = min(i__1,*nlat);
/*     set gaussian point nearest equator pointer */
    late = (*nlat + *nlat % 2) / 2;
/*     set number of grid points for analysis/synthesis */
    lat = *nlat;
    if (*isym != 0) {
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
    if (*lshagc < *nlat * ((l2 << 1) + l1 * 3 - 2) + l1 * 3 * (1 - l1) / 2 + *
	    nlon + 15) {
	return 0;
    }
    *ierror = 10;
/*     check temporary work space length */
    if (*isym == 0) {
/* Computing MAX */
	i__1 = l2 * 3;
	if (*lwork < *nlat * (*nlon * *nt + max(i__1,*nlon))) {
	    return 0;
	}
    } else {
/*     isym.ne.0 */
/* Computing MAX */
	i__1 = *nlat * 3;
	if (*lwork < l2 * (*nlon * *nt + max(i__1,*nlon))) {
	    return 0;
	}
    }
    *ierror = 0;
/*     starting address for gaussian wts in shigc and fft values */
    iwts = 1;
    ifft = *nlat + (*nlat << 1) * late + (l * (l - 1) / 2 + (*nlat - l) * (l 
	    - 1)) * 3 + 1;
/*     set pointers for internal storage of g and legendre polys */
    ipmn = lat * *nlon * *nt + 1;
    shagc1_(nlat, nlon, &l, &lat, isym, &g[g_offset], idg, jdg, nt, &a[
	    a_offset], &b[b_offset], mdab, ndab, &wshagc[1], &wshagc[iwts], &
	    wshagc[ifft], &late, &work[ipmn], &work[1]);
    return 0;
} /* shagc_ */

/* Subroutine */ int shagc1_(integer *nlat, integer *nlon, integer *l, 
	integer *lat, integer *mode, doublereal *gs, integer *idg, integer *
	jdg, integer *nt, doublereal *a, doublereal *b, integer *mdab, 
	integer *ndab, doublereal *w, doublereal *wts, doublereal *wfft, 
	integer *late, doublereal *pmn, doublereal *g)
{
    /* System generated locals */
    integer gs_dim1, gs_dim2, gs_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, g_dim1, g_dim2, g_offset, pmn_dim1, pmn_dim2, 
	    pmn_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer i__, j, k, m;
    doublereal t1, t2;
    integer km, is, ms, ns, lm1, nl2, lp1, mp1, np1, mp2;
    doublereal sfn;
    extern /* Subroutine */ int legin_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *), hrfftf_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *);

/*     set gs array internally in shagc1 */
    /* Parameter adjustments */
    --wts;
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
    --w;
    --wfft;
    pmn_dim1 = *nlat;
    pmn_dim2 = *late;
    pmn_offset = 1 + pmn_dim1 * (1 + pmn_dim2);
    pmn -= pmn_offset;

    /* Function Body */
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *lat;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		g[i__ + (j + k * g_dim2) * g_dim1] = gs[i__ + (j + k * 
			gs_dim2) * gs_dim1];
/* L100: */
	    }
	}
    }
/*     do fourier transform */
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	hrfftf_(lat, nlon, &g[(k * g_dim2 + 1) * g_dim1 + 1], lat, &wfft[1], &
		pmn[pmn_offset]);
/* L101: */
    }
/*     scale result */
    sfn = 2. / (doublereal) (*nlon);
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    i__1 = *lat;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		g[i__ + (j + k * g_dim2) * g_dim1] = sfn * g[i__ + (j + k * 
			g_dim2) * g_dim1];
/* L102: */
	    }
	}
    }
/*     compute using gaussian quadrature */
/*     a(n,m) = s (ga(theta,m)*pnm(theta)*sin(theta)*dtheta) */
/*     b(n,m) = s (gb(theta,m)*pnm(theta)*sin(theta)*dtheta) */
/*     here ga,gb are the cos(phi),sin(phi) coefficients of */
/*     the fourier expansion of g(theta,phi) in phi.  as a result */
/*     of the above fourier transform they are stored in array */
/*     g as follows: */
/*     for each theta(i) and k= l-1 */
/*     ga(0),ga(1),gb(1),ga(2),gb(2),...,ga(k-1),gb(k-1),ga(k) */
/*     correspond to (in the case nlon=l+l-2) */
/*     g(i,1),g(i,2),g(i,3),g(i,4),g(i,5),...,g(i,2l-4),g(i,2l-3),g(i,2l- */
/*     initialize coefficients to zero */
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlat;
	for (np1 = 1; np1 <= i__2; ++np1) {
	    i__3 = *l;
	    for (mp1 = 1; mp1 <= i__3; ++mp1) {
		a[mp1 + (np1 + k * a_dim2) * a_dim1] = 0.;
		b[mp1 + (np1 + k * b_dim2) * b_dim1] = 0.;
/* L103: */
	    }
	}
    }
/*     set m+1 limit on b(m+1) calculation */
    lm1 = *l;
    if (*nlon == *l + *l - 2) {
	lm1 = *l - 1;
    }
    if (*mode == 0) {
/*     for full sphere (mode=0) and even/odd reduction: */
/*     overwrite g(i) with (g(i)+g(nlat-i+1))*wts(i) */
/*     overwrite g(nlat-i+1) with (g(i)-g(nlat-i+1))*wts(i) */
	nl2 = *nlat / 2;
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		i__1 = nl2;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    is = *nlat - i__ + 1;
		    t1 = g[i__ + (j + k * g_dim2) * g_dim1];
		    t2 = g[is + (j + k * g_dim2) * g_dim1];
		    g[i__ + (j + k * g_dim2) * g_dim1] = wts[i__] * (t1 + t2);
		    g[is + (j + k * g_dim2) * g_dim1] = wts[i__] * (t1 - t2);
/* L105: */
		}
/*     adjust equator if necessary(nlat odd) */
		if (*nlat % 2 != 0) {
		    g[*late + (j + k * g_dim2) * g_dim1] = wts[*late] * g[*
			    late + (j + k * g_dim2) * g_dim1];
		}
/* L104: */
	    }
	}
/*     set m = 0 coefficients first */
	m = 0;
	legin_(mode, l, nlat, &m, &w[1], &pmn[pmn_offset], &km);
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *late;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		is = *nlat - i__ + 1;
		i__1 = *nlat;
		for (np1 = 1; np1 <= i__1; np1 += 2) {
/*     n even */
		    a[(np1 + k * a_dim2) * a_dim1 + 1] += g[i__ + (k * g_dim2 
			    + 1) * g_dim1] * pmn[np1 + (i__ + km * pmn_dim2) *
			     pmn_dim1];
/* L107: */
		}
		i__1 = *nlat;
		for (np1 = 2; np1 <= i__1; np1 += 2) {
/*     n odd */
		    a[(np1 + k * a_dim2) * a_dim1 + 1] += g[is + (k * g_dim2 
			    + 1) * g_dim1] * pmn[np1 + (i__ + km * pmn_dim2) *
			     pmn_dim1];
/* L108: */
		}
/* L106: */
	    }
	}
/*     compute coefficients for which b(m,n) is available */
	i__3 = lm1;
	for (mp1 = 2; mp1 <= i__3; ++mp1) {
	    m = mp1 - 1;
	    mp2 = m + 2;
/*     compute pmn for all i and n=m,...,l-1 */
	    legin_(mode, l, nlat, &m, &w[1], &pmn[pmn_offset], &km);
	    i__2 = *nt;
	    for (k = 1; k <= i__2; ++k) {
		i__1 = *late;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    is = *nlat - i__ + 1;
/*     n-m even */
		    i__4 = *nlat;
		    for (np1 = mp1; np1 <= i__4; np1 += 2) {
			a[mp1 + (np1 + k * a_dim2) * a_dim1] += g[i__ + ((m <<
				 1) + k * g_dim2) * g_dim1] * pmn[np1 + (i__ 
				+ km * pmn_dim2) * pmn_dim1];
			b[mp1 + (np1 + k * b_dim2) * b_dim1] += g[i__ + ((m <<
				 1) + 1 + k * g_dim2) * g_dim1] * pmn[np1 + (
				i__ + km * pmn_dim2) * pmn_dim1];
/* L112: */
		    }
/*     n-m odd */
		    i__4 = *nlat;
		    for (np1 = mp2; np1 <= i__4; np1 += 2) {
			a[mp1 + (np1 + k * a_dim2) * a_dim1] += g[is + ((m << 
				1) + k * g_dim2) * g_dim1] * pmn[np1 + (i__ + 
				km * pmn_dim2) * pmn_dim1];
			b[mp1 + (np1 + k * b_dim2) * b_dim1] += g[is + ((m << 
				1) + 1 + k * g_dim2) * g_dim1] * pmn[np1 + (
				i__ + km * pmn_dim2) * pmn_dim1];
/* L113: */
		    }
/* L111: */
		}
/* L110: */
	    }
/* L109: */
	}
	if (*nlon == *l + *l - 2) {
/*     compute a(l,np1) coefficients only */
	    m = *l - 1;
	    legin_(mode, l, nlat, &m, &w[1], &pmn[pmn_offset], &km);
	    i__3 = *nt;
	    for (k = 1; k <= i__3; ++k) {
		i__2 = *late;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    is = *nlat - i__ + 1;
/*     n-m even */
		    i__1 = *nlat;
		    for (np1 = *l; np1 <= i__1; np1 += 2) {
			a[*l + (np1 + k * a_dim2) * a_dim1] += g[i__ + (*nlon 
				+ k * g_dim2) * g_dim1] * .5 * pmn[np1 + (i__ 
				+ km * pmn_dim2) * pmn_dim1];
/* L124: */
		    }
		    lp1 = *l + 1;
/*     n-m odd */
		    i__1 = *nlat;
		    for (np1 = lp1; np1 <= i__1; np1 += 2) {
			a[*l + (np1 + k * a_dim2) * a_dim1] += g[is + (*nlon 
				+ k * g_dim2) * g_dim1] * .5 * pmn[np1 + (i__ 
				+ km * pmn_dim2) * pmn_dim1];
/* L125: */
		    }
/* L114: */
		}
	    }
	}
    } else {
/*     half sphere */
/*     overwrite g(i) with wts(i)*(g(i)+g(i)) for i=1,...,nlate/2 */
	nl2 = *nlat / 2;
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *nlon;
	    for (j = 1; j <= i__3; ++j) {
		i__1 = nl2;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    g[i__ + (j + k * g_dim2) * g_dim1] = wts[i__] * (g[i__ + (
			    j + k * g_dim2) * g_dim1] + g[i__ + (j + k * 
			    g_dim2) * g_dim1]);
/* L115: */
		}
/*     adjust equator separately if a grid point */
		if (nl2 < *late) {
		    g[*late + (j + k * g_dim2) * g_dim1] = wts[*late] * g[*
			    late + (j + k * g_dim2) * g_dim1];
		}
/* L116: */
	    }
	}
/*     set m = 0 coefficients first */
	m = 0;
	legin_(mode, l, nlat, &m, &w[1], &pmn[pmn_offset], &km);
	ms = 1;
	if (*mode == 1) {
	    ms = 2;
	}
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__2 = *late;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = *nlat;
		for (np1 = ms; np1 <= i__1; np1 += 2) {
		    a[(np1 + k * a_dim2) * a_dim1 + 1] += g[i__ + (k * g_dim2 
			    + 1) * g_dim1] * pmn[np1 + (i__ + km * pmn_dim2) *
			     pmn_dim1];
/* L117: */
		}
	    }
	}
/*     compute coefficients for which b(m,n) is available */
	i__1 = lm1;
	for (mp1 = 2; mp1 <= i__1; ++mp1) {
	    m = mp1 - 1;
	    ms = mp1;
	    if (*mode == 1) {
		ms = mp1 + 1;
	    }
/*     compute pmn for all i and n=m,...,nlat-1 */
	    legin_(mode, l, nlat, &m, &w[1], &pmn[pmn_offset], &km);
	    i__2 = *nt;
	    for (k = 1; k <= i__2; ++k) {
		i__3 = *late;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    i__4 = *nlat;
		    for (np1 = ms; np1 <= i__4; np1 += 2) {
			a[mp1 + (np1 + k * a_dim2) * a_dim1] += g[i__ + ((m <<
				 1) + k * g_dim2) * g_dim1] * pmn[np1 + (i__ 
				+ km * pmn_dim2) * pmn_dim1];
			b[mp1 + (np1 + k * b_dim2) * b_dim1] += g[i__ + ((m <<
				 1) + 1 + k * g_dim2) * g_dim1] * pmn[np1 + (
				i__ + km * pmn_dim2) * pmn_dim1];
/* L119: */
		    }
		}
	    }
/* L118: */
	}
	if (*nlon == *l + *l - 2) {
/*     compute coefficient a(l,np1) only */
	    m = *l - 1;
	    legin_(mode, l, nlat, &m, &w[1], &pmn[pmn_offset], &km);
	    ns = *l;
	    if (*mode == 1) {
		ns = *l + 1;
	    }
	    i__1 = *nt;
	    for (k = 1; k <= i__1; ++k) {
		i__4 = *late;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    i__3 = *nlat;
		    for (np1 = ns; np1 <= i__3; np1 += 2) {
			a[*l + (np1 + k * a_dim2) * a_dim1] += g[i__ + (*nlon 
				+ k * g_dim2) * g_dim1] * .5 * pmn[np1 + (i__ 
				+ km * pmn_dim2) * pmn_dim1];
/* L120: */
		    }
		}
	    }
	}
    }
    return 0;
} /* shagc1_ */

/* Subroutine */ int shagci_(integer *nlat, integer *nlon, doublereal *wshagc,
	 integer *lshagc, doublereal *dwork, integer *ldwork, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer l, i1, i2, i3, l1, l2, i4, i5, i6, i7, iw, late, idth, idwts;
    extern /* Subroutine */ int shagci1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);

/*     this subroutine must be called before calling shagc with */
/*     fixed nlat,nlon. it precomputes quantites such as the gaussian */
/*     points and weights, m=0,m=1 legendre polynomials, recursion */
/*     recursion coefficients. */
    /* Parameter adjustments */
    --wshagc;
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
    if (*lshagc < *nlat * ((l2 << 1) + l1 * 3 - 2) + l1 * 3 * (1 - l1) / 2 + *
	    nlon + 15) {
	return 0;
    }
    *ierror = 4;
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
    idwts = idth + *nlat;
    iw = idwts + *nlat;
    shagci1_(nlat, nlon, &l, &late, &wshagc[i1], &wshagc[i2], &wshagc[i3], &
	    wshagc[i4], &wshagc[i5], &wshagc[i6], &wshagc[i7], &dwork[idth], &
	    dwork[idwts], &dwork[iw], ierror);
    if (*ierror != 0) {
	*ierror = 5;
    }
    return 0;
} /* shagci_ */

/* Subroutine */ int shagci1_(integer *nlat, integer *nlon, integer *l, 
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
} /* shagci1_ */

