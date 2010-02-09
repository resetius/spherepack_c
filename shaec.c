/* shaec.f -- translated by f2c (version 20061008).
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

static integer c__2 = 2;
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


/* ... file shaec.f */

/*     this file contains code and documentation for subroutines */
/*     shaec and shaeci */

/* ... files which must be loaded with shaec.f */

/*     sphcom.f, hrfft.f */

/*     subroutine shaec(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab, */
/*    +                 wshaec,lshaec,work,lwork,ierror) */

/*     subroutine shaec performs the spherical harmonic analysis */
/*     on the array g and stores the result in the arrays a and b. */
/*     the analysis is performed on an equally spaced grid.  the */
/*     associated legendre functions are recomputed rather than stored */
/*     as they are in subroutine shaes.  the analysis is described */
/*     below at output parameters a,b. */

/*     input parameters */

/*     nlat   the number of colatitudes on the full sphere including the */
/*            poles. for example, nlat = 37 for a five degree grid. */
/*            nlat determines the grid increment in colatitude as */
/*            pi/(nlat-1).  if nlat is odd the equator is located at */
/*            grid point i=(nlat+1)/2. if nlat is even the equator is */
/*            located half way between points i=nlat/2 and i=nlat/2+1. */
/*            nlat must be at least 3. note: on the half sphere, the */
/*            number of grid points in the colatitudinal direction is */
/*            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd. */

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

/*     nt     the number of analyses.  in the program that calls shaec, */
/*            the arrays g,a and b can be three dimensional in which */
/*            case multiple analyses will be performed.  the third */
/*            index is the analysis index which assumes the values */
/*            k=1,...,nt.  for a single analysis set nt=1. the */
/*            discription of the remaining parameters is simplified */
/*            by assuming that nt=1 or that the arrays g,a and b */
/*            have only two dimensions. */

/*     g      a two or three dimensional array (see input parameter */
/*            nt) that contains the discrete function to be analyzed. */
/*            g(i,j) contains the value of the function at the colatitude */
/*            point theta(i) = (i-1)*pi/(nlat-1) and longitude point */
/*            phi(j) = (j-1)*2*pi/nlon. the index ranges are defined */
/*            above at the input parameter isym. */


/*     idg    the first dimension of the array g as it appears in the */
/*            program that calls shaec.  if isym equals zero then idg */
/*            must be at least nlat.  if isym is nonzero then idg */
/*            must be at least nlat/2 if nlat is even or at least */
/*            (nlat+1)/2 if nlat is odd. */

/*     jdg    the second dimension of the array g as it appears in the */
/*            program that calls shaec.  jdg must be at least nlon. */

/*     mdab   the first dimension of the arrays a and b as it appears */
/*            in the program that calls shaec. mdab must be at least */
/*            min0(nlat,(nlon+2)/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndab   the second dimension of the arrays a and b as it appears */
/*            in the program that calls shaec. ndab must be at least nlat */

/*     wshaec an array which must be initialized by subroutine shaeci. */
/*            once initialized, wshaec can be used repeatedly by shaec */
/*            as long as nlon and nlat remain unchanged.  wshaec must */
/*            not be altered between calls of shaec. */

/*     lshaec the dimension of the array wshaec as it appears in the */
/*            program that calls shaec. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshaec must be at least */

/*            2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls shaec. define */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            if isym is zero then lwork must be at least */

/*                    nlat*(nt*nlon+max0(3*l2,nlon)) */

/*            if isym is not zero then lwork must be at least */

/*                    l2*(nt*nlon+max0(3*nlat,nlon)) */

/*     ************************************************************** */

/*     output parameters */

/*     a,b    both a,b are two or three dimensional arrays (see input */
/*            parameter nt) that contain the spherical harmonic */
/*            coefficients in the representation of g(i,j) given in the */
/*            discription of subroutine shsec. for isym=0, a(m,n) and */
/*            b(m,n) are given by the equations listed below. symmetric */
/*            versions are used when isym is greater than zero. */



/*     definitions */

/*     1. the normalized associated legendre functions */

/*     pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m))) */
/*                       *sin(theta)**m/(2**n*factorial(n)) times the */
/*                       (n+m)th derivative of (x**2-1)**n with respect */
/*                       to x=cos(theta) */

/*     2. the normalized z functions for m even */

/*     zbar(m,n,theta) = 2/(nlat-1) times the sum from k=0 to k=nlat-1 of */
/*                       the integral from tau = 0 to tau = pi of */
/*                       cos(k*theta)*cos(k*tau)*pbar(m,n,tau)*sin(tau) */
/*                       (first and last terms in this sum are divided */
/*                       by 2) */

/*     3. the normalized z functions for m odd */

/*     zbar(m,n,theta) = 2/(nlat-1) times the sum from k=0 to k=nlat-1 of */
/*                       of the integral from tau = 0 to tau = pi of */
/*                       sin(k*theta)*sin(k*tau)*pbar(m,n,tau)*sin(tau) */

/*     4. the fourier transform of g(i,j). */

/*     c(m,i)          = 2/nlon times the sum from j=1 to j=nlon */
/*                       of g(i,j)*cos((m-1)*(j-1)*2*pi/nlon) */
/*                       (the first and last terms in this sum */
/*                       are divided by 2) */

/*     s(m,i)          = 2/nlon times the sum from j=2 to j=nlon */
/*                       of g(i,j)*sin((m-1)*(j-1)*2*pi/nlon) */

/*     5. the maximum (plus one) longitudinal wave number */

/*            mmax = min0(nlat,(nlon+2)/2) if nlon is even or */
/*            mmax = min0(nlat,(nlon+1)/2) if nlon is odd. */


/*     then for m=0,...,mmax-1 and n=m,...,nlat-1 the arrays a,b */
/*     are given by */

/*     a(m+1,n+1)      = the sum from i=1 to i=nlat of */
/*                       c(m+1,i)*zbar(m,n,theta(i)) */
/*                       (first and last terms in this sum are */
/*                       divided by 2) */

/*     b(m+1,n+1)      = the sum from i=1 to i=nlat of */
/*                       s(m+1,i)*zbar(m,n,theta(i)) */


/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of isym */
/*            = 4  error in the specification of nt */
/*            = 5  error in the specification of idg */
/*            = 6  error in the specification of jdg */
/*            = 7  error in the specification of mdab */
/*            = 8  error in the specification of ndab */
/*            = 9  error in the specification of lshaec */
/*            = 10 error in the specification of lwork */


/* **************************************************************** */
/*     subroutine shaeci(nlat,nlon,wshaec,lshaec,dwork,ldwork,ierror) */

/*     subroutine shaeci initializes the array wshaec which can then */
/*     be used repeatedly by subroutine shaec. */

/*     input parameters */

/*     nlat   the number of colatitudes on the full sphere including the */
/*            poles. for example, nlat = 37 for a five degree grid. */
/*            nlat determines the grid increment in colatitude as */
/*            pi/(nlat-1).  if nlat is odd the equator is located at */
/*            grid point i=(nlat+1)/2. if nlat is even the equator is */
/*            located half way between points i=nlat/2 and i=nlat/2+1. */
/*            nlat must be at least 3. note: on the half sphere, the */
/*            number of grid points in the colatitudinal direction is */
/*            nlat/2 if nlat is even or (nlat+1)/2 if nlat is odd. */

/*     nlon   the number of distinct londitude points.  nlon determines */
/*            the grid increment in longitude as 2*pi/nlon. for example */
/*            nlon = 72 for a five degree grid. nlon must be greater */
/*            than or equal to 4. the efficiency of the computation is */
/*            improved when nlon is a product of small prime numbers. */

/*     lshaec the dimension of the array wshaec as it appears in the */
/*            program that calls shaeci. the array wshaec is an output */
/*            parameter which is described below. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshaec must be at least */

/*            2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15 */

/*     dwork  a double precision dwork array that does not have to be saved. */

/*     ldwork the dimension of the array dwork as it appears in the */
/*            program that calls shaeci.  ldwork  must be at least */
/*            nlat+1. */


/*     output parameters */

/*     wshaec an array which is initialized for use by subroutine shaec. */
/*            once initialized, wshaec can be used repeatedly by shaec */
/*            as long as nlon and nlat remain unchanged.  wshaec must */
/*            not be altered between calls of shaec. */

/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of lshaec */
/*            = 4  error in the specification of ldwork */


/* ******************************************************************* */
/* Subroutine */ int shaec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a, doublereal *b, 
	integer *mdab, integer *ndab, doublereal *wshaec, integer *lshaec, doublereal *
	work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer g_dim1, g_dim2, g_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    static integer ls, iw1, nln, ist, lzz1, labc, imid, mmax;
    extern /* Subroutine */ int shaec1_(integer *, integer *, integer *, doublereal 
	    *, integer *, integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

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
    --wshaec;
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
    if (*nt < 0) {
	return 0;
    }
    *ierror = 5;
    if (*isym == 0 && *idg < *nlat || *isym != 0 && *idg < (*nlat + 1) / 2) {
	return 0;
    }
    *ierror = 6;
    if (*jdg < *nlon) {
	return 0;
    }
    *ierror = 7;
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mmax = min(i__1,i__2);
    if (*mdab < mmax) {
	return 0;
    }
    *ierror = 8;
    if (*ndab < *nlat) {
	return 0;
    }
    *ierror = 9;
    imid = (*nlat + 1) / 2;
    lzz1 = (*nlat << 1) * imid;
    labc = (mmax - 2) * (*nlat + *nlat - mmax - 1) * 3 / 2;
    if (*lshaec < lzz1 + labc + *nlon + 15) {
	return 0;
    }
    *ierror = 10;
    ls = *nlat;
    if (*isym > 0) {
	ls = imid;
    }
    nln = *nt * ls * *nlon;
/* Computing MAX */
    i__1 = ls * *nlon, i__2 = *nlat * 3 * imid;
    if (*lwork < nln + max(i__1,i__2)) {
	return 0;
    }
    *ierror = 0;
    ist = 0;
    if (*isym == 0) {
	ist = imid;
    }
    iw1 = lzz1 + labc + 1;
    shaec1_(nlat, isym, nt, &g[g_offset], idg, jdg, &a[a_offset], &b[b_offset]
	    , mdab, ndab, &imid, &ls, nlon, &work[1], &work[ist + 1], &work[
	    nln + 1], &work[nln + 1], &wshaec[1], &wshaec[iw1]);
    return 0;
} /* shaec_ */

/* Subroutine */ int shaec1_(integer *nlat, integer *isym, integer *nt, doublereal *
	g, integer *idgs, integer *jdgs, doublereal *a, doublereal *b, integer *mdab, 
	integer *ndab, integer *imid, integer *idg, integer *jdg, doublereal *ge, 
	doublereal *go, doublereal *work, doublereal *zb, doublereal *wzfin, doublereal *whrfft)
{
    /* System generated locals */
    integer g_dim1, g_dim2, g_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, ge_dim1, ge_dim2, ge_offset, go_dim1, go_dim2, 
	    go_offset, zb_dim1, zb_dim2, zb_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, m, i3, ls, mp1, np1, mp2, mdo, ndo;
    static doublereal fsn, tsn;
    static integer imm1, nlp1, modl, mmax, nlon;
    extern /* Subroutine */ int zfin_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *), hrfftf_(integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *);


/*     whrfft must have at least nlon+15 locations */
/*     wzfin must have 2*l*(nlat+1)/2 + ((l-3)*l+2)/2 locations */
/*     zb must have 3*l*(nlat+1)/2 locations */
/*     work must have ls*nlon locations */

    /* Parameter adjustments */
    g_dim1 = *idgs;
    g_dim2 = *jdgs;
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
    zb_dim1 = *imid;
    zb_dim2 = *nlat;
    zb_offset = 1 + zb_dim1 * (1 + zb_dim2);
    zb -= zb_offset;
    go_dim1 = *idg;
    go_dim2 = *jdg;
    go_offset = 1 + go_dim1 * (1 + go_dim2);
    go -= go_offset;
    ge_dim1 = *idg;
    ge_dim2 = *jdg;
    ge_offset = 1 + ge_dim1 * (1 + ge_dim2);
    ge -= ge_offset;
    --work;
    --wzfin;
    --whrfft;

    /* Function Body */
    ls = *idg;
    nlon = *jdg;
/* Computing MIN */
    i__1 = *nlat, i__2 = nlon / 2 + 1;
    mmax = min(i__1,i__2);
    mdo = mmax;
    if (mdo + mdo - 1 > nlon) {
	mdo = mmax - 1;
    }
    nlp1 = *nlat + 1;
    tsn = 2.f / nlon;
    fsn = 4.f / nlon;
    modl = *nlat % 2;
    imm1 = *imid;
    if (modl != 0) {
	imm1 = *imid - 1;
    }
    if (*isym != 0) {
	goto L15;
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = imm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = nlon;
	    for (j = 1; j <= i__3; ++j) {
		ge[i__ + (j + k * ge_dim2) * ge_dim1] = tsn * (g[i__ + (j + k 
			* g_dim2) * g_dim1] + g[nlp1 - i__ + (j + k * g_dim2) 
			* g_dim1]);
		go[i__ + (j + k * go_dim2) * go_dim1] = tsn * (g[i__ + (j + k 
			* g_dim2) * g_dim1] - g[nlp1 - i__ + (j + k * g_dim2) 
			* g_dim1]);
/* L5: */
	    }
	}
    }
    goto L30;
L15:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = imm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = nlon;
	    for (j = 1; j <= i__1; ++j) {
		ge[i__ + (j + k * ge_dim2) * ge_dim1] = fsn * g[i__ + (j + k *
			 g_dim2) * g_dim1];
/* L20: */
	    }
	}
    }
    if (*isym == 1) {
	goto L27;
    }
L30:
    if (modl == 0) {
	goto L27;
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = nlon;
	for (j = 1; j <= i__2; ++j) {
	    ge[*imid + (j + k * ge_dim2) * ge_dim1] = tsn * g[*imid + (j + k *
		     g_dim2) * g_dim1];
/* L25: */
	}
    }
L27:
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	hrfftf_(&ls, &nlon, &ge[(k * ge_dim2 + 1) * ge_dim1 + 1], &ls, &
		whrfft[1], &work[1]);
	if (nlon % 2 != 0) {
	    goto L35;
	}
	i__1 = ls;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ge[i__ + (nlon + k * ge_dim2) * ge_dim1] *= .5f;
/* L36: */
	}
L35:
	;
    }
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__1 = mmax;
	for (mp1 = 1; mp1 <= i__1; ++mp1) {
	    i__3 = *nlat;
	    for (np1 = mp1; np1 <= i__3; ++np1) {
		a[mp1 + (np1 + k * a_dim2) * a_dim1] = 0.f;
		b[mp1 + (np1 + k * b_dim2) * b_dim1] = 0.f;
/* L40: */
	    }
	}
    }
    if (*isym == 1) {
	goto L145;
    }
    zfin_(&c__2, nlat, &nlon, &c__0, &zb[zb_offset], &i3, &wzfin[1]);
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlat;
	    for (np1 = 1; np1 <= i__2; np1 += 2) {
		a[(np1 + k * a_dim2) * a_dim1 + 1] += zb[i__ + (np1 + i3 * 
			zb_dim2) * zb_dim1] * ge[i__ + (k * ge_dim2 + 1) * 
			ge_dim1];
/* L110: */
	    }
	}
    }
    ndo = *nlat;
    if (*nlat % 2 == 0) {
	ndo = *nlat - 1;
    }
    i__2 = mdo;
    for (mp1 = 2; mp1 <= i__2; ++mp1) {
	m = mp1 - 1;
	zfin_(&c__2, nlat, &nlon, &m, &zb[zb_offset], &i3, &wzfin[1]);
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = ndo;
		for (np1 = mp1; np1 <= i__4; np1 += 2) {
		    a[mp1 + (np1 + k * a_dim2) * a_dim1] += zb[i__ + (np1 + 
			    i3 * zb_dim2) * zb_dim1] * ge[i__ + ((mp1 << 1) - 
			    2 + k * ge_dim2) * ge_dim1];
		    b[mp1 + (np1 + k * b_dim2) * b_dim1] += zb[i__ + (np1 + 
			    i3 * zb_dim2) * zb_dim1] * ge[i__ + ((mp1 << 1) - 
			    1 + k * ge_dim2) * ge_dim1];
/* L120: */
		}
	    }
	}
    }
    if (mdo == mmax || mmax > ndo) {
	goto L135;
    }
    zfin_(&c__2, nlat, &nlon, &mdo, &zb[zb_offset], &i3, &wzfin[1]);
    i__4 = *nt;
    for (k = 1; k <= i__4; ++k) {
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__1 = ndo;
	    for (np1 = mmax; np1 <= i__1; np1 += 2) {
		a[mmax + (np1 + k * a_dim2) * a_dim1] += zb[i__ + (np1 + i3 * 
			zb_dim2) * zb_dim1] * ge[i__ + ((mmax << 1) - 2 + k * 
			ge_dim2) * ge_dim1];
/* L130: */
	    }
	}
    }
L135:
    if (*isym == 2) {
	return 0;
    }
L145:
    zfin_(&c__1, nlat, &nlon, &c__0, &zb[zb_offset], &i3, &wzfin[1]);
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__3 = imm1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = *nlat;
	    for (np1 = 2; np1 <= i__4; np1 += 2) {
		a[(np1 + k * a_dim2) * a_dim1 + 1] += zb[i__ + (np1 + i3 * 
			zb_dim2) * zb_dim1] * go[i__ + (k * go_dim2 + 1) * 
			go_dim1];
/* L150: */
	    }
	}
    }
    ndo = *nlat;
    if (*nlat % 2 != 0) {
	ndo = *nlat - 1;
    }
    i__4 = mdo;
    for (mp1 = 2; mp1 <= i__4; ++mp1) {
	m = mp1 - 1;
	mp2 = mp1 + 1;
	zfin_(&c__1, nlat, &nlon, &m, &zb[zb_offset], &i3, &wzfin[1]);
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = ndo;
		for (np1 = mp2; np1 <= i__2; np1 += 2) {
		    a[mp1 + (np1 + k * a_dim2) * a_dim1] += zb[i__ + (np1 + 
			    i3 * zb_dim2) * zb_dim1] * go[i__ + ((mp1 << 1) - 
			    2 + k * go_dim2) * go_dim1];
		    b[mp1 + (np1 + k * b_dim2) * b_dim1] += zb[i__ + (np1 + 
			    i3 * zb_dim2) * zb_dim1] * go[i__ + ((mp1 << 1) - 
			    1 + k * go_dim2) * go_dim1];
/* L160: */
		}
	    }
	}
    }
    mp2 = mmax + 1;
    if (mdo == mmax || mp2 > ndo) {
	return 0;
    }
    zfin_(&c__1, nlat, &nlon, &mdo, &zb[zb_offset], &i3, &wzfin[1]);
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__1 = imm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = ndo;
	    for (np1 = mp2; np1 <= i__3; np1 += 2) {
		a[mmax + (np1 + k * a_dim2) * a_dim1] += zb[i__ + (np1 + i3 * 
			zb_dim2) * zb_dim1] * go[i__ + ((mmax << 1) - 2 + k * 
			go_dim2) * go_dim1];
/* L170: */
	    }
	}
    }
    return 0;
} /* shaec1_ */

/* Subroutine */ int shaeci_(integer *nlat, integer *nlon, doublereal *wshaec, 
	integer *lshaec, doublereal *dwork, integer *ldwork, integer *ierror)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, lzz1, labc, imid, mmax;
    extern /* Subroutine */ int hrffti_(integer *, doublereal *), zfinit_(integer *,
	     integer *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --wshaec;
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
    *ierror = 3;
    imid = (*nlat + 1) / 2;
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mmax = min(i__1,i__2);
    lzz1 = (*nlat << 1) * imid;
    labc = (mmax - 2) * (*nlat + *nlat - mmax - 1) * 3 / 2;
    if (*lshaec < lzz1 + labc + *nlon + 15) {
	return 0;
    }
    *ierror = 4;
    if (*ldwork < *nlat + 1) {
	return 0;
    }
    *ierror = 0;
    zfinit_(nlat, nlon, &wshaec[1], &dwork[1]);
    iw1 = lzz1 + labc + 1;
    hrffti_(nlon, &wshaec[iw1]);
    return 0;
} /* shaeci_ */

