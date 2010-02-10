/* shses.f -- translated by f2c (version 20061008).
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
/*  .                         SPHEREPACK                       . */
/*  .                                                             . */
/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */



/* ... file shses.f */

/*     this file contains code and documentation for subroutines */
/*     shses and shsesi */

/* ... files which must be loaded with shses.f */

/*     sphcom.f, hrfft.f */

/*     subroutine shses(nlat,nlon,isym,nt,g,idg,jdg,a,b,mdab,ndab, */
/*    +                 wshses,lshses,work,lwork,ierror) */

/*     subroutine shses performs the spherical harmonic synthesis */
/*     on the arrays a and b and stores the result in the array g. */
/*     the synthesis is performed on an equally spaced grid.  the */
/*     associated legendre functions are stored rather than recomputed */
/*     as they are in subroutine shsec.  the synthesis is described */
/*     below at output parameter g. */

/* *** required files from spherepack2 */

/*     sphcom.f, hrfft.f */


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

/*     nt     the number of syntheses.  in the program that calls shses, */
/*            the arrays g,a and b can be three dimensional in which */
/*            case multiple syntheses will be performed.  the third */
/*            index is the synthesis index which assumes the values */
/*            k=1,...,nt.  for a single synthesis set nt=1. the */
/*            discription of the remaining parameters is simplified */
/*            by assuming that nt=1 or that the arrays g,a and b */
/*            have only two dimensions. */

/*     idg    the first dimension of the array g as it appears in the */
/*            program that calls shses.  if isym equals zero then idg */
/*            must be at least nlat.  if isym is nonzero then idg */
/*            must be at least nlat/2 if nlat is even or at least */
/*            (nlat+1)/2 if nlat is odd. */

/*     jdg    the second dimension of the array g as it appears in the */
/*            program that calls shses.  jdg must be at least nlon. */

/*     a,b    two or three dimensional arrays (see the input parameter */
/*            nt) that contain the coefficients in the spherical harmonic */
/*            expansion of g(i,j) given below at the definition of the */
/*            output parameter g.  a(m,n) and b(m,n) are defined for */
/*            indices m=1,...,mmax and n=m,...,nlat where mmax is the */
/*            maximum (plus one) longitudinal wave number given by */
/*            mmax = min0(nlat,(nlon+2)/2) if nlon is even or */
/*            mmax = min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     mdab   the first dimension of the arrays a and b as it appears */
/*            in the program that calls shses. mdab must be at least */
/*            min0(nlat,(nlon+2)/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndab   the second dimension of the arrays a and b as it appears */
/*            in the program that calls shses. ndab must be at least nlat */

/*     wshses an array which must be initialized by subroutine shsesi. */
/*            once initialized, wshses can be used repeatedly by shses */
/*            as long as nlon and nlat remain unchanged.  wshses must */
/*            not be altered between calls of shses. */

/*     lshses the dimension of the array wshses as it appears in the */
/*            program that calls shses. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshses must be at least */

/*               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15 */

/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls shses.  define */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            if isym is zero then lwork must be at least */

/*               (nt+1)*nlat*nlon */

/*            if isym is nonzero lwork must be at least */

/*               (nt+1)*l2*nlon. */

/*     ************************************************************** */

/*     output parameters */

/*     g      a two or three dimensional array (see input parameter */
/*            nt) that contains the spherical harmonic synthesis of */
/*            the arrays a and b at the colatitude point theta(i) = */
/*            (i-1)*pi/(nlat-1) and longitude point phi(j) = */
/*            (j-1)*2*pi/nlon. the index ranges are defined above at */
/*            at the input parameter isym.  for isym=0, g(i,j) is */
/*            given by the the equations listed below.  symmetric */
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
/*            = 9  error in the specification of lshses */
/*            = 10 error in the specification of lwork */


/* **************************************************************** */
/*     subroutine shsesi(nlat,nlon,wshses,lshses,work,lwork,dwork, */
/*    +                  ldwork,ierror) */

/*     subroutine shsesi initializes the array wshses which can then */
/*     be used repeatedly by subroutine shses. */

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

/*     lshses the dimension of the array wshses as it appears in the */
/*            program that calls shsesi. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshses must be at least */

/*               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15 */

/*     work   a real   work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in */
/*            the program that calls shsesi.  define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lwork must be at least */

/*               5*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2 */


/*     dwork  a double precision work array that does not have to be saved. */

/*     ldwork the dimension of the array dwork as it appears in the */
/*            program that calls shsesi.  ldwork must be at least nlat+1 */


/*     output parameters */

/*     wshses an array which is initialized for use by subroutine shses. */
/*            once initialized, wshses can be used repeatedly by shses */
/*            as long as nlon and nlat remain unchanged.  wshses must */
/*            not be altered between calls of shses. */

/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of lshses */
/*            = 4  error in the specification of lwork */
/*            = 5  error in the specification of ldwork */

/* **************************************************************** */
/* Subroutine */ int shses_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a,
	 doublereal *b, integer *mdab, integer *ndab, doublereal *wshses, 
	integer *lshses, doublereal *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer g_dim1, g_dim2, g_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    static integer ls, nln, ist, imid, mmax, lpimn;
    extern /* Subroutine */ int shses1_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *);

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
    --wshses;
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
    lpimn = imid * mmax * (*nlat + *nlat - mmax + 1) / 2;
    if (*lshses < lpimn + *nlon + 15) {
	return 0;
    }
    *ierror = 10;
    ls = *nlat;
    if (*isym > 0) {
	ls = imid;
    }
    nln = *nt * ls * *nlon;
    if (*lwork < nln + ls * *nlon) {
	return 0;
    }
    *ierror = 0;
    ist = 0;
    if (*isym == 0) {
	ist = imid;
    }
    shses1_(nlat, isym, nt, &g[g_offset], idg, jdg, &a[a_offset], &b[b_offset]
	    , mdab, ndab, &wshses[1], &imid, &ls, nlon, &work[1], &work[ist + 
	    1], &work[nln + 1], &wshses[lpimn + 1]);
    return 0;
} /* shses_ */

/* Subroutine */ int shses1_(integer *nlat, integer *isym, integer *nt, 
	doublereal *g, integer *idgs, integer *jdgs, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *p, integer *
	imid, integer *idg, integer *jdg, doublereal *ge, doublereal *go, 
	doublereal *work, doublereal *whrfft)
{
    /* System generated locals */
    integer g_dim1, g_dim2, g_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, p_dim1, p_offset, ge_dim1, ge_dim2, ge_offset, 
	    go_dim1, go_dim2, go_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, m, mb, mn, ls, mp1, np1, mp2, mdo, ndo, imm1, 
	    nlp1, modl, mmax, nlon;
    extern /* Subroutine */ int hrfftb_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *);

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
    p_dim1 = *imid;
    p_offset = 1 + p_dim1;
    p -= p_offset;
    go_dim1 = *idg;
    go_dim2 = *jdg;
    go_offset = 1 + go_dim1 * (1 + go_dim2);
    go -= go_offset;
    ge_dim1 = *idg;
    ge_dim2 = *jdg;
    ge_offset = 1 + ge_dim1 * (1 + ge_dim2);
    ge -= ge_offset;
    --work;
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
    modl = *nlat % 2;
    imm1 = *imid;
    if (modl != 0) {
	imm1 = *imid - 1;
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = nlon;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = ls;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		ge[i__ + (j + k * ge_dim2) * ge_dim1] = 0.;
/* L8000: */
/* L800: */
/* L80: */
	    }
	}
    }
    if (*isym == 1) {
	goto L125;
    }
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = *nlat;
	for (np1 = 1; np1 <= i__2; np1 += 2) {
	    i__1 = *imid;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ge[i__ + (k * ge_dim2 + 1) * ge_dim1] += a[(np1 + k * a_dim2) 
			* a_dim1 + 1] * p[i__ + np1 * p_dim1];
/* L100: */
	    }
	}
    }
    ndo = *nlat;
    if (*nlat % 2 == 0) {
	ndo = *nlat - 1;
    }
    i__1 = mdo;
    for (mp1 = 2; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	mb = m * (*nlat - 1) - m * (m - 1) / 2;
	i__2 = ndo;
	for (np1 = mp1; np1 <= i__2; np1 += 2) {
	    mn = mb + np1;
	    i__3 = *nt;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = *imid;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    ge[i__ + ((mp1 << 1) - 2 + k * ge_dim2) * ge_dim1] += a[
			    mp1 + (np1 + k * a_dim2) * a_dim1] * p[i__ + mn * 
			    p_dim1];
		    ge[i__ + ((mp1 << 1) - 1 + k * ge_dim2) * ge_dim1] += b[
			    mp1 + (np1 + k * b_dim2) * b_dim1] * p[i__ + mn * 
			    p_dim1];
/* L110: */
		}
	    }
	}
    }
    if (mdo == mmax || mmax > ndo) {
	goto L122;
    }
    mb = mdo * (*nlat - 1) - mdo * (mdo - 1) / 2;
    i__4 = ndo;
    for (np1 = mmax; np1 <= i__4; np1 += 2) {
	mn = mb + np1;
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__2 = *imid;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ge[i__ + ((mmax << 1) - 2 + k * ge_dim2) * ge_dim1] += a[mmax 
			+ (np1 + k * a_dim2) * a_dim1] * p[i__ + mn * p_dim1];
/* L120: */
	    }
	}
    }
L122:
    if (*isym == 2) {
	goto L155;
    }
L125:
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__3 = *nlat;
	for (np1 = 2; np1 <= i__3; np1 += 2) {
	    i__4 = imm1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		go[i__ + (k * go_dim2 + 1) * go_dim1] += a[(np1 + k * a_dim2) 
			* a_dim1 + 1] * p[i__ + np1 * p_dim1];
/* L140: */
	    }
	}
    }
    ndo = *nlat;
    if (*nlat % 2 != 0) {
	ndo = *nlat - 1;
    }
    i__4 = mdo;
    for (mp1 = 2; mp1 <= i__4; ++mp1) {
	mp2 = mp1 + 1;
	m = mp1 - 1;
	mb = m * (*nlat - 1) - m * (m - 1) / 2;
	i__3 = ndo;
	for (np1 = mp2; np1 <= i__3; np1 += 2) {
	    mn = mb + np1;
	    i__2 = *nt;
	    for (k = 1; k <= i__2; ++k) {
		i__1 = imm1;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    go[i__ + ((mp1 << 1) - 2 + k * go_dim2) * go_dim1] += a[
			    mp1 + (np1 + k * a_dim2) * a_dim1] * p[i__ + mn * 
			    p_dim1];
		    go[i__ + ((mp1 << 1) - 1 + k * go_dim2) * go_dim1] += b[
			    mp1 + (np1 + k * b_dim2) * b_dim1] * p[i__ + mn * 
			    p_dim1];
/* L150: */
		}
	    }
	}
    }
    mp2 = mmax + 1;
    if (mdo == mmax || mp2 > ndo) {
	goto L155;
    }
    mb = mdo * (*nlat - 1) - mdo * (mdo - 1) / 2;
    i__1 = ndo;
    for (np1 = mp2; np1 <= i__1; np1 += 2) {
	mn = mb + np1;
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		go[i__ + ((mmax << 1) - 2 + k * go_dim2) * go_dim1] += a[mmax 
			+ (np1 + k * a_dim2) * a_dim1] * p[i__ + mn * p_dim1];
/* L152: */
	    }
	}
    }
L155:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	if (nlon % 2 != 0) {
	    goto L157;
	}
	i__2 = ls;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ge[i__ + (nlon + k * ge_dim2) * ge_dim1] *= 2.;
/* L156: */
	}
L157:
	hrfftb_(&ls, &nlon, &ge[(k * ge_dim2 + 1) * ge_dim1 + 1], &ls, &
		whrfft[1], &work[1]);
/* L160: */
    }
    if (*isym != 0) {
	goto L180;
    }
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = nlon;
	for (j = 1; j <= i__2; ++j) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		g[i__ + (j + k * g_dim2) * g_dim1] = (ge[i__ + (j + k * 
			ge_dim2) * ge_dim1] + go[i__ + (j + k * go_dim2) * 
			go_dim1]) * .5;
		g[nlp1 - i__ + (j + k * g_dim2) * g_dim1] = (ge[i__ + (j + k *
			 ge_dim2) * ge_dim1] - go[i__ + (j + k * go_dim2) * 
			go_dim1]) * .5;
/* L175: */
	    }
	    if (modl == 0) {
		goto L170;
	    }
	    g[*imid + (j + k * g_dim2) * g_dim1] = ge[*imid + (j + k * 
		    ge_dim2) * ge_dim1] * .5;
L170:
	    ;
	}
    }
    return 0;
L180:
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__1 = nlon;
	    for (j = 1; j <= i__1; ++j) {
		g[i__ + (j + k * g_dim2) * g_dim1] = ge[i__ + (j + k * 
			ge_dim2) * ge_dim1] * .5;
/* L185: */
	    }
	}
    }
    return 0;
} /* shses1_ */

/* Subroutine */ int shsesi_(integer *nlat, integer *nlon, doublereal *wshses,
	 integer *lshses, doublereal *work, integer *lwork, doublereal *dwork,
	 integer *ldwork, integer *ierror)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1;
    extern /* Subroutine */ int ses1_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer labc, imid, mmax, lpimn;
    extern /* Subroutine */ int hrffti_(integer *, doublereal *);

    /* Parameter adjustments */
    --dwork;
    --work;
    --wshses;

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
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mmax = min(i__1,i__2);
    imid = (*nlat + 1) / 2;
    lpimn = imid * mmax * (*nlat + *nlat - mmax + 1) / 2;
    if (*lshses < lpimn + *nlon + 15) {
	return 0;
    }
    *ierror = 4;
    labc = (mmax - 2) * (*nlat + *nlat - mmax - 1) * 3 / 2;
    if (*lwork < *nlat * 5 * imid + labc) {
	return 0;
    }
    *ierror = 5;
    if (*ldwork < *nlat + 1) {
	return 0;
    }
    *ierror = 0;
    iw1 = *nlat * 3 * imid + 1;
    ses1_(nlat, nlon, &imid, &wshses[1], &work[1], &work[iw1], &dwork[1]);
    hrffti_(nlon, &wshses[lpimn + 1]);
    return 0;
} /* shsesi_ */

