/* slapgc.f -- translated by f2c (version 20061008).
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




/* ... file slapgc.f */

/*     this file includes documentation and code for */
/*     subroutine slapgc          i */

/* ... files which must be loaded with slapec.f */

/*     sphcom.f, hrfft.f, shagc.f, shsgc.f */



/*     subroutine slapgc(nlat,nlon,isym,nt,slap,ids,jds,a,b, */
/*    +mdab,ndab,wshsgc,lshsgc,work,lwork,ierror) */


/*     given the scalar spherical harmonic coefficients a and b, precomputed */
/*     by subroutine shagc for a scalar field sf, subroutine slapgc computes */
/*     the laplacian of sf in the scalar array slap.  slap(i,j) is the */
/*     laplacian of sf at the gaussian colatitude theta(i) (see nlat as */
/*     an input parameter) and east longitude lambda(j) = (j-1)*2*pi/nlon */
/*     on the sphere.  i.e. */

/*         slap(i,j) = */

/*                  2                2 */
/*         [1/sint*d (sf(i,j)/dlambda + d(sint*d(sf(i,j))/dtheta)/dtheta]/sint */


/*     where sint = sin(theta(i)).  the scalar laplacian in slap has the */
/*     same symmetry or absence of symmetry about the equator as the scalar */
/*     field sf.  the input parameters isym,nt,mdab,ndab must have the */
/*     same values used by shagc to compute a and b for sf. the associated */
/*     legendre functions are stored rather than recomputed as they are */
/*     in subroutine slapgc. */

/*     input parameters */

/*     nlat   the number of points in the gaussian colatitude grid on the */
/*            full sphere. these lie in the interval (0,pi) and are computed */
/*            in radians in theta(1) <...< theta(nlat) by subroutine gaqd. */
/*            if nlat is odd the equator will be included as the grid point */
/*            theta((nlat+1)/2).  if nlat is even the equator will be */
/*            excluded as a grid point and will lie half way between */
/*            theta(nlat/2) and theta(nlat/2+1). nlat must be at least 3. */
/*            note: on the half sphere, the number of grid points in the */
/*            colatitudinal direction is nlat/2 if nlat is even or */
/*            (nlat+1)/2 if nlat is odd. */

/*     nlon   the number of distinct longitude points.  nlon determines */
/*            the grid increment in longitude as 2*pi/nlon. for example */
/*            nlon = 72 for a five degree grid. nlon must be greater */
/*            than zero. the axisymmetric case corresponds to nlon=1. */
/*            the efficiency of the computation is improved when nlon */
/*            is a product of small prime numbers. */

/*     isym   this parameter should have the same value input to subroutine */
/*            shagc to compute the coefficients a and b for the scalar field */
/*            sf.  isym is set as follows: */

/*            = 0  no symmetries exist in sf about the equator. scalar */
/*                 synthesis is used to compute slap on the entire sphere. */
/*                 i.e., in the array slap(i,j) for i=1,...,nlat and */
/*                 j=1,...,nlon. */

/*           = 1  sf and slap are antisymmetric about the equator. the */
/*                synthesis used to compute slap is performed on the */
/*                northern hemisphere only.  if nlat is odd, slap(i,j) is */
/*                computed for i=1,...,(nlat+1)/2 and j=1,...,nlon.  if */
/*                nlat is even, slap(i,j) is computed for i=1,...,nlat/2 */
/*                and j=1,...,nlon. */


/*           = 2  sf and slap are symmetric about the equator. the */
/*                synthesis used to compute slap is performed on the */
/*                northern hemisphere only.  if nlat is odd, slap(i,j) is */
/*                computed for i=1,...,(nlat+1)/2 and j=1,...,nlon.  if */
/*                nlat is even, slap(i,j) is computed for i=1,...,nlat/2 */
/*                and j=1,...,nlon. */


/*     nt     the number of analyses.  in the program that calls slapgc */
/*            the arrays slap,a, and b can be three dimensional in which */
/*            case multiple synthesis will be performed.  the third index */
/*            is the synthesis index which assumes the values k=1,...,nt. */
/*            for a single analysis set nt=1. the description of the */
/*            remaining parameters is simplified by assuming that nt=1 */
/*            or that all the arrays are two dimensional. */

/*   ids      the first dimension of the array slap as it appears in the */
/*            program that calls slapgc.  if isym = 0 then ids must be at */
/*            least nlat.  if isym > 0 and nlat is even then ids must be */
/*            at least nlat/2. if isym > 0 and nlat is odd then ids must */
/*            be at least (nlat+1)/2. */

/*   jds      the second dimension of the array slap as it appears in the */
/*            program that calls slapgc. jds must be at least nlon. */


/*   a,b      two or three dimensional arrays (see input parameter nt) */
/*            that contain scalar spherical harmonic coefficients */
/*            of the scalar field sf as computed by subroutine shagc. */
/*     ***    a,b must be computed by shagc prior to calling slapgc. */


/*    mdab    the first dimension of the arrays a and b as it appears */
/*            in the program that calls slapgc.  mdab must be at */
/*            least min0(nlat,(nlon+2)/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*    ndab    the second dimension of the arrays a and b as it appears */
/*            in the program that calls slapgc. ndbc must be at least */
/*            least nlat. */

/*            mdab,ndab should have the same values input to shagc to */
/*            compute the coefficients a and b. */


/*    wshsgc  an array which must be initialized by subroutine shsgci. */
/*            once initialized, wshsgc */
/*            can be used repeatedly by slapgc as long as nlat and nlon */
/*            remain unchanged.  wshsgc must not be altered between calls */
/*            of slapgc. */

/*    lshsgc  the dimension of the array wshsgc as it appears in the */
/*            program that calls slapgc.  let */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshsgc must be at least */

/*               nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls slapgc. define */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */
/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            if isym = 0 let */

/*               lwkmin = nlat*(2*nt*nlon+max0(6*l2,nlon)+2*nt*l1+1. */

/*            if isym > 0 let */

/*               lwkmin = l2*(2*nt*nlon+max0(6*nlat,nlon))+nlat*(2*nt*l1+1) */


/*     then lwork must be greater than or equal to lwkmin (see ierror=10) */

/*     ************************************************************** */

/*     output parameters */


/*    slap    a two or three dimensional arrays (see input parameter nt) that */
/*            contain the scalar laplacian of the scalar field sf.  slap(i,j) */
/*            is the scalar laplacian at the gaussian colatitude theta(i) */
/*            and longitude lambda(j) = (j-1)*2*pi/nlon for i=1,...,nlat */
/*            and j=1,...,nlon. */


/*  ierror    a parameter which flags errors in input parameters as follows: */

/*            = 0  no errors detected */

/*            = 1  error in the specification of nlat */

/*            = 2  error in the specification of nlon */

/*            = 3  error in the specification of ityp */

/*            = 4  error in the specification of nt */

/*            = 5  error in the specification of ids */

/*            = 6  error in the specification of jds */

/*            = 7  error in the specification of mdbc */

/*            = 8  error in the specification of ndbc */

/*            = 9  error in the specification of lshsgc */

/*            = 10 error in the specification of lwork */


/* ********************************************************************** */

/*     end of documentation for slapgc */

/* ********************************************************************** */


/* Subroutine */ int slapgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *slap, integer *ids, integer *jds, doublereal 
	*a, doublereal *b, integer *mdab, integer *ndab, doublereal *wshsgc, 
	integer *lshsgc, doublereal *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer slap_dim1, slap_dim2, slap_offset, a_dim1, a_dim2, a_offset, 
	    b_dim1, b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    integer l1, l2, ia, ib, mn, ls, ifn, nln, iwk, lwk, imid, mmax, lwkmin;
    extern /* Subroutine */ int slapgc1_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);


/*     check input parameters */

    /* Parameter adjustments */
    slap_dim1 = *ids;
    slap_dim2 = *jds;
    slap_offset = 1 + slap_dim1 * (1 + slap_dim2);
    slap -= slap_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    --wshsgc;
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
    imid = (*nlat + 1) / 2;
    if (*isym == 0 && *ids < *nlat || *isym > 0 && *ids < imid) {
	return 0;
    }
    *ierror = 6;
    if (*jds < *nlon) {
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

/*     set and verify saved work space length */


/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    l1 = min(i__1,i__2);
    l2 = (*nlat + 1) / 2;
    if (*lshsgc < *nlat * ((l2 << 1) + l1 * 3 - 2) + l1 * 3 * (1 - l1) / 2 + *
	    nlon + 15) {
	return 0;
    }
    *ierror = 10;

/*     set and verify unsaved work space length */

    ls = *nlat;
    if (*isym > 0) {
	ls = imid;
    }
    nln = *nt * ls * *nlon;
    mn = mmax * *nlat * *nt;
/*     lwmin = nln+ls*nlon+2*mn+nlat */
/*     if (lwork .lt. lwmin) return */
    l2 = (*nlat + 1) / 2;
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    l1 = min(i__1,i__2);
    if (*isym == 0) {
/* Computing MAX */
	i__1 = l2 * 6;
	lwkmin = *nlat * ((*nt << 1) * *nlon + max(i__1,*nlon) + (*nt << 1) * 
		l1 + 1);
    } else {
/* Computing MAX */
	i__1 = *nlat * 6;
	lwkmin = l2 * ((*nt << 1) * *nlon + max(i__1,*nlon)) + *nlat * ((*nt 
		<< 1) * l1 + 1);
    }
    if (*lwork < lwkmin) {
	return 0;
    }
    *ierror = 0;

/*     set work space pointers */

    ia = 1;
    ib = ia + mn;
    ifn = ib + mn;
    iwk = ifn + *nlat;
    lwk = *lwork - (mn << 1) - *nlat;
    slapgc1_(nlat, nlon, isym, nt, &slap[slap_offset], ids, jds, &a[a_offset],
	     &b[b_offset], mdab, ndab, &work[ia], &work[ib], &mmax, &work[ifn]
	    , &wshsgc[1], lshsgc, &work[iwk], &lwk, ierror);
    return 0;
} /* slapgc_ */

/* Subroutine */ int slapgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *slap, integer *ids, integer *jds, doublereal 
	*a, doublereal *b, integer *mdab, integer *ndab, doublereal *alap, 
	doublereal *blap, integer *mmax, doublereal *fnn, doublereal *wsave, 
	integer *lsave, doublereal *wk, integer *lwk, integer *ierror)
{
    /* System generated locals */
    integer slap_dim1, slap_dim2, slap_offset, a_dim1, a_dim2, a_offset, 
	    b_dim1, b_dim2, b_offset, alap_dim1, alap_dim2, alap_offset, 
	    blap_dim1, blap_dim2, blap_offset, i__1, i__2, i__3;

    /* Local variables */
    integer k, m, n;
    doublereal fn;
    extern /* Subroutine */ int shsgc_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);


/*     set coefficient multiplyers */

    /* Parameter adjustments */
    --fnn;
    slap_dim1 = *ids;
    slap_dim2 = *jds;
    slap_offset = 1 + slap_dim1 * (1 + slap_dim2);
    slap -= slap_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    blap_dim1 = *mmax;
    blap_dim2 = *nlat;
    blap_offset = 1 + blap_dim1 * (1 + blap_dim2);
    blap -= blap_offset;
    alap_dim1 = *mmax;
    alap_dim2 = *nlat;
    alap_offset = 1 + alap_dim1 * (1 + alap_dim2);
    alap -= alap_offset;
    --wsave;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 2; n <= i__1; ++n) {
	fn = (doublereal) (n - 1);
	fnn[n] = fn * (fn + 1.);
/* L1: */
    }

/*     compute scalar laplacian coefficients for each vector field */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlat;
	for (n = 1; n <= i__2; ++n) {
	    i__3 = *mmax;
	    for (m = 1; m <= i__3; ++m) {
		alap[m + (n + k * alap_dim2) * alap_dim1] = 0.;
		blap[m + (n + k * blap_dim2) * blap_dim1] = 0.;
/* L4: */
	    }
/* L3: */
	}

/*     compute m=0 coefficients */

	i__2 = *nlat;
	for (n = 2; n <= i__2; ++n) {
	    alap[(n + k * alap_dim2) * alap_dim1 + 1] = -fnn[n] * a[(n + k * 
		    a_dim2) * a_dim1 + 1];
	    blap[(n + k * blap_dim2) * blap_dim1 + 1] = -fnn[n] * b[(n + k * 
		    b_dim2) * b_dim1 + 1];
/* L5: */
	}

/*     compute m>0 coefficients */

	i__2 = *mmax;
	for (m = 2; m <= i__2; ++m) {
	    i__3 = *nlat;
	    for (n = m; n <= i__3; ++n) {
		alap[m + (n + k * alap_dim2) * alap_dim1] = -fnn[n] * a[m + (
			n + k * a_dim2) * a_dim1];
		blap[m + (n + k * blap_dim2) * blap_dim1] = -fnn[n] * b[m + (
			n + k * b_dim2) * b_dim1];
/* L7: */
	    }
/* L6: */
	}
/* L2: */
    }

/*     synthesize alap,blap into slap */

    shsgc_(nlat, nlon, isym, nt, &slap[slap_offset], ids, jds, &alap[
	    alap_offset], &blap[blap_offset], mmax, nlat, &wsave[1], lsave, &
	    wk[1], lwk, ierror);
    return 0;
} /* slapgc1_ */

