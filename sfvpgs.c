/* sfvpgs.f -- translated by f2c (version 20061008).
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



/* ... file sfvpgs.f */

/*     this file includes documentation and code for */
/*     subroutine sfvpgs          i */

/* ... files which must be loaded with sfvpgs.f */

/*     sphcom.f, hrfft.f, vhags.f, shsgs.f, gaqd.f */


/*     subroutine sfvpgs(nlat,nlon,isym,nt,sf,vp,idv,jdv,br,bi,cr,ci, */
/*    +                   mdb,ndb,wshsgs,lshsgs,work,lwork,ierror) */

/*     given the vector spherical harmonic coefficients br,bi,cr,ci, */
/*     computed by subroutine vhags for a vector field (v,w), sfvpgs */
/*     computes a scalar stream function sf and scalar velocity potential */
/*     vp for (v,w).  (v,w) is expressed in terms of sf and vp by the */
/*     helmholtz relations (in mathematical spherical coordinates): */

/*          v = -1/sint*d(vp)/dlambda + d(st)/dtheta */

/*          w =  1/sint*d(st)/dlambda + d(vp)/dtheta */

/*     where sint = sin(theta).  w is the east longitudinal and v */
/*     is the colatitudinal component of the vector field from which */
/*     br,bi,cr,ci were precomputed.  required associated legendre */
/*     polynomials are stored rather than recomputed as they are in */
/*     subroutine sfvpgc. sf(i,j) and vp(i,j) are given at the i(th) */
/*     gaussian colatitude point theta(i) (see nlat description below) */
/*     and east longitude lambda(j) = (j-1)*2*pi/nlon on the sphere. */

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

/*     nlon   the number of distinct londitude points.  nlon determines */
/*            the grid increment in longitude as 2*pi/nlon. for example */
/*            nlon = 72 for a five degree grid. nlon must be greater than */
/*            3.  the efficiency of the computation is improved when nlon */
/*            is a product of small prime numbers. */


/*     isym   a parameter which determines whether the stream function and */
/*            velocity potential are computed on the full or half sphere */
/*            as follows: */

/*      = 0 */

/*            the symmetries/antsymmetries described in isym=1,2 below */
/*            do not exist in (v,w) about the equator.  in this case st */
/*            and vp are not necessarily symmetric or antisymmetric about */
/*            the equator.  sf and vp are computed on the entire sphere. */
/*            i.e., in arrays sf(i,j),vp(i,j) for i=1,...,nlat and */
/*            j=1,...,nlon. */

/*      = 1 */

/*            w is antisymmetric and v is symmetric about the equator. */
/*            in this case sf is symmetric and vp antisymmetric about */
/*            the equator and are computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the sf(i,j),vp(i,j) are computed */
/*            for i=1,...,(nlat+1)/2 and for j=1,...,nlon.  if nlat is */
/*            even then sf(i,j),vp(i,j) are computed for i=1,...,nlat/2 */
/*            and j=1,...,nlon. */

/*      = 2 */

/*            w is symmetric and v is antisymmetric about the equator. */
/*            in this case sf is antisymmetric and vp symmetric about */
/*            the equator and are computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the sf(i,j),vp(i,j) are computed */
/*            for i=1,...,(nlat+1)/2 and for j=1,...,nlon.  if nlat is */
/*            even then sf(i,j),vp(i,j) are computed for i=1,...,nlat/2 */
/*            and j=1,...,nlon. */

/*     nt     nt is the number of scalar and vector fields.  some */
/*            computational efficiency is obtained for multiple fields. arrays */
/*            can be three dimensional corresponding to an indexed multiple */
/*            vector field.  in this case multiple scalar synthesis will */
/*            be performed to compute sf,vp for each field.  the */
/*            third index is the synthesis index which assumes the values */
/*            k=1,...,nt.  for a single synthesis set nt = 1.  the */
/*            description of the remaining parameters is simplified by */
/*            assuming that nt=1 or that all the arrays are two dimensional. */

/*     idv    the first dimension of the arrays sf,vp as it appears in */
/*            the program that calls sfvpgs. if isym = 0 then idv */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then idv must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then idv must be at least (nlat+1)/2. */

/*     jdv    the second dimension of the arrays sf,vp as it appears in */
/*            the program that calls sfvpgs. jdv must be at least nlon. */

/*     br,bi, two or three dimensional arrays (see input parameter nt) */
/*     cr,ci  that contain vector spherical harmonic coefficients */
/*            of the vector field (v,w) as computed by subroutine vhags. */

/*     mdb    the first dimension of the arrays br,bi,cr,ci as it */
/*            appears in the program that calls sfvpgs. mdb must be at */
/*            least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndb    the second dimension of the arrays br,bi,cr,ci as it */
/*            appears in the program that calls sfvpgs. ndb must be at */
/*            least nlat. */

/*     wshsgs an array which must be initialized by subroutine shsgsi. */
/*            once initialized, wshsgs can be used repeatedly by sfvpgs */
/*            as long as nlon and nlat remain unchanged.  wshsgs must */
/*            not bel altered between calls of sfvpgs. */


/*     lshsgs the dimension of the array wshsgs as it appears in the */
/*            program that calls sfvpgs. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshsgs must be at least */

/*            nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls sfvpgs. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */

/*            if isym is zero then lwork must be at least */

/*                  nlat*nlon*(nt+1) + nlat*(2*l1*nt + 1) */

/*            if isym is nonzero then lwork must be at least */

/*                  l2*nlon*(nt+1) + nlat*(2*l1*nt + 1) */

/*     ************************************************************** */

/*     output parameters */

/*    sf,vp  two or three dimensional arrays (see input parameter nt) */
/*           that contains the stream function and velocity potential */
/*           of the vector field (v,w) whose coefficients br,bi,cr,ci */
/*           where precomputed by subroutine vhags.  sf(i,j),vp(i,j) */
/*           are given at the i(th) gaussian colatitude point theta(i) */
/*           and longitude point lambda(j) = (j-1)*2*pi/nlon.  the index */
/*           ranges are defined above at the input parameter isym. */


/*    ierror = 0  no errors */
/*           = 1  error in the specification of nlat */
/*           = 2  error in the specification of nlon */
/*           = 3  error in the specification of isym */
/*           = 4  error in the specification of nt */
/*           = 5  error in the specification of idv */
/*           = 6  error in the specification of jdv */
/*           = 7  error in the specification of mdb */
/*           = 8  error in the specification of ndb */
/*           = 9  error in the specification of lshsgs */
/*           = 10 error in the specification of lwork */
/* ********************************************************************** */

/* Subroutine */ int sfvpgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, doublereal *vp, integer *idv, integer *
	jdv, doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, 
	integer *mdb, integer *ndb, doublereal *wshsgs, integer *lshsgs, 
	doublereal *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer sf_dim1, sf_dim2, sf_offset, vp_dim1, vp_dim2, vp_offset, br_dim1,
	     br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, 
	    cr_dim2, cr_offset, ci_dim1, ci_dim2, ci_offset, i__1, i__2;

    /* Local variables */
    static integer l1, l2, ia, ib, mn, is, lp, ls, mab, lat, iwk, lwk, imid, 
	    late;
    extern /* Subroutine */ int stvpgs1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *);


/*     check input parameters */

    /* Parameter adjustments */
    vp_dim1 = *idv;
    vp_dim2 = *jdv;
    vp_offset = 1 + vp_dim1 * (1 + vp_dim2);
    vp -= vp_offset;
    sf_dim1 = *idv;
    sf_dim2 = *jdv;
    sf_offset = 1 + sf_dim1 * (1 + sf_dim2);
    sf -= sf_offset;
    ci_dim1 = *mdb;
    ci_dim2 = *ndb;
    ci_offset = 1 + ci_dim1 * (1 + ci_dim2);
    ci -= ci_offset;
    cr_dim1 = *mdb;
    cr_dim2 = *ndb;
    cr_offset = 1 + cr_dim1 * (1 + cr_dim2);
    cr -= cr_offset;
    bi_dim1 = *mdb;
    bi_dim2 = *ndb;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mdb;
    br_dim2 = *ndb;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
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
    if (*isym < 0 || *isym > 2) {
	return 0;
    }
    *ierror = 4;
    if (*nt < 0) {
	return 0;
    }
    *ierror = 5;
    imid = (*nlat + 1) / 2;
    if (*isym == 0 && *idv < *nlat || *isym > 0 && *idv < imid) {
	return 0;
    }
    *ierror = 6;
    if (*jdv < *nlon) {
	return 0;
    }
    *ierror = 7;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    if (*mdb < min(i__1,i__2)) {
	return 0;
    }
    *ierror = 8;
    if (*ndb < *nlat) {
	return 0;
    }
    *ierror = 9;
/* Computing MIN */
    i__1 = (*nlon + 2) / 2;
    l1 = min(i__1,*nlat);
    late = (*nlat + *nlat % 2) / 2;
    lat = *nlat;
    if (*isym != 0) {
	lat = late;
    }
    l2 = late;
/*     check permanent work space length */
    l2 = (*nlat + *nlat % 2) / 2;
/* Computing MIN */
    i__1 = (*nlon + 2) / 2;
    l1 = min(i__1,*nlat);
    lp = *nlat * ((l1 + l2) * 3 - 2) + (l1 - 1) * (l2 * ((*nlat << 1) - l1) - 
	    l1 * 3) / 2 + *nlon + 15;
    if (*lshsgs < lp) {
	return 0;
    }

/*     verify unsaved work space */

    *ierror = 10;
    ls = *nlat;
    if (*isym > 0) {
	ls = imid;
    }

/*     set first dimension for a,b (as requried by shsgs) */

/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mab = min(i__1,i__2);
    mn = mab * *nlat * *nt;
    if (*lwork < ls * *nlon + (*nt + 1) + *nlat * ((l1 << 1) * *nt + 1)) {
	return 0;
    }
    *ierror = 0;

/*     set work space pointers */

    ia = 1;
    ib = ia + mn;
    is = ib + mn;
    iwk = is + *nlat;
    lwk = *lwork - (mn << 1) - *nlat;
    stvpgs1_(nlat, nlon, isym, nt, &sf[sf_offset], &vp[vp_offset], idv, jdv, &
	    br[br_offset], &bi[bi_offset], &cr[cr_offset], &ci[ci_offset], 
	    mdb, ndb, &work[ia], &work[ib], &mab, &work[is], &wshsgs[1], 
	    lshsgs, &work[iwk], &lwk, ierror);
    return 0;
} /* sfvpgs_ */

/* Subroutine */ int stvpgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, doublereal *vp, integer *idv, integer *
	jdv, doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, 
	integer *mdb, integer *ndb, doublereal *a, doublereal *b, integer *
	mab, doublereal *fnn, doublereal *wshsgs, integer *lshsgs, doublereal 
	*wk, integer *lwk, integer *ierror)
{
    /* System generated locals */
    integer sf_dim1, sf_dim2, sf_offset, vp_dim1, vp_dim2, vp_offset, br_dim1,
	     br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, 
	    cr_dim2, cr_offset, ci_dim1, ci_dim2, ci_offset, a_dim1, a_dim2, 
	    a_offset, b_dim1, b_dim2, b_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, m, n, mmax;
    extern /* Subroutine */ int shsgs_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);


/*     set coefficient multiplyers */

    /* Parameter adjustments */
    --fnn;
    vp_dim1 = *idv;
    vp_dim2 = *jdv;
    vp_offset = 1 + vp_dim1 * (1 + vp_dim2);
    vp -= vp_offset;
    sf_dim1 = *idv;
    sf_dim2 = *jdv;
    sf_offset = 1 + sf_dim1 * (1 + sf_dim2);
    sf -= sf_offset;
    ci_dim1 = *mdb;
    ci_dim2 = *ndb;
    ci_offset = 1 + ci_dim1 * (1 + ci_dim2);
    ci -= ci_offset;
    cr_dim1 = *mdb;
    cr_dim2 = *ndb;
    cr_offset = 1 + cr_dim1 * (1 + cr_dim2);
    cr -= cr_offset;
    bi_dim1 = *mdb;
    bi_dim2 = *ndb;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mdb;
    br_dim2 = *ndb;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
    b_dim1 = *mab;
    b_dim2 = *nlat;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mab;
    a_dim2 = *nlat;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    --wshsgs;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 2; n <= i__1; ++n) {
	fnn[n] = 1. / sqrt((doublereal) (n * (n - 1)));
    }
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);

/*     compute st scalar coefficients from cr,ci */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlat;
	for (n = 1; n <= i__2; ++n) {
	    i__3 = *mab;
	    for (m = 1; m <= i__3; ++m) {
		a[m + (n + k * a_dim2) * a_dim1] = 0.;
		b[m + (n + k * b_dim2) * b_dim1] = 0.;
	    }
	}

/*     compute m=0 coefficients */

	i__2 = *nlat;
	for (n = 2; n <= i__2; ++n) {
	    a[(n + k * a_dim2) * a_dim1 + 1] = -fnn[n] * cr[(n + k * cr_dim2) 
		    * cr_dim1 + 1];
	    b[(n + k * b_dim2) * b_dim1 + 1] = -fnn[n] * ci[(n + k * ci_dim2) 
		    * ci_dim1 + 1];
	}

/*     compute m>0 coefficients using vector spherepack value for mmax */

	i__2 = mmax;
	for (m = 2; m <= i__2; ++m) {
	    i__3 = *nlat;
	    for (n = m; n <= i__3; ++n) {
		a[m + (n + k * a_dim2) * a_dim1] = -fnn[n] * cr[m + (n + k * 
			cr_dim2) * cr_dim1];
		b[m + (n + k * b_dim2) * b_dim1] = -fnn[n] * ci[m + (n + k * 
			ci_dim2) * ci_dim1];
	    }
	}
    }

/*     synthesize a,b into st */

    shsgs_(nlat, nlon, isym, nt, &sf[sf_offset], idv, jdv, &a[a_offset], &b[
	    b_offset], mab, nlat, &wshsgs[1], lshsgs, &wk[1], lwk, ierror);

/*    set coefficients for vp from br,bi */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlat;
	for (n = 1; n <= i__2; ++n) {
	    i__3 = *mab;
	    for (m = 1; m <= i__3; ++m) {
		a[m + (n + k * a_dim2) * a_dim1] = 0.;
		b[m + (n + k * b_dim2) * b_dim1] = 0.;
	    }
	}

/*     compute m=0 coefficients */

	i__2 = *nlat;
	for (n = 2; n <= i__2; ++n) {
	    a[(n + k * a_dim2) * a_dim1 + 1] = fnn[n] * br[(n + k * br_dim2) *
		     br_dim1 + 1];
	    b[(n + k * b_dim2) * b_dim1 + 1] = fnn[n] * bi[(n + k * bi_dim2) *
		     bi_dim1 + 1];
	}

/*     compute m>0 coefficients using vector spherepack value for mmax */

/* Computing MIN */
	i__2 = *nlat, i__3 = (*nlon + 1) / 2;
	mmax = min(i__2,i__3);
	i__2 = mmax;
	for (m = 2; m <= i__2; ++m) {
	    i__3 = *nlat;
	    for (n = m; n <= i__3; ++n) {
		a[m + (n + k * a_dim2) * a_dim1] = fnn[n] * br[m + (n + k * 
			br_dim2) * br_dim1];
		b[m + (n + k * b_dim2) * b_dim1] = fnn[n] * bi[m + (n + k * 
			bi_dim2) * bi_dim1];
	    }
	}
    }

/*     synthesize a,b into vp */

    shsgs_(nlat, nlon, isym, nt, &vp[vp_offset], idv, jdv, &a[a_offset], &b[
	    b_offset], mab, nlat, &wshsgs[1], lshsgs, &wk[1], lwk, ierror);
    return 0;
} /* stvpgs1_ */

