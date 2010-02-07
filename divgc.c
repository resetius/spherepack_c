/* divgc.f -- translated by f2c (version 20061008).
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



/* ... file divgc.f */

/*     this file includes documentation and code for */
/*     subroutine divgc          i */

/* ... files which must be loaded with divgc.f */

/*     sphcom.f, hrfft.f, vhagc.f, shsgc.f, gaqd.f */


/*     subroutine divgc(nlat,nlon,isym,nt,dv,idv,jdv,br,bi,mdb,ndb, */
/*    +                 wshsgc,lshsgc,work,lwork,ierror) */

/*     given the vector spherical harmonic coefficients br and bi, precomputed */
/*     by subroutine vhagc for a vector field (v,w), subroutine divgc */
/*     computes the divergence of the vector field in the scalar array dv. */
/*     dv(i,j) is the divergence at the gaussian colatitude point theta(i) */
/*     (see nlat as input parameter) and east longitude */

/*            lambda(j) = (j-1)*2*pi/nlon */

/*     on the sphere.  i.e. */

/*            dv(i,j) = 1/sint*[ d(sint*v(i,j))/dtheta + d(w(i,j))/dlambda ] */

/*     where sint = sin(theta(i)).  w is the east longitudinal and v */
/*     is the colatitudinal component of the vector field from which */
/*     br,bi were precomputed.  required associated legendre polynomials */
/*     are recomputed rather than stored as they are in subroutine divgs. */


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
/*            nlon = 72 for a five degree grid. nlon must be greater */
/*            than zero. the axisymmetric case corresponds to nlon=1. */
/*            the efficiency of the computation is improved when nlon */
/*            is a product of small prime numbers. */


/*     isym   a parameter which determines whether the divergence is */
/*            computed on the full or half sphere as follows: */

/*      = 0 */

/*            the symmetries/antsymmetries described in isym=1,2 below */
/*            do not exist in (v,w) about the equator.  in this case the */
/*            divergence is neither symmetric nor antisymmetric about */
/*            the equator.  the divergence is computed on the entire */
/*            sphere.  i.e., in the array dv(i,j) for i=1,...,nlat and */
/*            j=1,...,nlon. */

/*      = 1 */

/*            w is antisymmetric and v is symmetric about the equator. */
/*            in this case the divergence is antisymmetyric about */
/*            the equator and is computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the divergence is computed */
/*            in the array dv(i,j) for i=1,...,(nlat+1)/2 and for */
/*            j=1,...,nlon.  if nlat is even the divergence is computed */
/*            in the array dv(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*      = 2 */
/*            w is symmetric and v is antisymmetric about the equator */
/*            in this case the divergence is symmetyric about the */
/*            equator and is computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the divergence is computed */
/*            in the array dv(i,j) for i=1,...,(nlat+1)/2 and for */
/*            j=1,...,nlon.  if nlat is even the divergence is computed */
/*            in the array dv(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */


/*     nt     nt is the number of scalar and vector fields.  some */
/*            computational efficiency is obtained for multiple fields. */
/*            in the program that calls divgc, the arrays br,bi, and dv */
/*            can be three dimensional corresponding to an indexed multiple */
/*            vector field.  in this case multiple scalar synthesis will */
/*            be performed to compute the divergence for each field.  the */
/*            third index is the synthesis index which assumes the values */
/*            k=1,...,nt.  for a single synthesis set nt = 1.  the */
/*            description of the remaining parameters is simplified by */
/*            assuming that nt=1 or that all the arrays are two dimensional. */

/*     idv    the first dimension of the array dv as it appears in */
/*            the program that calls divgc. if isym = 0 then idv */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then idv must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then idv must be at least (nlat+1)/2. */

/*     jdv    the second dimension of the array dv as it appears in */
/*            the program that calls divgc. jdv must be at least nlon. */

/*     br,bi  two or three dimensional arrays (see input parameter nt) */
/*            that contain vector spherical harmonic coefficients */
/*            of the vector field (v,w) as computed by subroutine vhagc. */
/*     ***    br and bi must be computed by vhagc prior to calling */
/*            divgc. */

/*     mdb    the first dimension of the arrays br and bi as it */
/*            appears in the program that calls divgc. mdb must be at */
/*            least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndb    the second dimension of the arrays br and bi as it */
/*            appears in the program that calls divgc. ndb must be at */
/*            least nlat. */


/*     wshsgc an array which must be initialized by subroutine shsgci */
/*            once initialized, wshsgc can be used repeatedly by divgc */
/*            as long as nlon and nlat remain unchanged.  wshsgc must */
/*            not be altered between calls of divgc. */


/*     lshsgc the dimension of the array wshsgc as it appears in the */
/*            program that calls divgc. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshsgc must be at least */

/*               nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15 */

/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls divgc. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */


/*            if isym is zero then lwork must be at least */

/*               nlat*(nlon*nt+max0(3*l2,nlon) + 2*nt*l1+1) */

/*            if isym is not zero then lwork must be at least */

/*               l2*(nlon*nt+max0(3*nlat,nlon)) + nlat*(2*nt*l1+1) */


/*     ************************************************************** */

/*     output parameters */


/*    dv     a two or three dimensional array (see input parameter nt) */
/*           that contains the divergence of the vector field (v,w) */
/*           whose coefficients br,bi where computed by subroutine */
/*           vhagc.  dv(i,j) is the divergence at the gaussian colatitude */
/*           point theta(i) and longitude point lambda(j) = (j-1)*2*pi/nlon. */
/*           the index ranges are defined above at the input parameter */
/*           isym. */


/*    ierror = 0  no errors */
/*           = 1  error in the specification of nlat */
/*           = 2  error in the specification of nlon */
/*           = 3  error in the specification of isym */
/*           = 4  error in the specification of nt */
/*           = 5  error in the specification of idv */
/*           = 6  error in the specification of jdv */
/*           = 7  error in the specification of mdb */
/*           = 8  error in the specification of ndb */
/*           = 9  error in the specification of lshsgc */
/*           = 10 error in the specification of lwork */
/* ********************************************************************** */


/* Subroutine */ int divgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, real *dv, integer *idv, integer *jdv, real *br, real *bi,
	 integer *mdb, integer *ndb, real *wshsgc, integer *lshsgc, real *
	work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer dv_dim1, dv_dim2, dv_offset, br_dim1, br_dim2, br_offset, bi_dim1,
	     bi_dim2, bi_offset, i__1, i__2;

    /* Local variables */
    static integer l1, l2, ia, ib, mn, is, ls, mab, nln, iwk, lwk, imid, mmax,
	     lpimn;
    extern /* Subroutine */ int divgc1_(integer *, integer *, integer *, 
	    integer *, real *, integer *, integer *, real *, real *, integer *
	    , integer *, real *, real *, integer *, real *, real *, integer *,
	     real *, integer *, integer *);
    static integer lwkmin;


/*     check input parameters */

    /* Parameter adjustments */
    dv_dim1 = *idv;
    dv_dim2 = *jdv;
    dv_offset = 1 + dv_dim1 * (1 + dv_dim2);
    dv -= dv_offset;
    bi_dim1 = *mdb;
    bi_dim2 = *ndb;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mdb;
    br_dim2 = *ndb;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
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
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    mmax = min(i__1,i__2);
    *ierror = 8;
    if (*ndb < *nlat) {
	return 0;
    }
    *ierror = 9;
    imid = (*nlat + 1) / 2;
    lpimn = imid * mmax * (*nlat + *nlat - mmax + 1) / 2;
/*     check permanent work space length */
    l2 = (*nlat + 1) / 2;
/* Computing MIN */
    i__1 = (*nlon + 2) / 2;
    l1 = min(i__1,*nlat);
    if (*lshsgc < *nlat * ((l2 << 1) + l1 * 3 - 2) + l1 * 3 * (1 - l1) / 2 + *
	    nlon + 15) {
	return 0;
    }
    *ierror = 10;

/*     verify unsaved work space (add to what shsgc requires) */

    ls = *nlat;
    if (*isym > 0) {
	ls = imid;
    }
    nln = *nt * ls * *nlon;

/*     set first dimension for a,b (as requried by shsgc) */

/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mab = min(i__1,i__2);
    mn = mab * *nlat * *nt;
/*     if(lwork.lt. nln+ls*nlon+2*mn+nlat) return */
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    l1 = min(i__1,i__2);
    l2 = (*nlat + 1) / 2;
    if (*isym == 0) {
/* Computing MAX */
	i__1 = l2 * 3;
	lwkmin = *nlat * (*nt * *nlon + max(i__1,*nlon) + (*nt << 1) * l1 + 1)
		;
    } else {
/* Computing MAX */
	i__1 = *nlat * 3;
	lwkmin = l2 * (*nt * *nlon + max(i__1,*nlon)) + *nlat * ((*nt << 1) * 
		l1 + 1);
    }
    if (*lwork < lwkmin) {
	return 0;
    }
    *ierror = 0;

/*     set work space pointers */

    ia = 1;
    ib = ia + mn;
    is = ib + mn;
    iwk = is + *nlat;
    lwk = *lwork - (mn << 1) - *nlat;
    divgc1_(nlat, nlon, isym, nt, &dv[dv_offset], idv, jdv, &br[br_offset], &
	    bi[bi_offset], mdb, ndb, &work[ia], &work[ib], &mab, &work[is], &
	    wshsgc[1], lshsgc, &work[iwk], &lwk, ierror);
    return 0;
} /* divgc_ */

/* Subroutine */ int divgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, real *dv, integer *idv, integer *jdv, real *br, real *bi,
	 integer *mdb, integer *ndb, real *a, real *b, integer *mab, real *
	sqnn, real *wshsgc, integer *lshsgc, real *wk, integer *lwk, integer *
	ierror)
{
    /* System generated locals */
    integer dv_dim1, dv_dim2, dv_offset, br_dim1, br_dim2, br_offset, bi_dim1,
	     bi_dim2, bi_offset, a_dim1, a_dim2, a_offset, b_dim1, b_dim2, 
	    b_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, m, n;
    static real fn;
    static integer mmax;
    extern /* Subroutine */ int shsgc_(integer *, integer *, integer *, 
	    integer *, real *, integer *, integer *, real *, real *, integer *
	    , integer *, real *, integer *, real *, integer *, integer *);


/*     set coefficient multiplyers */

    /* Parameter adjustments */
    --sqnn;
    dv_dim1 = *idv;
    dv_dim2 = *jdv;
    dv_offset = 1 + dv_dim1 * (1 + dv_dim2);
    dv -= dv_offset;
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
    --wshsgc;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 2; n <= i__1; ++n) {
	fn = (real) (n - 1);
	sqnn[n] = sqrt(fn * (fn + 1.f));
/* L1: */
    }

/*     compute divergence scalar coefficients for each vector field */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlat;
	for (n = 1; n <= i__2; ++n) {
	    i__3 = *mab;
	    for (m = 1; m <= i__3; ++m) {
		a[m + (n + k * a_dim2) * a_dim1] = 0.f;
		b[m + (n + k * b_dim2) * b_dim1] = 0.f;
/* L4: */
	    }
/* L3: */
	}

/*     compute m=0 coefficients */

	i__2 = *nlat;
	for (n = 2; n <= i__2; ++n) {
	    a[(n + k * a_dim2) * a_dim1 + 1] = -sqnn[n] * br[(n + k * br_dim2)
		     * br_dim1 + 1];
	    b[(n + k * b_dim2) * b_dim1 + 1] = -sqnn[n] * bi[(n + k * bi_dim2)
		     * bi_dim1 + 1];
/* L5: */
	}

/*     compute m>0 coefficients using vector spherepack value for mmax */

/* Computing MIN */
	i__2 = *nlat, i__3 = (*nlon + 1) / 2;
	mmax = min(i__2,i__3);
	i__2 = mmax;
	for (m = 2; m <= i__2; ++m) {
	    i__3 = *nlat;
	    for (n = m; n <= i__3; ++n) {
		a[m + (n + k * a_dim2) * a_dim1] = -sqnn[n] * br[m + (n + k * 
			br_dim2) * br_dim1];
		b[m + (n + k * b_dim2) * b_dim1] = -sqnn[n] * bi[m + (n + k * 
			bi_dim2) * bi_dim1];
/* L7: */
	    }
/* L6: */
	}
/* L2: */
    }

/*     synthesize a,b into dv */

    shsgc_(nlat, nlon, isym, nt, &dv[dv_offset], idv, jdv, &a[a_offset], &b[
	    b_offset], mab, nlat, &wshsgc[1], lshsgc, &wk[1], lwk, ierror);
    return 0;
} /* divgc1_ */

