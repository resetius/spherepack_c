/* divec.f -- translated by f2c (version 20061008).
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



/* ... file divec.f */

/*     this file includes documentation and code for */
/*     subroutine divec          i */

/* ... files which must be loaded with divec.f */

/*     sphcom.f, hrfft.f, vhaec.f,shsec.f */


/*     subroutine divec(nlat,nlon,isym,nt,dv,idv,jdv,br,bi,mdb,ndb, */
/*    +                 wshsec,lshsec,work,lwork,ierror) */

/*     given the vector spherical harmonic coefficients br and bi, precomputed */
/*     by subroutine vhaec for a vector field (v,w), subroutine divec */
/*     computes the divergence of the vector field in the scalar array dv. */
/*     dv(i,j) is the divergence at the colatitude */

/*            theta(i) = (i-1)*pi/(nlat-1) */

/*     and east longitude */

/*            lambda(j) = (j-1)*2*pi/nlon */

/*     on the sphere.  i.e. */

/*            dv(i,j) = 1/sint*[ d(sint*v(i,j))/dtheta + d(w(i,j))/dlambda ] */

/*     where sint = sin(theta(i)).  w is the east longitudinal and v */
/*     is the colatitudinal component of the vector field from which */
/*     br,bi were precomputed.  required associated legendre polynomials */
/*     are recomputed rather than stored as they are in subroutine dives. */


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
/*            nlon = 72 for a five degree grid. nlon must be greater than */
/*            3.  the efficiency of the computation is improved when nlon */
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
/*            can be three dimensional corresponding to an indexed multiple */
/*            vector field.  in this case multiple scalar synthesis will */
/*            be performed to compute the divergence for each field.  the */
/*            third index is the synthesis index which assumes the values */
/*            k=1,...,nt.  for a single synthesis set nt = 1.  the */
/*            description of the remaining parameters is simplified by */
/*            assuming that nt=1 or that all the arrays are two dimensional. */

/*     idv    the first dimension of the array dv as it appears in */
/*            the program that calls divec. if isym = 0 then idv */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then idv must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then idv must be at least (nlat+1)/2. */

/*     jdv    the second dimension of the array dv as it appears in */
/*            the program that calls divec. jdv must be at least nlon. */

/*     br,bi  two or three dimensional arrays (see input parameter nt) */
/*            that contain vector spherical harmonic coefficients */
/*            of the vector field (v,w) as computed by subroutine vhaec. */
/*     ***    br and bi must be computed by vhaec prior to calling */
/*            divec. */

/*     mdb    the first dimension of the arrays br and bi as it */
/*            appears in the program that calls divec. mdb must be at */
/*            least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndb    the second dimension of the arrays br and bi as it */
/*            appears in the program that calls divec. ndb must be at */
/*            least nlat. */


/*     wshsec an array which must be initialized by subroutine shseci. */
/*            once initialized, */
/*            wshsec can be used repeatedly by divec as long as nlon */
/*            and nlat remain unchanged.  wshsec must not be altered */
/*            between calls of divec. */


/*     lshsec the dimension of the array wshsec as it appears in the */
/*            program that calls divec. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshsec must be at least */

/*            2*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls divec. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */

/*            if isym is zero then lwork must be at least */

/*               nlat*(nt*nlon+max0(3*l2,nlon)+2*nt*l1+1) */

/*            if isym is not zero then lwork must be at least */

/*               l2*(nt*nlon+max0(3*nlat,nlon)) + nlat*(2*nt*l1+1) */


/*     ************************************************************** */

/*     output parameters */


/*    dv     a two or three dimensional array (see input parameter nt) */
/*           that contains the divergence of the vector field (v,w) */
/*           whose coefficients br,bi where computed by subroutine */
/*           vhaec.  dv(i,j) is the divergence at the colatitude point */
/*           theta(i) = (i-1)*pi/(nlat-1) and longitude point */
/*           lambda(j) = (j-1)*2*pi/nlon. the index ranges are defined */
/*           above at the input parameter isym. */


/*    ierror = 0  no errors */
/*           = 1  error in the specification of nlat */
/*           = 2  error in the specification of nlon */
/*           = 3  error in the specification of isym */
/*           = 4  error in the specification of nt */
/*           = 5  error in the specification of idv */
/*           = 6  error in the specification of jdv */
/*           = 7  error in the specification of mdb */
/*           = 8  error in the specification of ndb */
/*           = 9  error in the specification of lshsec */
/*           = 10 error in the specification of lwork */
/* ********************************************************************** */

/* Subroutine */ int divec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *dv, integer *idv, integer *jdv, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *wshsec, integer *lshsec, doublereal *
	work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer dv_dim1, dv_dim2, dv_offset, br_dim1, br_dim2, br_offset, bi_dim1,
	     bi_dim2, bi_offset, i__1, i__2;

    /* Local variables */
    static integer l1, l2, ia, ib, mn, is, ls, mab, nln, iwk, lwk, lzz1, labc,
	     imid, mmax, lwmin;
    extern /* Subroutine */ int divec1_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, doublereal *, integer *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *);
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
    --wshsec;
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

/*     verify saved work space (same as shsec) */

    imid = (*nlat + 1) / 2;
    lzz1 = (*nlat << 1) * imid;
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) * 3 / 2;
    lwmin = lzz1 + labc + *nlon + 15;
    if (*lshsec < lwmin) {
	return 0;
    }

/*     verify unsaved work space (add to what shec requires) */

    *ierror = 10;
    ls = *nlat;
    if (*isym > 0) {
	ls = imid;
    }
    nln = *nt * ls * *nlon;

/*     set first dimension for a,b (as requried by shsec) */

/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mab = min(i__1,i__2);
    mn = mab * *nlat * *nt;
/*     if(lwork .lt. nln+max0(ls*nlon,3*nlat*imid)+2*mn+nlat) return */
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
    divec1_(nlat, nlon, isym, nt, &dv[dv_offset], idv, jdv, &br[br_offset], &
	    bi[bi_offset], mdb, ndb, &work[ia], &work[ib], &mab, &work[is], &
	    wshsec[1], lshsec, &work[iwk], &lwk, ierror);
    return 0;
} /* divec_ */

/* Subroutine */ int divec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *dv, integer *idv, integer *jdv, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *a, doublereal *b, integer *mab, doublereal *
	sqnn, doublereal *wshsec, integer *lshsec, doublereal *wk, integer *lwk, integer *
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
    static doublereal fn;
    static integer mmax;
    extern /* Subroutine */ int shsec_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, doublereal *, integer *
	    , integer *, doublereal *, integer *, doublereal *, integer *, integer *);


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
    --wshsec;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 2; n <= i__1; ++n) {
	fn = (doublereal) (n - 1);
	sqnn[n] = sqrt(fn * (fn + 1.));
/* L1: */
    }

/*     compute divergence scalar coefficients for each vector field */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlat;
	for (n = 1; n <= i__2; ++n) {
	    i__3 = *mab;
	    for (m = 1; m <= i__3; ++m) {
		a[m + (n + k * a_dim2) * a_dim1] = 0.;
		b[m + (n + k * b_dim2) * b_dim1] = 0.;
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

    shsec_(nlat, nlon, isym, nt, &dv[dv_offset], idv, jdv, &a[a_offset], &b[
	    b_offset], mab, nlat, &wshsec[1], lshsec, &wk[1], lwk, ierror);
    return 0;
} /* divec1_ */

