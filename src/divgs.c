/* divgs.f -- translated by f2c (version 20061008).
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




/* ... file divgs.f */

/*     this file includes documentation and code for */
/*     subroutine divgs          i */

/* ... files which must be loaded with divgs.f */

/*     sphcom.f, hrfft.f, vhags.f, shsgs.f, gaqd.f */


/*     subroutine divgs(nlat,nlon,isym,nt,divg,idiv,jdiv,br,bi,mdb,ndb, */
/*    +                 wshsgs,lshsgs,work,lwork,ierror) */

/*     given the vector spherical harmonic coefficients br and bi, precomputed */
/*     by subroutine vhags for a vector field (v,w), subroutine divgs */
/*     computes the divergence of the vector field in the scalar array divg. */
/*     divg(i,j) is the divergence at the gaussian colatitude point theta(i) */
/*     (see nlat as input parameter) and east longitude */

/*            lambda(j) = (j-1)*2*pi/nlon */

/*     on the sphere.  i.e. */

/*            dv(i,j) = 1/sint*[ d(sint*v(i,j))/dtheta + d(w(i,j))/dlambda ] */

/*     where sint = sin(theta(i)).  w is the east longitudinal and v */
/*     is the colatitudinal component of the vector field from which */
/*     br,bi were precomputed */


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
/*            sphere.  i.e., in the array divg(i,j) for i=1,...,nlat and */
/*            j=1,...,nlon. */

/*      = 1 */

/*            w is antisymmetric and v is symmetric about the equator. */
/*            in this case the divergence is antisymmetyric about */
/*            the equator and is computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the divergence is computed */
/*            in the array divg(i,j) for i=1,...,(nlat+1)/2 and for */
/*            j=1,...,nlon.  if nlat is even the divergence is computed */
/*            in the array divg(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*      = 2 */
/*            w is symmetric and v is antisymmetric about the equator */
/*            in this case the divergence is symmetyric about the */
/*            equator and is computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the divergence is computed */
/*            in the array divg(i,j) for i=1,...,(nlat+1)/2 and for */
/*            j=1,...,nlon.  if nlat is even the divergence is computed */
/*            in the array divg(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */


/*     nt     nt is the number of scalar and vector fields.  some */
/*            computational efficiency is obtained for multiple fields. */
/*            in the program that calls divgs, the arrays br,bi, and divg */
/*            can be three dimensional corresponding to an indexed multiple */
/*            vector field.  in this case multiple scalar synthesis will */
/*            be performed to compute the divergence for each field.  the */
/*            third index is the synthesis index which assumes the values */
/*            k=1,...,nt.  for a single synthesis set nt = 1.  the */
/*            description of the remaining parameters is simplified by */
/*            assuming that nt=1 or that all the arrays are two dimensional. */

/*     idiv   the first dimension of the array divg as it appears in */
/*            the program that calls divgs. if isym = 0 then idiv */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then idiv must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then idiv must be at least (nlat+1)/2. */

/*     jdiv   the second dimension of the array divg as it appears in */
/*            the program that calls divgs. jdiv must be at least nlon. */

/*     br,bi  two or three dimensional arrays (see input parameter nt) */
/*            that contain vector spherical harmonic coefficients */
/*            of the vector field (v,w) as computed by subroutine vhags. */
/*     ***    br and bi must be computed by vhags prior to calling */
/*            divgs. */

/*     mdb    the first dimension of the arrays br and bi as it */
/*            appears in the program that calls divgs. mdb must be at */
/*            least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndb    the second dimension of the arrays br and bi as it */
/*            appears in the program that calls divgs. ndb must be at */
/*            least nlat. */


/*     wshsgs an array which must be intialized by subroutine shsgsi. */
/*            once initialized, */
/*            wshsgs can be used repeatedly by divgs as long as nlon */
/*            and nlat remain unchanged.  wshsgs must not be altered */
/*            between calls of divgs. */


/*     lshsgs the dimension of the array wshsgs as it appears in the */
/*            program that calls divgs. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshsgs must be at least */

/*               nlat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*nlat-l1)-3*l1)/2+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls divgs. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */

/*            if isym = 0 then lwork must be at least */

/*               nlat*((nt+1)*nlon+2*nt*l1+1) */

/*            if isym > 0 then lwork must be at least */

/*               (nt+1)*l2*nlon+nlat*(2*nt*l1+1) */

/*     ************************************************************** */

/*     output parameters */


/*    divg   a two or three dimensional array (see input parameter nt) */
/*           that contains the divergence of the vector field (v,w) */
/*           whose coefficients br,bi where computed by subroutine */
/*           vhags.  divg(i,j) is the divergence at the gaussian colatitude */
/*           point theta(i) and longitude point lambda(j) = (j-1)*2*pi/nlon. */
/*           the index ranges are defined above at the input parameter */
/*           isym. */


/*    ierror = 0  no errors */
/*           = 1  error in the specification of nlat */
/*           = 2  error in the specification of nlon */
/*           = 3  error in the specification of isym */
/*           = 4  error in the specification of nt */
/*           = 5  error in the specification of idiv */
/*           = 6  error in the specification of jdiv */
/*           = 7  error in the specification of mdb */
/*           = 8  error in the specification of ndb */
/*           = 9  error in the specification of lshsgs */
/*           = 10 error in the specification of lwork */
/* ********************************************************************** */


/* Subroutine */ int divgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *divg, integer *idiv, integer *jdiv, 
	doublereal *br, doublereal *bi, integer *mdb, integer *ndb, 
	doublereal *wshsgs, integer *lshsgs, doublereal *work, integer *lwork,
	 integer *ierror)
{
    /* System generated locals */
    integer divg_dim1, divg_dim2, divg_offset, br_dim1, br_dim2, br_offset, 
	    bi_dim1, bi_dim2, bi_offset, i__1, i__2;

    /* Local variables */
    integer l1, l2, ia, ib, mn, lp, is, ls, mab, nln, iwk, lwk, imid, mmax, 
	    lpimn;
    extern /* Subroutine */ int divgs1_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);


/*     check input parameters */

    /* Parameter adjustments */
    divg_dim1 = *idiv;
    divg_dim2 = *jdiv;
    divg_offset = 1 + divg_dim1 * (1 + divg_dim2);
    divg -= divg_offset;
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
    if (*isym == 0 && *idiv < *nlat || *isym > 0 && *idiv < imid) {
	return 0;
    }
    *ierror = 6;
    if (*jdiv < *nlon) {
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
    lp = *nlat * ((l1 + l2) * 3 - 2) + (l1 - 1) * (l2 * ((*nlat << 1) - l1) - 
	    l1 * 3) / 2 + *nlon + 15;
    if (*lshsgs < lp) {
	return 0;
    }
    *ierror = 10;

/*     verify unsaved work space (add to what shses requires, file f3) */

    ls = *nlat;
    if (*isym > 0) {
	ls = imid;
    }
    nln = *nt * ls * *nlon;

/*     set first dimension for a,b (as requried by shses) */

/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mab = min(i__1,i__2);
    mn = mab * *nlat * *nt;
    if (*lwork < nln + ls * *nlon + (mn << 1) + *nlat) {
	return 0;
    }
    *ierror = 0;

/*     set work space pointers */

    ia = 1;
    ib = ia + mn;
    is = ib + mn;
    iwk = is + *nlat;
    lwk = *lwork - (mn << 1) - *nlat;
    divgs1_(nlat, nlon, isym, nt, &divg[divg_offset], idiv, jdiv, &br[
	    br_offset], &bi[bi_offset], mdb, ndb, &work[ia], &work[ib], &mab, 
	    &work[is], &wshsgs[1], lshsgs, &work[iwk], &lwk, ierror);
    return 0;
} /* divgs_ */

/* Subroutine */ int divgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *divg, integer *idiv, integer *jdiv, 
	doublereal *br, doublereal *bi, integer *mdb, integer *ndb, 
	doublereal *a, doublereal *b, integer *mab, doublereal *sqnn, 
	doublereal *wshsgs, integer *lshsgs, doublereal *wk, integer *lwk, 
	integer *ierror)
{
    /* System generated locals */
    integer divg_dim1, divg_dim2, divg_offset, br_dim1, br_dim2, br_offset, 
	    bi_dim1, bi_dim2, bi_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer k, m, n;
    doublereal fn;
    integer mmax;
    extern /* Subroutine */ int shsgs_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);


/*     set coefficient multiplyers */

    /* Parameter adjustments */
    --sqnn;
    divg_dim1 = *idiv;
    divg_dim2 = *jdiv;
    divg_offset = 1 + divg_dim1 * (1 + divg_dim2);
    divg -= divg_offset;
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

/*     synthesize a,b into divg */

    shsgs_(nlat, nlon, isym, nt, &divg[divg_offset], idiv, jdiv, &a[a_offset],
	     &b[b_offset], mab, nlat, &wshsgs[1], lshsgs, &wk[1], lwk, ierror)
	    ;
    return 0;
} /* divgs1_ */

