/* igrades.f -- translated by f2c (version 20061008).
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



/* ... file igrades.f */

/*     this file includes documentation and code for */
/*     subroutine igrades         i */

/* ... files which must be loaded with igradec.f */

/*     sphcom.f, hrfft.f, shses.f,vhaes.f */

/*     subroutine igrades(nlat,nlon,isym,nt,sf,isf,jsf,br,bi,mdb,ndb, */
/*    +                   wshses,lshses,work,lwork,ierror) */

/*     let br,bi,cr,ci be the vector spherical harmonic coefficients */
/*     precomputed by vhaes for a vector field (v,w).  let (v',w') be */
/*     the irrotational component of (v,w) (i.e., (v',w') is generated */
/*     by assuming cr,ci are zero and synthesizing br,bi with vhses). */
/*     then subroutine igrades computes a scalar field sf such that */

/*            gradient(sf) = (v',w'). */

/*     i.e., */

/*            v'(i,j) = d(sf(i,j))/dtheta          (colatitudinal component of */
/*                                                 the gradient) */
/*     and */

/*            w'(i,j) = 1/sint*d(sf(i,j))/dlambda  (east longitudinal component */
/*                                                 of the gradient) */

/*     at colatitude */

/*            theta(i) = (i-1)*pi/(nlat-1) */

/*     and longitude */

/*            lambda(j) = (j-1)*2*pi/nlon */

/*     where sint = sin(theta(i)).  required associated legendre polynomials */
/*     are stored rather than recomputed as they are in subroutine igradec. */

/*     note:  for an irrotational vector field (v,w), subroutine igrades */
/*     computes a scalar field whose gradient is (v,w).  in ay case, */
/*     subroutine igrades "inverts" the gradient subroutine grades. */


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


/*     isym   a parameter which determines whether the scalar field sf is */
/*            computed on the full or half sphere as follows: */

/*      = 0 */

/*            the symmetries/antsymmetries described in isym=1,2 below */
/*            do not exist in (v,w) about the equator.  in this case sf */
/*            is neither symmetric nor antisymmetric about the equator. */
/*            sf is computed on the entire sphere.  i.e., in the array */
/*            sf(i,j) for i=1,...,nlat and  j=1,...,nlon */

/*      = 1 */

/*            w is antisymmetric and v is symmetric about the equator. */
/*            in this case sf is antisymmetyric about the equator and */
/*            is computed for the northern hemisphere only.  i.e., */
/*            if nlat is odd sf is computed in the array sf(i,j) for */
/*            i=1,...,(nlat+1)/2 and for j=1,...,nlon.  if nlat is even */
/*            sf is computed in the array sf(i,j) for i=1,...,nlat/2 */
/*            and j=1,...,nlon. */

/*      = 2 */

/*            w is symmetric and v is antisymmetric about the equator. */
/*            in this case sf is symmetyric about the equator and */
/*            is computed for the northern hemisphere only.  i.e., */
/*            if nlat is odd sf is computed in the array sf(i,j) for */
/*            i=1,...,(nlat+1)/2 and for j=1,...,nlon.  if nlat is even */
/*            sf is computed in the array sf(i,j) for i=1,...,nlat/2 */
/*            and j=1,...,nlon. */


/*     nt     nt is the number of scalar and vector fields.  some */
/*            computational efficiency is obtained for multiple fields. */
/*            the arrays br,bi, and sf can be three dimensional corresponding */
/*            to an indexed multiple vector field (v,w).  in this case, */
/*            multiple scalar synthesis will be performed to compute each */
/*            scalar field.  the third index for br,bi, and sf is the synthesis */
/*            index which assumes the values k = 1,...,nt.  for a single */
/*            synthesis set nt = 1.  the description of the remaining */
/*            parameters is simplified by assuming that nt=1 or that br,bi, */
/*            and sf are two dimensional arrays. */

/*     isf    the first dimension of the array sf as it appears in */
/*            the program that calls igrades. if isym = 0 then isf */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then isf must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then isf must be at least (nlat+1)/2. */

/*     jsf    the second dimension of the array sf as it appears in */
/*            the program that calls igrades. jsf must be at least nlon. */

/*     br,bi  two or three dimensional arrays (see input parameter nt) */
/*            that contain vector spherical harmonic coefficients */
/*            of the vector field (v,w) as computed by subroutine vhaes. */
/*     ***    br,bi must be computed by vhaes prior to calling igrades. */

/*     mdb    the first dimension of the arrays br and bi as it appears in */
/*            the program that calls igrades (and vhaes). mdb must be at */
/*            least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndb    the second dimension of the arrays br and bi as it appears in */
/*            the program that calls igrades (and vhaes). ndb must be at */
/*            least nlat. */


/*  wshses    an array which must be initialized by subroutine igradesi */
/*            (or equivalently by subroutine shsesi).  once initialized, */
/*            wshses can be used repeatedly by igrades as long as nlon */
/*            and nlat remain unchanged.  wshses must not be altered */
/*            between calls of igrades. */


/*  lshses    the dimension of the array wshses as it appears in the */
/*            program that calls igrades. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd. */


/*            then lshses must be greater than or equal to */

/*               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15 */

/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls igrades. define */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */
/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            if isym = 0 lwork must be greater than or equal to */

/*               nlat*((nt+1)*nlon+2*nt*l1+1) */

/*            if isym > 0 lwork must be greater than or equal to */

/*               (nt+1)*l2*nlon+nlat*(2*nt*l1+1) */



/*     ************************************************************** */

/*     output parameters */


/*     sf    a two or three dimensional array (see input parameter nt) that */
/*           contain a scalar field whose gradient is the irrotational */
/*           component of the vector field (v,w).  the vector spherical */
/*           harmonic coefficients br,bi were precomputed by subroutine */
/*           vhaes.  sf(i,j) is given at the gaussian colatitude theta(i) */
/*           and longitude lambda(j) = (j-1)*2*pi/nlon.  the index ranges */
/*           are defined at input parameter isym. */


/*  ierror   = 0  no errors */
/*           = 1  error in the specification of nlat */
/*           = 2  error in the specification of nlon */
/*           = 3  error in the specification of isym */
/*           = 4  error in the specification of nt */
/*           = 5  error in the specification of isf */
/*           = 6  error in the specification of jsf */
/*           = 7  error in the specification of mdb */
/*           = 8  error in the specification of ndb */
/*           = 9  error in the specification of lshses */
/*           = 10 error in the specification of lwork */

/* ********************************************************************** */

/* Subroutine */ int igrades_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, real *sf, integer *isf, integer *jsf, real *br, real *bi,
	 integer *mdb, integer *ndb, real *wshses, integer *lshses, real *
	work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer sf_dim1, sf_dim2, sf_offset, br_dim1, br_dim2, br_offset, bi_dim1,
	     bi_dim2, bi_offset, i__1, i__2;

    /* Local variables */
    static integer ia, ib, mn, is, ls, mab, nln, iwk, imid, mmax, liwk, lpimn,
	     lwkmin;
    extern /* Subroutine */ int igrdes1_(integer *, integer *, integer *, 
	    integer *, real *, integer *, integer *, real *, real *, integer *
	    , real *, integer *, integer *, real *, real *, real *, integer *,
	     real *, integer *, integer *);


/*     check input parameters */

    /* Parameter adjustments */
    sf_dim1 = *isf;
    sf_dim2 = *jsf;
    sf_offset = 1 + sf_dim1 * (1 + sf_dim2);
    sf -= sf_offset;
    bi_dim1 = *mdb;
    bi_dim2 = *ndb;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mdb;
    br_dim2 = *ndb;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
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
    imid = (*nlat + 1) / 2;
    if (*isym == 0 && *isf < *nlat || *isym != 0 && *isf < imid) {
	return 0;
    }
    *ierror = 6;
    if (*jsf < *nlon) {
	return 0;
    }
    *ierror = 7;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    mmax = min(i__1,i__2);
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

/*     verify saved work space length */

    imid = (*nlat + 1) / 2;
    lpimn = imid * mmax * (*nlat + *nlat - mmax + 1) / 2;
    if (*lshses < lpimn + *nlon + 15) {
	return 0;
    }
    *ierror = 10;

/*     set minimum and verify unsaved work space */

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
    lwkmin = nln + ls * *nlon + (mn << 1) + *nlat;
    if (*lwork < lwkmin) {
	return 0;
    }
    *ierror = 0;

/*     set work space pointers */

    ia = 1;
    ib = ia + mn;
    is = ib + mn;
    iwk = is + *nlat;
    liwk = *lwork - (mn << 1) - *nlat;
    igrdes1_(nlat, nlon, isym, nt, &sf[sf_offset], isf, jsf, &work[ia], &work[
	    ib], &mab, &work[is], mdb, ndb, &br[br_offset], &bi[bi_offset], &
	    wshses[1], lshses, &work[iwk], &liwk, ierror);
    return 0;
} /* igrades_ */

/* Subroutine */ int igrdes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, real *sf, integer *isf, integer *jsf, real *a, real *b, 
	integer *mab, real *sqnn, integer *mdb, integer *ndb, real *br, real *
	bi, real *wshses, integer *lshses, real *wk, integer *lwk, integer *
	ierror)
{
    /* System generated locals */
    integer sf_dim1, sf_dim2, sf_offset, br_dim1, br_dim2, br_offset, bi_dim1,
	     bi_dim2, bi_offset, a_dim1, a_dim2, a_offset, b_dim1, b_dim2, 
	    b_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, m, n;
    static real fn;
    static integer mmax;
    extern /* Subroutine */ int shses_(integer *, integer *, integer *, 
	    integer *, real *, integer *, integer *, real *, real *, integer *
	    , integer *, real *, integer *, real *, integer *, integer *);


/*     preset coefficient multiplyers in vector */

    /* Parameter adjustments */
    --sqnn;
    sf_dim1 = *isf;
    sf_dim2 = *jsf;
    sf_offset = 1 + sf_dim1 * (1 + sf_dim2);
    sf -= sf_offset;
    b_dim1 = *mab;
    b_dim2 = *nlat;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mab;
    a_dim2 = *nlat;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    bi_dim1 = *mdb;
    bi_dim2 = *ndb;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mdb;
    br_dim2 = *ndb;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
    --wshses;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 2; n <= i__1; ++n) {
	fn = (real) (n - 1);
	sqnn[n] = 1.f / sqrt(fn * (fn + 1.f));
/* L1: */
    }

/*     set upper limit for vector m subscript */

/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);

/*     compute multiple scalar field coefficients */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {

/*     preset to 0.0 */

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
	    a[(n + k * a_dim2) * a_dim1 + 1] = br[(n + k * br_dim2) * br_dim1 
		    + 1] * sqnn[n];
	    b[(n + k * b_dim2) * b_dim1 + 1] = bi[(n + k * bi_dim2) * bi_dim1 
		    + 1] * sqnn[n];
/* L5: */
	}

/*     compute m>0 coefficients */

	i__2 = mmax;
	for (m = 2; m <= i__2; ++m) {
	    i__3 = *nlat;
	    for (n = m; n <= i__3; ++n) {
		a[m + (n + k * a_dim2) * a_dim1] = sqnn[n] * br[m + (n + k * 
			br_dim2) * br_dim1];
		b[m + (n + k * b_dim2) * b_dim1] = sqnn[n] * bi[m + (n + k * 
			bi_dim2) * bi_dim1];
/* L7: */
	    }
/* L6: */
	}
/* L2: */
    }

/*     scalar sythesize a,b into sf */

    shses_(nlat, nlon, isym, nt, &sf[sf_offset], isf, jsf, &a[a_offset], &b[
	    b_offset], mab, nlat, &wshses[1], lshses, &wk[1], lwk, ierror);
    return 0;
} /* igrdes1_ */

