/* islapes.f -- translated by f2c (version 20061008).
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


/* ... file islapes.f */

/*     this file includes documentation and code for */
/*     subroutine islapes         i */

/* ... files which must be loaded with islapes.f */

/*     sphcom.f, hrfft.f, shaes.f, shses.f */

/*     subroutine islapes(nlat,nlon,isym,nt,xlmbda,sf,ids,jds,a,b, */
/*    +mdab,ndab,wshses,lshses,work,lwork,pertrb,ierror) */

/*     islapes inverts the laplace or helmholz operator on an equally */
/*     spaced latitudinal grid using o(n**3) storage. given the */
/*     spherical harmonic coefficients a(m,n) and b(m,n) of the right */
/*     hand side slap(i,j), islapes computes a solution sf(i,j) to */
/*     the following helmhotz equation : */

/*           2                2 */
/*     [d(sf(i,j))/dlambda /sint + d(sint*d(sf(i,j))/dtheta)/dtheta]/sint */

/*                   - xlmbda * sf(i,j) = slap(i,j) */

/*      where sf(i,j) is computed at colatitude */

/*                 theta(i) = (i-1)*pi/(nlat-1) */

/*            and longitude */

/*                 lambda(j) = (j-1)*2*pi/nlon */

/*            for i=1,...,nlat and j=1,...,nlon. */


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

/*     nlon   the number of distinct longitude points.  nlon determines */
/*            the grid increment in longitude as 2*pi/nlon. for example */
/*            nlon = 72 for a five degree grid. nlon must be greater */
/*            than zero. the axisymmetric case corresponds to nlon=1. */
/*            the efficiency of the computation is improved when nlon */
/*            is a product of small prime numbers. */

/*     isym   this parameter should have the same value input to subroutine */
/*            shaes to compute the coefficients a and b for the scalar field */
/*            slap.  isym is set as follows: */

/*            = 0  no symmetries exist in slap about the equator. scalar */
/*                 synthesis is used to compute sf on the entire sphere. */
/*                 i.e., in the array sf(i,j) for i=1,...,nlat and */
/*                 j=1,...,nlon. */

/*           = 1  sf and slap are antisymmetric about the equator. the */
/*                synthesis used to compute sf is performed on the */
/*                northern hemisphere only.  if nlat is odd, sf(i,j) is */
/*                computed for i=1,...,(nlat+1)/2 and j=1,...,nlon.  if */
/*                nlat is even, sf(i,j) is computed for i=1,...,nlat/2 */
/*                and j=1,...,nlon. */


/*           = 2  sf and slap are symmetric about the equator. the */
/*                synthesis used to compute sf is performed on the */
/*                northern hemisphere only.  if nlat is odd, sf(i,j) is */
/*                computed for i=1,...,(nlat+1)/2 and j=1,...,nlon.  if */
/*                nlat is even, sf(i,j) is computed for i=1,...,nlat/2 */
/*                and j=1,...,nlon. */


/*   nt       the number of solutions. in the program that calls islapes */
/*            the arrays sf,a, and b can be three dimensional in which */
/*            case multiple solutions are computed. the third index */
/*            is the solution index with values k=1,...,nt. */
/*            for a single solution set nt=1. the description of the */
/*            remaining parameters is simplified by assuming that nt=1 */
/*            and sf,a,b are two dimensional. */

/*   xlmbda   a one dimensional array with nt elements. if xlmbda is */
/*            is identically zero islapes solves poisson's equation. */
/*            if xlmbda > 0.0 islapes solves the helmholtz equation. */
/*            if xlmbda < 0.0 the nonfatal error flag ierror=-1 is */
/*            returned. negative xlambda could result in a division */
/*            by zero. */

/*   ids      the first dimension of the array sf as it appears in the */
/*            program that calls islapes.  if isym = 0 then ids must be at */
/*            least nlat.  if isym > 0 and nlat is even then ids must be */
/*            at least nlat/2. if isym > 0 and nlat is odd then ids must */
/*            be at least (nlat+1)/2. */

/*   jds      the second dimension of the array sf as it appears in the */
/*            program that calls islapes. jds must be at least nlon. */


/*   a,b      two or three dimensional arrays (see input parameter nt) */
/*            that contain scalar spherical harmonic coefficients */
/*            of the scalar field slap. a,b must be computed by shaes */
/*            prior to calling islapes. */


/*   mdab     the first dimension of the arrays a and b as it appears */
/*            in the program that calls islapes.  mdab must be at */
/*            least min0(nlat,(nlon+2)/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*   ndab     the second dimension of the arrays a and b as it appears */
/*            in the program that calls islapes. ndbc must be at least */
/*            least nlat. */

/*            mdab,ndab should have the same values input to shaes to */
/*            compute the coefficients a and b. */


/*   wshses   an array which must be initialized by subroutine shsesi. */
/*            once initialized, wshses can be used repeatedly by */
/*            islapes as long as nlat and nlon  remain unchanged. */
/*            wshses must not be altered between calls of islapes. */

/*    lshses  the dimension of the array wshses as it appears in the */
/*            program that calls islapes.  let */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshses must be at least */

/*               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15 */

/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls islapes. define */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */
/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            if isym is zero then lwork must be at least */

/*               (nt+1)*nlat*nlon + nlat*(2*nt*l1+1) */

/*            if isym is nonzero lwork must be at least */

/*               (nt+1)*l2*nlon + nlat*(2*nt*l1+1) */


/*     ************************************************************** */

/*     output parameters */


/*    sf      a two or three dimensional arrays (see input parameter nt) that */
/*            inverts the scalar laplacian in slap - pertrb.  sf(i,j) is given */
/*            at the colatitude */

/*                 theta(i) = (i-1)*pi/(nlat-1) */

/*            and longitude */

/*                 lambda(j) = (j-1)*2*pi/nlon */

/*            for i=1,...,nlat and j=1,...,nlon. */

/*   pertrb  a one dimensional array with nt elements (see input */
/*           parameter nt). in the discription that follows we assume */
/*           that nt=1. if xlmbda > 0.0 then pertrb=0.0 is always */
/*           returned because the helmholtz operator is invertible. */
/*           if xlmbda = 0.0 then a solution exists only if a(1,1) */
/*           is zero. islapec sets a(1,1) to zero. the resulting */
/*           solution sf(i,j) solves poisson's equation with */
/*           pertrb = a(1,1)/(2.*sqrt(2.)) subtracted from the */
/*           right side slap(i,j). */

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

/*            = 9  error in the specification of lshses */

/*            = 10 error in the specification of lwork */


/* ********************************************************************** */

/*     end of documentation for islapes */

/* ********************************************************************** */

/* Subroutine */ int islapes_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, real *xlmbda, real *sf, integer *ids, integer *jds, real 
	*a, real *b, integer *mdab, integer *ndab, real *wshses, integer *
	lshses, real *work, integer *lwork, real *pertrb, integer *ierror)
{
    /* System generated locals */
    integer sf_dim1, sf_dim2, sf_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    static integer k, l1, l2, ia, ib, mn, ls, ifn, nln, iwk, lwk, imid, mmax, 
	    lpimn, lwkmin;
    extern /* Subroutine */ int islpes1_(integer *, integer *, integer *, 
	    integer *, real *, real *, integer *, integer *, real *, real *, 
	    integer *, integer *, real *, real *, integer *, real *, real *, 
	    integer *, real *, integer *, real *, integer *);


/*     check input parameters */

    /* Parameter adjustments */
    --pertrb;
    --xlmbda;
    sf_dim1 = *ids;
    sf_dim2 = *jds;
    sf_offset = 1 + sf_dim1 * (1 + sf_dim2);
    sf -= sf_offset;
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

    imid = (*nlat + 1) / 2;
    lpimn = imid * mmax * (*nlat + *nlat - mmax + 1) / 2;
    if (*lshses < lpimn + *nlon + 15) {
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
/*     lwkmin = nln+ls*nlon+2*mn+nlat */
/*     if (lwork .lt. lwkmin) return */
    l2 = (*nlat + 1) / 2;
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    l1 = min(i__1,i__2);
    if (*isym == 0) {
	lwkmin = (*nt + 1) * *nlat * *nlon + *nlat * ((*nt << 1) * l1 + 1);
    } else {
	lwkmin = (*nt + 1) * l2 * *nlon + *nlat * ((*nt << 1) * l1 + 1);
    }
    if (*lwork < lwkmin) {
	return 0;
    }
    *ierror = 0;

/*     check sign of xlmbda */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	if (xlmbda[k] < 0.f) {
	    *ierror = -1;
	}
    }

/*     set work space pointers */

    ia = 1;
    ib = ia + mn;
    ifn = ib + mn;
    iwk = ifn + *nlat;
    lwk = *lwork - (mn << 1) - *nlat;
    islpes1_(nlat, nlon, isym, nt, &xlmbda[1], &sf[sf_offset], ids, jds, &a[
	    a_offset], &b[b_offset], mdab, ndab, &work[ia], &work[ib], &mmax, 
	    &work[ifn], &wshses[1], lshses, &work[iwk], &lwk, &pertrb[1], 
	    ierror);
    return 0;
} /* islapes_ */

/* Subroutine */ int islpes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, real *xlmbda, real *sf, integer *ids, integer *jds, real 
	*a, real *b, integer *mdab, integer *ndab, real *as, real *bs, 
	integer *mmax, real *fnn, real *wshses, integer *lshses, real *wk, 
	integer *lwk, real *pertrb, integer *ierror)
{
    /* System generated locals */
    integer sf_dim1, sf_dim2, sf_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, as_dim1, as_dim2, as_offset, bs_dim1, bs_dim2, 
	    bs_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer k, m, n;
    static real fn;
    extern /* Subroutine */ int shses_(integer *, integer *, integer *, 
	    integer *, real *, integer *, integer *, real *, real *, integer *
	    , integer *, real *, integer *, real *, integer *, integer *);


/*     set multipliers and preset synthesis coefficients to zero */

    /* Parameter adjustments */
    --fnn;
    --pertrb;
    --xlmbda;
    sf_dim1 = *ids;
    sf_dim2 = *jds;
    sf_offset = 1 + sf_dim1 * (1 + sf_dim2);
    sf -= sf_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    bs_dim1 = *mmax;
    bs_dim2 = *nlat;
    bs_offset = 1 + bs_dim1 * (1 + bs_dim2);
    bs -= bs_offset;
    as_dim1 = *mmax;
    as_dim2 = *nlat;
    as_offset = 1 + as_dim1 * (1 + as_dim2);
    as -= as_offset;
    --wshses;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 1; n <= i__1; ++n) {
	fn = (real) (n - 1);
	fnn[n] = fn * (fn + 1.f);
	i__2 = *mmax;
	for (m = 1; m <= i__2; ++m) {
	    i__3 = *nt;
	    for (k = 1; k <= i__3; ++k) {
		as[m + (n + k * as_dim2) * as_dim1] = 0.f;
		bs[m + (n + k * bs_dim2) * bs_dim1] = 0.f;
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {

/*     compute synthesis coefficients for xlmbda zero or nonzero */

	if (xlmbda[k] == 0.f) {
	    i__2 = *nlat;
	    for (n = 2; n <= i__2; ++n) {
		as[(n + k * as_dim2) * as_dim1 + 1] = -a[(n + k * a_dim2) * 
			a_dim1 + 1] / fnn[n];
		bs[(n + k * bs_dim2) * bs_dim1 + 1] = -b[(n + k * b_dim2) * 
			b_dim1 + 1] / fnn[n];
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    as[m + (n + k * as_dim2) * as_dim1] = -a[m + (n + k * 
			    a_dim2) * a_dim1] / fnn[n];
		    bs[m + (n + k * bs_dim2) * bs_dim1] = -b[m + (n + k * 
			    b_dim2) * b_dim1] / fnn[n];
		}
	    }
	} else {

/*     xlmbda nonzero so operator invertible unless */
/*     -n*(n-1) = xlmbda(k) < 0.0  for some n */

	    pertrb[k] = 0.f;
	    i__2 = *nlat;
	    for (n = 1; n <= i__2; ++n) {
		as[(n + k * as_dim2) * as_dim1 + 1] = -a[(n + k * a_dim2) * 
			a_dim1 + 1] / (fnn[n] + xlmbda[k]);
		bs[(n + k * bs_dim2) * bs_dim1 + 1] = -b[(n + k * b_dim2) * 
			b_dim1 + 1] / (fnn[n] + xlmbda[k]);
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    as[m + (n + k * as_dim2) * as_dim1] = -a[m + (n + k * 
			    a_dim2) * a_dim1] / (fnn[n] + xlmbda[k]);
		    bs[m + (n + k * bs_dim2) * bs_dim1] = -b[m + (n + k * 
			    b_dim2) * b_dim1] / (fnn[n] + xlmbda[k]);
		}
	    }
	}
    }

/*     synthesize as,bs into sf */

    shses_(nlat, nlon, isym, nt, &sf[sf_offset], ids, jds, &as[as_offset], &
	    bs[bs_offset], mmax, nlat, &wshses[1], lshses, &wk[1], lwk, 
	    ierror);
    return 0;
} /* islpes1_ */

