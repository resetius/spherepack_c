/* isfvpec.f -- translated by f2c (version 20061008).
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




/* ... file isfvpec.f */

/*     this file includes documentation and code for */
/*     subroutine isfvpec          i */

/* ... files which must be loaded with isfvpec.f */

/*     sphcom.f, hrfft.f, vhsec.f,shaec.f */


/*     subroutine isfvpec(nlat,nlon,isym,nt,sf,vp,idv,jdv,as,bs,av,bv, */
/*    +                   mdb,ndb,wvhsec,lvhsec,work,lwork,ierror) */

/*     given the scalar spherical harmonic coefficients as,bs precomputed */
/*     by shaec for the scalar stream function sf and av,bv precomputed by */
/*     shaec for the scalar velocity potenital vp, subroutine isfvpec computes */
/*     the vector field (v,w) corresponding to sf and vp.  w is the east */
/*     longitudinal and v is the colatitudinal component of the vector field. */
/*     (v,w) is expressed in terms of sf,vp by the helmholtz relations (in */
/*     mathematical spherical coordinates): */

/*          v = -1/sin(theta)*d(vp)/dlambda + d(st)/dtheta */

/*          w =  1/sin(theta)*d(st)/dlambda + d(vp)/dtheta */

/*     required legendre functions are recomputed rather than stored as */
/*     they are in subroutine isfvpes.  v(i,j) and w(i,j) are given at */
/*     colatitude */

/*            theta(i) = (i-1)*pi/(nlat-1) */

/*     and east longitude */

/*            lambda(j) = (j-1)*2*pi/nlon */

/*     on the sphere (pi=4.0*atan(1.0)). */


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


/*     isym   a parameter which determines whether the vector field is */
/*            computed on the full or half sphere as follows: */

/*      = 0 */

/*            the symmetries/antsymmetries described in isym=1,2 below */
/*            do not exist in sf,vp about the equator.  in this case v */
/*            and w are not necessarily symmetric or antisymmetric about */
/*            equator.  v and w are computed on the entire sphere. */
/*            i.e., in arrays sf(i,j),vp(i,j) for i=1,...,nlat and */
/*            j=1,...,nlon. */

/*      = 1 */

/*            vp is antisymmetric and sf is symmetric about the equator. */
/*            in this case v is symmetric and w antisymmetric about */
/*            the equator and are computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the v(i,j),w(i,j) are computed */
/*            for i=1,...,(nlat+1)/2 and for j=1,...,nlon.  if nlat is */
/*            even then v(i,j),w(i,j) are computed for i=1,...,nlat/2 */
/*            and j=1,...,nlon. */

/*      = 2 */

/*            vp is symmetric and sf is antisymmetric about the equator. */
/*            in this case v is antisymmetric and w symmetric about */
/*            the equator and are computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the v(i,j),w(i,j) are computed */
/*            for i=1,...,(nlat+1)/2 and for j=1,...,nlon.  if nlat is */
/*            even then v(i,j),w(i,j) are computed for i=1,...,nlat/2 */
/*            and j=1,...,nlon. */

/*     nt     nt is the number of scalar and vector fields.  some */
/*            computational efficiency is obtained for multiple fields. arrays */
/*            can be three dimensional corresponding to an indexed multiple */
/*            vector field.  in this case multiple vector synthesis will */
/*            be performed to compute (v,w) for each field.  the */
/*            third index is the synthesis index which assumes the values */
/*            k=1,...,nt.  for a single synthesis set nt = 1.  the */
/*            description of the remaining parameters is simplified by */
/*            assuming that nt=1 or that all the arrays are two dimensional. */

/*     idv    the first dimension of the arrays v,w as it appears in */
/*            the program that calls isfvpec. if isym = 0 then idv */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then idv must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then idv must be at least (nlat+1)/2. */

/*     jdv    the second dimension of the arrays v,w as it appears in */
/*            the program that calls isfvpec. jdv must be at least nlon. */

/*     as,bs  two or three dimensional arrays (see input parameter nt) */
/*            that contain the spherical harmonic coefficients of */
/*            the scalar field sf as computed by subroutine shaec. */

/*     av,bv  two or three dimensional arrays (see input parameter nt) */
/*            that contain the spherical harmonic coefficients of */
/*            the scalar field vp as computed by subroutine shaec. */

/*     mdb    the first dimension of the arrays as,bs,av,bv as it */
/*            appears in the program that calls isfvpec. mdb must be at */
/*            least min0(nlat,(nlon+2)/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndb    the second dimension of the arrays as,bs,av,bv as it */
/*            appears in the program that calls isfvpec. ndb must be at */
/*            least nlat. */

/*     wvhsec an array which must be initialized by subroutine vhseci. */
/*            once initialized, wvhsec can be used repeatedly by isfvpec */
/*            as long as nlon and nlat remain unchanged.  wvhsec must */
/*            not bel altered between calls of isfvpec. */


/*     lvhsec the dimension of the array wvhsec as it appears in the */
/*            program that calls isfvpec. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lvhsec must be at least */


/*               4*nlat*l2+3*max0(l1-2,0)*(nlat+nlat-l1-1)+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls isfvpec. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */

/*            if isym = 0 then lwork must be at least */

/*               nlat*(2*nt*nlon+max0(6*l2,nlon)+4*l1*nt+1) */

/*            if isym = 1 or 2 then lwork must be at least */

/*               l2*(2*nt*nlon+max0(6*nlat,nlon))+nlat*(4*l1*nt+1) */

/*     ************************************************************** */

/*     output parameters */

/*    v,w    two or three dimensional arrays (see input parameter nt) */
/*           that contains the vector field corresponding to the stream */
/*           function sf and velocity potential vp whose coefficients, */
/*           as,bs (for sf) and av,bv (for vp), were precomputed by */
/*           subroutine shaec.  v(i,j) and w(i,j) are given at the */
/*           colatitude point */

/*                theta(i) = (i-1)*pi/(nlat-1) */

/*           and longitude point */

/*                lambda(j) = (j-1)*2*pi/nlon */

/*           the index ranges are defined above at the input parameter isym. */


/*    ierror = 0  no errors */
/*           = 1  error in the specification of nlat */
/*           = 2  error in the specification of nlon */
/*           = 3  error in the specification of isym */
/*           = 4  error in the specification of nt */
/*           = 5  error in the specification of idv */
/*           = 6  error in the specification of jdv */
/*           = 7  error in the specification of mdb */
/*           = 8  error in the specification of ndb */
/*           = 9  error in the specification of lvhsec */
/*           = 10 error in the specification of lwork */
/* ********************************************************************** */

/* Subroutine */ int isfvpec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idv, integer *jdv,
	 doublereal *as, doublereal *bs, doublereal *av, doublereal *bv, 
	integer *mdb, integer *ndb, doublereal *wvhsec, integer *lvhsec, 
	doublereal *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, as_dim1, 
	    as_dim2, as_offset, bs_dim1, bs_dim2, bs_offset, av_dim1, av_dim2,
	     av_offset, bv_dim1, bv_dim2, bv_offset, i__1, i__2;

    /* Local variables */
    static integer l1, l2, mn, is, ibi, ici, ibr, icr, iwk, lwk, lzz1, labc, 
	    mmax, lwmin;
    extern /* Subroutine */ int isfvpec1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);


/*     check input parameters */

    /* Parameter adjustments */
    w_dim1 = *idv;
    w_dim2 = *jdv;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    v_dim1 = *idv;
    v_dim2 = *jdv;
    v_offset = 1 + v_dim1 * (1 + v_dim2);
    v -= v_offset;
    bv_dim1 = *mdb;
    bv_dim2 = *ndb;
    bv_offset = 1 + bv_dim1 * (1 + bv_dim2);
    bv -= bv_offset;
    av_dim1 = *mdb;
    av_dim2 = *ndb;
    av_offset = 1 + av_dim1 * (1 + av_dim2);
    av -= av_offset;
    bs_dim1 = *mdb;
    bs_dim2 = *ndb;
    bs_offset = 1 + bs_dim1 * (1 + bs_dim2);
    bs -= bs_offset;
    as_dim1 = *mdb;
    as_dim2 = *ndb;
    as_offset = 1 + as_dim1 * (1 + as_dim2);
    as -= as_offset;
    --wvhsec;
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
    l2 = (*nlat + 1) / 2;
    if (*isym == 0 && *idv < *nlat || *isym > 0 && *idv < l2) {
	return 0;
    }
    *ierror = 6;
    if (*jdv < *nlon) {
	return 0;
    }
    *ierror = 7;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    l1 = min(i__1,i__2);
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    if (*mdb < min(i__1,i__2)) {
	return 0;
    }
    *ierror = 8;
    if (*ndb < *nlat) {
	return 0;
    }
    *ierror = 9;
    lzz1 = (*nlat << 1) * l2;
/* Computing MAX */
    i__1 = l1 - 2;
    labc = max(i__1,0) * (*nlat + *nlat - l1 - 1) * 3 / 2;
    if (*lvhsec < (lzz1 + labc << 1) + *nlon + 15) {
	return 0;
    }
    *ierror = 10;
    if (*isym == 0) {
/* Computing MAX */
	i__1 = l2 * 6;
	lwmin = *nlat * ((*nt << 1) * *nlon + max(i__1,*nlon) + (l1 << 2) * *
		nt + 1);
    } else {
/* Computing MAX */
	i__1 = *nlat * 6;
	lwmin = l2 * ((*nt << 1) * *nlon + max(i__1,*nlon)) + *nlat * ((l1 << 
		2) * *nt + 1);
    }
    if (*lwork < lwmin) {
	return 0;
    }

/*     set first dimension for br,bi,cr,ci (as requried by vhsec) */

/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
    mn = mmax * *nlat * *nt;
    *ierror = 0;

/*     set work space pointers */

    ibr = 1;
    ibi = ibr + mn;
    icr = ibi + mn;
    ici = icr + mn;
    is = ici + mn;
    iwk = is + *nlat;
    lwk = *lwork - (mn << 2) - *nlat;
    isfvpec1_(nlat, nlon, isym, nt, &v[v_offset], &w[w_offset], idv, jdv, &as[
	    as_offset], &bs[bs_offset], &av[av_offset], &bv[bv_offset], mdb, 
	    ndb, &work[ibr], &work[ibi], &work[icr], &work[ici], &l1, &work[
	    is], &wvhsec[1], lvhsec, &work[iwk], &lwk, ierror);
    return 0;
} /* isfvpec_ */

/* Subroutine */ int isfvpec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idv, integer *jdv,
	 doublereal *as, doublereal *bs, doublereal *av, doublereal *bv, 
	integer *mdb, integer *ndb, doublereal *br, doublereal *bi, 
	doublereal *cr, doublereal *ci, integer *mab, doublereal *fnn, 
	doublereal *wvhsec, integer *lvhsec, doublereal *wk, integer *lwk, 
	integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, as_dim1, 
	    as_dim2, as_offset, bs_dim1, bs_dim2, bs_offset, av_dim1, av_dim2,
	     av_offset, bv_dim1, bv_dim2, bv_offset, br_dim1, br_dim2, 
	    br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, cr_dim2, 
	    cr_offset, ci_dim1, ci_dim2, ci_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, m, n, mmax, ityp;
    extern /* Subroutine */ int vhsec_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);


/*     set coefficient multiplyers */

    /* Parameter adjustments */
    --fnn;
    w_dim1 = *idv;
    w_dim2 = *jdv;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    v_dim1 = *idv;
    v_dim2 = *jdv;
    v_offset = 1 + v_dim1 * (1 + v_dim2);
    v -= v_offset;
    bv_dim1 = *mdb;
    bv_dim2 = *ndb;
    bv_offset = 1 + bv_dim1 * (1 + bv_dim2);
    bv -= bv_offset;
    av_dim1 = *mdb;
    av_dim2 = *ndb;
    av_offset = 1 + av_dim1 * (1 + av_dim2);
    av -= av_offset;
    bs_dim1 = *mdb;
    bs_dim2 = *ndb;
    bs_offset = 1 + bs_dim1 * (1 + bs_dim2);
    bs -= bs_offset;
    as_dim1 = *mdb;
    as_dim2 = *ndb;
    as_offset = 1 + as_dim1 * (1 + as_dim2);
    as -= as_offset;
    ci_dim1 = *mab;
    ci_dim2 = *nlat;
    ci_offset = 1 + ci_dim1 * (1 + ci_dim2);
    ci -= ci_offset;
    cr_dim1 = *mab;
    cr_dim2 = *nlat;
    cr_offset = 1 + cr_dim1 * (1 + cr_dim2);
    cr -= cr_offset;
    bi_dim1 = *mab;
    bi_dim2 = *nlat;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mab;
    br_dim2 = *nlat;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
    --wvhsec;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 2; n <= i__1; ++n) {
	fnn[n] = -sqrt((doublereal) (n * (n - 1)));
    }
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);

/*     compute (v,w) coefficients from as,bs,av,bv */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlat;
	for (n = 1; n <= i__2; ++n) {
	    i__3 = *mab;
	    for (m = 1; m <= i__3; ++m) {
		br[m + (n + k * br_dim2) * br_dim1] = 0.;
		bi[m + (n + k * bi_dim2) * bi_dim1] = 0.;
		cr[m + (n + k * cr_dim2) * cr_dim1] = 0.;
		ci[m + (n + k * ci_dim2) * ci_dim1] = 0.;
	    }
	}

/*     compute m=0 coefficients */

	i__2 = *nlat;
	for (n = 2; n <= i__2; ++n) {
	    br[(n + k * br_dim2) * br_dim1 + 1] = -fnn[n] * av[(n + k * 
		    av_dim2) * av_dim1 + 1];
	    bi[(n + k * bi_dim2) * bi_dim1 + 1] = -fnn[n] * bv[(n + k * 
		    bv_dim2) * bv_dim1 + 1];
	    cr[(n + k * cr_dim2) * cr_dim1 + 1] = fnn[n] * as[(n + k * 
		    as_dim2) * as_dim1 + 1];
	    ci[(n + k * ci_dim2) * ci_dim1 + 1] = fnn[n] * bs[(n + k * 
		    bs_dim2) * bs_dim1 + 1];
	}

/*     compute m>0 coefficients using vector spherepack value for mmax */

	i__2 = mmax;
	for (m = 2; m <= i__2; ++m) {
	    i__3 = *nlat;
	    for (n = m; n <= i__3; ++n) {
		br[m + (n + k * br_dim2) * br_dim1] = -fnn[n] * av[m + (n + k 
			* av_dim2) * av_dim1];
		bi[m + (n + k * bi_dim2) * bi_dim1] = -fnn[n] * bv[m + (n + k 
			* bv_dim2) * bv_dim1];
		cr[m + (n + k * cr_dim2) * cr_dim1] = fnn[n] * as[m + (n + k *
			 as_dim2) * as_dim1];
		ci[m + (n + k * ci_dim2) * ci_dim1] = fnn[n] * bs[m + (n + k *
			 bs_dim2) * bs_dim1];
	    }
	}
    }

/*     synthesize br,bi,cr,ci into (v,w) */

    if (*isym == 0) {
	ityp = 0;
    } else if (*isym == 1) {
	ityp = 3;
    } else if (*isym == 2) {
	ityp = 6;
    }
    vhsec_(nlat, nlon, &ityp, nt, &v[v_offset], &w[w_offset], idv, jdv, &br[
	    br_offset], &bi[bi_offset], &cr[cr_offset], &ci[ci_offset], mab, 
	    nlat, &wvhsec[1], lvhsec, &wk[1], lwk, ierror);
    return 0;
} /* isfvpec1_ */

