/* idvtgs.f -- translated by f2c (version 20061008).
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



/* ... file idvtgs.f */

/*     this file includes documentation and code for */
/*     subroutine idvtgs         i */

/* ... files which must be loaded with idvtgs.f */

/*     sphcom.f, hrfft.f, vhsgs.f,shags.f, gaqd.f */


/*     subroutine idvtgs(nlat,nlon,isym,nt,v,w,idvw,jdvw,ad,bd,av,bv, */
/*    +mdab,ndab,wvhsgs,lvhsgs,work,lwork,pertbd,pertbv,ierror) */

/*     given the scalar spherical harmonic coefficients ad,bd precomputed */
/*     by subroutine shags for the scalar field divg and coefficients av,bv */
/*     precomputed by subroutine shags for the scalar field vort, subroutine */
/*     idvtgs computes a vector field (v,w) whose divergence is divg - pertbd */
/*     and whose vorticity is vort - pertbv.  w the is east longitude component */
/*     and v is the colatitudinal component of the velocity.  if nt=1 (see nt */
/*     below) pertrbd and pertbv are constants which must be subtracted from */
/*     divg and vort for (v,w) to exist (see the description of pertbd and */
/*     pertrbv below).  usually pertbd and pertbv are zero or small relative */
/*     to divg and vort.  w(i,j) and v(i,j) are the velocity components at */
/*     gaussian colatitude theta(i) (see nlat as input argument) and longitude */
/*     lambda(j) = (j-1)*2*pi/nlon */

/*     the */

/*            divergence(v(i,j),w(i,j)) */

/*         =  [d(sint*v)/dtheta + dw/dlambda]/sint */

/*         =  divg(i,j) - pertbd */

/*     and */

/*            vorticity(v(i,j),w(i,j)) */

/*         =  [-dv/dlambda + d(sint*w)/dtheta]/sint */

/*         =  vort(i,j) - pertbv */

/*     where */

/*            sint = cos(theta(i)). */


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
/*            than 3. the axisymmetric case corresponds to nlon=1. */
/*            the efficiency of the computation is improved when nlon */
/*            is a product of small prime numbers. */


/*     isym   isym determines whether (v,w) are computed on the full or half */
/*            sphere as follows: */

/*      = 0 */
/*            divg,vort are neither pairwise symmetric/antisymmetric nor */
/*            antisymmetric/symmetric about the equator as described for */
/*            isym = 1 or isym = 2  below.  in this case, the vector field */
/*            (v,w) is computed on the entire sphere.  i.e., in the arrays */
/*            w(i,j) and v(i,j) i=1,...,nlat and j=1,...,nlon. */

/*      = 1 */

/*            divg is antisymmetric and vort is symmetric about the equator. */
/*            in this case w is antisymmetric and v is symmetric about the */
/*            equator.  w and v are computed on the northern hemisphere only. */
/*            if nlat is odd they are computed for i=1,...,(nlat+1)/2 */
/*            and j=1,...,nlon.  if nlat is even they are computed for */
/*            i=1,...,nlat/2 and j=1,...,nlon. */

/*      = 2 */

/*            divg is symmetric and vort is antisymmetric about the equator. */
/*            in this case w is symmetric and v is antisymmetric about the */
/*            equator.  w and v are computed on the northern hemisphere only. */
/*            if nlat is odd they are computed for i=1,...,(nlat+1)/2 */
/*            and j=1,...,nlon.  if nlat is even they are computed for */
/*            i=1,...,nlat/2 and j=1,...,nlon. */


/*     nt     in the program that calls idvtgs, nt is the number of scalar */
/*            and vector fields.  some computational efficiency is obtained */
/*            for multiple fields.  the arrays ad,bd,av,bv,u, and v can be */
/*            three dimensional and pertbd,pertbv can be one dimensional */
/*            corresponding to indexed multiple arrays divg, vort.  in this */
/*            case, multiple synthesis will be performed to compute each */
/*            vector field.  the third index for ad,bd,av,bv,v,w and first */
/*            pertrbd,pertbv is the synthesis index which assumes the values */
/*            k=1,...,nt.  for a single synthesis set nt=1. the description of */
/*            remaining parameters is simplified by assuming that nt=1 or that */
/*            ad,bd,av,bv,v,w are two dimensional and pertbd,pertbv are */
/*            constants. */

/*     idvw   the first dimension of the arrays v,w as it appears in */
/*            the program that calls idvtgs. if isym = 0 then idvw */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then idvw must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then idvw must be at least (nlat+1)/2. */

/*     jdvw   the second dimension of the arrays v,w as it appears in */
/*            the program that calls idvtgs. jdvw must be at least nlon. */

/*     ad,bd  two or three dimensional arrays (see input parameter nt) */
/*            that contain scalar spherical harmonic coefficients */
/*            of the divergence array divg as computed by subroutine shags. */

/*     av,bv  two or three dimensional arrays (see input parameter nt) */
/*            that contain scalar spherical harmonic coefficients */
/*            of the vorticity array vort as computed by subroutine shags. */
/*     ***    ad,bd,av,bv must be computed by shags prior to calling idvtgs. */

/*     mdab   the first dimension of the arrays ad,bd,av,bv as it appears */
/*            in the program that calls idvtgs (and shags). mdab must be at */
/*            least min0(nlat,(nlon+2)/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndab   the second dimension of the arrays ad,bd,av,bv as it appears in */
/*            the program that calls idvtgs (and shags). ndab must be at */
/*            least nlat. */

/*  wvhsgs    an array which must be initialized by subroutine vhsgsi. */
/*            wvhsgs can be used repeatedly by idvtgs as long as nlon */
/*            and nlat remain unchanged.  wvhsgs must not be altered */
/*            between calls of idvtgs. */


/*  lvhsgs    the dimension of the array wvhsgs as it appears in the */
/*            program that calls idvtgs. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lvhsgs must be at least */

/*               l1*l2*(nlat+nlat-l1+1)+nlon+15+2*nlat */

/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls idvtgs. define */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */
/*               l1 = min(nlat,nlon/2+1)        if nlon is even or */
/*               l1 = min(nlat,(nlon+1)/2)      if nlon is odd */


/*            if isym = 0 then lwork must be at least */

/*                       nlat*((2*nt+1)*nlon+4*nt*l1+1) */

/*            if isym = 1 or 2 then lwork must be at least */

/*                       (2*nt+1)*l2*nlon+nlat*(4*nt*l1+1) */

/*     ************************************************************** */

/*     output parameters */


/*     v,w   two or three dimensional arrays (see input parameter nt) that */
/*           contain a vector field whose divergence is divg - pertbd and */
/*           whose vorticity is vort - pertbv.  w(i,j) is the east longitude */
/*           component and v(i,j) is the colatitudinal component of velocity */
/*           at the colatitude theta(i) = (i-1)*pi/(nlat-1) and longitude */
/*           lambda(j) = (j-1)*2*pi/nlon for i=1,...,nlat and j=1,...,nlon. */

/*   pertbd  a nt dimensional array (see input parameter nt and assume nt=1 */
/*           for the description that follows).  divg - pertbd is a scalar */
/*           field which can be the divergence of a vector field (v,w). */
/*           pertbd is related to the scalar harmonic coefficients ad,bd */
/*           of divg (computed by shags) by the formula */

/*                pertbd = ad(1,1)/(2.*sqrt(2.)) */

/*           an unperturbed divg can be the divergence of a vector field */
/*           only if ad(1,1) is zero.  if ad(1,1) is nonzero (flagged by */
/*           pertbd nonzero) then subtracting pertbd from divg yields a */
/*           scalar field for which ad(1,1) is zero.  usually pertbd is */
/*           zero or small relative to divg. */

/*   pertbv a nt dimensional array (see input parameter nt and assume nt=1 */
/*           for the description that follows).  vort - pertbv is a scalar */
/*           field which can be the vorticity of a vector field (v,w). */
/*           pertbv is related to the scalar harmonic coefficients av,bv */
/*           of vort (computed by shags) by the formula */

/*                pertbv = av(1,1)/(2.*sqrt(2.)) */

/*           an unperturbed vort can be the vorticity of a vector field */
/*           only if av(1,1) is zero.  if av(1,1) is nonzero (flagged by */
/*           pertbv nonzero) then subtracting pertbv from vort yields a */
/*           scalar field for which av(1,1) is zero.  usually pertbv is */
/*           zero or small relative to vort. */

/*    ierror = 0  no errors */
/*           = 1  error in the specification of nlat */
/*           = 2  error in the specification of nlon */
/*           = 3  error in the specification of isym */
/*           = 4  error in the specification of nt */
/*           = 5  error in the specification of idvw */
/*           = 6  error in the specification of jdvw */
/*           = 7  error in the specification of mdab */
/*           = 8  error in the specification of ndab */
/*           = 9  error in the specification of lvhsgs */
/*           = 10 error in the specification of lwork */
/* ********************************************************************** */


/* Subroutine */ int idvtgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *ad,
	 doublereal *bd, doublereal *av, doublereal *bv, integer *mdab, integer *ndab, doublereal *
	wvhsgs, integer *lvhsgs, doublereal *work, integer *lwork, doublereal *pertbd, 
	doublereal *pertbv, integer *ierror)
{
    /* System generated locals */
    integer w_dim1, w_dim2, w_offset, v_dim1, v_dim2, v_offset, ad_dim1, 
	    ad_dim2, ad_offset, bd_dim1, bd_dim2, bd_offset, av_dim1, av_dim2,
	     av_offset, bv_dim1, bv_dim2, bv_offset, i__1, i__2;

    /* Local variables */
    static integer mn, is, ibi, ici, ibr, icr, idz, iwk, imid, mmax, liwk, 
	    lzimn;
    extern /* Subroutine */ int idvtgs1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);


/*     check input parameters */

    /* Parameter adjustments */
    --pertbv;
    --pertbd;
    w_dim1 = *idvw;
    w_dim2 = *jdvw;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    v_dim1 = *idvw;
    v_dim2 = *jdvw;
    v_offset = 1 + v_dim1 * (1 + v_dim2);
    v -= v_offset;
    bv_dim1 = *mdab;
    bv_dim2 = *ndab;
    bv_offset = 1 + bv_dim1 * (1 + bv_dim2);
    bv -= bv_offset;
    av_dim1 = *mdab;
    av_dim2 = *ndab;
    av_offset = 1 + av_dim1 * (1 + av_dim2);
    av -= av_offset;
    bd_dim1 = *mdab;
    bd_dim2 = *ndab;
    bd_offset = 1 + bd_dim1 * (1 + bd_dim2);
    bd -= bd_offset;
    ad_dim1 = *mdab;
    ad_dim2 = *ndab;
    ad_offset = 1 + ad_dim1 * (1 + ad_dim2);
    ad -= ad_offset;
    --wvhsgs;
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
    if (*isym == 0 && *idvw < *nlat || *isym != 0 && *idvw < imid) {
	return 0;
    }
    *ierror = 6;
    if (*jdvw < *nlon) {
	return 0;
    }
    *ierror = 7;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    if (*mdab < min(i__1,i__2)) {
	return 0;
    }
    *ierror = 8;
    if (*ndab < *nlat) {
	return 0;
    }
    *ierror = 9;
    idz = mmax * (*nlat + *nlat - mmax + 1) / 2;
    lzimn = idz * imid;
    if (*lvhsgs < lzimn + lzimn + *nlon + 15) {
	return 0;
    }
    *ierror = 10;

/*     verify unsaved work space length */

    mn = mmax * *nlat * *nt;
    if (*isym != 0 && *lwork < ((*nt << 1) + 1) * imid * *nlon + (mn << 2) + *
	    nlat) {
	return 0;
    }
    if (*isym == 0 && *lwork < ((*nt << 1) + 1) * *nlat * *nlon + (mn << 2) + 
	    *nlat) {
	return 0;
    }
    *ierror = 0;

/*     set work space pointers */

    ibr = 1;
    ibi = ibr + mn;
    icr = ibi + mn;
    ici = icr + mn;
    is = ici + mn;
    iwk = is + *nlat;
    liwk = *lwork - (mn << 2) - *nlat;
    idvtgs1_(nlat, nlon, isym, nt, &v[v_offset], &w[w_offset], idvw, jdvw, &
	    work[ibr], &work[ibi], &work[icr], &work[ici], &mmax, &work[is], 
	    mdab, ndab, &ad[ad_offset], &bd[bd_offset], &av[av_offset], &bv[
	    bv_offset], &wvhsgs[1], lvhsgs, &work[iwk], &liwk, &pertbd[1], &
	    pertbv[1], ierror);
    return 0;
} /* idvtgs_ */

/* Subroutine */ int idvtgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mmax, doublereal *sqnn, integer *
	mdab, integer *ndab, doublereal *ad, doublereal *bd, doublereal *av, doublereal *bv, doublereal *
	wvhsgs, integer *lvhsgs, doublereal *wk, integer *lwk, doublereal *pertbd, doublereal *
	pertbv, integer *ierror)
{
    /* System generated locals */
    integer w_dim1, w_dim2, w_offset, v_dim1, v_dim2, v_offset, br_dim1, 
	    br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, cr_dim2,
	     cr_offset, ci_dim1, ci_dim2, ci_offset, ad_dim1, ad_dim2, 
	    ad_offset, bd_dim1, bd_dim2, bd_offset, av_dim1, av_dim2, 
	    av_offset, bv_dim1, bv_dim2, bv_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, m, n;
    static doublereal fn;
    static integer ityp;
    extern /* Subroutine */ int vhsgs_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);


/*     preset coefficient multiplyers in vector */

    /* Parameter adjustments */
    --sqnn;
    --pertbv;
    --pertbd;
    w_dim1 = *idvw;
    w_dim2 = *jdvw;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    v_dim1 = *idvw;
    v_dim2 = *jdvw;
    v_offset = 1 + v_dim1 * (1 + v_dim2);
    v -= v_offset;
    ci_dim1 = *mmax;
    ci_dim2 = *nlat;
    ci_offset = 1 + ci_dim1 * (1 + ci_dim2);
    ci -= ci_offset;
    cr_dim1 = *mmax;
    cr_dim2 = *nlat;
    cr_offset = 1 + cr_dim1 * (1 + cr_dim2);
    cr -= cr_offset;
    bi_dim1 = *mmax;
    bi_dim2 = *nlat;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mmax;
    br_dim2 = *nlat;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
    bv_dim1 = *mdab;
    bv_dim2 = *ndab;
    bv_offset = 1 + bv_dim1 * (1 + bv_dim2);
    bv -= bv_offset;
    av_dim1 = *mdab;
    av_dim2 = *ndab;
    av_offset = 1 + av_dim1 * (1 + av_dim2);
    av -= av_offset;
    bd_dim1 = *mdab;
    bd_dim2 = *ndab;
    bd_offset = 1 + bd_dim1 * (1 + bd_dim2);
    bd -= bd_offset;
    ad_dim1 = *mdab;
    ad_dim2 = *ndab;
    ad_offset = 1 + ad_dim1 * (1 + ad_dim2);
    ad -= ad_offset;
    --wvhsgs;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 2; n <= i__1; ++n) {
	fn = (doublereal) (n - 1);
	sqnn[n] = sqrt(fn * (fn + 1.));
/* L1: */
    }

/*     compute multiple vector fields coefficients */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {

/*     set divergence,vorticity perturbation constants */

	pertbd[k] = ad[(k * ad_dim2 + 1) * ad_dim1 + 1] / (sqrt(2.) * 2.);
	pertbv[k] = av[(k * av_dim2 + 1) * av_dim1 + 1] / (sqrt(2.) * 2.);

/*     preset br,bi,cr,ci to 0.0 */

	i__2 = *nlat;
	for (n = 1; n <= i__2; ++n) {
	    i__3 = *mmax;
	    for (m = 1; m <= i__3; ++m) {
		br[m + (n + k * br_dim2) * br_dim1] = 0.;
		bi[m + (n + k * bi_dim2) * bi_dim1] = 0.;
		cr[m + (n + k * cr_dim2) * cr_dim1] = 0.;
		ci[m + (n + k * ci_dim2) * ci_dim1] = 0.;
/* L4: */
	    }
/* L3: */
	}

/*     compute m=0 coefficients */

	i__2 = *nlat;
	for (n = 2; n <= i__2; ++n) {
	    br[(n + k * br_dim2) * br_dim1 + 1] = -ad[(n + k * ad_dim2) * 
		    ad_dim1 + 1] / sqnn[n];
	    bi[(n + k * bi_dim2) * bi_dim1 + 1] = -bd[(n + k * bd_dim2) * 
		    bd_dim1 + 1] / sqnn[n];
	    cr[(n + k * cr_dim2) * cr_dim1 + 1] = av[(n + k * av_dim2) * 
		    av_dim1 + 1] / sqnn[n];
	    ci[(n + k * ci_dim2) * ci_dim1 + 1] = bv[(n + k * bv_dim2) * 
		    bv_dim1 + 1] / sqnn[n];
/* L5: */
	}

/*     compute m>0 coefficients */

	i__2 = *mmax;
	for (m = 2; m <= i__2; ++m) {
	    i__3 = *nlat;
	    for (n = m; n <= i__3; ++n) {
		br[m + (n + k * br_dim2) * br_dim1] = -ad[m + (n + k * 
			ad_dim2) * ad_dim1] / sqnn[n];
		bi[m + (n + k * bi_dim2) * bi_dim1] = -bd[m + (n + k * 
			bd_dim2) * bd_dim1] / sqnn[n];
		cr[m + (n + k * cr_dim2) * cr_dim1] = av[m + (n + k * av_dim2)
			 * av_dim1] / sqnn[n];
		ci[m + (n + k * ci_dim2) * ci_dim1] = bv[m + (n + k * bv_dim2)
			 * bv_dim1] / sqnn[n];
/* L7: */
	    }
/* L6: */
	}
/* L2: */
    }

/*     set ityp for vector synthesis without assuming div=0 or curl=0 */

    if (*isym == 0) {
	ityp = 0;
    } else if (*isym == 1) {
	ityp = 3;
    } else if (*isym == 2) {
	ityp = 6;
    }

/*     sythesize br,bi,cr,ci into the vector field (v,w) */

    vhsgs_(nlat, nlon, &ityp, nt, &v[v_offset], &w[w_offset], idvw, jdvw, &br[
	    br_offset], &bi[bi_offset], &cr[cr_offset], &ci[ci_offset], mmax, 
	    nlat, &wvhsgs[1], lvhsgs, &wk[1], lwk, ierror);
    return 0;
} /* idvtgs1_ */

