/* ivlapes.f -- translated by f2c (version 20061008).
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




/* ... file ivlapes.f */

/*     this file includes documentation and code for */
/*     subroutine ivlapes */

/* ... files which must be loaded with ivlapes.f */

/*     sphcom.f, hrfft.f, vhaes.f, vhses.f */



/*     subroutine ivlapes(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci, */
/*    +mdbc,ndbc,wvhses,lvhses,work,lwork,ierror) */


/*     subroutine ivlapes computes a the vector field (v,w) whose vector */
/*     laplacian is (vlap,wlap).  w and wlap are east longitudinal */
/*     components of the vectors.  v and vlap are colatitudinal components */
/*     of the vectors.  br,bi,cr, and ci are the vector harmonic coefficients */
/*     of (vlap,wlap).  these must be precomputed by vhaes and are input */
/*     parameters to ivlapes.  (v,w) have the same symmetry or lack of */
/*     symmetry about the about the equator as (vlap,wlap).  the input */
/*     parameters ityp,nt,mdbc,ndbc must have the same values used by */
/*     vhaes to compute br,bi,cr, and ci for (vlap,wlap). */


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

/*     ityp   this parameter should have the same value input to subroutine */
/*            vhaes to compute the coefficients br,bi,cr, and ci for the */
/*            vector field (vlap,wlap).  ityp is set as follows: */

/*            = 0  no symmetries exist in (vlap,wlap) about the equator. (v,w) */
/*                 is computed and stored on the entire sphere in the arrays */
/*                 arrays v(i,j) and w(i,j) for i=1,...,nlat and j=1,...,nlon. */

/*            = 1  no symmetries exist in (vlap,wlap) about the equator. (v,w) */
/*                 is computed and stored on the entire sphere in the arrays */
/*                 v(i,j) and w(i,j) for i=1,...,nlat and j=1,...,nlon.  the */
/*                 vorticity of (vlap,wlap) is zero so the coefficients cr and */
/*                 ci are zero and are not used.  the vorticity of (v,w) is */
/*                 also zero. */


/*            = 2  no symmetries exist in (vlap,wlap) about the equator. (v,w) */
/*                 is computed and stored on the entire sphere in the arrays */
/*                 w(i,j) and v(i,j) for i=1,...,nlat and j=1,...,nlon.  the */
/*                 divergence of (vlap,wlap) is zero so the coefficients br and */
/*                 bi are zero and are not used.  the divergence of (v,w) is */
/*                 also zero. */

/*            = 3  wlap is antisymmetric and vlap is symmetric about the */
/*                 equator. consequently w is antisymmetric and v is symmetric. */
/*                 (v,w) is computed and stored on the northern hemisphere */
/*                 only.  if nlat is odd, storage is in the arrays v(i,j), */
/*                 w(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.  if nlat */
/*                 is even, storage is in the arrays v(i,j),w(i,j) for */
/*                 i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 4  wlap is antisymmetric and vlap is symmetric about the */
/*                 equator. consequently w is antisymmetric and v is symmetric. */
/*                 (v,w) is computed and stored on the northern hemisphere */
/*                 only.  if nlat is odd, storage is in the arrays v(i,j), */
/*                 w(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.  if nlat */
/*                 is even, storage is in the arrays v(i,j),w(i,j) for */
/*                 i=1,...,nlat/2 and j=1,...,nlon.  the vorticity of (vlap, */
/*                 wlap) is zero so the coefficients cr,ci are zero and */
/*                 are not used. the vorticity of (v,w) is also zero. */

/*            = 5  wlap is antisymmetric and vlap is symmetric about the */
/*                 equator. consequently w is antisymmetric and v is symmetric. */
/*                 (v,w) is computed and stored on the northern hemisphere */
/*                 only.  if nlat is odd, storage is in the arrays w(i,j), */
/*                 v(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.  if nlat */
/*                 is even, storage is in the arrays w(i,j),v(i,j) for */
/*                 i=1,...,nlat/2 and j=1,...,nlon.  the divergence of (vlap, */
/*                 wlap) is zero so the coefficients br,bi are zero and */
/*                 are not used. the divergence of (v,w) is also zero. */


/*            = 6  wlap is symmetric and vlap is antisymmetric about the */
/*                 equator. consequently w is symmetric and v is antisymmetric. */
/*                 (v,w) is computed and stored on the northern hemisphere */
/*                 only.  if nlat is odd, storage is in the arrays w(i,j), */
/*                 v(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.  if nlat */
/*                 is even, storage is in the arrays w(i,j),v(i,j) for */
/*                 i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 7  wlap is symmetric and vlap is antisymmetric about the */
/*                 equator. consequently w is symmetric and v is antisymmetric. */
/*                 (v,w) is computed and stored on the northern hemisphere */
/*                 only.  if nlat is odd, storage is in the arrays w(i,j), */
/*                 v(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.  if nlat */
/*                 is even, storage is in the arrays w(i,j),v(i,j) for */
/*                 i=1,...,nlat/2 and j=1,...,nlon.  the vorticity of (vlap, */
/*                 wlap) is zero so the coefficients cr,ci are zero and */
/*                 are not used. the vorticity of (v,w) is also zero. */

/*            = 8  wlap is symmetric and vlap is antisymmetric about the */
/*                 equator. consequently w is symmetric and v is antisymmetric. */
/*                 (v,w) is computed and stored on the northern hemisphere */
/*                 only.  if nlat is odd, storage is in the arrays w(i,j), */
/*                 v(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon.  if nlat */
/*                 is even, storage is in the arrays w(i,j),v(i,j) for */
/*                 i=1,...,nlat/2 and j=1,...,nlon.  the divergence of (vlap, */
/*                 wlap) is zero so the coefficients br,bi are zero and */
/*                 are not used. the divergence of (v,w) is also zero. */


/*     nt     nt is the number of vector fields (vlap,wlap). some computational */
/*            efficiency is obtained for multiple fields.  in the program */
/*            that calls ivlapes, the arrays v,w,br,bi,cr and ci can be */
/*            three dimensional corresponding to an indexed multiple vector */
/*            field.  in this case multiple vector synthesis will be performed */
/*            to compute the (v,w) for each field (vlap,wlap).  the third */
/*            index is the synthesis index which assumes the values k=1,...,nt. */
/*            for a single synthesis set nt=1.  the description of the */
/*            remaining parameters is simplified by assuming that nt=1 or */
/*            that all arrays are two dimensional. */

/*   idvw     the first dimension of the arrays w and v as it appears in */
/*            the program that calls ivlapes.  if ityp=0,1, or 2  then idvw */
/*            must be at least nlat.  if ityp > 2 and nlat is even then idvw */
/*            must be at least nlat/2. if ityp > 2 and nlat is odd then idvw */
/*            must be at least (nlat+1)/2. */

/*   jdvw     the second dimension of the arrays w and v as it appears in */
/*            the program that calls ivlapes. jdvw must be at least nlon. */


/*   br,bi    two or three dimensional arrays (see input parameter nt) */
/*   cr,ci    that contain vector spherical harmonic coefficients of the */
/*            vector field (vlap,wlap) as computed by subroutine vhaes. */
/*            br,bi,cr and ci must be computed by vhaes prior to calling */
/*            ivlapes.  if ityp=1,4, or 7 then cr,ci are not used and can */
/*            be dummy arguments.  if ityp=2,5, or 8 then br,bi are not */
/*            used and can be dummy arguments. */

/*    mdbc    the first dimension of the arrays br,bi,cr and ci as it */
/*            appears in the program that calls ivlapes.  mdbc must be */
/*            at least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*    ndbc    the second dimension of the arrays br,bi,cr and ci as it */
/*            appears in the program that calls ivlapes. ndbc must be at */
/*            least nlat. */

/*    wvhses  an array which must be initialized by subroutine vhsesi. */
/*            once initialized, wvhses */
/*            can be used repeatedly by ivlapes as long as nlat and nlon */
/*            remain unchanged.  wvhses must not be altered between calls */
/*            of ivlapes. */

/*    lvhses  the dimension of the array wvhses as it appears in the */
/*            program that calls ivlapes.  let */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd. */

/*            let */

/*               lsavmin = (l1*l2*(nlat+nlat-l1+1))/2+nlon+15 */

/*            then lvhses must be greater than or equal to lsavmin */
/*            (see ierror=9 below). */

/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls ivlapes. define */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */
/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            if ityp .le. 2 then */

/*               (2*nt+1)*nlat*nlon + nlat*(4*nt*l1+1) */

/*            or if ityp .gt. 2 then */

/*               (2*nt+1)*l2*nlon + nlat*(4*nt*l1+1) */

/*            will suffice as a length for lwork. */

/*     ************************************************************** */

/*     output parameters */


/*    v,w     two or three dimensional arrays (see input parameter nt) that */
/*            contain a vector field whose vector laplacian is (vlap,wlap). */
/*            w(i,j) is the east longitude and v(i,j) is the colatitudinal */
/*            component of the vector. v(i,j) and w(i,j) are given on the */
/*            sphere at the colatitude */

/*                 theta(i) = (i-1)*pi/(nlat-1) */

/*            for i=1,...,nlat and east longitude */

/*                 lambda(j) = (j-1)*2*pi/nlon */

/*            for j=1,...,nlon. */

/*            let cost and sint be the cosine and sine at colatitude theta. */
/*            let d( )/dlambda  and d( )/dtheta be the first order partial */
/*            derivatives in longitude and colatitude.  let sf be either v */
/*            or w.  define: */

/*                 del2s(sf) = [d(sint*d(sf)/dtheta)/dtheta + */
/*                               2            2 */
/*                              d (sf)/dlambda /sint]/sint */

/*            then the vector laplacian of (v,w) in (vlap,wlap) satisfies */

/*                 vlap = del2s(v) + (2*cost*dw/dlambda - v)/sint**2 */

/*            and */

/*                 wlap = del2s(w) - (2*cost*dv/dlambda + w)/sint**2 */


/*  ierror    a parameter which flags errors in input parameters as follows: */

/*            = 0  no errors detected */

/*            = 1  error in the specification of nlat */

/*            = 2  error in the specification of nlon */

/*            = 3  error in the specification of ityp */

/*            = 4  error in the specification of nt */

/*            = 5  error in the specification of idvw */

/*            = 6  error in the specification of jdvw */

/*            = 7  error in the specification of mdbc */

/*            = 8  error in the specification of ndbc */

/*            = 9  error in the specification of lvhses */

/*            = 10 error in the specification of lwork */


/* ********************************************************************** */

/*     end of documentation for ivlapes */

/* ********************************************************************** */

/* Subroutine */ int ivlapes_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *
	jdvw, doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, 
	integer *mdbc, integer *ndbc, doublereal *wvhses, integer *lvhses, 
	doublereal *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, br_dim1, 
	    br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, cr_dim2,
	     cr_offset, ci_dim1, ci_dim2, ci_offset, i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int ivlapes1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *);
    integer l1, l2, mn, ibi, ici, ibr, icr, ifn, idz, iwk, imid, mmax, liwk, 
	    lzimn, lwkmin, lsavmin;

    /* Parameter adjustments */
    w_dim1 = *idvw;
    w_dim2 = *jdvw;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    v_dim1 = *idvw;
    v_dim2 = *jdvw;
    v_offset = 1 + v_dim1 * (1 + v_dim2);
    v -= v_offset;
    ci_dim1 = *mdbc;
    ci_dim2 = *ndbc;
    ci_offset = 1 + ci_dim1 * (1 + ci_dim2);
    ci -= ci_offset;
    cr_dim1 = *mdbc;
    cr_dim2 = *ndbc;
    cr_offset = 1 + cr_dim1 * (1 + cr_dim2);
    cr -= cr_offset;
    bi_dim1 = *mdbc;
    bi_dim2 = *ndbc;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mdbc;
    br_dim2 = *ndbc;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
    --wvhses;
    --work;

    /* Function Body */
    *ierror = 1;
    if (*nlat < 3) {
	return 0;
    }
    *ierror = 2;
    if (*nlon < 1) {
	return 0;
    }
    *ierror = 3;
    if (*ityp < 0 || *ityp > 8) {
	return 0;
    }
    *ierror = 4;
    if (*nt < 0) {
	return 0;
    }
    *ierror = 5;
    imid = (*nlat + 1) / 2;
    if (*ityp <= 2 && *idvw < *nlat || *ityp > 2 && *idvw < imid) {
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
    if (*mdbc < mmax) {
	return 0;
    }
    *ierror = 8;
    if (*ndbc < *nlat) {
	return 0;
    }
    *ierror = 9;

/*     set minimum and verify saved workspace length */

    idz = mmax * (*nlat + *nlat - mmax + 1) / 2;
    lzimn = idz * imid;
    lsavmin = lzimn + lzimn + *nlon + 15;
    if (*lvhses < lsavmin) {
	return 0;
    }

/*     set minimum and verify unsaved work space length */

    mn = mmax * *nlat * *nt;
    l2 = (*nlat + 1) / 2;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    l1 = min(i__1,i__2);
    if (*ityp <= 2) {
	lwkmin = ((*nt << 1) + 1) * *nlat * *nlon + *nlat * ((*nt << 2) * l1 
		+ 1);
    } else {
	lwkmin = ((*nt << 1) + 1) * l2 * *nlon + *nlat * ((*nt << 2) * l1 + 1)
		;
    }
    if (*lwork < lwkmin) {
	return 0;
    }
    *ierror = 0;

/*     set work space pointers for vector laplacian coefficients */

    if (*ityp == 0 || *ityp == 3 || *ityp == 6) {
	ibr = 1;
	ibi = ibr + mn;
	icr = ibi + mn;
	ici = icr + mn;
    } else if (*ityp == 1 || *ityp == 4 || *ityp == 7) {
	ibr = 1;
	ibi = ibr + mn;
	icr = ibi + mn;
	ici = icr;
    } else {
	ibr = 1;
	ibi = 1;
	icr = ibi + mn;
	ici = icr + mn;
    }
    ifn = ici + mn;
    iwk = ifn + *nlat;
    if (*ityp == 0 || *ityp == 3 || *ityp == 6) {
	liwk = *lwork - (mn << 2) - *nlat;
    } else {
	liwk = *lwork - (mn << 1) - *nlat;
    }
    ivlapes1_(nlat, nlon, ityp, nt, &v[v_offset], &w[w_offset], idvw, jdvw, &
	    work[ibr], &work[ibi], &work[icr], &work[ici], &mmax, &work[ifn], 
	    mdbc, ndbc, &br[br_offset], &bi[bi_offset], &cr[cr_offset], &ci[
	    ci_offset], &wvhses[1], lvhses, &work[iwk], &liwk, ierror);
    return 0;
} /* ivlapes_ */

/* Subroutine */ int ivlapes1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *
	jdvw, doublereal *brvw, doublereal *bivw, doublereal *crvw, 
	doublereal *civw, integer *mmax, doublereal *fnn, integer *mdbc, 
	integer *ndbc, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, doublereal *wsave, integer *lsave, doublereal *wk, 
	integer *lwk, integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, brvw_dim1, 
	    brvw_dim2, brvw_offset, bivw_dim1, bivw_dim2, bivw_offset, 
	    crvw_dim1, crvw_dim2, crvw_offset, civw_dim1, civw_dim2, 
	    civw_offset, br_dim1, br_dim2, br_offset, bi_dim1, bi_dim2, 
	    bi_offset, cr_dim1, cr_dim2, cr_offset, ci_dim1, ci_dim2, 
	    ci_offset, i__1, i__2, i__3;

    /* Local variables */
    integer k, m, n;
    doublereal fn;
    extern /* Subroutine */ int vhses_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);


/*     preset coefficient multiplyers */

    /* Parameter adjustments */
    --fnn;
    w_dim1 = *idvw;
    w_dim2 = *jdvw;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    v_dim1 = *idvw;
    v_dim2 = *jdvw;
    v_offset = 1 + v_dim1 * (1 + v_dim2);
    v -= v_offset;
    civw_dim1 = *mmax;
    civw_dim2 = *nlat;
    civw_offset = 1 + civw_dim1 * (1 + civw_dim2);
    civw -= civw_offset;
    crvw_dim1 = *mmax;
    crvw_dim2 = *nlat;
    crvw_offset = 1 + crvw_dim1 * (1 + crvw_dim2);
    crvw -= crvw_offset;
    bivw_dim1 = *mmax;
    bivw_dim2 = *nlat;
    bivw_offset = 1 + bivw_dim1 * (1 + bivw_dim2);
    bivw -= bivw_offset;
    brvw_dim1 = *mmax;
    brvw_dim2 = *nlat;
    brvw_offset = 1 + brvw_dim1 * (1 + brvw_dim2);
    brvw -= brvw_offset;
    ci_dim1 = *mdbc;
    ci_dim2 = *ndbc;
    ci_offset = 1 + ci_dim1 * (1 + ci_dim2);
    ci -= ci_offset;
    cr_dim1 = *mdbc;
    cr_dim2 = *ndbc;
    cr_offset = 1 + cr_dim1 * (1 + cr_dim2);
    cr -= cr_offset;
    bi_dim1 = *mdbc;
    bi_dim2 = *ndbc;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mdbc;
    br_dim2 = *ndbc;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
    --wsave;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 2; n <= i__1; ++n) {
	fn = (doublereal) (n - 1);
	fnn[n] = -1. / (fn * (fn + 1.));
/* L1: */
    }

/*     set (u,v) coefficients from br,bi,cr,ci */

    if (*ityp == 0 || *ityp == 3 || *ityp == 6) {

/*     all coefficients needed */

	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *nlat;
	    for (n = 1; n <= i__2; ++n) {
		i__3 = *mmax;
		for (m = 1; m <= i__3; ++m) {
		    brvw[m + (n + k * brvw_dim2) * brvw_dim1] = 0.;
		    bivw[m + (n + k * bivw_dim2) * bivw_dim1] = 0.;
		    crvw[m + (n + k * crvw_dim2) * crvw_dim1] = 0.;
		    civw[m + (n + k * civw_dim2) * civw_dim1] = 0.;
/* L4: */
		}
/* L3: */
	    }
	    i__2 = *nlat;
	    for (n = 2; n <= i__2; ++n) {
		brvw[(n + k * brvw_dim2) * brvw_dim1 + 1] = fnn[n] * br[(n + 
			k * br_dim2) * br_dim1 + 1];
		bivw[(n + k * bivw_dim2) * bivw_dim1 + 1] = fnn[n] * bi[(n + 
			k * bi_dim2) * bi_dim1 + 1];
		crvw[(n + k * crvw_dim2) * crvw_dim1 + 1] = fnn[n] * cr[(n + 
			k * cr_dim2) * cr_dim1 + 1];
		civw[(n + k * civw_dim2) * civw_dim1 + 1] = fnn[n] * ci[(n + 
			k * ci_dim2) * ci_dim1 + 1];
/* L5: */
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    brvw[m + (n + k * brvw_dim2) * brvw_dim1] = fnn[n] * br[m 
			    + (n + k * br_dim2) * br_dim1];
		    bivw[m + (n + k * bivw_dim2) * bivw_dim1] = fnn[n] * bi[m 
			    + (n + k * bi_dim2) * bi_dim1];
		    crvw[m + (n + k * crvw_dim2) * crvw_dim1] = fnn[n] * cr[m 
			    + (n + k * cr_dim2) * cr_dim1];
		    civw[m + (n + k * civw_dim2) * civw_dim1] = fnn[n] * ci[m 
			    + (n + k * ci_dim2) * ci_dim1];
/* L7: */
		}
/* L6: */
	    }
/* L2: */
	}
    } else if (*ityp == 1 || *ityp == 4 || *ityp == 7) {

/*     vorticity is zero so cr,ci=0 not used */

	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *nlat;
	    for (n = 1; n <= i__2; ++n) {
		i__3 = *mmax;
		for (m = 1; m <= i__3; ++m) {
		    brvw[m + (n + k * brvw_dim2) * brvw_dim1] = 0.;
		    bivw[m + (n + k * bivw_dim2) * bivw_dim1] = 0.;
/* L14: */
		}
/* L13: */
	    }
	    i__2 = *nlat;
	    for (n = 2; n <= i__2; ++n) {
		brvw[(n + k * brvw_dim2) * brvw_dim1 + 1] = fnn[n] * br[(n + 
			k * br_dim2) * br_dim1 + 1];
		bivw[(n + k * bivw_dim2) * bivw_dim1 + 1] = fnn[n] * bi[(n + 
			k * bi_dim2) * bi_dim1 + 1];
/* L15: */
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    brvw[m + (n + k * brvw_dim2) * brvw_dim1] = fnn[n] * br[m 
			    + (n + k * br_dim2) * br_dim1];
		    bivw[m + (n + k * bivw_dim2) * bivw_dim1] = fnn[n] * bi[m 
			    + (n + k * bi_dim2) * bi_dim1];
/* L17: */
		}
/* L16: */
	    }
/* L12: */
	}
    } else {

/*     divergence is zero so br,bi=0 not used */

	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *nlat;
	    for (n = 1; n <= i__2; ++n) {
		i__3 = *mmax;
		for (m = 1; m <= i__3; ++m) {
		    crvw[m + (n + k * crvw_dim2) * crvw_dim1] = 0.;
		    civw[m + (n + k * civw_dim2) * civw_dim1] = 0.;
/* L24: */
		}
/* L23: */
	    }
	    i__2 = *nlat;
	    for (n = 2; n <= i__2; ++n) {
		crvw[(n + k * crvw_dim2) * crvw_dim1 + 1] = fnn[n] * cr[(n + 
			k * cr_dim2) * cr_dim1 + 1];
		civw[(n + k * civw_dim2) * civw_dim1 + 1] = fnn[n] * ci[(n + 
			k * ci_dim2) * ci_dim1 + 1];
/* L25: */
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    crvw[m + (n + k * crvw_dim2) * crvw_dim1] = fnn[n] * cr[m 
			    + (n + k * cr_dim2) * cr_dim1];
		    civw[m + (n + k * civw_dim2) * civw_dim1] = fnn[n] * ci[m 
			    + (n + k * ci_dim2) * ci_dim1];
/* L27: */
		}
/* L26: */
	    }
/* L22: */
	}
    }

/*     sythesize coefs into vector field (v,w) */

    vhses_(nlat, nlon, ityp, nt, &v[v_offset], &w[w_offset], idvw, jdvw, &
	    brvw[brvw_offset], &bivw[bivw_offset], &crvw[crvw_offset], &civw[
	    civw_offset], mmax, nlat, &wsave[1], lsave, &wk[1], lwk, ierror);
    return 0;
} /* ivlapes1_ */

