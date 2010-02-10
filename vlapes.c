/* vlapes.f -- translated by f2c (version 20061008).
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



/* ... file vlapes.f */

/*     this file includes documentation and code for */
/*     subroutine vlapes          i */

/* ... files which must be loaded with vlapes.f */

/*     sphcom.f, hrfft.f, vhaes.f, vhses.f */




/*     subroutine vlapes(nlat,nlon,ityp,nt,vlap,wlap,idvw,jdvw,br,bi,cr,ci, */
/*    +mdbc,ndbc,wvhses,lvhses,work,lwork,ierror) */


/*     subroutine vlapes computes the vector laplacian of the vector field */
/*     (v,w) in (vlap,wlap) (see the definition of the vector laplacian at */
/*     the output parameter description of vlap,wlap below).  w and wlap */
/*     are east longitudinal components of the vectors.  v and vlap are */
/*     colatitudinal components of the vectors.  br,bi,cr, and ci are the */
/*     vector harmonic coefficients of (v,w).  these must be precomputed by */
/*     vhaes and are input parameters to vlapes.  the laplacian components */
/*     in (vlap,wlap) have the same symmetry or lack of symmetry about the */
/*     equator as (v,w).  the input parameters ityp,nt,mdbc,nbdc must have */
/*     the same values used by vhaes to compute br,bi,cr, and ci for (v,w). */
/*     vlap(i,j) and wlap(i,j) are given on the sphere at the colatitude */

/*            theta(i) = (i-1)*pi/(nlat-1) */

/*     for i=1,...,nlat and east longitude */

/*            lambda(j) = (j-1)*2*pi/nlon */

/*     for j=1,...,nlon. */


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
/*            vector field (v,w).  ityp is set as follows: */

/*            = 0  no symmetries exist in (v,w) about the equator. (vlap,wlap) */
/*                 is computed and stored on the entire sphere in the arrays */
/*                 vlap(i,j) and wlap(i,j) for i=1,...,nlat and j=1,...,nlon. */


/*            = 1  no symmetries exist in (v,w) about the equator. (vlap,wlap) */
/*                 is computed and stored on the entire sphere in the arrays */
/*                 vlap(i,j) and wlap(i,j) for i=1,...,nlat and j=1,...,nlon. */
/*                 the vorticity of (v,w) is zero so the coefficients cr and */
/*                 ci are zero and are not used.  the vorticity of (vlap,wlap) */
/*                 is also zero. */


/*            = 2  no symmetries exist in (v,w) about the equator. (vlap,wlap) */
/*                 is computed and stored on the entire sphere in the arrays */
/*                 vlap(i,j) and wlap(i,j) for i=1,...,nlat and j=1,...,nlon. */
/*                 the divergence of (v,w) is zero so the coefficients br and */
/*                 bi are zero and are not used.  the divergence of (vlap,wlap) */
/*                 is also zero. */

/*            = 3  w is antisymmetric and v is symmetric about the equator. */
/*                 consequently wlap is antisymmetric and vlap is symmetric. */
/*                 (vlap,wlap) is computed and stored on the northern */
/*                 hemisphere only.  if nlat is odd, storage is in the arrays */
/*                 vlap(i,j),wlap(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even, storage is in the arrays vlap(i,j), */
/*                 wlap(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 4  w is antisymmetric and v is symmetric about the equator. */
/*                 consequently wlap is antisymmetric and vlap is symmetric. */
/*                 (vlap,wlap) is computed and stored on the northern */
/*                 hemisphere only.  if nlat is odd, storage is in the arrays */
/*                 vlap(i,j),wlap(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even, storage is in the arrays vlap(i,j), */
/*                 wlap(i,j) for i=1,...,nlat/2 and j=1,...,nlon.  the */
/*                 vorticity of (v,w) is zero so the coefficients cr,ci are */
/*                 zero and are not used. the vorticity of (vlap,wlap) is */
/*                 also zero. */

/*            = 5  w is antisymmetric and v is symmetric about the equator. */
/*                 consequently wlap is antisymmetric and vlap is symmetric. */
/*                 (vlap,wlap) is computed and stored on the northern */
/*                 hemisphere only.  if nlat is odd, storage is in the arrays */
/*                 vlap(i,j),wlap(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even, storage is in the arrays vlap(i,j), */
/*                 wlap(i,j) for i=1,...,nlat/2 and j=1,...,nlon.  the */
/*                 divergence of (v,w) is zero so the coefficients br,bi */
/*                 are zero and are not used. the divergence of (vlap,wlap) */
/*                 is also zero. */


/*            = 6  w is symmetric and v is antisymmetric about the equator. */
/*                 consequently wlap is symmetric and vlap is antisymmetric. */
/*                 (vlap,wlap) is computed and stored on the northern */
/*                 hemisphere only.  if nlat is odd, storage is in the arrays */
/*                 vlap(i,j),wlap(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even, storage is in the arrays vlap(i,j), */
/*                 wlap(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 7  w is symmetric and v is antisymmetric about the equator. */
/*                 consequently wlap is symmetric and vlap is antisymmetric. */
/*                 (vlap,wlap) is computed and stored on the northern */
/*                 hemisphere only.  if nlat is odd, storage is in the arrays */
/*                 vlap(i,j),wlap(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even, storage is in the arrays vlap(i,j), */
/*                 wlap(i,j) for i=1,...,nlat/2 and j=1,...,nlon.  the */
/*                 vorticity of (v,w) is zero so the coefficients cr,ci are */
/*                 zero and are not used. the vorticity of (vlap,wlap) is */
/*                 also zero. */

/*            = 8  w is symmetric and v is antisymmetric about the equator. */
/*                 consequently wlap is symmetric and vlap is antisymmetric. */
/*                 (vlap,wlap) is computed and stored on the northern */
/*                 hemisphere only.  if nlat is odd, storage is in the arrays */
/*                 vlap(i,j),wlap(i,j) for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even, storage is in the arrays vlap(i,j), */
/*                 wlap(i,j) for i=1,...,nlat/2 and j=1,...,nlon.  the */
/*                 divergence of (v,w) is zero so the coefficients br,bi */
/*                 are zero and are not used. the divergence of (vlap,wlap) */
/*                 is also zero. */


/*     nt     nt is the number of vector fields (v,w).  some computational */
/*            efficiency is obtained for multiple fields.  in the program */
/*            that calls vlapes, the arrays vlap,wlap,br,bi,cr and ci */
/*            can be three dimensional corresponding to an indexed multiple */
/*            vector field.  in this case multiple vector synthesis will */
/*            be performed to compute the vector laplacian for each field. */
/*            the third index is the synthesis index which assumes the values */
/*            k=1,...,nt.  for a single synthesis set nt=1.  the description */
/*            of the remaining parameters is simplified by assuming that nt=1 */
/*            or that all arrays are two dimensional. */

/*   idvw     the first dimension of the arrays vlap and wlap as it appears */
/*            in the program that calls vlapes.  if ityp=0,1, or 2  then idvw */
/*            must be at least nlat.  if ityp > 2 and nlat is even then idvw */
/*            must be at least nlat/2. if ityp > 2 and nlat is odd then idvw */
/*            must be at least (nlat+1)/2. */

/*   jdvw     the second dimension of the arrays vlap and wlap as it appears */
/*            in the program that calls vlapes. jdvw must be at least nlon. */


/*   br,bi    two or three dimensional arrays (see input parameter nt) */
/*   cr,ci    that contain vector spherical harmonic coefficients */
/*            of the vector field (v,w) as computed by subroutine vhaes. */
/*            br,bi,cr and ci must be computed by vhaes prior to calling */
/*            vlapes.  if ityp=1,4, or 7 then cr,ci are not used and can */
/*            be dummy arguments.  if ityp=2,5, or 8 then br,bi are not */
/*            used and can be dummy arguments. */

/*    mdbc    the first dimension of the arrays br,bi,cr and ci as it */
/*            appears in the program that calls vlapes.  mdbc must be */
/*            at least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*    ndbc    the second dimension of the arrays br,bi,cr and ci as it */
/*            appears in the program that calls vlapes. ndbc must be at */
/*            least nlat. */

/*    wvhses  an array which must be initialized by subroutine vhsesi. */
/*            once initialized, vhses */
/*            can be used repeatedly by vlapes as long as nlat and nlon */
/*            remain unchanged.  wvhses must not be altered between calls */
/*            of vlapes. */

/*    lvhses  the dimension of the array wvhses as it appears in the */
/*            program that calls vlapes.  let */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd. */

/*            then lvhses must be greater than or equal */

/*               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15 */

/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls vlapes. define */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */
/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            if ityp .le. 2 then */

/*               (2*nt+1)*nlat*nlon + nlat*(4*nt*l1+1) */

/*            or if ityp .gt. 2 then */

/*               (2*nt+1)*l2*nlon + nlat*(4*nt*l1+1) */

/*            will suffice as a length for lwork. */

/*     ************************************************************** */

/*     output parameters */


/*    vlap,   two or three dimensional arrays (see input parameter nt) that */
/*    wlap    contain the vector laplacian of the field (v,w).  wlap(i,j) is */
/*            the east longitude component and vlap(i,j) is the colatitudinal */
/*            component of the vector laplacian.  the definition of the */
/*            vector laplacian follows: */

/*            let cost and sint be the cosine and sine at colatitude theta. */
/*            let d( )/dlambda  and d( )/dtheta be the first order partial */
/*            derivatives in longitude and colatitude.  let del2 be the scalar */
/*            laplacian operator */

/*                 del2(s) = [d(sint*d(s)/dtheta)/dtheta + */
/*                             2            2 */
/*                            d (s)/dlambda /sint]/sint */

/*            then the vector laplacian opeator */

/*                 dvel2(v,w) = (vlap,wlap) */

/*            is defined by */

/*                 vlap = del2(v) - (2*cost*dw/dlambda + v)/sint**2 */

/*                 wlap = del2(w) + (2*cost*dv/dlambda - w)/sint**2 */

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

/*     end of documentation for vlapes */

/* ********************************************************************** */

/* Subroutine */ int vlapes_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vlap, doublereal *wlap, integer *idvw, integer *jdvw, 
	doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, integer *mdbc, integer *ndbc, 
	doublereal *wvhses, integer *lvhses, doublereal *work, integer *lwork, integer *
	ierror)
{
    /* System generated locals */
    integer vlap_dim1, vlap_dim2, vlap_offset, wlap_dim1, wlap_dim2, 
	    wlap_offset, br_dim1, br_dim2, br_offset, bi_dim1, bi_dim2, 
	    bi_offset, cr_dim1, cr_dim2, cr_offset, ci_dim1, ci_dim2, 
	    ci_offset, i__1, i__2;

    /* Local variables */
    static integer l1, l2, mn, ibi, ici, ibr, icr, ifn, idz, iwk, imid, mmax, 
	    liwk, lzimn, lwkmin;
    extern /* Subroutine */ int vlapes1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    static integer lsavmin;

    /* Parameter adjustments */
    wlap_dim1 = *idvw;
    wlap_dim2 = *jdvw;
    wlap_offset = 1 + wlap_dim1 * (1 + wlap_dim2);
    wlap -= wlap_offset;
    vlap_dim1 = *idvw;
    vlap_dim2 = *jdvw;
    vlap_offset = 1 + vlap_dim1 * (1 + vlap_dim2);
    vlap -= vlap_offset;
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
    idz = mmax * (*nlat + *nlat - mmax + 1) / 2;
    lzimn = idz * imid;
    lsavmin = lzimn + lzimn + *nlon + 15;
    if (*lvhses < lsavmin) {
	return 0;
    }

/*     verify unsaved work space length */

    mn = mmax * *nlat * *nt;
    l2 = (*nlat + 1) / 2;
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
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
    liwk = *lwork - (mn << 2) - *nlat;
    vlapes1_(nlat, nlon, ityp, nt, &vlap[vlap_offset], &wlap[wlap_offset], 
	    idvw, jdvw, &work[ibr], &work[ibi], &work[icr], &work[ici], &mmax,
	     &work[ifn], mdbc, ndbc, &br[br_offset], &bi[bi_offset], &cr[
	    cr_offset], &ci[ci_offset], &wvhses[1], lvhses, &work[iwk], &liwk,
	     ierror);
    return 0;
} /* vlapes_ */

/* Subroutine */ int vlapes1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vlap, doublereal *wlap, integer *idvw, integer *jdvw, 
	doublereal *brlap, doublereal *bilap, doublereal *crlap, doublereal *cilap, integer *mmax, 
	doublereal *fnn, integer *mdb, integer *ndb, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, doublereal *wsave, integer *lsave, doublereal *wk, integer *lwk, 
	integer *ierror)
{
    /* System generated locals */
    integer vlap_dim1, vlap_dim2, vlap_offset, wlap_dim1, wlap_dim2, 
	    wlap_offset, brlap_dim1, brlap_dim2, brlap_offset, bilap_dim1, 
	    bilap_dim2, bilap_offset, crlap_dim1, crlap_dim2, crlap_offset, 
	    cilap_dim1, cilap_dim2, cilap_offset, br_dim1, br_dim2, br_offset,
	     bi_dim1, bi_dim2, bi_offset, cr_dim1, cr_dim2, cr_offset, 
	    ci_dim1, ci_dim2, ci_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer k, m, n;
    static doublereal fn;
    extern /* Subroutine */ int vhses_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);


/*     preset coefficient multiplyers */

    /* Parameter adjustments */
    --fnn;
    wlap_dim1 = *idvw;
    wlap_dim2 = *jdvw;
    wlap_offset = 1 + wlap_dim1 * (1 + wlap_dim2);
    wlap -= wlap_offset;
    vlap_dim1 = *idvw;
    vlap_dim2 = *jdvw;
    vlap_offset = 1 + vlap_dim1 * (1 + vlap_dim2);
    vlap -= vlap_offset;
    cilap_dim1 = *mmax;
    cilap_dim2 = *nlat;
    cilap_offset = 1 + cilap_dim1 * (1 + cilap_dim2);
    cilap -= cilap_offset;
    crlap_dim1 = *mmax;
    crlap_dim2 = *nlat;
    crlap_offset = 1 + crlap_dim1 * (1 + crlap_dim2);
    crlap -= crlap_offset;
    bilap_dim1 = *mmax;
    bilap_dim2 = *nlat;
    bilap_offset = 1 + bilap_dim1 * (1 + bilap_dim2);
    bilap -= bilap_offset;
    brlap_dim1 = *mmax;
    brlap_dim2 = *nlat;
    brlap_offset = 1 + brlap_dim1 * (1 + brlap_dim2);
    brlap -= brlap_offset;
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
    --wsave;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 2; n <= i__1; ++n) {
	fn = (doublereal) (n - 1);
	fnn[n] = -fn * (fn + 1.);
/* L1: */
    }

/*     set laplacian coefficients from br,bi,cr,ci */

    if (*ityp == 0 || *ityp == 3 || *ityp == 6) {

/*     all coefficients needed */

	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = *nlat;
	    for (n = 1; n <= i__2; ++n) {
		i__3 = *mmax;
		for (m = 1; m <= i__3; ++m) {
		    brlap[m + (n + k * brlap_dim2) * brlap_dim1] = 0.;
		    bilap[m + (n + k * bilap_dim2) * bilap_dim1] = 0.;
		    crlap[m + (n + k * crlap_dim2) * crlap_dim1] = 0.;
		    cilap[m + (n + k * cilap_dim2) * cilap_dim1] = 0.;
/* L4: */
		}
/* L3: */
	    }
	    i__2 = *nlat;
	    for (n = 2; n <= i__2; ++n) {
		brlap[(n + k * brlap_dim2) * brlap_dim1 + 1] = fnn[n] * br[(n 
			+ k * br_dim2) * br_dim1 + 1];
		bilap[(n + k * bilap_dim2) * bilap_dim1 + 1] = fnn[n] * bi[(n 
			+ k * bi_dim2) * bi_dim1 + 1];
		crlap[(n + k * crlap_dim2) * crlap_dim1 + 1] = fnn[n] * cr[(n 
			+ k * cr_dim2) * cr_dim1 + 1];
		cilap[(n + k * cilap_dim2) * cilap_dim1 + 1] = fnn[n] * ci[(n 
			+ k * ci_dim2) * ci_dim1 + 1];
/* L5: */
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    brlap[m + (n + k * brlap_dim2) * brlap_dim1] = fnn[n] * 
			    br[m + (n + k * br_dim2) * br_dim1];
		    bilap[m + (n + k * bilap_dim2) * bilap_dim1] = fnn[n] * 
			    bi[m + (n + k * bi_dim2) * bi_dim1];
		    crlap[m + (n + k * crlap_dim2) * crlap_dim1] = fnn[n] * 
			    cr[m + (n + k * cr_dim2) * cr_dim1];
		    cilap[m + (n + k * cilap_dim2) * cilap_dim1] = fnn[n] * 
			    ci[m + (n + k * ci_dim2) * ci_dim1];
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
		    brlap[m + (n + k * brlap_dim2) * brlap_dim1] = 0.;
		    bilap[m + (n + k * bilap_dim2) * bilap_dim1] = 0.;
/* L14: */
		}
/* L13: */
	    }
	    i__2 = *nlat;
	    for (n = 2; n <= i__2; ++n) {
		brlap[(n + k * brlap_dim2) * brlap_dim1 + 1] = fnn[n] * br[(n 
			+ k * br_dim2) * br_dim1 + 1];
		bilap[(n + k * bilap_dim2) * bilap_dim1 + 1] = fnn[n] * bi[(n 
			+ k * bi_dim2) * bi_dim1 + 1];
/* L15: */
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    brlap[m + (n + k * brlap_dim2) * brlap_dim1] = fnn[n] * 
			    br[m + (n + k * br_dim2) * br_dim1];
		    bilap[m + (n + k * bilap_dim2) * bilap_dim1] = fnn[n] * 
			    bi[m + (n + k * bi_dim2) * bi_dim1];
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
		    crlap[m + (n + k * crlap_dim2) * crlap_dim1] = 0.;
		    cilap[m + (n + k * cilap_dim2) * cilap_dim1] = 0.;
/* L24: */
		}
/* L23: */
	    }
	    i__2 = *nlat;
	    for (n = 2; n <= i__2; ++n) {
		crlap[(n + k * crlap_dim2) * crlap_dim1 + 1] = fnn[n] * cr[(n 
			+ k * cr_dim2) * cr_dim1 + 1];
		cilap[(n + k * cilap_dim2) * cilap_dim1 + 1] = fnn[n] * ci[(n 
			+ k * ci_dim2) * ci_dim1 + 1];
/* L25: */
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    crlap[m + (n + k * crlap_dim2) * crlap_dim1] = fnn[n] * 
			    cr[m + (n + k * cr_dim2) * cr_dim1];
		    cilap[m + (n + k * cilap_dim2) * cilap_dim1] = fnn[n] * 
			    ci[m + (n + k * ci_dim2) * ci_dim1];
/* L27: */
		}
/* L26: */
	    }
/* L22: */
	}
    }

/*     sythesize coefs into vector field (vlap,wlap) */

    vhses_(nlat, nlon, ityp, nt, &vlap[vlap_offset], &wlap[wlap_offset], idvw,
	     jdvw, &brlap[brlap_offset], &bilap[bilap_offset], &crlap[
	    crlap_offset], &cilap[cilap_offset], mmax, nlat, &wsave[1], lsave,
	     &wk[1], lwk, ierror);
    return 0;
} /* vlapes1_ */

