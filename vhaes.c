/* vhaes.f -- translated by f2c (version 20061008).
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

/* Table of constant values */

static integer c__0 = 0;


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



/* ... file vhaes.f */

/*     this file contains code and documentation for subroutines */
/*     vhaes and vhaesi */

/* ... files which must be loaded with vhaes.f */

/*     sphcom.f, hrfft.f */


/*     subroutine vhaes(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci, */
/*    +                 mdab,ndab,wvhaes,lvhaes,work,lwork,ierror) */

/*     subroutine vhaes performs the vector spherical harmonic analysis */
/*     on the vector field (v,w) and stores the result in the arrays */
/*     br, bi, cr, and ci. v(i,j) and w(i,j) are the colatitudinal */
/*     (measured from the north pole) and east longitudinal components */
/*     respectively, located at colatitude theta(i) = (i-1)*pi/(nlat-1) */
/*     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral */
/*     representation of (v,w) is given at output parameters v,w in */
/*     subroutine vhses. */

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
/*            nlon = 72 for a five degree grid. nlon must be greater */
/*            than zero. the axisymmetric case corresponds to nlon=1. */
/*            the efficiency of the computation is improved when nlon */
/*            is a product of small prime numbers. */

/*     ityp   = 0  no symmetries exist about the equator. the analysis */
/*                 is performed on the entire sphere.  i.e. on the */
/*                 arrays v(i,j),w(i,j) for i=1,...,nlat and */
/*                 j=1,...,nlon. */

/*            = 1  no symmetries exist about the equator. the analysis */
/*                 is performed on the entire sphere.  i.e. on the */
/*                 arrays v(i,j),w(i,j) for i=1,...,nlat and */
/*                 j=1,...,nlon. the curl of (v,w) is zero. that is, */
/*                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. */
/*                 the coefficients cr and ci are zero. */

/*            = 2  no symmetries exist about the equator. the analysis */
/*                 is performed on the entire sphere.  i.e. on the */
/*                 arrays v(i,j),w(i,j) for i=1,...,nlat and */
/*                 j=1,...,nlon. the divergence of (v,w) is zero. i.e., */
/*                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. */
/*                 the coefficients br and bi are zero. */

/*            = 3  v is symmetric and w is antisymmetric about the */
/*                 equator. the analysis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the analysis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the analysis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 4  v is symmetric and w is antisymmetric about the */
/*                 equator. the analysis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the analysis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the analysis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */
/*                 the curl of (v,w) is zero. that is, */
/*                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. */
/*                 the coefficients cr and ci are zero. */

/*            = 5  v is symmetric and w is antisymmetric about the */
/*                 equator. the analysis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the analysis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the analysis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */
/*                 the divergence of (v,w) is zero. i.e., */
/*                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. */
/*                 the coefficients br and bi are zero. */

/*            = 6  v is antisymmetric and w is symmetric about the */
/*                 equator. the analysis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the analysis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the analysis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 7  v is antisymmetric and w is symmetric about the */
/*                 equator. the analysis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the analysis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the analysis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */
/*                 the curl of (v,w) is zero. that is, */
/*                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. */
/*                 the coefficients cr and ci are zero. */

/*            = 8  v is antisymmetric and w is symmetric about the */
/*                 equator. the analysis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the analysis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the analysis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */
/*                 the divergence of (v,w) is zero. i.e., */
/*                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. */
/*                 the coefficients br and bi are zero. */


/*     nt     the number of analyses.  in the program that calls vhaes, */
/*            the arrays v,w,br,bi,cr, and ci can be three dimensional */
/*            in which case multiple analyses will be performed. */
/*            the third index is the analysis index which assumes the */
/*            values k=1,...,nt.  for a single analysis set nt=1. the */
/*            discription of the remaining parameters is simplified */
/*            by assuming that nt=1 or that all the arrays are two */
/*            dimensional. */

/*     v,w    two or three dimensional arrays (see input parameter nt) */
/*            that contain the vector function to be analyzed. */
/*            v is the colatitudnal component and w is the east */
/*            longitudinal component. v(i,j),w(i,j) contain the */
/*            components at colatitude theta(i) = (i-1)*pi/(nlat-1) */
/*            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges */
/*            are defined above at the input parameter ityp. */

/*     idvw   the first dimension of the arrays v,w as it appears in */
/*            the program that calls vhaes. if ityp .le. 2 then idvw */
/*            must be at least nlat.  if ityp .gt. 2 and nlat is */
/*            even then idvw must be at least nlat/2. if ityp .gt. 2 */
/*            and nlat is odd then idvw must be at least (nlat+1)/2. */

/*     jdvw   the second dimension of the arrays v,w as it appears in */
/*            the program that calls vhaes. jdvw must be at least nlon. */

/*     mdab   the first dimension of the arrays br,bi,cr, and ci as it */
/*            appears in the program that calls vhaes. mdab must be at */
/*            least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndab   the second dimension of the arrays br,bi,cr, and ci as it */
/*            appears in the program that calls vhaes. ndab must be at */
/*            least nlat. */

/*     lvhaes an array which must be initialized by subroutine vhaesi. */
/*            once initialized, wvhaes can be used repeatedly by vhaes */
/*            as long as nlon and nlat remain unchanged.  wvhaes must */
/*            not be altered between calls of vhaes. */

/*     lvhaes the dimension of the array wvhaes as it appears in the */
/*            program that calls vhaes. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lvhaes must be at least */

/*            l1*l2(nlat+nlat-l1+1)+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls vhaes. define */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            if ityp .le. 2 then lwork must be at least */

/*                       (2*nt+1)*nlat*nlon */

/*            if ityp .gt. 2 then lwork must be at least */

/*                        (2*nt+1)*l2*nlon */

/*     ************************************************************** */

/*     output parameters */

/*     br,bi  two or three dimensional arrays (see input parameter nt) */
/*     cr,ci  that contain the vector spherical harmonic coefficients */
/*            in the spectral representation of v(i,j) and w(i,j) given */
/*            in the discription of subroutine vhses. br(mp1,np1), */
/*            bi(mp1,np1),cr(mp1,np1), and ci(mp1,np1) are computed */
/*            for mp1=1,...,mmax and np1=mp1,...,nlat except for np1=nlat */
/*            and odd mp1. mmax=min0(nlat,nlon/2) if nlon is even or */
/*            mmax=min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of ityp */
/*            = 4  error in the specification of nt */
/*            = 5  error in the specification of idvw */
/*            = 6  error in the specification of jdvw */
/*            = 7  error in the specification of mdab */
/*            = 8  error in the specification of ndab */
/*            = 9  error in the specification of lvhaes */
/*            = 10 error in the specification of lwork */

/* ******************************************************** */

/*     subroutine vhaesi(nlat,nlon,wvhaes,lvhaes,work,lwork,dwork, */
/*    +                  ldwork,ierror) */

/*     subroutine vhaesi initializes the array wvhaes which can then be */
/*     used repeatedly by subroutine vhaes until nlat or nlon is changed. */

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
/*            nlon = 72 for a five degree grid. nlon must be greater */
/*            than zero. the axisymmetric case corresponds to nlon=1. */
/*            the efficiency of the computation is improved when nlon */
/*            is a product of small prime numbers. */

/*     lvhaes the dimension of the array wvhaes as it appears in the */
/*            program that calls vhaes. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lvhaes must be at least */

/*               l1*l2*(nlat+nlat-l1+1)+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls vhaes. lwork must be at least */

/*              3*(max0(l1-2,0)*(nlat+nlat-l1-1))/2+5*l2*nlat */

/*     dwork  an unsaved double precision work space */

/*     ldwork the length of the array dwork as it appears in the */
/*            program that calls vhaesi.  ldwork must be at least */
/*            2*(nlat+1) */


/*     ************************************************************** */

/*     output parameters */

/*     wvhaes an array which is initialized for use by subroutine vhaes. */
/*            once initialized, wvhaes can be used repeatedly by vhaes */
/*            as long as nlat or nlon remain unchanged.  wvhaes must not */
/*            be altered between calls of vhaes. */


/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of lvhaes */
/*            = 4  error in the specification of lwork */
/*            = 5  error in the specification of ldwork */


/* Subroutine */ int vhaes_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *
	jdvw, doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, 
	integer *mdab, integer *ndab, doublereal *wvhaes, integer *lvhaes, 
	doublereal *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, br_dim1, 
	    br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, cr_dim2,
	     cr_offset, ci_dim1, ci_dim2, ci_offset, i__1, i__2;

    /* Local variables */
    integer iw1, iw2, iw3, iw4, jw1, jw2, idv, lnl, idz, ist, imid, mmax, 
	    lzimn;
    extern /* Subroutine */ int vhaes1_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, doublereal *
	    , doublereal *, doublereal *);

    /* Parameter adjustments */
    w_dim1 = *idvw;
    w_dim2 = *jdvw;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    v_dim1 = *idvw;
    v_dim2 = *jdvw;
    v_offset = 1 + v_dim1 * (1 + v_dim2);
    v -= v_offset;
    ci_dim1 = *mdab;
    ci_dim2 = *ndab;
    ci_offset = 1 + ci_dim1 * (1 + ci_dim2);
    ci -= ci_offset;
    cr_dim1 = *mdab;
    cr_dim2 = *ndab;
    cr_offset = 1 + cr_dim1 * (1 + cr_dim2);
    cr -= cr_offset;
    bi_dim1 = *mdab;
    bi_dim2 = *ndab;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mdab;
    br_dim2 = *ndab;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
    --wvhaes;
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
    if (*mdab < mmax) {
	return 0;
    }
    *ierror = 8;
    if (*ndab < *nlat) {
	return 0;
    }
    *ierror = 9;
    idz = mmax * (*nlat + *nlat - mmax + 1) / 2;
    lzimn = idz * imid;
    if (*lvhaes < lzimn + lzimn + *nlon + 15) {
	return 0;
    }
    *ierror = 10;
    idv = *nlat;
    if (*ityp > 2) {
	idv = imid;
    }
    lnl = *nt * idv * *nlon;
    if (*lwork < lnl + lnl + idv * *nlon) {
	return 0;
    }
    *ierror = 0;
    ist = 0;
    if (*ityp <= 2) {
	ist = imid;
    }
    iw1 = ist + 1;
    iw2 = lnl + 1;
    iw3 = iw2 + ist;
    iw4 = iw2 + lnl;
    jw1 = lzimn + 1;
    jw2 = jw1 + lzimn;
    vhaes1_(nlat, nlon, ityp, nt, &imid, idvw, jdvw, &v[v_offset], &w[
	    w_offset], mdab, ndab, &br[br_offset], &bi[bi_offset], &cr[
	    cr_offset], &ci[ci_offset], &idv, &work[1], &work[iw1], &work[iw2]
	    , &work[iw3], &work[iw4], &idz, &wvhaes[1], &wvhaes[jw1], &wvhaes[
	    jw2]);
    return 0;
} /* vhaes_ */

/* Subroutine */ int vhaes1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *
	v, doublereal *w, integer *mdab, integer *ndab, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci, integer *idv, 
	doublereal *ve, doublereal *vo, doublereal *we, doublereal *wo, 
	doublereal *work, integer *idz, doublereal *zv, doublereal *zw, 
	doublereal *wrfft)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, br_dim1, 
	    br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, cr_dim2,
	     cr_offset, ci_dim1, ci_dim2, ci_offset, ve_dim1, ve_dim2, 
	    ve_offset, vo_dim1, vo_dim2, vo_offset, we_dim1, we_dim2, 
	    we_offset, wo_dim1, wo_dim2, wo_offset, zv_dim1, zv_offset, 
	    zw_dim1, zw_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer i__, j, k, m, mb, mp1, np1, mp2;
    doublereal fsn, tsn;
    integer ndo1, ndo2, imm1, nlp1, mlat, mmax, mlon, itypp;
    extern /* Subroutine */ int hrfftf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *);

    /* Parameter adjustments */
    w_dim1 = *idvw;
    w_dim2 = *jdvw;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    v_dim1 = *idvw;
    v_dim2 = *jdvw;
    v_offset = 1 + v_dim1 * (1 + v_dim2);
    v -= v_offset;
    ci_dim1 = *mdab;
    ci_dim2 = *ndab;
    ci_offset = 1 + ci_dim1 * (1 + ci_dim2);
    ci -= ci_offset;
    cr_dim1 = *mdab;
    cr_dim2 = *ndab;
    cr_offset = 1 + cr_dim1 * (1 + cr_dim2);
    cr -= cr_offset;
    bi_dim1 = *mdab;
    bi_dim2 = *ndab;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mdab;
    br_dim2 = *ndab;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
    wo_dim1 = *idv;
    wo_dim2 = *nlon;
    wo_offset = 1 + wo_dim1 * (1 + wo_dim2);
    wo -= wo_offset;
    we_dim1 = *idv;
    we_dim2 = *nlon;
    we_offset = 1 + we_dim1 * (1 + we_dim2);
    we -= we_offset;
    vo_dim1 = *idv;
    vo_dim2 = *nlon;
    vo_offset = 1 + vo_dim1 * (1 + vo_dim2);
    vo -= vo_offset;
    ve_dim1 = *idv;
    ve_dim2 = *nlon;
    ve_offset = 1 + ve_dim1 * (1 + ve_dim2);
    ve -= ve_offset;
    --work;
    zw_dim1 = *idz;
    zw_offset = 1 + zw_dim1;
    zw -= zw_offset;
    zv_dim1 = *idz;
    zv_offset = 1 + zv_dim1;
    zv -= zv_offset;
    --wrfft;

    /* Function Body */
    nlp1 = *nlat + 1;
    tsn = 2. / *nlon;
    fsn = 4. / *nlon;
    mlat = *nlat % 2;
    mlon = *nlon % 2;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
    imm1 = *imid;
    if (mlat != 0) {
	imm1 = *imid - 1;
    }
    if (*ityp > 2) {
	goto L3;
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = imm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *nlon;
	    for (j = 1; j <= i__3; ++j) {
		ve[i__ + (j + k * ve_dim2) * ve_dim1] = tsn * (v[i__ + (j + k 
			* v_dim2) * v_dim1] + v[nlp1 - i__ + (j + k * v_dim2) 
			* v_dim1]);
		vo[i__ + (j + k * vo_dim2) * vo_dim1] = tsn * (v[i__ + (j + k 
			* v_dim2) * v_dim1] - v[nlp1 - i__ + (j + k * v_dim2) 
			* v_dim1]);
		we[i__ + (j + k * we_dim2) * we_dim1] = tsn * (w[i__ + (j + k 
			* w_dim2) * w_dim1] + w[nlp1 - i__ + (j + k * w_dim2) 
			* w_dim1]);
		wo[i__ + (j + k * wo_dim2) * wo_dim1] = tsn * (w[i__ + (j + k 
			* w_dim2) * w_dim1] - w[nlp1 - i__ + (j + k * w_dim2) 
			* w_dim1]);
/* L5: */
	    }
	}
    }
    goto L2;
L3:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = imm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = *nlon;
	    for (j = 1; j <= i__1; ++j) {
		ve[i__ + (j + k * ve_dim2) * ve_dim1] = fsn * v[i__ + (j + k *
			 v_dim2) * v_dim1];
		vo[i__ + (j + k * vo_dim2) * vo_dim1] = fsn * v[i__ + (j + k *
			 v_dim2) * v_dim1];
		we[i__ + (j + k * we_dim2) * we_dim1] = fsn * w[i__ + (j + k *
			 w_dim2) * w_dim1];
		wo[i__ + (j + k * wo_dim2) * wo_dim1] = fsn * w[i__ + (j + k *
			 w_dim2) * w_dim1];
/* L8: */
	    }
	}
    }
L2:
    if (mlat == 0) {
	goto L7;
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    ve[*imid + (j + k * ve_dim2) * ve_dim1] = tsn * v[*imid + (j + k *
		     v_dim2) * v_dim1];
	    we[*imid + (j + k * we_dim2) * we_dim1] = tsn * w[*imid + (j + k *
		     w_dim2) * w_dim1];
/* L6: */
	}
    }
L7:
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	hrfftf_(idv, nlon, &ve[(k * ve_dim2 + 1) * ve_dim1 + 1], idv, &wrfft[
		1], &work[1]);
	hrfftf_(idv, nlon, &we[(k * we_dim2 + 1) * we_dim1 + 1], idv, &wrfft[
		1], &work[1]);
/* L9: */
    }
    ndo1 = *nlat;
    ndo2 = *nlat;
    if (mlat != 0) {
	ndo1 = *nlat - 1;
    }
    if (mlat == 0) {
	ndo2 = *nlat - 1;
    }
    if (*ityp == 2 || *ityp == 5 || *ityp == 8) {
	goto L11;
    }
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__1 = mmax;
	for (mp1 = 1; mp1 <= i__1; ++mp1) {
	    i__3 = *nlat;
	    for (np1 = mp1; np1 <= i__3; ++np1) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] = 0.;
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = 0.;
/* L10: */
	    }
	}
    }
L11:
    if (*ityp == 1 || *ityp == 4 || *ityp == 7) {
	goto L13;
    }
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__1 = mmax;
	for (mp1 = 1; mp1 <= i__1; ++mp1) {
	    i__2 = *nlat;
	    for (np1 = mp1; np1 <= i__2; ++np1) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = 0.;
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = 0.;
/* L12: */
	    }
	}
    }
L13:
    itypp = *ityp + 1;
    switch (itypp) {
	case 1:  goto L1;
	case 2:  goto L100;
	case 3:  goto L200;
	case 4:  goto L300;
	case 5:  goto L400;
	case 6:  goto L500;
	case 7:  goto L600;
	case 8:  goto L700;
	case 9:  goto L800;
    }

/*     case ityp=0 ,  no symmetries */

/*     case m=0 */

L1:
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = ndo2;
	    for (np1 = 2; np1 <= i__3; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += zv[np1 + i__ * 
			zv_dim1] * ve[i__ + (k * ve_dim2 + 1) * ve_dim1];
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= zv[np1 + i__ * 
			zv_dim1] * we[i__ + (k * we_dim2 + 1) * we_dim1];
/* L15: */
	    }
	}
    }
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__1 = imm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = ndo1;
	    for (np1 = 3; np1 <= i__2; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += zv[np1 + i__ * 
			zv_dim1] * vo[i__ + (k * vo_dim2 + 1) * vo_dim1];
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= zv[np1 + i__ * 
			zv_dim1] * wo[i__ + (k * wo_dim2 + 1) * wo_dim1];
/* L16: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	return 0;
    }
    i__2 = mmax;
    for (mp1 = 2; mp1 <= i__2; ++mp1) {
	m = mp1 - 1;
	mb = m * (*nlat - 1) - m * (m - 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L17;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = ndo1;
		for (np1 = mp1; np1 <= i__4; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2)
			     * vo_dim1] + zw[np1 + mb + i__ * zw_dim1] * we[
			    i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1];
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2)
			     * vo_dim1] - zw[np1 + mb + i__ * zw_dim1] * we[
			    i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1];
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2)
			     * wo_dim1] + zw[np1 + mb + i__ * zw_dim1] * ve[
			    i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1];
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2)
			     * wo_dim1] - zw[np1 + mb + i__ * zw_dim1] * ve[
			    i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1];
/* L23: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L17;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__3 = ndo1;
	    for (np1 = mp1; np1 <= i__3; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += zw[np1 + mb + *
			imid * zw_dim1] * we[*imid + ((mp1 << 1) - 1 + k * 
			we_dim2) * we_dim1];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] -= zw[np1 + mb + *
			imid * zw_dim1] * we[*imid + ((mp1 << 1) - 2 + k * 
			we_dim2) * we_dim1];
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] += zw[np1 + mb + *
			imid * zw_dim1] * ve[*imid + ((mp1 << 1) - 1 + k * 
			ve_dim2) * ve_dim1];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= zw[np1 + mb + *
			imid * zw_dim1] * ve[*imid + ((mp1 << 1) - 2 + k * 
			ve_dim2) * ve_dim1];
/* L24: */
	    }
	}
L17:
	if (mp2 > ndo2) {
	    goto L20;
	}
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__4 = imm1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__1 = ndo2;
		for (np1 = mp2; np1 <= i__1; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2)
			     * ve_dim1] + zw[np1 + mb + i__ * zw_dim1] * wo[
			    i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1];
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2)
			     * ve_dim1] - zw[np1 + mb + i__ * zw_dim1] * wo[
			    i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1];
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * we[i__ + ((mp1 << 1) - 2 + k * we_dim2)
			     * we_dim1] + zw[np1 + mb + i__ * zw_dim1] * vo[
			    i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1];
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * we[i__ + ((mp1 << 1) - 1 + k * we_dim2)
			     * we_dim1] - zw[np1 + mb + i__ * zw_dim1] * vo[
			    i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1];
/* L21: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L20;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__4 = ndo2;
	    for (np1 = mp2; np1 <= i__4; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += zv[np1 + mb + *
			imid * zv_dim1] * ve[*imid + ((mp1 << 1) - 2 + k * 
			ve_dim2) * ve_dim1];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] += zv[np1 + mb + *
			imid * zv_dim1] * ve[*imid + ((mp1 << 1) - 1 + k * 
			ve_dim2) * ve_dim1];
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] -= zv[np1 + mb + *
			imid * zv_dim1] * we[*imid + ((mp1 << 1) - 2 + k * 
			we_dim2) * we_dim1];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= zv[np1 + mb + *
			imid * zv_dim1] * we[*imid + ((mp1 << 1) - 1 + k * 
			we_dim2) * we_dim1];
/* L22: */
	    }
	}
L20:
	;
    }
    return 0;

/*     case ityp=1 ,  no symmetries but cr and ci equal zero */

/*     case m=0 */

L100:
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__4 = *imid;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    i__1 = ndo2;
	    for (np1 = 2; np1 <= i__1; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += zv[np1 + i__ * 
			zv_dim1] * ve[i__ + (k * ve_dim2 + 1) * ve_dim1];
/* L115: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__4 = imm1;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    i__2 = ndo1;
	    for (np1 = 3; np1 <= i__2; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += zv[np1 + i__ * 
			zv_dim1] * vo[i__ + (k * vo_dim2 + 1) * vo_dim1];
/* L116: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	return 0;
    }
    i__2 = mmax;
    for (mp1 = 2; mp1 <= i__2; ++mp1) {
	m = mp1 - 1;
	mb = m * (*nlat - 1) - m * (m - 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L117;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__3 = ndo1;
		for (np1 = mp1; np1 <= i__3; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2)
			     * vo_dim1] + zw[np1 + mb + i__ * zw_dim1] * we[
			    i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1];
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2)
			     * vo_dim1] - zw[np1 + mb + i__ * zw_dim1] * we[
			    i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1];
/* L123: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L117;
	}
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += zw[np1 + mb + *
			imid * zw_dim1] * we[*imid + ((mp1 << 1) - 1 + k * 
			we_dim2) * we_dim1];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] -= zw[np1 + mb + *
			imid * zw_dim1] * we[*imid + ((mp1 << 1) - 2 + k * 
			we_dim2) * we_dim1];
/* L124: */
	    }
	}
L117:
	if (mp2 > ndo2) {
	    goto L120;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = ndo2;
		for (np1 = mp2; np1 <= i__4; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2)
			     * ve_dim1] + zw[np1 + mb + i__ * zw_dim1] * wo[
			    i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1];
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2)
			     * ve_dim1] - zw[np1 + mb + i__ * zw_dim1] * wo[
			    i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1];
/* L121: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L120;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__3 = ndo2;
	    for (np1 = mp2; np1 <= i__3; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += zv[np1 + mb + *
			imid * zv_dim1] * ve[*imid + ((mp1 << 1) - 2 + k * 
			ve_dim2) * ve_dim1];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] += zv[np1 + mb + *
			imid * zv_dim1] * ve[*imid + ((mp1 << 1) - 1 + k * 
			ve_dim2) * ve_dim1];
/* L122: */
	    }
	}
L120:
	;
    }
    return 0;

/*     case ityp=2 ,  no symmetries but br and bi equal zero */

/*     case m=0 */

L200:
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = ndo2;
	    for (np1 = 2; np1 <= i__4; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= zv[np1 + i__ * 
			zv_dim1] * we[i__ + (k * we_dim2 + 1) * we_dim1];
/* L215: */
	    }
	}
    }
    i__4 = *nt;
    for (k = 1; k <= i__4; ++k) {
	i__3 = imm1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__2 = ndo1;
	    for (np1 = 3; np1 <= i__2; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= zv[np1 + i__ * 
			zv_dim1] * wo[i__ + (k * wo_dim2 + 1) * wo_dim1];
/* L216: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	return 0;
    }
    i__2 = mmax;
    for (mp1 = 2; mp1 <= i__2; ++mp1) {
	m = mp1 - 1;
	mb = m * (*nlat - 1) - m * (m - 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L217;
	}
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__4 = imm1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__1 = ndo1;
		for (np1 = mp1; np1 <= i__1; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2)
			     * wo_dim1] + zw[np1 + mb + i__ * zw_dim1] * ve[
			    i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1];
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2)
			     * wo_dim1] - zw[np1 + mb + i__ * zw_dim1] * ve[
			    i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1];
/* L223: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L217;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__4 = ndo1;
	    for (np1 = mp1; np1 <= i__4; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] += zw[np1 + mb + *
			imid * zw_dim1] * ve[*imid + ((mp1 << 1) - 1 + k * 
			ve_dim2) * ve_dim1];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= zw[np1 + mb + *
			imid * zw_dim1] * ve[*imid + ((mp1 << 1) - 2 + k * 
			ve_dim2) * ve_dim1];
/* L224: */
	    }
	}
L217:
	if (mp2 > ndo2) {
	    goto L220;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__3 = ndo2;
		for (np1 = mp2; np1 <= i__3; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * we[i__ + ((mp1 << 1) - 2 + k * we_dim2)
			     * we_dim1] + zw[np1 + mb + i__ * zw_dim1] * vo[
			    i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1];
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * we[i__ + ((mp1 << 1) - 1 + k * we_dim2)
			     * we_dim1] - zw[np1 + mb + i__ * zw_dim1] * vo[
			    i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1];
/* L221: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L220;
	}
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__1 = ndo2;
	    for (np1 = mp2; np1 <= i__1; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] -= zv[np1 + mb + *
			imid * zv_dim1] * we[*imid + ((mp1 << 1) - 2 + k * 
			we_dim2) * we_dim1];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= zv[np1 + mb + *
			imid * zv_dim1] * we[*imid + ((mp1 << 1) - 1 + k * 
			we_dim2) * we_dim1];
/* L222: */
	    }
	}
L220:
	;
    }
    return 0;

/*     case ityp=3 ,  v even , w odd */

/*     case m=0 */

L300:
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = ndo2;
	    for (np1 = 2; np1 <= i__3; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += zv[np1 + i__ * 
			zv_dim1] * ve[i__ + (k * ve_dim2 + 1) * ve_dim1];
/* L315: */
	    }
	}
    }
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__1 = imm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = ndo1;
	    for (np1 = 3; np1 <= i__2; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= zv[np1 + i__ * 
			zv_dim1] * wo[i__ + (k * wo_dim2 + 1) * wo_dim1];
/* L316: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	return 0;
    }
    i__2 = mmax;
    for (mp1 = 2; mp1 <= i__2; ++mp1) {
	m = mp1 - 1;
	mb = m * (*nlat - 1) - m * (m - 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L317;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = ndo1;
		for (np1 = mp1; np1 <= i__4; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2)
			     * wo_dim1] + zw[np1 + mb + i__ * zw_dim1] * ve[
			    i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1];
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2)
			     * wo_dim1] - zw[np1 + mb + i__ * zw_dim1] * ve[
			    i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1];
/* L323: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L317;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__3 = ndo1;
	    for (np1 = mp1; np1 <= i__3; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] += zw[np1 + mb + *
			imid * zw_dim1] * ve[*imid + ((mp1 << 1) - 1 + k * 
			ve_dim2) * ve_dim1];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= zw[np1 + mb + *
			imid * zw_dim1] * ve[*imid + ((mp1 << 1) - 2 + k * 
			ve_dim2) * ve_dim1];
/* L324: */
	    }
	}
L317:
	if (mp2 > ndo2) {
	    goto L320;
	}
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__4 = imm1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__1 = ndo2;
		for (np1 = mp2; np1 <= i__1; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2)
			     * ve_dim1] + zw[np1 + mb + i__ * zw_dim1] * wo[
			    i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1];
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2)
			     * ve_dim1] - zw[np1 + mb + i__ * zw_dim1] * wo[
			    i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1];
/* L321: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L320;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__4 = ndo2;
	    for (np1 = mp2; np1 <= i__4; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += zv[np1 + mb + *
			imid * zv_dim1] * ve[*imid + ((mp1 << 1) - 2 + k * 
			ve_dim2) * ve_dim1];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] += zv[np1 + mb + *
			imid * zv_dim1] * ve[*imid + ((mp1 << 1) - 1 + k * 
			ve_dim2) * ve_dim1];
/* L322: */
	    }
	}
L320:
	;
    }
    return 0;

/*     case ityp=4 ,  v even, w odd, and cr and ci equal 0. */

/*     case m=0 */

L400:
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__4 = *imid;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    i__1 = ndo2;
	    for (np1 = 2; np1 <= i__1; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += zv[np1 + i__ * 
			zv_dim1] * ve[i__ + (k * ve_dim2 + 1) * ve_dim1];
/* L415: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	return 0;
    }
    i__1 = mmax;
    for (mp1 = 2; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	mb = m * (*nlat - 1) - m * (m - 1) / 2;
	mp2 = mp1 + 1;
	if (mp2 > ndo2) {
	    goto L420;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__2 = imm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = ndo2;
		for (np1 = mp2; np1 <= i__3; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2)
			     * ve_dim1] + zw[np1 + mb + i__ * zw_dim1] * wo[
			    i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1];
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2)
			     * ve_dim1] - zw[np1 + mb + i__ * zw_dim1] * wo[
			    i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1];
/* L421: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L420;
	}
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__2 = ndo2;
	    for (np1 = mp2; np1 <= i__2; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += zv[np1 + mb + *
			imid * zv_dim1] * ve[*imid + ((mp1 << 1) - 2 + k * 
			ve_dim2) * ve_dim1];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] += zv[np1 + mb + *
			imid * zv_dim1] * ve[*imid + ((mp1 << 1) - 1 + k * 
			ve_dim2) * ve_dim1];
/* L422: */
	    }
	}
L420:
	;
    }
    return 0;

/*     case ityp=5   v even, w odd, and br and bi equal zero */

/*     case m=0 */

L500:
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = imm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = ndo1;
	    for (np1 = 3; np1 <= i__3; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= zv[np1 + i__ * 
			zv_dim1] * wo[i__ + (k * wo_dim2 + 1) * wo_dim1];
/* L516: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	return 0;
    }
    i__3 = mmax;
    for (mp1 = 2; mp1 <= i__3; ++mp1) {
	m = mp1 - 1;
	mb = m * (*nlat - 1) - m * (m - 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L520;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__4 = ndo1;
		for (np1 = mp1; np1 <= i__4; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2)
			     * wo_dim1] + zw[np1 + mb + i__ * zw_dim1] * ve[
			    i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1];
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2)
			     * wo_dim1] - zw[np1 + mb + i__ * zw_dim1] * ve[
			    i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1];
/* L523: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L520;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] += zw[np1 + mb + *
			imid * zw_dim1] * ve[*imid + ((mp1 << 1) - 1 + k * 
			ve_dim2) * ve_dim1];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= zw[np1 + mb + *
			imid * zw_dim1] * ve[*imid + ((mp1 << 1) - 2 + k * 
			ve_dim2) * ve_dim1];
/* L524: */
	    }
	}
L520:
	;
    }
    return 0;

/*     case ityp=6 ,  v odd , w even */

/*     case m=0 */

L600:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__4 = ndo2;
	    for (np1 = 2; np1 <= i__4; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= zv[np1 + i__ * 
			zv_dim1] * we[i__ + (k * we_dim2 + 1) * we_dim1];
/* L615: */
	    }
	}
    }
    i__4 = *nt;
    for (k = 1; k <= i__4; ++k) {
	i__1 = imm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = ndo1;
	    for (np1 = 3; np1 <= i__3; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += zv[np1 + i__ * 
			zv_dim1] * vo[i__ + (k * vo_dim2 + 1) * vo_dim1];
/* L616: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	return 0;
    }
    i__3 = mmax;
    for (mp1 = 2; mp1 <= i__3; ++mp1) {
	m = mp1 - 1;
	mb = m * (*nlat - 1) - m * (m - 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L617;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__4 = imm1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__2 = ndo1;
		for (np1 = mp1; np1 <= i__2; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2)
			     * vo_dim1] + zw[np1 + mb + i__ * zw_dim1] * we[
			    i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1];
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2)
			     * vo_dim1] - zw[np1 + mb + i__ * zw_dim1] * we[
			    i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1];
/* L623: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L617;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__4 = ndo1;
	    for (np1 = mp1; np1 <= i__4; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += zw[np1 + mb + *
			imid * zw_dim1] * we[*imid + ((mp1 << 1) - 1 + k * 
			we_dim2) * we_dim1];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] -= zw[np1 + mb + *
			imid * zw_dim1] * we[*imid + ((mp1 << 1) - 2 + k * 
			we_dim2) * we_dim1];
/* L624: */
	    }
	}
L617:
	if (mp2 > ndo2) {
	    goto L620;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__2 = imm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = ndo2;
		for (np1 = mp2; np1 <= i__1; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * we[i__ + ((mp1 << 1) - 2 + k * we_dim2)
			     * we_dim1] + zw[np1 + mb + i__ * zw_dim1] * vo[
			    i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1];
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * we[i__ + ((mp1 << 1) - 1 + k * we_dim2)
			     * we_dim1] - zw[np1 + mb + i__ * zw_dim1] * vo[
			    i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1];
/* L621: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L620;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = ndo2;
	    for (np1 = mp2; np1 <= i__2; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] -= zv[np1 + mb + *
			imid * zv_dim1] * we[*imid + ((mp1 << 1) - 2 + k * 
			we_dim2) * we_dim1];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= zv[np1 + mb + *
			imid * zv_dim1] * we[*imid + ((mp1 << 1) - 1 + k * 
			we_dim2) * we_dim1];
/* L622: */
	    }
	}
L620:
	;
    }
    return 0;

/*     case ityp=7   v odd, w even, and cr and ci equal zero */

/*     case m=0 */

L700:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = imm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__1 = ndo1;
	    for (np1 = 3; np1 <= i__1; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += zv[np1 + i__ * 
			zv_dim1] * vo[i__ + (k * vo_dim2 + 1) * vo_dim1];
/* L716: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	return 0;
    }
    i__1 = mmax;
    for (mp1 = 2; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	mb = m * (*nlat - 1) - m * (m - 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L720;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		i__4 = ndo1;
		for (np1 = mp1; np1 <= i__4; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2)
			     * vo_dim1] + zw[np1 + mb + i__ * zw_dim1] * we[
			    i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1];
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + zv[np1 + mb + i__ * 
			    zv_dim1] * vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2)
			     * vo_dim1] - zw[np1 + mb + i__ * zw_dim1] * we[
			    i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1];
/* L723: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L720;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__3 = ndo1;
	    for (np1 = mp1; np1 <= i__3; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += zw[np1 + mb + *
			imid * zw_dim1] * we[*imid + ((mp1 << 1) - 1 + k * 
			we_dim2) * we_dim1];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] -= zw[np1 + mb + *
			imid * zw_dim1] * we[*imid + ((mp1 << 1) - 2 + k * 
			we_dim2) * we_dim1];
/* L724: */
	    }
	}
L720:
	;
    }
    return 0;

/*     case ityp=8   v odd, w even, and both br and bi equal zero */

/*     case m=0 */

L800:
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__4 = ndo2;
	    for (np1 = 2; np1 <= i__4; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= zv[np1 + i__ * 
			zv_dim1] * we[i__ + (k * we_dim2 + 1) * we_dim1];
/* L815: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	return 0;
    }
    i__4 = mmax;
    for (mp1 = 2; mp1 <= i__4; ++mp1) {
	m = mp1 - 1;
	mb = m * (*nlat - 1) - m * (m - 1) / 2;
	mp2 = mp1 + 1;
	if (mp2 > ndo2) {
	    goto L820;
	}
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = ndo2;
		for (np1 = mp2; np1 <= i__2; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * we[i__ + ((mp1 << 1) - 2 + k * we_dim2)
			     * we_dim1] + zw[np1 + mb + i__ * zw_dim1] * vo[
			    i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1];
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - zv[np1 + mb + i__ * 
			    zv_dim1] * we[i__ + ((mp1 << 1) - 1 + k * we_dim2)
			     * we_dim1] - zw[np1 + mb + i__ * zw_dim1] * vo[
			    i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1];
/* L821: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L820;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo2;
	    for (np1 = mp2; np1 <= i__1; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] -= zv[np1 + mb + *
			imid * zv_dim1] * we[*imid + ((mp1 << 1) - 2 + k * 
			we_dim2) * we_dim1];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= zv[np1 + mb + *
			imid * zv_dim1] * we[*imid + ((mp1 << 1) - 1 + k * 
			we_dim2) * we_dim1];
/* L822: */
	    }
	}
L820:
	;
    }
    return 0;
} /* vhaes1_ */


/*     dwork must be of length at least 2*(nlat+1) */

/* Subroutine */ int vhaesi_(integer *nlat, integer *nlon, doublereal *wvhaes,
	 integer *lvhaes, doublereal *work, integer *lwork, doublereal *dwork,
	 integer *ldwork, integer *ierror)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer iw1, idz;
    extern /* Subroutine */ int vea1_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *);
    integer labc, imid, mmax, lzimn;
    extern /* Subroutine */ int hrffti_(integer *, doublereal *);

    /* Parameter adjustments */
    --wvhaes;
    --work;
    --dwork;

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
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
    imid = (*nlat + 1) / 2;
    lzimn = imid * mmax * (*nlat + *nlat - mmax + 1) / 2;
    if (*lvhaes < lzimn + lzimn + *nlon + 15) {
	return 0;
    }
    *ierror = 4;
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) * 3 / 2;
    if (*lwork < *nlat * 5 * imid + labc) {
	return 0;
    }
    *ierror = 5;
    if (*ldwork < *nlat + 1 << 1) {
	return 0;
    }
    *ierror = 0;
    iw1 = *nlat * 3 * imid + 1;
    idz = mmax * (*nlat + *nlat - mmax + 1) / 2;
    vea1_(nlat, nlon, &imid, &wvhaes[1], &wvhaes[lzimn + 1], &idz, &work[1], &
	    work[iw1], &dwork[1]);
    hrffti_(nlon, &wvhaes[(lzimn << 1) + 1]);
    return 0;
} /* vhaesi_ */

/* Subroutine */ int vea1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *zv, doublereal *zw, integer *idz, doublereal *zin, 
	doublereal *wzvin, doublereal *dwork)
{
    /* System generated locals */
    integer zv_dim1, zv_offset, zw_dim1, zw_offset, zin_dim1, zin_dim2, 
	    zin_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, m, i3, mn, mp1, np1, mmax;
    extern /* Subroutine */ int zvin_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *), zwin_(integer *
	    , integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *), zvinit_(integer *, integer *, doublereal *, 
	    doublereal *), zwinit_(integer *, integer *, doublereal *, 
	    doublereal *);

    /* Parameter adjustments */
    zin_dim1 = *imid;
    zin_dim2 = *nlat;
    zin_offset = 1 + zin_dim1 * (1 + zin_dim2);
    zin -= zin_offset;
    zw_dim1 = *idz;
    zw_offset = 1 + zw_dim1;
    zw -= zw_offset;
    zv_dim1 = *idz;
    zv_offset = 1 + zv_dim1;
    zv -= zv_offset;
    --wzvin;
    --dwork;

    /* Function Body */
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
    zvinit_(nlat, nlon, &wzvin[1], &dwork[1]);
    i__1 = mmax;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	zvin_(&c__0, nlat, nlon, &m, &zin[zin_offset], &i3, &wzvin[1]);
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    mn = m * (*nlat - 1) - m * (m - 1) / 2 + np1;
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		zv[mn + i__ * zv_dim1] = zin[i__ + (np1 + i3 * zin_dim2) * 
			zin_dim1];
/* L33: */
	    }
	}
    }
    zwinit_(nlat, nlon, &wzvin[1], &dwork[1]);
    i__3 = mmax;
    for (mp1 = 1; mp1 <= i__3; ++mp1) {
	m = mp1 - 1;
	zwin_(&c__0, nlat, nlon, &m, &zin[zin_offset], &i3, &wzvin[1]);
	i__2 = *nlat;
	for (np1 = mp1; np1 <= i__2; ++np1) {
	    mn = m * (*nlat - 1) - m * (m - 1) / 2 + np1;
	    i__1 = *imid;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		zw[mn + i__ * zw_dim1] = zin[i__ + (np1 + i3 * zin_dim2) * 
			zin_dim1];
/* L34: */
	    }
	}
    }
    return 0;
} /* vea1_ */

