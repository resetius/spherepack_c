/* vhagc.f -- translated by f2c (version 20061008).
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
static integer c__1 = 1;
static integer c__2 = 2;


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



/* ... file vhagc.f */

/*     this file contains code and documentation for subroutines */
/*     vhagc and vhagci */

/* ... files which must be loaded with vhagc.f */

/*     sphcom.f, hrfft.f, gaqd.f */


/*     subroutine vhagc(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci, */
/*    +                 mdab,ndab,wvhagc,lvhagc,work,lwork,ierror) */

/*     subroutine vhagc performs the vector spherical harmonic analysis */
/*     on the vector field (v,w) and stores the result in the arrays */
/*     br,bi,cr, and ci. v(i,j) and w(i,j) are the colatitudinal */
/*     (measured from the north pole) and east longitudinal components */
/*     respectively, located at the gaussian colatitude point theta(i) */
/*     and longitude phi(j) = (j-1)*2*pi/nlon. the spectral */
/*     representation of (v,w) is given at output parameters v,w in */
/*     subroutine vhsec. */

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


/*     nt     the number of analyses.  in the program that calls vhagc, */
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
/*            the program that calls vhagc. if ityp .le. 2 then idvw */
/*            must be at least nlat.  if ityp .gt. 2 and nlat is */
/*            even then idvw must be at least nlat/2. if ityp .gt. 2 */
/*            and nlat is odd then idvw must be at least (nlat+1)/2. */

/*     jdvw   the second dimension of the arrays v,w as it appears in */
/*            the program that calls vhagc. jdvw must be at least nlon. */

/*     mdab   the first dimension of the arrays br,bi,cr, and ci as it */
/*            appears in the program that calls vhagc. mdab must be at */
/*            least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndab   the second dimension of the arrays br,bi,cr, and ci as it */
/*            appears in the program that calls vhagc. ndab must be at */
/*            least nlat. */

/*     wvhagc an array which must be initialized by subroutine vhagci. */
/*            once initialized, wvhagc can be used repeatedly by vhagc */
/*            as long as nlon and nlat remain unchanged.  wvhagc must */
/*            not be altered between calls of vhagc. */

/*     lvhagc the dimension of the array wvhagc as it appears in the */
/*            program that calls vhagc. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lvhagc must be at least */

/*               4*nlat*l2+3*max0(l1-2,0)*(2*nlat-l1-1)+nlon+l2+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls vhagc. define */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            if ityp .le. 2 then lwork must be at least */

/*               2*nlat*(2*nlon*nt+3*l2) */

/*            if ityp .gt. 2 then lwork must be at least */

/*               2*l2*(2*nlon*nt+3*nlat) */



/*     ************************************************************** */

/*     output parameters */

/*     br,bi  two or three dimensional arrays (see input parameter nt) */
/*     cr,ci  that contain the vector spherical harmonic coefficients */
/*            in the spectral representation of v(i,j) and w(i,j) given */
/*            in the discription of subroutine vhsec. br(mp1,np1), */
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
/*            = 9  error in the specification of lvhagc */
/*            = 10 error in the specification of lwork */

/* **************************************************************** */

/*     subroutine vhagci(nlat,nlon,wvhagc,lvhagc,dwork,ldwork,ierror) */

/*     subroutine vhagci initializes the array wvhagc which can then be */
/*     used repeatedly by subroutine vhagc until nlat or nlon is changed. */

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

/*     lvhagc the dimension of the array wvhagc as it appears in the */
/*            program that calls vhagci.  define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lvhagc must be at least */

/*               4*nlat*l2+3*max0(l1-2,0)*(2*nlat-l1-1)+nlon+l2+15 */


/*     dwork  a double precision work array that does not need to be saved */

/*     ldwork the dimension of the array dwork as it appears in the */
/*            program that calls vhagci. ldwork must be at least */

/*               2*nlat*(nlat+1)+1 */


/*     ************************************************************** */

/*     output parameters */

/*     wvhagc an array which is initialized for use by subroutine vhagc. */
/*            once initialized, wvhagc can be used repeatedly by vhagc */
/*            as long as nlat and nlon remain unchanged.  wvhagc must not */
/*            be altered between calls of vhagc. */


/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of lvhagc */
/*            = 4  error in the specification of lwork */

/* Subroutine */ int vhagc_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, real *v, real *w, integer *idvw, integer *jdvw, real *br,
	 real *bi, real *cr, real *ci, integer *mdab, integer *ndab, real *
	wvhagc, integer *lvhagc, real *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, br_dim1, 
	    br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, cr_dim2,
	     cr_offset, ci_dim1, ci_dim2, ci_offset, i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, iw3, iw4, iw5, jw1, jw2, jw3, idv, lnl, ist, 
	    lzz1, labc, imid, mmax;
    extern /* Subroutine */ int vhagc1_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, real *, real *, 
	    integer *, integer *, real *, real *, real *, real *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *);
    static integer lwzvin;

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
    --wvhagc;
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
    lzz1 = (*nlat << 1) * imid;
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) * 3 / 2;
    if (*lvhagc < (lzz1 + labc << 1) + *nlon + imid + 15) {
	return 0;
    }
    *ierror = 10;
    if (*ityp <= 2 && *lwork < *nlat * ((*nlon << 2) * *nt + imid * 6)) {
	return 0;
    }
    if (*ityp > 2 && *lwork < imid * ((*nlon << 2) * *nt + *nlat * 6)) {
	return 0;
    }
    *ierror = 0;
    idv = *nlat;
    if (*ityp > 2) {
	idv = imid;
    }
    lnl = *nt * idv * *nlon;
    ist = 0;
    if (*ityp <= 2) {
	ist = imid;
    }
    iw1 = ist + 1;
    iw2 = lnl + 1;
    iw3 = iw2 + ist;
    iw4 = iw2 + lnl;
    iw5 = iw4 + imid * 3 * *nlat;
    lwzvin = lzz1 + labc;
    jw1 = (*nlat + 1) / 2 + 1;
    jw2 = jw1 + lwzvin;
    jw3 = jw2 + lwzvin;
    vhagc1_(nlat, nlon, ityp, nt, &imid, idvw, jdvw, &v[v_offset], &w[
	    w_offset], mdab, ndab, &br[br_offset], &bi[bi_offset], &cr[
	    cr_offset], &ci[ci_offset], &idv, &work[1], &work[iw1], &work[iw2]
	    , &work[iw3], &work[iw4], &work[iw5], &wvhagc[1], &wvhagc[jw1], &
	    wvhagc[jw2], &wvhagc[jw3]);
    return 0;
} /* vhagc_ */

/* Subroutine */ int vhagc1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, real *v, 
	real *w, integer *mdab, integer *ndab, real *br, real *bi, real *cr, 
	real *ci, integer *idv, real *ve, real *vo, real *we, real *wo, real *
	vb, real *wb, real *wts, real *wvbin, real *wwbin, real *wrfft)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, br_dim1, 
	    br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, cr_dim2,
	     cr_offset, ci_dim1, ci_dim2, ci_offset, ve_dim1, ve_dim2, 
	    ve_offset, vo_dim1, vo_dim2, vo_offset, we_dim1, we_dim2, 
	    we_offset, wo_dim1, wo_dim2, wo_offset, vb_dim1, vb_dim2, 
	    vb_offset, wb_dim1, wb_dim2, wb_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, m, iv, iw;
    static real tv, tw;
    static integer mp1, np1, mp2;
    static real fsn, tsn;
    static integer ndo1, ndo2, imm1, nlp1;
    static real tve1, tve2, twe1, twe2, tvo1, tvo2, two1, two2;
    static integer mlat;
    extern /* Subroutine */ int vbin_(integer *, integer *, integer *, 
	    integer *, real *, integer *, real *), wbin_(integer *, integer *,
	     integer *, integer *, real *, integer *, real *);
    static integer mmax, mlon, itypp;
    extern /* Subroutine */ int hrfftf_(integer *, integer *, real *, integer 
	    *, real *, real *);

    /* Parameter adjustments */
    wb_dim1 = *imid;
    wb_dim2 = *nlat;
    wb_offset = 1 + wb_dim1 * (1 + wb_dim2);
    wb -= wb_offset;
    vb_dim1 = *imid;
    vb_dim2 = *nlat;
    vb_offset = 1 + vb_dim1 * (1 + vb_dim2);
    vb -= vb_offset;
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
    --wts;
    --wvbin;
    --wwbin;
    --wrfft;

    /* Function Body */
    nlp1 = *nlat + 1;
    tsn = 2.f / *nlon;
    fsn = 4.f / *nlon;
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
		1], &vb[vb_offset]);
	hrfftf_(idv, nlon, &we[(k * we_dim2 + 1) * we_dim1 + 1], idv, &wrfft[
		1], &vb[vb_offset]);
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
		br[mp1 + (np1 + k * br_dim2) * br_dim1] = 0.f;
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = 0.f;
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
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = 0.f;
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = 0.f;
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

L1:
    vbin_(&c__0, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m=0 */

    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tv = ve[i__ + (k * ve_dim2 + 1) * ve_dim1] * wts[i__];
	    tw = we[i__ + (k * we_dim2 + 1) * we_dim1] * wts[i__];
	    i__3 = ndo2;
	    for (np1 = 2; np1 <= i__3; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tv;
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tw;
/* L10015: */
	    }
/* L1015: */
	}
/* L15: */
    }
    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__1 = imm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tv = vo[i__ + (k * vo_dim2 + 1) * vo_dim1] * wts[i__];
	    tw = wo[i__ + (k * wo_dim2 + 1) * wo_dim1] * wts[i__];
	    i__3 = ndo1;
	    for (np1 = 3; np1 <= i__3; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tv;
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tw;
/* L10016: */
	    }
/* L1016: */
	}
/* L16: */
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	return 0;
    }
    i__2 = mmax;
    for (mp1 = 2; mp1 <= i__2; ++mp1) {
	m = mp1 - 1;
	mp2 = mp1 + 1;
	vbin_(&c__0, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__0, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L17;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {

/*     set temps to optimize quadrature */

		tvo1 = vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		tvo2 = vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		tve1 = ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		tve2 = ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		two1 = wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		two2 = wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		twe1 = we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] * 
			wts[i__];
		twe2 = we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] * 
			wts[i__];
		i__4 = ndo1;
		for (np1 = mp1; np1 <= i__4; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tvo2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * twe1;
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tvo1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * twe2;
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * two2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tve1;
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * two1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tve2;
/* L10023: */
		}
/* L1023: */
	    }
/* L23: */
	}
	if (mlat == 0) {
	    goto L17;
	}
	i__ = *imid;
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__3 = ndo1;
	    for (np1 = mp1; np1 <= i__3; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * we[i__ + ((mp1 << 1) - 1 + k *
			 we_dim2) * we_dim1] * wts[i__];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] -= wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * we[i__ + ((mp1 << 1) - 2 + k *
			 we_dim2) * we_dim1] * wts[i__];
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] += wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * ve[i__ + ((mp1 << 1) - 1 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * ve[i__ + ((mp1 << 1) - 2 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
/* L1024: */
	    }
/* L24: */
	}
L17:
	if (mp2 > ndo2) {
	    goto L20;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		tvo1 = vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		tvo2 = vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		tve1 = ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		tve2 = ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		two1 = wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		two2 = wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		twe1 = we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] * 
			wts[i__];
		twe2 = we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] * 
			wts[i__];
		i__4 = ndo2;
		for (np1 = mp2; np1 <= i__4; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tve2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * two1;
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tve1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * two2;
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * twe2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tvo1;
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * twe1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tvo2;
/* L10021: */
		}
/* L1021: */
	    }
/* L21: */
	}
	if (mlat == 0) {
	    goto L20;
	}
	i__ = *imid;
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__3 = ndo2;
	    for (np1 = mp2; np1 <= i__3; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * ve[i__ + ((mp1 << 1) - 2 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] += vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * ve[i__ + ((mp1 << 1) - 1 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] -= vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * we[i__ + ((mp1 << 1) - 2 + k *
			 we_dim2) * we_dim1] * wts[i__];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * we[i__ + ((mp1 << 1) - 1 + k *
			 we_dim2) * we_dim1] * wts[i__];
/* L1022: */
	    }
/* L22: */
	}
L20:
	;
    }
    return 0;

/*     case ityp=1 ,  no symmetries but cr and ci equal zero */

L100:
    vbin_(&c__0, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m=0 */

    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tv = ve[i__ + (k * ve_dim2 + 1) * ve_dim1] * wts[i__];
	    i__3 = ndo2;
	    for (np1 = 2; np1 <= i__3; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tv;
/* L115: */
	    }
	}
    }
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__1 = imm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tv = vo[i__ + (k * vo_dim2 + 1) * vo_dim1] * wts[i__];
	    i__2 = ndo1;
	    for (np1 = 3; np1 <= i__2; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tv;
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
	mp2 = mp1 + 1;
	vbin_(&c__0, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__0, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L117;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		tvo1 = vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		tvo2 = vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		twe1 = we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] * 
			wts[i__];
		twe2 = we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] * 
			wts[i__];
		i__4 = ndo1;
		for (np1 = mp1; np1 <= i__4; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tvo2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * twe1;
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tvo1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * twe2;
/* L123: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L117;
	}
	i__ = *imid;
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__3 = ndo1;
	    for (np1 = mp1; np1 <= i__3; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * we[i__ + ((mp1 << 1) - 1 + k *
			 we_dim2) * we_dim1] * wts[i__];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] -= wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * we[i__ + ((mp1 << 1) - 2 + k *
			 we_dim2) * we_dim1] * wts[i__];
/* L124: */
	    }
	}
L117:
	if (mp2 > ndo2) {
	    goto L120;
	}
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__4 = imm1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		tve1 = ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		tve2 = ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		two1 = wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		two2 = wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		i__1 = ndo2;
		for (np1 = mp2; np1 <= i__1; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tve2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * two1;
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tve1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * two2;
/* L121: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L120;
	}
	i__ = *imid;
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__4 = ndo2;
	    for (np1 = mp2; np1 <= i__4; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * ve[i__ + ((mp1 << 1) - 2 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] += vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * ve[i__ + ((mp1 << 1) - 1 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
/* L122: */
	    }
	}
L120:
	;
    }
    return 0;

/*     case ityp=2 ,  no symmetries but br and bi equal zero */

L200:
    vbin_(&c__0, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m=0 */

    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__4 = *imid;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    tw = we[i__ + (k * we_dim2 + 1) * we_dim1] * wts[i__];
	    i__1 = ndo2;
	    for (np1 = 2; np1 <= i__1; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tw;
/* L215: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__4 = imm1;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    tw = wo[i__ + (k * wo_dim2 + 1) * wo_dim1] * wts[i__];
	    i__2 = ndo1;
	    for (np1 = 3; np1 <= i__2; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tw;
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
	mp2 = mp1 + 1;
	vbin_(&c__0, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__0, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L217;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		tve1 = ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		tve2 = ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		two1 = wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		two2 = wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		i__3 = ndo1;
		for (np1 = mp1; np1 <= i__3; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * two2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tve1;
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * two1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tve2;
/* L223: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L217;
	}
	i__ = *imid;
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] += wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * ve[i__ + ((mp1 << 1) - 1 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * ve[i__ + ((mp1 << 1) - 2 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
/* L224: */
	    }
	}
L217:
	if (mp2 > ndo2) {
	    goto L220;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		twe1 = we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] * 
			wts[i__];
		twe2 = we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] * 
			wts[i__];
		tvo1 = vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		tvo2 = vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		i__4 = ndo2;
		for (np1 = mp2; np1 <= i__4; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * twe2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tvo1;
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * twe1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tvo2;
/* L221: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L220;
	}
	i__ = *imid;
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__3 = ndo2;
	    for (np1 = mp2; np1 <= i__3; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] -= vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * we[i__ + ((mp1 << 1) - 2 + k *
			 we_dim2) * we_dim1] * wts[i__];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * we[i__ + ((mp1 << 1) - 1 + k *
			 we_dim2) * we_dim1] * wts[i__];
/* L222: */
	    }
	}
L220:
	;
    }
    return 0;

/*     case ityp=3 ,  v even , w odd */

L300:
    vbin_(&c__0, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m=0 */

    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    tv = ve[i__ + (k * ve_dim2 + 1) * ve_dim1] * wts[i__];
	    i__4 = ndo2;
	    for (np1 = 2; np1 <= i__4; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tv;
/* L315: */
	    }
	}
    }
    i__4 = *nt;
    for (k = 1; k <= i__4; ++k) {
	i__3 = imm1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    tw = wo[i__ + (k * wo_dim2 + 1) * wo_dim1] * wts[i__];
	    i__2 = ndo1;
	    for (np1 = 3; np1 <= i__2; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tw;
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
	mp2 = mp1 + 1;
	vbin_(&c__0, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__0, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L317;
	}
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__4 = imm1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		two1 = wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		two2 = wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		tve1 = ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		tve2 = ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		i__1 = ndo1;
		for (np1 = mp1; np1 <= i__1; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * two2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tve1;
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * two1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tve2;
/* L323: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L317;
	}
	i__ = *imid;
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__4 = ndo1;
	    for (np1 = mp1; np1 <= i__4; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] += wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * ve[i__ + ((mp1 << 1) - 1 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * ve[i__ + ((mp1 << 1) - 2 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
/* L324: */
	    }
	}
L317:
	if (mp2 > ndo2) {
	    goto L320;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		two1 = wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		two2 = wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		tve1 = ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		tve2 = ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		i__3 = ndo2;
		for (np1 = mp2; np1 <= i__3; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tve2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * two1;
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tve1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * two2;
/* L321: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L320;
	}
	i__ = *imid;
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__1 = ndo2;
	    for (np1 = mp2; np1 <= i__1; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * ve[i__ + ((mp1 << 1) - 2 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] += vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * ve[i__ + ((mp1 << 1) - 1 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
/* L322: */
	    }
	}
L320:
	;
    }
    return 0;

/*     case ityp=4 ,  v even, w odd, and cr and ci equal 0. */

L400:
    vbin_(&c__1, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m=0 */

    i__2 = *nt;
    for (k = 1; k <= i__2; ++k) {
	i__1 = *imid;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tv = ve[i__ + (k * ve_dim2 + 1) * ve_dim1] * wts[i__];
	    i__3 = ndo2;
	    for (np1 = 2; np1 <= i__3; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tv;
/* L415: */
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
	mp2 = mp1 + 1;
	vbin_(&c__1, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__1, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp2 > ndo2) {
	    goto L420;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = imm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		two1 = wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		two2 = wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		tve1 = ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		tve2 = ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		i__4 = ndo2;
		for (np1 = mp2; np1 <= i__4; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tve2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * two1;
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tve1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * two2;
/* L421: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L420;
	}
	i__ = *imid;
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__2 = ndo2;
	    for (np1 = mp2; np1 <= i__2; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * ve[i__ + ((mp1 << 1) - 2 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] += vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * ve[i__ + ((mp1 << 1) - 1 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
/* L422: */
	    }
	}
L420:
	;
    }
    return 0;

/*     case ityp=5   v even, w odd, and br and bi equal zero */

L500:
    vbin_(&c__2, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m=0 */

    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = imm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    tw = wo[i__ + (k * wo_dim2 + 1) * wo_dim1] * wts[i__];
	    i__4 = ndo1;
	    for (np1 = 3; np1 <= i__4; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tw;
/* L516: */
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
	mp2 = mp1 + 1;
	vbin_(&c__2, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__2, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L520;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		two1 = wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		two2 = wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] * 
			wts[i__];
		tve1 = ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		tve2 = ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] * 
			wts[i__];
		i__1 = ndo1;
		for (np1 = mp1; np1 <= i__1; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * two2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tve1;
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * two1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tve2;
/* L523: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L520;
	}
	i__ = *imid;
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__3 = ndo1;
	    for (np1 = mp1; np1 <= i__3; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] += wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * ve[i__ + ((mp1 << 1) - 1 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * ve[i__ + ((mp1 << 1) - 2 + k *
			 ve_dim2) * ve_dim1] * wts[i__];
/* L524: */
	    }
	}
L520:
	;
    }
    return 0;

/*     case ityp=6 ,  v odd , w even */

L600:
    vbin_(&c__0, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m=0 */

    i__4 = *nt;
    for (k = 1; k <= i__4; ++k) {
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    tw = we[i__ + (k * we_dim2 + 1) * we_dim1] * wts[i__];
	    i__1 = ndo2;
	    for (np1 = 2; np1 <= i__1; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tw;
/* L615: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__3 = imm1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    tv = vo[i__ + (k * vo_dim2 + 1) * vo_dim1] * wts[i__];
	    i__4 = ndo1;
	    for (np1 = 3; np1 <= i__4; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tv;
/* L616: */
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
	mp2 = mp1 + 1;
	vbin_(&c__0, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__0, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L617;
	}
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		twe1 = we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] * 
			wts[i__];
		twe2 = we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] * 
			wts[i__];
		tvo1 = vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		tvo2 = vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		i__2 = ndo1;
		for (np1 = mp1; np1 <= i__2; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tvo2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * twe1;
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tvo1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * twe2;
/* L623: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L617;
	}
	i__ = *imid;
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * we[i__ + ((mp1 << 1) - 1 + k *
			 we_dim2) * we_dim1] * wts[i__];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] -= wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * we[i__ + ((mp1 << 1) - 2 + k *
			 we_dim2) * we_dim1] * wts[i__];
/* L624: */
	    }
	}
L617:
	if (mp2 > ndo2) {
	    goto L620;
	}
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__2 = imm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		twe1 = we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] * 
			wts[i__];
		twe2 = we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] * 
			wts[i__];
		tvo1 = vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		tvo2 = vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		i__3 = ndo2;
		for (np1 = mp2; np1 <= i__3; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * twe2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tvo1;
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * twe1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tvo2;
/* L621: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L620;
	}
	i__ = *imid;
	i__3 = *nt;
	for (k = 1; k <= i__3; ++k) {
	    i__2 = ndo2;
	    for (np1 = mp2; np1 <= i__2; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] -= vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * we[i__ + ((mp1 << 1) - 2 + k *
			 we_dim2) * we_dim1] * wts[i__];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * we[i__ + ((mp1 << 1) - 1 + k *
			 we_dim2) * we_dim1] * wts[i__];
/* L622: */
	    }
	}
L620:
	;
    }
    return 0;

/*     case ityp=7   v odd, w even, and cr and ci equal zero */

L700:
    vbin_(&c__2, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m=0 */

    i__4 = *nt;
    for (k = 1; k <= i__4; ++k) {
	i__2 = imm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    tv = vo[i__ + (k * vo_dim2 + 1) * vo_dim1] * wts[i__];
	    i__3 = ndo1;
	    for (np1 = 3; np1 <= i__3; np1 += 2) {
		br[(np1 + k * br_dim2) * br_dim1 + 1] += vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tv;
/* L716: */
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
	mp2 = mp1 + 1;
	vbin_(&c__2, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__2, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L720;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__4 = imm1;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		twe1 = we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] * 
			wts[i__];
		twe2 = we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] * 
			wts[i__];
		tvo1 = vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		tvo2 = vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		i__1 = ndo1;
		for (np1 = mp1; np1 <= i__1; np1 += 2) {
		    br[mp1 + (np1 + k * br_dim2) * br_dim1] = br[mp1 + (np1 + 
			    k * br_dim2) * br_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tvo2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * twe1;
		    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] = bi[mp1 + (np1 + 
			    k * bi_dim2) * bi_dim1] + vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * tvo1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * twe2;
/* L723: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L720;
	}
	i__ = *imid;
	i__1 = *nt;
	for (k = 1; k <= i__1; ++k) {
	    i__4 = ndo1;
	    for (np1 = mp1; np1 <= i__4; np1 += 2) {
		br[mp1 + (np1 + k * br_dim2) * br_dim1] += wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * we[i__ + ((mp1 << 1) - 1 + k *
			 we_dim2) * we_dim1] * wts[i__];
		bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] -= wb[i__ + (np1 + iw 
			* wb_dim2) * wb_dim1] * we[i__ + ((mp1 << 1) - 2 + k *
			 we_dim2) * we_dim1] * wts[i__];
/* L724: */
	    }
	}
L720:
	;
    }
    return 0;

/*     case ityp=8   v odd, w even, and both br and bi equal zero */

L800:
    vbin_(&c__1, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m=0 */

    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__4 = *imid;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    tw = we[i__ + (k * we_dim2 + 1) * we_dim1] * wts[i__];
	    i__1 = ndo2;
	    for (np1 = 2; np1 <= i__1; np1 += 2) {
		cr[(np1 + k * cr_dim2) * cr_dim1 + 1] -= vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1] * tw;
/* L815: */
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
	mp2 = mp1 + 1;
	vbin_(&c__1, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__1, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp2 > ndo2) {
	    goto L820;
	}
	i__4 = *nt;
	for (k = 1; k <= i__4; ++k) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		twe1 = we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] * 
			wts[i__];
		twe2 = we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] * 
			wts[i__];
		tvo1 = vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		tvo2 = vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] * 
			wts[i__];
		i__2 = ndo2;
		for (np1 = mp2; np1 <= i__2; np1 += 2) {
		    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] = cr[mp1 + (np1 + 
			    k * cr_dim2) * cr_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * twe2 + wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tvo1;
		    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] = ci[mp1 + (np1 + 
			    k * ci_dim2) * ci_dim1] - vb[i__ + (np1 + iv * 
			    vb_dim2) * vb_dim1] * twe1 - wb[i__ + (np1 + iw * 
			    wb_dim2) * wb_dim1] * tvo2;
/* L821: */
		}
	    }
	}
	if (mlat == 0) {
	    goto L820;
	}
	i__ = *imid;
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = ndo2;
	    for (np1 = mp2; np1 <= i__3; np1 += 2) {
		cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] -= vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * we[i__ + ((mp1 << 1) - 2 + k *
			 we_dim2) * we_dim1] * wts[i__];
		ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] -= vb[i__ + (np1 + iv 
			* vb_dim2) * vb_dim1] * we[i__ + ((mp1 << 1) - 1 + k *
			 we_dim2) * we_dim1] * wts[i__];
/* L822: */
	    }
	}
L820:
	;
    }
    return 0;
} /* vhagc1_ */

/* Subroutine */ int vhagci_(integer *nlat, integer *nlon, real *wvhagc, 
	integer *lvhagc, doublereal *dwork, integer *ldwork, integer *ierror)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, jw1, jw2, jw3, iw2, iw3, lwk, lzz1, labc;
    extern /* Subroutine */ int gaqd_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    static integer imid, mmax, iwrk;
    extern /* Subroutine */ int hrffti_(integer *, real *), vbgint_(integer *,
	     integer *, doublereal *, real *, doublereal *), wbgint_(integer *
	    , integer *, doublereal *, real *, doublereal *);
    static integer lwvbin;
    extern /* Subroutine */ int setwts_(integer *, doublereal *, real *);

    /* Parameter adjustments */
    --dwork;
    --wvhagc;

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
    imid = (*nlat + 1) / 2;
    lzz1 = (*nlat << 1) * imid;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) * 3 / 2;
    imid = (*nlat + 1) / 2;
    if (*lvhagc < (lzz1 + labc << 1) + *nlon + imid + 15) {
	return 0;
    }
    *ierror = 4;
    if (*ldwork < (*nlat << 1) * (*nlat + 1) + 1) {
	return 0;
    }
    *ierror = 0;

/*     compute gaussian points in first nlat+1 words of dwork */
/*     double precision */

    lwk = *nlat * (*nlat + 2);
    jw1 = 1;
/*     jw2 = jw1+nlat+nlat */
/*     jw3 = jw2+nlat+nlat */
    jw2 = jw1 + *nlat;
    jw3 = jw2 + *nlat;
    gaqd_(nlat, &dwork[jw1], &dwork[jw2], &dwork[jw3], &lwk, ierror);
    imid = (*nlat + 1) / 2;

/*     set first imid words of double precision weights in dwork */
/*     as single precision in first imid words of wvhagc */

    setwts_(&imid, &dwork[*nlat + 1], &wvhagc[1]);

/*     first nlat+1 words of dwork contain  double theta */

/*     iwrk = nlat+2 */
    iwrk = (*nlat + 1) / 2 + 1;
    iw1 = imid + 1;
    vbgint_(nlat, nlon, &dwork[1], &wvhagc[iw1], &dwork[iwrk]);
    lwvbin = lzz1 + labc;
    iw2 = iw1 + lwvbin;
    wbgint_(nlat, nlon, &dwork[1], &wvhagc[iw2], &dwork[iwrk]);
    iw3 = iw2 + lwvbin;
    hrffti_(nlon, &wvhagc[iw3]);
    return 0;
} /* vhagci_ */

/* Subroutine */ int setwts_(integer *imid, doublereal *dwts, real *wts)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;


/*     set first imid =(nlat+1)/2 of double precision weights in dwts */
/*     as single precision in wts */

    /* Parameter adjustments */
    --wts;
    --dwts;

    /* Function Body */
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wts[i__] = dwts[i__];
/* L1: */
    }
    return 0;
} /* setwts_ */

