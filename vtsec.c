/* vtsec.f -- translated by f2c (version 20061008).
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
/*  .                         SPHEREPACK                          . */
/*  .                                                             . */
/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */



/* ... file vtsec.f */

/*     this file includes documentation and code for */
/*     subroutines vtsec and vtseci */

/* ... files which must be loaded with vtsec.f */

/*     sphcom.f, hrfft.f, vhaec.f, vhsec.f */

/*     subroutine vtsec(nlat,nlon,ityp,nt,vt,wt,idvw,jdvw,br,bi,cr,ci, */
/*    +                 mdab,ndab,wvts,lwvts,work,lwork,ierror) */

/*     given the vector harmonic analysis br,bi,cr, and ci (computed */
/*     by subroutine vhaec) of some vector function (v,w), this */
/*     subroutine computes the vector function (vt,wt) which is */
/*     the derivative of (v,w) with respect to colatitude theta. vtsec */
/*     is similar to vhsec except the vector harmonics are replaced by */
/*     their derivative with respect to colatitude with the result that */
/*     (vt,wt) is computed instead of (v,w). vt(i,j) is the derivative */
/*     of the colatitudinal component v(i,j) at the point theta(i) = */
/*     (i-1)*pi/(nlat-1) and longitude phi(j) = (j-1)*2*pi/nlon. the */
/*     spectral representation of (vt,wt) is given below at output */
/*     parameters vt,wt. */

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

/*     ityp   = 0  no symmetries exist about the equator. the synthesis */
/*                 is performed on the entire sphere. i.e. the arrays */
/*                 vt(i,j),wt(i,j) are computed for i=1,...,nlat and */
/*                 j=1,...,nlon. */

/*            = 1  no symmetries exist about the equator however the */
/*                 the coefficients cr and ci are zero. the synthesis */
/*                 is performed on the entire sphere.  i.e. the arrays */
/*                 vt(i,j),wt(i,j) are computed for i=1,...,nlat and */
/*                 j=1,...,nlon. */

/*            = 2  no symmetries exist about the equator however the */
/*                 the coefficients br and bi are zero. the synthesis */
/*                 is performed on the entire sphere.  i.e. the arrays */
/*                 vt(i,j),wt(i,j) are computed for i=1,...,nlat and */
/*                 j=1,...,nlon. */

/*            = 3  vt is odd and wt is even about the equator. the */
/*                 synthesis is performed on the northern hemisphere */
/*                 only.  i.e., if nlat is odd the arrays vt(i,j),wt(i,j) */
/*                 are computed for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even the arrays vt(i,j),wt(i,j) are computed */
/*                 for i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 4  vt is odd and wt is even about the equator and the */
/*                 coefficients cr and ci are zero. the synthesis is */
/*                 performed on the northern hemisphere only. i.e. if */
/*                 nlat is odd the arrays vt(i,j),wt(i,j) are computed */
/*                 for i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the arrays vt(i,j),wt(i,j) are computed for */
/*                 i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 5  vt is odd and wt is even about the equator and the */
/*                 coefficients br and bi are zero. the synthesis is */
/*                 performed on the northern hemisphere only. i.e. if */
/*                 nlat is odd the arrays vt(i,j),wt(i,j) are computed */
/*                 for i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the arrays vt(i,j),wt(i,j) are computed for */
/*                 i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 6  vt is even and wt is odd about the equator. the */
/*                 synthesis is performed on the northern hemisphere */
/*                 only.  i.e., if nlat is odd the arrays vt(i,j),wt(i,j) */
/*                 are computed for i=1,...,(nlat+1)/2 and j=1,...,nlon. */
/*                 if nlat is even the arrays vt(i,j),wt(i,j) are computed */
/*                 for i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 7  vt is even and wt is odd about the equator and the */
/*                 coefficients cr and ci are zero. the synthesis is */
/*                 performed on the northern hemisphere only. i.e. if */
/*                 nlat is odd the arrays vt(i,j),wt(i,j) are computed */
/*                 for i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the arrays vt(i,j),wt(i,j) are computed for */
/*                 i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 8  vt is even and wt is odd about the equator and the */
/*                 coefficients br and bi are zero. the synthesis is */
/*                 performed on the northern hemisphere only. i.e. if */
/*                 nlat is odd the arrays vt(i,j),wt(i,j) are computed */
/*                 for i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the arrays vt(i,j),wt(i,j) are computed for */
/*                 i=1,...,nlat/2 and j=1,...,nlon. */

/*     nt     the number of syntheses.  in the program that calls vtsec, */
/*            the arrays vt,wt,br,bi,cr, and ci can be three dimensional */
/*            in which case multiple syntheses will be performed. */
/*            the third index is the synthesis index which assumes the */
/*            values k=1,...,nt.  for a single synthesis set nt=1. the */
/*            discription of the remaining parameters is simplified */
/*            by assuming that nt=1 or that all the arrays are two */
/*            dimensional. */

/*     idvw   the first dimension of the arrays vt,wt as it appears in */
/*            the program that calls vtsec. if ityp .le. 2 then idvw */
/*            must be at least nlat.  if ityp .gt. 2 and nlat is */
/*            even then idvw must be at least nlat/2. if ityp .gt. 2 */
/*            and nlat is odd then idvw must be at least (nlat+1)/2. */

/*     jdvw   the second dimension of the arrays vt,wt as it appears in */
/*            the program that calls vtsec. jdvw must be at least nlon. */

/*     br,bi  two or three dimensional arrays (see input parameter nt) */
/*     cr,ci  that contain the vector spherical harmonic coefficients */
/*            of (v,w) as computed by subroutine vhaec. */

/*     mdab   the first dimension of the arrays br,bi,cr, and ci as it */
/*            appears in the program that calls vtsec. mdab must be at */
/*            least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndab   the second dimension of the arrays br,bi,cr, and ci as it */
/*            appears in the program that calls vtsec. ndab must be at */
/*            least nlat. */

/*     wvts   an array which must be initialized by subroutine vtseci. */
/*            once initialized, wvts can be used repeatedly by vtsec */
/*            as long as nlon and nlat remain unchanged.  wvts must */
/*            not be altered between calls of vtsec. */

/*     lwvts  the dimension of the array wvts as it appears in the */
/*            program that calls vtsec. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lwvts must be at least */

/*            4*nlat*l2+3*max0(l1-2,0)*(nlat+nlat-l1-1)+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls vtsec. define */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            if ityp .le. 2 then lwork must be at least */

/*                    nlat*(2*nt*nlon+max0(6*l2,nlon)) */

/*            if ityp .gt. 2 then lwork must be at least */

/*                    l2*(2*nt*nlon+max0(6*nlat,nlon)) */

/*     ************************************************************** */

/*     output parameters */

/*     vt,wt  two or three dimensional arrays (see input parameter nt) */
/*            in which the derivative of (v,w) with respect to */
/*            colatitude theta is stored. vt(i,j),wt(i,j) contain the */
/*            derivatives at colatitude theta(i) = (i-1)*pi/(nlat-1) */
/*            and longitude phi(j) = (j-1)*2*pi/nlon. the index ranges */
/*            are defined above at the input parameter ityp. vt and wt */
/*            are computed from the formulas for v and w given in */
/*            subroutine vhsec but with vbar and wbar replaced with */
/*            their derivatives with respect to colatitude. these */
/*            derivatives are denoted by vtbar and wtbar. */

/*   in terms of doublereal variables this expansion takes the form */

/*             for i=1,...,nlat and  j=1,...,nlon */

/*     vt(i,j) = the sum from n=1 to n=nlat-1 of */

/*               .5*br(1,n+1)*vtbar(0,n,theta(i)) */

/*     plus the sum from m=1 to m=mmax-1 of the sum from n=m to */
/*     n=nlat-1 of the doublereal part of */

/*       (br(m+1,n+1)*vtbar(m,n,theta(i)) */
/*                   -ci(m+1,n+1)*wtbar(m,n,theta(i)))*cos(m*phi(j)) */
/*      -(bi(m+1,n+1)*vtbar(m,n,theta(i)) */
/*                   +cr(m+1,n+1)*wtbar(m,n,theta(i)))*sin(m*phi(j)) */

/*    and for i=1,...,nlat and  j=1,...,nlon */

/*     wt(i,j) = the sum from n=1 to n=nlat-1 of */

/*              -.5*cr(1,n+1)*vtbar(0,n,theta(i)) */

/*     plus the sum from m=1 to m=mmax-1 of the sum from n=m to */
/*     n=nlat-1 of the doublereal part of */

/*      -(cr(m+1,n+1)*vtbar(m,n,theta(i)) */
/*                   +bi(m+1,n+1)*wtbar(m,n,theta(i)))*cos(m*phi(j)) */
/*      +(ci(m+1,n+1)*vtbar(m,n,theta(i)) */
/*                   -br(m+1,n+1)*wtbar(m,n,theta(i)))*sin(m*phi(j)) */


/*      br(m+1,nlat),bi(m+1,nlat),cr(m+1,nlat), and ci(m+1,nlat) are */
/*      assumed zero for m even. */


/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of ityp */
/*            = 4  error in the specification of nt */
/*            = 5  error in the specification of idvw */
/*            = 6  error in the specification of jdvw */
/*            = 7  error in the specification of mdab */
/*            = 8  error in the specification of ndab */
/*            = 9  error in the specification of lwvts */
/*            = 10 error in the specification of lwork */


/* ******************************************************************* */

/*     subroutine vtseci(nlat,nlon,wvts,lwvts,dwork,ldwork,ierror) */

/*     subroutine vtseci initializes the array wvts which can then be */
/*     used repeatedly by subroutine vtsec until nlat or nlon is changed. */

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

/*     lwvts  the dimension of the array wvts as it appears in the */
/*            program that calls vtsec. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lwvts must be at least */

/*            4*nlat*l2+3*max0(l1-2,0)*(nlat+nlat-l1-1)+nlon+15 */


/*     dwork  a double precision work array that does not have to be saved. */

/*     ldwork the dimension of the array work as it appears in the */
/*            program that calls vtsec. lwork must be at least */
/*            2*(nlat+1) */

/*     ************************************************************** */

/*     output parameters */

/*     wvts   an array which is initialized for use by subroutine vtsec. */
/*            once initialized, wvts can be used repeatedly by vtsec */
/*            as long as nlat or nlon remain unchanged.  wvts must not */
/*            be altered between calls of vtsec. */


/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of lwvts */
/*            = 4  error in the specification of ldwork */

/* ********************************************************************** */

/* Subroutine */ int vtsec_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vt, doublereal *wt, integer *idvw, integer *jdvw, doublereal *
	br, doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvts, integer *lwvts, doublereal *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer vt_dim1, vt_dim2, vt_offset, wt_dim1, wt_dim2, wt_offset, br_dim1,
	     br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, 
	    cr_dim2, cr_offset, ci_dim1, ci_dim2, ci_offset, i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, iw3, iw4, iw5, jw1, jw2, idv, lnl, ist, lzz1, 
	    labc, imid, mmax;
    extern /* Subroutine */ int vtsec1_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static integer lwzvin;


    /* Parameter adjustments */
    wt_dim1 = *idvw;
    wt_dim2 = *jdvw;
    wt_offset = 1 + wt_dim1 * (1 + wt_dim2);
    wt -= wt_offset;
    vt_dim1 = *idvw;
    vt_dim2 = *jdvw;
    vt_offset = 1 + vt_dim1 * (1 + vt_dim2);
    vt -= vt_offset;
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
    --wvts;
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
    if (*lwvts < (lzz1 + labc << 1) + *nlon + 15) {
	return 0;
    }
    *ierror = 10;
/* Computing MAX */
    i__1 = imid * 6;
    if (*ityp <= 2 && *lwork < *nlat * ((*nt << 1) * *nlon + max(i__1,*nlon)))
	     {
	return 0;
    }
/* Computing MAX */
    i__1 = *nlat * 6;
    if (*ityp > 2 && *lwork < imid * ((*nt << 1) * *nlon + max(i__1,*nlon))) {
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
    lzz1 = (*nlat << 1) * imid;
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) * 3 / 2;
    lwzvin = lzz1 + labc;
    jw1 = lwzvin + 1;
    jw2 = jw1 + lwzvin;
    vtsec1_(nlat, nlon, ityp, nt, &imid, idvw, jdvw, &vt[vt_offset], &wt[
	    wt_offset], mdab, ndab, &br[br_offset], &bi[bi_offset], &cr[
	    cr_offset], &ci[ci_offset], &idv, &work[1], &work[iw1], &work[iw2]
	    , &work[iw3], &work[iw4], &work[iw5], &wvts[1], &wvts[jw1], &wvts[
	    jw2]);
    return 0;
} /* vtsec_ */

/* Subroutine */ int vtsec1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *vt, 
	doublereal *wt, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *vte, doublereal *vto, doublereal *wte, doublereal *wto, 
	doublereal *vb, doublereal *wb, doublereal *wvbin, doublereal *wwbin, doublereal *wrfft)
{
    /* System generated locals */
    integer vt_dim1, vt_dim2, vt_offset, wt_dim1, wt_dim2, wt_offset, br_dim1,
	     br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, 
	    cr_dim2, cr_offset, ci_dim1, ci_dim2, ci_offset, vte_dim1, 
	    vte_dim2, vte_offset, vto_dim1, vto_dim2, vto_offset, wte_dim1, 
	    wte_dim2, wte_offset, wto_dim1, wto_dim2, wto_offset, vb_dim1, 
	    vb_dim2, vb_offset, wb_dim1, wb_dim2, wb_offset, i__1, i__2, i__3,
	     i__4;

    /* Local variables */
    static integer i__, j, k, m, iv, iw, mp1, np1, mp2, ndo1, ndo2, imm1, 
	    nlp1, mlat;
    extern /* Subroutine */ int vbin_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *), wbin_(integer *, integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *);
    static integer mmax, mlon, itypp;
    extern /* Subroutine */ int hrfftb_(integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *);

    /* Parameter adjustments */
    wb_dim1 = *imid;
    wb_dim2 = *nlat;
    wb_offset = 1 + wb_dim1 * (1 + wb_dim2);
    wb -= wb_offset;
    vb_dim1 = *imid;
    vb_dim2 = *nlat;
    vb_offset = 1 + vb_dim1 * (1 + vb_dim2);
    vb -= vb_offset;
    wt_dim1 = *idvw;
    wt_dim2 = *jdvw;
    wt_offset = 1 + wt_dim1 * (1 + wt_dim2);
    wt -= wt_offset;
    vt_dim1 = *idvw;
    vt_dim2 = *jdvw;
    vt_offset = 1 + vt_dim1 * (1 + vt_dim2);
    vt -= vt_offset;
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
    wto_dim1 = *idv;
    wto_dim2 = *nlon;
    wto_offset = 1 + wto_dim1 * (1 + wto_dim2);
    wto -= wto_offset;
    wte_dim1 = *idv;
    wte_dim2 = *nlon;
    wte_offset = 1 + wte_dim1 * (1 + wte_dim2);
    wte -= wte_offset;
    vto_dim1 = *idv;
    vto_dim2 = *nlon;
    vto_offset = 1 + vto_dim1 * (1 + vto_dim2);
    vto -= vto_offset;
    vte_dim1 = *idv;
    vte_dim2 = *nlon;
    vte_offset = 1 + vte_dim1 * (1 + vte_dim2);
    vte -= vte_offset;
    --wvbin;
    --wwbin;
    --wrfft;

    /* Function Body */
    nlp1 = *nlat + 1;
    mlat = *nlat % 2;
    mlon = *nlon % 2;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
    imm1 = *imid;
    if (mlat != 0) {
	imm1 = *imid - 1;
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = *idv;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		vte[i__ + (j + k * vte_dim2) * vte_dim1] = 0.f;
		wte[i__ + (j + k * wte_dim2) * wte_dim1] = 0.f;
/* L10: */
	    }
	}
    }
    ndo1 = *nlat;
    ndo2 = *nlat;
    if (mlat != 0) {
	ndo1 = *nlat - 1;
    }
    if (mlat == 0) {
	ndo2 = *nlat - 1;
    }
/* L18: */
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

/*     case ityp=0   no symmetries */

L1:
    vbin_(&c__0, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m = 0 */

    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		vto[i__ + (k * vto_dim2 + 1) * vto_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
		wto[i__ + (k * wto_dim2 + 1) * wto_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L15: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		vte[i__ + (k * vte_dim2 + 1) * vte_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
		wte[i__ + (k * wte_dim2 + 1) * wte_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L16: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	goto L950;
    }
    i__3 = mmax;
    for (mp1 = 2; mp1 <= i__3; ++mp1) {
	m = mp1 - 1;
	mp2 = mp1 + 1;
	vbin_(&c__0, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__0, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L26;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vte[i__ + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    vto[i__ + ((mp1 << 1) - 2 + k * vto_dim2) * vto_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    vte[i__ + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    vto[i__ + ((mp1 << 1) - 1 + k * vto_dim2) * vto_dim1] += 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wte[i__ + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wto[i__ + ((mp1 << 1) - 2 + k * wto_dim2) * wto_dim1] -= 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wte[i__ + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wto[i__ + ((mp1 << 1) - 1 + k * wto_dim2) * wto_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
/* L23: */
		}
		if (mlat == 0) {
		    goto L24;
		}
		vte[*imid + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
		vte[*imid + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
		wte[*imid + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
		wte[*imid + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
L24:
		;
	    }
/* L25: */
	}
L26:
	if (mp2 > ndo2) {
	    goto L30;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo2;
	    for (np1 = mp2; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vto[i__ + ((mp1 << 1) - 2 + k * vto_dim2) * vto_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    vte[i__ + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    vto[i__ + ((mp1 << 1) - 1 + k * vto_dim2) * vto_dim1] += 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    vte[i__ + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wto[i__ + ((mp1 << 1) - 2 + k * wto_dim2) * wto_dim1] -= 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wte[i__ + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wto[i__ + ((mp1 << 1) - 1 + k * wto_dim2) * wto_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wte[i__ + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
/* L27: */
		}
		if (mlat == 0) {
		    goto L28;
		}
		vte[*imid + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
		vte[*imid + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
		wte[*imid + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
		wte[*imid + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
L28:
		;
	    }
/* L29: */
	}
L30:
	;
    }
    goto L950;

/*     case ityp=1   no symmetries,  cr and ci equal zero */

L100:
    vbin_(&c__0, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m = 0 */

    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		vto[i__ + (k * vto_dim2 + 1) * vto_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L115: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		vte[i__ + (k * vte_dim2 + 1) * vte_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L116: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	goto L950;
    }
    i__3 = mmax;
    for (mp1 = 2; mp1 <= i__3; ++mp1) {
	m = mp1 - 1;
	mp2 = mp1 + 1;
	vbin_(&c__0, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__0, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L126;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vte[i__ + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    vte[i__ + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wto[i__ + ((mp1 << 1) - 2 + k * wto_dim2) * wto_dim1] -= 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wto[i__ + ((mp1 << 1) - 1 + k * wto_dim2) * wto_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
/* L123: */
		}
		if (mlat == 0) {
		    goto L124;
		}
		vte[*imid + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
		vte[*imid + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
L124:
		;
	    }
/* L125: */
	}
L126:
	if (mp2 > ndo2) {
	    goto L130;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo2;
	    for (np1 = mp2; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vto[i__ + ((mp1 << 1) - 2 + k * vto_dim2) * vto_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    vto[i__ + ((mp1 << 1) - 1 + k * vto_dim2) * vto_dim1] += 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wte[i__ + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wte[i__ + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
/* L127: */
		}
		if (mlat == 0) {
		    goto L128;
		}
		wte[*imid + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
		wte[*imid + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
L128:
		;
	    }
/* L129: */
	}
L130:
	;
    }
    goto L950;

/*     case ityp=2   no symmetries,  br and bi are equal to zero */

L200:
    vbin_(&c__0, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m = 0 */

    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		wto[i__ + (k * wto_dim2 + 1) * wto_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L215: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		wte[i__ + (k * wte_dim2 + 1) * wte_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L216: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	goto L950;
    }
    i__3 = mmax;
    for (mp1 = 2; mp1 <= i__3; ++mp1) {
	m = mp1 - 1;
	mp2 = mp1 + 1;
	vbin_(&c__0, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__0, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L226;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vto[i__ + ((mp1 << 1) - 2 + k * vto_dim2) * vto_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    vto[i__ + ((mp1 << 1) - 1 + k * vto_dim2) * vto_dim1] += 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wte[i__ + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wte[i__ + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
/* L223: */
		}
		if (mlat == 0) {
		    goto L224;
		}
		wte[*imid + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
		wte[*imid + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
L224:
		;
	    }
/* L225: */
	}
L226:
	if (mp2 > ndo2) {
	    goto L230;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo2;
	    for (np1 = mp2; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vte[i__ + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    vte[i__ + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wto[i__ + ((mp1 << 1) - 2 + k * wto_dim2) * wto_dim1] -= 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wto[i__ + ((mp1 << 1) - 1 + k * wto_dim2) * wto_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
/* L227: */
		}
		if (mlat == 0) {
		    goto L228;
		}
		vte[*imid + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
		vte[*imid + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
L228:
		;
	    }
/* L229: */
	}
L230:
	;
    }
    goto L950;

/*     case ityp=3   v odd,  w even */

L300:
    vbin_(&c__0, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m = 0 */

    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		vto[i__ + (k * vto_dim2 + 1) * vto_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L315: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		wte[i__ + (k * wte_dim2 + 1) * wte_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L316: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	goto L950;
    }
    i__3 = mmax;
    for (mp1 = 2; mp1 <= i__3; ++mp1) {
	m = mp1 - 1;
	mp2 = mp1 + 1;
	vbin_(&c__0, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__0, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L326;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vto[i__ + ((mp1 << 1) - 2 + k * vto_dim2) * vto_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    vto[i__ + ((mp1 << 1) - 1 + k * vto_dim2) * vto_dim1] += 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wte[i__ + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wte[i__ + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
/* L323: */
		}
		if (mlat == 0) {
		    goto L324;
		}
		wte[*imid + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
		wte[*imid + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
L324:
		;
	    }
/* L325: */
	}
L326:
	if (mp2 > ndo2) {
	    goto L330;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo2;
	    for (np1 = mp2; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vto[i__ + ((mp1 << 1) - 2 + k * vto_dim2) * vto_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    vto[i__ + ((mp1 << 1) - 1 + k * vto_dim2) * vto_dim1] += 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wte[i__ + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wte[i__ + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
/* L327: */
		}
		if (mlat == 0) {
		    goto L328;
		}
		wte[*imid + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
		wte[*imid + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
L328:
		;
	    }
/* L329: */
	}
L330:
	;
    }
    goto L950;

/*     case ityp=4   v odd,  w even, and both cr and ci equal zero */

L400:
    vbin_(&c__1, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m = 0 */

    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		vto[i__ + (k * vto_dim2 + 1) * vto_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L415: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	goto L950;
    }
    i__1 = mmax;
    for (mp1 = 2; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	mp2 = mp1 + 1;
	vbin_(&c__1, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__1, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp2 > ndo2) {
	    goto L430;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = ndo2;
	    for (np1 = mp2; np1 <= i__3; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vto[i__ + ((mp1 << 1) - 2 + k * vto_dim2) * vto_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    vto[i__ + ((mp1 << 1) - 1 + k * vto_dim2) * vto_dim1] += 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wte[i__ + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wte[i__ + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
/* L427: */
		}
		if (mlat == 0) {
		    goto L428;
		}
		wte[*imid + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
		wte[*imid + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
L428:
		;
	    }
/* L429: */
	}
L430:
	;
    }
    goto L950;

/*     case ityp=5   v odd,  w even,     br and bi equal zero */

L500:
    vbin_(&c__2, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m = 0 */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		wte[i__ + (k * wte_dim2 + 1) * wte_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L516: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	goto L950;
    }
    i__3 = mmax;
    for (mp1 = 2; mp1 <= i__3; ++mp1) {
	m = mp1 - 1;
	mp2 = mp1 + 1;
	vbin_(&c__2, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__2, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L530;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vto[i__ + ((mp1 << 1) - 2 + k * vto_dim2) * vto_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    vto[i__ + ((mp1 << 1) - 1 + k * vto_dim2) * vto_dim1] += 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wte[i__ + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wte[i__ + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
/* L523: */
		}
		if (mlat == 0) {
		    goto L524;
		}
		wte[*imid + ((mp1 << 1) - 2 + k * wte_dim2) * wte_dim1] -= cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
		wte[*imid + ((mp1 << 1) - 1 + k * wte_dim2) * wte_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
L524:
		;
	    }
/* L525: */
	}
L530:
	;
    }
    goto L950;

/*     case ityp=6   v even  ,  w odd */

L600:
    vbin_(&c__0, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m = 0 */

    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		wto[i__ + (k * wto_dim2 + 1) * wto_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L615: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		vte[i__ + (k * vte_dim2 + 1) * vte_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L616: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	goto L950;
    }
    i__3 = mmax;
    for (mp1 = 2; mp1 <= i__3; ++mp1) {
	m = mp1 - 1;
	mp2 = mp1 + 1;
	vbin_(&c__0, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__0, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L626;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vte[i__ + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    vte[i__ + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wto[i__ + ((mp1 << 1) - 2 + k * wto_dim2) * wto_dim1] -= 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wto[i__ + ((mp1 << 1) - 1 + k * wto_dim2) * wto_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
/* L623: */
		}
		if (mlat == 0) {
		    goto L624;
		}
		vte[*imid + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
		vte[*imid + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
L624:
		;
	    }
/* L625: */
	}
L626:
	if (mp2 > ndo2) {
	    goto L630;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo2;
	    for (np1 = mp2; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vte[i__ + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    vte[i__ + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wto[i__ + ((mp1 << 1) - 2 + k * wto_dim2) * wto_dim1] -= 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wto[i__ + ((mp1 << 1) - 1 + k * wto_dim2) * wto_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
/* L627: */
		}
		if (mlat == 0) {
		    goto L628;
		}
		vte[*imid + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
		vte[*imid + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
L628:
		;
	    }
/* L629: */
	}
L630:
	;
    }
    goto L950;

/*     case ityp=7   v even, w odd   cr and ci equal zero */

L700:
    vbin_(&c__2, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m = 0 */

    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__1 = *imid;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		vte[i__ + (k * vte_dim2 + 1) * vte_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L716: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	goto L950;
    }
    i__1 = mmax;
    for (mp1 = 2; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	mp2 = mp1 + 1;
	vbin_(&c__2, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__2, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp1 > ndo1) {
	    goto L730;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = ndo1;
	    for (np1 = mp1; np1 <= i__3; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vte[i__ + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    vte[i__ + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wto[i__ + ((mp1 << 1) - 2 + k * wto_dim2) * wto_dim1] -= 
			    bi[mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wto[i__ + ((mp1 << 1) - 1 + k * wto_dim2) * wto_dim1] += 
			    br[mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
/* L723: */
		}
		if (mlat == 0) {
		    goto L724;
		}
		vte[*imid + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
		vte[*imid + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[*imid + (
			np1 + iv * vb_dim2) * vb_dim1];
L724:
		;
	    }
/* L725: */
	}
L730:
	;
    }
    goto L950;

/*     case ityp=8   v even,  w odd   br and bi equal zero */

L800:
    vbin_(&c__1, nlat, nlon, &c__0, &vb[vb_offset], &iv, &wvbin[1]);

/*     case m = 0 */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		wto[i__ + (k * wto_dim2 + 1) * wto_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + (np1 + iv * 
			vb_dim2) * vb_dim1];
/* L815: */
	    }
	}
    }

/*     case m = 1 through nlat-1 */

    if (mmax < 2) {
	goto L950;
    }
    i__3 = mmax;
    for (mp1 = 2; mp1 <= i__3; ++mp1) {
	m = mp1 - 1;
	mp2 = mp1 + 1;
	vbin_(&c__1, nlat, nlon, &m, &vb[vb_offset], &iv, &wvbin[1]);
	wbin_(&c__1, nlat, nlon, &m, &wb[wb_offset], &iw, &wwbin[1]);
	if (mp2 > ndo2) {
	    goto L830;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo2;
	    for (np1 = mp2; np1 <= i__1; np1 += 2) {
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vte[i__ + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    vte[i__ + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ 
			    + (np1 + iw * wb_dim2) * wb_dim1];
		    wto[i__ + ((mp1 << 1) - 2 + k * wto_dim2) * wto_dim1] -= 
			    cr[mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
		    wto[i__ + ((mp1 << 1) - 1 + k * wto_dim2) * wto_dim1] -= 
			    ci[mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ 
			    + (np1 + iv * vb_dim2) * vb_dim1];
/* L827: */
		}
		if (mlat == 0) {
		    goto L828;
		}
		vte[*imid + ((mp1 << 1) - 2 + k * vte_dim2) * vte_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
		vte[*imid + ((mp1 << 1) - 1 + k * vte_dim2) * vte_dim1] += cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[*imid + (
			np1 + iw * wb_dim2) * wb_dim1];
L828:
		;
	    }
/* L829: */
	}
L830:
	;
    }
L950:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	hrfftb_(idv, nlon, &vte[(k * vte_dim2 + 1) * vte_dim1 + 1], idv, &
		wrfft[1], &vb[vb_offset]);
	hrfftb_(idv, nlon, &wte[(k * wte_dim2 + 1) * wte_dim1 + 1], idv, &
		wrfft[1], &vb[vb_offset]);
/* L14: */
    }
    if (*ityp > 2) {
	goto L12;
    }
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		vt[i__ + (j + k * vt_dim2) * vt_dim1] = (vte[i__ + (j + k * 
			vte_dim2) * vte_dim1] + vto[i__ + (j + k * vto_dim2) *
			 vto_dim1]) * .5f;
		wt[i__ + (j + k * wt_dim2) * wt_dim1] = (wte[i__ + (j + k * 
			wte_dim2) * wte_dim1] + wto[i__ + (j + k * wto_dim2) *
			 wto_dim1]) * .5f;
		vt[nlp1 - i__ + (j + k * vt_dim2) * vt_dim1] = (vte[i__ + (j 
			+ k * vte_dim2) * vte_dim1] - vto[i__ + (j + k * 
			vto_dim2) * vto_dim1]) * .5f;
		wt[nlp1 - i__ + (j + k * wt_dim2) * wt_dim1] = (wte[i__ + (j 
			+ k * wte_dim2) * wte_dim1] - wto[i__ + (j + k * 
			wto_dim2) * wto_dim1]) * .5f;
/* L60: */
	    }
	}
    }
    goto L13;
L12:
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		vt[i__ + (j + k * vt_dim2) * vt_dim1] = vte[i__ + (j + k * 
			vte_dim2) * vte_dim1] * .5f;
		wt[i__ + (j + k * wt_dim2) * wt_dim1] = wte[i__ + (j + k * 
			wte_dim2) * wte_dim1] * .5f;
/* L11: */
	    }
	}
    }
L13:
    if (mlat == 0) {
	return 0;
    }
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    vt[*imid + (j + k * vt_dim2) * vt_dim1] = vte[*imid + (j + k * 
		    vte_dim2) * vte_dim1] * .5f;
	    wt[*imid + (j + k * wt_dim2) * wt_dim1] = wte[*imid + (j + k * 
		    wte_dim2) * wte_dim1] * .5f;
/* L65: */
	}
    }
    return 0;
} /* vtsec1_ */

/* Subroutine */ int vtseci_(integer *nlat, integer *nlon, doublereal *wvts, 
	integer *lwvts, doublereal *dwork, integer *ldwork, integer *ierror)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer iw1, iw2, lzz1, labc, imid, mmax;
    extern /* Subroutine */ int hrffti_(integer *, doublereal *);
    static integer lwvbin;
    extern /* Subroutine */ int vtinit_(integer *, integer *, doublereal *, 
	    doublereal *), wtinit_(integer *, integer *, doublereal *, doublereal *)
	    ;


    /* Parameter adjustments */
    --wvts;
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
    imid = (*nlat + 1) / 2;
    lzz1 = (*nlat << 1) * imid;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    mmax = min(i__1,i__2);
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) * 3 / 2;
    if (*lwvts < (lzz1 + labc << 1) + *nlon + 15) {
	return 0;
    }
    *ierror = 4;
    if (*ldwork < (*nlat << 1) + 2) {
	return 0;
    }
    *ierror = 0;
    vtinit_(nlat, nlon, &wvts[1], &dwork[1]);
    lwvbin = lzz1 + labc;
    iw1 = lwvbin + 1;
    wtinit_(nlat, nlon, &wvts[iw1], &dwork[1]);
    iw2 = iw1 + lwvbin;
    hrffti_(nlon, &wvts[iw2]);
    return 0;
} /* vtseci_ */

