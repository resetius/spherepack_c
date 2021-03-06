/* vhsgs.f -- translated by f2c (version 20061008).
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



/* ... file vhsgs.f */

/*     this file contains code and documentation for subroutines */
/*     vhsgs and vhsgsi */

/* ... files which must be loaded with vhsgs.f */

/*     sphcom.f, hrfft.f, gaqd.f */

/*     subroutine vhsgs(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci, */
/*    +                 mdab,ndab,wvhsgs,lvhsgs,work,lwork,ierror) */


/*     subroutine vhsgs performs the vector spherical harmonic synthesis */
/*     of the arrays br, bi, cr, and ci and stores the result in the */
/*     arrays v and w.  the synthesis is performed on an equally spaced */
/*     longitude grid and a gaussian colatitude grid (measured from */
/*     the north pole). v(i,j) and w(i,j) are the colatitudinal and */
/*     east longitudinal components respectively, located at the i(th) */
/*     colatitude gaussian point (see nlat below) and longitude */
/*     phi(j) = (j-1)*2*pi/nlon.  the spectral respresentation of (v,w) */
/*     is given below at output parameters v,w. */

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

/*     ityp   = 0  no symmetries exist about the equator. the synthesis */
/*                 is performed on the entire sphere.  i.e. on the */
/*                 arrays v(i,j),w(i,j) for i=1,...,nlat and */
/*                 j=1,...,nlon. */

/*            = 1  no symmetries exist about the equator. the synthesis */
/*                 is performed on the entire sphere.  i.e. on the */
/*                 arrays v(i,j),w(i,j) for i=1,...,nlat and */
/*                 j=1,...,nlon. the curl of (v,w) is zero. that is, */
/*                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. */
/*                 the coefficients cr and ci are zero. */

/*            = 2  no symmetries exist about the equator. the synthesis */
/*                 is performed on the entire sphere.  i.e. on the */
/*                 arrays v(i,j),w(i,j) for i=1,...,nlat and */
/*                 j=1,...,nlon. the divergence of (v,w) is zero. i.e., */
/*                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. */
/*                 the coefficients br and bi are zero. */

/*            = 3  v is symmetric and w is antisymmetric about the */
/*                 equator. the synthesis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the synthesis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the synthesis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 4  v is symmetric and w is antisymmetric about the */
/*                 equator. the synthesis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the synthesis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the synthesis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */
/*                 the curl of (v,w) is zero. that is, */
/*                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. */
/*                 the coefficients cr and ci are zero. */

/*            = 5  v is symmetric and w is antisymmetric about the */
/*                 equator. the synthesis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the synthesis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the synthesis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */
/*                 the divergence of (v,w) is zero. i.e., */
/*                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. */
/*                 the coefficients br and bi are zero. */

/*            = 6  v is antisymmetric and w is symmetric about the */
/*                 equator. the synthesis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the synthesis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the synthesis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*            = 7  v is antisymmetric and w is symmetric about the */
/*                 equator. the synthesis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the synthesis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the synthesis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */
/*                 the curl of (v,w) is zero. that is, */
/*                 (d/dtheta (sin(theta) w) - dv/dphi)/sin(theta) = 0. */
/*                 the coefficients cr and ci are zero. */

/*            = 8  v is antisymmetric and w is symmetric about the */
/*                 equator. the synthesis is performed on the northern */
/*                 hemisphere only.  i.e., if nlat is odd the synthesis */
/*                 is performed on the arrays v(i,j),w(i,j) for */
/*                 i=1,...,(nlat+1)/2 and j=1,...,nlon. if nlat is */
/*                 even the synthesis is performed on the the arrays */
/*                 v(i,j),w(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */
/*                 the divergence of (v,w) is zero. i.e., */
/*                 (d/dtheta (sin(theta) v) + dw/dphi)/sin(theta) = 0. */
/*                 the coefficients br and bi are zero. */


/*     nt     the number of syntheses.  in the program that calls vhsgs, */
/*            the arrays v,w,br,bi,cr, and ci can be three dimensional */
/*            in which case multiple syntheses will be performed. */
/*            the third index is the synthesis index which assumes the */
/*            values k=1,...,nt.  for a single synthesis set nt=1. the */
/*            discription of the remaining parameters is simplified */
/*            by assuming that nt=1 or that all the arrays are two */
/*            dimensional. */

/*     idvw   the first dimension of the arrays v,w as it appears in */
/*            the program that calls vhags. if ityp .le. 2 then idvw */
/*            must be at least nlat.  if ityp .gt. 2 and nlat is */
/*            even then idvw must be at least nlat/2. if ityp .gt. 2 */
/*            and nlat is odd then idvw must be at least (nlat+1)/2. */

/*     jdvw   the second dimension of the arrays v,w as it appears in */
/*            the program that calls vhsgs. jdvw must be at least nlon. */

/*     br,bi  two or three dimensional arrays (see input parameter nt) */
/*     cr,ci  that contain the vector spherical harmonic coefficients */
/*            in the spectral representation of v(i,j) and w(i,j) given */
/*            below at the discription of output parameters v and w. */

/*     mdab   the first dimension of the arrays br,bi,cr, and ci as it */
/*            appears in the program that calls vhsgs. mdab must be at */
/*            least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndab   the second dimension of the arrays br,bi,cr, and ci as it */
/*            appears in the program that calls vhsgs. ndab must be at */
/*            least nlat. */

/*     wvhsgs an array which must be initialized by subroutine vhsgsi. */
/*            once initialized, wvhsgs can be used repeatedly by vhsgs */
/*            as long as nlon and nlat remain unchanged.  wvhsgs must */
/*            not be altered between calls of vhsgs. */

/*     lvhsgs the dimension of the array wvhsgs as it appears in the */
/*            program that calls vhsgs. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lvhsgs must be at least */

/*                 l1*l2*(nlat+nlat-l1+1)+nlon+15+2*nlat */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls vhsgs. define */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            if ityp .le. 2 then lwork must be at least */

/*                       (2*nt+1)*nlat*nlon */

/*            if ityp .gt. 2 then lwork must be at least */

/*                        (2*nt+1)*l2*nlon */

/*     ************************************************************** */

/*     output parameters */

/*     v,w    two or three dimensional arrays (see input parameter nt) */
/*            in which the synthesis is stored. v is the colatitudinal */
/*            component and w is the east longitudinal component. */
/*            v(i,j),w(i,j) contain the components at the guassian colatitude */
/*            point theta(i) and longitude phi(j) = (j-1)*2*pi/nlon. */
/*            the index ranges are defined above at the input parameter */
/*            ityp. v and w are computed from the formulas given below. */


/*     define */

/*     1.  theta is colatitude and phi is east longitude */

/*     2.  the normalized associated legendre funnctions */

/*         pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m) */
/*                        /(2*factorial(n+m)))*sin(theta)**m/(2**n* */
/*                        factorial(n)) times the (n+m)th derivative */
/*                        of (x**2-1)**n with respect to x=cos(theta) */

/*     3.  vbar(m,n,theta) = the derivative of pbar(m,n,theta) with */
/*                           respect to theta divided by the square */
/*                           root of n(n+1). */

/*         vbar(m,n,theta) is more easily computed in the form */

/*         vbar(m,n,theta) = (sqrt((n+m)*(n-m+1))*pbar(m-1,n,theta) */
/*         -sqrt((n-m)*(n+m+1))*pbar(m+1,n,theta))/(2*sqrt(n*(n+1))) */

/*     4.  wbar(m,n,theta) = m/(sin(theta))*pbar(m,n,theta) divided */
/*                           by the square root of n(n+1). */

/*         wbar(m,n,theta) is more easily computed in the form */

/*         wbar(m,n,theta) = sqrt((2n+1)/(2n-1))*(sqrt((n+m)*(n+m-1)) */
/*         *pbar(m-1,n-1,theta)+sqrt((n-m)*(n-m-1))*pbar(m+1,n-1,theta)) */
/*         /(2*sqrt(n*(n+1))) */


/*    the colatitudnal dependence of the normalized surface vector */
/*                spherical harmonics are defined by */

/*     5.    bbar(m,n,theta) = (vbar(m,n,theta),i*wbar(m,n,theta)) */

/*     6.    cbar(m,n,theta) = (i*wbar(m,n,theta),-vbar(m,n,theta)) */


/*    the coordinate to index mappings */

/*     7.   theta(i) = i(th) gaussian grid point and phi(j) = (j-1)*2*pi/nlon */


/*     the maximum (plus one) longitudinal wave number */

/*     8.     mmax = min0(nlat,nlon/2) if nlon is even or */
/*            mmax = min0(nlat,(nlon+1)/2) if nlon is odd. */

/*    if we further define the output vector as */

/*     9.    h(i,j) = (v(i,j),w(i,j)) */

/*    and the complex coefficients */

/*     10.   b(m,n) = cmplx(br(m+1,n+1),bi(m+1,n+1)) */

/*     11.   c(m,n) = cmplx(cr(m+1,n+1),ci(m+1,n+1)) */


/*    then for i=1,...,nlat and  j=1,...,nlon */

/*        the expansion for real h(i,j) takes the form */

/*     h(i,j) = the sum from n=1 to n=nlat-1 of the real part of */

/*         .5*(b(0,n)*bbar(0,n,theta(i))+c(0,n)*cbar(0,n,theta(i))) */

/*     plus the sum from m=1 to m=mmax-1 of the sum from n=m to */
/*     n=nlat-1 of the real part of */

/*              b(m,n)*bbar(m,n,theta(i))*exp(i*m*phi(j)) */
/*             +c(m,n)*cbar(m,n,theta(i))*exp(i*m*phi(j)) */

/*   ************************************************************* */

/*   in terms of real variables this expansion takes the form */

/*             for i=1,...,nlat and  j=1,...,nlon */

/*     v(i,j) = the sum from n=1 to n=nlat-1 of */

/*               .5*br(1,n+1)*vbar(0,n,theta(i)) */

/*     plus the sum from m=1 to m=mmax-1 of the sum from n=m to */
/*     n=nlat-1 of the real part of */

/*       (br(m+1,n+1)*vbar(m,n,theta(i))-ci(m+1,n+1)*wbar(m,n,theta(i))) */
/*                                          *cos(m*phi(j)) */
/*      -(bi(m+1,n+1)*vbar(m,n,theta(i))+cr(m+1,n+1)*wbar(m,n,theta(i))) */
/*                                          *sin(m*phi(j)) */

/*    and for i=1,...,nlat and  j=1,...,nlon */

/*     w(i,j) = the sum from n=1 to n=nlat-1 of */

/*              -.5*cr(1,n+1)*vbar(0,n,theta(i)) */

/*     plus the sum from m=1 to m=mmax-1 of the sum from n=m to */
/*     n=nlat-1 of the real part of */

/*      -(cr(m+1,n+1)*vbar(m,n,theta(i))+bi(m+1,n+1)*wbar(m,n,theta(i))) */
/*                                          *cos(m*phi(j)) */
/*      +(ci(m+1,n+1)*vbar(m,n,theta(i))-br(m+1,n+1)*wbar(m,n,theta(i))) */
/*                                          *sin(m*phi(j)) */


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
/*            = 9  error in the specification of lvhsgs */
/*            = 10 error in the specification of lwork */


/*     subroutine vhsgsi(nlat,nlon,wvhsgs,lvhsgs,dwork,ldwork,ierror) */

/*     subroutine vhsgsi initializes the array wvhsgs which can then be */
/*     used repeatedly by subroutine vhsgs until nlat or nlon is changed. */

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

/*     lvhsgs the dimension of the array wvhsgs as it appears in the */
/*            program that calls vhsgs. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lvhsgs must be at least */

/*                 l1*l2*(nlat+nlat-l1+1)+nlon+15+2*nlat */

/*     dwork a double precision work array that does not need to be saved */

/*     ldwork the dimension of the array dwork as it appears in the */
/*            program that calls vhsgsi. ldwork must be at least */

/*                 (3*nlat*(nlat+3)+2)/2 */

/*     ************************************************************** */

/*     output parameters */

/*     wvhsgs an array which is initialized for use by subroutine vhsgs. */
/*            once initialized, wvhsgs can be used repeatedly by vhsgs */
/*            as long as nlat and nlon remain unchanged.  wvhsgs must not */
/*            be altered between calls of vhsgs. */


/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of lvhsgs */
/*            = 4  error in the specification of lwork */

/* Subroutine */ int vhsgs_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *
	jdvw, doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, 
	integer *mdab, integer *ndab, doublereal *wvhsgs, integer *lvhsgs, 
	doublereal *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, br_dim1, 
	    br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, cr_dim2,
	     cr_offset, ci_dim1, ci_dim2, ci_offset, i__1, i__2;

    /* Local variables */
    integer iw1, jw1, jw2, jw3, iw2, iw3, iw4, idv, lnl, idz, lmn, ist, imid, 
	    mmax, lzimn;
    extern /* Subroutine */ int vhsgs1_(integer *, integer *, integer *, 
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
    --wvhsgs;
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
    if (*lvhsgs < lzimn + lzimn + *nlon + 15) {
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

/*     set wvhsgs pointers */

    lmn = *nlat * (*nlat + 1) / 2;
    jw1 = 1;
    jw2 = jw1 + imid * lmn;
    jw3 = jw2 + imid * lmn;

/*     set work pointers */

    iw1 = ist + 1;
    iw2 = lnl + 1;
    iw3 = iw2 + ist;
    iw4 = iw2 + lnl;
    vhsgs1_(nlat, nlon, ityp, nt, &imid, idvw, jdvw, &v[v_offset], &w[
	    w_offset], mdab, ndab, &br[br_offset], &bi[bi_offset], &cr[
	    cr_offset], &ci[ci_offset], &idv, &work[1], &work[iw1], &work[iw2]
	    , &work[iw3], &work[iw4], &idz, &wvhsgs[jw1], &wvhsgs[jw2], &
	    wvhsgs[jw3]);
    return 0;
} /* vhsgs_ */

/* Subroutine */ int vhsgs1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *
	v, doublereal *w, integer *mdab, integer *ndab, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci, integer *idv, 
	doublereal *ve, doublereal *vo, doublereal *we, doublereal *wo, 
	doublereal *work, integer *idz, doublereal *vb, doublereal *wb, 
	doublereal *wrfft)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, br_dim1, 
	    br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, cr_dim1, cr_dim2,
	     cr_offset, ci_dim1, ci_dim2, ci_offset, ve_dim1, ve_dim2, 
	    ve_offset, vo_dim1, vo_dim2, vo_offset, we_dim1, we_dim2, 
	    we_offset, wo_dim1, wo_dim2, wo_offset, vb_dim1, vb_offset, 
	    wb_dim1, wb_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    integer i__, j, k, m, mb, mn, mp1, np1, mp2, ndo1, ndo2, imm1, nlp1, mlat,
	     mmax, mlon, itypp;
    extern /* Subroutine */ int hrfftb_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *);

    /* Parameter adjustments */
    wb_dim1 = *imid;
    wb_offset = 1 + wb_dim1;
    wb -= wb_offset;
    vb_dim1 = *imid;
    vb_offset = 1 + vb_dim1;
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
    --work;
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
		ve[i__ + (j + k * ve_dim2) * ve_dim1] = 0.;
		we[i__ + (j + k * we_dim2) * we_dim1] = 0.;
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

/*     case m = 0 */

L1:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = *imid;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ve[i__ + (k * ve_dim2 + 1) * ve_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + np1 * vb_dim1];
		we[i__ + (k * we_dim2 + 1) * we_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + np1 * vb_dim1];
/* L15: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		vo[i__ + (k * vo_dim2 + 1) * vo_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + np1 * vb_dim1];
		wo[i__ + (k * wo_dim2 + 1) * wo_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + np1 * vb_dim1];
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
/*     mb = m*(nlat-1)-(m*(m-1))/2 */
	mb = m * *nlat - m * (m + 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L26;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] += bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] -= cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ + 
			    mn * wb_dim1];
/* L23: */
		}
		if (mlat == 0) {
		    goto L24;
		}
		ve[*imid + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[*imid + mn *
			 wb_dim1];
		ve[*imid + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[*imid + mn *
			 wb_dim1];
		we[*imid + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[*imid + mn *
			 wb_dim1];
		we[*imid + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * wb[*imid + mn *
			 wb_dim1];
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
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] += cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] -= bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ + 
			    mn * wb_dim1];
/* L27: */
		}
		if (mlat == 0) {
		    goto L28;
		}
		ve[*imid + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * vb[*imid + mn *
			 vb_dim1];
		ve[*imid + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[*imid + mn *
			 vb_dim1];
		we[*imid + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[*imid + mn *
			 vb_dim1];
		we[*imid + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[*imid + mn *
			 vb_dim1];
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

/*     case m = 0 */

L100:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = *imid;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ve[i__ + (k * ve_dim2 + 1) * ve_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + np1 * vb_dim1];
/* L115: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		vo[i__ + (k * vo_dim2 + 1) * vo_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + np1 * vb_dim1];
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
/*     mb = m*(nlat-1)-(m*(m-1))/2 */
	mb = m * *nlat - m * (m + 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L126;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] += bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ + 
			    mn * wb_dim1];
/* L123: */
		}
		if (mlat == 0) {
		    goto L124;
		}
		we[*imid + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[*imid + mn *
			 wb_dim1];
		we[*imid + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * wb[*imid + mn *
			 wb_dim1];
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
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] -= bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ + 
			    mn * wb_dim1];
/* L127: */
		}
		if (mlat == 0) {
		    goto L128;
		}
		ve[*imid + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * vb[*imid + mn *
			 vb_dim1];
		ve[*imid + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[*imid + mn *
			 vb_dim1];
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

/*     case m = 0 */

L200:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = *imid;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		we[i__ + (k * we_dim2 + 1) * we_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + np1 * vb_dim1];
/* L215: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		wo[i__ + (k * wo_dim2 + 1) * wo_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + np1 * vb_dim1];
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
/*     mb = m*(nlat-1)-(m*(m-1))/2 */
	mb = m * *nlat - m * (m + 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L226;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] -= cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ + 
			    mn * vb_dim1];
/* L223: */
		}
		if (mlat == 0) {
		    goto L224;
		}
		ve[*imid + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[*imid + mn *
			 wb_dim1];
		ve[*imid + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[*imid + mn *
			 wb_dim1];
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
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] += cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ + 
			    mn * vb_dim1];
/* L227: */
		}
		if (mlat == 0) {
		    goto L228;
		}
		we[*imid + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[*imid + mn *
			 vb_dim1];
		we[*imid + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[*imid + mn *
			 vb_dim1];
L228:
		;
	    }
/* L229: */
	}
L230:
	;
    }
    goto L950;

/*     case ityp=3   v even,  w odd */

/*     case m = 0 */

L300:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = *imid;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ve[i__ + (k * ve_dim2 + 1) * ve_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + np1 * vb_dim1];
/* L315: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		wo[i__ + (k * wo_dim2 + 1) * wo_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + np1 * vb_dim1];
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
/*     mb = m*(nlat-1)-(m*(m-1))/2 */
	mb = m * *nlat - m * (m + 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L326;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] -= cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ + 
			    mn * vb_dim1];
/* L323: */
		}
		if (mlat == 0) {
		    goto L324;
		}
		ve[*imid + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[*imid + mn *
			 wb_dim1];
		ve[*imid + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[*imid + mn *
			 wb_dim1];
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
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] -= bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ + 
			    mn * wb_dim1];
/* L327: */
		}
		if (mlat == 0) {
		    goto L328;
		}
		ve[*imid + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * vb[*imid + mn *
			 vb_dim1];
		ve[*imid + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[*imid + mn *
			 vb_dim1];
L328:
		;
	    }
/* L329: */
	}
L330:
	;
    }
    goto L950;

/*     case ityp=4   v even,  w odd, and both cr and ci equal zero */

/*     case m = 0 */

L400:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = *imid;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ve[i__ + (k * ve_dim2 + 1) * ve_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + np1 * vb_dim1];
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
/*     mb = m*(nlat-1)-(m*(m-1))/2 */
	mb = m * *nlat - m * (m + 1) / 2;
	mp2 = mp1 + 1;
	if (mp2 > ndo2) {
	    goto L430;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = ndo2;
	    for (np1 = mp2; np1 <= i__3; np1 += 2) {
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] -= bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ + 
			    mn * wb_dim1];
/* L427: */
		}
		if (mlat == 0) {
		    goto L428;
		}
		ve[*imid + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * vb[*imid + mn *
			 vb_dim1];
		ve[*imid + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[*imid + mn *
			 vb_dim1];
L428:
		;
	    }
/* L429: */
	}
L430:
	;
    }
    goto L950;

/*     case ityp=5   v even,  w odd,     br and bi equal zero */

/*     case m = 0 */

L500:
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		wo[i__ + (k * wo_dim2 + 1) * wo_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + np1 * vb_dim1];
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
/*     mb = m*(nlat-1)-(m*(m-1))/2 */
	mb = m * *nlat - m * (m + 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L530;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    ve[i__ + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    ve[i__ + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    wo[i__ + ((mp1 << 1) - 2 + k * wo_dim2) * wo_dim1] -= cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    wo[i__ + ((mp1 << 1) - 1 + k * wo_dim2) * wo_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ + 
			    mn * vb_dim1];
/* L523: */
		}
		if (mlat == 0) {
		    goto L524;
		}
		ve[*imid + ((mp1 << 1) - 2 + k * ve_dim2) * ve_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[*imid + mn *
			 wb_dim1];
		ve[*imid + ((mp1 << 1) - 1 + k * ve_dim2) * ve_dim1] += cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[*imid + mn *
			 wb_dim1];
L524:
		;
	    }
/* L525: */
	}
L530:
	;
    }
    goto L950;

/*     case ityp=6   v odd  ,  w even */

/*     case m = 0 */

L600:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__1 = *imid;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		we[i__ + (k * we_dim2 + 1) * we_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + np1 * vb_dim1];
/* L615: */
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__3 = imm1;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		vo[i__ + (k * vo_dim2 + 1) * vo_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + np1 * vb_dim1];
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
/*     mb = m*(nlat-1)-(m*(m-1))/2 */
	mb = m * *nlat - m * (m + 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L626;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo1;
	    for (np1 = mp1; np1 <= i__1; np1 += 2) {
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] += bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ + 
			    mn * wb_dim1];
/* L623: */
		}
		if (mlat == 0) {
		    goto L624;
		}
		we[*imid + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[*imid + mn *
			 wb_dim1];
		we[*imid + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * wb[*imid + mn *
			 wb_dim1];
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
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] += cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ + 
			    mn * vb_dim1];
/* L627: */
		}
		if (mlat == 0) {
		    goto L628;
		}
		we[*imid + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[*imid + mn *
			 vb_dim1];
		we[*imid + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[*imid + mn *
			 vb_dim1];
L628:
		;
	    }
/* L629: */
	}
L630:
	;
    }
    goto L950;

/*     case ityp=7   v odd, w even   cr and ci equal zero */

/*     case m = 0 */

L700:
    i__3 = *nt;
    for (k = 1; k <= i__3; ++k) {
	i__2 = ndo1;
	for (np1 = 3; np1 <= i__2; np1 += 2) {
	    i__1 = imm1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		vo[i__ + (k * vo_dim2 + 1) * vo_dim1] += br[(np1 + k * 
			br_dim2) * br_dim1 + 1] * vb[i__ + np1 * vb_dim1];
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
/*     mb = m*(nlat-1)-(m*(m-1))/2 */
	mb = m * *nlat - m * (m + 1) / 2;
	mp2 = mp1 + 1;
	if (mp1 > ndo1) {
	    goto L730;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = ndo1;
	    for (np1 = mp1; np1 <= i__3; np1 += 2) {
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] += bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= bi[
			    mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] += br[
			    mp1 + (np1 + k * br_dim2) * br_dim1] * wb[i__ + 
			    mn * wb_dim1];
/* L723: */
		}
		if (mlat == 0) {
		    goto L724;
		}
		we[*imid + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= bi[
			mp1 + (np1 + k * bi_dim2) * bi_dim1] * wb[*imid + mn *
			 wb_dim1];
		we[*imid + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] += br[
			mp1 + (np1 + k * br_dim2) * br_dim1] * wb[*imid + mn *
			 wb_dim1];
L724:
		;
	    }
/* L725: */
	}
L730:
	;
    }
    goto L950;

/*     case ityp=8   v odd,  w even   br and bi equal zero */

/*     case m = 0 */

L800:
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = ndo2;
	for (np1 = 2; np1 <= i__2; np1 += 2) {
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		we[i__ + (k * we_dim2 + 1) * we_dim1] -= cr[(np1 + k * 
			cr_dim2) * cr_dim1 + 1] * vb[i__ + np1 * vb_dim1];
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
/*     mb = m*(nlat-1)-(m*(m-1))/2 */
	mb = m * *nlat - m * (m + 1) / 2;
	mp2 = mp1 + 1;
	if (mp2 > ndo2) {
	    goto L830;
	}
	i__2 = *nt;
	for (k = 1; k <= i__2; ++k) {
	    i__1 = ndo2;
	    for (np1 = mp2; np1 <= i__1; np1 += 2) {
		mn = mb + np1;
		i__4 = imm1;
		for (i__ = 1; i__ <= i__4; ++i__) {
		    vo[i__ + ((mp1 << 1) - 2 + k * vo_dim2) * vo_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    vo[i__ + ((mp1 << 1) - 1 + k * vo_dim2) * vo_dim1] += cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * wb[i__ + 
			    mn * wb_dim1];
		    we[i__ + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= cr[
			    mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[i__ + 
			    mn * vb_dim1];
		    we[i__ + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] -= ci[
			    mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[i__ + 
			    mn * vb_dim1];
/* L827: */
		}
		if (mlat == 0) {
		    goto L828;
		}
		we[*imid + ((mp1 << 1) - 2 + k * we_dim2) * we_dim1] -= cr[
			mp1 + (np1 + k * cr_dim2) * cr_dim1] * vb[*imid + mn *
			 vb_dim1];
		we[*imid + ((mp1 << 1) - 1 + k * we_dim2) * we_dim1] -= ci[
			mp1 + (np1 + k * ci_dim2) * ci_dim1] * vb[*imid + mn *
			 vb_dim1];
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
	hrfftb_(idv, nlon, &ve[(k * ve_dim2 + 1) * ve_dim1 + 1], idv, &wrfft[
		1], &work[1]);
	hrfftb_(idv, nlon, &we[(k * we_dim2 + 1) * we_dim1 + 1], idv, &wrfft[
		1], &work[1]);
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
		v[i__ + (j + k * v_dim2) * v_dim1] = (ve[i__ + (j + k * 
			ve_dim2) * ve_dim1] + vo[i__ + (j + k * vo_dim2) * 
			vo_dim1]) * .5;
		w[i__ + (j + k * w_dim2) * w_dim1] = (we[i__ + (j + k * 
			we_dim2) * we_dim1] + wo[i__ + (j + k * wo_dim2) * 
			wo_dim1]) * .5;
		v[nlp1 - i__ + (j + k * v_dim2) * v_dim1] = (ve[i__ + (j + k *
			 ve_dim2) * ve_dim1] - vo[i__ + (j + k * vo_dim2) * 
			vo_dim1]) * .5;
		w[nlp1 - i__ + (j + k * w_dim2) * w_dim1] = (we[i__ + (j + k *
			 we_dim2) * we_dim1] - wo[i__ + (j + k * wo_dim2) * 
			wo_dim1]) * .5;
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
		v[i__ + (j + k * v_dim2) * v_dim1] = ve[i__ + (j + k * 
			ve_dim2) * ve_dim1] * .5;
		w[i__ + (j + k * w_dim2) * w_dim1] = we[i__ + (j + k * 
			we_dim2) * we_dim1] * .5;
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
	    v[*imid + (j + k * v_dim2) * v_dim1] = ve[*imid + (j + k * 
		    ve_dim2) * ve_dim1] * .5;
	    w[*imid + (j + k * w_dim2) * w_dim1] = we[*imid + (j + k * 
		    we_dim2) * we_dim1] * .5;
/* L65: */
	}
    }
    return 0;
} /* vhsgs1_ */

/* Subroutine */ int vhsgsi_(integer *nlat, integer *nlon, doublereal *wvhsgs,
	 integer *lvhsgs, doublereal *dwork, integer *ldwork, integer *ierror)
{
    integer iw1, jw1, jw2, jw3, iw2, iw3, iw4, lmn, imid;
    extern /* Subroutine */ int vhgsi1_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), hrffti_(integer *, doublereal *);


/*     subroutine vhsfsi computes the gaussian points theta, gauss */
/*     weights wts, and the components vb and wb of the vector */
/*     harmonics. all quantities are computed internally in double */
/*     precision but returned in single precision and are therfore */
/*     accurate to single precision. */

/*     set imid = (nlat+1)/2 and lmn=(nlat*(nlat+1))/2 then */
/*     wvhsgs must have 2*(imid*lmn+nlat)+nlon+15 locations */

/*     double precision array dwork must have */
/*       3*nlat*(nlat+1)+5*nlat+1 = nlat*(3*nlat+8)+1 */
/*     locations which is determined by the size of dthet, */
/*     dwts, dwork, and dpbar in vhsgs1 */

    /* Parameter adjustments */
    --dwork;
    --wvhsgs;

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
    lmn = *nlat * (*nlat + 1) / 2;
    if (*lvhsgs < (imid * lmn << 1) + *nlon + 15) {
	return 0;
    }
    *ierror = 4;
    if (*ldwork < (*nlat * 3 * (*nlat + 3) + 2) / 2) {
	return 0;
    }
    *ierror = 0;

/*     set saved work space pointers */

    jw1 = 1;
    jw2 = jw1 + imid * lmn;
    jw3 = jw2 + imid * lmn;

/*     set unsaved work space pointers */

    iw1 = 1;
    iw2 = iw1 + *nlat;
    iw3 = iw2 + *nlat;
    iw4 = iw3 + imid * 3 * *nlat;
/*     iw2 = iw1+nlat+nlat */
/*     iw3 = iw2+nlat+nlat */
/*     iw4 = iw3+6*imid*nlat */
    vhgsi1_(nlat, &imid, &wvhsgs[jw1], &wvhsgs[jw2], &dwork[iw1], &dwork[iw2],
	     &dwork[iw3], &dwork[iw4]);
    hrffti_(nlon, &wvhsgs[jw3]);
    return 0;
} /* vhsgsi_ */

/* Subroutine */ int vhgsi1_(integer *nlat, integer *imid, doublereal *vb, 
	doublereal *wb, doublereal *dthet, doublereal *dwts, doublereal *
	dpbar, doublereal *work)
{
    /* System generated locals */
    integer vb_dim1, vb_offset, wb_dim1, wb_offset, dpbar_dim1, dpbar_dim2, 
	    dpbar_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, m, n, id, nm, mn, np, ix, iy, nz;
    doublereal dcf;
    integer lwk;
    doublereal abel, bbel, cbel;
    extern /* Subroutine */ int gaqd_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);
    extern integer indx_(integer *, integer *, integer *);
    doublereal ssqr2;
    extern /* Subroutine */ int dnlfk_(integer *, integer *, doublereal *), 
	    dnlft_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *);
    integer ierror;


/*     compute gauss points and weights */
/*     use dpbar (length 3*nnlat*(nnlat+1)) as work space for gaqd */

    /* Parameter adjustments */
    dpbar_dim1 = *imid;
    dpbar_dim2 = *nlat;
    dpbar_offset = 1 + dpbar_dim1 * (1 + dpbar_dim2);
    dpbar -= dpbar_offset;
    wb_dim1 = *imid;
    wb_offset = 1 + wb_dim1;
    wb -= wb_offset;
    vb_dim1 = *imid;
    vb_offset = 1 + vb_dim1;
    vb -= vb_offset;
    --dthet;
    --dwts;
    --work;

    /* Function Body */
    lwk = *nlat * (*nlat + 2);
    gaqd_(nlat, &dthet[1], &dwts[1], &dpbar[dpbar_offset], &lwk, &ierror);

/*     compute associated legendre functions */

/*     compute m=n=0 legendre polynomials for all theta(i) */

    ssqr2 = 1. / sqrt(2.);
    i__1 = *imid;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dpbar[i__ + (dpbar_dim2 + 1) * dpbar_dim1] = ssqr2;
	vb[i__ + vb_dim1] = 0.;
	wb[i__ + wb_dim1] = 0.;
/* L90: */
    }

/*     main loop for remaining vb, and wb */

    i__1 = *nlat - 1;
    for (n = 1; n <= i__1; ++n) {
	nm = (n - 2) % 3 + 1;
	nz = (n - 1) % 3 + 1;
	np = n % 3 + 1;

/*     compute dpbar for m=0 */

	dnlfk_(&c__0, &n, &work[1]);
	mn = indx_(&c__0, &n, nlat);
	i__2 = *imid;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dnlft_(&c__0, &n, &dthet[i__], &work[1], &dpbar[i__ + (np * 
		    dpbar_dim2 + 1) * dpbar_dim1]);
/*      pbar(i,mn) = dpbar(i,1,np) */
/* L105: */
	}

/*     compute dpbar for m=1 */

	dnlfk_(&c__1, &n, &work[1]);
	mn = indx_(&c__1, &n, nlat);
	i__2 = *imid;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    dnlft_(&c__1, &n, &dthet[i__], &work[1], &dpbar[i__ + (np * 
		    dpbar_dim2 + 2) * dpbar_dim1]);
/*      pbar(i,mn) = dpbar(i,2,np) */
/* L106: */
	}
/* L104: */

/*     compute and store dpbar for m=2,n */

	if (n < 2) {
	    goto L108;
	}
	i__2 = n;
	for (m = 2; m <= i__2; ++m) {
	    abel = sqrt((doublereal) (((n << 1) + 1) * (m + n - 2) * (m + n - 
		    3)) / (doublereal) (((n << 1) - 3) * (m + n - 1) * (m + n)
		    ));
	    bbel = sqrt((doublereal) (((n << 1) + 1) * (n - m - 1) * (n - m)) 
		    / (doublereal) (((n << 1) - 3) * (m + n - 1) * (m + n)));
	    cbel = sqrt((doublereal) ((n - m + 1) * (n - m + 2)) / (
		    doublereal) ((m + n - 1) * (m + n)));
	    id = indx_(&m, &n, nlat);
	    if (m >= n - 1) {
		goto L102;
	    }
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		dpbar[i__ + (m + 1 + np * dpbar_dim2) * dpbar_dim1] = abel * 
			dpbar[i__ + (m - 1 + nm * dpbar_dim2) * dpbar_dim1] + 
			bbel * dpbar[i__ + (m + 1 + nm * dpbar_dim2) * 
			dpbar_dim1] - cbel * dpbar[i__ + (m - 1 + np * 
			dpbar_dim2) * dpbar_dim1];
/*      pbar(i,id) = dpbar(i,m+1,np) */
/* L103: */
	    }
	    goto L107;
L102:
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		dpbar[i__ + (m + 1 + np * dpbar_dim2) * dpbar_dim1] = abel * 
			dpbar[i__ + (m - 1 + nm * dpbar_dim2) * dpbar_dim1] - 
			cbel * dpbar[i__ + (m - 1 + np * dpbar_dim2) * 
			dpbar_dim1];
/*      pbar(i,id) = dpbar(i,m+1,np) */
/* L101: */
	    }
L107:
	    ;
	}

/*     compute the derivative of the functions */

L108:
	ix = indx_(&c__0, &n, nlat);
	iy = indx_(&n, &n, nlat);
	i__2 = *imid;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    vb[i__ + ix * vb_dim1] = -dpbar[i__ + (np * dpbar_dim2 + 2) * 
		    dpbar_dim1];
	    vb[i__ + iy * vb_dim1] = dpbar[i__ + (n + np * dpbar_dim2) * 
		    dpbar_dim1] / sqrt((doublereal) (n + 1 << 1));
/* L125: */
	}

	if (n == 1) {
	    goto L131;
	}
	dcf = sqrt((doublereal) ((n << 2) * (n + 1)));
	i__2 = n - 1;
	for (m = 1; m <= i__2; ++m) {
	    ix = indx_(&m, &n, nlat);
	    abel = sqrt((doublereal) ((n + m) * (n - m + 1))) / dcf;
	    bbel = sqrt((doublereal) ((n - m) * (n + m + 1))) / dcf;
	    i__3 = *imid;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		vb[i__ + ix * vb_dim1] = abel * dpbar[i__ + (m + np * 
			dpbar_dim2) * dpbar_dim1] - bbel * dpbar[i__ + (m + 2 
			+ np * dpbar_dim2) * dpbar_dim1];
/* L130: */
	    }
	}

/*     compute the vector harmonic w(theta) = m*pbar/cos(theta) */

/*     set wb=0 for m=0 */

L131:
	ix = indx_(&c__0, &n, nlat);
	i__3 = *imid;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    wb[i__ + ix * wb_dim1] = 0.;
/* L220: */
	}

/*     compute wb for m=1,n */

	dcf = sqrt((doublereal) (n + n + 1) / (doublereal) ((n << 2) * (n + 1)
		 * (n + n - 1)));
	i__3 = n;
	for (m = 1; m <= i__3; ++m) {
	    ix = indx_(&m, &n, nlat);
	    abel = dcf * sqrt((doublereal) ((n + m) * (n + m - 1)));
	    bbel = dcf * sqrt((doublereal) ((n - m) * (n - m - 1)));
	    if (m >= n - 1) {
		goto L231;
	    }
	    i__2 = *imid;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		wb[i__ + ix * wb_dim1] = abel * dpbar[i__ + (m + nz * 
			dpbar_dim2) * dpbar_dim1] + bbel * dpbar[i__ + (m + 2 
			+ nz * dpbar_dim2) * dpbar_dim1];
/* L229: */
	    }
	    goto L230;
L231:
	    i__2 = *imid;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		wb[i__ + ix * wb_dim1] = abel * dpbar[i__ + (m + nz * 
			dpbar_dim2) * dpbar_dim1];
/* L228: */
	    }
L230:
	    ;
	}
/* L100: */
    }
    return 0;
} /* vhgsi1_ */

