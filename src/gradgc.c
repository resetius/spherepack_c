/* gradgc.f -- translated by f2c (version 20061008).
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


/* ... file gradgc.f */

/*     this file includes documentation and code for */
/*     subroutine gradgc         i */

/* ... files which must be loaded with gradgc.f */

/*     sphcom.f, hrfft.f, shagc.f,vhsgc.f */

/*     subroutine gradgc(nlat,nlon,isym,nt,v,w,idvw,jdvw,a,b,mdab,ndab, */
/*    +                  wvhsgc,lvhsgc,work,lwork,ierror) */

/*     given the scalar spherical harmonic coefficients a and b, precomputed */
/*     by subroutine shagc for a scalar field sf, subroutine gradgc computes */
/*     an irrotational vector field (v,w) such that */

/*           gradient(sf) = (v,w). */

/*     v is the colatitudinal and w is the east longitudinal component */
/*     of the gradient.  i.e., */

/*            v(i,j) = d(sf(i,j))/dtheta */

/*     and */

/*            w(i,j) = 1/sint*d(sf(i,j))/dlambda */

/*     at the gaussian colatitude point theta(i) (see nlat as input */
/*     parameter) and longitude lambda(j) = (j-1)*2*pi/nlon where */
/*     where sint = sin(theta(i)).  required associated legendre polynomials */
/*     are recomputed rather than stored as they are in subroutine gradgs. this */
/*     saves storage (compare lsav with lsav in gradgs) but increases */
/*     computational requirements. */


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
/*            nlon = 72 for a five degree grid. nlon must be greater than */
/*            3.  the efficiency of the computation is improved when nlon */
/*            is a product of small prime numbers. */


/*     isym   this has the same value as the isym that was input to */
/*            subroutine shagc to compute the arrays a and b from the */
/*            scalar field sf.  isym determines whether (v,w) are */
/*            computed on the full or half sphere as follows: */

/*      = 0 */

/*           sf is not symmetric about the equator. in this case */
/*           the vector field (v,w) is computed on the entire sphere. */
/*           i.e., in the arrays  v(i,j),w(i,j) for i=1,...,nlat and */
/*           j=1,...,nlon. */

/*      = 1 */

/*           sf is antisymmetric about the equator. in this case w is */
/*           antisymmetric and v is symmetric about the equator. w */
/*           and v are computed on the northern hemisphere only.  i.e., */
/*           if nlat is odd they are computed for i=1,...,(nlat+1)/2 */
/*           and j=1,...,nlon.  if nlat is even they are computed for */
/*           i=1,...,nlat/2 and j=1,...,nlon. */

/*      = 2 */

/*           sf is symmetric about the equator. in this case w is */
/*           symmetric and v is antisymmetric about the equator. w */
/*           and v are computed on the northern hemisphere only.  i.e., */
/*           if nlat is odd they are computed for i=1,...,(nlat+1)/2 */
/*           and j=1,...,nlon.  if nlat is even they are computed for */
/*           i=1,...,nlat/2 and j=1,...,nlon. */


/*     nt     nt is the number of scalar and vector fields.  some */
/*            computational efficiency is obtained for multiple fields. */
/*            the arrays a,b,v, and w can be three dimensional corresponding */
/*            to an indexed multiple array sf.  in this case, multiple */
/*            vector synthesis will be performed to compute each vector */
/*            field.  the third index for a,b,v, and w is the synthesis */
/*            index which assumes the values k = 1,...,nt.  for a single */
/*            synthesis set nt = 1.  the description of the remaining */
/*            parameters is simplified by assuming that nt=1 or that a,b,v, */
/*            and w are two dimensional arrays. */

/*     idvw   the first dimension of the arrays v,w as it appears in */
/*            the program that calls gradgc. if isym = 0 then idvw */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then idvw must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then idvw must be at least (nlat+1)/2. */

/*     jdvw   the second dimension of the arrays v,w as it appears in */
/*            the program that calls gradgc. jdvw must be at least nlon. */

/*     a,b    two or three dimensional arrays (see input parameter nt) */
/*            that contain scalar spherical harmonic coefficients */
/*            of the scalar field array sf as computed by subroutine shagc. */
/*     ***    a,b must be computed by shagc prior to calling gradgc. */

/*     mdab   the first dimension of the arrays a and b as it appears in */
/*            the program that calls gradgc (and shagc). mdab must be at */
/*            least min0(nlat,(nlon+2)/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndab   the second dimension of the arrays a and b as it appears in */
/*            the program that calls gradgc (and shagc). ndab must be at */
/*            least nlat. */


/*     wvhsgc an array which must be initialized by subroutine vhsgci. */
/*            once initialized, */
/*            wvhsgc can be used repeatedly by gradgc as long as nlon */
/*            and nlat remain unchanged.  wvhsgc must not be altered */
/*            between calls of gradgc. */


/*     lvhsgc the dimension of the array wvhsgc as it appears in the */
/*            program that calls gradgc. Let */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2                  if nlat is even or */
/*               l2 = (nlat+1)/2              if nlat is odd */

/*            lvhsgc must be at least */

/*               4*nlat*l2+3*max0(l1-2,0)*(nlat+nlat-l1-1)+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls gradgc. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2                  if nlat is even or */
/*               l2 = (nlat+1)/2              if nlat is odd */

/*            if isym = 0 then lwork must be at least */

/*                nlat*(2*nt*nlon+max0(6*l2,nlon)) + nlat*(2*l1*nt+1) */

/*            if isym = 1 or 2 then lwork must be at least */

/*                l2*(2*nt*nlon+max0(6*nlat,nlon)) + nlat*(2*l1*nt+1) */



/*     ************************************************************** */

/*     output parameters */


/*     v,w   two or three dimensional arrays (see input parameter nt) that */
/*           contain an irrotational vector field such that the gradient of */
/*           the scalar field sf is (v,w).  w(i,j) is the east longitude */
/*           component and v(i,j) is the colatitudinal component of velocity */
/*           at gaussian colatitude and longitude lambda(j) = (j-1)*2*pi/nlon */
/*           the indices for v and w are defined at the input parameter */
/*           isym.  the vorticity of (v,w) is zero.  note that any nonzero */
/*           vector field on the sphere will be multiple valued at the poles */
/*           [reference swarztrauber]. */


/*  ierror   = 0  no errors */
/*           = 1  error in the specification of nlat */
/*           = 2  error in the specification of nlon */
/*           = 3  error in the specification of isym */
/*           = 4  error in the specification of nt */
/*           = 5  error in the specification of idvw */
/*           = 6  error in the specification of jdvw */
/*           = 7  error in the specification of mdab */
/*           = 8  error in the specification of ndab */
/*           = 9  error in the specification of lvhsgc */
/*           = 10 error in the specification of lwork */
/* ********************************************************************** */


/* Subroutine */ int gradgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *
	jdvw, doublereal *a, doublereal *b, integer *mdab, integer *ndab, 
	doublereal *wvhsgc, integer *lvhsgc, doublereal *work, integer *lwork,
	 integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, a_dim1, 
	    a_dim2, a_offset, b_dim1, b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    integer l1, l2, mn, is, ibi, ibr, iwk, imid, mmax, liwk, lwmin, lwkmin;
    extern /* Subroutine */ int gradgc1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);


/*     check input parameters */

    /* Parameter adjustments */
    w_dim1 = *idvw;
    w_dim2 = *jdvw;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    v_dim1 = *idvw;
    v_dim2 = *jdvw;
    v_offset = 1 + v_dim1 * (1 + v_dim2);
    v -= v_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    --wvhsgc;
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

/*     verify minimum saved work space length */

/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    l1 = min(i__1,i__2);
    l2 = (*nlat + 1) / 2;
/* Computing MAX */
    i__1 = l1 - 2;
    lwmin = (*nlat << 2) * l2 + max(i__1,0) * 3 * ((*nlat << 1) - l1 - 1) + *
	    nlon + 15;
    if (*lvhsgc < lwmin) {
	return 0;
    }
    *ierror = 10;

/*     verify minimum unsaved work space length */

    if (*isym == 0) {
/* Computing MAX */
	i__1 = l2 * 6;
	lwkmin = *nlat * ((*nt << 1) * *nlon + max(i__1,*nlon) + (l1 << 1) * *
		nt + 1);
    } else {
/* Computing MAX */
	i__1 = *nlat * 6;
	lwkmin = l2 * ((*nt << 1) * *nlon + max(i__1,*nlon)) + *nlat * ((l1 <<
		 1) * *nt + 1);
    }
    if (*lwork < lwkmin) {
	return 0;
    }
    *ierror = 0;

/*     set work space pointers */

    mn = mmax * *nlat * *nt;
    ibr = 1;
    ibi = ibr + mn;
    is = ibi + mn;
    iwk = is + *nlat;
    liwk = *lwork - (mn << 1) - *nlat;
    gradgc1_(nlat, nlon, isym, nt, &v[v_offset], &w[w_offset], idvw, jdvw, &
	    work[ibr], &work[ibi], &mmax, &work[is], mdab, ndab, &a[a_offset],
	     &b[b_offset], &wvhsgc[1], lvhsgc, &work[iwk], &liwk, ierror);
    return 0;
} /* gradgc_ */

/* Subroutine */ int gradgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *
	jdvw, doublereal *br, doublereal *bi, integer *mmax, doublereal *sqnn,
	 integer *mdab, integer *ndab, doublereal *a, doublereal *b, 
	doublereal *wvhsgc, integer *lvhsgc, doublereal *wk, integer *lwk, 
	integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, br_dim1, 
	    br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, a_dim1, a_dim2, 
	    a_offset, b_dim1, b_dim2, b_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer k, m, n;
    doublereal ci, fn, cr;
    integer ityp;
    extern /* Subroutine */ int vhsgc_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);


/*     preset coefficient multiplyers in vector */

    /* Parameter adjustments */
    --sqnn;
    w_dim1 = *idvw;
    w_dim2 = *jdvw;
    w_offset = 1 + w_dim1 * (1 + w_dim2);
    w -= w_offset;
    v_dim1 = *idvw;
    v_dim2 = *jdvw;
    v_offset = 1 + v_dim1 * (1 + v_dim2);
    v -= v_offset;
    bi_dim1 = *mmax;
    bi_dim2 = *nlat;
    bi_offset = 1 + bi_dim1 * (1 + bi_dim2);
    bi -= bi_offset;
    br_dim1 = *mmax;
    br_dim2 = *nlat;
    br_offset = 1 + br_dim1 * (1 + br_dim2);
    br -= br_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    --wvhsgc;
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

/*     preset br,bi to 0.0 */

	i__2 = *nlat;
	for (n = 1; n <= i__2; ++n) {
	    i__3 = *mmax;
	    for (m = 1; m <= i__3; ++m) {
		br[m + (n + k * br_dim2) * br_dim1] = 0.;
		bi[m + (n + k * bi_dim2) * bi_dim1] = 0.;
/* L4: */
	    }
/* L3: */
	}

/*     compute m=0 coefficients */

	i__2 = *nlat;
	for (n = 2; n <= i__2; ++n) {
	    br[(n + k * br_dim2) * br_dim1 + 1] = sqnn[n] * a[(n + k * a_dim2)
		     * a_dim1 + 1];
	    bi[(n + k * bi_dim2) * bi_dim1 + 1] = sqnn[n] * b[(n + k * b_dim2)
		     * b_dim1 + 1];
/* L5: */
	}

/*     compute m>0 coefficients */

	i__2 = *mmax;
	for (m = 2; m <= i__2; ++m) {
	    i__3 = *nlat;
	    for (n = m; n <= i__3; ++n) {
		br[m + (n + k * br_dim2) * br_dim1] = sqnn[n] * a[m + (n + k *
			 a_dim2) * a_dim1];
		bi[m + (n + k * bi_dim2) * bi_dim1] = sqnn[n] * b[m + (n + k *
			 b_dim2) * b_dim1];
/* L7: */
	    }
/* L6: */
	}
/* L2: */
    }

/*     set ityp for irrotational vector synthesis to compute gradient */

    if (*isym == 0) {
	ityp = 1;
    } else if (*isym == 1) {
	ityp = 4;
    } else if (*isym == 2) {
	ityp = 7;
    }

/*     vector sythesize br,bi into (v,w) (cr,ci are dummy variables) */

    vhsgc_(nlat, nlon, &ityp, nt, &v[v_offset], &w[w_offset], idvw, jdvw, &br[
	    br_offset], &bi[bi_offset], &cr, &ci, mmax, nlat, &wvhsgc[1], 
	    lvhsgc, &wk[1], lwk, ierror);
    return 0;
} /* gradgc1_ */

