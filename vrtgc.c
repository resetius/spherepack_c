/* vrtgc.f -- translated by f2c (version 20061008).
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



/* ... file vrtgc.f */

/*     this file includes documentation and code for */
/*     subroutine divec          i */

/* ... files which must be loaded with vrtgc.f */

/*     sphcom.f, hrfft.f, vhagc.f, shsgc.f, gaqd.f */

/*     subroutine vrtgc(nlat,nlon,isym,nt,vort,ivrt,jvrt,cr,ci,mdc,ndc, */
/*    +                 wshsgc,lshsgc,work,lwork,ierror) */

/*     given the vector spherical harmonic coefficients cr and ci, precomputed */
/*     by subroutine vhagc for a vector field (v,w), subroutine vrtgc */
/*     computes the vorticity of the vector field in the scalar array */
/*     vort.  vort(i,j) is the vorticity at the gaussian colatitude */
/*     theta(i) (see nlat as input parameter) and longitude */
/*     lambda(j) = (j-1)*2*pi/nlon on the sphere.  i.e., */

/*            vort(i,j) =  [-dv/dlambda + d(sint*w)/dtheta]/sint */

/*     where sint = sin(theta(i)).  w is the east longitudinal and v */
/*     is the colatitudinal component of the vector field from which */
/*     cr,ci were precomputed.  required associated legendre polynomials */
/*     are recomputed rather than stored as they are in subroutine vrtgs. */


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


/*     isym   a parameter which determines whether the vorticity is */
/*            computed on the full or half sphere as follows: */

/*      = 0 */
/*            the symmetries/antsymmetries described in isym=1,2 below */
/*            do not exist in (v,w) about the equator.  in this case the */
/*            vorticity is neither symmetric nor antisymmetric about */
/*            the equator.  the vorticity is computed on the entire */
/*            sphere.  i.e., in the array vort(i,j) for i=1,...,nlat and */
/*            j=1,...,nlon. */

/*      = 1 */
/*            w is antisymmetric and v is symmetric about the equator. */
/*            in this case the vorticity is symmetyric about the */
/*            equator and is computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the vorticity is computed */
/*            in the array vort(i,j) for i=1,...,(nlat+1)/2 and for */
/*            j=1,...,nlon.  if nlat is even the vorticity is computed */
/*            in the array vort(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*      = 2 */
/*            w is symmetric and v is antisymmetric about the equator */
/*            in this case the vorticity is antisymmetric about the */
/*            equator and is computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the vorticity is computed */
/*            in the array vort(i,j) for i=1,...,(nlat+1)/2 and for */
/*            j=1,...,nlon.  if nlat is even the vorticity is computed */
/*            in the array vort(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */


/*      nt    nt is the number of scalar and vector fields.  some */
/*            computational efficiency is obtained for multiple fields. */
/*            in the program that calls vrtgc, the arrays cr,ci, and vort */
/*            can be three dimensional corresponding to an indexed multiple */
/*            vector field.  in this case multiple scalar synthesis will */
/*            be performed to compute the vorticity for each field.  the */
/*            third index is the synthesis index which assumes the values */
/*            k=1,...,nt.  for a single synthesis set nt = 1.  the */
/*            description of the remaining parameters is simplified by */
/*            assuming that nt=1 or that all the arrays are two dimensional. */

/*     ivrt   the first dimension of the array vort as it appears in */
/*            the program that calls vrtgc. if isym = 0 then ivrt */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then ivrt must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then ivrt must be at least (nlat+1)/2. */

/*     jvrt   the second dimension of the array vort as it appears in */
/*            the program that calls vrtgc. jvrt must be at least nlon. */

/*    cr,ci   two or three dimensional arrays (see input parameter nt) */
/*            that contain vector spherical harmonic coefficients */
/*            of the vector field (v,w) as computed by subroutine vhagc. */
/*     ***    cr and ci must be computed by vhagc prior to calling */
/*            vrtgc. */

/*      mdc   the first dimension of the arrays cr and ci as it */
/*            appears in the program that calls vrtgc. mdc must be at */
/*            least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*      ndc   the second dimension of the arrays cr and ci as it */
/*            appears in the program that calls vrtgc. ndc must be at */
/*            least nlat. */

/*   wshsgc   an array which must be initialized by subroutine shsgci. */
/*            once initialized, */
/*            wshsgc can be used repeatedly by vrtgc as long as nlon */
/*            and nlat remain unchanged.  wshsgc must not be altered */
/*            between calls of vrtgc */

/*   lshsgc   the dimension of the array wshsgc   as it appears in the */
/*            program that calls vrtgc. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshsgc must be at least */

/*               nlat*(2*l2+3*l1-2)+3*l1*(1-l1)/2+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*    lwork   the dimension of the array work as it appears in the */
/*            program that calls vrtgc. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd. */

/*            if isym is zero then lwork must be at least */

/*               nlat*(nlon*nt+max0(3*l2,nlon) + 2*nt*l1+1) */

/*            if isym is not zero then lwork must be at least */

/*               l2*(nlon*nt+max0(3*nlat,nlon)) + nlat*(2*nt*l1+1) */


/*     ************************************************************** */

/*     output parameters */


/*     vort   a two or three dimensional array (see input parameter nt) */
/*            that contains the vorticity of the vector field (v,w) */
/*            whose coefficients cr,ci where computed by subroutine vhagc. */
/*            vort(i,j) is the vorticity at the gaussian colatitude point */
/*            theta(i) and longitude point lambda(j) = (j-1)*2*pi/nlon. */
/*            the index ranges are defined above at the input parameter */
/*            isym. */


/*   ierror   an error parameter which indicates fatal errors with input */
/*            parameters when returned positive. */
/*          = 0  no errors */
/*          = 1  error in the specification of nlat */
/*          = 2  error in the specification of nlon */
/*          = 3  error in the specification of isym */
/*          = 4  error in the specification of nt */
/*          = 5  error in the specification of ivrt */
/*          = 6  error in the specification of jvrt */
/*          = 7  error in the specification of mdc */
/*          = 8  error in the specification of ndc */
/*          = 9  error in the specification of lshsgc */
/*          = 10 error in the specification of lwork */
/* ********************************************************************** */


/* Subroutine */ int vrtgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *vort, integer *ivrt, integer *jvrt, doublereal *cr, doublereal 
	*ci, integer *mdc, integer *ndc, doublereal *wshsgc, integer *lshsgc, doublereal *
	work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer vort_dim1, vort_dim2, vort_offset, cr_dim1, cr_dim2, cr_offset, 
	    ci_dim1, ci_dim2, ci_offset, i__1, i__2;

    /* Local variables */
    static integer l1, l2, ia, ib, mn, is, ls, mab, nln, iwk, lwk, imid, mmax,
	     lpimn, lwmin;
    extern /* Subroutine */ int vrtgc1_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, doublereal *, integer *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *);
    static integer lwkmin;


/*     check input parameters */

    /* Parameter adjustments */
    vort_dim1 = *ivrt;
    vort_dim2 = *jvrt;
    vort_offset = 1 + vort_dim1 * (1 + vort_dim2);
    vort -= vort_offset;
    ci_dim1 = *mdc;
    ci_dim2 = *ndc;
    ci_offset = 1 + ci_dim1 * (1 + ci_dim2);
    ci -= ci_offset;
    cr_dim1 = *mdc;
    cr_dim2 = *ndc;
    cr_offset = 1 + cr_dim1 * (1 + cr_dim2);
    cr -= cr_offset;
    --wshsgc;
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
    if (*isym == 0 && *ivrt < *nlat || *isym > 0 && *ivrt < imid) {
	return 0;
    }
    *ierror = 6;
    if (*jvrt < *nlon) {
	return 0;
    }
    *ierror = 7;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    if (*mdc < min(i__1,i__2)) {
	return 0;
    }
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    mmax = min(i__1,i__2);
    *ierror = 8;
    if (*ndc < *nlat) {
	return 0;
    }
    *ierror = 9;
    imid = (*nlat + 1) / 2;
    lpimn = imid * mmax * (*nlat + *nlat - mmax + 1) / 2;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    l1 = min(i__1,i__2);
    l2 = (*nlat + 1) / 2;
    lwmin = *nlat * ((l2 << 1) + l1 * 3 - 2) + l1 * 3 * (1 - l1) / 2 + *nlon 
	    + 15;
    if (*lshsgc < lwmin) {
	return 0;
    }
    *ierror = 10;

/*     verify unsaved work space (add to what shses requires, file f3) */


/*     set first dimension for a,b (as requried by shses) */

/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mab = min(i__1,i__2);
    mn = mab * *nlat * *nt;
    ls = *nlat;
    if (*isym > 0) {
	ls = imid;
    }
    nln = *nt * ls * *nlon;
/*     if(lwork.lt. nln+ls*nlon+2*mn+nlat) return */
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    l1 = min(i__1,i__2);
    l2 = (*nlat + 1) / 2;
    if (*isym == 0) {
/* Computing MAX */
	i__1 = l2 * 3;
	lwkmin = *nlat * (*nt * *nlon + max(i__1,*nlon) + (*nt << 1) * l1 + 1)
		;
    } else {
/* Computing MAX */
	i__1 = *nlat * 3;
	lwkmin = l2 * (*nt * *nlon + max(i__1,*nlon)) + *nlat * ((*nt << 1) * 
		l1 + 1);
    }
    if (*lwork < lwkmin) {
	return 0;
    }
    *ierror = 0;

/*     set work space pointers */

    ia = 1;
    ib = ia + mn;
    is = ib + mn;
    iwk = is + *nlat;
    lwk = *lwork - (mn << 1) - *nlat;
    vrtgc1_(nlat, nlon, isym, nt, &vort[vort_offset], ivrt, jvrt, &cr[
	    cr_offset], &ci[ci_offset], mdc, ndc, &work[ia], &work[ib], &mab, 
	    &work[is], &wshsgc[1], lshsgc, &work[iwk], &lwk, ierror);
    return 0;
} /* vrtgc_ */

/* Subroutine */ int vrtgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *vort, integer *ivrt, integer *jvrt, doublereal *cr, doublereal 
	*ci, integer *mdc, integer *ndc, doublereal *a, doublereal *b, integer *mab, doublereal 
	*sqnn, doublereal *wsav, integer *lwsav, doublereal *wk, integer *lwk, integer *
	ierror)
{
    /* System generated locals */
    integer vort_dim1, vort_dim2, vort_offset, cr_dim1, cr_dim2, cr_offset, 
	    ci_dim1, ci_dim2, ci_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, m, n;
    static doublereal fn;
    static integer mmax;
    extern /* Subroutine */ int shsgc_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, doublereal *, integer *
	    , integer *, doublereal *, integer *, doublereal *, integer *, integer *);


/*     set coefficient multiplyers */

    /* Parameter adjustments */
    --sqnn;
    vort_dim1 = *ivrt;
    vort_dim2 = *jvrt;
    vort_offset = 1 + vort_dim1 * (1 + vort_dim2);
    vort -= vort_offset;
    ci_dim1 = *mdc;
    ci_dim2 = *ndc;
    ci_offset = 1 + ci_dim1 * (1 + ci_dim2);
    ci -= ci_offset;
    cr_dim1 = *mdc;
    cr_dim2 = *ndc;
    cr_offset = 1 + cr_dim1 * (1 + cr_dim2);
    cr -= cr_offset;
    b_dim1 = *mab;
    b_dim2 = *nlat;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mab;
    a_dim2 = *nlat;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    --wsav;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 2; n <= i__1; ++n) {
	fn = (doublereal) (n - 1);
	sqnn[n] = sqrt(fn * (fn + 1.));
/* L1: */
    }

/*     compute divergence scalar coefficients for each vector field */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlat;
	for (n = 1; n <= i__2; ++n) {
	    i__3 = *mab;
	    for (m = 1; m <= i__3; ++m) {
		a[m + (n + k * a_dim2) * a_dim1] = 0.;
		b[m + (n + k * b_dim2) * b_dim1] = 0.;
/* L4: */
	    }
/* L3: */
	}

/*     compute m=0 coefficients */

	i__2 = *nlat;
	for (n = 2; n <= i__2; ++n) {
	    a[(n + k * a_dim2) * a_dim1 + 1] = sqnn[n] * cr[(n + k * cr_dim2) 
		    * cr_dim1 + 1];
	    b[(n + k * b_dim2) * b_dim1 + 1] = sqnn[n] * ci[(n + k * ci_dim2) 
		    * ci_dim1 + 1];
/* L5: */
	}

/*     compute m>0 coefficients */

/* Computing MIN */
	i__2 = *nlat, i__3 = (*nlon + 1) / 2;
	mmax = min(i__2,i__3);
	i__2 = mmax;
	for (m = 2; m <= i__2; ++m) {
	    i__3 = *nlat;
	    for (n = m; n <= i__3; ++n) {
		a[m + (n + k * a_dim2) * a_dim1] = sqnn[n] * cr[m + (n + k * 
			cr_dim2) * cr_dim1];
		b[m + (n + k * b_dim2) * b_dim1] = sqnn[n] * ci[m + (n + k * 
			ci_dim2) * ci_dim1];
/* L7: */
	    }
/* L6: */
	}
/* L2: */
    }

/*     synthesize a,b into vort */

    shsgc_(nlat, nlon, isym, nt, &vort[vort_offset], ivrt, jvrt, &a[a_offset],
	     &b[b_offset], mab, nlat, &wsav[1], lwsav, &wk[1], lwk, ierror);
    return 0;
} /* vrtgc1_ */

