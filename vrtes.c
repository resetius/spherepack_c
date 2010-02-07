/* vrtes.f -- translated by f2c (version 20061008).
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



/* ... file vrtes.f */

/*     this file includes documentation and code for */
/*     subroutine divec          i */

/* ... files which must be loaded with vrtes.f */

/*     sphcom.f, hrfft.f, vhaes.f,shses.f */

/*     subroutine vrtes(nlat,nlon,isym,nt,vt,ivrt,jvrt,cr,ci,mdc,ndc, */
/*    +                 wshses,lshses,work,lwork,ierror) */

/*     given the vector spherical harmonic coefficients cr and ci, precomputed */
/*     by subroutine vhaes for a vector field (v,w), subroutine vrtes */
/*     computes the vorticity of the vector field in the scalar array */
/*     vt.  vt(i,j) is the vorticity at the colatitude */

/*            theta(i) = (i-1)*pi/(nlat-1) */

/*     and longitude */

/*            lambda(j) = (j-1)*2*pi/nlon */

/*     on the sphere.  i.e., */

/*            vt(i,j) =  [-dv/dlambda + d(sint*w)/dtheta]/sint */

/*     where sint = sin(theta(i)).  w is the east longitudinal and v */
/*     is the colatitudinal component of the vector field from which */
/*     cr,ci were precomputed.  required associated legendre polynomials */
/*     are stored rather than recomputed  as they are in subroutine vrtec. */


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
/*            sphere.  i.e., in the array vt(i,j) for i=1,...,nlat and */
/*            j=1,...,nlon. */

/*      = 1 */
/*            w is antisymmetric and v is symmetric about the equator. */
/*            in this case the vorticity is symmetyric about the */
/*            equator and is computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the vorticity is computed */
/*            in the array vt(i,j) for i=1,...,(nlat+1)/2 and for */
/*            j=1,...,nlon.  if nlat is even the vorticity is computed */
/*            in the array vt(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */

/*      = 2 */
/*            w is symmetric and v is antisymmetric about the equator */
/*            in this case the vorticity is antisymmetric about the */
/*            equator and is computed for the northern hemisphere */
/*            only.  i.e., if nlat is odd the vorticity is computed */
/*            in the array vt(i,j) for i=1,...,(nlat+1)/2 and for */
/*            j=1,...,nlon.  if nlat is even the vorticity is computed */
/*            in the array vt(i,j) for i=1,...,nlat/2 and j=1,...,nlon. */


/*      nt    nt is the number of scalar and vector fields.  some */
/*            computational efficiency is obtained for multiple fields. */
/*            in the program that calls vrtes, the arrays cr,ci, and vort */
/*            can be three dimensional corresponding to an indexed multiple */
/*            vector field.  in this case multiple scalar synthesis will */
/*            be performed to compute the vorticity for each field.  the */
/*            third index is the synthesis index which assumes the values */
/*            k=1,...,nt.  for a single synthesis set nt = 1.  the */
/*            description of the remaining parameters is simplified by */
/*            assuming that nt=1 or that all the arrays are two dimensional. */

/*     ivrt   the first dimension of the array vt as it appears in */
/*            the program that calls vrtes. if isym = 0 then ivrt */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then ivrt must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then ivrt must be at least (nlat+1)/2. */

/*     jvrt   the second dimension of the array vt as it appears in */
/*            the program that calls vrtes. jvrt must be at least nlon. */

/*    cr,ci   two or three dimensional arrays (see input parameter nt) */
/*            that contain vector spherical harmonic coefficients */
/*            of the vector field (v,w) as computed by subroutine vhaes. */
/*     ***    cr and ci must be computed by vhaes prior to calling */
/*            vrtes. */

/*      mdc   the first dimension of the arrays cr and ci as it */
/*            appears in the program that calls vrtes. mdc must be at */
/*            least min0(nlat,nlon/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*      ndc   the second dimension of the arrays cr and ci as it */
/*            appears in the program that calls vrtes. ndc must be at */
/*            least nlat. */

/*   wshses   an array which must be initialized by subroutine shsesi. */
/*            once initialized, */
/*            wshses can be used repeatedly by vrtes as long as nlon */
/*            and nlat remain unchanged.  wshses must not be altered */
/*            between calls of vrtes */

/*   lshses   the dimension of the array wshses as it appears in the */
/*            program that calls vrtes. define */

/*               l1 = min0(nlat,(nlon+2)/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lshses must be at least */

/*               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15 */

/*     work   a work array that does not have to be saved. */

/*    lwork   the dimension of the array work as it appears in the */
/*            program that calls vrtes. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd. */

/*            if isym = 0 then lwork must be at least */

/*               nlat*((nt+1)*nlon+2*nt*l1+1) */

/*            if isym > 0 then lwork must be at least */

/*               (nt+1)*l2*nlon+nlat*(2*nt*l1+1) */


/*     ************************************************************** */

/*     output parameters */


/*     vt     a two or three dimensional array (see input parameter nt) */
/*            that contains the vorticity of the vector field (v,w) */
/*            whose coefficients cr,ci where computed by subroutine vhaes. */
/*            vt(i,j) is the vorticity at the colatitude point theta(i) = */
/*            (i-1)*pi/(nlat-1) and longitude point lambda(j) = */
/*            (j-1)*2*pi/nlon. the index ranges are defined above at the */
/*            input parameter isym. */


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
/*          = 9  error in the specification of lshses */
/*          = 10 error in the specification of lwork */
/* ********************************************************************** */


/* Subroutine */ int vrtes_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, real *vort, integer *ivrt, integer *jvrt, real *cr, real 
	*ci, integer *mdc, integer *ndc, real *wshses, integer *lshses, real *
	work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer vort_dim1, vort_dim2, vort_offset, cr_dim1, cr_dim2, cr_offset, 
	    ci_dim1, ci_dim2, ci_offset, i__1, i__2;

    /* Local variables */
    static integer ia, ib, mn, is, ls, mab, nln, iwk, lwk, imid, mmax, lpimn;
    extern /* Subroutine */ int vrtes1_(integer *, integer *, integer *, 
	    integer *, real *, integer *, integer *, real *, real *, integer *
	    , integer *, real *, real *, integer *, real *, real *, integer *,
	     real *, integer *, integer *);


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
    --wshses;
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
    if (*lshses < lpimn + *nlon + 15) {
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
    if (*lwork < nln + ls * *nlon + (mn << 1) + *nlat) {
	return 0;
    }
    *ierror = 0;

/*     set work space pointers */

    ia = 1;
    ib = ia + mn;
    is = ib + mn;
    iwk = is + *nlat;
    lwk = *lwork - (mn << 1) - *nlat;
    vrtes1_(nlat, nlon, isym, nt, &vort[vort_offset], ivrt, jvrt, &cr[
	    cr_offset], &ci[ci_offset], mdc, ndc, &work[ia], &work[ib], &mab, 
	    &work[is], &wshses[1], lshses, &work[iwk], &lwk, ierror);
    return 0;
} /* vrtes_ */

/* Subroutine */ int vrtes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, real *vort, integer *ivrt, integer *jvrt, real *cr, real 
	*ci, integer *mdc, integer *ndc, real *a, real *b, integer *mab, real 
	*sqnn, real *wsav, integer *lwsav, real *wk, integer *lwk, integer *
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
    static real fn;
    static integer mmax;
    extern /* Subroutine */ int shses_(integer *, integer *, integer *, 
	    integer *, real *, integer *, integer *, real *, real *, integer *
	    , integer *, real *, integer *, real *, integer *, integer *);


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
	fn = (real) (n - 1);
	sqnn[n] = sqrt(fn * (fn + 1.f));
/* L1: */
    }

/*     compute divergence scalar coefficients for each vector field */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *nlat;
	for (n = 1; n <= i__2; ++n) {
	    i__3 = *mab;
	    for (m = 1; m <= i__3; ++m) {
		a[m + (n + k * a_dim2) * a_dim1] = 0.f;
		b[m + (n + k * b_dim2) * b_dim1] = 0.f;
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

    shses_(nlat, nlon, isym, nt, &vort[vort_offset], ivrt, jvrt, &a[a_offset],
	     &b[b_offset], mab, nlat, &wsav[1], lwsav, &wk[1], lwk, ierror);
    return 0;
} /* vrtes1_ */

