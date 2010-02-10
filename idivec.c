/* idivec.f -- translated by f2c (version 20061008).
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



/* ... file idivec.f */

/*     this file includes documentation and code for */
/*     subroutine idivec          i */

/* ... files which must be loaded with idivec.f */

/*     sphcom.f, hrfft.f, vhsec.f,shaec.f */



/*     subroutine idivec(nlat,nlon,isym,nt,v,w,idvw,jdvw,a,b,mdab,ndab, */
/*    +                  wvhsec,lvhsec,work,lwork,pertrb,ierror) */

/*     given the scalar spherical harmonic coefficients a and b, precomputed */
/*     by subroutine shaec for a scalar array dv, subroutine idivec computes */
/*     an irrotational vector field (v,w) whose divergence is dv - pertrb. */
/*     w is the east longitude component and v is the colatitudinal component. */
/*     pertrb is a constant which must be subtracted from dv for (v,w) to */
/*     exist (see the description of pertrb below).  usually pertrb is zero */
/*     or small relative to dv.  the vorticity of (v,w), as computed by */
/*     vortec, is the zero scalar field.  v(i,j) and w(i,j) are the */
/*     velocity components at colatitude */

/*            theta(i) = (i-1)*pi/(nlat-1) */

/*     and longitude */

/*            lambda(j) = (j-1)*2*pi/nlon. */

/*     the */

/*            divergence[v(i,j),w(i,j)] */

/*          = [d(w(i,j)/dlambda + d(sint*v(i,j))/dtheta]/sint */

/*          = dv(i,j) - pertrb */

/*     and */

/*            vorticity(v(i,j),w(i,j)) */

/*         =  [dv/dlambda - d(sint*w)/dtheta]/sint */

/*         =  0.0 */

/*     where sint = sin(theta(i)).  required associated legendre polynomials */
/*     are recomputed rather than stored as they are in subroutine idives. */

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
/*            nlon = 72 for a five degree grid. nlon must be greater than */
/*            3.  the efficiency of the computation is improved when nlon */
/*            is a product of small prime numbers. */


/*     isym   this has the same value as the isym that was input to */
/*            subroutine shaec to compute the arrays a and b from the */
/*            scalar field dv.  isym determines whether (v,w) are */
/*            computed on the full or half sphere as follows: */

/*      = 0 */

/*           dv is not symmetric about the equator. in this case */
/*           the vector field (v,w) is computed on the entire sphere. */
/*           i.e., in the arrays  v(i,j),w(i,j) for i=1,...,nlat and */
/*           j=1,...,nlon. */

/*      = 1 */

/*           dv is antisymmetric about the equator. in this case w is */
/*           antisymmetric and v is symmetric about the equator. w */
/*           and v are computed on the northern hemisphere only.  i.e., */
/*           if nlat is odd they are computed for i=1,...,(nlat+1)/2 */
/*           and j=1,...,nlon.  if nlat is even they are computed for */
/*           i=1,...,nlat/2 and j=1,...,nlon. */

/*      = 2 */

/*           dv is symmetric about the equator. in this case w is */
/*           symmetric and v is antisymmetric about the equator. w */
/*           and v are computed on the northern hemisphere only.  i.e., */
/*           if nlat is odd they are computed for i=1,...,(nlat+1)/2 */
/*           and j=1,...,nlon.  if nlat is even they are computed for */
/*           i=1,...,nlat/2 and j=1,...,nlon. */


/*     nt     nt is the number of divergence and vector fields.  some */
/*            computational efficiency is obtained for multiple fields. */
/*            the arrays a,b,v, and w can be three dimensional and pertrb */
/*            can be one dimensional corresponding to an indexed multiple */
/*            array dv.  in this case, multiple vector synthesis will be */
/*            performed to compute each vector field.  the third index for */
/*            a,b,v,w and first for pertrb is the synthesis index which */
/*            assumes the values k = 1,...,nt.  for a single synthesis set */
/*            nt = 1.  the description of the remaining parameters is */
/*            simplified by assuming that nt=1 or that a,b,v,w are two */
/*            dimensional and pertrb is a constant. */

/*     idvw   the first dimension of the arrays v,w as it appears in */
/*            the program that calls idivec. if isym = 0 then idvw */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then idvw must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then idvw must be at least (nlat+1)/2. */

/*     jdvw   the second dimension of the arrays v,w as it appears in */
/*            the program that calls idivec. jdvw must be at least nlon. */

/*     a,b    two or three dimensional arrays (see input parameter nt) */
/*            that contain scalar spherical harmonic coefficients */
/*            of the divergence array dv as computed by subroutine shaec. */
/*     ***    a,b must be computed by shaec prior to calling idivec. */

/*     mdab   the first dimension of the arrays a and b as it appears in */
/*            the program that calls idivec (and shaec). mdab must be at */
/*            least min0(nlat,(nlon+2)/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndab   the second dimension of the arrays a and b as it appears in */
/*            the program that calls idivec (and shaec). ndab must be at */
/*            least nlat. */


/*  wvhsec    an array which must be initialized by subroutine vhseci. */
/*            once initialized, */
/*            wvhsec can be used repeatedly by idivec as long as nlon */
/*            and nlat remain unchanged.  wvhsec must not be altered */
/*            between calls of idivec. */


/*  lvhsec    the dimension of the array wvhsec as it appears in the */
/*            program that calls idivec. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lvhsec must be at least */

/*            4*nlat*l2+3*max0(l1-2,0)*(nlat+nlat-l1-1)+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls idivec. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2                  if nlat is even or */
/*               l2 = (nlat+1)/2              if nlat is odd */

/*            if isym = 0 then lwork must be at least */

/*               nlat*(2*nt*nlon+max0(6*l2,nlon) + 2*nt*l1 + 1) */

/*            if isym = 1 or 2 then lwork must be at least */

/*               l2*(2*nt*nlon+max0(6*nlat,nlon)) + nlat*(2*l1*nt+1) */

/*     ************************************************************** */

/*     output parameters */


/*     v,w   two or three dimensional arrays (see input parameter nt) that */
/*           contain an irrotational vector field whose divergence is */
/*           dv-pertrb at the colatitude point theta(i)=(i-1)*pi/(nlat-1) */
/*           and longitude point lambda(j)=(j-1)*2*pi/nlon.  w is the east */
/*           longitude component and v is the colatitudinal component.  the */
/*           indices for w and v are defined at the input parameter isym. */
/*           the curl or vorticity of (v,w) is the zero vector field.  note */
/*           that any nonzero vector field on the sphere will be multiple */
/*           valued at the poles [reference swarztrauber]. */

/*   pertrb  a nt dimensional array (see input parameter nt and assume nt=1 */
/*           for the description that follows).  dv - pertrb is a scalar */
/*           field which can be the divergence of a vector field (v,w). */
/*           pertrb is related to the scalar harmonic coefficients a,b */
/*           of dv (computed by shaec) by the formula */

/*                pertrb = a(1,1)/(2.*sqrt(2.)) */



/*           the unperturbed scalar field dv can be the divergence of a */
/*           vector field only if a(1,1) is zero.  if a(1,1) is nonzero */
/*           (flagged by pertrb nonzero) then subtracting pertrb from */
/*           dv yields a scalar field for which a(1,1) is zero. */

/*    ierror = 0  no errors */
/*           = 1  error in the specification of nlat */
/*           = 2  error in the specification of nlon */
/*           = 3  error in the specification of isym */
/*           = 4  error in the specification of nt */
/*           = 5  error in the specification of idvw */
/*           = 6  error in the specification of jdvw */
/*           = 7  error in the specification of mdab */
/*           = 8  error in the specification of ndab */
/*           = 9  error in the specification of lvhsec */
/*           = 10 error in the specification of lwork */
/* ********************************************************************** */


/* Subroutine */ int idivec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhsec, integer *lvhsec, 
	doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, a_dim1, 
	    a_dim2, a_offset, b_dim1, b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    static integer l1, l2, mn, is, ibi, ibr, idz, iwk, imid, mmax, liwk, 
	    lwmin, lzimn;
    extern /* Subroutine */ int idvec1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);


/*     check input parameters */

    /* Parameter adjustments */
    --pertrb;
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
    --wvhsec;
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
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 1) / 2;
    l1 = min(i__1,i__2);
    l2 = (*nlat + 1) / 2;
/* Computing MAX */
    i__1 = l1 - 2;
    lwmin = (*nlat << 2) * l2 + max(i__1,0) * 3 * (*nlat + *nlat - l1 - 1) + *
	    nlon + 15;
    if (*lvhsec < lwmin) {
	return 0;
    }
    *ierror = 10;

/*     verify unsaved work space length */

    mn = mmax * *nlat * *nt;
/* Computing MAX */
    i__1 = imid * 6;
    if (*isym != 0 && *lwork < *nlat * ((*nt << 1) * *nlon + max(i__1,*nlon)) 
	    + (mn << 1) + *nlat) {
	return 0;
    }
/* Computing MAX */
    i__1 = *nlat * 6;
    if (*isym == 0 && *lwork < imid * ((*nt << 1) * *nlon + max(i__1,*nlon)) 
	    + (mn << 1) + *nlat) {
	return 0;
    }
    *ierror = 0;

/*     set work space pointers */

    ibr = 1;
    ibi = ibr + mn;
    is = ibi + mn;
    iwk = is + *nlat;
    liwk = *lwork - (mn << 1) - *nlat;
    idvec1_(nlat, nlon, isym, nt, &v[v_offset], &w[w_offset], idvw, jdvw, &
	    work[ibr], &work[ibi], &mmax, &work[is], mdab, ndab, &a[a_offset],
	     &b[b_offset], &wvhsec[1], lvhsec, &work[iwk], &liwk, &pertrb[1], 
	    ierror);
    return 0;
} /* idivec_ */

/* Subroutine */ int idvec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wvhsec, integer *lvhsec, doublereal *wk, integer *
	lwk, doublereal *pertrb, integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, br_dim1, 
	    br_dim2, br_offset, bi_dim1, bi_dim2, bi_offset, a_dim1, a_dim2, 
	    a_offset, b_dim1, b_dim2, b_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, m, n;
    static doublereal ci, fn, cr;
    static integer ityp;
    extern /* Subroutine */ int vhsec_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);

/*     preset coefficient multiplyers in vector */

    /* Parameter adjustments */
    --sqnn;
    --pertrb;
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
    --wvhsec;
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

/*     set divergence field perturbation adjustment */

	pertrb[k] = a[(k * a_dim2 + 1) * a_dim1 + 1] / (sqrt(2.) * 2.);

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
	    br[(n + k * br_dim2) * br_dim1 + 1] = -a[(n + k * a_dim2) * 
		    a_dim1 + 1] / sqnn[n];
	    bi[(n + k * bi_dim2) * bi_dim1 + 1] = -b[(n + k * b_dim2) * 
		    b_dim1 + 1] / sqnn[n];
/* L5: */
	}

/*     compute m>0 coefficients */

	i__2 = *mmax;
	for (m = 2; m <= i__2; ++m) {
	    i__3 = *nlat;
	    for (n = m; n <= i__3; ++n) {
		br[m + (n + k * br_dim2) * br_dim1] = -a[m + (n + k * a_dim2) 
			* a_dim1] / sqnn[n];
		bi[m + (n + k * bi_dim2) * bi_dim1] = -b[m + (n + k * b_dim2) 
			* b_dim1] / sqnn[n];
/* L7: */
	    }
/* L6: */
	}
/* L2: */
    }

/*     set ityp for vector synthesis with curl=0 */

    if (*isym == 0) {
	ityp = 1;
    } else if (*isym == 1) {
	ityp = 4;
    } else if (*isym == 2) {
	ityp = 7;
    }

/*     vector sythesize br,bi into irrotational (v,w) */

    vhsec_(nlat, nlon, &ityp, nt, &v[v_offset], &w[w_offset], idvw, jdvw, &br[
	    br_offset], &bi[bi_offset], &cr, &ci, mmax, nlat, &wvhsec[1], 
	    lvhsec, &wk[1], lwk, ierror);
    return 0;
} /* idvec1_ */

