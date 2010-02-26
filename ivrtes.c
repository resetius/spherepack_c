/* ivrtes.f -- translated by f2c (version 20061008).
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
/*  .                         SPHEREPACK                       . */
/*  .                                                             . */
/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */



/* ... file ivrtes.f */

/*     this file includes documentation and code for */
/*     subroutine ivrtes         i */

/* ... files which must be loaded with ivrtes.f */

/*     sphcom.f, hrfft.f, vhses.f,shaes.f */


/*     subroutine ivrtes(nlat,nlon,isym,nt,v,w,idvw,jdvw,a,b,mdab,ndab, */
/*    +                  wvhses,lvhses,work,lwork,pertrb,ierror) */

/*     given the scalar spherical harmonic coefficients a and b, precomputed */
/*     by subroutine shaes for a scalar array vort, subroutine ivrtes computes */
/*     a divergence free vector field (v,w) whose vorticity is vort - pertrb. */
/*     w is the east longitude component and v is the colatitudinal component. */
/*     pertrb is a constant which must be subtracted from vort for (v,w) to */
/*     exist (see the description of pertrb below).  usually pertrb is zero */
/*     or small relative to vort.  the divergence of (v,w), as computed by */
/*     ivrtes, is the zero scalar field.  i.e., v(i,j) and w(i,j) are the */
/*     colaatitudinal and east longitude velocity components at colatitude */

/*            theta(i) = (i-1)*pi/(nlat-1) */

/*     and longitude */

/*            lambda(j) = (j-1)*2*pi/nlon. */

/*     the */

/*            vorticity(v(i,j),w(i,j)) */

/*         =  [-dv/dlambda + d(sint*w)/dtheta]/sint */

/*         =  vort(i,j) - pertrb */

/*     and */

/*            divergence(v(i,j),w(i,j)) */

/*         =  [d(sint*v)/dtheta + dw/dlambda]/sint */

/*         =  0.0 */

/*     where sint = sin(theta(i)).  required associated legendre polynomials */
/*     are stored rather than recomputed as they are in subroutine ivrtec. */


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


/*     isym   this has the same value as the isym that was input to */
/*            subroutine shaes to compute the arrays a and b.  isym */
/*            determines whether (v,w) are computed on the full or half */
/*            sphere as follows: */

/*      = 0 */
/*            vort is not symmetric about the equator. in this case */
/*            the vector field (v,w) is computed on the entire sphere. */
/*            i.e., in the arrays  v(i,j),w(i,j) for i=1,...,nlat and */
/*            j=1,...,nlon. */

/*      = 1 */

/*            vort is symmetric about the equator. in this case w is */
/*            antiymmetric and v is symmetric about the equator. v */
/*            and w are computed on the northern hemisphere only.  i.e., */
/*            if nlat is odd they are computed for i=1,...,(nlat+1)/2 */
/*            and j=1,...,nlon.  if nlat is even they are computed for */
/*            i=1,...,nlat/2 and j=1,...,nlon. */

/*       = 2 */

/*            vort is antisymmetric about the equator. in this case w is */
/*            symmetric and v is antisymmetric about the equator. w */
/*            and v are computed on the northern hemisphere only.  i.e., */
/*            if nlat is odd they are computed for i=1,...,(nlat+1)/2 */
/*            and j=1,...,nlon.  if nlat is even they are computed for */
/*            i=1,...,nlat/2 and j=1,...,nlon. */


/*     nt     in the program that calls ivrtes, nt is the number of vorticity */
/*            and vector fields.  some computational efficiency is obtained */
/*            for multiple fields.  the arrays a,b,v, and w can be three */
/*            dimensional and pertrb can be one dimensional corresponding */
/*            to an indexed multiple array vort.  in this case, multiple vector */
/*            synthesis will be performed to compute each vector field.  the */
/*            third index for a,b,v,w and first for pertrb is the synthesis */
/*            index which assumes the values k=1,...,nt.  for a single */
/*            synthesis set nt=1. the description of the remaining parameters */
/*            is simplified by assuming that nt=1 or that a,b,v,w are two */
/*            dimensional and pertrb is a constant. */

/*     idvw   the first dimension of the arrays v,w as it appears in */
/*            the program that calls ivrtes. if isym = 0 then idvw */
/*            must be at least nlat.  if isym = 1 or 2 and nlat is */
/*            even then idvw must be at least nlat/2. if isym = 1 or 2 */
/*            and nlat is odd then idvw must be at least (nlat+1)/2. */

/*     jdvw   the second dimension of the arrays v,w as it appears in */
/*            the program that calls ivrtes. jdvw must be at least nlon. */

/*     a,b    two or three dimensional arrays (see input parameter nt) */
/*            that contain scalar spherical harmonic coefficients */
/*            of the vorticity array vort as computed by subroutine shaes. */
/*     ***    a,b must be computed by shaes prior to calling ivrtes. */

/*     mdab   the first dimension of the arrays a and b as it appears in */
/*            the program that calls ivrtes (and shaes). mdab must be at */
/*            least min0(nlat,(nlon+2)/2) if nlon is even or at least */
/*            min0(nlat,(nlon+1)/2) if nlon is odd. */

/*     ndab   the second dimension of the arrays a and b as it appears in */
/*            the program that calls ivrtes (and shaes). ndab must be at */
/*            least nlat. */


/*  wvhses    an array which must be initialized by subroutine vhsesi. */
/*            once initialized */
/*            wvhses can be used repeatedly by ivrtes as long as nlon */
/*            and nlat remain unchanged.  wvhses must not be altered */
/*            between calls of ivrtes. */


/*  lvhses    the dimension of the array wvhses as it appears in the */
/*            program that calls ivrtes. define */

/*               l1 = min0(nlat,nlon/2) if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2) if nlon is odd */

/*            and */

/*               l2 = nlat/2        if nlat is even or */
/*               l2 = (nlat+1)/2    if nlat is odd */

/*            then lvhses must be at least */

/*               (l1*l2*(nlat+nlat-l1+1))/2+nlon+15 */


/*     work   a work array that does not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls ivrtes. define */

/*               l2 = nlat/2                    if nlat is even or */
/*               l2 = (nlat+1)/2                if nlat is odd */
/*               l1 = min0(nlat,nlon/2)         if nlon is even or */
/*               l1 = min0(nlat,(nlon+1)/2)     if nlon is odd */

/*            if isym = 0 then lwork must be at least */

/*               nlat*((2*nt+1)*nlon+2*l1*nt+1) */

/*            if isym = 1 or 2 then lwork must be at least */

/*               (2*nt+1)*l2*nlon+nlat*(2*nt*l1+1) */

/*     ************************************************************** */

/*     output parameters */


/*     v,w   two or three dimensional arrays (see input parameter nt) that */
/*           contain a divergence free vector field whose vorticity is */
/*           vort - pertrb at the lattitude point theta(i)=pi/2-(i-1)*pi/(nlat-1) */
/*           and longitude point lambda(j)=(j-1)*2*pi/nlon.  w is the east */
/*           longitude component and v is the colatitudinal component.  the */
/*           indices for v and w are defined at the input parameter isym. */
/*           the divergence of (v,w) is the zero scalar field. */

/*   pertrb  a nt dimensional array (see input parameter nt and assume nt=1 */
/*           for the description that follows).  vort - pertrb is a scalar */
/*           field which can be the vorticity of a vector field (v,w). */
/*           pertrb is related to the scalar harmonic coefficients a,b */
/*           of vort (computed by shaes) by the formula */

/*                pertrb = a(1,1)/(2.*sqrt(2.)) */

/*           an unperturbed vort can be the vorticity of a vector field */
/*           only if a(1,1) is zero.  if a(1,1) is nonzero (flagged by */
/*           pertrb nonzero) then subtracting pertrb from vort yields a */
/*           scalar field for which a(1,1) is zero. */

/*    ierror = 0  no errors */
/*           = 1  error in the specification of nlat */
/*           = 2  error in the specification of nlon */
/*           = 3  error in the specification of isym */
/*           = 4  error in the specification of nt */
/*           = 5  error in the specification of idvw */
/*           = 6  error in the specification of jdvw */
/*           = 7  error in the specification of mdab */
/*           = 8  error in the specification of ndab */
/*           = 9  error in the specification of lvhses */
/*           = 10 error in the specification of lwork */
/* ********************************************************************** */


/* Subroutine */ int ivrtes_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *
	jdvw, doublereal *a, doublereal *b, integer *mdab, integer *ndab, 
	doublereal *wvhses, integer *lvhses, doublereal *work, integer *lwork,
	 doublereal *pertrb, integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, a_dim1, 
	    a_dim2, a_offset, b_dim1, b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    integer l1, l2, mn, is, ici, icr, iwk, lzz1, labc, imid, mmax, liwk, 
	    lwmin;
    extern /* Subroutine */ int ivtes1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);


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
    --wvhses;
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
    lzz1 = (*nlat << 1) * imid;
/* Computing MAX */
    i__1 = mmax - 2;
    labc = max(i__1,0) * (*nlat + *nlat - mmax - 1) * 3 / 2;
/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    l1 = min(i__1,i__2);
    l2 = (*nlat + 1) / 2;
    lwmin = l1 * l2 * (*nlat + *nlat - l1 + 1) / 2 + *nlon + 15;
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

    icr = 1;
    ici = icr + mn;
    is = ici + mn;
    iwk = is + *nlat;
    liwk = *lwork - (mn << 1) - *nlat;
    ivtes1_(nlat, nlon, isym, nt, &v[v_offset], &w[w_offset], idvw, jdvw, &
	    work[icr], &work[ici], &mmax, &work[is], mdab, ndab, &a[a_offset],
	     &b[b_offset], &wvhses[1], lvhses, &work[iwk], &liwk, &pertrb[1], 
	    ierror);
    return 0;
} /* ivrtes_ */

/* Subroutine */ int ivtes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *
	jdvw, doublereal *cr, doublereal *ci, integer *mmax, doublereal *sqnn,
	 integer *mdab, integer *ndab, doublereal *a, doublereal *b, 
	doublereal *wsav, integer *lwsav, doublereal *wk, integer *lwk, 
	doublereal *pertrb, integer *ierror)
{
    /* System generated locals */
    integer v_dim1, v_dim2, v_offset, w_dim1, w_dim2, w_offset, cr_dim1, 
	    cr_dim2, cr_offset, ci_dim1, ci_dim2, ci_offset, a_dim1, a_dim2, 
	    a_offset, b_dim1, b_dim2, b_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer k, m, n;
    doublereal bi, fn, br;
    integer ityp;
    extern /* Subroutine */ int vhses_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);


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
    ci_dim1 = *mmax;
    ci_dim2 = *nlat;
    ci_offset = 1 + ci_dim1 * (1 + ci_dim2);
    ci -= ci_offset;
    cr_dim1 = *mmax;
    cr_dim2 = *nlat;
    cr_offset = 1 + cr_dim1 * (1 + cr_dim2);
    cr -= cr_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
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

/*     compute multiple vector fields coefficients */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {

/*     set vorticity field perturbation adjustment */

	pertrb[k] = a[(k * a_dim2 + 1) * a_dim1 + 1] / (sqrt(2.) * 2.);

/*     preset br,bi to 0.0 */

	i__2 = *nlat;
	for (n = 1; n <= i__2; ++n) {
	    i__3 = *mmax;
	    for (m = 1; m <= i__3; ++m) {
		cr[m + (n + k * cr_dim2) * cr_dim1] = 0.;
		ci[m + (n + k * ci_dim2) * ci_dim1] = 0.;
/* L4: */
	    }
/* L3: */
	}

/*     compute m=0 coefficients */

	i__2 = *nlat;
	for (n = 2; n <= i__2; ++n) {
	    cr[(n + k * cr_dim2) * cr_dim1 + 1] = a[(n + k * a_dim2) * a_dim1 
		    + 1] / sqnn[n];
	    ci[(n + k * ci_dim2) * ci_dim1 + 1] = b[(n + k * b_dim2) * b_dim1 
		    + 1] / sqnn[n];
/* L5: */
	}

/*     compute m>0 coefficients */

	i__2 = *mmax;
	for (m = 2; m <= i__2; ++m) {
	    i__3 = *nlat;
	    for (n = m; n <= i__3; ++n) {
		cr[m + (n + k * cr_dim2) * cr_dim1] = a[m + (n + k * a_dim2) *
			 a_dim1] / sqnn[n];
		ci[m + (n + k * ci_dim2) * ci_dim1] = b[m + (n + k * b_dim2) *
			 b_dim1] / sqnn[n];
/* L7: */
	    }
/* L6: */
	}
/* L2: */
    }

/*     set ityp for vector synthesis with divergence=0 */

    if (*isym == 0) {
	ityp = 2;
    } else if (*isym == 1) {
	ityp = 5;
    } else if (*isym == 2) {
	ityp = 8;
    }

/*     vector sythesize cr,ci into divergence free vector field (v,w) */

    vhses_(nlat, nlon, &ityp, nt, &v[v_offset], &w[w_offset], idvw, jdvw, &br,
	     &bi, &cr[cr_offset], &ci[ci_offset], mmax, nlat, &wsav[1], lwsav,
	     &wk[1], lwk, ierror);
    return 0;
} /* ivtes1_ */

