
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
/* Modifications */

//-----------------------------------------  (lapl^3 - b)u = f -------------------//

/* Subroutine */ int islapec_3_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *xlmbda, doublereal *sf, integer *ids, 
	integer *jds, doublereal *a, doublereal *b, integer *mdab, integer *
	ndab, doublereal *wshsec, integer *lshsec, doublereal *work, integer *
	lwork, doublereal *pertrb, integer *ierror)
{
    /* System generated locals */
    integer sf_dim1, sf_dim2, sf_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    integer k, l1, l2, ia, ib, mn, ls, ifn, nln, iwk, lwk, imid, mmax, lwmin, 
	    lwkmin;
    extern /* Subroutine */ int islpec1_3_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);


/*     check input parameters */

    /* Parameter adjustments */
    --pertrb;
    --xlmbda;
    sf_dim1 = *ids;
    sf_dim2 = *jds;
    sf_offset = 1 + sf_dim1 * (1 + sf_dim2);
    sf -= sf_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    --wshsec;
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
    if (*isym == 0 && *ids < *nlat || *isym > 0 && *ids < imid) {
	return 0;
    }
    *ierror = 6;
    if (*jds < *nlon) {
	return 0;
    }
    *ierror = 7;
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mmax = min(i__1,i__2);
    if (*mdab < mmax) {
	return 0;
    }
    *ierror = 8;
    if (*ndab < *nlat) {
	return 0;
    }
    *ierror = 9;

/*     set and verify saved work space length */


/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    l1 = min(i__1,i__2);
    l2 = (*nlat + 1) / 2;
    lwmin = (*nlat << 1) * l2 + (l1 - 2) * (*nlat + *nlat - l1 - 1) * 3 / 2 + 
	    *nlon + 15;
    if (*lshsec < lwmin) {
	return 0;
    }
    *ierror = 10;

/*     set and verify unsaved work space length */

    ls = *nlat;
    if (*isym > 0) {
	ls = imid;
    }
    nln = *nt * ls * *nlon;
    mn = mmax * *nlat * *nt;
/*     lwmin = nln+ls*nlon+2*mn+nlat */
/*     if (lwork .lt. lwmin) return */
    l2 = (*nlat + 1) / 2;
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    l1 = min(i__1,i__2);
    if (*isym == 0) {
/* Computing MAX */
	i__1 = l2 * 6;
	lwkmin = *nlat * ((*nt << 1) * *nlon + max(i__1,*nlon) + (*nt << 1) * 
		l1 + 1);
    } else {
/* Computing MAX */
	i__1 = *nlat * 6;
	lwkmin = l2 * ((*nt << 1) * *nlon + max(i__1,*nlon)) + *nlat * ((*nt 
		<< 1) * l1 + 1);
    }
    if (*lwork < lwkmin) {
	return 0;
    }
    *ierror = 0;

/*     check sign of xlmbda */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	if (xlmbda[k] < 0.) {
	    *ierror = -1;
	}
    }

/*     set work space pointers */

    ia = 1;
    ib = ia + mn;
    ifn = ib + mn;
    iwk = ifn + *nlat;
    lwk = *lwork - (mn << 1) - *nlat;
    islpec1_3_(nlat, nlon, isym, nt, &xlmbda[1], &sf[sf_offset], ids, jds, &a[
	    a_offset], &b[b_offset], mdab, ndab, &work[ia], &work[ib], &mmax, 
	    &work[ifn], &wshsec[1], lshsec, &work[iwk], &lwk, &pertrb[1], 
	    ierror);
    return 0;
} /* islapec_3_ */

/* Subroutine */ int islpec1_3_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *xlmbda, doublereal *sf, integer *ids, 
	integer *jds, doublereal *a, doublereal *b, integer *mdab, integer *
	ndab, doublereal *as, doublereal *bs, integer *mmax, doublereal *fnn, 
	doublereal *wshsec, integer *lshsec, doublereal *wk, integer *lwk, 
	doublereal *pertrb, integer *ierror)
{
    /* System generated locals */
    integer sf_dim1, sf_dim2, sf_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, as_dim1, as_dim2, as_offset, bs_dim1, bs_dim2, 
	    bs_offset, i__1, i__2, i__3;

    /* Local variables */
    integer k, m, n;
    doublereal fn;
    extern /* Subroutine */ int shsec_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);


/*     set multipliers and preset synthesis coefficients to zero */

    /* Parameter adjustments */
    --fnn;
    --pertrb;
    --xlmbda;
    sf_dim1 = *ids;
    sf_dim2 = *jds;
    sf_offset = 1 + sf_dim1 * (1 + sf_dim2);
    sf -= sf_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    bs_dim1 = *mmax;
    bs_dim2 = *nlat;
    bs_offset = 1 + bs_dim1 * (1 + bs_dim2);
    bs -= bs_offset;
    as_dim1 = *mmax;
    as_dim2 = *nlat;
    as_offset = 1 + as_dim1 * (1 + as_dim2);
    as -= as_offset;
    --wshsec;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 1; n <= i__1; ++n) {
	fn = (doublereal) (n - 1);
	fnn[n] = fn * (fn + 1.);
	fnn[n] *= fn * (fn +1.) * fn * (fn +1.);
	i__2 = *mmax;
	for (m = 1; m <= i__2; ++m) {
	    i__3 = *nt;
	    for (k = 1; k <= i__3; ++k) {
		as[m + (n + k * as_dim2) * as_dim1] = 0.;
		bs[m + (n + k * bs_dim2) * bs_dim1] = 0.;
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {

/*     compute synthesis coefficients for xlmbda zero or nonzero */

	if (xlmbda[k] == 0.) {
	    i__2 = *nlat;
	    for (n = 2; n <= i__2; ++n) {
		as[(n + k * as_dim2) * as_dim1 + 1] = -a[(n + k * a_dim2) * 
			a_dim1 + 1] / fnn[n];
		bs[(n + k * bs_dim2) * bs_dim1 + 1] = -b[(n + k * b_dim2) * 
			b_dim1 + 1] / fnn[n];
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    as[m + (n + k * as_dim2) * as_dim1] = -a[m + (n + k * 
			    a_dim2) * a_dim1] / fnn[n];
		    bs[m + (n + k * bs_dim2) * bs_dim1] = -b[m + (n + k * 
			    b_dim2) * b_dim1] / fnn[n];
		}
	    }
	} else {

/*     xlmbda nonzero so operator invertible unless */
/*     -n*(n-1) = xlmbda(k) < 0.0  for some n */

	    pertrb[k] = 0.;
	    i__2 = *nlat;
	    for (n = 1; n <= i__2; ++n) {
		as[(n + k * as_dim2) * as_dim1 + 1] = -a[(n + k * a_dim2) * 
			a_dim1 + 1] / (fnn[n] + xlmbda[k]);
		bs[(n + k * bs_dim2) * bs_dim1 + 1] = -b[(n + k * b_dim2) * 
			b_dim1 + 1] / (fnn[n] + xlmbda[k]);
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    as[m + (n + k * as_dim2) * as_dim1] = -a[m + (n + k * 
			    a_dim2) * a_dim1] / (fnn[n] + xlmbda[k]);
		    bs[m + (n + k * bs_dim2) * bs_dim1] = -b[m + (n + k * 
			    b_dim2) * b_dim1] / (fnn[n] + xlmbda[k]);
		}
	    }
	}
    }

/*     synthesize as,bs into sf */

    shsec_(nlat, nlon, isym, nt, &sf[sf_offset], ids, jds, &as[as_offset], &
	    bs[bs_offset], mmax, nlat, &wshsec[1], lshsec, &wk[1], lwk, 
	    ierror);
    return 0;
} /* islpec1_3_ */




//----------------------------------------- (mu2 lapl^3 + mu lapl -  b)u = f  ---------------//

/* Subroutine */ int islapec_l3_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *mu2, doublereal *mu, 
	doublereal *xlmbda, doublereal *sf, integer *ids, 
	integer *jds, doublereal *a, doublereal *b, integer *mdab, integer *
	ndab, doublereal *wshsec, integer *lshsec, doublereal *work, integer *
	lwork, doublereal *pertrb, integer *ierror)
{
    /* System generated locals */
    integer sf_dim1, sf_dim2, sf_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, i__1, i__2;

    /* Local variables */
    integer k, l1, l2, ia, ib, mn, ls, ifn, nln, iwk, lwk, imid, mmax, lwmin, 
	    lwkmin;
    extern /* Subroutine */ int islpec1_l3_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	    integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);


/*     check input parameters */

    /* Parameter adjustments */
    --pertrb;
    --xlmbda;
    sf_dim1 = *ids;
    sf_dim2 = *jds;
    sf_offset = 1 + sf_dim1 * (1 + sf_dim2);
    sf -= sf_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    --wshsec;
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
    if (*isym == 0 && *ids < *nlat || *isym > 0 && *ids < imid) {
	return 0;
    }
    *ierror = 6;
    if (*jds < *nlon) {
	return 0;
    }
    *ierror = 7;
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    mmax = min(i__1,i__2);
    if (*mdab < mmax) {
	return 0;
    }
    *ierror = 8;
    if (*ndab < *nlat) {
	return 0;
    }
    *ierror = 9;

/*     set and verify saved work space length */


/* Computing MIN */
    i__1 = *nlat, i__2 = (*nlon + 2) / 2;
    l1 = min(i__1,i__2);
    l2 = (*nlat + 1) / 2;
    lwmin = (*nlat << 1) * l2 + (l1 - 2) * (*nlat + *nlat - l1 - 1) * 3 / 2 + 
	    *nlon + 15;
    if (*lshsec < lwmin) {
	return 0;
    }
    *ierror = 10;

/*     set and verify unsaved work space length */

    ls = *nlat;
    if (*isym > 0) {
	ls = imid;
    }
    nln = *nt * ls * *nlon;
    mn = mmax * *nlat * *nt;
/*     lwmin = nln+ls*nlon+2*mn+nlat */
/*     if (lwork .lt. lwmin) return */
    l2 = (*nlat + 1) / 2;
/* Computing MIN */
    i__1 = *nlat, i__2 = *nlon / 2 + 1;
    l1 = min(i__1,i__2);
    if (*isym == 0) {
/* Computing MAX */
	i__1 = l2 * 6;
	lwkmin = *nlat * ((*nt << 1) * *nlon + max(i__1,*nlon) + (*nt << 1) * 
		l1 + 1);
    } else {
/* Computing MAX */
	i__1 = *nlat * 6;
	lwkmin = l2 * ((*nt << 1) * *nlon + max(i__1,*nlon)) + *nlat * ((*nt 
		<< 1) * l1 + 1);
    }
    if (*lwork < lwkmin) {
	return 0;
    }
    *ierror = 0;

/*     check sign of xlmbda */

    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {
	if (xlmbda[k] < 0.) {
	    *ierror = -1;
	}
    }

/*     set work space pointers */

    ia = 1;
    ib = ia + mn;
    ifn = ib + mn;
    iwk = ifn + *nlat;
    lwk = *lwork - (mn << 1) - *nlat;
    islpec1_l3_(nlat, nlon, isym, nt, mu2, mu, &xlmbda[1], &sf[sf_offset], ids, jds, &a[
	    a_offset], &b[b_offset], mdab, ndab, &work[ia], &work[ib], &mmax, 
	    &work[ifn], &wshsec[1], lshsec, &work[iwk], &lwk, &pertrb[1], 
	    ierror);
    return 0;
} /* islapec_l3_ */

/* Subroutine */ int islpec1_l3_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal * mu2, doublereal *mu, 
	doublereal *xlmbda, doublereal *sf, integer *ids, 
	integer *jds, doublereal *a, doublereal *b, integer *mdab, integer *
	ndab, doublereal *as, doublereal *bs, integer *mmax, doublereal *fnn, 
	doublereal *wshsec, integer *lshsec, doublereal *wk, integer *lwk, 
	doublereal *pertrb, integer *ierror)
{
    /* System generated locals */
    integer sf_dim1, sf_dim2, sf_offset, a_dim1, a_dim2, a_offset, b_dim1, 
	    b_dim2, b_offset, as_dim1, as_dim2, as_offset, bs_dim1, bs_dim2, 
	    bs_offset, i__1, i__2, i__3;

    /* Local variables */
    integer k, m, n;
    doublereal fn;
    extern /* Subroutine */ int shsec_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);


/*     set multipliers and preset synthesis coefficients to zero */

    /* Parameter adjustments */
    --fnn;
    --pertrb;
    --xlmbda;
    sf_dim1 = *ids;
    sf_dim2 = *jds;
    sf_offset = 1 + sf_dim1 * (1 + sf_dim2);
    sf -= sf_offset;
    b_dim1 = *mdab;
    b_dim2 = *ndab;
    b_offset = 1 + b_dim1 * (1 + b_dim2);
    b -= b_offset;
    a_dim1 = *mdab;
    a_dim2 = *ndab;
    a_offset = 1 + a_dim1 * (1 + a_dim2);
    a -= a_offset;
    bs_dim1 = *mmax;
    bs_dim2 = *nlat;
    bs_offset = 1 + bs_dim1 * (1 + bs_dim2);
    bs -= bs_offset;
    as_dim1 = *mmax;
    as_dim2 = *nlat;
    as_offset = 1 + as_dim1 * (1 + as_dim2);
    as -= as_offset;
    --wshsec;
    --wk;

    /* Function Body */
    i__1 = *nlat;
    for (n = 1; n <= i__1; ++n) {
	fn = (doublereal) (n - 1);
	fnn[n] = fn * (fn + 1.);
	i__2 = *mmax;
	for (m = 1; m <= i__2; ++m) {
	    i__3 = *nt;
	    for (k = 1; k <= i__3; ++k) {
		as[m + (n + k * as_dim2) * as_dim1] = 0.;
		bs[m + (n + k * bs_dim2) * bs_dim1] = 0.;
	    }
	}
    }
    i__1 = *nt;
    for (k = 1; k <= i__1; ++k) {

/*     compute synthesis coefficients for xlmbda zero or nonzero */

	if (xlmbda[k] == 0.) {
	    i__2 = *nlat;
	    for (n = 2; n <= i__2; ++n) {
		as[(n + k * as_dim2) * as_dim1 + 1] = -a[(n + k * a_dim2) * 
			a_dim1 + 1] / 
	((*mu2)*fnn[n]*fnn[n]*fnn[n] + (*mu)*fnn[n]);
		bs[(n + k * bs_dim2) * bs_dim1 + 1] = -b[(n + k * b_dim2) * 
			b_dim1 + 1] / 
	((*mu2)*fnn[n]*fnn[n]*fnn[n] + (*mu)*fnn[n]);
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    as[m + (n + k * as_dim2) * as_dim1] = -a[m + (n + k * 
			    a_dim2) * a_dim1] / 
	((*mu2)*fnn[n]*fnn[n]*fnn[n] + (*mu)*fnn[n]);
		    bs[m + (n + k * bs_dim2) * bs_dim1] = -b[m + (n + k * 
			    b_dim2) * b_dim1] /
	((*mu2)*fnn[n]*fnn[n]*fnn[n] + (*mu)*fnn[n]);
		}
	    }
	} else {

/*     xlmbda nonzero so operator invertible unless */
/*     -n*(n-1) = xlmbda(k) < 0.0  for some n */

	    pertrb[k] = 0.;
	    i__2 = *nlat;
	    for (n = 1; n <= i__2; ++n) {
		as[(n + k * as_dim2) * as_dim1 + 1] = -a[(n + k * a_dim2) * 
			a_dim1 + 1] / 
	((*mu2)*fnn[n]*fnn[n]*fnn[n] + (*mu)*fnn[n] + xlmbda[k]);
		bs[(n + k * bs_dim2) * bs_dim1 + 1] = -b[(n + k * b_dim2) * 
			b_dim1 + 1] / 
	((*mu2)*fnn[n]*fnn[n]*fnn[n] + (*mu)*fnn[n] + xlmbda[k]);
	    }
	    i__2 = *mmax;
	    for (m = 2; m <= i__2; ++m) {
		i__3 = *nlat;
		for (n = m; n <= i__3; ++n) {
		    as[m + (n + k * as_dim2) * as_dim1] = -a[m + (n + k * 
			    a_dim2) * a_dim1] / 
	((*mu2)*fnn[n]*fnn[n]*fnn[n] + (*mu)*fnn[n] + xlmbda[k]);
		    bs[m + (n + k * bs_dim2) * bs_dim1] = -b[m + (n + k * 
			    b_dim2) * b_dim1] / 
	((*mu2)*fnn[n]*fnn[n]*fnn[n] + (*mu)*fnn[n] + xlmbda[k]);
		}
	    }
	}
    }

/*     synthesize as,bs into sf */

    shsec_(nlat, nlon, isym, nt, &sf[sf_offset], ids, jds, &as[as_offset], &
	    bs[bs_offset], mmax, nlat, &wshsec[1], lshsec, &wk[1], lwk, 
	    ierror);
    return 0;
} /* islpec1_l3_ */



