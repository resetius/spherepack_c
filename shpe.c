/* shpe.f -- translated by f2c (version 20061008).
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

static integer c__10 = 10;
static integer c__0 = 0;
static doublereal c_b54 = 2.;
static integer c__1 = 1;
static doublereal c_b113 = -1.;
static doublereal c_b191 = 1.;


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

/*                            August 2003 */

/* ... file shpe.f */

/*     this file contains code and documentation for subroutines */
/*     shpei and shpe. */

/* ... files which must be loaded with shpe.f */

/*     hrfft.f */

/*     subroutine shpei initializes arrays wshp and iwshp for */
/*     subsequent repeated use by subroutine shpe, which */
/*     performs the harmonic projection equivalent to a */
/*     harmonic analysis followed by harmonic synthesis */
/*     but faster and with less memory. (see description of */
/*     subroutine shpe below) */

/*     subroutine shpei(nlat,nlon,isym,mtrunc,wshp,lwshp,iwshp, */
/*    1 liwshp,work,lwork,ierror) */

/*     shpei initializes arrays wshp and iwshp for repeated use */
/*     by subroutine shpe .... */

/*     input parameters */

/*     nlat   the number of colatitudes on the full sphere including the */
/*            poles. for example, nlat = 37 for a five degree grid. */
/*            nlat determines the grid increment in colatitude as */
/*            pi/(nlat-1).  if nlat is odd the equator is located at */
/*            grid point i=(nlat+1)/2. if nlat is even the equator is */
/*            located half way between points i=nlat/2 and i=nlat/2+1. */
/*            nlat must be at least 3. */

/*     nlon   the number of distinct londitude points.  nlon determines */
/*            the grid increment in longitude as 2*pi/nlon. for example */
/*            nlon = 72 for a five degree grid. nlon must be greater */
/*            than or equal to 4. the efficiency of the computation is */
/*            improved when nlon is a product of small prime numbers. */
/*            nlon must beat least 4. */

/*     isym   currently not used. */

/*     mtrunc the highest longitudinal wave number retained in the */
/*            projection. It must be less than or equal to */
/*            the minimum of nlat-1 and nlon/2. The first wave */
/*            number is zero. For example, if wave numbers 0 and */
/*            1 are desired then mtrunc = 1. */

/*     lwshp  the dimension of the array wshp as it appears in the */
/*            program that calls shpei. It must be at least */
/*            2*(nlat+1)**2+nlon+log2(nlon) */

/*     liwshp the dimension of the array iwshp as it appears in the */
/*            program that calls shpei. It must be at least */
/*            4*(nlat+1). */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls shpei. It must be at least */
/*            1.25*(nlat+1)**2+7*nlat+8. */

/*     ************************************************************** */

/*     output parameters */

/*     wshp   a single precision array that must be saved for */
/*            repeated use by subroutine shpe. */

/*     iwshp  an integer array that must be saved for repeated */
/*            use by subroutine shpe. */

/*     work   a double precision work array that does */
/*            not have to be saved. */

/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of isym */
/*            = 4  error in the specification of mtrunc */
/*            = 5  error in the specification of lwshp */
/*            = 6  error in the specification of liwshp */
/*            = 7  error in the specification of lwork */

/* Subroutine */ int shpei_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, doublereal *wshp, integer *lwshp, integer *iwshp, integer *
	liwshp, doublereal *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer iw1, iw2, iw3, lw1, iw4, jw1, jw2, jw3, jw4, kw1, kw2, kw3,
	     kw4, kw5, kw6, kw7, kw8, kw9, kw10, kw11, kw12, nte, kw13, mmax, 
	    mlwk, nloc1, nloc2, log2n;
    extern /* Subroutine */ int shpei1_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), hrffti_(integer *, doublereal *);


    /* Parameter adjustments */
    --work;
    --iwshp;
    --wshp;

    /* Function Body */
    *ierror = 1;
    if (*nlat < 3) {
	return 0;
    }
    *ierror = 2;
    if (*nlon < 4) {
	return 0;
    }
/*      ierror = 3 */
/*      if(isym.lt.0 .or. isym.gt.2) return */
    *ierror = 4;
/* Computing MIN */
    i__1 = *nlat - 1, i__2 = *nlon / 2;
    mmax = min(i__1,i__2);
    if (*mtrunc < 0 || *mtrunc > mmax) {
	return 0;
    }
    *ierror = 5;
/* Computing 2nd power */
    i__1 = *nlat + 1;
    lw1 = i__1 * i__1 << 1;
    log2n = log((doublereal) (*nlon)) / log(2.);
    if (*lwshp < lw1 + *nlon + log2n) {
	return 0;
    }
    *ierror = 6;
    if (*liwshp < *nlat + 1 << 2) {
	return 0;
    }
    *ierror = 7;
/* Computing 2nd power */
    i__1 = *nlat + 1;
    mlwk = i__1 * i__1 * 1.25 + *nlat * 7 + 8;
    if (*lwork < mlwk) {
	return 0;
    }
    *ierror = 0;

    hrffti_(nlon, &wshp[lw1 + 1]);

    nte = (*nlat + 1) / 2;
    nloc1 = (nte << 1) * nte;
    nloc2 = *nlat + 1;
    iw1 = 1;
    iw2 = iw1 + nloc1;
    iw3 = iw2 + nloc1;
    iw4 = iw3 + nloc1;
    jw1 = 1;
    jw2 = jw1 + nloc2;
    jw3 = jw2 + nloc2;
    jw4 = jw3 + nloc2;
    kw1 = 1;
    kw2 = kw1 + nte;
    kw3 = kw2 + nte;
    kw4 = kw3 + nte;
    kw5 = kw4 + nte + 1;
    kw6 = kw5 + nte;
    kw7 = kw6 + nte;
    kw8 = kw7 + nte;
    kw9 = kw8 + nte;
    kw10 = kw9 + nloc2 + nloc2;
    kw11 = kw10 + nloc2;
    kw12 = kw11 + nloc1;
    kw13 = kw12 + nloc1;

    shpei1_(nlat, nlon, isym, mtrunc, &nte, ierror, &wshp[iw1], &wshp[iw2], &
	    wshp[iw3], &wshp[iw4], &iwshp[jw1], &iwshp[jw2], &iwshp[jw3], &
	    iwshp[jw4], &work[kw1], &work[kw2], &work[kw3], &work[kw4], &work[
	    kw5], &work[kw6], &work[kw7], &work[kw8], &work[kw9], &work[kw10],
	     &work[kw11], &work[kw12], &work[kw11], &work[kw12], &work[kw13]);
    return 0;
} /* shpei_ */

/* Subroutine */ int shpei1_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, integer *idp, integer *ierror, doublereal *pe, doublereal *po, 
	doublereal *ze, doublereal *zo, integer *ipse, integer *jzse, integer *ipso, 
	integer *jzso, doublereal *cp, doublereal *work, doublereal *wx, 
	doublereal *s, doublereal *e, doublereal *thet, doublereal *xx, 
	doublereal *z__, doublereal *a, doublereal *b, doublereal *we, 
	doublereal *ped, doublereal *wo, doublereal *pod, doublereal *u)
{
    /* System generated locals */
    integer u_dim1, u_offset, we_dim1, we_dim2, we_offset, ped_dim1, ped_dim2,
	     ped_offset, wo_dim1, wo_dim2, wo_offset, pod_dim1, pod_dim2, 
	    pod_offset, pe_dim1, pe_dim2, pe_offset, po_dim1, po_dim2, 
	    po_offset, ze_dim1, ze_dim2, ze_offset, zo_dim1, zo_dim2, 
	    zo_offset, ipse_dim1, ipse_offset, jzse_dim1, jzse_offset, 
	    ipso_dim1, ipso_offset, jzso_dim1, jzso_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, m, n;
    static doublereal v, a1, b1, c1, pi;
    static integer ip;
    extern /* Subroutine */ int gs_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    static integer it, js, mp1, ms2, ns2;
    static doublereal dfn;
    static integer nem, nte;
    static doublereal toe;
    static integer nom, nto;
    static doublereal sum;
    extern doublereal rand_(integer *);
    static integer lock, info, nshe[2], modn, nsho[2];
    static doublereal tusl;
    static integer mxtr;
    extern /* Subroutine */ int dlfkp_(integer *, integer *, doublereal *), 
	    dsvdc_(doublereal *, integer *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);
    static doublereal dthet;
    static integer mrank;
    extern /* Subroutine */ int dlftp_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static integer nrank;
    extern /* Subroutine */ int trunc_(integer *, integer *, integer *, 
	    doublereal *, integer *, integer *), normal_(integer *, 
	    doublereal *, integer *, doublereal *);




    /* Parameter adjustments */
    u_dim1 = *idp;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    pod_dim1 = *idp;
    pod_dim2 = *idp;
    pod_offset = 1 + pod_dim1 * (1 + pod_dim2);
    pod -= pod_offset;
    wo_dim1 = *idp;
    wo_dim2 = *idp;
    wo_offset = 1 + wo_dim1 * (1 + wo_dim2);
    wo -= wo_offset;
    ped_dim1 = *idp;
    ped_dim2 = *idp;
    ped_offset = 1 + ped_dim1 * (1 + ped_dim2);
    ped -= ped_offset;
    we_dim1 = *idp;
    we_dim2 = *idp;
    we_offset = 1 + we_dim1 * (1 + we_dim2);
    we -= we_offset;
    --b;
    --a;
    --z__;
    --xx;
    --thet;
    --e;
    --s;
    --wx;
    --work;
    --cp;
    jzso_dim1 = *idp;
    jzso_offset = 1 + jzso_dim1;
    jzso -= jzso_offset;
    ipso_dim1 = *idp;
    ipso_offset = 1 + ipso_dim1;
    ipso -= ipso_offset;
    jzse_dim1 = *idp;
    jzse_offset = 1 + jzse_dim1;
    jzse -= jzse_offset;
    ipse_dim1 = *idp;
    ipse_offset = 1 + ipse_dim1;
    ipse -= ipse_offset;
    zo_dim1 = *idp;
    zo_dim2 = *idp;
    zo_offset = 1 + zo_dim1 * (1 + zo_dim2);
    zo -= zo_offset;
    ze_dim1 = *idp;
    ze_dim2 = *idp;
    ze_offset = 1 + ze_dim1 * (1 + ze_dim2);
    ze -= ze_offset;
    po_dim1 = *idp;
    po_dim2 = *idp;
    po_offset = 1 + po_dim1 * (1 + po_dim2);
    po -= po_offset;
    pe_dim1 = *idp;
    pe_dim2 = *idp;
    pe_offset = 1 + pe_dim1 * (1 + pe_dim2);
    pe -= pe_offset;

    /* Function Body */
    ns2 = *nlat / 2;
    modn = *nlat - ns2 - ns2;
    nte = (*nlat + 1) / 2;
    nto = *nlat - nte;
    tusl = 0.;
    toe = 0.;

/*     compute grid distribution */

    pi = atan(1.) * 4.;
    dthet = pi / (*nlat - 1);
    i__1 = nte;
    for (i__ = 1; i__ <= i__1; ++i__) {
	thet[i__] = (i__ - 1) * dthet;
    }

/*     compute weight matrices for even functions */

    for (mp1 = 1; mp1 <= 2; ++mp1) {
	m = mp1 - 1;
	mrank = *nlat - m - m;
	nem = (mrank + 1) / 2;
	i__1 = nem;
	for (j = 1; j <= i__1; ++j) {
	    n = j + j + m - 2;
	    dlfkp_(&m, &n, &cp[1]);
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dlftp_(&m, &n, &thet[i__], &cp[1], &ped[i__ + (j + mp1 * 
			ped_dim2) * ped_dim1]);
	    }
	    if (m > 0) {
		ped[(j + mp1 * ped_dim2) * ped_dim1 + 1] = 0.;
	    }
	}
	dsvdc_(&ped[m + 1 + (mp1 * ped_dim2 + 1) * ped_dim1], idp, &nem, &nem,
		 &s[1], &e[1], &u[u_offset], idp, &v, idp, &work[1], &c__10, &
		info);

	i__1 = nem;
	for (j = 1; j <= i__1; ++j) {
	    s[j] = 1. / (s[j] * s[j]);
	}

/*     compute weight matrix as u  s sup -2 u transpose */

	i__1 = nte;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		we[i__ + (j + mp1 * we_dim2) * we_dim1] = 0.;
	    }
	}
	i__1 = nem;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = nem;
	    for (j = 1; j <= i__2; ++j) {
		sum = 0.;
		i__3 = nem;
		for (k = 1; k <= i__3; ++k) {
		    sum += s[k] * u[i__ + k * u_dim1] * u[j + k * u_dim1];
		}
		we[i__ + m + (j + m + mp1 * we_dim2) * we_dim1] = sum;
	    }
	}
/* L40: */
    }
    we[((we_dim2 << 1) + 1) * we_dim1 + 1] = 1.;

/*     compute n**2 basis (even functions) */

    i__1 = *nlat + *nlat - 2;
    for (n = 1; n <= i__1; ++n) {
	dfn = (doublereal) n;
	a[n] = sqrt(dfn * (dfn + 1.));
    }
    i__1 = *nlat - 1;
    for (n = 1; n <= i__1; ++n) {
	dfn = (doublereal) n;
	b[n] = sqrt((dfn + dfn + 3.) / (dfn + dfn - 1.));
    }

/* Computing MIN */
    i__1 = *nlat - 1, i__2 = *nlon / 2, i__1 = min(i__1,i__2);
    mxtr = min(i__1,*mtrunc);
    ip = 2;
    i__1 = mxtr + 1;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	m = mp1 - 1;
	ip = 3 - ip;
	ms2 = mp1 / 2;
	nrank = ms2 + ms2;
	mrank = *nlat - nrank;
	nem = (mrank + 1) / 2;

/*     compute associated legendre functions */

	if (m <= 1) {
	    i__2 = nem;
	    for (j = 1; j <= i__2; ++j) {
		n = j + j + m - 2;
		dlfkp_(&m, &n, &cp[1]);
		i__3 = nte;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    dlftp_(&m, &n, &thet[i__], &cp[1], &ped[i__ + (j + ms2 + 
			    ip * ped_dim2) * ped_dim1]);
		}
/* L202: */
		if (m > 0) {
		    ped[(j + ms2 + ip * ped_dim2) * ped_dim1 + 1] = 0.;
		}
/* L205: */
	    }

	} else {

	    i__2 = nem;
	    for (j = 1; j <= i__2; ++j) {
		n = j + j + m - 2;
		if (m > 1 && n > mxtr) {
		    i__3 = nte;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			u[i__ + (j + ms2) * u_dim1] = ped[i__ + (j + ms2 + ip 
				* ped_dim2) * ped_dim1];
		    }
		    goto L207;
		}
		a1 = b[n - 1] * a[n + m - 3] / a[n + m - 1];
		b1 = a[n - m + 1] / a[n + m - 1];
		if (n - m <= 1) {
		    i__3 = nte;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			u[i__ + (j + ms2) * u_dim1] = a1 * ped[i__ + (j + ms2 
				- 1 + ip * ped_dim2) * ped_dim1] - b1 * ped[
				i__ + (j + ms2 + ip * ped_dim2) * ped_dim1];
		    }
		} else {
		    c1 = b[n - 1] * a[n - m - 1] / a[n + m - 1];
		    i__3 = nte;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			u[i__ + (j + ms2) * u_dim1] = a1 * ped[i__ + (j + ms2 
				- 1 + ip * ped_dim2) * ped_dim1] - b1 * ped[
				i__ + (j + ms2 + ip * ped_dim2) * ped_dim1] + 
				c1 * u[i__ + (j + ms2 - 1) * u_dim1];
		    }
		}
L207:
		;
	    }
	    i__2 = nem;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = nte;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    ped[i__ + (j + ms2 + ip * ped_dim2) * ped_dim1] = u[i__ + 
			    (j + ms2) * u_dim1];
		}
	    }
	}

	if (ms2 <= 0 || ms2 >= nte) {
	    goto L200;
	}
	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xx[i__] = rand_(&c__0);
	}
	it = 0;
L201:
	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] = 0.;
	    wx[i__] = 0.;
	    i__3 = nte;
	    for (j = 1; j <= i__3; ++j) {
		wx[i__] += we[i__ + (j + ip * we_dim2) * we_dim1] * xx[j];
	    }
	}
	i__2 = nte;
	for (j = 1; j <= i__2; ++j) {
	    if (j == ms2) {
		goto L220;
	    }
	    gs_(&nte, &wx[1], &ped[(j + ip * ped_dim2) * ped_dim1 + 1], &z__[
		    1]);
L220:
	    ;
	}

	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xx[i__] -= z__[i__];
	}
	normal_(&nte, &xx[1], idp, &we[(ip * we_dim2 + 1) * we_dim1 + 1]);
	++it;
	if (it <= 2) {
	    goto L201;
	}
	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ped[i__ + (ms2 + ip * ped_dim2) * ped_dim1] = xx[i__];
	}
L200:
	;
    }

/*     reorder if mtrunc is less than nlat-1 */
/*         case of even functions */

    if (modn == 0) {
	nshe[0] = (*nlat - *mtrunc - 1) / 2;
	nshe[1] = (*nlat - *mtrunc - 2) / 2;
    } else {
	nshe[0] = (*nlat - *mtrunc) / 2;
	nshe[1] = (*nlat - *mtrunc - 1) / 2;
    }

    for (mp1 = 1; mp1 <= 2; ++mp1) {
	i__1 = nte;
	for (j = 1; j <= i__1; ++j) {
	    js = j + nshe[mp1 - 1];
	    if (js > nte) {
		js -= nte;
	    }
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		u[i__ + js * u_dim1] = ped[i__ + (j + mp1 * ped_dim2) * 
			ped_dim1];
	    }
	}
	i__1 = nte;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ped[i__ + (j + mp1 * ped_dim2) * ped_dim1] = u[i__ + j * 
			u_dim1];
	    }
	}
/* L210: */
    }

    trunc_(&c__0, &nte, idp, &ped[(ped_dim2 + 1) * ped_dim1 + 1], &nte, &ipse[
	    ipse_dim1 + 1]);
    trunc_(&c__0, &nte, idp, &ped[((ped_dim2 << 1) + 1) * ped_dim1 + 1], &nte,
	     &ipse[(ipse_dim1 << 1) + 1]);

/*     compute the analysis matrices */

    for (ip = 1; ip <= 2; ++ip) {
	i__1 = nte;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lock = 0;
	    i__2 = nte;
	    for (j = 1; j <= i__2; ++j) {
		sum = 0.;
		i__3 = nte;
		for (k = 1; k <= i__3; ++k) {
		    sum += ped[k + (i__ + ip * ped_dim2) * ped_dim1] * we[k + 
			    (j + ip * we_dim2) * we_dim1];
		}
		pe[i__ + (j + ip * pe_dim2) * pe_dim1] = ped[i__ + (j + ip * 
			ped_dim2) * ped_dim1];
		ze[j + (i__ + ip * ze_dim2) * ze_dim1] = sum;
		if (abs(sum) > 5e-8 && lock == 0) {
		    lock = 1;
		    jzse[i__ + ip * jzse_dim1] = j;
		}
	    }
	}
/* L250: */
    }

/*     compute weight matrices for odd functions */

    for (mp1 = 1; mp1 <= 2; ++mp1) {
	m = mp1 - 1;
	mrank = *nlat - m - m;
	nem = (mrank + 1) / 2;
	nom = mrank - nem;
	i__1 = nom;
	for (j = 1; j <= i__1; ++j) {
	    n = j + j + m - 1;
	    dlfkp_(&m, &n, &cp[1]);
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		dlftp_(&m, &n, &thet[i__], &cp[1], &pod[i__ + (j + mp1 * 
			pod_dim2) * pod_dim1]);
	    }
	    if (modn == 1) {
		pod[nte + (j + mp1 * pod_dim2) * pod_dim1] = 0.;
	    }
	}
	dsvdc_(&pod[m + 1 + (mp1 * pod_dim2 + 1) * pod_dim1], idp, &nom, &nom,
		 &s[1], &e[1], &u[u_offset], idp, &v, idp, &work[1], &c__10, &
		info);

	i__1 = nom;
	for (j = 1; j <= i__1; ++j) {
	    s[j] = 1. / (s[j] * s[j]);
	}

/*     compute weight matrix as u  s sup -2 u transpose */

	i__1 = nte;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		wo[i__ + (j + mp1 * wo_dim2) * wo_dim1] = 0.;
	    }
	}
	i__1 = nom;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = nom;
	    for (j = 1; j <= i__2; ++j) {
		sum = 0.;
		i__3 = nom;
		for (k = 1; k <= i__3; ++k) {
		    sum += s[k] * u[i__ + k * u_dim1] * u[j + k * u_dim1];
		}
		wo[i__ + m + (j + m + mp1 * wo_dim2) * wo_dim1] = sum;
	    }
	}
/* L50: */
    }
    wo[((wo_dim2 << 1) + 1) * wo_dim1 + 1] = 1.;
    if (modn == 1) {
	wo[nte + (nte + wo_dim2) * wo_dim1] = 1.;
	wo[nte + (nte + (wo_dim2 << 1)) * wo_dim1] = 1.;
    }

/*     compute n**2 basis (odd functions) */

    ip = 2;
    i__1 = mxtr + 1;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	ip = 3 - ip;
	m = mp1 - 1;
	ms2 = mp1 / 2;
	nrank = ms2 + ms2;
	mrank = *nlat - nrank;
	nem = (mrank + 1) / 2;
	nom = mrank - nem;

/*     compute associated legendre functions */

	if (m <= 1) {
	    i__2 = nom;
	    for (j = 1; j <= i__2; ++j) {
		n = j + j + m - 1;
		dlfkp_(&m, &n, &cp[1]);
		i__3 = nte;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    dlftp_(&m, &n, &thet[i__], &cp[1], &pod[i__ + (j + ms2 + 
			    ip * pod_dim2) * pod_dim1]);
		}
/* L302: */
		if (modn == 1) {
		    pod[nte + (j + ms2 + ip * pod_dim2) * pod_dim1] = 0.;
		}
		if (m > 0) {
		    pod[(j + ms2 + ip * pod_dim2) * pod_dim1 + 1] = 0.;
		}
/* L305: */
	    }

	} else {

	    i__2 = nom;
	    for (j = 1; j <= i__2; ++j) {
		n = j + j + m - 1;
		if (m > 1 && n > mxtr) {
		    i__3 = nte;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			u[i__ + (j + ms2) * u_dim1] = pod[i__ + (j + ms2 + ip 
				* pod_dim2) * pod_dim1];
		    }
		    goto L304;
		}
		a1 = b[n - 1] * a[n + m - 3] / a[n + m - 1];
		b1 = a[n - m + 1] / a[n + m - 1];
		if (n - m <= 1) {
		    i__3 = nte;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			u[i__ + (j + ms2) * u_dim1] = a1 * pod[i__ + (j + ms2 
				- 1 + ip * pod_dim2) * pod_dim1] - b1 * pod[
				i__ + (j + ms2 + ip * pod_dim2) * pod_dim1];
		    }
		} else {
		    c1 = b[n - 1] * a[n - m - 1] / a[n + m - 1];
		    i__3 = nte;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			u[i__ + (j + ms2) * u_dim1] = a1 * pod[i__ + (j + ms2 
				- 1 + ip * pod_dim2) * pod_dim1] - b1 * pod[
				i__ + (j + ms2 + ip * pod_dim2) * pod_dim1] + 
				c1 * u[i__ + (j + ms2 - 1) * u_dim1];
		    }
		}
L304:
		if (modn == 1) {
		    u[nte + (j + ms2) * u_dim1] = 0.;
		}
/* L307: */
	    }
	    i__2 = nom;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = nte;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    pod[i__ + (j + ms2 + ip * pod_dim2) * pod_dim1] = u[i__ + 
			    (j + ms2) * u_dim1];
		}
	    }
	}

	if (ms2 <= 0 || ms2 >= nto) {
	    goto L300;
	}
	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xx[i__] = rand_(&c__0);
	}
	if (modn == 1) {
	    xx[nte] = 0.;
	}
	it = 0;
L306:
	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] = 0.;
	    wx[i__] = 0.;
	    i__3 = nto;
	    for (j = 1; j <= i__3; ++j) {
		wx[i__] += wo[i__ + (j + ip * wo_dim2) * wo_dim1] * xx[j];
	    }
	}
	i__2 = nto;
	for (j = 1; j <= i__2; ++j) {
	    if (j == ms2) {
		goto L330;
	    }
	    gs_(&nte, &wx[1], &pod[(j + ip * pod_dim2) * pod_dim1 + 1], &z__[
		    1]);
L330:
	    ;
	}

	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xx[i__] -= z__[i__];
	}
	normal_(&nte, &xx[1], idp, &wo[(ip * wo_dim2 + 1) * wo_dim1 + 1]);
	++it;
	if (it <= 2) {
	    goto L306;
	}
	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    pod[i__ + (ms2 + ip * pod_dim2) * pod_dim1] = xx[i__];
	}
	if (modn == 1) {
	    pod[nte + (ms2 + ip * pod_dim2) * pod_dim1] = 0.;
	}
L300:
	;
    }

/*     reorder if mtrunc is less than nlat-1 */
/*        case of odd functions */

    if (modn == 0) {
	nsho[0] = (*nlat - *mtrunc) / 2;
	nsho[1] = (*nlat - *mtrunc - 1) / 2;
    } else {
	nsho[0] = (*nlat - *mtrunc - 1) / 2;
	nsho[1] = (*nlat - *mtrunc - 2) / 2;
    }

    for (mp1 = 1; mp1 <= 2; ++mp1) {
	i__1 = nto;
	for (j = 1; j <= i__1; ++j) {
	    js = j + nsho[mp1 - 1];
	    if (js > nto) {
		js -= nto;
	    }
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		u[i__ + js * u_dim1] = pod[i__ + (j + mp1 * pod_dim2) * 
			pod_dim1];
	    }
	}
	i__1 = nto;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		pod[i__ + (j + mp1 * pod_dim2) * pod_dim1] = u[i__ + j * 
			u_dim1];
	    }
	}
/* L310: */
    }

    trunc_(&c__0, &nte, idp, &pod[(pod_dim2 + 1) * pod_dim1 + 1], &nto, &ipso[
	    ipso_dim1 + 1]);
    trunc_(&c__0, &nte, idp, &pod[((pod_dim2 << 1) + 1) * pod_dim1 + 1], &nto,
	     &ipso[(ipso_dim1 << 1) + 1]);

/*     compute the analysis matrices (odd functions) */

    for (ip = 1; ip <= 2; ++ip) {
	i__1 = nto;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lock = 0;
	    i__2 = nto;
	    for (j = 1; j <= i__2; ++j) {
		sum = 0.;
		i__3 = nte;
		for (k = 1; k <= i__3; ++k) {
		    sum += pod[k + (i__ + ip * pod_dim2) * pod_dim1] * wo[k + 
			    (j + ip * wo_dim2) * wo_dim1];
		}
		po[i__ + (j + ip * po_dim2) * po_dim1] = pod[i__ + (j + ip * 
			pod_dim2) * pod_dim1];
		zo[j + (i__ + ip * zo_dim2) * zo_dim1] = sum;
		if (abs(sum) > 5e-8 && lock == 0) {
		    lock = 1;
		    jzso[i__ + ip * jzso_dim1] = j;
		}
	    }
	}
    }
    return 0;
} /* shpei1_ */



/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
/*  .                                                             . */
/*  .                  copyright (c) 1998 by UCAR                 . */
/*  .                                                             . */
/*  .       University Corporation for Atmospheric Research       . */
/*  .                                                             . */
/*  .                      all rights reserved                    . */
/*  .                                                             . */
/*  .                                                             . */
/*  .                         SPHEREPACK3.0                       . */
/*  .                                                             . */
/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* ... file shpe.f */

/* ... files which must be loaded with shpe.f */

/*     sphcom.f, hrfft.f */

/*     the n**2 projection with complement, odd/even */
/*     factorization and zero truncation on an */
/*     equally spaced grid as defined in the JCP paper */
/*     "Generalized discrete spherical harmonic transforms" */
/*     by Paul N. Swarztrauber and William F. Spotz */
/*     It is equivalent to a harmonic analysis followed */
/*     by a synthesis except faster and requires less memory. */

/*     subroutine shpe(nlat,nlon,isym,mtrunc,x,y,idxy, */
/*    1        wshp,lwshp,iwshp,liwshp,work,lwork,ierror) */

/*     shpe projects the array x onto the set of functions represented */
/*     by a discrete set of spherical harmonics. */

/*     input parameters */

/*     nlat   the number of colatitudes on the full sphere including the */
/*            poles. for example, nlat = 37 for a five degree grid. */
/*            nlat determines the grid increment in colatitude as */
/*            pi/(nlat-1).  if nlat is odd the equator is located at */
/*            grid point i=(nlat+1)/2. if nlat is even the equator is */
/*            located half way between points i=nlat/2 and i=nlat/2+1. */
/*            nlat must be at least 3. */

/*     nlon   the number of distinct londitude points.  nlon determines */
/*            the grid increment in longitude as 2*pi/nlon. for example */
/*            nlon = 72 for a five degree grid. nlon must be greater */
/*            than or equal to 4. the efficiency of the computation is */
/*            improved when nlon is a product of small prime numbers. */
/*            nlon must beat least 4. */

/*     isym   currently not used. */

/*     mtrunc the highest longitudinal wave number retained in the */
/*            projection. It must be less than or equal to */
/*            the minimum of nlat-1 and nlon/2. The first wave */
/*            number is zero. For example, if wave numbers 0 and */
/*            1 are desired then mtrunc = 1. */
/*            zero. */

/*     x      a two dimensional array that contains the the nlat */
/*            by nlon array x(i,j) defined at the colatitude point */
/*            theta(i) = (i-1)*pi/(nlat-1) and longitude point phi(j) = */
/*            (j-1)*2*pi/nlon. */

/*     idxy   the first dimension of the arrays x and y as they */
/*            appear in the program that calls shpe. It must be */
/*            at least nlat. */

/*     wshp   a single precision array that must be saved for */
/*            repeated use by subroutine shpe. */

/*     lwshp  the dimension of the array wshp as it appears in the */
/*            program that calls shpei. It must be at least */
/*            2*(nlat+1)**2+nlon+log2(nlon) */

/*     iwshp  an integer array that must be saved for repeated */
/*            use by subroutine shpe. */


/*     liwshp the dimension of the array iwshp as it appears in the */
/*            program that calls shpei. It must be at least */
/*            4*(nlat+1). */

/*     work   a single precision work array that does */
/*            not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls shpe. It must be at least */
/*            max(nlat*nlon,4*(nlat+1)). */

/*     ************************************************************** */

/*     output parameters */

/*     y      an nlat by nlon single precision array that contains */
/*            the projection of x onto the set of functions that */
/*            can be represented by the discrete set of spherical */
/*            harmonics. The arrays x(i,j) and y(i,j) are located */
/*            at colatitude point theta(i) = (i-1)*pi/(nlat-1) and */
/*            longitude point phi(j) = (j-1)*2*pi/nlon. */

/*     ierror = 0  no errors */
/*            = 1  error in the specification of nlat */
/*            = 2  error in the specification of nlon */
/*            = 3  error in the specification of isym */
/*            = 4  error in the specification of mtrunc */
/*            = 5  error in the specification of lwshp */
/*            = 6  error in the specification of liwshp */
/*            = 7  error in the specification of lwork */

/* Subroutine */ int shpe_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, doublereal *x, doublereal *y, integer *idxy, doublereal *wshp, integer 
	*lwshp, integer *iwshp, integer *liwshp, doublereal *work, integer *lwork, 
	integer *ierror)
{
    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, i__1, i__2;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal sn;
    static integer iw1, iw2, iw3, lw1, iw4, jw1, jw2, jw3, jw4, nte, mmax, 
	    mwrk, nloc1, nloc2;
    extern /* Subroutine */ int shpe1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *);
    static integer log2n;
    extern /* Subroutine */ int hrfftb_(integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *), hrfftf_(integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *);



    /* Parameter adjustments */
    y_dim1 = *idxy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    x_dim1 = *idxy;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --wshp;
    --iwshp;
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
/*      ierror = 3 */
/*      if(isym.lt.0 .or. isym.gt.2) return */
    *ierror = 4;
/* Computing MIN */
    i__1 = *nlat - 1, i__2 = *nlon / 2;
    mmax = min(i__1,i__2);
    if (*mtrunc < 0 || *mtrunc > mmax) {
	return 0;
    }
    *ierror = 5;
    log2n = log((doublereal) (*nlon)) / log(2.);
/* Computing 2nd power */
    i__1 = *nlat + 1;
    lw1 = i__1 * i__1 << 1;
    if (*lwshp < lw1 + *nlon + log2n) {
	return 0;
    }
    *ierror = 6;
    if (*liwshp < *nlat + 1 << 2) {
	return 0;
    }
    *ierror = 7;
/* Computing MAX */
    i__1 = *nlat * *nlon, i__2 = *nlat + 1 << 2;
    mwrk = max(i__1,i__2);
    if (*lwork < mwrk) {
	return 0;
    }
    *ierror = 0;

    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ + j * y_dim1] = x[i__ + j * x_dim1];
	}
    }
    hrfftf_(nlat, nlon, &y[y_offset], idxy, &wshp[lw1 + 1], &work[1]);

    nte = (*nlat + 1) / 2;
    nloc1 = (nte << 1) * nte;
    nloc2 = *nlat + 1;
    iw1 = 1;
    iw2 = iw1 + nloc1;
    iw3 = iw2 + nloc1;
    iw4 = iw3 + nloc1;
    jw1 = 1;
    jw2 = jw1 + nloc2;
    jw3 = jw2 + nloc2;
    jw4 = jw3 + nloc2;

    shpe1_(nlat, nlon, isym, mtrunc, &y[y_offset], &y[y_offset], idxy, ierror,
	     &nte, &wshp[iw1], &wshp[iw2], &wshp[iw3], &wshp[iw4], &iwshp[jw1]
	    , &iwshp[jw2], &iwshp[jw3], &iwshp[jw4], &work[jw1], &work[jw2], &
	    work[jw3], &work[jw4]);

    hrfftb_(nlat, nlon, &y[y_offset], idxy, &wshp[lw1 + 1], &work[1]);

    sn = 1. / *nlon;
    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ + j * y_dim1] = sn * y[i__ + j * y_dim1];
	}
    }
    return 0;
} /* shpe_ */

/* Subroutine */ int shpe1_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, doublereal *sx, doublereal *sy, integer *idxy, integer *ierror, 
	integer *idp, doublereal *pe, doublereal *po, doublereal *ze, doublereal *zo, integer *ipse, 
	integer *jzse, integer *ipso, integer *jzso, doublereal *xe, doublereal *xo, doublereal 
	*ye, doublereal *yo)
{
    /* System generated locals */
    integer sx_dim1, sx_offset, sy_dim1, sy_offset, pe_dim1, pe_dim2, 
	    pe_offset, po_dim1, po_dim2, po_offset, ze_dim1, ze_dim2, 
	    ze_offset, zo_dim1, zo_dim2, zo_offset, ipse_dim1, ipse_offset, 
	    jzse_dim1, jzse_offset, ipso_dim1, ipso_offset, jzso_dim1, 
	    jzso_offset, xe_dim1, xe_offset, xo_dim1, xo_offset, ye_dim1, 
	    ye_offset, yo_dim1, yo_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, m, ip, js, mp1, ms2, ns2, nec, nem, noc, nte, mpm, 
	    nom, nto, nshe[2], modn, nsho[2], mxtr, mrank, nrank;
    extern /* Subroutine */ int tmxmx_(integer *, integer *, integer *, doublereal *
	    , integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *);



    /* Parameter adjustments */
    sy_dim1 = *idxy;
    sy_offset = 1 + sy_dim1;
    sy -= sy_offset;
    sx_dim1 = *idxy;
    sx_offset = 1 + sx_dim1;
    sx -= sx_offset;
    yo_dim1 = *idp;
    yo_offset = 1 + yo_dim1;
    yo -= yo_offset;
    ye_dim1 = *idp;
    ye_offset = 1 + ye_dim1;
    ye -= ye_offset;
    xo_dim1 = *idp;
    xo_offset = 1 + xo_dim1;
    xo -= xo_offset;
    xe_dim1 = *idp;
    xe_offset = 1 + xe_dim1;
    xe -= xe_offset;
    jzso_dim1 = *idp;
    jzso_offset = 1 + jzso_dim1;
    jzso -= jzso_offset;
    ipso_dim1 = *idp;
    ipso_offset = 1 + ipso_dim1;
    ipso -= ipso_offset;
    jzse_dim1 = *idp;
    jzse_offset = 1 + jzse_dim1;
    jzse -= jzse_offset;
    ipse_dim1 = *idp;
    ipse_offset = 1 + ipse_dim1;
    ipse -= ipse_offset;
    zo_dim1 = *idp;
    zo_dim2 = *idp;
    zo_offset = 1 + zo_dim1 * (1 + zo_dim2);
    zo -= zo_offset;
    ze_dim1 = *idp;
    ze_dim2 = *idp;
    ze_offset = 1 + ze_dim1 * (1 + ze_dim2);
    ze -= ze_offset;
    po_dim1 = *idp;
    po_dim2 = *idp;
    po_offset = 1 + po_dim1 * (1 + po_dim2);
    po -= po_offset;
    pe_dim1 = *idp;
    pe_dim2 = *idp;
    pe_offset = 1 + pe_dim1 * (1 + pe_dim2);
    pe -= pe_offset;

    /* Function Body */
    ns2 = *nlat / 2;
    modn = *nlat - ns2 - ns2;
    nte = (*nlat + 1) / 2;
    nto = *nlat - nte;

    if (modn == 0) {
	nshe[0] = (*nlat - *mtrunc - 1) / 2;
	nshe[1] = (*nlat - *mtrunc - 2) / 2;
	nsho[0] = (*nlat - *mtrunc) / 2;
	nsho[1] = (*nlat - *mtrunc - 1) / 2;
    } else {
	nshe[0] = (*nlat - *mtrunc) / 2;
	nshe[1] = (*nlat - *mtrunc - 1) / 2;
	nsho[0] = (*nlat - *mtrunc - 1) / 2;
	nsho[1] = (*nlat - *mtrunc - 2) / 2;
    }
/* Computing MIN */
    i__1 = *nlat - 1, i__2 = *nlon / 2, i__1 = min(i__1,i__2);
    mxtr = min(i__1,*mtrunc);
    ip = 2;
    i__1 = mxtr + 1;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	ip = 3 - ip;
	if (mxtr == *nlat - 1 && mp1 <= 2) {
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sy[i__ + mp1 * sy_dim1] = sx[i__ + mp1 * sx_dim1];
	    }
	    if (mp1 == 2) {
		sy[(sy_dim1 << 1) + 1] = 0.;
		sy[*nlat + (sy_dim1 << 1)] = 0.;
	    }
	    if (*nlon >= 3) {
		sy[sy_dim1 * 3 + 1] = 0.;
		sy[*nlat + sy_dim1 * 3] = 0.;
		i__2 = *nlat - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    sy[i__ + sy_dim1 * 3] = sx[i__ + sx_dim1 * 3];
		}
	    }
	    goto L100;
	}
	m = mp1 - 1;
/* Computing MAX */
	i__2 = 1, i__3 = m + m;
	mpm = max(i__2,i__3);
	ms2 = mp1 / 2;
/* Computing MIN */
	i__2 = *nlat - m, i__3 = *nlat - ms2 - ms2;
	mrank = min(i__2,i__3);
/*      mrank = mxtr+1-ms2-ms2 */
	nrank = *nlat - mrank;
	nem = (mrank + 1) / 2 - nshe[ip - 1];
	nom = mrank - (mrank + 1) / 2 - nsho[ip - 1];
	nec = nte - nem;
	noc = nto - nom;

	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xe[i__ + xe_dim1] = (sx[i__ + mpm * sx_dim1] + sx[*nlat + 1 - i__ 
		    + mpm * sx_dim1]) * .5;
	    xo[i__ + xo_dim1] = (sx[i__ + mpm * sx_dim1] - sx[*nlat + 1 - i__ 
		    + mpm * sx_dim1]) * .5;
	}
	if (mpm < *nlon) {
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		xe[i__ + (xe_dim1 << 1)] = (sx[i__ + (mpm + 1) * sx_dim1] + 
			sx[*nlat + 1 - i__ + (mpm + 1) * sx_dim1]) * .5;
		xo[i__ + (xo_dim1 << 1)] = (sx[i__ + (mpm + 1) * sx_dim1] - 
			sx[*nlat + 1 - i__ + (mpm + 1) * sx_dim1]) * .5;
	    }
	}
	if (nec * 3 < nem << 1 || nem == 0) {
	    tmxmx_(&nte, &nec, idp, &pe[(ip * pe_dim2 + 1) * pe_dim1 + 1], &
		    nte, idp, &ze[(ip * ze_dim2 + 1) * ze_dim1 + 1], &xe[
		    xe_offset], &ye[ye_offset], &ipse[ip * ipse_dim1 + 1], &
		    jzse[ip * jzse_dim1 + 1]);
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ye[i__ + ye_dim1] = xe[i__ + xe_dim1] - ye[i__ + ye_dim1];
	    }
	    if (mpm < *nlon && m != 0) {
		i__2 = nte;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ye[i__ + (ye_dim1 << 1)] = xe[i__ + (xe_dim1 << 1)] - ye[
			    i__ + (ye_dim1 << 1)];
		}
	    }
	} else {
	    tmxmx_(&nte, &nem, idp, &pe[(nec + 1 + ip * pe_dim2) * pe_dim1 + 
		    1], &nte, idp, &ze[(nec + 1 + ip * ze_dim2) * ze_dim1 + 1]
		    , &xe[xe_offset], &ye[ye_offset], &ipse[nec + 1 + ip * 
		    ipse_dim1], &jzse[nec + 1 + ip * jzse_dim1]);
	}
	if (noc * 3 < nom << 1 || nom == 0) {
	    tmxmx_(&nto, &noc, idp, &po[(ip * po_dim2 + 1) * po_dim1 + 1], &
		    nto, idp, &zo[(ip * zo_dim2 + 1) * zo_dim1 + 1], &xo[
		    xo_offset], &yo[yo_offset], &ipso[ip * ipso_dim1 + 1], &
		    jzso[ip * jzso_dim1 + 1]);
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		yo[i__ + yo_dim1] = xo[i__ + xo_dim1] - yo[i__ + yo_dim1];
	    }
	    if (mpm < *nlon && m != 0) {
		i__2 = nte;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    yo[i__ + (yo_dim1 << 1)] = xo[i__ + (xo_dim1 << 1)] - yo[
			    i__ + (yo_dim1 << 1)];
		}
	    }
	} else {
	    tmxmx_(&nto, &nom, idp, &po[(noc + 1 + ip * po_dim2) * po_dim1 + 
		    1], &nto, idp, &zo[(noc + 1 + ip * zo_dim2) * zo_dim1 + 1]
		    , &xo[xo_offset], &yo[yo_offset], &ipso[noc + 1 + ip * 
		    ipso_dim1], &jzso[noc + 1 + ip * jzso_dim1]);
	}
	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sy[i__ + mpm * sy_dim1] = ye[i__ + ye_dim1] + yo[i__ + yo_dim1];
	    sy[*nlat + 1 - i__ + mpm * sy_dim1] = ye[i__ + ye_dim1] - yo[i__ 
		    + yo_dim1];
	}
	if (mpm < *nlon && m != 0) {
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sy[i__ + (mpm + 1) * sy_dim1] = ye[i__ + (ye_dim1 << 1)] + yo[
			i__ + (yo_dim1 << 1)];
		sy[*nlat + 1 - i__ + (mpm + 1) * sy_dim1] = ye[i__ + (ye_dim1 
			<< 1)] - yo[i__ + (yo_dim1 << 1)];
	    }
	}
L100:
	;
    }
    js = mxtr + mxtr + 2;
    i__1 = *nlon;
    for (j = js; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sy[i__ + j * sy_dim1] = 0.;
	}
    }
    return 0;
} /* shpe1_ */

/* Subroutine */ static int mxm_(integer *lr, integer *lc, integer *ld, doublereal *
	a, integer *mc, integer *md, doublereal *b, integer *nd, doublereal *
	c__)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, k;

    /* Parameter adjustments */
    a_dim1 = *ld;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *md;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *nd;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    i__1 = *lr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *mc;
	for (j = 1; j <= i__2; ++j) {
	    c__[i__ + j * c_dim1] = 0.;
	    i__3 = *lc;
	    for (k = 1; k <= i__3; ++k) {
		c__[i__ + j * c_dim1] += a[i__ + k * a_dim1] * b[k + j * 
			b_dim1];
	    }
	}
    }
    return 0;
} /* mxm_ */

/* Subroutine */ int smxm_(integer *lr, integer *lc, integer *ld, doublereal *a, 
	integer *mc, integer *md, doublereal *b, integer *nd, doublereal *c__)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer i__, j, k;

    /* Parameter adjustments */
    a_dim1 = *ld;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *md;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    c_dim1 = *nd;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    i__1 = *lr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *mc;
	for (j = 1; j <= i__2; ++j) {
	    c__[i__ + j * c_dim1] = 0.;
	    i__3 = *lc;
	    for (k = 1; k <= i__3; ++k) {
		c__[i__ + j * c_dim1] += a[i__ + k * a_dim1] * b[k + j * 
			b_dim1];
	    }
	}
    }
    return 0;
} /* smxm_ */

/* Subroutine */ int mxmx_(integer *lr, integer *lc, integer *ld, doublereal *a, 
	integer *mc, integer *md, doublereal *b, doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset, y_dim1, 
	    y_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal sum1, sum2;

    /* Parameter adjustments */
    y_dim1 = *ld;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    x_dim1 = *ld;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    a_dim1 = *ld;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *md;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    i__1 = *lr;
    for (k = 1; k <= i__1; ++k) {
	y[k + y_dim1] = 0.;
	y[k + (y_dim1 << 1)] = 0.;
    }

    if (*lc <= 0) {
	return 0;
    }
    i__1 = *lc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum1 = 0.;
	sum2 = 0.;
	i__2 = *mc;
	for (j = 1; j <= i__2; ++j) {
	    sum1 += b[i__ + j * b_dim1] * x[j + x_dim1];
	    sum2 += b[i__ + j * b_dim1] * x[j + (x_dim1 << 1)];
	}
	i__2 = *lr;
	for (k = 1; k <= i__2; ++k) {
	    y[k + y_dim1] += sum1 * a[k + i__ * a_dim1];
	    y[k + (y_dim1 << 1)] += sum2 * a[k + i__ * a_dim1];
	}
    }
    return 0;
} /* mxmx_ */

/* Subroutine */ int dmxmx_(integer *lr, integer *lc, integer *ld, doublereal 
	*a, integer *mc, integer *md, doublereal *b, doublereal *x, 
	doublereal *y)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset, y_dim1, 
	    y_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal sum1, sum2;

    /* Parameter adjustments */
    y_dim1 = *ld;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    x_dim1 = *ld;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    a_dim1 = *ld;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *md;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    i__1 = *lr;
    for (k = 1; k <= i__1; ++k) {
	y[k + y_dim1] = 0.;
	y[k + (y_dim1 << 1)] = 0.;
    }

    if (*lc <= 0) {
	return 0;
    }
    i__1 = *lc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum1 = 0.;
	sum2 = 0.;
	i__2 = *mc;
	for (j = 1; j <= i__2; ++j) {
	    sum1 += b[i__ + j * b_dim1] * x[j + x_dim1];
	    sum2 += b[i__ + j * b_dim1] * x[j + (x_dim1 << 1)];
	}
	i__2 = *lr;
	for (k = 1; k <= i__2; ++k) {
	    y[k + y_dim1] += sum1 * a[k + i__ * a_dim1];
	    y[k + (y_dim1 << 1)] += sum2 * a[k + i__ * a_dim1];
	}
    }
    return 0;
} /* dmxmx_ */

/* Subroutine */ int tmxmx_(integer *lr, integer *lc, integer *ld, doublereal *a, 
	integer *mc, integer *md, doublereal *b, doublereal *x, doublereal *y, integer *is, 
	integer *js)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset, y_dim1, 
	    y_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, kmx;
    static doublereal sum1, sum2;


    /* Parameter adjustments */
    y_dim1 = *ld;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    x_dim1 = *ld;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    a_dim1 = *ld;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    b_dim1 = *md;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --is;
    --js;

    /* Function Body */
/* Computing MIN */
    i__1 = *lr + 1;
    kmx = min(i__1,*ld);
    i__1 = kmx;
    for (k = 1; k <= i__1; ++k) {
	y[k + y_dim1] = 0.;
	y[k + (y_dim1 << 1)] = 0.;
    }
    if (*lc <= 0) {
	return 0;
    }

    i__1 = *lc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum1 = 0.;
	sum2 = 0.;
	i__2 = *mc;
	for (j = js[i__]; j <= i__2; ++j) {
	    sum1 += b[j + i__ * b_dim1] * x[j + x_dim1];
	    sum2 += b[j + i__ * b_dim1] * x[j + (x_dim1 << 1)];
	}
	i__2 = *lr;
	for (k = is[i__]; k <= i__2; ++k) {
	    y[k + y_dim1] += sum1 * a[k + i__ * a_dim1];
	    y[k + (y_dim1 << 1)] += sum2 * a[k + i__ * a_dim1];
	}
    }
    return 0;
} /* tmxmx_ */

/* Subroutine */ int trunc_(integer *irc, integer *n, integer *idp, 
	doublereal *a, integer *nrc, integer *ijs)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;


/*     irc = 0 for columns , or irc = 1 for rows */

    /* Parameter adjustments */
    --ijs;
    a_dim1 = *idp;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    if (*irc != 0) {
	goto L30;
    }
    i__1 = *nrc;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ijs[j] = i__;
	    if ((d__1 = a[i__ + j * a_dim1], abs(d__1)) > 5e-8) {
		goto L20;
	    }
	}
L20:
	;
    }
    return 0;
L30:
    i__1 = *nrc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    ijs[i__] = j;
	    if ((d__1 = a[i__ + j * a_dim1], abs(d__1)) > 5e-8) {
		goto L50;
	    }
	}
L50:
	;
    }
    return 0;
} /* trunc_ */

/* Subroutine */ int gs_(integer *n, doublereal *x, doublereal *y, doublereal 
	*z__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal sum;


/*     accumulate innerproducts of x with respect to y. */

    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    sum = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum += x[i__] * y[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__] += sum * y[i__];
    }
    return 0;
} /* gs_ */

/* Subroutine */ int normal_(integer *n, doublereal *x, integer *id, 
	doublereal *q)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal sum, sqs;


/*     normalize x */

    /* Parameter adjustments */
    --x;
    q_dim1 = *id;
    q_offset = 1 + q_dim1;
    q -= q_offset;

    /* Function Body */
    sqs = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    sum += q[i__ + j * q_dim1] * x[j];
	}
	sqs += sum * x[i__];
    }

    sqs = sqrt(sqs);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] /= sqs;
    }
    return 0;
} /* normal_ */

/* Subroutine */ int coe_(integer *moe, integer *n, doublereal *x, doublereal 
	*dmax__)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, nh;

    /* Parameter adjustments */
    --x;

    /* Function Body */
    nh = (*n + 1) / 2;
    *dmax__ = 0.;
    if (*moe != 0) {
	goto L1;
    }
    i__1 = nh;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = *dmax__, d__3 = (d__1 = x[i__] - x[*n - i__ + 1], abs(d__1));
	*dmax__ = max(d__2,d__3);
	x[i__] = (x[i__] + x[*n - i__ + 1]) * .5;
	x[*n - i__ + 1] = x[i__];
    }
    return 0;
L1:
    i__1 = nh;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = *dmax__, d__3 = (d__1 = x[i__] + x[*n - i__ + 1], abs(d__1));
	*dmax__ = max(d__2,d__3);
	x[i__] = (x[i__] - x[*n - i__ + 1]) * .5;
	x[*n - i__ + 1] = -x[i__];
    }
    if (*n % 2 != 0) {
	x[nh] = 0.;
    }
    return 0;
} /* coe_ */

/*     subroutine dlfkp(m,n,cp) */

/*     subroutine dlfkp computes the coefficients in the trigonometric */
/*     expansion of the normalized associated legendre functions: */

/*     pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m))) */
/*                        *sin(theta)**m/(2**n*factorial(n)) times the */
/*                        (n+m)th derivative of (x**2-1)**n with respect */
/*                        to x=cos(theta) */

/*     where theta is colatitude. */

/*     subroutine dlfkp computes the coefficients cp(k) in the */
/*     following trigonometric expansion of pbar(m,n,theta). */

/*            1) for n even and m even, pbar(m,n,theta) = */
/*               .5*cp(1) plus the sum from k=1 to k=n/2 */
/*               of cp(k)*cos(2*k*th) */

/*            2) for n even and m odd, pbar(m,n,theta) = */
/*               the sum from k=1 to k=n/2 of */
/*               cp(k)*sin(2*k*th) */

/*            3) for n odd and m even, pbar(m,n,theta) = */
/*               the sum from k=1 to k=(n+1)/2 of */
/*               cp(k)*cos((2*k-1)*th) */

/*            4) for n odd and m odd,  pbar(m,n,theta) = */
/*               the sum from k=1 to k=(n+1)/2 of */
/*               cp(k)*sin((2*k-1)*th) */

/*     input parameters */

/*     m      is the order of pbar(n,m,theta). m can be any integer */
/*            however pbar(n,m,theta) = 0  if abs(m) is greater than */
/*            n and pbar(n,m,theta) = (-1)**m*pbar(n,-m,theta) for */
/*            negative m. */

/*     n      nonnegative integer specifying the degree of */
/*            pbar(n,m,theta) */

/*     output parameters */

/*     cp     a double precision array that contains the fourier */
/*            coefficients for pbar(m,n,theta). the length of the */
/*            array depends on the parity of m and n */

/*                  parity            length of cp */

/*               n even m even           n/2+1 */
/*               n even m odd             n/2 */
/*               n odd  m even          (n+1)/2 */
/*               n odd  m odd           (n+1)/2 */


/* **************************************************************** */
/* Subroutine */ int dlfkp_(integer *m, integer *n, doublereal *cp)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal), pow_di(doublereal *, integer *);

    /* Local variables */
    static integer i__, l;
    static doublereal a1, b1, c1, t1, t2;
    static integer ma;
    static doublereal fk, cp2, pm1;
    static integer nex;
    static doublereal fden, fnmh, fnum, fnnp1;
    static integer nmms2;
    static doublereal fnmsq;



    /* Parameter adjustments */
    --cp;

    /* Function Body */
    cp[1] = 0.;
    ma = abs(*m);
    if (ma > *n) {
	return 0;
    }
    if ((i__1 = *n - 1) < 0) {
	goto L2;
    } else if (i__1 == 0) {
	goto L3;
    } else {
	goto L5;
    }
L2:
    cp[1] = sqrt(2.);
    return 0;
L3:
    if (ma != 0) {
	goto L4;
    }
    cp[1] = sqrt(1.5);
    return 0;
L4:
    cp[1] = sqrt(.75);
    if (*m == -1) {
	cp[1] = -cp[1];
    }
    return 0;
L5:
    if ((*n + ma) % 2 != 0) {
	goto L10;
    }
    nmms2 = (*n - ma) / 2;
    fnum = (doublereal) (*n + ma + 1);
    fnmh = (doublereal) (*n - ma + 1);
    pm1 = 1.;
    goto L15;
L10:
    nmms2 = (*n - ma - 1) / 2;
    fnum = (doublereal) (*n + ma + 2);
    fnmh = (doublereal) (*n - ma + 2);
    pm1 = -1.;
/*      t1 = 1. */
/*      t1 = 2.d0**(n-1) */
/*      t1 = 1.d0/t1 */
L15:
    t1 = 9.5367431640625e-7;
    nex = 20;
    fden = 2.;
    if (nmms2 < 1) {
	goto L20;
    }
    i__1 = nmms2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t1 = fnum * t1 / fden;
	if (t1 > 1048576.) {
	    t1 /= 1099511627776.;
	    nex += 40;
	}
	fnum += 2.;
	fden += 2.;
/* L18: */
    }
L20:
    i__1 = *n - 1 - nex;
    t1 /= pow_di(&c_b54, &i__1);
    if (ma / 2 % 2 != 0) {
	t1 = -t1;
    }
    t2 = 1.;
    if (ma == 0) {
	goto L26;
    }
    i__1 = ma;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t2 = fnmh * t2 / (fnmh + pm1);
	fnmh += 2.;
/* L25: */
    }
L26:
    cp2 = t1 * sqrt((*n + .5) * t2);
    fnnp1 = (doublereal) (*n * (*n + 1));
    fnmsq = fnnp1 - ma * 2. * ma;
    l = (*n + 1) / 2;
    if (*n % 2 == 0 && ma % 2 == 0) {
	++l;
    }
    cp[l] = cp2;
    if (*m >= 0) {
	goto L29;
    }
    if (ma % 2 != 0) {
	cp[l] = -cp[l];
    }
L29:
    if (l <= 1) {
	return 0;
    }
    fk = (doublereal) (*n);
    a1 = (fk - 2.) * (fk - 1.) - fnnp1;
    b1 = (fk * fk - fnmsq) * 2.;
    cp[l - 1] = b1 * cp[l] / a1;
L30:
    --l;
    if (l <= 1) {
	return 0;
    }
    fk += -2.;
    a1 = (fk - 2.) * (fk - 1.) - fnnp1;
    b1 = (fk * fk - fnmsq) * -2.;
    c1 = (fk + 1.) * (fk + 2.) - fnnp1;
    cp[l - 1] = -(b1 * cp[l] + c1 * cp[l + 1]) / a1;
    goto L30;
} /* dlfkp_ */

/* Subroutine */ int dlftp_(integer *m, integer *n, doublereal *theta, 
	doublereal *cp, doublereal *pb)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal chh, cdt;
    static integer kdo;
    static doublereal cth, sdt, sth;
    static integer mmod, nmod;

    /* Parameter adjustments */
    --cp;

    /* Function Body */
    cdt = cos(*theta + *theta);
    sdt = sin(*theta + *theta);
    nmod = *n % 2;
    mmod = abs(*m) % 2;
    if (nmod <= 0) {
	goto L1;
    } else {
	goto L2;
    }
L1:
    if (mmod <= 0) {
	goto L3;
    } else {
	goto L4;
    }

/*     n even, m even */

L3:
    kdo = *n / 2;
    *pb = cp[1] * .5;
    if (*n == 0) {
	return 0;
    }
    cth = cdt;
    sth = sdt;
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k+1)*dcos(2*k*theta) */
	*pb += cp[k + 1] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L170: */
    }
    return 0;

/*     n even, m odd */

L4:
    kdo = *n / 2;
    *pb = 0.;
    cth = cdt;
    sth = sdt;
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k)*dsin(2*k*theta) */
	*pb += cp[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L180: */
    }
    return 0;
L2:
    if (mmod <= 0) {
	goto L13;
    } else {
	goto L14;
    }

/*     n odd, m even */

L13:
    kdo = (*n + 1) / 2;
    *pb = 0.;
    cth = cos(*theta);
    sth = sin(*theta);
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k)*dcos((2*k-1)*theta) */
	*pb += cp[k] * cth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L190: */
    }
    return 0;

/*     n odd, m odd */

L14:
    kdo = (*n + 1) / 2;
    *pb = 0.;
    cth = cos(*theta);
    sth = sin(*theta);
    i__1 = kdo;
    for (k = 1; k <= i__1; ++k) {
/*     pb = pb+cp(k)*dsin((2*k-1)*theta) */
	*pb += cp[k] * sth;
	chh = cdt * cth - sdt * sth;
	sth = sdt * cth + cdt * sth;
	cth = chh;
/* L200: */
    }
    return 0;
} /* dlftp_ */


/* Subroutine */ int dsvdc_(doublereal *x, integer *ldx, integer *n, integer *
	p, doublereal *s, doublereal *e, doublereal *u, integer *ldu, 
	doublereal *v, integer *ldv, doublereal *work, integer *job, integer *
	info)
{
    /* System generated locals */
    integer x_dim1, x_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal b, c__, f, g;
    static integer i__, j, k, l, m;
    static doublereal t, t1, el;
    static integer kk;
    static doublereal cs;
    static integer ll, mm, ls;
    static doublereal sl;
    static integer lu;
    static doublereal sm, sn;
    static integer lm1, mm1, lp1, mp1, nct, ncu, lls, nrt;
    static doublereal emm1, smm1;
    static integer kase;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer jobu, iter;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal test;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer nctp1, nrtp1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale, shift;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), drotg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer maxit;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical wantu, wantv;
    static doublereal ztest;



/*     dsvdc is a subroutine to reduce a double precision nxp matrix x */
/*     by orthogonal transformations u and v to diagonal form.  the */
/*     diagonal elements s(i) are the singular values of x.  the */
/*     columns of u are the corresponding left singular vectors, */
/*     and the columns of v the right singular vectors. */

/*     on entry */

/*         x         double precision(ldx,p), where ldx.ge.n. */
/*                   x contains the matrix whose singular value */
/*                   decomposition is to be computed.  x is */
/*                   destroyed by dsvdc. */

/*         ldx       integer. */
/*                   ldx is the leading dimension of the array x. */

/*         n         integer. */
/*                   n is the number of rows of the matrix x. */

/*         p         integer. */
/*                   p is the number of columns of the matrix x. */

/*         ldu       integer. */
/*                   ldu is the leading dimension of the array u. */
/*                   (see below). */

/*         ldv       integer. */
/*                   ldv is the leading dimension of the array v. */
/*                   (see below). */

/*         work      double precision(n). */
/*                   work is a scratch array. */

/*         job       integer. */
/*                   job controls the computation of the singular */
/*                   vectors.  it has the decimal expansion ab */
/*                   with the following meaning */

/*                        a.eq.0    do not compute the left singular */
/*                                  vectors. */
/*                        a.eq.1    return the n left singular vectors */
/*                                  in u. */
/*                        a.ge.2    return the first min(n,p) singular */
/*                                  vectors in u. */
/*                        b.eq.0    do not compute the right singular */
/*                                  vectors. */
/*                        b.eq.1    return the right singular vectors */
/*                                  in v. */

/*     on return */

/*         s         double precision(mm), where mm=min(n+1,p). */
/*                   the first min(n,p) entries of s contain the */
/*                   singular values of x arranged in descending */
/*                   order of magnitude. */

/*         e         double precision(p), */
/*                   e ordinarily contains zeros.  however see the */
/*                   discussion of info for exceptions. */

/*         u         double precision(ldu,k), where ldu.ge.n.  if */
/*                                   joba.eq.1 then k.eq.n, if joba.ge.2 */
/*                                   then k.eq.min(n,p). */
/*                   u contains the matrix of left singular vectors. */
/*                   u is not referenced if joba.eq.0.  if n.le.p */
/*                   or if joba.eq.2, then u may be identified with x */
/*                   in the subroutine call. */

/*         v         double precision(ldv,p), where ldv.ge.p. */
/*                   v contains the matrix of right singular vectors. */
/*                   v is not referenced if job.eq.0.  if p.le.n, */
/*                   then v may be identified with x in the */
/*                   subroutine call. */

/*         info      integer. */
/*                   the singular values (and their corresponding */
/*                   singular vectors) s(info+1),s(info+2),...,s(m) */
/*                   are correct (here m=min(n,p)).  thus if */
/*                   info.eq.0, all the singular values and their */
/*                   vectors are correct.  in any event, the matrix */
/*                   b = trans(u)*x*v is the bidiagonal matrix */
/*                   with the elements of s on its diagonal and the */
/*                   elements of e on its super-diagonal (trans(u) */
/*                   is the transpose of u).  thus the singular */
/*                   values of x and b are the same. */

/*     linpack. this version dated 08/14/78 . */
/*              correction made to shift 2/84. */
/*     g.w. stewart, university of maryland, argonne national lab. */

/*     dsvdc uses the following functions and subprograms. */

/*     external drot */
/*     blas daxpy,ddot,dscal,dswap,dnrm2,drotg */
/*     fortran dabs,dmax1,max0,min0,mod,dsqrt */

/*     internal variables */

/*      double precision ddot,t,r */


/*     set the maximum number of iterations. */

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --s;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --work;

    /* Function Body */
    maxit = 30;

/*     determine what is to be computed. */

    wantu = FALSE_;
    wantv = FALSE_;
    jobu = *job % 100 / 10;
    ncu = *n;
    if (jobu > 1) {
	ncu = min(*n,*p);
    }
    if (jobu != 0) {
	wantu = TRUE_;
    }
    if (*job % 10 != 0) {
	wantv = TRUE_;
    }

/*     reduce x to bidiagonal form, storing the diagonal elements */
/*     in s and the super-diagonal elements in e. */

    *info = 0;
/* Computing MIN */
    i__1 = *n - 1;
    nct = min(i__1,*p);
/* Computing MAX */
/* Computing MIN */
    i__3 = *p - 2;
    i__1 = 0, i__2 = min(i__3,*n);
    nrt = max(i__1,i__2);
    lu = max(nct,nrt);
    if (lu < 1) {
	goto L170;
    }
    i__1 = lu;
    for (l = 1; l <= i__1; ++l) {
	lp1 = l + 1;
	if (l > nct) {
	    goto L20;
	}

/*           compute the transformation for the l-th column and */
/*           place the l-th diagonal in s(l). */

	i__2 = *n - l + 1;
	s[l] = dnrm2_(&i__2, &x[l + l * x_dim1], &c__1);
	if (s[l] == 0.) {
	    goto L10;
	}
	if (x[l + l * x_dim1] != 0.) {
	    s[l] = d_sign(&s[l], &x[l + l * x_dim1]);
	}
	i__2 = *n - l + 1;
	d__1 = 1. / s[l];
	dscal_(&i__2, &d__1, &x[l + l * x_dim1], &c__1);
	x[l + l * x_dim1] += 1.;
L10:
	s[l] = -s[l];
L20:
	if (*p < lp1) {
	    goto L50;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    if (l > nct) {
		goto L30;
	    }
	    if (s[l] == 0.) {
		goto L30;
	    }

/*              apply the transformation. */

	    i__3 = *n - l + 1;
	    t = -ddot_(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1) / x[l + l * x_dim1];
	    i__3 = *n - l + 1;
	    daxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
L30:

/*           place the l-th row of x into  e for the */
/*           subsequent calculation of the row transformation. */

	    e[j] = x[l + j * x_dim1];
/* L40: */
	}
L50:
	if (! wantu || l > nct) {
	    goto L70;
	}

/*           place the transformation in u for subsequent back */
/*           multiplication. */

	i__2 = *n;
	for (i__ = l; i__ <= i__2; ++i__) {
	    u[i__ + l * u_dim1] = x[i__ + l * x_dim1];
/* L60: */
	}
L70:
	if (l > nrt) {
	    goto L150;
	}

/*           compute the l-th row transformation and place the */
/*           l-th super-diagonal in e(l). */

	i__2 = *p - l;
	e[l] = dnrm2_(&i__2, &e[lp1], &c__1);
	if (e[l] == 0.) {
	    goto L80;
	}
	if (e[lp1] != 0.) {
	    e[l] = d_sign(&e[l], &e[lp1]);
	}
	i__2 = *p - l;
	d__1 = 1. / e[l];
	dscal_(&i__2, &d__1, &e[lp1], &c__1);
	e[lp1] += 1.;
L80:
	e[l] = -e[l];
	if (lp1 > *n || e[l] == 0.) {
	    goto L120;
	}

/*              apply the transformation. */

	i__2 = *n;
	for (i__ = lp1; i__ <= i__2; ++i__) {
	    work[i__] = 0.;
/* L90: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    daxpy_(&i__3, &e[j], &x[lp1 + j * x_dim1], &c__1, &work[lp1], &
		    c__1);
/* L100: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    d__1 = -e[j] / e[lp1];
	    daxpy_(&i__3, &d__1, &work[lp1], &c__1, &x[lp1 + j * x_dim1], &
		    c__1);
/* L110: */
	}
L120:
	if (! wantv) {
	    goto L140;
	}

/*              place the transformation in v for subsequent */
/*              back multiplication. */

	i__2 = *p;
	for (i__ = lp1; i__ <= i__2; ++i__) {
	    v[i__ + l * v_dim1] = e[i__];
/* L130: */
	}
L140:
L150:
/* L160: */
	;
    }
L170:

/*     set up the final bidiagonal matrix or order m. */

/* Computing MIN */
    i__1 = *p, i__2 = *n + 1;
    m = min(i__1,i__2);
    nctp1 = nct + 1;
    nrtp1 = nrt + 1;
    if (nct < *p) {
	s[nctp1] = x[nctp1 + nctp1 * x_dim1];
    }
    if (*n < m) {
	s[m] = 0.;
    }
    if (nrtp1 < m) {
	e[nrtp1] = x[nrtp1 + m * x_dim1];
    }
    e[m] = 0.;

/*     if required, generate u. */

    if (! wantu) {
	goto L300;
    }
    if (ncu < nctp1) {
	goto L200;
    }
    i__1 = ncu;
    for (j = nctp1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    u[i__ + j * u_dim1] = 0.;
/* L180: */
	}
	u[j + j * u_dim1] = 1.;
/* L190: */
    }
L200:
    if (nct < 1) {
	goto L290;
    }
    i__1 = nct;
    for (ll = 1; ll <= i__1; ++ll) {
	l = nct - ll + 1;
	if (s[l] == 0.) {
	    goto L250;
	}
	lp1 = l + 1;
	if (ncu < lp1) {
	    goto L220;
	}
	i__2 = ncu;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    t = -ddot_(&i__3, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &
		    c__1) / u[l + l * u_dim1];
	    i__3 = *n - l + 1;
	    daxpy_(&i__3, &t, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &
		    c__1);
/* L210: */
	}
L220:
	i__2 = *n - l + 1;
	dscal_(&i__2, &c_b113, &u[l + l * u_dim1], &c__1);
	u[l + l * u_dim1] += 1.;
	lm1 = l - 1;
	if (lm1 < 1) {
	    goto L240;
	}
	i__2 = lm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    u[i__ + l * u_dim1] = 0.;
/* L230: */
	}
L240:
	goto L270;
L250:
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    u[i__ + l * u_dim1] = 0.;
/* L260: */
	}
	u[l + l * u_dim1] = 1.;
L270:
/* L280: */
	;
    }
L290:
L300:

/*     if it is required, generate v. */

    if (! wantv) {
	goto L350;
    }
    i__1 = *p;
    for (ll = 1; ll <= i__1; ++ll) {
	l = *p - ll + 1;
	lp1 = l + 1;
	if (l > nrt) {
	    goto L320;
	}
	if (e[l] == 0.) {
	    goto L320;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *p - l;
	    t = -ddot_(&i__3, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1) / v[lp1 + l * v_dim1];
	    i__3 = *p - l;
	    daxpy_(&i__3, &t, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1);
/* L310: */
	}
L320:
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    v[i__ + l * v_dim1] = 0.;
/* L330: */
	}
	v[l + l * v_dim1] = 1.;
/* L340: */
    }
L350:

/*     main iteration loop for the singular values. */

    mm = m;
    iter = 0;
L360:

/*        quit if all the singular values have been found. */

/*     ...exit */
    if (m == 0) {
	goto L620;
    }

/*        if too many iterations have been performed, set */
/*        flag and return. */

    if (iter < maxit) {
	goto L370;
    }
    *info = m;
/*     ......exit */
    goto L620;
L370:

/*        this section of the program inspects for */
/*        negligible elements in the s and e arrays.  on */
/*        completion the variables kase and l are set as follows. */

/*           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m */
/*           kase = 2     if s(l) is negligible and l.lt.m */
/*           kase = 3     if e(l-1) is negligible, l.lt.m, and */
/*                        s(l), ..., s(m) are not negligible (qr step). */
/*           kase = 4     if e(m-1) is negligible (convergence). */

    i__1 = m;
    for (ll = 1; ll <= i__1; ++ll) {
	l = m - ll;
/*        ...exit */
	if (l == 0) {
	    goto L400;
	}
	test = (d__1 = s[l], abs(d__1)) + (d__2 = s[l + 1], abs(d__2));
	ztest = test + (d__1 = e[l], abs(d__1));
	if (ztest != test) {
	    goto L380;
	}
	e[l] = 0.;
/*        ......exit */
	goto L400;
L380:
/* L390: */
	;
    }
L400:
    if (l != m - 1) {
	goto L410;
    }
    kase = 4;
    goto L480;
L410:
    lp1 = l + 1;
    mp1 = m + 1;
    i__1 = mp1;
    for (lls = lp1; lls <= i__1; ++lls) {
	ls = m - lls + lp1;
/*           ...exit */
	if (ls == l) {
	    goto L440;
	}
	test = 0.;
	if (ls != m) {
	    test += (d__1 = e[ls], abs(d__1));
	}
	if (ls != l + 1) {
	    test += (d__1 = e[ls - 1], abs(d__1));
	}
	ztest = test + (d__1 = s[ls], abs(d__1));
	if (ztest != test) {
	    goto L420;
	}
	s[ls] = 0.;
/*           ......exit */
	goto L440;
L420:
/* L430: */
	;
    }
L440:
    if (ls != l) {
	goto L450;
    }
    kase = 3;
    goto L470;
L450:
    if (ls != m) {
	goto L460;
    }
    kase = 1;
    goto L470;
L460:
    kase = 2;
    l = ls;
L470:
L480:
    ++l;

/*        perform the task indicated by kase. */

    switch (kase) {
	case 1:  goto L490;
	case 2:  goto L520;
	case 3:  goto L540;
	case 4:  goto L570;
    }

/*        deflate negligible s(m). */

L490:
    mm1 = m - 1;
    f = e[m - 1];
    e[m - 1] = 0.;
    i__1 = mm1;
    for (kk = l; kk <= i__1; ++kk) {
	k = mm1 - kk + l;
	t1 = s[k];
	drotg_(&t1, &f, &cs, &sn);
	s[k] = t1;
	if (k == l) {
	    goto L500;
	}
	f = -sn * e[k - 1];
	e[k - 1] = cs * e[k - 1];
L500:
	if (wantv) {
	    drot_(p, &v[k * v_dim1 + 1], &c__1, &v[m * v_dim1 + 1], &c__1, &
		    cs, &sn);
	}
/* L510: */
    }
    goto L610;

/*        split at negligible s(l). */

L520:
    f = e[l - 1];
    e[l - 1] = 0.;
    i__1 = m;
    for (k = l; k <= i__1; ++k) {
	t1 = s[k];
	drotg_(&t1, &f, &cs, &sn);
	s[k] = t1;
	f = -sn * e[k];
	e[k] = cs * e[k];
	if (wantu) {
	    drot_(n, &u[k * u_dim1 + 1], &c__1, &u[(l - 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L530: */
    }
    goto L610;

/*        perform one qr step. */

L540:

/*           calculate the shift. */

/* Computing MAX */
    d__6 = (d__1 = s[m], abs(d__1)), d__7 = (d__2 = s[m - 1], abs(d__2)), 
	    d__6 = max(d__6,d__7), d__7 = (d__3 = e[m - 1], abs(d__3)), d__6 =
	     max(d__6,d__7), d__7 = (d__4 = s[l], abs(d__4)), d__6 = max(d__6,
	    d__7), d__7 = (d__5 = e[l], abs(d__5));
    scale = max(d__6,d__7);
    sm = s[m] / scale;
    smm1 = s[m - 1] / scale;
    emm1 = e[m - 1] / scale;
    sl = s[l] / scale;
    el = e[l] / scale;
/* Computing 2nd power */
    d__1 = emm1;
    b = ((smm1 + sm) * (smm1 - sm) + d__1 * d__1) / 2.;
/* Computing 2nd power */
    d__1 = sm * emm1;
    c__ = d__1 * d__1;
    shift = 0.;
    if (b == 0. && c__ == 0.) {
	goto L550;
    }
/* Computing 2nd power */
    d__1 = b;
    shift = sqrt(d__1 * d__1 + c__);
    if (b < 0.) {
	shift = -shift;
    }
    shift = c__ / (b + shift);
L550:
    f = (sl + sm) * (sl - sm) + shift;
    g = sl * el;

/*           chase zeros. */

    mm1 = m - 1;
    i__1 = mm1;
    for (k = l; k <= i__1; ++k) {
	drotg_(&f, &g, &cs, &sn);
	if (k != l) {
	    e[k - 1] = f;
	}
	f = cs * s[k] + sn * e[k];
	e[k] = cs * e[k] - sn * s[k];
	g = sn * s[k + 1];
	s[k + 1] = cs * s[k + 1];
	if (wantv) {
	    drot_(p, &v[k * v_dim1 + 1], &c__1, &v[(k + 1) * v_dim1 + 1], &
		    c__1, &cs, &sn);
	}
	drotg_(&f, &g, &cs, &sn);
	s[k] = f;
	f = cs * e[k] + sn * s[k + 1];
	s[k + 1] = -sn * e[k] + cs * s[k + 1];
	g = sn * e[k + 1];
	e[k + 1] = cs * e[k + 1];
	if (wantu && k < *n) {
	    drot_(n, &u[k * u_dim1 + 1], &c__1, &u[(k + 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L560: */
    }
    e[m - 1] = f;
    ++iter;
    goto L610;

/*        convergence. */

L570:

/*           make the singular value  positive. */

    if (s[l] >= 0.) {
	goto L580;
    }
    s[l] = -s[l];
    if (wantv) {
	dscal_(p, &c_b113, &v[l * v_dim1 + 1], &c__1);
    }
L580:

/*           order the singular value. */

L590:
    if (l == mm) {
	goto L600;
    }
/*           ...exit */
    if (s[l] >= s[l + 1]) {
	goto L600;
    }
    t = s[l];
    s[l] = s[l + 1];
    s[l + 1] = t;
    if (wantv && l < *p) {
	dswap_(p, &v[l * v_dim1 + 1], &c__1, &v[(l + 1) * v_dim1 + 1], &c__1);
    }
    if (wantu && l < *n) {
	dswap_(n, &u[l * u_dim1 + 1], &c__1, &u[(l + 1) * u_dim1 + 1], &c__1);
    }
    ++l;
    goto L590;
L600:
    iter = 0;
    --m;
L610:
    goto L360;
L620:
    return 0;
} /* dsvdc_ */

/* Subroutine */ int daxpy_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx, doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] += *da * dx[i__];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	dy[i__] += *da * dx[i__];
	dy[i__ + 1] += *da * dx[i__ + 1];
	dy[i__ + 2] += *da * dx[i__ + 2];
	dy[i__ + 3] += *da * dx[i__ + 3];
/* L50: */
    }
    return 0;
} /* daxpy_ */

doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static doublereal dtemp;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[i__] * dy[i__];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
		i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 
		4] * dy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */

doublereal dnrm2_(integer *n, doublereal *x, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer ix;
    static doublereal ssq, norm, scale, absxi;

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. */

/*  DNRM2 returns the euclidean norm of a vector via the function */
/*  name, so that */

/*     DNRM2 := sqrt( x'*x ) */



/*  -- This version written on 25-October-1982. */
/*     Modified on 14-October-1993 to inline the call to DLASSQ. */
/*     Sven Hammarling, Nag Ltd. */


/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (*n < 1 || *incx < 1) {
	norm = 0.;
    } else if (*n == 1) {
	norm = abs(x[1]);
    } else {
	scale = 0.;
	ssq = 1.;
/*        The following loop is equivalent to this call to the LAPACK */
/*        auxiliary routine: */
/*        CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */

	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	    if (x[ix] != 0.) {
		absxi = (d__1 = x[ix], abs(d__1));
		if (scale < absxi) {
/* Computing 2nd power */
		    d__1 = scale / absxi;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / scale;
		    ssq += d__1 * d__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

/*     End of DNRM2. */

} /* dnrm2_ */

/* Subroutine */ int drot_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ix, iy;
    static doublereal dtemp;


/*     applies a plane rotation. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = *c__ * dx[ix] + *s * dy[iy];
	dy[iy] = *c__ * dy[iy] - *s * dx[ix];
	dx[ix] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = *c__ * dx[i__] + *s * dy[i__];
	dy[i__] = *c__ * dy[i__] - *s * dx[i__];
	dx[i__] = dtemp;
/* L30: */
    }
    return 0;
} /* drot_ */

/* Subroutine */ int drotg_(doublereal *da, doublereal *db, doublereal *c__, 
	doublereal *s)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal r__, z__, roe, scale;


/*     construct givens plane rotation. */
/*     jack dongarra, linpack, 3/11/78. */


    roe = *db;
    if (abs(*da) > abs(*db)) {
	roe = *da;
    }
    scale = abs(*da) + abs(*db);
    if (scale != 0.) {
	goto L10;
    }
    *c__ = 1.;
    *s = 0.;
    r__ = 0.;
    z__ = 0.;
    goto L20;
L10:
/* Computing 2nd power */
    d__1 = *da / scale;
/* Computing 2nd power */
    d__2 = *db / scale;
    r__ = scale * sqrt(d__1 * d__1 + d__2 * d__2);
    r__ = d_sign(&c_b191, &roe) * r__;
    *c__ = *da / r__;
    *s = *db / r__;
    z__ = 1.;
    if (abs(*da) > abs(*db)) {
	z__ = *s;
    }
    if (abs(*db) >= abs(*da) && *c__ != 0.) {
	z__ = 1. / *c__;
    }
L20:
    *da = r__;
    *db = z__;
    return 0;
} /* drotg_ */

/* Subroutine */ int dscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, mp1, nincx;


/*     scales a vector by a constant. */
/*     uses unrolled loops for increment equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dx[i__] = *da * dx[i__];
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	dx[i__] = *da * dx[i__];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 5) {
	dx[i__] = *da * dx[i__];
	dx[i__ + 1] = *da * dx[i__ + 1];
	dx[i__ + 2] = *da * dx[i__ + 2];
	dx[i__ + 3] = *da * dx[i__ + 3];
	dx[i__ + 4] = *da * dx[i__ + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */

/* Subroutine */ int dswap_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static doublereal dtemp;


/*     interchanges two vectors. */
/*     uses unrolled loops for increments equal one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 12/3/93, array(1) declarations changed to array(*) */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal */
/*         to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = dx[ix];
	dx[ix] = dy[iy];
	dy[iy] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       code for both increments equal to 1 */


/*       clean-up loop */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = dx[i__];
	dx[i__] = dy[i__];
	dy[i__] = dtemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 3) {
	dtemp = dx[i__];
	dx[i__] = dy[i__];
	dy[i__] = dtemp;
	dtemp = dx[i__ + 1];
	dx[i__ + 1] = dy[i__ + 1];
	dy[i__ + 1] = dtemp;
	dtemp = dx[i__ + 2];
	dx[i__ + 2] = dy[i__ + 2];
	dy[i__ + 2] = dtemp;
/* L50: */
    }
    return 0;
} /* dswap_ */

