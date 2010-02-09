/* shpg.f -- translated by f2c (version 20061008).
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

static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b54 = 2.;
static doublereal c_b71 = 1.;


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

/*                           August 2003 */

/* ... in file shpg.f */

/*     this file contains code and documentation for subroutines */
/*     shpgi and shpg. */

/* ... files which must be loaded with shpg.f */

/*     hrfft.f */

/*     shpgi initializes the arrays wshp and iwshp for subsequent */
/*     use in subroutine shpg, which performs the harmonic projection */
/*     which is equivalent to a harmonic analysis followed by */
/*     harmonic synthesis but faster and with less memory. */
/*     (see description of subroutine shpg below). */

/*     subroutine shpgi(nlat,nlon,isym,mtrunc,wshp,lwshp,iwshp, */
/*    1 liwshp,work,lwork,ierror) */

/*     shpgi initializes arrays wshp and iwshp for repeated use */
/*     by subroutine shpg .... */

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
/*            nlon must be at least 4. */

/*     isym   currently not used, no equatorial symmetries assumed, */
/*            only whole sphere computations. */

/*     mtrunc the highest longitudinal wave number retained in the */
/*            projection. It must be less than or equal to */
/*            the minimum of nlat-1 and nlon/2. The first wave */
/*            number is zero. For example, if wave numbers 0 and */
/*            1 are desired then mtrunc = 1. */

/*     lwshp  the dimension of the array wshp as it appears in the */
/*            program that calls shpgi. It must be at least */
/*            2*(nlat+1)**2+nlon+log2(nlon) */

/*     liwshp the dimension of the array iwshp as it appears in the */
/*            program that calls shpgi. It must be at least */
/*            4*(nlat+1). */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls shpgi. It must be at least */
/*            1.25*(nlat+1)**2+7*nlat+8. */

/*     ************************************************************** */

/*     output parameters */

/*     wshp   a single precision array that must be saved for */
/*            repeated use by subroutine shpg. */

/*     iwshp  an integer array that must be saved for repeated */
/*            use by subroutine shpg. */

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

/* Subroutine */ int shpgi_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, doublereal *wshp, integer *lwshp, integer *iwshp, integer *
	liwshp, doublereal *work, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer iw1, iw2, iw3, lw1, iw4, jw1, jw2, jw3, jw4, kw1, kw2, kw3,
	     kw4, kw5, kw6, kw7, kw8, kw9, kw10, kw11, nte, mmax, mlwk, ktot, 
	    nloc1, nloc2, log2n;
    extern /* Subroutine */ int shpgi1_(integer *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), hrffti_(integer *, doublereal *);


    /* Parameter adjustments */
    --work;
    --iwshp;
    --wshp;

    /* Function Body */
    *ierror = 1;
    if (*nlat < 1) {
	return 0;
    }
    *ierror = 2;
    if (*nlon < 1) {
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
    log2n = log((doublereal) (*nlon)) / log(2.f);
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
    mlwk = i__1 * i__1 * 1.25f + *nlat * 7 + 8;
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
    kw4 = kw3 + (nte << 1);
    kw5 = kw4 + (nte << 1);
    kw6 = kw5 + nte;
    kw7 = kw6 + nte;
    kw8 = kw7 + (nte << 2);
    kw9 = kw8 + (nte << 1);
    kw10 = kw9 + nloc1;
    kw11 = kw10 + nloc1;
    ktot = kw11 + nte * nte;

    shpgi1_(nlat, nlon, isym, mtrunc, &nte, ierror, &wshp[iw1], &wshp[iw2], &
	    wshp[iw3], &wshp[iw4], &iwshp[jw1], &iwshp[jw2], &iwshp[jw3], &
	    iwshp[jw4], &work[kw1], &work[kw2], &work[kw3], &work[kw4], &work[
	    kw5], &work[kw6], &work[kw7], &work[kw8], &work[kw9], &work[kw10],
	     &work[kw11]);
    return 0;
} /* shpgi_ */

/* Subroutine */ int shpgi1_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, integer *idp, integer *ierror, doublereal *pe, doublereal *po, 
	doublereal *ze, doublereal *zo, integer *ipse, integer *jzse, integer *ipso, 
	integer *jzso, doublereal *cp, doublereal *wx, doublereal *thet, 
	doublereal *gwts, doublereal *xx, doublereal *z__, doublereal *a, 
	doublereal *b, doublereal *ped, doublereal *pod, doublereal *u)
{
    /* Format strings */
    static char fmt_160[] = "(\002 error in gaqd =\002,i5)";

    /* System generated locals */
    integer ped_dim1, ped_dim2, ped_offset, pod_dim1, pod_dim2, pod_offset, 
	    u_dim1, u_offset, pe_dim1, pe_dim2, pe_offset, po_dim1, po_dim2, 
	    po_offset, ze_dim1, ze_dim2, ze_offset, zo_dim1, zo_dim2, 
	    zo_offset, ipse_dim1, ipse_offset, jzse_dim1, jzse_offset, 
	    ipso_dim1, ipso_offset, jzso_dim1, jzso_offset, i__1, i__2, i__3;
    doublereal r__1, r__2, r__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, m, n;
    static doublereal a1, b1, c1;
    static integer ip;
    extern /* Subroutine */ int gs_(integer *, doublereal *, doublereal *, 
	    doublereal *);
    static integer it, js, mp1, ms2, ns2, nec;
    static doublereal dfn;
    static integer nem, noc, nte;
    static doublereal toe;
    static integer nom, nto, nmx;
    static doublereal sum;
    static doublereal sum1;
    extern doublereal rand_(integer *);
    static integer lock;
    static doublereal dmax__;
    static integer nshe[2], modn, ierr, nsho[2];
    static doublereal work;
    static doublereal tusl;
    static integer mxtr;
    static doublereal zort[8192]	/* was [64][64][2] */;
    extern /* Subroutine */ int dlfkg_(integer *, integer *, doublereal *), 
	    gaqdp_(integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *), dlftg_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), trunc_(integer *, integer *, integer 
	    *, doublereal *, integer *, integer *);
    static integer lwork;
    extern /* Subroutine */ int normal_(integer *, doublereal *, integer *, 
	    doublereal *);

    /* Fortran I/O blocks */
    static cilist io___37 = { 0, 6, 0, fmt_160, 0 };





    /* Parameter adjustments */
    --gwts;
    --thet;
    u_dim1 = *idp;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    pod_dim1 = *idp;
    pod_dim2 = *idp;
    pod_offset = 1 + pod_dim1 * (1 + pod_dim2);
    pod -= pod_offset;
    ped_dim1 = *idp;
    ped_dim2 = *idp;
    ped_offset = 1 + ped_dim1 * (1 + ped_dim2);
    ped -= ped_offset;
    --b;
    --a;
    --z__;
    --xx;
    --wx;
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
    tusl = 0.f;
    toe = 0.f;

/*     compute gauss grid distribution */

    lwork = *nlat + 1;
    gaqdp_(nlat, &thet[1], &gwts[1], &work, &lwork, &ierr);
    if (ierr != 0) {
	s_wsfe(&io___37);
	do_fio(&c__1, (char *)&ierr, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    i__1 = nto;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gwts[i__] += gwts[i__];
    }

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
	nem = (*nlat - m + 1) / 2;
	nec = nte - nem;

/*     compute associated legendre functions */

	if (m <= 1) {
	    i__2 = nem;
	    for (j = 1; j <= i__2; ++j) {
		n = j + j + m - 2;
		dlfkg_(&m, &n, &cp[1]);
		i__3 = nte;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    dlftg_(&m, &n, &thet[i__], &cp[1], &ped[i__ + (j + nec + 
			    ip * ped_dim2) * ped_dim1]);
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
			u[i__ + (j + nec) * u_dim1] = ped[i__ + (j + nec + ip 
				* ped_dim2) * ped_dim1];
		    }
		    goto L207;
		}
		a1 = b[n - 1] * a[n + m - 3] / a[n + m - 1];
		b1 = a[n - m + 1] / a[n + m - 1];
		if (n - m <= 1) {
		    i__3 = nte;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			u[i__ + (j + nec) * u_dim1] = a1 * ped[i__ + (j + nec 
				- 1 + ip * ped_dim2) * ped_dim1] - b1 * ped[
				i__ + (j + nec + ip * ped_dim2) * ped_dim1];
		    }
		} else {
		    c1 = b[n - 1] * a[n - m - 1] / a[n + m - 1];
		    i__3 = nte;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			u[i__ + (j + nec) * u_dim1] = a1 * ped[i__ + (j + nec 
				- 1 + ip * ped_dim2) * ped_dim1] - b1 * ped[
				i__ + (j + nec + ip * ped_dim2) * ped_dim1] + 
				c1 * u[i__ + (j + nec - 1) * u_dim1];
		    }
		}
L207:
		;
	    }
	    i__2 = nem;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = nte;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    ped[i__ + (j + nec + ip * ped_dim2) * ped_dim1] = u[i__ + 
			    (j + nec) * u_dim1];
		}
	    }
	}
	if (nec <= 0) {
	    goto L200;
	}

/*     generate orthogonal vector */

	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xx[i__] = rand_(&c__0);
	}

	it = 0;
L201:
	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] = 0.;
	    wx[i__] = gwts[i__] * xx[i__];
	}
	i__2 = nte;
	for (j = 1; j <= i__2; ++j) {
	    if (j == nec) {
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
	normal_(&nte, &xx[1], idp, &gwts[1]);
	++it;
	if (it <= 2) {
	    goto L201;
	}
	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ped[i__ + (nec + ip * ped_dim2) * ped_dim1] = xx[i__];
	}
L200:
	;
    }

/*     reorder if mtrunc is less than nlat-1 */
/*         case of even functions */

    nmx = *nlat - mxtr;
    if (modn == 1) {
	nshe[0] = nmx / 2;
	nshe[1] = (nmx - 1) / 2;
    } else {
	nshe[0] = (nmx - 1) / 2;
	nshe[1] = nmx / 2;
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
		sum = ped[j + (i__ + ip * ped_dim2) * ped_dim1] * gwts[j];
		ze[j + (i__ + ip * ze_dim2) * ze_dim1] = sum;
		pe[i__ + (j + ip * pe_dim2) * pe_dim1] = ped[i__ + (j + ip * 
			ped_dim2) * ped_dim1];
		if (abs(sum) > 5e-8 && lock == 0) {
		    lock = 1;
		    jzse[i__ + ip * jzse_dim1] = j;
		}
	    }
	}
/* L250: */
    }

/*     check orthogonality of pe(i,j,mp1)  mp1=1,2 */

    for (ip = 1; ip <= 2; ++ip) {
	dmax__ = 0.f;
	i__1 = nte;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = nte;
	    for (j = 1; j <= i__2; ++j) {
		sum1 = 0.f;
		i__3 = nte;
		for (k = 1; k <= i__3; ++k) {
		    sum1 += ze[k + (i__ + ip * ze_dim2) * ze_dim1] * pe[k + (
			    j + ip * pe_dim2) * pe_dim1];
		}
		zo[i__ + (j + ip * zo_dim2) * zo_dim1] = sum1;
		if (i__ != j) {
/* Computing MAX */
		    r__1 = dmax__, r__2 = dabs(sum1);
		    dmax__ = dmax(r__1,r__2);
		} else {
/* Computing MAX */
		    r__2 = dmax__, r__3 = (r__1 = sum1 - 1.f, dabs(r__1));
		    dmax__ = dmax(r__2,r__3);
		}
	    }
	}
    }

/*     compute n**2 basis (odd functions) */

    ip = 2;
    i__1 = mxtr + 1;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	ip = 3 - ip;
	m = mp1 - 1;
	ms2 = mp1 / 2;
	nem = (*nlat - m + 1) / 2;
	nom = *nlat - m - nem;
	noc = nto - nom;

/*     compute associated legendre functions */

	if (m <= 1) {
	    i__2 = nom;
	    for (j = 1; j <= i__2; ++j) {
		n = j + j + m - 1;
		dlfkg_(&m, &n, &cp[1]);
		i__3 = nte;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    dlftg_(&m, &n, &thet[i__], &cp[1], &pod[i__ + (j + noc + 
			    ip * pod_dim2) * pod_dim1]);
		}
		if (modn > 0) {
		    pod[nte + (j + noc + ip * pod_dim2) * pod_dim1] = 0.;
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
			u[i__ + (j + noc) * u_dim1] = pod[i__ + (j + noc + ip 
				* pod_dim2) * pod_dim1];
		    }
		    goto L304;
		}
		a1 = b[n - 1] * a[n + m - 3] / a[n + m - 1];
		b1 = a[n - m + 1] / a[n + m - 1];
		if (n - m <= 1) {
		    i__3 = nte;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			u[i__ + (j + noc) * u_dim1] = a1 * pod[i__ + (j + noc 
				- 1 + ip * pod_dim2) * pod_dim1] - b1 * pod[
				i__ + (j + noc + ip * pod_dim2) * pod_dim1];
		    }
		} else {
		    c1 = b[n - 1] * a[n - m - 1] / a[n + m - 1];
		    i__3 = nte;
		    for (i__ = 1; i__ <= i__3; ++i__) {
			u[i__ + (j + noc) * u_dim1] = a1 * pod[i__ + (j + noc 
				- 1 + ip * pod_dim2) * pod_dim1] - b1 * pod[
				i__ + (j + noc + ip * pod_dim2) * pod_dim1] + 
				c1 * u[i__ + (j + noc - 1) * u_dim1];
		    }
		}
L304:
		if (modn == 1) {
		    u[nte + (j + noc) * u_dim1] = 0.;
		}
/* L307: */
	    }
	    i__2 = nom;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = nte;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    pod[i__ + (j + noc + ip * pod_dim2) * pod_dim1] = u[i__ + 
			    (j + noc) * u_dim1];
		}
	    }
	}

	if (noc <= 0) {
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
	    z__[i__] = 0.f;
	    wx[i__] = gwts[i__] * xx[i__];
	}
	i__2 = nto;
	for (j = 1; j <= i__2; ++j) {
	    if (j == noc) {
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
	normal_(&nte, &xx[1], idp, &gwts[1]);
	++it;
	if (it <= 2) {
	    goto L306;
	}
	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    pod[i__ + (noc + ip * pod_dim2) * pod_dim1] = xx[i__];
	}
	if (modn == 1) {
	    pod[nte + (noc + ip * pod_dim2) * pod_dim1] = 0.;
	}
L300:
	;
    }

    nmx = *nlat - mxtr;
    if (modn == 1) {
	nsho[0] = (nmx - 1) / 2;
	nsho[1] = nmx / 2;
    } else {
	nsho[0] = nmx / 2;
	nsho[1] = (nmx - 1) / 2;
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
		sum = pod[j + (i__ + ip * pod_dim2) * pod_dim1] * gwts[j];
		zo[j + (i__ + ip * zo_dim2) * zo_dim1] = sum;
		po[i__ + (j + ip * po_dim2) * po_dim1] = pod[i__ + (j + ip * 
			pod_dim2) * pod_dim1];
		if (abs(sum) > 5e-8 && lock == 0) {
		    lock = 1;
		    jzso[i__ + ip * jzso_dim1] = j;
		}
	    }
	}
    }

/*     check orthogonality of po(i,j,mp1)  mp1=1,2 */

    for (ip = 1; ip <= 2; ++ip) {
	dmax__ = 0.f;
	i__1 = nto;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = nto;
	    for (j = 1; j <= i__2; ++j) {
		sum1 = 0.f;
		i__3 = nto;
		for (k = 1; k <= i__3; ++k) {
		    sum1 += zo[k + (i__ + ip * zo_dim2) * zo_dim1] * po[k + (
			    j + ip * po_dim2) * po_dim1];
		}
		zort[i__ + (j + (ip << 6) << 6) - 4161] = sum1;
		if (i__ != j) {
/* Computing MAX */
		    r__1 = dmax__, r__2 = dabs(sum1);
		    dmax__ = dmax(r__1,r__2);
		} else {
/* Computing MAX */
		    r__2 = dmax__, r__3 = (r__1 = sum1 - 1.f, dabs(r__1));
		    dmax__ = dmax(r__2,r__3);
		}
	    }
	}
    }
    return 0;
} /* shpgi1_ */



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

/* ... file shpg.f */

/* ... files which must be loaded with shpg.f */

/*     hrfft.f */

/*     shpg computes the harmonic projection, which is */
/*     equivalent to a harmonic analysis (forward) followed */
/*     by a harmonic synthesis (backward transform). */
/*     shpg uses the n**2 projection or complement when appropriate */
/*     as well as  odd/even factorization and zero truncation on an */
/*     on a Gaussian distributed grid as defined in the JCP paper */
/*     "Generalized discrete spherical harmonic transforms" */
/*     by Paul N. Swarztrauber and William F. Spotz */
/*     J. Comp. Phys., 159(2000) pp. 213-230. */

/*     subroutine shpg(nlat,nlon,isym,mtrunc,x,y,idxy, */
/*    1        wshp,lwshp,iwshp,liwshp,work,lwork,ierror) */

/*     shpg projects the array x onto the set of functions represented */
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
/*            nlon must be at least 4. */

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
/*            appear in the program that calls shpg. It must be */
/*            at least nlat. */

/*     wshp   a single precision array that must be saved for */
/*            repeated use by subroutine shpg. */

/*     lwshp  the dimension of the array wshp as it appears in the */
/*            program that calls shpgi. It must be at least */
/*            2*(nlat+1)**2+nlon+log2(nlon) */

/*     iwshp  an integer array that must be saved for repeated */
/*            use by subroutine shpg. */


/*     liwshp the dimension of the array iwshp as it appears in the */
/*            program that calls shpgi. It must be at least */
/*            4*(nlat+1). */

/*     work   a single precision work array that does */
/*            not have to be saved. */

/*     lwork  the dimension of the array work as it appears in the */
/*            program that calls shpg. It must be at least */
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

/* Subroutine */ int shpg_(integer *nlat, integer *nlon, integer *isym, 
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
	    mwrk, nloc1, nloc2, log2n;
    extern /* Subroutine */ int shpg1_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *), hrfftb_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *), hrfftf_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *);



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
    if (*nlat < 1) {
	return 0;
    }
    *ierror = 2;
    if (*nlon < 1) {
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
    log2n = log((doublereal) (*nlon)) / log(2.f);
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

    shpg1_(nlat, nlon, isym, mtrunc, &y[y_offset], &y[y_offset], idxy, ierror,
	     &nte, &wshp[iw1], &wshp[iw2], &wshp[iw3], &wshp[iw4], &iwshp[jw1]
	    , &iwshp[jw2], &iwshp[jw3], &iwshp[jw4], &work[jw1], &work[jw2], &
	    work[jw3], &work[jw4]);

    hrfftb_(nlat, nlon, &y[y_offset], idxy, &wshp[lw1 + 1], &work[1]);

    sn = 1.f / *nlon;
    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__ + j * y_dim1] = sn * y[i__ + j * y_dim1];
	}
    }
    return 0;
} /* shpg_ */

/* Subroutine */ int shpg1_(integer *nlat, integer *nlon, integer *isym, 
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
    static integer i__, j, m, ip, js, mp1, ms2, ns2, lag, nec, nem, noc, nte, 
	    mpm, nom, nto, nmx, nshe[2], modn, nsho[2], mxtr;
    extern /* Subroutine */ int tmxmx_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *);



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

/* Computing MIN */
    i__1 = *nlat - 1, i__2 = *nlon / 2, i__1 = min(i__1,i__2);
    mxtr = min(i__1,*mtrunc);
    nmx = *nlat - mxtr;
    if (modn == 1) {
	nshe[0] = nmx / 2;
	nshe[1] = (nmx - 1) / 2;
	nsho[0] = (nmx - 1) / 2;
	nsho[1] = nmx / 2;
    } else {
	nshe[0] = (nmx - 1) / 2;
	nshe[1] = nmx / 2;
	nsho[0] = nmx / 2;
	nsho[1] = (nmx - 1) / 2;
    }

    ip = 2;
    i__1 = mxtr + 1;
    for (mp1 = 1; mp1 <= i__1; ++mp1) {
	ip = 3 - ip;
	if (mxtr == *nlat - 1 && mp1 == 1) {
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sy[i__ + mp1 * sy_dim1] = sx[i__ + mp1 * sx_dim1];
	    }
/*      if(mp1.eq.2) then */
/*      sy(1,2) = 0. */
/*      sy(nlat,2) = 0. */
/*      end if */
/*      if(nlon.ge.3) then */
/*      sy(1,3) = 0. */
/*      sy(nlat,3) = 0. */
/*      do i=2,nlat-1 */
/*      sy(i,3) = sx(i,3) */
/*      end do */
/*      end if */
	    goto L100;
	}
	m = mp1 - 1;
/* Computing MAX */
	i__2 = 1, i__3 = m + m;
	mpm = max(i__2,i__3);
	ms2 = mp1 / 2;
/*      mrank = min(nlat-m,nlat-ms2-ms2) */
/*      nrank = nlat-mrank */
/*      nem = (mrank+1)/2-nshe(ip) */
/*      nom = mrank-(mrank+1)/2-nsho(ip) */
	nem = (*nlat - m + 1) / 2 - nshe[ip - 1];
	nom = (*nlat - m) / 2 - nsho[ip - 1];
	nec = nte - nem;
	noc = nto - nom;
	i__2 = nte;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    xe[i__ + xe_dim1] = (sx[i__ + mpm * sx_dim1] + sx[*nlat + 1 - i__ 
		    + mpm * sx_dim1]) * .5f;
	    xo[i__ + xo_dim1] = (sx[i__ + mpm * sx_dim1] - sx[*nlat + 1 - i__ 
		    + mpm * sx_dim1]) * .5f;
	}
/*      if(modn.eq.1) then */
/*      xe(nte,1) = sx(nte,mpm) */
/*      xo(nte,1) = 0. */
/*      end if */
	if (mpm < *nlon) {
	    i__2 = nte;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		xe[i__ + (xe_dim1 << 1)] = (sx[i__ + (mpm + 1) * sx_dim1] + 
			sx[*nlat + 1 - i__ + (mpm + 1) * sx_dim1]) * .5f;
		xo[i__ + (xo_dim1 << 1)] = (sx[i__ + (mpm + 1) * sx_dim1] - 
			sx[*nlat + 1 - i__ + (mpm + 1) * sx_dim1]) * .5f;
	    }
/*      if(modn.eq.1) then */
/*      xe(nte,2) = sx(nte,mpm+1) */
/*      xo(nte,2) = 0. */
/*      end if */
	}
	lag = 0;
	if (m == 0 || mpm == *nlon) {
	    lag = 1;
	}
	if (nec * 3 < nem << 1 || nem == 0) {
	    tmxmx_(&lag, &nte, &nec, idp, &pe[(ip * pe_dim2 + 1) * pe_dim1 + 
		    1], &nte, idp, &ze[(ip * ze_dim2 + 1) * ze_dim1 + 1], &xe[
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
	    tmxmx_(&lag, &nte, &nem, idp, &pe[(nec + 1 + ip * pe_dim2) * 
		    pe_dim1 + 1], &nte, idp, &ze[(nec + 1 + ip * ze_dim2) * 
		    ze_dim1 + 1], &xe[xe_offset], &ye[ye_offset], &ipse[nec + 
		    1 + ip * ipse_dim1], &jzse[nec + 1 + ip * jzse_dim1]);
	}
	if (noc * 3 < nom << 1 || nom == 0) {
	    tmxmx_(&lag, &nto, &noc, idp, &po[(ip * po_dim2 + 1) * po_dim1 + 
		    1], &nto, idp, &zo[(ip * zo_dim2 + 1) * zo_dim1 + 1], &xo[
		    xo_offset], &yo[yo_offset], &ipso[ip * ipso_dim1 + 1], &
		    jzso[ip * jzso_dim1 + 1]);
	    i__2 = nto;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		yo[i__ + yo_dim1] = xo[i__ + xo_dim1] - yo[i__ + yo_dim1];
	    }
	    if (mpm < *nlon && m != 0) {
		i__2 = nto;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    yo[i__ + (yo_dim1 << 1)] = xo[i__ + (xo_dim1 << 1)] - yo[
			    i__ + (yo_dim1 << 1)];
		}
	    }
	} else {
	    tmxmx_(&lag, &nto, &nom, idp, &po[(noc + 1 + ip * po_dim2) * 
		    po_dim1 + 1], &nto, idp, &zo[(noc + 1 + ip * zo_dim2) * 
		    zo_dim1 + 1], &xo[xo_offset], &yo[yo_offset], &ipso[noc + 
		    1 + ip * ipso_dim1], &jzso[noc + 1 + ip * jzso_dim1]);
	}
	i__2 = nto;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sy[i__ + mpm * sy_dim1] = ye[i__ + ye_dim1] + yo[i__ + yo_dim1];
	    sy[*nlat + 1 - i__ + mpm * sy_dim1] = ye[i__ + ye_dim1] - yo[i__ 
		    + yo_dim1];
	}
	if (nte > nto) {
	    sy[nte + mpm * sy_dim1] = ye[nte + ye_dim1];
	}
	if (mpm < *nlon && m != 0) {
	    i__2 = nto;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sy[i__ + (mpm + 1) * sy_dim1] = ye[i__ + (ye_dim1 << 1)] + yo[
			i__ + (yo_dim1 << 1)];
		sy[*nlat + 1 - i__ + (mpm + 1) * sy_dim1] = ye[i__ + (ye_dim1 
			<< 1)] - yo[i__ + (yo_dim1 << 1)];
	    }
	    if (nte > nto) {
		sy[nte + (mpm + 1) * sy_dim1] = ye[nte + (ye_dim1 << 1)];
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
	    sy[i__ + j * sy_dim1] = 0.f;
	}
    }
    return 0;
} /* shpg1_ */

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
	    c__[i__ + j * c_dim1] = 0.f;
	    i__3 = *lc;
	    for (k = 1; k <= i__3; ++k) {
		c__[i__ + j * c_dim1] += a[i__ + k * a_dim1] * b[k + j * 
			b_dim1];
	    }
	}
    }
    return 0;
} /* mxm_ */

/* Subroutine */ static int smxm_(integer *lr, integer *lc, integer *ld, doublereal *a, 
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
	    c__[i__ + j * c_dim1] = 0.f;
	    i__3 = *lc;
	    for (k = 1; k <= i__3; ++k) {
		c__[i__ + j * c_dim1] += a[i__ + k * a_dim1] * b[k + j * 
			b_dim1];
	    }
	}
    }
    return 0;
} /* smxm_ */

/* Subroutine */ static int mxmx_(integer *lr, integer *lc, integer *ld, doublereal *a, 
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
	y[k + y_dim1] = 0.f;
	y[k + (y_dim1 << 1)] = 0.f;
    }

    if (*lc <= 0) {
	return 0;
    }
    i__1 = *lc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum1 = 0.f;
	sum2 = 0.f;
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

/* Subroutine */ static int dmxmx_(integer *lr, integer *lc, integer *ld, doublereal 
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
	y[k + y_dim1] = 0.f;
	y[k + (y_dim1 << 1)] = 0.f;
    }

    if (*lc <= 0) {
	return 0;
    }
    i__1 = *lc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum1 = 0.f;
	sum2 = 0.f;
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

/* Subroutine */ static int tmxmx_(integer *lag, integer *lr, integer *lc, integer *
	ld, doublereal *a, integer *mc, integer *md, doublereal *b, doublereal *x, doublereal *y, 
	integer *is, integer *js)
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
    if (*lag == 1) {
	i__1 = kmx;
	for (k = 1; k <= i__1; ++k) {
	    y[k + y_dim1] = 0.f;
	}
/*      if(lc.eq.0) then */
/*      do k=1,lr */
/*      y(k,1) = x(k,1) */
/*      end do */
/*      return */
/*      end if */
	if (*lc <= 0) {
	    return 0;
	}
	i__1 = *lc;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sum1 = 0.f;
	    i__2 = *mc;
	    for (j = js[i__]; j <= i__2; ++j) {
		sum1 += b[j + i__ * b_dim1] * x[j + x_dim1];
	    }
	    i__2 = *lr;
	    for (k = is[i__]; k <= i__2; ++k) {
		y[k + y_dim1] += sum1 * a[k + i__ * a_dim1];
	    }
	}
	return 0;
    }
    i__1 = kmx;
    for (k = 1; k <= i__1; ++k) {
	y[k + y_dim1] = 0.f;
	y[k + (y_dim1 << 1)] = 0.f;
    }
    if (*lc <= 0) {
	return 0;
    }

    i__1 = *lc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum1 = 0.f;
	sum2 = 0.f;
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

/* Subroutine */ static int trunc_(integer *irc, integer *n, integer *idp, 
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

/* Subroutine */ static int gs_(integer *n, doublereal *x, doublereal *y, doublereal 
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
    sum = 0.f;
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

/* Subroutine */ static int normal_(integer *n, doublereal *x, integer *id, 
	doublereal *q)
{
    /* Format strings */
    static char fmt_3[] = "(\002 norm of z is zero in subroutine normal\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal sqs;

    /* Fortran I/O blocks */
    static cilist io___132 = { 0, 6, 0, fmt_3, 0 };



/*     normalize x */

    /* Parameter adjustments */
    --q;
    --x;

    /* Function Body */
    sqs = 0.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*      sum = 0. */
/*      do j=1,n */
/*      sum = sum+q(i,j)*x(j) */
/*      end do */
/*      sqs = sqs+sum*x(i) */
	sqs += q[i__] * x[i__] * x[i__];
    }

    if (sqs != 0.) {
	goto L4;
    }
    s_wsfe(&io___132);
    e_wsfe();
    return 0;
L4:
    sqs = sqrt(sqs);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] /= sqs;
    }
    return 0;
} /* normal_ */

/* Subroutine */ static int coe_(integer *moe, integer *n, doublereal *x, doublereal 
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
    *dmax__ = 0.f;
    if (*moe != 0) {
	goto L1;
    }
    i__1 = nh;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = *dmax__, d__3 = (d__1 = x[i__] - x[*n - i__ + 1], abs(d__1));
	*dmax__ = max(d__2,d__3);
	x[i__] = (x[i__] + x[*n - i__ + 1]) * .5f;
	x[*n - i__ + 1] = x[i__];
    }
    return 0;
L1:
    i__1 = nh;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__2 = *dmax__, d__3 = (d__1 = x[i__] + x[*n - i__ + 1], abs(d__1));
	*dmax__ = max(d__2,d__3);
	x[i__] = (x[i__] - x[*n - i__ + 1]) * .5f;
	x[*n - i__ + 1] = -x[i__];
    }
    if (*n % 2 != 0) {
	x[nh] = 0.f;
    }
    return 0;
} /* coe_ */

/*     subroutine dlfkg(m,n,cp) */

/*     subroutine dlfkg computes the coefficients in the trigonometric */
/*     expansion of the normalized associated legendre functions: */

/*     pbar(m,n,theta) = sqrt((2*n+1)*factorial(n-m)/(2*factorial(n+m))) */
/*                        *sin(theta)**m/(2**n*factorial(n)) times the */
/*                        (n+m)th derivative of (x**2-1)**n with respect */
/*                        to x=cos(theta) */

/*     where theta is colatitude. */

/*     subroutine dlfkg computes the coefficients cp(k) in the */
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
/* Subroutine */ int dlfkg_(integer *m, integer *n, doublereal *cp)
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
    cp[1] = 0.f;
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
	fnum += 2.f;
	fden += 2.f;
/* L18: */
    }
L20:
    i__1 = *n - 1 - nex;
    t1 /= pow_di(&c_b54, &i__1);
    if (ma / 2 % 2 != 0) {
	t1 = -t1;
    }
    t2 = 1.f;
    if (ma == 0) {
	goto L26;
    }
    i__1 = ma;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t2 = fnmh * t2 / (fnmh + pm1);
	fnmh += 2.f;
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
    a1 = (fk - 2.f) * (fk - 1.f) - fnnp1;
    b1 = (fk * fk - fnmsq) * 2.f;
    cp[l - 1] = b1 * cp[l] / a1;
L30:
    --l;
    if (l <= 1) {
	return 0;
    }
    fk += -2.f;
    a1 = (fk - 2.f) * (fk - 1.f) - fnnp1;
    b1 = (fk * fk - fnmsq) * -2.f;
    c1 = (fk + 1.f) * (fk + 2.f) - fnnp1;
    cp[l - 1] = -(b1 * cp[l] + c1 * cp[l + 1]) / a1;
    goto L30;
} /* dlfkg_ */

/* Subroutine */ int dlftg_(integer *m, integer *n, doublereal *theta, 
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
    *pb = cp[1] * .5f;
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
    *pb = 0.f;
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
    *pb = 0.f;
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
    *pb = 0.f;
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
} /* dlftg_ */


/* Subroutine */ int gaqdp_(integer *nlat, doublereal *theta, doublereal *wts,
	 doublereal *w, integer *lwork, integer *ierror)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double acos(doublereal), sqrt(doublereal), atan(doublereal), cos(
	    doublereal), sin(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal x, pb, pi, cz;
    static integer it, ns2;
    static doublereal dpb;
    static integer idx;
    static doublereal eps;
    static integer nix;
    static doublereal sum, pis2, dcor, cmax;
    static doublereal sgnd;
    static doublereal zero;
    extern /* Subroutine */ int cpdp1_(integer *, doublereal *, doublereal *, 
	    doublereal *), tpdp1_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer nhalf, mnlat;
    static doublereal zhold;
    extern doublereal dzepp_(doublereal *);
    static doublereal zlast, zprev, dthalf, dtheta;


/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
/*  .                                                             . */
/*  .                  copyright (c) 2001 by ucar                 . */
/*  .                                                             . */
/*  .       university corporation for atmospheric research       . */
/*  .                                                             . */
/*  .                      all rights reserved                    . */
/*  .                                                             . */
/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/*                             April 2002 */

/*     gauss points and weights are computed using the fourier-newton */
/*     described in "on computing the points and weights for */
/*     gauss-legendre quadrature", paul n. swarztrauber, siam journal */
/*     on scientific computing that has been accepted for publication. */
/*     This routine is faster and more accurate than older program */
/*     with the same name. */

/*     subroutine gaqdp computes the nlat gaussian colatitudes and weights */
/*     in double precision. the colatitudes are in radians and lie in the */
/*     in the interval (0,pi). */

/*     input parameters */

/*     nlat    the number of gaussian colatitudes in the interval (0,pi) */
/*             (between the two poles).  nlat must be greater than zero. */

/*     w       unused variable that permits a simple exchange with the */
/*             old routine with the same name in spherepack. */

/*     lwork   unused variable that permits a simple exchange with the */
/*             old routine with the same name in spherepack. */

/*     output parameters */

/*     theta   a double precision array with length nlat */
/*             containing the gaussian colatitudes in */
/*             increasing radians on the interval (0,pi). */

/*     wts     a double precision array with lenght nlat */
/*             containing the gaussian weights. */

/*     ierror = 0 no errors */
/*            = 1 if nlat.le.0 */

/*  ***************************************************************** */


/*     check work space length */

    /* Parameter adjustments */
    --wts;
    --theta;

    /* Function Body */
    *ierror = 1;
    if (*nlat <= 0) {
	return 0;
    }
    *ierror = 0;

/*     compute weights and points analytically when nlat=1,2 */

    if (*nlat == 1) {
	theta[1] = acos(0.);
	wts[1] = 2.;
	return 0;
    }
    if (*nlat == 2) {
	x = sqrt(.33333333333333331);
	theta[1] = acos(x);
	theta[2] = acos(-x);
	wts[1] = 1.;
	wts[2] = 1.;
	return 0;
    }
    eps = sqrt(dzepp_(&c_b71));
    eps *= sqrt(eps);
    pis2 = atan(1.) * 2.;
    pi = pis2 + pis2;
    mnlat = *nlat % 2;
    ns2 = *nlat / 2;
    nhalf = (*nlat + 1) / 2;
    idx = ns2 + 2;

    cpdp1_(nlat, &cz, &theta[ns2 + 1], &wts[ns2 + 1]);

    dtheta = pis2 / nhalf;
    dthalf = dtheta / 2.;
    cmax = dtheta * .2;

/*     estimate first point next to theta = pi/2 */

    if (mnlat != 0) {
	zero = pis2 - dtheta;
	zprev = pis2;
	nix = nhalf - 1;
    } else {
	zero = pis2 - dthalf;
	nix = nhalf;
    }
L9:
    it = 0;
L10:
    ++it;
    zlast = zero;

/*     newton iterations */

    tpdp1_(nlat, &zero, &cz, &theta[ns2 + 1], &wts[ns2 + 1], &pb, &dpb);
    dcor = pb / dpb;
    sgnd = 1.f;
    if (dcor != 0.) {
	sgnd = dcor / abs(dcor);
    }
/* Computing MIN */
    d__1 = abs(dcor);
    dcor = sgnd * min(d__1,cmax);
    zero -= dcor;
    if ((d__1 = zero - zlast, abs(d__1)) > eps * abs(zero)) {
	goto L10;
    }
    theta[nix] = zero;
    zhold = zero;
/*      wts(nix) = (nlat+nlat+1)/(dpb*dpb) */

/*     yakimiw's formula permits using old pb and dpb */

/* Computing 2nd power */
    d__1 = dpb + pb * cos(zlast) / sin(zlast);
    wts[nix] = (*nlat + *nlat + 1) / (d__1 * d__1);
    --nix;
    if (nix == 0) {
	goto L30;
    }
    if (nix == nhalf - 1) {
	zero = zero * 3.f - pi;
    }
    if (nix < nhalf - 1) {
	zero = zero + zero - zprev;
    }
    zprev = zhold;
    goto L9;

/*     extend points and weights via symmetries */

L30:
    if (mnlat != 0) {
	theta[nhalf] = pis2;
	tpdp1_(nlat, &pis2, &cz, &theta[ns2 + 1], &wts[ns2 + 1], &pb, &dpb);
	wts[nhalf] = (*nlat + *nlat + 1) / (dpb * dpb);
    }
    i__1 = ns2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wts[*nlat - i__ + 1] = wts[i__];
	theta[*nlat - i__ + 1] = pi - theta[i__];
    }
    sum = 0.;
    i__1 = *nlat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum += wts[i__];
    }
    i__1 = *nlat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wts[i__] = wts[i__] * 2. / sum;
    }
    return 0;
} /* gaqdp_ */

/* Subroutine */ int cpdp1_(integer *n, doublereal *cz, doublereal *cp, 
	doublereal *dcp)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static doublereal t1, t2, t3, t4;
    static integer ncp;


/*     computes the fourier coefficients of the legendre */
/*     polynomial p_n^0 and its derivative. */
/*     n is the degree and n/2 or (n+1)/2 */
/*     coefficients are returned in cp depending on whether */
/*     n is even or odd. The same number of coefficients */
/*     are returned in dcp. For n even the constant */
/*     coefficient is returned in cz. */

    /* Parameter adjustments */
    --dcp;
    --cp;

    /* Function Body */
    ncp = (*n + 1) / 2;
    t1 = -1.;
    t2 = *n + 1.;
    t3 = 0.;
    t4 = *n + *n + 1.;
    if (*n % 2 == 0) {
	cp[ncp] = 1.;
	for (j = ncp; j >= 2; --j) {
	    t1 += 2.;
	    t2 += -1.;
	    t3 += 1.;
	    t4 += -2.;
	    cp[j - 1] = t1 * t2 / (t3 * t4) * cp[j];
	}
	t1 += 2.;
	t2 += -1.;
	t3 += 1.;
	t4 += -2.;
	*cz = t1 * t2 / (t3 * t4) * cp[1];
	i__1 = ncp;
	for (j = 1; j <= i__1; ++j) {
	    dcp[j] = (j + j) * cp[j];
	}
    } else {
	cp[ncp] = 1.;
	for (j = ncp - 1; j >= 1; --j) {
	    t1 += 2.;
	    t2 += -1.;
	    t3 += 1.;
	    t4 += -2.;
	    cp[j] = t1 * t2 / (t3 * t4) * cp[j + 1];
	}
	i__1 = ncp;
	for (j = 1; j <= i__1; ++j) {
	    dcp[j] = (j + j - 1) * cp[j];
	}
    }
    return 0;
} /* cpdp1_ */

/* Subroutine */ int tpdp1_(integer *n, doublereal *theta, doublereal *cz, 
	doublereal *cp, doublereal *dcp, doublereal *pb, doublereal *dpb)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer k;
    static doublereal fn, chh, cdt;
    static integer kdo;
    static doublereal cth, sdt, sth;


/*     computes pn(theta) and its derivative dpb(theta) with */
/*     respect to theta */


    /* Parameter adjustments */
    --dcp;
    --cp;

    /* Function Body */
    fn = (doublereal) (*n);
    cdt = cos(*theta + *theta);
    sdt = sin(*theta + *theta);
    if (*n % 2 == 0) {

/*     n even */

	kdo = *n / 2;
	*pb = *cz * .5;
	*dpb = 0.;
	if (*n > 0) {
	    cth = cdt;
	    sth = sdt;
	    i__1 = kdo;
	    for (k = 1; k <= i__1; ++k) {
/*      pb = pb+cp(k)*cos(2*k*theta) */
		*pb += cp[k] * cth;
/*      dpb = dpb-(k+k)*cp(k)*sin(2*k*theta) */
		*dpb -= dcp[k] * sth;
		chh = cdt * cth - sdt * sth;
		sth = sdt * cth + cdt * sth;
		cth = chh;
/* L170: */
	    }
	}
    } else {

/*     n odd */

	kdo = (*n + 1) / 2;
	*pb = 0.;
	*dpb = 0.;
	cth = cos(*theta);
	sth = sin(*theta);
	i__1 = kdo;
	for (k = 1; k <= i__1; ++k) {
/*      pb = pb+cp(k)*cos((2*k-1)*theta) */
	    *pb += cp[k] * cth;
/*      dpb = dpb-(k+k-1)*cp(k)*sin((2*k-1)*theta) */
	    *dpb -= dcp[k] * sth;
	    chh = cdt * cth - sdt * sth;
	    sth = sdt * cth + cdt * sth;
	    cth = chh;
/* L190: */
	}
    }
    return 0;
} /* tpdp1_ */

doublereal dzepp_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;
    doublereal d__1;

    /* Local variables */
    static doublereal a, b, c__, eps;


/*     estimate unit roundoff in quantities of size x. */


/*     this program should function properly on all systems */
/*     satisfying the following two assumptions, */
/*        1.  the base used in representing floating point */
/*            numbers is not a power of three. */
/*        2.  the quantity  a  in statement 10 is represented to */
/*            the accuracy used in floating point variables */
/*            that are stored in memory. */
/*     the statement number 10 and the go to 10 are intended to */
/*     force optimizing compilers to generate code satisfying */
/*     assumption 2. */
/*     under these assumptions, it should be true that, */
/*            a  is not exactly equal to four-thirds, */
/*            b  has a zero for its last bit or digit, */
/*            c  is not exactly equal to one, */
/*            eps  measures the separation of 1.0 from */
/*                 the next larger floating point number. */
/*     the developers of eispack would appreciate being informed */
/*     about any systems where these assumptions do not hold. */

/*     this version dated 4/6/83. */

    a = 1.3333333333333333;
L10:
    b = a - 1.;
    c__ = b + b + b;
    eps = (d__1 = c__ - 1., abs(d__1));
    if (eps == 0.) {
	goto L10;
    }
    ret_val = eps * abs(*x);
    return ret_val;
} /* dzepp_ */

