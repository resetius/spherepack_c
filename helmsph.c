/* helmsph.f -- translated by f2c (version 20061008).
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


/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */
/*  .                                                             . */
/*  .                  Copyright (C) 1998 by UCAR                 . */
/*  .                                                             . */
/*  .       University Corporation for Atmospheric Research       . */
/*  .                                                             . */
/*  .                      All Rights Reserved                    . */
/*  .                                                             . */
/*  .                                                             . */
/*  .                          SPHEREPACK                         . */
/*  .                                                             . */
/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */



/* ... file helmsph.f */

/*     this file contains a program for solving the Helmholtz */
/*     equation with constant 1.0 on a ten degree grid on the full sphere */

/* ... required spherepack3.0 files */

/*     islapec.f, shaec.f, shsec.f, sphcom.f, hrfft.f */

/* ... description */

/*     let theta be latitude and phi be east longitude in radians. */
/*     and let */


/*       x = cos(theta)*sin(phi) */
/*       y = cos(theta)*cos(phi) */
/*       z = sint(theta) */

/*     be the cartesian coordinates corresponding to theta and phi. */
/*     on the unit sphere.  The exact solution */

/*        ue(theta,phi) = (1.+x*y)*exp(z) */

/*     is used to set the right hand side and compute error. */


/* ********************************************************************** */

/* OUTPUT FROM EXECUTING THE PROGRAM BELOW */
/* WITH 32 AND 64 BIT FLOATING POINT ARITHMETIC */

/* Helmholtz approximation on a ten degree grid */
/* nlat = 19   nlon = 36 */
/* xlmbda =  1.00   pertrb =  0.000E+00 */
/* maximum error =  0.715E-06 *** (32 BIT) */
/* maximum error =  0.114E-12 *** (64 BIT) */

/* *********************************************** */
/* *********************************************** */
/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_200[] = "(\002 shaeci, ierror = \002,i2)";
    static char fmt_201[] = "(\002 shseci, ierror = \002,i2)";
    static char fmt_202[] = "(\002 shaec , ierror = \002,i2)";
    static char fmt_100[] = "(\002 helmholtz approximation on a ten degree g"
	    "rid\002/\002 nlat = \002,i3,2x,\002 nlon = \002,i3)";
    static char fmt_103[] = "(\002 islapec, ierror = \002,i2)";
    static char fmt_204[] = "(\002 xlmbda = \002,f5.2,2x,\002 pertrb = \002,"
	    "e10.3,/\002 maximum error = \002,e10.3)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double atan(doublereal), sin(doublereal), cos(doublereal), exp(doublereal)
	    ;
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal a[361]	/* was [19][19] */, b[361]	/* was [19][
	    19] */;
    static integer i__, j;
    static doublereal r__[684]	/* was [19][36] */, u[684]	/* was [19][
	    36] */, x, y, z__, pi, ue, ez;
    static integer nt;
    static doublereal phi, dlat, dlon;
    static integer nlat;
    static doublereal cosp[36], errm;
    static integer nlon;
    static doublereal cost[19], sinp[36];
    extern /* Subroutine */ int exit_(integer *);
    static doublereal sint[19];
    static integer isym;
    static doublereal work[3249];
    extern /* Subroutine */ int shaec_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *);
    static doublereal theta, dwork[20];
    static integer lwork;
    extern /* Subroutine */ int shaeci_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    static integer lshaec;
    static doublereal xlmbda, wshaec[1451];
    extern /* Subroutine */ int shseci_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    static integer lshsec;
    static doublereal wshsec[1451], pertrb;
    static integer ierror, ldwork;
    extern /* Subroutine */ int islapec_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static cilist io___26 = { 0, 6, 0, fmt_200, 0 };
    static cilist io___28 = { 0, 6, 0, fmt_201, 0 };
    static cilist io___34 = { 0, 6, 0, fmt_202, 0 };
    static cilist io___35 = { 0, 6, 0, fmt_100, 0 };
    static cilist io___38 = { 0, 6, 0, fmt_103, 0 };
    static cilist io___42 = { 0, 6, 0, fmt_204, 0 };



/*     set grid size with parameter statements */


/*     set saved and unsaved work space lengths in terms of nnlat,nnlon */
/*     (see documentation for shaec,shsec,islapec) */


/*     set double precision work space length for initializations */


/*     dimension arrays */

    pi = atan(1.) * 4.;

/*     set helmholtz constant */

    xlmbda = 1.;

/*     set work space length arguments */

    lwork = 3249;
    ldwork = 20;
    lshaec = 1451;
    lshsec = 1451;

/*     set grid size arguments */

    nlat = 19;
    nlon = 36;

/*     set sine and cosine vectors */

    dlat = pi / (nlat - 1);
    dlon = (pi + pi) / nlon;
    i__1 = nlat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	theta = pi * -.5 + (i__ - 1) * dlat;
	sint[i__ - 1] = sin(theta);
	cost[i__ - 1] = cos(theta);
    }
    i__1 = nlon;
    for (j = 1; j <= i__1; ++j) {
	phi = (j - 1) * dlon;
	sinp[j - 1] = sin(phi);
	cosp[j - 1] = cos(phi);
    }

/*     set right hand side as helmholtz operator */
/*     applied to ue = (1.+x*y)*exp(z) */

    i__1 = nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x = cost[i__ - 1] * cosp[j - 1];
	    y = cost[i__ - 1] * sinp[j - 1];
	    z__ = sint[i__ - 1];
	    r__[i__ + j * 19 - 20] = -(x * y * (z__ * z__ + (z__ + 1.) * 6.) 
		    + z__ * (z__ + 2.)) * exp(z__);
	}
    }

/*     initialize saved work space arrays for scalar harmonic */
/*     analysis and Helmholtz inversion of r */

    shaeci_(&nlat, &nlon, wshaec, &lshaec, dwork, &ldwork, &ierror);
    if (ierror > 0) {
	s_wsfe(&io___26);
	do_fio(&c__1, (char *)&ierror, (ftnlen)sizeof(integer));
	e_wsfe();
	exit_(&c__0);
    }
    shseci_(&nlat, &nlon, wshsec, &lshsec, dwork, &ldwork, &ierror);
    if (ierror > 0) {
	s_wsfe(&io___28);
	do_fio(&c__1, (char *)&ierror, (ftnlen)sizeof(integer));
	e_wsfe();
	exit_(&c__0);
    }

/*     set no symmetry and one array */

    isym = 0;
    nt = 1;

/*     compute coefficients of r for input to islapec */

    shaec_(&nlat, &nlon, &isym, &nt, r__, &nlat, &nlon, a, b, &nlat, &nlat, 
	    wshaec, &lshaec, work, &lwork, &ierror);
    if (ierror > 0) {
	s_wsfe(&io___34);
	do_fio(&c__1, (char *)&ierror, (ftnlen)sizeof(integer));
	e_wsfe();
	exit_(&c__0);
    }

/*     solve Helmholtz equation on the sphere in u */

    s_wsfe(&io___35);
    do_fio(&c__1, (char *)&nlat, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nlon, (ftnlen)sizeof(integer));
    e_wsfe();
    islapec_(&nlat, &nlon, &isym, &nt, &xlmbda, u, &nlat, &nlon, a, b, &nlat, 
	    &nlat, wshsec, &lshsec, work, &lwork, &pertrb, &ierror);
    if (ierror != 0) {
	s_wsfe(&io___38);
	do_fio(&c__1, (char *)&ierror, (ftnlen)sizeof(integer));
	e_wsfe();
	if (ierror > 0) {
	    exit_(&c__0);
	}
    }

/*     compute and print maximum error in u */

    errm = 0.;
    i__1 = nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    x = cost[i__ - 1] * cosp[j - 1];
	    y = cost[i__ - 1] * sinp[j - 1];
	    z__ = sint[i__ - 1];
	    ez = exp(z__);
	    ue = (x * y + 1.) * ez;
/* Computing MAX */
	    d__2 = errm, d__3 = (d__1 = u[i__ + j * 19 - 20] - ue, abs(d__1));
	    errm = max(d__2,d__3);
	}
    }
    s_wsfe(&io___42);
    do_fio(&c__1, (char *)&xlmbda, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&pertrb, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&errm, (ftnlen)sizeof(doublereal));
    e_wsfe();
    return 0;
} /* MAIN__ */

/* Main program alias */ int helmsph_ () { MAIN__ (); return 0; }
