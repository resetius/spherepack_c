/* sshifte.f -- translated by f2c (version 20061008).
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


/* ... file sshifte.f contains code and documentation for subroutine sshifte */
/*     and its' initialization subroutine sshifti */

/* ... required files off spherepack3.0 */

/*     hrfft.f */

/*     subroutine sshifte(ioff,nlon,nlat,goff,greg,wsav,lsav,work,lwork,ier) */

/* *** purpose */

/*     subroutine sshifte does a highly accurate 1/2 grid increment shift */
/*     in both longitude and latitude of equally spaced data on the sphere. */
/*     data is transferred between the nlon by nlat "offset grid" in goff */
/*     (which excludes poles) and the nlon by nlat+1 "regular grid" in greg */
/*     (which includes poles).  the transfer can go from goff to greg or from */
/*     greg to goff (see ioff).  the grids which underly goff and greg are */
/*     described below.  the north and south poles are at latitude 0.5*pi and */
/*     -0.5*pi radians respectively where pi = 4.*atan(1.). */

/* *** grid descriptions */

/*     let dlon = (pi+pi)/nlon and dlat = pi/nlat be the uniform grid */
/*     increments in longitude and latitude */

/*     offset grid */

/*     the "1/2 increment offset" grid (long(j),lat(i)) on which goff(j,i) */
/*     is given (ioff=0) or generated (ioff=1) is */

/*          long(j) =0.5*dlon + (j-1)*dlon  (j=1,...,nlon) */

/*     and */

/*          lat(i) = -0.5*pi + 0.5*dlat + (i-1)*dlat (i=1,...,nlat) */

/*     the data in goff is "shifted" one half a grid increment in longitude */
/*     and latitude and excludes the poles.  each goff(j,1) is given at */
/*     latitude -0.5*pi+0.5*dlat and goff(j,nlat) is given at 0.5*pi-0.5*dlat */
/*     (1/2 a grid increment away from the poles).  each goff(1,i),goff(nlon,i) */
/*     is given at longitude 0.5*dlon and 2.*pi-0.5*dlon. */

/*     regular grid */

/*     let dlat,dlon be as above.  then the nlon by nlat+1 grid on which */
/*     greg(j,i) is generated (ioff=0) or given (ioff=1) is given by */

/*          lone(j) = (j-1)*dlon (j=1,...,nlon) */

/*      and */

/*          late(i) = -0.5*pi + (i-1)*dlat (i=1,...,nlat+1) */

/*     values in greg include the poles and start at zero degrees longitude. */

/* *** remark */

/*     subroutine sshifte can be used in conjunction with subroutine trssph */
/*     when transferring data from an equally spaced "1/2 increment offset" */
/*     grid to a gaussian or equally spaced grid (which includes poles) of */
/*     any resolution.  this problem (personal communication with dennis */
/*     shea) is encountered in geophysical modeling and data analysis. */

/* *** method */

/*     fast fourier transform software from spherepack2 and trigonometric */
/*     identities are used to accurately "shift" periodic vectors half a */
/*     grid increment in latitude and longitude.  latitudinal shifts are */
/*     accomplished by setting periodic 2*nlat vectors over the pole for each */
/*     longitude.  when nlon is odd, this requires an additional longitude */
/*     shift.  longitudinal shifts are then executed for each shifted latitude. */
/*     when necessary (ioff=0) poles are obtained by averaging the nlon */
/*     shifted polar values. */

/* *** required files from spherepack3.0 */

/*     hrfft.f */

/* *** argument description */

/* ... ioff */

/*     ioff = 0 if values on the offset grid in goff are given and values */
/*              on the regular grid in greg are to be generated. */

/*     ioff = 1 if values on the regular grid in greg are given and values */
/*              on the offset grid in goff are to be generated. */

/* ... nlon */

/*     the number of longitude points on both the "offset" and "regular" */
/*     uniform grid in longitude (see "grid description" above).  nlon */
/*     is also the first dimension of array goff and greg.  nlon determines */
/*     the grid increment in longitude as dlon = 2.*pi/nlon.  for example, */
/*     nlon = 144 for a 2.5 degree grid.  nlon can be even or odd and must */
/*     be greater than or equal to 4.  the efficiency of the computation */
/*     is improved when nlon is a product of small primes. */

/* ... nlat */

/*     the number of latitude points on the "offset" uniform grid.  nlat+1 */
/*     is the number of latitude points on the "regular" uniform grid (see */
/*     "grid description" above).  nlat is the second dimension of array goff. */
/*     nlat+1 must be the second dimension of the array greg in the program */
/*     calling sshifte.  nlat determines the grid in latitude as pi/nlat. */
/*     for example, nlat = 36 for a five degree grid.  nlat must be at least 3. */

/* ... goff */

/*     a nlon by nlat array that contains data on the offset grid */
/*     described above.  goff is a given input argument if ioff=0. */
/*     goff is a generated output argument if ioff=1. */

/* ... greg */

/*     a nlon by nlat+1 array that contains data on the regular grid */
/*     described above.  greg is a given input argument if ioff=1. */
/*     greg is a generated output argument if ioff=0. */

/* ... wsav */

/*     a doublereal saved work space array that must be initialized by calling */
/*     subroutine sshifti(ioff,nlon,nlat,wsav,ier) before calling sshifte. */
/*     wsav can then be used repeatedly by sshifte as long as ioff, nlon, */
/*     and nlat do not change.  this bypasses redundant computations and */
/*     saves time.  undetectable errors will result if sshifte is called */
/*     without initializing wsav whenever ioff, nlon, or nlat change. */

/* ... lsav */

/*     the length of the saved work space wsav in the routine calling sshifte */
/*     and sshifti.  lsave must be greater than or equal to 2*(2*nlat+nlon+16). */

/* ... work */

/*     a doublereal unsaved work space */

/* ... lwork */

/*     the length of the unsaved work space in the routine calling sshifte */
/*     lwork must be greater than or equal to 2*nlon*(nlat+1) if nlon is even. */
/*     lwork must be greater than or equal to nlon*(5*nlat+1) if nlon is odd. */

/* ... ier */

/*     indicates errors in input parameters */

/*     = 0 if no errors are detected */

/*     = 1 if ioff is not equal to 0 or 1 */

/*     = 1 if nlon < 4 */

/*     = 2 if nlat < 3 */

/*     = 3 if lsave < 2*(nlon+2*nlat+16) */

/*     = 4 if lwork < 2*nlon*(nlat+1) for nlon even or */
/*            lwork < nlon*(5*nlat+1) for nlon odd */

/* *** end of sshifte documentation */

/*     subroutine sshifti(ioff,nlon,nlat,lsav,wsav,ier) */

/*     subroutine sshifti initializes the saved work space wsav */
/*     for ioff and nlon and nlat (see documentation for sshifte). */
/*     sshifti must be called before sshifte whenever ioff or nlon */
/*     or nlat change. */

/* ... ier */

/*     = 0 if no errors with input arguments */

/*     = 1 if ioff is not 0 or 1 */

/*     = 2 if nlon < 4 */

/*     = 3 if nlat < 3 */

/*     = 4 if lsav < 2*(2*nlat+nlon+16) */

/* *** end of sshifti documentation */

/* Subroutine */ int sshifte_(integer *ioff, integer *nlon, integer *nlat, 
	doublereal *goff, doublereal *greg, doublereal *wsav, integer *lsav, doublereal *wrk, integer 
	*lwrk, integer *ier)
{
    /* System generated locals */
    integer goff_dim1, goff_offset, greg_dim1, greg_offset;

    /* Local variables */
    static integer i1, i2, n2, nr, nlat2, nlatp1;
    extern /* Subroutine */ int shftoff_(integer *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, doublereal *), shftreg_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *);


/*     check input parameters */

    /* Parameter adjustments */
    greg_dim1 = *nlon;
    greg_offset = 1 + greg_dim1;
    greg -= greg_offset;
    goff_dim1 = *nlon;
    goff_offset = 1 + goff_dim1;
    goff -= goff_offset;
    --wsav;
    --wrk;

    /* Function Body */
    *ier = 1;
    if (*ioff * (*ioff - 1) != 0) {
	return 0;
    }
    *ier = 2;
    if (*nlon < 4) {
	return 0;
    }
    *ier = 3;
    if (*nlat < 3) {
	return 0;
    }
    *ier = 4;
    if (*lsav < (*nlat << 1) + *nlon + 16 << 1) {
	return 0;
    }
    *ier = 5;
    n2 = (*nlon + 1) / 2;
    if (n2 << 1 == *nlon) {
	if (*lwrk < (*nlon << 1) * (*nlat + 1)) {
	    return 0;
	}
	i1 = 1;
	nr = n2;
    } else {
	if (*lwrk < *nlon * (*nlat * 5 + 1)) {
	    return 0;
	}
	i1 = (*nlat << 1) * *nlon + 1;
	nr = *nlon;
    }
    *ier = 0;
    nlat2 = *nlat + *nlat;
    i2 = i1 + (*nlat + 1) * *nlon;
    if (*ioff == 0) {
	shftoff_(nlon, nlat, &goff[goff_offset], &greg[greg_offset], &wsav[1],
		 &nr, &nlat2, &wrk[1], &wrk[i1], &wrk[i2]);
    } else {
	nlatp1 = *nlat + 1;
	shftreg_(nlon, nlat, &goff[goff_offset], &greg[greg_offset], &wsav[1],
		 &nr, &nlat2, &nlatp1, &wrk[1], &wrk[i1], &wrk[i2]);
    }
    return 0;
} /* sshifte_ */

/* Subroutine */ int shftoff_(integer *nlon, integer *nlat, doublereal *goff, doublereal *
	greg, doublereal *wsav, integer *nr, integer *nlat2, doublereal *rlat, doublereal *rlon,
	 doublereal *wrk)
{
    /* System generated locals */
    integer goff_dim1, goff_offset, greg_dim1, greg_offset, rlat_dim1, 
	    rlat_offset, rlon_dim1, rlon_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, n2, js, isav;
    extern /* Subroutine */ int shifth_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal gnorth, gsouth;


/*     shift offset grid to regular grid, i.e., */
/*     goff is given, greg is to be generated */

    /* Parameter adjustments */
    rlon_dim1 = *nlat;
    rlon_offset = 1 + rlon_dim1;
    rlon -= rlon_offset;
    greg_dim1 = *nlon;
    greg_offset = 1 + greg_dim1;
    greg -= greg_offset;
    goff_dim1 = *nlon;
    goff_offset = 1 + goff_dim1;
    goff -= goff_offset;
    --wsav;
    rlat_dim1 = *nr;
    rlat_offset = 1 + rlat_dim1;
    rlat -= rlat_offset;
    --wrk;

    /* Function Body */
    isav = (*nlat << 2) + 17;
    n2 = (*nlon + 1) / 2;

/*     execute full circle latitude shifts for nlon odd or even */

    if (n2 << 1 > *nlon) {

/*     odd number of longitudes */

	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		rlon[i__ + j * rlon_dim1] = goff[j + i__ * goff_dim1];
	    }
	}

/*       half shift in longitude */

	shifth_(nlat, nlon, &rlon[rlon_offset], &wsav[isav], &wrk[1]);

/*       set full 2*nlat circles in rlat using shifted values in rlon */

	i__1 = n2 - 1;
	for (j = 1; j <= i__1; ++j) {
	    js = j + n2;
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		rlat[j + i__ * rlat_dim1] = goff[j + i__ * goff_dim1];
		rlat[j + (*nlat + i__) * rlat_dim1] = rlon[*nlat + 1 - i__ + 
			js * rlon_dim1];
	    }
	}
	i__1 = *nlon;
	for (j = n2; j <= i__1; ++j) {
	    js = j - n2 + 1;
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		rlat[j + i__ * rlat_dim1] = goff[j + i__ * goff_dim1];
		rlat[j + (*nlat + i__) * rlat_dim1] = rlon[*nlat + 1 - i__ + 
			js * rlon_dim1];
	    }
	}

/*       shift the nlon rlat vectors one half latitude grid */

	shifth_(nlon, nlat2, &rlat[rlat_offset], &wsav[1], &wrk[1]);

/*       set nonpole values in greg and average for poles */

	gnorth = 0.;
	gsouth = 0.;
	i__1 = *nlon;
	for (j = 1; j <= i__1; ++j) {
	    gnorth += rlat[j + rlat_dim1];
	    gsouth += rlat[j + (*nlat + 1) * rlat_dim1];
	    i__2 = *nlat;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		greg[j + i__ * greg_dim1] = rlat[j + i__ * rlat_dim1];
	    }
	}
	gnorth /= *nlon;
	gsouth /= *nlon;
    } else {

/*     even number of longitudes (no initial longitude shift necessary) */
/*     set full 2*nlat circles (over poles) for each longitude pair (j,js) */

	i__1 = n2;
	for (j = 1; j <= i__1; ++j) {
	    js = n2 + j;
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		rlat[j + i__ * rlat_dim1] = goff[j + i__ * goff_dim1];
		rlat[j + (*nlat + i__) * rlat_dim1] = goff[js + (*nlat + 1 - 
			i__) * goff_dim1];
	    }
	}

/*       shift the n2=(nlon+1)/2 rlat vectors one half latitude grid */

	shifth_(&n2, nlat2, &rlat[rlat_offset], &wsav[1], &wrk[1]);

/*       set nonpole values in greg and average poles */

	gnorth = 0.;
	gsouth = 0.;
	i__1 = n2;
	for (j = 1; j <= i__1; ++j) {
	    js = n2 + j;
	    gnorth += rlat[j + rlat_dim1];
	    gsouth += rlat[j + (*nlat + 1) * rlat_dim1];
	    i__2 = *nlat;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		greg[j + i__ * greg_dim1] = rlat[j + i__ * rlat_dim1];
		greg[js + i__ * greg_dim1] = rlat[j + (*nlat2 - i__ + 2) * 
			rlat_dim1];
	    }
	}
	gnorth /= n2;
	gsouth /= n2;
    }

/*     set poles */

    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	greg[j + greg_dim1] = gnorth;
	greg[j + (*nlat + 1) * greg_dim1] = gsouth;
    }

/*     execute full circle longitude shift */

    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    rlon[i__ + j * rlon_dim1] = greg[j + i__ * greg_dim1];
	}
    }
    shifth_(nlat, nlon, &rlon[rlon_offset], &wsav[isav], &wrk[1]);
    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    greg[j + i__ * greg_dim1] = rlon[i__ + j * rlon_dim1];
	}
    }
    return 0;
} /* shftoff_ */

/* Subroutine */ int shftreg_(integer *nlon, integer *nlat, doublereal *goff, doublereal *
	greg, doublereal *wsav, integer *nr, integer *nlat2, integer *nlatp1, doublereal *
	rlat, doublereal *rlon, doublereal *wrk)
{
    /* System generated locals */
    integer goff_dim1, goff_offset, greg_dim1, greg_offset, rlat_dim1, 
	    rlat_offset, rlon_dim1, rlon_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, n2, js, isav;
    extern /* Subroutine */ int shifth_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *);


/*     shift regular grid to offset grid, i.e., */
/*     greg is given, goff is to be generated */

    /* Parameter adjustments */
    goff_dim1 = *nlon;
    goff_offset = 1 + goff_dim1;
    goff -= goff_offset;
    --wsav;
    rlat_dim1 = *nr;
    rlat_offset = 1 + rlat_dim1;
    rlat -= rlat_offset;
    rlon_dim1 = *nlatp1;
    rlon_offset = 1 + rlon_dim1;
    rlon -= rlon_offset;
    greg_dim1 = *nlon;
    greg_offset = 1 + greg_dim1;
    greg -= greg_offset;
    --wrk;

    /* Function Body */
    isav = (*nlat << 2) + 17;
    n2 = (*nlon + 1) / 2;

/*     execute full circle latitude shifts for nlon odd or even */

    if (n2 << 1 > *nlon) {

/*     odd number of longitudes */

	i__1 = *nlat + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		rlon[i__ + j * rlon_dim1] = greg[j + i__ * greg_dim1];
	    }
	}

/*       half shift in longitude in rlon */

	i__1 = *nlat + 1;
	shifth_(&i__1, nlon, &rlon[rlon_offset], &wsav[isav], &wrk[1]);

/*       set full 2*nlat circles in rlat using shifted values in rlon */

	i__1 = n2;
	for (j = 1; j <= i__1; ++j) {
	    js = j + n2 - 1;
	    rlat[j + rlat_dim1] = greg[j + greg_dim1];
	    i__2 = *nlat;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		rlat[j + i__ * rlat_dim1] = greg[j + i__ * greg_dim1];
		rlat[j + (*nlat + i__) * rlat_dim1] = rlon[*nlat + 2 - i__ + 
			js * rlon_dim1];
	    }
	    rlat[j + (*nlat + 1) * rlat_dim1] = greg[j + (*nlat + 1) * 
		    greg_dim1];
	}
	i__1 = *nlon;
	for (j = n2 + 1; j <= i__1; ++j) {
	    js = j - n2;
	    rlat[j + rlat_dim1] = greg[j + greg_dim1];
	    i__2 = *nlat;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		rlat[j + i__ * rlat_dim1] = greg[j + i__ * greg_dim1];
		rlat[j + (*nlat + i__) * rlat_dim1] = rlon[*nlat + 2 - i__ + 
			js * rlon_dim1];
	    }
	    rlat[j + (*nlat + 1) * rlat_dim1] = greg[j + (*nlat + 1) * 
		    greg_dim1];
	}

/*       shift the nlon rlat vectors one halflatitude grid */

	shifth_(nlon, nlat2, &rlat[rlat_offset], &wsav[1], &wrk[1]);

/*       set values in goff */

	i__1 = *nlon;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		goff[j + i__ * goff_dim1] = rlat[j + i__ * rlat_dim1];
	    }
	}
    } else {

/*     even number of longitudes (no initial longitude shift necessary) */
/*     set full 2*nlat circles (over poles) for each longitude pair (j,js) */

	i__1 = n2;
	for (j = 1; j <= i__1; ++j) {
	    js = n2 + j;
	    rlat[j + rlat_dim1] = greg[j + greg_dim1];
	    i__2 = *nlat;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		rlat[j + i__ * rlat_dim1] = greg[j + i__ * greg_dim1];
		rlat[j + (*nlat + i__) * rlat_dim1] = greg[js + (*nlat + 2 - 
			i__) * greg_dim1];
	    }
	    rlat[j + (*nlat + 1) * rlat_dim1] = greg[j + (*nlat + 1) * 
		    greg_dim1];
	}

/*       shift the n2=(nlon+1)/2 rlat vectors one half latitude grid */

	shifth_(&n2, nlat2, &rlat[rlat_offset], &wsav[1], &wrk[1]);

/*       set values in goff */

	i__1 = n2;
	for (j = 1; j <= i__1; ++j) {
	    js = n2 + j;
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		goff[j + i__ * goff_dim1] = rlat[j + i__ * rlat_dim1];
		goff[js + i__ * goff_dim1] = rlat[j + (*nlat2 + 1 - i__) * 
			rlat_dim1];
	    }
	}
    }

/*     execute full circle longitude shift for all latitude circles */

    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    rlon[i__ + j * rlon_dim1] = goff[j + i__ * goff_dim1];
	}
    }
    i__1 = *nlat + 1;
    shifth_(&i__1, nlon, &rlon[rlon_offset], &wsav[isav], &wrk[1]);
    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    goff[j + i__ * goff_dim1] = rlon[i__ + j * rlon_dim1];
	}
    }
    return 0;
} /* shftreg_ */

/* Subroutine */ int sshifti_(integer *ioff, integer *nlon, integer *nlat, 
	integer *lsav, doublereal *wsav, integer *ier)
{
    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static doublereal dp, pi, dlat, dlon;
    static integer isav, nlat2;
    extern /* Subroutine */ int shifthi_(integer *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --wsav;

    /* Function Body */
    *ier = 1;
    if (*ioff * (*ioff - 1) != 0) {
	return 0;
    }
    *ier = 2;
    if (*nlon < 4) {
	return 0;
    }
    *ier = 3;
    if (*nlat < 3) {
	return 0;
    }
    *ier = 4;
    if (*lsav < (*nlat << 1) + *nlon + 16 << 1) {
	return 0;
    }
    *ier = 0;
    pi = atan(1.) * 4.;

/*     set lat,long increments */

    dlat = pi / *nlat;
    dlon = (pi + pi) / *nlon;

/*     initialize wsav for left or right latitude shifts */

    if (*ioff == 0) {
	dp = dlat * -.5;
    } else {
	dp = dlat * .5;
    }
    nlat2 = *nlat + *nlat;
    shifthi_(&nlat2, &dp, &wsav[1]);

/*     initialize wsav for left or right longitude shifts */

    if (*ioff == 0) {
	dp = dlon * -.5;
    } else {
	dp = dlon * .5;
    }
    isav = (*nlat << 2) + 17;
    shifthi_(nlon, &dp, &wsav[isav]);
    return 0;
} /* sshifti_ */

/* Subroutine */ int shifth_(integer *m, integer *n, doublereal *r__, doublereal *wsav, 
	doublereal *work)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    static integer k, l, n2;
    static doublereal r2km1, r2km2;
    extern /* Subroutine */ int hrfftb_(integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *), hrfftf_(integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *);

    /* Parameter adjustments */
    r_dim1 = *m;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --wsav;
    --work;

    /* Function Body */
    n2 = (*n + 1) / 2;

/*     compute fourier coefficients for r on shifted grid */

    hrfftf_(m, n, &r__[r_offset], m, &wsav[*n + 2], &work[1]);
    i__1 = *m;
    for (l = 1; l <= i__1; ++l) {
	i__2 = n2;
	for (k = 2; k <= i__2; ++k) {
	    r2km2 = r__[l + (k + k - 2) * r_dim1];
	    r2km1 = r__[l + (k + k - 1) * r_dim1];
	    r__[l + (k + k - 2) * r_dim1] = r2km2 * wsav[n2 + k] - r2km1 * 
		    wsav[k];
	    r__[l + (k + k - 1) * r_dim1] = r2km2 * wsav[k] + r2km1 * wsav[n2 
		    + k];
	}
    }

/*     shift r with fourier synthesis and normalization */

    hrfftb_(m, n, &r__[r_offset], m, &wsav[*n + 2], &work[1]);
    i__1 = *m;
    for (l = 1; l <= i__1; ++l) {
	i__2 = *n;
	for (k = 1; k <= i__2; ++k) {
	    r__[l + k * r_dim1] /= *n;
	}
    }
    return 0;
} /* shifth_ */

/* Subroutine */ int shifthi_(integer *n, doublereal *dp, doublereal *wsav)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static integer k, n2;
    extern /* Subroutine */ int hrffti_(integer *, doublereal *);


/*     initialize wsav for subroutine shifth */

    /* Parameter adjustments */
    --wsav;

    /* Function Body */
    n2 = (*n + 1) / 2;
    i__1 = n2;
    for (k = 2; k <= i__1; ++k) {
	wsav[k] = sin((k - 1) * *dp);
	wsav[k + n2] = cos((k - 1) * *dp);
    }
    hrffti_(n, &wsav[*n + 2]);
    return 0;
} /* shifthi_ */

