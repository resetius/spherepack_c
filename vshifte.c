/* vshifte.f -- translated by f2c (version 20061008).
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



/* ... file vshifte.f contains code and documentation for subroutine vshifte */
/*     and its initialization subroutine vshifti */

/* ... required files */

/*     hrfft.f */

/*     subroutine vshifte(ioff,nlon,nlat,uoff,voff,ureg,vreg, */
/*    +                   wsave,lsave,work,lwork,ierror) */

/* *** purpose */

/*     subroutine vshifte does a highly accurate 1/2 grid increment shift */
/*     in both longitude and latitude of equally spaced vector data on the */
/*     sphere. data is transferred between the nlon by nlat "offset grid" */
/*     in (uoff,voff) (which excludes poles) and the nlon by nlat+1 "regular */
/*     grid" in (ureg,vreg) (which includes poles).  the transfer can go from */
/*     (uoff,voff) to (ureg,vreg) or vice versa (see ioff).  the grids which */
/*     underly the vector fields are described below.  the north and south */
/*     pole are at 0.5*pi and-0.5*pi radians respectively (pi=4.*atan(1.)). */
/*     uoff and ureg are the east longitudinal vector data components.  voff */
/*     and vreg are the latitudinal vector data components. */

/*     subroutine sshifte can be used to shift scalar data on the sphere. */
/*     notice that scalar and vector quantities are fundamentally different */
/*     on the sphere.  for example, vectors are discontinuous and multiple */
/*     valued at the poles.  scalars are continuous and single valued at the */
/*     poles. erroneous results would be produced if one attempted to shift */
/*     vector fields with subroutine sshifte applied to each component of */
/*     of the vector. */

/* *** grid descriptions */

/*     let dlon = (pi+pi)/nlon and dlat = pi/nlat be the uniform grid */
/*     increments in longitude and latitude */

/*     offset grid */

/*     the "1/2 increment offset" grid (long(j),lat(i)) on which uoff(j,i) */
/*     and voff(j,i) are given (ioff=0) or generated (ioff=1) is */

/*          long(j) =0.5*dlon + (j-1)*dlon  (j=1,...,nlon) */

/*     and */

/*          lat(i) = -0.5*pi + 0.5*dlat + (i-1)*dlat (i=1,...,nlat) */

/*     the data in (uoff,voff) is "shifted" one half a grid increment in both */
/*     longitude and latitude and excludes the poles.  each uoff(j,1),voff(j,1) */
/*     is given at latitude -pi/2+dlat/2.  uoff(j,nlat),voff(j,nlat) is */
/*     given at pi/2-dlat/2 (1/2 a grid increment away from the poles). */
/*     uoff(1,i),voff(1,i) is given at longitude dlon/2.  each uoff(nlon,i), */
/*     voff(nlon,i) is given at longitude 2*pi-dlon/2. */

/*     regular grid */

/*     let dlat,dlon be as above.  then the nlon by nlat+1 grid on which */
/*     ureg(j,i),vreg(j,i) are generated (ioff=0) or given (ioff=1) is */

/*          lone(j) = (j-1)*dlon (j=1,...,nlon) */

/*      and */

/*          late(i) = -0.5*pi + (i-1)*dlat (i=1,...,nlat+1) */

/*     values in ureg,vreg include the poles and start at zero degrees */
/*     longitude and at the south pole this is the "usual" equally spaced */
/*     grid in geophysical coordinates. */

/* *** remark */

/*     subroutine vshifte can be used in conjunction with subroutine trvsph */
/*     when transferring vector data from an equally spaced "1/2 increment */
/*     offset" grid to a gaussian or equally spaced grid (which includes poles) */
/*     of any resolution.  this problem (personal communication with dennis */
/*     shea) is encountered in geophysical modeling and data analysis. */

/* *** method */

/*     fast fourier transform software from spherepack2 and trigonometric */
/*     identities are used to accurately "shift" periodic vectors half a */
/*     grid increment in latitude and longitude.  latitudinal shifts are */
/*     accomplished by setting periodic 2*nlat vectors over the pole for each */
/*     longitude.  vector values must be negated on one side of the pole */
/*     to maintain periodicity prior to the 2*nlat shift over the poles. */
/*     when nlon is odd, the 2*nlat latitudinal shift requires an additional */
/*     longitude shift to obtain symmetry necessary for full circle shifts */
/*     over the poles.  finally longitudinal shifts are executed for each */
/*     shifted latitude. */

/* *** argument description */

/* ... ioff */

/*     ioff = 0 if values on the offset grid in (uoff,voff) are given and */
/*              values on the regular grid in (ureg,vreg) are to be generated. */

/*     ioff = 1 if values on the regular grid in (ureg,vreg) are given and */
/*              values on the offset grid in (uoff,voff) are to be generated. */

/* ... nlon */

/*     the number of longitude points on both the "offset" and "regular" */
/*     uniform grid in longitude (see "grid description" above).  nlon */
/*     is also the first dimension of uoff,voff,ureg,vreg.  nlon determines */
/*     the grid increment in longitude as dlon = 2.*pi/nlon.  for example, */
/*     nlon = 144 for a 2.5 degree grid.  nlon can be even or odd and must */
/*     be greater than or equal to 4.  the efficiency of the computation */
/*     is improved when nlon is a product of small primes. */

/* ... nlat */

/*     the number of latitude points on the "offset" uniform grid.  nlat+1 */
/*     is the number of latitude points on the "regular" uniform grid (see */
/*     "grid description" above).  nlat is the second dimension of uoff,voff. */
/*     nlat+1 must be the second dimension of ureg,vreg in the program */
/*     calling vshifte.  nlat determines the grid in latitude as pi/nlat. */
/*     for example, nlat = 36 for a five degree grid.  nlat must be at least 3. */

/* ... uoff */

/*     a nlon by nlat array that contains the east longitudinal vector */
/*     data component on the offset grid described above.  uoff is a */
/*     given input argument if ioff=0.  uoff is a generated output */
/*     argument if ioff=1. */

/* ... voff */

/*     a nlon by nlat array that contains the latitudinal vector data */
/*     component on the offset grid described above.  voff is a given */
/*     input argument if ioff=0.  voff is a generated output argument */
/*     if ioff=1. */

/* ... ureg */

/*     a nlon by nlat+1 array that contains the east longitudinal vector */
/*     data component on the regular grid described above.  ureg is a given */
/*     input argument if ioff=1.  ureg is a generated output argument */
/*     if ioff=0. */

/* ... vreg */

/*     a nlon by nlat+1 array that contains the latitudinal vector data */
/*     component on the regular grid described above.  vreg is a given */
/*     input argument if ioff=1.  vreg is a generated output argument */
/*     if ioff=0. */

/* ... wsav */

/*     a real saved work space array that must be initialized by calling */
/*     subroutine vshifti(ioff,nlon,nlat,wsav,ier) before calling vshifte. */
/*     wsav can then be used repeatedly by vshifte as long as ioff, nlon, */
/*     and nlat do not change.  this bypasses redundant computations and */
/*     saves time.  undetectable errors will result if vshifte is called */
/*     without initializing wsav whenever ioff, nlon, or nlat change. */

/* ... lsav */

/*     the length of the saved work space wsav in the routine calling vshifte */
/*     and sshifti.  lsave must be greater than or equal to 2*(2*nlat+nlon+16). */

/* ... work */

/*     a real unsaved work space */

/* ... lwork */

/*     the length of the unsaved work space in the routine calling vshifte */
/*     if nlon is even then lwork must be greater than or equal to */

/*          2*nlon*(nlat+1) */

/*     if nlon is odd then lwork must be greater than or equal to */

/*          nlon*(5*nlat+1) */

/* ... ier */

/*     indicates errors in input parameters */

/*     = 0 if no errors are detected */

/*     = 1 if ioff is not equal to 0 or 1 */

/*     = 2 if nlon < 4 */

/*     = 3 if nlat < 3 */

/*     = 4 if lsave < 2*(nlon+2*nlat)+32 */

/*     = 5 if lwork < 2*nlon*(nlat+1) for nlon even or */
/*            lwork < nlon*(5*nlat+1) for nlon odd */

/* *** end of vshifte documentation */

/*     subroutine vshifti(ioff,nlon,nlat,lsav,wsav,ier) */

/*     subroutine vshifti initializes the saved work space wsav */
/*     for ioff and nlon and nlat (see documentation for vshifte). */
/*     vshifti must be called before vshifte whenever ioff or nlon */
/*     or nlat change. */

/* ... ier */

/*     = 0 if no errors with input arguments */

/*     = 1 if ioff is not 0 or 1 */

/*     = 2 if nlon < 4 */

/*     = 3 if nlat < 3 */

/*     = 4 if lsav < 2*(2*nlat+nlon+16) */

/* *** end of vshifti documentation */

/* Subroutine */ int vshifte_(integer *ioff, integer *nlon, integer *nlat, 
	doublereal *uoff, doublereal *voff, doublereal *ureg, doublereal *
	vreg, doublereal *wsav, integer *lsav, doublereal *wrk, integer *lwrk,
	 integer *ier)
{
    /* System generated locals */
    integer uoff_dim1, uoff_offset, voff_dim1, voff_offset, ureg_dim1, 
	    ureg_offset, vreg_dim1, vreg_offset;

    /* Local variables */
    static integer i1, i2, i3, n2, nr, nlat2, nlatp1;
    extern /* Subroutine */ int vhftoff_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), vhftreg_(
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);


/*     check input parameters */

    /* Parameter adjustments */
    vreg_dim1 = *nlon;
    vreg_offset = 1 + vreg_dim1;
    vreg -= vreg_offset;
    ureg_dim1 = *nlon;
    ureg_offset = 1 + ureg_dim1;
    ureg -= ureg_offset;
    voff_dim1 = *nlon;
    voff_offset = 1 + voff_dim1;
    voff -= voff_offset;
    uoff_dim1 = *nlon;
    uoff_offset = 1 + uoff_dim1;
    uoff -= uoff_offset;
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
    nlat2 = *nlat + *nlat;
    nlatp1 = *nlat + 1;
    n2 = (*nlon + 1) / 2;
    *ier = 5;
    if (n2 << 1 == *nlon) {
	if (*lwrk < (*nlon << 1) * (*nlat + 1)) {
	    return 0;
	}
	nr = n2;
	i1 = 1;
	i2 = 1;
	i3 = i2 + *nlon * nlatp1;
    } else {
	if (*lwrk < *nlon * (*nlat * 5 + 1)) {
	    return 0;
	}
	nr = *nlon;
	i1 = 1;
	i2 = i1 + nlat2 * *nlon;
	i3 = i2 + nlatp1 * *nlon;
    }
    *ier = 0;
    if (*ioff == 0) {

/*     shift (uoff,voff) to (ureg,vreg) */

	vhftoff_(nlon, nlat, &uoff[uoff_offset], &ureg[ureg_offset], &wsav[1],
		 &nr, &nlat2, &nlatp1, &wrk[i1], &wrk[i2], &wrk[i2], &wrk[i3])
		;
	vhftoff_(nlon, nlat, &voff[voff_offset], &vreg[vreg_offset], &wsav[1],
		 &nr, &nlat2, &nlatp1, &wrk[i1], &wrk[i2], &wrk[i2], &wrk[i3])
		;
    } else {

/*     shift (ureg,vreg) to (uoff,voff) */

	vhftreg_(nlon, nlat, &uoff[uoff_offset], &ureg[ureg_offset], &wsav[1],
		 &nr, &nlat2, &nlatp1, &wrk[i1], &wrk[i2], &wrk[i2], &wrk[i3])
		;
	vhftreg_(nlon, nlat, &voff[voff_offset], &vreg[vreg_offset], &wsav[1],
		 &nr, &nlat2, &nlatp1, &wrk[i1], &wrk[i2], &wrk[i2], &wrk[i3])
		;
    }
    return 0;
} /* vshifte_ */

/* Subroutine */ int vhftoff_(integer *nlon, integer *nlat, doublereal *uoff, 
	doublereal *ureg, doublereal *wsav, integer *nr, integer *nlat2, 
	integer *nlatp1, doublereal *rlatu, doublereal *rlonu, doublereal *
	rlou, doublereal *wrk)
{
    /* System generated locals */
    integer uoff_dim1, uoff_offset, ureg_dim1, ureg_offset, rlatu_dim1, 
	    rlatu_offset, rlonu_dim1, rlonu_offset, rlou_dim1, rlou_offset, 
	    i__1, i__2;

    /* Local variables */
    static integer i__, j, n2, js, isav;
    extern /* Subroutine */ int vhifth_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);


/*     generate ureg from uoff (a vector component!) */

    /* Parameter adjustments */
    rlou_dim1 = *nlat;
    rlou_offset = 1 + rlou_dim1;
    rlou -= rlou_offset;
    uoff_dim1 = *nlon;
    uoff_offset = 1 + uoff_dim1;
    uoff -= uoff_offset;
    --wsav;
    rlatu_dim1 = *nr;
    rlatu_offset = 1 + rlatu_dim1;
    rlatu -= rlatu_offset;
    rlonu_dim1 = *nlatp1;
    rlonu_offset = 1 + rlonu_dim1;
    rlonu -= rlonu_offset;
    ureg_dim1 = *nlon;
    ureg_offset = 1 + ureg_dim1;
    ureg -= ureg_offset;
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
		rlou[i__ + j * rlou_dim1] = uoff[j + i__ * uoff_dim1];
	    }
	}

/*       half shift in longitude */

	vhifth_(nlat, nlon, &rlou[rlou_offset], &wsav[isav], &wrk[1]);

/*       set full 2*nlat circles in rlatu using shifted values in rlonu */

	i__1 = n2 - 1;
	for (j = 1; j <= i__1; ++j) {
	    js = j + n2;
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		rlatu[j + i__ * rlatu_dim1] = uoff[j + i__ * uoff_dim1];
		rlatu[j + (*nlat + i__) * rlatu_dim1] = -rlou[*nlat + 1 - i__ 
			+ js * rlou_dim1];
	    }
	}
	i__1 = *nlon;
	for (j = n2; j <= i__1; ++j) {
	    js = j - n2 + 1;
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		rlatu[j + i__ * rlatu_dim1] = uoff[j + i__ * uoff_dim1];
		rlatu[j + (*nlat + i__) * rlatu_dim1] = -rlou[*nlat + 1 - i__ 
			+ js * rlou_dim1];
	    }
	}

/*       shift the nlon rlat vectors one half latitude grid */

	vhifth_(nlon, nlat2, &rlatu[rlatu_offset], &wsav[1], &wrk[1]);

/*       set in ureg */

	i__1 = *nlon;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *nlat + 1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		ureg[j + i__ * ureg_dim1] = rlatu[j + i__ * rlatu_dim1];
	    }
	}
    } else {

/*     even number of longitudes (no initial longitude shift necessary) */
/*     set full 2*nlat circles (over poles) for each longitude pair (j,js) */
/*     negating js vector side for periodicity */

	i__1 = n2;
	for (j = 1; j <= i__1; ++j) {
	    js = n2 + j;
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		rlatu[j + i__ * rlatu_dim1] = uoff[j + i__ * uoff_dim1];
		rlatu[j + (*nlat + i__) * rlatu_dim1] = -uoff[js + (*nlatp1 - 
			i__) * uoff_dim1];
	    }
	}

/*       shift the n2=(nlon+1)/2 rlat vectors one half latitude grid */

	vhifth_(&n2, nlat2, &rlatu[rlatu_offset], &wsav[1], &wrk[1]);

/*       set ureg,vreg shifted in latitude */

	i__1 = n2;
	for (j = 1; j <= i__1; ++j) {
	    js = n2 + j;
	    ureg[j + ureg_dim1] = rlatu[j + rlatu_dim1];
	    ureg[js + ureg_dim1] = -rlatu[j + rlatu_dim1];
	    i__2 = *nlatp1;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		ureg[j + i__ * ureg_dim1] = rlatu[j + i__ * rlatu_dim1];
		ureg[js + i__ * ureg_dim1] = -rlatu[j + (*nlat2 - i__ + 2) * 
			rlatu_dim1];
	    }
	}
    }

/*     execute full circle longitude shift */

    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlatp1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    rlonu[i__ + j * rlonu_dim1] = ureg[j + i__ * ureg_dim1];
	}
    }
    vhifth_(nlatp1, nlon, &rlonu[rlonu_offset], &wsav[isav], &wrk[1]);
    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlatp1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ureg[j + i__ * ureg_dim1] = rlonu[i__ + j * rlonu_dim1];
	}
    }
    return 0;
} /* vhftoff_ */

/* Subroutine */ int vhftreg_(integer *nlon, integer *nlat, doublereal *uoff, 
	doublereal *ureg, doublereal *wsav, integer *nr, integer *nlat2, 
	integer *nlatp1, doublereal *rlatu, doublereal *rlonu, doublereal *
	rlou, doublereal *wrk)
{
    /* System generated locals */
    integer uoff_dim1, uoff_offset, ureg_dim1, ureg_offset, rlatu_dim1, 
	    rlatu_offset, rlonu_dim1, rlonu_offset, rlou_dim1, rlou_offset, 
	    i__1, i__2;

    /* Local variables */
    static integer i__, j, n2, js, isav;
    extern /* Subroutine */ int vhifth_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);


/*     generate uoff vector component from ureg */

    /* Parameter adjustments */
    rlou_dim1 = *nlat;
    rlou_offset = 1 + rlou_dim1;
    rlou -= rlou_offset;
    uoff_dim1 = *nlon;
    uoff_offset = 1 + uoff_dim1;
    uoff -= uoff_offset;
    --wsav;
    rlatu_dim1 = *nr;
    rlatu_offset = 1 + rlatu_dim1;
    rlatu -= rlatu_offset;
    rlonu_dim1 = *nlatp1;
    rlonu_offset = 1 + rlonu_dim1;
    rlonu -= rlonu_offset;
    ureg_dim1 = *nlon;
    ureg_offset = 1 + ureg_dim1;
    ureg -= ureg_offset;
    --wrk;

    /* Function Body */
    isav = (*nlat << 2) + 17;
    n2 = (*nlon + 1) / 2;

/*     execute full circle latitude shifts for nlon odd or even */

    if (n2 << 1 > *nlon) {

/*     odd number of longitudes */

	i__1 = *nlatp1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		rlonu[i__ + j * rlonu_dim1] = ureg[j + i__ * ureg_dim1];
	    }
	}

/*       half shift in longitude in rlon */

	vhifth_(nlatp1, nlon, &rlonu[rlonu_offset], &wsav[isav], &wrk[1]);

/*       set full 2*nlat circles in rlat using shifted values in rlon */

	i__1 = n2;
	for (j = 1; j <= i__1; ++j) {
	    js = j + n2 - 1;
	    rlatu[j + rlatu_dim1] = ureg[j + ureg_dim1];
	    i__2 = *nlat;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		rlatu[j + i__ * rlatu_dim1] = ureg[j + i__ * ureg_dim1];
		rlatu[j + (*nlat + i__) * rlatu_dim1] = -rlonu[*nlat + 2 - 
			i__ + js * rlonu_dim1];
	    }
	    rlatu[j + (*nlat + 1) * rlatu_dim1] = ureg[j + (*nlat + 1) * 
		    ureg_dim1];
	}
	i__1 = *nlon;
	for (j = n2 + 1; j <= i__1; ++j) {
	    js = j - n2;
	    rlatu[j + rlatu_dim1] = ureg[j + ureg_dim1];
	    i__2 = *nlat;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		rlatu[j + i__ * rlatu_dim1] = ureg[j + i__ * ureg_dim1];
		rlatu[j + (*nlat + i__) * rlatu_dim1] = -rlonu[*nlat + 2 - 
			i__ + js * rlonu_dim1];
	    }
	    rlatu[j + (*nlat + 1) * rlatu_dim1] = ureg[j + (*nlat + 1) * 
		    ureg_dim1];
	}

/*       shift the nlon rlat vectors one halflatitude grid */

	vhifth_(nlon, nlat2, &rlatu[rlatu_offset], &wsav[1], &wrk[1]);

/*       set values in uoff */

	i__1 = *nlon;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		uoff[j + i__ * uoff_dim1] = rlatu[j + i__ * rlatu_dim1];
	    }
	}
    } else {

/*     even number of longitudes (no initial longitude shift necessary) */
/*     set full 2*nlat circles (over poles) for each longitude pair (j,js) */

	i__1 = n2;
	for (j = 1; j <= i__1; ++j) {
	    js = n2 + j;
	    rlatu[j + rlatu_dim1] = ureg[j + ureg_dim1];
	    i__2 = *nlat;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		rlatu[j + i__ * rlatu_dim1] = ureg[j + i__ * ureg_dim1];
		rlatu[j + (*nlat + i__) * rlatu_dim1] = -ureg[js + (*nlat + 2 
			- i__) * ureg_dim1];
	    }
	    rlatu[j + (*nlat + 1) * rlatu_dim1] = ureg[j + (*nlat + 1) * 
		    ureg_dim1];
	}

/*       shift the n2=(nlon+1)/2 rlat vectors one half latitude grid */

	vhifth_(&n2, nlat2, &rlatu[rlatu_offset], &wsav[1], &wrk[1]);

/*       set values in uoff */

	i__1 = n2;
	for (j = 1; j <= i__1; ++j) {
	    js = n2 + j;
	    i__2 = *nlat;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		uoff[j + i__ * uoff_dim1] = rlatu[j + i__ * rlatu_dim1];
		uoff[js + i__ * uoff_dim1] = -rlatu[j + (*nlat2 + 1 - i__) * 
			rlatu_dim1];
	    }
	}
    }

/*     execute full circle longitude shift for all latitude circles */

    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    rlou[i__ + j * rlou_dim1] = uoff[j + i__ * uoff_dim1];
	}
    }
    vhifth_(nlat, nlon, &rlou[rlou_offset], &wsav[isav], &wrk[1]);
    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    uoff[j + i__ * uoff_dim1] = rlou[i__ + j * rlou_dim1];
	}
    }
    return 0;
} /* vhftreg_ */

/* Subroutine */ int vshifti_(integer *ioff, integer *nlon, integer *nlat, 
	integer *lsav, doublereal *wsav, integer *ier)
{
    /* Builtin functions */
    double atan(doublereal);

    /* Local variables */
    static doublereal dp, pi, dlat, dlon;
    static integer isav, nlat2;
    extern /* Subroutine */ int vhifthi_(integer *, doublereal *, doublereal *
	    );


/*     initialize wsav for vshifte */

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

/*     set left or right latitude shifts */

    if (*ioff == 0) {
	dp = dlat * -.5;
    } else {
	dp = dlat * .5;
    }
    nlat2 = *nlat + *nlat;
    vhifthi_(&nlat2, &dp, &wsav[1]);

/*     set left or right longitude shifts */

    if (*ioff == 0) {
	dp = dlon * -.5;
    } else {
	dp = dlon * .5;
    }
    isav = (*nlat << 2) + 17;
    vhifthi_(nlon, &dp, &wsav[isav]);
    return 0;
} /* vshifti_ */

/* Subroutine */ int vhifth_(integer *m, integer *n, doublereal *r__, 
	doublereal *wsav, doublereal *work)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    static integer k, l, n2;
    static doublereal r2km1, r2km2;
    extern /* Subroutine */ int hrfftb_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *), hrfftf_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *);

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
} /* vhifth_ */

/* Subroutine */ int vhifthi_(integer *n, doublereal *dp, doublereal *wsav)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal);

    /* Local variables */
    static integer k, n2;
    extern /* Subroutine */ int hrffti_(integer *, doublereal *);


/*     initialize wsav for subroutine vhifth */

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
} /* vhifthi_ */

