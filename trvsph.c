/* trvsph.f -- translated by f2c (version 20061008).
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


/* ... file trvsph.f */

/*     this file contains documentation and code for subroutine trvsph */

/* ... required files */

/*     sphcom.f, hrfft.f, gaqd.f, vhaec.f, vhsec.f, vhagc.f, vhsgc.f */

/*     subroutine trvsph (intl,igrida,nlona,nlata,iveca,ua,va, */
/*    +igridb,nlonb,nlatb,ivecb,ub,vb,wsave,lsave,lsvmin,work, */
/*    +lwork,lwkmin,dwork,ldwork,ier) */

/* *** author */

/*     John C. Adams (NCAR 1997), email: johnad@ncar.ucar.edu */

/* *** purpose */

/*     subroutine trvsph transfers vector data given in (ua,va) on a grid on */
/*     the full sphere to vector data in (ub,vb) on a grid on the full sphere. */
/*     the grids on which (ua,va) is given and (ub,vb) is generated can be */
/*     specified independently of each other (see the input arguments igrida, */
/*     igridb,iveca,ivecb).  ua and ub are the east longitudinal components of */
/*     the given and transformed vector fields.  va is either the latitudinal */
/*     or colatitudinal component of the given vector field (see iveca). */
/*     vb is either the latitudinal or colatitudinal component of the */
/*     transformed vector field (see ivecb).  for transferring scalar data */
/*     on the sphere, use subroutine trssph. */

/* *   notice that scalar and vector quantities are fundamentally different */
/*     on the sphere.  for example, vectors are discontinuous and multiple */
/*     valued at the poles.  scalars are continuous and single valued at the */
/*     poles. erroneous results would be produced if one attempted to transfer */
/*     vector fields between grids with subroutine trssph applied to each */
/*     component of the vector. */

/* *** underlying grid assumptions and a description */

/*     discussions with the ncar scd data support group and others indicate */
/*     there is no standard grid for storing observational or model generated */
/*     data on the sphere.  subroutine trvsph was designed to handle most */
/*     cases likely to be encountered when moving data from one grid format */
/*     to another. */

/*     the grid on which (ua,va) is given must be equally spaced in longitude */
/*     and either equally spaced or gaussian in latitude (or colatitude). */
/*     longitude, which can be either the first or second dimension of ua,va */
/*     subdivides [0,2pi) excluding the periodic point 2pi.  (co)latitude, */
/*     which can be the second or first dimension of ua,va, has south */
/*     to north or north to south orientation with increasing subscript */
/*     value in ua,va (see the argument igrida). */

/*     the grid on which ub,vb is generated must be equally spaced in longitude */
/*     and either equally spaced or gaussian in latitude (or colatitude). */
/*     longitude, which can be either the first or second dimension of ub,vb */
/*     subdivides [0,2pi) excluding the periodic point 2pi.  (co)latitude, */
/*     which can be the second or first dimension of ub,vb, has south */
/*     to north or north to south orientation with increasing subscript */
/*     value in db (see the argument igridb). */

/*     let nlon be either nlona or nlonb (the number of grid points in */
/*     longitude.  the longitude grid subdivides [0,2pi) into nlon spaced */
/*     points */

/*          (j-1)*2.*pi/nlon  (j=1,...,nlon). */

/*     it is not necessary to communicate to subroutine trvsph whether the */
/*     underlying grids are in latitude or colatitude.  it is only necessary */
/*     to communicate whether they run south to north or north to south with */
/*     increasing subscripts.  a brief discussion of latitude and colatitude */
/*     follows.  equally spaced latitude grids are assumed to subdivide */
/*     [-pi/2,pi/2] with the south pole at -pi/2 and north pole at pi/2. */
/*     equally spaced colatitude grids subdivide [0,pi] with the north pole */
/*     at 0 and south pole at pi.  equally spaced partitions on the sphere */
/*     include both poles.  gaussian latitude grids subdivide (-pi/2,pi/2) */
/*     and gaussian colatitude grids subdivide (0,pi).  gaussian grids do not */
/*     include the poles.  the gaussian grid points are uniquely determined by */
/*     the size of the partition.  they can be computed in colatitude in */
/*     (0,pi) (north to south) in double precision by the spherepack subroutine */
/*     gaqd.  let nlat be nlata or nlatb if either the ua,va or ub,vb grid is */
/*     gaussian.  let */

/*        north pole                             south pole */
/*        ----------                             ---------- */
/*           0.0    <  cth(1) < ... < cth(nlat)  <   pi */


/*     be nlat gaussian colatitude points in the interval (0,pi) and let */

/*        south pole                        north pole */
/*        ----------                        ---------- */
/*           -pi/2  < th(1) < ... < th(nlat) < pi/2 */

/*     be nlat gaussian latitude points in the open interval (-pi/2,pi/2). */
/*     these are related by */

/*          th(i) = -pi/2 + cth(i)  (i=1,...,nlat) */

/*     if the (ua,va) or (ub,vb) grid is equally spaced in (co)latitude then */

/*          ctht(i) = (i-1)*pi/(nlat-1) */
/*                                               (i=1,...,nlat) */
/*          tht(i) = -pi/2 + (i-1)*pi/(nlat-1) */

/*     define the equally spaced (north to south) colatitude and (south to */
/*     north) latitude grids. */

/* *** method (simplified description) */

/*    (1) */

/*     the vector field (ua,va) is reformated to a vector field in mathematical */
/*     spherical coordinates using array transpositions, subscript reordering */
/*     and negation of va as necessary (see arguments igrida,iveca). */

/*     (2) */

/*     a vector harmonic analysis is performed on the result from (1) */

/*     (3) */

/*     a vector harmonic synthesis is performed on the (ub,vb) grid */
/*     using as many coefficients from (2) as possible (i.e., as */
/*     as is consistent with the size of the ub,vb grid). */

/*     (4) */

/*     the vector field generated in (3) is transformed from mathematical */
/*     spherical coordinates to the form flagged by ivecb and igridb in */
/*     (ub,vb) using array transpositions, subscript reordering and negation */
/*     as necessary */


/* *** advantages */

/*     the use of vector spherical harmonics to transfer vector data is */
/*     highly accurate and preserves properties of vectors on the sphere. */
/*     the method produces a weighted least squares fit to vector data in */
/*     which waves are resolved uniformly on the full sphere.  high frequencies */
/*     induced by closeness of grid points near the poles (due to computational */
/*     or observational errors) are smoothed.  the method is consistent with */
/*     methods used to generate vector data in numerical spectral models based */
/*     on spherical harmonics.  for more discussion of these and related issues, */
/*     see "on the spectral approximation of discrete scalar and vector */
/*     functions on the sphere," siam j. numer. anal., vol. 16, december 1979, */
/*     pp. 934-949, by paul swarztrauber. */


/* *** comment */

/*     on a nlon by nlat or nlat by nlon grid (gaussian or equally spaced) */
/*     spherical harmonic analysis generates and synthesis utilizes */
/*     min0(nlat,(nlon+2)/2)) by nlat coefficients.  consequently, for */
/*     ua,va and ub,vb,  if either */

/*             min0(nlatb,(nlonb+2)/2) < min0(nlata,(nlona+2)/2) */

/*     or if */

/*             nlatb < nlata */

/*     then all the coefficients generated by an analysis of ua,va cannot be */
/*     used in the synthesis which generates ub,vb.  in this case "information" */
/*     can be lost in generating ub,vb.  more precisely, information will be */
/*     lost if the analysis of ua,va yields nonzero coefficients which are */
/*     outside the coefficient bounds determined by the ub,vb grid. still */
/*     transference with vector spherical harmonics will yield results */
/*     consistent with grid resolution and is highly accurate. */

/* *** input arguments */

/* ... intl */

/*     an initialization argument which should be zero on an initial call to */
/*     trvsph.  intl should be one if trvsph is being recalled and */

/*          igrida,nlona,nlata,iveca,igridb,nlonb,nlatb,ivecb */

/*     have not changed from the previous call.  if any of these arguments have */
/*     changed intl=0 must be used to avoid undetectable errors.  when allowed, */
/*     calls with intl=1 bypass redundant computation and save time.  it can */
/*     be used when transferring multiple vector data sets with the same */
/*     underlying grids. */

/* ... igrida */

/*     an integer vector dimensioned two which identifies the underlying grid */
/*     on the full sphere for the given vector data (ua,va) as follows: */

/*     igrida(1) */

/*     = -1 */
/*     if the latitude (or colatitude) grid for ua,va is an equally spaced */
/*     partition of [-pi/2,pi/2] ( or [0,pi] ) including the poles which */
/*     runs north to south with increasing subscript value */

/*     = +1 */
/*     if the latitude (or colatitude) grid for ua,va is an equally spaced */
/*     partition of [-pi/2,pi/2] ( or [0,pi] ) including the poles which */
/*     runs south to north with increasing subscript value */

/*     = -2 */
/*     if the latitude (or colatitude) grid for ua,va is a gaussian partition */
/*     of (-pi/2,pi/2) ( or (0,pi) ) excluding the poles which runs north */
/*     to south with increasing subscript value */

/*     = +2 */
/*     if the latitude (or colatitude) grid for ua,va is a gaussian partition */
/*     of (-pi/2,pi/2) ( or (0,pi) ) excluding the poles which runs south */
/*     north with increasing subscript value */

/*     igrida(2) */

/*     = 0 if the underlying grid for ua,va is a nlona by nlata */

/*     = 1 if the underlying grid for ua,va is a nlata by nlona */


/* ... nlona */

/*     the number of longitude points on the uniform grid which partitions */
/*     [0,2pi) for the given vector (ua,va).  nlona is also the first or second */
/*     dimension of ua,va (see igrida(2)) in the program which calls trvsph. */
/*     nlona determines the grid increment in longitude as 2*pi/nlona. for */
/*     example nlona = 72 for a five degree grid.  nlona must be greater than */
/*     or equal to 4.  the efficiency of the computation is improved when */
/*     nlona is a product of small prime numbers */

/* ... nlata */

/*     the number of points in the latitude (or colatitude) grid for the */
/*     given vector (ua,va).  nlata is also the first or second dimension */
/*     of ua and va (see igrida(2)) in the program which calls trvsph. */
/*     if nlata is odd then the equator will be located at the (nlata+1)/2 */
/*     gaussian grid point.  if nlata is even then the equator will be */
/*     located half way between the nlata/2 and nlata/2+1 grid points. */

/* ... iveca */

/*     if iveca=0 is input then va is the latitudinal component of the */
/*     given vector field. if iveca=1 then va is the colatitudinal */
/*     compoenent of the given vector field.  in either case, ua must */
/*     be the east longitudinal component of the given vector field. */

/* *** note: */
/*     igrida(1)=-1 or igrida(1)=-2, igrida(2)=1, and iveca=1 corresponds */
/*     to the "usual" mathematical spherical coordinate system required */
/*     by most of the drivers in spherepack2.  igrida(1)=1 or igrida(1)=2, */
/*     igrida(2)=0, and iveca=0 corresponds to the "usual" geophysical */
/*     spherical coordinate system. */


/* ... ua */

/*     ua is the east longitudinal component of the given vector field. */
/*     ua must be dimensioned nlona by nlata in the program calling trvsph if */
/*     igrida(2) = 0.  ua must be dimensioned nlata by nlona in the program */
/*     calling trvsph if igrida(2) = 1.  if ua is not properly dimensioned */
/*     and if the latitude (colatitude) values do not run south to north or */
/*     north to south as flagged by igrida(1) (this cannot be checked!) then */
/*     incorrect results will be produced. */


/* ... va */

/*     va is either the latitudinal or colatitudinal componenet of the */
/*     given vector field (see iveca).  va must be dimensioned nlona by */
/*     nlata in the program calling trvsph if igrida(2)=0.  va must be */
/*     dimensioned nlata by nlona in the program calling trvsph if */
/*     igrida(2)=1.  if va is not properly dimensioned or if the latitude */
/*     (colatitude) values do not run south to north or north to south */
/*     as flagged by igrida(1) (this cannot be checked!) then incorrect */
/*     results will be produced. */

/* ... igridb */

/*     an integer vector dimensioned two which identifies the underlying grid */
/*     on the full sphere for the transformed vector (ub,vb) as follows: */

/*     igridb(1) */

/*     = -1 */
/*     if the latitude (or colatitude) grid for ub,vb is an equally spaced */
/*     partition of [-pi/2,pi/2] ( or [0,pi] ) including the poles which */
/*     north to south */

/*     = +1 */
/*     if the latitude (or colatitude) grid for ub,vb is an equally spaced */
/*     partition of [-pi/2,pi/2] ( or [0,pi] ) including the poles which */
/*     south to north */

/*     = -2 */
/*     if the latitude (or colatitude) grid for ub,vb is a gaussian partition */
/*     of (-pi/2,pi/2) ( or (0,pi) ) excluding the poles which runs north to */
/*     south */

/*     = +2 */
/*     if the latitude (or colatitude) grid for ub,vb is a gaussian partition */
/*     of (-pi/2,pi/2) ( or (0,pi) ) excluding the poles which runs south to */
/*     north */

/*     igridb(2) */

/*     = 0 if the underlying grid for ub,vb is a nlonb by nlatb */

/*     = 1 if the underlying grid for ub,vb is a nlatb by nlonb */


/* ... nlonb */

/*     the number of longitude points on the uniform grid which partitions */
/*     [0,2pi) for the transformed vector (ub,vb).  nlonb is also the first or */
/*     second dimension of ub and vb (see igridb(2)) in the program which calls */
/*     trvsph.  nlonb determines the grid increment in longitude as 2*pi/nlonb. */
/*     for example nlonb = 72 for a five degree grid.  nlonb must be greater */
/*     than or equal to 4.  the efficiency of the computation is improved when */
/*     nlonb is a product of small prime numbers */

/* ... nlatb */

/*     the number of points in the latitude (or colatitude) grid for the */
/*     transformed vector (ub,vb).  nlatb is also the first or second dimension */
/*     of ub and vb (see igridb(2)) in the program which calls trvsph. */
/*     if nlatb is odd then the equator will be located at the (nlatb+1)/2 */
/*     gaussian grid point.  if nlatb is even then the equator will be */
/*     located half way between the nlatb/2 and nlatb/2+1 grid points. */

/* ... ivecb */

/*     if ivecb=0 is input then vb is the latitudinal component of the */
/*     given vector field. if ivecb=1 then vb is the colatitudinal */
/*     compoenent of the given vector field.  in either case, ub must */
/*     be the east longitudinal component of the given vector field. */

/* *** note: */
/*     igridb(1)=-1 or igridb(1)=-2, igridb(2)=1, and ivecb=1 corresponds */
/*     to the "usual" mathematical spherical coordinate system required */
/*     by most of the drivers in spherepack2.  igridb(1)=1 or igridb(1)=2, */
/*     igridb(2)=0, and ivecb=0 corresponds to the "usual" geophysical */
/*     spherical coordinate system. */

/* ... wsave */

/*     a saved work space array that can be utilized repeatedly by trvsph */
/*     as long as the arguments nlata,nlona,nlatb,nlonb remain unchanged. */
/*     wsave is set by a intl=0 call to trvsph.  wsave must not be altered */
/*     when trvsph is being recalled with intl=1. */

/* ... lsave */

/*     the dimension of the work space wsave as it appears in the program */
/*     that calls trvsph.  the minimum required value of lsave for the */
/*     current set of input arguments is set in the output argument lsvmin. */
/*     it can be determined by calling trvsph with lsave=0 and printing lsvmin. */

/*          la1 = min0(nlata,(nlona+1)/2), la2 = (nlata+1)/2 */

/*          lb1 = min0(nlatb,(nlonb+1)/2), lb2 = (nlatb+1)/2 */

/*          lwa = 4*nlata*la2+3*max0(la1-2,0)*(2*nlata-la1-1)+la2+nlona+15 */

/*          lwb = 4*nlatb*lb2+3*max0(lb1-2,0)*(2*nlatb-lb1-1)+nlonb+15 */

/*      then */

/*          lsvmin = lwa + lwb */

/*      is the minimal required work space length of wsave */


/* ... work */

/*     a work array that does not have to be preserved */

/* ... lwork */

/*     the dimension of the array work as it appears in the program that */
/*     calls trvsph. the minimum required value of lwork for the current */
/*     set of input arguments is set in the output argument lwkmin. */
/*     it can be determined by calling trvsph with lwork=0 and printing */
/*     lwkmin.  an estimate for lwork follows.  let nlat = max0(nlata,nlatb), */
/*     nlon = max0(nlona,nlonb) and l1 = min0(nlat,(nlon+2)/2).  with these */
/*     these definitions, the quantity */

/*            2*nlat*(8*l1 + 4*nlon + 3) */

/*     will suffice as a length for the unsaved work space.  this formula */
/*     may overestimate the required minimum value for lwork.  the exact */
/*     minimum value can be predetermined by calling trvsph wtih lwork=0 */
/*     and printout of lwkmin. */

/* ... dwork */

/*     a double precision work array that does not have to be preserved. */

/* ... ldwork */

/*     the length of dwork in the routine calling trvsph */
/*     Let */

/*       nlat = max0(nlata,nlatb) */

/*     ldwork must be at least 2*nlat*(nlat+1)+1 */


/* *** output arguments */


/* ... ub */

/*     a two dimensional array that contains the east longitudinal component */
/*     of the transformed vector data.  ub */
/*     must be dimensioned nlonb by nlatb in the program calling trvsph if */
/*     igridb(2)=0.  ub must be dimensioned nlatb by nlonb in the program */
/*     calling trvsph if igridb(2)=1.  if ub is not properly dimensioned */
/*     and if the latitude (colatitude) values do not run south to north or */
/*     north to south as flagged by igrdb(1) (this cannot be checked!) then */
/*     incorrect results will be produced. */


/* ... vb */

/*     a two dimensional array that contains the latitudinal or colatitudinal */
/*     component of the transformed vector data (see ivecb). */
/*     vb must be dimensioned nlonb by nlatb in the program calling trvsph if */
/*     igridb(2)=0.  vb must be dimensioned nlatb by nlonb in the program */
/*     calling trvsph if igridb(2)=1.  if vb is not properly dimensioned */
/*     and if the latitude (colatitude) values do not run south to north or */
/*     north to south as flagged by igrdb(1) (this cannot be checked!) then */
/*     incorrect results will be produced. */

/* ... lsvmin */

/*     the minimum length of the saved work space in wsave. */
/*     lsvmin is computed even if lsave < lsvmin (ier = 10). */

/* ... lwkmin */

/*     the minimum length of the unsaved work space in work. */
/*     lwkmin is computed even if lwork < lwkmin (ier = 11). */


/* *** error argument */

/* ... ier = 0  if no errors are detected */

/*         = 1  if intl is not 0 or 1 */

/*         = 2  if igrida(1) is not -1 or +1 or -2 or +2 */

/*         = 3  if igrida(2) is not 0 or 1 */

/*         = 4  if nlona is less than 4 */

/*         = 5  if nlata is less than 3 */

/*         = 6  if iveca is not 0 or 1 */

/*         = 7  if igridb(1) is not -1 or +1 or -2 or +2 */

/*         = 8  if igridb(2) is not 0 or 1 */

/*         = 9  if nlonb is less than 4 */

/*         =10  if nlatb is less than 3 */

/*         =11  if ivecb is not 0 or 1 */

/*         =12  if there is insufficient saved work space (lsave < lsvmin) */

/*         =13  if there is insufficient unsaved work space (lwork < lwkmin) */

/*         =14  indicates failure in an eigenvalue routine which computes */
/*              gaussian weights and points */

/*         =15  if ldwork is too small (insufficient double precision */
/*              unsaved work space) */

/* ***************************************************** */
/* ***************************************************** */

/*     end of argument description ... code follows */

/* ***************************************************** */
/* ***************************************************** */

/* Subroutine */ int trvsph_(integer *intl, integer *igrida, integer *nlona, 
	integer *nlata, integer *iveca, doublereal *ua, doublereal *va, 
	integer *igridb, integer *nlonb, integer *nlatb, integer *ivecb, 
	doublereal *ub, doublereal *vb, doublereal *wsave, integer *lsave, 
	integer *lsvmin, doublereal *work, integer *lwork, integer *lwkmin, 
	doublereal *dwork, integer *ldwork, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer jb, ig, iw, nt, lw, la1, la2, lb1, lb2, lwa, lwb, lwk1, 
	    lwk2, iabi, iaci, ibbi, ibci, iabr, iacr, ibbr, ibcr, nlat;
    extern /* Subroutine */ int negv_(integer *, integer *, doublereal *);
    static integer ityp, igrda, igrdb;
    extern /* Subroutine */ int vhagc_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), vhaec_(integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *), 
	    vhsec_(integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), vhsgc_(integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *), trvab_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    ), vhaeci_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *), vhagci_(integer *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *), 
	    vhseci_(integer *, integer *, doublereal *, integer *, doublereal 
	    *, integer *, integer *), vhsgci_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *), 
	    covlat_(integer *, integer *, doublereal *), trvplat_(integer *, 
	    integer *, doublereal *, doublereal *);


/*     include a save statement to ensure local variables in trvsph, set during */
/*     an intl=0 call, are preserved if trvsph is recalled with intl=1 */


/*     check input arguments */

    /* Parameter adjustments */
    --dwork;
    --work;
    --wsave;
    --vb;
    --ub;
    --igridb;
    --va;
    --ua;
    --igrida;

    /* Function Body */
    *ier = 1;
    if (*intl * (*intl - 1) != 0) {
	return 0;
    }
    *ier = 2;
    ig = igrida[1];
    if ((ig - 1) * (ig + 1) * (ig - 2) * (ig + 2) != 0) {
	return 0;
    }
    *ier = 3;
    ig = igrida[2];
    if (ig * (ig - 1) != 0) {
	return 0;
    }
    *ier = 4;
    if (*nlona < 4) {
	return 0;
    }
    *ier = 5;
    if (*nlata < 3) {
	return 0;
    }
    *ier = 6;
    if (*iveca * (*iveca - 1) != 0) {
	return 0;
    }
    *ier = 7;
    ig = igridb[1];
    if ((ig - 1) * (ig + 1) * (ig - 2) * (ig + 2) != 0) {
	return 0;
    }
    *ier = 8;
    ig = igridb[2];
    if (ig * (ig - 1) != 0) {
	return 0;
    }
    *ier = 9;
    if (*nlonb < 4) {
	return 0;
    }
    *ier = 10;
    if (*nlatb < 3) {
	return 0;
    }
    *ier = 11;
    if (*ivecb * (*ivecb - 1) != 0) {
	return 0;
    }
    *ier = 0;
    igrda = abs(igrida[1]);
    igrdb = abs(igridb[1]);
    if (*intl == 0) {
/* Computing MIN */
	i__1 = *nlata, i__2 = (*nlona + 1) / 2;
	la1 = min(i__1,i__2);
	la2 = (*nlata + 1) / 2;
/* Computing MIN */
	i__1 = *nlatb, i__2 = (*nlonb + 1) / 2;
	lb1 = min(i__1,i__2);
	lb2 = (*nlatb + 1) / 2;

/*     saved space for analysis on a grid */

/* Computing MAX */
	i__1 = la1 - 2;
	lwa = (*nlata << 2) * la2 + max(i__1,0) * 3 * ((*nlata << 1) - la1 - 
		1) + la2 + *nlona + 15;

/*     set saved work space length for synthesis on b grid */

/* Computing MAX */
	i__1 = lb1 - 2;
	lwb = (*nlatb << 2) * lb2 + max(i__1,0) * 3 * ((*nlatb << 1) - lb1 - 
		1) + *nlonb + 15;

/*     set minimum required saved work space length */

	*lsvmin = lwa + lwb;

/*     set wsave pointer */

	jb = lwa + 1;

/*     set pointers for vector spherical harmonic coefs in work */

	iabr = 1;
	iabi = iabr + la1 * *nlata;
	iacr = iabi + la1 * *nlata;
	iaci = iacr + la1 * *nlata;
	ibbr = iaci + la1 * *nlata;
	ibbi = ibbr + lb1 * *nlatb;
	ibcr = ibbi + lb1 * *nlatb;
	ibci = ibcr + lb1 * *nlatb;

/*     set pointers for remaining work */

	iw = ibci + lb1 * *nlatb;

/*     set remaining work space length in lw */

	lw = *lwork - iw;

/*     compute unsaved space for analysis and synthesis */

/* Computing MAX */
	i__1 = la2 * 6;
	lwk1 = (*nlata << 1) * ((*nlona << 1) + max(i__1,*nlona));
/* Computing MAX */
	i__1 = lb2 * 6;
	lwk2 = (*nlatb << 1) * ((*nlonb << 1) + max(i__1,*nlonb));

/*     set minimum unsaved work space required by trvsph */

	*lwkmin = iw + max(lwk1,lwk2);

/*     set error flags if saved or unsaved work space is insufficient */

	*ier = 12;
	if (*lsave < *lsvmin) {
	    return 0;
	}
	*ier = 13;
	if (*lwork < *lwkmin) {
	    return 0;
	}
	*ier = 15;
	nlat = max(*nlata,*nlatb);
	if (*ldwork < (nlat << 1) * (nlat + 1) + 1) {
	    return 0;
	}
	*ier = 0;
	if (igrda == 1) {

/*     initialize wsave for equally spaced analysis */

	    vhaeci_(nlata, nlona, &wsave[1], &lwa, &dwork[1], ldwork, ier);
	} else {

/*     initialize wsave for gaussian analysis */

	    vhagci_(nlata, nlona, &wsave[1], &lwa, &dwork[1], ldwork, ier);
	    if (*ier != 0) {

/*     flag failure in spherepack gaussian software */

		*ier = 14;
		return 0;
	    }
	}
	if (igrdb == 2) {

/*     initialize wsave for gaussian synthesis */

	    vhsgci_(nlatb, nlonb, &wsave[jb], &lwb, &dwork[1], ldwork, ier);
	    if (*ier != 0) {

/*     flag failure in spherepack gaussian software */

		*ier = 14;
		return 0;
	    }
	} else {

/*     initialize wsave for equally spaced synthesis */

	    vhseci_(nlatb, nlonb, &wsave[jb], &lwb, &dwork[1], ldwork, ier);
	}

/*     end of initialization (intl=0) call */

    }

/*     convert the vector field (ua,va) to mathematical spherical coordinates */

    if (igrida[2] == 0) {
	trvplat_(nlona, nlata, &ua[1], &work[1]);
	trvplat_(nlona, nlata, &va[1], &work[1]);
    }
    if (igrida[1] > 0) {
	covlat_(nlata, nlona, &ua[1]);
	covlat_(nlata, nlona, &va[1]);
    }
    if (*iveca == 0) {
	negv_(nlata, nlona, &va[1]);
    }
    nt = 1;
    ityp = 0;

/*     analyze vector field */

    if (igrda == 2) {
	vhagc_(nlata, nlona, &ityp, &nt, &va[1], &ua[1], nlata, nlona, &work[
		iabr], &work[iabi], &work[iacr], &work[iaci], &la1, nlata, &
		wsave[1], &lwa, &work[iw], &lw, ier);
    } else {
	vhaec_(nlata, nlona, &ityp, &nt, &va[1], &ua[1], nlata, nlona, &work[
		iabr], &work[iabi], &work[iacr], &work[iaci], &la1, nlata, &
		wsave[1], &lwa, &work[iw], &lw, ier);
    }

/*     transfer a grid coefficients to b grid coefficients */

    trvab_(&la1, nlata, &work[iabr], &work[iabi], &work[iacr], &work[iaci], &
	    lb1, nlatb, &work[ibbr], &work[ibbi], &work[ibcr], &work[ibci]);

/*     synthesize on b grid */

    if (igrdb == 1) {
	vhsec_(nlatb, nlonb, &ityp, &nt, &vb[1], &ub[1], nlatb, nlonb, &work[
		ibbr], &work[ibbi], &work[ibcr], &work[ibci], &lb1, nlatb, &
		wsave[jb], &lwb, &work[iw], &lw, ier);
    } else {
	vhsgc_(nlatb, nlonb, &ityp, &nt, &vb[1], &ub[1], nlatb, nlonb, &work[
		ibbr], &work[ibbi], &work[ibcr], &work[ibci], &lb1, nlatb, &
		wsave[jb], &lwb, &work[iw], &lw, ier);
    }

/*     restore a grid and b grid vector fields (now in math coordinates) to */
/*     agree with grid flags in igrida,iveca,igridb,ivecb */

    if (*iveca == 0) {
	negv_(nlata, nlona, &va[1]);
    }
    if (*ivecb == 0) {
	negv_(nlatb, nlonb, &vb[1]);
    }
    if (igrida[1] > 0) {
	covlat_(nlata, nlona, &ua[1]);
	covlat_(nlata, nlona, &va[1]);
    }
    if (igridb[1] > 0) {
	covlat_(nlatb, nlonb, &ub[1]);
	covlat_(nlatb, nlonb, &vb[1]);
    }
    if (igrida[2] == 0) {
	trvplat_(nlata, nlona, &ua[1], &work[1]);
	trvplat_(nlata, nlona, &va[1], &work[1]);
    }
    if (igridb[2] == 0) {
	trvplat_(nlatb, nlonb, &ub[1], &work[1]);
	trvplat_(nlatb, nlonb, &vb[1], &work[1]);
    }
    return 0;
} /* trvsph_ */

/* Subroutine */ int negv_(integer *nlat, integer *nlon, doublereal *v)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;


/*     negate (co)latitudinal vector componenet */

    /* Parameter adjustments */
    v_dim1 = *nlat;
    v_offset = 1 + v_dim1;
    v -= v_offset;

    /* Function Body */
    i__1 = *nlon;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nlat;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    v[i__ + j * v_dim1] = -v[i__ + j * v_dim1];
	}
    }
    return 0;
} /* negv_ */

/* Subroutine */ int trvab_(integer *ma, integer *na, doublereal *abr, 
	doublereal *abi, doublereal *acr, doublereal *aci, integer *mb, 
	integer *nb, doublereal *bbr, doublereal *bbi, doublereal *bcr, 
	doublereal *bci)
{
    /* System generated locals */
    integer abr_dim1, abr_offset, abi_dim1, abi_offset, acr_dim1, acr_offset, 
	    aci_dim1, aci_offset, bbr_dim1, bbr_offset, bbi_dim1, bbi_offset, 
	    bcr_dim1, bcr_offset, bci_dim1, bci_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, m, n;


/*     set coefficients for b grid from coefficients for a grid */

    /* Parameter adjustments */
    aci_dim1 = *ma;
    aci_offset = 1 + aci_dim1;
    aci -= aci_offset;
    acr_dim1 = *ma;
    acr_offset = 1 + acr_dim1;
    acr -= acr_offset;
    abi_dim1 = *ma;
    abi_offset = 1 + abi_dim1;
    abi -= abi_offset;
    abr_dim1 = *ma;
    abr_offset = 1 + abr_dim1;
    abr -= abr_offset;
    bci_dim1 = *mb;
    bci_offset = 1 + bci_dim1;
    bci -= bci_offset;
    bcr_dim1 = *mb;
    bcr_offset = 1 + bcr_dim1;
    bcr -= bcr_offset;
    bbi_dim1 = *mb;
    bbi_offset = 1 + bbi_dim1;
    bbi -= bbi_offset;
    bbr_dim1 = *mb;
    bbr_offset = 1 + bbr_dim1;
    bbr -= bbr_offset;

    /* Function Body */
    m = min(*ma,*mb);
    n = min(*na,*nb);
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    bbr[i__ + j * bbr_dim1] = abr[i__ + j * abr_dim1];
	    bbi[i__ + j * bbi_dim1] = abi[i__ + j * abi_dim1];
	    bcr[i__ + j * bcr_dim1] = acr[i__ + j * acr_dim1];
	    bci[i__ + j * bci_dim1] = aci[i__ + j * aci_dim1];
	}
    }

/*     set coefs outside triangle to zero */

    i__1 = *mb;
    for (i__ = m + 1; i__ <= i__1; ++i__) {
	i__2 = *nb;
	for (j = 1; j <= i__2; ++j) {
	    bbr[i__ + j * bbr_dim1] = 0.;
	    bbi[i__ + j * bbi_dim1] = 0.;
	    bcr[i__ + j * bcr_dim1] = 0.;
	    bci[i__ + j * bci_dim1] = 0.;
	}
    }
    i__1 = *nb;
    for (j = n + 1; j <= i__1; ++j) {
	i__2 = *mb;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    bbr[i__ + j * bbr_dim1] = 0.;
	    bbi[i__ + j * bbi_dim1] = 0.;
	    bcr[i__ + j * bcr_dim1] = 0.;
	    bci[i__ + j * bci_dim1] = 0.;
	}
    }
    return 0;
} /* trvab_ */

/* Subroutine */ int trvplat_(integer *n, integer *m, doublereal *data, 
	doublereal *work)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, j, ij, ji;


/*     transpose the n by m array data to a m by n array data */
/*     work must be at least n*m words long */

    /* Parameter adjustments */
    --work;
    --data;

    /* Function Body */
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ij = (j - 1) * *n + i__;
	    work[ij] = data[ij];
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    ji = (i__ - 1) * *m + j;
	    ij = (j - 1) * *n + i__;
	    data[ji] = work[ij];
	}
    }
    return 0;
} /* trvplat_ */

/* Subroutine */ int covlat_(integer *nlat, integer *nlon, doublereal *data)
{
    /* System generated locals */
    integer data_dim1, data_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, ib;
    doublereal temp;
    integer nlat2;


/*     reverse order of latitude (colatitude) grids */

    /* Parameter adjustments */
    data_dim1 = *nlat;
    data_offset = 1 + data_dim1;
    data -= data_offset;

    /* Function Body */
    nlat2 = *nlat / 2;
    i__1 = nlat2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ib = *nlat - i__ + 1;
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    temp = data[i__ + j * data_dim1];
	    data[i__ + j * data_dim1] = data[ib + j * data_dim1];
	    data[ib + j * data_dim1] = temp;
	}
    }
    return 0;
} /* covlat_ */

