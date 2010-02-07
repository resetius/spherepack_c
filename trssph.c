/* trssph.f -- translated by f2c (version 20061008).
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


/* ... file trssph.f */

/*     contains documentation and code for subroutine trssph */

/* ... required files */

/*     sphcom.f, hrfft.f, gaqd.f, shaec.f, shsec.f, shagc.f, shsgc.f */


/*     subroutine trssph(intl,igrida,nlona,nlata,da,igridb,nlonb,nlatb, */
/*    +db,wsave,lsave,lsvmin,work,lwork,lwkmin,dwork,ldwork,ier) */

/* *** purpose */

/*     subroutine trssph transfers data given in array da on a grid on the */
/*     full sphere to data in array db on a grid on the full sphere.  the */
/*     grids on which da is given and db is generated can be specified */
/*     independently of each other (see description below and the arguments */
/*     igrida,igridb).  for transferring vector data on the sphere, use */
/*     subroutine trvsph. */
/*     notice that scalar and vector quantities are fundamentally different */
/*     on the sphere.  for example, vectors are discontinuous and multiple */
/*     valued at the poles.  scalars are continuous and single valued at the */
/*     poles. erroneous results would be produced if one attempted to transfer */
/*     vector fields between grids with subroutine trssph applied to each */
/*     component of the vector. */


/* *** underlying grid assumptions and a description */

/*     discussions with the ncar scd data support group and others indicate */
/*     there is no standard grid for storing observational or model generated */
/*     data on the sphere.  subroutine trssph was designed to handle most */
/*     cases likely to be encountered when moving data from one grid format */
/*     to another. */

/*     the grid on which da is given must be equally spaced in longitude */
/*     and either equally spaced or gaussian in latitude (or colatitude). */
/*     longitude, which can be either the first or second dimension of da, */
/*     subdivides [0,2pi) excluding the periodic point 2pi.  (co)latitude, */
/*     which can be the second or first dimension of da, has south */
/*     to north or north to south orientation with increasing subscript */
/*     value in da (see the argument igrida). */

/*     the grid on which db is generated must be equally spaced in longitude */
/*     and either equally spaced or gaussian in latitude (or colatitude). */
/*     longitude, which can be either the first or second dimension of db, */
/*     subdivides [0,2pi) excluding the periodic point 2pi.  (co)latitude, */
/*     which can be the second or first dimension of db, has south */
/*     to north or north to south orientation with increasing subscript */
/*     value in db (see the argument igridb). */

/*     let nlon be either nlona or nlonb (the number of grid points in */
/*     longitude.  the longitude grid subdivides [0,2pi) into nlon spaced */
/*     points */

/*          (j-1)*2.*pi/nlon  (j=1,...,nlon). */

/*     it is not necessary to communicate to subroutine trssph whether the */
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
/*     gaqd.  let nlat be nlata or nlatb if either the da or db grid is */
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

/*     if the da or db grid is equally spaced in (co)latitude then */

/*          ctht(i) = (i-1)*pi/(nlat-1) */
/*                                               (i=1,...,nlat) */
/*          tht(i) = -pi/2 + (i-1)*pi/(nlat-1) */

/*     define the equally spaced (north to south) colatitude and (south to */
/*     north) latitude grids. */


/* *** method (simplified description) */

/*     for simplicity, assume da is a nlat by nlon data tabulation and da(i,j) */
/*     is the value at latitude theta(i) and longitude phi(j).  then */
/*     coefficients a(m,n) and b(m,n) can be determined so that da(i,j) is */
/*     approximated by the sum */

/*         l-1  n */
/*     (a) sum sum pbar(m,n,theta(i))*(a(m,n)*cos(m*phi(j)+b(m,n)*sin(m*phi(j)) */
/*         n=0 m=0 */

/*     here pbar(n,m,theta) are the normalized associated legendre functions */
/*     and l = min0(nlat,(nlon+2)/2).  the determination of a(m,n) and b(m,n) */
/*     is called spherical harmonic analysis. a sum of this form can then be */
/*     used to regenerate the data in db on the new grid with the known */
/*     a(m,n) and b(m,n).  this is referred to spherical harmonic synthesis. */
/*     analysis and synthesis subroutines from the software package spherepack, */
/*     are used for these purposes. */

/*     if da or db is not in mathematical spherical coordinates then array */
/*     transposition and/or subscript reordering is used prior to harmonic */
/*     analysis and after harmonic synthesis. */

/* *** advantages */

/*     the use of surface spherical harmonics to transfer spherical grid data */
/*     has advantages over pointwise grid interpolation schemes on the sphere. */
/*     it is highly accurate.  if p(x,y,z) is any polynomial of degree n or */
/*     less in x,y,z cartesian coordinates which is restricted to the surface */
/*     of the sphere, then p is exactly represented by sums of the form (a) */
/*     whenever n = mino(nlat,nlon/2) (i.e., transfers with spherical harmonics */
/*     have n(th) order accuracy.  by way of contrast, bilinear interpolation */
/*     schemes are exact for polynomials of degree one.  bicubic interpolation */
/*     is exact only for polynomials of degree three or less.  the method */
/*     also produces a weighted least squares fit to the data in which waves */
/*     are resolved uniformly on the full sphere.  high frequencies, induced */
/*     by closeness of grid points near the poles (due to computational */
/*     or observational errors) are smoothed.  finally, the method is */
/*     consistent with methods used to generate data in numerical spectral */
/*     models based on spherical harmonics.  for more discussion of these and */
/*     related issues,  see the article: "on the spectral approximation of */
/*     discrete scalar and vector functions on the sphere," siam j. numer. */
/*     anal., vol 16. dec 1979, pp. 934-949, by paul swarztrauber. */


/* *** comment */

/*     on a nlon by nlat or nlat by nlon grid (gaussian or equally spaced) */
/*     spherical harmonic analysis generates and synthesis utilizes */
/*     min0(nlat,(nlon+2)/2)) by nlat coefficients.  consequently, for */
/*     da and db,  if either */

/*             min0(nlatb,(nlonb+2)/2) < min0(nlata,(nlona+2)/2) */

/*     or if */

/*             nlatb < nlata */

/*     then all the coefficients generated by an analysis of da cannot be used */
/*     in the synthesis which generates db.  in this case "information" can be */
/*     lost in generating db.  more precisely, information will be lost if the */
/*     analysis of da yields nonzero coefficients which are outside the bounds */
/*     determined by the db grid.  nevertheless, transference of values with */
/*     spherical harmonics will yield results consistent with grid resolution */
/*     and is highly accurate. */


/* *** input arguments */

/* ... intl */

/*     an initialization argument which should be zero on an initial call to */
/*     trssph.  intl should be one if trssph is being recalled and */

/*          igrida,nlona,nlata,igridb,nlonb,nlatb */

/*     have not changed from the previous call.  if any of these arguments */
/*     have changed, intl=0 must be used to avoid undetectable errors.  calls */
/*     with intl=1 bypass redundant computation and save time.  it can be used */
/*     when transferring multiple data sets with the same underlying grids. */


/* ... igrida */

/*     an integer vector dimensioned two which identifies the underlying grid */
/*     on the full sphere for the given data array da as follows: */

/*     igrida(1) */

/*     = -1 */
/*     if the latitude (or colatitude) grid for da is an equally spaced */
/*     partition of [-pi/2,pi/2] ( or [0,pi] ) including the poles which */
/*     runs north to south */

/*     = +1 */
/*     if the latitude (or colatitude) grid for da is an equally spaced */
/*     partition of [-pi/2,pi/2] ( or [0,pi] ) including the poles which */
/*     runs south to north */

/*     = -2 */
/*     if the latitude (or colatitude) grid for da is a gaussian partition */
/*     of (-pi/2,pi/2) ( or (0,pi) ) excluding the poles which runs north */
/*     to south */

/*     = +2 */
/*     if the latitude (or colatitude) grid for da is a gaussian partition */
/*     of (-pi/2,pi/2) ( or (0,pi) ) excluding the poles which runs south */
/*     north */

/*     igrida(2) */

/*     = 0 if the underlying grid for da is a nlona by nlata */

/*     = 1 if the underlying grid for da is a nlata by nlona */


/* ... nlona */

/*     the number of longitude points on the uniform grid which partitions */
/*     [0,2pi) for the given data array da.  nlona is also the first or second */
/*     dimension of da (see igrida(2)) in the program which calls trssph. */
/*     nlona determines the grid increment in longitude as 2*pi/nlona. for */
/*     example nlona = 72 for a five degree grid.  nlona must be greater than */
/*     or equal to 4.  the efficiency of the computation is improved when */
/*     nlona is a product of small prime numbers */

/* ... nlata */

/*     the number of points in the latitude (or colatitude) grid */
/*     for the given data array da.  nlata is also the first or second */
/*     dimension of da (see igrida(2)) in the program which calls trssph. */
/*     if nlata is odd then the equator will be located at the (nlata+1)/2 */
/*     gaussian grid point.  if nlata is even then the equator will be */
/*     located half way between the nlata/2 and nlata/2+1 grid points. */

/* *** note: */
/*     igrida(1)=-1 or igrida(1)=-2 and igrida(2)=1 corresponds to */
/*     the "usual" mathematical spherical coordinate system required */
/*     by most of the drivers in spherepack2.  igrida(1)=1 or igrida(1)=2 */
/*     and igrida(2)=0 corresponds to the "usual" geophysical spherical */
/*     coordinate system. */

/* ... da */

/*     a two dimensional array that contains the data to be transferred. */
/*     da must be dimensioned nlona by nlata in the program calling trssph if */
/*     igrida(2) = 0.  da must be dimensioned nlata by nlona in the program */
/*     calling trssph if igrida(2) = 1.  if da is not properly dimensioned */
/*     and if the latitude (colatitude) values do not run south to north or */
/*     north to south as flagged by igrida(1) (this cannot be checked!) then */
/*     incorrect results will be produced. */

/* ... igridb */

/*     an integer vector dimensioned two which identifies the underlying grid */
/*     on the full sphere for the transformed data array db as follows: */

/*     igridb(1) */

/*     = -1 */
/*     if the latitude (or colatitude) grid for db is an equally spaced */
/*     partition of [-pi/2,pi/2] ( or [0,pi] ) including the poles which */
/*     north to south */

/*     = +1 */
/*     if the latitude (or colatitude) grid for db is an equally spaced */
/*     partition of [-pi/2,pi/2] ( or [0,pi] ) including the poles which */
/*     south to north */

/*     = -2 */
/*     if the latitude (or colatitude) grid for db is a gaussian partition */
/*     of (-pi/2,pi/2) ( or (0,pi) ) excluding the poles which runs north to */
/*     south */

/*     = +2 */
/*     if the latitude (or colatitude) grid for db is a gaussian partition */
/*     of (-pi/2,pi/2) ( or (0,pi) ) excluding the poles which runs south to */
/*     north */


/*     igridb(2) */

/*     = 0 if the underlying grid for db is a nlonb by nlatb */

/*     = 1 if the underlying grid for db is a nlatb by nlonb */


/* ... nlonb */

/*     the number of longitude points on the uniform grid which partitions */
/*     [0,2pi) for the transformed data array db.  nlonb is also the first or */
/*     second dimension of db (see igridb(2)) in the program which calls */
/*     trssph.  nlonb determines the grid increment in longitude as 2*pi/nlonb. */
/*     for example nlonb = 72 for a five degree grid.  nlonb must be greater */
/*     than or equal to 4.  the efficiency of the computation is improved when */
/*     nlonb is a product of small prime numbers */

/* ... nlatb */

/*     the number of points in the latitude (or colatitude) grid */
/*     for the transformed data array db.  nlatb is also the first or second */
/*     dimension of db (see igridb(2)) in the program which calls trssph. */
/*     if nlatb is odd then the equator will be located at the (nlatb+1)/2 */
/*     gaussian grid point.  if nlatb is even then the equator will be */
/*     located half way between the nlatb/2 and nlatb/2+1 grid points. */

/* ... wsave */

/*     a saved work space array that can be utilized repeatedly by trssph */
/*     as long as the arguments nlata,nlona,nlatb,nlonb remain unchanged. */
/*     wsave is set by a intl=0 call to trssph.  wsave must not be altered */
/*     when trssph is being recalled with intl=1. */

/* ... lsave */

/*     the dimension of the work space wsave as it appears in the program */
/*     that calls trssph.  the minimum required value of lsave for the */
/*     current set of input arguments is set in the output argument lsvmin. */
/*     it can be determined by calling trssph with lsave=0 and printing lsvmin. */
/*     let */

/*          lwa =  2*nlata*la2+3*((la1-2)*(nlata+nlata-la1-1))/2+nlona+15 */

/*     if the grid for da is equally spaced in (co)latitude.  let */

/*          lwa = nlata*(2*la2+3*la1-2)+3*la1*(1-la1)/2+nlona+15 */

/*     if the grid for da is gaussian in (co)latitude. */
/*     let */

/*          lwb = nlatb*(2*lb2+3*lb1-2)+3*lb1*(1-lb1)/2+nlonb+15 */

/*     if the grid for db is gaussian in (co)latitude.  let */

/*          lwb = 2*nlatb*lb2+3*((lb1-2)*(nlatb+nlatb-lb1-1))/2+nlonb+15 */

/*     if the grid for db is equally spaced in (co)latitude.  then */
/*     the quantity */

/*          lwa + lwb */

/*     is the minimum required length of wsave.  this value is returned */
/*     in the output argument lsvmin even if lsave is to small (ierror=10) */

/* ... work */

/*     a real work array that does not have to be preserved */

/* ... lwork */

/*     the dimension of the array work as it appears in the program */
/*     calling trssph. the minimum required value of lwork for the current */
/*     set of input arguments is set in the output argument lwkmin. */
/*     it can be determined by calling trssph with lwork=0 and printing */
/*     lwkmin.  an estimate for lwork follows.  let nlat,nlon,l1,l2 be */
/*     defined by */

/*       nlat = max0(nlata,nlatb), nlon = nax0(nlona,nlonb), */
/*       l1 = min0(nlat,(nlon+2)/2), l2 = (nlat+1)/2 */

/*     then the quantity */

/*          nlat*(4*l1+nlon+2*nlat+4)+3*((l1-2)*2*(2*nlat-l1-1))/2 */

/*     will suffice as a length for the unsaved work space. */

/*  *  both of the formulas above for lsave and lwork may overestimate the */
/*     required minimum values.  they can be predetermined by calling trssph */
/*     with lsave=lwork=0 and printout of lsvmin and lwkmin. */

/* ... dwork */

/*     a double precision work array that does not have to be preserved. */

/* ... ldwork */

/*     The length of dwork in the routine calling trssph. */
/*     Let */

/*       nlat = max0(nlata,nlatb) */

/*     ldwork must be at least nlat*(nlat+4) */

/* *** output arguments */


/* ... db */

/*     a two dimensional array that contains the transformed data.  db */
/*     must be dimensioned nlonb by nlatb in the program calling trssph if */
/*     igridb(2) = 0 or 1.  db must be dimensioned nlatb by nlonb in the */
/*     program calling trssph if igridb(2) = 1.  if db is not properly */
/*     dimensioned and if the latitude (colatitude) values do not run south */
/*     north or north to south as flagged by igrdb(1) (this cannot be checked!) */
/*     then incorrect results will be produced. */

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

/*         = 6  if igridb(1) is not -1 or +1 or -2 or +2 */

/*         = 7  if igridb(2) is not 0 or 1 */

/*         = 8  if nlonb is less than 4 */

/*         = 9  if nlatb is less than 3 */
/*         =10  if there is insufficient saved work space (lsave < lsvmin) */

/*         =11  if there is insufficient unsaved work space (lwork < lwkmin) */

/*         =12  indicates failure in an eigenvalue routine which computes */
/*              gaussian weights and points */

/*         =13  if ldwork is too small (insufficient unsaved double precision */
/*              work space) */

/* ***************************************************** */
/* ***************************************************** */

/*     end of argument description ... code follows */

/* ***************************************************** */
/* ***************************************************** */

/* Subroutine */ int trssph_(integer *intl, integer *igrida, integer *nlona, 
	integer *nlata, real *da, integer *igridb, integer *nlonb, integer *
	nlatb, real *db, real *wsave, integer *lsave, integer *lsvmin, real *
	work, integer *lwork, integer *lwkmin, doublereal *dwork, integer *
	ldwork, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer jb, ig, iw, nt, lw, la1, la2, lb1, lb2, iaa, iab, iba, ibb,
	     lwa, lwb, lwk3, lwk4;
    extern /* Subroutine */ int trab_(integer *, integer *, real *, real *, 
	    integer *, integer *, real *, real *);
    static integer nlat, isym;
    extern /* Subroutine */ int shaec_(integer *, integer *, integer *, 
	    integer *, real *, integer *, integer *, real *, real *, integer *
	    , integer *, real *, integer *, real *, integer *, integer *), 
	    shagc_(integer *, integer *, integer *, integer *, real *, 
	    integer *, integer *, real *, real *, integer *, integer *, real *
	    , integer *, real *, integer *, integer *);
    static integer igrda, igrdb;
    extern /* Subroutine */ int shsec_(integer *, integer *, integer *, 
	    integer *, real *, integer *, integer *, real *, real *, integer *
	    , integer *, real *, integer *, real *, integer *, integer *), 
	    shsgc_(integer *, integer *, integer *, integer *, real *, 
	    integer *, integer *, real *, real *, integer *, integer *, real *
	    , integer *, real *, integer *, integer *), shaeci_(integer *, 
	    integer *, real *, integer *, doublereal *, integer *, integer *),
	     shagci_(integer *, integer *, real *, integer *, doublereal *, 
	    integer *, integer *), shseci_(integer *, integer *, real *, 
	    integer *, doublereal *, integer *, integer *), shsgci_(integer *,
	     integer *, real *, integer *, doublereal *, integer *, integer *)
	    , convlat_(integer *, integer *, real *), trsplat_(integer *, 
	    integer *, real *, real *);


/*     include a save statement to ensure local variables in trssph, set during */
/*     an intl=0 call, are preserved if trssph is recalled with intl=1 */


/*     check input arguments */

    /* Parameter adjustments */
    --dwork;
    --work;
    --wsave;
    --db;
    --igridb;
    --da;
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
    ig = igridb[1];
    if ((ig - 1) * (ig + 1) * (ig - 2) * (ig + 2) != 0) {
	return 0;
    }
    *ier = 7;
    ig = igridb[2];
    if (ig * (ig - 1) != 0) {
	return 0;
    }
    *ier = 8;
    if (*nlonb < 4) {
	return 0;
    }
    *ier = 9;
    if (*nlatb < 3) {
	return 0;
    }
    *ier = 0;
    igrda = abs(igrida[1]);
    igrdb = abs(igridb[1]);
    if (*intl == 0) {
/* Computing MIN */
	i__1 = *nlata, i__2 = (*nlona + 2) / 2;
	la1 = min(i__1,i__2);
	la2 = (*nlata + 1) / 2;
/* Computing MIN */
	i__1 = *nlatb, i__2 = (*nlonb + 2) / 2;
	lb1 = min(i__1,i__2);
	lb2 = (*nlatb + 1) / 2;

/*     set saved work space length for analysis */

	if (igrda == 1) {

/*     saved space for analysis on  equally spaced grid */

	    lwa = (*nlata << 1) * la2 + (la1 - 2) * (*nlata + *nlata - la1 - 
		    1) * 3 / 2 + *nlona + 15;
	} else {

/*     saved space for analysis on gaussian grid */

	    lwa = *nlata * ((la2 << 1) + la1 * 3 - 2) + la1 * 3 * (1 - la1) / 
		    2 + *nlona + 15;
	}

/*     set wsave pointer */

	jb = lwa + 1;

/*     set pointers for spherical harmonic coefs */

	iaa = 1;
	iba = iaa + la1 * *nlata;
	iab = iba + la1 * *nlata;
	if (igrdb == 2) {

/*     set saved work space length for gaussian synthesis */

	    lwb = *nlatb * ((lb2 << 1) + lb1 * 3 - 2) + lb1 * 3 * (1 - lb1) / 
		    2 + *nlonb + 15;
	} else {

/*     set saved work space length for equally spaced synthesis */

	    lwb = (*nlatb << 1) * lb2 + (lb1 - 2) * (*nlatb + *nlatb - lb1 - 
		    1) * 3 / 2 + *nlonb + 15;
	}

/*     set minimum saved work space length */

	*lsvmin = lwa + lwb;

/*     set remaining harmonic pointer */

	ibb = iab + lb1 * *nlatb;

/*     set pointers for remaining work */

	iw = ibb + lb1 * *nlatb;

/*     set remaining work space length in lw */

	lw = *lwork - iw;
	lwk3 = *nlata * *nlona << 1;
	lwk4 = *nlatb * *nlonb << 1;

/*     set minimum unsaved work space required by trssph */

	*lwkmin = iw + max(lwk3,lwk4);

/*     set error flags if saved or unsaved work spaces are insufficient */

	*ier = 10;
	if (*lsave < *lsvmin) {
	    return 0;
	}
	*ier = 11;
	if (*lwork < *lwkmin) {
	    return 0;
	}
	*ier = 13;
	nlat = max(*nlata,*nlatb);
	if (*ldwork < nlat * (nlat + 4)) {
	    return 0;
	}
	*ier = 0;
	if (igrda == 1) {

/*     initialize wsave for equally spaced analysis */

	    shaeci_(nlata, nlona, &wsave[1], &lwa, &dwork[1], ldwork, ier);
	} else {

/*     initialize wsave for gaussian analysis */

	    shagci_(nlata, nlona, &wsave[1], &lwa, &dwork[1], ldwork, ier);
	    if (*ier != 0) {

/*     flag failure in spherepack gaussian software */

		*ier = 12;
		return 0;
	    }
	}
	if (igrdb == 2) {

/*     initialize wsave for gaussian synthesis */

	    shsgci_(nlatb, nlonb, &wsave[jb], &lwb, &dwork[1], ldwork, ier);
	    if (*ier != 0) {

/*     flag failure in spherepack gaussian software */

		*ier = 12;
		return 0;
	    }
	} else {

/*     initialize wsave for equally spaced synthesis */

	    shseci_(nlatb, nlonb, &wsave[jb], &lwb, &dwork[1], ldwork, ier);
	}

/*     end of initialization (intl=0) call */

    }

/*     transpose and/or reorder (co)latitude if necessary for da */
/*     (arrays must have latitude (colatitude) as the first dimension */
/*     and run north to south for spherepack software) */

    if (igrida[2] == 0) {
	trsplat_(nlona, nlata, &da[1], &work[1]);
    }
    if (igrida[1] > 0) {
	convlat_(nlata, nlona, &da[1]);
    }
    nt = 1;
    isym = 0;
    if (igrda == 2) {

/*     do spherical harmonic analysis of "adjusted" da on gaussian grid */

	shagc_(nlata, nlona, &isym, &nt, &da[1], nlata, nlona, &work[iaa], &
		work[iba], &la1, nlata, &wsave[1], &lwa, &work[iw], &lw, ier);
    } else {

/*     do spherical harmonic analysis of "adjusted" da on equally spaced grid */

	shaec_(nlata, nlona, &isym, &nt, &da[1], nlata, nlona, &work[iaa], &
		work[iba], &la1, nlata, &wsave[1], &lwa, &work[iw], &lw, ier);
    }

/*     transfer da grid coefficients to db grid coefficients */
/*     truncating to zero as necessary */

    trab_(&la1, nlata, &work[iaa], &work[iba], &lb1, nlatb, &work[iab], &work[
	    ibb]);
    if (igrdb == 1) {

/*     do spherical harmonic synthesis on nlatb by nlonb equally spaced grid */

	shsec_(nlatb, nlonb, &isym, &nt, &db[1], nlatb, nlonb, &work[iab], &
		work[ibb], &lb1, nlatb, &wsave[jb], &lwb, &work[iw], &lw, ier)
		;
    } else {

/*     do spherical harmonic synthesis on nlatb by nlonb gaussian grid */

	shsgc_(nlatb, nlonb, &isym, &nt, &db[1], nlatb, nlonb, &work[iab], &
		work[ibb], &lb1, nlatb, &wsave[jb], &lwb, &work[iw], &lw, ier)
		;
    }

/*     both da,db are currently latitude by longitude north to south arrays */
/*     restore da and set db to agree with flags in igrida and igridb */

    if (igrida[1] > 0) {
	convlat_(nlata, nlona, &da[1]);
    }
    if (igridb[1] > 0) {
	convlat_(nlatb, nlonb, &db[1]);
    }
    if (igrida[2] == 0) {
	trsplat_(nlata, nlona, &da[1], &work[1]);
    }
    if (igridb[2] == 0) {
	trsplat_(nlatb, nlonb, &db[1], &work[1]);
    }
    return 0;
} /* trssph_ */

/* Subroutine */ int trab_(integer *ma, integer *na, real *aa, real *ba, 
	integer *mb, integer *nb, real *ab, real *bb)
{
    /* System generated locals */
    integer aa_dim1, aa_offset, ba_dim1, ba_offset, ab_dim1, ab_offset, 
	    bb_dim1, bb_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, m, n;


/*     set coefficients for b grid from coefficients for a grid */

    /* Parameter adjustments */
    ba_dim1 = *ma;
    ba_offset = 1 + ba_dim1;
    ba -= ba_offset;
    aa_dim1 = *ma;
    aa_offset = 1 + aa_dim1;
    aa -= aa_offset;
    bb_dim1 = *mb;
    bb_offset = 1 + bb_dim1;
    bb -= bb_offset;
    ab_dim1 = *mb;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;

    /* Function Body */
    m = min(*ma,*mb);
    n = min(*na,*nb);
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ab[i__ + j * ab_dim1] = aa[i__ + j * aa_dim1];
	    bb[i__ + j * bb_dim1] = ba[i__ + j * ba_dim1];
	}
    }

/*     set coefs outside triangle to zero */

    i__1 = *mb;
    for (i__ = m + 1; i__ <= i__1; ++i__) {
	i__2 = *nb;
	for (j = 1; j <= i__2; ++j) {
	    ab[i__ + j * ab_dim1] = 0.f;
	    bb[i__ + j * bb_dim1] = 0.f;
	}
    }
    i__1 = *nb;
    for (j = n + 1; j <= i__1; ++j) {
	i__2 = *mb;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ab[i__ + j * ab_dim1] = 0.f;
	    bb[i__ + j * bb_dim1] = 0.f;
	}
    }
    return 0;
} /* trab_ */

/* Subroutine */ int trsplat_(integer *n, integer *m, real *data, real *work)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, ij, ji;


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
} /* trsplat_ */

/* Subroutine */ int convlat_(integer *nlat, integer *nlon, real *data)
{
    /* System generated locals */
    integer data_dim1, data_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ib;
    static real temp;
    static integer nlat2;


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
} /* convlat_ */

