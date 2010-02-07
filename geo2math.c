/* geo2math.f -- translated by f2c (version 20061008).
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
/*  .                          SPHEREPACK                         . */
/*  .                                                             . */
/*  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */



/* ... geo2math.f */

/*     file geo2math.f contains subroutines for converting scalar and */
/*     vector fields between geophysical and mathematical spherical */
/*     coordinates.  The latter is required when using most spherepack */
/*     software.  The four main subroutines in geo2math.f are described */
/*     as follows: */

/*     (1) subroutine geo2maths(ig,nlon,nlat,sg,sm,work) */

/*         converts the nlon by nlat scalar field sg given in */
/*         geophysical coordinates to the nlat by nlon scalar */
/*         field sm given in mathematical coordinates. sg and sm */
/*         can be identical in the program calling geo2maths. */

/*     (2) subroutine math2geos(ig,nlat,nlon,sm,sg,work) */

/*         converts the nlat by nlon scalar field sm given in */
/*         mathematical coordinates to the nlon by nlat scalar */
/*         field sg given in geophysical coordinates.  sm and */
/*         sg can be identical in the program calling math2geos. */

/*     (3) subroutine geo2mathv(ig,nlon,nlat,ug,vg,vm,wm,work) */

/*         converts the nlon by nlat vector field (ug,vg) given */
/*         in geophysical coordinates to the nlat by nlon vector */
/*         field (vm,wm) in mathematical coordinates.  ug and wm */
/*         can be identical in the program calling geo2mathv.  vg */
/*         and vm can be identical in the program calling geo2mathv. */

/*     (4) subroutine math2geov(ig,nlat,nlon,vm,wm,ug,vg,work) */

/*         converts the nlat by nlon vector field (vm,wm) given */
/*         in mathematical coordinates to the nlon by nlat vector */
/*         field (ug,vg) in spherical coordinates.  vm and vg can */
/*         be identical in the program calling math2geov.  wm and */
/*         ug can be identical in the program calling math2geov. */

/* *** (1),(2),(3),(4) argument description. */

/* ... ig */

/*     = 0 if the latitude values in the geophysical arrays sg,ug,vg are */
/*         ordered south to north with increasing latitude subscript */
/*         i=1,2,...,nlat. */

/*     = 1 if the latitude values in the geophysical arrays sg,ug,vg are */
/*         ordered north to south with increasing latitude subscript */
/*         i=1,2,...,nlat. */

/* ... nlon */

/*     the number of distinct londitude points.  nlon determines */
/*     the grid increment in longitude as 2*pi/nlon. nlon is the first */
/*     dimension of the geophysical arrays sg,ug,vg and the second */
/*     dimension of the mathematical arrays sm,vm,wm.  The longitude */
/*     grid is given by phi(j) = (j-1)*2*pi/nlon j=1,...,nlon. */

/* ... nlat */

/*     the number of distinct latitude and colatitude points and the */
/*     first dimension of sm,vm,wm and second dimension of sg,ug,vg. */
/*     If the (co)latitude grid is equally spaced then the grid increment */
/*     is dlat=pi/(nlat-1).  In this case the south to north latitude grid */
/*     underlying is */

/*        lat(i) = -0.5*pi + (i-1)*dlat (i=1,...,nlat) */

/*     and the north to south colatitude grid underlying sm,vm,wm is */

/*        colat(i) = (i-1)*dlat  (i=1,...,nlat) */

/*     If the grid is Gaussian let thetag(i) be the north to south colatitude */
/*     grid (as computed by the spherepack routine gaqd).  In this case */

/*        colat(i) = thetag(i) (i=1,...,nlat) */

/*     and */

/*        lat(i) = -0.5*pi + thetag(i) (i=1,...,nlat) */

/*     In either case lat(i) = colat(nlat-i+1) for all i. */
/*     If nlat is odd the equator is located at the (nlat+1)/2 */
/*     latitude or colatitude point.  If nlat is even the equator is */
/*     half way between the nlat/2 and nlat/2+1 latitude or colatitude */
/*     points.  The equally spaced (co)latitude grid includes the poles. */
/*     The Gaussian grid excludes the poles. */

/* ... sg,sm */

/*     In (1),(2) sg is a nlon by nlat array containing the scalar field */
/*     in geophysical coordinates.  Latitude values in sg are ordered from */
/*     the southern to the northern hemisphere with increasing latitude */
/*     subscript if ig = 0 or ordered from the northern hemisphere to the */
/*     southern hemisphere if ig = 1.  sm is a nlat by nlon array containing */
/*     the scalar field in mathematical coordinates.  Colatitude values in sm */
/*     are ordered from the north to the south hemisphere  with increasing */
/*     colatitude subscript (i=1,...,nlat).  The (co)latitude grid for sg and */
/*     sm can be equally spaced or Gaussian.  sg and sm can be equivalenced or */
/*     be identical in the routine calling geo2maths or math2geos.  sg and */
/*     sm are related by */

/*          sm(nlat-i+1,j) = sg(j,i)    (if ig = 0) */

/*     or */

/*          sm(i,j) = sg(j,i)            (if ig = 1) */

/*     for i=1,...,nlat and j=1,...,nlon. This formula is not used because */
/*     the two arrays can be equivalenced or identical arguments in the */
/*     program calling geo2maths or math2geos. */

/* ... ug,vg,vm,wm */

/*     In (3),(4) ug is a nlon by nlat array containing the longitudinal */
/*     vector component.  vg is a nlon by nlat array containing the */
/*     latitudinal vector component.  Values in (ug,vg) are ordered */
/*     from the southern to the northern hemisphere with increasing */
/*     latitude subscript if ig = 0 or from the northern to southern */
/*     hemisphere if ig = 1.  vm is a nlat by nlon array containing the */
/*     the colatitudinal vector component.  wm is a nlat by nlon array */
/*     containing the east longitudinal vector component.  Values in */
/*     (vm,wm) are ordered from the northern to the southern hemisphere */
/*     with increasing colatitude subscript.  The (co)latitude grid for */
/*     both vector fields can be equally spaced or Gaussian.  ug,wm and */
/*     vg,vm can be equivalenced or be identical in the program calling */
/*     geo2mathv or math2geov.  They are related by */

/*          ug(j,nlat-i+1) =  wm(i,j) */
/*                                      (ig = 0) */
/*          vg(j,nlat-i+1) = -vm(i,j) */

/*     or */

/*          ug(j,i) =  wm(i,j) */
/*                                      (ig = 1) */
/*          vg(j,i) = -vm(i,j) */


/*     for i=1,...,nlat and j=1,...,nlon.  These formulas are not */
/*     used because ug,wm and vg,vm can be equivalenced or identical */
/*     arguments in the program calling math2geov or geo2mathv. */

/*     Let ib = nlat-i+1 for i=1,...,nlat.  Summarizing: */
/*     sg(j,i) or ug(j,i),vg(j,i) are values at (phi(j),lat(i)) if ig = 0 */
/*     sg(j,i) or ug(j,i),vg(j,i) are values at (phi(j),lat(ib)) if ig = 1 */
/*     sm(i,j) or vm(i,j),wm(i,j) are values at (colat(i),phi(j)) */

/* ... work is an unsaved real work space of length at least nlon*nlat */
/*     in the routine calling (1),(2),(3), or (4).  It is used to simplify */
/*     a nonsquare array transposition in case it is required. */

/* *** example (1) */

/*     suppose you wish to compute the divergence of (ug,vg) on a Gaussian */
/*     grid in geophysical coordinates using the stored Legendre polynomial */
/*     routines on SPHEREPACK 2.0. */

/*     (1) call geo2mathv to set vm,wm from ug,vg */

/*     (2) call vhags to compute the vector harmonic coefficients of vm,wm */

/*     (3) call divgs with the coefficients from (2) to compute the divergence */
/*         dv in mathematical spherical coordinates on the UNIT sphere. */

/*     (4) call math2geos to convert the scalar divergence dv back to */
/*         geophysical spherical coordinates. */

/*     (5) divide dv by R (the radius of the earth) to compute divergence */
/*         on the earth (scaling from unit sphere computation in (3)). */

/* *** example (2) */

/*     suppose you wish to compute a vector field (ug,vg) corresponding */
/*     to  a given divergence dvg and vorticity vtg (all in geophysical */
/*     coordinates) on an equally spaced (co)latitude grid using the */
/*     computed Legendre polynomial software. */

/*     (1) call geo2maths to set dvm from dvg */

/*     (2) call geo2maths to set vtm from vts */

/*     (3) call shaec to compute the scalar harmonic coefficients of dvm */

/*     (4) call shaec to compute the scalar harmonic coefficients of vtm */

/*     (5) call idvtec to compute (vm,wm) using the coefficients from (3),(4). */

/*     (6) call math2geov to set (ug,vg) from (vm,wm) */

/*     (7) multiply (ug,vg) by the earth's radius R for scaling */
/*         from the unit sphere computation in (5) */

/* *** END OF DOCUMENTATION ... CODE FOLLOWS: */


/* Subroutine */ int geo2maths_(integer *ig, integer *nlon, integer *nlat, 
	real *sg, real *sm, real *work)
{
    /* System generated locals */
    integer sg_dim1, sg_offset, sm_dim1, sm_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ij;


/*     transpose sg into sm and reverse colatitude subscript order */
/*     if necessary */

    /* Parameter adjustments */
    sm_dim1 = *nlat;
    sm_offset = 1 + sm_dim1;
    sm -= sm_offset;
    sg_dim1 = *nlon;
    sg_offset = 1 + sg_dim1;
    sg -= sg_offset;
    --work;

    /* Function Body */
    i__1 = *nlat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    ij = (j - 1) * *nlat + i__;
	    work[ij] = sg[j + i__ * sg_dim1];
	}
    }
    if (*ig == 0) {
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		sm[*nlat - i__ + 1 + j * sm_dim1] = work[ij];
	    }
	}
    } else {
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		sm[i__ + j * sm_dim1] = work[ij];
	    }
	}
    }
    return 0;
} /* geo2maths_ */

/* Subroutine */ int math2geos_(integer *ig, integer *nlat, integer *nlon, 
	real *sm, real *sg, real *work)
{
    /* System generated locals */
    integer sm_dim1, sm_offset, sg_dim1, sg_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ij;


/*     transpose sm into sg and reverse colatitude subscript order */
/*     if necessary */

    /* Parameter adjustments */
    sg_dim1 = *nlon;
    sg_offset = 1 + sg_dim1;
    sg -= sg_offset;
    sm_dim1 = *nlat;
    sm_offset = 1 + sm_dim1;
    sm -= sm_offset;
    --work;

    /* Function Body */
    i__1 = *nlat;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *nlon;
	for (j = 1; j <= i__2; ++j) {
	    ij = (j - 1) * *nlat + i__;
	    work[ij] = sm[i__ + j * sm_dim1];
	}
    }
    if (*ig == 0) {
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		sg[j + (*nlat - i__ + 1) * sg_dim1] = work[ij];
	    }
	}
    } else {
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		sg[j + i__ * sg_dim1] = work[ij];
	    }
	}
    }
    return 0;
} /* math2geos_ */

/* Subroutine */ int geo2mathv_(integer *ig, integer *nlon, integer *nlat, 
	real *ug, real *vg, real *vm, real *wm, real *work)
{
    /* System generated locals */
    integer ug_dim1, ug_offset, vg_dim1, vg_offset, vm_dim1, vm_offset, 
	    wm_dim1, wm_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ij;


/*     convert vg to vm, ug to wm */

    /* Parameter adjustments */
    wm_dim1 = *nlat;
    wm_offset = 1 + wm_dim1;
    wm -= wm_offset;
    vm_dim1 = *nlat;
    vm_offset = 1 + vm_dim1;
    vm -= vm_offset;
    vg_dim1 = *nlon;
    vg_offset = 1 + vg_dim1;
    vg -= vg_offset;
    ug_dim1 = *nlon;
    ug_offset = 1 + ug_dim1;
    ug -= ug_offset;
    --work;

    /* Function Body */
    if (*ig == 0) {
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		work[ij] = vg[j + i__ * vg_dim1];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		vm[*nlat - i__ + 1 + j * vm_dim1] = -work[ij];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		work[ij] = ug[j + i__ * ug_dim1];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		wm[*nlat - i__ + 1 + j * wm_dim1] = work[ij];
	    }
	}
    } else {
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		work[ij] = vg[j + i__ * vg_dim1];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		vm[i__ + j * vm_dim1] = -work[ij];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		work[ij] = ug[j + i__ * ug_dim1];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		wm[i__ + j * wm_dim1] = work[ij];
	    }
	}
    }
    return 0;
} /* geo2mathv_ */

/* Subroutine */ int math2geov_(integer *ig, integer *nlat, integer *nlon, 
	real *vm, real *wm, real *ug, real *vg, real *work)
{
    /* System generated locals */
    integer vm_dim1, vm_offset, wm_dim1, wm_offset, ug_dim1, ug_offset, 
	    vg_dim1, vg_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ij;


/*     convert vm to vg, wm to ug */

    /* Parameter adjustments */
    vg_dim1 = *nlon;
    vg_offset = 1 + vg_dim1;
    vg -= vg_offset;
    ug_dim1 = *nlon;
    ug_offset = 1 + ug_dim1;
    ug -= ug_offset;
    wm_dim1 = *nlat;
    wm_offset = 1 + wm_dim1;
    wm -= wm_offset;
    vm_dim1 = *nlat;
    vm_offset = 1 + vm_dim1;
    vm -= vm_offset;
    --work;

    /* Function Body */
    if (*ig == 0) {
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		work[ij] = vm[i__ + j * vm_dim1];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		vg[j + (*nlat - i__ + 1) * vg_dim1] = -work[ij];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		work[ij] = wm[i__ + j * wm_dim1];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		ug[j + (*nlat - i__ + 1) * ug_dim1] = work[ij];
	    }
	}
    } else {
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		work[ij] = vm[i__ + j * vm_dim1];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		vg[j + i__ * vg_dim1] = -work[ij];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		work[ij] = wm[i__ + j * wm_dim1];
	    }
	}
	i__1 = *nlat;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = *nlon;
	    for (j = 1; j <= i__2; ++j) {
		ij = (j - 1) * *nlat + i__;
		ug[j + i__ * ug_dim1] = work[ij];
	    }
	}
    }
    return 0;
} /* math2geov_ */

