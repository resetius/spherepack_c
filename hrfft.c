/* hrfft.f -- translated by f2c (version 20061008).
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

/* Common Block Declarations */

struct {
    doublereal tfft;
} hrf_;

#define hrf_1 hrf_


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


/* ... file hrfft.f */

/*     this file contains a multiple fft package for spherepack3.0. */
/*     it includes code and documentation for performing fast fourier */
/*     transforms (see subroutines hrffti,hrfftf and hrfftb) */

/* ********************************************************************** */

/*     subroutine hrffti(n,wsave) */

/*     subroutine hrffti initializes the array wsave which is used in */
/*     both hrfftf and hrfftb. the prime factorization of n together */
/*     with a tabulation of the trigonometric functions are computed and */
/*     stored in wsave. */

/*     input parameter */

/*     n       the length of the sequence to be transformed. */

/*     output parameter */

/*     wsave   a work array which must be dimensioned at least 2*n+15. */
/*             the same work array can be used for both hrfftf and */
/*             hrfftb as long as n remains unchanged. different wsave */
/*             arrays are required for different values of n. the */
/*             contents of wsave must not be changed between calls */
/*             of hrfftf or hrfftb. */

/* ********************************************************************** */

/*     subroutine hrfftf(m,n,r,mdimr,wsave,work) */

/*     subroutine hrfftf computes the fourier coefficients of m doublereal */
/*     perodic sequences (fourier analysis); i.e. hrfftf computes the */
/*     doublereal fft of m sequences each with length n. the transform is */
/*     defined below at output parameter r. */

/*     input parameters */

/*     m       the number of sequences. */

/*     n       the length of all m sequences.  the method is most */
/*             efficient when n is a product of small primes. n may */
/*             change as long as different work arrays are provided */

/*     r       r(m,n) is a two dimensional doublereal array that contains m */
/*             sequences each with length n. */

/*     mdimr   the first dimension of the r array as it appears */
/*             in the program that calls hrfftf. mdimr must be */
/*             greater than or equal to m. */


/*     wsave   a work array with at least least 2*n+15 locations */
/*             in the program that calls hrfftf. the wsave array must be */
/*             initialized by calling subroutine hrffti(n,wsave) and a */
/*             different wsave array must be used for each different */
/*             value of n. this initialization does not have to be */
/*             repeated so long as n remains unchanged thus subsequent */
/*             transforms can be obtained faster than the first. */
/*             the same wsave array can be used by hrfftf and hrfftb. */

/*     work    a doublereal work array with m*n locations. */


/*     output parameters */

/*     r      for all j=1,...,m */

/*             r(j,1) = the sum from i=1 to i=n of r(j,i) */

/*             if n is even set l =n/2   , if n is odd set l = (n+1)/2 */

/*               then for k = 2,...,l */

/*                  r(j,2*k-2) = the sum from i = 1 to i = n of */

/*                       r(j,i)*cos((k-1)*(i-1)*2*pi/n) */

/*                  r(j,2*k-1) = the sum from i = 1 to i = n of */

/*                      -r(j,i)*sin((k-1)*(i-1)*2*pi/n) */

/*             if n is even */

/*                  r(j,n) = the sum from i = 1 to i = n of */

/*                       (-1)**(i-1)*r(j,i) */

/*      *****  note */
/*                  this transform is unnormalized since a call of hrfftf */
/*                  followed by a call of hrfftb will multiply the input */
/*                  sequence by n. */

/*     wsave   contains results which must not be destroyed between */
/*             calls of hrfftf or hrfftb. */

/*     work    a doublereal work array with m*n locations that does */
/*             not have to be saved. */

/* ********************************************************************** */

/*     subroutine hrfftb(m,n,r,mdimr,wsave,work) */

/*     subroutine hrfftb computes the doublereal perodic sequence of m */
/*     sequences from their fourier coefficients (fourier synthesis). */
/*     the transform is defined below at output parameter r. */

/*     input parameters */

/*     m       the number of sequences. */

/*     n       the length of all m sequences.  the method is most */
/*             efficient when n is a product of small primes. n may */
/*             change as long as different work arrays are provided */

/*     r       r(m,n) is a two dimensional doublereal array that contains */
/*             the fourier coefficients of m sequences each with */
/*             length n. */

/*     mdimr   the first dimension of the r array as it appears */
/*             in the program that calls hrfftb. mdimr must be */
/*             greater than or equal to m. */

/*     wsave   a work array which must be dimensioned at least 2*n+15. */
/*             in the program that calls hrfftb. the wsave array must be */
/*             initialized by calling subroutine hrffti(n,wsave) and a */
/*             different wsave array must be used for each different */
/*             value of n. this initialization does not have to be */
/*             repeated so long as n remains unchanged thus subsequent */
/*             transforms can be obtained faster than the first. */
/*             the same wsave array can be used by hrfftf and hrfftb. */

/*     work    a doublereal work array with m*n locations. */


/*     output parameters */

/*     r      for all j=1,...,m */

/*             for n even and for i = 1,...,n */

/*                  r(j,i) = r(j,1)+(-1)**(i-1)*r(j,n) */

/*                       plus the sum from k=2 to k=n/2 of */

/*                        2.*r(j,2*k-2)*cos((k-1)*(i-1)*2*pi/n) */

/*                       -2.*r(j,2*k-1)*sin((k-1)*(i-1)*2*pi/n) */

/*             for n odd and for i = 1,...,n */

/*                  r(j,i) = r(j,1) plus the sum from k=2 to k=(n+1)/2 of */

/*                       2.*r(j,2*k-2)*cos((k-1)*(i-1)*2*pi/n) */

/*                      -2.*r(j,2*k-1)*sin((k-1)*(i-1)*2*pi/n) */

/*      *****  note */
/*                  this transform is unnormalized since a call of hrfftf */
/*                  followed by a call of hrfftb will multiply the input */
/*                  sequence by n. */

/*     wsave   contains results which must not be destroyed between */
/*             calls of hrfftb or hrfftf. */

/*     work    a doublereal work array with m*n locations that does not */
/*             have to be saved */

/* ********************************************************************** */



/* Subroutine */ int hrffti_(integer *n, doublereal *wsave)
{
    extern /* Subroutine */ int hrfti1_(integer *, doublereal *, doublereal *);

    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    hrf_1.tfft = 0.f;
    if (*n == 1) {
	return 0;
    }
    hrfti1_(n, &wsave[1], &wsave[*n + 1]);
    return 0;
} /* hrffti_ */

/* Subroutine */ int hrfti1_(integer *n, doublereal *wa, doublereal *fac)
{
    /* Initialized data */

    static integer ntryh[4] = { 4,2,3,5 };

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double atan(doublereal), cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, k1, l1, l2, ib;
    static doublereal fi;
    static integer ld, ii, nf, ip, nl, is, nq, nr;
    static doublereal arg;
    static integer ido, ipm;
    static doublereal tpi;
    static integer nfm1;
    static doublereal argh;
    static integer ntry;
    static doublereal argld;


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa;
    --fac;

    /* Function Body */
    nl = *n;
    nf = 0;
    j = 0;
L101:
    ++j;
    if (j - 4 <= 0) {
	goto L102;
    } else {
	goto L103;
    }
L102:
    ntry = ntryh[j - 1];
    goto L104;
L103:
    ntry += 2;
L104:
    nq = nl / ntry;
    nr = nl - ntry * nq;
    if (nr != 0) {
	goto L101;
    } else {
	goto L105;
    }
L105:
    ++nf;
    fac[nf + 2] = (doublereal) ntry;
    nl = nq;
    if (ntry != 2) {
	goto L107;
    }
    if (nf == 1) {
	goto L107;
    }
    i__1 = nf;
    for (i__ = 2; i__ <= i__1; ++i__) {
	ib = nf - i__ + 2;
	fac[ib + 2] = fac[ib + 1];
/* L106: */
    }
    fac[3] = 2.f;
L107:
    if (nl != 1) {
	goto L104;
    }
    fac[1] = (doublereal) (*n);
    fac[2] = (doublereal) nf;
    tpi = atan(1.) * 8.;
    argh = tpi / (doublereal) (*n);
    is = 0;
    nfm1 = nf - 1;
    l1 = 1;
    if (nfm1 == 0) {
	return 0;
    }
    i__1 = nfm1;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = fac[k1 + 2];
	ld = 0;
	l2 = l1 * ip;
	ido = *n / l2;
	ipm = ip - 1;
	i__2 = ipm;
	for (j = 1; j <= i__2; ++j) {
	    ld += l1;
	    i__ = is;
	    argld = (doublereal) ld * argh;
	    fi = 0.f;
	    i__3 = ido;
	    for (ii = 3; ii <= i__3; ii += 2) {
		i__ += 2;
		fi += 1.f;
		arg = fi * argld;
		wa[i__ - 1] = cos(arg);
		wa[i__] = sin(arg);
/* L108: */
	    }
	    is += ido;
/* L109: */
	}
	l1 = l2;
/* L110: */
    }
    return 0;
} /* hrfti1_ */

/* Subroutine */ int hrfftf_(integer *m, integer *n, doublereal *r__, integer *
	mdimr, doublereal *whrfft, doublereal *work)
{
    /* System generated locals */
    integer r_dim1, r_offset;

    /* Local variables */
    extern /* Subroutine */ int hrftf1_(integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *, doublereal *);


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --whrfft;
    r_dim1 = *mdimr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --work;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
/*     tstart = second(dum) */
    hrftf1_(m, n, &r__[r_offset], mdimr, &work[1], &whrfft[1], &whrfft[*n + 1]
	    );
/*     tfft = tfft+second(dum)-tstart */
    return 0;
} /* hrfftf_ */

/* Subroutine */ int hrftf1_(integer *m, integer *n, doublereal *c__, integer *
	mdimc, doublereal *ch, doublereal *wa, doublereal *fac)
{
    /* System generated locals */
    integer ch_dim1, ch_offset, c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k1, l1, l2, na, kh, nf, ip, iw, ix2, ix3, ix4, ido,
	     idl1;
    extern /* Subroutine */ int hradf2_(integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *), hradf3_(integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *
	    , doublereal *), hradf4_(integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *), hradf5_(
	    integer *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *), hradfg_(integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, doublereal *, doublereal *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *);


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa;
    ch_dim1 = *m;
    ch_offset = 1 + ch_dim1;
    ch -= ch_offset;
    c_dim1 = *mdimc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --fac;

    /* Function Body */
    nf = fac[2];
    na = 1;
    l2 = *n;
    iw = *n;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	kh = nf - k1;
	ip = fac[kh + 3];
	l1 = l2 / ip;
	ido = *n / l2;
	idl1 = ido * l1;
	iw -= (ip - 1) * ido;
	na = 1 - na;
	if (ip != 4) {
	    goto L102;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	if (na != 0) {
	    goto L101;
	}
	hradf4_(m, &ido, &l1, &c__[c_offset], mdimc, &ch[ch_offset], m, &wa[
		iw], &wa[ix2], &wa[ix3]);
	goto L110;
L101:
	hradf4_(m, &ido, &l1, &ch[ch_offset], m, &c__[c_offset], mdimc, &wa[
		iw], &wa[ix2], &wa[ix3]);
	goto L110;
L102:
	if (ip != 2) {
	    goto L104;
	}
	if (na != 0) {
	    goto L103;
	}
	hradf2_(m, &ido, &l1, &c__[c_offset], mdimc, &ch[ch_offset], m, &wa[
		iw]);
	goto L110;
L103:
	hradf2_(m, &ido, &l1, &ch[ch_offset], m, &c__[c_offset], mdimc, &wa[
		iw]);
	goto L110;
L104:
	if (ip != 3) {
	    goto L106;
	}
	ix2 = iw + ido;
	if (na != 0) {
	    goto L105;
	}
	hradf3_(m, &ido, &l1, &c__[c_offset], mdimc, &ch[ch_offset], m, &wa[
		iw], &wa[ix2]);
	goto L110;
L105:
	hradf3_(m, &ido, &l1, &ch[ch_offset], m, &c__[c_offset], mdimc, &wa[
		iw], &wa[ix2]);
	goto L110;
L106:
	if (ip != 5) {
	    goto L108;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	ix4 = ix3 + ido;
	if (na != 0) {
	    goto L107;
	}
	hradf5_(m, &ido, &l1, &c__[c_offset], mdimc, &ch[ch_offset], m, &wa[
		iw], &wa[ix2], &wa[ix3], &wa[ix4]);
	goto L110;
L107:
	hradf5_(m, &ido, &l1, &ch[ch_offset], m, &c__[c_offset], mdimc, &wa[
		iw], &wa[ix2], &wa[ix3], &wa[ix4]);
	goto L110;
L108:
	if (ido == 1) {
	    na = 1 - na;
	}
	if (na != 0) {
	    goto L109;
	}
	hradfg_(m, &ido, &ip, &l1, &idl1, &c__[c_offset], &c__[c_offset], &
		c__[c_offset], mdimc, &ch[ch_offset], &ch[ch_offset], m, &wa[
		iw]);
	na = 1;
	goto L110;
L109:
	hradfg_(m, &ido, &ip, &l1, &idl1, &ch[ch_offset], &ch[ch_offset], &ch[
		ch_offset], m, &c__[c_offset], &c__[c_offset], mdimc, &wa[iw])
		;
	na = 0;
L110:
	l2 = l1;
/* L111: */
    }
    if (na == 1) {
	return 0;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    c__[i__ + j * c_dim1] = ch[i__ + j * ch_dim1];
/* L112: */
	}
    }
    return 0;
} /* hrftf1_ */

/* Subroutine */ int hradf4_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2, 
	doublereal *wa3)
{
    /* System generated locals */
    integer cc_dim1, cc_dim2, cc_dim3, cc_offset, ch_dim1, ch_dim2, ch_offset,
	     i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, k, m, ic, idp2;
    static doublereal hsqt2;


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa3;
    --wa2;
    --wa1;
    cc_dim1 = *mdimcc;
    cc_dim2 = *ido;
    cc_dim3 = *l1;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * (1 + cc_dim3));
    cc -= cc_offset;
    ch_dim1 = *mdimch;
    ch_dim2 = *ido;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * 5);
    ch -= ch_offset;

    /* Function Body */
    hsqt2 = sqrt(2.f) / 2.f;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + (((k << 2) + 1) * ch_dim2 + 1) * ch_dim1] = cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + (
		    cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] + (cc[m + ((k + 
		    cc_dim3) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + cc_dim3 
		    * 3) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + (*ido + ((k << 2) + 4) * ch_dim2) * ch_dim1] = cc[m + ((k 
		    + cc_dim3) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + 
		    cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1] - (cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + (
		    cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + (*ido + ((k << 2) + 2) * ch_dim2) * ch_dim1] = cc[m + ((k 
		    + cc_dim3) * cc_dim2 + 1) * cc_dim1] - cc[m + ((k + 
		    cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1];
	    ch[m + (((k << 2) + 3) * ch_dim2 + 1) * ch_dim1] = cc[m + ((k + (
		    cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] - cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1];
/* L1001: */
	}
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + (i__ - 1 + ((k << 2) + 1) * ch_dim2) * ch_dim1] = wa1[
			i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k 
			+ (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa3[i__ - 
			2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2)
			 * cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1]) + (cc[m + (i__ 
			- 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + (wa2[i__ 
			- 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) 
			* cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k + 
			cc_dim3 * 3) * cc_dim2) * cc_dim1]));
		ch[m + (ic - 1 + ((k << 2) + 4) * ch_dim2) * ch_dim1] = cc[m 
			+ (i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + (
			wa2[i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1]) - (wa1[i__ - 2] 
			* cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1] + (wa3[i__ - 2] * cc[m + 
			(i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1] 
			+ wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 << 2)) * 
			cc_dim2) * cc_dim1]));
		ch[m + (i__ + ((k << 2) + 1) * ch_dim2) * ch_dim1] = wa1[i__ 
			- 2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) *
			 cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa3[i__ - 2] *
			 cc[m + (i__ + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] - wa3[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1]) + (cc[m + (i__ 
			+ (k + cc_dim3) * cc_dim2) * cc_dim1] + (wa2[i__ - 2] 
			* cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) * 
			cc_dim1] - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 3) * cc_dim2) * cc_dim1]));
		ch[m + (ic + ((k << 2) + 4) * ch_dim2) * ch_dim1] = wa1[i__ - 
			2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa3[i__ - 2] *
			 cc[m + (i__ + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] - wa3[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1]) - (cc[m + (i__ 
			+ (k + cc_dim3) * cc_dim2) * cc_dim1] + (wa2[i__ - 2] 
			* cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) * 
			cc_dim1] - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 3) * cc_dim2) * cc_dim1]));
		ch[m + (i__ - 1 + ((k << 2) + 3) * ch_dim2) * ch_dim1] = wa1[
			i__ - 2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa3[
			i__ - 2] * cc[m + (i__ + (k + (cc_dim3 << 2)) * 
			cc_dim2) * cc_dim1] - wa3[i__ - 1] * cc[m + (i__ - 1 
			+ (k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1]) + (cc[m 
			+ (i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] - (
			wa2[i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1]));
		ch[m + (ic - 1 + ((k << 2) + 2) * ch_dim2) * ch_dim1] = cc[m 
			+ (i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] - (
			wa2[i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1]) - (wa1[i__ - 2] 
			* cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa3[i__ - 2] *
			 cc[m + (i__ + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] - wa3[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1]));
		ch[m + (i__ + ((k << 2) + 3) * ch_dim2) * ch_dim1] = wa3[i__ 
			- 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * 
			cc_dim2) * cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k 
			+ (cc_dim3 << 2)) * cc_dim2) * cc_dim1] - (wa1[i__ - 
			2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2)
			 * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1]) + (cc[m + (i__ 
			+ (k + cc_dim3) * cc_dim2) * cc_dim1] - (wa2[i__ - 2] 
			* cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) * 
			cc_dim1] - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 3) * cc_dim2) * cc_dim1]));
		ch[m + (ic + ((k << 2) + 2) * ch_dim2) * ch_dim1] = wa3[i__ - 
			2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2)
			 * cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1] - (wa1[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1]) - (cc[m + (i__ + (k + 
			cc_dim3) * cc_dim2) * cc_dim1] - (wa2[i__ - 2] * cc[m 
			+ (i__ + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] - 
			wa2[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1]));
/* L1003: */
	    }
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + (*ido + ((k << 2) + 1) * ch_dim2) * ch_dim1] = hsqt2 * (cc[
		    m + (*ido + (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] - 
		    cc[m + (*ido + (k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1])
		     + cc[m + (*ido + (k + cc_dim3) * cc_dim2) * cc_dim1];
	    ch[m + (*ido + ((k << 2) + 3) * ch_dim2) * ch_dim1] = cc[m + (*
		    ido + (k + cc_dim3) * cc_dim2) * cc_dim1] - hsqt2 * (cc[m 
		    + (*ido + (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] - cc[
		    m + (*ido + (k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1]);
	    ch[m + (((k << 2) + 2) * ch_dim2 + 1) * ch_dim1] = -hsqt2 * (cc[m 
		    + (*ido + (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + cc[
		    m + (*ido + (k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1]) - 
		    cc[m + (*ido + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1];
	    ch[m + (((k << 2) + 4) * ch_dim2 + 1) * ch_dim1] = -hsqt2 * (cc[m 
		    + (*ido + (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + cc[
		    m + (*ido + (k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1]) + 
		    cc[m + (*ido + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1];
/* L1006: */
	}
/* L106: */
    }
L107:
    return 0;
} /* hradf4_ */

/* Subroutine */ int hradf2_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1)
{
    /* System generated locals */
    integer ch_dim1, ch_dim2, ch_offset, cc_dim1, cc_dim2, cc_dim3, cc_offset,
	     i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, m, ic, idp2;


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa1;
    cc_dim1 = *mdimcc;
    cc_dim2 = *ido;
    cc_dim3 = *l1;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * (1 + cc_dim3));
    cc -= cc_offset;
    ch_dim1 = *mdimch;
    ch_dim2 = *ido;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * 3);
    ch -= ch_offset;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + (((k << 1) + 1) * ch_dim2 + 1) * ch_dim1] = cc[m + ((k + 
		    cc_dim3) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1];
	    ch[m + (*ido + ((k << 1) + 2) * ch_dim2) * ch_dim1] = cc[m + ((k 
		    + cc_dim3) * cc_dim2 + 1) * cc_dim1] - cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1];
/* L1001: */
	}
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + (i__ + ((k << 1) + 1) * ch_dim2) * ch_dim1] = cc[m + (
			i__ + (k + cc_dim3) * cc_dim2) * cc_dim1] + (wa1[i__ 
			- 2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) *
			 cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1]);
		ch[m + (ic + ((k << 1) + 2) * ch_dim2) * ch_dim1] = wa1[i__ - 
			2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] - cc[m + (i__ + (
			k + cc_dim3) * cc_dim2) * cc_dim1];
		ch[m + (i__ - 1 + ((k << 1) + 1) * ch_dim2) * ch_dim1] = cc[m 
			+ (i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + (
			wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) 
			* cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (
			k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1]);
		ch[m + (ic - 1 + ((k << 1) + 2) * ch_dim2) * ch_dim1] = cc[m 
			+ (i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] - (
			wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) 
			* cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (
			k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1]);
/* L1003: */
	    }
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + (((k << 1) + 2) * ch_dim2 + 1) * ch_dim1] = -cc[m + (*ido 
		    + (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1];
	    ch[m + (*ido + ((k << 1) + 1) * ch_dim2) * ch_dim1] = cc[m + (*
		    ido + (k + cc_dim3) * cc_dim2) * cc_dim1];
/* L1006: */
	}
/* L106: */
    }
L107:
    return 0;
} /* hradf2_ */

/* Subroutine */ int hradf3_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2)
{
    /* System generated locals */
    integer ch_dim1, ch_dim2, ch_offset, cc_dim1, cc_dim2, cc_dim3, cc_offset,
	     i__1, i__2, i__3;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, k, m, ic;
    static doublereal arg;
    static integer idp2;
    static doublereal taui, taur;
    extern doublereal pimach_(void);


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa2;
    --wa1;
    cc_dim1 = *mdimcc;
    cc_dim2 = *ido;
    cc_dim3 = *l1;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * (1 + cc_dim3));
    cc -= cc_offset;
    ch_dim1 = *mdimch;
    ch_dim2 = *ido;
    ch_offset = 1 + ch_dim1 * (1 + (ch_dim2 << 2));
    ch -= ch_offset;

    /* Function Body */
    arg = pimach_() * 2.f / 3.f;
    taur = cos(arg);
    taui = sin(arg);
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + ((k * 3 + 1) * ch_dim2 + 1) * ch_dim1] = cc[m + ((k + 
		    cc_dim3) * cc_dim2 + 1) * cc_dim1] + (cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + 
		    cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + ((k * 3 + 3) * ch_dim2 + 1) * ch_dim1] = taui * (cc[m + ((
		    k + cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1] - cc[m + ((k + 
		    (cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + (*ido + (k * 3 + 2) * ch_dim2) * ch_dim1] = cc[m + ((k + 
		    cc_dim3) * cc_dim2 + 1) * cc_dim1] + taur * (cc[m + ((k + 
		    (cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + 
		    cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
/* L1001: */
	}
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + (i__ - 1 + (k * 3 + 1) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + (wa1[
			i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k 
			+ (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa2[i__ - 
			2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * 
			cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k + cc_dim3 *
			 3) * cc_dim2) * cc_dim1]));
		ch[m + (i__ + (k * 3 + 1) * ch_dim2) * ch_dim1] = cc[m + (i__ 
			+ (k + cc_dim3) * cc_dim2) * cc_dim1] + (wa1[i__ - 2] 
			* cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa2[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1]
			 - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) 
			* cc_dim2) * cc_dim1]));
		ch[m + (i__ - 1 + (k * 3 + 3) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + taur *
			 (wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)
			) * cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1])) + taui * (wa1[
			i__ - 2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa2[
			i__ - 2] * cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) 
			* cc_dim1] - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 3) * cc_dim2) * cc_dim1]));
		ch[m + (ic - 1 + (k * 3 + 2) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + taur *
			 (wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)
			) * cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1])) - taui * (wa1[
			i__ - 2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa2[
			i__ - 2] * cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) 
			* cc_dim1] - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 3) * cc_dim2) * cc_dim1]));
		ch[m + (i__ + (k * 3 + 3) * ch_dim2) * ch_dim1] = cc[m + (i__ 
			+ (k + cc_dim3) * cc_dim2) * cc_dim1] + taur * (wa1[
			i__ - 2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa2[
			i__ - 2] * cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) 
			* cc_dim1] - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 3) * cc_dim2) * cc_dim1])) + taui * (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] - (wa1[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1]));
		ch[m + (ic + (k * 3 + 2) * ch_dim2) * ch_dim1] = taui * (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] - (wa1[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1])) - (cc[m + (i__ + (k + 
			cc_dim3) * cc_dim2) * cc_dim1] + taur * (wa1[i__ - 2] 
			* cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa2[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1]
			 - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) 
			* cc_dim2) * cc_dim1])));
/* L1002: */
	    }
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* hradf3_ */

/* Subroutine */ int hradf5_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2, 
	doublereal *wa3, doublereal *wa4)
{
    /* System generated locals */
    integer cc_dim1, cc_dim2, cc_dim3, cc_offset, ch_dim1, ch_dim2, ch_offset,
	     i__1, i__2, i__3;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, k, m, ic;
    static doublereal arg, ti11, ti12, tr11, tr12;
    static integer idp2;
    extern doublereal pimach_(void);


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa4;
    --wa3;
    --wa2;
    --wa1;
    cc_dim1 = *mdimcc;
    cc_dim2 = *ido;
    cc_dim3 = *l1;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * (1 + cc_dim3));
    cc -= cc_offset;
    ch_dim1 = *mdimch;
    ch_dim2 = *ido;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * 6);
    ch -= ch_offset;

    /* Function Body */
    arg = pimach_() * 2.f / 5.f;
    tr11 = cos(arg);
    ti11 = sin(arg);
    tr12 = cos(arg * 2.f);
    ti12 = sin(arg * 2.f);
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + ((k * 5 + 1) * ch_dim2 + 1) * ch_dim1] = cc[m + ((k + 
		    cc_dim3) * cc_dim2 + 1) * cc_dim1] + (cc[m + ((k + 
		    cc_dim3 * 5) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1]) + (cc[m + ((k + 
		    (cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + 
		    cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + (*ido + (k * 5 + 2) * ch_dim2) * ch_dim1] = cc[m + ((k + 
		    cc_dim3) * cc_dim2 + 1) * cc_dim1] + tr11 * (cc[m + ((k + 
		    cc_dim3 * 5) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1]) + tr12 * (cc[m 
		    + ((k + (cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] + cc[m 
		    + ((k + cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + ((k * 5 + 3) * ch_dim2 + 1) * ch_dim1] = ti11 * (cc[m + ((
		    k + cc_dim3 * 5) * cc_dim2 + 1) * cc_dim1] - cc[m + ((k + 
		    (cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1]) + ti12 * (cc[m 
		    + ((k + (cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] - cc[m 
		    + ((k + cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + (*ido + (k * 5 + 4) * ch_dim2) * ch_dim1] = cc[m + ((k + 
		    cc_dim3) * cc_dim2 + 1) * cc_dim1] + tr12 * (cc[m + ((k + 
		    cc_dim3 * 5) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1]) + tr11 * (cc[m 
		    + ((k + (cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] + cc[m 
		    + ((k + cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + ((k * 5 + 5) * ch_dim2 + 1) * ch_dim1] = ti12 * (cc[m + ((
		    k + cc_dim3 * 5) * cc_dim2 + 1) * cc_dim1] - cc[m + ((k + 
		    (cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1]) - ti11 * (cc[m 
		    + ((k + (cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] - cc[m 
		    + ((k + cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
/* L1001: */
	}
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + (i__ - 1 + (k * 5 + 1) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + (wa1[
			i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k 
			+ (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[i__ - 
			2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * cc_dim2) * 
			cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k + cc_dim3 *
			 5) * cc_dim2) * cc_dim1])) + (wa2[i__ - 2] * cc[m + (
			i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] + 
			wa2[i__ - 1] * cc[m + (i__ + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + (wa3[i__ - 2] * cc[m + (i__ - 1 
			+ (k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1] + wa3[
			i__ - 1] * cc[m + (i__ + (k + (cc_dim3 << 2)) * 
			cc_dim2) * cc_dim1]));
		ch[m + (i__ + (k * 5 + 1) * ch_dim2) * ch_dim1] = cc[m + (i__ 
			+ (k + cc_dim3) * cc_dim2) * cc_dim1] + (wa1[i__ - 2] 
			* cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) * cc_dim1]
			 - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) 
			* cc_dim2) * cc_dim1])) + (wa2[i__ - 2] * cc[m + (i__ 
			+ (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] - wa2[i__ - 
			1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * 
			cc_dim1] + (wa3[i__ - 2] * cc[m + (i__ + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 1] * 
			cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1]));
		ch[m + (i__ - 1 + (k * 5 + 3) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + tr11 *
			 (wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)
			) * cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1]) + tr12 * (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] + wa3[i__ - 2] * 
			cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1]) + ti11 * (wa1[i__ - 2] *
			 cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa4[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) * cc_dim1]
			 - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) 
			* cc_dim2) * cc_dim1])) + ti12 * (wa2[i__ - 2] * cc[m 
			+ (i__ + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] - 
			wa2[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] - (wa3[i__ - 2] * cc[m + (i__ + (
			k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 
			1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2)
			 * cc_dim1]));
		ch[m + (ic - 1 + (k * 5 + 2) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + tr11 *
			 (wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)
			) * cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1]) + tr12 * (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] + wa3[i__ - 2] * 
			cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1]) - (ti11 * (wa1[i__ - 2] 
			* cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa4[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) * cc_dim1]
			 - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) 
			* cc_dim2) * cc_dim1])) + ti12 * (wa2[i__ - 2] * cc[m 
			+ (i__ + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] - 
			wa2[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] - (wa3[i__ - 2] * cc[m + (i__ + (
			k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 
			1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2)
			 * cc_dim1])));
		ch[m + (i__ + (k * 5 + 3) * ch_dim2) * ch_dim1] = cc[m + (i__ 
			+ (k + cc_dim3) * cc_dim2) * cc_dim1] + tr11 * (wa1[
			i__ - 2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[
			i__ - 2] * cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) 
			* cc_dim1] - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 5) * cc_dim2) * cc_dim1])) + tr12 * (wa2[
			i__ - 2] * cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) 
			* cc_dim1] - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 3) * cc_dim2) * cc_dim1] + (wa3[i__ - 2] * 
			cc[m + (i__ + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] - wa3[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1])) + (ti11 * (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1] - (wa1[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1])) + ti12 * (wa3[i__ - 2] 
			* cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1] - (wa2[i__ - 2] * cc[m + 
			(i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] + 
			wa2[i__ - 1] * cc[m + (i__ + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1])));
		ch[m + (ic + (k * 5 + 2) * ch_dim2) * ch_dim1] = ti11 * (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1] - (wa1[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1])) + ti12 * (wa3[i__ - 2] 
			* cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1] - (wa2[i__ - 2] * cc[m + 
			(i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] + 
			wa2[i__ - 1] * cc[m + (i__ + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1])) - (cc[m + (i__ + (k + cc_dim3) *
			 cc_dim2) * cc_dim1] + tr11 * (wa1[i__ - 2] * cc[m + (
			i__ + (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] - 
			wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) 
			* cc_dim2) * cc_dim1] + (wa4[i__ - 2] * cc[m + (i__ + 
			(k + cc_dim3 * 5) * cc_dim2) * cc_dim1] - wa4[i__ - 1]
			 * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * cc_dim2) * 
			cc_dim1])) + tr12 * (wa2[i__ - 2] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] - wa2[i__ - 1] * 
			cc[m + (i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * 
			cc_dim1] + (wa3[i__ - 2] * cc[m + (i__ + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 1] * 
			cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1])));
		ch[m + (i__ - 1 + (k * 5 + 5) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + tr12 *
			 (wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)
			) * cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1])) + tr11 * (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] + (wa3[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1])) + (ti12 * (wa1[i__ - 2]
			 * cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa4[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) * cc_dim1]
			 - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) 
			* cc_dim2) * cc_dim1])) - ti11 * (wa2[i__ - 2] * cc[m 
			+ (i__ + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] - 
			wa2[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] - (wa3[i__ - 2] * cc[m + (i__ + (
			k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 
			1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2)
			 * cc_dim1])));
		ch[m + (ic - 1 + (k * 5 + 4) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + tr12 *
			 (wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)
			) * cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1])) + tr11 * (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] + (wa3[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1])) - (ti12 * (wa1[i__ - 2]
			 * cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa4[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) * cc_dim1]
			 - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) 
			* cc_dim2) * cc_dim1])) - ti11 * (wa2[i__ - 2] * cc[m 
			+ (i__ + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] - 
			wa2[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] - (wa3[i__ - 2] * cc[m + (i__ + (
			k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 
			1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2)
			 * cc_dim1])));
		ch[m + (i__ + (k * 5 + 5) * ch_dim2) * ch_dim1] = cc[m + (i__ 
			+ (k + cc_dim3) * cc_dim2) * cc_dim1] + tr12 * (wa1[
			i__ - 2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[
			i__ - 2] * cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) 
			* cc_dim1] - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 5) * cc_dim2) * cc_dim1])) + tr11 * (wa2[
			i__ - 2] * cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) 
			* cc_dim1] - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 3) * cc_dim2) * cc_dim1] + (wa3[i__ - 2] * 
			cc[m + (i__ + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] - wa3[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1])) + (ti12 * (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1] - (wa1[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1])) - ti11 * (wa3[i__ - 2] 
			* cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1] - (wa2[i__ - 2] * cc[m + 
			(i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] + 
			wa2[i__ - 1] * cc[m + (i__ + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1])));
		ch[m + (ic + (k * 5 + 4) * ch_dim2) * ch_dim1] = ti12 * (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1] - (wa1[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1])) - ti11 * (wa3[i__ - 2] 
			* cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1] - (wa2[i__ - 2] * cc[m + 
			(i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] + 
			wa2[i__ - 1] * cc[m + (i__ + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1])) - (cc[m + (i__ + (k + cc_dim3) *
			 cc_dim2) * cc_dim1] + tr12 * (wa1[i__ - 2] * cc[m + (
			i__ + (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] - 
			wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) 
			* cc_dim2) * cc_dim1] + (wa4[i__ - 2] * cc[m + (i__ + 
			(k + cc_dim3 * 5) * cc_dim2) * cc_dim1] - wa4[i__ - 1]
			 * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * cc_dim2) * 
			cc_dim1])) + tr11 * (wa2[i__ - 2] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] - wa2[i__ - 1] * 
			cc[m + (i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * 
			cc_dim1] + (wa3[i__ - 2] * cc[m + (i__ + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 1] * 
			cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1])));
/* L1002: */
	    }
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* hradf5_ */

/* Subroutine */ int hradfg_(integer *mp, integer *ido, integer *ip, integer *
	l1, integer *idl1, doublereal *cc, doublereal *c1, doublereal *c2, integer *mdimcc, 
	doublereal *ch, doublereal *ch2, integer *mdimch, doublereal *wa)
{
    /* System generated locals */
    integer ch_dim1, ch_dim2, ch_dim3, ch_offset, cc_dim1, cc_dim2, cc_dim3, 
	    cc_offset, c1_dim1, c1_dim2, c1_dim3, c1_offset, c2_dim1, c2_dim2,
	     c2_offset, ch2_dim1, ch2_dim2, ch2_offset, i__1, i__2, i__3, 
	    i__4;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, k, l, m, j2, ic, jc, lc, ik, is;
    static doublereal dc2, ai1, ai2, ar1, ar2, ds2;
    static integer nbd;
    static doublereal dcp, arg, dsp, tpi, ar1h, ar2h;
    static integer idp2, ipp2, idij, ipph;
    extern doublereal pimach_(void);


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa;
    c2_dim1 = *mdimcc;
    c2_dim2 = *idl1;
    c2_offset = 1 + c2_dim1 * (1 + c2_dim2);
    c2 -= c2_offset;
    c1_dim1 = *mdimcc;
    c1_dim2 = *ido;
    c1_dim3 = *l1;
    c1_offset = 1 + c1_dim1 * (1 + c1_dim2 * (1 + c1_dim3));
    c1 -= c1_offset;
    cc_dim1 = *mdimcc;
    cc_dim2 = *ido;
    cc_dim3 = *ip;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * (1 + cc_dim3));
    cc -= cc_offset;
    ch2_dim1 = *mdimch;
    ch2_dim2 = *idl1;
    ch2_offset = 1 + ch2_dim1 * (1 + ch2_dim2);
    ch2 -= ch2_offset;
    ch_dim1 = *mdimch;
    ch_dim2 = *ido;
    ch_dim3 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * (1 + ch_dim3));
    ch -= ch_offset;

    /* Function Body */
    tpi = pimach_() * 2.f;
    arg = tpi / (doublereal) (*ip);
    dcp = cos(arg);
    dsp = sin(arg);
    ipph = (*ip + 1) / 2;
    ipp2 = *ip + 2;
    idp2 = *ido + 2;
    nbd = (*ido - 1) / 2;
    if (*ido == 1) {
	goto L119;
    }
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch2[m + (ik + ch2_dim2) * ch2_dim1] = c2[m + (ik + c2_dim2) * 
		    c2_dim1];
/* L1001: */
	}
/* L101: */
    }
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1] = c1[m + (
			(k + j * c1_dim3) * c1_dim2 + 1) * c1_dim1];
/* L1002: */
	    }
/* L102: */
	}
/* L103: */
    }
    if (nbd > *l1) {
	goto L107;
    }
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	idij = is;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) * ch_dim1] 
			    = wa[idij - 1] * c1[m + (i__ - 1 + (k + j * 
			    c1_dim3) * c1_dim2) * c1_dim1] + wa[idij] * c1[m 
			    + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
		    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1] = 
			    wa[idij - 1] * c1[m + (i__ + (k + j * c1_dim3) * 
			    c1_dim2) * c1_dim1] - wa[idij] * c1[m + (i__ - 1 
			    + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
/* L1004: */
		}
/* L104: */
	    }
/* L105: */
	}
/* L106: */
    }
    goto L111;
L107:
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = is;
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		idij += 2;
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) * ch_dim1] 
			    = wa[idij - 1] * c1[m + (i__ - 1 + (k + j * 
			    c1_dim3) * c1_dim2) * c1_dim1] + wa[idij] * c1[m 
			    + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
		    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1] = 
			    wa[idij - 1] * c1[m + (i__ + (k + j * c1_dim3) * 
			    c1_dim2) * c1_dim1] - wa[idij] * c1[m + (i__ - 1 
			    + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
/* L1008: */
		}
/* L108: */
	    }
/* L109: */
	}
/* L110: */
    }
L111:
    if (nbd < *l1) {
	goto L115;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    c1[m + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) * c1_dim1] 
			    = ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) *
			     ch_dim1] + ch[m + (i__ - 1 + (k + jc * ch_dim3) *
			     ch_dim2) * ch_dim1];
		    c1[m + (i__ - 1 + (k + jc * c1_dim3) * c1_dim2) * c1_dim1]
			     = ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    c1[m + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] + ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    c1[m + (i__ + (k + jc * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m + (i__ - 1 + (k + jc * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1];
/* L1012: */
		}
/* L112: */
	    }
/* L113: */
	}
/* L114: */
    }
    goto L121;
L115:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    c1[m + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) * c1_dim1] 
			    = ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) *
			     ch_dim1] + ch[m + (i__ - 1 + (k + jc * ch_dim3) *
			     ch_dim2) * ch_dim1];
		    c1[m + (i__ - 1 + (k + jc * c1_dim3) * c1_dim2) * c1_dim1]
			     = ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    c1[m + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] + ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    c1[m + (i__ + (k + jc * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m + (i__ - 1 + (k + jc * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1];
/* L1016: */
		}
/* L116: */
	    }
/* L117: */
	}
/* L118: */
    }
    goto L121;
L119:
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    c2[m + (ik + c2_dim2) * c2_dim1] = ch2[m + (ik + ch2_dim2) * 
		    ch2_dim1];
/* L1020: */
	}
/* L120: */
    }
L121:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		c1[m + ((k + j * c1_dim3) * c1_dim2 + 1) * c1_dim1] = ch[m + (
			(k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1] + ch[m + (
			(k + jc * ch_dim3) * ch_dim2 + 1) * ch_dim1];
		c1[m + ((k + jc * c1_dim3) * c1_dim2 + 1) * c1_dim1] = ch[m + 
			((k + jc * ch_dim3) * ch_dim2 + 1) * ch_dim1] - ch[m 
			+ ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1];
/* L1022: */
	    }
/* L122: */
	}
/* L123: */
    }

    ar1 = 1.f;
    ai1 = 0.f;
    i__1 = ipph;
    for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	ar1h = dcp * ar1 - dsp * ai1;
	ai1 = dcp * ai1 + dsp * ar1;
	ar1 = ar1h;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch2[m + (ik + l * ch2_dim2) * ch2_dim1] = c2[m + (ik + 
			c2_dim2) * c2_dim1] + ar1 * c2[m + (ik + (c2_dim2 << 
			1)) * c2_dim1];
		ch2[m + (ik + lc * ch2_dim2) * ch2_dim1] = ai1 * c2[m + (ik + 
			*ip * c2_dim2) * c2_dim1];
/* L1024: */
	    }
/* L124: */
	}
	dc2 = ar1;
	ds2 = ai1;
	ar2 = ar1;
	ai2 = ai1;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    ar2h = dc2 * ar2 - ds2 * ai2;
	    ai2 = dc2 * ai2 + ds2 * ar2;
	    ar2 = ar2h;
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    ch2[m + (ik + l * ch2_dim2) * ch2_dim1] += ar2 * c2[m + (
			    ik + j * c2_dim2) * c2_dim1];
		    ch2[m + (ik + lc * ch2_dim2) * ch2_dim1] += ai2 * c2[m + (
			    ik + jc * c2_dim2) * c2_dim1];
/* L1025: */
		}
/* L125: */
	    }
/* L126: */
	}
/* L127: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch2[m + (ik + ch2_dim2) * ch2_dim1] += c2[m + (ik + j * 
			c2_dim2) * c2_dim1];
/* L1028: */
	    }
/* L128: */
	}
/* L129: */
    }

    if (*ido < *l1) {
	goto L132;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		cc[m + (i__ + (k * cc_dim3 + 1) * cc_dim2) * cc_dim1] = ch[m 
			+ (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1];
/* L1030: */
	    }
/* L130: */
	}
/* L131: */
    }
    goto L135;
L132:
    i__1 = *ido;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		cc[m + (i__ + (k * cc_dim3 + 1) * cc_dim2) * cc_dim1] = ch[m 
			+ (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1];
/* L1033: */
	    }
/* L133: */
	}
/* L134: */
    }
L135:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		cc[m + (*ido + (j2 - 2 + k * cc_dim3) * cc_dim2) * cc_dim1] = 
			ch[m + ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1];
		cc[m + ((j2 - 1 + k * cc_dim3) * cc_dim2 + 1) * cc_dim1] = ch[
			m + ((k + jc * ch_dim3) * ch_dim2 + 1) * ch_dim1];
/* L1036: */
	    }
/* L136: */
	}
/* L137: */
    }
    if (*ido == 1) {
	return 0;
    }
    if (nbd < *l1) {
	goto L141;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		ic = idp2 - i__;
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    cc[m + (i__ - 1 + (j2 - 1 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] + ch[m + (i__ - 1 + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m + (ic - 1 + (j2 - 2 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] - ch[m + (i__ - 1 + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m + (i__ + (j2 - 1 + k * cc_dim3) * cc_dim2) * cc_dim1]
			     = ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] + ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    cc[m + (ic + (j2 - 2 + k * cc_dim3) * cc_dim2) * cc_dim1] 
			    = ch[m + (i__ + (k + jc * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1];
/* L1038: */
		}
/* L138: */
	    }
/* L139: */
	}
/* L140: */
    }
    return 0;
L141:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    cc[m + (i__ - 1 + (j2 - 1 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] + ch[m + (i__ - 1 + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m + (ic - 1 + (j2 - 2 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] - ch[m + (i__ - 1 + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m + (i__ + (j2 - 1 + k * cc_dim3) * cc_dim2) * cc_dim1]
			     = ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] + ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    cc[m + (ic + (j2 - 2 + k * cc_dim3) * cc_dim2) * cc_dim1] 
			    = ch[m + (i__ + (k + jc * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1];
/* L1042: */
		}
/* L142: */
	    }
/* L143: */
	}
/* L144: */
    }
    return 0;
} /* hradfg_ */

doublereal pimach_(void)
{
    /* System generated locals */
    doublereal ret_val;

    ret_val = 3.14159265358979f;
    return ret_val;
} /* pimach_ */

/* Subroutine */ int hrfftb_(integer *m, integer *n, doublereal *r__, integer *
	mdimr, doublereal *whrfft, doublereal *work)
{
    /* System generated locals */
    integer r_dim1, r_offset;

    /* Local variables */
    extern /* Subroutine */ int hrftb1_(integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *, doublereal *);


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --whrfft;
    r_dim1 = *mdimr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --work;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
/*     tstart = second(dum) */
    hrftb1_(m, n, &r__[r_offset], mdimr, &work[1], &whrfft[1], &whrfft[*n + 1]
	    );
/*     tfft = tfft+second(dum)-tstart */
    return 0;
} /* hrfftb_ */

/* Subroutine */ int hrftb1_(integer *m, integer *n, doublereal *c__, integer *
	mdimc, doublereal *ch, doublereal *wa, doublereal *fac)
{
    /* System generated locals */
    integer ch_dim1, ch_offset, c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k1, l1, l2, na, nf, ip, iw, ix2, ix3, ix4, ido, 
	    idl1;
    extern /* Subroutine */ int hradb2_(integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *), hradb3_(integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, integer *, doublereal *
	    , doublereal *), hradb4_(integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, doublereal *), hradb5_(
	    integer *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *), hradbg_(integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, doublereal *, doublereal *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *);


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa;
    ch_dim1 = *m;
    ch_offset = 1 + ch_dim1;
    ch -= ch_offset;
    c_dim1 = *mdimc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --fac;

    /* Function Body */
    nf = fac[2];
    na = 0;
    l1 = 1;
    iw = 1;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	ip = fac[k1 + 2];
	l2 = ip * l1;
	ido = *n / l2;
	idl1 = ido * l1;
	if (ip != 4) {
	    goto L103;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	if (na != 0) {
	    goto L101;
	}
	hradb4_(m, &ido, &l1, &c__[c_offset], mdimc, &ch[ch_offset], m, &wa[
		iw], &wa[ix2], &wa[ix3]);
	goto L102;
L101:
	hradb4_(m, &ido, &l1, &ch[ch_offset], m, &c__[c_offset], mdimc, &wa[
		iw], &wa[ix2], &wa[ix3]);
L102:
	na = 1 - na;
	goto L115;
L103:
	if (ip != 2) {
	    goto L106;
	}
	if (na != 0) {
	    goto L104;
	}
	hradb2_(m, &ido, &l1, &c__[c_offset], mdimc, &ch[ch_offset], m, &wa[
		iw]);
	goto L105;
L104:
	hradb2_(m, &ido, &l1, &ch[ch_offset], m, &c__[c_offset], mdimc, &wa[
		iw]);
L105:
	na = 1 - na;
	goto L115;
L106:
	if (ip != 3) {
	    goto L109;
	}
	ix2 = iw + ido;
	if (na != 0) {
	    goto L107;
	}
	hradb3_(m, &ido, &l1, &c__[c_offset], mdimc, &ch[ch_offset], m, &wa[
		iw], &wa[ix2]);
	goto L108;
L107:
	hradb3_(m, &ido, &l1, &ch[ch_offset], m, &c__[c_offset], mdimc, &wa[
		iw], &wa[ix2]);
L108:
	na = 1 - na;
	goto L115;
L109:
	if (ip != 5) {
	    goto L112;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	ix4 = ix3 + ido;
	if (na != 0) {
	    goto L110;
	}
	hradb5_(m, &ido, &l1, &c__[c_offset], mdimc, &ch[ch_offset], m, &wa[
		iw], &wa[ix2], &wa[ix3], &wa[ix4]);
	goto L111;
L110:
	hradb5_(m, &ido, &l1, &ch[ch_offset], m, &c__[c_offset], mdimc, &wa[
		iw], &wa[ix2], &wa[ix3], &wa[ix4]);
L111:
	na = 1 - na;
	goto L115;
L112:
	if (na != 0) {
	    goto L113;
	}
	hradbg_(m, &ido, &ip, &l1, &idl1, &c__[c_offset], &c__[c_offset], &
		c__[c_offset], mdimc, &ch[ch_offset], &ch[ch_offset], m, &wa[
		iw]);
	goto L114;
L113:
	hradbg_(m, &ido, &ip, &l1, &idl1, &ch[ch_offset], &ch[ch_offset], &ch[
		ch_offset], m, &c__[c_offset], &c__[c_offset], mdimc, &wa[iw])
		;
L114:
	if (ido == 1) {
	    na = 1 - na;
	}
L115:
	l1 = l2;
	iw += (ip - 1) * ido;
/* L116: */
    }
    if (na == 0) {
	return 0;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    c__[i__ + j * c_dim1] = ch[i__ + j * ch_dim1];
/* L117: */
	}
    }
    return 0;
} /* hrftb1_ */

/* Subroutine */ int hradbg_(integer *mp, integer *ido, integer *ip, integer *
	l1, integer *idl1, doublereal *cc, doublereal *c1, doublereal *c2, integer *mdimcc, 
	doublereal *ch, doublereal *ch2, integer *mdimch, doublereal *wa)
{
    /* System generated locals */
    integer ch_dim1, ch_dim2, ch_dim3, ch_offset, cc_dim1, cc_dim2, cc_dim3, 
	    cc_offset, c1_dim1, c1_dim2, c1_dim3, c1_offset, c2_dim1, c2_dim2,
	     c2_offset, ch2_dim1, ch2_dim2, ch2_offset, i__1, i__2, i__3, 
	    i__4;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, k, l, m, j2, ic, jc, lc, ik, is;
    static doublereal dc2, ai1, ai2, ar1, ar2, ds2;
    static integer nbd;
    static doublereal dcp, arg, dsp, tpi, ar1h, ar2h;
    static integer idp2, ipp2, idij, ipph;
    extern doublereal pimach_(void);


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa;
    c2_dim1 = *mdimcc;
    c2_dim2 = *idl1;
    c2_offset = 1 + c2_dim1 * (1 + c2_dim2);
    c2 -= c2_offset;
    c1_dim1 = *mdimcc;
    c1_dim2 = *ido;
    c1_dim3 = *l1;
    c1_offset = 1 + c1_dim1 * (1 + c1_dim2 * (1 + c1_dim3));
    c1 -= c1_offset;
    cc_dim1 = *mdimcc;
    cc_dim2 = *ido;
    cc_dim3 = *ip;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * (1 + cc_dim3));
    cc -= cc_offset;
    ch2_dim1 = *mdimch;
    ch2_dim2 = *idl1;
    ch2_offset = 1 + ch2_dim1 * (1 + ch2_dim2);
    ch2 -= ch2_offset;
    ch_dim1 = *mdimch;
    ch_dim2 = *ido;
    ch_dim3 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * (1 + ch_dim3));
    ch -= ch_offset;

    /* Function Body */
    tpi = pimach_() * 2.f;
    arg = tpi / (doublereal) (*ip);
    dcp = cos(arg);
    dsp = sin(arg);
    idp2 = *ido + 2;
    nbd = (*ido - 1) / 2;
    ipp2 = *ip + 2;
    ipph = (*ip + 1) / 2;
    if (*ido < *l1) {
	goto L103;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m + (
			i__ + (k * cc_dim3 + 1) * cc_dim2) * cc_dim1];
/* L1001: */
	    }
/* L101: */
	}
/* L102: */
    }
    goto L106;
L103:
    i__1 = *ido;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m + (
			i__ + (k * cc_dim3 + 1) * cc_dim2) * cc_dim1];
/* L1004: */
	    }
/* L104: */
	}
/* L105: */
    }
L106:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1] = cc[m + (
			*ido + (j2 - 2 + k * cc_dim3) * cc_dim2) * cc_dim1] + 
			cc[m + (*ido + (j2 - 2 + k * cc_dim3) * cc_dim2) * 
			cc_dim1];
		ch[m + ((k + jc * ch_dim3) * ch_dim2 + 1) * ch_dim1] = cc[m + 
			((j2 - 1 + k * cc_dim3) * cc_dim2 + 1) * cc_dim1] + 
			cc[m + ((j2 - 1 + k * cc_dim3) * cc_dim2 + 1) * 
			cc_dim1];
/* L1007: */
	    }
/* L107: */
	}
/* L108: */
    }
    if (*ido == 1) {
	goto L116;
    }
    if (nbd < *l1) {
	goto L112;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		ic = idp2 - i__;
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) * ch_dim1] 
			    = cc[m + (i__ - 1 + ((j << 1) - 1 + k * cc_dim3) *
			     cc_dim2) * cc_dim1] + cc[m + (ic - 1 + ((j << 1) 
			    - 2 + k * cc_dim3) * cc_dim2) * cc_dim1];
		    ch[m + (i__ - 1 + (k + jc * ch_dim3) * ch_dim2) * ch_dim1]
			     = cc[m + (i__ - 1 + ((j << 1) - 1 + k * cc_dim3) 
			    * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + ((j << 1)
			     - 2 + k * cc_dim3) * cc_dim2) * cc_dim1];
		    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1] = 
			    cc[m + (i__ + ((j << 1) - 1 + k * cc_dim3) * 
			    cc_dim2) * cc_dim1] - cc[m + (ic + ((j << 1) - 2 
			    + k * cc_dim3) * cc_dim2) * cc_dim1];
		    ch[m + (i__ + (k + jc * ch_dim3) * ch_dim2) * ch_dim1] = 
			    cc[m + (i__ + ((j << 1) - 1 + k * cc_dim3) * 
			    cc_dim2) * cc_dim1] + cc[m + (ic + ((j << 1) - 2 
			    + k * cc_dim3) * cc_dim2) * cc_dim1];
/* L1009: */
		}
/* L109: */
	    }
/* L110: */
	}
/* L111: */
    }
    goto L116;
L112:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) * ch_dim1] 
			    = cc[m + (i__ - 1 + ((j << 1) - 1 + k * cc_dim3) *
			     cc_dim2) * cc_dim1] + cc[m + (ic - 1 + ((j << 1) 
			    - 2 + k * cc_dim3) * cc_dim2) * cc_dim1];
		    ch[m + (i__ - 1 + (k + jc * ch_dim3) * ch_dim2) * ch_dim1]
			     = cc[m + (i__ - 1 + ((j << 1) - 1 + k * cc_dim3) 
			    * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + ((j << 1)
			     - 2 + k * cc_dim3) * cc_dim2) * cc_dim1];
		    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1] = 
			    cc[m + (i__ + ((j << 1) - 1 + k * cc_dim3) * 
			    cc_dim2) * cc_dim1] - cc[m + (ic + ((j << 1) - 2 
			    + k * cc_dim3) * cc_dim2) * cc_dim1];
		    ch[m + (i__ + (k + jc * ch_dim3) * ch_dim2) * ch_dim1] = 
			    cc[m + (i__ + ((j << 1) - 1 + k * cc_dim3) * 
			    cc_dim2) * cc_dim1] + cc[m + (ic + ((j << 1) - 2 
			    + k * cc_dim3) * cc_dim2) * cc_dim1];
/* L1013: */
		}
/* L113: */
	    }
/* L114: */
	}
/* L115: */
    }
L116:
    ar1 = 1.f;
    ai1 = 0.f;
    i__1 = ipph;
    for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	ar1h = dcp * ar1 - dsp * ai1;
	ai1 = dcp * ai1 + dsp * ar1;
	ar1 = ar1h;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		c2[m + (ik + l * c2_dim2) * c2_dim1] = ch2[m + (ik + ch2_dim2)
			 * ch2_dim1] + ar1 * ch2[m + (ik + (ch2_dim2 << 1)) * 
			ch2_dim1];
		c2[m + (ik + lc * c2_dim2) * c2_dim1] = ai1 * ch2[m + (ik + *
			ip * ch2_dim2) * ch2_dim1];
/* L1017: */
	    }
/* L117: */
	}
	dc2 = ar1;
	ds2 = ai1;
	ar2 = ar1;
	ai2 = ai1;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    ar2h = dc2 * ar2 - ds2 * ai2;
	    ai2 = dc2 * ai2 + ds2 * ar2;
	    ar2 = ar2h;
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    c2[m + (ik + l * c2_dim2) * c2_dim1] += ar2 * ch2[m + (ik 
			    + j * ch2_dim2) * ch2_dim1];
		    c2[m + (ik + lc * c2_dim2) * c2_dim1] += ai2 * ch2[m + (
			    ik + jc * ch2_dim2) * ch2_dim1];
/* L1018: */
		}
/* L118: */
	    }
/* L119: */
	}
/* L120: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch2[m + (ik + ch2_dim2) * ch2_dim1] += ch2[m + (ik + j * 
			ch2_dim2) * ch2_dim1];
/* L1021: */
	    }
/* L121: */
	}
/* L122: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1] = c1[m + (
			(k + j * c1_dim3) * c1_dim2 + 1) * c1_dim1] - c1[m + (
			(k + jc * c1_dim3) * c1_dim2 + 1) * c1_dim1];
		ch[m + ((k + jc * ch_dim3) * ch_dim2 + 1) * ch_dim1] = c1[m + 
			((k + j * c1_dim3) * c1_dim2 + 1) * c1_dim1] + c1[m + 
			((k + jc * c1_dim3) * c1_dim2 + 1) * c1_dim1];
/* L1023: */
	    }
/* L123: */
	}
/* L124: */
    }
    if (*ido == 1) {
	goto L132;
    }
    if (nbd < *l1) {
	goto L128;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) * ch_dim1] 
			    = c1[m + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) *
			     c1_dim1] - c1[m + (i__ + (k + jc * c1_dim3) * 
			    c1_dim2) * c1_dim1];
		    ch[m + (i__ - 1 + (k + jc * ch_dim3) * ch_dim2) * ch_dim1]
			     = c1[m + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) 
			    * c1_dim1] + c1[m + (i__ + (k + jc * c1_dim3) * 
			    c1_dim2) * c1_dim1];
		    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1] = 
			    c1[m + (i__ + (k + j * c1_dim3) * c1_dim2) * 
			    c1_dim1] + c1[m + (i__ - 1 + (k + jc * c1_dim3) * 
			    c1_dim2) * c1_dim1];
		    ch[m + (i__ + (k + jc * ch_dim3) * ch_dim2) * ch_dim1] = 
			    c1[m + (i__ + (k + j * c1_dim3) * c1_dim2) * 
			    c1_dim1] - c1[m + (i__ - 1 + (k + jc * c1_dim3) * 
			    c1_dim2) * c1_dim1];
/* L1025: */
		}
/* L125: */
	    }
/* L126: */
	}
/* L127: */
    }
    goto L132;
L128:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) * ch_dim1] 
			    = c1[m + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) *
			     c1_dim1] - c1[m + (i__ + (k + jc * c1_dim3) * 
			    c1_dim2) * c1_dim1];
		    ch[m + (i__ - 1 + (k + jc * ch_dim3) * ch_dim2) * ch_dim1]
			     = c1[m + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) 
			    * c1_dim1] + c1[m + (i__ + (k + jc * c1_dim3) * 
			    c1_dim2) * c1_dim1];
		    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1] = 
			    c1[m + (i__ + (k + j * c1_dim3) * c1_dim2) * 
			    c1_dim1] + c1[m + (i__ - 1 + (k + jc * c1_dim3) * 
			    c1_dim2) * c1_dim1];
		    ch[m + (i__ + (k + jc * ch_dim3) * ch_dim2) * ch_dim1] = 
			    c1[m + (i__ + (k + j * c1_dim3) * c1_dim2) * 
			    c1_dim1] - c1[m + (i__ - 1 + (k + jc * c1_dim3) * 
			    c1_dim2) * c1_dim1];
/* L1029: */
		}
/* L129: */
	    }
/* L130: */
	}
/* L131: */
    }
L132:
    if (*ido == 1) {
	return 0;
    }
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    c2[m + (ik + c2_dim2) * c2_dim1] = ch2[m + (ik + ch2_dim2) * 
		    ch2_dim1];
/* L1033: */
	}
/* L133: */
    }
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		c1[m + ((k + j * c1_dim3) * c1_dim2 + 1) * c1_dim1] = ch[m + (
			(k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1];
/* L1034: */
	    }
/* L134: */
	}
/* L135: */
    }
    if (nbd > *l1) {
	goto L139;
    }
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	idij = is;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    c1[m + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) * c1_dim1] 
			    = wa[idij - 1] * ch[m + (i__ - 1 + (k + j * 
			    ch_dim3) * ch_dim2) * ch_dim1] - wa[idij] * ch[m 
			    + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1];
		    c1[m + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1] = 
			    wa[idij - 1] * ch[m + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] + wa[idij] * ch[m + (i__ - 1 
			    + (k + j * ch_dim3) * ch_dim2) * ch_dim1];
/* L1036: */
		}
/* L136: */
	    }
/* L137: */
	}
/* L138: */
    }
    goto L143;
L139:
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = is;
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		idij += 2;
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    c1[m + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) * c1_dim1] 
			    = wa[idij - 1] * ch[m + (i__ - 1 + (k + j * 
			    ch_dim3) * ch_dim2) * ch_dim1] - wa[idij] * ch[m 
			    + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1];
		    c1[m + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1] = 
			    wa[idij - 1] * ch[m + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] + wa[idij] * ch[m + (i__ - 1 
			    + (k + j * ch_dim3) * ch_dim2) * ch_dim1];
/* L1040: */
		}
/* L140: */
	    }
/* L141: */
	}
/* L142: */
    }
L143:
    return 0;
} /* hradbg_ */

/* Subroutine */ int hradb4_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2, 
	doublereal *wa3)
{
    /* System generated locals */
    integer cc_dim1, cc_dim2, cc_offset, ch_dim1, ch_dim2, ch_dim3, ch_offset,
	     i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, k, m, ic, idp2;
    static doublereal sqrt2;


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa3;
    --wa2;
    --wa1;
    cc_dim1 = *mdimcc;
    cc_dim2 = *ido;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * 5);
    cc -= cc_offset;
    ch_dim1 = *mdimch;
    ch_dim2 = *ido;
    ch_dim3 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * (1 + ch_dim3));
    ch -= ch_offset;

    /* Function Body */
    sqrt2 = sqrt(2.f);
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + ((k + ch_dim3 * 3) * ch_dim2 + 1) * ch_dim1] = cc[m + (((k 
		    << 2) + 1) * cc_dim2 + 1) * cc_dim1] + cc[m + (*ido + ((k 
		    << 2) + 4) * cc_dim2) * cc_dim1] - (cc[m + (*ido + ((k << 
		    2) + 2) * cc_dim2) * cc_dim1] + cc[m + (*ido + ((k << 2) 
		    + 2) * cc_dim2) * cc_dim1]);
	    ch[m + ((k + ch_dim3) * ch_dim2 + 1) * ch_dim1] = cc[m + (((k << 
		    2) + 1) * cc_dim2 + 1) * cc_dim1] + cc[m + (*ido + ((k << 
		    2) + 4) * cc_dim2) * cc_dim1] + (cc[m + (*ido + ((k << 2) 
		    + 2) * cc_dim2) * cc_dim1] + cc[m + (*ido + ((k << 2) + 2)
		     * cc_dim2) * cc_dim1]);
	    ch[m + ((k + (ch_dim3 << 2)) * ch_dim2 + 1) * ch_dim1] = cc[m + ((
		    (k << 2) + 1) * cc_dim2 + 1) * cc_dim1] - cc[m + (*ido + (
		    (k << 2) + 4) * cc_dim2) * cc_dim1] + (cc[m + (((k << 2) 
		    + 3) * cc_dim2 + 1) * cc_dim1] + cc[m + (((k << 2) + 3) * 
		    cc_dim2 + 1) * cc_dim1]);
	    ch[m + ((k + (ch_dim3 << 1)) * ch_dim2 + 1) * ch_dim1] = cc[m + ((
		    (k << 2) + 1) * cc_dim2 + 1) * cc_dim1] - cc[m + (*ido + (
		    (k << 2) + 4) * cc_dim2) * cc_dim1] - (cc[m + (((k << 2) 
		    + 3) * cc_dim2 + 1) * cc_dim1] + cc[m + (((k << 2) + 3) * 
		    cc_dim2 + 1) * cc_dim1]);
/* L1001: */
	}
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + (i__ - 1 + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m 
			+ (i__ - 1 + ((k << 2) + 1) * cc_dim2) * cc_dim1] + 
			cc[m + (ic - 1 + ((k << 2) + 4) * cc_dim2) * cc_dim1] 
			+ (cc[m + (i__ - 1 + ((k << 2) + 3) * cc_dim2) * 
			cc_dim1] + cc[m + (ic - 1 + ((k << 2) + 2) * cc_dim2) 
			* cc_dim1]);
		ch[m + (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m + (
			i__ + ((k << 2) + 1) * cc_dim2) * cc_dim1] - cc[m + (
			ic + ((k << 2) + 4) * cc_dim2) * cc_dim1] + (cc[m + (
			i__ + ((k << 2) + 3) * cc_dim2) * cc_dim1] - cc[m + (
			ic + ((k << 2) + 2) * cc_dim2) * cc_dim1]);
		ch[m + (i__ - 1 + (k + (ch_dim3 << 1)) * ch_dim2) * ch_dim1] =
			 wa1[i__ - 2] * (cc[m + (i__ - 1 + ((k << 2) + 1) * 
			cc_dim2) * cc_dim1] - cc[m + (ic - 1 + ((k << 2) + 4) 
			* cc_dim2) * cc_dim1] - (cc[m + (i__ + ((k << 2) + 3) 
			* cc_dim2) * cc_dim1] + cc[m + (ic + ((k << 2) + 2) * 
			cc_dim2) * cc_dim1])) - wa1[i__ - 1] * (cc[m + (i__ + 
			((k << 2) + 1) * cc_dim2) * cc_dim1] + cc[m + (ic + ((
			k << 2) + 4) * cc_dim2) * cc_dim1] + (cc[m + (i__ - 1 
			+ ((k << 2) + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 
			1 + ((k << 2) + 2) * cc_dim2) * cc_dim1]));
		ch[m + (i__ + (k + (ch_dim3 << 1)) * ch_dim2) * ch_dim1] = 
			wa1[i__ - 2] * (cc[m + (i__ + ((k << 2) + 1) * 
			cc_dim2) * cc_dim1] + cc[m + (ic + ((k << 2) + 4) * 
			cc_dim2) * cc_dim1] + (cc[m + (i__ - 1 + ((k << 2) + 
			3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + ((k << 2) 
			+ 2) * cc_dim2) * cc_dim1])) + wa1[i__ - 1] * (cc[m + 
			(i__ - 1 + ((k << 2) + 1) * cc_dim2) * cc_dim1] - cc[
			m + (ic - 1 + ((k << 2) + 4) * cc_dim2) * cc_dim1] - (
			cc[m + (i__ + ((k << 2) + 3) * cc_dim2) * cc_dim1] + 
			cc[m + (ic + ((k << 2) + 2) * cc_dim2) * cc_dim1]));
		ch[m + (i__ - 1 + (k + ch_dim3 * 3) * ch_dim2) * ch_dim1] = 
			wa2[i__ - 2] * (cc[m + (i__ - 1 + ((k << 2) + 1) * 
			cc_dim2) * cc_dim1] + cc[m + (ic - 1 + ((k << 2) + 4) 
			* cc_dim2) * cc_dim1] - (cc[m + (i__ - 1 + ((k << 2) 
			+ 3) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + ((k << 
			2) + 2) * cc_dim2) * cc_dim1])) - wa2[i__ - 1] * (cc[
			m + (i__ + ((k << 2) + 1) * cc_dim2) * cc_dim1] - cc[
			m + (ic + ((k << 2) + 4) * cc_dim2) * cc_dim1] - (cc[
			m + (i__ + ((k << 2) + 3) * cc_dim2) * cc_dim1] - cc[
			m + (ic + ((k << 2) + 2) * cc_dim2) * cc_dim1]));
		ch[m + (i__ + (k + ch_dim3 * 3) * ch_dim2) * ch_dim1] = wa2[
			i__ - 2] * (cc[m + (i__ + ((k << 2) + 1) * cc_dim2) * 
			cc_dim1] - cc[m + (ic + ((k << 2) + 4) * cc_dim2) * 
			cc_dim1] - (cc[m + (i__ + ((k << 2) + 3) * cc_dim2) * 
			cc_dim1] - cc[m + (ic + ((k << 2) + 2) * cc_dim2) * 
			cc_dim1])) + wa2[i__ - 1] * (cc[m + (i__ - 1 + ((k << 
			2) + 1) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + ((k 
			<< 2) + 4) * cc_dim2) * cc_dim1] - (cc[m + (i__ - 1 + 
			((k << 2) + 3) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 
			+ ((k << 2) + 2) * cc_dim2) * cc_dim1]));
		ch[m + (i__ - 1 + (k + (ch_dim3 << 2)) * ch_dim2) * ch_dim1] =
			 wa3[i__ - 2] * (cc[m + (i__ - 1 + ((k << 2) + 1) * 
			cc_dim2) * cc_dim1] - cc[m + (ic - 1 + ((k << 2) + 4) 
			* cc_dim2) * cc_dim1] + (cc[m + (i__ + ((k << 2) + 3) 
			* cc_dim2) * cc_dim1] + cc[m + (ic + ((k << 2) + 2) * 
			cc_dim2) * cc_dim1])) - wa3[i__ - 1] * (cc[m + (i__ + 
			((k << 2) + 1) * cc_dim2) * cc_dim1] + cc[m + (ic + ((
			k << 2) + 4) * cc_dim2) * cc_dim1] - (cc[m + (i__ - 1 
			+ ((k << 2) + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 
			1 + ((k << 2) + 2) * cc_dim2) * cc_dim1]));
		ch[m + (i__ + (k + (ch_dim3 << 2)) * ch_dim2) * ch_dim1] = 
			wa3[i__ - 2] * (cc[m + (i__ + ((k << 2) + 1) * 
			cc_dim2) * cc_dim1] + cc[m + (ic + ((k << 2) + 4) * 
			cc_dim2) * cc_dim1] - (cc[m + (i__ - 1 + ((k << 2) + 
			3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + ((k << 2) 
			+ 2) * cc_dim2) * cc_dim1])) + wa3[i__ - 1] * (cc[m + 
			(i__ - 1 + ((k << 2) + 1) * cc_dim2) * cc_dim1] - cc[
			m + (ic - 1 + ((k << 2) + 4) * cc_dim2) * cc_dim1] + (
			cc[m + (i__ + ((k << 2) + 3) * cc_dim2) * cc_dim1] + 
			cc[m + (ic + ((k << 2) + 2) * cc_dim2) * cc_dim1]));
/* L1002: */
	    }
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + (*ido + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m + (*ido 
		    + ((k << 2) + 1) * cc_dim2) * cc_dim1] + cc[m + (*ido + ((
		    k << 2) + 3) * cc_dim2) * cc_dim1] + (cc[m + (*ido + ((k 
		    << 2) + 1) * cc_dim2) * cc_dim1] + cc[m + (*ido + ((k << 
		    2) + 3) * cc_dim2) * cc_dim1]);
	    ch[m + (*ido + (k + (ch_dim3 << 1)) * ch_dim2) * ch_dim1] = sqrt2 
		    * (cc[m + (*ido + ((k << 2) + 1) * cc_dim2) * cc_dim1] - 
		    cc[m + (*ido + ((k << 2) + 3) * cc_dim2) * cc_dim1] - (cc[
		    m + (((k << 2) + 2) * cc_dim2 + 1) * cc_dim1] + cc[m + (((
		    k << 2) + 4) * cc_dim2 + 1) * cc_dim1]));
	    ch[m + (*ido + (k + ch_dim3 * 3) * ch_dim2) * ch_dim1] = cc[m + ((
		    (k << 2) + 4) * cc_dim2 + 1) * cc_dim1] - cc[m + (((k << 
		    2) + 2) * cc_dim2 + 1) * cc_dim1] + (cc[m + (((k << 2) + 
		    4) * cc_dim2 + 1) * cc_dim1] - cc[m + (((k << 2) + 2) * 
		    cc_dim2 + 1) * cc_dim1]);
	    ch[m + (*ido + (k + (ch_dim3 << 2)) * ch_dim2) * ch_dim1] = 
		    -sqrt2 * (cc[m + (*ido + ((k << 2) + 1) * cc_dim2) * 
		    cc_dim1] - cc[m + (*ido + ((k << 2) + 3) * cc_dim2) * 
		    cc_dim1] + (cc[m + (((k << 2) + 2) * cc_dim2 + 1) * 
		    cc_dim1] + cc[m + (((k << 2) + 4) * cc_dim2 + 1) * 
		    cc_dim1]));
/* L1003: */
	}
/* L106: */
    }
L107:
    return 0;
} /* hradb4_ */

/* Subroutine */ int hradb2_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1)
{
    /* System generated locals */
    integer cc_dim1, cc_dim2, cc_offset, ch_dim1, ch_dim2, ch_dim3, ch_offset,
	     i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, m, ic, idp2;


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa1;
    cc_dim1 = *mdimcc;
    cc_dim2 = *ido;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * 3);
    cc -= cc_offset;
    ch_dim1 = *mdimch;
    ch_dim2 = *ido;
    ch_dim3 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * (1 + ch_dim3));
    ch -= ch_offset;

    /* Function Body */
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + ((k + ch_dim3) * ch_dim2 + 1) * ch_dim1] = cc[m + (((k << 
		    1) + 1) * cc_dim2 + 1) * cc_dim1] + cc[m + (*ido + ((k << 
		    1) + 2) * cc_dim2) * cc_dim1];
	    ch[m + ((k + (ch_dim3 << 1)) * ch_dim2 + 1) * ch_dim1] = cc[m + ((
		    (k << 1) + 1) * cc_dim2 + 1) * cc_dim1] - cc[m + (*ido + (
		    (k << 1) + 2) * cc_dim2) * cc_dim1];
/* L1001: */
	}
/* L101: */
    }
    if ((i__1 = *ido - 2) < 0) {
	goto L107;
    } else if (i__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + (i__ - 1 + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m 
			+ (i__ - 1 + ((k << 1) + 1) * cc_dim2) * cc_dim1] + 
			cc[m + (ic - 1 + ((k << 1) + 2) * cc_dim2) * cc_dim1];
		ch[m + (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m + (
			i__ + ((k << 1) + 1) * cc_dim2) * cc_dim1] - cc[m + (
			ic + ((k << 1) + 2) * cc_dim2) * cc_dim1];
		ch[m + (i__ - 1 + (k + (ch_dim3 << 1)) * ch_dim2) * ch_dim1] =
			 wa1[i__ - 2] * (cc[m + (i__ - 1 + ((k << 1) + 1) * 
			cc_dim2) * cc_dim1] - cc[m + (ic - 1 + ((k << 1) + 2) 
			* cc_dim2) * cc_dim1]) - wa1[i__ - 1] * (cc[m + (i__ 
			+ ((k << 1) + 1) * cc_dim2) * cc_dim1] + cc[m + (ic + 
			((k << 1) + 2) * cc_dim2) * cc_dim1]);
		ch[m + (i__ + (k + (ch_dim3 << 1)) * ch_dim2) * ch_dim1] = 
			wa1[i__ - 2] * (cc[m + (i__ + ((k << 1) + 1) * 
			cc_dim2) * cc_dim1] + cc[m + (ic + ((k << 1) + 2) * 
			cc_dim2) * cc_dim1]) + wa1[i__ - 1] * (cc[m + (i__ - 
			1 + ((k << 1) + 1) * cc_dim2) * cc_dim1] - cc[m + (ic 
			- 1 + ((k << 1) + 2) * cc_dim2) * cc_dim1]);
/* L1002: */
	    }
/* L103: */
	}
/* L104: */
    }
    if (*ido % 2 == 1) {
	return 0;
    }
L105:
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + (*ido + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m + (*ido 
		    + ((k << 1) + 1) * cc_dim2) * cc_dim1] + cc[m + (*ido + ((
		    k << 1) + 1) * cc_dim2) * cc_dim1];
	    ch[m + (*ido + (k + (ch_dim3 << 1)) * ch_dim2) * ch_dim1] = -(cc[
		    m + (((k << 1) + 2) * cc_dim2 + 1) * cc_dim1] + cc[m + (((
		    k << 1) + 2) * cc_dim2 + 1) * cc_dim1]);
/* L1003: */
	}
/* L106: */
    }
L107:
    return 0;
} /* hradb2_ */

/* Subroutine */ int hradb3_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2)
{
    /* System generated locals */
    integer cc_dim1, cc_dim2, cc_offset, ch_dim1, ch_dim2, ch_dim3, ch_offset,
	     i__1, i__2, i__3;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, k, m, ic;
    static doublereal arg;
    static integer idp2;
    static doublereal taui, taur;
    extern doublereal pimach_(void);


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa2;
    --wa1;
    cc_dim1 = *mdimcc;
    cc_dim2 = *ido;
    cc_offset = 1 + cc_dim1 * (1 + (cc_dim2 << 2));
    cc -= cc_offset;
    ch_dim1 = *mdimch;
    ch_dim2 = *ido;
    ch_dim3 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * (1 + ch_dim3));
    ch -= ch_offset;

    /* Function Body */
    arg = pimach_() * 2.f / 3.f;
    taur = cos(arg);
    taui = sin(arg);
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + ((k + ch_dim3) * ch_dim2 + 1) * ch_dim1] = cc[m + ((k * 3 
		    + 1) * cc_dim2 + 1) * cc_dim1] + cc[m + (*ido + (k * 3 + 
		    2) * cc_dim2) * cc_dim1] * 2.f;
	    ch[m + ((k + (ch_dim3 << 1)) * ch_dim2 + 1) * ch_dim1] = cc[m + ((
		    k * 3 + 1) * cc_dim2 + 1) * cc_dim1] + taur * 2.f * cc[m 
		    + (*ido + (k * 3 + 2) * cc_dim2) * cc_dim1] - taui * 2.f *
		     cc[m + ((k * 3 + 3) * cc_dim2 + 1) * cc_dim1];
	    ch[m + ((k + ch_dim3 * 3) * ch_dim2 + 1) * ch_dim1] = cc[m + ((k *
		     3 + 1) * cc_dim2 + 1) * cc_dim1] + taur * 2.f * cc[m + (*
		    ido + (k * 3 + 2) * cc_dim2) * cc_dim1] + taui * 2.f * cc[
		    m + ((k * 3 + 3) * cc_dim2 + 1) * cc_dim1];
/* L1001: */
	}
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + (i__ - 1 + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m 
			+ (i__ - 1 + (k * 3 + 1) * cc_dim2) * cc_dim1] + (cc[
			m + (i__ - 1 + (k * 3 + 3) * cc_dim2) * cc_dim1] + cc[
			m + (ic - 1 + (k * 3 + 2) * cc_dim2) * cc_dim1]);
		ch[m + (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m + (
			i__ + (k * 3 + 1) * cc_dim2) * cc_dim1] + (cc[m + (
			i__ + (k * 3 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic 
			+ (k * 3 + 2) * cc_dim2) * cc_dim1]);
		ch[m + (i__ - 1 + (k + (ch_dim3 << 1)) * ch_dim2) * ch_dim1] =
			 wa1[i__ - 2] * (cc[m + (i__ - 1 + (k * 3 + 1) * 
			cc_dim2) * cc_dim1] + taur * (cc[m + (i__ - 1 + (k * 
			3 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 
			3 + 2) * cc_dim2) * cc_dim1]) - taui * (cc[m + (i__ + 
			(k * 3 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic + (k * 
			3 + 2) * cc_dim2) * cc_dim1])) - wa1[i__ - 1] * (cc[m 
			+ (i__ + (k * 3 + 1) * cc_dim2) * cc_dim1] + taur * (
			cc[m + (i__ + (k * 3 + 3) * cc_dim2) * cc_dim1] - cc[
			m + (ic + (k * 3 + 2) * cc_dim2) * cc_dim1]) + taui * 
			(cc[m + (i__ - 1 + (k * 3 + 3) * cc_dim2) * cc_dim1] 
			- cc[m + (ic - 1 + (k * 3 + 2) * cc_dim2) * cc_dim1]))
			;
		ch[m + (i__ + (k + (ch_dim3 << 1)) * ch_dim2) * ch_dim1] = 
			wa1[i__ - 2] * (cc[m + (i__ + (k * 3 + 1) * cc_dim2) *
			 cc_dim1] + taur * (cc[m + (i__ + (k * 3 + 3) * 
			cc_dim2) * cc_dim1] - cc[m + (ic + (k * 3 + 2) * 
			cc_dim2) * cc_dim1]) + taui * (cc[m + (i__ - 1 + (k * 
			3 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + (k * 
			3 + 2) * cc_dim2) * cc_dim1])) + wa1[i__ - 1] * (cc[m 
			+ (i__ - 1 + (k * 3 + 1) * cc_dim2) * cc_dim1] + taur 
			* (cc[m + (i__ - 1 + (k * 3 + 3) * cc_dim2) * cc_dim1]
			 + cc[m + (ic - 1 + (k * 3 + 2) * cc_dim2) * cc_dim1])
			 - taui * (cc[m + (i__ + (k * 3 + 3) * cc_dim2) * 
			cc_dim1] + cc[m + (ic + (k * 3 + 2) * cc_dim2) * 
			cc_dim1]));
		ch[m + (i__ - 1 + (k + ch_dim3 * 3) * ch_dim2) * ch_dim1] = 
			wa2[i__ - 2] * (cc[m + (i__ - 1 + (k * 3 + 1) * 
			cc_dim2) * cc_dim1] + taur * (cc[m + (i__ - 1 + (k * 
			3 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 
			3 + 2) * cc_dim2) * cc_dim1]) + taui * (cc[m + (i__ + 
			(k * 3 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic + (k * 
			3 + 2) * cc_dim2) * cc_dim1])) - wa2[i__ - 1] * (cc[m 
			+ (i__ + (k * 3 + 1) * cc_dim2) * cc_dim1] + taur * (
			cc[m + (i__ + (k * 3 + 3) * cc_dim2) * cc_dim1] - cc[
			m + (ic + (k * 3 + 2) * cc_dim2) * cc_dim1]) - taui * 
			(cc[m + (i__ - 1 + (k * 3 + 3) * cc_dim2) * cc_dim1] 
			- cc[m + (ic - 1 + (k * 3 + 2) * cc_dim2) * cc_dim1]))
			;
		ch[m + (i__ + (k + ch_dim3 * 3) * ch_dim2) * ch_dim1] = wa2[
			i__ - 2] * (cc[m + (i__ + (k * 3 + 1) * cc_dim2) * 
			cc_dim1] + taur * (cc[m + (i__ + (k * 3 + 3) * 
			cc_dim2) * cc_dim1] - cc[m + (ic + (k * 3 + 2) * 
			cc_dim2) * cc_dim1]) - taui * (cc[m + (i__ - 1 + (k * 
			3 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + (k * 
			3 + 2) * cc_dim2) * cc_dim1])) + wa2[i__ - 1] * (cc[m 
			+ (i__ - 1 + (k * 3 + 1) * cc_dim2) * cc_dim1] + taur 
			* (cc[m + (i__ - 1 + (k * 3 + 3) * cc_dim2) * cc_dim1]
			 + cc[m + (ic - 1 + (k * 3 + 2) * cc_dim2) * cc_dim1])
			 + taui * (cc[m + (i__ + (k * 3 + 3) * cc_dim2) * 
			cc_dim1] + cc[m + (ic + (k * 3 + 2) * cc_dim2) * 
			cc_dim1]));
/* L1002: */
	    }
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* hradb3_ */

/* Subroutine */ int hradb5_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2, 
	doublereal *wa3, doublereal *wa4)
{
    /* System generated locals */
    integer cc_dim1, cc_dim2, cc_offset, ch_dim1, ch_dim2, ch_dim3, ch_offset,
	     i__1, i__2, i__3;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, k, m, ic;
    static doublereal arg, ti11, ti12, tr11, tr12;
    static integer idp2;
    extern doublereal pimach_(void);


/*     a multiple fft package for spherepack */

    /* Parameter adjustments */
    --wa4;
    --wa3;
    --wa2;
    --wa1;
    cc_dim1 = *mdimcc;
    cc_dim2 = *ido;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * 6);
    cc -= cc_offset;
    ch_dim1 = *mdimch;
    ch_dim2 = *ido;
    ch_dim3 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * (1 + ch_dim3));
    ch -= ch_offset;

    /* Function Body */
    arg = pimach_() * 2.f / 5.f;
    tr11 = cos(arg);
    ti11 = sin(arg);
    tr12 = cos(arg * 2.f);
    ti12 = sin(arg * 2.f);
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + ((k + ch_dim3) * ch_dim2 + 1) * ch_dim1] = cc[m + ((k * 5 
		    + 1) * cc_dim2 + 1) * cc_dim1] + cc[m + (*ido + (k * 5 + 
		    2) * cc_dim2) * cc_dim1] * 2.f + cc[m + (*ido + (k * 5 + 
		    4) * cc_dim2) * cc_dim1] * 2.f;
	    ch[m + ((k + (ch_dim3 << 1)) * ch_dim2 + 1) * ch_dim1] = cc[m + ((
		    k * 5 + 1) * cc_dim2 + 1) * cc_dim1] + tr11 * 2.f * cc[m 
		    + (*ido + (k * 5 + 2) * cc_dim2) * cc_dim1] + tr12 * 2.f *
		     cc[m + (*ido + (k * 5 + 4) * cc_dim2) * cc_dim1] - (ti11 
		    * 2.f * cc[m + ((k * 5 + 3) * cc_dim2 + 1) * cc_dim1] + 
		    ti12 * 2.f * cc[m + ((k * 5 + 5) * cc_dim2 + 1) * cc_dim1]
		    );
	    ch[m + ((k + ch_dim3 * 3) * ch_dim2 + 1) * ch_dim1] = cc[m + ((k *
		     5 + 1) * cc_dim2 + 1) * cc_dim1] + tr12 * 2.f * cc[m + (*
		    ido + (k * 5 + 2) * cc_dim2) * cc_dim1] + tr11 * 2.f * cc[
		    m + (*ido + (k * 5 + 4) * cc_dim2) * cc_dim1] - (ti12 * 
		    2.f * cc[m + ((k * 5 + 3) * cc_dim2 + 1) * cc_dim1] - 
		    ti11 * 2.f * cc[m + ((k * 5 + 5) * cc_dim2 + 1) * cc_dim1]
		    );
	    ch[m + ((k + (ch_dim3 << 2)) * ch_dim2 + 1) * ch_dim1] = cc[m + ((
		    k * 5 + 1) * cc_dim2 + 1) * cc_dim1] + tr12 * 2.f * cc[m 
		    + (*ido + (k * 5 + 2) * cc_dim2) * cc_dim1] + tr11 * 2.f *
		     cc[m + (*ido + (k * 5 + 4) * cc_dim2) * cc_dim1] + (ti12 
		    * 2.f * cc[m + ((k * 5 + 3) * cc_dim2 + 1) * cc_dim1] - 
		    ti11 * 2.f * cc[m + ((k * 5 + 5) * cc_dim2 + 1) * cc_dim1]
		    );
	    ch[m + ((k + ch_dim3 * 5) * ch_dim2 + 1) * ch_dim1] = cc[m + ((k *
		     5 + 1) * cc_dim2 + 1) * cc_dim1] + tr11 * 2.f * cc[m + (*
		    ido + (k * 5 + 2) * cc_dim2) * cc_dim1] + tr12 * 2.f * cc[
		    m + (*ido + (k * 5 + 4) * cc_dim2) * cc_dim1] + (ti11 * 
		    2.f * cc[m + ((k * 5 + 3) * cc_dim2 + 1) * cc_dim1] + 
		    ti12 * 2.f * cc[m + ((k * 5 + 5) * cc_dim2 + 1) * cc_dim1]
		    );
/* L1001: */
	}
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + (i__ - 1 + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m 
			+ (i__ - 1 + (k * 5 + 1) * cc_dim2) * cc_dim1] + (cc[
			m + (i__ - 1 + (k * 5 + 3) * cc_dim2) * cc_dim1] + cc[
			m + (ic - 1 + (k * 5 + 2) * cc_dim2) * cc_dim1]) + (
			cc[m + (i__ - 1 + (k * 5 + 5) * cc_dim2) * cc_dim1] + 
			cc[m + (ic - 1 + (k * 5 + 4) * cc_dim2) * cc_dim1]);
		ch[m + (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1] = cc[m + (
			i__ + (k * 5 + 1) * cc_dim2) * cc_dim1] + (cc[m + (
			i__ + (k * 5 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic 
			+ (k * 5 + 2) * cc_dim2) * cc_dim1]) + (cc[m + (i__ + 
			(k * 5 + 5) * cc_dim2) * cc_dim1] - cc[m + (ic + (k * 
			5 + 4) * cc_dim2) * cc_dim1]);
		ch[m + (i__ - 1 + (k + (ch_dim3 << 1)) * ch_dim2) * ch_dim1] =
			 wa1[i__ - 2] * (cc[m + (i__ - 1 + (k * 5 + 1) * 
			cc_dim2) * cc_dim1] + tr11 * (cc[m + (i__ - 1 + (k * 
			5 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) + tr12 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1]) - (ti11 * (cc[
			m + (i__ + (k * 5 + 3) * cc_dim2) * cc_dim1] + cc[m + 
			(ic + (k * 5 + 2) * cc_dim2) * cc_dim1]) + ti12 * (cc[
			m + (i__ + (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + 
			(ic + (k * 5 + 4) * cc_dim2) * cc_dim1]))) - wa1[i__ 
			- 1] * (cc[m + (i__ + (k * 5 + 1) * cc_dim2) * 
			cc_dim1] + tr11 * (cc[m + (i__ + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr12 * (cc[m + (i__ + (k * 5 + 
			5) * cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 4) * 
			cc_dim2) * cc_dim1]) + (ti11 * (cc[m + (i__ - 1 + (k *
			 5 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) + ti12 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] - cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1])));
		ch[m + (i__ + (k + (ch_dim3 << 1)) * ch_dim2) * ch_dim1] = 
			wa1[i__ - 2] * (cc[m + (i__ + (k * 5 + 1) * cc_dim2) *
			 cc_dim1] + tr11 * (cc[m + (i__ + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr12 * (cc[m + (i__ + (k * 5 + 
			5) * cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 4) * 
			cc_dim2) * cc_dim1]) + (ti11 * (cc[m + (i__ - 1 + (k *
			 5 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) + ti12 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] - cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1]))) + wa1[i__ - 
			1] * (cc[m + (i__ - 1 + (k * 5 + 1) * cc_dim2) * 
			cc_dim1] + tr11 * (cc[m + (i__ - 1 + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr12 * (cc[m + (i__ - 1 + (k * 
			5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 
			5 + 4) * cc_dim2) * cc_dim1]) - (ti11 * (cc[m + (i__ 
			+ (k * 5 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic + (k 
			* 5 + 2) * cc_dim2) * cc_dim1]) + ti12 * (cc[m + (i__ 
			+ (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic + (k 
			* 5 + 4) * cc_dim2) * cc_dim1])));
		ch[m + (i__ - 1 + (k + ch_dim3 * 3) * ch_dim2) * ch_dim1] = 
			wa2[i__ - 2] * (cc[m + (i__ - 1 + (k * 5 + 1) * 
			cc_dim2) * cc_dim1] + tr12 * (cc[m + (i__ - 1 + (k * 
			5 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) + tr11 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1]) - (ti12 * (cc[
			m + (i__ + (k * 5 + 3) * cc_dim2) * cc_dim1] + cc[m + 
			(ic + (k * 5 + 2) * cc_dim2) * cc_dim1]) - ti11 * (cc[
			m + (i__ + (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + 
			(ic + (k * 5 + 4) * cc_dim2) * cc_dim1]))) - wa2[i__ 
			- 1] * (cc[m + (i__ + (k * 5 + 1) * cc_dim2) * 
			cc_dim1] + tr12 * (cc[m + (i__ + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr11 * (cc[m + (i__ + (k * 5 + 
			5) * cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 4) * 
			cc_dim2) * cc_dim1]) + (ti12 * (cc[m + (i__ - 1 + (k *
			 5 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) - ti11 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] - cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1])));
		ch[m + (i__ + (k + ch_dim3 * 3) * ch_dim2) * ch_dim1] = wa2[
			i__ - 2] * (cc[m + (i__ + (k * 5 + 1) * cc_dim2) * 
			cc_dim1] + tr12 * (cc[m + (i__ + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr11 * (cc[m + (i__ + (k * 5 + 
			5) * cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 4) * 
			cc_dim2) * cc_dim1]) + (ti12 * (cc[m + (i__ - 1 + (k *
			 5 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) - ti11 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] - cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1]))) + wa2[i__ - 
			1] * (cc[m + (i__ - 1 + (k * 5 + 1) * cc_dim2) * 
			cc_dim1] + tr12 * (cc[m + (i__ - 1 + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr11 * (cc[m + (i__ - 1 + (k * 
			5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 
			5 + 4) * cc_dim2) * cc_dim1]) - (ti12 * (cc[m + (i__ 
			+ (k * 5 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic + (k 
			* 5 + 2) * cc_dim2) * cc_dim1]) - ti11 * (cc[m + (i__ 
			+ (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic + (k 
			* 5 + 4) * cc_dim2) * cc_dim1])));
		ch[m + (i__ - 1 + (k + (ch_dim3 << 2)) * ch_dim2) * ch_dim1] =
			 wa3[i__ - 2] * (cc[m + (i__ - 1 + (k * 5 + 1) * 
			cc_dim2) * cc_dim1] + tr12 * (cc[m + (i__ - 1 + (k * 
			5 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) + tr11 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1]) + (ti12 * (cc[
			m + (i__ + (k * 5 + 3) * cc_dim2) * cc_dim1] + cc[m + 
			(ic + (k * 5 + 2) * cc_dim2) * cc_dim1]) - ti11 * (cc[
			m + (i__ + (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + 
			(ic + (k * 5 + 4) * cc_dim2) * cc_dim1]))) - wa3[i__ 
			- 1] * (cc[m + (i__ + (k * 5 + 1) * cc_dim2) * 
			cc_dim1] + tr12 * (cc[m + (i__ + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr11 * (cc[m + (i__ + (k * 5 + 
			5) * cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 4) * 
			cc_dim2) * cc_dim1]) - (ti12 * (cc[m + (i__ - 1 + (k *
			 5 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) - ti11 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] - cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1])));
		ch[m + (i__ + (k + (ch_dim3 << 2)) * ch_dim2) * ch_dim1] = 
			wa3[i__ - 2] * (cc[m + (i__ + (k * 5 + 1) * cc_dim2) *
			 cc_dim1] + tr12 * (cc[m + (i__ + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr11 * (cc[m + (i__ + (k * 5 + 
			5) * cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 4) * 
			cc_dim2) * cc_dim1]) - (ti12 * (cc[m + (i__ - 1 + (k *
			 5 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) - ti11 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] - cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1]))) + wa3[i__ - 
			1] * (cc[m + (i__ - 1 + (k * 5 + 1) * cc_dim2) * 
			cc_dim1] + tr12 * (cc[m + (i__ - 1 + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr11 * (cc[m + (i__ - 1 + (k * 
			5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 
			5 + 4) * cc_dim2) * cc_dim1]) + (ti12 * (cc[m + (i__ 
			+ (k * 5 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic + (k 
			* 5 + 2) * cc_dim2) * cc_dim1]) - ti11 * (cc[m + (i__ 
			+ (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic + (k 
			* 5 + 4) * cc_dim2) * cc_dim1])));
		ch[m + (i__ - 1 + (k + ch_dim3 * 5) * ch_dim2) * ch_dim1] = 
			wa4[i__ - 2] * (cc[m + (i__ - 1 + (k * 5 + 1) * 
			cc_dim2) * cc_dim1] + tr11 * (cc[m + (i__ - 1 + (k * 
			5 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) + tr12 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1]) + (ti11 * (cc[
			m + (i__ + (k * 5 + 3) * cc_dim2) * cc_dim1] + cc[m + 
			(ic + (k * 5 + 2) * cc_dim2) * cc_dim1]) + ti12 * (cc[
			m + (i__ + (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + 
			(ic + (k * 5 + 4) * cc_dim2) * cc_dim1]))) - wa4[i__ 
			- 1] * (cc[m + (i__ + (k * 5 + 1) * cc_dim2) * 
			cc_dim1] + tr11 * (cc[m + (i__ + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr12 * (cc[m + (i__ + (k * 5 + 
			5) * cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 4) * 
			cc_dim2) * cc_dim1]) - (ti11 * (cc[m + (i__ - 1 + (k *
			 5 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) + ti12 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] - cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1])));
		ch[m + (i__ + (k + ch_dim3 * 5) * ch_dim2) * ch_dim1] = wa4[
			i__ - 2] * (cc[m + (i__ + (k * 5 + 1) * cc_dim2) * 
			cc_dim1] + tr11 * (cc[m + (i__ + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr12 * (cc[m + (i__ + (k * 5 + 
			5) * cc_dim2) * cc_dim1] - cc[m + (ic + (k * 5 + 4) * 
			cc_dim2) * cc_dim1]) - (ti11 * (cc[m + (i__ - 1 + (k *
			 5 + 3) * cc_dim2) * cc_dim1] - cc[m + (ic - 1 + (k * 
			5 + 2) * cc_dim2) * cc_dim1]) + ti12 * (cc[m + (i__ - 
			1 + (k * 5 + 5) * cc_dim2) * cc_dim1] - cc[m + (ic - 
			1 + (k * 5 + 4) * cc_dim2) * cc_dim1]))) + wa4[i__ - 
			1] * (cc[m + (i__ - 1 + (k * 5 + 1) * cc_dim2) * 
			cc_dim1] + tr11 * (cc[m + (i__ - 1 + (k * 5 + 3) * 
			cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 5 + 2) * 
			cc_dim2) * cc_dim1]) + tr12 * (cc[m + (i__ - 1 + (k * 
			5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic - 1 + (k * 
			5 + 4) * cc_dim2) * cc_dim1]) + (ti11 * (cc[m + (i__ 
			+ (k * 5 + 3) * cc_dim2) * cc_dim1] + cc[m + (ic + (k 
			* 5 + 2) * cc_dim2) * cc_dim1]) + ti12 * (cc[m + (i__ 
			+ (k * 5 + 5) * cc_dim2) * cc_dim1] + cc[m + (ic + (k 
			* 5 + 4) * cc_dim2) * cc_dim1])));
/* L1002: */
	    }
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* hradb5_ */

