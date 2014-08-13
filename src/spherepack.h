
#ifndef SPHEREPACK_H
#define SPHEREPACK_H

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif
/* Subroutine */ int alfk_(integer *n, integer *m, doublereal *cp);

/* Subroutine */ int lfim_(integer *init, doublereal *theta, integer *l, integer *n,
	 integer *nm, doublereal *pb, integer *id, doublereal *wlfim);

/* Subroutine */ int lfim1_(integer *init, doublereal *theta, integer *l, integer *
	n, integer *nm, integer *id, doublereal *p3, doublereal *phz, doublereal *ph1, doublereal *p1,
	 doublereal *p2, doublereal *cp);

/* Subroutine */ int lfin_(integer *init, doublereal *theta, integer *l, integer *m,
	 integer *nm, doublereal *pb, integer *id, doublereal *wlfin);

/* Subroutine */ int lfin1_(integer *init, doublereal *theta, integer *l, integer *
	m, integer *nm, integer *id, doublereal *p3, doublereal *phz, doublereal *ph1, doublereal *p1,
	 doublereal *p2, doublereal *cp);

/* Subroutine */ int lfpt_(integer *n, integer *m, doublereal *theta, doublereal *cp, 
	doublereal *pb);

/* Subroutine */ int divec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *dv, integer *idv, integer *jdv, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *wshsec, integer *lshsec, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int divec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *dv, integer *idv, integer *jdv, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *a, doublereal *b, integer *mab, doublereal *
	sqnn, doublereal *wshsec, integer *lshsec, doublereal *wk, integer *lwk, integer *
	ierror);

/* Subroutine */ int dives_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *dv, integer *idv, integer *jdv, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *wshses, integer *lshses, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int dives1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *dv, integer *idv, integer *jdv, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *a, doublereal *b, integer *mab, doublereal *
	sqnn, doublereal *wshses, integer *lshses, doublereal *wk, integer *lwk, integer *
	ierror);

/* Subroutine */ int divgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *dv, integer *idv, integer *jdv, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *wshsgc, integer *lshsgc, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int divgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *dv, integer *idv, integer *jdv, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *a, doublereal *b, integer *mab, doublereal *
	sqnn, doublereal *wshsgc, integer *lshsgc, doublereal *wk, integer *lwk, integer *
	ierror);

/* Subroutine */ int divgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *divg, integer *idiv, integer *jdiv, doublereal *br, doublereal 
	*bi, integer *mdb, integer *ndb, doublereal *wshsgs, integer *lshsgs, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int divgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *divg, integer *idiv, integer *jdiv, doublereal *br, doublereal 
	*bi, integer *mdb, integer *ndb, doublereal *a, doublereal *b, integer *mab, doublereal 
	*sqnn, doublereal *wshsgs, integer *lshsgs, doublereal *wk, integer *lwk, integer 
	*ierror);

/* Subroutine */ int gaqd_(integer *nlat, doublereal *theta, doublereal *wts, 
	doublereal *w, integer *lwork, integer *ierror);

/* Subroutine */ int cpdp_(integer *n, doublereal *cz, doublereal *cp, 
	doublereal *dcp);

/* Subroutine */ int tpdp_(integer *n, doublereal *theta, doublereal *cz, 
	doublereal *cp, doublereal *dcp, doublereal *pb, doublereal *dpb);

/* Subroutine */ int geo2maths_(integer *ig, integer *nlon, integer *nlat, 
	doublereal *sg, doublereal *sm, doublereal *work);

/* Subroutine */ int math2geos_(integer *ig, integer *nlat, integer *nlon, 
	doublereal *sm, doublereal *sg, doublereal *work);

/* Subroutine */ int geo2mathv_(integer *ig, integer *nlon, integer *nlat, 
	doublereal *ug, doublereal *vg, doublereal *vm, doublereal *wm, doublereal *work);

/* Subroutine */ int math2geov_(integer *ig, integer *nlat, integer *nlon, 
	doublereal *vm, doublereal *wm, doublereal *ug, doublereal *vg, doublereal *work);

/* Subroutine */ int gradec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhsec, integer *lvhsec, 
	doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int gradec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wvhsec, integer *lvhsec, doublereal *wk, integer *
	lwk, integer *ierror);

/* Subroutine */ int grades_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhses, integer *lvhses, 
	doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int grades1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wvhses, integer *lvhses, doublereal *wk, integer *
	lwk, integer *ierror);

/* Subroutine */ int gradgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhsgc, integer *lvhsgc, 
	doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int gradgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wvhsgc, integer *lvhsgc, doublereal *wk, integer *
	lwk, integer *ierror);

/* Subroutine */ int gradgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhsgs, integer *lvhsgs, 
	doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int gradgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wvhsgs, integer *lvhsgs, doublereal *wk, integer *
	lwk, integer *ierror);

/* Subroutine */ int hrffti_(integer *n, doublereal *wsave);

/* Subroutine */ int hrfti1_(integer *n, doublereal *wa, doublereal *fac);

/* Subroutine */ int hrfftf_(integer *m, integer *n, doublereal *r__, integer *
	mdimr, doublereal *whrfft, doublereal *work);

/* Subroutine */ int hrftf1_(integer *m, integer *n, doublereal *c__, integer *
	mdimc, doublereal *ch, doublereal *wa, doublereal *fac);

/* Subroutine */ int hradf4_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2, 
	doublereal *wa3);

/* Subroutine */ int hradf2_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1);

/* Subroutine */ int hradf3_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2);

/* Subroutine */ int hradf5_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2, 
	doublereal *wa3, doublereal *wa4);

/* Subroutine */ int hradfg_(integer *mp, integer *ido, integer *ip, integer *
	l1, integer *idl1, doublereal *cc, doublereal *c1, doublereal *c2, integer *mdimcc, 
	doublereal *ch, doublereal *ch2, integer *mdimch, doublereal *wa);

/* Subroutine */ int hrfftb_(integer *m, integer *n, doublereal *r__, integer *
	mdimr, doublereal *whrfft, doublereal *work);

/* Subroutine */ int hrftb1_(integer *m, integer *n, doublereal *c__, integer *
	mdimc, doublereal *ch, doublereal *wa, doublereal *fac);

/* Subroutine */ int hradbg_(integer *mp, integer *ido, integer *ip, integer *
	l1, integer *idl1, doublereal *cc, doublereal *c1, doublereal *c2, integer *mdimcc, 
	doublereal *ch, doublereal *ch2, integer *mdimch, doublereal *wa);

/* Subroutine */ int hradb4_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2, 
	doublereal *wa3);

/* Subroutine */ int hradb2_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1);

/* Subroutine */ int hradb3_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2);

/* Subroutine */ int hradb5_(integer *mp, integer *ido, integer *l1, doublereal *cc,
	 integer *mdimcc, doublereal *ch, integer *mdimch, doublereal *wa1, doublereal *wa2, 
	doublereal *wa3, doublereal *wa4);

/* Subroutine */ int idivec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhsec, integer *lvhsec, 
	doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int idvec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wvhsec, integer *lvhsec, doublereal *wk, integer *
	lwk, doublereal *pertrb, integer *ierror);

/* Subroutine */ int idives_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhses, integer *lvhses, 
	doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int idves1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wsav, integer *lwsav, doublereal *wk, integer *lwk, 
	doublereal *pertrb, integer *ierror);

/* Subroutine */ int idivgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhsgc, integer *lvhsgc, 
	doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int idvgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wsav, integer *lwsav, doublereal *wk, integer *lwk, 
	doublereal *pertrb, integer *ierror);

/* Subroutine */ int idivgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhsgs, integer *lvhsgs, 
	doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int idvgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wsav, integer *lwsav, doublereal *wk, integer *lwk, 
	doublereal *pertrb, integer *ierror);

/* Subroutine */ int idvtec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *ad,
	 doublereal *bd, doublereal *av, doublereal *bv, integer *mdab, integer *ndab, doublereal *
	wvhsec, integer *lvhsec, doublereal *work, integer *lwork, doublereal *pertbd, 
	doublereal *pertbv, integer *ierror);

/* Subroutine */ int idvtec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mmax, doublereal *sqnn, integer *
	mdab, integer *ndab, doublereal *ad, doublereal *bd, doublereal *av, doublereal *bv, doublereal *
	wvhsec, integer *lvhsec, doublereal *wk, integer *lwk, doublereal *pertbd, doublereal *
	pertbv, integer *ierror);

/* Subroutine */ int idvtes_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *ad,
	 doublereal *bd, doublereal *av, doublereal *bv, integer *mdab, integer *ndab, doublereal *
	wvhses, integer *lvhses, doublereal *work, integer *lwork, doublereal *pertbd, 
	doublereal *pertbv, integer *ierror);

/* Subroutine */ int idvtes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mmax, doublereal *sqnn, integer *
	mdab, integer *ndab, doublereal *ad, doublereal *bd, doublereal *av, doublereal *bv, doublereal *
	widvtes, integer *lidvtes, doublereal *wk, integer *lwk, doublereal *pertbd, doublereal 
	*pertbv, integer *ierror);

/* Subroutine */ int idvtgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *ad,
	 doublereal *bd, doublereal *av, doublereal *bv, integer *mdab, integer *ndab, doublereal *
	wvhsgc, integer *lvhsgc, doublereal *work, integer *lwork, doublereal *pertbd, 
	doublereal *pertbv, integer *ierror);

/* Subroutine */ int idvtgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mmax, doublereal *sqnn, integer *
	mdab, integer *ndab, doublereal *ad, doublereal *bd, doublereal *av, doublereal *bv, doublereal *
	wvhsgc, integer *lvhsgc, doublereal *wk, integer *lwk, doublereal *pertbd, doublereal *
	pertbv, integer *ierror);

/* Subroutine */ int idvtgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *ad,
	 doublereal *bd, doublereal *av, doublereal *bv, integer *mdab, integer *ndab, doublereal *
	wvhsgs, integer *lvhsgs, doublereal *work, integer *lwork, doublereal *pertbd, 
	doublereal *pertbv, integer *ierror);

/* Subroutine */ int idvtgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mmax, doublereal *sqnn, integer *
	mdab, integer *ndab, doublereal *ad, doublereal *bd, doublereal *av, doublereal *bv, doublereal *
	wvhsgs, integer *lvhsgs, doublereal *wk, integer *lwk, doublereal *pertbd, doublereal *
	pertbv, integer *ierror);

/* Subroutine */ int igradec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, integer *isf, integer *jsf, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *wshsec, integer *lshsec, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int igrdec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, integer *isf, integer *jsf, doublereal *a, doublereal *b, 
	integer *mab, doublereal *sqnn, integer *mdb, integer *ndb, doublereal *br, doublereal *
	bi, doublereal *wshsec, integer *lshsec, doublereal *wk, integer *lwk, integer *
	ierror);

/* Subroutine */ int igrades_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, integer *isf, integer *jsf, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *wshses, integer *lshses, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int igrdes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, integer *isf, integer *jsf, doublereal *a, doublereal *b, 
	integer *mab, doublereal *sqnn, integer *mdb, integer *ndb, doublereal *br, doublereal *
	bi, doublereal *wshses, integer *lshses, doublereal *wk, integer *lwk, integer *
	ierror);

/* Subroutine */ int igradgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, integer *isf, integer *jsf, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *wshsgc, integer *lshsgc, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int igrdgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, integer *isf, integer *jsf, doublereal *a, doublereal *b, 
	integer *mab, doublereal *sqnn, integer *mdb, integer *ndb, doublereal *br, doublereal *
	bi, doublereal *wsav, integer *lsav, doublereal *wk, integer *lwk, integer *
	ierror);

/* Subroutine */ int igradgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, integer *isf, integer *jsf, doublereal *br, doublereal *bi,
	 integer *mdb, integer *ndb, doublereal *wshsgs, integer *lshsgs, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int igrdgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, integer *isf, integer *jsf, doublereal *a, doublereal *b, 
	integer *mab, doublereal *sqnn, integer *mdb, integer *ndb, doublereal *br, doublereal *
	bi, doublereal *wsav, integer *lsav, doublereal *wk, integer *lwk, integer *
	ierror);

/* Subroutine */ int ihgeod_(integer *m, integer *idp, integer *jdp, doublereal *x, 
	doublereal *y, doublereal *z__);

/* Subroutine */ int ctos_(doublereal *x, doublereal *y, doublereal *z__, doublereal *r__, doublereal *
	theta, doublereal *phi);

/* Subroutine */ int stoc_(doublereal *r__, doublereal *theta, doublereal *phi, doublereal *x, doublereal *
	y, doublereal *z__);

/* Subroutine */ int isfvpec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idv, integer *jdv, doublereal *as, 
	doublereal *bs, doublereal *av, doublereal *bv, integer *mdb, integer *ndb, doublereal *
	wvhsec, integer *lvhsec, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int isfvpec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idv, integer *jdv, doublereal *as, 
	doublereal *bs, doublereal *av, doublereal *bv, integer *mdb, integer *ndb, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci, integer *mab, doublereal *fnn, doublereal *wvhsec, 
	integer *lvhsec, doublereal *wk, integer *lwk, integer *ierror);

/* Subroutine */ int isfvpes_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idv, integer *jdv, doublereal *as, 
	doublereal *bs, doublereal *av, doublereal *bv, integer *mdb, integer *ndb, doublereal *
	wvhses, integer *lvhses, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int isfvpes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idv, integer *jdv, doublereal *as, 
	doublereal *bs, doublereal *av, doublereal *bv, integer *mdb, integer *ndb, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci, integer *mab, doublereal *fnn, doublereal *wvhses, 
	integer *lvhses, doublereal *wk, integer *lwk, integer *ierror);

/* Subroutine */ int isfvpgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idv, integer *jdv, doublereal *as, 
	doublereal *bs, doublereal *av, doublereal *bv, integer *mdb, integer *ndb, doublereal *
	wvhsgc, integer *lvhsgc, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int isfvpgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idv, integer *jdv, doublereal *as, 
	doublereal *bs, doublereal *av, doublereal *bv, integer *mdb, integer *ndb, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci, integer *mab, doublereal *fnn, doublereal *wvhsgc, 
	integer *lvhsgc, doublereal *wk, integer *lwk, integer *ierror);

/* Subroutine */ int isfvpgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idv, integer *jdv, doublereal *as, 
	doublereal *bs, doublereal *av, doublereal *bv, integer *mdb, integer *ndb, doublereal *
	wvhsgs, integer *lvhsgs, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int isfvpgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idv, integer *jdv, doublereal *as, 
	doublereal *bs, doublereal *av, doublereal *bv, integer *mdb, integer *ndb, doublereal *br, 
	doublereal *bi, doublereal *cr, doublereal *ci, integer *mab, doublereal *fnn, doublereal *wvhsgs, 
	integer *lvhsgs, doublereal *wk, integer *lwk, integer *ierror);

/* Subroutine */ int islapec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *xlmbda, doublereal *sf, integer *ids, integer *jds, doublereal 
	*a, doublereal *b, integer *mdab, integer *ndab, doublereal *wshsec, integer *
	lshsec, doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int islpec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *xlmbda, doublereal *sf, integer *ids, integer *jds, doublereal 
	*a, doublereal *b, integer *mdab, integer *ndab, doublereal *as, doublereal *bs, 
	integer *mmax, doublereal *fnn, doublereal *wshsec, integer *lshsec, doublereal *wk, 
	integer *lwk, doublereal *pertrb, integer *ierror);

/* Subroutine */ int islapec_l3_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *mu2, doublereal *mu, 
	doublereal *xlmbda, doublereal *sf, integer *ids, 
	integer *jds, doublereal *a, doublereal *b, integer *mdab, integer *
	ndab, doublereal *wshsec, integer *lshsec, doublereal *work, integer *
	lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int islpec1_l3_(integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	    integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);

/* Subroutine */ int islapes_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *xlmbda, doublereal *sf, integer *ids, integer *jds, doublereal 
	*a, doublereal *b, integer *mdab, integer *ndab, doublereal *wshses, integer *
	lshses, doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int islpes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *xlmbda, doublereal *sf, integer *ids, integer *jds, doublereal 
	*a, doublereal *b, integer *mdab, integer *ndab, doublereal *as, doublereal *bs, 
	integer *mmax, doublereal *fnn, doublereal *wshses, integer *lshses, doublereal *wk, 
	integer *lwk, doublereal *pertrb, integer *ierror);

/* Subroutine */ int islapgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *xlmbda, doublereal *sf, integer *ids, integer *jds, doublereal 
	*a, doublereal *b, integer *mdab, integer *ndab, doublereal *wshsgc, integer *
	lshsgc, doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int islpgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *xlmbda, doublereal *sf, integer *ids, integer *jds, doublereal 
	*a, doublereal *b, integer *mdab, integer *ndab, doublereal *as, doublereal *bs, 
	integer *mmax, doublereal *fnn, doublereal *wsav, integer *lsav, doublereal *wk, 
	integer *lwk, doublereal *pertrb, integer *ierror);

/* Subroutine */ int islapgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *xlmbda, doublereal *sf, integer *ids, integer *jds, doublereal 
	*a, doublereal *b, integer *mdab, integer *ndab, doublereal *wshsgs, integer *
	lshsgs, doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int islpgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *xlmbda, doublereal *sf, integer *ids, integer *jds, doublereal 
	*a, doublereal *b, integer *mdab, integer *ndab, doublereal *as, doublereal *bs, 
	integer *mmax, doublereal *fnn, doublereal *wsav, integer *lsav, doublereal *wk, 
	integer *lwk, doublereal *pertrb, integer *ierror);

/* Subroutine */ int ivlapec_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdbc, integer *ndbc, doublereal *
	wvhsec, integer *lvhsec, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int ivlapec1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *
	brvw, doublereal *bivw, doublereal *crvw, doublereal *civw, integer *mmax, doublereal *fnn, 
	integer *mdbc, integer *ndbc, doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, 
	doublereal *wsave, integer *lwsav, doublereal *wk, integer *lwk, integer *ierror);

/* Subroutine */ int ivlapes_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdbc, integer *ndbc, doublereal *
	wvhses, integer *lvhses, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int ivlapes1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *
	brvw, doublereal *bivw, doublereal *crvw, doublereal *civw, integer *mmax, doublereal *fnn, 
	integer *mdbc, integer *ndbc, doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, 
	doublereal *wsave, integer *lsave, doublereal *wk, integer *lwk, integer *ierror);

/* Subroutine */ int ivlapgc_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdbc, integer *ndbc, doublereal *
	wvhsgc, integer *lvhsgc, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int ivlapgc1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *
	brvw, doublereal *bivw, doublereal *crvw, doublereal *civw, integer *mmax, doublereal *fnn, 
	integer *mdbc, integer *ndbc, doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, 
	doublereal *wvhsgc, integer *lvhsgc, doublereal *wk, integer *lwk, integer *
	ierror);

/* Subroutine */ int ivlapgs_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdbc, integer *ndbc, doublereal *
	wvhsgs, integer *lvhsgs, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int ivlapgs1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *
	brvw, doublereal *bivw, doublereal *crvw, doublereal *civw, integer *mmax, doublereal *fnn, 
	integer *mdbc, integer *ndbc, doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, 
	doublereal *wsave, integer *lsave, doublereal *wk, integer *lwk, integer *ierror);

/* Subroutine */ int ivrtec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhsec, integer *lvhsec, 
	doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int ivtec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *cr,
	 doublereal *ci, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wsav, integer *lwsav, doublereal *wk, integer *lwk, 
	doublereal *pertrb, integer *ierror);

/* Subroutine */ int ivrtes_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhses, integer *lvhses, 
	doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int ivtes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *cr,
	 doublereal *ci, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wsav, integer *lwsav, doublereal *wk, integer *lwk, 
	doublereal *pertrb, integer *ierror);

/* Subroutine */ int ivrtgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhsgc, integer *lvhsgc, 
	doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int ivtgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *cr,
	 doublereal *ci, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wsav, integer *lsav, doublereal *wk, integer *lwk, 
	doublereal *pertrb, integer *ierror);

/* Subroutine */ int ivrtgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *a, 
	doublereal *b, integer *mdab, integer *ndab, doublereal *wvhsgs, integer *lvhsgs, 
	doublereal *work, integer *lwork, doublereal *pertrb, integer *ierror);

/* Subroutine */ int ivtgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *cr,
	 doublereal *ci, integer *mmax, doublereal *sqnn, integer *mdab, integer *ndab, 
	doublereal *a, doublereal *b, doublereal *wsav, integer *lsav, doublereal *wk, integer *lwk, 
	doublereal *pertrb, integer *ierror);

/* Subroutine */ int lfim_(integer *init, doublereal *theta, integer *l, integer *n,
	 integer *nm, doublereal *pb, integer *id, doublereal *wlfim);

/* Subroutine */ int lfim1_(integer *init, doublereal *theta, integer *l, integer *
	n, integer *nm, integer *id, doublereal *p3, doublereal *phz, doublereal *ph1, doublereal *p1,
	 doublereal *p2, doublereal *cp);

/* Subroutine */ int lfin_(integer *init, doublereal *theta, integer *l, integer *m,
	 integer *nm, doublereal *pb, integer *id, doublereal *wlfin);

/* Subroutine */ int lfin1_(integer *init, doublereal *theta, integer *l, integer *
	m, integer *nm, integer *id, doublereal *p3, doublereal *phz, doublereal *ph1, doublereal *p1,
	 doublereal *p2, doublereal *cp);

/* Subroutine */ int lfp_(integer *init, integer *n, integer *m, integer *l, 
	doublereal *cp, doublereal *pb, doublereal *w);

/* Subroutine */ int lfp1_(integer *init, integer *n, integer *m, integer *l, 
	doublereal *cp, doublereal *p, doublereal *wsave1, doublereal *wsave2, doublereal *wsave3);

/* Subroutine */ int lfpt_(integer *n, integer *m, doublereal *theta, doublereal *cp, 
	doublereal *pb);

/* Subroutine */ int sfvpec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, doublereal *vp, integer *idv, integer *jdv, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdb, integer *ndb, doublereal *
	wshsec, integer *lshsec, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int sfvpec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, doublereal *vp, integer *idv, integer *jdv, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdb, integer *ndb, doublereal *a, 
	doublereal *b, integer *mab, doublereal *fnn, doublereal *wshsec, integer *lshsec, doublereal 
	*wk, integer *lwk, integer *ierror);

/* Subroutine */ int sfvpes_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, doublereal *vp, integer *idv, integer *jdv, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdb, integer *ndb, doublereal *
	wshses, integer *lshses, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int sfvpes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, doublereal *vp, integer *idv, integer *jdv, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdb, integer *ndb, doublereal *a, 
	doublereal *b, integer *mab, doublereal *fnn, doublereal *wshses, integer *lshses, doublereal 
	*wk, integer *lwk, integer *ierror);

/* Subroutine */ int sfvpgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, doublereal *vp, integer *idv, integer *jdv, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdb, integer *ndb, doublereal *
	wshsgc, integer *lshsgc, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int sfvpgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, doublereal *vp, integer *idv, integer *jdv, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdb, integer *ndb, doublereal *a, 
	doublereal *b, integer *mab, doublereal *fnn, doublereal *wshsgc, integer *lshsgc, doublereal 
	*wk, integer *lwk, integer *ierror);

/* Subroutine */ int sfvpgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, doublereal *vp, integer *idv, integer *jdv, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdb, integer *ndb, doublereal *
	wshsgs, integer *lshsgs, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int stvpgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *sf, doublereal *vp, integer *idv, integer *jdv, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdb, integer *ndb, doublereal *a, 
	doublereal *b, integer *mab, doublereal *fnn, doublereal *wshsgs, integer *lshsgs, doublereal 
	*wk, integer *lwk, integer *ierror);

/* Subroutine */ int shaec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a, doublereal *b, 
	integer *mdab, integer *ndab, doublereal *wshaec, integer *lshaec, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int shaec1_(integer *nlat, integer *isym, integer *nt, doublereal *
	g, integer *idgs, integer *jdgs, doublereal *a, doublereal *b, integer *mdab, 
	integer *ndab, integer *imid, integer *idg, integer *jdg, doublereal *ge, 
	doublereal *go, doublereal *work, doublereal *zb, doublereal *wzfin, doublereal *whrfft);

/* Subroutine */ int shaeci_(integer *nlat, integer *nlon, doublereal *wshaec, 
	integer *lshaec, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int shaes_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a, doublereal *b, 
	integer *mdab, integer *ndab, doublereal *wshaes, integer *lshaes, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int shaes1_(integer *nlat, integer *isym, integer *nt, doublereal *
	g, integer *idgs, integer *jdgs, doublereal *a, doublereal *b, integer *mdab, 
	integer *ndab, doublereal *z__, integer *idz, integer *idg, integer *jdg, 
	doublereal *ge, doublereal *go, doublereal *work, doublereal *whrfft);

/* Subroutine */ int shaesi_(integer *nlat, integer *nlon, doublereal *wshaes, 
	integer *lshaes, doublereal *work, integer *lwork, doublereal *dwork, 
	integer *ldwork, integer *ierror);

/* Subroutine */ int shagc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a, doublereal *b, 
	integer *mdab, integer *ndab, doublereal *wshagc, integer *lshagc, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int shagc1_(integer *nlat, integer *nlon, integer *l, 
	integer *lat, integer *mode, doublereal *gs, integer *idg, integer *jdg, 
	integer *nt, doublereal *a, doublereal *b, integer *mdab, integer *ndab, doublereal *w, 
	doublereal *wts, doublereal *wfft, integer *late, doublereal *pmn, doublereal *g);

/* Subroutine */ int shagci_(integer *nlat, integer *nlon, doublereal *wshagc, 
	integer *lshagc, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int shagci1_(integer *nlat, integer *nlon, integer *l, 
	integer *late, doublereal *wts, doublereal *p0n, doublereal *p1n, doublereal *abel, doublereal *
	bbel, doublereal *cbel, doublereal *wfft, doublereal *dtheta, doublereal *dwts, 
	doublereal *work, integer *ier);

/* Subroutine */ int shags_(integer *nlat, integer *nlon, integer *mode, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a, doublereal *b, 
	integer *mdab, integer *ndab, doublereal *wshags, integer *lshags, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int shags1_(integer *nlat, integer *nlon, integer *l, 
	integer *lat, integer *mode, doublereal *gs, integer *idg, integer *jdg, 
	integer *nt, doublereal *a, doublereal *b, integer *mdab, integer *ndab, doublereal *
	wts, doublereal *wfft, doublereal *pmn, integer *late, doublereal *g, doublereal *work);

/* Subroutine */ int shagsi_(integer *nlat, integer *nlon, doublereal *wshags, 
	integer *lshags, doublereal *work, integer *lwork, doublereal *dwork, 
	integer *ldwork, integer *ierror);

/* Subroutine */ int shagss1_(integer *nlat, integer *l, integer *late, doublereal *
	w, doublereal *pmn, doublereal *pmnf);

/* Subroutine */ int shagsp_(integer *nlat, integer *nlon, doublereal *wshags, 
	integer *lshags, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int shagsp1_(integer *nlat, integer *nlon, integer *l, 
	integer *late, doublereal *wts, doublereal *p0n, doublereal *p1n, doublereal *abel, doublereal *
	bbel, doublereal *cbel, doublereal *wfft, doublereal *dtheta, doublereal *dwts, 
	doublereal *work, integer *ier);

/* Subroutine */ int shigc_(integer *nlat, integer *nlon, doublereal *wshigc, 
	integer *lshigc, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int shigc1_(integer *nlat, integer *nlon, integer *l, 
	integer *late, doublereal *wts, doublereal *p0n, doublereal *p1n, doublereal *abel, doublereal *
	bbel, doublereal *cbel, doublereal *wfft, doublereal *dtheta, doublereal *dwts, 
	doublereal *work, integer *ier);

/* Subroutine */ int shigs_(integer *nlat, integer *nlon, doublereal *wshigs, 
	integer *lshigs, doublereal *work, integer *lwork, doublereal *dwork, 
	integer *ldwork, integer *ierror);

/* Subroutine */ int shigss1_(integer *nlat, integer *l, integer *late, doublereal *
	w, doublereal *pmn, doublereal *pmnf);

/* Subroutine */ int shigsp_(integer *nlat, integer *nlon, doublereal *wshigs, 
	integer *lshigs, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int shigsp1_(integer *nlat, integer *nlon, integer *l, 
	integer *late, doublereal *wts, doublereal *p0n, doublereal *p1n, doublereal *abel, doublereal *
	bbel, doublereal *cbel, doublereal *wfft, doublereal *dtheta, doublereal *dwts, 
	doublereal *work, integer *ier);

/* Subroutine */ int shpei_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, doublereal *wshp, integer *lwshp, integer *iwshp, integer *
	liwshp, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int shpei1_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, integer *idp, integer *ierror, doublereal *pe, doublereal *po, 
	doublereal *ze, doublereal *zo, integer *ipse, integer *jzse, integer *ipso, 
	integer *jzso, doublereal *cp, doublereal *work, doublereal *wx, 
	doublereal *s, doublereal *e, doublereal *thet, doublereal *xx, 
	doublereal *z__, doublereal *a, doublereal *b, doublereal *we, 
	doublereal *ped, doublereal *wo, doublereal *pod, doublereal *u);

/* Subroutine */ int shpe_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, doublereal *x, doublereal *y, integer *idxy, doublereal *wshp, integer 
	*lwshp, integer *iwshp, integer *liwshp, doublereal *work, integer *lwork, 
	integer *ierror);

/* Subroutine */ int shpe1_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, doublereal *sx, doublereal *sy, integer *idxy, integer *ierror, 
	integer *idp, doublereal *pe, doublereal *po, doublereal *ze, doublereal *zo, integer *ipse, 
	integer *jzse, integer *ipso, integer *jzso, doublereal *xe, doublereal *xo, doublereal 
	*ye, doublereal *yo);

/* Subroutine */ int mxm_(integer *lr, integer *lc, integer *ld, doublereal *
	a, integer *mc, integer *md, doublereal *b, integer *nd, doublereal *
	c__);

/* Subroutine */ int smxm_(integer *lr, integer *lc, integer *ld, doublereal *a, 
	integer *mc, integer *md, doublereal *b, integer *nd, doublereal *c__);

/* Subroutine */ int mxmx_(integer *lr, integer *lc, integer *ld, doublereal *a, 
	integer *mc, integer *md, doublereal *b, doublereal *x, doublereal *y);

/* Subroutine */ int dmxmx_(integer *lr, integer *lc, integer *ld, doublereal 
	*a, integer *mc, integer *md, doublereal *b, doublereal *x, 
	doublereal *y);

/* Subroutine */ int trunc_(integer *irc, integer *n, integer *idp, 
	doublereal *a, integer *nrc, integer *ijs);

/* Subroutine */ int gs_(integer *n, doublereal *x, doublereal *y, doublereal 
	*z__);

/* Subroutine */ int normal_(integer *n, doublereal *x, integer *id, 
	doublereal *q);

/* Subroutine */ int coe_(integer *moe, integer *n, doublereal *x, doublereal 
	*dmax__);

/* Subroutine */ int dlfkp_(integer *m, integer *n, doublereal *cp);

/* Subroutine */ int dlftp_(integer *m, integer *n, doublereal *theta, 
	doublereal *cp, doublereal *pb);

/* Subroutine */ int dsvdc_(doublereal *x, integer *ldx, integer *n, integer *
	p, doublereal *s, doublereal *e, doublereal *u, integer *ldu, 
	doublereal *v, integer *ldv, doublereal *work, integer *job, integer *
	info);

/* Subroutine */ int daxpy_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx, doublereal *dy, integer *incy);

/* Subroutine */ int drot_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *c__, doublereal *s);

/* Subroutine */ int drotg_(doublereal *da, doublereal *db, doublereal *c__, 
	doublereal *s);

/* Subroutine */ int dscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx);

/* Subroutine */ int dswap_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy);

/* Subroutine */ int shpgi_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, doublereal *wshp, integer *lwshp, integer *iwshp, integer *
	liwshp, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int shpgi1_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, integer *idp, integer *ierror, doublereal *pe, doublereal *po, 
	doublereal *ze, doublereal *zo, integer *ipse, integer *jzse, integer *ipso, 
	integer *jzso, doublereal *cp, doublereal *wx, doublereal *thet, 
	doublereal *gwts, doublereal *xx, doublereal *z__, doublereal *a, 
	doublereal *b, doublereal *ped, doublereal *pod, doublereal *u);

/* Subroutine */ int shpg_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, doublereal *x, doublereal *y, integer *idxy, doublereal *wshp, integer 
	*lwshp, integer *iwshp, integer *liwshp, doublereal *work, integer *lwork, 
	integer *ierror);

/* Subroutine */ int shpg1_(integer *nlat, integer *nlon, integer *isym, 
	integer *mtrunc, doublereal *sx, doublereal *sy, integer *idxy, integer *ierror, 
	integer *idp, doublereal *pe, doublereal *po, doublereal *ze, doublereal *zo, integer *ipse, 
	integer *jzse, integer *ipso, integer *jzso, doublereal *xe, doublereal *xo, doublereal 
	*ye, doublereal *yo);

/* Subroutine */ int mxm_(integer *lr, integer *lc, integer *ld, doublereal *
	a, integer *mc, integer *md, doublereal *b, integer *nd, doublereal *
	c__);

/* Subroutine */ int smxm_(integer *lr, integer *lc, integer *ld, doublereal *a, 
	integer *mc, integer *md, doublereal *b, integer *nd, doublereal *c__);

/* Subroutine */ int mxmx_(integer *lr, integer *lc, integer *ld, doublereal *a, 
	integer *mc, integer *md, doublereal *b, doublereal *x, doublereal *y);

/* Subroutine */ int dmxmx_(integer *lr, integer *lc, integer *ld, doublereal 
	*a, integer *mc, integer *md, doublereal *b, doublereal *x, 
	doublereal *y);

/* Subroutine */ int trunc_(integer *irc, integer *n, integer *idp, 
	doublereal *a, integer *nrc, integer *ijs);

/* Subroutine */ int gs_(integer *n, doublereal *x, doublereal *y, doublereal 
	*z__);

/* Subroutine */ int normal_(integer *n, doublereal *x, integer *id, 
	doublereal *q);

/* Subroutine */ int coe_(integer *moe, integer *n, doublereal *x, doublereal 
	*dmax__);

/* Subroutine */ int dlfkg_(integer *m, integer *n, doublereal *cp);

/* Subroutine */ int dlftg_(integer *m, integer *n, doublereal *theta, 
	doublereal *cp, doublereal *pb);

/* Subroutine */ int gaqdp_(integer *nlat, doublereal *theta, doublereal *wts,
	 doublereal *w, integer *lwork, integer *ierror);

/* Subroutine */ int cpdp1_(integer *n, doublereal *cz, doublereal *cp, 
	doublereal *dcp);

/* Subroutine */ int tpdp1_(integer *n, doublereal *theta, doublereal *cz, 
	doublereal *cp, doublereal *dcp, doublereal *pb, doublereal *dpb);

/* Subroutine */ int shsec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a, doublereal *b, 
	integer *mdab, integer *ndab, doublereal *wshsec, integer *lshsec, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int shsec1_(integer *nlat, integer *isym, integer *nt, doublereal *
	g, integer *idgs, integer *jdgs, doublereal *a, doublereal *b, integer *mdab, 
	integer *ndab, integer *imid, integer *idg, integer *jdg, doublereal *ge, 
	doublereal *go, doublereal *work, doublereal *pb, doublereal *walin, doublereal *whrfft);

/* Subroutine */ int shseci_(integer *nlat, integer *nlon, doublereal *wshsec, 
	integer *lshsec, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int shses_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a, doublereal *b, 
	integer *mdab, integer *ndab, doublereal *wshses, integer *lshses, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int shses1_(integer *nlat, integer *isym, integer *nt, doublereal *
	g, integer *idgs, integer *jdgs, doublereal *a, doublereal *b, integer *mdab, 
	integer *ndab, doublereal *p, integer *imid, integer *idg, integer *jdg, 
	doublereal *ge, doublereal *go, doublereal *work, doublereal *whrfft);

/* Subroutine */ int shsesi_(integer *nlat, integer *nlon, doublereal *wshses, 
	integer *lshses, doublereal *work, integer *lwork, doublereal *dwork, 
	integer *ldwork, integer *ierror);

/* Subroutine */ int shsgc_(integer *nlat, integer *nlon, integer *mode, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a, doublereal *b, 
	integer *mdab, integer *ndab, doublereal *wshsgc, integer *lshsgc, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int shsgc1_(integer *nlat, integer *nlon, integer *l, 
	integer *lat, integer *mode, doublereal *gs, integer *idg, integer *jdg, 
	integer *nt, doublereal *a, doublereal *b, integer *mdab, integer *ndab, doublereal *w, 
	doublereal *wfft, integer *late, doublereal *pmn, doublereal *g);

/* Subroutine */ int shsgci_(integer *nlat, integer *nlon, doublereal *wshsgc, 
	integer *lshsgc, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int shsgci1_(integer *nlat, integer *nlon, integer *l, 
	integer *late, doublereal *wts, doublereal *p0n, doublereal *p1n, doublereal *abel, doublereal *
	bbel, doublereal *cbel, doublereal *wfft, doublereal *dtheta, doublereal *dwts, 
	doublereal *work, integer *ier);

/* Subroutine */ int shsgs_(integer *nlat, integer *nlon, integer *mode, 
	integer *nt, doublereal *g, integer *idg, integer *jdg, doublereal *a, doublereal *b, 
	integer *mdab, integer *ndab, doublereal *wshsgs, integer *lshsgs, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int shsgs1_(integer *nlat, integer *nlon, integer *l, 
	integer *lat, integer *mode, doublereal *gs, integer *idg, integer *jdg, 
	integer *nt, doublereal *a, doublereal *b, integer *mdab, integer *ndab, doublereal *
	wfft, doublereal *pmn, integer *late, doublereal *g, doublereal *work);

/* Subroutine */ int shsgsi_(integer *nlat, integer *nlon, doublereal *wshsgs, 
	integer *lshsgs, doublereal *work, integer *lwork, doublereal *dwork, 
	integer *ldwork, integer *ierror);

/* Subroutine */ int shsgss1_(integer *nlat, integer *l, integer *late, doublereal *
	w, doublereal *pmn, doublereal *pmnf);

/* Subroutine */ int shsgsp_(integer *nlat, integer *nlon, doublereal *wshsgs, 
	integer *lshsgs, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int shsgsp1_(integer *nlat, integer *nlon, integer *l, 
	integer *late, doublereal *wts, doublereal *p0n, doublereal *p1n, doublereal *abel, doublereal *
	bbel, doublereal *cbel, doublereal *wfft, doublereal *dtheta, doublereal *dwts, 
	doublereal *work, integer *ier);

/* Subroutine */ int slapec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *slap, integer *ids, integer *jds, doublereal *a, doublereal *b,
	 integer *mdab, integer *ndab, doublereal *wshsec, integer *lshsec, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int slapec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *slap, integer *ids, integer *jds, doublereal *a, doublereal *b,
	 integer *mdab, integer *ndab, doublereal *alap, doublereal *blap, integer *mmax, 
	doublereal *fnn, doublereal *wshsec, integer *lshsec, doublereal *wk, integer *lwk, 
	integer *ierror);

/* Subroutine */ int slapes_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *slap, integer *ids, integer *jds, doublereal *a, doublereal *b,
	 integer *mdab, integer *ndab, doublereal *wshses, integer *lshses, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int slapes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *slap, integer *ids, integer *jds, doublereal *a, doublereal *b,
	 integer *mdab, integer *ndab, doublereal *alap, doublereal *blap, integer *mmax, 
	doublereal *fnn, doublereal *wsave, integer *lsave, doublereal *wk, integer *lwk, 
	integer *ierror);

/* Subroutine */ int slapgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *slap, integer *ids, integer *jds, doublereal *a, doublereal *b,
	 integer *mdab, integer *ndab, doublereal *wshsgc, integer *lshsgc, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int slapgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *slap, integer *ids, integer *jds, doublereal *a, doublereal *b,
	 integer *mdab, integer *ndab, doublereal *alap, doublereal *blap, integer *mmax, 
	doublereal *fnn, doublereal *wsave, integer *lsave, doublereal *wk, integer *lwk, 
	integer *ierror);

/* Subroutine */ int slapgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *slap, integer *ids, integer *jds, doublereal *a, doublereal *b,
	 integer *mdab, integer *ndab, doublereal *wshsgs, integer *lshsgs, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int slapgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *slap, integer *ids, integer *jds, doublereal *a, doublereal *b,
	 integer *mdab, integer *ndab, doublereal *alap, doublereal *blap, integer *mmax, 
	doublereal *fnn, doublereal *wsave, integer *lsave, doublereal *wk, integer *lwk, 
	integer *ierror);

/* Subroutine */ int dnlfk_(integer *m, integer *n, doublereal *cp);

/* Subroutine */ int dnlft_(integer *m, integer *n, doublereal *theta, 
	doublereal *cp, doublereal *pb);

/* Subroutine */ int dnlftd_(integer *m, integer *n, doublereal *theta, 
	doublereal *cp, doublereal *pb);

/* Subroutine */ int legin_(integer *mode, integer *l, integer *nlat, integer 
	*m, doublereal *w, doublereal *pmn, integer *km);

/* Subroutine */ int legin1_(integer *mode, integer *l, integer *nlat, 
	integer *late, integer *m, doublereal *p0n, doublereal *p1n, doublereal *abel, doublereal *
	bbel, doublereal *cbel, doublereal *pmn, integer *km);

/* Subroutine */ int zfin_(integer *isym, integer *nlat, integer *nlon, 
	integer *m, doublereal *z__, integer *i3, doublereal *wzfin);

/* Subroutine */ int zfin1_(integer *isym, integer *nlat, integer *m, doublereal *
	z__, integer *imid, integer *i3, doublereal *zz, doublereal *z1, doublereal *a, doublereal *b,
	 doublereal *c__);

/* Subroutine */ int zfinit_(integer *nlat, integer *nlon, doublereal *wzfin, 
	doublereal *dwork);

/* Subroutine */ int zfini1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *z__, doublereal *abc, doublereal *cz, doublereal *work);

/* Subroutine */ int dnzfk_(integer *nlat, integer *m, integer *n, doublereal 
	*cz, doublereal *work);

/* Subroutine */ int dnzft_(integer *nlat, integer *m, integer *n, doublereal 
	*th, doublereal *cz, doublereal *zh);

/* Subroutine */ int alin_(integer *isym, integer *nlat, integer *nlon, 
	integer *m, doublereal *p, integer *i3, doublereal *walin);

/* Subroutine */ int alin1_(integer *isym, integer *nlat, integer *m, doublereal *p,
	 integer *imid, integer *i3, doublereal *pz, doublereal *p1, doublereal *a, doublereal *b, 
	doublereal *c__);

/* Subroutine */ int alinit_(integer *nlat, integer *nlon, doublereal *walin, 
	doublereal *dwork);

/* Subroutine */ int alini1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *p, doublereal *abc, doublereal *cp);

/* Subroutine */ int rabcp_(integer *nlat, integer *nlon, doublereal *abc);

/* Subroutine */ int rabcp1_(integer *nlat, integer *nlon, doublereal *a, doublereal *b, 
	doublereal *c__);

/* Subroutine */ int sea1_(integer *nlat, integer *nlon, integer *imid, doublereal *
	z__, integer *idz, doublereal *zin, doublereal *wzfin, doublereal *dwork);

/* Subroutine */ int ses1_(integer *nlat, integer *nlon, integer *imid, doublereal *
	p, doublereal *pin, doublereal *walin, doublereal *dwork);

/* Subroutine */ int zvinit_(integer *nlat, integer *nlon, doublereal *wzvin, 
	doublereal *dwork);

/* Subroutine */ int zvini1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *zv, doublereal *abc, doublereal *czv, doublereal *work);

/* Subroutine */ int zwinit_(integer *nlat, integer *nlon, doublereal *wzwin, 
	doublereal *dwork);

/* Subroutine */ int zwini1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *zw, doublereal *abc, doublereal *czw, doublereal *work);

/* Subroutine */ int zvin_(integer *ityp, integer *nlat, integer *nlon, 
	integer *m, doublereal *zv, integer *i3, doublereal *wzvin);

/* Subroutine */ int zvin1_(integer *ityp, integer *nlat, integer *m, doublereal *
	zv, integer *imid, integer *i3, doublereal *zvz, doublereal *zv1, doublereal *a, doublereal *
	b, doublereal *c__);

/* Subroutine */ int zwin_(integer *ityp, integer *nlat, integer *nlon, 
	integer *m, doublereal *zw, integer *i3, doublereal *wzwin);

/* Subroutine */ int zwin1_(integer *ityp, integer *nlat, integer *m, doublereal *
	zw, integer *imid, integer *i3, doublereal *zw1, doublereal *zw2, doublereal *a, doublereal *
	b, doublereal *c__);

/* Subroutine */ int vbinit_(integer *nlat, integer *nlon, doublereal *wvbin, 
	doublereal *dwork);

/* Subroutine */ int vbini1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *vb, doublereal *abc, doublereal *cvb, doublereal *work);

/* Subroutine */ int wbinit_(integer *nlat, integer *nlon, doublereal *wwbin, 
	doublereal *dwork);

/* Subroutine */ int wbini1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *wb, doublereal *abc, doublereal *cwb, doublereal *work);

/* Subroutine */ int vbin_(integer *ityp, integer *nlat, integer *nlon, 
	integer *m, doublereal *vb, integer *i3, doublereal *wvbin);

/* Subroutine */ int vbin1_(integer *ityp, integer *nlat, integer *m, doublereal *
	vb, integer *imid, integer *i3, doublereal *vbz, doublereal *vb1, doublereal *a, doublereal *
	b, doublereal *c__);

/* Subroutine */ int wbin_(integer *ityp, integer *nlat, integer *nlon, 
	integer *m, doublereal *wb, integer *i3, doublereal *wwbin);

/* Subroutine */ int wbin1_(integer *ityp, integer *nlat, integer *m, doublereal *
	wb, integer *imid, integer *i3, doublereal *wb1, doublereal *wb2, doublereal *a, doublereal *
	b, doublereal *c__);

/* Subroutine */ int dzvk_(integer *nlat, integer *m, integer *n, doublereal *
	czv, doublereal *work);

/* Subroutine */ int dzvt_(integer *nlat, integer *m, integer *n, doublereal *
	th, doublereal *czv, doublereal *zvh);

/* Subroutine */ int dzwk_(integer *nlat, integer *m, integer *n, doublereal *
	czw, doublereal *work);

/* Subroutine */ int dzwt_(integer *nlat, integer *m, integer *n, doublereal *
	th, doublereal *czw, doublereal *zwh);

/* Subroutine */ int dvbk_(integer *m, integer *n, doublereal *cv, doublereal 
	*work);

/* Subroutine */ int dwbk_(integer *m, integer *n, doublereal *cw, doublereal 
	*work);

/* Subroutine */ int dvbt_(integer *m, integer *n, doublereal *theta, 
	doublereal *cv, doublereal *vh);

/* Subroutine */ int dwbt_(integer *m, integer *n, doublereal *theta, 
	doublereal *cw, doublereal *wh);

/* Subroutine */ int rabcv_(integer *nlat, integer *nlon, doublereal *abc);

/* Subroutine */ int rabcv1_(integer *nlat, integer *nlon, doublereal *a, doublereal *b, 
	doublereal *c__);

/* Subroutine */ int rabcw_(integer *nlat, integer *nlon, doublereal *abc);

/* Subroutine */ int rabcw1_(integer *nlat, integer *nlon, doublereal *a, doublereal *b, 
	doublereal *c__);

/* Subroutine */ int vtinit_(integer *nlat, integer *nlon, doublereal *wvbin, 
	doublereal *dwork);

/* Subroutine */ int vtini1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *vb, doublereal *abc, doublereal *cvb, doublereal *work);

/* Subroutine */ int wtinit_(integer *nlat, integer *nlon, doublereal *wwbin, 
	doublereal *dwork);

/* Subroutine */ int wtini1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *wb, doublereal *abc, doublereal *cwb, doublereal *work);

/* Subroutine */ int vtgint_(integer *nlat, integer *nlon, doublereal *theta, 
	doublereal *wvbin, doublereal *work);

/* Subroutine */ int vtgit1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *theta, doublereal *vb, doublereal *abc, doublereal *cvb, doublereal *
	work);

/* Subroutine */ int wtgint_(integer *nlat, integer *nlon, doublereal *theta, 
	doublereal *wwbin, doublereal *work);

/* Subroutine */ int wtgit1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *theta, doublereal *wb, doublereal *abc, doublereal *cwb, doublereal *
	work);

/* Subroutine */ int dvtk_(integer *m, integer *n, doublereal *cv, doublereal 
	*work);

/* Subroutine */ int dwtk_(integer *m, integer *n, doublereal *cw, doublereal 
	*work);

/* Subroutine */ int dvtt_(integer *m, integer *n, doublereal *theta, 
	doublereal *cv, doublereal *vh);

/* Subroutine */ int dwtt_(integer *m, integer *n, doublereal *theta, 
	doublereal *cw, doublereal *wh);

/* Subroutine */ int vbgint_(integer *nlat, integer *nlon, doublereal *theta, 
	doublereal *wvbin, doublereal *work);

/* Subroutine */ int vbgit1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *theta, doublereal *vb, doublereal *abc, doublereal *cvb, doublereal *
	work);

/* Subroutine */ int wbgint_(integer *nlat, integer *nlon, doublereal *theta, 
	doublereal *wwbin, doublereal *work);

/* Subroutine */ int wbgit1_(integer *nlat, integer *nlon, integer *imid, 
	doublereal *theta, doublereal *wb, doublereal *abc, doublereal *cwb, doublereal *
	work);

/* Subroutine */ int sshifte_(integer *ioff, integer *nlon, integer *nlat, 
	doublereal *goff, doublereal *greg, doublereal *wsav, integer *lsav, doublereal *wrk, integer 
	*lwrk, integer *ier);

/* Subroutine */ int shftoff_(integer *nlon, integer *nlat, doublereal *goff, doublereal *
	greg, doublereal *wsav, integer *nr, integer *nlat2, doublereal *rlat, doublereal *rlon,
	 doublereal *wrk);

/* Subroutine */ int shftreg_(integer *nlon, integer *nlat, doublereal *goff, doublereal *
	greg, doublereal *wsav, integer *nr, integer *nlat2, integer *nlatp1, doublereal *
	rlat, doublereal *rlon, doublereal *wrk);

/* Subroutine */ int sshifti_(integer *ioff, integer *nlon, integer *nlat, 
	integer *lsav, doublereal *wsav, integer *ier);

/* Subroutine */ int shifth_(integer *m, integer *n, doublereal *r__, doublereal *wsav, 
	doublereal *work);

/* Subroutine */ int shifthi_(integer *n, doublereal *dp, doublereal *wsav);

/* Subroutine */ int trssph_(integer *intl, integer *igrida, integer *nlona, 
	integer *nlata, doublereal *da, integer *igridb, integer *nlonb, integer *
	nlatb, doublereal *db, doublereal *wsave, integer *lsave, integer *lsvmin, doublereal *
	work, integer *lwork, integer *lwkmin, doublereal *dwork, integer *
	ldwork, integer *ier);

/* Subroutine */ int trab_(integer *ma, integer *na, doublereal *aa, doublereal *ba, 
	integer *mb, integer *nb, doublereal *ab, doublereal *bb);

/* Subroutine */ int trsplat_(integer *n, integer *m, doublereal *data, doublereal *work);

/* Subroutine */ int convlat_(integer *nlat, integer *nlon, doublereal *data);

/* Subroutine */ int trvsph_(integer *intl, integer *igrida, integer *nlona, 
	integer *nlata, integer *iveca, doublereal *ua, doublereal *va, integer *igridb, 
	integer *nlonb, integer *nlatb, integer *ivecb, doublereal *ub, doublereal *vb, 
	doublereal *wsave, integer *lsave, integer *lsvmin, doublereal *work, integer *
	lwork, integer *lwkmin, doublereal *dwork, integer *ldwork, integer *
	ier);

/* Subroutine */ int negv_(integer *nlat, integer *nlon, doublereal *v);

/* Subroutine */ int trvab_(integer *ma, integer *na, doublereal *abr, doublereal *abi, 
	doublereal *acr, doublereal *aci, integer *mb, integer *nb, doublereal *bbr, doublereal *bbi, 
	doublereal *bcr, doublereal *bci);

/* Subroutine */ int trvplat_(integer *n, integer *m, doublereal *data, doublereal *work);

/* Subroutine */ int covlat_(integer *nlat, integer *nlon, doublereal *data);

/* Subroutine */ int vhaec_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvhaec, integer *lvhaec, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vhaec1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *v, 
	doublereal *w, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *ve, doublereal *vo, doublereal *we, doublereal *wo, doublereal *
	zv, doublereal *zw, doublereal *wzvin, doublereal *wzwin, doublereal *wrfft);

/* Subroutine */ int vhaeci_(integer *nlat, integer *nlon, doublereal *wvhaec, 
	integer *lvhaec, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int vhaes_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvhaes, integer *lvhaes, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vhaes1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *v, 
	doublereal *w, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *ve, doublereal *vo, doublereal *we, doublereal *wo, doublereal *
	work, integer *idz, doublereal *zv, doublereal *zw, doublereal *wrfft);

/* Subroutine */ int vhaesi_(integer *nlat, integer *nlon, doublereal *wvhaes, 
	integer *lvhaes, doublereal *work, integer *lwork, doublereal *dwork, 
	integer *ldwork, integer *ierror);

/* Subroutine */ int vea1_(integer *nlat, integer *nlon, integer *imid, doublereal *
	zv, doublereal *zw, integer *idz, doublereal *zin, doublereal *wzvin, doublereal *dwork);

/* Subroutine */ int vhagc_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvhagc, integer *lvhagc, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vhagc1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *v, 
	doublereal *w, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *ve, doublereal *vo, doublereal *we, doublereal *wo, doublereal *
	vb, doublereal *wb, doublereal *wts, doublereal *wvbin, doublereal *wwbin, doublereal *wrfft);

/* Subroutine */ int vhagci_(integer *nlat, integer *nlon, doublereal *wvhagc, 
	integer *lvhagc, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int setwts_(integer *imid, doublereal *dwts, doublereal *wts);

/* Subroutine */ int vhags_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvhags, integer *lvhags, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vhags1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *v, 
	doublereal *w, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *ve, doublereal *vo, doublereal *we, doublereal *wo, doublereal *
	work, integer *idz, doublereal *vb, doublereal *wb, doublereal *wrfft);

/* Subroutine */ int vhagsi_(integer *nlat, integer *nlon, doublereal *wvhags, 
	integer *lvhags, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int vhgai1_(integer *nlat, integer *imid, doublereal *vb, doublereal *wb,
	 doublereal *dthet, doublereal *dwts, doublereal *dpbar, doublereal *
	work);

/* Subroutine */ int vhsec_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvhsec, integer *lvhsec, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vhsec1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *v, 
	doublereal *w, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *ve, doublereal *vo, doublereal *we, doublereal *wo, doublereal *
	vb, doublereal *wb, doublereal *wvbin, doublereal *wwbin, doublereal *wrfft);

/* Subroutine */ int vhseci_(integer *nlat, integer *nlon, doublereal *wvhsec, 
	integer *lvhsec, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int vhses_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvhses, integer *lvhses, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vhses1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *v, 
	doublereal *w, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *ve, doublereal *vo, doublereal *we, doublereal *wo, doublereal *
	work, integer *idz, doublereal *vb, doublereal *wb, doublereal *wrfft);

/* Subroutine */ int vhsesi_(integer *nlat, integer *nlon, doublereal *wvhses, 
	integer *lvhses, doublereal *work, integer *lwork, doublereal *dwork, 
	integer *ldwork, integer *ierror);

/* Subroutine */ int ves1_(integer *nlat, integer *nlon, integer *imid, doublereal *
	vb, doublereal *wb, integer *idz, doublereal *vin, doublereal *wzvin, doublereal *dwork);

/* Subroutine */ int vhsgc_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvhsgc, integer *lvhsgc, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vhsgc1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *v, 
	doublereal *w, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *ve, doublereal *vo, doublereal *we, doublereal *wo, doublereal *
	vb, doublereal *wb, doublereal *wvbin, doublereal *wwbin, doublereal *wrfft);

/* Subroutine */ int vhsgci_(integer *nlat, integer *nlon, doublereal *wvhsgc, 
	integer *lvhsgc, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int vhsgs_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *v, doublereal *w, integer *idvw, integer *jdvw, doublereal *br,
	 doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvhsgs, integer *lvhsgs, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vhsgs1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *v, 
	doublereal *w, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *ve, doublereal *vo, doublereal *we, doublereal *wo, doublereal *
	work, integer *idz, doublereal *vb, doublereal *wb, doublereal *wrfft);

/* Subroutine */ int vhsgsi_(integer *nlat, integer *nlon, doublereal *wvhsgs, 
	integer *lvhsgs, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int vhgsi1_(integer *nlat, integer *imid, doublereal *vb, doublereal *wb,
	 doublereal *dthet, doublereal *dwts, doublereal *dpbar, doublereal *
	work);

/* Subroutine */ int visequ_(integer *nlat, integer *nlon, doublereal *h__, integer 
	*len, doublereal *eyer, doublereal *eyelat, doublereal *eyelon, doublereal *wk, integer *lwk, 
	integer *iwk, integer *liwk, integer *ierror);

/* Subroutine */ int interp_(doublereal *h__, integer *len, integer *m, integer *n, 
	doublereal *w1, doublereal *w2, integer *iflag);

/* Subroutine */ int sptc_(doublereal *r__, integer *len, integer *m, integer *n, 
	doublereal *clat, doublereal *slat, doublereal *x, doublereal *y, doublereal *z__);

/* Subroutine */ int diag_(integer *m, integer *n, doublereal *xp, doublereal *yp, 
	integer *iflag);

/* Subroutine */ int stride_(integer *m, integer *n, integer *mst, integer *
	mfac);

/* Subroutine */ int triang_(integer *m, integer *n, doublereal *x, doublereal *y, doublereal *
	z__, integer *itri, doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, 
	doublereal *z2, doublereal *x3, doublereal *y3, doublereal *z3, integer *ityp, integer *iflag,
	 integer *mst);

/* Subroutine */ int vsurf_(doublereal *xeye, doublereal *yeye, doublereal *zeye, integer *ntri,
	 doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, doublereal *z2, doublereal *x3,
	 doublereal *y3, doublereal *z3, integer *itype, doublereal *work, integer *iwork);

/* Subroutine */ int vsurf1_(doublereal *xeye, doublereal *yeye, doublereal *zeye, integer *
	ntri, doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, doublereal *z2, 
	doublereal *x3, doublereal *y3, doublereal *z3, integer *itype, doublereal *px1, doublereal *py1, 
	doublereal *px2, doublereal *py2, doublereal *px3, doublereal *py3, doublereal *vx1, doublereal *vy1, 
	doublereal *vx2, doublereal *vy2, doublereal *vx3, doublereal *vy3, doublereal *tl, doublereal *tr, 
	integer *kh, integer *next, integer *istart, integer *ifinal);

/* Subroutine */ int prjct_(integer *init, doublereal *xeye, doublereal *yeye, doublereal *zeye,
	 doublereal *x, doublereal *y, doublereal *z__, doublereal *px, doublereal *py);

/* Subroutine */ int box_(integer *isd, integer *istart, integer *next, 
	integer *l, integer *list);

/* Subroutine */ int projct_(integer *m, integer *n, doublereal *xeye, doublereal *yeye, 
	doublereal *zeye, doublereal *x, doublereal *y, doublereal *z__, doublereal *px, doublereal *py);

/* Subroutine */ int visgau_(integer *nlat, integer *nlon, doublereal *h__, integer 
	*len, doublereal *eyer, doublereal *eyelat, doublereal *eyelon, doublereal *theta, doublereal 
	*wk, integer *lwk, integer *iwk, integer *liwk, integer *ierror);

/* Subroutine */ int embed_(integer *nlat, integer *nlon, doublereal *h__, integer *
	len, doublereal *hg);

/* Subroutine */ int intrpg_(doublereal *h__, integer *m, integer *n, doublereal *w1, 
	doublereal *w2, integer *iflag);

/* Subroutine */ int sptcg_(doublereal *r__, integer *m, integer *n, doublereal *
	theta, doublereal *clat, doublereal *slat, doublereal *x, doublereal *y, doublereal *z__);

/* Subroutine */ int diag_(integer *m, integer *n, doublereal *xp, doublereal *yp, 
	integer *iflag);

/* Subroutine */ int stride_(integer *m, integer *n, integer *mst, integer *
	mfac);

/* Subroutine */ int trigau_(integer *m, integer *n, doublereal *x, doublereal *y, doublereal *
	z__, integer *itri, doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, 
	doublereal *z2, doublereal *x3, doublereal *y3, doublereal *z3, integer *ityp, integer *iflag,
	 integer *mst);

/* Subroutine */ int vsurf_(doublereal *xeye, doublereal *yeye, doublereal *zeye, integer *ntri,
	 doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, doublereal *z2, doublereal *x3,
	 doublereal *y3, doublereal *z3, integer *itype, doublereal *work, integer *iwork);

/* Subroutine */ int vsurf1_(doublereal *xeye, doublereal *yeye, doublereal *zeye, integer *
	ntri, doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, doublereal *z2, 
	doublereal *x3, doublereal *y3, doublereal *z3, integer *itype, doublereal *px1, doublereal *py1, 
	doublereal *px2, doublereal *py2, doublereal *px3, doublereal *py3, doublereal *vx1, doublereal *vy1, 
	doublereal *vx2, doublereal *vy2, doublereal *vx3, doublereal *vy3, doublereal *tl, doublereal *tr, 
	integer *kh, integer *next, integer *istart, integer *ifinal);

/* Subroutine */ int prjct_(integer *init, doublereal *xeye, doublereal *yeye, doublereal *zeye,
	 doublereal *x, doublereal *y, doublereal *z__, doublereal *px, doublereal *py);

/* Subroutine */ int box_(integer *isd, integer *istart, integer *next, 
	integer *l, integer *list);

/* Subroutine */ int projct_(integer *m, integer *n, doublereal *xeye, doublereal *yeye, 
	doublereal *zeye, doublereal *x, doublereal *y, doublereal *z__, doublereal *px, doublereal *py);

/* Subroutine */ int visgeo_(integer *m, integer *idp, integer *jdp, doublereal *x, 
	doublereal *y, doublereal *z__, doublereal *h__, doublereal *eyer, doublereal *eyelat, doublereal *eyelon,
	 doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer 
	*ierror);

/* Subroutine */ int visgeo1_(integer *m, integer *idp, integer *jdp, doublereal *
	h__, doublereal *eyer, doublereal *eyelat, doublereal *eyelon, doublereal *xi, doublereal *yi, doublereal 
	*zi, doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, doublereal *z2, doublereal 
	*x3, doublereal *y3, doublereal *z3, integer *itype, doublereal *work, doublereal *x, doublereal *y,
	 doublereal *z__, integer *iwork);

/* Subroutine */ int ctos_(doublereal *x, doublereal *y, doublereal *z__, doublereal *r__, doublereal *
	theta, doublereal *phi);

/* Subroutine */ int stoc_(doublereal *r__, doublereal *theta, doublereal *phi, doublereal *x, doublereal *
	y, doublereal *z__);

/* Subroutine */ int vsurf_(doublereal *xeye, doublereal *yeye, doublereal *zeye, integer *ntri,
	 doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, doublereal *z2, doublereal *x3,
	 doublereal *y3, doublereal *z3, integer *itype, doublereal *work, integer *iwork);

/* Subroutine */ int vsurf1_(doublereal *xeye, doublereal *yeye, doublereal *zeye, integer *
	ntri, doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, doublereal *z2, 
	doublereal *x3, doublereal *y3, doublereal *z3, integer *itype, doublereal *px1, doublereal *py1, 
	doublereal *px2, doublereal *py2, doublereal *px3, doublereal *py3, doublereal *vx1, doublereal *vy1, 
	doublereal *vx2, doublereal *vy2, doublereal *vx3, doublereal *vy3, doublereal *tl, doublereal *tr, 
	integer *kh, integer *next, integer *istart, integer *ifinal);

/* Subroutine */ int prjct_(integer *init, doublereal *xeye, doublereal *yeye, doublereal *zeye,
	 doublereal *x, doublereal *y, doublereal *z__, doublereal *px, doublereal *py);

/* Subroutine */ int box_(integer *isd, integer *istart, integer *next, 
	integer *l, integer *list);

/* Subroutine */ int vlapec_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vlap, doublereal *wlap, integer *idvw, integer *jdvw, 
	doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, integer *mdbc, integer *ndbc, 
	doublereal *wvhsec, integer *lvhsec, doublereal *work, integer *lwork, integer *
	ierror);

/* Subroutine */ int vlapec1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vlap, doublereal *wlap, integer *idvw, integer *jdvw, 
	doublereal *brlap, doublereal *bilap, doublereal *crlap, doublereal *cilap, integer *mmax, 
	doublereal *fnn, integer *mdb, integer *ndb, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, doublereal *wsave, integer *lwsav, doublereal *wk, integer *lwk, 
	integer *ierror);

/* Subroutine */ int vlapes_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vlap, doublereal *wlap, integer *idvw, integer *jdvw, 
	doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, integer *mdbc, integer *ndbc, 
	doublereal *wvhses, integer *lvhses, doublereal *work, integer *lwork, integer *
	ierror);

/* Subroutine */ int vlapes1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vlap, doublereal *wlap, integer *idvw, integer *jdvw, 
	doublereal *brlap, doublereal *bilap, doublereal *crlap, doublereal *cilap, integer *mmax, 
	doublereal *fnn, integer *mdb, integer *ndb, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, doublereal *wsave, integer *lsave, doublereal *wk, integer *lwk, 
	integer *ierror);

/* Subroutine */ int vlapgc_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vlap, doublereal *wlap, integer *idvw, integer *jdvw, 
	doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, integer *mdbc, integer *ndbc, 
	doublereal *wvhsgc, integer *lvhsgc, doublereal *work, integer *lwork, integer *
	ierror);

/* Subroutine */ int vlapgc1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vlap, doublereal *wlap, integer *idvw, integer *jdvw, 
	doublereal *brlap, doublereal *bilap, doublereal *crlap, doublereal *cilap, integer *mmax, 
	doublereal *fnn, integer *mdb, integer *ndb, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, doublereal *wsave, integer *lwsav, doublereal *wk, integer *lwk, 
	integer *ierror);

/* Subroutine */ int vlapgs_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vlap, doublereal *wlap, integer *idvw, integer *jdvw, 
	doublereal *br, doublereal *bi, doublereal *cr, doublereal *ci, integer *mdbc, integer *ndbc, 
	doublereal *wvhsgs, integer *lvhsgs, doublereal *work, integer *lwork, integer *
	ierror);

/* Subroutine */ int vlapgs1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vlap, doublereal *wlap, integer *idvw, integer *jdvw, 
	doublereal *brlap, doublereal *bilap, doublereal *crlap, doublereal *cilap, integer *mmax, 
	doublereal *fnn, integer *mdb, integer *ndb, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, doublereal *wsave, integer *lsave, doublereal *wk, integer *lwk, 
	integer *ierror);

/* Subroutine */ int vrtec_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *vort, integer *ivrt, integer *jvrt, doublereal *cr, doublereal 
	*ci, integer *mdc, integer *ndc, doublereal *wshsec, integer *lshsec, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int vrtec1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *vort, integer *ivrt, integer *jvrt, doublereal *cr, doublereal 
	*ci, integer *mdc, integer *ndc, doublereal *a, doublereal *b, integer *mab, doublereal 
	*sqnn, doublereal *wshsec, integer *lshsec, doublereal *wk, integer *lwk, integer 
	*ierror);

/* Subroutine */ int vrtes_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *vort, integer *ivrt, integer *jvrt, doublereal *cr, doublereal 
	*ci, integer *mdc, integer *ndc, doublereal *wshses, integer *lshses, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int vrtes1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *vort, integer *ivrt, integer *jvrt, doublereal *cr, doublereal 
	*ci, integer *mdc, integer *ndc, doublereal *a, doublereal *b, integer *mab, doublereal 
	*sqnn, doublereal *wsav, integer *lwsav, doublereal *wk, integer *lwk, integer *
	ierror);

/* Subroutine */ int vrtgc_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *vort, integer *ivrt, integer *jvrt, doublereal *cr, doublereal 
	*ci, integer *mdc, integer *ndc, doublereal *wshsgc, integer *lshsgc, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int vrtgc1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *vort, integer *ivrt, integer *jvrt, doublereal *cr, doublereal 
	*ci, integer *mdc, integer *ndc, doublereal *a, doublereal *b, integer *mab, doublereal 
	*sqnn, doublereal *wsav, integer *lwsav, doublereal *wk, integer *lwk, integer *
	ierror);

/* Subroutine */ int vrtgs_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *vort, integer *ivrt, integer *jvrt, doublereal *cr, doublereal 
	*ci, integer *mdc, integer *ndc, doublereal *wshsgs, integer *lshsgs, doublereal *
	work, integer *lwork, integer *ierror);

/* Subroutine */ int vrtgs1_(integer *nlat, integer *nlon, integer *isym, 
	integer *nt, doublereal *vort, integer *ivrt, integer *jvrt, doublereal *cr, doublereal 
	*ci, integer *mdc, integer *ndc, doublereal *a, doublereal *b, integer *mab, doublereal 
	*sqnn, doublereal *wsav, integer *lwsav, doublereal *wk, integer *lwk, integer *
	ierror);

/* Subroutine */ int vshifte_(integer *ioff, integer *nlon, integer *nlat, 
	doublereal *uoff, doublereal *voff, doublereal *ureg, doublereal *vreg, doublereal *wsav, integer *
	lsav, doublereal *wrk, integer *lwrk, integer *ier);

/* Subroutine */ int vhftoff_(integer *nlon, integer *nlat, doublereal *uoff, doublereal *
	ureg, doublereal *wsav, integer *nr, integer *nlat2, integer *nlatp1, doublereal *
	rlatu, doublereal *rlonu, doublereal *rlou, doublereal *wrk);

/* Subroutine */ int vhftreg_(integer *nlon, integer *nlat, doublereal *uoff, doublereal *
	ureg, doublereal *wsav, integer *nr, integer *nlat2, integer *nlatp1, doublereal *
	rlatu, doublereal *rlonu, doublereal *rlou, doublereal *wrk);

/* Subroutine */ int vshifti_(integer *ioff, integer *nlon, integer *nlat, 
	integer *lsav, doublereal *wsav, integer *ier);

/* Subroutine */ int vhifth_(integer *m, integer *n, doublereal *r__, doublereal *wsav, 
	doublereal *work);

/* Subroutine */ int vhifthi_(integer *n, doublereal *dp, doublereal *wsav);

/* Subroutine */ int vsurf_(doublereal *xeye, doublereal *yeye, doublereal *zeye, integer *ntri,
	 doublereal *x1, doublereal *y1, doublereal *z1, doublereal *x2, doublereal *y2, doublereal *z2, doublereal *x3,
	 doublereal *y3, doublereal *z3, integer *itype, doublereal *work, integer *iwork);

/* Subroutine */ int vtsec_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vt, doublereal *wt, integer *idvw, integer *jdvw, doublereal *
	br, doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvts, integer *lwvts, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vtsec1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *vt, 
	doublereal *wt, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *vte, doublereal *vto, doublereal *wte, doublereal *wto, 
	doublereal *vb, doublereal *wb, doublereal *wvbin, doublereal *wwbin, doublereal *wrfft);

/* Subroutine */ int vtseci_(integer *nlat, integer *nlon, doublereal *wvts, 
	integer *lwvts, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int vtses_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vt, doublereal *wt, integer *idvw, integer *jdvw, doublereal *
	br, doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvts, integer *lwvts, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vtses1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *vt, 
	doublereal *wt, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *vte, doublereal *vto, doublereal *wte, doublereal *wto, 
	doublereal *work, integer *idz, doublereal *vb, doublereal *wb, doublereal *wrfft);

/* Subroutine */ int vtsesi_(integer *nlat, integer *nlon, doublereal *wvts, 
	integer *lwvts, doublereal *work, integer *lwork, doublereal *dwork, 
	integer *ldwork, integer *ierror);

/* Subroutine */ int vet1_(integer *nlat, integer *nlon, integer *imid, doublereal *
	vb, doublereal *wb, integer *idz, doublereal *vin, doublereal *wzvin, doublereal *dwork);

/* Subroutine */ int vtsgc_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vt, doublereal *wt, integer *idvw, integer *jdvw, doublereal *
	br, doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvts, integer *lwvts, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vtsgc1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *vt, 
	doublereal *wt, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *vte, doublereal *vto, doublereal *wte, doublereal *wto, 
	doublereal *vb, doublereal *wb, doublereal *wvbin, doublereal *wwbin, doublereal *wrfft);

/* Subroutine */ int vtsgci_(integer *nlat, integer *nlon, doublereal *wvts, 
	integer *lwvts, doublereal *dwork, integer *ldwork, integer *ierror);

/* Subroutine */ int vtsgs_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, doublereal *vt, doublereal *wt, integer *idvw, integer *jdvw, doublereal *
	br, doublereal *bi, doublereal *cr, doublereal *ci, integer *mdab, integer *ndab, doublereal *
	wvts, integer *lwvts, doublereal *work, integer *lwork, integer *ierror);

/* Subroutine */ int vtsgs1_(integer *nlat, integer *nlon, integer *ityp, 
	integer *nt, integer *imid, integer *idvw, integer *jdvw, doublereal *vt, 
	doublereal *wt, integer *mdab, integer *ndab, doublereal *br, doublereal *bi, doublereal *cr, 
	doublereal *ci, integer *idv, doublereal *vte, doublereal *vto, doublereal *wte, doublereal *wto, 
	doublereal *work, integer *idz, doublereal *vb, doublereal *wb, doublereal *wrfft);

/* Subroutine */ int vtsgsi_(integer *nlat, integer *nlon, doublereal *wvts, 
	integer *lwvts, doublereal *work, integer *lwork, doublereal *dwork, 
	integer *ldwork, integer *ierror);

/* Subroutine */ int vetg1_(integer *nlat, integer *nlon, integer *imid, doublereal 
	*vb, doublereal *wb, doublereal *vin, doublereal *wvbin, doublereal *theta, doublereal *
	wts, doublereal *dwork, integer *ierror);


#undef abs
#undef dabs
#undef min
#undef max
#undef dmin
#undef dmax
#undef bit_test
#undef bit_clear
#undef bit_set

#ifdef __cplusplus
}
#endif

#endif /*SPHEREPACK_H*/
