/* pushp.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

Extern struct {
    doublereal aw0, xkx, xky, wcoefz[3], xlamd, fbess0, delaw, awd, awx, awy, 
	    gamma0, delgam, rxbeam, rybeam, alphax, alphay, emitx, emity, 
	    xbeam, ybeam, pxbeam, pybeam, cuttail, curpeak, conditx, condity, 
	    bunch, bunchphase, emod, emodphase, xlamds, prad0, zrayl, zwaist, 
	    rmax0, zsep, delz, zstop, quadf, quadd, fl, dl, drl, f1st, qfdx, 
	    qfdy, sl, solen, curlen, shotnoise, svar, dgrid, eloss, version, 
	    ibfield, imagl, idril;
    integer iseed, nwig, nsec, npart, ncar, lbc, nscr, nscz, nptr, ildgam, 
	    ildpsi, ildx, ildy, ildpx, ildpy, itgaus, nbins, iphsty, ishsty, 
	    ippart, ispart, ipradi, isradi, iertyp, iwityp, idump, iotail, 
	    nharm, magin, magout, lout[35], ffspec, ntail, nslice, iall, itdp,
	     ipseed, iscan, nscan, isntyp, isravg, isrsig, iorb, ndcut, 
	    idmppar, idmpfld, ilog, igamgaus, convharm, alignradf, offsetradf,
	     multconv;
    char beamfile[30], fieldfile[30], maginfile[30], magoutfile[30], 
	    outputfile[30], inputfile[30], scan[30], distfile[30], partfile[
	    30], filetype[30], radfile[30];
} inputcom_;

#define inputcom_1 inputcom_

Extern struct {
    //doublereal xpart[70001], ypart[70001], px[70001], py[70001], gamma[70001],
	   //  theta[70001], xporb[70001], yporb[70001], btpar[70001], btper[
	   // 70001], ez[70001], wx[70001], wy[70001], xcuren, dedz, tdmin, 
	   // tdmax, delcharge, dtd, charge;
    //integer lostid[70001], lost, losttot, ipos[280004]	/* was [4][70001] */;
	doublereal *xpart, *ypart, *px, *py, *gamma,
			*theta, *xporb, *yporb, *btpar, *btper, *ez, *wx, *wy, 
			xcuren, dedz, tdmin, tdmax, delcharge, dtd, charge;
	integer *lostid, lost, losttot, *ipos;
} beamcom_;

#define beamcom_1 beamcom_

Extern struct {
    //doublecomplex crwork3[116964], cpart1[70001], cpart2[70001];
    //doublereal k2gg[70001], k2pp[70001], k3gg[70001], k3pp[70001], p1[70001], 
	   // p2[70001];
	doublecomplex *crwork3, *cpart1, *cpart2;
	doublereal *k2gg, *k2pp, *k3gg, *k3pp, *p1, *p2;
} workspace_;

#define workspace_1 workspace_

/* Table of constant values */

static integer c__1 = 1;
static integer c_n16 = -16;

/* Subroutine */ int pushp_(istepz, xkper0)
integer *istepz;
doublereal *xkper0;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int chk_loss__();
    static doublereal stpz;
    extern /* Subroutine */ int partsorc_();
    static integer n;
    extern /* Subroutine */ int track_(), partsim_();

/*     ================================================================== */
/*     advance the particles, using runge-kutta (fourth order) method */
/*     ------------------------------------------------------------------ */





/*     error codes */

/* genesis version */
/* platform */
/* indicator for original fil */
/* indicator for sdds filetyp */
/* # of particles */
/* # of integration steps */
/* # of slices */
/* maximum of harmonics */
/* maximum of particle i */
/* <> 0 keeps distribution in */
/* energy units (mc^2) in ev */
/* vacuum impedence in ohms */
/* speed of light * electron */
/* pi */
/* pi/2 */
/* 2*pi */
/* check i for precission */
/* check ii for precission */
/* number of radial points fo */
/* # of gridpoints of cartesi */


/*     function prototypes */


/*     ------------------------------------------------------ */
/*     all input variables */

/*     wiggler */
/*     electron beam */
/*     radiation */
/*     grid-quantities */
/*     control */
/*     strong focusing */
/*     loading */
/*     output */
/*     external files */
/*     time-dependency */
/*     scan */
/*     extension */


/*     ------------------------------------------------------------------ */
/*     electron beam */





/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







    track_(istepz, xkper0);

/* advance particle transversly h */
    partsorc_(istepz);

/*     first step */
/*     ------------------------------------------------------------------ */

/* get source term at z0+delz/2 */
    i__1 = inputcom_1.npart;
    for (n = 1; n <= i__1; ++n) {
/* clear work arrays */
	workspace_1.k2gg[n - 1] = (float)0.;
	workspace_1.k2pp[n - 1] = (float)0.;
    }

/*     ode at z0 */

/* n */
    partsim_(beamcom_1.gamma, beamcom_1.theta, workspace_1.k2gg, 
	    workspace_1.k2pp);

/*     second step */
/*     ------------------------------------------------------------------ */

    stpz = inputcom_1.delz * (float).5 * 6.28318530717958;
    i__1 = inputcom_1.npart;
    for (n = 1; n <= i__1; ++n) {
/* k1 = d/2 k2 + k1 */
	beamcom_1.gamma[n - 1] = stpz * workspace_1.k2gg[n - 1] + 
		beamcom_1.gamma[n - 1];
	beamcom_1.theta[n - 1] = stpz * workspace_1.k2pp[n - 1] + 
		beamcom_1.theta[n - 1];
/* k3 = k2 */
	workspace_1.k3gg[n - 1] = workspace_1.k2gg[n - 1];
	workspace_1.k3pp[n - 1] = workspace_1.k2pp[n - 1];
/* k2 = 0 */
	workspace_1.k2gg[n - 1] = (float)0.;
	workspace_1.k2pp[n - 1] = (float)0.;
    }
/* ode at z0+delz/2 */
/* n */
    partsim_(beamcom_1.gamma, beamcom_1.theta, workspace_1.k2gg, 
	    workspace_1.k2pp);

/*     third step */
/*     ------------------------------------------------------------------ */

    i__1 = inputcom_1.npart;
    for (n = 1; n <= i__1; ++n) {
/* k1 = d/2 k2 + k1 */
	beamcom_1.gamma[n - 1] = stpz * workspace_1.k2gg[n - 1] + 
		beamcom_1.gamma[n - 1];
	beamcom_1.theta[n - 1] = stpz * workspace_1.k2pp[n - 1] + 
		beamcom_1.theta[n - 1];
/* k1 = -d/2 k3 + k1 */
	beamcom_1.gamma[n - 1] = -stpz * workspace_1.k3gg[n - 1] + 
		beamcom_1.gamma[n - 1];
	beamcom_1.theta[n - 1] = -stpz * workspace_1.k3pp[n - 1] + 
		beamcom_1.theta[n - 1];
/* k3 = k2/6 */
	workspace_1.k3gg[n - 1] /= (float)6.;
	workspace_1.k3pp[n - 1] /= (float)6.;
/* k2 = -k2/2 */
	workspace_1.k2gg[n - 1] = -workspace_1.k2gg[n - 1] / (float)2.;
	workspace_1.k2pp[n - 1] = -workspace_1.k2pp[n - 1] / (float)2.;
    }

/* ode at z0+delz/2 */

/* n */
    partsim_(beamcom_1.gamma, beamcom_1.theta, workspace_1.k2gg, 
	    workspace_1.k2pp);

/*     fourth step */
/*     ------------------------------------------------------------------ */

    stpz = inputcom_1.delz * 6.28318530717958;
    i__1 = inputcom_1.npart;
    for (n = 1; n <= i__1; ++n) {
/* k1 = d k2 + k1 */
	beamcom_1.gamma[n - 1] = stpz * workspace_1.k2gg[n - 1] + 
		beamcom_1.gamma[n - 1];
	beamcom_1.theta[n - 1] = stpz * workspace_1.k2pp[n - 1] + 
		beamcom_1.theta[n - 1];
/* k3 = -k2 + k3 */
	workspace_1.k3gg[n - 1] = -workspace_1.k2gg[n - 1] + workspace_1.k3gg[
		n - 1];
	workspace_1.k3pp[n - 1] = -workspace_1.k2pp[n - 1] + workspace_1.k3pp[
		n - 1];
/* k2 = -2 k2 */
	workspace_1.k2gg[n - 1] *= (float)2.;
	workspace_1.k2pp[n - 1] *= (float)2.;
    }

/* ode at z0+delz */

/* n */
    partsim_(beamcom_1.gamma, beamcom_1.theta, workspace_1.k2gg, 
	    workspace_1.k2pp);

    i__1 = inputcom_1.npart;
    for (n = 1; n <= i__1; ++n) {
	beamcom_1.gamma[n - 1] += stpz * (workspace_1.k3gg[n - 1] + 
		workspace_1.k2gg[n - 1] / (float)6.);
	beamcom_1.theta[n - 1] += stpz * (workspace_1.k3pp[n - 1] + 
		workspace_1.k2pp[n - 1] / (float)6.);
    }

/* n */
    track_(istepz, xkper0);

/*     wakefields or other energy losses if selected */

/* advance particle transversly hal */
    if (beamcom_1.dedz != (float)0.) {
	i__1 = inputcom_1.npart;
	for (n = 1; n <= i__1; ++n) {
	    beamcom_1.gamma[n - 1] += beamcom_1.dedz;
	}
    }

/*     check for particles outside the gridd */

    chk_loss__();

    return 0;
} /* pushp_ */




/* pushp */
/* Subroutine */ int chk_loss__()
{
    /* Format strings */
    static char fmt_100[] = "(f4.0)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi();

    /* Local variables */
    static integer idel;
    static doublereal rtmp;
    extern integer printerr_();
    static integer i__, j, k;
    static char closs[30];
    static integer delip[70001], mpart;

    /* Fortran I/O blocks */
    static icilist io___11 = { 0, closs, 0, fmt_100, 30, 1 };


/*     ======================================================================== */
/*     checks for lost particles, reorganizing the particle arrays */
/*     ------------------------------------------------------------------------ */





/*     error codes */

/* genesis version */
/* platform */
/* indicator for original fil */
/* indicator for sdds filetyp */
/* # of particles */
/* # of integration steps */
/* # of slices */
/* maximum of harmonics */
/* maximum of particle i */
/* <> 0 keeps distribution in */
/* energy units (mc^2) in ev */
/* vacuum impedence in ohms */
/* speed of light * electron */
/* pi */
/* pi/2 */
/* 2*pi */
/* check i for precission */
/* check ii for precission */
/* number of radial points fo */
/* # of gridpoints of cartesi */


/*     function prototypes */


/*     ------------------------------------------------------ */
/*     all input variables */

/*     wiggler */
/*     electron beam */
/*     radiation */
/*     grid-quantities */
/*     control */
/*     strong focusing */
/*     loading */
/*     output */
/*     external files */
/*     time-dependency */
/*     scan */
/*     extension */



/*     ------------------------------------------------------------------ */
/*     electron beam */





    if (beamcom_1.lost <= 0) {
	return 0;
    }

/* no loss-initialized in cut-tail */
    mpart = inputcom_1.npart / inputcom_1.nbins;
    j = 0;

    i__1 = inputcom_1.nbins - 2;
    for (k = 0; k <= i__1; ++k) {
	i__2 = mpart;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* run over one set of mirro */
	    if (beamcom_1.lostid[i__ + k * mpart - 1] != 0) {
/* lost ? */
		++j;
/* count & */
		delip[j - 1] = i__;
/* get index */
	    }
	}
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* make sure that the found */
	    beamcom_1.lostid[delip[i__ - 1] + (k + 1) * mpart - 1] = 0;
/* are not countet in next s */
	}
    }

    i__1 = mpart;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* search last set */
	if (beamcom_1.lostid[i__ + (inputcom_1.nbins - 1) * mpart - 1] != 0) {
	    ++j;
	    delip[j - 1] = i__;
	}
    }

    i__1 = j;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = inputcom_1.nbins - 1;
	for (k = 0; k <= i__2; ++k) {
	    beamcom_1.gamma[delip[i__ - 1] + k * mpart - 1] = (float)-1.;
	}
    }

    idel = 0;
    i__1 = inputcom_1.npart;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (beamcom_1.gamma[i__ - 1] > (float)0.) {
	    beamcom_1.gamma[i__ - idel - 1] = beamcom_1.gamma[i__ - 1];
	    beamcom_1.theta[i__ - idel - 1] = beamcom_1.theta[i__ - 1];
	    beamcom_1.xpart[i__ - idel - 1] = beamcom_1.xpart[i__ - 1];
	    beamcom_1.ypart[i__ - idel - 1] = beamcom_1.ypart[i__ - 1];
	    beamcom_1.px[i__ - idel - 1] = beamcom_1.px[i__ - 1];
	    beamcom_1.py[i__ - idel - 1] = beamcom_1.py[i__ - 1];
	} else {
	    ++idel;
	}
	beamcom_1.lostid[i__ - 1] = 0;
/* clear flags of lost particles */
    }

    beamcom_1.lost = 0;

/*     get numbers right */

    inputcom_1.npart -= inputcom_1.nbins * j;
    beamcom_1.xcuren = beamcom_1.xcuren * (real) inputcom_1.npart / (real) (
	    inputcom_1.npart + inputcom_1.nbins * j);
    rtmp = 1. - (doublereal) inputcom_1.npart / (doublereal) (
	    inputcom_1.npart + inputcom_1.nbins * j);
    if (rtmp > (float).01) {
	s_wsfi(&io___11);
	d__1 = rtmp * (float)100.;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfi();
	i__ = printerr_(&c_n16, closs, (ftnlen)30);
    }

    return 0;
} /* chk_loss__ */

