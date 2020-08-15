/* rpos.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

Extern struct {
    //doublecomplex crfield[29241], crsource[29241], crmatc[171], cstep, crhm[29241], cwet[171], cbet[171];
	doublecomplex *crfield, *crsource, crmatc[513], cstep, *crhm, cwet[513], cbet[513];
    doublereal dxy, xks, radphase;
} cartcom_;

#define cartcom_1 cartcom_

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
    //doublereal awz[10001], awdz[10000], solz[10000], awerx[10000], awery[
	   // 10000], qfld[10001], fbess, magversion, unitlength, dqfx[10001], 
	   // dqfy[10001], awdx[10001], awdy[10001];
	doublereal *awz, *awdz, *solz, *awerx, *awery, *qfld, 
			fbess, magversion, unitlength, *dqfx, *dqfy, *awdx, *awdy;
    integer iran, nstepz, itap;
} wigcom_;

#define wigcom_1 wigcom_

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
    doublereal xkw0, xkper0, sval, svalout, gamma0_in__;
    integer npart0, inorun;
} simcom_;

#define simcom_1 simcom_

/* Subroutine */ int rpos_(istepz, xx, yy)
integer *istepz;
doublereal *xx, *yy;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(), sin(), cos();

    /* Local variables */
    static doublereal x, awtmp, wxlow, wylow;
    static integer ip, ix1, iy1, ix2, iy2;
    extern doublereal faw2_();

/*     ================================================================== */
/*     locates the position of the electron on its actual trajectory */
/*     apply orbit correction to account for wiggle motion. */
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

/*     -------------------------------------------------------------------- */
/*     cartesian mesh */




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
/*     wiggler parameters */




/*     ------------------------------------------------------------------ */
/*     electron beam */





/*     simulation control and normalisation parameter */



    /* Parameter adjustments */
    --yy;
    --xx;

    /* Function Body */
    x = ((doublereal) (*istepz) + (float).5) * inputcom_1.delz * 
	    6.28318530717958;

/*     orbit correction ? */
/*     ------------------------------------------------------------------- */
    if (inputcom_1.iorb == 0) {
	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
/* no orbit correction */
	    beamcom_1.xporb[ip - 1] = xx[ip];
	    beamcom_1.yporb[ip - 1] = yy[ip];
	}
    } else {
	if (wigcom_1.fbess < (float)1.) {
	    i__1 = inputcom_1.npart;
	    for (ip = 1; ip <= i__1; ++ip) {
/* planar undulat */
		awtmp = sqrt(faw2_(istepz, &xx[ip], &yy[ip]));
/* aw at particle */
		beamcom_1.xporb[ip - 1] = xx[ip] - awtmp * sin(x) / 
			inputcom_1.gamma0;
		beamcom_1.yporb[ip - 1] = yy[ip];
	    }
	} else {
	    i__1 = inputcom_1.npart;
	    for (ip = 1; ip <= i__1; ++ip) {
/* helical undula */
		awtmp = sqrt(faw2_(istepz, &xx[ip], &yy[ip]));
/* aw at particle */
		beamcom_1.xporb[ip - 1] = xx[ip] - awtmp * sin(x) / 
			inputcom_1.gamma0;
		beamcom_1.yporb[ip - 1] = yy[ip] - awtmp * cos(x) / 
			inputcom_1.gamma0;
	    }
	}
    }

/*     linear interpolation */
/*     ----------------------------------------------------------- */
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	ix1 = (integer) ((beamcom_1.xporb[ip - 1] + inputcom_1.dgrid * 
		simcom_1.xkw0) / cartcom_1.dxy) + 1;
/* index in cart */
	iy1 = (integer) ((beamcom_1.yporb[ip - 1] + inputcom_1.dgrid * 
		simcom_1.xkw0) / cartcom_1.dxy) + 1;
	if (beamcom_1.xporb[ip - 1] < -inputcom_1.dgrid * simcom_1.xkw0) {
	    ix1 = 1;
	    ix2 = 1;
	    wxlow = 0.;
	    ++beamcom_1.lost;
	    beamcom_1.lostid[ip - 1] = 1;
	} else {
	    if (beamcom_1.xporb[ip - 1] >= inputcom_1.dgrid * simcom_1.xkw0) {
		ix1 = inputcom_1.ncar;
		ix2 = inputcom_1.ncar;
		wxlow = 1.;
		++beamcom_1.lost;
		beamcom_1.lostid[ip - 1] = 1;
	    } else {
		ix2 = ix1 + 1;
		wxlow = beamcom_1.xporb[ip - 1] + inputcom_1.dgrid * 
			simcom_1.xkw0 - cartcom_1.dxy * (real) (ix1 - 1);
		wxlow = 1. - wxlow / cartcom_1.dxy;
	    }
	}
	if (beamcom_1.yporb[ip - 1] < -inputcom_1.dgrid * simcom_1.xkw0) {
	    iy1 = 1;
	    iy2 = 1;
	    wylow = 0.;
	    ++beamcom_1.lost;
	    beamcom_1.lostid[ip - 1] = 1;
	} else {
	    if (beamcom_1.yporb[ip - 1] >= inputcom_1.dgrid * simcom_1.xkw0) {
		iy1 = inputcom_1.ncar;
		iy2 = inputcom_1.ncar;
		wylow = 1.;
		++beamcom_1.lost;
		beamcom_1.lostid[ip - 1] = 1;
	    } else {
		iy2 = iy1 + 1;
		wylow = beamcom_1.yporb[ip - 1] + inputcom_1.dgrid * 
			simcom_1.xkw0 - cartcom_1.dxy * (real) (iy1 - 1);
		wylow = 1. - wylow / cartcom_1.dxy;
	    }
	}
	beamcom_1.ipos[(ip << 2) - 4] = (iy1 - 1) * inputcom_1.ncar + ix1;
	beamcom_1.ipos[(ip << 2) - 3] = (iy2 - 1) * inputcom_1.ncar + ix1;
	beamcom_1.ipos[(ip << 2) - 2] = (iy1 - 1) * inputcom_1.ncar + ix2;
	beamcom_1.ipos[(ip << 2) - 1] = (iy2 - 1) * inputcom_1.ncar + ix2;
	beamcom_1.wx[ip - 1] = wxlow;
	beamcom_1.wy[ip - 1] = wylow;
    }

/* ip */
    return 0;
} /* rpos_ */


/* Subroutine */ int getpsi_(psi)
doublereal *psi;
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(), atan2();

    /* Local variables */
    static integer ip;
    static doublecomplex clocal;
    static doublereal philoc, wei;
    static integer idx;

/*     ================================================================== */
/*     calculates the total phase psi as the sum of the radiation phase phi */
/*     and the particla phase theta. */
/*     getpsi is only called by outpart to get a non moving bucket. */
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

/*     -------------------------------------------------------------------- */
/*     cartesian mesh */




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





    /* Parameter adjustments */
    --psi;

    /* Function Body */
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	wei = beamcom_1.wx[ip - 1] * beamcom_1.wy[ip - 1];
	idx = beamcom_1.ipos[(ip << 2) - 4];
	i__2 = idx - 1;
	z__1.r = wei * cartcom_1.crfield[i__2].r, z__1.i = wei * 
		cartcom_1.crfield[i__2].i;
	clocal.r = z__1.r, clocal.i = z__1.i;
	wei = beamcom_1.wx[ip - 1] * (1. - beamcom_1.wy[ip - 1]);
	idx = beamcom_1.ipos[(ip << 2) - 3];
	i__2 = idx - 1;
	z__2.r = wei * cartcom_1.crfield[i__2].r, z__2.i = wei * 
		cartcom_1.crfield[i__2].i;
	z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	clocal.r = z__1.r, clocal.i = z__1.i;
	wei = (1. - beamcom_1.wx[ip - 1]) * beamcom_1.wy[ip - 1];
	idx = beamcom_1.ipos[(ip << 2) - 2];
	i__2 = idx - 1;
	z__2.r = wei * cartcom_1.crfield[i__2].r, z__2.i = wei * 
		cartcom_1.crfield[i__2].i;
	z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	clocal.r = z__1.r, clocal.i = z__1.i;
	wei = (1. - beamcom_1.wx[ip - 1]) * (1. - beamcom_1.wy[ip - 1]);
	idx = beamcom_1.ipos[(ip << 2) - 1];
	i__2 = idx - 1;
	z__2.r = wei * cartcom_1.crfield[i__2].r, z__2.i = wei * 
		cartcom_1.crfield[i__2].i;
	z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	clocal.r = z__1.r, clocal.i = z__1.i;
	philoc = atan2(d_imag(&clocal), (doublereal) clocal.r);
	psi[ip] = philoc + beamcom_1.theta[ip - 1];
    }
/* ip */
    return 0;
} /* getpsi_ */

