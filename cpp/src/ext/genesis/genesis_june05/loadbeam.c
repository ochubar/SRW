/* loadbeam.f -- translated by f2c (version 20000118).
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
    doublereal distversion, distrev;
    integer iout[24], nout, nfld, npar, ndump, firstout, nfin, irecpar, 
	    irecfld, kout, ndmp2, npin, nlog, ndis, ncoldis, iconv2t, iconv2g,
	     icolpar[10], ndistsize, iconv2px, iconv2py, nprobe, ftype, 
	    ftdist, ftpart, ftfield;
} iocom_;

#define iocom_1 iocom_

Extern struct {
    //doublecomplex crwork3[116964], cpart1[70001], cpart2[70001];
    //doublereal k2gg[70001], k2pp[70001], k3gg[70001], k3pp[70001], p1[70001], 
	   // p2[70001];
	doublecomplex *crwork3, *cpart1, *cpart2;
	doublereal *k2gg, *k2pp, *k3gg, *k3pp, *p1, *p2;
} workspace_;

#define workspace_1 workspace_

/* Table of constant values */

static integer c_n23 = -23;

/* Subroutine */ int loadbeam_(islice, xkper0)
integer *islice;
doublereal *xkper0;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sin(), sqrt();

    /* Local variables */
    extern integer readpart_();
    extern /* Subroutine */ int loaddist_();
    extern integer printerr_();
    static integer i__;
    extern /* Subroutine */ int shotnoise_penman__();
    extern doublereal hammv_();
    static integer mpart;
    extern /* Subroutine */ int shotnoise_fawley__(), loadquiet_();
    static integer ip;

/*     =================================================================== */
/*     this routine fills the phase space for one slice of phase space */
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





/*     ------------------------------------------------------------------ */
/*     input/output control */





/*     initialize particle loss */

    beamcom_1.lost = 0;
    i__1 = inputcom_1.npart;
    for (i__ = 1; i__ <= i__1; ++i__) {
	beamcom_1.lostid[i__ - 1] = 0;
    }

    if (inputcom_1.iall != 0) {
	inputcom_1.ildpsi = -abs(inputcom_1.ildpsi);
/* reinitialize all hammersley sequences */
    }

/*     fill phase */

    mpart = inputcom_1.npart / inputcom_1.nbins;

/* particles per bin */
    i__1 = mpart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.theta[ip - 1] = hammv_(&inputcom_1.ildpsi) * (float)2. * 
		3.14159265358979 / (doublereal) inputcom_1.nbins - 
		3.14159265358979;
/* load in first bin */
    }

/*     branches for different loading methods */

    if (iocom_1.npin > 0) {
	i__ = readpart_(islice);
/* if reading from file */
	if (inputcom_1.convharm > 1 && inputcom_1.multconv != 0) {
	    shotnoise_fawley__();
	}
	return 0;
/* skip the rest (loading is done) */
    }

/*     btpar is calculated in readpart */

    if (iocom_1.ndis > 0) {
	loaddist_(&mpart, islice);
/* load from distribution file */
    } else {
	loadquiet_(&mpart);
/* internal load (quiet start) */
    }

/*     normalized transverse position */

    i__1 = mpart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.xpart[ip - 1] *= *xkper0;
	beamcom_1.ypart[ip - 1] *= *xkper0;
    }

/*     mirror particle in remaining bins */

    i__1 = inputcom_1.nbins - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = mpart;
	for (ip = 1; ip <= i__2; ++ip) {
	    beamcom_1.xpart[ip + i__ * mpart - 1] = beamcom_1.xpart[ip - 1];
	    beamcom_1.ypart[ip + i__ * mpart - 1] = beamcom_1.ypart[ip - 1];
	    beamcom_1.px[ip + i__ * mpart - 1] = beamcom_1.px[ip - 1];
	    beamcom_1.py[ip + i__ * mpart - 1] = beamcom_1.py[ip - 1];
	    beamcom_1.gamma[ip + i__ * mpart - 1] = beamcom_1.gamma[ip - 1];
	    beamcom_1.theta[ip + i__ * mpart - 1] = beamcom_1.theta[ip - 1] + 
		    (doublereal) i__ * (float)2. * 3.14159265358979 / (
		    doublereal) inputcom_1.nbins;
	    beamcom_1.lostid[ip + i__ * mpart - 1] = beamcom_1.lostid[ip - 1];
	}
    }

/*     add shotnoise */

    if (beamcom_1.xcuren < 0.) {
	i__ = printerr_(&c_n23, "xcuren <=0 in loadbeam", (ftnlen)22);
	beamcom_1.xcuren = 0.;
    } else {
	if (inputcom_1.isntyp == 0) {
	    shotnoise_fawley__();
	} else {
	    shotnoise_penman__();
	}
    }

/*     add longitudinal correlations (energy modulation, prebunching) */

    if (inputcom_1.bunch != (float)0. || inputcom_1.emod != (float)0.) {
	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
	    beamcom_1.gamma[ip - 1] -= inputcom_1.emod * sin(beamcom_1.theta[
		    ip - 1] - inputcom_1.emodphase);
	    beamcom_1.theta[ip - 1] -= inputcom_1.bunch * (float)2. * sin(
		    beamcom_1.theta[ip - 1] - inputcom_1.bunchphase);
/*              bunching is J_1(2.*bunch), approximately = bunch */
	}
    }

/*     calculate init. parallel velocity (needed in first call of track) */

    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
/* Computing 2nd power */
	d__1 = beamcom_1.px[ip - 1];
/* Computing 2nd power */
	d__2 = beamcom_1.py[ip - 1];
/* Computing 2nd power */
	d__3 = beamcom_1.gamma[ip - 1];
	beamcom_1.btpar[ip - 1] = sqrt(1. - (d__1 * d__1 + d__2 * d__2 + (
		float)1.) / (d__3 * d__3));
/* parallel velocity */
    }

    return 0;
} /* loadbeam_ */




/* loadbeam */
/* Subroutine */ int shotnoise_fawley__()
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double log(), sqrt(), sin();

    /* Local variables */
    static doublereal phin, enum__;
    static integer i__, j;
    static doublereal ecorr;
    static integer mpart, iharm;
    static doublereal an;
    static integer jj, ip;
    extern doublereal ran1_();

/*     ================================================================== */
/*     shotnoise algortihm following fawley */
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







    if (inputcom_1.shotnoise < 1e-25) {
	return 0;
    }

    mpart = inputcom_1.npart / inputcom_1.nbins;
    enum__ = beamcom_1.xcuren * inputcom_1.xlamds * inputcom_1.zsep / 
	    4.803302e-11 / (real) mpart;
/* #electron in slice */
    if (enum__ < 10.) {
	enum__ = 10.;
    }
/* catch low current error */
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	workspace_1.p1[ip - 1] = (float)0.;
    }

    i__1 = (inputcom_1.nbins - 1) / 2;
    for (iharm = 1; iharm <= i__1; ++iharm) {

/*       adjust shotnoise for imported distribution */

	ecorr = (float)1.;
	if (inputcom_1.convharm > 1 && iharm * inputcom_1.convharm <= 
		inputcom_1.nharm) {
	    ecorr = ((real) (inputcom_1.convharm * inputcom_1.convharm) - (
		    float)1.) / (real) (inputcom_1.convharm * 
		    inputcom_1.convharm);
	}
	i__2 = mpart;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    phin = ran1_(&inputcom_1.ipseed) * 6.28318530717958;
	    an = sqrt(-log(ran1_(&inputcom_1.ipseed)) * ecorr / enum__) * (
		    float)2. / (real) iharm;
	    i__3 = inputcom_1.nbins;
	    for (j = 1; j <= i__3; ++j) {
		jj = (j - 1) * mpart + i__;
		workspace_1.p1[jj - 1] -= an * sin((real) iharm * 
			beamcom_1.theta[jj - 1] + phin);
	    }
	}
    }

    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.theta[ip - 1] += workspace_1.p1[ip - 1] * 
		inputcom_1.shotnoise;
    }

    return 0;
} /* shotnoise_fawley__ */


/* of shotnoise_fawley */
/* Subroutine */ int shotnoise_penman__()
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(), sin();

    /* Local variables */
    static doublereal enum__, ratio;
    static integer ip;
    static doublereal snoise, sn1, sn2;
    extern doublereal ran1_();

/*     ================================================================== */
/*     shotnoise algorithm following Penman */
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





    if (inputcom_1.shotnoise < 1e-25) {
	return 0;
    }

    enum__ = beamcom_1.xcuren * inputcom_1.xlamds * inputcom_1.zsep / 
	    4.803302e-11;
/* #electron in */
    if (enum__ < 10.) {
	enum__ = 10.;
    }
/* catch low current error */
    ratio = (real) inputcom_1.npart / enum__;
    snoise = sqrt((real) inputcom_1.npart * 3. / enum__);

/*     for npart ~ enum the approximation is not precise enough */

/* shot noise parameter */
    if (ratio > (float).02) {
	sn1 = snoise;
	sn2 = 1.570796326794895;
L1:
	snoise = (sn1 + sn2) * (float).5;
/* Computing 2nd power */
	d__1 = sin(snoise);
	ratio = d__1 * d__1 / snoise / snoise - (float)1. + (real) 
		inputcom_1.npart / enum__;
	if (ratio < 0.) {
	    sn2 = snoise;
	} else {
	    sn1 = snoise;
	}
	if (sn2 - sn1 > (float).002) {
	    goto L1;
	}
    }

    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.theta[ip - 1] += snoise * (1. - ran1_(&inputcom_1.ipseed) * 
		2.);
    }

    return 0;
} /* shotnoise_penman__ */


/* of shotnoise_penman */
/* Subroutine */ int loadquiet_(mpart)
integer *mpart;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double log(), sqrt();

    /* Local variables */
    extern /* Subroutine */ int last_(), cut_tail__();
    static doublereal betaxinv, betayinv;
    extern integer printerr_();
    static doublereal x;
    extern doublereal hammv_();
    static doublereal y, r2;
    static integer ip;
    static doublereal xd;
    extern doublereal gasham_();
    static doublereal yd, xy, xx, yy, r2p, ampl_x__, ampl_y__, fac, pxd, pyd, 
	    xpd, ypd;
    extern /* Subroutine */ int compmom_();

/*     ================================================================== */
/*     do quiet loading of transverse phase space */
/*     ------------------------------------------------------------------- */





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

/*     ------------------------------------------------------------------ */
/*     electron beam */






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



/*     fill energy with standard diviation 1 */

    if (inputcom_1.igamgaus != 0) {
	i__1 = *mpart;
	for (ip = 1; ip <= i__1; ++ip) {
	    beamcom_1.gamma[ip - 1] = gasham_(&inputcom_1.ildgam);
/* gauss distribution */
	}
    } else {
	i__1 = *mpart;
	for (ip = 1; ip <= i__1; ++ip) {
	    beamcom_1.gamma[ip - 1] = hammv_(&inputcom_1.ildgam) * 2. - 1.;
/* uniform distribution */
	}
    }

    *mpart /= 2;

/*     fill x,y,px,py parabolic between [-1,1] */

/* only half particles per one bin, others symmetr */
    ip = 0;
L10:
    ++ip;
L20:
    beamcom_1.xpart[ip - 1] = hammv_(&inputcom_1.ildx) * 2. - 1.;
    beamcom_1.ypart[ip - 1] = hammv_(&inputcom_1.ildy) * 2. - 1.;
    beamcom_1.px[ip - 1] = hammv_(&inputcom_1.ildpx) * 2. - 1.;
    beamcom_1.py[ip - 1] = hammv_(&inputcom_1.ildpy) * 2. - 1.;
/* Computing 2nd power */
    d__1 = beamcom_1.xpart[ip - 1];
/* Computing 2nd power */
    d__2 = beamcom_1.ypart[ip - 1];
    r2 = d__1 * d__1 + d__2 * d__2;
    if (r2 >= (float)1.) {
	goto L20;
    }
/* Computing 2nd power */
    d__1 = beamcom_1.px[ip - 1];
/* Computing 2nd power */
    d__2 = beamcom_1.py[ip - 1];
    r2p = d__1 * d__1 + d__2 * d__2 + r2;
    if (r2p >= (float)1.) {
	goto L20;
    }
    if (ip < *mpart) {
	goto L10;
    }

/*     enforce profile */

    i__1 = *mpart;
    for (ip = 1; ip <= i__1; ++ip) {
/* Computing 2nd power */
	d__1 = beamcom_1.xpart[ip - 1];
/* Computing 2nd power */
	d__2 = beamcom_1.ypart[ip - 1];
/* Computing 2nd power */
	d__3 = beamcom_1.px[ip - 1];
/* Computing 2nd power */
	d__4 = beamcom_1.py[ip - 1];
	r2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4;
	if (r2 > 1e-25) {
	    fac = 1.;
/* parabolic */
	    if (inputcom_1.itgaus == 1) {
		fac = sqrt(log(1. - r2) * -1. / r2);
	    }
/* gaussian */
	    if (inputcom_1.itgaus == 2) {
		fac = 1. / sqrt(r2);
	    }
/* step */
	    beamcom_1.xpart[ip - 1] *= fac;
	    beamcom_1.ypart[ip - 1] *= fac;
	    beamcom_1.px[ip - 1] *= fac;
	    beamcom_1.py[ip - 1] *= fac;
	}
    }

/*     correct phasespace loading */

    i__1 = *mpart;
    for (ip = 1; ip <= i__1; ++ip) {
/* <x>=<y>=0 locally */
	beamcom_1.xpart[ip + *mpart - 1] = -beamcom_1.xpart[ip - 1];
	beamcom_1.ypart[ip + *mpart - 1] = -beamcom_1.ypart[ip - 1];
	beamcom_1.px[ip + *mpart - 1] = -beamcom_1.px[ip - 1];
	beamcom_1.py[ip + *mpart - 1] = -beamcom_1.py[ip - 1];
    }

    *mpart <<= 1;

/*     remove any correlation between x,px and y,py */

    compmom_(mpart, &x, &xx, &y, &yy, &xy, beamcom_1.xpart, beamcom_1.px);
    xd = (xy - x * y) / (xx - x * x);
    compmom_(mpart, &x, &xx, &y, &yy, &xy, beamcom_1.ypart, beamcom_1.py);
    yd = (xy - x * y) / (xx - x * x);
    i__1 = *mpart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.px[ip - 1] -= xd * beamcom_1.xpart[ip - 1];
	beamcom_1.py[ip - 1] -= yd * beamcom_1.ypart[ip - 1];
    }

/*     normalize distritution -> rms value is 1. */

    compmom_(mpart, &x, &xx, &y, &yy, &xy, beamcom_1.xpart, beamcom_1.ypart);
/* compute moments */
    xd = (float)1. / sqrt(xx);
    yd = (float)1. / sqrt(yy);
    i__1 = *mpart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.xpart[ip - 1] *= xd;
	beamcom_1.ypart[ip - 1] *= yd;
    }

    compmom_(mpart, &x, &xx, &y, &yy, &xy, beamcom_1.px, beamcom_1.py);
/* compute moments */
    xd = (float)1. / sqrt(xx);
    yd = (float)1. / sqrt(yy);
    i__1 = *mpart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.px[ip - 1] *= xd;
	beamcom_1.py[ip - 1] *= yd;
    }

    cut_tail__(mpart);

    compmom_(mpart, &x, &xx, &y, &yy, &xy, beamcom_1.gamma, beamcom_1.theta);
/* compute moments */
    xx -= x * x;
    xd = (float)1. / sqrt(xx);
    if (inputcom_1.igamgaus == 0) {
	xd /= sqrt((float)3.);
    }
    i__1 = *mpart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.gamma[ip - 1] = (beamcom_1.gamma[ip - 1] - x) * xd;
    }

/*     scale and shift distribution */

    ip = 0;
    if (inputcom_1.rxbeam <= 0.) {
	ip = printerr_(&c_n23, "rxbeam in loadquiet", (ftnlen)19);
    }
    if (inputcom_1.rybeam <= 0.) {
	ip = printerr_(&c_n23, "rybeam in loadquiet", (ftnlen)19);
    }
    if (inputcom_1.gamma0 <= 1.) {
	ip = printerr_(&c_n23, "gamma0 in loadquiet", (ftnlen)19);
    }
    if (ip < 0) {
	last_();
    }

/* Computing 2nd power */
    d__1 = inputcom_1.rxbeam;
    betaxinv = inputcom_1.emitx / (d__1 * d__1);
/* calculation of 1/betax => */
/* Computing 2nd power */
    d__1 = inputcom_1.rybeam;
    betayinv = inputcom_1.emity / (d__1 * d__1);
/* no sigularity for emit=0 */
    pxd = inputcom_1.rxbeam * betaxinv;
/* rms px at x=0 */
    pyd = inputcom_1.rybeam * betayinv;
/* rms py at y=0 */
    xpd = -inputcom_1.alphax * betaxinv;
/* slope of phase space elli */
    ypd = -inputcom_1.alphay * betayinv;
    i__1 = *mpart;
    for (ip = 1; ip <= i__1; ++ip) {
/*        define amplitudes, simplest before tilt and rescale is applied */
/* Computing 2nd power */
	d__1 = beamcom_1.px[ip - 1];
/* Computing 2nd power */
	d__2 = beamcom_1.xpart[ip - 1];
	ampl_x__ = inputcom_1.emitx * (float).5 * (d__1 * d__1 + d__2 * d__2);
/* Computing 2nd power */
	d__1 = beamcom_1.py[ip - 1];
/* Computing 2nd power */
	d__2 = beamcom_1.ypart[ip - 1];
	ampl_y__ = inputcom_1.emity * (float).5 * (d__1 * d__1 + d__2 * d__2);
	beamcom_1.px[ip - 1] = xpd * inputcom_1.rxbeam * beamcom_1.xpart[ip - 
		1] + beamcom_1.px[ip - 1] * pxd;
/* scale+shift px */
	beamcom_1.py[ip - 1] = ypd * inputcom_1.rybeam * beamcom_1.ypart[ip - 
		1] + beamcom_1.py[ip - 1] * pyd;
/* scale+shift py */
	beamcom_1.xpart[ip - 1] *= inputcom_1.rxbeam;
/* scale x to right size */
	beamcom_1.ypart[ip - 1] *= inputcom_1.rybeam;
/* scale y to right size */
	beamcom_1.gamma[ip - 1] = inputcom_1.gamma0 + inputcom_1.delgam * 
		beamcom_1.gamma[ip - 1];
	beamcom_1.gamma[ip - 1] = beamcom_1.gamma[ip - 1] + 
		inputcom_1.conditx * ampl_x__ + inputcom_1.condity * ampl_y__;
    }

/*     shift center in phase space */

    i__1 = *mpart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.xpart[ip - 1] += inputcom_1.xbeam;
	beamcom_1.ypart[ip - 1] += inputcom_1.ybeam;
	beamcom_1.px[ip - 1] += inputcom_1.pxbeam;
	beamcom_1.py[ip - 1] += inputcom_1.pybeam;
    }
    return 0;
} /* loadquiet_ */



/* Subroutine */ int loaddist_(mpart, islice)
integer *mpart, *islice;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer mget;
    extern /* Subroutine */ int neighbor_();
    static integer i__, n1, n2;
    static doublereal t0, t1;
    extern /* Subroutine */ int readslice_(), scaledist_();
    static doublereal ga, xa, ya, gs, xs, ys;
    extern /* Subroutine */ int switch_();
    static doublereal pxa, pya, pxs, pys;
    extern /* Subroutine */ int compmom_();
    extern doublereal ran1_();
    static doublereal tmp1, tmp2, tmp3;

/*     ================================================================= */
/*     load slice from distribution */
/*     ----------------------------------------------------------------- */





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





    t0 = beamcom_1.tdmax - (inputcom_1.ntail + *islice - 1) * inputcom_1.zsep 
	    * inputcom_1.xlamds / 3e8 - beamcom_1.dtd * (float).5;
/* note sign for */
    t1 = t0 + beamcom_1.dtd;
    readslice_(&mget, beamcom_1.xpart, beamcom_1.px, beamcom_1.ypart, 
	    beamcom_1.py, beamcom_1.gamma, &t0, &t1);
    beamcom_1.xcuren = beamcom_1.delcharge * (doublereal) mget / 
	    beamcom_1.dtd;

    if (mget < 5) {
	i__1 = *mpart;
	for (i__ = mget + 1; i__ <= i__1; ++i__) {
	    beamcom_1.xpart[i__ - 1] = 0.;
	    beamcom_1.ypart[i__ - 1] = (float)0.;
	    beamcom_1.px[i__ - 1] = (float)0.;
	    beamcom_1.py[i__ - 1] = (float)0.;
	    beamcom_1.gamma[i__ - 1] = inputcom_1.gamma0;
	}
	return 0;
    }

/*     adjust and normalize distirbution */

    compmom_(&mget, &xa, &xs, &pxa, &pxs, &tmp3, beamcom_1.xpart, 
	    beamcom_1.px);
    compmom_(&mget, &ya, &ys, &pya, &pys, &tmp3, beamcom_1.ypart, 
	    beamcom_1.py);
    compmom_(&mget, &ga, &gs, &tmp1, &tmp2, &tmp3, beamcom_1.gamma, 
	    beamcom_1.xpart);

/*     the absolute function should be redundent but sometime a round off */
/*     error can cause a negative number. the error was pointed out by Gregg Penn */

    xs = sqrt((d__1 = xs - xa * xa, abs(d__1)));
    pxs = sqrt((d__1 = pxs - pxa * pxa, abs(d__1)));
    ys = sqrt((d__1 = ys - ya * ya, abs(d__1)));
    pys = sqrt((d__1 = pys - pya * pya, abs(d__1)));
    gs = sqrt((d__1 = gs - ga * ga, abs(d__1)));
    if (xs < 1e-25) {
	xs = 1.;
    }
    if (pxs < 1e-25) {
	pxs = 1.;
    }
    if (ys < 1e-25) {
	ys = 1.;
    }
    if (pys < 1e-25) {
	pys = 1.;
    }
    if (gs < 1e-25) {
	gs = 1.;
    }
    d__1 = -xa;
    d__2 = 1. / xs;
    scaledist_(beamcom_1.xpart, &mget, &d__1, &d__2);
    d__1 = -pxa;
    d__2 = 1. / pxs;
    scaledist_(beamcom_1.px, &mget, &d__1, &d__2);
    d__1 = -ya;
    d__2 = 1. / ys;
    scaledist_(beamcom_1.ypart, &mget, &d__1, &d__2);
    d__1 = -pya;
    d__2 = 1. / pys;
    scaledist_(beamcom_1.py, &mget, &d__1, &d__2);
    d__1 = -ga;
    d__2 = 1. / gs;
    scaledist_(beamcom_1.gamma, &mget, &d__1, &d__2);

/*     massage proto-distribution */

    if (mget >= *mpart) {

/*         remove particles */

	i__1 = mget - 1;
	for (i__ = *mpart; i__ <= i__1; ++i__) {
	    n1 = (integer) (mget * ran1_(&inputcom_1.iseed)) + 1;
	    switch_(beamcom_1.xpart, beamcom_1.px, beamcom_1.ypart, 
		    beamcom_1.py, beamcom_1.gamma, &n1, &mget);
	    --mget;
	}
    } else {

/*        add particle */

L1:
	n1 = (integer) ((doublereal) mget * ran1_(&inputcom_1.iseed)) + 1;
	neighbor_(beamcom_1.xpart, beamcom_1.px, beamcom_1.ypart, 
		beamcom_1.py, beamcom_1.gamma, &n1, &n2, &mget, &
		inputcom_1.iseed);
	beamcom_1.xpart[mget] = beamcom_1.xpart[n1 - 1] + (beamcom_1.xpart[n2 
		- 1] - beamcom_1.xpart[n1 - 1]) * ((ran1_(&inputcom_1.iseed) 
		+ ran1_(&inputcom_1.iseed)) * (float).5);
	beamcom_1.px[mget] = beamcom_1.px[n1 - 1] + (beamcom_1.px[n2 - 1] - 
		beamcom_1.px[n1 - 1]) * ((ran1_(&inputcom_1.iseed) + ran1_(&
		inputcom_1.iseed)) * (float).5);
	beamcom_1.ypart[mget] = beamcom_1.ypart[n1 - 1] + (beamcom_1.ypart[n2 
		- 1] - beamcom_1.ypart[n1 - 1]) * ((ran1_(&inputcom_1.iseed) 
		+ ran1_(&inputcom_1.iseed)) * (float).5);
	beamcom_1.py[mget] = beamcom_1.py[n1 - 1] + (beamcom_1.py[n2 - 1] - 
		beamcom_1.py[n1 - 1]) * ((ran1_(&inputcom_1.iseed) + ran1_(&
		inputcom_1.iseed)) * (float).5);
	beamcom_1.gamma[mget] = beamcom_1.gamma[n1 - 1] + (beamcom_1.gamma[n2 
		- 1] - beamcom_1.gamma[n1 - 1]) * ((ran1_(&inputcom_1.iseed) 
		+ ran1_(&inputcom_1.iseed)) * (float).5);
	++mget;
	if (mget < *mpart) {
	    goto L1;
	}
    }

/*     scale back */

    d__1 = xa / xs;
    scaledist_(beamcom_1.xpart, mpart, &d__1, &xs);
    d__1 = pxa / pxs;
    scaledist_(beamcom_1.px, mpart, &d__1, &pxs);
    d__1 = ya / ys;
    scaledist_(beamcom_1.ypart, mpart, &d__1, &ys);
    d__1 = pya / pys;
    scaledist_(beamcom_1.py, mpart, &d__1, &pys);
    d__1 = ga / gs;
    scaledist_(beamcom_1.gamma, mpart, &d__1, &gs);

    return 0;
} /* loaddist_ */



/* Subroutine */ int compmom_(n, x, x2, y, y2, xy, xa, ya)
integer *n;
doublereal *x, *x2, *y, *y2, *xy, *xa, *ya;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal dd;

/*     ================================================================== */
/*     compute moments */
/*     ------------------------------------------------------------------ */


/* >=n */
    /* Parameter adjustments */
    --ya;
    --xa;

    /* Function Body */
    *x = 0.;
    *x2 = 0.;
    *y = 0.;
    *y2 = 0.;
    *xy = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*x += xa[i__];
	*y += ya[i__];
	*x2 += xa[i__] * xa[i__];
	*y2 += ya[i__] * ya[i__];
	*xy += xa[i__] * ya[i__];
    }
    dd = 1. / (doublereal) (*n);
/* enforce double precission */
    *x *= dd;
    *y *= dd;
    *x2 *= dd;
    *y2 *= dd;
    *xy *= dd;

    return 0;
} /* compmom_ */



/* compmom */
/* Subroutine */ int switch_(x, px, y, py, g, n1, n2)
doublereal *x, *px, *y, *py, *g;
integer *n1, *n2;
{
    static doublereal tmp;

/*     ================================================================= */
/*     switch two particles */
/*     ----------------------------------------------------------------- */


    /* Parameter adjustments */
    --g;
    --py;
    --y;
    --px;
    --x;

    /* Function Body */
    tmp = x[*n1];
    x[*n1] = x[*n2];
    x[*n2] = tmp;
    tmp = px[*n1];
    px[*n1] = px[*n2];
    px[*n2] = tmp;
    tmp = y[*n1];
    y[*n1] = y[*n2];
    y[*n2] = tmp;
    tmp = py[*n1];
    py[*n1] = py[*n2];
    py[*n2] = tmp;
    tmp = g[*n1];
    g[*n1] = g[*n2];
    g[*n2] = tmp;
    return 0;
} /* switch_ */

/* Subroutine */ int scaledist_(x, n, a, b)
doublereal *x;
integer *n;
doublereal *a, *b;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/*     ============================================== */
/*     scales distribution to x -> b*(x+a) */
/*     ---------------------------------------------- */


    /* Parameter adjustments */
    --x;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = *b * (x[i__] + *a);
    }
    return 0;
} /* scaledist_ */


/* Subroutine */ int neighbor_(x, px, y, py, g, n1, n2, n, iseed)
doublereal *x, *px, *y, *py, *g;
integer *n1, *n2, *n, *iseed;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal rmin;
    static integer i__;
    static doublereal r__;
    extern doublereal ran1_();

/*     ================================================================= */
/*     search for particle n2 closest to n1 */
/*     ----------------------------------------------------------------- */


    /* Parameter adjustments */
    --g;
    --py;
    --y;
    --px;
    --x;

    /* Function Body */
    rmin = (float)-1.;
    *n2 = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = x[*n1] - x[i__];
/* Computing 2nd power */
	d__2 = px[*n1] - px[i__];
	r__ = d__1 * d__1 * ran1_(iseed) + d__2 * d__2 * ran1_(iseed);
/* Computing 2nd power */
	d__1 = y[*n1] - y[i__];
/* Computing 2nd power */
	d__2 = py[*n1] - py[i__];
	r__ = r__ + d__1 * d__1 * ran1_(iseed) + d__2 * d__2 * ran1_(iseed);
/* Computing 2nd power */
	d__1 = g[*n1] - g[i__];
	r__ += d__1 * d__1 * ran1_(iseed);
	if (i__ != *n1) {
	    if (rmin < 0. || r__ < rmin) {
		*n2 = i__;
		rmin = r__;
	    }
	}
    }
    if (*n2 == 0) {
	*n2 = *n1;
/* copy himself if only one particle is present */
    }
    return 0;
} /* neighbor_ */


/* Subroutine */ int cut_tail__(mp)
integer *mp;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal r__;

/*     ======================================================= */
/*     collimation of the transverse tails */
/*     ---------------------------------------------------------- */





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

/*     ------------------------------------------------------------------ */
/*     electron beam */






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



    if (inputcom_1.cuttail <= 0.) {
	return 0;
    }

    i__1 = *mp;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = beamcom_1.xpart[i__ - 1];
/* Computing 2nd power */
	d__2 = beamcom_1.px[i__ - 1];
	r__ = d__1 * d__1 + d__2 * d__2;
	if (r__ > inputcom_1.cuttail) {
	    beamcom_1.lost = 1;
	    beamcom_1.lostid[i__ - 1] = 1;
	}
/* Computing 2nd power */
	d__1 = beamcom_1.ypart[i__ - 1];
/* Computing 2nd power */
	d__2 = beamcom_1.py[i__ - 1];
	r__ = d__1 * d__1 + d__2 * d__2;
	if (r__ > inputcom_1.cuttail) {
	    beamcom_1.lost = 1;
	    beamcom_1.lostid[i__ - 1] = 1;
	}
    }

    return 0;
} /* cut_tail__ */

