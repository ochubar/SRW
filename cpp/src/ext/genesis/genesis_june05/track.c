/* track.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

Extern struct {
	doublereal aw0,	xkx, xky, wcoefz[3], xlamd,	fbess0,	delaw, awd,	awx, awy, 
		gamma0,	delgam,	rxbeam,	rybeam,	alphax,	alphay,	emitx, emity, 
		xbeam, ybeam, pxbeam, pybeam, cuttail, curpeak,	conditx, condity, 
		bunch, bunchphase, emod, emodphase,	xlamds,	prad0, zrayl, zwaist, 
		rmax0, zsep, delz, zstop, quadf, quadd,	fl,	dl,	drl, f1st, qfdx, 
		qfdy, sl, solen, curlen, shotnoise,	svar, dgrid, eloss,	version, 
		ibfield, imagl,	idril;
	integer	iseed, nwig, nsec, npart, ncar,	lbc, nscr, nscz, nptr, ildgam, 
		ildpsi,	ildx, ildy,	ildpx, ildpy, itgaus, nbins, iphsty, ishsty, 
		ippart,	ispart,	ipradi,	isradi,	iertyp,	iwityp,	idump, iotail, 
		nharm, magin, magout, lout[35],	ffspec,	ntail, nslice, iall, itdp,
		ipseed,	iscan,	nscan, isntyp, isravg, isrsig, iorb, ndcut,	
		idmppar, idmpfld, ilog,	igamgaus, convharm,	alignradf, offsetradf,
		multconv;
	char beamfile[30], fieldfile[30], maginfile[30], magoutfile[30], 
		outputfile[30], inputfile[30], scan[30], distfile[30], partfile[30], filetype[30], radfile[30];
} inputcom_;

#define inputcom_1 inputcom_

Extern struct {
	doublereal *awz, *awdz, *solz, *awerx, *awery, *qfld, 
		fbess, magversion, unitlength, *dqfx, *dqfy, *awdx, *awdy;
	integer iran, nstepz, itap;
} wigcom_;

#define wigcom_1 wigcom_

Extern struct {
	doublereal *xpart, *ypart, *px, *py, *gamma,
		*theta, *xporb, *yporb, *btpar, *btper, *ez, *wx, *wy, 
		xcuren, dedz, tdmin, tdmax, delcharge, dtd, charge;
	integer *lostid, lost, losttot, *ipos;
} beamcom_;

#define beamcom_1 beamcom_

Extern struct {
	doublecomplex *crwork3, *cpart1, *cpart2;
	doublereal *k2gg, *k2pp, *k3gg, *k3pp, *p1, *p2;
} workspace_;

#define workspace_1 workspace_

/* Subroutine */ int track_(istepz, xkper0)
integer *istepz;
doublereal *xkper0;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(), cos(), sin(), cosh(), sinh();

    /* Local variables */
    static doublereal xoff, yoff, a1, a2, a3;
    static integer ip;
    static doublereal qx, qy, betpar0, foc, omg;

/*     ================================================================== */
/*     calculates exact soultion for transverse motion */
/*     this subroutine is call before and after runge-kutta integration */
/*     of phase and energy */
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
/*     wiggler parameters */




/*     ------------------------------------------------------------------ */
/*     electron beam */





/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







    inputcom_1.delz *= (float).5;

/* Computing 2nd power */
    d__1 = wigcom_1.awz[*istepz];
    betpar0 = d__1 * d__1 + (float)1.;
/* add weak focusing */
    betpar0 = sqrt(1. - betpar0 / inputcom_1.gamma0 / inputcom_1.gamma0);
/* Computing 2nd power */
    d__1 = *xkper0;
/* Computing 2nd power */
    d__2 = wigcom_1.awz[*istepz];
    qx = wigcom_1.qfld[*istepz] + inputcom_1.xkx * (d__1 * d__1) * (d__2 * 
	    d__2) / inputcom_1.gamma0 / betpar0;
/* Computing 2nd power */
    d__1 = *xkper0;
/* Computing 2nd power */
    d__2 = wigcom_1.awz[*istepz];
    qy = -wigcom_1.qfld[*istepz] + inputcom_1.xky * (d__1 * d__1) * (d__2 * 
	    d__2) / inputcom_1.gamma0 / betpar0;

/*     calculate magnetic center of quadrupole */

    xoff = 0.;
    yoff = 0.;

/* sven   the extra factor xkper0 comes from the normalization of x and y */

    if (abs(qx) > 1e-7) {
/* Computing 2nd power */
	d__1 = wigcom_1.awz[*istepz];
	xoff = wigcom_1.awdx[*istepz] * inputcom_1.xkx * *xkper0 * *xkper0 * (
		d__1 * d__1) / inputcom_1.gamma0 / betpar0;
	xoff = (xoff + wigcom_1.qfld[*istepz] * wigcom_1.dqfx[*istepz] * *
		xkper0) / qx;
    }
    if (abs(qy) > 1e-7) {
/* Computing 2nd power */
	d__1 = wigcom_1.awz[*istepz];
	yoff = wigcom_1.awdy[*istepz] * inputcom_1.xky * *xkper0 * *xkper0 * (
		d__1 * d__1) / inputcom_1.gamma0 / betpar0;
	yoff = (yoff - wigcom_1.qfld[*istepz] * wigcom_1.dqfy[*istepz] * *
		xkper0) / qy;
    }
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	workspace_1.p1[ip - 1] = beamcom_1.xpart[ip - 1] - xoff;
/* position relative to magnetic center */
	workspace_1.p2[ip - 1] = beamcom_1.ypart[ip - 1] - yoff;
/* of quadrupole field */
    }

    if (qx == (float)0.) {
	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
	    beamcom_1.xpart[ip - 1] += beamcom_1.px[ip - 1] * 
		    6.28318530717958 * inputcom_1.delz / beamcom_1.gamma[ip - 
		    1] / beamcom_1.btpar[ip - 1];
	}
    } else {
	if (qx > (float)0.) {
	    i__1 = inputcom_1.npart;
	    for (ip = 1; ip <= i__1; ++ip) {
		foc = sqrt(abs(qx) / beamcom_1.gamma[ip - 1] / 
			beamcom_1.btpar[ip - 1]);
		omg = foc * inputcom_1.delz * inputcom_1.xlamd;
		a1 = cos(omg);
		a2 = sin(omg) / foc;
		a3 = -a2 * foc * foc;
		beamcom_1.xpart[ip - 1] = a1 * workspace_1.p1[ip - 1] + a2 * 
			beamcom_1.px[ip - 1] / beamcom_1.gamma[ip - 1] / 
			beamcom_1.btpar[ip - 1] * *xkper0;
		beamcom_1.xpart[ip - 1] += xoff;
		beamcom_1.px[ip - 1] = a3 * workspace_1.p1[ip - 1] * 
			beamcom_1.gamma[ip - 1] * beamcom_1.btpar[ip - 1] / *
			xkper0 + a1 * beamcom_1.px[ip - 1];
	    }
	} else {
	    i__1 = inputcom_1.npart;
	    for (ip = 1; ip <= i__1; ++ip) {
		foc = sqrt(abs(qx) / beamcom_1.gamma[ip - 1] / 
			beamcom_1.btpar[ip - 1]);
		omg = foc * inputcom_1.delz * inputcom_1.xlamd;
		a1 = cosh(omg);
		a2 = sinh(omg) / foc;
		a3 = a2 * foc * foc;
		beamcom_1.xpart[ip - 1] = a1 * workspace_1.p1[ip - 1] + a2 * 
			beamcom_1.px[ip - 1] / beamcom_1.gamma[ip - 1] / 
			beamcom_1.btpar[ip - 1] * *xkper0;
		beamcom_1.xpart[ip - 1] += xoff;
		beamcom_1.px[ip - 1] = a3 * workspace_1.p1[ip - 1] * 
			beamcom_1.gamma[ip - 1] * beamcom_1.btpar[ip - 1] / *
			xkper0 + a1 * beamcom_1.px[ip - 1];
	    }
	}
    }

/*     field error */

    if ((d__1 = wigcom_1.awerx[*istepz - 1], abs(d__1)) > 1e-25) {
	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
	    beamcom_1.px[ip - 1] += wigcom_1.awerx[*istepz - 1] * 
		    inputcom_1.delz * 6.28318530717958;
/* sven            xpart(ip)=xpart(ip)+                         !kick at 0.5*delz */
/* sven     +                awerx(istepz)*0.5*delz*twopi/gamma(ip)/btpar(ip) */
	}
    }

/*     solenoid field */

    if ((d__1 = wigcom_1.solz[*istepz - 1], abs(d__1)) > 1e-25) {
	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
	    a1 = wigcom_1.solz[*istepz - 1] * inputcom_1.delz * 
		    inputcom_1.xlamd / beamcom_1.gamma[ip - 1];
	    beamcom_1.px[ip - 1] += a1 * beamcom_1.py[ip - 1];
	    beamcom_1.xpart[ip - 1] += a1 * (float).5 / beamcom_1.gamma[ip - 
		    1] / beamcom_1.btpar[ip - 1];
/* kick at 0.5* */
	}
    }

/*     and now for the y-plane */

    if (qy == (float)0.) {
	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
	    beamcom_1.ypart[ip - 1] += beamcom_1.py[ip - 1] * inputcom_1.delz 
		    * 6.28318530717958 / beamcom_1.gamma[ip - 1] / 
		    beamcom_1.btpar[ip - 1];
	}
    } else {
	if (qy > (float)0.) {
	    i__1 = inputcom_1.npart;
	    for (ip = 1; ip <= i__1; ++ip) {
		foc = sqrt(abs(qy) / beamcom_1.gamma[ip - 1] / 
			beamcom_1.btpar[ip - 1]);
		omg = foc * inputcom_1.delz * inputcom_1.xlamd;
		a1 = cos(omg);
		a2 = sin(omg) / foc;
		a3 = -a2 * foc * foc;
		beamcom_1.ypart[ip - 1] = a1 * workspace_1.p2[ip - 1] + a2 * 
			beamcom_1.py[ip - 1] / beamcom_1.gamma[ip - 1] / 
			beamcom_1.btpar[ip - 1] * *xkper0;
		beamcom_1.ypart[ip - 1] += yoff;
		beamcom_1.py[ip - 1] = a3 * workspace_1.p2[ip - 1] * 
			beamcom_1.gamma[ip - 1] * beamcom_1.btpar[ip - 1] / *
			xkper0 + a1 * beamcom_1.py[ip - 1];
	    }
	} else {
	    i__1 = inputcom_1.npart;
	    for (ip = 1; ip <= i__1; ++ip) {
		foc = sqrt(abs(qy) / beamcom_1.gamma[ip - 1] / 
			beamcom_1.btpar[ip - 1]);
		omg = foc * inputcom_1.delz * inputcom_1.xlamd;
		a1 = cosh(omg);
		a2 = sinh(omg) / foc;
		a3 = a2 * foc * foc;
		beamcom_1.ypart[ip - 1] = a1 * workspace_1.p2[ip - 1] + a2 * 
			beamcom_1.py[ip - 1] / beamcom_1.gamma[ip - 1] / 
			beamcom_1.btpar[ip - 1] * *xkper0;
		beamcom_1.ypart[ip - 1] += yoff;
		beamcom_1.py[ip - 1] = a3 * workspace_1.p2[ip - 1] * 
			beamcom_1.gamma[ip - 1] * beamcom_1.btpar[ip - 1] / *
			xkper0 + a1 * beamcom_1.py[ip - 1];
	    }
	}
    }

    if ((d__1 = wigcom_1.awery[*istepz - 1], abs(d__1)) > 1e-25) {
	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
	    beamcom_1.py[ip - 1] += wigcom_1.awery[*istepz - 1] * 
		    inputcom_1.delz * 6.28318530717958;
/* sven            ypart(ip)=ypart(ip) */
/* sven     +               +awery(istepz)*0.5*delz*twopi/gamma(ip)/btpar(ip) */
	}
    }


/*     solenoid field */

    if ((d__1 = wigcom_1.solz[*istepz - 1], abs(d__1)) > 1e-25) {
	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
	    a1 = wigcom_1.solz[*istepz - 1] * inputcom_1.delz * 
		    inputcom_1.xlamd / beamcom_1.gamma[ip - 1];
	    beamcom_1.py[ip - 1] -= a1 * beamcom_1.px[ip - 1];
	    beamcom_1.ypart[ip - 1] -= a1 * (float).5 / beamcom_1.gamma[ip - 
		    1] / beamcom_1.btpar[ip - 1];
/* kick at 0.5* */
	}
    }

    inputcom_1.delz *= (float)2.;

    return 0;
} /* track_ */

