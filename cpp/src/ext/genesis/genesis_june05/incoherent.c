/* incoherent.f -- translated by f2c (version 20000118).
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

/* Subroutine */ int incoherent_(awz)
doublereal *awz;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer k, mpart, ip;
    static doublereal dgamavg, dgamsig, gam0;
    extern doublereal ran1_();
    static doublereal tmp2, xkw0;

/*     ================================================================== */
/*     compute elnergy lost and spread due to synchrtron radiation */
/*     only correct for planar undulators (almost correct for helical) */
/*     from  draft by s. reichle, modified by p. elleaume  august 1999 */
/*     the effect ofenergy loss and spread  has been tested successfully */
/*     with respect to e. saldin nima381 (1996) p 545-547 */
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





    if (*awz < 1e-25) {
	return 0;
    }
/* drift section! */
    if (inputcom_1.isravg == 0 && inputcom_1.isrsig == 0) {
	return 0;
    }

    xkw0 = 6.28318530717958 / inputcom_1.xlamd;
    gam0 = (float)0.;
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	gam0 += beamcom_1.gamma[ip - 1];
    }
    gam0 /= (doublereal) inputcom_1.npart;

/*     increase of energy spread */

/* Computing 2nd power */
    d__1 = xkw0 * *awz;
/* Computing 2nd power */
    d__2 = *awz;
    dgamsig = d__1 * d__1 * 1.015e-27 * (*awz * (float)1.697 + (float)1. / (*
	    awz * (float)1.88 + (float)1. + d__2 * d__2 * (float).8));
    if (inputcom_1.iwityp != 0) {
/* helical undula */
/* Computing 2nd power */
	d__1 = xkw0 * *awz;
/* Computing 2nd power */
	d__2 = *awz;
	dgamsig = d__1 * d__1 * 1.015e-27 * (*awz * (float)1.42 + (float)1. / 
		(*awz * (float)1.5 + (float)1. + d__2 * d__2 * (float).95));
    }

/* Computing 4th power */
    d__1 = gam0, d__1 *= d__1;
    dgamsig = sqrt(dgamsig * (d__1 * d__1) * xkw0 * inputcom_1.xlamd * 
	    inputcom_1.delz) * sqrt((float)3.) * (doublereal) 
	    inputcom_1.isrsig;

/*     average energy loss */

/* sqrt(3) uniform distribution */
/* Computing 2nd power */
    d__1 = xkw0 * gam0 * *awz;
    dgamavg = d__1 * d__1 * 1.88e-15 * inputcom_1.delz * inputcom_1.xlamd * (
	    doublereal) inputcom_1.isravg;

    mpart = (integer) (inputcom_1.npart / (doublereal) inputcom_1.nbins);
    gam0 = (float)0.;
    if (inputcom_1.isrsig != 0) {
	i__1 = mpart;
	for (ip = 1; ip <= i__1; ++ip) {
	    tmp2 = (ran1_(&inputcom_1.iseed) * (float)2. - (float)1.) * 
		    dgamsig;
	    gam0 += tmp2;
	    i__2 = inputcom_1.nbins - 1;
	    for (k = 0; k <= i__2; ++k) {
		beamcom_1.gamma[ip + k * mpart - 1] += tmp2;
	    }
	}
    }
    gam0 /= (doublereal) mpart;
    if (gam0 + dgamavg != (float)0.) {
	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
	    beamcom_1.gamma[ip - 1] = beamcom_1.gamma[ip - 1] - gam0 - 
		    dgamavg;
	}
    }

    return 0;
} /* incoherent_ */

