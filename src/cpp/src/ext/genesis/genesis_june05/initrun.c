/* initrun.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
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
	integer *lostid, lost, losttot, *ipos;	/* was [4][70001] */

	//double *pTest;

} beamcom_;

#define beamcom_1 beamcom_

Extern struct {
    doublereal xkw0, xkper0, sval, svalout, gamma0_in__;
    integer npart0, inorun;
} simcom_;

#define simcom_1 simcom_

Extern struct {
    //doublecomplex crfield[29241], crsource[29241], crmatc[171], cstep, crhm[29241], cwet[171], cbet[171];
	doublecomplex *crfield, *crsource, crmatc[513], cstep, *crhm, cwet[513], cbet[513];
    doublereal dxy, xks, radphase;
} cartcom_;

#define cartcom_1 cartcom_

Extern struct {
    //doublereal tgam0[30000], tdgam[30000], temitx[30000], temity[30000], 
	   // txrms[30000], tyrms[30000], txpos[30000], typos[30000], tpxpos[
	   // 30000], tpypos[30000], talphx[30000], talphy[30000], tcurrent[
	   // 30000], tpos[30000], tloss[30000], distgam[250000], distx[250000],
	   //  disty[250000], distpx[250000], distpy[250000], distt[250000], 
	   // tradpos[30000], tzrayl[30000], tzwaist[30000], tprad0[30000], 
	   // tradphase[30000];
	doublereal *tgam0, *tdgam, *temitx, *temity, *txrms, *tyrms, *txpos, 
			*typos, *tpxpos, *tpypos, *talphx, *talphy, *tcurrent, *tpos, *tloss, 
			*distgam, *distx, *disty, *distpx, *distpy, *distt,
			*tradpos, *tzrayl, *tzwaist, *tprad0, *tradphase;
    integer ndata, nsep, nslp, ndist, nraddata;
} tbunchcom_;

#define tbunchcom_1 tbunchcom_

Extern struct {
    //doublereal awz[10001], awdz[10000], solz[10000], awerx[10000], awery[
	   // 10000], qfld[10001], fbess, magversion, unitlength, dqfx[10001], 
	   // dqfy[10001], awdx[10001], awdy[10001];
	doublereal *awz, *awdz, *solz, *awerx, *awery, *qfld, 
			fbess, magversion, unitlength, *dqfx, *dqfy, *awdx, *awdy;
    integer iran, nstepz, itap;
} wigcom_;

#define wigcom_1 wigcom_

/* Table of constant values */

static integer c__1 = 1;
static integer c_n24 = -24;

/* Subroutine */ int initrun_()
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    extern /* Subroutine */ int magfield_(), last_(), scaninit_();
    extern integer printerr_();
    extern doublereal bessj0_(), bessj1_();
    static integer ip;
    static doublereal xi;
    extern /* Subroutine */ int loadslpfld_(), getdiag_();
    extern doublereal ran1_();

/*     ================================================================== */
/*     initialize the run by setting up */
/*     the precalculated matrices for the field solver and */
/*     claculating/normalizing some auxiliary variables */

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




/*     simulation control and normalisation parameter */


/*     -------------------------------------------------------------------- */
/*     cartesian mesh */



/*     -------------------------------------------------------------------- */
/*     time dependency of bunch + slippage field */





/*     ------------------------------------------------------------------ */
/*     wiggler parameters */





/*     seeding of random number generator */

    xi = ran1_(&inputcom_1.iseed);

/*     calculate coupling factor if necessary */

/* init ran1 */
    if (inputcom_1.fbess0 == (float)0.) {
	if (inputcom_1.iwityp == 0) {
/* Computing 2nd power */
	    d__1 = inputcom_1.aw0;
/* Computing 2nd power */
	    d__2 = inputcom_1.aw0;
	    xi = d__1 * d__1 / (d__2 * d__2 + 1.) / 2;
	    inputcom_1.fbess0 = bessj0_(&xi) - bessj1_(&xi);
	} else {
	    inputcom_1.fbess0 = (float)1.;
	}
    }

    beamcom_1.dedz = inputcom_1.eloss * inputcom_1.delz * inputcom_1.xlamd / 510999.06;
    beamcom_1.xcuren = inputcom_1.curpeak;
    simcom_1.npart0 = inputcom_1.npart;

/*     normalizations */
/*     ------------------------------------------------------------------ */

/* save total number or part */
    simcom_1.xkw0 = 6.28318530717958 / inputcom_1.xlamd;
/* wiggler wavenumber */
    simcom_1.xkper0 = simcom_1.xkw0;

/*     magnetic field */
/*     ------------------------------------------------------------------ */

/* transverse normalisation */
    magfield_(&simcom_1.xkper0, &c__1);

/* magnetic field descript */
    if (simcom_1.inorun != 0) {
	ip = printerr_(&c_n24, "Termination enforced by user", (ftnlen)28);
	last_();
    }

/*     slipping length */

    tbunchcom_1.nsep = (integer) (inputcom_1.zsep / inputcom_1.delz);
/* steps between field slip */
    tbunchcom_1.nslp = wigcom_1.nstepz / tbunchcom_1.nsep;
/* total slippage steps */
    if (wigcom_1.nstepz % tbunchcom_1.nsep != 0) {
	++tbunchcom_1.nslp;
    }
/* would be shorter */

/*     contruct grid properties (grid spacing, precalculated matrices) */

/* if not added the effecti */
	
    cartcom_1.dxy = simcom_1.xkw0 * 2. * inputcom_1.dgrid / (real) (inputcom_1.ncar - 1);
    cartcom_1.xks = 6.28318530717958 / inputcom_1.xlamds;

/*     time dependencies */

    loadslpfld_(&tbunchcom_1.nslp);

/*     scanning */

/* input field for first slice and seed */
    scaninit_();

/*     matrix initialization */

/* initialize scanning */
    d__1 = inputcom_1.delz * inputcom_1.xlamd;
    d__2 = cartcom_1.dxy / simcom_1.xkper0;
    getdiag_(&d__1, &d__2, &cartcom_1.xks);

/*     clear space charge field for case that space charge is disabled */

    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
/* clear space charge term */
	beamcom_1.ez[ip - 1] = 0.;
    }


/* ip */
    return 0;
} /* initrun_ */

#ifdef __cplusplus
}
#endif

