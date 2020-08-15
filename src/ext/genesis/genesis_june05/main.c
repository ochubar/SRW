/* main.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

Extern struct {
    doublereal xkw0, xkper0, sval, svalout, gamma0_in__;
    integer npart0, inorun;
} simcom_;

#define simcom_1 simcom_

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
    //doublecomplex crfield[29241], crsource[29241], crmatc[171], cstep, crhm[29241], cwet[171], cbet[171];
	doublecomplex *crfield, *crsource, crmatc[513], cstep, *crhm, cwet[513], cbet[513];
    doublereal dxy, xks, radphase;
} cartcom_;

#define cartcom_1 cartcom_

Extern struct {
    //doublecomplex crwork3[116964], cpart1[70001], cpart2[70001];
    //doublereal k2gg[70001], k2pp[70001], k3gg[70001], k3pp[70001], p1[70001], 
	   // p2[70001];
	doublecomplex *crwork3, *cpart1, *cpart2;
	doublereal *k2gg, *k2pp, *k3gg, *k3pp, *p1, *p2;
} workspace_;

#define workspace_1 workspace_

Extern struct {
    //doublereal awz[10001], awdz[10000], solz[10000], awerx[10000], awery[
	   // 10000], qfld[10001], fbess, magversion, unitlength, dqfx[10001], 
	   // dqfy[10001], awdx[10001], awdy[10001];
	doublereal *awz, *awdz, *solz, *awerx, *awery, *qfld, 
			fbess, magversion, unitlength, *dqfx, *dqfy, *awdx, *awdy;
    integer iran, nstepz, itap;
} wigcom_;

#define wigcom_1 wigcom_

/*     ############################################################# */

/*             deusque dixit fiat lux et facta est lux */

/*     ############################################################# */

/*     genesis 1.3 is a time dependent fel code. the very basic structure */
/*     and the naming of most variables are taken from its precessor */
/*     tda3d. */

/*     particle integration:  runge-kutta 4th order */
/*     field integration:     alternating direction implicit-methode */
/*     space charge field:    inversion of tridiagonal matrix */

/*     genesis 1.3 is a cpu and memory expensive program and migth */
/*     exceed the requirement on older platforms. it can be partially */
/*     reduced by excluding time-dependent simulation or, as an */
/*     alternative, using a scratch file for most of the stored */
/*     information. */


/*     author:     sven reiche     (main algorithm, version 1.0) */
/*     co-author:  bart faatz      (external modules, version 1.0) */
/*                 pascal elleaume (quantum fluctuation) */
/*                 micheal borland (sdds) */
/*                 robert soliday  (sdds) */
/*                 Greg Penn       (conditioning, HGHG) */
/*                 Ati Meseck      (HGHG) */

/*     ############################################################### */

/*     genesis 1.3 was written in 1997 at desy, hamburg as a part */
/*     of my  ph.d.-thesis. i intended as an open-source project. */
/*     i'm willing to discuss with others about modifications, */
/*     found bugs and extensions which migth become official in the */
/*     release of the next version of genesis 1.3 */

/*     the current contact is */
/*              reiche@physics.ucla.edu */


/*     sven reiche, ucla, 08/24/04 */

/*     ################################################################# */
/*     ------------------------------------------------------------------ */
/*     main unit */
/*     ------------------------------------------------------------------ */
/* Main program */ MAIN__()
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    /* Subroutine */ int s_stop();

    /* Local variables */
    static integer islp, isep;
    extern /* Subroutine */ int loadbeam_(), chk_loss__(), last_(), stepz_(), 
	    swapfield_(), doscan_();
    static integer islice, lstepz, istepz;
    extern /* Subroutine */ int initio_(), dotime_(), output_(), loadrad_(), 
	    outglob_(), initrun_(), outhist_(), outdump_();


/* insert common blocks from file */




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

/*     simulation control and normalisation parameter */



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


/*     -------------------------------------------------------------------- */
/*     time dependency of bunch + slippage field */




/*     -------------------------------------------------------------------- */
/*     cartesian mesh */



/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







/*     ------------------------------------------------------------------ */
/*     wiggler parameters */





    initio_(); /* open files, read in input */
    initrun_(); /* initialize simulation */
    outglob_();

/*     start loop for each slice */

/* output of global parameter (t-independent) */
    i__1 = inputcom_1.nslice;
    for (islice = 1; islice <= i__1; ++islice) {

/*     initial loading */

/* loop for each slice */
	istepz = 0;

	doscan_(&islice);
/* update scan value */
	dotime_(&islice);

/* calculate time-dependent paramete */
	loadrad_(&islice);
/* radiation field loading */
	loadbeam_(&islice, &simcom_1.xkw0);
/* particle loading */
	chk_loss__();

/* remove cut particle */
	output_(&istepz, &islice, &simcom_1.xkw0);

/*       propagate beam for nu wiggler periods */

	i__2 = tbunchcom_1.nslp;
	for (islp = 1; islp <= i__2; ++islp) {
/* loop over slippage (advance field */
	    lstepz = tbunchcom_1.nsep;
	    if (islp == tbunchcom_1.nslp) {
		lstepz = wigcom_1.nstepz - (islp - 1) * tbunchcom_1.nsep;
	    }

/* correct for l */
	    i__3 = lstepz;
	    for (isep = 1; isep <= i__3; ++isep) {
/* loop 'steady state' simu */
		++istepz;
		stepz_(&istepz, &simcom_1.xkw0);
/* advance one step in z */
		output_(&istepz, &islice, &simcom_1.xkw0);
	    }

	    if (inputcom_1.itdp != 0 && islp < tbunchcom_1.nslp) {
		swapfield_(&islp);
/* advance field in time dep. s */
	    }

	}
	outhist_(&islice);
	outdump_(&islice);
/* dump rad + part if neede */
    }

    last_();
/* end run */
    s_stop("", (ftnlen)0);
} /* MAIN__ */

/* Main program alias */ int genesis_ () { MAIN__ (); return 0;}
