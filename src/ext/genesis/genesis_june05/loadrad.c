/* loadrad.f -- translated by f2c (version 20000118).
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
    //doublecomplex crwork3[116964], cpart1[70001], cpart2[70001];
    //doublereal k2gg[70001], k2pp[70001], k3gg[70001], k3pp[70001], p1[70001], 
	   // p2[70001];
	doublecomplex *crwork3, *cpart1, *cpart2;
	doublereal *k2gg, *k2pp, *k3gg, *k3pp, *p1, *p2;
} workspace_;

#define workspace_1 workspace_

Extern struct {
    doublereal distversion, distrev;
    integer iout[24], nout, nfld, npar, ndump, firstout, nfin, irecpar, 
	    irecfld, kout, ndmp2, npin, nlog, ndis, ncoldis, iconv2t, iconv2g,
	     icolpar[10], ndistsize, iconv2px, iconv2py, nprobe, ftype, 
	    ftdist, ftpart, ftfield;
} iocom_;

#define iocom_1 iocom_

Extern struct {
    //doublereal error[10000], gain[10000], phimid[10000], whalf[10000], logp[
	   // 10000], powmid[10000], xrms[10000], yrms[10000], xpos[10000], 
	   // ypos[10000], pmodhist[50000]	/* was [5][10000] */, gamhist[
	   // 10000], diver[10000], pradol, pinit, bunphase[50000]	/* 
	   // was [5][10000] */, dgamhist[10000], ffield[10000];
	doublereal *error, *gain, *phimid, *whalf, *logp, *powmid, *xrms, *yrms, *xpos, 
			*ypos, *pmodhist, *gamhist, *diver, pradol, pinit, 
			*bunphase, *dgamhist, *ffield;
	integer ihist;
} diagcom_;

#define diagcom_1 diagcom_

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
    doublereal xkw0, xkper0, sval, svalout, gamma0_in__;
    integer npart0, inorun;
} simcom_;

#define simcom_1 simcom_

Extern struct {
    //doublecomplex crtime[15960700];
	doublecomplex *crtime;
    integer ntmp;
} tslipcom_;

#define tslipcom_1 tslipcom_

/* Table of constant values */

static integer c_n23 = -23;
static integer c_n20 = -20;
static integer c_n22 = -22;

/* Subroutine */ int loadrad_(islice)
integer *islice;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg();

    /* Local variables */
    static integer irec, ierr;
    extern /* Subroutine */ int last_();
    extern integer readfield_();
    extern /* Subroutine */ int gauss_hermite__();
    static integer ix;
    static doublereal pradin;

/*     ========================================================= */
/*     fills the array crfield with initial field */
/*     for start up from noise this should be small */
/*     --------------------------------------------------------- */





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


/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */






/*     ------------------------------------------------------------------ */
/*     input/output control */




/*     ------------------------------------------------------------------ */
/*     diagnostic */


/*     -------------------------------------------------------------------- */
/*     time dependency of bunch + slippage field */





/*     simulation control and normalisation parameter */



    diagcom_1.pradol = inputcom_1.prad0;

/*     load initial field */
/*     ------------------------------------------------------------------ */

/* first halfstep no gain (see diagno) */
    if (iocom_1.nfin <= 0) {
	gauss_hermite__(cartcom_1.crfield, &inputcom_1.prad0, &
		inputcom_1.zrayl, &inputcom_1.zwaist, &cartcom_1.xks, &
		cartcom_1.radphase);
/* loa */
    } else {
	irec = tbunchcom_1.nslp - 1 + *islice;
/* get record number */
	if (inputcom_1.alignradf != 0) {
	    irec = inputcom_1.offsetradf + *islice;
	}
/* add offset, when se */
	if (inputcom_1.itdp == 0) {
	    irec = 1;
	}
/* scan+ss -> use 1sr record */
	if (irec > 0) {
/* physical record? */
	    ierr = readfield_(cartcom_1.crfield, &irec);
/* get field from file */
	    if (ierr < 0) {
		last_();
	    }
/* stop if error occured */
	} else {
	    gauss_hermite__(cartcom_1.crfield, &inputcom_1.prad0, &
		    inputcom_1.zrayl, &inputcom_1.zwaist, &cartcom_1.xks, &
		    cartcom_1.radphase);
	}
	pradin = (float)0.;
	i__1 = inputcom_1.ncar * inputcom_1.ncar;
	for (ix = 1; ix <= i__1; ++ix) {
/* copy to crfield */
	    i__2 = ix - 1;
	    d_cnjg(&z__2, &cartcom_1.crfield[ix - 1]);
	    z__1.r = cartcom_1.crfield[i__2].r * z__2.r - cartcom_1.crfield[
		    i__2].i * z__2.i, z__1.i = cartcom_1.crfield[i__2].r * 
		    z__2.i + cartcom_1.crfield[i__2].i * z__2.r;
	    pradin += z__1.r;
	}
/* Computing 2nd power */
	d__1 = cartcom_1.dxy * 510999.06 * simcom_1.xkper0 / cartcom_1.xks;
	inputcom_1.prad0 = pradin * (d__1 * d__1) / 376.73;
	diagcom_1.pradol = inputcom_1.prad0;
    }
    return 0;
} /* loadrad_ */



/* of loadrad */
/* Subroutine */ int gauss_hermite__(cfld, power, zr, zw, rks, phase)
doublecomplex *cfld;
doublereal *power, *zr, *zw, *rks, *phase;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    void z_div();
    double sqrt();
    void z_exp();

    /* Local variables */
    extern /* Subroutine */ int last_();
    extern integer printerr_();
    static doublereal zscal;
    static integer ix, iy;
    static doublecomplex cgauss;
    static integer idx;
    static doublereal xcr, ycr, rcr2;

/*     ======================================================= */
/*     fills array cfld with the fundamental gauss-hermite mode */
/*     using the total power, wavenumber rks, rayleigh length zr */
/*     and waist position zw. */
/*     -------------------------------------------------------- */





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



/*     check for unphysical parameters */

    /* Parameter adjustments */
    --cfld;

    /* Function Body */
    idx = 0;
    if (*zr <= 0.) {
	idx = printerr_(&c_n23, "zrayl in gauss_hermite", (ftnlen)22);
    }
    if (*rks <= 0.) {
	idx = printerr_(&c_n23, "xks in gauss_hermite", (ftnlen)20);
    }
    if (*power < 0.) {
	idx = printerr_(&c_n23, "prad0 in gauss_hermite", (ftnlen)22);
    }
    if (idx < 0) {
	last_();
    }

    d__1 = *rks * .5;
    z__2.r = d__1, z__2.i = 0.;
    d__2 = -(*zw);
    z__3.r = *zr, z__3.i = d__2;
    z_div(&z__1, &z__2, &z__3);
    cgauss.r = z__1.r, cgauss.i = z__1.i;
/* see siegman */
/* Computing 2nd power */
    d__1 = simcom_1.xkper0;
    zscal = sqrt(*power * 753.46000000000004 / 3.14159265358979 * cgauss.r) * 
	    *rks / (d__1 * d__1) / 510999.06;
/* ? */
/* L1: */
    i__1 = inputcom_1.ncar;
    for (iy = 1; iy <= i__1; ++iy) {
	i__2 = inputcom_1.ncar;
	for (ix = 1; ix <= i__2; ++ix) {
	    idx = (iy - 1) * inputcom_1.ncar + ix;
	    xcr = cartcom_1.dxy * (real) (ix - 1) / simcom_1.xkw0 - 
		    inputcom_1.dgrid;
	    ycr = cartcom_1.dxy * (real) (iy - 1) / simcom_1.xkw0 - 
		    inputcom_1.dgrid;
	    rcr2 = xcr * xcr + ycr * ycr;
	    i__3 = idx;
	    z__5.r = -cgauss.r, z__5.i = -cgauss.i;
	    z__4.r = rcr2 * z__5.r, z__4.i = rcr2 * z__5.i;
	    z__6.r = *phase * 0., z__6.i = *phase * 1.;
	    z__3.r = z__4.r + z__6.r, z__3.i = z__4.i + z__6.i;
	    z_exp(&z__2, &z__3);
	    z__1.r = zscal * z__2.r, z__1.i = zscal * z__2.i;
	    cfld[i__3].r = z__1.r, cfld[i__3].i = z__1.i;
/* gauss */
	}
/* ix */
    }
/* iy */
    return 0;
} /* gauss_hermite__ */


/* Subroutine */ int loadslpfld_(nslp)
integer *nslp;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer irec, islp, ierr;
    extern /* Subroutine */ int last_();
    extern integer printerr_(), readfield_();
    extern /* Subroutine */ int dotimerad_(), gauss_hermite__(), pushtimerec_(
	    );

/*     ========================================================= */
/*     fills the array crtime with a seeding field for the first */
/*     slice. */
/*     --------------------------------------------------------- */





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


/*     -------------------------------------------------------------------- */
/*     cartesian mesh */



/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */






/*     ------------------------------------------------------------------ */
/*     input/output control */





/*     file:   timerec.com */

/*     time dependent common block - huge!!!! */
/*     ------------------------------------------------------------------ */




/*     -------------------------------------------------------------------- */
/*     slippage field */

/* must be > than ncar*ncar*nslp */
/* <>0 -> slippage is stored on disk */




    if (inputcom_1.itdp == 0) {
	return 0;
    }

/*     check for limitation in timerecord */
/*     --------------------------------------------------------------- */

    if (inputcom_1.nslice < *nslp * (1 - inputcom_1.iotail)) {
	ierr = printerr_(&c_n20, "no output - ntail too small", (ftnlen)27);
    }
    if (TRUE_) {
	if (*nslp * inputcom_1.ncar * inputcom_1.ncar > 15960700) {
	    ierr = printerr_(&c_n22, " ", (ftnlen)1);
	    last_();
	}
    }

/*     seding of the random number generator */
/*     ---------------------------------------------------------------- */
    i__1 = *nslp - 1;
    for (islp = 1; islp <= i__1; ++islp) {
	if (iocom_1.nfin > 0) {
	    irec = islp;
	    if (inputcom_1.alignradf != 0) {
		irec = irec - *nslp + 1 + inputcom_1.offsetradf;
	    }
	    if (irec > 0) {
		ierr = readfield_(workspace_1.crwork3, &irec);
/* get field from file (record */
		if (ierr < 0) {
		    last_();
		}
/* record nslp is loaded with */
	    } else {
		gauss_hermite__(workspace_1.crwork3, &inputcom_1.prad0, &
			inputcom_1.zrayl, &inputcom_1.zwaist, &cartcom_1.xks, 
			&cartcom_1.radphase);
	    }
	} else {
	    i__2 = 1 - islp;
	    dotimerad_(&i__2);
/* get time-dependence of slippage */
	    gauss_hermite__(workspace_1.crwork3, &inputcom_1.prad0, &
		    inputcom_1.zrayl, &inputcom_1.zwaist, &cartcom_1.xks, &
		    cartcom_1.radphase);
	}
	i__2 = *nslp - islp;
	pushtimerec_(workspace_1.crwork3, &inputcom_1.ncar, &i__2);
    }
    return 0;
} /* loadslpfld_ */



/* Subroutine */ int swapfield_(islp)
integer *islp;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer it;
    extern /* Subroutine */ int pulltimerec_(), pushtimerec_();

/*     ======================================== */
/*     swap current field with then time-record */
/*     ---------------------------------------- */





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


/*     -------------------------------------------------------------------- */
/*     cartesian mesh */




/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







    i__1 = inputcom_1.ncar * inputcom_1.ncar;
    for (it = 1; it <= i__1; ++it) {
	i__2 = it - 1;
	i__3 = it - 1;
	workspace_1.crwork3[i__2].r = cartcom_1.crfield[i__3].r, 
		workspace_1.crwork3[i__2].i = cartcom_1.crfield[i__3].i;
    }
    pulltimerec_(cartcom_1.crfield, &inputcom_1.ncar, islp);
    pushtimerec_(workspace_1.crwork3, &inputcom_1.ncar, islp);
    return 0;
} /* swapfield_ */

