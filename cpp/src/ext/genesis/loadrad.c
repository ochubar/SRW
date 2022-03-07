/* loadrad.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

Extern struct {
    //doublecomplex crfield[476847], crsource[68121], crmatc[1827], cstep[7], crhm[68121], cwet[1827], cbet[1827];
    doublecomplex *crfield, *crsource, *crmatc, cstep[7], *crhm, *cwet, *cbet;
    doublereal dxy, xks, radphase, besselcoupling[7];
    integer nhloop, hloop[7];
} cartcom_;

#define cartcom_1 cartcom_

Extern struct {
    doublereal aw0, xkx, xky, wcoefz[3], xlamd, fbess0, delaw, awd, awx, awy, 
	    gamma0, delgam, rxbeam, rybeam, alphax, alphay, emitx, emity, 
	    xbeam, ybeam, pxbeam, pybeam, cuttail, curpeak, conditx, condity, 
	    bunch, bunchphase, emod, emodphase, xlamds, prad0, zrayl, zwaist, 
	    rmax0, zsep, delz, zstop, quadf, quadd, fl, dl, drl, f1st, qfdx, 
	    qfdy, sl, solen, curlen, shotnoise, svar, dgrid, eloss, version, 
	    ibfield, imagl, idril, igamref, pradh0, itram11, itram12, itram13,
	     itram14, itram15, itram16, itram21, itram22, itram23, itram24, 
	    itram25, itram26, itram31, itram32, itram33, itram34, itram35, 
	    itram36, itram41, itram42, itram43, itram44, itram45, itram46, 
	    itram51, itram52, itram53, itram54, itram55, itram56, itram61, 
	    itram62, itram63, itram64, itram65, itram66, rmax0sc;
    integer iseed, nwig, nsec, npart, ncar, lbc, nscr, nscz, nptr, ildgam, 
	    ildpsi, ildx, ildy, ildpx, ildpy, itgaus, nbins, iphsty, ishsty, 
	    ippart, ispart, ipradi, isradi, iertyp, iwityp, idump, iotail, 
	    nharm, iallharm, iharmsc, magin, magout, lout[40], ffspec, ntail, 
	    nslice, iall, itdp, ipseed, iscan, nscan, isntyp, isravg, isrsig, 
	    iorb, ndcut, idmppar, idmpfld, ilog, igamgaus, convharm, 
	    alignradf, offsetradf, multconv, trama, iscrkup, inverfc;
    char beamfile[30], fieldfile[30], maginfile[30], magoutfile[30], 
	    outputfile[30], inputfile[30], scan[30], distfile[30], partfile[
	    30], filetype[30], radfile[30];
} inputcom_;

#define inputcom_1 inputcom_

Extern struct {
    //doublecomplex crwork3[1907388], cpart1[7000007], cpart2[1000001], cpart3[1000001];
    //doublereal k2gg[1000001], k2pp[1000001], k3gg[1000001], k3pp[1000001], p1[1000001], p2[1000001];
    //integer iwork[1000001];
    doublecomplex *crwork3, *cpart1, *cpart2, *cpart3;
    doublereal *k2gg, *k2pp, *k3gg, *k3pp, *p1, *p2;
    integer *iwork;
} workspace_;

#define workspace_1 workspace_

Extern struct {
    doublereal distversion, distrev;
    integer iout[39], nout, nfld, npar, ndump, firstout, nfin, irecpar, 
	    irecfld, kout, ndmp2, npin, nlog, ndis, ncoldis, iconv2t, iconv2g,
	     icolpar[10], ndistsize, iconv2px, iconv2py, nprobe, ftype, 
	    ftdist, ftpart, ftfield, ndumph[6], nfldh[6];
} iocom_;

#define iocom_1 iocom_

Extern struct {
    //doublereal error[10000], gain[10000], phimid[70000]	/* was [7][10000] */, 
	   // whalf[10000], logp[10000], powmid, xrms[10000], yrms[10000], xpos[
	   // 10000], ypos[10000], pmodhist[70000]	/* was [7][10000] */, 
	   // gamhist[10000], diver[10000], pradol, pinit, bunphase[70000]	
	   // /* was [7][10000] */, dgamhist[10000], ffield[70000]	/* 
	   // was [7][10000] */, pradoln[7], pgainhist[70000]	/* was [7][
	   // 10000] */, pmidhist[70000]	/* was [7][10000] */;
    doublereal *error, *gain, *phimid, *whalf, *logp, powmid, *xrms, *yrms, *xpos, 
		*ypos, *pmodhist, *gamhist, *diver, pradol, pinit, 
		*bunphase, *dgamhist, *ffield, pradoln[7], 
		*pgainhist, *pmidhist;
    integer ihist;
} diagcom_;

#define diagcom_1 diagcom_

Extern struct {
    //doublereal tgam0[50000], tdgam[50000], temitx[50000], temity[50000], 
	   // txrms[50000], tyrms[50000], txpos[50000], typos[50000], tpxpos[
	   // 50000], tpypos[50000], talphx[50000], talphy[50000], tcurrent[
	   // 50000], tpos[50000], tloss[50000], distgam[1250000], distx[
	   // 1250000], disty[1250000], distpx[1250000], distpy[1250000], distt[
	   // 1250000], tradpos[50000], tzrayl[50000], tzwaist[50000], tprad0[
	   // 50000], tradphase[50000];
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
    //doublecomplex crtime[79803500];
    doublecomplex *crtime;
    integer ntmp;
} tslipcom_;

#define tslipcom_1 tslipcom_

Extern struct {
    integer mpi_id__, mpi_err__, mpi_size__, mpi_loop__, nfldmpi, nparmpi, nfldhmpi[6];
} mpicom_;

#define mpicom_1 mpicom_

/* Table of constant values */

static integer c__1 = 1;
static integer c_n23 = -23;
static integer c_n20 = -20;
static integer c_n22 = -22;

/* Subroutine */ int loadrad_(islice)
integer *islice;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg();

    /* Local variables */
    static integer irec, ierr;
    extern /* Subroutine */ int last_();
    static integer n;
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
/* maximum of particle in imp */
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
/*     transfermatrix */


/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */






/*     ------------------------------------------------------------------ */
/*     input/output control */




/*     ------------------------------------------------------------------ */
/*     diagnostic */



/*     -------------------------------------------------------------------- */
/*     time dependency of bunch + slippage field */





/*     simulation control and normalisation parameter */



/*     initialize field */
/*     --------------------------- */
    i__1 = inputcom_1.ncar * inputcom_1.ncar * cartcom_1.nhloop;
    for (ix = 1; ix <= i__1; ++ix) {
	i__2 = ix - 1;
	cartcom_1.crfield[i__2].r = 0., cartcom_1.crfield[i__2].i = 0.;
    }

    diagcom_1.pradoln[0] = inputcom_1.prad0;
/* first halfstep no gain (see diagno) */
    for (n = 2; n <= 7; ++n) {
	diagcom_1.pradoln[n - 1] = 0.;
/* kg */
	if (n == inputcom_1.nharm) {
	    diagcom_1.pradoln[n - 1] = inputcom_1.pradh0;
	}
    }

/*     load initial field */
/*     ------------------------------------------------------------------ */


    if (iocom_1.nfin <= 0) {
	gauss_hermite__(cartcom_1.crfield, &inputcom_1.prad0, &
		inputcom_1.zrayl, &inputcom_1.zwaist, &cartcom_1.xks, &
		cartcom_1.radphase, &c__1);
/* load gauss-hermite mode for all harmonic */
	if (inputcom_1.nharm > 1 && inputcom_1.pradh0 > 0.) {
	    d__1 = inputcom_1.zrayl * (doublereal) inputcom_1.nharm;
	    d__2 = cartcom_1.xks * (doublereal) inputcom_1.nharm;
	    gauss_hermite__(cartcom_1.crfield, &inputcom_1.pradh0, &d__1, &
		    inputcom_1.zwaist, &d__2, &cartcom_1.radphase, &
		    inputcom_1.nharm);
/* load gauss-hermite mod */
	}
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
/* get fundamental field from */
	    if (ierr < 0) {
		last_();
	    }
/* stop if error occured */
	} else {
	    d__1 = inputcom_1.zrayl * (doublereal) inputcom_1.nharm;
	    d__2 = cartcom_1.xks * (doublereal) inputcom_1.nharm;
	    gauss_hermite__(cartcom_1.crfield, &inputcom_1.prad0, &d__1, &
		    inputcom_1.zwaist, &d__2, &cartcom_1.radphase, &c__1);
	    if (inputcom_1.nharm > 1 && inputcom_1.pradh0 > 0.) {
		d__1 = inputcom_1.zrayl * (doublereal) inputcom_1.nharm;
		d__2 = cartcom_1.xks * (doublereal) inputcom_1.nharm;
		gauss_hermite__(cartcom_1.crfield, &inputcom_1.pradh0, &d__1, 
			&inputcom_1.zwaist, &d__2, &cartcom_1.radphase, &
			inputcom_1.nharm);
/* load gauss-hermite m */
	    }
	}
	pradin = 0.;
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
	diagcom_1.pradoln[0] = inputcom_1.prad0;
    }
    return 0;
} /* loadrad_ */



/* of loadrad */
/* Subroutine */ int gauss_hermite__(cfld, power, zr, zw, rks, phase, harm)
doublecomplex *cfld;
doublereal *power, *zr, *zw, *rks, *phase;
integer *harm;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    void z_div();
    double sqrt();
    void z_exp(), d_cnjg();

    /* Local variables */
    static integer ioff;
    static doublereal dump;
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

/*     Note - only the fundamental is loaded. harmonics are set to zero */
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
/* maximum of particle in imp */
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
/*     transfermatrix */


/*     ------------------------------------------------------------------ */
/*     diagnostic */




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
	idx = printerr_(&c_n23, "power in gauss_hermite", (ftnlen)22);
    }
    if (idx < 0) {
	last_();
    }

    ioff = inputcom_1.ncar * inputcom_1.ncar * (*harm - 1);
/* offset for harmonics */
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
    dump = 0.;
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
	    i__3 = idx + ioff;
	    z__5.r = -cgauss.r, z__5.i = -cgauss.i;
	    z__4.r = rcr2 * z__5.r, z__4.i = rcr2 * z__5.i;
	    z__6.r = *phase * 0., z__6.i = *phase * 1.;
	    z__3.r = z__4.r + z__6.r, z__3.i = z__4.i + z__6.i;
	    z_exp(&z__2, &z__3);
	    z__1.r = zscal * z__2.r, z__1.i = zscal * z__2.i;
	    cfld[i__3].r = z__1.r, cfld[i__3].i = z__1.i;
/* gaussian */
	    i__3 = idx + ioff;
	    d_cnjg(&z__2, &cfld[idx + ioff]);
	    z__1.r = cfld[i__3].r * z__2.r - cfld[i__3].i * z__2.i, z__1.i = 
		    cfld[i__3].r * z__2.i + cfld[i__3].i * z__2.r;
	    dump += z__1.r;
/* =sum */
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
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer irec, islp, ierr;
    extern integer printerr_();
    extern /* Subroutine */ int last_();
    static integer i__;
    extern integer readfield_();
    extern /* Subroutine */ int gauss_hermite__(), dotimerad_();
    static integer ix;
    extern /* Subroutine */ int pushtimerec_();

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
/* maximum of particle in imp */
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
/*     transfermatrix */


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

/* must be > than ncar*ncar*nslp(98) */
/* <>0 -> slippage is stored on disk */




/*     ------------------------------------------------------------------ */
/*     diagnostic */




    if (inputcom_1.itdp == 0) {
	return 0;
    }

/*     check for limitation in timerecord */
/*     --------------------------------------------------------------- */

    if (inputcom_1.nslice < *nslp * (1 - inputcom_1.iotail)) {
	ierr = printerr_(&c_n20, "no output - ntail too small", (ftnlen)27);
    }
    if (TRUE_) {
	if (*nslp * inputcom_1.ncar * inputcom_1.ncar * cartcom_1.nhloop > 
		79803500) {
	    ierr = printerr_(&c_n22, " ", (ftnlen)1);
	    last_();
	}
    }

/*     load initial slippage field from file or internal generation */
/*     ---------------------------------------------------------------- */

    i__1 = *nslp - 1;
    for (islp = 1; islp <= i__1; ++islp) {

	i__2 = inputcom_1.ncar * inputcom_1.ncar * cartcom_1.nhloop;
	for (ix = 1; ix <= i__2; ++ix) {
	    i__3 = ix - 1;
	    workspace_1.crwork3[i__3].r = 0., workspace_1.crwork3[i__3].i = 
		    0.;
/* initialize the radiation field */
	}

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
		diagcom_1.pradoln[0] = inputcom_1.prad0;
		for (i__ = 2; i__ <= 7; ++i__) {
		    diagcom_1.pradoln[i__ - 1] = (float)0.;
		    if (i__ == inputcom_1.nharm) {
			diagcom_1.pradoln[i__ - 1] = inputcom_1.pradh0;
		    }
		}
		gauss_hermite__(workspace_1.crwork3, &inputcom_1.prad0, &
			inputcom_1.zrayl, &inputcom_1.zwaist, &cartcom_1.xks, 
			&cartcom_1.radphase, &c__1);
		if (inputcom_1.nharm > 1 && inputcom_1.pradh0 > 0.) {
		    d__1 = inputcom_1.zrayl * (doublereal) inputcom_1.nharm;
		    d__2 = cartcom_1.xks * (doublereal) inputcom_1.nharm;
		    gauss_hermite__(workspace_1.crwork3, &inputcom_1.pradh0, &
			    d__1, &inputcom_1.zwaist, &d__2, &
			    cartcom_1.radphase, &inputcom_1.nharm);
/* lo */
		}
	    }
	} else {
	    i__2 = islp + 1 - *nslp;
	    dotimerad_(&i__2);
/* get time-dependence of sli */
	    gauss_hermite__(workspace_1.crwork3, &inputcom_1.prad0, &
		    inputcom_1.zrayl, &inputcom_1.zwaist, &cartcom_1.xks, &
		    cartcom_1.radphase, &c__1);
	    if (inputcom_1.nharm > 1 && inputcom_1.pradh0 > 0.) {
		d__1 = inputcom_1.zrayl * (doublereal) inputcom_1.nharm;
		d__2 = cartcom_1.xks * (doublereal) inputcom_1.nharm;
		gauss_hermite__(workspace_1.crwork3, &inputcom_1.pradh0, &
			d__1, &inputcom_1.zwaist, &d__2, &cartcom_1.radphase, 
			&inputcom_1.nharm);
/* load */
	    }
	    diagcom_1.pradoln[0] = inputcom_1.prad0;
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
    extern /* Subroutine */ int mpi_send__(), mpi_recv__();
    static integer it, status[1], mpi_bot__, mpi_top__, memsize;
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
/* maximum of particle in imp */
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
/*     common block for mpi implementation for genesis */





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
/*     transfermatrix */


/*     -------------------------------------------------------------------- */
/*     cartesian mesh */




/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







    memsize = inputcom_1.ncar * inputcom_1.ncar * cartcom_1.nhloop;
    if (mpicom_1.mpi_loop__ > 1) {

	i__1 = memsize;
	for (it = 1; it <= i__1; ++it) {
	    i__2 = it - 1;
	    i__3 = it - 1;
	    workspace_1.crwork3[i__2].r = cartcom_1.crfield[i__3].r, 
		    workspace_1.crwork3[i__2].i = cartcom_1.crfield[i__3].i;
	}

	mpi_top__ = mpicom_1.mpi_id__ + 1;
	if (mpi_top__ >= mpicom_1.mpi_loop__) {
	    mpi_top__ = 0;
	}
	mpi_bot__ = mpicom_1.mpi_id__ - 1;
	if (mpi_bot__ < 0) {
	    mpi_bot__ = mpicom_1.mpi_loop__ - 1;
	}

	if (mpicom_1.mpi_id__ % 2 == 0) {
	    mpi_send__(workspace_1.crwork3, &memsize, &c__1, &mpi_top__, &
		    mpicom_1.mpi_id__, &c__1, &mpicom_1.mpi_err__);
	    mpi_recv__(cartcom_1.crfield, &memsize, &c__1, &mpi_bot__, &
		    mpi_bot__, &c__1, status, &mpicom_1.mpi_err__);
	} else {
	    mpi_recv__(cartcom_1.crfield, &memsize, &c__1, &mpi_bot__, &
		    mpi_bot__, &c__1, status, &mpicom_1.mpi_err__);
	    mpi_send__(workspace_1.crwork3, &memsize, &c__1, &mpi_top__, &
		    mpicom_1.mpi_id__, &c__1, &mpicom_1.mpi_err__);
	}
    }
    if (mpicom_1.mpi_id__ > 0) {
	return 0;
    }
    i__1 = memsize;
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

