/* main.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

Extern struct {
    integer mpi_id__, mpi_err__, mpi_size__, mpi_loop__, nfldmpi, nparmpi, nfldhmpi[6];
} mpicom_;

#define mpicom_1 mpicom_

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
    //doublecomplex crfield[476847], crsource[68121], crmatc[1827], cstep[7], crhm[68121], cwet[1827], cbet[1827];
    doublecomplex *crfield, *crsource, *crmatc, cstep[7], *crhm, *cwet, *cbet;
    doublereal dxy, xks, radphase, besselcoupling[7];
    integer nhloop, hloop[7];
} cartcom_;

#define cartcom_1 cartcom_

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
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer islp, isep;
    extern /* Subroutine */ int mpi_init__(), loadbeam_(), chk_loss__(), 
	    outdumpslippage_(), last_(), mpi_comm_rank__(), mpi_comm_size__(),
	     openoutputbinmpi_(), stepz_(), swapfield_(), mpi_merge__(), 
	    doscan_();
    static integer lstepz, istepz, islice;
    extern /* Subroutine */ int initio_(), initrun_(), outglob_(), dotime_(), 
	    loadrad_(), output_(), outhist_(), outdump_(), closeoutputbinmpi_(
	    );







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
/*     transfermatrix */


/*     -------------------------------------------------------------------- */
/*     time dependency of bunch + slippage field */




/*     -------------------------------------------------------------------- */
/*     cartesian mesh */



/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







/*     ------------------------------------------------------------------ */
/*     wiggler parameters */





    mpi_init__(&mpicom_1.mpi_err__);
    mpi_comm_rank__(&c__1, &mpicom_1.mpi_id__, &mpicom_1.mpi_err__);
    mpi_comm_size__(&c__1, &mpicom_1.mpi_size__, &mpicom_1.mpi_err__);

    initio_();

/* open files, read in input */
    initrun_();
/* initialize simulation */
    outglob_();

/*     temporary file for debugging */

/*      open(69,file='debug.dat',status='unknown',access='direct', */
/*     +    recl=16*nptr) */

/*     start loop for each slice */

/* output of global parameter (t-independent) */
    i__1 = inputcom_1.nslice;
    i__2 = mpicom_1.mpi_size__;
    for (islice = mpicom_1.mpi_id__ + 1; i__2 < 0 ? islice >= i__1 : islice <=
	     i__1; islice += i__2) {

/* loop for each slice */
	mpicom_1.mpi_loop__ = mpicom_1.mpi_size__;
	if (islice - mpicom_1.mpi_id__ >= inputcom_1.nslice - 
		mpicom_1.mpi_size__ + 1) {
	    mpicom_1.mpi_loop__ = inputcom_1.nslice % mpicom_1.mpi_size__;
	}
	if (mpicom_1.mpi_loop__ == 0) {
	    mpicom_1.mpi_loop__ = mpicom_1.mpi_size__;
	}

/*     initial loading */

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
	openoutputbinmpi_(&islice);

/* open binary outputfile for mpi fe */
	output_(&istepz, &islice, &simcom_1.xkw0);

/*       propagate beam for nu wiggler periods */

	i__3 = tbunchcom_1.nslp;
	for (islp = 1; islp <= i__3; ++islp) {
/* loop over slippage (advance field */
	    lstepz = tbunchcom_1.nsep;
	    if (islp == tbunchcom_1.nslp) {
		lstepz = wigcom_1.nstepz - (islp - 1) * tbunchcom_1.nsep;
	    }

/* correct for l */
	    i__4 = lstepz;
	    for (isep = 1; isep <= i__4; ++isep) {
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
	closeoutputbinmpi_();
/* close binary files for m */
    }

/*     merge has to be done first so that the field dump is already */
/*     created and the records from the slippage field */
/*     are written at the right position */


    mpi_merge__();

/* merge single outputfiles into one */
    outdumpslippage_();

/*     temporary for debugging */

/*      close(69) */

/* write slippage field, which escaped */
    last_();

/* end run */
} /* MAIN__ */

/* Main program alias */ int genesis_ () { MAIN__ (); return 0;} //OC port
