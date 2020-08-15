/* tdepend.f -- translated by f2c (version 20000118).
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
    //doublereal xpart[1000001], ypart[1000001], px[1000001], py[1000001], 
	   // gamma[1000001], theta[1000001], xporb[1000001], yporb[1000001], 
	   // btpar[1000001], btper[1000001], ez[1000001], wx[1000001], wy[
	   // 1000001], xcuren, dedz, tdmin, tdmax, delcharge, dtd, charge;
    //integer lostid[1000001], lost, losttot, ipos[4000004]	/* was [4][1000001] */;
    doublereal *xpart, *ypart, *px, *py, *gamma, 
		*theta, *xporb, *yporb, *btpar, *btper, *ez, *wx, *wy, 
		xcuren, dedz, tdmin, tdmax, delcharge, dtd, charge;
    integer *lostid, lost, losttot, *ipos;
} beamcom_;

#define beamcom_1 beamcom_

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
    //doublecomplex crfield[476847], crsource[68121], crmatc[1827], cstep[7], crhm[68121], cwet[1827], cbet[1827];
    doublecomplex *crfield, *crsource, *crmatc, cstep[7], *crhm, *cwet, *cbet;
    doublereal dxy, xks, radphase, besselcoupling[7];
    integer nhloop, hloop[7];
} cartcom_;

#define cartcom_1 cartcom_

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c_n13 = -13;

/* Subroutine */ int dotime_(islice)
integer *islice;
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double exp();
    integer s_wsli(), do_lio(), e_wsli();

    /* Local variables */
    static doublereal zpos;
    extern integer printerr_();
    static integer i__;
    static char cdiff[30];
    static doublereal w1, w2;
    extern /* Subroutine */ int dotimerad_();
    extern integer luf_();
    static integer idx;
    static doublereal invcur;
    extern doublereal ran1_();

    /* Fortran I/O blocks */
    static icilist io___6 = { 0, cdiff, 0, 0, 30, 1 };
    static icilist io___8 = { 0, cdiff, 0, 0, 30, 1 };


/*     =================================================================== */
/*     set the beam parameter for the case of time-dependence */
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


/*     ------------------------------------------------------------------ */
/*     electron beam */




/*     -------------------------------------------------------------------- */
/*     time dependency of bunch + slippage field */





/*     simulation control and normalisation parameter */


    if (inputcom_1.itdp == 0) {
	return 0;
    }

    if (*islice == 1) {
	simcom_1.npart0 = inputcom_1.npart;
/* save #particles */
	w1 = ran1_(&inputcom_1.ipseed);
/* init ran1 */
    }

    inputcom_1.npart = simcom_1.npart0;

/* compensate former particle losses */
    if (tbunchcom_1.ndata <= 1) {
/* internal generation of time-dependence */
	if (inputcom_1.curlen <= (float)0.) {
	    invcur = (float)0.;
	} else {
	    invcur = (float)1. / inputcom_1.curlen;
	}
	zpos = (doublereal) (inputcom_1.ntail + *islice - 1) * 
		inputcom_1.zsep * inputcom_1.xlamds * invcur;
/* normalized t-pos */
	beamcom_1.xcuren = inputcom_1.curpeak * exp(zpos * -.5 * zpos);
/* beam current */
	dotimerad_(islice);
	return 0;
    }

    zpos = (inputcom_1.ntail + *islice - 1) * inputcom_1.zsep * 
	    inputcom_1.xlamds;

/* position in m */
    idx = luf_(&zpos, tbunchcom_1.tpos, &tbunchcom_1.ndata) - 1;
/* find position in array */
    if (idx <= 0) {
	idx = 1;
    }
    if (idx >= tbunchcom_1.ndata) {
	idx = tbunchcom_1.ndata - 1;
    }

    if (zpos < tbunchcom_1.tpos[0]) {
	s_wsli(&io___6);
	d__1 = zpos - tbunchcom_1.tpos[0];
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsli();
	i__ = printerr_(&c_n13, cdiff, (ftnlen)30);
    }
    if (zpos > tbunchcom_1.tpos[tbunchcom_1.ndata - 1]) {
	s_wsli(&io___8);
	d__1 = zpos - tbunchcom_1.tpos[tbunchcom_1.ndata - 1];
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsli();
	i__ = printerr_(&c_n13, cdiff, (ftnlen)30);
    }

    w2 = (zpos - tbunchcom_1.tpos[idx - 1]) / (tbunchcom_1.tpos[idx] - 
	    tbunchcom_1.tpos[idx - 1]);
/* weight of higher inde */
    w1 = 1. - w2;
/* weight of lower index */
    inputcom_1.gamma0 = w1 * tbunchcom_1.tgam0[idx - 1] + w2 * 
	    tbunchcom_1.tgam0[idx];
/* interpolation */
    inputcom_1.delgam = w1 * tbunchcom_1.tdgam[idx - 1] + w2 * 
	    tbunchcom_1.tdgam[idx];
/* temporary stored */
    inputcom_1.rxbeam = w1 * tbunchcom_1.txrms[idx - 1] + w2 * 
	    tbunchcom_1.txrms[idx];
/* into working arays */
    inputcom_1.rybeam = w1 * tbunchcom_1.tyrms[idx - 1] + w2 * 
	    tbunchcom_1.tyrms[idx];
    inputcom_1.xbeam = w1 * tbunchcom_1.txpos[idx - 1] + w2 * 
	    tbunchcom_1.txpos[idx];
    inputcom_1.ybeam = w1 * tbunchcom_1.typos[idx - 1] + w2 * 
	    tbunchcom_1.typos[idx];
    inputcom_1.emitx = w1 * tbunchcom_1.temitx[idx - 1] + w2 * 
	    tbunchcom_1.temitx[idx];
    inputcom_1.emity = w1 * tbunchcom_1.temity[idx - 1] + w2 * 
	    tbunchcom_1.temity[idx];
    inputcom_1.pxbeam = w1 * tbunchcom_1.tpxpos[idx - 1] + w2 * 
	    tbunchcom_1.tpxpos[idx];
    inputcom_1.pybeam = w1 * tbunchcom_1.tpypos[idx - 1] + w2 * 
	    tbunchcom_1.tpypos[idx];
    inputcom_1.alphax = w1 * tbunchcom_1.talphx[idx - 1] + w2 * 
	    tbunchcom_1.talphx[idx];
    inputcom_1.alphay = w1 * tbunchcom_1.talphy[idx - 1] + w2 * 
	    tbunchcom_1.talphy[idx];
    beamcom_1.xcuren = w1 * tbunchcom_1.tcurrent[idx - 1] + w2 * 
	    tbunchcom_1.tcurrent[idx];
    beamcom_1.dedz = (w1 * tbunchcom_1.tloss[idx - 1] + w2 * 
	    tbunchcom_1.tloss[idx]) * inputcom_1.delz * inputcom_1.xlamd / 
	    510999.06;
    dotimerad_(islice);
    return 0;
} /* dotime_ */

/* of dotime */
/* Subroutine */ int dotimerad_(islice)
integer *islice;
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    integer s_wsli(), do_lio(), e_wsli();

    /* Local variables */
    static doublereal zpos;
    extern integer printerr_();
    static integer i__;
    static char cdiff[30];
    static doublereal w1, w2;
    extern integer luf_();
    static integer idx;

    /* Fortran I/O blocks */
    static icilist io___13 = { 0, cdiff, 0, 0, 30, 1 };
    static icilist io___15 = { 0, cdiff, 0, 0, 30, 1 };


/*     =================================================================== */
/*     set the beam parameter for the case of time-dependence */
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



/*     -------------------------------------------------------------------- */
/*     time dependency of bunch + slippage field */





/*     simulation control and normalisation parameter */


    if (tbunchcom_1.nraddata <= 1) {
	return 0;
    }
    zpos = (inputcom_1.ntail + *islice - 1) * inputcom_1.zsep * 
	    inputcom_1.xlamds;
/* position in m */
    idx = luf_(&zpos, tbunchcom_1.tradpos, &tbunchcom_1.nraddata) - 1;
/* find position in array */
    if (idx <= 0) {
	s_wsli(&io___13);
	d__1 = zpos - tbunchcom_1.tradpos[0];
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsli();
	i__ = printerr_(&c_n13, cdiff, (ftnlen)30);
	idx = 1;
    }
    if (idx >= tbunchcom_1.nraddata) {
	s_wsli(&io___15);
	d__1 = zpos - tbunchcom_1.tradpos[tbunchcom_1.ndata - 1];
	do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsli();
	i__ = printerr_(&c_n13, cdiff, (ftnlen)30);
	idx = tbunchcom_1.nraddata - 1;
    }
    w2 = (zpos - tbunchcom_1.tradpos[idx - 1]) / (tbunchcom_1.tradpos[idx] - 
	    tbunchcom_1.tradpos[idx - 1]);
/* weight of hi */
    w1 = 1. - w2;
/* weight of lower index */
    inputcom_1.prad0 = w1 * tbunchcom_1.tprad0[idx - 1] + w2 * 
	    tbunchcom_1.tprad0[idx];
/* interpolation */
    inputcom_1.zrayl = w1 * tbunchcom_1.tzrayl[idx - 1] + w2 * 
	    tbunchcom_1.tzrayl[idx];
/* temporary stored */
    inputcom_1.zwaist = w1 * tbunchcom_1.tzwaist[idx - 1] + w2 * 
	    tbunchcom_1.tzwaist[idx];
/* into working arays */
    cartcom_1.radphase = w1 * tbunchcom_1.tradphase[idx - 1] + w2 * 
	    tbunchcom_1.tradphase[idx];
    if (inputcom_1.prad0 < 0.) {
	inputcom_1.prad0 = 0.;
    }
    return 0;
} /* dotimerad_ */

