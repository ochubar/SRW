/* partsim.f -- translated by f2c (version 20000118).
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
    //doublecomplex crwork3[1907388], cpart1[7000007], cpart2[1000001], cpart3[1000001];
    //doublereal k2gg[1000001], k2pp[1000001], k3gg[1000001], k3pp[1000001], p1[1000001], p2[1000001];
    //integer iwork[1000001];
    doublecomplex *crwork3, *cpart1, *cpart2, *cpart3;
    doublereal *k2gg, *k2pp, *k3gg, *k3pp, *p1, *p2;
    integer *iwork;
} workspace_;

#define workspace_1 workspace_

Extern struct {
    //doublecomplex crfield[476847], crsource[68121], crmatc[1827], cstep[7], crhm[68121], cwet[1827], cbet[1827];
    doublecomplex *crfield, *crsource, *crmatc, cstep[7], *crhm, *cwet, *cbet;
    doublereal dxy, xks, radphase, besselcoupling[7];
    integer nhloop, hloop[7];
} cartcom_;

#define cartcom_1 cartcom_

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
    doublereal xkw0, xkper0, sval, svalout, gamma0_in__;
    integer npart0, inorun;
} simcom_;

#define simcom_1 simcom_

/* Table of constant values */

static integer c_n1 = -1;

/* Subroutine */ int partsim_(tgam, tthet, dgam, dthet, istepz)
doublereal *tgam, *tthet, *dgam, *dthet;
integer *istepz;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, ddd_oc; //OC port
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sin(), cos(), sqrt(), d_imag();

    /* Local variables */
    static doublecomplex ctmp;
    static integer nharmpart;
    static doublereal btper0;
    static integer ih;
    static doublereal ztemp1, ztemp2;
    static integer ip;

/*     ================================================================== */
/*     define the system of ode of the canonic variables */
/*     t... are the values of the variables at current position */
/*     d... are the value of the differentioan equation */
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





/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







/*     diff-eq. for longitudinal variables: */
/*        - energy (gamma) */
/*        - phase (theta) */
/*     ------------------------------------------------------------------ */


    /* Parameter adjustments */
    --dthet;
    --dgam;
    --tthet;
    --tgam;

    /* Function Body */
    ztemp1 = inputcom_1.xlamds * -2. / inputcom_1.xlamd;
    ztemp2 = inputcom_1.xlamd / inputcom_1.xlamds;

    nharmpart = 1;
/* default - only fundamentalacts on electron */
    if (inputcom_1.iharmsc != 0) {
	nharmpart = inputcom_1.nharm;
/* self-consistent equation of motion */
    }

/*      nharmpart=1 */

    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {

	ctmp.r = 0., ctmp.i = 0.;
	i__2 = nharmpart;
	for (ih = 1; ih <= i__2; ++ih) {
/* loop over harmonics if self-consitent metho */
	    i__3 = ip + (ih - 1) * inputcom_1.npart - 1;
	    d__1 = cos((doublereal) ih * tthet[ip]);
	    d__2 = -sin((doublereal) ih * tthet[ip]);
	    z__3.r = d__1, z__3.i = d__2;
	    z__2.r = workspace_1.cpart1[i__3].r * z__3.r - workspace_1.cpart1[
		    i__3].i * z__3.i, z__2.i = workspace_1.cpart1[i__3].r * 
		    z__3.i + workspace_1.cpart1[i__3].i * z__3.r;
	    z__1.r = ctmp.r + z__2.r, z__1.i = ctmp.i + z__2.i;
	    ctmp.r = z__1.r, ctmp.i = z__1.i;
	}

/* ih */
	btper0 = beamcom_1.btper[ip - 1] + ztemp1 * ctmp.r;
/* perpendicular velocity */
/* Computing 2nd power */
	d__1 = tgam[ip];

	//OC port
	//beamcom_1.btpar[ip - 1] = sqrt(1. - btper0 / (d__1 * d__1));
	ddd_oc = 1. - btper0 / (d__1 * d__1);
	if(ddd_oc <= 0) ddd_oc = 0.000000001;
	beamcom_1.btpar[ip - 1] = sqrt(ddd_oc);
	//OC port end

/* parallel velocity */
	dthet[ip] = dthet[ip] + ztemp2 * ((float)1. - (float)1. / 
		beamcom_1.btpar[ip - 1]) + (float)1.;
/* dtheta/dz */
	dgam[ip] = dgam[ip] + d_imag(&ctmp) / beamcom_1.btpar[ip - 1] / tgam[
		ip] - beamcom_1.ez[ip - 1];
/* dgam */
    }
/* ip */
    return 0;
} /* partsim_ */


/* partsim */
/* Subroutine */ int partsorc_(istepz)
integer *istepz;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    double sqrt();
    void d_cnjg();

    /* Local variables */
    static doublereal rtmp;
    extern /* Subroutine */ int rpos_();
    static integer nharmpart, ih, ip;
    static doublecomplex clocal;
    static doublereal xi;
    static integer ioffset;
    extern doublereal faw2_();
    static doublereal aw2;
    static integer idx1, idx2, idx3, idx4;
    static doublereal wei1, wei2, wei3, wei4;

/*     ================================================================== */
/*     calculates source term for gamma-theta integration */
/*     when higher harmonic coupling is considered */
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


/*     ------------------------------------------------------------------ */
/*     wiggler parameters */




/*     ------------------------------------------------------------------ */
/*     electron beam */




/*     simulation control and normalisation parameter */



/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







/*     call space charge routine to calculate the local field Ez */


    nharmpart = 1;
/* default - only fundamentalacts on electron */
    if (inputcom_1.iharmsc != 0) {
	nharmpart = inputcom_1.nharm;
/* self-consistent equation of motion */
    }

    rpos_(istepz, beamcom_1.xpart, beamcom_1.ypart);

/* position o */
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	aw2 = faw2_(istepz, &beamcom_1.xporb[ip - 1], &beamcom_1.yporb[ip - 1]
		);
/* square of wiggler fiel */
	aw2 += wigcom_1.awdz[*istepz - 1] * wigcom_1.awdz[*istepz - 1];
/* artificial delay */
	beamcom_1.btper[ip - 1] = aw2 + beamcom_1.px[ip - 1] * beamcom_1.px[
		ip - 1] + beamcom_1.py[ip - 1] * beamcom_1.py[ip - 1] + (
		float)1.;
	rtmp = sqrt(aw2) * inputcom_1.xlamds / inputcom_1.xlamd;

/*       interpolation to the grid (index and weight of the 4 surrounding grid points) */

/* effective K-paramete */
	wei1 = beamcom_1.wx[ip - 1] * beamcom_1.wy[ip - 1];
	idx1 = beamcom_1.ipos[(ip << 2) - 4];
	wei2 = beamcom_1.wx[ip - 1] * (1. - beamcom_1.wy[ip - 1]);
	idx2 = beamcom_1.ipos[(ip << 2) - 3];
	wei3 = (1. - beamcom_1.wx[ip - 1]) * beamcom_1.wy[ip - 1];
	idx3 = beamcom_1.ipos[(ip << 2) - 2];
	wei4 = (1. - beamcom_1.wx[ip - 1]) * (1. - beamcom_1.wy[ip - 1]);
	idx4 = beamcom_1.ipos[(ip << 2) - 1];
	i__2 = nharmpart;
	for (ih = 1; ih <= i__2; ih += 2) {
/* sum over odd */
	    ioffset = (ih - 1) * inputcom_1.ncar * inputcom_1.ncar;
	    i__3 = idx1 + ioffset - 1;
	    z__1.r = wei1 * cartcom_1.crfield[i__3].r, z__1.i = wei1 * 
		    cartcom_1.crfield[i__3].i;
	    clocal.r = z__1.r, clocal.i = z__1.i;
/* note Genesis */
	    i__3 = idx2 + ioffset - 1;
	    z__2.r = wei2 * cartcom_1.crfield[i__3].r, z__2.i = wei2 * 
		    cartcom_1.crfield[i__3].i;
	    z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	    clocal.r = z__1.r, clocal.i = z__1.i;
/* u = k*n/ k_u */
	    i__3 = idx3 + ioffset - 1;
	    z__2.r = wei3 * cartcom_1.crfield[i__3].r, z__2.i = wei3 * 
		    cartcom_1.crfield[i__3].i;
	    z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	    clocal.r = z__1.r, clocal.i = z__1.i;
	    i__3 = idx4 + ioffset - 1;
	    z__2.r = wei4 * cartcom_1.crfield[i__3].r, z__2.i = wei4 * 
		    cartcom_1.crfield[i__3].i;
	    z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	    clocal.r = z__1.r, clocal.i = z__1.i;
	    i__3 = ip + (ih - 1) * inputcom_1.npart - 1;
	    d_cnjg(&z__4, &clocal);
	    z__3.r = rtmp * z__4.r, z__3.i = rtmp * z__4.i;
	    i__4 = ih - 1;
	    z__2.r = cartcom_1.besselcoupling[i__4] * z__3.r, z__2.i = 
		    cartcom_1.besselcoupling[i__4] * z__3.i;
	    d__1 = (doublereal) ih;
	    z__1.r = z__2.r / d__1, z__1.i = z__2.i / d__1;
	    workspace_1.cpart1[i__3].r = z__1.r, workspace_1.cpart1[i__3].i = 
		    z__1.i;
	}
	i__2 = nharmpart;
	for (ih = 2; ih <= i__2; ih += 2) {
/* sum over even */
	    xi = (doublereal) ih * inputcom_1.xlamd / inputcom_1.xlamds / 
		    beamcom_1.gamma[ip - 1] / beamcom_1.gamma[ip - 1] * sqrt(
		    aw2);
	    ioffset = (ih - 1) * inputcom_1.ncar * inputcom_1.ncar;
	    i__3 = idx1 + ioffset - 1;
	    z__1.r = wei1 * cartcom_1.crfield[i__3].r, z__1.i = wei1 * 
		    cartcom_1.crfield[i__3].i;
	    clocal.r = z__1.r, clocal.i = z__1.i;
/* note Genesis */
	    i__3 = idx2 + ioffset - 1;
	    z__2.r = wei2 * cartcom_1.crfield[i__3].r, z__2.i = wei2 * 
		    cartcom_1.crfield[i__3].i;
	    z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	    clocal.r = z__1.r, clocal.i = z__1.i;
/* u = k*n/ k_u */
	    i__3 = idx3 + ioffset - 1;
	    z__2.r = wei3 * cartcom_1.crfield[i__3].r, z__2.i = wei3 * 
		    cartcom_1.crfield[i__3].i;
	    z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	    clocal.r = z__1.r, clocal.i = z__1.i;
	    i__3 = idx4 + ioffset - 1;
	    z__2.r = wei4 * cartcom_1.crfield[i__3].r, z__2.i = wei4 * 
		    cartcom_1.crfield[i__3].i;
	    z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	    clocal.r = z__1.r, clocal.i = z__1.i;
	    i__3 = ip + (ih - 1) * inputcom_1.npart - 1;
	    d_cnjg(&z__6, &clocal);
	    z__5.r = rtmp * z__6.r, z__5.i = rtmp * z__6.i;
	    i__4 = ih - 1;
	    z__4.r = cartcom_1.besselcoupling[i__4] * z__5.r, z__4.i = 
		    cartcom_1.besselcoupling[i__4] * z__5.i;
	    d__1 = (doublereal) ih;
	    z__3.r = z__4.r / d__1, z__3.i = z__4.i / d__1;
	    z__2.r = xi * z__3.r, z__2.i = xi * z__3.i;
	    d__2 = sqrt((float)2.) * beamcom_1.px[ip - 1];
	    z__7.r = 0., z__7.i = d__2;
	    z__1.r = z__2.r * z__7.r - z__2.i * z__7.i, z__1.i = z__2.r * 
		    z__7.i + z__2.i * z__7.r;
	    workspace_1.cpart1[i__3].r = z__1.r, workspace_1.cpart1[i__3].i = 
		    z__1.i;
/* missing px he */
	}
    }

    return 0;
} /* partsorc_ */



/* partsim */
/* Subroutine */ int harmcoupling_(awloc)
doublereal *awloc;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer pow_ii();

    /* Local variables */
    extern doublereal bessj_();
    static integer ih;
    static doublereal xi;

/*     ============================================================ */
/*     routine to calculate the coupling to higher modes */
/*     ------------------------------------------------------ */





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



    if (*awloc < 1e-25) {
/* drift -> set coupling to zero */
	i__1 = inputcom_1.nharm;
	for (ih = 1; ih <= i__1; ++ih) {
	    cartcom_1.besselcoupling[ih - 1] = (float)0.;
	}
	return 0;
    }

/* Computing 2nd power */
    d__1 = *awloc;
/* Computing 2nd power */
    d__2 = *awloc;
    xi = d__1 * d__1 / (d__2 * d__2 + 1.) / 2;

    if (inputcom_1.iwityp != 0) {
/* helical undulator */
	cartcom_1.besselcoupling[0] = (float)1.;
	i__1 = inputcom_1.nharm;
	for (ih = 2; ih <= i__1; ++ih) {
/* even harmonic is disabled due to t */
	    cartcom_1.besselcoupling[ih - 1] = (float)0.;
/* break in symmetry */
	}
    } else {
	i__1 = inputcom_1.nharm;
	for (ih = 1; ih <= i__1; ++ih) {
	    if (ih % 2 == 1) {
		i__2 = (ih - 1) / 2;
		d__1 = xi * (doublereal) ih;
		i__3 = (ih + 1) / 2;
		d__2 = xi * (doublereal) ih;
		i__4 = (ih - 1) / 2;
		cartcom_1.besselcoupling[ih - 1] = (bessj_(&i__2, &d__1) - 
			bessj_(&i__3, &d__2)) * pow_ii(&c_n1, &i__4);
	    } else {
		i__2 = (ih - 2) / 2;
		d__1 = xi * (doublereal) ih;
		i__3 = (ih + 2) / 2;
		d__2 = xi * (doublereal) ih;
		i__4 = (ih - 2) / 2;
		cartcom_1.besselcoupling[ih - 1] = (bessj_(&i__2, &d__1) - 
			bessj_(&i__3, &d__2)) * (float).5 * pow_ii(&c_n1, &
			i__4);
	    }
/* mulitplied with the transverse momentum */
/* note that for full coupling it has to */
	}
    }
    return 0;
} /* harmcoupling_ */

