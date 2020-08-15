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
    //doublecomplex crwork3[116964], cpart1[70001], cpart2[70001];
    //doublereal k2gg[70001], k2pp[70001], k3gg[70001], k3pp[70001], p1[70001], 
	   // p2[70001];
	doublecomplex *crwork3, *cpart1, *cpart2;
	doublereal *k2gg, *k2pp, *k3gg, *k3pp, *p1, *p2;
} workspace_;

#define workspace_1 workspace_

Extern struct {
    //doublecomplex crfield[29241], crsource[29241], crmatc[171], cstep, crhm[29241], cwet[171], cbet[171];
	doublecomplex *crfield, *crsource, crmatc[513], cstep, *crhm, cwet[513], cbet[513];
    doublereal dxy, xks, radphase;
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

/* Subroutine */ int partsim_(tgam, tthet, dgam, dthet)
doublereal *tgam, *tthet, *dgam, *dthet;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, ddd_debug;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double sin(), cos(), sqrt(), d_imag();

    /* Local variables */
    static doublecomplex ctmp;
    static doublereal btper0, ztemp1, ztemp2;
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
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {

	i__2 = ip - 1;
	d__1 = cos(tthet[ip]);
	d__2 = -sin(tthet[ip]);
	z__2.r = d__1, z__2.i = d__2;
	z__1.r = workspace_1.cpart1[i__2].r * z__2.r - workspace_1.cpart1[
		i__2].i * z__2.i, z__1.i = workspace_1.cpart1[i__2].r * 
		z__2.i + workspace_1.cpart1[i__2].i * z__2.r;
	ctmp.r = z__1.r, ctmp.i = z__1.i;

	btper0 = beamcom_1.btper[ip - 1] + ztemp1 * ctmp.r;
/* perpendicular velocity */
/* Computing 2nd power */
	d__1 = tgam[ip];

	//OC TEST
	ddd_debug = 1. - btper0 / (d__1 * d__1);
	if(ddd_debug < 0) ddd_debug = 0.000001;
	//beamcom_1.btpar[ip - 1] = sqrt(1. - btper0 / (d__1 * d__1));
	beamcom_1.btpar[ip - 1] = sqrt(ddd_debug);

/* parallel velocity */
	dthet[ip] = dthet[ip] + ztemp2 * ((float)1. - (float)1. / 
		beamcom_1.btpar[ip - 1]) + (float)1.;
/* dtheta/dz */
	dgam[ip] = dgam[ip] + d_imag(&ctmp) / beamcom_1.btpar[ip - 1] / tgam[
		ip] - beamcom_1.ez[ip - 1];
/* dgamma/dz */
    }
/* ip */
    return 0;
} /* partsim_ */


/* partsim */
/* Subroutine */ int partsorc_(istepz)
integer *istepz;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg();
    double sqrt();

    /* Local variables */
    extern /* Subroutine */ int rpos_();
    static integer ip;
    static doublecomplex clocal;
    static doublereal zconst, aw2;
    static integer idx;
    static doublereal wei;
    extern doublereal faw2_();

/*     ================================================================== */
/*     calculates source term for gamma-theta integration */
/*     the particle transverse position and the radiation field are taken */
/*     at the middle of the integration step. */
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


/*     ------------------------------------------------------------------ */
/*     wiggler parameters */




/*     ------------------------------------------------------------------ */
/*     electron beam */




/*     simulation control and normalisation parameter */



/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







    rpos_(istepz, beamcom_1.xpart, beamcom_1.ypart);

/* position o */
    zconst = inputcom_1.xlamds / inputcom_1.xlamd * wigcom_1.fbess;
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	wei = beamcom_1.wx[ip - 1] * beamcom_1.wy[ip - 1];
	idx = beamcom_1.ipos[(ip << 2) - 4];
	i__2 = idx - 1;
	z__1.r = wei * cartcom_1.crfield[i__2].r, z__1.i = wei * 
		cartcom_1.crfield[i__2].i;
	clocal.r = z__1.r, clocal.i = z__1.i;
	wei = beamcom_1.wx[ip - 1] * (1. - beamcom_1.wy[ip - 1]);
	idx = beamcom_1.ipos[(ip << 2) - 3];
	i__2 = idx - 1;
	z__2.r = wei * cartcom_1.crfield[i__2].r, z__2.i = wei * 
		cartcom_1.crfield[i__2].i;
	z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	clocal.r = z__1.r, clocal.i = z__1.i;
	wei = (1. - beamcom_1.wx[ip - 1]) * beamcom_1.wy[ip - 1];
	idx = beamcom_1.ipos[(ip << 2) - 2];
	i__2 = idx - 1;
	z__2.r = wei * cartcom_1.crfield[i__2].r, z__2.i = wei * 
		cartcom_1.crfield[i__2].i;
	z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	clocal.r = z__1.r, clocal.i = z__1.i;
	wei = (1. - beamcom_1.wx[ip - 1]) * (1. - beamcom_1.wy[ip - 1]);
	idx = beamcom_1.ipos[(ip << 2) - 1];
	i__2 = idx - 1;
	z__2.r = wei * cartcom_1.crfield[i__2].r, z__2.i = wei * 
		cartcom_1.crfield[i__2].i;
	z__1.r = clocal.r + z__2.r, z__1.i = clocal.i + z__2.i;
	clocal.r = z__1.r, clocal.i = z__1.i;

/* field at particle position */
	aw2 = faw2_(istepz, &beamcom_1.xporb[ip - 1], &beamcom_1.yporb[ip - 1]
		);
/* square of wiggler fiel */
	aw2 += wigcom_1.awdz[*istepz - 1] * wigcom_1.awdz[*istepz - 1];
/* artificial delay */
	beamcom_1.btper[ip - 1] = aw2 + beamcom_1.px[ip - 1] * beamcom_1.px[
		ip - 1] + beamcom_1.py[ip - 1] * beamcom_1.py[ip - 1] + (
		float)1.;
	i__2 = ip - 1;
	d_cnjg(&z__3, &clocal);
	z__2.r = zconst * z__3.r, z__2.i = zconst * z__3.i;
	d__1 = sqrt(aw2);
	z__1.r = d__1 * z__2.r, z__1.i = d__1 * z__2.i;
	workspace_1.cpart1[i__2].r = z__1.r, workspace_1.cpart1[i__2].i = 
		z__1.i;
    }
/* ip */
    return 0;
} /* partsorc_ */

