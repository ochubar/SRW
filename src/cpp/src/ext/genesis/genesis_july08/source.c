/* source.f -- translated by f2c (version 20000118).
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
    //doublereal awz[10001], awdz[10000], solz[10000], awerx[10000], awery[
	   // 10000], qfld[10001], fbess, magversion, unitlength, dqfx[10001], 
	   // dqfy[10001], awdx[10001], awdy[10001];
    doublereal *awz, *awdz, *solz, *awerx, *awery, *qfld, 
		fbess, magversion, unitlength, *dqfx, *dqfy, *awdx, *awdy;
    integer iran, nstepz, itap;
} wigcom_;

#define wigcom_1 wigcom_

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

/* Subroutine */ int source_(istepz, i__)
integer *istepz, *i__;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double sqrt(), sin(), cos();

    /* Local variables */
    extern /* Subroutine */ int rpos_();
    static doublereal evencoupling;
    static integer j;
    static doublecomplex ctemp;
    static doublereal stemp, awloc;
    static integer ip, idx;
    static doublereal wei;
    extern doublereal faw2_();

/*     ================================================================== */
/*     construct the source for the wave equation */
/*     = radiation of the electron beam */
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





/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */








    i__1 = inputcom_1.ncar * inputcom_1.ncar;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j - 1;
	cartcom_1.crsource[i__2].r = 0., cartcom_1.crsource[i__2].i = 0.;
/* clear 2d array */
    }

/* j */
    if (wigcom_1.awz[*istepz] < 1e-25) {
	return 0;
    }

/* drift !!! */
    stemp = beamcom_1.xcuren * 3.6862103034005585e-4 / (real) 
	    inputcom_1.npart * inputcom_1.xlamd / inputcom_1.xlamds;
    stemp = stemp * inputcom_1.delz * 6.28318530717958 / cartcom_1.dxy / 
	    cartcom_1.dxy / 2.;

/*     debugged - adding harmonic number to the source term */

/* constant factor 2 because so */
    stemp *= (doublereal) (*i__);



    rpos_(istepz, beamcom_1.xpart, beamcom_1.ypart);

/* get particle position on grid */
    if (*i__ % 2 != 0) {

/*     calculating the odd harmonics */

	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
/* get undulator field at particle */
	    workspace_1.p1[ip - 1] = sqrt(faw2_(istepz, &beamcom_1.xporb[ip - 
		    1], &beamcom_1.yporb[ip - 1])) * stemp / beamcom_1.gamma[
		    ip - 1] / beamcom_1.btpar[ip - 1] * 
		    cartcom_1.besselcoupling[*i__ - 1];
	    i__2 = ip - 1;
	    i__3 = ip - 1;
	    z__1.r = 0., z__1.i = workspace_1.p1[i__3];
	    workspace_1.cpart1[i__2].r = z__1.r, workspace_1.cpart1[i__2].i = 
		    z__1.i;
	}

/* ip */
    } else {

/*     calculate the even harmonics, given by i */

	evencoupling = (doublereal) (*i__) * inputcom_1.xlamd / 
		inputcom_1.xlamds;
	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
	    awloc = sqrt(faw2_(istepz, &beamcom_1.xporb[ip - 1], &
		    beamcom_1.yporb[ip - 1]));
	    workspace_1.p1[ip - 1] = awloc * stemp / beamcom_1.gamma[ip - 1] /
		     beamcom_1.btpar[ip - 1] * cartcom_1.besselcoupling[*i__ 
		    - 1];
	    workspace_1.p1[ip - 1] = workspace_1.p1[ip - 1] / beamcom_1.gamma[
		    ip - 1] / beamcom_1.gamma[ip - 1] * evencoupling * awloc;
/* nk */
	    i__2 = ip - 1;
	    d__1 = -workspace_1.p1[ip - 1];
	    d__2 = sqrt((float)2.) * beamcom_1.px[ip - 1];
	    z__2.r = d__2, z__2.i = 0.;
	    z__1.r = d__1 * z__2.r, z__1.i = d__1 * z__2.i;
	    workspace_1.cpart1[i__2].r = z__1.r, workspace_1.cpart1[i__2].i = 
		    z__1.i;
/* pi/2 phas */
	}
    }

/*     load source term with local bunching factor */

    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	i__2 = ip - 1;
	d__1 = cos((doublereal) (*i__) * beamcom_1.theta[ip - 1]);
	d__2 = -sin((doublereal) (*i__) * beamcom_1.theta[ip - 1]);
	z__2.r = d__1, z__2.i = d__2;
	z__1.r = workspace_1.cpart1[i__2].r * z__2.r - workspace_1.cpart1[
		i__2].i * z__2.i, z__1.i = workspace_1.cpart1[i__2].r * 
		z__2.i + workspace_1.cpart1[i__2].i * z__2.r;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	wei = beamcom_1.wx[ip - 1] * beamcom_1.wy[ip - 1];
	idx = beamcom_1.ipos[(ip << 2) - 4];
	i__2 = idx - 1;
	i__3 = idx - 1;
	z__2.r = wei * ctemp.r, z__2.i = wei * ctemp.i;
	z__1.r = cartcom_1.crsource[i__3].r + z__2.r, z__1.i = 
		cartcom_1.crsource[i__3].i + z__2.i;
	cartcom_1.crsource[i__2].r = z__1.r, cartcom_1.crsource[i__2].i = 
		z__1.i;
	wei = beamcom_1.wx[ip - 1] * (1. - beamcom_1.wy[ip - 1]);
	idx = beamcom_1.ipos[(ip << 2) - 3];
	i__2 = idx - 1;
	i__3 = idx - 1;
	z__2.r = wei * ctemp.r, z__2.i = wei * ctemp.i;
	z__1.r = cartcom_1.crsource[i__3].r + z__2.r, z__1.i = 
		cartcom_1.crsource[i__3].i + z__2.i;
	cartcom_1.crsource[i__2].r = z__1.r, cartcom_1.crsource[i__2].i = 
		z__1.i;
	wei = (1. - beamcom_1.wx[ip - 1]) * beamcom_1.wy[ip - 1];
	idx = beamcom_1.ipos[(ip << 2) - 2];
	i__2 = idx - 1;
	i__3 = idx - 1;
	z__2.r = wei * ctemp.r, z__2.i = wei * ctemp.i;
	z__1.r = cartcom_1.crsource[i__3].r + z__2.r, z__1.i = 
		cartcom_1.crsource[i__3].i + z__2.i;
	cartcom_1.crsource[i__2].r = z__1.r, cartcom_1.crsource[i__2].i = 
		z__1.i;
	wei = (1. - beamcom_1.wx[ip - 1]) * (1. - beamcom_1.wy[ip - 1]);
	idx = beamcom_1.ipos[(ip << 2) - 1];
	i__2 = idx - 1;
	i__3 = idx - 1;
	z__2.r = wei * ctemp.r, z__2.i = wei * ctemp.i;
	z__1.r = cartcom_1.crsource[i__3].r + z__2.r, z__1.i = 
		cartcom_1.crsource[i__3].i + z__2.i;
	cartcom_1.crsource[i__2].r = z__1.r, cartcom_1.crsource[i__2].i = 
		z__1.i;
    }

    return 0;
} /* source_ */

