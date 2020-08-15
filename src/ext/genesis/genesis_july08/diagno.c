/* diagno.f -- translated by f2c (version 20000118).
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
    //doublecomplex crfield[476847], crsource[68121], crmatc[1827], cstep[7], crhm[68121], cwet[1827], cbet[1827];
    doublecomplex *crfield, *crsource, *crmatc, cstep[7], *crhm, *cwet, *cbet;
    doublereal dxy, xks, radphase, besselcoupling[7];
    integer nhloop, hloop[7];
} cartcom_;

#define cartcom_1 cartcom_

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

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int diagno_(istepz)
integer *istepz;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void d_cnjg();
    double log(), sqrt(), d_imag(), atan2(), sin(), cos();
    integer pow_ii();

    /* Local variables */
    static integer ioff;
    static doublereal xavg, yavg, ptot, wwcr;
    static doublecomplex ctmp;
    static integer i__, n, i0, i1, nctmp;
    static doublereal tpsin;
    static integer ip, ix, iy, nn[2];
    static doublereal tpcos, gainavg, xxsum, yysum, cr2, crsum, pradn;
    extern /* Subroutine */ int fourn_();

/*     ================================================================== */
/*     some diagnostics: */
/*     the radiation power must be calculated for each integration step */
/*     otherwise error will be wrong. */
/*     all calculation are stored in a history arrays which will be */
/*     written to a file ad the end of the run. */
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
/*     cartesian mesh */



/*     ------------------------------------------------------------------ */
/*     electron beam */




/*     ------------------------------------------------------------------ */
/*     diagnostic */



/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







/* unofficial */
/*     ------------------------------------------------------------------ */
/*     wiggler parameters */





    if (inputcom_1.iphsty <= 0) {
	return 0;
    }
/* no output at all */
    if (*istepz % inputcom_1.iphsty != 0) {
	return 0;
    }
/* no evaluation at this step */
    if (*istepz == 0) {
	diagcom_1.ihist = 1;
/* first initialization of step counter */
    } else {
	++diagcom_1.ihist;
/* advancing step counter for history c */
    }

/*     diagnostic: radiation field */
/*     ----------------------------------------------------------------- */

/*     radiation power */

    i__1 = cartcom_1.nhloop;
    for (n = 1; n <= i__1; ++n) {
/* looping over harmonics */
	crsum = 0.;
	ctmp.r = 0., ctmp.i = 0.;
	ioff = (n - 1) * inputcom_1.ncar * inputcom_1.ncar;
	i__2 = inputcom_1.ncar * inputcom_1.ncar + ioff;
	for (i__ = ioff + 1; i__ <= i__2; ++i__) {
	    i__3 = i__ - 1;
	    d_cnjg(&z__2, &cartcom_1.crfield[i__ - 1]);
	    z__1.r = cartcom_1.crfield[i__3].r * z__2.r - cartcom_1.crfield[
		    i__3].i * z__2.i, z__1.i = cartcom_1.crfield[i__3].r * 
		    z__2.i + cartcom_1.crfield[i__3].i * z__2.r;
	    wwcr = z__1.r;
/* =sum of |aij|^2 */
	    crsum += wwcr;
	    i__3 = i__ - 1;
	    z__1.r = ctmp.r + cartcom_1.crfield[i__3].r, z__1.i = ctmp.i + 
		    cartcom_1.crfield[i__3].i;
	    ctmp.r = z__1.r, ctmp.i = z__1.i;
	}
/* Computing 2nd power */
	d__1 = cartcom_1.dxy * 510999.06 * simcom_1.xkper0 / cartcom_1.xks / 
		cartcom_1.hloop[n - 1];
	pradn = crsum * (d__1 * d__1) / 376.73;
/* = current */
	if (n == 1) {
	    if (diagcom_1.pradoln[n - 1] > 0. && pradn > 0.) {
		diagcom_1.logp[diagcom_1.ihist - 1] = log(pradn / 
			diagcom_1.pradoln[n - 1]) / (inputcom_1.delz * 
			inputcom_1.xlamd);
/* log deriva */
	    } else {
		diagcom_1.logp[diagcom_1.ihist - 1] = (float)0.;
	    }
	}
	gainavg = (pradn + diagcom_1.pradoln[n - 1]) * .5;
/* average with old po */
	diagcom_1.pgainhist[cartcom_1.hloop[n - 1] + diagcom_1.ihist * 7 - 8] = pradn;
	diagcom_1.pradoln[n - 1] = pradn;

/*       on-axis far field intensity */

/* store actual value */
	z__5.r = ctmp.r * 510999.06, z__5.i = ctmp.i * 510999.06;
/* Computing 2nd power */
	d__2 = simcom_1.xkper0;
	d__1 = d__2 * d__2;
	z__4.r = d__1 * z__5.r, z__4.i = d__1 * z__5.i;
	z__3.r = z__4.r / cartcom_1.xks, z__3.i = z__4.i / cartcom_1.xks;
	i__2 = n - 1;
	d__3 = (doublereal) cartcom_1.hloop[i__2];
	z__2.r = z__3.r / d__3, z__2.i = z__3.i / d__3;
	d__4 = sqrt(376.73);
	z__1.r = z__2.r / d__4, z__1.i = z__2.i / d__4;
	ctmp.r = z__1.r, ctmp.i = z__1.i;
/* scale it to */
	d_cnjg(&z__2, &ctmp);
	z__1.r = ctmp.r * z__2.r - ctmp.i * z__2.i, z__1.i = ctmp.r * z__2.i 
		+ ctmp.i * z__2.r;
	diagcom_1.ffield[cartcom_1.hloop[n - 1] + diagcom_1.ihist * 7 - 8] = 
		z__1.r;

/*       on-axis near field intensity */

/* far field on-axi */
	i__ = inputcom_1.ncar * (inputcom_1.ncar - 1) / 2 + (inputcom_1.ncar 
		+ 1) / 2 + (n - 1) * inputcom_1.ncar * inputcom_1.ncar;
	i__2 = i__ - 1;
	d_cnjg(&z__2, &cartcom_1.crfield[i__ - 1]);
	z__1.r = cartcom_1.crfield[i__2].r * z__2.r - cartcom_1.crfield[i__2]
		.i * z__2.i, z__1.i = cartcom_1.crfield[i__2].r * z__2.i + 
		cartcom_1.crfield[i__2].i * z__2.r;
/* Computing 2nd power */
	d__1 = cartcom_1.dxy * 510999.06 * simcom_1.xkper0 / cartcom_1.xks / 
		cartcom_1.hloop[n - 1];
	diagcom_1.pmidhist[cartcom_1.hloop[n - 1] + diagcom_1.ihist * 7 - 8] =
		 z__1.r * (d__1 * d__1) / 376.73;
/* kg */
	diagcom_1.phimid[cartcom_1.hloop[n - 1] + diagcom_1.ihist * 7 - 8] = 
		atan2(d_imag(&cartcom_1.crfield[i__ - 1]), (doublereal) 
		cartcom_1.crfield[i__ - 1].r);
	if (inputcom_1.ffspec < 0) {
	    diagcom_1.phimid[cartcom_1.hloop[n - 1] + diagcom_1.ihist * 7 - 8]
		     = atan2(d_imag(&ctmp), (doublereal) ctmp.r);
	    diagcom_1.pmidhist[cartcom_1.hloop[n - 1] + diagcom_1.ihist * 7 - 
		    8] = diagcom_1.ffield[cartcom_1.hloop[n - 1] + 
		    diagcom_1.ihist * 7 - 8];
	}
	if (inputcom_1.ffspec > 0) {
	    diagcom_1.pmidhist[cartcom_1.hloop[n - 1] + diagcom_1.ihist * 7 - 
		    8] = diagcom_1.pgainhist[cartcom_1.hloop[n - 1] + 
		    diagcom_1.ihist * 7 - 8];
/* kg */
	}
    }

/*     radiation size  of fundamental */

/* end of loop over harmonics */
    if (diagcom_1.pradoln[0] > 0.) {
	crsum = (float)0.;
	xavg = 0.;
	yavg = 0.;
	cr2 = 0.;
	i__1 = inputcom_1.ncar;
	for (iy = 1; iy <= i__1; ++iy) {
	    i__2 = inputcom_1.ncar;
	    for (ix = 1; ix <= i__2; ++ix) {
		i__ = (iy - 1) * inputcom_1.ncar + ix;
		i__3 = i__ - 1;
		d_cnjg(&z__2, &cartcom_1.crfield[i__ - 1]);
		z__1.r = cartcom_1.crfield[i__3].r * z__2.r - 
			cartcom_1.crfield[i__3].i * z__2.i, z__1.i = 
			cartcom_1.crfield[i__3].r * z__2.i + 
			cartcom_1.crfield[i__3].i * z__2.r;
		wwcr = z__1.r;
		crsum += wwcr;
		xavg += wwcr * (doublereal) ix;
/* sum up positio */
		yavg += wwcr * (doublereal) iy;
		cr2 += wwcr * (doublereal) (ix * ix + iy * iy);
	    }
	}
	xavg /= crsum;
/* center of radiat */
	yavg /= crsum;
	cr2 = cr2 / crsum - xavg * xavg - yavg * yavg;
    } else {
	cr2 = 0.;
    }
    diagcom_1.whalf[diagcom_1.ihist - 1] = sqrt(cr2) * cartcom_1.dxy / 
	    simcom_1.xkper0;

/*     diagnostic: electron beam */
/*     ------------------------------------------------------------------ */

/*     energy */

/* rms radiation size */
    diagcom_1.gamhist[diagcom_1.ihist - 1] = (float)0.;
    diagcom_1.dgamhist[diagcom_1.ihist - 1] = (float)0.;
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	diagcom_1.gamhist[diagcom_1.ihist - 1] += beamcom_1.gamma[ip - 1] - 
		simcom_1.gamma0_in__;
/* Computing 2nd power */
	d__1 = beamcom_1.gamma[ip - 1] - simcom_1.gamma0_in__;
	diagcom_1.dgamhist[diagcom_1.ihist - 1] += d__1 * d__1;
    }
    if (inputcom_1.npart > 0) {
	diagcom_1.gamhist[diagcom_1.ihist - 1] /= (doublereal) 
		inputcom_1.npart;
/* Computing 2nd power */
	d__1 = diagcom_1.gamhist[diagcom_1.ihist - 1];
	diagcom_1.dgamhist[diagcom_1.ihist - 1] = sqrt(diagcom_1.dgamhist[
		diagcom_1.ihist - 1] / (doublereal) inputcom_1.npart - d__1 * 
		d__1);

    }
    ptot = diagcom_1.pradoln[0] + beamcom_1.xcuren * 510999.06 * 
	    diagcom_1.gamhist[diagcom_1.ihist - 1];
/* 1st part of total */
    if (*istepz == 0) {
	diagcom_1.pinit = ptot;
    }

/*     bunching at nharm harmonics */

    i__1 = cartcom_1.nhloop;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tpsin = 0.;
	tpcos = 0.;
	i__2 = inputcom_1.npart;
	for (ip = 1; ip <= i__2; ++ip) {
	    tpsin += sin(beamcom_1.theta[ip - 1] * (doublereal) i__);
/* add up phases */
	    tpcos += cos(beamcom_1.theta[ip - 1] * (doublereal) i__);
	}
/* ip */
	if (inputcom_1.npart > 0) {
	    tpsin /= (doublereal) inputcom_1.npart;
	    tpcos /= (doublereal) inputcom_1.npart;
	}
/* Computing 2nd power */
	d__1 = tpsin;
/* Computing 2nd power */
	d__2 = tpcos;
	diagcom_1.pmodhist[i__ + diagcom_1.ihist * 7 - 8] = sqrt(d__1 * d__1 + d__2 * d__2);
/* bunching factor */
	diagcom_1.bunphase[i__ + diagcom_1.ihist * 7 - 8] = atan2(tpsin, tpcos);
/* bunching phase */
    }
    for (i__ = cartcom_1.nhloop + 1; i__ <= 7; ++i__) {
	diagcom_1.pmodhist[i__ + diagcom_1.ihist * 7 - 8] = (float)0.;
	diagcom_1.bunphase[i__ + diagcom_1.ihist * 7 - 8] = (float)0.;
    }

/*     beam radius & energy spread */

    xxsum = 0.;
/* reset counter */
    yysum = 0.;
    diagcom_1.xpos[diagcom_1.ihist - 1] = 0.;
    diagcom_1.ypos[diagcom_1.ihist - 1] = 0.;

    i__1 = inputcom_1.npart;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = beamcom_1.xpart[i__ - 1];
	xxsum += d__1 * d__1;
/* sum of radii squared */
/* Computing 2nd power */
	d__1 = beamcom_1.ypart[i__ - 1];
	yysum += d__1 * d__1;
	diagcom_1.xpos[diagcom_1.ihist - 1] += beamcom_1.xpart[i__ - 1];
/* sum of radii */
	diagcom_1.ypos[diagcom_1.ihist - 1] += beamcom_1.ypart[i__ - 1];
    }
/* i */
    if (inputcom_1.npart > 0) {
	diagcom_1.xpos[diagcom_1.ihist - 1] = diagcom_1.xpos[diagcom_1.ihist 
		- 1] / (doublereal) inputcom_1.npart / simcom_1.xkper0;
/* mean value */
	diagcom_1.ypos[diagcom_1.ihist - 1] = diagcom_1.ypos[diagcom_1.ihist 
		- 1] / (doublereal) inputcom_1.npart / simcom_1.xkper0;
	xxsum = xxsum / (doublereal) inputcom_1.npart / simcom_1.xkper0 / 
		simcom_1.xkper0;
/* mean square value */
	yysum = yysum / (doublereal) inputcom_1.npart / simcom_1.xkper0 / 
		simcom_1.xkper0;
    }
/* Computing 2nd power */
    d__1 = diagcom_1.xpos[diagcom_1.ihist - 1];
    diagcom_1.xrms[diagcom_1.ihist - 1] = sqrt(xxsum - d__1 * d__1);
/* Computing 2nd power */
    d__1 = diagcom_1.ypos[diagcom_1.ihist - 1];
    diagcom_1.yrms[diagcom_1.ihist - 1] = sqrt(yysum - d__1 * d__1);

/*     energy conservation */
/*     ---------------------------------------------------------------- */

    diagcom_1.error[diagcom_1.ihist - 1] = (float)0.;
    if (gainavg != (float)0.) {
	diagcom_1.error[diagcom_1.ihist - 1] = (ptot / gainavg - 
		diagcom_1.pinit / gainavg) * 100.;
    }

    if (inputcom_1.lout[5] == 0) {
	return 0;
    }

/*     diffraction angle of radiation field */
/*     --------------------------------------------------------------- */

    if (diagcom_1.pradoln[0] == 0.) {
/* no radiation field */
	diagcom_1.diver[diagcom_1.ihist - 1] = (float)0.;
	return 0;
    }

    i__1 = (integer) (log((real) inputcom_1.ncar) / log((float)2.)) + 1;
    nctmp = pow_ii(&c__2, &i__1);
/* with nctmp> ncar */
    i__1 = nctmp * nctmp;
    for (i1 = 1; i1 <= i__1; ++i1) {
	i__2 = i1 - 1;
	workspace_1.crwork3[i__2].r = 0., workspace_1.crwork3[i__2].i = 0.;
/* clear working */
    }
    i__ = (nctmp - inputcom_1.ncar) / 2;
/* first index in bigge */
    i__1 = inputcom_1.ncar;
    for (ix = 1; ix <= i__1; ++ix) {
	i__2 = inputcom_1.ncar;
	for (iy = 1; iy <= i__2; ++iy) {
	    i0 = (ix - 1) * inputcom_1.ncar + iy;
	    i1 = (ix - 1 + i__) * nctmp + iy + i__;
	    i__3 = i1 - 1;
	    i__4 = i0 - 1;
	    workspace_1.crwork3[i__3].r = cartcom_1.crfield[i__4].r, 
		    workspace_1.crwork3[i__3].i = cartcom_1.crfield[i__4].i;
/* copy field around mi */
	}
    }
    nn[0] = nctmp;
/* size of mesh */
    nn[1] = nctmp;

/*     debug */

/*      crsum=0.0 */
/*      do i=1,ncar*ncar */
/*         =crsum+dble(crwork3(i)*conjg(crwork3(i))) */
/*      enddo */
/*      write (*,*) 'field power before fft',crsum */
    fourn_(workspace_1.crwork3, nn, &c__2, &c__1);

/*     debug */

/*      crsum=0.0 */
/*      do i=1,ncar*ncar */
/*         crsum=crsum+dble(crwork3(i)*conjg(crwork3(i))) */
/*      enddo */
/*      write (*,*) 'field power before fft',crsum */
/*      pause */

/* 2d fft with complex */
    i__1 = nctmp / 2;
    for (ix = 1; ix <= i__1; ++ix) {
/* rearrange fft output */
	i__2 = nctmp / 2;
	for (iy = 1; iy <= i__2; ++iy) {
	    i0 = (ix - 1) * nctmp + iy;
	    i1 = (ix - 1 + nctmp / 2) * nctmp + iy + nctmp / 2;
	    i__3 = i1 - 1;
	    ctmp.r = workspace_1.crwork3[i__3].r, ctmp.i = 
		    workspace_1.crwork3[i__3].i;
	    i__3 = i1 - 1;
	    i__4 = i0 - 1;
	    workspace_1.crwork3[i__3].r = workspace_1.crwork3[i__4].r, 
		    workspace_1.crwork3[i__3].i = workspace_1.crwork3[i__4].i;
	    i__3 = i0 - 1;
	    workspace_1.crwork3[i__3].r = ctmp.r, workspace_1.crwork3[i__3].i 
		    = ctmp.i;
	    i0 = (ix - 1) * nctmp + iy + nctmp / 2;
	    i1 = (ix + nctmp / 2 - 1) * nctmp + iy;
	    i__3 = i1 - 1;
	    ctmp.r = workspace_1.crwork3[i__3].r, ctmp.i = 
		    workspace_1.crwork3[i__3].i;
	    i__3 = i1 - 1;
	    i__4 = i0 - 1;
	    workspace_1.crwork3[i__3].r = workspace_1.crwork3[i__4].r, 
		    workspace_1.crwork3[i__3].i = workspace_1.crwork3[i__4].i;
	    i__3 = i0 - 1;
	    workspace_1.crwork3[i__3].r = ctmp.r, workspace_1.crwork3[i__3].i 
		    = ctmp.i;
	}
    }

    xavg = 0.;
    yavg = 0.;
    crsum = 0.;
    i__1 = nctmp;
    for (iy = 1; iy <= i__1; ++iy) {
	i__2 = nctmp;
	for (ix = 1; ix <= i__2; ++ix) {
	    i__ = (iy - 1) * nctmp + ix;
	    i__3 = i__ - 1;
	    d_cnjg(&z__2, &workspace_1.crwork3[i__ - 1]);
	    z__1.r = workspace_1.crwork3[i__3].r * z__2.r - 
		    workspace_1.crwork3[i__3].i * z__2.i, z__1.i = 
		    workspace_1.crwork3[i__3].r * z__2.i + 
		    workspace_1.crwork3[i__3].i * z__2.r;
	    wwcr = z__1.r;
	    xavg += wwcr * (doublereal) ix;
/* sum up spatial f */
	    yavg += wwcr * (doublereal) iy;
	    crsum += wwcr;
	}
    }
    if (crsum > (float)0.) {
	xavg /= crsum;
/* center of spat */
	yavg /= crsum;
    } else {
	xavg = (float)0.;
	yavg = (float)0.;
    }

    cr2 = 0.;
    i__1 = nctmp;
    for (iy = 1; iy <= i__1; ++iy) {
	i__2 = nctmp;
	for (ix = 1; ix <= i__2; ++ix) {
	    i__ = (iy - 1) * nctmp + ix;
	    i__3 = i__ - 1;
	    d_cnjg(&z__2, &workspace_1.crwork3[i__ - 1]);
	    z__1.r = workspace_1.crwork3[i__3].r * z__2.r - 
		    workspace_1.crwork3[i__3].i * z__2.i, z__1.i = 
		    workspace_1.crwork3[i__3].r * z__2.i + 
		    workspace_1.crwork3[i__3].i * z__2.r;
	    wwcr = z__1.r;
/* Computing 2nd power */
	    d__1 = (doublereal) ix - xavg;
/* Computing 2nd power */
	    d__2 = (doublereal) iy - yavg;
	    cr2 += wwcr * (d__1 * d__1 + d__2 * d__2);
	}
    }

    if (crsum <= (float)0.) {
	cr2 = (float)0.;
	crsum = (float)1.;
    }

    diagcom_1.diver[diagcom_1.ihist - 1] = sqrt(cr2 / crsum) * 
	    simcom_1.xkper0 * inputcom_1.xlamds / cartcom_1.dxy / (doublereal)
	     nctmp;

/* rms */
    return 0;
} /* diagno_ */

