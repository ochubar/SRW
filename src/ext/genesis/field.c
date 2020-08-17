/* field.f -- translated by f2c (version 20000118).
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

/* Table of constant values */

static doublecomplex c_b5 = {1.,0.};

/* Subroutine */ int field_(ihloop)
integer *ihloop;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9, i__10;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Local variables */
    static integer ioff, ix, idx;
    extern /* Subroutine */ int tridagx_(), tridagy_();

/*     ================================================================== */
/*     integrate the wave equation one step in z for cartesian mesh */
/*     using adi - methode (alternationg direction implicit) */
/*        1: u(n+1/2)=u(n)+alpha/2(dx u(n+1/2)+dy u(n)) */
/*        2: u(n+1)=u(n+1/2)+alpha/2(dx u(n+1/2)+dy u(n+1)) */
/*        to use tridiag methode transpose data array: */
/*        (1)->u(n)->u(i,j)->u(j,i)->u(n)->(2)->transpose again */
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



    ioff = (*ihloop - 1) * inputcom_1.ncar * inputcom_1.ncar;

/*     homogenious part to right hand side of diff equation (1) */
/*     ------------------------------------------------------------------ */

/* note that ihloop is the harmonic cou */
    i__1 = inputcom_1.ncar;
    for (ix = 1; ix <= i__1; ++ix) {
	i__2 = ix - 1;
	i__3 = ix - 1;
	i__4 = ix + ioff - 1;
	z__2.r = cartcom_1.crsource[i__3].r + cartcom_1.crfield[i__4].r, 
		z__2.i = cartcom_1.crsource[i__3].i + cartcom_1.crfield[i__4]
		.i;
	i__5 = *ihloop - 1;
	i__6 = ix + ioff + inputcom_1.ncar - 1;
	i__7 = ix + ioff - 1;
	z__5.r = cartcom_1.crfield[i__7].r * (float)2., z__5.i = 
		cartcom_1.crfield[i__7].i * (float)2.;
	z__4.r = cartcom_1.crfield[i__6].r - z__5.r, z__4.i = 
		cartcom_1.crfield[i__6].i - z__5.i;
	z__3.r = cartcom_1.cstep[i__5].r * z__4.r - cartcom_1.cstep[i__5].i * 
		z__4.i, z__3.i = cartcom_1.cstep[i__5].r * z__4.i + 
		cartcom_1.cstep[i__5].i * z__4.r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	cartcom_1.crhm[i__2].r = z__1.r, cartcom_1.crhm[i__2].i = z__1.i;
/* boundary :fie */
    }
    i__1 = inputcom_1.ncar * (inputcom_1.ncar - 1);
    for (idx = inputcom_1.ncar + 1; idx <= i__1; ++idx) {
	i__2 = idx - 1;
	i__3 = idx - 1;
	i__4 = idx + ioff - 1;
	z__2.r = cartcom_1.crsource[i__3].r + cartcom_1.crfield[i__4].r, 
		z__2.i = cartcom_1.crsource[i__3].i + cartcom_1.crfield[i__4]
		.i;
	i__5 = *ihloop - 1;
	i__6 = idx + inputcom_1.ncar + ioff - 1;
	i__7 = idx - inputcom_1.ncar + ioff - 1;
	z__5.r = cartcom_1.crfield[i__6].r + cartcom_1.crfield[i__7].r, 
		z__5.i = cartcom_1.crfield[i__6].i + cartcom_1.crfield[i__7]
		.i;
	i__8 = idx + ioff - 1;
	z__6.r = cartcom_1.crfield[i__8].r * (float)2., z__6.i = 
		cartcom_1.crfield[i__8].i * (float)2.;
	z__4.r = z__5.r - z__6.r, z__4.i = z__5.i - z__6.i;
	z__3.r = cartcom_1.cstep[i__5].r * z__4.r - cartcom_1.cstep[i__5].i * 
		z__4.i, z__3.i = cartcom_1.cstep[i__5].r * z__4.i + 
		cartcom_1.cstep[i__5].i * z__4.r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	cartcom_1.crhm[i__2].r = z__1.r, cartcom_1.crhm[i__2].i = z__1.i;
    }
    i__1 = inputcom_1.ncar * inputcom_1.ncar;
    for (idx = inputcom_1.ncar * (inputcom_1.ncar - 1) + 1; idx <= i__1; 
	    ++idx) {
	i__2 = idx - 1;
	i__3 = idx - 1;
	i__4 = idx + ioff - 1;
	z__2.r = cartcom_1.crsource[i__3].r + cartcom_1.crfield[i__4].r, 
		z__2.i = cartcom_1.crsource[i__3].i + cartcom_1.crfield[i__4]
		.i;
	i__5 = *ihloop - 1;
	i__6 = idx - inputcom_1.ncar + ioff - 1;
	i__7 = idx + ioff - 1;
	z__5.r = cartcom_1.crfield[i__7].r * (float)2., z__5.i = 
		cartcom_1.crfield[i__7].i * (float)2.;
	z__4.r = cartcom_1.crfield[i__6].r - z__5.r, z__4.i = 
		cartcom_1.crfield[i__6].i - z__5.i;
	z__3.r = cartcom_1.cstep[i__5].r * z__4.r - cartcom_1.cstep[i__5].i * 
		z__4.i, z__3.i = cartcom_1.cstep[i__5].r * z__4.i + 
		cartcom_1.cstep[i__5].i * z__4.r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	cartcom_1.crhm[i__2].r = z__1.r, cartcom_1.crhm[i__2].i = z__1.i;
/* bounda */
    }


/*     neumann boundary condition */
/*     ------------------------------------------------------------------ */

    if (inputcom_1.lbc != 0) {
	idx = inputcom_1.ncar * (inputcom_1.ncar - 1);
	i__1 = inputcom_1.ncar;
	for (ix = 1; ix <= i__1; ++ix) {
	    i__2 = ix - 1;
	    i__3 = ix - 1;
	    i__4 = *ihloop - 1;
	    i__5 = ix + inputcom_1.ncar + ioff - 1;
	    z__2.r = cartcom_1.cstep[i__4].r * cartcom_1.crfield[i__5].r - 
		    cartcom_1.cstep[i__4].i * cartcom_1.crfield[i__5].i, 
		    z__2.i = cartcom_1.cstep[i__4].r * cartcom_1.crfield[i__5]
		    .i + cartcom_1.cstep[i__4].i * cartcom_1.crfield[i__5].r;
	    z__1.r = cartcom_1.crhm[i__3].r + z__2.r, z__1.i = cartcom_1.crhm[
		    i__3].i + z__2.i;
	    cartcom_1.crhm[i__2].r = z__1.r, cartcom_1.crhm[i__2].i = z__1.i;
	    i__2 = idx + ix - 1;
	    i__3 = idx + ix - 1;
	    i__4 = *ihloop - 1;
	    i__5 = idx + ix - inputcom_1.ncar + ioff - 1;
	    z__2.r = cartcom_1.cstep[i__4].r * cartcom_1.crfield[i__5].r - 
		    cartcom_1.cstep[i__4].i * cartcom_1.crfield[i__5].i, 
		    z__2.i = cartcom_1.cstep[i__4].r * cartcom_1.crfield[i__5]
		    .i + cartcom_1.cstep[i__4].i * cartcom_1.crfield[i__5].r;
	    z__1.r = cartcom_1.crhm[i__3].r + z__2.r, z__1.i = cartcom_1.crhm[
		    i__3].i + z__2.i;
	    cartcom_1.crhm[i__2].r = z__1.r, cartcom_1.crhm[i__2].i = z__1.i;
	}
    }

/*     solve the tridiagonal system 1 */
/*     ------------------------------------------------------------------ */
    tridagx_(cartcom_1.crmatc, cartcom_1.crhm, cartcom_1.crfield, ihloop);


/*     homogenious part to right hand side of diff equation (2) */
/*     ------------------------------------------------------------------ */


    i__1 = inputcom_1.ncar * (inputcom_1.ncar - 1) + 1;
    i__2 = inputcom_1.ncar;
    for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	i__3 = ix - 1;
	i__4 = ix - 1;
	i__5 = ix + ioff - 1;
	z__2.r = cartcom_1.crsource[i__4].r + cartcom_1.crfield[i__5].r, 
		z__2.i = cartcom_1.crsource[i__4].i + cartcom_1.crfield[i__5]
		.i;
	i__6 = *ihloop - 1;
	i__7 = ix + 1 + ioff - 1;
	i__8 = ix + ioff - 1;
	z__5.r = cartcom_1.crfield[i__8].r * (float)2., z__5.i = 
		cartcom_1.crfield[i__8].i * (float)2.;
	z__4.r = cartcom_1.crfield[i__7].r - z__5.r, z__4.i = 
		cartcom_1.crfield[i__7].i - z__5.i;
	z__3.r = cartcom_1.cstep[i__6].r * z__4.r - cartcom_1.cstep[i__6].i * 
		z__4.i, z__3.i = cartcom_1.cstep[i__6].r * z__4.i + 
		cartcom_1.cstep[i__6].i * z__4.r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	cartcom_1.crhm[i__3].r = z__1.r, cartcom_1.crhm[i__3].i = z__1.i;
	i__3 = ix + inputcom_1.ncar - 2;
	for (idx = ix + 1; idx <= i__3; ++idx) {
	    i__4 = idx - 1;
	    i__5 = idx - 1;
	    i__6 = idx + ioff - 1;
	    z__2.r = cartcom_1.crsource[i__5].r + cartcom_1.crfield[i__6].r, 
		    z__2.i = cartcom_1.crsource[i__5].i + cartcom_1.crfield[
		    i__6].i;
	    i__7 = *ihloop - 1;
	    i__8 = idx + 1 + ioff - 1;
	    i__9 = idx - 1 + ioff - 1;
	    z__5.r = cartcom_1.crfield[i__8].r + cartcom_1.crfield[i__9].r, 
		    z__5.i = cartcom_1.crfield[i__8].i + cartcom_1.crfield[
		    i__9].i;
	    i__10 = idx + ioff - 1;
	    z__6.r = cartcom_1.crfield[i__10].r * (float)2., z__6.i = 
		    cartcom_1.crfield[i__10].i * (float)2.;
	    z__4.r = z__5.r - z__6.r, z__4.i = z__5.i - z__6.i;
	    z__3.r = cartcom_1.cstep[i__7].r * z__4.r - cartcom_1.cstep[i__7]
		    .i * z__4.i, z__3.i = cartcom_1.cstep[i__7].r * z__4.i + 
		    cartcom_1.cstep[i__7].i * z__4.r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    cartcom_1.crhm[i__4].r = z__1.r, cartcom_1.crhm[i__4].i = z__1.i;
	}
	idx = ix + inputcom_1.ncar - 1;
	i__3 = idx - 1;
	i__4 = idx - 1;
	i__5 = idx + ioff - 1;
	z__2.r = cartcom_1.crsource[i__4].r + cartcom_1.crfield[i__5].r, 
		z__2.i = cartcom_1.crsource[i__4].i + cartcom_1.crfield[i__5]
		.i;
	i__6 = *ihloop - 1;
	i__7 = idx - 1 + ioff - 1;
	i__8 = idx + ioff - 1;
	z__5.r = cartcom_1.crfield[i__8].r * (float)2., z__5.i = 
		cartcom_1.crfield[i__8].i * (float)2.;
	z__4.r = cartcom_1.crfield[i__7].r - z__5.r, z__4.i = 
		cartcom_1.crfield[i__7].i - z__5.i;
	z__3.r = cartcom_1.cstep[i__6].r * z__4.r - cartcom_1.cstep[i__6].i * 
		z__4.i, z__3.i = cartcom_1.cstep[i__6].r * z__4.i + 
		cartcom_1.cstep[i__6].i * z__4.r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	cartcom_1.crhm[i__3].r = z__1.r, cartcom_1.crhm[i__3].i = z__1.i;
    }


/*     neumann boundary condition */
/*     ------------------------------------------------------------------ */

    if (inputcom_1.lbc != 0) {
	i__2 = inputcom_1.ncar;
	for (ix = 1; ix <= i__2; ++ix) {
	    idx = inputcom_1.ncar * (ix - 1) + 1;
	    i__1 = idx - 1;
	    i__3 = idx - 1;
	    i__4 = *ihloop - 1;
	    i__5 = idx + 1 + ioff - 1;
	    z__2.r = cartcom_1.cstep[i__4].r * cartcom_1.crfield[i__5].r - 
		    cartcom_1.cstep[i__4].i * cartcom_1.crfield[i__5].i, 
		    z__2.i = cartcom_1.cstep[i__4].r * cartcom_1.crfield[i__5]
		    .i + cartcom_1.cstep[i__4].i * cartcom_1.crfield[i__5].r;
	    z__1.r = cartcom_1.crhm[i__3].r + z__2.r, z__1.i = cartcom_1.crhm[
		    i__3].i + z__2.i;
	    cartcom_1.crhm[i__1].r = z__1.r, cartcom_1.crhm[i__1].i = z__1.i;
	    idx = idx + inputcom_1.ncar - 1;
	    i__1 = idx - 1;
	    i__3 = idx - 1;
	    i__4 = *ihloop - 1;
	    i__5 = idx - 1 + ioff - 1;
	    z__2.r = cartcom_1.cstep[i__4].r * cartcom_1.crfield[i__5].r - 
		    cartcom_1.cstep[i__4].i * cartcom_1.crfield[i__5].i, 
		    z__2.i = cartcom_1.cstep[i__4].r * cartcom_1.crfield[i__5]
		    .i + cartcom_1.cstep[i__4].i * cartcom_1.crfield[i__5].r;
	    z__1.r = cartcom_1.crhm[i__3].r + z__2.r, z__1.i = cartcom_1.crhm[
		    i__3].i + z__2.i;
	    cartcom_1.crhm[i__1].r = z__1.r, cartcom_1.crhm[i__1].i = z__1.i;
	}
    }

/*     solve the tridiagonal system 2 */
/*     ------------------------------------------------------------------ */

    tridagy_(cartcom_1.crmatc, cartcom_1.crhm, cartcom_1.crfield, ihloop);

    return 0;
} /* field_ */



/* fieldcar */
/* Subroutine */ int tridagx_(c__, r__, u, h__)
doublecomplex *c__, *r__, *u;
integer *h__;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer ioff1, ioff2, k, i__;

/*     ================================================================== */
/*     solve a tridiagonal system for cartesian mesh in x direction */
/*     cbet and cwet are precalculated in auxval */
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


/*     -------------------------------------------------------------------- */
/*     cartesian mesh */




    /* Parameter adjustments */
    --u;
    --r__;
    --c__;

    /* Function Body */
    ioff1 = inputcom_1.ncar * (*h__ - 1);
    ioff2 = inputcom_1.ncar * ioff1;

    i__1 = inputcom_1.ncar * (inputcom_1.ncar - 1);
    i__2 = inputcom_1.ncar;
    for (i__ = 0; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	i__3 = i__ + 1 + ioff2;
	i__4 = i__ + 1;
	i__5 = ioff1;
	z__1.r = r__[i__4].r * cartcom_1.cbet[i__5].r - r__[i__4].i * 
		cartcom_1.cbet[i__5].i, z__1.i = r__[i__4].r * cartcom_1.cbet[
		i__5].i + r__[i__4].i * cartcom_1.cbet[i__5].r;
	u[i__3].r = z__1.r, u[i__3].i = z__1.i;
	i__3 = inputcom_1.ncar;
	for (k = 2; k <= i__3; ++k) {
	    i__4 = k + i__ + ioff2;
	    i__5 = k + i__;
	    i__6 = ioff1 + k;
	    i__7 = k + i__ - 1 + ioff2;
	    z__3.r = c__[i__6].r * u[i__7].r - c__[i__6].i * u[i__7].i, 
		    z__3.i = c__[i__6].r * u[i__7].i + c__[i__6].i * u[i__7]
		    .r;
	    z__2.r = r__[i__5].r - z__3.r, z__2.i = r__[i__5].i - z__3.i;
	    i__8 = k + ioff1 - 1;
	    z__1.r = z__2.r * cartcom_1.cbet[i__8].r - z__2.i * 
		    cartcom_1.cbet[i__8].i, z__1.i = z__2.r * cartcom_1.cbet[
		    i__8].i + z__2.i * cartcom_1.cbet[i__8].r;
	    u[i__4].r = z__1.r, u[i__4].i = z__1.i;
	}
/* k */
	for (k = inputcom_1.ncar - 1; k >= 1; --k) {
	    i__3 = k + i__ + ioff2;
	    i__4 = k + i__ + ioff2;
	    i__5 = k + 1 + ioff1 - 1;
	    i__6 = k + i__ + 1 + ioff2;
	    z__2.r = cartcom_1.cwet[i__5].r * u[i__6].r - cartcom_1.cwet[i__5]
		    .i * u[i__6].i, z__2.i = cartcom_1.cwet[i__5].r * u[i__6]
		    .i + cartcom_1.cwet[i__5].i * u[i__6].r;
	    z__1.r = u[i__4].r - z__2.r, z__1.i = u[i__4].i - z__2.i;
	    u[i__3].r = z__1.r, u[i__3].i = z__1.i;
	}
    }

    return 0;
} /* tridagx_ */



/* tridag */
/* Subroutine */ int tridagy_(c__, r__, u, h__)
doublecomplex *c__, *r__, *u;
integer *h__;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer ioff1, ioff2, k, i__, n;

/*     tridagy(crmatc,crhm,crfield,ihloop) */
/*     ================================================================== */
/*     solve a tridiagonal system for cartesian mesh in y direction */
/*     cbet and cwet are precalculated in auxval */
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



    /* Parameter adjustments */
    --u;
    --r__;
    --c__;

    /* Function Body */
    ioff1 = inputcom_1.ncar * (*h__ - 1);
    ioff2 = inputcom_1.ncar * ioff1;

    i__1 = inputcom_1.ncar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + ioff2;
	i__3 = i__;
	i__4 = ioff1;
	z__1.r = r__[i__3].r * cartcom_1.cbet[i__4].r - r__[i__3].i * 
		cartcom_1.cbet[i__4].i, z__1.i = r__[i__3].r * cartcom_1.cbet[
		i__4].i + r__[i__3].i * cartcom_1.cbet[i__4].r;
	u[i__2].r = z__1.r, u[i__2].i = z__1.i;
    }
    i__1 = inputcom_1.ncar;
    for (k = 2; k <= i__1; ++k) {
	n = k * inputcom_1.ncar - inputcom_1.ncar;
	i__2 = inputcom_1.ncar;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = n + i__ + ioff2;
	    i__4 = n + i__;
	    i__5 = ioff1 + k;
	    i__6 = n + i__ - inputcom_1.ncar + ioff2;
	    z__3.r = c__[i__5].r * u[i__6].r - c__[i__5].i * u[i__6].i, 
		    z__3.i = c__[i__5].r * u[i__6].i + c__[i__5].i * u[i__6]
		    .r;
	    z__2.r = r__[i__4].r - z__3.r, z__2.i = r__[i__4].i - z__3.i;
	    i__7 = k + ioff1 - 1;
	    z__1.r = z__2.r * cartcom_1.cbet[i__7].r - z__2.i * 
		    cartcom_1.cbet[i__7].i, z__1.i = z__2.r * cartcom_1.cbet[
		    i__7].i + z__2.i * cartcom_1.cbet[i__7].r;
	    u[i__3].r = z__1.r, u[i__3].i = z__1.i;
	}
    }
    for (k = inputcom_1.ncar - 1; k >= 1; --k) {
	n = k * inputcom_1.ncar - inputcom_1.ncar;
	i__1 = inputcom_1.ncar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = n + i__ + ioff2;
	    i__3 = n + i__ + ioff2;
	    i__4 = k + 1 + ioff1 - 1;
	    i__5 = n + i__ + inputcom_1.ncar + ioff2;
	    z__2.r = cartcom_1.cwet[i__4].r * u[i__5].r - cartcom_1.cwet[i__4]
		    .i * u[i__5].i, z__2.i = cartcom_1.cwet[i__4].r * u[i__5]
		    .i + cartcom_1.cwet[i__4].i * u[i__5].r;
	    z__1.r = u[i__3].r - z__2.r, z__1.i = u[i__3].i - z__2.i;
	    u[i__2].r = z__1.r, u[i__2].i = z__1.i;
	}
    }

    return 0;
} /* tridagy_ */




/* tridag */
/* Subroutine */ int getdiag_(stepsize, gridsize, wavenumber)
doublereal *stepsize, *gridsize, *wavenumber;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div();

    /* Local variables */
    static integer icar;
    static doublereal mmid[261], mlow[261], mupp[261], rtmp;
    static doublecomplex cwrk1[261], cwrk2[261];
    static integer ix, ihloop;

/*     ====================================================================== */
/*     construct the diagonal matrix for field equation */
/*     do some precalculation for field solver */
/*     ---------------------------------------------------------------------- */





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







/*     construction of the diagonal maxtrix for cartesian mesh */
/*     -------------------------------------------------------------------- */
    i__1 = cartcom_1.nhloop;
    for (ihloop = 1; ihloop <= i__1; ++ihloop) {
	rtmp = *stepsize * .25 / (*wavenumber * (doublereal) cartcom_1.hloop[
		ihloop - 1]);
	rtmp /= *gridsize * *gridsize;
/* factor dz/(4 ks dx */
	i__2 = ihloop - 1;
	z__1.r = 0., z__1.i = rtmp;
	cartcom_1.cstep[i__2].r = z__1.r, cartcom_1.cstep[i__2].i = z__1.i;
/* complex value - see field eq */
	if (inputcom_1.lbc != 0) {
	    inputcom_1.lbc = 1;
	}
/* boundary condition */
	mupp[0] = rtmp;
/* one edge of mesh */
	mmid[0] = -((doublereal) (2 - inputcom_1.lbc)) * rtmp;
/* boundary condition a=0 or da/dz=0 */
	mlow[0] = 0.;
	i__2 = inputcom_1.ncar - 1;
	for (ix = 2; ix <= i__2; ++ix) {
	    mupp[ix - 1] = rtmp;
/* inside of mesh -> 2nd derivation possibl */
	    mmid[ix - 1] = rtmp * -2.;
	    mlow[ix - 1] = rtmp;
	}
	mupp[inputcom_1.ncar - 1] = 0.;
/* other edge of mesh */
	mmid[inputcom_1.ncar - 1] = -((doublereal) (2 - inputcom_1.lbc)) * 
		rtmp;
	mlow[inputcom_1.ncar - 1] = rtmp;

/*     construct complex matrix crmat=(i-im) for */
/*     field equation  (i-im)*a(t+1)=(i+im)a(t) */
/*     ------------------------------------------------------------------------- */

	i__2 = inputcom_1.ncar;
	for (icar = 1; icar <= i__2; ++icar) {
	    i__3 = icar - 1;
	    i__4 = icar - 1;
	    z__2.r = 0., z__2.i = mupp[i__4];
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    cwrk1[i__3].r = z__1.r, cwrk1[i__3].i = z__1.i;
/* store value t */
	    i__3 = icar - 1;
	    i__4 = icar - 1;
	    z__2.r = 0., z__2.i = mmid[i__4];
	    z__1.r = 1. - z__2.r, z__1.i = 0. - z__2.i;
	    cwrk2[i__3].r = z__1.r, cwrk2[i__3].i = z__1.i;
/* same here */
	    i__3 = inputcom_1.ncar * (ihloop - 1) + icar - 1;
	    i__4 = icar - 1;
	    z__2.r = 0., z__2.i = mlow[i__4];
	    z__1.r = -z__2.r, z__1.i = -z__2.i;
	    cartcom_1.crmatc[i__3].r = z__1.r, cartcom_1.crmatc[i__3].i = 
		    z__1.i;
	}

/*     precalculated constants for tridiag subroutine */
/*     ------------------------------------------------------------------------ */
	i__2 = (ihloop - 1) * inputcom_1.ncar;
	z_div(&z__1, &c_b5, cwrk2);
	cartcom_1.cbet[i__2].r = z__1.r, cartcom_1.cbet[i__2].i = z__1.i;
	i__2 = (ihloop - 1) * inputcom_1.ncar;
	cartcom_1.cwet[i__2].r = (float)0., cartcom_1.cwet[i__2].i = (float)
		0.;
	i__2 = inputcom_1.ncar;
	for (icar = 2; icar <= i__2; ++icar) {
	    i__3 = (ihloop - 1) * inputcom_1.ncar + icar - 1;
	    i__4 = icar - 2;
	    i__5 = (ihloop - 1) * inputcom_1.ncar + icar - 2;
	    z__1.r = cwrk1[i__4].r * cartcom_1.cbet[i__5].r - cwrk1[i__4].i * 
		    cartcom_1.cbet[i__5].i, z__1.i = cwrk1[i__4].r * 
		    cartcom_1.cbet[i__5].i + cwrk1[i__4].i * cartcom_1.cbet[
		    i__5].r;
	    cartcom_1.cwet[i__3].r = z__1.r, cartcom_1.cwet[i__3].i = z__1.i;
	    i__3 = (ihloop - 1) * inputcom_1.ncar + icar - 1;
	    i__4 = icar - 1;
	    i__5 = inputcom_1.ncar * (ihloop - 1) + icar - 1;
	    i__6 = (ihloop - 1) * inputcom_1.ncar + icar - 1;
	    z__3.r = cartcom_1.crmatc[i__5].r * cartcom_1.cwet[i__6].r - 
		    cartcom_1.crmatc[i__5].i * cartcom_1.cwet[i__6].i, z__3.i 
		    = cartcom_1.crmatc[i__5].r * cartcom_1.cwet[i__6].i + 
		    cartcom_1.crmatc[i__5].i * cartcom_1.cwet[i__6].r;
	    z__2.r = cwrk2[i__4].r - z__3.r, z__2.i = cwrk2[i__4].i - z__3.i;
	    z_div(&z__1, &c_b5, &z__2);
	    cartcom_1.cbet[i__3].r = z__1.r, cartcom_1.cbet[i__3].i = z__1.i;
	}
    }
    return 0;
} /* getdiag_ */

