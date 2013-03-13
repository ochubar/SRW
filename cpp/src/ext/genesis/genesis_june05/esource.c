/* esource.f -- translated by f2c (version 20000118).
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

/* Subroutine */ int esource_(istepz)
integer *istepz;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(), log();
    void pow_zi();
    double sin(), cos();
    void z_div();

    /* Local variables */
    static doublecomplex coef;
    static doublereal drsc, rmid[100], rdig[100], xmid, ymid, rlog[100];
    static integer j, m;
    static doublecomplex cscsource[100];
    static doublereal ezmax;
    static doublecomplex crtmp1[100], crtmp2[100];
    static integer ip, ir;
    static doublecomplex vn;
    static doublereal econst, rscmax;
    extern /* Subroutine */ int trirad_();
    static doublecomplex cma[100], cmb[100], cmc[100];
    static doublereal vol[100], xks, xkw0;

/*     ================================================================== */
/*     calculates the space charge field. */
/*     all particle are discretized on a radial mesh for all */
/*     selected azimutal and longitudinal fourier modes. */
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








    if (inputcom_1.nscz <= 0) {
	return 0;
    }

/* no space charge selec */
    xks = 6.28318530717958 / inputcom_1.xlamds;
    xkw0 = 6.28318530717958 / inputcom_1.xlamd;
    xmid = 0.;
/* get centroid position */
    ymid = 0.;
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	xmid += beamcom_1.xpart[ip - 1];
	ymid += beamcom_1.ypart[ip - 1];
    }
/* i */
    if (inputcom_1.npart > 0) {
	xmid /= (doublereal) inputcom_1.npart;
/* mean value */
	ymid /= (doublereal) inputcom_1.npart;
    }

    rscmax = (float)0.;
/* clear max radius */
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.ez[ip - 1] = 0.;
/* clear old space charg */
/* Computing 2nd power */
	d__1 = beamcom_1.xpart[ip - 1] - xmid;
/* Computing 2nd power */
	d__2 = beamcom_1.ypart[ip - 1] - ymid;
	workspace_1.p1[ip - 1] = sqrt(d__1 * d__1 + d__2 * d__2) / xkw0;
/* ge */
	if (workspace_1.p1[ip - 1] > rscmax) {
	    rscmax = workspace_1.p1[ip - 1];
	}
/* look for maximum radi */
	if (workspace_1.p1[ip - 1] == (float)0.) {
	    i__2 = ip - 1;
	    workspace_1.cpart1[i__2].r = 1., workspace_1.cpart1[i__2].i = 0.;
/* particle in origin */
	} else {
	    i__2 = ip - 1;
	    d__1 = beamcom_1.ypart[ip - 1] - ymid;
	    d__2 = -beamcom_1.xpart[ip - 1] + xmid;
	    z__2.r = d__1, z__2.i = d__2;
	    i__3 = ip - 1;
	    z__1.r = z__2.r / workspace_1.p1[i__3], z__1.i = z__2.i / 
		    workspace_1.p1[i__3];
	    workspace_1.cpart1[i__2].r = z__1.r, workspace_1.cpart1[i__2].i = 
		    z__1.i;
/* az */
	}
    }

    rscmax *= (float)2.;
/* maximum of radial mes */
    if (rscmax <= 0.) {
	rscmax = inputcom_1.dgrid * xkw0;
    }
/*                                                 !set to radiation grid size */

/* catch error of zero s */
    drsc = rscmax / (real) (inputcom_1.nptr - 1);
/* gridpoint distance */
    i__1 = inputcom_1.nptr;
    for (ir = 2; ir <= i__1; ++ir) {
	vol[ir - 1] = drsc * 3.14159265358979 * drsc * (ir * (float)2. - 1);
/* 2d volume around grid */
	rlog[ir - 1] = log((real) ir / (real) (ir - 1));
/* log term of higher mo */
	rdig[ir - 1] = (real) (ir - 1) * 6.28318530717958;
/* diagonal terms above/ */
    }
    vol[0] = drsc * 3.14159265358979 * drsc;
/* volume of origin */
    rlog[0] = log((float).5);
/* shielding radius */
    rdig[0] = (float)0.;
/* no lower element at o */
    rdig[inputcom_1.nptr] = (float)0.;

/* no upper element at b */
    econst = beamcom_1.xcuren * 7.3724206068011169e-4 / (real) 
	    inputcom_1.npart / (xks + xkw0);
/* source term nor */
/* Computing 2nd power */
    d__2 = xks;
/* Computing 2nd power */
    d__3 = xks + xkw0;
    d__1 = 1. / (d__2 * d__2 - d__3 * d__3);
    z__1.r = d__1, z__1.i = 0.;
    coef.r = z__1.r, coef.i = z__1.i;

/* matrix normalisatio */
    i__1 = inputcom_1.nscr;
    for (m = -inputcom_1.nscr; m <= i__1; ++m) {

/* count over azimutal m */
	i__2 = inputcom_1.npart;
	for (ip = 1; ip <= i__2; ++ip) {
/* get right azimutal co */
	    i__3 = ip - 1;
	    pow_zi(&z__1, &workspace_1.cpart1[ip - 1], &m);
	    workspace_1.cpart2[i__3].r = z__1.r, workspace_1.cpart2[i__3].i = 
		    z__1.i;
	}

	i__2 = inputcom_1.nptr;
	for (ir = 1; ir <= i__2; ++ir) {
	    rmid[ir - 1] = -rdig[ir - 1] - rdig[ir] - (real) (m * m) * 
		    6.28318530717958 * rlog[ir - 1];
/* main */
	}
	rmid[inputcom_1.nptr - 1] -= (real) inputcom_1.nptr * 
		6.28318530717958;

/* direchlet boundary */
	i__2 = inputcom_1.nscz;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = inputcom_1.nptr;
	    for (ir = 1; ir <= i__3; ++ir) {
		i__4 = ir - 1;
		cscsource[i__4].r = 0., cscsource[i__4].i = 0.;
/* clear source */
	    }
	    i__3 = inputcom_1.npart;
	    for (ip = 1; ip <= i__3; ++ip) {
		ir = (integer) (workspace_1.p1[ip - 1] / drsc) + 1;
/* radial position */
		i__4 = ir - 1;
		i__5 = ir - 1;
		i__6 = ip - 1;
		d__1 = cos(beamcom_1.theta[ip - 1] * j);
		d__2 = -sin(beamcom_1.theta[ip - 1] * j);
		z__3.r = d__1, z__3.i = d__2;
		z__2.r = workspace_1.cpart2[i__6].r * z__3.r - 
			workspace_1.cpart2[i__6].i * z__3.i, z__2.i = 
			workspace_1.cpart2[i__6].r * z__3.i + 
			workspace_1.cpart2[i__6].i * z__3.r;
		z__1.r = cscsource[i__5].r + z__2.r, z__1.i = cscsource[i__5]
			.i + z__2.i;
		cscsource[i__4].r = z__1.r, cscsource[i__4].i = z__1.i;
/* add long and azi */
	    }
	    i__3 = inputcom_1.nptr;
	    for (ir = 1; ir <= i__3; ++ir) {
		d__1 = econst / (real) j / vol[ir - 1];
		z__1.r = 0., z__1.i = d__1;
		vn.r = z__1.r, vn.i = z__1.i;
/* complex norm. te */
		i__4 = ir - 1;
		i__5 = ir - 1;
		z__1.r = vn.r * cscsource[i__5].r - vn.i * cscsource[i__5].i, 
			z__1.i = vn.r * cscsource[i__5].i + vn.i * cscsource[
			i__5].r;
		cscsource[i__4].r = z__1.r, cscsource[i__4].i = z__1.i;
/* scale source t */
		i__4 = ir - 1;
		d__1 = rdig[ir - 1] / j / j / vol[ir - 1];
		z__2.r = d__1, z__2.i = 0.;
		z__1.r = coef.r * z__2.r - coef.i * z__2.i, z__1.i = coef.r * 
			z__2.i + coef.i * z__2.r;
		cma[i__4].r = z__1.r, cma[i__4].i = z__1.i;
/* construc */
		i__4 = ir - 1;
		d__1 = rmid[ir - 1] / j / j / vol[ir - 1];
		z__3.r = d__1, z__3.i = 0.;
		z__2.r = coef.r * z__3.r - coef.i * z__3.i, z__2.i = coef.r * 
			z__3.i + coef.i * z__3.r;
		z__1.r = z__2.r + 1., z__1.i = z__2.i + 0.;
		cmb[i__4].r = z__1.r, cmb[i__4].i = z__1.i;
		i__4 = ir - 1;
		d__1 = rdig[ir] / j / j / vol[ir - 1];
		z__2.r = d__1, z__2.i = 0.;
		z__1.r = coef.r * z__2.r - coef.i * z__2.i, z__1.i = coef.r * 
			z__2.i + coef.i * z__2.r;
		cmc[i__4].r = z__1.r, cmc[i__4].i = z__1.i;
	    }

	    trirad_(cma, cmb, cmc, cscsource, crtmp1, crtmp2, &
		    inputcom_1.nptr);

/* solve */
	    i__3 = inputcom_1.npart;
	    for (ip = 1; ip <= i__3; ++ip) {
/* sum up fourier coeffi */
		ir = (integer) (workspace_1.p1[ip - 1] / drsc) + 1;
		z_div(&z__2, &crtmp1[ir - 1], &workspace_1.cpart2[ip - 1]);
		d__1 = cos(beamcom_1.theta[ip - 1] * j);
		d__2 = sin(beamcom_1.theta[ip - 1] * j);
		z__3.r = d__1, z__3.i = d__2;
		z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * 
			z__3.i + z__2.i * z__3.r;
		beamcom_1.ez[ip - 1] += z__1.r * (float)2.;
	    }
	}
    }
    ezmax = (float)0.;
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.ez[ip - 1] /= xkw0;
/* scale due to normalized z */
	if (ezmax < beamcom_1.ez[ip - 1]) {
	    ezmax = beamcom_1.ez[ip - 1];
	}
    }
    return 0;
} /* esource_ */

/* Subroutine */ int trirad_(a, b, c__, r__, u, w, n)
doublecomplex *a, *b, *c__, *r__, *u, *w;
integer *n;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div();

    /* Local variables */
    static integer k;
    static doublecomplex bet;

/*     ================================================================== */
/*     solve a tridiagonal system for radial mesh */
/*     only called by esource for space charge calculation */
/*     ------------------------------------------------------------------ */


    /* Parameter adjustments */
    --w;
    --u;
    --r__;
    --c__;
    --b;
    --a;

    /* Function Body */
    bet.r = b[1].r, bet.i = b[1].i;
    z_div(&z__1, &r__[1], &bet);
    u[1].r = z__1.r, u[1].i = z__1.i;
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	i__2 = k;
	z_div(&z__1, &c__[k - 1], &bet);
	w[i__2].r = z__1.r, w[i__2].i = z__1.i;
	i__2 = k;
	i__3 = k;
	i__4 = k;
	z__2.r = a[i__3].r * w[i__4].r - a[i__3].i * w[i__4].i, z__2.i = a[
		i__3].r * w[i__4].i + a[i__3].i * w[i__4].r;
	z__1.r = b[i__2].r - z__2.r, z__1.i = b[i__2].i - z__2.i;
	bet.r = z__1.r, bet.i = z__1.i;
	i__2 = k;
	i__3 = k;
	i__4 = k;
	i__5 = k - 1;
	z__3.r = a[i__4].r * u[i__5].r - a[i__4].i * u[i__5].i, z__3.i = a[
		i__4].r * u[i__5].i + a[i__4].i * u[i__5].r;
	z__2.r = r__[i__3].r - z__3.r, z__2.i = r__[i__3].i - z__3.i;
	z_div(&z__1, &z__2, &bet);
	u[i__2].r = z__1.r, u[i__2].i = z__1.i;
    }
/* k */
    for (k = *n - 1; k >= 1; --k) {
	i__1 = k;
	i__2 = k;
	i__3 = k + 1;
	i__4 = k + 1;
	z__2.r = w[i__3].r * u[i__4].r - w[i__3].i * u[i__4].i, z__2.i = w[
		i__3].r * u[i__4].i + w[i__3].i * u[i__4].r;
	z__1.r = u[i__2].r - z__2.r, z__1.i = u[i__2].i - z__2.i;
	u[i__1].r = z__1.r, u[i__1].i = z__1.i;
    }

/* k */
    return 0;
} /* trirad_ */

