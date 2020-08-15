/* magfield.f -- translated by f2c (version 20000118).
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
    //doublereal awz[10001], awdz[10000], solz[10000], awerx[10000], awery[
	   // 10000], qfld[10001], fbess, magversion, unitlength, dqfx[10001], 
	   // dqfy[10001], awdx[10001], awdy[10001];
	doublereal *awz, *awdz, *solz, *awerx, *awery, *qfld, 
			fbess, magversion, unitlength, *dqfx, *dqfy, *awdx, *awdy;
    integer iran, nstepz, itap;
} wigcom_;

#define wigcom_1 wigcom_

/* Table of constant values */

static integer c_n11 = -11;
static integer c_n8 = -8;
static real c_b14 = (float)-1.;
static integer c__8 = 8;
static integer c__1 = 1;
static integer c_n9 = -9;
static integer c_n10 = -10;
static integer c_n12 = -12;
static integer c_n2 = -2;
static integer c_n17 = -17;

/* Subroutine */ int magfield_(xkw0, isup)
doublereal *xkw0;
integer *isup;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer i_dnnt();
    double sqrt(), pow_ri();

    /* Local variables */
    static integer ndrl, nlat, itmp;
    static doublereal norm, rold, atmp, corx[10000], cory[10000];
    extern /* Subroutine */ int last_(), magwrite_();
    extern integer printerr_();
    static integer i__, i1, i2, iskip, nwigz, iserr, nsecl;
    static doublereal rn;
    extern /* Subroutine */ int chk_maglen__();
    extern doublereal gasran_();
    static integer nfodol, ndl, nfl;
    extern /* Subroutine */ int magread_();
    static integer nsl, imz[11];
    static doublereal rnx, rny;
    extern doublereal faw0_(), ran1_();
    static integer n1st;

/*     ================================================================== */
/*     calculates the magnetic structure of the undulator */
/*     outputs: awz   - wiggler field on axis */
/*              awdz  - artificial k-parameter for matching of drifts */
/*              awerr - kick of momentum due to field errors */
/*              qx    - focusing strength in x (>0 -> focusing) */
/*              qy    - focusing strength in y */

/*     working arrays: p1 - readin quadrupole field */
/*                     k2gg - quad offset in x */
/*                     k3gg - quad offset in y */
/*                     k2pp - correction in x */
/*                     k3pp - correction in x */

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
/*     wiggler parameters */





/* ======================================================================= */
/* set default version number for external magnetic file. */

    wigcom_1.magversion = (float).1;
/* version of magnetic file */
    wigcom_1.unitlength = 0.;

/* ======================================================================= */

/*     clear arrays */

/* unit length in magfile header */
    for (i__ = 1; i__ <= 11; ++i__) {
	imz[i__ - 1] = 0;
    }
    for (i__ = 1; i__ <= 10000; ++i__) {
	wigcom_1.awz[i__] = (float)0.;
	wigcom_1.awdz[i__ - 1] = (float)0.;
	wigcom_1.qfld[i__] = (float)0.;
	wigcom_1.solz[i__ - 1] = (float)0.;
	wigcom_1.awerx[i__ - 1] = (float)0.;
	wigcom_1.awery[i__ - 1] = (float)0.;
	wigcom_1.dqfx[i__] = (float)0.;
	wigcom_1.dqfy[i__] = (float)0.;
	corx[i__ - 1] = (float)0.;
	cory[i__ - 1] = (float)0.;
	wigcom_1.awdx[i__] = (float)0.;
	wigcom_1.awdy[i__] = (float)0.;
    }

/*     read external file */

    if (inputcom_1.magin != 0) {
	magread_(corx, cory, imz, isup);
    }

/*     estimate length of the undulator */

    wigcom_1.nstepz = 0;
    nlat = 1;
    d__1 = (inputcom_1.fl + inputcom_1.dl + inputcom_1.drl * (float)2.) / 
	    inputcom_1.delz;
    nfodol = i_dnnt(&d__1);
/* steps of internal focusing la */
    if (nfodol <= 0) {
	nfodol = 1;
    }
    d__1 = inputcom_1.nwig / inputcom_1.delz;
    nwigz = i_dnnt(&d__1);
/* steps of internal undulator s */
    if (imz[3] > 0) {
	nfodol = imz[3];
    }
/* structure defined in input fi */
    if (imz[0] > 0) {
	inputcom_1.nsec = 1;
	nwigz = imz[0];
	wigcom_1.nstepz = imz[0];
/* aw0 in input file determines */
    }

    if (wigcom_1.nstepz <= 0) {
/* length yet not defined */
	if (nfodol >= nwigz) {
	    nlat = 1;
/* # of lattice per section */
	} else {
	    nlat = nwigz / nfodol;
	    if (nwigz % nfodol > 0) {
		++nlat;
	    }
	}
	wigcom_1.nstepz = inputcom_1.nsec * nlat * nfodol;
	//wigcom_1.nstepz = inputcom_1.nsec*nwigz + (inputcom_1.nsec - 1)*nfodol; //OC test

    }
    if (inputcom_1.zstop <= 0.) {
	inputcom_1.zstop = wigcom_1.nstepz * inputcom_1.delz * 
		inputcom_1.xlamd;
    }
/* bart      if (zstop.le.(nstepz*delz*xlamd)) nstepz=zstop/xlamd/delz */
    if (inputcom_1.zstop <= wigcom_1.nstepz * inputcom_1.delz * 
	    inputcom_1.xlamd) {
	d__1 = inputcom_1.zstop / inputcom_1.xlamd / inputcom_1.delz;
	wigcom_1.nstepz = i_dnnt(&d__1);
    }
    if (wigcom_1.nstepz > 10000) {
	i__ = printerr_(&c_n11, "undulator field", (ftnlen)15);
	wigcom_1.nstepz = 10000;
	inputcom_1.zstop = wigcom_1.nstepz * inputcom_1.delz * 
		inputcom_1.xlamd;
    }

    if (inputcom_1.magin != 0) {
	chk_maglen__(imz, &wigcom_1.nstepz);
    }

/*     check for inconsistent input */

/* check whether inp */
    itmp = 0;
    if (imz[3] <= 0) {
	if (inputcom_1.fl / inputcom_1.delz - (integer) (inputcom_1.fl / 
		inputcom_1.delz) * (float)1. > 1e-7) {
	    itmp = printerr_(&c_n8, "FL not a multiple of DELZ", (ftnlen)25);
	}
	if (inputcom_1.dl / inputcom_1.delz - (integer) (inputcom_1.dl / 
		inputcom_1.delz) * (float)1. > 1e-7) {
	    itmp = printerr_(&c_n8, "DL not a multiple of DELZ", (ftnlen)25);
	}
	if (inputcom_1.drl / inputcom_1.delz - (integer) (inputcom_1.drl / 
		inputcom_1.delz) * (float)1. > 1e-7) {
	    itmp = printerr_(&c_n8, "DRL not a multiple of DELZ", (ftnlen)26);
	}
	if (inputcom_1.f1st / inputcom_1.delz - (integer) (inputcom_1.f1st / 
		inputcom_1.delz) * (float)1. > 1e-7) {
	    itmp = printerr_(&c_n8, "F1ST not a multiple of DELZ", (ftnlen)27)
		    ;
	}
    }
    if (imz[8] <= 0) {
	if (inputcom_1.sl / inputcom_1.delz - (integer) (inputcom_1.sl / 
		inputcom_1.delz) * (float)1. > 1e-7) {
	    itmp = printerr_(&c_n8, "SL not a multiple of DELZ", (ftnlen)25);
	}
    }
    if (itmp < 0) {
	last_();
    }

/*     generatic magnetic field */

    d__1 = inputcom_1.fl / inputcom_1.delz;
    nfl = i_dnnt(&d__1);
    d__1 = inputcom_1.dl / inputcom_1.delz;
    ndl = i_dnnt(&d__1);
    d__1 = inputcom_1.drl / inputcom_1.delz;
    ndrl = i_dnnt(&d__1);
    d__1 = inputcom_1.f1st / inputcom_1.delz;
    n1st = i_dnnt(&d__1) + 1;
    d__1 = inputcom_1.sl / inputcom_1.delz;
    nsl = i_dnnt(&d__1);
    if (imz[8] > 0) {
	nsl = imz[8];
    }

/*     check for extreme long fodo cells */

    nsecl = nfodol * nlat;
    if (nfodol > nwigz << 1 && nfl == ndl) {
/* if fulfilled nl */
	nsecl /= 2;
/* nsecl is even because fodo cell is symmetric */
    }

    i__1 = inputcom_1.nsec;
    for (i1 = 1; i1 <= i__1; ++i1) {
	rnx = inputcom_1.awx * (ran1_(&inputcom_1.iseed) * 2. - 1.);
/* module offset in x */
	rny = inputcom_1.awy * (ran1_(&inputcom_1.iseed) * 2. - 1.);
/* module offset in y */
	i__2 = nsecl;
	for (i2 = 1; i2 <= i__2; ++i2) {
	    i__ = i2 + (i1 - 1) * nsecl;
	    if (i__ <= wigcom_1.nstepz) {

/*   if field not externally defined or if no taper defined */

/* check for array boundar */
		if (imz[0] <= 0) {
/* main undulator field not defined */
		    if (i2 <= nwigz) {
			wigcom_1.awz[i__] = inputcom_1.aw0;
			if (imz[9] <= 0) {
			    wigcom_1.awdx[i__] = rnx;
			}
/* offset in x */
			if (imz[10] <= 0) {
			    wigcom_1.awdy[i__] = rny;
			}
/* offset in y */
		    }

/* sven      undulator offset are only generated when the main undulator */
/* sven      is NOT defined in the external file. This is in difference to */
/* sven      quadrupole offset, where errors can be added later. Reason is that */
/* sven      genesis cannot distinguish modules from taper or errors! */
/* sven      hopefully the C++ version of genesis might solve this problem, because */
/* sven      magnetic elements are used in linked lists, thus distingishable! */

/* sven      Although accepted by Genesis, defining offsets in the external files */
/* sven      while the main undulator is generated internally is quite illogical!! */

		}

		if (imz[1] <= 0) {
/* drif section */
		    if (wigcom_1.awz[i__] < 1e-25) {
			wigcom_1.awdz[i__ - 1] = inputcom_1.awd;
		    }
		}

		if (n1st > nfodol) {
		    n1st = 1;
		}
/* check for overflow of c */
		if (imz[3] <= 0) {
/* quadrupole field */
		    if (n1st <= nfl) {
			wigcom_1.qfld[i__] = inputcom_1.quadf;
		    }
		    if (n1st > nfl + ndrl && n1st <= nfl + ndrl + ndl) {
			wigcom_1.qfld[i__] = inputcom_1.quadd;
		    }
		}
		++n1st;

		if (imz[8] <= 0) {
/* solenoid field */
		    if (i2 <= nsl) {
			wigcom_1.solz[i__ - 1] = inputcom_1.solen;
		    }
		}
	    }
	}
    }

/*     field taper: apply to existing lattice of awz(i) - either defined externally or aw0 */

    i__1 = wigcom_1.nstepz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	atmp = wigcom_1.awz[i__];
	wigcom_1.awz[i__] = faw0_(&i__, inputcom_1.wcoefz, &inputcom_1.delz, &
		inputcom_1.xlamd, &atmp, &wigcom_1.nstepz);
    }


/*     field errors: skip this part if iertyp=0 or if delaw=0 (no errors) */

/* bart      if((iertyp.eq.0).or.(delaw.lt.small)) then   ! no/ignore errors in external file */
    if (inputcom_1.iertyp != 0 && inputcom_1.delaw > 1e-7) {
/* no/ignore errors i */
	rn = (float)0.;
	iserr = 0;
	iskip = 1;
	if (inputcom_1.iertyp < 0) {
	    iskip = 2;
	}
/* correlated error */
	i__1 = wigcom_1.nstepz;
	i__2 = iskip;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	    if (abs(inputcom_1.iertyp) == 1) {
		rn = inputcom_1.delaw * (ran1_(&inputcom_1.iseed) * 2. - 1.) *
			 sqrt(3.);
/* uniform */
	    } else {
		rn = inputcom_1.delaw * gasran_(&inputcom_1.iseed) * sqrt(2.);
/* gaussian */
	    }
	    wigcom_1.awz[i__] *= rn + (float)1.;
	}
	wigcom_1.awz[0] = wigcom_1.awz[1];
	if (iskip == 2) {
/* no/ignore errors in external file */
	    i__2 = wigcom_1.nstepz - 2;
	    i__1 = iskip;
	    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
		if (wigcom_1.awz[i__] * wigcom_1.awz[i__ + 2] > 1e-7) {
		    wigcom_1.awz[i__ + 1] = (wigcom_1.awz[i__] + wigcom_1.awz[
			    i__ + 2]) / (float)2.;
		} else {
		    if (wigcom_1.awz[i__ + 1] > 1e-7) {
			wigcom_1.awz[i__ + 1] = wigcom_1.awz[i__];
		    }
		}
	    }
	}
/*  correlated errors */
    }
/* iertyp = 0 or delaw = 0. */
    wigcom_1.awz[0] = wigcom_1.awz[1];

/* 0 */
    norm = sqrt(2.);
/* normalization for the kick due to field */
    if (inputcom_1.iwityp != 0) {
	norm = (float)1.;
    }
/* bart      if (imz(2).le.0) then */
    i__1 = wigcom_1.nstepz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (wigcom_1.awz[i__ - 1] > 1e-25 && wigcom_1.awz[i__] > 1e-25) {
	    wigcom_1.awerx[i__ - 1] += pow_ri(&c_b14, &i__) * (wigcom_1.awz[
		    i__] - wigcom_1.awz[i__ - 1]) * norm / 3.14159265358979;
	}
    }
/* bart      endif */

    i__1 = wigcom_1.nstepz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wigcom_1.awery[i__ - 1] = wigcom_1.awerx[i__ - 1] * (doublereal) 
		inputcom_1.iwityp;
/* is nonzero for helical undulato */
    }

/*     scale quadrupole & solenoid field */
    i__1 = wigcom_1.nstepz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wigcom_1.qfld[i__] *= (float)586.;
	wigcom_1.solz[i__ - 1] *= (float)586.;
    }

/*     qfdx and qfdy can both be zero and the loop won't do anything. */
/*     the only check is, whether quadrupole errors are define in the */
/*     external input files. */

    rold = (float)0.;
    rnx = (float)0.;
    rny = (float)0.;
    if (imz[4] == 0 && imz[5] == 0) {
	i__1 = wigcom_1.nstepz;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if ((d__1 = wigcom_1.qfld[i__], abs(d__1)) > 1e-25 && abs(rold) < 
		    1e-25) {
/* b */
		if (inputcom_1.qfdx > 1e-7) {
		    rnx = inputcom_1.qfdx * (ran1_(&inputcom_1.iseed) * (
			    float)2. - 1);
		}
/* offset */
		if (inputcom_1.qfdy > 1e-7) {
		    rny = inputcom_1.qfdy * (ran1_(&inputcom_1.iseed) * (
			    float)2. - 1);
		}
	    }
	    rold = wigcom_1.qfld[i__];
	    if ((d__1 = wigcom_1.qfld[i__], abs(d__1)) > 1e-25) {
		wigcom_1.dqfx[i__] = rnx;
/* doen't matter whether */
		wigcom_1.dqfy[i__] = rny;
	    }
	}
    }

/*     write file with field description. */

    if (inputcom_1.magout != 0) {
	magwrite_(corx, cory, xkw0);
    }

/*     combine kicks off both planes */

    i__1 = wigcom_1.nstepz;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* sven         awerx(i)=awerx(i)+(dqfx(i)*qfld(i)+corx(i))/xkw0 !add kick of quad offset */
/* sven         awery(i)=awery(i)+(dqfy(i)*qfld(i)+cory(i))/xkw0 !to kick of field errors */
	wigcom_1.awerx[i__ - 1] += corx[i__ - 1] / *xkw0;
/* add kick to field errors */
	wigcom_1.awery[i__ - 1] += cory[i__ - 1] / *xkw0;
/* quadrupole offsets are now trea */
	wigcom_1.awdx[i__] *= *xkw0;
	wigcom_1.awdy[i__] *= *xkw0;
    }

/*     needed for output before 1st step */

    wigcom_1.qfld[0] = wigcom_1.qfld[1];
    wigcom_1.dqfx[0] = wigcom_1.dqfx[1];
    wigcom_1.dqfy[0] = wigcom_1.dqfy[1];
    wigcom_1.awdx[0] = wigcom_1.awdx[1];
    wigcom_1.awdy[0] = wigcom_1.awdy[1];

    return 0;
} /* magfield_ */



doublereal faw0_(i__, wcoefz, delz, xlamd, aw0, nstepz)
integer *i__;
doublereal *wcoefz, *delz, *xlamd, *aw0;
integer *nstepz;
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal z__, pctap, taplen, fac;

/*     ================================================================== */
/*     the wiggler amplitude profile along the undulator */
/*     the passed variable is the integration step istepz */
/*     this routine is called only once at the beginning */

/*     the tapering depends on wcoefz: */
/*             wcoefz(1): taper start location in z */
/*             wcoefz(2): gradient of field taper */
/*             wcoefz(3): the type of taper */
/*                             =0 -> no taper */
/*                             =1 -> linear taper */
/*                             =2 -> quadratic taper */
/*     ------------------------------------------------------------------ */



    /* Parameter adjustments */
    --wcoefz;

    /* Function Body */
    z__ = (doublereal) (*i__) * *delz * *xlamd - wcoefz[1];
/* position relative to taper sta */
    ret_val = *aw0;
    if (z__ <= (float)0. || wcoefz[3] == (float)0.) {
	return ret_val;
    }

/*    if external magnetic file includes taper, don't do anything */

/* before taper or no tap */
    pctap = wcoefz[2];
/* taper gradient */
    taplen = (doublereal) (*nstepz) * *delz * *xlamd - wcoefz[1];

/*     taper models */
/*     ------------------------------------------------------------------ */

/* taper length */
    if (wcoefz[3] < 1.1) {
	fac = 1. - pctap * z__ / taplen;
/* linear taper */
    } else {
	if (wcoefz[3] < 2.1) {
	    fac = 1. - pctap * z__ * z__ / taplen / taplen;
/* quadratic taper */
	} else {
	    fac = 1.;
	}
    }
    ret_val = fac * ret_val;

    return ret_val;
} /* faw0_ */



/* faw0 */
doublereal faw2_(i__, x, y)
integer *i__;
doublereal *x, *y;
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal xt, yt;

/*     ================================================================== */
/*     calculation of the square of the off-axis wiggler field at step i */
/*     the dependency on x and y is given by the wiggler type */
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

/*     ------------------------------------------------------------------ */
/*     wiggler parameters */






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



    xt = *x - wigcom_1.awdx[*i__];
/* include offset */
    yt = *y - wigcom_1.awdy[*i__];
/* in x and y */
    ret_val = wigcom_1.awz[*i__] * wigcom_1.awz[*i__] * (inputcom_1.xkx * xt *
	     xt + 1. + inputcom_1.xky * yt * yt);

    return ret_val;
} /* faw2_ */



/* faw2 */
/* Subroutine */ int magwrite_(corx, cory, xkw0)
doublereal *corx, *cory, *xkw0;
{
    /* Format strings */
    static char fmt_50[] = "(a)";
    static char fmt_40[] = "(a,1x,1f4.2,1x,a)";
    static char fmt_45[] = "(a,1x,1f7.5,1x,a)";
    static char fmt_60[] = "(a3,1x,1pe14.4,1x,i4,1x,i4)";
    static char fmt_30[] = "(a3,1x,1pe14.4,1x,1pe14.4,1x,i4)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy();
    integer s_wsfe(), do_fio(), e_wsfe(), f_clos();

    /* Local variables */
    static doublereal rold, rcur;
    extern integer opentextfile_();
    static integer i__, k, nmout, ic, nr;
    static char cid[3*11];

    /* Fortran I/O blocks */
    static cilist io___33 = { 0, 0, 0, fmt_50, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_45, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_60, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_30, 0 };


/*     =================================================================== */
/*     output of the used magentic field */
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

/*     ------------------------------------------------------------------ */
/*     wiggler parameters */






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



    /* Parameter adjustments */
    --cory;
    --corx;

    /* Function Body */
    if (inputcom_1.magin == 0) {
	wigcom_1.magversion = (float)1.;
    }
/* write in new format unless spec. */
    s_copy(cid, "AW ", (ftnlen)3, (ftnlen)3);
    s_copy(cid + 3, "DP ", (ftnlen)3, (ftnlen)3);
/* not used in output -> can be recalculated from aw */
    s_copy(cid + 6, "QF ", (ftnlen)3, (ftnlen)3);
    s_copy(cid + 9, "QX ", (ftnlen)3, (ftnlen)3);
    s_copy(cid + 12, "QY ", (ftnlen)3, (ftnlen)3);
    s_copy(cid + 15, "AD ", (ftnlen)3, (ftnlen)3);
    s_copy(cid + 18, "SL ", (ftnlen)3, (ftnlen)3);
    s_copy(cid + 21, "CX ", (ftnlen)3, (ftnlen)3);
    s_copy(cid + 24, "CY ", (ftnlen)3, (ftnlen)3);
    s_copy(cid + 27, "AX", (ftnlen)3, (ftnlen)2);
    s_copy(cid + 30, "AY", (ftnlen)3, (ftnlen)2);

    nmout = opentextfile_(inputcom_1.magoutfile, "unknown", &c__8, (ftnlen)30,
	     (ftnlen)7);
    if (nmout < 0) {
	return 0;
    }

/*  write header of the file */

    if (wigcom_1.magversion > (float).12) {
	io___33.ciunit = nmout;
	s_wsfe(&io___33);
	do_fio(&c__1, "# header is included", (ftnlen)20);
	e_wsfe();
	io___34.ciunit = nmout;
	s_wsfe(&io___34);
	do_fio(&c__1, "? VERSION=", (ftnlen)10);
	do_fio(&c__1, (char *)&wigcom_1.magversion, (ftnlen)sizeof(doublereal)
		);
	do_fio(&c__1, " including new format", (ftnlen)21);
	e_wsfe();
	io___35.ciunit = nmout;
	s_wsfe(&io___35);
	do_fio(&c__1, "? UNITLENGTH=", (ftnlen)13);
	d__1 = inputcom_1.xlamd * inputcom_1.delz;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, ":unit length in header", (ftnlen)22);
	e_wsfe();
    }

    for (k = 1; k <= 11; ++k) {
	i__ = 1;
	ic = 1;
	nr = 0;
	rold = wigcom_1.awz[i__];
	if (k == 2) {
	    rold = 0.;
	}
	if (k == 3) {
	    rold = wigcom_1.qfld[i__] / (float)586.;
	}
	if (k == 4) {
	    rold = wigcom_1.dqfx[i__];
	}
	if (k == 5) {
	    rold = wigcom_1.dqfy[i__];
	}
	if (k == 6) {
	    rold = wigcom_1.awdz[i__ - 1];
	}
	if (k == 7) {
	    rold = wigcom_1.solz[i__ - 1] / (float)586.;
	}
	if (k == 8) {
	    rold = corx[i__];
	}
	if (k == 9) {
	    rold = cory[i__];
	}
	if (k == 10) {
	    rold = wigcom_1.awdx[i__];
	}
	if (k == 11) {
	    rold = wigcom_1.awdy[i__];
	}
L1:
	++i__;
	++ic;
	rcur = wigcom_1.awz[i__];
	if (k == 2) {
	    rcur = 0.;
	}
	if (k == 3) {
	    rcur = wigcom_1.qfld[i__] / (float)586.;
	}
	if (k == 4) {
	    rcur = wigcom_1.dqfx[i__];
	}
	if (k == 5) {
	    rcur = wigcom_1.dqfy[i__];
	}
	if (k == 6) {
	    rcur = wigcom_1.awdz[i__ - 1];
	}
	if (k == 7) {
	    rcur = wigcom_1.solz[i__ - 1] / (float)586.;
	}
	if (k == 8) {
	    rcur = corx[i__];
	}
	if (k == 9) {
	    rcur = cory[i__];
	}
	if (k == 10) {
	    rcur = wigcom_1.awdx[i__];
	}
	if (k == 11) {
	    rcur = wigcom_1.awdy[i__];
	}
	if (i__ > wigcom_1.nstepz) {
/* bart            if ((ic.gt.2).and.(k.ne.2)) then */
	    if (ic >= 2 && k != 2) {
		if (wigcom_1.magversion > (float).12) {
		    io___42.ciunit = nmout;
		    s_wsfe(&io___42);
		    do_fio(&c__1, cid + (k - 1) * 3, (ftnlen)3);
		    do_fio(&c__1, (char *)&rold, (ftnlen)sizeof(doublereal));
		    i__1 = ic - 1;
		    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&nr, (ftnlen)sizeof(integer));
		    e_wsfe();
		    nr = 0;
		} else {
		    io___43.ciunit = nmout;
		    s_wsfe(&io___43);
		    do_fio(&c__1, cid + (k - 1) * 3, (ftnlen)3);
		    do_fio(&c__1, (char *)&rold, (ftnlen)sizeof(doublereal));
		    d__1 = inputcom_1.delz * inputcom_1.xlamd;
		    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		    i__1 = ic - 1;
		    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		    e_wsfe();
		}
	    }
	    goto L2;
	}
	if ((d__1 = rcur - rold, abs(d__1)) > 1e-25 && k != 2) {
	    if (wigcom_1.magversion > (float).12) {
		if (abs(rold) < 1e-25) {
		    nr = ic - 1;
		} else {
		    io___44.ciunit = nmout;
		    s_wsfe(&io___44);
		    do_fio(&c__1, cid + (k - 1) * 3, (ftnlen)3);
		    do_fio(&c__1, (char *)&rold, (ftnlen)sizeof(doublereal));
		    i__1 = ic - 1;
		    do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		    do_fio(&c__1, (char *)&nr, (ftnlen)sizeof(integer));
		    e_wsfe();
		    nr = 0;
		}
	    } else {
		io___45.ciunit = nmout;
		s_wsfe(&io___45);
		do_fio(&c__1, cid + (k - 1) * 3, (ftnlen)3);
		do_fio(&c__1, (char *)&rold, (ftnlen)sizeof(doublereal));
		d__1 = inputcom_1.delz * inputcom_1.xlamd;
		do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
		i__1 = ic - 1;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    rold = rcur;
	    ic = 1;
	}
	goto L1;
L2:
	;
    }

    cl__1.cerr = 0;
    cl__1.cunit = nmout;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;

/*     format statements */

} /* magwrite_ */




/* Subroutine */ int magread_(corx, cory, imz, isup)
doublereal *corx, *cory;
integer *imz, *isup;
{
    /* Format strings */
    static char fmt_1000[] = "(a)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy();
    integer i_len(), s_rsfe(), do_fio(), e_rsfe(), i_indx(), f_clos(), s_cmp()
	    , i_dnnt();

    /* Local variables */
    static char line[255], cmagtype[30*11];
    static integer imzb[11], idum, loop, nmin, ncol, ierr, nlen;
    extern /* Subroutine */ int last_(), getfirstchar_();
    extern integer printerr_(), opentextfile_();
    static integer i__, j, k;
    static doublereal r1, r2, r3;
    static integer ntemp, nloop, ninfo;
    extern /* Subroutine */ int closefile_();
    extern integer extractnumber_();
    static integer nr;
    static doublereal values[4];
    extern /* Subroutine */ int getmagfileinfo_();
    static integer nr2;
    extern integer extractval_();
    static char cin[30];
    static doublereal val;
    static integer idx, loopcnt;
    extern /* Subroutine */ int touppercase_();
    static integer int_version__;

    /* Fortran I/O blocks */
    static cilist io___55 = { 1, 0, 1, fmt_1000, 0 };


/*     ======================================================================== */
/*     input of the magnetic field. */
/*     only those parameters are replace which are included in the list */
/*     format : indicator - strength - length - number */
/*     inticator: aw - wiggler field */
/*                ad - phase matching between drift section */
/*                dp - kick due to field errors */
/*                qf - quadrupole field */
/*                qx - quadrupole offset in x */
/*                qy - quadrupole offset in y */
/*                cx - corrector strength in x */
/*                cy - corrector strength in y */
/*                sl - solenoid field */
/* 				 ax - undulator offset in x */
/* 				 ay - undulator offset in y */

/*     note: dp,co qx are combined into awerx. output is only dp */
/*     -------------------------------------------------------------------------- */





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

/*     ------------------------------------------------------------------ */
/*     wiggler parameters */






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



    /* Parameter adjustments */
    --imz;
    --cory;
    --corx;

    /* Function Body */
    if (*isup != 0) {
	if (inputcom_1.wcoefz[1] != 0.) {
	    ierr = printerr_(&c_n9, "Taper defined in namelist", (ftnlen)25);
	}
	if (inputcom_1.iertyp != 0 && abs(inputcom_1.delaw) > 1e-7) {
	    ierr = printerr_(&c_n9, "Wiggler errors defined in namelist", (
		    ftnlen)34);
	}
	if (abs(inputcom_1.qfdx) > 1e-7) {
	    ierr = printerr_(&c_n9, "Random quad offset errors in x defined \
in namelist", (ftnlen)50);
	}
	if (abs(inputcom_1.qfdy) > 1e-7) {
	    ierr = printerr_(&c_n9, "Random quad offset errors in y defined \
in namelist", (ftnlen)50);
	}
	if (abs(inputcom_1.awx) > 1e-7) {
	    ierr = printerr_(&c_n9, "Random undulator offset errors in x def\
ined in namelist", (ftnlen)55);
	}
	if (abs(inputcom_1.awy) > 1e-7) {
	    ierr = printerr_(&c_n9, "Random undulator offset errors in y def\
ined in namelist", (ftnlen)55);
	}
    }

    s_copy(cmagtype, "undulator field", (ftnlen)30, (ftnlen)15);
    s_copy(cmagtype + 30, "drift section", (ftnlen)30, (ftnlen)13);
    s_copy(cmagtype + 60, "field errors", (ftnlen)30, (ftnlen)12);
/* not used anymore */
    s_copy(cmagtype + 90, "quadrupole field", (ftnlen)30, (ftnlen)16);
    s_copy(cmagtype + 120, "quadrupole offset in x", (ftnlen)30, (ftnlen)22);
    s_copy(cmagtype + 150, "quadrupole offset in y", (ftnlen)30, (ftnlen)22);
    s_copy(cmagtype + 180, "orbit correction in x", (ftnlen)30, (ftnlen)21);
    s_copy(cmagtype + 210, "orbit correction in y", (ftnlen)30, (ftnlen)21);
    s_copy(cmagtype + 240, "solenoid field", (ftnlen)30, (ftnlen)14);
    s_copy(cmagtype + 270, "undulator offset in x", (ftnlen)30, (ftnlen)21);
    s_copy(cmagtype + 300, "undulator offset in y", (ftnlen)30, (ftnlen)21);
    wigcom_1.unitlength = (float)0.;

/* unit length has to be checked for new version */
    nlen = i_len(line, (ftnlen)255);
/* line size can be easily changed */
    loop = 0;
    loopcnt = 0;
    int_version__ = 2;
/* check if version set to 0.1 */
    ncol = 3;

/*     open file */

/* default older version to read 3 numbers after th */
    nmin = opentextfile_(inputcom_1.maginfile, "old", &c__8, (ftnlen)30, (
	    ftnlen)3);
    if (nmin < 0) {
	return 0;
    }

L1:
    io___55.ciunit = nmin;
    i__1 = s_rsfe(&io___55);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_fio(&c__1, line, (ftnlen)255);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = e_rsfe();
L100001:
    if (i__1 < 0) {
	goto L2;
    }
    if (i__1 > 0) {
	goto L100;
    }

/* read line */
    touppercase_(line, (ftnlen)255);
/* convert to upper case */
    getfirstchar_(line, &idx, (ftnlen)255);

/* get index of first non-spac */
    if (idx == 0) {
	goto L1;
    }
/* empty line */
    if (*(unsigned char *)&line[idx - 1] == '#') {
	goto L1;
    }
/* comment line */
    if (*(unsigned char *)&line[idx - 1] == '?') {
/* information line */
	i__1 = idx;
	getmagfileinfo_(line + i__1, nlen - i__1);
	if (wigcom_1.magversion < (float).11) {
	    ncol = 3;
	    int_version__ = 1;
	} else {
	    ncol = 3;
/* read 4 number at version >= */
	}
	goto L1;
    }
    if (*(unsigned char *)&line[idx - 1] == '!') {
/* start/end loop structure */
	i__1 = idx;
	s_copy(line, line + i__1, (ftnlen)255, nlen - i__1);
	getfirstchar_(line, &idx, (ftnlen)255);
	if (idx == 0) {
	    i__ = printerr_(&c_n9, "Empty line behind \"!\": ignored", (
		    ftnlen)30);
	    goto L1;
	}

/* check for start of loop: set loop counter */

	nloop = i_indx(line, "LOOP", (ftnlen)255, (ftnlen)4);
	if (nloop == 0) {
	    ierr = printerr_(&c_n9, "Undefined command in maginfile", (ftnlen)
		    30);
	    goto L1;
	}
	if (nloop == idx) {
	    if (loop == 1) {
/* start next loop before ending previous */
		nloop = printerr_(&c_n9, "Illegal ending of loop", (ftnlen)22)
			;
/* for */
		cl__1.cerr = 0;
		cl__1.cunit = nmin;
		cl__1.csta = 0;
		f_clos(&cl__1);
		last_();
	    }
	    i__1 = idx + 3;
	    ierr = extractnumber_(line + i__1, &val, nlen - i__1);
	    if (ierr < 0) {
		i__ = printerr_(&c_n9, "Undefined loop argument", (ftnlen)23);
		goto L1;
	    } else {
		loopcnt = (integer) val;
	    }
	    if (loopcnt > 1) {
		loop = 1;
		for (k = 1; k <= 9; ++k) {
		    imzb[k - 1] = imz[k];
		}
	    }
	    goto L1;
	}

/* check for loop ending: reset loopcounter and copy loop */

	nloop = i_indx(line, "ENDLOOP", (ftnlen)255, (ftnlen)7);
	if (nloop == idx && loop == 1) {
/* a loop is ende */
	    i__1 = loopcnt;
	    for (j = 2; j <= i__1; ++j) {
/* copy filed loopcnt-1 times */
		for (k = 1; k <= 11; ++k) {
		    ntemp = imz[k] - imzb[k - 1];
		    if (imz[1] + ntemp > 10000) {
			idum = printerr_(&c_n11, cmagtype + (k - 1) * 30, (
				ftnlen)30);
		    } else {
			i__2 = ntemp;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    if (k == 1) {
				wigcom_1.awz[imz[1] + i__] = wigcom_1.awz[
					imzb[0] + i__];
			    }
			    if (k == 2) {
				wigcom_1.awdz[imz[2] + i__ - 1] = 
					wigcom_1.awdz[imzb[1] + i__ - 1];
			    }
			    if (k == 4) {
				wigcom_1.qfld[imz[4] + i__] = wigcom_1.qfld[
					imzb[3] + i__];
			    }
			    if (k == 5) {
				wigcom_1.dqfx[imz[5] + i__] = wigcom_1.dqfx[
					imzb[4] + i__];
			    }
			    if (k == 6) {
				wigcom_1.dqfy[imz[6] + i__] = wigcom_1.dqfy[
					imzb[5] + i__];
			    }
			    if (k == 7) {
				corx[imz[7] + i__] = corx[imzb[6] + i__];
			    }
			    if (k == 8) {
				cory[imz[8] + i__] = cory[imzb[7] + i__];
			    }
			    if (k == 9) {
				wigcom_1.solz[imz[9] + i__ - 1] = 
					wigcom_1.solz[imzb[8] + i__ - 1];
			    }
			    if (k == 10) {
				wigcom_1.awdx[imz[10] + i__] = wigcom_1.awdx[
					imzb[9] + i__];
			    }
			    if (k == 11) {
				wigcom_1.awdy[imz[11] + i__] = wigcom_1.awdy[
					imzb[10] + i__];
			    }
			}
			imz[k] += ntemp;
			imzb[k - 1] += ntemp;
		    }
		}
	    }
	    loop = 0;
/* finish the round of copying */
	    loopcnt = 0;
	    goto L1;
	}
    }

/*     process input line */

    s_copy(cin, line + (idx - 1), (ftnlen)30, (ftnlen)2);
/* get identifier */
    k = -1;
    if (s_cmp(cin, "AW", (ftnlen)2, (ftnlen)2) == 0) {
	k = 1;
    }
/* main wiggler field */
    if (s_cmp(cin, "AD", (ftnlen)2, (ftnlen)2) == 0) {
	k = 2;
    }
/* art. wiggler field between section */
    if (s_cmp(cin, "QF", (ftnlen)2, (ftnlen)2) == 0) {
	k = 4;
    }
/* bart      if (cin(1:2).eq.'QD') k=5 !quadrupole offset in x */
/* quadrupole field strength */
    if (s_cmp(cin, "QX", (ftnlen)2, (ftnlen)2) == 0) {
	k = 5;
    }
/* quadrupole offset in x */
    if (s_cmp(cin, "QY", (ftnlen)2, (ftnlen)2) == 0) {
	k = 6;
    }
/* quadrupole offset in y */
    if (s_cmp(cin, "CX", (ftnlen)2, (ftnlen)2) == 0) {
	k = 7;
    }
/* correcter strength in x */
    if (s_cmp(cin, "CY", (ftnlen)2, (ftnlen)2) == 0) {
	k = 8;
    }
/* corrector strength in y */
    if (s_cmp(cin, "SL", (ftnlen)2, (ftnlen)2) == 0) {
	k = 9;
    }
/* solenoid strength */
    if (s_cmp(cin, "AX", (ftnlen)2, (ftnlen)2) == 0) {
	k = 10;
    }
/* undulator offset in x */
    if (s_cmp(cin, "AY", (ftnlen)2, (ftnlen)2) == 0) {
	k = 11;
    }
/* undulator offset in y */
    if (k < 0) {
	ierr = printerr_(&c_n10, cin, (ftnlen)2);
	goto L1;
    }

    i__1 = idx + 1;
    s_copy(line, line + i__1, (ftnlen)255, nlen - i__1);
    getfirstchar_(line, &idx, (ftnlen)255);
/* eliminate spaces in the beginning */
    if (idx > 0) {
	s_copy(line, line + (idx - 1), (ftnlen)255, nlen - (idx - 1));
    }
    ierr = extractval_(line, values, &ncol, (ftnlen)255);
    if (ierr < 0) {
	ierr = printerr_(&c_n10, line, (ftnlen)255);
	goto L1;
    }
    r1 = values[0];

/*     filling the arrays */

/* strength/offset */
    if (wigcom_1.magversion >= (float).12) {
	nr = i_dnnt(&values[1]);
/* length in unit length */
	nr2 = i_dnnt(&values[2]);
/* separation to prev. element */
	if (nr < 0 || nr2 < 0) {
	    idum = printerr_(&c_n8, "Negative element length/distance in mag\
infile", (ftnlen)45);
	    last_();
	}
	r2 = wigcom_1.unitlength;
/* unit length */
	r3 = r2 * nr2;
	d__1 = r3 / inputcom_1.delz / inputcom_1.xlamd;
	nr2 = i_dnnt(&d__1);
	if ((d__1 = nr2 - r3 / inputcom_1.delz / inputcom_1.xlamd, abs(d__1)) 
		> 1e-7) {
	    i__ = printerr_(&c_n12, line, (ftnlen)255);
/* different error needed ???????? */
	}
	if (imz[k] + nr2 > 10000) {
	    idum = printerr_(&c_n11, cmagtype + (k - 1) * 30, (ftnlen)30);
	} else {
	    imz[k] += nr2;
	}
    } else {
	r2 = values[1];
/* unit length */
	nr = i_dnnt(&values[2]);
/* length in unit length */
	if (nr < 0 || r2 < -1e-7) {
	    idum = printerr_(&c_n8, "negative element length/distance in mag\
infile", (ftnlen)45);
	    last_();
	}
    }

    r2 *= nr;
/* full length of this section */
    d__1 = r2 / inputcom_1.delz / inputcom_1.xlamd;
    nr = i_dnnt(&d__1);
/* # int. steps for this section */
    if ((d__1 = nr - r2 / inputcom_1.delz / inputcom_1.xlamd, abs(d__1)) > 
	    1e-7) {
	i__ = printerr_(&c_n12, line, (ftnlen)255);
    }
    if (imz[k] + nr > 10000) {
	idum = printerr_(&c_n11, cmagtype + (k - 1) * 30, (ftnlen)30);
    } else {
	i__1 = nr;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (k == 1) {
		wigcom_1.awz[imz[1] + i__] = r1;
	    }
	    if (k == 2) {
		wigcom_1.awdz[imz[2] + i__ - 1] = r1;
	    }
	    if (k == 4) {
		wigcom_1.qfld[imz[4] + i__] = r1;
	    }
	    if (k == 5) {
		wigcom_1.dqfx[imz[5] + i__] = r1;
	    }
	    if (k == 6) {
		wigcom_1.dqfy[imz[6] + i__] = r1;
	    }
	    if (k == 7) {
		corx[imz[7] + i__] = r1;
	    }
	    if (k == 8) {
		cory[imz[8] + i__] = r1;
	    }
	    if (k == 9) {
		wigcom_1.solz[imz[9] + i__ - 1] = r1;
	    }
	    if (k == 10) {
		wigcom_1.awdx[imz[10] + i__] = r1;
	    }
	    if (k == 11) {
		wigcom_1.awdy[imz[11] + i__] = r1;
	    }
	}
	imz[k] += nr;
    }
    goto L1;

/*     final processing + closing files */

/* read new line */
L2:
    if (int_version__ == 2) {
	wigcom_1.magversion = (float)1.;
    }
    if (loop == 1) {
	ninfo = printerr_(&c_n9, "\"LOOP\" not terminated: no loop", (ftnlen)
		30);
    }
    closefile_(&nmin);

    return 0;

/*     error */

L100:
    idum = printerr_(&c_n2, inputcom_1.maginfile, (ftnlen)30);
    goto L1;

/*     format statements */


} /* magread_ */


/* Subroutine */ int chk_maglen__(imz, nstepz)
integer *imz, *nstepz;
{
    /* Builtin functions */
    /* Subroutine */ int s_copy();

    /* Local variables */
    static char cmagtype[30*11];
    extern integer printerr_();
    static integer i__, j;

/*     =================================================================== */
/*     checks whether the user supplied file for the description of the */
/*     magnetic fields is incomplete */
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


    /* Parameter adjustments */
    --imz;

    /* Function Body */
    s_copy(cmagtype, "undulator field", (ftnlen)30, (ftnlen)15);
    s_copy(cmagtype + 30, "drift section", (ftnlen)30, (ftnlen)13);
    s_copy(cmagtype + 90, "quadrupole field", (ftnlen)30, (ftnlen)16);
    s_copy(cmagtype + 120, "quadrupole offset in x", (ftnlen)30, (ftnlen)22);
    s_copy(cmagtype + 150, "quadrupole offset in y", (ftnlen)30, (ftnlen)22);
    s_copy(cmagtype + 180, "orbit correction in x", (ftnlen)30, (ftnlen)21);
    s_copy(cmagtype + 210, "orbit correction in y", (ftnlen)30, (ftnlen)21);
    s_copy(cmagtype + 240, "solenoid field", (ftnlen)30, (ftnlen)14);
    s_copy(cmagtype + 270, "undulator offset in x", (ftnlen)30, (ftnlen)21);
    s_copy(cmagtype + 300, "undulator offset in y", (ftnlen)30, (ftnlen)21);

    for (i__ = 1; i__ <= 11; ++i__) {
	if (imz[i__] > 1 && imz[i__] < *nstepz) {
	    j = printerr_(&c_n17, cmagtype + (i__ - 1) * 30, (ftnlen)30);
	}
    }
    return 0;
} /* chk_maglen__ */



/* of chk_maglen */
/* Subroutine */ int getmagfileinfo_(line, line_len)
char *line;
ftnlen line_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(), i_indx();

    /* Local variables */
    static integer ierr;
    extern /* Subroutine */ int getfirstchar_();
    extern integer printerr_();
    static integer i__, n;
    extern integer extractnumber_();
    static doublereal val;
    static integer idx;

/*     ================================================================= */
/*     extract information from beamfile */
/*     ----------------------------------------------------------------- */





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


/*     ------------------------------------------------------------------ */
/*     wiggler parameters */





    n = i_len(line, line_len);
    getfirstchar_(line, &idx, line_len);

/*     version number */

/* get first character should be identi */
    i__ = i_indx(line, "VERSION", line_len, (ftnlen)7);
/* check for version number */
    if (i__ > 0 && i__ == idx) {
	i__1 = i__ + 6;
	ierr = extractnumber_(line + i__1, &val, n - i__1);
	if (ierr < 0) {
	    i__ = printerr_(&c_n9, "Unrecognized information line in maginfi\
le", (ftnlen)42);
	} else {
	    wigcom_1.magversion = val;
	}
	return 0;
    }


    i__ = i_indx(line, "UNITLENGTH", line_len, (ftnlen)10);
/* check for unit length */
    if (i__ > 0 && i__ == idx) {
	i__1 = i__ + 6;
	ierr = extractnumber_(line + i__1, &val, n - i__1);
	if (ierr < 0) {
	    i__ = printerr_(&c_n9, "Unrecognized information line in maginfi\
le", (ftnlen)42);
	} else {
	    wigcom_1.unitlength = val;
	}
	return 0;
    }

/*     unrecognized */

    i__ = printerr_(&c_n9, "Unrecognized information line in beamfile", (
	    ftnlen)41);
    return 0;
} /* getmagfileinfo_ */

