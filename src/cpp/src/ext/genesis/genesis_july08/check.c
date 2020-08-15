/* check.f -- translated by f2c (version 20000118).
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
    doublereal xkw0, xkper0, sval, svalout, gamma0_in__;
    integer npart0, inorun;
} simcom_;

#define simcom_1 simcom_

Extern struct {
    integer mpi_id__, mpi_err__, mpi_size__, mpi_loop__, nfldmpi, nparmpi, nfldhmpi[6];
} mpicom_;

#define mpicom_1 mpicom_

Extern struct {
    //doublecomplex crfield[476847], crsource[68121], crmatc[1827], cstep[7], crhm[68121], cwet[1827], cbet[1827];
    doublecomplex *crfield, *crsource, *crmatc, cstep[7], *crhm, *cwet, *cbet;
    doublereal dxy, xks, radphase, besselcoupling[7];
    integer nhloop, hloop[7];
} cartcom_;

#define cartcom_1 cartcom_

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
    doublereal distversion, distrev;
    integer iout[39], nout, nfld, npar, ndump, firstout, nfin, irecpar, 
	    irecfld, kout, ndmp2, npin, nlog, ndis, ncoldis, iconv2t, iconv2g,
	     icolpar[10], ndistsize, iconv2px, iconv2py, nprobe, ftype, 
	    ftdist, ftpart, ftfield, ndumph[6], nfldh[6];
} iocom_;

#define iocom_1 iocom_

/* Table of constant values */

static integer c_n9 = -9;
static integer c_n14 = -14;
static integer c_n15 = -15;
static integer c_n5 = -5;
static integer c_n8 = -8;
static integer c_n21 = -21;
static integer c__1 = 1;
static integer c_n18 = -18;

/* Subroutine */ int chk_input__()
{
    /* Format strings */
    static char fmt_100[] = "(\002Please enter magnetic input file name\002)";
    static char fmt_200[] = "(a30)";
    static char fmt_110[] = "(\002Please enter magnetic output file name\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt();
    integer i_indx(), s_wsfe(), e_wsfe(), s_rsfe(), do_fio(), e_rsfe();

    /* Local variables */
    extern integer chk_scan__();
    extern /* Subroutine */ int last_();
    extern integer printerr_();
    static integer i__, ix;
    static doublereal rw0;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, fmt_100, 0 };
    static cilist io___4 = { 0, 5, 0, fmt_200, 0 };
    static cilist io___5 = { 0, 6, 0, fmt_110, 0 };
    static cilist io___6 = { 0, 5, 0, fmt_200, 0 };


/*     ============================================================ */
/*     check for inconsistencies of correlated input parameter */
/*     such as magin and maginfile. guerantee compability for */
/*     older versions of genesis 1.3 */
/*     ------------------------------------------------------------ */





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


/*     simulation control and normalisation parameter */


/*     ------------------------------------------------------------------ */
/*     common block for mpi implementation for genesis */





/*     -------------------------------------------------------------------- */
/*     cartesian mesh */




/*     set seeds for random number generator to negative value */

    inputcom_1.iseed = -abs(inputcom_1.iseed);
    inputcom_1.ipseed = -abs(inputcom_1.ipseed) - mpicom_1.mpi_id__;
    inputcom_1.quadf = abs(inputcom_1.quadf);
    inputcom_1.quadd = -abs(inputcom_1.quadd);

/*     save the input value for gamma0 for diagnostic output and advanced scan */

    simcom_1.gamma0_in__ = inputcom_1.gamma0;

/*     set flags to either 0 or 1 + adjust values */

    if (inputcom_1.isravg != 0) {
	inputcom_1.isravg = 1;
    }
    if (inputcom_1.isrsig != 0) {
	inputcom_1.isrsig = 1;
    }
    if (inputcom_1.lbc != 0) {
	inputcom_1.lbc = 1;
    }
    if (inputcom_1.iorb != 0) {
	inputcom_1.iorb = 1;
    }
    if (inputcom_1.magin != 0) {
	inputcom_1.magin = 1;
    }
    simcom_1.inorun = 0;
/* do not stop after initialization */
    if (inputcom_1.magout < 0) {
	simcom_1.inorun = 1;
    }
/* enforce termination after initializa */
    if (inputcom_1.magout != 0) {
	inputcom_1.magout = 1;
    }
    if (inputcom_1.idump != 0) {
	inputcom_1.idump = 1;
    }
    if (inputcom_1.idmppar != 0) {
	inputcom_1.idmppar = 1;
    }
    if (inputcom_1.idmpfld != 0) {
	inputcom_1.idmpfld = 1;
    }
    if (inputcom_1.iotail != 0) {
	inputcom_1.iotail = 1;
    }
    if (inputcom_1.ilog != 0) {
	inputcom_1.ilog = 1;
    }
    if (inputcom_1.iall != 0) {
	inputcom_1.iall = 1;
    }
    if (inputcom_1.itdp != 0) {
	inputcom_1.itdp = 1;
    }
    if (inputcom_1.wcoefz[1] > 1.) {
	inputcom_1.wcoefz[1] = 1.;
    }
    if (inputcom_1.iwityp != 0) {
	inputcom_1.iwityp = 1;
    }
    if (inputcom_1.delaw < 1e-7) {
	inputcom_1.iertyp = 0;
    }
    if (inputcom_1.iertyp == 0) {
	inputcom_1.delaw = 0.;
    }
    if (inputcom_1.ffspec != 0) {
	inputcom_1.ffspec /= abs(inputcom_1.ffspec);
    }
    if (inputcom_1.isntyp != 0) {
	inputcom_1.isntyp = 1;
    }
    if (inputcom_1.dgrid <= 1e-7 && inputcom_1.zrayl > (float)0.) {
/* grid size determine */
/* Computing 2nd power */
	d__1 = inputcom_1.zwaist / inputcom_1.zrayl;
	rw0 = sqrt(inputcom_1.zrayl * inputcom_1.xlamds / 3.14159265358979 * (
		d__1 * d__1 + 1.));
/* Computing 2nd power */
	d__1 = inputcom_1.rxbeam;
/* Computing 2nd power */
	d__2 = inputcom_1.rybeam;
	inputcom_1.dgrid = inputcom_1.rmax0 * (rw0 + sqrt(d__1 * d__1 + d__2 *
		 d__2)) / 2.;
    }

/*     define harmonic content */

    cartcom_1.nhloop = 1;
    cartcom_1.hloop[0] = 1;
    if (inputcom_1.nharm > 1) {
/* changed to avoid nharm=0 settin */
	if (inputcom_1.iallharm != 0) {
/* denotes which harmonics shall */
	    cartcom_1.nhloop = inputcom_1.nharm;
	    i__1 = inputcom_1.nharm;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* fill array with harmonic numbers */
		cartcom_1.hloop[i__ - 1] = i__;
	    }
	} else {
	    cartcom_1.nhloop = 2;
/* else only fundamental + one ha */
	    cartcom_1.hloop[1] = inputcom_1.nharm;
	}
    } else {
	inputcom_1.pradh0 = (float)0.;
    }

/*     check for version specific changes */

    if (inputcom_1.idmppar == 0) {
	inputcom_1.idmppar = inputcom_1.idump;
    }

    if (inputcom_1.version < (float)1.) {
	inputcom_1.ffspec = 0;
/* phase in near field */
	inputcom_1.isntyp = 1;
/* shotnoise with Penman algorithm */
	simcom_1.inorun = 0;
    } else {
	if (inputcom_1.idmpfld == 0) {
	    inputcom_1.idmpfld = inputcom_1.idump;
	}
    }

/*     check case if time dependent code is selected */

    if (inputcom_1.itdp == 0) {
	inputcom_1.nslice = 1;
/* one slice */
	inputcom_1.zsep = inputcom_1.delz;
/* one bucket */
	inputcom_1.ntail = 0;
/* middle of the beam */
	inputcom_1.curlen = -1.;
/* step profile in z */
	inputcom_1.shotnoise = 0.;
/* no phase fluctuation */
	inputcom_1.iall = 1;
/* reference loading for scan */
	inputcom_1.ishsty = 1;
/* make sure for output */
	inputcom_1.isradi = 1;
/* -- " -- */
	inputcom_1.ispart = 1;
/* -- " -- */
	inputcom_1.iotail = 0;
/* cut tails */
    }

/*     check for scan function */

    if (i_indx(inputcom_1.scan, " ", (ftnlen)30, (ftnlen)1) != 1) {
	inputcom_1.iscan = chk_scan__(inputcom_1.scan, &inputcom_1.iscan, (
		ftnlen)30);
    }

    if (inputcom_1.iscan > 0) {
	if (inputcom_1.nscan <= 1) {
	    i__ = printerr_(&c_n9, "NSCAN too small for scan", (ftnlen)24);
	    inputcom_1.iscan = 0;
	} else {
	    if (inputcom_1.magin + inputcom_1.magout > 0 && (inputcom_1.iscan 
		    == 5 || inputcom_1.iscan == 6 || inputcom_1.iscan == 13 ||
		     inputcom_1.iscan == 14)) {
		if (inputcom_1.magin != 0) {
		    i__ = printerr_(&c_n14, inputcom_1.maginfile, (ftnlen)30);
		    inputcom_1.iscan = 0;
		} else {
		    i__ = printerr_(&c_n15, inputcom_1.magoutfile, (ftnlen)30)
			    ;
		    inputcom_1.iscan = 0;
		}
	    } else {
		inputcom_1.nslice = inputcom_1.nscan;
		inputcom_1.iall = 1;
		inputcom_1.zsep = inputcom_1.delz;
		inputcom_1.curlen = -1.;
		inputcom_1.ntail = 0;
		inputcom_1.iotail = 0;
		inputcom_1.shotnoise = (float)0.;
		if (inputcom_1.itdp != 0) {
		    i__ = printerr_(&c_n5, " ", (ftnlen)1);
		    inputcom_1.itdp = 0;
		}
		if (inputcom_1.iscan > 22 && i_indx(inputcom_1.beamfile, 
			" ", (ftnlen)30, (ftnlen)1) == 0) {
		    i__ = printerr_(&c_n8, "scan feature requires BEAMFILE", (
			    ftnlen)30);
		    last_();
		}
	    }
	}
    }

/*     check for magnet field input & output */

    if (i_indx(inputcom_1.maginfile, " ", (ftnlen)30, (ftnlen)1) != 1) {
	inputcom_1.magin = 1;
    } else {
	if (inputcom_1.magin != 0) {
	    if (inputcom_1.ilog != 0) {
		i__ = printerr_(&c_n21, "MAGINFILE", (ftnlen)9);
/* interaction required */
		last_();
	    }
L1:
	    s_wsfe(&io___3);
	    e_wsfe();
	    s_rsfe(&io___4);
	    do_fio(&c__1, inputcom_1.maginfile, (ftnlen)30);
	    e_rsfe();
/* get magnetic input file name */
	    if (i_indx(inputcom_1.maginfile, " ", (ftnlen)30, (ftnlen)1) == 1)
		     {
		i__ = printerr_(&c_n18, inputcom_1.maginfile, (ftnlen)30);
		goto L1;
	    }
	}
    }
    if (i_indx(inputcom_1.magoutfile, " ", (ftnlen)30, (ftnlen)1) != 1) {
	inputcom_1.magout = 1;
    } else {
	if (inputcom_1.magout != 0) {
	    if (inputcom_1.ilog != 0) {
		i__ = printerr_(&c_n21, "MAGOUTFILE", (ftnlen)10);
/* interaction required */
		last_();
	    }
L2:
	    s_wsfe(&io___5);
	    e_wsfe();
	    s_rsfe(&io___6);
	    do_fio(&c__1, inputcom_1.magoutfile, (ftnlen)30);
	    e_rsfe();
/* get magnetic input file name */
	    if (i_indx(inputcom_1.magoutfile, " ", (ftnlen)30, (ftnlen)1) == 
		    1) {
		i__ = printerr_(&c_n18, inputcom_1.magoutfile, (ftnlen)30);
		goto L2;
	    }
	}
    }

/*     check for output parameters */

    if (inputcom_1.ishsty < 1) {
	inputcom_1.ishsty = 1;
    }
    if (inputcom_1.ispart < 1) {
	inputcom_1.ispart = 1;
    }
    if (inputcom_1.isradi < 1) {
	inputcom_1.isradi = 1;
    }
    ix = 0;
    for (i__ = 1; i__ <= 39; ++i__) {
	if (inputcom_1.lout[i__ - 1] != 0) {
	    ++ix;
	    inputcom_1.lout[i__ - 1] = 1;
	}
    }
    if (ix == 0) {
	inputcom_1.iphsty = 0;
    }
    return 0;

/*     format statements */


} /* chk_input__ */

/* chk_input */
integer chk_bnd__()
{
    /* System generated locals */
    integer ret_val;
    doublereal d__1;

    /* Local variables */
    static integer ibas[7], itmp;
    extern /* Subroutine */ int last_();
    extern integer printerr_();
    static integer i__, i1, i2;

/*     ================================================================== */
/*     checks some boundaries of the input file. */
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
/*     time dependency of bunch + slippage field */





/*     ------------------------------------------------------------------ */
/*     input/output control */






/*     check for the case that convharm is set if no partfile is defined */

    if (iocom_1.npin <= 0) {
	inputcom_1.convharm = 1;
    }

    ret_val = 0;
    itmp = 0;

/*     case if nslice is smaller than 1  (auto-adjustment) */

    if (inputcom_1.nslice <= 0) {
	if (inputcom_1.curlen < 0.) {
/* step profile */
	    inputcom_1.ntail = 0;
	    inputcom_1.nslice = (integer) (abs(inputcom_1.curlen) / 
		    inputcom_1.xlamds / inputcom_1.zsep);
	} else {
/* gaussian */
	    inputcom_1.ntail = -((integer) (inputcom_1.curlen * 3. / 
		    inputcom_1.xlamds / inputcom_1.zsep));
	    inputcom_1.nslice = (integer) (inputcom_1.curlen * 6. / 
		    inputcom_1.xlamds / inputcom_1.zsep);
	}
    }

/*     adjustment for nslice if beamfile determines the scan */

    if (inputcom_1.iscan > 22) {
	if (tbunchcom_1.ndata <= 0) {
	    i__ = printerr_(&c_n8, "BEAMFILE for scan not defined", (ftnlen)
		    29);
	    last_();
	}
	inputcom_1.nslice = tbunchcom_1.ndata;
	inputcom_1.nscan = tbunchcom_1.ndata;
    }
    if (inputcom_1.aw0 <= 0. && inputcom_1.magin == 0) {
	itmp = printerr_(&c_n8, "No resonable wiggler field defined", (ftnlen)
		34);
/* a */
    }
    if (inputcom_1.nwig <= 0) {
	itmp = printerr_(&c_n8, "NWIG must be positive and non-zero", (ftnlen)
		34);
    }
    if (inputcom_1.delz <= 0.) {
	itmp = printerr_(&c_n8, "DELZ must be positive and non-zero", (ftnlen)
		34);
    }
    if (inputcom_1.zsep < 1. && inputcom_1.itdp == 1) {
	itmp = printerr_(&c_n8, "ZSEP must be al least 1", (ftnlen)23);
    }
    if (inputcom_1.xlamd <= 0.) {
	itmp = printerr_(&c_n8, "XLAMD must be positive and non-zero", (
		ftnlen)35);
    }
    if (inputcom_1.gamma0 - abs(inputcom_1.delgam) * 4 < 1.) {
	itmp = printerr_(&c_n8, "energy GAMMA0 too small", (ftnlen)23);
/* abort */
    }
    if (inputcom_1.npart > 1000001) {
	i__ = printerr_(&c_n9, "NPART > NPMAX - setting NPART=NPMAX", (ftnlen)
		35);
    }
    if (inputcom_1.nbins < 4) {
	i__ = printerr_(&c_n9, "NBINS too small - setting NBINS=4", (ftnlen)
		33);
	inputcom_1.nbins = 4;
    }
    if (inputcom_1.npart % (inputcom_1.nbins << 2) != 0) {
	itmp = printerr_(&c_n8, "NPART not a multiple of 4*NBINS", (ftnlen)31)
		;
/* abor */
    }
    for (i1 = inputcom_1.nharm + 1; i1 <= 7; ++i1) {
	if (inputcom_1.lout[i1 + 13] != 0) {
	    i__ = printerr_(&c_n9, "no harmonic output above NHARM", (ftnlen)
		    30);
	}
    }
    if (inputcom_1.idmppar > inputcom_1.nharm) {
	i__ = printerr_(&c_n9, "No dump possible (IDMPPAR > IHARM)", (ftnlen)
		34);
	inputcom_1.idmppar = 0;
    }
    if (inputcom_1.xlamds <= 0.) {
	itmp = printerr_(&c_n8, "XLAMDS must be positive", (ftnlen)23);
/* abort */
    }
    if (inputcom_1.prad0 < 0.) {
	itmp = printerr_(&c_n8, "PRAD0 must not be negative", (ftnlen)26);
/* abort */
    }
    if (inputcom_1.pradh0 < 0.) {
	itmp = printerr_(&c_n8, "PRADH0 must not be negative", (ftnlen)27);
/* abort */
    }
    ibas[0] = inputcom_1.ildpsi;
    ibas[1] = inputcom_1.ildx;
    ibas[2] = inputcom_1.ildy;
    ibas[3] = inputcom_1.ildpx;
    ibas[4] = inputcom_1.ildpy;
    ibas[5] = inputcom_1.ildgam;
    ibas[6] = inputcom_1.ildgam + 1;
    for (i1 = 1; i1 <= 7; ++i1) {
	for (i2 = i1 + 1; i2 <= 7; ++i2) {
	    if (ibas[i1 - 1] == ibas[i2 - 1]) {
/* no abort */
		i__ = printerr_(&c_n9, "Identical bases in Hammersley sequen\
ces", (ftnlen)39);
	    }
	}
    }
    if (abs(inputcom_1.iertyp) != 0 && (d__1 = inputcom_1.delz - .5, abs(d__1)
	    ) > 1e-7) {
	itmp = printerr_(&c_n8, "DELZ must be 0.5 for field errors", (ftnlen)
		33);
    }
    if (inputcom_1.iscan > 25) {
	i__ = printerr_(&c_n9, "Invalid scan parameter - setting ISCAN=0", (
		ftnlen)40);
	inputcom_1.iscan = 0;
    }
    if (inputcom_1.iscan > 0 && (inputcom_1.nscan <= 1 || abs(inputcom_1.svar)
	     < 1e-7)) {
	itmp = printerr_(&c_n8, "Invalid scan range (NSCAN,SVAL)", (ftnlen)31)
		;
    }
    if (inputcom_1.zsep / inputcom_1.delz - (integer) (inputcom_1.zsep / 
	    inputcom_1.delz) * (float)1. > 1e-7) {
	itmp = printerr_(&c_n8, "ZSEP not a multiple of DELZ", (ftnlen)27);
    }
    if (inputcom_1.nslice > 50000) {
	itmp = printerr_(&c_n8, "Too many slices (NSLICE>NSMAX)", (ftnlen)30);
    }
    if (inputcom_1.nslice <= 0) {
	itmp = printerr_(&c_n8, "NSLICE < 1", (ftnlen)10);
    }
    if (inputcom_1.zrayl <= (float)0.) {
	itmp = printerr_(&c_n8, "ZRAYL must be larger than 0", (ftnlen)27);
    }
    if (inputcom_1.ncar % 2 == 0) {
	itmp = printerr_(&c_n8, "NCAR not an odd integer", (ftnlen)23);
    }
    if (inputcom_1.ncar > 261) {
	i__ = printerr_(&c_n9, "NCAR too large - setting NCAR=NCMAX", (ftnlen)
		35);
    }
    if (inputcom_1.nptr > 999) {
	i__ = printerr_(&c_n9, "NPTR too large - setting NPTR=NRGRID", (
		ftnlen)36);
    }
    if (inputcom_1.nptr < 2) {
	i__ = printerr_(&c_n9, "NPTR too small - disabling space charge", (
		ftnlen)39);
	inputcom_1.nscz = 0;
	inputcom_1.nscr = 0;
    }
    if (inputcom_1.nscz >= inputcom_1.nbins / 2 + 1) {
/* somehow empirical boundary */
	i__ = printerr_(&c_n9, "NSCZ too large - setting NSCZ=2", (ftnlen)31);
	inputcom_1.nscz = 2;
    }
    if (inputcom_1.nharm > inputcom_1.nbins / 2 + 1) {
	i__ = printerr_(&c_n9, "Higher harmonics are inaccurate (NHARM)", (
		ftnlen)39);
    }

    ret_val = itmp;
    return ret_val;
} /* chk_bnd__ */


/* chk_bnd */
integer chk_scan__(c0, iscn, c0_len)
char *c0;
integer *iscn;
ftnlen c0_len;
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_cmp();

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int touppercase_();

/*     ============================================================ */
/*     check for string in input scan */
/*     ------------------------------------------------------------ */


    ret_val = *iscn;
/* if not found use value if iscan */
    touppercase_(c0, (ftnlen)30);
    j = 1;
    for (i__ = 1; i__ <= 30; ++i__) {
/* conversion to uppercase */
	if (*(unsigned char *)&c0[i__ - 1] > ' ') {
	    j = i__;
	}
/* end of string? */
    }
    if (s_cmp(c0, "GAMMA0", j, (ftnlen)6) == 0) {
	ret_val = 1;
    }
    if (s_cmp(c0, "DELGAM", j, (ftnlen)6) == 0) {
	ret_val = 2;
    }
    if (s_cmp(c0, "CURPEAK", j, (ftnlen)7) == 0) {
	ret_val = 3;
    }
    if (s_cmp(c0, "XLAMDS", j, (ftnlen)6) == 0) {
	ret_val = 4;
    }
    if (s_cmp(c0, "AW0", j, (ftnlen)3) == 0) {
	ret_val = 5;
    }
    if (s_cmp(c0, "ISEED", j, (ftnlen)5) == 0) {
	ret_val = 6;
    }
    if (s_cmp(c0, "PXBEAM", j, (ftnlen)6) == 0) {
	ret_val = 7;
    }
    if (s_cmp(c0, "PYBEAM", j, (ftnlen)6) == 0) {
	ret_val = 8;
    }
    if (s_cmp(c0, "RXBEAM", j, (ftnlen)6) == 0) {
	ret_val = 11;
    }
    if (s_cmp(c0, "RYBEAM", j, (ftnlen)6) == 0) {
	ret_val = 12;
    }
    if (s_cmp(c0, "XBEAM", j, (ftnlen)5) == 0) {
	ret_val = 9;
    }
    if (s_cmp(c0, "YBEAM", j, (ftnlen)5) == 0) {
	ret_val = 10;
    }
    if (s_cmp(c0, "XLAMD", j, (ftnlen)5) == 0) {
	ret_val = 13;
    }
    if (s_cmp(c0, "DELAW", j, (ftnlen)5) == 0) {
	ret_val = 14;
    }
    if (s_cmp(c0, "ALPHAX", j, (ftnlen)6) == 0) {
	ret_val = 15;
    }
    if (s_cmp(c0, "ALPHAY", j, (ftnlen)6) == 0) {
	ret_val = 16;
    }
    if (s_cmp(c0, "EMITX", j, (ftnlen)5) == 0) {
	ret_val = 17;
    }
    if (s_cmp(c0, "EMITY", j, (ftnlen)5) == 0) {
	ret_val = 18;
    }
    if (s_cmp(c0, "PRAD0", j, (ftnlen)5) == 0) {
	ret_val = 19;
    }
    if (s_cmp(c0, "ZRAYL", j, (ftnlen)5) == 0) {
	ret_val = 20;
    }
    if (s_cmp(c0, "ZWAIST", j, (ftnlen)6) == 0) {
	ret_val = 21;
    }
    if (s_cmp(c0, "AWD", j, (ftnlen)3) == 0) {
	ret_val = 22;
    }
    if (s_cmp(c0, "BEAMFILE", j, (ftnlen)8) == 0) {
	ret_val = 23;
    }
    if (s_cmp(c0, "BEAMOPT", j, (ftnlen)7) == 0) {
	ret_val = 24;
    }
    if (s_cmp(c0, "BEAMGAM", j, (ftnlen)7) == 0) {
	ret_val = 25;
    }
    return ret_val;
} /* chk_scan__ */

