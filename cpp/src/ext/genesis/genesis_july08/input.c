/* input.f -- translated by f2c (version 20000118).
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
    doublereal distversion, distrev;
    integer iout[39], nout, nfld, npar, ndump, firstout, nfin, irecpar, 
	    irecfld, kout, ndmp2, npin, nlog, ndis, ncoldis, iconv2t, iconv2g,
	     icolpar[10], ndistsize, iconv2px, iconv2py, nprobe, ftype, 
	    ftdist, ftpart, ftfield, ndumph[6], nfldh[6];
} iocom_;

#define iocom_1 iocom_

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
    doublereal xkw0, xkper0, sval, svalout, gamma0_in__;
    integer npart0, inorun;
} simcom_;

#define simcom_1 simcom_

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

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__9 = 9;
static integer c_n1 = -1;
static integer c__10 = 10;
static integer c__11 = 11;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__0 = 0;
static integer c__17 = 17;
static integer c_n2 = -2;
static integer c__30 = 30;
//static integer c_n3 = -3; //OC030110
//static integer c_n4 = -4;
static integer c_n6 = -6;
static integer c__8 = 8;
static integer c_n8 = -8;
static integer c_n9 = -9;
static integer c__2 = 2;
static integer c_n19 = -19;

/* Subroutine */ int initio_()
{
    /* Format strings */
    static char fmt_2[] = "(\002.node\002,i6.6)";
    static char fmt_100[] = "(a3,i1.1)";

    /* System generated locals */
    address a__1[3];
    integer i__1[3], i__2, i__3, i__4;
    olist o__1;

    /* Builtin functions */
    integer i_indx();
    /* Subroutine */ int s_copy();
    integer s_wsfi(), do_fio(), e_wsfi();
    /* Subroutine */ int s_cat();
    integer s_wsle(), do_lio(), e_wsle(), f_open();

    /* Local variables */
    static char file[34];
    extern integer readdistfile_();
    static char file_ext__[4];
    extern /* Subroutine */ int last_(), readbeamfile_();
    static integer ierr1, ierr2;
    extern integer printerr_(), openbininput_();
    static integer i__;
    extern /* Subroutine */ int first_(), chk_input__();
    static integer ih;
    extern integer readin_(), strlen_(), chk_bnd__(), openoutputfile_();
    static char file_id__[11];
    extern /* Subroutine */ int readradfile_();
    extern integer openbinfile_();
    extern /* Subroutine */ int opentimerec_(), touppercase_();

    /* Fortran I/O blocks */
    static icilist io___5 = { 0, file_id__, 0, fmt_2, 11, 1 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };
    static icilist io___10 = { 0, file_ext__, 0, fmt_100, 4, 1 };
    static icilist io___11 = { 0, file_ext__, 0, fmt_100, 4, 1 };


/*     ================================================================== */
/*     manages all initial input for genesis */
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
/*     input/output control */




/*     ------------------------------------------------------------------ */
/*     common block for mpi implementation for genesis */




/*     -------------------------------------------------------------------- */
/*     cartesian mesh */




    iocom_1.nprobe = 30;
/* filenumber for filetype probin */
    iocom_1.nlog = 6;
    first_();
    ierr1 = readin_();

/* read namelist */
    touppercase_(inputcom_1.filetype, (ftnlen)30);
/* convert to upper case let */
    iocom_1.ftype = 1;
    i__ = i_indx(inputcom_1.filetype, "SDDS", (ftnlen)30, (ftnlen)4);
    if (i__ > 0) {
	iocom_1.ftype = 2;
    }

    if (inputcom_1.ilog != 0) {
	if (i_indx(inputcom_1.outputfile, " ", (ftnlen)30, (ftnlen)1) == 1) {
	    s_copy(file, "log-file", (ftnlen)34, (ftnlen)8);
	    goto L5;
	} else {
	    s_copy(file_id__, "", (ftnlen)11, (ftnlen)0);
	    if (mpicom_1.mpi_size__ > 1) {
		s_wsfi(&io___5);
		do_fio(&c__1, (char *)&mpicom_1.mpi_id__, (ftnlen)sizeof(
			integer));
		e_wsfi();
/* apppend process ID if run i */
	    }
/* Writing concatenation */
	    i__1[0] = strlen_(inputcom_1.outputfile, (ftnlen)30), a__1[0] = 
		    inputcom_1.outputfile;
	    i__1[1] = 4, a__1[1] = ".log";
	    i__1[2] = 11, a__1[2] = file_id__;
	    s_cat(file, a__1, i__1, &c__3, (ftnlen)34);
	    s_wsle(&io___6);
	    do_lio(&c__9, &c__1, "logfile: ", (ftnlen)9);
	    do_lio(&c__9, &c__1, file, (ftnlen)34);
	    e_wsle();
	}
	o__1.oerr = 1;
	o__1.ounit = 16;
	o__1.ofnmlen = 34;
	o__1.ofnm = file;
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__2 = f_open(&o__1);
	if (i__2 != 0) {
	    goto L5;
	}
	iocom_1.nlog = 16;
L5:
	if (iocom_1.nlog == 6) {
	    i__ = printerr_(&c_n1, file, (ftnlen)34);
	}
    }

    ierr2 = openoutputfile_(&ierr1, "template.in", (ftnlen)11);
/* open outputfile+ */
    if (ierr2 + ierr1 < 0) {
	last_();
    }

/* error occured? */
    chk_input__();

/*     open addition output files */

/* make input consistent */
    if (inputcom_1.ipradi > 0) {
/* is output of radiation field selected? */
	iocom_1.irecfld = 1;
	i__2 = (inputcom_1.ncar << 3) * inputcom_1.ncar;
	iocom_1.nfld = openbinfile_(inputcom_1.outputfile, "fld ", &c__10, &
		i__2, (ftnlen)30, (ftnlen)4);
/* real and imag seprated */
	i__2 = cartcom_1.nhloop;
	for (ih = 2; ih <= i__2; ++ih) {
	    s_wsfi(&io___10);
	    do_fio(&c__1, "fld", (ftnlen)3);
	    do_fio(&c__1, (char *)&cartcom_1.hloop[ih - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfi();
	    i__3 = ih + 20;
	    i__4 = (inputcom_1.ncar << 3) * inputcom_1.ncar;
	    iocom_1.nfldh[ih - 2] = openbinfile_(inputcom_1.outputfile, 
		    file_ext__, &i__3, &i__4, (ftnlen)30, (ftnlen)4);
	}
    } else {
	iocom_1.nfld = -1;
    }

    if (inputcom_1.ippart > 0) {
/* is output of particle distribution desired */
	iocom_1.irecpar = 1;
	i__2 = inputcom_1.npart << 3;
	iocom_1.npar = openbinfile_(inputcom_1.outputfile, "par ", &c__11, &
		i__2, (ftnlen)30, (ftnlen)4);
    } else {
	iocom_1.npar = -1;
    }
    if (inputcom_1.idmpfld != 0) {
/* dumped radiation field? */
	i__2 = (inputcom_1.ncar << 4) * inputcom_1.ncar;
	iocom_1.ndump = openbinfile_(inputcom_1.outputfile, "dfl ", &c__12, &
		i__2, (ftnlen)30, (ftnlen)4);
	i__2 = cartcom_1.nhloop;
	for (ih = 2; ih <= i__2; ++ih) {
	    s_wsfi(&io___11);
	    do_fio(&c__1, "dfl", (ftnlen)3);
	    do_fio(&c__1, (char *)&cartcom_1.hloop[ih - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfi();
	    i__3 = ih + 30;
	    i__4 = (inputcom_1.ncar << 4) * inputcom_1.ncar;
	    iocom_1.ndumph[ih - 2] = openbinfile_(inputcom_1.outputfile, 
		    file_ext__, &i__3, &i__4, (ftnlen)30, (ftnlen)4);
	}
    } else {
	iocom_1.ndump = -1;
    }

    if (inputcom_1.idmppar != 0) {
/* should the radiation field be dumped at */
	i__2 = inputcom_1.npart << 3;
	iocom_1.ndmp2 = openbinfile_(inputcom_1.outputfile, "dpa ", &c__13, &
		i__2, (ftnlen)30, (ftnlen)4);
    } else {
	iocom_1.ndmp2 = -1;
    }

    if (inputcom_1.itdp != 0) {
	opentimerec_(&inputcom_1.ncar);
    }

/*     open additional input files (maginfile is opened in magfield.f) */

/* prepare the crtime */
    readbeamfile_(inputcom_1.beamfile, (ftnlen)30);
/* read external description file for */
    readradfile_(inputcom_1.radfile, (ftnlen)30);

/* read external description file for rad */
    i__2 = inputcom_1.ncar * inputcom_1.ncar << 4;
    iocom_1.nfin = openbininput_(inputcom_1.fieldfile, &c__14, &i__2, &c__1, (
	    ftnlen)30);
/* returns -1 if no n */
    i__2 = inputcom_1.npart << 3;
    iocom_1.npin = openbininput_(inputcom_1.partfile, &c__15, &i__2, &c__0, (
	    ftnlen)30);
/* ---------- " ----- */
    iocom_1.ndis = readdistfile_(inputcom_1.distfile, &c__17, (ftnlen)30);

/*     check for boundary violation or unphysical input parameters */

/* ---------- " ----- */
    ierr1 = chk_bnd__();
/* check for boundary violation */
    if (ierr1 < 0) {
	last_();
    }

    return 0;

} /* initio_ */


/* of ioinit */
integer openbininput_(file, nio, size, isfield, file_len)
char *file;
integer *nio, *size, *isfield;
ftnlen file_len;
{
    /* System generated locals */
    integer ret_val, i__1;
    olist o__1;

    /* Builtin functions */
    integer i_indx(), f_open();

    /* Local variables */
    extern /* Subroutine */ int last_();
    extern integer printerr_();
    static integer ft;
    extern integer detectfiletype_();

/*     ================================================================= */
/*     opens binary input files (filed and part files) and checks for */
/*     the filetype */
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
/*     input/output control */




    ret_val = -1;
    if (i_indx(file, " ", file_len, (ftnlen)1) == 1) {
	return ret_val;
    }
/* no file selected */
    ret_val = *nio;
    ft = detectfiletype_(file, file_len);
    if (*isfield == 1) {
	iocom_1.ftfield = ft;
/* is field */
    } else {
	iocom_1.ftpart = ft;
/* is particle */
    }

    o__1.oerr = 1;
    o__1.ounit = *nio;
    o__1.ofnmlen = file_len;
    o__1.ofnm = file;
    o__1.orl = *size;
    o__1.osta = "old";
    o__1.oacc = "direct";
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L100;
    }
    return ret_val;

L100:
    ret_val = printerr_(&c_n1, file, file_len);
    last_();
    return ret_val;
} /* openbininput_ */



integer detectfiletype_(file, file_len)
char *file;
ftnlen file_len;
{
    /* Format strings */
    static char fmt_200[] = "(a)";

    /* System generated locals */
    integer ret_val, i__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), f_clos(), i_indx();

    /* Local variables */
    static char line[80];
    extern /* Subroutine */ int last_();
    extern integer printerr_();
    static integer i__;
    extern /* Subroutine */ int touppercase_();

    /* Fortran I/O blocks */
    static cilist io___13 = { 1, 0, 0, fmt_200, 0 };


/*     ================================================================= */
/*     the routine tries to read the beginning of the file */
/*     if the first line contains sdds then it returns the constant */
/*     sdds, otherwise it returns the constant original */
/*     ---------------------------------------------------------------- */





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
/*     input/output control */





    ret_val = 1;
    o__1.oerr = 1;
    o__1.ounit = iocom_1.nprobe;
    o__1.ofnmlen = file_len;
    o__1.ofnm = file;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L100;
    }
    io___13.ciunit = iocom_1.nprobe;
    i__1 = s_rsfe(&io___13);
    if (i__1 != 0) {
	goto L110;
    }
    i__1 = do_fio(&c__1, line, (ftnlen)80);
    if (i__1 != 0) {
	goto L110;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L110;
    }
    cl__1.cerr = 0;
    cl__1.cunit = iocom_1.nprobe;
    cl__1.csta = 0;
    f_clos(&cl__1);
    touppercase_(line, (ftnlen)80);
    i__ = i_indx(line, "SDDS", (ftnlen)80, (ftnlen)4);
    if (i__ > 0) {
	ret_val = 2;
    }
    return ret_val;

L100:
    ret_val = printerr_(&c_n1, file, file_len);
    last_();
    return ret_val;
L110:
    ret_val = printerr_(&c_n2, file, file_len);
    last_();
    return ret_val;

} /* detectfiletype_ */



/* of detectfiletype */
integer readin_()
{
    /* Format strings */
	//OC030110
    //static char fmt_100[] = "(\002Please enter input file name \002)";
    //static char fmt_110[] = "(a30)";

    /* System generated locals */
    integer ret_val; //, i__1; //OC port
    olist o__1;
    //cllist cl__1; //OC port

    /* Builtin functions */
    integer s_wsfe(), e_wsfe(), s_rsfe(), do_fio(), e_rsfe(), f_open(), 
	    s_rsne(), f_clos();

    /* Local variables */
    extern /* Subroutine */ int last_();
    extern integer printerr_();
    extern /* Subroutine */ int mpi_bcast__(), preset_();
    static integer nin;

    /* Namelist stuff */
/**  //OC030110
    static Vardesc iall_dv = { "IALL", (char *)&inputcom_1.iall, (ftnlen *)0, 
	    3 };
    static Vardesc ncar_dv = { "NCAR", (char *)&inputcom_1.ncar, (ftnlen *)0, 
	    3 };
    static Vardesc emod_dv = { "EMOD", (char *)&inputcom_1.emod, (ftnlen *)0, 
	    5 };
    static Vardesc beamfile_dv = { "BEAMFILE", inputcom_1.beamfile, (ftnlen *)
	    0, -30 };
    static Vardesc scan_dv = { "SCAN", inputcom_1.scan, (ftnlen *)0, -30 };
    static Vardesc nsec_dv = { "NSEC", (char *)&inputcom_1.nsec, (ftnlen *)0, 
	    3 };
    static Vardesc ilog_dv = { "ILOG", (char *)&inputcom_1.ilog, (ftnlen *)0, 
	    3 };
    static Vardesc iorb_dv = { "IORB", (char *)&inputcom_1.iorb, (ftnlen *)0, 
	    3 };
    static Vardesc delz_dv = { "DELZ", (char *)&inputcom_1.delz, (ftnlen *)0, 
	    5 };
    static Vardesc ildx_dv = { "ILDX", (char *)&inputcom_1.ildx, (ftnlen *)0, 
	    3 };
    static Vardesc ildy_dv = { "ILDY", (char *)&inputcom_1.ildy, (ftnlen *)0, 
	    3 };
    static Vardesc qfdx_dv = { "QFDX", (char *)&inputcom_1.qfdx, (ftnlen *)0, 
	    5 };
    static Vardesc qfdy_dv = { "QFDY", (char *)&inputcom_1.qfdy, (ftnlen *)0, 
	    5 };
    static Vardesc nwig_dv = { "NWIG", (char *)&inputcom_1.nwig, (ftnlen *)0, 
	    3 };
    static Vardesc nscr_dv = { "NSCR", (char *)&inputcom_1.nscr, (ftnlen *)0, 
	    3 };
    static Vardesc itdp_dv = { "ITDP", (char *)&inputcom_1.itdp, (ftnlen *)0, 
	    3 };
    static Vardesc iallharm_dv = { "IALLHARM", (char *)&inputcom_1.iallharm, (
	    ftnlen *)0, 3 };
    static Vardesc svar_dv = { "SVAR", (char *)&inputcom_1.svar, (ftnlen *)0, 
	    5 };
    static Vardesc igamgaus_dv = { "IGAMGAUS", (char *)&inputcom_1.igamgaus, (
	    ftnlen *)0, 3 };
    static Vardesc nscz_dv = { "NSCZ", (char *)&inputcom_1.nscz, (ftnlen *)0, 
	    3 };
    static Vardesc zsep_dv = { "ZSEP", (char *)&inputcom_1.zsep, (ftnlen *)0, 
	    5 };
    static Vardesc nptr_dv = { "NPTR", (char *)&inputcom_1.nptr, (ftnlen *)0, 
	    3 };
    static ftnlen lout_dims[] = { 1, 40, 1 };
    static Vardesc lout_dv = { "LOUT", (char *)inputcom_1.lout, lout_dims, 3 }
	    ;
    static Vardesc partfile_dv = { "PARTFILE", inputcom_1.partfile, (ftnlen *)
	    0, -30 };
    static Vardesc distfile_dv = { "DISTFILE", inputcom_1.distfile, (ftnlen *)
	    0, -30 };
    static Vardesc convharm_dv = { "CONVHARM", (char *)&inputcom_1.convharm, (
	    ftnlen *)0, 3 };
    static Vardesc filetype_dv = { "FILETYPE", inputcom_1.filetype, (ftnlen *)
	    0, -30 };
    static Vardesc prad0_dv = { "PRAD0", (char *)&inputcom_1.prad0, (ftnlen *)
	    0, 5 };
    static Vardesc multconv_dv = { "MULTCONV", (char *)&inputcom_1.multconv, (
	    ftnlen *)0, 3 };
    static Vardesc rmax0_dv = { "RMAX0", (char *)&inputcom_1.rmax0, (ftnlen *)
	    0, 5 };
    static Vardesc iseed_dv = { "ISEED", (char *)&inputcom_1.iseed, (ftnlen *)
	    0, 3 };
    static Vardesc dgrid_dv = { "DGRID", (char *)&inputcom_1.dgrid, (ftnlen *)
	    0, 5 };
    static Vardesc magin_dv = { "MAGIN", (char *)&inputcom_1.magin, (ftnlen *)
	    0, 3 };
    static Vardesc delaw_dv = { "DELAW", (char *)&inputcom_1.delaw, (ftnlen *)
	    0, 5 };
    static Vardesc xbeam_dv = { "XBEAM", (char *)&inputcom_1.xbeam, (ftnlen *)
	    0, 5 };
    static Vardesc ybeam_dv = { "YBEAM", (char *)&inputcom_1.ybeam, (ftnlen *)
	    0, 5 };
    static Vardesc bunch_dv = { "BUNCH", (char *)&inputcom_1.bunch, (ftnlen *)
	    0, 5 };
    static Vardesc quadf_dv = { "QUADF", (char *)&inputcom_1.quadf, (ftnlen *)
	    0, 5 };
    static Vardesc quadd_dv = { "QUADD", (char *)&inputcom_1.quadd, (ftnlen *)
	    0, 5 };
    static Vardesc fieldfile_dv = { "FIELDFILE", inputcom_1.fieldfile, (
	    ftnlen *)0, -30 };
    static Vardesc iscan_dv = { "ISCAN", (char *)&inputcom_1.iscan, (ftnlen *)
	    0, 3 };
    static Vardesc xlamd_dv = { "XLAMD", (char *)&inputcom_1.xlamd, (ftnlen *)
	    0, 5 };
    static Vardesc nbins_dv = { "NBINS", (char *)&inputcom_1.nbins, (ftnlen *)
	    0, 3 };
    static Vardesc nharm_dv = { "NHARM", (char *)&inputcom_1.nharm, (ftnlen *)
	    0, 3 };
    static Vardesc idump_dv = { "IDUMP", (char *)&inputcom_1.idump, (ftnlen *)
	    0, 3 };
    static Vardesc solen_dv = { "SOLEN", (char *)&inputcom_1.solen, (ftnlen *)
	    0, 5 };
    static Vardesc ildpx_dv = { "ILDPX", (char *)&inputcom_1.ildpx, (ftnlen *)
	    0, 3 };
    static Vardesc ildpy_dv = { "ILDPY", (char *)&inputcom_1.ildpy, (ftnlen *)
	    0, 3 };
    static Vardesc emodphase_dv = { "EMODPHASE", (char *)&
	    inputcom_1.emodphase, (ftnlen *)0, 5 };
    static Vardesc emitx_dv = { "EMITX", (char *)&inputcom_1.emitx, (ftnlen *)
	    0, 5 };
    static Vardesc emity_dv = { "EMITY", (char *)&inputcom_1.emity, (ftnlen *)
	    0, 5 };
    static Vardesc npart_dv = { "NPART", (char *)&inputcom_1.npart, (ftnlen *)
	    0, 3 };
    static Vardesc maginfile_dv = { "MAGINFILE", inputcom_1.maginfile, (
	    ftnlen *)0, -30 };
    static Vardesc ntail_dv = { "NTAIL", (char *)&inputcom_1.ntail, (ftnlen *)
	    0, 3 };
    static Vardesc ndcut_dv = { "NDCUT", (char *)&inputcom_1.ndcut, (ftnlen *)
	    0, 3 };
    static Vardesc gamma0_dv = { "GAMMA0", (char *)&inputcom_1.gamma0, (
	    ftnlen *)0, 5 };
    static Vardesc zrayl_dv = { "ZRAYL", (char *)&inputcom_1.zrayl, (ftnlen *)
	    0, 5 };
    static Vardesc nscan_dv = { "NSCAN", (char *)&inputcom_1.nscan, (ftnlen *)
	    0, 3 };
    static Vardesc alignradf_dv = { "ALIGNRADF", (char *)&
	    inputcom_1.alignradf, (ftnlen *)0, 3 };
    static Vardesc eloss_dv = { "ELOSS", (char *)&inputcom_1.eloss, (ftnlen *)
	    0, 5 };
    static Vardesc imagl_dv = { "IMAGL", (char *)&inputcom_1.imagl, (ftnlen *)
	    0, 5 };
    static Vardesc idril_dv = { "IDRIL", (char *)&inputcom_1.idril, (ftnlen *)
	    0, 5 };
    static Vardesc trama_dv = { "TRAMA", (char *)&inputcom_1.trama, (ftnlen *)
	    0, 3 };
    static Vardesc pradh0_dv = { "PRADH0", (char *)&inputcom_1.pradh0, (
	    ftnlen *)0, 5 };
    static Vardesc zstop_dv = { "ZSTOP", (char *)&inputcom_1.zstop, (ftnlen *)
	    0, 5 };
    static Vardesc fbess0_dv = { "FBESS0", (char *)&inputcom_1.fbess0, (
	    ftnlen *)0, 5 };
    static Vardesc shotnoise_dv = { "SHOTNOISE", (char *)&
	    inputcom_1.shotnoise, (ftnlen *)0, 5 };
    static Vardesc dl_dv = { "DL", (char *)&inputcom_1.dl, (ftnlen *)0, 5 };
    static Vardesc fl_dv = { "FL", (char *)&inputcom_1.fl, (ftnlen *)0, 5 };
    static Vardesc delgam_dv = { "DELGAM", (char *)&inputcom_1.delgam, (
	    ftnlen *)0, 5 };
    static Vardesc ildgam_dv = { "ILDGAM", (char *)&inputcom_1.ildgam, (
	    ftnlen *)0, 3 };
    static Vardesc sl_dv = { "SL", (char *)&inputcom_1.sl, (ftnlen *)0, 5 };
    static Vardesc ffspec_dv = { "FFSPEC", (char *)&inputcom_1.ffspec, (
	    ftnlen *)0, 3 };
    static Vardesc ipradi_dv = { "IPRADI", (char *)&inputcom_1.ipradi, (
	    ftnlen *)0, 3 };
    static Vardesc ipseed_dv = { "IPSEED", (char *)&inputcom_1.ipseed, (
	    ftnlen *)0, 3 };
    static Vardesc pxbeam_dv = { "PXBEAM", (char *)&inputcom_1.pxbeam, (
	    ftnlen *)0, 5 };
    static Vardesc alphax_dv = { "ALPHAX", (char *)&inputcom_1.alphax, (
	    ftnlen *)0, 5 };
    static Vardesc rxbeam_dv = { "RXBEAM", (char *)&inputcom_1.rxbeam, (
	    ftnlen *)0, 5 };
    static Vardesc rybeam_dv = { "RYBEAM", (char *)&inputcom_1.rybeam, (
	    ftnlen *)0, 5 };
    static Vardesc alphay_dv = { "ALPHAY", (char *)&inputcom_1.alphay, (
	    ftnlen *)0, 5 };
    static Vardesc pybeam_dv = { "PYBEAM", (char *)&inputcom_1.pybeam, (
	    ftnlen *)0, 5 };
    static Vardesc isradi_dv = { "ISRADI", (char *)&inputcom_1.isradi, (
	    ftnlen *)0, 3 };
    static Vardesc ildpsi_dv = { "ILDPSI", (char *)&inputcom_1.ildpsi, (
	    ftnlen *)0, 3 };
    static Vardesc iotail_dv = { "IOTAIL", (char *)&inputcom_1.iotail, (
	    ftnlen *)0, 3 };
    static Vardesc xlamds_dv = { "XLAMDS", (char *)&inputcom_1.xlamds, (
	    ftnlen *)0, 5 };
    static Vardesc curlen_dv = { "CURLEN", (char *)&inputcom_1.curlen, (
	    ftnlen *)0, 5 };
    static Vardesc nslice_dv = { "NSLICE", (char *)&inputcom_1.nslice, (
	    ftnlen *)0, 3 };
    static Vardesc itgaus_dv = { "ITGAUS", (char *)&inputcom_1.itgaus, (
	    ftnlen *)0, 3 };
    static ftnlen wcoefz_dims[] = { 1, 3, 1 };
    static Vardesc wcoefz_dv = { "WCOEFZ", (char *)inputcom_1.wcoefz, 
	    wcoefz_dims, 5 };
    static Vardesc magout_dv = { "MAGOUT", (char *)&inputcom_1.magout, (
	    ftnlen *)0, 3 };
    static Vardesc bunchphase_dv = { "BUNCHPHASE", (char *)&
	    inputcom_1.bunchphase, (ftnlen *)0, 5 };
    static Vardesc ippart_dv = { "IPPART", (char *)&inputcom_1.ippart, (
	    ftnlen *)0, 3 };
    static Vardesc isravg_dv = { "ISRAVG", (char *)&inputcom_1.isravg, (
	    ftnlen *)0, 3 };
    static Vardesc ispart_dv = { "ISPART", (char *)&inputcom_1.ispart, (
	    ftnlen *)0, 3 };
    static Vardesc isrsig_dv = { "ISRSIG", (char *)&inputcom_1.isrsig, (
	    ftnlen *)0, 3 };
    static Vardesc offsetradf_dv = { "OFFSETRADF", (char *)&
	    inputcom_1.offsetradf, (ftnlen *)0, 3 };
    static Vardesc itram11_dv = { "ITRAM11", (char *)&inputcom_1.itram11, (
	    ftnlen *)0, 5 };
    static Vardesc itram12_dv = { "ITRAM12", (char *)&inputcom_1.itram12, (
	    ftnlen *)0, 5 };
    static Vardesc aw0_dv = { "AW0", (char *)&inputcom_1.aw0, (ftnlen *)0, 5 }
	    ;
    static Vardesc itram13_dv = { "ITRAM13", (char *)&inputcom_1.itram13, (
	    ftnlen *)0, 5 };
    static Vardesc itram14_dv = { "ITRAM14", (char *)&inputcom_1.itram14, (
	    ftnlen *)0, 5 };
    static Vardesc magoutfile_dv = { "MAGOUTFILE", inputcom_1.magoutfile, (
	    ftnlen *)0, -30 };
    static Vardesc iertyp_dv = { "IERTYP", (char *)&inputcom_1.iertyp, (
	    ftnlen *)0, 3 };
    static Vardesc itram15_dv = { "ITRAM15", (char *)&inputcom_1.itram15, (
	    ftnlen *)0, 5 };
    static Vardesc itram16_dv = { "ITRAM16", (char *)&inputcom_1.itram16, (
	    ftnlen *)0, 5 };
    static Vardesc itram21_dv = { "ITRAM21", (char *)&inputcom_1.itram21, (
	    ftnlen *)0, 5 };
    static Vardesc iphsty_dv = { "IPHSTY", (char *)&inputcom_1.iphsty, (
	    ftnlen *)0, 3 };
    static Vardesc zwaist_dv = { "ZWAIST", (char *)&inputcom_1.zwaist, (
	    ftnlen *)0, 5 };
    static Vardesc itram22_dv = { "ITRAM22", (char *)&inputcom_1.itram22, (
	    ftnlen *)0, 5 };
    static Vardesc ishsty_dv = { "ISHSTY", (char *)&inputcom_1.ishsty, (
	    ftnlen *)0, 3 };
    static Vardesc itram23_dv = { "ITRAM23", (char *)&inputcom_1.itram23, (
	    ftnlen *)0, 5 };
    static Vardesc iwityp_dv = { "IWITYP", (char *)&inputcom_1.iwityp, (
	    ftnlen *)0, 3 };
    static Vardesc isntyp_dv = { "ISNTYP", (char *)&inputcom_1.isntyp, (
	    ftnlen *)0, 3 };
    static Vardesc itram24_dv = { "ITRAM24", (char *)&inputcom_1.itram24, (
	    ftnlen *)0, 5 };
    static Vardesc itram25_dv = { "ITRAM25", (char *)&inputcom_1.itram25, (
	    ftnlen *)0, 5 };
    static Vardesc itram26_dv = { "ITRAM26", (char *)&inputcom_1.itram26, (
	    ftnlen *)0, 5 };
    static Vardesc itram31_dv = { "ITRAM31", (char *)&inputcom_1.itram31, (
	    ftnlen *)0, 5 };
    static Vardesc itram32_dv = { "ITRAM32", (char *)&inputcom_1.itram32, (
	    ftnlen *)0, 5 };
    static Vardesc itram33_dv = { "ITRAM33", (char *)&inputcom_1.itram33, (
	    ftnlen *)0, 5 };
    static Vardesc itram34_dv = { "ITRAM34", (char *)&inputcom_1.itram34, (
	    ftnlen *)0, 5 };
    static Vardesc itram35_dv = { "ITRAM35", (char *)&inputcom_1.itram35, (
	    ftnlen *)0, 5 };
    static Vardesc itram36_dv = { "ITRAM36", (char *)&inputcom_1.itram36, (
	    ftnlen *)0, 5 };
    static Vardesc itram41_dv = { "ITRAM41", (char *)&inputcom_1.itram41, (
	    ftnlen *)0, 5 };
    static Vardesc itram42_dv = { "ITRAM42", (char *)&inputcom_1.itram42, (
	    ftnlen *)0, 5 };
    static Vardesc itram43_dv = { "ITRAM43", (char *)&inputcom_1.itram43, (
	    ftnlen *)0, 5 };
    static Vardesc itram44_dv = { "ITRAM44", (char *)&inputcom_1.itram44, (
	    ftnlen *)0, 5 };
    static Vardesc itram45_dv = { "ITRAM45", (char *)&inputcom_1.itram45, (
	    ftnlen *)0, 5 };
    static Vardesc itram46_dv = { "ITRAM46", (char *)&inputcom_1.itram46, (
	    ftnlen *)0, 5 };
    static Vardesc itram51_dv = { "ITRAM51", (char *)&inputcom_1.itram51, (
	    ftnlen *)0, 5 };
    static Vardesc itram52_dv = { "ITRAM52", (char *)&inputcom_1.itram52, (
	    ftnlen *)0, 5 };
    static Vardesc itram53_dv = { "ITRAM53", (char *)&inputcom_1.itram53, (
	    ftnlen *)0, 5 };
    static Vardesc rmax0sc_dv = { "RMAX0SC", (char *)&inputcom_1.rmax0sc, (
	    ftnlen *)0, 5 };
    static Vardesc itram54_dv = { "ITRAM54", (char *)&inputcom_1.itram54, (
	    ftnlen *)0, 5 };
    static Vardesc outputfile_dv = { "OUTPUTFILE", inputcom_1.outputfile, (
	    ftnlen *)0, -30 };
    static Vardesc lbc_dv = { "LBC", (char *)&inputcom_1.lbc, (ftnlen *)0, 3 }
	    ;
    static Vardesc itram55_dv = { "ITRAM55", (char *)&inputcom_1.itram55, (
	    ftnlen *)0, 5 };
    static Vardesc itram56_dv = { "ITRAM56", (char *)&inputcom_1.itram56, (
	    ftnlen *)0, 5 };
    static Vardesc itram61_dv = { "ITRAM61", (char *)&inputcom_1.itram61, (
	    ftnlen *)0, 5 };
    static Vardesc itram62_dv = { "ITRAM62", (char *)&inputcom_1.itram62, (
	    ftnlen *)0, 5 };
    static Vardesc itram63_dv = { "ITRAM63", (char *)&inputcom_1.itram63, (
	    ftnlen *)0, 5 };
    static Vardesc itram64_dv = { "ITRAM64", (char *)&inputcom_1.itram64, (
	    ftnlen *)0, 5 };
    static Vardesc itram65_dv = { "ITRAM65", (char *)&inputcom_1.itram65, (
	    ftnlen *)0, 5 };
    static Vardesc itram66_dv = { "ITRAM66", (char *)&inputcom_1.itram66, (
	    ftnlen *)0, 5 };
    static Vardesc awd_dv = { "AWD", (char *)&inputcom_1.awd, (ftnlen *)0, 5 }
	    ;
    static Vardesc ibfield_dv = { "IBFIELD", (char *)&inputcom_1.ibfield, (
	    ftnlen *)0, 5 };
    static Vardesc drl_dv = { "DRL", (char *)&inputcom_1.drl, (ftnlen *)0, 5 }
	    ;
    static Vardesc radfile_dv = { "RADFILE", inputcom_1.radfile, (ftnlen *)0, 
	    -30 };
    static Vardesc igamref_dv = { "IGAMREF", (char *)&inputcom_1.igamref, (
	    ftnlen *)0, 5 };
    static Vardesc idmpfld_dv = { "IDMPFLD", (char *)&inputcom_1.idmpfld, (
	    ftnlen *)0, 3 };
    static Vardesc awx_dv = { "AWX", (char *)&inputcom_1.awx, (ftnlen *)0, 5 }
	    ;
    static Vardesc awy_dv = { "AWY", (char *)&inputcom_1.awy, (ftnlen *)0, 5 }
	    ;
    static Vardesc iharmsc_dv = { "IHARMSC", (char *)&inputcom_1.iharmsc, (
	    ftnlen *)0, 3 };
    static Vardesc curpeak_dv = { "CURPEAK", (char *)&inputcom_1.curpeak, (
	    ftnlen *)0, 5 };
    static Vardesc xkx_dv = { "XKX", (char *)&inputcom_1.xkx, (ftnlen *)0, 5 }
	    ;
    static Vardesc xky_dv = { "XKY", (char *)&inputcom_1.xky, (ftnlen *)0, 5 }
	    ;
    static Vardesc inverfc_dv = { "INVERFC", (char *)&inputcom_1.inverfc, (
	    ftnlen *)0, 3 };
    static Vardesc idmppar_dv = { "IDMPPAR", (char *)&inputcom_1.idmppar, (
	    ftnlen *)0, 3 };
    static Vardesc cuttail_dv = { "CUTTAIL", (char *)&inputcom_1.cuttail, (
	    ftnlen *)0, 5 };
    static Vardesc conditx_dv = { "CONDITX", (char *)&inputcom_1.conditx, (
	    ftnlen *)0, 5 };
    static Vardesc condity_dv = { "CONDITY", (char *)&inputcom_1.condity, (
	    ftnlen *)0, 5 };
    static Vardesc iscrkup_dv = { "ISCRKUP", (char *)&inputcom_1.iscrkup, (
	    ftnlen *)0, 3 };
    static Vardesc version_dv = { "VERSION", (char *)&inputcom_1.version, (
	    ftnlen *)0, 5 };
    static Vardesc f1st_dv = { "F1ST", (char *)&inputcom_1.f1st, (ftnlen *)0, 
	    5 };
**/
/**  //OC030110
    static Vardesc *newrun_vl[] = { &aw0_dv, &xkx_dv, &xky_dv, &wcoefz_dv, &
	    xlamd_dv, &fbess0_dv, &delaw_dv, &iertyp_dv, &iwityp_dv, &awd_dv, 
	    &iseed_dv, &awx_dv, &awy_dv, &npart_dv, &gamma0_dv, &delgam_dv, &
	    rxbeam_dv, &rybeam_dv, &alphax_dv, &alphay_dv, &emitx_dv, &
	    emity_dv, &xbeam_dv, &ybeam_dv, &pxbeam_dv, &pybeam_dv, &
	    cuttail_dv, &curpeak_dv, &conditx_dv, &condity_dv, &bunch_dv, &
	    bunchphase_dv, &emod_dv, &emodphase_dv, &xlamds_dv, &prad0_dv, &
	    zrayl_dv, &zwaist_dv, &pradh0_dv, &ncar_dv, &lbc_dv, &rmax0_dv, &
	    dgrid_dv, &nscr_dv, &nscz_dv, &nptr_dv, &rmax0sc_dv, &iscrkup_dv, 
	    &nwig_dv, &delz_dv, &zsep_dv, &nsec_dv, &iorb_dv, &zstop_dv, &
	    magin_dv, &magout_dv, &nbins_dv, &version_dv, &quadf_dv, &
	    quadd_dv, &fl_dv, &dl_dv, &drl_dv, &f1st_dv, &qfdx_dv, &qfdy_dv, &
	    sl_dv, &solen_dv, &ildgam_dv, &ildpsi_dv, &ildx_dv, &ildy_dv, &
	    ildpx_dv, &ildpy_dv, &itgaus_dv, &lout_dv, &igamgaus_dv, &
	    inverfc_dv, &iphsty_dv, &ishsty_dv, &ippart_dv, &ispart_dv, &
	    ipradi_dv, &isradi_dv, &idump_dv, &iotail_dv, &nharm_dv, &
	    iallharm_dv, &iharmsc_dv, &idmppar_dv, &idmpfld_dv, &ilog_dv, &
	    ffspec_dv, &beamfile_dv, &fieldfile_dv, &maginfile_dv, &
	    magoutfile_dv, &outputfile_dv, &partfile_dv, &distfile_dv, &
	    filetype_dv, &radfile_dv, &curlen_dv, &ntail_dv, &nslice_dv, &
	    shotnoise_dv, &iall_dv, &itdp_dv, &ipseed_dv, &isntyp_dv, &
	    iscan_dv, &nscan_dv, &svar_dv, &scan_dv, &isravg_dv, &isrsig_dv, &
	    eloss_dv, &ndcut_dv, &ibfield_dv, &imagl_dv, &idril_dv, &
	    convharm_dv, &alignradf_dv, &offsetradf_dv, &multconv_dv, &
	    igamref_dv, &trama_dv, &itram11_dv, &itram12_dv, &itram13_dv, &
	    itram14_dv, &itram15_dv, &itram16_dv, &itram21_dv, &itram22_dv, &
	    itram23_dv, &itram24_dv, &itram25_dv, &itram26_dv, &itram31_dv, &
	    itram32_dv, &itram33_dv, &itram34_dv, &itram35_dv, &itram36_dv, &
	    itram41_dv, &itram42_dv, &itram43_dv, &itram44_dv, &itram45_dv, &
	    itram46_dv, &itram51_dv, &itram52_dv, &itram53_dv, &itram54_dv, &
	    itram55_dv, &itram56_dv, &itram61_dv, &itram62_dv, &itram63_dv, &
	    itram64_dv, &itram65_dv, &itram66_dv };
**/
    //static Namelist newrun = { "NEWRUN", newrun_vl, 163 }; //OC030110

    /* Fortran I/O blocks */
    //static cilist io___18 = { 0, 6, 0, fmt_100, 0 }; //OC030110
    //static cilist io___19 = { 0, 5, 0, fmt_110, 0 };
    //static cilist io___20 = { 1, 0, 1, (char *)&newrun, 0 };


/*     ================================================================= */
/*     this routine reads in the user input files. */
/*     it assumes the standard fortran namelist format. */
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

/*     ------------------------------------------------------------------ */
/*     common block for mpi implementation for genesis */






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
/*     temporary included parameter */

/*     initialize input/output */

    nin = 8;
    ret_val = 0;
    preset_();
    if (mpicom_1.mpi_id__ == 0) {

	//s_wsfe(&io___18); //OC port
	//e_wsfe(); //OC port

/* initialize input p */
	//s_rsfe(&io___19); //OC port
	//do_fio(&c__1, inputcom_1.inputfile, (ftnlen)30); //OC port
	//e_rsfe(); //OC port

/* get input filename */
    }
    if (mpicom_1.mpi_size__ > 1) {
	mpi_bcast__(inputcom_1.inputfile, &c__30, &c__1, &c__0, &c__1, &
		mpicom_1.mpi_err__, (ftnlen)30);
    }

    o__1.oerr = 1;
    o__1.ounit = nin;
    o__1.ofnmlen = 30;
    o__1.ofnm = inputcom_1.inputfile;
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;

	//OC port
    //i__1 = f_open(&o__1);
    //if (i__1 != 0) {
	//goto L10;
    //}
	//OC port end

/* open input file. */
	//OC port
    //io___20.ciunit = nin;
    //i__1 = s_rsne(&io___20);
    //if (i__1 != 0) {
	//goto L20;
    //}
	//OC port end

/* read in namel */
	//OC port
    //cl__1.cerr = 0;
    //cl__1.cunit = nin;
    //cl__1.csta = 0;
    //f_clos(&cl__1);
	//OC port end

/* close file */
    return ret_val;

//OC port
//L10:
//    ret_val = printerr_(&c_n3, inputcom_1.inputfile, (ftnlen)30);
//    return ret_val;
//L20:
//    ret_val = printerr_(&c_n4, inputcom_1.inputfile, (ftnlen)30);
//    last_();
//    return ret_val;
//OC port end

/*     format statements */

} /* readin_ */


/* readin */
/* Subroutine */ int readbeamfile_(file, file_len)
char *file;
ftnlen file_len;
{
    /* Format strings */
    static char fmt_100[] = "(a)";
    static char fmt_110[] = "(\002Auto-adjustment of time window:\002,/,\002\
nslice=\002,i6,/\002ntail =\002,i6)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    cllist cl__1;

    /* Builtin functions */
    integer i_indx(), s_rsfe(), do_fio(), e_rsfe(), s_wsli(), do_lio(), 
	    e_wsli(), i_dnnt();
    double sqrt();
    integer f_clos(), s_wsfe(), e_wsfe();

    /* Local variables */
    static integer ncol, ipar[15], itmp;
    static char cerr__[50], line[511];
    static doublereal tmin, tmax;
    extern /* Subroutine */ int last_(), getfirstchar_();
    extern integer printerr_(), opentextfile_();
    static integer i__, j, idata, ft, ix;
    extern integer extractval_(), detectfiletype_();
    static integer nin;
    static doublereal ver, values[15], reverse;
    extern /* Subroutine */ int getbeamfileinfo_(), touppercase_();

    /* Fortran I/O blocks */
    static cilist io___30 = { 1, 0, 1, fmt_100, 0 };
    static icilist io___35 = { 0, cerr__, 0, 0, 50, 1 };
    static icilist io___36 = { 0, cerr__, 0, 0, 50, 1 };
    static cilist io___40 = { 0, 0, 0, fmt_110, 0 };


/*     ============================================================= */
/*     read the file for external description of the electron beam */
/*     ------------------------------------------------------------- */





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





    itmp = 0;
    tbunchcom_1.ndata = -1;
    if (i_indx(file, " ", file_len, (ftnlen)1) == 1) {
	return 0;
    }

    if (inputcom_1.iscan > 0 && inputcom_1.iscan < 23) {
	i__ = printerr_(&c_n6, file, file_len);
	return 0;
    }

/*     read file */

    ft = detectfiletype_(file, file_len);

/* check for filetype */
    nin = opentextfile_(file, "old", &c__8, file_len, (ftnlen)3);
    if (nin < 0) {
	last_();
    }

/* stop program on error */
    tbunchcom_1.ndata = -1;
/* # of rows not defined */
    idata = 0;
/* # of rows read */
    ncol = 15;
/* # of elements per line */
    reverse = (float)1.;
/* tail for ZPOS and head for TPOS comes first i */
    ver = (float).1;
    i__1 = ncol;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipar[i__ - 1] = i__;
/* basic order of input */
    }

L1:
    io___30.ciunit = nin;
    i__1 = s_rsfe(&io___30);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = do_fio(&c__1, line, (ftnlen)511);
    if (i__1 != 0) {
	goto L100001;
    }
    i__1 = e_rsfe();
L100001:
    if (i__1 < 0) {
	goto L50;
    }
    if (i__1 > 0) {
	goto L20;
    }

/*     processing line */

    touppercase_(line, (ftnlen)511);
    getfirstchar_(line, &ix, (ftnlen)511);

    if (ix == 0) {
	goto L1;
    }
/* empty line */
    if (*(unsigned char *)&line[ix - 1] == '#') {
	goto L1;
    }
/* no comment used */
    if (*(unsigned char *)&line[ix - 1] == '?') {
	getbeamfileinfo_(line, ipar, &ncol, &tbunchcom_1.ndata, &ver, &
		reverse, (ftnlen)511);
/* check */
	goto L1;
    }

    if (tbunchcom_1.ndata < 0 && ver < (float)1.) {
/* old version */
	i__ = extractval_(line, values, &c__1, (ftnlen)511);
	if (i__ < 0) {
	    s_wsli(&io___35);
	    do_lio(&c__9, &c__1, "Line number of BEAMFILE cannot be determin\
ed", (ftnlen)44);
	    e_wsli();
	    i__ = printerr_(&c_n8, cerr__, (ftnlen)50);
	}
	tbunchcom_1.ndata = i_dnnt(values);
	goto L1;
    }

    i__ = extractval_(line, values, &ncol, (ftnlen)511);
    if (i__ < 0) {
	s_wsli(&io___36);
	do_lio(&c__9, &c__1, "BEAMFILE data line ", (ftnlen)19);
	i__1 = idata + 1;
	do_lio(&c__3, &c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " has bad format", (ftnlen)15);
	e_wsli();
	i__ = printerr_(&c_n8, cerr__, (ftnlen)50);
	last_();
    }

    ++idata;

/*     set default values */

    tbunchcom_1.tpos[idata - 1] = inputcom_1.xlamds * inputcom_1.zsep * idata;
/* can we make this the default? */
    tbunchcom_1.tgam0[idata - 1] = inputcom_1.gamma0;
    tbunchcom_1.tdgam[idata - 1] = inputcom_1.delgam;
    tbunchcom_1.temitx[idata - 1] = inputcom_1.emitx;
    tbunchcom_1.temity[idata - 1] = inputcom_1.emity;
    tbunchcom_1.txrms[idata - 1] = inputcom_1.rxbeam;
    tbunchcom_1.tyrms[idata - 1] = inputcom_1.rybeam;
    tbunchcom_1.txpos[idata - 1] = inputcom_1.xbeam;
    tbunchcom_1.typos[idata - 1] = inputcom_1.ybeam;
    tbunchcom_1.tpxpos[idata - 1] = inputcom_1.pxbeam;
    tbunchcom_1.tpypos[idata - 1] = inputcom_1.pybeam;
    tbunchcom_1.talphx[idata - 1] = inputcom_1.alphax;
    tbunchcom_1.talphy[idata - 1] = inputcom_1.alphay;
    tbunchcom_1.tcurrent[idata - 1] = inputcom_1.curpeak;
    tbunchcom_1.tloss[idata - 1] = inputcom_1.eloss;

/*     write over with input data */

    i__1 = ncol;
    for (j = 1; j <= i__1; ++j) {
	if (ipar[j - 1] == 1) {
	    tbunchcom_1.tpos[idata - 1] = reverse * values[j - 1];
	}
	if (ipar[j - 1] == -1) {
	    tbunchcom_1.tpos[idata - 1] = -reverse * values[j - 1] * (float)
		    3e8;
	}
/* ti */
	if (ipar[j - 1] == 2) {
	    tbunchcom_1.tgam0[idata - 1] = values[j - 1];
	}
	if (ipar[j - 1] == 3) {
	    tbunchcom_1.tdgam[idata - 1] = values[j - 1];
	}
	if (ipar[j - 1] == 4) {
	    tbunchcom_1.temitx[idata - 1] = values[j - 1];
	}
	if (ipar[j - 1] == 5) {
	    tbunchcom_1.temity[idata - 1] = values[j - 1];
	}
	if ((i__2 = ipar[j - 1], abs(i__2)) == 6) {
	    tbunchcom_1.txrms[idata - 1] = values[j - 1];
	}
/* save for d */
	if ((i__2 = ipar[j - 1], abs(i__2)) == 7) {
	    tbunchcom_1.tyrms[idata - 1] = values[j - 1];
	}
/* RXBEAM/BET */
	if (ipar[j - 1] == 8) {
	    tbunchcom_1.txpos[idata - 1] = values[j - 1];
	}
	if (ipar[j - 1] == 9) {
	    tbunchcom_1.typos[idata - 1] = values[j - 1];
	}
	if ((i__2 = ipar[j - 1], abs(i__2)) == 10) {
	    tbunchcom_1.tpxpos[idata - 1] = values[j - 1];
	}
/* save for d */
	if ((i__2 = ipar[j - 1], abs(i__2)) == 11) {
	    tbunchcom_1.tpypos[idata - 1] = values[j - 1];
	}
/* PXBEAM/XPR */
	if (ipar[j - 1] == 12) {
	    tbunchcom_1.talphx[idata - 1] = values[j - 1];
	}
	if (ipar[j - 1] == 13) {
	    tbunchcom_1.talphy[idata - 1] = values[j - 1];
	}
	if (ipar[j - 1] == 14) {
	    tbunchcom_1.tcurrent[idata - 1] = values[j - 1];
	}
	if (ipar[j - 1] == 15) {
	    tbunchcom_1.tloss[idata - 1] = values[j - 1];
	}
    }

/*     check for unphysical parameters */

    if (tbunchcom_1.tgam0[idata - 1] - (d__1 = tbunchcom_1.tdgam[idata - 1], 
	    abs(d__1)) * 4 < 1.) {
	itmp = printerr_(&c_n8, "Energy GAMMA0 too small in BEAMFILE", (
		ftnlen)35);
    } else {
	i__1 = ncol;
	for (j = 1; j <= i__1; ++j) {
/* calculate beam sizes (avods floating poi */
	    if (ipar[j - 1] == -6) {
		tbunchcom_1.txrms[idata - 1] = sqrt(tbunchcom_1.txrms[idata - 
			1] * tbunchcom_1.temitx[idata - 1] / 
			tbunchcom_1.tgam0[idata - 1]);
	    }
	    if (ipar[j - 1] == -7) {
		tbunchcom_1.tyrms[idata - 1] = sqrt(tbunchcom_1.tyrms[idata - 
			1] * tbunchcom_1.temity[idata - 1] / 
			tbunchcom_1.tgam0[idata - 1]);
	    }
	    if (ipar[j - 1] == -10) {
		tbunchcom_1.tpxpos[idata - 1] *= tbunchcom_1.tgam0[idata - 1];
	    }
	    if (ipar[j - 1] == -11) {
		tbunchcom_1.tpypos[idata - 1] *= tbunchcom_1.tgam0[idata - 1];
	    }
	}
    }

    if (tbunchcom_1.tcurrent[idata - 1] < 0.) {
	itmp = printerr_(&c_n8, "Current negative in BEAMFILE", (ftnlen)28);
/* abort */
    }

    if (tbunchcom_1.temitx[idata - 1] < 0.) {
	itmp = printerr_(&c_n8, "EMITX negative in BEAMFILE", (ftnlen)26);
/* abort */
    }

    if (tbunchcom_1.temity[idata - 1] < 0.) {
	itmp = printerr_(&c_n8, "EMITY negative in BEAMFILE", (ftnlen)26);
/* abort */
    }

    if (idata == 1) {
	tmin = tbunchcom_1.tpos[idata - 1];
	tmax = tbunchcom_1.tpos[idata - 1];
    } else {
	if (tbunchcom_1.tpos[idata - 1] > tmax) {
	    tmax = tbunchcom_1.tpos[idata - 1];
	}
	if (tbunchcom_1.tpos[idata - 1] < tmin) {
	    tmin = tbunchcom_1.tpos[idata - 1];
	}
    }

    goto L1;

L50:
    cl__1.cerr = 0;
    cl__1.cunit = nin;
    cl__1.csta = 0;
    f_clos(&cl__1);

    if (tbunchcom_1.ndata >= 0 && idata != tbunchcom_1.ndata) {
	i__ = printerr_(&c_n9, "BEAMFILE has fewer lines than defined", (
		ftnlen)37);
    }
    tbunchcom_1.ndata = idata;
    if (idata < 2) {
	i__ = printerr_(&c_n8, "BEAMFILE contains less than 2 valid lines", (
		ftnlen)41);
	last_();
    }
    if (itmp != 0) {
	last_();
    }

    if (ver >= (float)1.) {
	i__1 = tbunchcom_1.ndata;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tbunchcom_1.tpos[i__ - 1] -= tmin;
/* set time window to zero */
	}
    }

    if (inputcom_1.nslice <= 0) {
	inputcom_1.nslice = (integer) ((tmax - tmin) / inputcom_1.xlamds / 
		inputcom_1.zsep);
	if (ver >= (float)1.) {
	    inputcom_1.ntail = 0;
	} else {
	    inputcom_1.ntail = (integer) (tmin / inputcom_1.xlamds / 
		    inputcom_1.zsep);
	}
	io___40.ciunit = inputcom_1.ilog;
	s_wsfe(&io___40);
	do_fio(&c__1, (char *)&inputcom_1.nslice, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&inputcom_1.ntail, (ftnlen)sizeof(integer));
	e_wsfe();
    }

    return 0;

L20:
    i__ = printerr_(&c_n2, file, file_len);
    cl__1.cerr = 0;
    cl__1.cunit = nin;
    cl__1.csta = 0;
    f_clos(&cl__1);
    last_();
    return 0;

/*     format statement */


} /* readbeamfile_ */


/* readbeamfile */
/* Subroutine */ int getbeamfileinfo_(line, ipar, ncol, ndata, ver, reverse, 
	line_len)
char *line;
integer *ipar, *ncol, *ndata;
doublereal *ver, *reverse;
ftnlen line_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_indx(), i_len(), i_dnnt();
    /* Subroutine */ int s_copy();

    /* Local variables */
    static integer iarg, ierr;
    extern /* Subroutine */ int getfirstchar_();
    extern integer printerr_();
    static integer i__, j, n;
    static char cline[511];
    extern integer extractnumber_();
    static integer ix1, ix2;
    static doublereal val;
    static integer haszpos;

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


/*     version number */

    /* Parameter adjustments */
    --ipar;

    /* Function Body */
    i__ = i_indx(line, "VERSION", line_len, (ftnlen)7);
/* check for version number */
    n = i_len(line, line_len);
    if (i__ > 0) {
	i__1 = i__ + 6;
	ierr = extractnumber_(line + i__1, &val, n - i__1);
	if (ierr < 0) {
	    i__ = printerr_(&c_n9, "Unrecognized information line in BEAMFILE"
		    , (ftnlen)41);
	} else {
	    *ver = val;
	}
	return 0;
    }

/*     reverse order */

    i__ = i_indx(line, "REVERSE", line_len, (ftnlen)7);
    n = i_len(line, line_len);
    if (i__ > 0) {
	*reverse = (float)-1.;
	return 0;
    }

/*     line numbers */

    i__ = i_indx(line, "SIZE", line_len, (ftnlen)4);
/* check for size argument (aka ndata) */
    if (i__ > 0) {
	i__1 = i__ + 3;
	ierr = extractnumber_(line + i__1, &val, n - i__1);
	if (ierr < 0) {
	    i__ = printerr_(&c_n9, "Unrecognized information line in BEAMFILE"
		    , (ftnlen)41);
	} else {
	    *ndata = i_dnnt(&val);
	}
	return 0;
    }

/*     colum order */

    i__ = i_indx(line, "COLUMNS", line_len, (ftnlen)7);
/* check for colums headers */
    if (i__ > 0) {
	for (j = 1; j <= 15; ++j) {
	    ipar[j] = 0;
	}
	*ncol = 0;
	i__1 = i__ + 6;
	s_copy(cline, line + i__1, (ftnlen)511, n - i__1);
	haszpos = 0;
L1:
	getfirstchar_(cline, &ix1, (ftnlen)511);
	if (ix1 > 0) {
	    ix2 = 255;
/* search backwards */
	    i__1 = ix1 + 1;
	    for (j = 255; j >= i__1; --j) {
/* for first space after ix1 */
		if (*(unsigned char *)&cline[j - 1] == ' ') {
		    ix2 = j;
		}
	    }
	    iarg = 0;
	    if (i_indx(cline + (ix1 - 1), "ZPOS", ix2 - (ix1 - 1), (ftnlen)4) 
		    != 0) {
		iarg = 1;
	    }
	    if (i_indx(cline + (ix1 - 1), "TPOS", ix2 - (ix1 - 1), (ftnlen)4) 
		    != 0) {
		iarg = -1;
	    }
	    if (i_indx(cline + (ix1 - 1), "GAMMA0", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = 2;
	    }
	    if (i_indx(cline + (ix1 - 1), "DELGAM", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = 3;
	    }
	    if (i_indx(cline + (ix1 - 1), "EMITX", ix2 - (ix1 - 1), (ftnlen)5)
		     != 0) {
		iarg = 4;
	    }
	    if (i_indx(cline + (ix1 - 1), "EMITY", ix2 - (ix1 - 1), (ftnlen)5)
		     != 0) {
		iarg = 5;
	    }
	    if (i_indx(cline + (ix1 - 1), "XBEAM", ix2 - (ix1 - 1), (ftnlen)5)
		     != 0) {
		iarg = 8;
	    }
/* XBEAM can */
	    if (i_indx(cline + (ix1 - 1), "YBEAM", ix2 - (ix1 - 1), (ftnlen)5)
		     != 0) {
		iarg = 9;
	    }
/* Thus compa */
	    if (i_indx(cline + (ix1 - 1), "RXBEAM", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = 6;
	    }
	    if (i_indx(cline + (ix1 - 1), "RYBEAM", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = 7;
	    }
	    if (i_indx(cline + (ix1 - 1), "BETAX", ix2 - (ix1 - 1), (ftnlen)5)
		     != 0) {
		iarg = -6;
	    }
	    if (i_indx(cline + (ix1 - 1), "BETAY", ix2 - (ix1 - 1), (ftnlen)5)
		     != 0) {
		iarg = -7;
	    }
	    if (i_indx(cline + (ix1 - 1), "PXBEAM", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = 10;
	    }
	    if (i_indx(cline + (ix1 - 1), "PYBEAM", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = 11;
	    }
	    if (i_indx(cline + (ix1 - 1), "XPRIME", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = -10;
	    }
	    if (i_indx(cline + (ix1 - 1), "YPRIME", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = -11;
	    }
	    if (i_indx(cline + (ix1 - 1), "ALPHAX", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = 12;
	    }
	    if (i_indx(cline + (ix1 - 1), "ALPHAY", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = 13;
	    }
	    if (i_indx(cline + (ix1 - 1), "CURPEAK", ix2 - (ix1 - 1), (ftnlen)
		    7) != 0) {
		iarg = 14;
	    }
	    if (i_indx(cline + (ix1 - 1), "ELOSS", ix2 - (ix1 - 1), (ftnlen)5)
		     != 0) {
		iarg = 15;
	    }

	    if (iarg == 0) {
		for (j = 1; j <= 15; ++j) {
		    ipar[j] = j;
		}
		*ncol = 15;
		j = printerr_(&c_n9, "Unrecognized information line in BEAMF\
ILE", (ftnlen)41);
		return 0;
	    } else {
		++(*ncol);
		ipar[*ncol] = iarg;
		if (abs(iarg) == 1) {
		    haszpos = 1;
		}
		i__1 = ix2;
		s_copy(cline, cline + i__1, (ftnlen)511, 255 - i__1);
		if (*ncol < 15) {
		    goto L1;
		}
	    }
	}
	if (haszpos < 1) {
	    for (j = 1; j <= 15; ++j) {
		ipar[j] = j;
	    }
	    *ncol = 15;
	    j = printerr_(&c_n9, "ZPOS/TPOS column not specified in BEAMFILE",
		     (ftnlen)42);
	    return 0;
	}
	return 0;
    }
    i__ = printerr_(&c_n9, "Unrecognized information line in BEAMFILE", (
	    ftnlen)41);
    return 0;
} /* getbeamfileinfo_ */


integer readfield_(cin, irec)
doublecomplex *cin;
integer *irec;
{
    /* System generated locals */
    integer ret_val, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    double sqrt();
    integer s_rdue(), do_uio(), e_rdue();

    /* Local variables */
    extern integer printerr_();
    static integer ix;
    static doublereal scltmp;

    /* Fortran I/O blocks */
    static cilist io___52 = { 1, 0, 0, 0, 0 };


/*     ============================================== */
/*     read field from input file */
/*     ---------------------------------------------- */





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
/*     input/output control */





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




/*     simulation control and normalisation parameter */



    /* Parameter adjustments */
    --cin;

    /* Function Body */
    scltmp = cartcom_1.xks * sqrt(376.73) / (cartcom_1.dxy * 510999.06 * 
	    simcom_1.xkper0);

    ret_val = 0;
    io___52.ciunit = iocom_1.nfin;
    io___52.cirec = *irec;
    i__1 = s_rdue(&io___52);
    if (i__1 != 0) {
	goto L1;
    }
    i__2 = inputcom_1.ncar * inputcom_1.ncar;
    for (ix = 1; ix <= i__2; ++ix) {
	i__1 = do_uio(&c__2, (char *)&cin[ix], (ftnlen)sizeof(doublereal));
	if (i__1 != 0) {
	    goto L1;
	}
    }
    i__1 = e_rdue();
    if (i__1 != 0) {
	goto L1;
    }
    i__1 = inputcom_1.ncar * inputcom_1.ncar;
    for (ix = 1; ix <= i__1; ++ix) {
	i__2 = ix;
	i__3 = ix;
	z__1.r = scltmp * cin[i__3].r, z__1.i = scltmp * cin[i__3].i;
	cin[i__2].r = z__1.r, cin[i__2].i = z__1.i;
    }

    return ret_val;
L1:
    ret_val = printerr_(&c_n19, inputcom_1.fieldfile, (ftnlen)30);
    return ret_val;
} /* readfield_ */



/* Subroutine */ int importdispersion_()
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12, d__13, d__14, d__15;

    /* Builtin functions */
    double asin(), tan(), sin(), cos();

    /* Local variables */
    static integer ierr;
    extern integer printerr_();
    static integer i__;
    static doublereal iarho, iaphi, imagl_old__, ypart_old__, mam, py_old__, 
	    ma12, ma33, ma34, ma43, ma56;

/*     ================================================================= */
/*     apply dispersion to imported beam file from readpart */
/*     subroutine supplied by Atoosa Meseck from Bessy */
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





    if (inputcom_1.ibfield == 0.) {
	return 0;
    }
    if (inputcom_1.idril < 0.) {
	ierr = printerr_(&c_n9, "IDRIL<0-NO DISP.SECTION", (ftnlen)23);
	return 0;
    }

    if (inputcom_1.igamref <= 0.) {
	inputcom_1.igamref = inputcom_1.gamma0;
    }

    inputcom_1.ibfield = abs(inputcom_1.ibfield);

    imagl_old__ = inputcom_1.imagl;
    iarho = inputcom_1.igamref * (float)5.11e-4 / (inputcom_1.ibfield * (
	    float).299793);
    iaphi = asin(inputcom_1.imagl / iarho);
    mam = tan(iaphi) / iarho;
    inputcom_1.imagl = iaphi * iarho;
    ma12 = inputcom_1.idril * 3 + iarho * 4 * sin(iaphi) * cos(iaphi) + 
	    inputcom_1.idril * 2 * cos(iaphi) * cos(iaphi);
/* Computing 3rd power */
    d__1 = mam;
/* Computing 2nd power */
    d__2 = inputcom_1.imagl;
/* Computing 2nd power */
    d__3 = inputcom_1.idril;
/* Computing 3rd power */
    d__4 = inputcom_1.idril;
/* Computing 4th power */
    d__5 = mam, d__5 *= d__5;
/* Computing 3rd power */
    d__6 = inputcom_1.idril;
/* Computing 2nd power */
    d__7 = inputcom_1.idril * 2 * inputcom_1.imagl;
/* Computing 4th power */
    d__8 = inputcom_1.idril, d__8 *= d__8;
    ma33 = mam * (inputcom_1.idril * -10 - inputcom_1.imagl * 8) + mam * mam *
	     (inputcom_1.imagl * 26 * inputcom_1.idril + inputcom_1.idril * 
	    15 * inputcom_1.idril + inputcom_1.imagl * 8 * inputcom_1.imagl) 
	    + d__1 * (d__1 * d__1) * (inputcom_1.idril * -12 * (d__2 * d__2) 
	    - d__3 * d__3 * 20 * inputcom_1.imagl - d__4 * (d__4 * d__4) * 7) 
	    + 1 + d__5 * d__5 * (d__6 * (d__6 * d__6) * 4 * inputcom_1.imagl 
	    + d__7 * d__7 + d__8 * d__8);
/* ma44=ma33 */
/* Computing 2nd power */
    d__1 = inputcom_1.imagl;
/* Computing 2nd power */
    d__2 = inputcom_1.idril;
/* Computing 2nd power */
    d__3 = inputcom_1.idril;
/* Computing 3rd power */
    d__4 = inputcom_1.idril;
/* Computing 2nd power */
    d__5 = inputcom_1.imagl;
/* Computing 3rd power */
    d__6 = mam;
/* Computing 2nd power */
    d__7 = inputcom_1.idril * inputcom_1.imagl;
/* Computing 3rd power */
    d__8 = inputcom_1.idril;
/* Computing 4th power */
    d__9 = inputcom_1.idril, d__9 *= d__9;
/* Computing 4th power */
    d__10 = mam, d__10 *= d__10;
/* Computing 3rd power */
    d__11 = inputcom_1.idril;
/* Computing 2nd power */
    d__12 = inputcom_1.imagl;
/* Computing 4th power */
    d__13 = inputcom_1.idril, d__13 *= d__13;
/* Computing 5th power */
    d__14 = inputcom_1.idril, d__15 = d__14, d__14 *= d__14;
    ma34 = mam * (inputcom_1.idril * -28 * inputcom_1.imagl - d__1 * d__1 * 8 
	    - d__2 * d__2 * 20) + mam * mam * (inputcom_1.imagl * 44 * (d__3 *
	     d__3) + d__4 * (d__4 * d__4) * 21 + inputcom_1.idril * 20 * (
	    d__5 * d__5)) + d__6 * (d__6 * d__6) * (d__7 * d__7 * -16 - d__8 *
	     (d__8 * d__8) * 24 * inputcom_1.imagl - d__9 * d__9 * 8) + d__10 
	    * d__10 * (d__11 * (d__11 * d__11) * 4 * (d__12 * d__12) + d__13 *
	     d__13 * 4 * inputcom_1.imagl + d__15 * (d__14 * d__14)) + 
	    inputcom_1.idril * 5 + inputcom_1.imagl * 4;
/* Computing 2nd power */
    d__1 = iarho;
/* Computing 2nd power */
    d__2 = tan(iaphi) * inputcom_1.idril;
/* Computing 4th power */
    d__3 = iarho, d__3 *= d__3;
    ma43 = tan(iaphi) * (iarho * -2 + tan(iaphi) * 2 * inputcom_1.imagl + tan(
	    iaphi) * inputcom_1.idril) * (d__1 * d__1 * 2 - tan(iaphi) * 4 * 
	    inputcom_1.imagl * iarho - tan(iaphi) * 4 * inputcom_1.idril * 
	    iarho + inputcom_1.idril * 2 * tan(iaphi) * tan(iaphi) * 
	    inputcom_1.imagl + d__2 * d__2) / (d__3 * d__3);
/* Computing 2nd power */
    d__1 = inputcom_1.igamref;
    ma56 = iarho * 8 * sin(iaphi) - iarho * 4 * sin(iaphi) * cos(iaphi) + 
	    inputcom_1.idril * 2 - inputcom_1.idril * 2 * cos(iaphi) * cos(
	    iaphi) - inputcom_1.imagl * 4 + (inputcom_1.idril * 5 + 
	    inputcom_1.imagl * 4) / (d__1 * d__1);

    inputcom_1.imagl = imagl_old__;

    i__1 = inputcom_1.npart;
    for (i__ = 1; i__ <= i__1; ++i__) {
	beamcom_1.theta[i__ - 1] += ma56 * (beamcom_1.gamma[i__ - 1] - 
		inputcom_1.igamref) / inputcom_1.igamref * 6.28318530717958 / 
		inputcom_1.xlamds / inputcom_1.convharm;
	beamcom_1.xpart[i__ - 1] += ma12 * beamcom_1.px[i__ - 1] / 
		beamcom_1.gamma[i__ - 1];
	ypart_old__ = beamcom_1.ypart[i__ - 1];
	py_old__ = beamcom_1.py[i__ - 1];
	beamcom_1.ypart[i__ - 1] = ma33 * ypart_old__ + ma34 * py_old__ / 
		beamcom_1.gamma[i__ - 1];
	beamcom_1.py[i__ - 1] = ma43 * ypart_old__ * beamcom_1.gamma[i__ - 1] 
		+ ma33 * py_old__;
    }
    return 0;
} /* importdispersion_ */



/* of import dispersion */
/* Subroutine */ int importtransfer_()
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal ypart_old__, xpart_old__, gamma_old__, theta_old__, 
	    py_old__, px_old__;

/*     ================================================================= */
/*     Transfer matrix calculation supplied by A. Meseck. */
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

/*     ------------------------------------------------------------------ */
/*     input/output control */





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



/*     simulation control and normalisation parameter */



/*     ------------------------------------------------------------------ */
/*     electron beam */





    if ((doublereal) inputcom_1.trama == 0.) {
	return 0;
    }
    if (inputcom_1.igamref <= 0.) {
	inputcom_1.igamref = inputcom_1.gamma0;
    }
    i__1 = inputcom_1.npart;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*  Denormalize */
	beamcom_1.px[i__ - 1] /= beamcom_1.gamma[i__ - 1];
	beamcom_1.py[i__ - 1] /= beamcom_1.gamma[i__ - 1];
/* c */
	xpart_old__ = beamcom_1.xpart[i__ - 1];
	ypart_old__ = beamcom_1.ypart[i__ - 1];
	px_old__ = beamcom_1.px[i__ - 1];
	py_old__ = beamcom_1.py[i__ - 1];
	theta_old__ = beamcom_1.theta[i__ - 1];
	gamma_old__ = beamcom_1.gamma[i__ - 1];
	beamcom_1.xpart[i__ - 1] = inputcom_1.itram11 * xpart_old__ + 
		inputcom_1.itram12 * px_old__ + inputcom_1.itram13 * 
		ypart_old__ + inputcom_1.itram14 * py_old__ + 
		inputcom_1.itram15 * beamcom_1.theta[i__ - 1] * 
		inputcom_1.xlamds * inputcom_1.convharm / 6.28318530717958 + 
		inputcom_1.itram16 * (beamcom_1.gamma[i__ - 1] - 
		inputcom_1.igamref) / inputcom_1.igamref;
	beamcom_1.px[i__ - 1] = inputcom_1.itram21 * xpart_old__ + 
		inputcom_1.itram22 * px_old__ + inputcom_1.itram23 * 
		ypart_old__ + inputcom_1.itram24 * py_old__ + 
		inputcom_1.itram25 * theta_old__ * inputcom_1.xlamds * 
		inputcom_1.convharm / 6.28318530717958 + inputcom_1.itram26 * 
		(beamcom_1.gamma[i__ - 1] - inputcom_1.igamref) / 
		inputcom_1.igamref;
	beamcom_1.ypart[i__ - 1] = inputcom_1.itram31 * xpart_old__ + 
		inputcom_1.itram32 * px_old__ + inputcom_1.itram33 * 
		ypart_old__ + inputcom_1.itram34 * py_old__ + 
		inputcom_1.itram35 * theta_old__ * inputcom_1.xlamds * 
		inputcom_1.convharm / 6.28318530717958 + inputcom_1.itram36 * 
		(beamcom_1.gamma[i__ - 1] - inputcom_1.igamref) / 
		inputcom_1.igamref;
	beamcom_1.py[i__ - 1] = inputcom_1.itram41 * xpart_old__ + 
		inputcom_1.itram42 * px_old__ + inputcom_1.itram43 * 
		ypart_old__ + inputcom_1.itram44 * py_old__ + 
		inputcom_1.itram45 * theta_old__ * inputcom_1.xlamds * 
		inputcom_1.convharm / 6.28318530717958 + inputcom_1.itram46 * 
		(beamcom_1.gamma[i__ - 1] - inputcom_1.igamref) / 
		inputcom_1.igamref;
	beamcom_1.theta[i__ - 1] = inputcom_1.itram55 * theta_old__ + 
		inputcom_1.itram56 * ((beamcom_1.gamma[i__ - 1] - 
		inputcom_1.igamref) / inputcom_1.igamref) * 6.28318530717958 /
		 inputcom_1.xlamds / inputcom_1.convharm + (
		inputcom_1.itram51 * xpart_old__ + inputcom_1.itram52 * 
		px_old__ + inputcom_1.itram53 * ypart_old__ + 
		inputcom_1.itram54 * py_old__) * 6.28318530717958 / 
		inputcom_1.xlamds / inputcom_1.convharm;
	beamcom_1.gamma[i__ - 1] = (inputcom_1.itram61 * xpart_old__ + 
		inputcom_1.itram62 * px_old__ + inputcom_1.itram63 * 
		ypart_old__ + inputcom_1.itram64 * py_old__ + 
		inputcom_1.itram65 * theta_old__ * inputcom_1.xlamds * 
		inputcom_1.convharm / 6.28318530717958) * inputcom_1.igamref 
		+ inputcom_1.itram66 * (beamcom_1.gamma[i__ - 1] - 
		inputcom_1.igamref) + inputcom_1.igamref;
/* normalization */
	beamcom_1.px[i__ - 1] *= beamcom_1.gamma[i__ - 1];
	beamcom_1.py[i__ - 1] *= beamcom_1.gamma[i__ - 1];
/* c */
    }
    return 0;
} /* importtransfer_ */



/* of import transfermatrix */
integer readpart_(islice)
integer *islice;
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_rdue(), do_uio(), e_rdue();
    double sqrt();

    /* Local variables */
    static integer idel;
    extern /* Subroutine */ int last_();
    extern integer printerr_();
    static integer i__, j;
    extern /* Subroutine */ int importdispersion_(), importtransfer_();

    /* Fortran I/O blocks */
    static cilist io___75 = { 1, 0, 0, 0, 0 };
    static cilist io___77 = { 1, 0, 0, 0, 0 };
    static cilist io___78 = { 1, 0, 0, 0, 0 };
    static cilist io___79 = { 1, 0, 0, 0, 0 };
    static cilist io___80 = { 1, 0, 0, 0, 0 };
    static cilist io___81 = { 1, 0, 0, 0, 0 };


/*     ================================================================= */
/*     load complete set of particle from file */
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




/*     ------------------------------------------------------------------ */
/*     input/output control */





/*     simulation control and normalisation parameter */




    ret_val = 0;

    if (inputcom_1.multconv == 0) {
	j = (*islice - 1) * 6 + 1;
/* reading every slice */
    } else {
	j = (*islice - 1) / inputcom_1.convharm * 6 + 1;
/* reading every (1/convharm)th sli */
    }

    io___75.ciunit = iocom_1.npin;
    io___75.cirec = j;
    i__1 = s_rdue(&io___75);
    if (i__1 != 0) {
	goto L100;
    }
    i__2 = simcom_1.npart0;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = do_uio(&c__1, (char *)&beamcom_1.gamma[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L100;
	}
    }
    i__1 = e_rdue();
    if (i__1 != 0) {
	goto L100;
    }
    io___77.ciunit = iocom_1.npin;
    io___77.cirec = j + 1;
    i__1 = s_rdue(&io___77);
    if (i__1 != 0) {
	goto L100;
    }
    i__2 = simcom_1.npart0;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = do_uio(&c__1, (char *)&beamcom_1.theta[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L100;
	}
    }
    i__1 = e_rdue();
    if (i__1 != 0) {
	goto L100;
    }
    io___78.ciunit = iocom_1.npin;
    io___78.cirec = j + 2;
    i__1 = s_rdue(&io___78);
    if (i__1 != 0) {
	goto L100;
    }
    i__2 = simcom_1.npart0;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = do_uio(&c__1, (char *)&beamcom_1.xpart[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L100;
	}
    }
    i__1 = e_rdue();
    if (i__1 != 0) {
	goto L100;
    }
    io___79.ciunit = iocom_1.npin;
    io___79.cirec = j + 3;
    i__1 = s_rdue(&io___79);
    if (i__1 != 0) {
	goto L100;
    }
    i__2 = simcom_1.npart0;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = do_uio(&c__1, (char *)&beamcom_1.ypart[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (i__1 != 0) {
	    goto L100;
	}
    }
    i__1 = e_rdue();
    if (i__1 != 0) {
	goto L100;
    }
    io___80.ciunit = iocom_1.npin;
    io___80.cirec = j + 4;
    i__1 = s_rdue(&io___80);
    if (i__1 != 0) {
	goto L100;
    }
    i__2 = simcom_1.npart0;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = do_uio(&c__1, (char *)&beamcom_1.px[i__ - 1], (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L100;
	}
    }
    i__1 = e_rdue();
    if (i__1 != 0) {
	goto L100;
    }
    io___81.ciunit = iocom_1.npin;
    io___81.cirec = j + 5;
    i__1 = s_rdue(&io___81);
    if (i__1 != 0) {
	goto L100;
    }
    i__2 = simcom_1.npart0;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = do_uio(&c__1, (char *)&beamcom_1.py[i__ - 1], (ftnlen)sizeof(
		doublereal));
	if (i__1 != 0) {
	    goto L100;
	}
    }
    i__1 = e_rdue();
    if (i__1 != 0) {
	goto L100;
    }

/*     apply transfermatrix to particle distribution */
    importtransfer_();

/*     dispersive section */
    importdispersion_();


/*     convert to higher harmonic */

    if (inputcom_1.convharm > 1) {
	i__1 = inputcom_1.npart;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    beamcom_1.theta[i__ - 1] = (real) inputcom_1.convharm * 
		    beamcom_1.theta[i__ - 1];
	}
    }

/*     calculate init. perpendicular velocity (needed in first call of track) */

    i__1 = simcom_1.npart0;
    for (i__ = 1; i__ <= i__1; ++i__) {
	beamcom_1.xpart[i__ - 1] *= simcom_1.xkper0;
	beamcom_1.ypart[i__ - 1] *= simcom_1.xkper0;
/* Computing 2nd power */
	d__1 = beamcom_1.px[i__ - 1];
/* Computing 2nd power */
	d__2 = beamcom_1.py[i__ - 1];
/* Computing 2nd power */
	d__3 = beamcom_1.gamma[i__ - 1];
	beamcom_1.btpar[i__ - 1] = sqrt(1. - (d__1 * d__1 + d__2 * d__2 + (
		float)1.) / (d__3 * d__3));
/* parallel velocity */
    }

/*     check for particle losses from previous run */

    idel = 0;
    i__1 = inputcom_1.npart;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (beamcom_1.gamma[i__ - 1] > (float)0.) {
	    beamcom_1.gamma[i__ - idel - 1] = beamcom_1.gamma[i__ - 1];
	    beamcom_1.theta[i__ - idel - 1] = beamcom_1.theta[i__ - 1];
	    beamcom_1.xpart[i__ - idel - 1] = beamcom_1.xpart[i__ - 1];
	    beamcom_1.ypart[i__ - idel - 1] = beamcom_1.ypart[i__ - 1];
	    beamcom_1.px[i__ - idel - 1] = beamcom_1.px[i__ - 1];
	    beamcom_1.py[i__ - idel - 1] = beamcom_1.py[i__ - 1];
	} else {
	    ++idel;
	}
    }
    inputcom_1.npart = simcom_1.npart0 - idel;
    beamcom_1.xcuren = beamcom_1.xcuren * (real) inputcom_1.npart / (real) 
	    simcom_1.npart0;

    return ret_val;
L100:
    ret_val = printerr_(&c_n2, inputcom_1.partfile, (ftnlen)30);
    last_();
    return ret_val;
} /* readpart_ */


/* Subroutine */ int readslice_(nget, x, px, y, py, g, tmin, tmax)
integer *nget;
doublereal *x, *px, *y, *py, *g, *tmin, *tmax;
{
    /* Format strings */
    static char fmt_200[] = "(a)";

    /* System generated locals */
    integer i__1;
    alist al__1;

    /* Builtin functions */
    integer s_rsfe(), do_fio(), e_rsfe(), s_wsli(), do_lio(), e_wsli(), 
	    i_dnnt();
    double sqrt();
    integer f_rew();

    /* Local variables */
    static char line[255], cerr__[255];
    extern /* Subroutine */ int last_(), getfirstchar_();
    extern integer printerr_();
    static integer i__, n0, ip;
    static doublereal tg;
    static integer ix;
    static doublereal tt, tx, ty, values[10];
    extern integer extractval_();
    static doublereal tpx, tpy;
    extern /* Subroutine */ int touppercase_();

    /* Fortran I/O blocks */
    static cilist io___85 = { 1, 0, 1, fmt_200, 0 };
    static icilist io___91 = { 0, cerr__, 0, 0, 255, 1 };
    static icilist io___92 = { 0, cerr__, 0, 0, 255, 1 };


/*     ================================================================= */
/*     read slice [tmin,tmax] of particle from distribution file */
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
/*     time dependency of bunch + slippage field */




/*     ------------------------------------------------------------------ */
/*     input/output control */





    /* Parameter adjustments */
    --g;
    --py;
    --y;
    --px;
    --x;

    /* Function Body */
    *nget = 0;

    if (TRUE_) {
	i__1 = tbunchcom_1.ndist;
	for (ip = 1; ip <= i__1; ++ip) {
	    if (tbunchcom_1.distt[ip - 1] >= *tmin && tbunchcom_1.distt[ip - 
		    1] <= *tmax) {
/* in sli */
		++(*nget);
		x[*nget] = tbunchcom_1.distx[ip - 1];
/* add to raw distr */
		px[*nget] = tbunchcom_1.distpx[ip - 1];
		y[*nget] = tbunchcom_1.disty[ip - 1];
		py[*nget] = tbunchcom_1.distpy[ip - 1];
		g[*nget] = tbunchcom_1.distgam[ip - 1];
	    }
	}
	return 0;
    }

    n0 = -1;

L1:
    io___85.ciunit = iocom_1.ndis;
    i__1 = s_rsfe(&io___85);
    if (i__1 != 0) {
	goto L100002;
    }
    i__1 = do_fio(&c__1, line, (ftnlen)255);
    if (i__1 != 0) {
	goto L100002;
    }
    i__1 = e_rsfe();
L100002:
    if (i__1 < 0) {
	goto L50;
    }
    if (i__1 > 0) {
	goto L100;
    }

/*     processing line */

    touppercase_(line, (ftnlen)255);
    getfirstchar_(line, &ix, (ftnlen)255);
    if (ix == 0) {
	goto L1;
    }
/* empty line */
    if (*(unsigned char *)&line[ix - 1] == '#') {
	goto L1;
    }
/* comment line */
    if (*(unsigned char *)&line[ix - 1] == '?') {
	goto L1;
    }

/* information allready processed */
    if (n0 < 0 && iocom_1.distversion < (float)1.) {
	i__ = extractval_(line, values, &c__1, (ftnlen)255);
	if (i__ < 0) {
	    s_wsli(&io___91);
	    do_lio(&c__9, &c__1, "DISTFILE has invalid input line: ", (ftnlen)
		    33);
	    do_lio(&c__9, &c__1, line, (ftnlen)255);
	    e_wsli();
	    i__ = printerr_(&c_n8, cerr__, (ftnlen)255);
	    last_();
	}
	n0 = i_dnnt(values);
	goto L1;
    }

/*     get record-tuple */

    i__ = extractval_(line, values, &iocom_1.ncoldis, (ftnlen)255);
    if (i__ < 0) {
	s_wsli(&io___92);
	do_lio(&c__9, &c__1, "DISTFILE has invalid input line: ", (ftnlen)33);
	do_lio(&c__9, &c__1, line, (ftnlen)255);
	e_wsli();
	i__ = printerr_(&c_n8, cerr__, (ftnlen)255);
	last_();
    }
    i__1 = iocom_1.ncoldis;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iocom_1.icolpar[i__ - 1] == 1) {
	    tx = values[i__ - 1];
	}
	if (iocom_1.icolpar[i__ - 1] == 2) {
	    tpx = values[i__ - 1];
	}
	if (iocom_1.icolpar[i__ - 1] == 3) {
	    ty = values[i__ - 1];
	}
	if (iocom_1.icolpar[i__ - 1] == 4) {
	    tpy = values[i__ - 1];
	}
	if (iocom_1.icolpar[i__ - 1] == 5) {
	    tt = iocom_1.distrev * values[i__ - 1];
	}
	if (iocom_1.icolpar[i__ - 1] == 6) {
	    tg = values[i__ - 1];
	}
    }
    if (iocom_1.iconv2t != 0) {
	tt = -tt / (float)3e8;
    }
/* convert from z to t */
    if (iocom_1.iconv2g != 0) {
	tg = sqrt(tg * tg + 1.);
    }
/* convert from p to gamma */
    if (iocom_1.iconv2px != 0) {
	tpx *= tg;
    }
/* convert from x' to px */
    if (iocom_1.iconv2py != 0) {
	tpy *= tg;
    }

/* convert from y' to py */
    if (tt >= *tmin && tt <= *tmax) {
/* in slice */
	++(*nget);
	x[*nget] = tx;
/* add to raw distribution */
	px[*nget] = tpx;
	y[*nget] = ty;
	py[*nget] = tpy;
	g[*nget] = tg;
    }
    goto L1;

L50:
    al__1.aerr = 0;
    al__1.aunit = iocom_1.ndis;
    f_rew(&al__1);
/* set file pointer back */
    return 0;

L100:
    i__ = printerr_(&c_n2, "DISTFILE", (ftnlen)8);
    last_();
    return 0;


} /* readslice_ */


integer readdistfile_(file, nio, file_len)
char *file;
integer *nio;
ftnlen file_len;
{
    /* Format strings */
    static char fmt_200[] = "(a)";

    /* System generated locals */
    address a__1[2];
    integer ret_val, i__1, i__2[2];
    cllist cl__1;
    alist al__1;

    /* Builtin functions */
    integer i_indx(), s_rsfe(), do_fio(), e_rsfe();
    /* Subroutine */ int s_cat();
    integer i_dnnt(), f_rew(), f_clos();
    double sqrt();

    /* Local variables */
    static char cerr__[255], line[255];
    static integer nget;
    extern /* Subroutine */ int last_(), getfirstchar_();
    extern integer printerr_(), opentextfile_();
    static integer i__, ip, ix;
    extern integer extractval_(), detectfiletype_();
    static integer niotmp;
    static doublereal tt, values[10];
    extern /* Subroutine */ int touppercase_(), getdistfileinfo_();

    /* Fortran I/O blocks */
    static cilist io___102 = { 1, 0, 1, fmt_200, 0 };


/*     ================================================================== */
/*     open an external file containing the distribution. read the file */
/*     geting the parameter range etc. */
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

/*     ------------------------------------------------------------------ */
/*     electron beam */





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





    iocom_1.ndistsize = -1;
    nget = 0;

/*     default settings */

    iocom_1.iconv2g = 1;
/* needs to convert from p to gamma */
    iocom_1.iconv2t = 0;
/* long coordinate is time */
    iocom_1.iconv2px = 0;
/* convert from x' to px */
    iocom_1.iconv2py = 0;
/* convert from y' to py */
    iocom_1.distversion = (float).1;
/* first number is number of particles */
    iocom_1.distrev = (float)1.;
/* normal order of long. position */
    iocom_1.ncoldis = 6;
/* 6 dimension */
    for (i__ = 1; i__ <= 10; ++i__) {
	iocom_1.icolpar[i__ - 1] = i__;
/* x,px,y,py,t,gamma,id order */
    }

    ret_val = -1;
    if (i_indx(file, " ", file_len, (ftnlen)1) == 1) {
	return ret_val;
    }

/* no file selected */
    iocom_1.ftdist = detectfiletype_(file, file_len);

/* sven some compiler gives an error here (segmentation fault) */

/* check for filetype */
    niotmp = opentextfile_(file, "old", nio, file_len, (ftnlen)3);
    ret_val = niotmp;
    if (niotmp < 0) {
	last_();
    }

L1:
    io___102.ciunit = *nio;
    i__1 = s_rsfe(&io___102);
    if (i__1 != 0) {
	goto L100003;
    }
    i__1 = do_fio(&c__1, line, (ftnlen)255);
    if (i__1 != 0) {
	goto L100003;
    }
    i__1 = e_rsfe();
L100003:
    if (i__1 < 0) {
	goto L50;
    }
    if (i__1 > 0) {
	goto L100;
    }

/*     processing line */

    touppercase_(line, (ftnlen)255);
    getfirstchar_(line, &ix, (ftnlen)255);
    if (ix == 0) {
	goto L1;
    }
/* empty line */
    if (*(unsigned char *)&line[ix - 1] == '#') {
	goto L1;
    }
/* comment line */
    if (*(unsigned char *)&line[ix - 1] == '?') {
/* read information line */
	getdistfileinfo_(line, (ftnlen)255);
	goto L1;
    }

    if (iocom_1.ndistsize < 0 && iocom_1.distversion < (float)1.) {
	i__ = extractval_(line, values, &c__1, (ftnlen)255);
	if (i__ < 0) {
/* Writing concatenation */
	    i__2[0] = 32, a__1[0] = "DISTFILE has invalid input line:";
	    i__2[1] = 255, a__1[1] = line;
	    s_cat(cerr__, a__1, i__2, &c__2, (ftnlen)255);
	    i__ = printerr_(&c_n8, cerr__, (ftnlen)255);
	    last_();
	}
	iocom_1.ndistsize = i_dnnt(values);
	goto L1;
    }

/*     get record-tuple */

    i__ = extractval_(line, values, &iocom_1.ncoldis, (ftnlen)255);
    if (i__ < 0) {
/* Writing concatenation */
	i__2[0] = 32, a__1[0] = "distfile has invalid input line:";
	i__2[1] = 255, a__1[1] = line;
	s_cat(cerr__, a__1, i__2, &c__2, (ftnlen)255);
	i__ = printerr_(&c_n8, cerr__, (ftnlen)255);
	last_();
    }

    ++nget;
    if (TRUE_) {
/* save variabls if kept in memory */
	if (nget > 1250000) {
/* check for overflow */
	    i__ = printerr_(&c_n8, "DISTFILE size exceeds NDMAX", (ftnlen)27);
	    last_();
	}
	i__1 = iocom_1.ncoldis;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (iocom_1.icolpar[i__ - 1] == 1) {
		tbunchcom_1.distx[nget - 1] = values[i__ - 1];
	    }
	    if (iocom_1.icolpar[i__ - 1] == 2) {
		tbunchcom_1.distpx[nget - 1] = values[i__ - 1];
	    }
	    if (iocom_1.icolpar[i__ - 1] == 3) {
		tbunchcom_1.disty[nget - 1] = values[i__ - 1];
	    }
	    if (iocom_1.icolpar[i__ - 1] == 4) {
		tbunchcom_1.distpy[nget - 1] = values[i__ - 1];
	    }
	    if (iocom_1.icolpar[i__ - 1] == 5) {
		tbunchcom_1.distt[nget - 1] = iocom_1.distrev * values[i__ - 
			1];
	    }
	    if (iocom_1.icolpar[i__ - 1] == 6) {
		tbunchcom_1.distgam[nget - 1] = values[i__ - 1];
	    }
	}
    }
    i__1 = iocom_1.ncoldis;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iocom_1.icolpar[i__ - 1] == 5) {
	    tt = iocom_1.distrev * values[i__ - 1];
	}
/* catch time value */
    }
    if (iocom_1.iconv2t != 0) {
	tt = -tt / (float)3e8;
    }

/* convert from z to t */
    if (nget == 1) {
	beamcom_1.tdmin = tt;
	beamcom_1.tdmax = tt;
    }
    if (tt < beamcom_1.tdmin) {
	beamcom_1.tdmin = tt;
    }
/* adjust min and max */
    if (tt > beamcom_1.tdmax) {
	beamcom_1.tdmax = tt;
    }
    goto L1;

L50:
    al__1.aerr = 0;
    al__1.aunit = *nio;
    f_rew(&al__1);
/* go back to file beginning */
    if (TRUE_) {
	cl__1.cerr = 0;
	cl__1.cunit = *nio;
	cl__1.csta = 0;
	f_clos(&cl__1);
	tbunchcom_1.ndist = nget;
	if (iocom_1.iconv2t != 0) {
	    i__1 = nget;
	    for (ip = 1; ip <= i__1; ++ip) {
		tbunchcom_1.distt[ip - 1] = -tbunchcom_1.distt[ip - 1] / (
			float)3e8;
/* convert from space coordinate */
	    }
	}
	if (iocom_1.iconv2g != 0) {
	    i__1 = nget;
	    for (ip = 1; ip <= i__1; ++ip) {
		tbunchcom_1.distgam[ip - 1] = sqrt(tbunchcom_1.distgam[ip - 1]
			 * tbunchcom_1.distgam[ip - 1] + 1.);
/* convert */
	    }
	}
	if (iocom_1.iconv2px != 0) {
	    i__1 = nget;
	    for (ip = 1; ip <= i__1; ++ip) {
		tbunchcom_1.distpx[ip - 1] *= tbunchcom_1.distgam[ip - 1];
/* convert from x' to px */
	    }
	}
	if (iocom_1.iconv2py != 0) {
	    i__1 = nget;
	    for (ip = 1; ip <= i__1; ++ip) {
		tbunchcom_1.distpy[ip - 1] *= tbunchcom_1.distgam[ip - 1];
/* convert from y' to py */
	    }
	}
    }

/*     set time window */

    if (inputcom_1.nslice <= 0) {
	inputcom_1.nslice = (integer) ((beamcom_1.tdmax - beamcom_1.tdmin) * 
		3e8 / inputcom_1.xlamds / inputcom_1.zsep);
	inputcom_1.ntail = 0;
    }

    if (iocom_1.ndistsize >= 0 && nget != iocom_1.ndistsize) {
	i__ = printerr_(&c_n9, "DISTFILE has fewer lines than defined", (
		ftnlen)37);
    }

    if (inputcom_1.ndcut <= 0) {
	inputcom_1.ndcut = nget * inputcom_1.nharm / inputcom_1.npart;
    }
/* self optimizing */
    if (inputcom_1.ndcut <= 0) {
	inputcom_1.ndcut = 1;
    }

    if (beamcom_1.charge <= 0.) {
	i__ = printerr_(&c_n8, "CHARGE for DISTFILE is not defined", (ftnlen)
		34);
	last_();
    }

    beamcom_1.delcharge = beamcom_1.charge / (doublereal) nget;
    beamcom_1.dtd = (beamcom_1.tdmax - beamcom_1.tdmin) / (doublereal) 
	    inputcom_1.ndcut;

    return ret_val;

L100:
    ret_val = printerr_(&c_n2, file, file_len);
    last_();
    return ret_val;


} /* readdistfile_ */


/* Subroutine */ int getdistfileinfo_(line, line_len)
char *line;
ftnlen line_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len(), i_indx(), i_dnnt();
    /* Subroutine */ int s_copy();

    /* Local variables */
    static integer iarg, narg, ierr;
    extern /* Subroutine */ int last_(), getfirstchar_();
    extern integer printerr_();
    static integer i__, j, n;
    static char cline[255];
    extern integer extractnumber_();
    static integer ix1, ix2, ncount;
    static doublereal val;

/*     ================================================================= */
/*     extract information from distfile */
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





/*     ------------------------------------------------------------------ */
/*     input/output control */





    narg = 0;

/* some compile initialize it to unity. */
    n = i_len(line, line_len);
/*      call getfirstchar(line,idx) ! get first character should be identical to index !! */

/*     version number */

    i__ = i_indx(line, "VERSION", line_len, (ftnlen)7);
/* check for version number */
    if (i__ > 0) {
	i__1 = i__ + 6;
	ierr = extractnumber_(line + i__1, &val, n - i__1);
	if (ierr < 0) {
	    i__ = printerr_(&c_n9, "Unrecognized information line in distfile"
		    , (ftnlen)41);
	} else {
	    iocom_1.distversion = val;
	}
	return 0;
    }

/*     order of distributon */

    i__ = i_indx(line, "REVERSE", line_len, (ftnlen)7);
    if (i__ > 0) {
	iocom_1.distrev = (float)-1.;
	return 0;
    }

/*     beam charge */

    i__ = i_indx(line, "CHARGE", line_len, (ftnlen)6);
    if (i__ > 0) {
	i__1 = i__ + 5;
	ierr = extractnumber_(line + i__1, &val, n - i__1);
	if (ierr < 0) {
	    i__ = printerr_(&c_n9, "Unrecognized information line in distfile"
		    , (ftnlen)41);
	} else {
	    beamcom_1.charge = val;
	}
	return 0;
    }

/*     cuts in long. phase space */

/*     disabled - no program specific instruction should be present */


/*      i=index(line,'ndcut') */
/*      if (i.gt.0) then */
/*        ierr=extractnumber(line(i+5:n),val) */
/*        if (ierr.lt.0) then */
/*           i=printerr(errinwarn, */
/*     c               'unrecognized information line in distfile') */
/*        else */
/*           ndcut=nint(val) */
/*        endif */
/*        return */
/*      endif */

/*     size of record */

    i__ = i_indx(line, "SIZE", line_len, (ftnlen)4);
    if (i__ > 0) {
	i__1 = i__ + 3;
	ierr = extractnumber_(line + i__1, &val, n - i__1);
	if (ierr < 0) {
	    i__ = printerr_(&c_n9, "Unrecognized information line in distfile"
		    , (ftnlen)41);
	} else {
	    iocom_1.ndistsize = i_dnnt(&val);
	}
	return 0;
    }

/*     get order and dimension */

    i__ = i_indx(line, "COLUMNS", line_len, (ftnlen)7);
/* check for colums headers */
    if (i__ > 0) {
	for (j = 1; j <= 10; ++j) {
	    iocom_1.icolpar[j - 1] = 0;
	}
	iocom_1.ncoldis = 0;
	ncount = 0;
	i__1 = i__ + 6;
	s_copy(cline, line + i__1, (ftnlen)255, i_len(line, line_len) - i__1);
L1:
	getfirstchar_(cline, &ix1, (ftnlen)255);
	if (ix1 > 0) {
	    ix2 = 255;
/* search backwards */
	    i__1 = ix1 + 1;
	    for (j = 255; j >= i__1; --j) {
/* for first space after ix1 */
		if (*(unsigned char *)&cline[j - 1] == ' ') {
		    ix2 = j;
		}
	    }
	    iarg = 0;
	    if (i_indx(cline + (ix1 - 1), "P", ix2 - (ix1 - 1), (ftnlen)1) != 
		    0) {
		iarg = 6;
	    }
	    if (i_indx(cline + (ix1 - 1), "GAMMA", ix2 - (ix1 - 1), (ftnlen)5)
		     != 0) {
		iarg = -6;
	    }
	    if (i_indx(cline + (ix1 - 1), "X", ix2 - (ix1 - 1), (ftnlen)1) != 
		    0) {
		iarg = 1;
	    }
	    if (i_indx(cline + (ix1 - 1), "PX", ix2 - (ix1 - 1), (ftnlen)2) !=
		     0) {
		iarg = 2;
	    }
	    if (i_indx(cline + (ix1 - 1), "XPRIME", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = -2;
	    }
	    if (i_indx(cline + (ix1 - 1), "Y", ix2 - (ix1 - 1), (ftnlen)1) != 
		    0) {
		iarg = 3;
	    }
	    if (i_indx(cline + (ix1 - 1), "PY", ix2 - (ix1 - 1), (ftnlen)2) !=
		     0) {
		iarg = 4;
	    }
	    if (i_indx(cline + (ix1 - 1), "YPRIME", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = -4;
	    }
	    if (i_indx(cline + (ix1 - 1), "T", ix2 - (ix1 - 1), (ftnlen)1) != 
		    0) {
		iarg = 5;
	    }
	    if (i_indx(cline + (ix1 - 1), "Z", ix2 - (ix1 - 1), (ftnlen)1) != 
		    0) {
		iarg = -5;
	    }

	    ++ncount;
	    narg += abs(iarg);
	    ++iocom_1.ncoldis;
	    iocom_1.icolpar[iocom_1.ncoldis - 1] = abs(iarg);
	    if (iarg == 2) {
		iocom_1.iconv2px = 0;
	    }
	    if (iarg == -2) {
		iocom_1.iconv2px = 1;
	    }
	    if (iarg == 4) {
		iocom_1.iconv2py = 0;
	    }
	    if (iarg == -4) {
		iocom_1.iconv2py = 1;
	    }
	    if (iarg == 5) {
		iocom_1.iconv2t = 0;
	    }
	    if (iarg == -5) {
		iocom_1.iconv2t = 1;
	    }
	    if (iarg == 6) {
		iocom_1.iconv2g = 1;
	    }
	    if (iarg == -6) {
		iocom_1.iconv2g = 0;
	    }
	    i__1 = ix2;
	    s_copy(cline, cline + i__1, (ftnlen)255, 255 - i__1);
	    if (ncount < 10) {
		goto L1;
	    }
	}
	if (narg != 21) {
	    j = printerr_(&c_n9, "Not all dimensions defined in DISTFILE", (
		    ftnlen)38);
	    last_();
	}
	return 0;
    }

/*     unrecognized */

    i__ = printerr_(&c_n9, "Unrecognized information line in distfile", (
	    ftnlen)41);
    return 0;
} /* getdistfileinfo_ */


/* Subroutine */ int preset_()
{
    /* Builtin functions */
    /* Subroutine */ int s_copy();

    /* Local variables */
    static integer i__;

/*     ================================================================== */
/*     sets default values of program inputs. */
/*     default parameter modelled after the pegasus fel */
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




/*     wiggler parameters: */

    inputcom_1.aw0 = (float).735;
/* dimensionless wiggler amplitude (rms). */
    inputcom_1.xkx = 0.;
/* weak focusing strength in x-plane. */
    inputcom_1.xky = 1.;
/* weak focusing strength in x-plane. */
    inputcom_1.wcoefz[0] = 0.;
/* field taper start location in z. */
    inputcom_1.wcoefz[1] = 0.;
/* field taper gradient */
    inputcom_1.wcoefz[2] = 0.;
/* field taper model */
    inputcom_1.delaw = 0.;
/* relative spread in aw */
    inputcom_1.iertyp = 0;
/* field error distribution type (<0 = corre */
    inputcom_1.iwityp = 0;
/* wiggler type (0=planar, 1= wiggler) */
    inputcom_1.awd = (float).735;
/* virtual wiggler field for drift space */
    inputcom_1.iseed = -1;
/* initial seed for wiggler error generation */
    inputcom_1.fbess0 = 0.;
/* beam-radiation coupling */
    inputcom_1.xlamd = .0205;
/* wiggler period  (m) */
    inputcom_1.awx = (float)0.;
/* max offset in x for undulator misalignment */
    inputcom_1.awy = (float)0.;

/*     electron beam parameters: */

/* max offset in y for undulator misalignment */
    inputcom_1.npart = 8192;
/* number of simulation particles. */
    inputcom_1.gamma0 = 35.2;
/* electron beam lorentz factor (energy) */
    inputcom_1.delgam = .005;
/* rms energy (gamma) spread */
    inputcom_1.rxbeam = 1.121e-4;
/* rms beam size in x-plane */
    inputcom_1.rybeam = 1.121e-4;
/* rms beam size in y-plane */
    inputcom_1.alphax = 0.;
/* twiss alpha parameter in x */
    inputcom_1.alphay = 0.;
/* twiss alpha parameter in y */
    inputcom_1.emitx = 2e-6;
/* normalized emittance x-plane (pi m-rad) */
    inputcom_1.emity = 2e-6;
/* normalized emittance y-plane (pi m-rad) */
    inputcom_1.xbeam = 0.;
/* center in x (m) */
    inputcom_1.ybeam = 0.;
/* center in y (m) */
    inputcom_1.pxbeam = 0.;
/* center in px */
    inputcom_1.pybeam = 0.;
/* center in py */
    inputcom_1.cuttail = (float)-1.;
/* no collimation transverse tail */
    inputcom_1.curpeak = 250.;
/* peak current */
    inputcom_1.conditx = 0.;
/* conditioning in x plane (1/m) */
    inputcom_1.condity = 0.;
/* conditioning in y plane (1/m) */
    inputcom_1.bunch = 0.;
/* prebunching, fraction */
    inputcom_1.bunchphase = 0.;
/* phase for prebunching */
    inputcom_1.emod = 0.;
/* energy modulation (in gamma) */
    inputcom_1.emodphase = 0.;

/*     radiation: */

/* phase for energy modulation */
    inputcom_1.xlamds = 1.2852e-5;
/* output radiation wavelength */
    inputcom_1.prad0 = 10.;
/* input power */
    inputcom_1.pradh0 = 0.;
/* input power of higher harmonics. */
    inputcom_1.zrayl = .5;
/* rayleigh range (of radiation) */
    inputcom_1.zwaist = 0.;
/* z location of (radiation) waist. */
    cartcom_1.radphase = 0.;

/*     numerical control parameters: */

/* although not an input parameter, here is the best */
    inputcom_1.ildgam = 5;
/* energy loading parameter. */
    inputcom_1.ildpsi = 7;
/* phase loading parameter. */
    inputcom_1.ildx = 1;
/* x-plane loading parameter. */
    inputcom_1.ildy = 2;
/* y-plane loading parameter. */
    inputcom_1.ildpx = 3;
/* x momentum loading parameter. */
    inputcom_1.ildpy = 4;
/* y momentum loading parameter. */
    inputcom_1.itgaus = 1;
/* gaussian (<>0) or uniform (=0) loading */
    inputcom_1.nbins = 4;
/* # of bins in the phase coordinate */
    inputcom_1.igamgaus = 1;
/* gaussian (<>0) or uniform (=0) enegy load */
    inputcom_1.inverfc = 0;

/*     mesh discretization: */

/* <>0 uses inverse error function to load G */
    inputcom_1.ncar = 151;
/* mesh points in one dimension (xy-grid). */
    inputcom_1.lbc = 0;
/* boundary condition (0=diriqlet,<>0 neuman */
    inputcom_1.rmax0 = 9.;
/* maximum edge of grid. */
    inputcom_1.nscr = 0;
/* # radial modes for space charge */
    inputcom_1.nscz = 0;
/* # longitudinal modes for space charge */
    inputcom_1.nptr = 40;
/* radial grid points for space charge */
    inputcom_1.dgrid = 0.;
/* grid size(-dgrid to dgrid) if dgrid > 0 */
    inputcom_1.rmax0sc = 0.;
/* if <0 the grid for space charge calculati */
    inputcom_1.iscrkup = 0;

/*     integration control parameter: */

/* if <> 0 then space charge is calculated a */
    inputcom_1.nwig = 98;
/* number of undulator periods per module */
    inputcom_1.zsep = (float)1.;
/* seperation of slices in units of xlamds */
    inputcom_1.delz = (float)1.;
/* integration step in units of xlamd */
    inputcom_1.nsec = 1;
/* number of sections */
    inputcom_1.iorb = 0;
/* orbit correction flag (<>0 -> use term). */
    inputcom_1.zstop = (float)-1.;
/* stop simulation @ z=zstop */
    inputcom_1.magin = 0;
/* read magnetic field description if magin */
    inputcom_1.magout = 0;
/* write magnetic field description if magou */
    inputcom_1.version = (float).1;

/*     scan control parameter */

/* assume oldest version */
    inputcom_1.iscan = 0;
/* >0 -> iscanth parameter selected */
    inputcom_1.nscan = 3;
/* number of scans */
    inputcom_1.svar = (float).01;
/* rel. variation [-sval,sval] */
    s_copy(inputcom_1.scan, " ", (ftnlen)30, (ftnlen)1);

/*     output control parameters: */

/* scan defined by name */
    inputcom_1.iphsty = 1;
/* history - steps in z */
    inputcom_1.ipradi = 0;
/* radiation field - steps in z */
    inputcom_1.ippart = 0;
/* particles - steps in z */
    inputcom_1.ishsty = 1;
/* history - steps in t */
    inputcom_1.isradi = 0;
/* radiation field - steps in t */
    inputcom_1.ispart = 0;
/* particles - steps in t */
    for (i__ = 1; i__ <= 11; ++i__) {
	inputcom_1.lout[i__ - 1] = 1;
/* flags for output */
    }
    inputcom_1.lout[5] = 0;
/* disable diffraction calculation */
    for (i__ = 12; i__ <= 21; ++i__) {
	inputcom_1.lout[i__ - 1] = 0;
    }
    inputcom_1.iotail = 0;
/* <>0 => output include also slippage */
    inputcom_1.idump = 0;
/* <>0 => dump complete radiation field. */
    inputcom_1.nharm = 1;
/* # or harmonics in the bunching factor */
    inputcom_1.iallharm = 0;
/* <>0 => all higher harmonics are calculate */
    inputcom_1.iharmsc = 0;
/* <>0 => include effects on electron motion */
    inputcom_1.idmpfld = 0;
/* <>0 => dump radiation field */
    inputcom_1.idmppar = 0;
/* <>0 => dump particle distribution, with > */
    inputcom_1.ilog = 0;
/* <>0 => terminal output written to file */
    inputcom_1.ffspec = 0;
/* <0 => on-axis intensity in far field */
/* =0 => on-axis intensity in near field */

/*     file names */

/* >0 => total power in near field */
    s_copy(inputcom_1.beamfile, " ", (ftnlen)30, (ftnlen)1);
/* beam description fiel */
    s_copy(inputcom_1.fieldfile, " ", (ftnlen)30, (ftnlen)1);
/* input radiation field */
    s_copy(inputcom_1.maginfile, " ", (ftnlen)30, (ftnlen)1);
/* input magnetic field */
    s_copy(inputcom_1.magoutfile, " ", (ftnlen)30, (ftnlen)1);
/* output magnetic field */
    s_copy(inputcom_1.outputfile, "template.out", (ftnlen)30, (ftnlen)12);
/* output file */
    s_copy(inputcom_1.partfile, " ", (ftnlen)30, (ftnlen)1);
/* input particle file */
    s_copy(inputcom_1.distfile, " ", (ftnlen)30, (ftnlen)1);
/* input distribution file */
    s_copy(inputcom_1.radfile, " ", (ftnlen)30, (ftnlen)1);
/* input radiation file */
    s_copy(inputcom_1.filetype, "ORIGINAL", (ftnlen)30, (ftnlen)8);

/*     focussing: */

/* filetype for output files (sdds,xml) */
    inputcom_1.quadf = 1.23;
/* quad focus strength */
    inputcom_1.quadd = 0.;
/* quad defocus strength */
    inputcom_1.qfdx = (float)0.;
/* quadrupole offset in x */
    inputcom_1.qfdy = (float)0.;
/* quadrupole offset in y */
    inputcom_1.fl = (float)98.;
/* focus section length */
    inputcom_1.dl = (float)0.;
/* defocus section length */
    inputcom_1.drl = (float)0.;
/* drift length */
    inputcom_1.f1st = (float)0.;
/* start fodo at this point */
    inputcom_1.solen = (float)0.;
/* strength of solenoid field */
    inputcom_1.sl = (float)0.;

/*     time dependency parameters: */

/* solenoid length */
    inputcom_1.curlen = .001;
/* rms bunch length in m (<0 -> uniform beam */
    inputcom_1.ntail = -253;
/* starting slice in beam relative to curren */
    inputcom_1.nslice = 408;
/* #slices */
    inputcom_1.itdp = 0;
/* <>0 => time dependent code */
    inputcom_1.shotnoise = 1.;
/* scaling of the shotnoise meanvalue */
    inputcom_1.iall = 0;
/* <>0 => load phasespace with same loading */
    inputcom_1.ipseed = -1;
/* seeding of the random numbers for shotnoi */
    inputcom_1.ndcut = -1;
/* =<0 self optimized binning of ext. dist. */
    inputcom_1.isntyp = 0;

/*     extensions */

/* =0 -> fawley algorithm <>0 -> Penman algo */
    inputcom_1.isravg = 0;
/* <>0 enables energy loss by incorerent rad */
    inputcom_1.isrsig = 0;
/* <>0 enables growth of energy spread by in */
    inputcom_1.eloss = (float)0.;
/* energy loss per meter */
    inputcom_1.convharm = 1;
/* the harmonic to convert */
    inputcom_1.multconv = 0;
/* <>0 imported + converted slice is used mu */
    inputcom_1.ibfield = (float)0.;
/* field strength of magnetic chicane */
    inputcom_1.imagl = (float)0.;
/* length of bending magnet of chicane */
    inputcom_1.idril = (float)0.;
/* length of drift between bending magnets */
    inputcom_1.igamref = (float)0.;
/* length of reference energy of compressor */
    inputcom_1.alignradf = 0;
/* <>0 imported radfile is aligned to electron beam */
    inputcom_1.offsetradf = 0;

/*    transfermatrix */

/* if aligned, number of slices to skip */
    inputcom_1.trama = 0;
/* <>0 allows to transform 6D phase space */
    inputcom_1.itram11 = 1.;
/* transfer matrix set equal to identity ma */
    inputcom_1.itram12 = 0.;
    inputcom_1.itram13 = 0.;
    inputcom_1.itram14 = 0.;
    inputcom_1.itram15 = 0.;
    inputcom_1.itram16 = 0.;
    inputcom_1.itram21 = 0.;
    inputcom_1.itram22 = 1.;
    inputcom_1.itram23 = 0.;
    inputcom_1.itram24 = 0.;
    inputcom_1.itram25 = 0.;
    inputcom_1.itram26 = 0.;
    inputcom_1.itram31 = 0.;
    inputcom_1.itram32 = 0.;
    inputcom_1.itram33 = 1.;
    inputcom_1.itram34 = 0.;
    inputcom_1.itram35 = 0.;
    inputcom_1.itram36 = 0.;
    inputcom_1.itram41 = 0.;
    inputcom_1.itram42 = 0.;
    inputcom_1.itram43 = 0.;
    inputcom_1.itram44 = 1.;
    inputcom_1.itram45 = 0.;
    inputcom_1.itram46 = 0.;
    inputcom_1.itram51 = 0.;
    inputcom_1.itram52 = 0.;
    inputcom_1.itram53 = 0.;
    inputcom_1.itram54 = 0.;
    inputcom_1.itram55 = 1.;
    inputcom_1.itram56 = 0.;
    inputcom_1.itram61 = 0.;
    inputcom_1.itram62 = 0.;
    inputcom_1.itram63 = 0.;
    inputcom_1.itram64 = 0.;
    inputcom_1.itram65 = 0.;
    inputcom_1.itram66 = 1.;

    return 0;
} /* preset_ */





/* preset */
/* Subroutine */ int readradfile_(file, file_len)
char *file;
ftnlen file_len;
{
    /* Format strings */
    static char fmt_100[] = "(a)";
    static char fmt_110[] = "(\002Auto-adjustment of time window:\002,/,\002\
nslice=\002,i6,/\002ntail =\002,i6)";

    /* System generated locals */
    integer i__1;
    cllist cl__1;

    /* Builtin functions */
    integer i_indx(), s_rsfe(), do_fio(), e_rsfe(), s_wsli(), do_lio(), 
	    e_wsli(), i_dnnt(), f_clos(), s_wsfe(), e_wsfe();

    /* Local variables */
    static integer ncol, ipar[5], itmp;
    static char cerr__[50], line[511];
    static doublereal tmin, tmax;
    extern /* Subroutine */ int last_(), getfirstchar_();
    extern integer printerr_(), opentextfile_();
    static integer i__, j, idata, ft, ix;
    extern integer extractval_(), detectfiletype_();
    static integer nin;
    static doublereal ver, values[5], reverse, zoffset;
    extern /* Subroutine */ int getradfileinfo_(), touppercase_();

    /* Fortran I/O blocks */
    static cilist io___131 = { 1, 0, 1, fmt_100, 0 };
    static icilist io___136 = { 0, cerr__, 0, 0, 50, 1 };
    static icilist io___137 = { 0, cerr__, 0, 0, 50, 1 };
    static cilist io___141 = { 0, 0, 0, fmt_110, 0 };


/*     ============================================================= */
/*     read the file for external description of the radiation beam */
/*     ------------------------------------------------------------- */





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





    itmp = 0;
    tbunchcom_1.nraddata = -1;
    if (i_indx(file, " ", file_len, (ftnlen)1) == 1) {
	return 0;
    }

    if (inputcom_1.iscan > 0 && inputcom_1.iscan < 25) {
	i__ = printerr_(&c_n6, file, file_len);
	return 0;
    }

/*     read file */

    ft = detectfiletype_(file, file_len);

/* check for filetype */
    nin = opentextfile_(file, "old", &c__8, file_len, (ftnlen)3);
    if (nin < 0) {
	last_();
    }

/* stop program on error */
    tbunchcom_1.nraddata = -1;
/* # of rows not defined */
    idata = 0;
/* # of rows read */
    ncol = 5;
/* # of elements per line */
    reverse = (float)1.;
/* tail for ZPOS and head for TPOS comes first i */
    zoffset = (float)0.;
/* offset of shifting the radiation profile */
    ver = (float).1;
    i__1 = ncol;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ipar[i__ - 1] = i__;
/* basic order of input */
    }

L1:
    io___131.ciunit = nin;
    i__1 = s_rsfe(&io___131);
    if (i__1 != 0) {
	goto L100004;
    }
    i__1 = do_fio(&c__1, line, (ftnlen)511);
    if (i__1 != 0) {
	goto L100004;
    }
    i__1 = e_rsfe();
L100004:
    if (i__1 < 0) {
	goto L50;
    }
    if (i__1 > 0) {
	goto L20;
    }

/*     processing line */

    touppercase_(line, (ftnlen)511);
    getfirstchar_(line, &ix, (ftnlen)511);

    if (ix == 0) {
	goto L1;
    }
/* empty line */
    if (*(unsigned char *)&line[ix - 1] == '#') {
	goto L1;
    }
/* no comment used */
    if (*(unsigned char *)&line[ix - 1] == '?') {
	getradfileinfo_(line, ipar, &ncol, &tbunchcom_1.nraddata, &ver, &
		reverse, &zoffset, (ftnlen)511);
	goto L1;
    }

    if (tbunchcom_1.nraddata < 0 && ver < (float)1.) {
/* old version */
	i__ = extractval_(line, values, &c__1, (ftnlen)511);
	if (i__ < 0) {
	    s_wsli(&io___136);
	    do_lio(&c__9, &c__1, "Line number of RADFILE cannot be determined"
		    , (ftnlen)43);
	    e_wsli();
	    i__ = printerr_(&c_n8, cerr__, (ftnlen)50);
	}
	tbunchcom_1.nraddata = i_dnnt(values);
	goto L1;
    }

    i__ = extractval_(line, values, &ncol, (ftnlen)511);
    if (i__ < 0) {
	s_wsli(&io___137);
	do_lio(&c__9, &c__1, "RADFILE data line ", (ftnlen)18);
	i__1 = idata + 1;
	do_lio(&c__3, &c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " has bad format", (ftnlen)15);
	e_wsli();
	i__ = printerr_(&c_n8, cerr__, (ftnlen)50);
	last_();
    }

    ++idata;

/*     set default values */

    tbunchcom_1.tradpos[idata - 1] = inputcom_1.xlamds * inputcom_1.zsep * 
	    idata;
/* can we make this the default? */
    tbunchcom_1.tzrayl[idata - 1] = inputcom_1.zrayl;
    tbunchcom_1.tzwaist[idata - 1] = inputcom_1.zwaist;
    tbunchcom_1.tprad0[idata - 1] = inputcom_1.prad0;
    tbunchcom_1.tradphase[idata - 1] = (float)0.;

/*     write over with input data */

/* default no detuning or chirp. */
    i__1 = ncol;
    for (j = 1; j <= i__1; ++j) {
	if (ipar[j - 1] == 1) {
	    tbunchcom_1.tradpos[idata - 1] = reverse * values[j - 1];
	}
	if (ipar[j - 1] == -1) {
	    tbunchcom_1.tradpos[idata - 1] = -reverse * values[j - 1] * (
		    float)3e8;
	}
	if (ipar[j - 1] == 2) {
	    tbunchcom_1.tprad0[idata - 1] = values[j - 1];
	}
	if (ipar[j - 1] == 3) {
	    tbunchcom_1.tzrayl[idata - 1] = values[j - 1];
	}
	if (ipar[j - 1] == 4) {
	    tbunchcom_1.tzwaist[idata - 1] = values[j - 1];
	}
	if (ipar[j - 1] == 5) {
	    tbunchcom_1.tradphase[idata - 1] = values[j - 1];
	}
    }

/*     check for unphysical parameters */

    if (tbunchcom_1.tprad0[idata - 1] < 0.) {
	itmp = printerr_(&c_n8, "Radiation power negative in RADFILE", (
		ftnlen)35);
    }

    if (tbunchcom_1.tzrayl[idata - 1] < 0.) {
	itmp = printerr_(&c_n8, "ZRAYL negative in RADFILE", (ftnlen)25);
/* abort */
    }

    if (idata == 1) {
	tmin = tbunchcom_1.tradpos[idata - 1];
	tmax = tbunchcom_1.tradpos[idata - 1];
    } else {
	if (tbunchcom_1.tradpos[idata - 1] > tmax) {
	    tmax = tbunchcom_1.tradpos[idata - 1];
	}
	if (tbunchcom_1.tradpos[idata - 1] < tmin) {
	    tmin = tbunchcom_1.tradpos[idata - 1];
	}
    }

    goto L1;

L50:
    cl__1.cerr = 0;
    cl__1.cunit = nin;
    cl__1.csta = 0;
    f_clos(&cl__1);

    if (tbunchcom_1.nraddata >= 0 && idata != tbunchcom_1.nraddata) {
	i__ = printerr_(&c_n9, "RADFILE has fewer lines than defined", (
		ftnlen)36);
    }
    tbunchcom_1.nraddata = idata;
    if (idata < 2) {
	i__ = printerr_(&c_n8, "RADFILE contains less than 2 valid lines", (
		ftnlen)40);
	last_();
    }
    if (itmp != 0) {
	last_();
    }

    if (ver >= (float)1.) {
	i__1 = tbunchcom_1.nraddata;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    tbunchcom_1.tradpos[i__ - 1] -= tmin;
/* set time window to zero */
	}
    }

    i__1 = tbunchcom_1.nraddata;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tbunchcom_1.tradpos[i__ - 1] += zoffset;
/* set time window to zero */
    }

    if (inputcom_1.nslice <= 0) {
	inputcom_1.nslice = (integer) ((tmax - tmin) / inputcom_1.xlamds / 
		inputcom_1.zsep);
	if (ver >= (float)1.) {
	    inputcom_1.ntail = (integer) (zoffset / inputcom_1.xlamds / 
		    inputcom_1.zsep);
	} else {
	    inputcom_1.ntail = (integer) ((tmin + zoffset) / 
		    inputcom_1.xlamds / inputcom_1.zsep);
	}
	io___141.ciunit = inputcom_1.ilog;
	s_wsfe(&io___141);
	do_fio(&c__1, (char *)&inputcom_1.nslice, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&inputcom_1.ntail, (ftnlen)sizeof(integer));
	e_wsfe();
    }

    return 0;

L20:
    i__ = printerr_(&c_n2, file, file_len);
    cl__1.cerr = 0;
    cl__1.cunit = nin;
    cl__1.csta = 0;
    f_clos(&cl__1);
    last_();
    return 0;

/*     format statement */


} /* readradfile_ */


/* readradfile */
/* Subroutine */ int getradfileinfo_(line, ipar, ncol, ndata, ver, reverse, 
	zoffset, line_len)
char *line;
integer *ipar, *ncol, *ndata;
doublereal *ver, *reverse, *zoffset;
ftnlen line_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_indx(), i_len(), i_dnnt();
    /* Subroutine */ int s_copy();

    /* Local variables */
    static integer iarg, ierr;
    extern /* Subroutine */ int getfirstchar_();
    extern integer printerr_();
    static integer i__, j, n;
    static char cline[511];
    extern integer extractnumber_();
    static integer ix1, ix2;
    static doublereal val;
    static integer haszpos;

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


/*     version number */

    /* Parameter adjustments */
    --ipar;

    /* Function Body */
    i__ = i_indx(line, "VERSION", line_len, (ftnlen)7);
/* check for version number */
    n = i_len(line, line_len);
    if (i__ > 0) {
	i__1 = i__ + 6;
	ierr = extractnumber_(line + i__1, &val, n - i__1);
	if (ierr < 0) {
	    i__ = printerr_(&c_n9, "Unrecognized information line in RADFILE",
		     (ftnlen)40);
	} else {
	    *ver = val;
	}
	return 0;
    }

/*     reverse order */

    i__ = i_indx(line, "REVERSE", line_len, (ftnlen)7);
    n = i_len(line, line_len);
    if (i__ > 0) {
	*reverse = (float)-1.;
	return 0;
    }

/*     line numbers */

    i__ = i_indx(line, "SIZE", line_len, (ftnlen)4);
/* check for size argument (aka ndata) */
    if (i__ > 0) {
	i__1 = i__ + 3;
	ierr = extractnumber_(line + i__1, &val, n - i__1);
	if (ierr < 0) {
	    i__ = printerr_(&c_n9, "Unrecognized information line in RADFILE",
		     (ftnlen)40);
	} else {
	    *ndata = i_dnnt(&val);
	}
	return 0;
    }

/*     longitudinal offset */

    i__ = i_indx(line, "OFFSET", line_len, (ftnlen)6);
/* check for size argument (aka ndata) */
    if (i__ > 0) {
	i__1 = i__ + 5;
	ierr = extractnumber_(line + i__1, &val, n - i__1);
	if (ierr < 0) {
	    i__ = printerr_(&c_n9, "Unrecognized information line in RADFILE",
		     (ftnlen)40);
	} else {
	    *zoffset = val;
	}
	return 0;
    }

/*     column order */

    i__ = i_indx(line, "COLUMNS", line_len, (ftnlen)7);
/* check for colums headers */
    if (i__ > 0) {
	for (j = 1; j <= 5; ++j) {
	    ipar[j] = 0;
	}
	*ncol = 0;
	i__1 = i__ + 6;
	s_copy(cline, line + i__1, (ftnlen)511, n - i__1);
	haszpos = 0;
L1:
	getfirstchar_(cline, &ix1, (ftnlen)511);
	if (ix1 > 0) {
	    ix2 = 255;
/* search backwards */
	    i__1 = ix1 + 1;
	    for (j = 255; j >= i__1; --j) {
/* for first space after ix1 */
		if (*(unsigned char *)&cline[j - 1] == ' ') {
		    ix2 = j;
		}
	    }
	    iarg = 0;
	    if (i_indx(cline + (ix1 - 1), "ZPOS", ix2 - (ix1 - 1), (ftnlen)4) 
		    != 0) {
		iarg = 1;
	    }
	    if (i_indx(cline + (ix1 - 1), "TPOS", ix2 - (ix1 - 1), (ftnlen)4) 
		    != 0) {
		iarg = -1;
	    }
	    if (i_indx(cline + (ix1 - 1), "PRAD0", ix2 - (ix1 - 1), (ftnlen)5)
		     != 0) {
		iarg = 2;
	    }
	    if (i_indx(cline + (ix1 - 1), "ZRAYL", ix2 - (ix1 - 1), (ftnlen)5)
		     != 0) {
		iarg = 3;
	    }
	    if (i_indx(cline + (ix1 - 1), "ZWAIST", ix2 - (ix1 - 1), (ftnlen)
		    6) != 0) {
		iarg = 4;
	    }
	    if (i_indx(cline + (ix1 - 1), "PHASE", ix2 - (ix1 - 1), (ftnlen)5)
		     != 0) {
		iarg = 5;
	    }

	    if (iarg == 0) {
		for (j = 1; j <= 5; ++j) {
		    ipar[j] = j;
		}
		*ncol = 5;
		j = printerr_(&c_n9, "Unrecognized information line in RADFI\
LE", (ftnlen)40);
		return 0;
	    } else {
		++(*ncol);
		ipar[*ncol] = iarg;
		if (abs(iarg) == 1) {
		    haszpos = 1;
		}
		i__1 = ix2;
		s_copy(cline, cline + i__1, (ftnlen)511, 255 - i__1);
		if (*ncol < 5) {
		    goto L1;
		}
	    }
	}
	if (haszpos < 1) {
	    for (j = 1; j <= 5; ++j) {
		ipar[j] = j;
	    }
	    *ncol = 5;
	    j = printerr_(&c_n9, "ZPOS/TPOS column not specified in RADFILE", 
		    (ftnlen)41);
	    return 0;
	}
	return 0;
    }
    i__ = printerr_(&c_n9, "Unrecognized information line in RADFILE", (
	    ftnlen)40);
    return 0;
} /* getradfileinfo_ */

