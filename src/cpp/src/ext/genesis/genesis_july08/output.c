/* output.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

Extern struct {
    doublereal distversion, distrev;
    integer iout[39], nout, nfld, npar, ndump, firstout, nfin, irecpar, 
	    irecfld, kout, ndmp2, npin, nlog, ndis, ncoldis, iconv2t, iconv2g,
	     icolpar[10], ndistsize, iconv2px, iconv2py, nprobe, ftype, 
	    ftdist, ftpart, ftfield, ndumph[6], nfldh[6];
} iocom_;

#define iocom_1 iocom_

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
    integer mpi_id__, mpi_err__, mpi_size__, mpi_loop__, nfldmpi, nparmpi, nfldhmpi[6];
} mpicom_;

#define mpicom_1 mpicom_

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
    //doublecomplex crfield[476847], crsource[68121], crmatc[1827], cstep[7], crhm[68121], cwet[1827], cbet[1827];
    doublecomplex *crfield, *crsource, *crmatc, cstep[7], *crhm, *cwet, *cbet;
    doublereal dxy, xks, radphase, besselcoupling[7];
    integer nhloop, hloop[7];
} cartcom_;

#define cartcom_1 cartcom_

Extern struct {
    doublereal xkw0, xkper0, sval, svalout, gamma0_in__;
    integer npart0, inorun;
} simcom_;

#define simcom_1 simcom_

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
    //doublecomplex crtime[79803500];
    doublecomplex *crtime;
    integer ntmp;
} tslipcom_;

#define tslipcom_1 tslipcom_

/* Table of constant values */

static integer c_n21 = -21;
static integer c__1 = 1;
static integer c_n18 = -18;
static integer c__30 = 30;
static integer c__0 = 0;
static doublereal c_b22 = 1.;
static integer c_n1 = -1;
static integer c_n20 = -20;
static integer c__2 = 2;
static integer c__3 = 3;
static doublereal c_b399 = 10.;
static integer c__4 = 4;

/* Subroutine */ int output_(istepz, islice, xkw0)
integer *istepz, *islice;
doublereal *xkw0;
{
    extern /* Subroutine */ int outfield_(), diagno_(), status_(), outpart_();

/*     ============================================= */
/*     calls all output function */
/*     --------------------------------------------- */





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





/*     ------------------------------------------------------------------ */
/*     diagnostic */




    diagno_(istepz);
    status_(istepz, islice);
    if (*islice <= iocom_1.firstout) {
	return 0;
    }
    outfield_(istepz, islice, xkw0);
    outpart_(istepz, islice);
    return 0;
} /* output_ */

integer openoutputfile_(ierror, templatefilename, templatefilename_len)
integer *ierror;
char *templatefilename;
ftnlen templatefilename_len;
{
    /* Format strings */
    static char fmt_100[] = "(\002Please enter output file name \002)";
    static char fmt_105[] = "(a30)";
    static char fmt_110[] = "(\002- - - - - - - - - - - - - - - - - - - - - \
- - ------------------------------------------------\002,/,\002Genesis 1.3 o\
utput start\002,/,\002(Version \002,f3.1,\002 \002,a,\002)\002,/,/,\002Input\
 file name:  \002,a30,/)";
    static char fmt_1000[] = "(1x,\002$newrun\002,/,1x,\002aw0   =\002,1pd14\
.6,/,1x,\002xkx   =\002,1pd14.6,/,1x,\002xky   =\002,1pd14.6,/,1x,\002wcoefz=\
\002,3(1pd14.6,1x),/,1x,\002xlamd =\002,1pd14.6,/,1x,\002fbess0=\002,1pd14.6\
,/,1x,\002delaw =\002,1pd14.6,/,1x,\002iertyp=\002,i5,/,1x,\002iwityp=\002,i\
5,/,1x,\002awd   =\002,1pd14.6,/,1x,\002awx   =\002,1pd14.6,/,1x,\002awy   \
=\002,1pd14.6,/,1x,\002iseed =\002,i5,/,1x,\002npart =\002,i7,/,1x,\002gamma\
0=\002,1pd14.6,/,1x,\002delgam=\002,1pd14.6,/,1x,\002rxbeam=\002,1pd14.6,/,1\
x,\002rybeam=\002,1pd14.6,/,1x,\002alphax=\002,1pd14.6,/,1x,\002alphay=\002,\
1pd14.6,/,1x,\002emitx =\002,1pd14.6,/,1x,\002emity =\002,1pd14.6,/,1x,\002x\
beam =\002,1pd14.6,/,1x,\002ybeam =\002,1pd14.6,/,1x,\002pxbeam=\002,1pd14.6\
,/,1x,\002pybeam=\002,1pd14.6,/,1x,\002conditx =\002,1pd14.6,/,1x,\002condit\
y =\002,1pd14.6,/,1x,\002bunch =\002,1pd14.6,/,1x,\002bunchphase =\002,1pd14\
.6,/,1x,\002emod =\002,1pd14.6,/,1x,\002emodphase =\002,1pd14.6,/,1x,\002xla\
mds=\002,1pd14.6,/,1x,\002prad0 =\002,1pd14.6,/,1x,\002pradh0=\002,1pd14.6,/\
,1x,\002zrayl =\002,1pd14.6,/,1x,\002zwaist=\002,1pd14.6)";
    static char fmt_1001[] = "(1x,\002ncar  =\002,i5,/,1x,\002lbc   =\002,i5\
,/,1x,\002rmax0 =\002,1pd14.6,/,1x,\002dgrid =\002,1pd14.6,/,1x,\002nscr  \
=\002,i5,/,1x,\002nscz  =\002,i5,/,1x,\002nptr  =\002,i5,/,1x,\002nwig  =\
\002,i5,/,1x,\002zsep  =\002,1pd14.6,/,1x,\002delz  =\002,1pd14.6,/,1x,\002n\
sec  =\002,i5,/,1x,\002iorb  =\002,i5,/,1x,\002zstop =\002,1pd14.6,/,1x,\002\
magin =\002,i5,/,1x,\002magout=\002,i5,/,1x,\002quadf =\002,1pd14.6,/,1x,\
\002quadd =\002,1pd14.6,/,1x,\002fl    =\002,1pd14.6,/,1x,\002dl    =\002,1p\
d14.6,/,1x,\002drl   =\002,1pd14.6,/,1x,\002f1st  =\002,1pd14.6,/,1x,\002qfd\
x  =\002,1pd14.6,/,1x,\002qfdy  =\002,1pd14.6,/,1x,\002solen =\002,1pd14.6,/\
,1x,\002sl    =\002,1pd14.6,/,1x,\002ildgam=\002,i5,/,1x,\002ildpsi=\002,i5,\
/,1x,\002ildx  =\002,i5,/,1x,\002ildy  =\002,i5,/,1x,\002ildpx =\002,i5,/,1x,\
\002ildpy =\002,i5,/,1x,\002itgaus=\002,i5,/,1x,\002nbins =\002,i5,/,1x,\002\
igamgaus =\002,i5,/,1x,\002inverfc =\002,i5,/,1x,\002lout  =\002,19(1x,i1))";
    static char fmt_1002[] = "(1x,\002iphsty=\002,i5,/,1x,\002ishsty=\002,i5\
,/,1x,\002ippart=\002,i5,/,1x,\002ispart=\002,i5,/,1x,\002ipradi=\002,i5,/,1\
x,\002isradi=\002,i5,/,1x,\002idump =\002,i5,/,1x,\002iotail=\002,i5,/,1x\
,\002nharm =\002,i5,/,1x,\002iallharm =\002,i5,/,1x,\002iharmsc =\002,i5,/,1\
x,\002curpeak=\002,1pd14.6,/,1x,\002curlen=\002,1pd14.6,/,1x,\002ntail =\002\
,i5,/,1x,\002nslice=\002,i5,/,1x,\002shotnoise=\002,1pd14.6,/,1x,\002isntyp\
=\002,i5,/,1x,\002iall  =\002,i5,/,1x,\002itdp  =\002,i5,/,1x,\002ipseed=\
\002,i5,/,1x,\002iscan =\002,i5,/,1x,\002nscan =\002,i5,/,1x,\002svar  =\002\
,1pd14.6,/,1x,\002isravg=\002,i5,/,1x,\002isrsig=\002,i5,/,1x,\002cuttail\
=\002,1pd14.6,/,1x,\002eloss =\002,1pd14.6,/,1x,\002version=\002,1pd14.6,/,1\
x,\002ndcut =\002,i5)";
    static char fmt_1003[] = "(1x,\002idmpfld=\002,i5,/,1x,\002idmppar=\002,\
i5,/,1x,\002ilog  =\002,i5,/,1x,\002ffspec=\002,i5,/,1x,\002convharm=\002,i5\
,/,1x,\002ibfield=\002,1pd14.6,/,1x,\002imagl=  \002,1pd14.6,/,1x,\002idril=\
  \002,1pd14.6,/,1x,\002alignradf=\002,i5,/,1x,\002offsetradf=\002,i5,/,1x\
,\002multconv=\002,i5,/,1x,\002igamref=\002,1pd14.6,/,1x,\002rmax0sc=\002,1p\
d14.6,/,1x,\002iscrkup=\002,i5)";
    static char fmt_1004[] = "(1x,\002trama=\002,i5,/,1x,\002itram11=\002,1p\
d14.6,/,1x,\002itram12=\002,1pd14.6,/,1x,\002itram13=\002,1pd14.6,/,1x,\002i\
tram14=\002,1pd14.6,/,1x,\002itram15=\002,1pd14.6,/,1x,\002itram16=\002,1pd1\
4.6,/,1x,\002itram21=\002,1pd14.6,/,1x,\002itram22=\002,1pd14.6,/,1x,\002itr\
am23=\002,1pd14.6,/,1x,\002itram24=\002,1pd14.6,/,1x,\002itram25=\002,1pd14.\
6,/,1x,\002itram26=\002,1pd14.6,/,1x,\002itram31=\002,1pd14.6,/,1x,\002itram\
32=\002,1pd14.6,/,1x,\002itram33=\002,1pd14.6,/,1x,\002itram34=\002,1pd14.6,\
/,1x,\002itram35=\002,1pd14.6,/,1x,\002itram36=\002,1pd14.6,/,1x,\002itram41=\
\002,1pd14.6,/,1x,\002itram42=\002,1pd14.6,/,1x,\002itram43=\002,1pd14.6,/,1\
x,\002itram44=\002,1pd14.6,/,1x,\002itram45=\002,1pd14.6,/,1x,\002itram46\
=\002,1pd14.6,/,1x,\002itram51=\002,1pd14.6,/,1x,\002itram52=\002,1pd14.6,/,\
1x,\002itram53=\002,1pd14.6,/,1x,\002itram54=\002,1pd14.6,/,1x,\002itram55\
=\002,1pd14.6,/,1x,\002itram56=\002,1pd14.6,/,1x,\002itram61=\002,1pd14.6,/,\
1x,\002itram62=\002,1pd14.6,/,1x,\002itram63=\002,1pd14.6,/,1x,\002itram64\
=\002,1pd14.6,/,1x,\002itram65=\002,1pd14.6,/,1x,\002itram66=\002,1pd14.6)";
    static char fmt_1010[] = "(1x,a,\002'\002,a,\002'\002)";
    static char fmt_1020[] = "(1x,a4)";

    /* System generated locals */
    integer ret_val, i__1;
    olist o__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy();
    integer i_indx(), s_wsfe(), e_wsfe(), s_rsfe(), do_fio(), e_rsfe(), 
	    f_open();

    /* Local variables */
    static integer istr;
    extern /* Subroutine */ int last_();
    extern integer printerr_();
    static integer i__, isdefined;
    extern /* Subroutine */ int mpi_bcast__();
    extern integer strlen_();

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, fmt_100, 0 };
    static cilist io___4 = { 0, 5, 0, fmt_105, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___6 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___7 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_1002, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_1003, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_1004, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_1020, 0 };


/*     ================================================================== */
/*     initial output of genesis. */
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



/*     ------------------------------------------------------------------ */
/*     input/output control */





    isdefined = 0;
    if (*ierror == -3) {
	s_copy(inputcom_1.outputfile, templatefilename, (ftnlen)30, 
		templatefilename_len);
/* input file not file -> generate t */
    } else {
	isdefined = 1;
	if (i_indx(inputcom_1.outputfile, " ", (ftnlen)30, (ftnlen)1) == 1) {
/* no outputfile defined */
	    isdefined = 0;
	    if (inputcom_1.ilog != 0) {
		i__ = printerr_(&c_n21, "OUTPUTFILE", (ftnlen)10);
/* interaction require */
		last_();
	    }
	    if (mpicom_1.mpi_id__ == 0) {
L1:
		s_wsfe(&io___3);
		e_wsfe();
		s_rsfe(&io___4);
		do_fio(&c__1, inputcom_1.outputfile, (ftnlen)30);
		e_rsfe();
/* output file name not defined */
		if (i_indx(inputcom_1.outputfile, " ", (ftnlen)30, (ftnlen)1) 
			== 1) {
		    i__ = printerr_(&c_n18, inputcom_1.outputfile, (ftnlen)30)
			    ;
		    goto L1;
		}
	    }
	    if (mpicom_1.mpi_size__ > 1) {
		mpi_bcast__(inputcom_1.outputfile, &c__30, &c__1, &c__0, &
			c__1, &mpicom_1.mpi_err__, (ftnlen)30);
/* send outputfile name to a */
	    }
	}
/* ask for output filename */
    }

    iocom_1.nout = 9;
    ret_val = 0;
    if (mpicom_1.mpi_id__ > 0) {
	return ret_val;
/* only one template file or header file has to be writ */
    }

    o__1.oerr = 1;
    o__1.ounit = iocom_1.nout;
    o__1.ofnmlen = 30;
    o__1.ofnm = inputcom_1.outputfile;
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L10;
    }

/* open input fil */
    if (*ierror == 0) {
	io___5.ciunit = iocom_1.nout;
	s_wsfe(&io___5);
	do_fio(&c__1, (char *)&c_b22, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, "Unix", (ftnlen)4);
	do_fio(&c__1, inputcom_1.inputfile, (ftnlen)30);
	e_wsfe();
    } else {
	s_copy(inputcom_1.outputfile, "", (ftnlen)30, (ftnlen)0);
/* template is written -> prevents outputfile='temp */
    }

/*     dump input parmeters */

    io___6.ciunit = iocom_1.nout;
    s_wsfe(&io___6);
    do_fio(&c__1, (char *)&inputcom_1.aw0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.xkx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.xky, (ftnlen)sizeof(doublereal));
    for (i__ = 1; i__ <= 3; ++i__) {
	do_fio(&c__1, (char *)&inputcom_1.wcoefz[i__ - 1], (ftnlen)sizeof(
		doublereal));
    }
    do_fio(&c__1, (char *)&inputcom_1.xlamd, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.fbess0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.delaw, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.iertyp, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.iwityp, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.awd, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.awx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.awy, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.iseed, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.npart, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.gamma0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.delgam, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.rxbeam, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.rybeam, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.alphax, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.alphay, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.emitx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.emity, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.xbeam, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.ybeam, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.pxbeam, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.pybeam, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.conditx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.condity, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.bunch, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.bunchphase, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.emod, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.emodphase, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.xlamds, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.prad0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.pradh0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.zrayl, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.zwaist, (ftnlen)sizeof(doublereal));
    e_wsfe();
    io___7.ciunit = iocom_1.nout;
    s_wsfe(&io___7);
    do_fio(&c__1, (char *)&inputcom_1.ncar, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.lbc, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.rmax0, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.dgrid, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.nscr, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.nscz, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.nptr, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.nwig, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.zsep, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.delz, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.nsec, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.iorb, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.zstop, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.magin, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.magout, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.quadf, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.quadd, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.fl, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.dl, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.drl, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.f1st, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.qfdx, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.qfdy, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.solen, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.sl, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.ildgam, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ildpsi, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ildx, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ildy, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ildpx, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ildpy, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.itgaus, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.nbins, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.igamgaus, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.inverfc, (ftnlen)sizeof(integer));
    for (i__ = 1; i__ <= 21; ++i__) {
	do_fio(&c__1, (char *)&inputcom_1.lout[i__ - 1], (ftnlen)sizeof(
		integer));
    }
    e_wsfe();
    io___8.ciunit = iocom_1.nout;
    s_wsfe(&io___8);
    do_fio(&c__1, (char *)&inputcom_1.iphsty, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ishsty, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ippart, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ispart, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ipradi, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.isradi, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.idump, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.iotail, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.nharm, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.iallharm, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.iharmsc, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.curpeak, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.curlen, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.ntail, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.nslice, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.shotnoise, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.isntyp, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.iall, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.itdp, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ipseed, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.iscan, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.nscan, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.svar, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.isravg, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.isrsig, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.cuttail, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.eloss, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.version, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.ndcut, (ftnlen)sizeof(integer));
    e_wsfe();
    io___9.ciunit = iocom_1.nout;
    s_wsfe(&io___9);
    do_fio(&c__1, (char *)&inputcom_1.idmpfld, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.idmppar, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ilog, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ffspec, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.convharm, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.ibfield, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.imagl, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.idril, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.alignradf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.offsetradf, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.multconv, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.igamref, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.rmax0sc, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.iscrkup, (ftnlen)sizeof(integer));
    e_wsfe();
    io___10.ciunit = iocom_1.nout;
    s_wsfe(&io___10);
    do_fio(&c__1, (char *)&inputcom_1.trama, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&inputcom_1.itram11, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram12, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram13, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram14, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram15, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram16, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram21, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram22, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram23, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram24, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram25, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram26, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram31, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram32, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram33, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram34, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram35, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram36, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram41, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram42, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram43, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram44, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram45, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram46, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram51, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram52, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram53, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram54, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram55, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram56, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram61, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram62, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram63, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram64, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram65, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&inputcom_1.itram66, (ftnlen)sizeof(doublereal));
    e_wsfe();

    istr = strlen_(inputcom_1.beamfile, (ftnlen)30);
    if (istr > 1) {
	io___12.ciunit = iocom_1.nout;
	s_wsfe(&io___12);
	do_fio(&c__1, "beamfile =", (ftnlen)10);
	do_fio(&c__1, inputcom_1.beamfile, istr);
	e_wsfe();
    }
    istr = strlen_(inputcom_1.fieldfile, (ftnlen)30);
    if (istr > 1) {
	io___13.ciunit = iocom_1.nout;
	s_wsfe(&io___13);
	do_fio(&c__1, "fieldfile =", (ftnlen)11);
	do_fio(&c__1, inputcom_1.fieldfile, istr);
	e_wsfe();
    }
    istr = strlen_(inputcom_1.scan, (ftnlen)30);
    if (istr > 1) {
	io___14.ciunit = iocom_1.nout;
	s_wsfe(&io___14);
	do_fio(&c__1, "scan =", (ftnlen)6);
	do_fio(&c__1, inputcom_1.scan, istr);
	e_wsfe();
    }
    istr = strlen_(inputcom_1.outputfile, (ftnlen)30);
    if (istr > 1) {
	io___15.ciunit = iocom_1.nout;
	s_wsfe(&io___15);
	do_fio(&c__1, "outputfile =", (ftnlen)12);
	do_fio(&c__1, inputcom_1.outputfile, istr);
	e_wsfe();
    }
    istr = strlen_(inputcom_1.maginfile, (ftnlen)30);
    if (istr > 1) {
	io___16.ciunit = iocom_1.nout;
	s_wsfe(&io___16);
	do_fio(&c__1, "maginfile =", (ftnlen)11);
	do_fio(&c__1, inputcom_1.maginfile, istr);
	e_wsfe();
    }
    istr = strlen_(inputcom_1.magoutfile, (ftnlen)30);
    if (istr > 1) {
	io___17.ciunit = iocom_1.nout;
	s_wsfe(&io___17);
	do_fio(&c__1, "magoutfile =", (ftnlen)12);
	do_fio(&c__1, inputcom_1.magoutfile, istr);
	e_wsfe();
    }
    istr = strlen_(inputcom_1.partfile, (ftnlen)30);
    if (istr > 1) {
	io___18.ciunit = iocom_1.nout;
	s_wsfe(&io___18);
	do_fio(&c__1, "partfile =", (ftnlen)10);
	do_fio(&c__1, inputcom_1.partfile, istr);
	e_wsfe();
    }
    istr = strlen_(inputcom_1.distfile, (ftnlen)30);
    if (istr > 1) {
	io___19.ciunit = iocom_1.nout;
	s_wsfe(&io___19);
	do_fio(&c__1, "distfile =", (ftnlen)10);
	do_fio(&c__1, inputcom_1.distfile, istr);
	e_wsfe();
    }
    istr = strlen_(inputcom_1.radfile, (ftnlen)30);
    if (istr > 1) {
	io___20.ciunit = iocom_1.nout;
	s_wsfe(&io___20);
	do_fio(&c__1, "radfile =", (ftnlen)9);
	do_fio(&c__1, inputcom_1.radfile, istr);
	e_wsfe();
    }
    istr = strlen_(inputcom_1.filetype, (ftnlen)30);
    if (istr > 1) {
	io___21.ciunit = iocom_1.nout;
	s_wsfe(&io___21);
	do_fio(&c__1, "filetype =", (ftnlen)10);
	do_fio(&c__1, inputcom_1.filetype, istr);
	e_wsfe();
    }
    io___22.ciunit = iocom_1.nout;
    s_wsfe(&io___22);
    do_fio(&c__1, "$end", (ftnlen)4);
    e_wsfe();

    return ret_val;

L10:
    ret_val = printerr_(&c_n1, inputcom_1.outputfile, (ftnlen)30);
    return ret_val;

/*     format statement */


} /* openoutputfile_ */



/* of outheader */
/* Subroutine */ int outglob_()
{
    /* Format strings */
    static char fmt_10[] = "(a)";
    static char fmt_13[] = "(30i2)";
    static char fmt_11[] = "(i5,a)";
    static char fmt_12[] = "(e14.4,a)";
    static char fmt_14[] = "(i7,a)";
    static char fmt_50[] = "(\002Size of \002,a,\002 file [Mbytes]:\002,i4)";
    static char fmt_20[] = "((3(a14)))";
    static char fmt_21[] = "((3(1pe14.4)))";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe(), s_wsfi(), e_wsfi();
    /* Subroutine */ int s_copy();
    integer f_clos();

    /* Local variables */
    extern integer printerr_();
    static integer itmp, i1, i2, i3;
    static char titel[14*3], cwarn[40];
    static integer iz;

    /* Fortran I/O blocks */
    static cilist io___24 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_13, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_12, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_12, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_12, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_14, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_11, 0 };
    static icilist io___46 = { 0, cwarn, 0, fmt_50, 40, 1 };
    static icilist io___47 = { 0, cwarn, 0, fmt_50, 40, 1 };
    static icilist io___48 = { 0, cwarn, 0, fmt_50, 40, 1 };
    static cilist io___50 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___51 = { 0, 0, 0, fmt_21, 0 };


/*     ================================================================== */
/*     output of global parameter (t-independend): */
/*     z, wiggler field */
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
/*     wiggler parameters */




/*     ------------------------------------------------------------------ */
/*     input/output control */




/*     -------------------------------------------------------------------- */
/*     cartesian mesh */



/*     simulation control and normalisation parameter */


/*     -------------------------------------------------------------------- */
/*     time dependency of bunch + slippage field */





/*     ------------------------------------------------------------------ */
/*     common block for mpi implementation for genesis */





    if (mpicom_1.mpi_id__ > 0) {
	return 0;
    }

/*     ------------------------------------------------------------------ */
/*     output of t independent values */

/*     parameter for idl about file size etc. */

/* global variables only written by head no */
    for (iz = 16; iz <= 21; ++iz) {
/* starting at 16 lout indicates output for */
	inputcom_1.lout[iz - 1] <<= 2;
	if (inputcom_1.iallharm == 0 && iz - 14 != inputcom_1.nharm) {
	    inputcom_1.lout[iz - 1] = 0;
	}
    }

    io___24.ciunit = iocom_1.nout;
    s_wsfe(&io___24);
    do_fio(&c__1, "flags for output parameter", (ftnlen)26);
    e_wsfe();
    io___25.ciunit = iocom_1.nout;
    s_wsfe(&io___25);
    for (iz = 1; iz <= 21; ++iz) {
	do_fio(&c__1, (char *)&inputcom_1.lout[iz - 1], (ftnlen)sizeof(
		integer));
    }
    e_wsfe();

    itmp = 0;
    if (inputcom_1.iphsty > 0) {
	itmp = wigcom_1.nstepz / inputcom_1.iphsty + 1;
    }
    io___27.ciunit = iocom_1.nout;
    s_wsfe(&io___27);
    do_fio(&c__1, (char *)&itmp, (ftnlen)sizeof(integer));
    do_fio(&c__1, " entries per record", (ftnlen)19);
    e_wsfe();
    i1 = itmp;

    itmp = inputcom_1.nslice;
    iocom_1.firstout = tbunchcom_1.nslp * (1 - inputcom_1.iotail) * 
	    inputcom_1.itdp;
/* =0 for scan, steady state and io */
    itmp = inputcom_1.nslice - iocom_1.firstout;
    itmp /= inputcom_1.ishsty;
    io___29.ciunit = iocom_1.nout;
    s_wsfe(&io___29);
    do_fio(&c__1, (char *)&itmp, (ftnlen)sizeof(integer));
    do_fio(&c__1, " history records", (ftnlen)16);
    e_wsfe();

    io___30.ciunit = iocom_1.nout;
    s_wsfe(&io___30);
    do_fio(&c__1, (char *)&inputcom_1.xlamds, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, " wavelength", (ftnlen)11);
    e_wsfe();
    io___31.ciunit = iocom_1.nout;
    s_wsfe(&io___31);
    d__1 = (doublereal) inputcom_1.ishsty * inputcom_1.zsep * 
	    inputcom_1.xlamds;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, " seperation of output slices", (ftnlen)28);
    e_wsfe();

    io___32.ciunit = iocom_1.nout;
    s_wsfe(&io___32);
    do_fio(&c__1, (char *)&inputcom_1.ncar, (ftnlen)sizeof(integer));
    do_fio(&c__1, " number of gridpoints", (ftnlen)21);
    e_wsfe();
    io___33.ciunit = iocom_1.nout;
    s_wsfe(&io___33);
    d__1 = cartcom_1.dxy / simcom_1.xkper0;
    do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, " meshsize", (ftnlen)9);
    e_wsfe();
    io___34.ciunit = iocom_1.nout;
    s_wsfe(&io___34);
    do_fio(&c__1, (char *)&inputcom_1.npart, (ftnlen)sizeof(integer));
    do_fio(&c__1, " number of particles", (ftnlen)20);
    e_wsfe();

    if (inputcom_1.ippart > 0) {
	io___35.ciunit = iocom_1.nout;
	s_wsfe(&io___35);
	i__1 = wigcom_1.nstepz / inputcom_1.ippart + 1;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, " particle: records in z", (ftnlen)23);
	e_wsfe();
	io___36.ciunit = iocom_1.nout;
	s_wsfe(&io___36);
	i__1 = itmp / inputcom_1.ispart;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, " particle: records in t", (ftnlen)23);
	e_wsfe();
    } else {
	io___37.ciunit = iocom_1.nout;
	s_wsfe(&io___37);
	do_fio(&c__1, (char *)&c__0, (ftnlen)sizeof(integer));
	do_fio(&c__1, " particle: records in z", (ftnlen)23);
	e_wsfe();
	io___38.ciunit = iocom_1.nout;
	s_wsfe(&io___38);
	do_fio(&c__1, (char *)&c__0, (ftnlen)sizeof(integer));
	do_fio(&c__1, " particle: records in t", (ftnlen)23);
	e_wsfe();
    }
    if (inputcom_1.ipradi > 0) {
	io___39.ciunit = iocom_1.nout;
	s_wsfe(&io___39);
	i__1 = wigcom_1.nstepz / inputcom_1.ipradi + 1;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, " field: records in z", (ftnlen)20);
	e_wsfe();
	io___40.ciunit = iocom_1.nout;
	s_wsfe(&io___40);
	i__1 = itmp / inputcom_1.isradi;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	do_fio(&c__1, " field: records in t", (ftnlen)20);
	e_wsfe();
    } else {
	io___41.ciunit = iocom_1.nout;
	s_wsfe(&io___41);
	do_fio(&c__1, (char *)&c__0, (ftnlen)sizeof(integer));
	do_fio(&c__1, " field: records in z", (ftnlen)20);
	e_wsfe();
	io___42.ciunit = iocom_1.nout;
	s_wsfe(&io___42);
	do_fio(&c__1, (char *)&c__0, (ftnlen)sizeof(integer));
	do_fio(&c__1, " field: records in t", (ftnlen)20);
	e_wsfe();
    }

/*     calculation of output filesizes. */

    if (inputcom_1.itdp != 0 || inputcom_1.iscan > 0) {
	i2 = 0;
	for (i3 = 1; i3 <= 21; ++i3) {
	    i2 += inputcom_1.lout[i3 - 1];
/* count number of output parameters */
	}
	if (inputcom_1.iscan > 0) {
	    itmp = inputcom_1.nscan;
	}
	i2 *= 14;
/* each entry has 14 character */
	i2 = i2 * i1 * itmp / 1024 / 1024;
/* mutiply by #records and lines per rec */
	if (i2 > 1) {
	    s_wsfi(&io___46);
	    do_fio(&c__1, "history", (ftnlen)7);
	    do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
	    e_wsfi();
	    i2 = printerr_(&c_n20, cwarn, (ftnlen)40);
	}
    }
    if (inputcom_1.ippart > 0) {
	i2 = (wigcom_1.nstepz / inputcom_1.ippart + 1) * itmp / 
		inputcom_1.ispart;
/* number of records and entries */
	i2 = i2 * 48 * inputcom_1.npart / 1024 / 1024;
/* 6 variables (real*8) per part */
	if (i2 > 1) {
	    s_wsfi(&io___47);
	    do_fio(&c__1, "particle", (ftnlen)8);
	    do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
	    e_wsfi();
	    i2 = printerr_(&c_n20, cwarn, (ftnlen)40);
	}
    }
    if (inputcom_1.ipradi > 0) {
	i2 = (wigcom_1.nstepz / inputcom_1.ipradi + 1) * itmp / 
		inputcom_1.isradi;
/* number of records and entries */
	i2 = (i2 << 4) * inputcom_1.ncar * inputcom_1.ncar / 1024 / 1024;
/* ncar**2 grid points (complex* */
	if (i2 > 1) {
	    s_wsfi(&io___48);
	    do_fio(&c__1, "radiation", (ftnlen)9);
	    do_fio(&c__1, (char *)&i2, (ftnlen)sizeof(integer));
	    e_wsfi();
	    i2 = printerr_(&c_n20, cwarn, (ftnlen)40);
	}
    }

/*     output of global parameter (undulator field) */

    if (inputcom_1.iphsty <= 0) {
	return 0;
    }
    s_copy(titel, "    z[m]", (ftnlen)14, (ftnlen)8);
    s_copy(titel + 14, "    aw", (ftnlen)14, (ftnlen)6);
    s_copy(titel + 28, "    qfld", (ftnlen)14, (ftnlen)8);
    io___50.ciunit = iocom_1.nout;
    s_wsfe(&io___50);
    for (iz = 1; iz <= 3; ++iz) {
	do_fio(&c__1, titel + (iz - 1) * 14, (ftnlen)14);
    }
    e_wsfe();
/* table heading */
    i__1 = wigcom_1.nstepz;
    i__2 = inputcom_1.iphsty;
    for (iz = 0; i__2 < 0 ? iz >= i__1 : iz <= i__1; iz += i__2) {
	io___51.ciunit = iocom_1.nout;
	s_wsfe(&io___51);
	d__1 = iz * inputcom_1.delz * inputcom_1.xlamd;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&wigcom_1.awz[iz], (ftnlen)sizeof(doublereal));
	d__2 = wigcom_1.qfld[iz] / (float)586.;
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (mpicom_1.mpi_size__ > 1) {
/* redirect output if mpi is running */
	cl__1.cerr = 0;
	cl__1.cunit = iocom_1.nout;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

    return 0;


} /* outglob_ */


/* Subroutine */ int outhist_(islice)
integer *islice;
{
    /* Format strings */
    static char fmt_2[] = "(\002.slice\002,i6.6)";
    static char fmt_30[] = "((50(1pe14.4)))";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2, i__3;
    char ch__1[42];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi();
    /* Subroutine */ int s_cat();
    integer f_open(), s_wsfe(), e_wsfe(), f_clos();

    /* Local variables */
    static doublereal vout[39];
    static integer n;
    extern /* Subroutine */ int outhistheader_();
    static integer ih, il;
    extern integer strlen_();
    static integer ill;
    static char file_id__[12];

    /* Fortran I/O blocks */
    static icilist io___53 = { 0, file_id__, 0, fmt_2, 12, 1 };
    static cilist io___59 = { 0, 0, 0, fmt_30, 0 };


/*     ================================================================== */
/*     output calculation results */
/*        - history record */
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
/*     diagnostic */



/*     -------------------------------------------------------------------- */
/*     cartesian mesh */



/*     ------------------------------------------------------------------ */
/*     input/output control */





/*     ------------------------------------------------------------------ */
/*     common block for mpi implementation for genesis */





    if (inputcom_1.iphsty <= 0) {
	return 0;
    }
/* no output at all */
    if (*islice % inputcom_1.ishsty != 0) {
	return 0;
    }
/* output every ishstyt */
    if (*islice <= iocom_1.firstout) {
	return 0;
    }

/* no output if in first slip */
    if (mpicom_1.mpi_size__ > 1) {
/* creata temporary files */
	s_wsfi(&io___53);
	do_fio(&c__1, (char *)&(*islice), (ftnlen)sizeof(integer));
	e_wsfi();
	o__1.oerr = 0;
	o__1.ounit = iocom_1.nout;
	o__1.ofnmlen = strlen_(inputcom_1.outputfile, (ftnlen)30) + 12;
/* Writing concatenation */
	i__1[0] = strlen_(inputcom_1.outputfile, (ftnlen)30), a__1[0] = 
		inputcom_1.outputfile;
	i__1[1] = 12, a__1[1] = file_id__;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)42);
	o__1.ofnm = ch__1;
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }

    outhistheader_(islice);


    i__2 = diagcom_1.ihist;
    for (ih = 1; ih <= i__2; ++ih) {
	n = 0;
	i__3 = iocom_1.kout;
	for (il = 1; il <= i__3; ++il) {
	    if (iocom_1.iout[il - 1] == 1) {
		vout[il - 1] = diagcom_1.pgainhist[ih * 7 - 7];
	    }
	    if (iocom_1.iout[il - 1] == 2) {
		vout[il - 1] = diagcom_1.logp[ih - 1];
	    }
	    if (iocom_1.iout[il - 1] == 3) {
		vout[il - 1] = diagcom_1.pmidhist[ih * 7 - 7];
	    }
	    if (iocom_1.iout[il - 1] == 4) {
		vout[il - 1] = diagcom_1.phimid[ih * 7 - 7];
	    }
	    if (iocom_1.iout[il - 1] == 5) {
		vout[il - 1] = diagcom_1.whalf[ih - 1];
	    }
	    if (iocom_1.iout[il - 1] == 6) {
		vout[il - 1] = diagcom_1.diver[ih - 1];
	    }
	    if (iocom_1.iout[il - 1] == 7) {
		vout[il - 1] = diagcom_1.gamhist[ih - 1];
	    }
	    if (iocom_1.iout[il - 1] == 8) {
		vout[il - 1] = diagcom_1.pmodhist[ih * 7 - 7];
	    }
	    if (iocom_1.iout[il - 1] == 9) {
		vout[il - 1] = diagcom_1.xrms[ih - 1];
	    }
	    if (iocom_1.iout[il - 1] == 10) {
		vout[il - 1] = diagcom_1.yrms[ih - 1];
	    }
	    if (iocom_1.iout[il - 1] == 11) {
		vout[il - 1] = diagcom_1.error[ih - 1];
	    }
	    if (iocom_1.iout[il - 1] == 12) {
		vout[il - 1] = diagcom_1.xpos[ih - 1];
	    }
	    if (iocom_1.iout[il - 1] == 13) {
		vout[il - 1] = diagcom_1.ypos[ih - 1];
	    }
	    if (iocom_1.iout[il - 1] == 14) {
		vout[il - 1] = diagcom_1.dgamhist[ih - 1];
	    }
	    if (iocom_1.iout[il - 1] == 15) {
		vout[il - 1] = diagcom_1.ffield[ih * 7 - 7];
	    }
/*     output of harmonic content 16 -> 2nd harmonc, 17-> 3rd harm etc. */
/*     one entry in lout indicates to print four output parameters: */

	    for (ill = 2; ill <= 7; ++ill) {
		if (iocom_1.iout[il - 1] == ill + 14) {
		    vout[il + n - 1] = diagcom_1.pmodhist[ill + ih * 7 - 8];
		    vout[il + n] = diagcom_1.pgainhist[ill + ih * 7 - 8];
		    vout[il + n + 1] = diagcom_1.bunphase[ill + ih * 7 - 8];
		    vout[il + n + 2] = diagcom_1.pmidhist[ill + ih * 7 - 8];
		    n += 3;
		}
	    }
	}
	io___59.ciunit = iocom_1.nout;
	s_wsfe(&io___59);
	i__3 = iocom_1.kout;
	for (il = 1; il <= i__3; ++il) {
	    do_fio(&c__1, (char *)&vout[il - 1], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
    }

/*     format statements */


    if (mpicom_1.mpi_size__ > 1) {
	cl__1.cerr = 0;
	cl__1.cunit = iocom_1.nout;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }

    return 0;
} /* outhist_ */


/* output */
/* Subroutine */ int outhistheader_(islice)
integer *islice;
{
    /* Format strings */
    static char fmt_50[] = "(i1.1)";
    static char fmt_5[] = "(\002***  writing history record for slice \002,i\
5)";
    static char fmt_10[] = "(/\002********** output: slice \002,i5/5x,\002  \
    =================\002/1x,1pe14.4,\002 current\002//)";
    static char fmt_11[] = "(/\002********** output: slice \002,i5/5x,\002  \
    =================\002/1x,1pe14.4,\002 scan value\002//)";
    static char fmt_20[] = "((50(a14)))";

    /* System generated locals */
    address a__1[3];
    integer i__1[3], i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy();
    integer s_wsfi(), do_fio(), e_wsfi();
    /* Subroutine */ int s_cat();
    integer s_wsfe(), e_wsfe();

    /* Local variables */
    static char carh[1];
    static integer m;
    static char titel[14*39];
    static integer ih, iz;

    /* Fortran I/O blocks */
    static icilist io___63 = { 0, carh, 0, fmt_50, 1, 1 };
    static cilist io___66 = { 0, 0, 0, fmt_5, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___68 = { 0, 0, 0, fmt_11, 0 };
    static cilist io___69 = { 0, 0, 0, fmt_20, 0 };


/*     ================================================================== */
/*     output calculation results */
/*        - history record */
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



/*     ------------------------------------------------------------------ */
/*     electron beam */




/*     ------------------------------------------------------------------ */
/*     input/output control */





/*     simulation control and normalisation parameter */



    if (inputcom_1.iphsty <= 0) {
	return 0;
    }
/* no output at all */
    if (*islice % inputcom_1.ishsty != 0) {
	return 0;
    }

/*     ----------------------------------------------------------------- */
/*     create output array for optional output */

/* output ishstyth slic */
    s_copy(titel, "    power", (ftnlen)14, (ftnlen)9);
    s_copy(titel + 14, "    increment", (ftnlen)14, (ftnlen)13);
    s_copy(titel + 28, "    p_mid", (ftnlen)14, (ftnlen)9);
    s_copy(titel + 42, "    phi_mid", (ftnlen)14, (ftnlen)11);
    s_copy(titel + 56, "    r_size", (ftnlen)14, (ftnlen)10);
    s_copy(titel + 70, "    angle", (ftnlen)14, (ftnlen)9);
    s_copy(titel + 84, "    energy", (ftnlen)14, (ftnlen)10);
    s_copy(titel + 98, "    bunching", (ftnlen)14, (ftnlen)12);
    s_copy(titel + 112, "    xrms", (ftnlen)14, (ftnlen)8);
    s_copy(titel + 126, "    yrms", (ftnlen)14, (ftnlen)8);
    s_copy(titel + 140, "    error", (ftnlen)14, (ftnlen)9);
    s_copy(titel + 154, "    <x>", (ftnlen)14, (ftnlen)7);
    s_copy(titel + 168, "    <y>", (ftnlen)14, (ftnlen)7);
    s_copy(titel + 182, "    e-spread", (ftnlen)14, (ftnlen)12);
    s_copy(titel + 196, "    far_field", (ftnlen)14, (ftnlen)13);
/*     KG ---------------------------------------------------- */
    s_copy(titel + 210, "    2nd_bunching", (ftnlen)14, (ftnlen)16);
    s_copy(titel + 224, "    2nd_power", (ftnlen)14, (ftnlen)13);
    s_copy(titel + 238, "    2nd_phase", (ftnlen)14, (ftnlen)13);
    s_copy(titel + 252, "    2nd_p-mid", (ftnlen)14, (ftnlen)13);
    s_copy(titel + 266, "    3rd_bunching", (ftnlen)14, (ftnlen)16);
    s_copy(titel + 280, "    3rd_power", (ftnlen)14, (ftnlen)13);
    s_copy(titel + 294, "    3rd_phase", (ftnlen)14, (ftnlen)13);
    s_copy(titel + 308, "    3rd_p-mid", (ftnlen)14, (ftnlen)13);
    for (ih = 4; ih <= 7; ++ih) {
	s_wsfi(&io___63);
	do_fio(&c__1, (char *)&ih, (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 3, a__1[0] = "   ";
	i__1[1] = 1, a__1[1] = carh;
	i__1[2] = 11, a__1[2] = "th_bunching";
	//s_cat(titel + ((ih - 4 << 2) + 23) * 14, a__1, i__1, &c__3, (ftnlen)14);
	s_cat(titel + (((ih - 4) << 2) + 23) * 14, a__1, i__1, &c__3, (ftnlen)14); //OC port
/* Writing concatenation */
	i__1[0] = 3, a__1[0] = "   ";
	i__1[1] = 1, a__1[1] = carh;
	i__1[2] = 8, a__1[2] = "th_power";
	//s_cat(titel + ((ih - 4 << 2) + 24) * 14, a__1, i__1, &c__3, (ftnlen)14);
	s_cat(titel + (((ih - 4) << 2) + 24) * 14, a__1, i__1, &c__3, (ftnlen)14); //OC port
/* Writing concatenation */
	i__1[0] = 3, a__1[0] = "   ";
	i__1[1] = 1, a__1[1] = carh;
	i__1[2] = 8, a__1[2] = "th_phase";
	//s_cat(titel + ((ih - 4 << 2) + 25) * 14, a__1, i__1, &c__3, (ftnlen)14);
	s_cat(titel + (((ih - 4) << 2) + 25) * 14, a__1, i__1, &c__3, (ftnlen)14); //OC port
/* Writing concatenation */
	i__1[0] = 3, a__1[0] = "   ";
	i__1[1] = 1, a__1[1] = carh;
	i__1[2] = 8, a__1[2] = "th_p-mid";
	//s_cat(titel + ((ih - 4 << 2) + 26) * 14, a__1, i__1, &c__3, (ftnlen)14);
	s_cat(titel + (((ih - 4) << 2) + 26) * 14, a__1, i__1, &c__3, (ftnlen)14); //OC port
    }

/*     -------------------------------------------------------------------- */
/*     select output */

    iocom_1.kout = 0;
    m = 0;
    for (iz = 1; iz <= 15; ++iz) {
	iocom_1.iout[iz - 1] = 0;
	if (inputcom_1.lout[iz - 1] != 0) {
	    ++iocom_1.kout;
	    iocom_1.iout[iocom_1.kout - 1] = iz;
	}
    }

    for (iz = 16; iz <= 21; ++iz) {
	if (inputcom_1.lout[iz - 1] != 0) {
	    iocom_1.iout[iocom_1.kout] = m + iz;
	    iocom_1.iout[iocom_1.kout + 1] = m + iz + 1;
	    iocom_1.iout[iocom_1.kout + 2] = m + iz + 2;
	    iocom_1.iout[iocom_1.kout + 3] = m + iz + 3;
	    iocom_1.kout += 4;
	}
	m += 3;
    }


/*     --------------------------------------------------------------------- */
/*     write header */

    io___66.ciunit = iocom_1.nlog;
    s_wsfe(&io___66);
    do_fio(&c__1, (char *)&(*islice), (ftnlen)sizeof(integer));
    e_wsfe();

    if (inputcom_1.iscan <= 0 || inputcom_1.iscan > 22) {
	io___67.ciunit = iocom_1.nout;
	s_wsfe(&io___67);
	do_fio(&c__1, (char *)&(*islice), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&beamcom_1.xcuren, (ftnlen)sizeof(doublereal));
	e_wsfe();
/* time dependence */
    } else {
	io___68.ciunit = iocom_1.nout;
	s_wsfe(&io___68);
	do_fio(&c__1, (char *)&(*islice), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&simcom_1.svalout, (ftnlen)sizeof(doublereal));
	e_wsfe();
/* scan parameter */
    }
    io___69.ciunit = iocom_1.nout;
    s_wsfe(&io___69);
    i__2 = iocom_1.kout;
    for (iz = 1; iz <= i__2; ++iz) {
	do_fio(&c__1, titel + (iocom_1.iout[iz - 1] - 1) * 14, (ftnlen)14);
    }
    e_wsfe();

/*     format statements */
/*     ------------------------------------------------------------------ */
/* table heading */

    return 0;
} /* outhistheader_ */




/* outputheader */
/* Subroutine */ int outpart_(istepz, islice)
integer *istepz, *islice;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wdue(), do_uio(), e_wdue();

    /* Local variables */
    extern /* Subroutine */ int rpos_();
    static integer iz;
    extern /* Subroutine */ int getpsi_();

    /* Fortran I/O blocks */
    static cilist io___71 = { 0, 0, 0, 0, 0 };
    static cilist io___72 = { 0, 0, 0, 0, 0 };
    static cilist io___73 = { 0, 0, 0, 0, 0 };
    static cilist io___74 = { 0, 0, 0, 0, 0 };
    static cilist io___75 = { 0, 0, 0, 0, 0 };
    static cilist io___76 = { 0, 0, 0, 0, 0 };


/*     ================================================================== */
/*     output of global parameter (t-independend): */
/*     z, wiggler field */
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


/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







/*     simulation control and normalisation parameter */



/*     ------------------------------------------------------------------ */
/*     output of t independent values with first slice */

    if (inputcom_1.ippart <= 0 || inputcom_1.ispart <= 0) {
	return 0;
    }
/* no output at all */
    if (*istepz % inputcom_1.ippart != 0) {
	return 0;
    }
/* output ippartth step */
    if (*islice % inputcom_1.ispart != 0) {
	return 0;
    }

/* output ispartth slic */
    if (*istepz == 0) {
	rpos_(&c__0, beamcom_1.xpart, beamcom_1.ypart);
    }
    getpsi_(workspace_1.p1);

    if (inputcom_1.npart < simcom_1.npart0) {
/* check for particle loss */
	i__1 = simcom_1.npart0;
	for (iz = inputcom_1.npart + 1; iz <= i__1; ++iz) {
/* indicate lost particles with neg. */
	    beamcom_1.gamma[iz - 1] = (float)-1.;
	}
    }

    io___71.ciunit = iocom_1.npar;
    io___71.cirec = iocom_1.irecpar;
    s_wdue(&io___71);
    i__1 = simcom_1.npart0;
    for (iz = 1; iz <= i__1; ++iz) {
	do_uio(&c__1, (char *)&beamcom_1.gamma[iz - 1], (ftnlen)sizeof(
		doublereal));
    }
    e_wdue();
    io___72.ciunit = iocom_1.npar;
    io___72.cirec = iocom_1.irecpar + 1;
    s_wdue(&io___72);
    i__1 = simcom_1.npart0;
    for (iz = 1; iz <= i__1; ++iz) {
	do_uio(&c__1, (char *)&workspace_1.p1[iz - 1], (ftnlen)sizeof(
		doublereal));
    }
    e_wdue();
/*      write(npar,rec=irecpar+1) (theta(iz),iz=1,npart0) */
    io___73.ciunit = iocom_1.npar;
    io___73.cirec = iocom_1.irecpar + 2;
    s_wdue(&io___73);
    i__1 = simcom_1.npart0;
    for (iz = 1; iz <= i__1; ++iz) {
	d__1 = beamcom_1.xpart[iz - 1] / simcom_1.xkper0;
	do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    }
    e_wdue();
    io___74.ciunit = iocom_1.npar;
    io___74.cirec = iocom_1.irecpar + 3;
    s_wdue(&io___74);
    i__1 = simcom_1.npart0;
    for (iz = 1; iz <= i__1; ++iz) {
	d__1 = beamcom_1.ypart[iz - 1] / simcom_1.xkper0;
	do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    }
    e_wdue();
    io___75.ciunit = iocom_1.npar;
    io___75.cirec = iocom_1.irecpar + 4;
    s_wdue(&io___75);
    i__1 = simcom_1.npart0;
    for (iz = 1; iz <= i__1; ++iz) {
	do_uio(&c__1, (char *)&beamcom_1.px[iz - 1], (ftnlen)sizeof(
		doublereal));
    }
    e_wdue();
    io___76.ciunit = iocom_1.npar;
    io___76.cirec = iocom_1.irecpar + 5;
    s_wdue(&io___76);
    i__1 = simcom_1.npart0;
    for (iz = 1; iz <= i__1; ++iz) {
	do_uio(&c__1, (char *)&beamcom_1.py[iz - 1], (ftnlen)sizeof(
		doublereal));
    }
    e_wdue();
    iocom_1.irecpar += 6;
    return 0;
} /* outpart_ */




/* Subroutine */ int status_(istepz, islice)
integer *istepz, *islice;
{
    /* Format strings */
    static char fmt_20[] = "(\002Slice \002,i5,\002: Simulation \002,i3,\002\
% completed.\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double d_mod();
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static doublereal xper, yper;

    /* Fortran I/O blocks */
    static cilist io___79 = { 0, 0, 0, fmt_20, 0 };


/*     ================================================================== */
/*     let user know % complete at every 10%. */
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
/*     wiggler parameters */





/*     ------------------------------------------------------------------ */
/*     input/output control */





    xper = (real) (*istepz) * 100. / (real) wigcom_1.nstepz;
    yper = (real) (*istepz - 1) * 100. / (real) wigcom_1.nstepz;
    if (d_mod(&xper, &c_b399) < d_mod(&yper, &c_b399)) {
	io___79.ciunit = iocom_1.nlog;
	s_wsfe(&io___79);
	do_fio(&c__1, (char *)&(*islice), (ftnlen)sizeof(integer));
	i__1 = (integer) xper;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsfe();
    }

    return 0;
} /* status_ */




/* status */
/* Subroutine */ int outfield_(istepz, islice, xkper0)
integer *istepz, *islice;
doublereal *xkper0;
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();
    integer s_wdue(), do_uio(), e_wdue();
    double d_imag();

    /* Local variables */
    static integer i__, ifile, ih;
    static doublereal scltmp;
    static integer ioffset;

    /* Fortran I/O blocks */
    static cilist io___81 = { 0, 0, 0, 0, 0 };
    static cilist io___83 = { 0, 0, 0, 0, 0 };
    static cilist io___87 = { 0, 0, 0, 0, 0 };
    static cilist io___88 = { 0, 0, 0, 0, 0 };


/*     ================================================================== */
/*     dump fieldarray */
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





/*     -------------------------------------------------------------------- */
/*     cartesian mesh */




    if (inputcom_1.ipradi <= 0 || inputcom_1.isradi <= 0) {
	return 0;
    }
/* no output at all */
    if (*istepz % inputcom_1.ipradi != 0) {
	return 0;
    }
/* output ipradith step */
    if (*islice % inputcom_1.isradi != 0) {
	return 0;
    }

/* output isradith slic */
    scltmp = cartcom_1.dxy * 510999.06 * *xkper0 / cartcom_1.xks / sqrt(
	    376.73);

    io___81.ciunit = iocom_1.nfld;
    io___81.cirec = iocom_1.irecfld;
    s_wdue(&io___81);
    i__1 = inputcom_1.ncar * inputcom_1.ncar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ - 1;
	d__1 = scltmp * cartcom_1.crfield[i__2].r;
	do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    }
    e_wdue();
    io___83.ciunit = iocom_1.nfld;
    io___83.cirec = iocom_1.irecfld + 1;
    s_wdue(&io___83);
    i__2 = inputcom_1.ncar * inputcom_1.ncar;
    for (i__ = 1; i__ <= i__2; ++i__) {
	d__1 = scltmp * d_imag(&cartcom_1.crfield[i__ - 1]);
	do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
    }
    e_wdue();
    i__2 = cartcom_1.nhloop;
    for (ih = 2; ih <= i__2; ++ih) {
	ioffset = (ih - 1) * inputcom_1.ncar * inputcom_1.ncar;
	ifile = iocom_1.nfldh[ih - 2];
	io___87.ciunit = ifile;
	io___87.cirec = iocom_1.irecfld;
	s_wdue(&io___87);
	i__1 = inputcom_1.ncar * inputcom_1.ncar + ioffset;
	for (i__ = ioffset + 1; i__ <= i__1; ++i__) {
	    i__3 = i__ - 1;
	    d__1 = scltmp / cartcom_1.hloop[ih - 1] * cartcom_1.crfield[i__3]
		    .r;
	    do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	}
	e_wdue();
	io___88.ciunit = ifile;
	io___88.cirec = iocom_1.irecfld + 1;
	s_wdue(&io___88);
	i__3 = inputcom_1.ncar * inputcom_1.ncar + ioffset;
	for (i__ = ioffset + 1; i__ <= i__3; ++i__) {
	    d__1 = scltmp / cartcom_1.hloop[ih - 1] * d_imag(&
		    cartcom_1.crfield[i__ - 1]);
	    do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	}
	e_wdue();
    }

    iocom_1.irecfld += 2;
    return 0;
} /* outfield_ */



/* Subroutine */ int outdump_(islice)
integer *islice;
{
    /* Format strings */
    static char fmt_2[] = "(\002.slice\002,i6.6)";
    static char fmt_3[] = "(i1.1)";

    /* System generated locals */
    address a__1[3], a__2[4];
    integer i__1, i__2[3], i__3, i__4[4], i__5, i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;
    char ch__1[46], ch__2[47];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi();
    /* Subroutine */ int s_cat();
    integer f_open(), s_wdue(), do_uio(), e_wdue(), f_clos();
    double sqrt();

    /* Local variables */
    static integer ih, ioffset;
    static doublereal scltmp;
    static integer ndmp2tmp;
    extern integer printerr_();
    static char file_id__[12], harm_id__[1];
    static integer i__, j, ifile;
    extern integer strlen_();

    /* Fortran I/O blocks */
    static icilist io___90 = { 0, file_id__, 0, fmt_2, 12, 1 };
    static cilist io___94 = { 0, 0, 0, 0, 0 };
    static cilist io___95 = { 0, 0, 0, 0, 0 };
    static cilist io___96 = { 0, 0, 0, 0, 0 };
    static cilist io___97 = { 0, 0, 0, 0, 0 };
    static cilist io___98 = { 0, 0, 0, 0, 0 };
    static cilist io___99 = { 0, 0, 0, 0, 0 };
    static cilist io___101 = { 0, 0, 0, 0, 0 };
    static icilist io___106 = { 0, harm_id__, 0, fmt_3, 1, 1 };
    static cilist io___107 = { 0, 0, 0, 0, 0 };


/*     ==================================================================== */
/*     dump complete field array for future use */
/*     -------------------------------------------------------------------- */





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



/*     file:   timerec.com */

/*     time dependent common block - huge!!!! */
/*     ------------------------------------------------------------------ */




/*     -------------------------------------------------------------------- */
/*     slippage field */

/* must be > than ncar*ncar*nslp(98) */
/* <>0 -> slippage is stored on disk */



/*     ------------------------------------------------------------------ */
/*     input/output control */




/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







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
/*     electron beam */




/*     simulation control and normalisation parameter */



/*     ------------------------------------------------------------------ */
/*     common block for mpi implementation for genesis */





    if (*islice <= iocom_1.firstout) {
	return 0;
    }

/* suppress for IOTAIL=0 */
    s_wsfi(&io___90);
    do_fio(&c__1, (char *)&(*islice), (ftnlen)sizeof(integer));
    e_wsfi();

/*     particle distribution */

    if (inputcom_1.idmppar != 0) {
	if (inputcom_1.npart < simcom_1.npart0) {
	    i__1 = simcom_1.npart0;
	    for (i__ = inputcom_1.npart + 1; i__ <= i__1; ++i__) {
/* check for particle losses */
		beamcom_1.gamma[i__ - 1] = (float)-1.;
/* indicate lost particles wit */
	    }
	}

/*       writing the record */

	j = (*islice - iocom_1.firstout - 1) * 6 + 1;

	if (mpicom_1.mpi_size__ > 1) {
/* creata temporary files */
	    ndmp2tmp = iocom_1.ndmp2;
	    iocom_1.ndmp2 += 30;
	    o__1.oerr = 1;
	    o__1.ounit = iocom_1.ndmp2;
	    o__1.ofnmlen = strlen_(inputcom_1.outputfile, (ftnlen)30) + 16;
/* Writing concatenation */
	    i__2[0] = strlen_(inputcom_1.outputfile, (ftnlen)30), a__1[0] = 
		    inputcom_1.outputfile;
	    i__2[1] = 4, a__1[1] = ".dpa";
	    i__2[2] = 12, a__1[2] = file_id__;
	    s_cat(ch__1, a__1, i__2, &c__3, (ftnlen)46);
	    o__1.ofnm = ch__1;
	    o__1.orl = inputcom_1.npart << 3;
	    o__1.osta = "unknown";
	    o__1.oacc = "direct";
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L100;
	    }
	    j = 1;
	}

	io___94.ciunit = iocom_1.ndmp2;
	io___94.cirec = j;
	s_wdue(&io___94);
	i__1 = simcom_1.npart0;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_uio(&c__1, (char *)&beamcom_1.gamma[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wdue();
	io___95.ciunit = iocom_1.ndmp2;
	io___95.cirec = j + 1;
	s_wdue(&io___95);
	i__1 = simcom_1.npart0;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_uio(&c__1, (char *)&beamcom_1.theta[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wdue();
	io___96.ciunit = iocom_1.ndmp2;
	io___96.cirec = j + 2;
	s_wdue(&io___96);
	i__1 = simcom_1.npart0;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__1 = beamcom_1.xpart[i__ - 1] / simcom_1.xkper0;
	    do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	}
	e_wdue();
	io___97.ciunit = iocom_1.ndmp2;
	io___97.cirec = j + 3;
	s_wdue(&io___97);
	i__1 = simcom_1.npart0;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    d__1 = beamcom_1.ypart[i__ - 1] / simcom_1.xkper0;
	    do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	}
	e_wdue();
	io___98.ciunit = iocom_1.ndmp2;
	io___98.cirec = j + 4;
	s_wdue(&io___98);
	i__1 = simcom_1.npart0;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_uio(&c__1, (char *)&beamcom_1.px[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wdue();
	io___99.ciunit = iocom_1.ndmp2;
	io___99.cirec = j + 5;
	s_wdue(&io___99);
	i__1 = simcom_1.npart0;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_uio(&c__1, (char *)&beamcom_1.py[i__ - 1], (ftnlen)sizeof(
		    doublereal));
	}
	e_wdue();

	if (mpicom_1.mpi_size__ > 1) {
	    cl__1.cerr = 0;
	    cl__1.cunit = iocom_1.ndmp2;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    iocom_1.ndmp2 = ndmp2tmp;
	}
    }

/*     field distribution */

    if (inputcom_1.idmpfld == 0) {
	return 0;
    }

/*     problems arise if the dump is used for another run */
/*     any change in undulator period. the field has twice a scaling */
/*     with xkper0 - 1. eikonal equation + 2. normalization of dxy */

    scltmp = cartcom_1.dxy * 510999.06 * simcom_1.xkper0 / cartcom_1.xks / 
	    sqrt(376.73);


    j = *islice - iocom_1.firstout;
    if (mpicom_1.mpi_size__ > 1) {
/* creata temporary files */
	ndmp2tmp = iocom_1.ndump;
	iocom_1.ndump += 30;
	o__1.oerr = 1;
	o__1.ounit = iocom_1.ndump;
	o__1.ofnmlen = strlen_(inputcom_1.outputfile, (ftnlen)30) + 16;
/* Writing concatenation */
	i__2[0] = strlen_(inputcom_1.outputfile, (ftnlen)30), a__1[0] = 
		inputcom_1.outputfile;
	i__2[1] = 4, a__1[1] = ".dfl";
	i__2[2] = 12, a__1[2] = file_id__;
	s_cat(ch__1, a__1, i__2, &c__3, (ftnlen)46);
	o__1.ofnm = ch__1;
	o__1.orl = (inputcom_1.ncar << 4) * inputcom_1.ncar;
	o__1.osta = "unknown";
	o__1.oacc = "direct";
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L100;
	}
	j = 1;
    }

    io___101.ciunit = iocom_1.ndump;
    io___101.cirec = j;
    s_wdue(&io___101);
    i__1 = inputcom_1.ncar * inputcom_1.ncar;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__3 = i__ - 1;
	z__2.r = scltmp * cartcom_1.crfield[i__3].r, z__2.i = scltmp * 
		cartcom_1.crfield[i__3].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	do_uio(&c__2, (char *)&z__1, (ftnlen)sizeof(doublereal));
    }
    e_wdue();

    if (mpicom_1.mpi_size__ > 1) {
	cl__1.cerr = 0;
	cl__1.cunit = iocom_1.ndump;
	cl__1.csta = 0;
	f_clos(&cl__1);
	iocom_1.ndump = ndmp2tmp;
    }

/*     harmonics */

/* L3: */
    i__3 = cartcom_1.nhloop;
    for (ih = 2; ih <= i__3; ++ih) {
	j = *islice - iocom_1.firstout;
	ifile = iocom_1.ndumph[ih - 2];
	ioffset = inputcom_1.ncar * inputcom_1.ncar * (ih - 1);
	if (mpicom_1.mpi_size__ > 1) {
	    s_wsfi(&io___106);
	    do_fio(&c__1, (char *)&cartcom_1.hloop[ih - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfi();
	    ifile += 30;
	    o__1.oerr = 1;
	    o__1.ounit = ifile;
	    o__1.ofnmlen = strlen_(inputcom_1.outputfile, (ftnlen)30) + 17;
/* Writing concatenation */
	    i__4[0] = strlen_(inputcom_1.outputfile, (ftnlen)30), a__2[0] = 
		    inputcom_1.outputfile;
	    i__4[1] = 4, a__2[1] = ".dfl";
	    i__4[2] = 1, a__2[2] = harm_id__;
	    i__4[3] = 12, a__2[3] = file_id__;
	    s_cat(ch__2, a__2, i__4, &c__4, (ftnlen)47);
	    o__1.ofnm = ch__2;
	    o__1.orl = (inputcom_1.ncar << 4) * inputcom_1.ncar;
	    o__1.osta = "unknown";
	    o__1.oacc = "direct";
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__1 = f_open(&o__1);
	    if (i__1 != 0) {
		goto L100;
	    }
	    j = 1;
	}

	io___107.ciunit = ifile;
	io___107.cirec = j;
	s_wdue(&io___107);
	i__1 = inputcom_1.ncar * inputcom_1.ncar + ioffset;
	for (i__ = ioffset + 1; i__ <= i__1; ++i__) {
	    i__5 = i__ - 1;
	    z__3.r = scltmp * cartcom_1.crfield[i__5].r, z__3.i = scltmp * 
		    cartcom_1.crfield[i__5].i;
	    i__6 = ih - 1;
	    d__1 = (doublereal) cartcom_1.hloop[i__6];
	    z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    do_uio(&c__2, (char *)&z__1, (ftnlen)sizeof(doublereal));
	}
	e_wdue();

	if (mpicom_1.mpi_size__ > 1) {
	    cl__1.cerr = 0;
	    cl__1.cunit = ifile;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	}

    }

    return 0;

L100:
    i__ = printerr_(&c_n1, "MPI Binary Temp Files", (ftnlen)21);
    return 0;
} /* outdump_ */


/* Subroutine */ int outdumpslippage_()
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt();
    integer s_wdue(), do_uio(), e_wdue();

    /* Local variables */
    static integer i__, j, ifile, i0, ih, ioffset;
    static doublereal scltmp;
    extern /* Subroutine */ int pulltimerec_();

    /* Fortran I/O blocks */
    static cilist io___111 = { 0, 0, 0, 0, 0 };
    static cilist io___116 = { 0, 0, 0, 0, 0 };


/*     ================================== */
/*     dumps the escaped slippage field ahead of the bunch */
/*     --------------------------------- */




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




/*     ------------------------------------------------------------------ */
/*     input/output control */




/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */







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



    if (inputcom_1.nslice <= iocom_1.firstout) {
	return 0;
    }
/* suppress for IOTAIL=0 */
    if (inputcom_1.itdp == 0) {
	return 0;
    }
    if (inputcom_1.idmpfld == 0) {
	return 0;
    }
    if (mpicom_1.mpi_id__ > 0) {
	return 0;
    }

    scltmp = cartcom_1.dxy * 510999.06 * simcom_1.xkper0 / cartcom_1.xks / 
	    sqrt(376.73);


    for (j = tbunchcom_1.nslp - 1; j >= 1; --j) {
/* dump field , escaping beam */
	pulltimerec_(workspace_1.crwork3, &inputcom_1.ncar, &j);
	i0 = inputcom_1.nslice + tbunchcom_1.nslp - j - iocom_1.firstout;
	io___111.ciunit = iocom_1.ndump;
	io___111.cirec = i0;
	s_wdue(&io___111);
	i__1 = inputcom_1.ncar * inputcom_1.ncar;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = i__ - 1;
	    z__2.r = scltmp * workspace_1.crwork3[i__2].r, z__2.i = scltmp * 
		    workspace_1.crwork3[i__2].i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    do_uio(&c__2, (char *)&z__1, (ftnlen)sizeof(doublereal));
	}
	e_wdue();
	i__2 = cartcom_1.nhloop;
	for (ih = 2; ih <= i__2; ++ih) {
	    ifile = iocom_1.ndumph[ih - 2];
	    ioffset = (ih - 1) * inputcom_1.ncar * inputcom_1.ncar;
	    io___116.ciunit = ifile;
	    io___116.cirec = i0;
	    s_wdue(&io___116);
	    i__1 = inputcom_1.ncar * inputcom_1.ncar + ioffset;
	    for (i__ = ioffset + 1; i__ <= i__1; ++i__) {
		i__3 = i__ - 1;
		z__3.r = scltmp * workspace_1.crwork3[i__3].r, z__3.i = 
			scltmp * workspace_1.crwork3[i__3].i;
		i__4 = ih - 1;
		d__1 = (doublereal) cartcom_1.hloop[i__4];
		z__2.r = z__3.r / d__1, z__2.i = z__3.i / d__1;
		z__1.r = z__2.r, z__1.i = z__2.i;
		do_uio(&c__2, (char *)&z__1, (ftnlen)sizeof(doublereal));
	    }
	    e_wdue();
	}
    }
    return 0;
} /* outdumpslippage_ */


/* Subroutine */ int closefile_(nio)
integer *nio;
{
    /* System generated locals */
    cllist cl__1;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(), f_clos();

    /* Local variables */
    static logical isop;

/*     ================================================================= */
/*     closing file */
/*     --------------------------------------------------------------- */


    if (*nio > 6) {
	ioin__1.inerr = 0;
	ioin__1.inunit = *nio;
	ioin__1.infile = 0;
	ioin__1.inex = 0;
	ioin__1.inopen = &isop;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
	if (isop) {
	    cl__1.cerr = 0;
	    cl__1.cunit = *nio;
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	}
/* close history file */
    }
    return 0;
} /* closefile_ */


/* of closefile */
integer opentextfile_(file, status, nio, file_len, status_len)
char *file, *status;
integer *nio;
ftnlen file_len;
ftnlen status_len;
{
    /* System generated locals */
    integer ret_val, i__1;
    olist o__1;

    /* Builtin functions */
    integer f_open();

    /* Local variables */
    extern integer printerr_();

/*     ================================================================== */
/*     open ascii file (sequential access) */
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

    ret_val = *nio;
    o__1.oerr = 1;
    o__1.ounit = *nio;
    o__1.ofnmlen = file_len;
    o__1.ofnm = file;
    o__1.orl = 0;
    o__1.osta = status;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L100;
    }
    return ret_val;
L100:
    ret_val = printerr_(&c_n1, file, file_len);
    return ret_val;
} /* opentextfile_ */


/* of opentextfile */
integer openbinfile_(root, extension, nio, nsize, root_len, extension_len)
char *root, *extension;
integer *nio, *nsize;
ftnlen root_len;
ftnlen extension_len;
{
    /* System generated locals */
    address a__1[3];
    integer ret_val, i__1[3], i__2;
    olist o__1;

    /* Builtin functions */
    integer i_indx();
    /* Subroutine */ int s_cat();
    integer f_open();

    /* Local variables */
    static char filename[36];
    extern integer printerr_();
    static integer j, jj;

/*     ================================================================== */
/*     open binary file (direct access) as addition output file */
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

    ret_val = *nio;
    j = i_indx(root, " ", (ftnlen)30, (ftnlen)1);
    if (j == 0) {
	j = 31;
    }
    --j;
    jj = i_indx(extension, " ", (ftnlen)4, (ftnlen)1);
    if (jj == 0) {
	jj = 5;
    }
    --jj;
/* Writing concatenation */
    i__1[0] = j, a__1[0] = root;
    i__1[1] = 1, a__1[1] = ".";
    i__1[2] = jj, a__1[2] = extension;
    s_cat(filename, a__1, i__1, &c__3, (ftnlen)36);
    o__1.oerr = 1;
    o__1.ounit = *nio;
    o__1.ofnmlen = 36;
    o__1.ofnm = filename;
    o__1.orl = *nsize;
    o__1.osta = "unknown";
    o__1.oacc = "direct";
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__2 = f_open(&o__1);
    if (i__2 != 0) {
	goto L100;
    }
    return ret_val;
L100:
    ret_val = printerr_(&c_n1, filename, (ftnlen)36);
    return ret_val;
} /* openbinfile_ */



/* of openbinfile */
/* Subroutine */ int openoutputbinmpi_(islice)
integer *islice;
{
    /* Format strings */
    static char fmt_2[] = "(\002.slice\002,i6.6)";
    static char fmt_5[] = "(i1.1)";

    /* System generated locals */
    address a__1[3], a__2[4];
    integer i__1, i__2[3], i__3, i__4[4];
    char ch__1[46], ch__2[47];
    olist o__1;

    /* Builtin functions */
    integer s_wsfi(), do_fio(), e_wsfi(), i_indx();
    /* Subroutine */ int s_cat();
    integer f_open();

    /* Local variables */
    static integer iopenerr;
    extern integer printerr_();
    static integer i__, j;
    static char file_harm__[1];
    extern integer strlen_();
    static char file_id__[12];

    /* Fortran I/O blocks */
    static icilist io___122 = { 0, file_id__, 0, fmt_2, 12, 1 };
    static icilist io___126 = { 0, file_harm__, 0, fmt_5, 1, 1 };


/*     =========================================== */
/*     open binary file for field, particle, dump field and dump particle */
/*     ------------------------------------------ */





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




    if (mpicom_1.mpi_size__ <= 1) {
	return 0;
    }

/* no mpi operation */
    s_wsfi(&io___122);
    do_fio(&c__1, (char *)&(*islice), (ftnlen)sizeof(integer));
    e_wsfi();

    j = i_indx(inputcom_1.outputfile, " ", (ftnlen)30, (ftnlen)1);
    if (j == 0) {
	j = strlen_(inputcom_1.outputfile, (ftnlen)30) + 1;
    }
    --j;

    if (iocom_1.nfld > 0) {
	iocom_1.irecfld = 1;
/* reset record counter */
	mpicom_1.nfldmpi = iocom_1.nfld;
/* save original file ID */
	iocom_1.nfld += 30;
	o__1.oerr = 1;
	o__1.ounit = iocom_1.nfld;
	o__1.ofnmlen = j + 16;
/* Writing concatenation */
	i__2[0] = j, a__1[0] = inputcom_1.outputfile;
	i__2[1] = 4, a__1[1] = ".fld";
	i__2[2] = 12, a__1[2] = file_id__;
	s_cat(ch__1, a__1, i__2, &c__3, (ftnlen)46);
	o__1.ofnm = ch__1;
	o__1.orl = inputcom_1.ncar * inputcom_1.ncar << 4;
	o__1.osta = "unknown";
	o__1.oacc = "direct";
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L100;
	}
	i__1 = cartcom_1.nhloop;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    mpicom_1.nfldhmpi[i__ - 2] = iocom_1.nfldh[i__ - 2];
	    iocom_1.nfldh[i__ - 2] += 30;
	    s_wsfi(&io___126);
	    do_fio(&c__1, (char *)&cartcom_1.hloop[i__ - 1], (ftnlen)sizeof(
		    integer));
	    e_wsfi();
	    o__1.oerr = 1;
	    o__1.ounit = iocom_1.nfldh[i__ - 2];
	    o__1.ofnmlen = j + 17;
/* Writing concatenation */
	    i__4[0] = j, a__2[0] = inputcom_1.outputfile;
	    i__4[1] = 4, a__2[1] = ".fld";
	    i__4[2] = 1, a__2[2] = file_harm__;
	    i__4[3] = 12, a__2[3] = file_id__;
	    s_cat(ch__2, a__2, i__4, &c__4, (ftnlen)47);
	    o__1.ofnm = ch__2;
	    o__1.orl = inputcom_1.ncar * inputcom_1.ncar << 4;
	    o__1.osta = "unknown";
	    o__1.oacc = "direct";
	    o__1.ofm = 0;
	    o__1.oblnk = 0;
	    i__3 = f_open(&o__1);
	    if (i__3 != 0) {
		goto L100;
	    }
	}
    }

    if (iocom_1.npar > 0) {
	iocom_1.irecpar = 1;
	mpicom_1.nparmpi = iocom_1.npar;
	iocom_1.npar += 30;
	o__1.oerr = 1;
	o__1.ounit = iocom_1.npar;
	o__1.ofnmlen = j + 16;
/* Writing concatenation */
	i__2[0] = j, a__1[0] = inputcom_1.outputfile;
	i__2[1] = 4, a__1[1] = ".par";
	i__2[2] = 12, a__1[2] = file_id__;
	s_cat(ch__1, a__1, i__2, &c__3, (ftnlen)46);
	o__1.ofnm = ch__1;
	o__1.orl = inputcom_1.npart << 3;
	o__1.osta = "unknown";
	o__1.oacc = "direct";
	o__1.ofm = 0;
	o__1.oblnk = 0;
	i__1 = f_open(&o__1);
	if (i__1 != 0) {
	    goto L100;
	}
    }

    return 0;
L100:
    iopenerr = printerr_(&c_n1, "MPI Binary Temp Files", (ftnlen)21);
    return 0;
} /* openoutputbinmpi_ */


/* Subroutine */ int closeoutputbinmpi_()
{
    /* System generated locals */
    integer i__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_clos();

    /* Local variables */
    static integer i__;

/*     =========================================== */
/*     close binary file for field, particle, dump field and dump particle */
/*     ------------------------------------------ */





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




/*     ------------------------------------------------------------------ */
/*     input/output control */




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



    if (mpicom_1.mpi_size__ <= 1) {
	return 0;
    }

/* no mpi operation */
    if (iocom_1.nfld > 0) {
	cl__1.cerr = 0;
	cl__1.cunit = iocom_1.nfld;
	cl__1.csta = 0;
	f_clos(&cl__1);
	iocom_1.nfld = mpicom_1.nfldmpi;
	i__1 = cartcom_1.nhloop;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    cl__1.cerr = 0;
	    cl__1.cunit = iocom_1.nfldh[i__ - 2];
	    cl__1.csta = 0;
	    f_clos(&cl__1);
	    iocom_1.nfldh[i__ - 2] = mpicom_1.nfldhmpi[i__ - 2];
	}
    }

    if (iocom_1.npar > 0) {
	cl__1.cerr = 0;
	cl__1.cunit = iocom_1.npar;
	cl__1.csta = 0;
	f_clos(&cl__1);
	iocom_1.npar = mpicom_1.nparmpi;
    }

    return 0;
} /* closeoutputbinmpi_ */


/* Subroutine */ int first_()
{
    /* Format strings */
    static char fmt_100[] = "(\002-------------------------------\002,/,\002\
Genesis 1.3 has begun execution\002,/,\002(Version \002,f3.1,\002 \002,a,\
\002)\002,/)";

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Fortran I/O blocks */
    static cilist io___129 = { 0, 0, 0, fmt_100, 0 };


/*     ============================================ */
/*     initial information for user */
/*     -------------------------------------------- */





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




    io___129.ciunit = iocom_1.nlog;
    s_wsfe(&io___129);
    do_fio(&c__1, (char *)&c_b22, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, "Unix", (ftnlen)4);
    e_wsfe();
    return 0;
} /* first_ */

/* Subroutine */ int last_()
{
    /* Format strings */
    static char fmt_100[] = "(\002***  closing files\002)";
    static char fmt_200[] = "(/,\002Genesis run has finished\002,/,\002-----\
-------------------\002)";

    /* Builtin functions */
    integer s_wsfe(), e_wsfe();
    /* Subroutine */ int s_stop();

    /* Local variables */
    extern /* Subroutine */ int mpi_finalize__(), closetimerec_(), closefile_(
	    );
    static integer ih;

    /* Fortran I/O blocks */
    static cilist io___130 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___132 = { 0, 0, 0, fmt_200, 0 };


/*     ================================================================== */
/*     called at end of run. */
/*     closes all files, which must stay open during the run */
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





/*     ------------------------------------------------------------------ */
/*     common block for mpi implementation for genesis */





    io___130.ciunit = iocom_1.nlog;
    s_wsfe(&io___130);
    e_wsfe();
    closefile_(&iocom_1.nout);
/* standard output */
    closefile_(&iocom_1.nfld);
/* field output */
    for (ih = 2; ih <= 7; ++ih) {
	closefile_(&iocom_1.nfldh[ih - 1]);
/* harmonic field output */
	closefile_(&iocom_1.ndumph[ih - 1]);
/* dumped harmonic field */
    }
    closefile_(&iocom_1.npar);
/* particle output */
    closefile_(&iocom_1.nfin);
/* field input */
    closefile_(&iocom_1.npin);
/* particle input */
    closefile_(&iocom_1.ndump);
/* dumped field */
    closefile_(&iocom_1.ndmp2);
/* dumped particle */
    closefile_(&iocom_1.ndis);
/* input distribution */
    closetimerec_();

    io___132.ciunit = iocom_1.nlog;
    s_wsfe(&io___132);
    e_wsfe();

/* genesis has finished */
    if (iocom_1.nlog != 6) {
	closefile_(&iocom_1.nlog);
    }
/* log file */
    mpi_finalize__(&mpicom_1.mpi_err__);
    s_stop("", (ftnlen)0);
	
	return 0; //OC port
} /* last_ */


/* last */
integer printerr_(ierr, text, text_len)
integer *ierr;
char *text;
ftnlen text_len;
{
    /* Format strings */
    static char fmt_100[] = "(\002***  File-error: \002,a,/,\002***  cannot \
be opened\002)";
    static char fmt_101[] = "(\002***  File-error: \002,a,/,\002***  cannot \
be accessed\002)";
    static char fmt_102[] = "(\002***  File-error: \002,a,/,\002***  cannot \
be opened\002,/,\002***  creating template file: template.in\002)";
    static char fmt_103[] = "(\002***  File-error: \002,a,/,\002***  error i\
n namelise $newrun\002)";
    static char fmt_104[] = "(\002***  Scan-warning: conflict with ITDP\002,\
/,\002***  using scan-feature\002,a)";
    static char fmt_105[] = "(\002***  Scan-warning: conflict with BEAMFIL\
E\002,/,\002***  ignoring BEAMFILE: \002,a)";
    static char fmt_106[] = "(\002***  Beamfile-warning: size exceeds NSMA\
X\002,/,\002***  \002,a)";
    static char fmt_107[] = "(\002***  Input-error: \002,a)";
    static char fmt_108[] = "(\002***  Input-warning: \002,a)";
    static char fmt_109[] = "(\002***  Input-error: cannot convert to indivi\
diual input\002,/,\002***  \002,a)";
    static char fmt_110[] = "(\002***  Numerical-error: boundary exceeded o\
f\002,a,/,\002***  ignoring exceeding elements\002)";
    static char fmt_111[] = "(\002***  Round-warning: section not multiple o\
f XLAMD\002,/,\002***  MAGINFILE:\002,a)";
    static char fmt_112[] = "(\002***  Extrapolation-warning: exceeding time\
 window of\002,/,\002***  BEAMFILE by:\002,a)";
    static char fmt_113[] = "(\002***  Scan-error: conflict with MAGINFILE\
:\002,a,/,\002***  disabling scan-feature\002)";
    static char fmt_114[] = "(\002***  Scan-error: conflict with MAGOUTFILE\
:\002,a,/,\002***  disabling scan-feature\002)";
    static char fmt_115[] = "(\002***  Warning: particle loss of \002,a\
,\002%\002)";
    static char fmt_116[] = "(\002***  Warning: external magnet definition t\
oo short for \002,a)";
    static char fmt_117[] = "(\002***  Error: invalid filename:\002,a)";
    static char fmt_118[] = "(\002***  File-error: cannot read from FIELDFIL\
E:\002,a)";
    static char fmt_119[] = "(\002***  Warning: \002,a)";
    static char fmt_120[] = "(\002***  Error: cannot run in background mode\
.\002,/,\002***  information needed for \002,a)";
    static char fmt_121[] = "(\002***  Error: CRTIME cannot hold slippage fi\
eld.\002,/,\002***  see manual for allocating more memory\002,a)";
    static char fmt_122[] = "(\002***  Error: unphysical parameter for loadi\
ng\002,/,\002***  \002,a)";
    static char fmt_123[] = "(\002***  \002,a)";

    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Fortran I/O blocks */
    static cilist io___133 = { 0, 0, 0, fmt_100, 0 };
    static cilist io___134 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___135 = { 0, 0, 0, fmt_102, 0 };
    static cilist io___136 = { 0, 0, 0, fmt_103, 0 };
    static cilist io___137 = { 0, 0, 0, fmt_104, 0 };
    static cilist io___138 = { 0, 0, 0, fmt_105, 0 };
    static cilist io___139 = { 0, 0, 0, fmt_106, 0 };
    static cilist io___140 = { 0, 0, 0, fmt_107, 0 };
    static cilist io___141 = { 0, 0, 0, fmt_108, 0 };
    static cilist io___142 = { 0, 0, 0, fmt_109, 0 };
    static cilist io___143 = { 0, 0, 0, fmt_110, 0 };
    static cilist io___144 = { 0, 0, 0, fmt_111, 0 };
    static cilist io___145 = { 0, 0, 0, fmt_112, 0 };
    static cilist io___146 = { 0, 0, 0, fmt_113, 0 };
    static cilist io___147 = { 0, 0, 0, fmt_114, 0 };
    static cilist io___148 = { 0, 0, 0, fmt_115, 0 };
    static cilist io___149 = { 0, 0, 0, fmt_116, 0 };
    static cilist io___150 = { 0, 0, 0, fmt_117, 0 };
    static cilist io___151 = { 0, 0, 0, fmt_118, 0 };
    static cilist io___152 = { 0, 0, 0, fmt_119, 0 };
    static cilist io___153 = { 0, 0, 0, fmt_120, 0 };
    static cilist io___154 = { 0, 0, 0, fmt_121, 0 };
    static cilist io___155 = { 0, 0, 0, fmt_122, 0 };
    static cilist io___156 = { 0, 0, 0, fmt_123, 0 };


/*     ======================================================== */
/*     print error messages */
/*     -------------------------------------------------------- */




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




    ret_val = *ierr;
    if (*ierr >= 0) {
	return ret_val;
    }
    switch ((int)(abs(*ierr))) {
	case 1:  goto L10;
	case 2:  goto L11;
	case 3:  goto L12;
	case 4:  goto L13;
	case 5:  goto L14;
	case 6:  goto L15;
	case 7:  goto L16;
	case 8:  goto L17;
	case 9:  goto L18;
	case 10:  goto L19;
	case 11:  goto L20;
	case 12:  goto L21;
	case 13:  goto L22;
	case 14:  goto L23;
	case 15:  goto L24;
	case 16:  goto L25;
	case 17:  goto L26;
	case 18:  goto L27;
	case 19:  goto L28;
	case 20:  goto L29;
	case 21:  goto L30;
	case 22:  goto L31;
	case 23:  goto L32;
	case 24:  goto L33;
    }
L10:
    io___133.ciunit = iocom_1.nlog;
    s_wsfe(&io___133);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L11:
    io___134.ciunit = iocom_1.nlog;
    s_wsfe(&io___134);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L12:
    io___135.ciunit = iocom_1.nlog;
    s_wsfe(&io___135);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L13:
    io___136.ciunit = iocom_1.nlog;
    s_wsfe(&io___136);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L14:
    io___137.ciunit = iocom_1.nlog;
    s_wsfe(&io___137);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L15:
    io___138.ciunit = iocom_1.nlog;
    s_wsfe(&io___138);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L16:
    io___139.ciunit = iocom_1.nlog;
    s_wsfe(&io___139);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L17:
    io___140.ciunit = iocom_1.nlog;
    s_wsfe(&io___140);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L18:
    io___141.ciunit = iocom_1.nlog;
    s_wsfe(&io___141);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L19:
    io___142.ciunit = iocom_1.nlog;
    s_wsfe(&io___142);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L20:
    io___143.ciunit = iocom_1.nlog;
    s_wsfe(&io___143);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L21:
    io___144.ciunit = iocom_1.nlog;
    s_wsfe(&io___144);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L22:
    io___145.ciunit = iocom_1.nlog;
    s_wsfe(&io___145);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L23:
    io___146.ciunit = iocom_1.nlog;
    s_wsfe(&io___146);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L24:
    io___147.ciunit = iocom_1.nlog;
    s_wsfe(&io___147);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L25:
    io___148.ciunit = iocom_1.nlog;
    s_wsfe(&io___148);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L26:
    io___149.ciunit = iocom_1.nlog;
    s_wsfe(&io___149);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L27:
    io___150.ciunit = iocom_1.nlog;
    s_wsfe(&io___150);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L28:
    io___151.ciunit = iocom_1.nlog;
    s_wsfe(&io___151);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L29:
    io___152.ciunit = iocom_1.nlog;
    s_wsfe(&io___152);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L30:
    io___153.ciunit = iocom_1.nlog;
    s_wsfe(&io___153);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L31:
    io___154.ciunit = iocom_1.nlog;
    s_wsfe(&io___154);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L32:
    io___155.ciunit = iocom_1.nlog;
    s_wsfe(&io___155);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;
L33:
    io___156.ciunit = iocom_1.nlog;
    s_wsfe(&io___156);
    do_fio(&c__1, text, text_len);
    e_wsfe();
    return ret_val;

/*     format statements */

} /* printerr_ */

