/************************************************************************//**
 * File: srsase.cpp
 * Description: Wrapper class for calculation of SASE using GENESIS (F2C-ed), with input and output data in SRW formats
 * Project: Synchrotron Radiation Workshop
 * First release: 2000
 *
 * Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * All Rights Reserved
 *
 * @author O.Chubar, P.Elleaume
 * @version 1.0
 ***************************************************************************/

#include "srsase.h"
#include "srmlttsk.h"
#include "gmfft.h"
#include "sroptdrf.h"

/* F2C: IMT 10Sep95  Declare jump buffer used to recover from exception exits & aborts */
#if defined(__MWERKS__) || defined(__MAC__) || defined(TPM_F2C) || defined(SPM_F2C) || defined(CW_F2C_MAC) || defined(CW_F2C_WIN32) 
	#include <setjmp.h>
	#ifdef __cplusplus
	extern "C" {
	#endif
	jmp_buf gRecoverToConsole;
	#ifdef __cplusplus
	}
	#endif
#endif /* Macintosh C compilers and MW Win32 */

//*************************************************************************

extern srTIntVect gVectWarnNos;

//*************************************************************************

extern "C" {

	int readin_();
	int chk_input__();
	int chk_bnd__();
	int outglob_();
	int preset_();
	int initio_();
	int initrun_();
	int doscan_(f2c_integer*);
	int dotime_(f2c_integer*);
	int loadrad_(f2c_integer*);
	int loadbeam_(f2c_integer*, f2c_doublereal*);
	int chk_loss__();
	int output_(f2c_integer*, f2c_integer*, f2c_doublereal*);
	int stepz_(f2c_integer*, f2c_doublereal*);
	int swapfield_(f2c_integer*);
	int diagno_(f2c_integer*);
	int rpos_(f2c_doublereal*, f2c_doublereal*, f2c_doublereal*);
	int getpsi_(f2c_doublereal*);
	int shotnoise_penman__();
	int loadquiet_(f2c_integer*);
	int loaddist_(f2c_integer*, f2c_integer*);
	int shotnoise_fawley__();
	f2c_doublereal hammv_(f2c_integer*);
	int importdispersion_();
	//int loadslpfld_(f2c_integer*);
	int getdiag_(f2c_doublereal*, f2c_doublereal*, f2c_doublereal*);
	int magfield_(f2c_doublereal*, f2c_integer*);
	int scaninit_();
	f2c_doublereal bessj0_(f2c_doublereal*);
	f2c_doublereal bessj1_(f2c_doublereal*);
	f2c_doublereal ran1_(f2c_integer*);
    //int last_();
	//int gauss_hermite__(f2c_doublecomplex* cfld, f2c_doublereal* power, f2c_doublereal* zr, f2c_doublereal* zw, f2c_doublereal* rks, f2c_doublereal* phase);
	int gauss_hermite__(f2c_doublecomplex* cfld, f2c_doublereal* power, f2c_doublereal* zr, f2c_doublereal* zw, f2c_doublereal* rks, f2c_doublereal* phase, f2c_integer* harm);
	int dotimerad_(f2c_integer* islice);
	int pushtimerec_(f2c_doublecomplex* cpush, f2c_integer* n, f2c_integer* irec);
	void d_cnjg(f2c_doublecomplex*, f2c_doublecomplex*);
	int mpi_init__(f2c_integer*);
	int mpi_comm_rank__(f2c_integer*, f2c_integer*, f2c_integer*);
	int mpi_comm_size__(f2c_integer*, f2c_integer*, f2c_integer*);
	int mpi_merge__();
	int importtransfer_();

	//int genesis_();
	//int auxval_();
	//int loading_();

	//extern struct {
	//	f2c_doublereal aw0,	xkx, xky, wcoefz[3], xlamd,	fbess0,	delaw, awd,	awx, awy, 
	//		gamma0,	delgam,	rxbeam,	rybeam,	alphax,	alphay,	emitx, emity, 
	//		xbeam, ybeam, pxbeam, pybeam, cuttail, curpeak,	conditx, condity, 
	//		bunch, bunchphase, emod, emodphase,	xlamds,	prad0, zrayl, zwaist, 
	//		rmax0, zsep, delz, zstop, quadf, quadd,	fl,	dl,	drl, f1st, qfdx, 
	//		qfdy, sl, solen, curlen, shotnoise,	svar, dgrid, eloss,	version, 
	//		ibfield, imagl,	idril;
	//	f2c_integer	iseed, nwig, nsec, npart, ncar,	lbc, nscr, nscz, nptr, ildgam, 
	//		ildpsi,	ildx, ildy,	ildpx, ildpy, itgaus, nbins, iphsty, ishsty, 
	//		ippart,	ispart,	ipradi,	isradi,	iertyp,	iwityp,	idump, iotail, 
	//		nharm, magin, magout, lout[35],	ffspec,	ntail, nslice, iall, itdp,
	//		ipseed,	iscan,	nscan, isntyp, isravg, isrsig, iorb, ndcut,	
	//		idmppar, idmpfld, ilog,	igamgaus, convharm,	alignradf, offsetradf,
	//		multconv;
	//	char beamfile[30], fieldfile[30], maginfile[30], magoutfile[30], 
	//		outputfile[30], inputfile[30], scan[30], distfile[30], partfile[30], filetype[30], radfile[30];
	//} inputcom_;
	extern struct {
		f2c_doublereal aw0, xkx, xky, wcoefz[3], xlamd, fbess0, delaw, awd, awx, awy, 
			gamma0, delgam, rxbeam, rybeam, alphax, alphay, emitx, emity, 
			xbeam, ybeam, pxbeam, pybeam, cuttail, curpeak, conditx, condity, 
			bunch, bunchphase, emod, emodphase, xlamds, prad0, zrayl, zwaist, 
			rmax0, zsep, delz, zstop, quadf, quadd, fl, dl, drl, f1st, qfdx, 
			qfdy, sl, solen, curlen, shotnoise, svar, dgrid, eloss, version, 
			ibfield, imagl, idril, igamref, pradh0, 
			itram11, itram12, itram13, itram14, itram15, itram16, 
			itram21, itram22, itram23, itram24, itram25, itram26, 
			itram31, itram32, itram33, itram34, itram35, itram36, 
			itram41, itram42, itram43, itram44, itram45, itram46, 
			itram51, itram52, itram53, itram54, itram55, itram56, 
			itram61, itram62, itram63, itram64, itram65, itram66, 
			rmax0sc;
		f2c_integer iseed, nwig, nsec, npart, ncar, lbc, nscr, nscz, nptr, ildgam, 
			ildpsi, ildx, ildy, ildpx, ildpy, itgaus, nbins, iphsty, ishsty, 
			ippart, ispart, ipradi, isradi, iertyp, iwityp, idump, iotail, 
			nharm, iallharm, iharmsc, magin, magout, lout[40], ffspec, ntail, 
			nslice, iall, itdp, ipseed, iscan, nscan, isntyp, isravg, isrsig, 
			iorb, ndcut, idmppar, idmpfld, ilog, igamgaus, convharm, 
			alignradf, offsetradf, multconv, trama, iscrkup, inverfc;
		char beamfile[30], fieldfile[30], maginfile[30], magoutfile[30], 
			outputfile[30], inputfile[30], scan[30], distfile[30], partfile[30], filetype[30], radfile[30];
	} inputcom_;

	//extern struct {
	//	f2c_doublereal distversion, distrev;
	//	f2c_integer iout[24], nout, nfld, npar, ndump, firstout, nfin, irecpar, 
	//		irecfld, kout, ndmp2, npin, nlog, ndis, ncoldis, iconv2t, iconv2g,
	//		icolpar[10], ndistsize, iconv2px, iconv2py, nprobe, ftype, 
	//		ftdist, ftpart, ftfield;
	//} iocom_;
	extern struct {
		f2c_doublereal distversion, distrev;
		f2c_integer iout[39], nout, nfld, npar, ndump, firstout, nfin, irecpar, 
			irecfld, kout, ndmp2, npin, nlog, ndis, ncoldis, iconv2t, iconv2g,
			icolpar[10], ndistsize, iconv2px, iconv2py, nprobe, ftype, 
			ftdist, ftpart, ftfield, ndumph[6], nfldh[6];
	} iocom_;

	//extern struct {
	//	//f2c_doublereal tgam0[30000], tdgam[30000], temitx[30000], temity[30000], txrms[30000], tyrms[30000], txpos[30000], 
	//	//	typos[30000], tpxpos[30000], tpypos[30000], talphx[30000], talphy[30000], tcurrent[30000], tpos[30000], tloss[30000], 
	//	//	distgam[250000], distx[250000],	disty[250000], distpx[250000], distpy[250000], distt[250000], 
	//	//	tradpos[30000], tzrayl[30000], tzwaist[30000], tprad0[30000], tradphase[30000];
	//	f2c_doublereal *tgam0, *tdgam, *temitx, *temity, *txrms, *tyrms, *txpos, 
	//		*typos, *tpxpos, *tpypos, *talphx, *talphy, *tcurrent, *tpos, *tloss, 
	//		*distgam, *distx, *disty, *distpx, *distpy, *distt,
	//		*tradpos, *tzrayl, *tzwaist, *tprad0, *tradphase;
	//	f2c_integer ndata, nsep, nslp, ndist, nraddata;
	//} tbunchcom_;
	extern struct {
		//doublereal tgam0[50000], tdgam[50000], temitx[50000], temity[50000], 
		// txrms[50000], tyrms[50000], txpos[50000], typos[50000], tpxpos[
		// 50000], tpypos[50000], talphx[50000], talphy[50000], tcurrent[
		// 50000], tpos[50000], tloss[50000], distgam[1250000], distx[
		// 1250000], disty[1250000], distpx[1250000], distpy[1250000], distt[
		// 1250000], tradpos[50000], tzrayl[50000], tzwaist[50000], tprad0[
		// 50000], tradphase[50000];
		f2c_doublereal *tgam0, *tdgam, *temitx, *temity, *txrms, *tyrms, *txpos, 
			*typos, *tpxpos, *tpypos, *talphx, *talphy, *tcurrent, *tpos, *tloss, 
			*distgam, *distx, *disty, *distpx, *distpy, *distt, 
			*tradpos, *tzrayl, *tzwaist, *tprad0, *tradphase;
		f2c_integer ndata, nsep, nslp, ndist, nraddata;
	} tbunchcom_;

	//extern struct {
	//	//f2c_doublecomplex crfield[29241], crsource[29241], crmatc[171], cstep, crhm[29241], cwet[171], cbet[171];
	//	f2c_doublecomplex *crfield, *crsource, crmatc[513], cstep, *crhm, cwet[513], cbet[513];
	//	f2c_doublereal dxy, xks, radphase;
	//} cartcom_;
	extern struct {
		//doublecomplex crfield[476847], crsource[68121], crmatc[1827], cstep[7], crhm[68121], cwet[1827], cbet[1827];
		f2c_doublecomplex *crfield, *crsource, *crmatc, cstep[7], *crhm, *cwet, *cbet;
		f2c_doublereal dxy, xks, radphase, besselcoupling[7];
		f2c_integer nhloop, hloop[7];
	} cartcom_;

	//extern struct {
	//	f2c_doublereal xkw0, xkper0, sval, svalout, gamma0_in__;
	//	f2c_integer npart0, inorun;
	//} simcom_;
	extern struct {
		f2c_doublereal xkw0, xkper0, sval, svalout, gamma0_in__;
		f2c_integer npart0, inorun;
	} simcom_;

	//extern struct {
	//	//f2c_doublereal xpart[70001], ypart[70001], px[70001], py[70001], gamma[70001],
	//	//	theta[70001], xporb[70001], yporb[70001], btpar[70001], btper[70001], ez[70001], wx[70001], wy[70001], 
	//	//	xcuren, dedz, tdmin, tdmax, delcharge, dtd, charge;
	//	//f2c_integer lostid[70001], lost, losttot, ipos[280004];	/* was [4][70001] */
	//	f2c_doublereal *xpart, *ypart, *px, *py, *gamma,
	//		*theta, *xporb, *yporb, *btpar, *btper, *ez, *wx, *wy, 
	//		xcuren, dedz, tdmin, tdmax, delcharge, dtd, charge;
	//	f2c_integer *lostid, lost, losttot, *ipos;
	//} beamcom_;
	extern struct {
		//doublereal xpart[1000001], ypart[1000001], px[1000001], py[1000001], 
		// gamma[1000001], theta[1000001], xporb[1000001], yporb[1000001], 
		// btpar[1000001], btper[1000001], ez[1000001], wx[1000001], wy[
		// 1000001], xcuren, dedz, tdmin, tdmax, delcharge, dtd, charge;
		//integer lostid[1000001], lost, losttot, ipos[4000004]	/* was [4][1000001] */;
		f2c_doublereal *xpart, *ypart, *px, *py, *gamma, 
			*theta, *xporb, *yporb, *btpar, *btper, *ez, *wx, *wy, 
			xcuren, dedz, tdmin, tdmax, delcharge, dtd, charge;
		f2c_integer *lostid, lost, losttot, *ipos;
	} beamcom_;

	//extern struct {
	//	//f2c_doublereal awz[10001], awdz[10000], solz[10000], awerx[10000], awery[10000], qfld[10001], 
	//	//	fbess, magversion, unitlength, dqfx[10001], dqfy[10001], awdx[10001], awdy[10001];
	//	f2c_doublereal *awz, *awdz, *solz, *awerx, *awery, *qfld, 
	//		fbess, magversion, unitlength, *dqfx, *dqfy, *awdx, *awdy;
	//	f2c_integer iran, nstepz, itap;
	//} wigcom_;
	extern struct {
		//doublereal awz[10001], awdz[10000], solz[10000], awerx[10000], awery[
		// 10000], qfld[10001], fbess, magversion, unitlength, dqfx[10001], 
		// dqfy[10001], awdx[10001], awdy[10001];
		f2c_doublereal *awz, *awdz, *solz, *awerx, *awery, *qfld, 
			fbess, magversion, unitlength, *dqfx, *dqfy, *awdx, *awdy;
		f2c_integer iran, nstepz, itap;
	} wigcom_;

	//extern struct {
	//	//f2c_doublereal error[10000], gain[10000], phimid[10000], whalf[10000], logp[10000], powmid[10000], xrms[10000], yrms[10000], xpos[10000], 
	//	//	ypos[10000], pmodhist[50000] /* was [5][10000] */, gamhist[10000], diver[10000], pradol, pinit, 
	//	//	bunphase[50000]	/* was [5][10000] */, dgamhist[10000], ffield[10000];
	//	f2c_doublereal *error, *gain, *phimid, *whalf, *logp, *powmid, *xrms, *yrms, *xpos, 
	//		*ypos, *pmodhist, *gamhist, *diver, pradol, pinit, 
	//		*bunphase, *dgamhist, *ffield;
	//	f2c_integer ihist;
	//} diagcom_;
	extern struct {
		//doublereal error[10000], gain[10000], phimid[70000]	/* was [7][10000] */, 
		// whalf[10000], logp[10000], powmid, xrms[10000], yrms[10000], xpos[
		// 10000], ypos[10000], pmodhist[70000]	/* was [7][10000] */, 
		// gamhist[10000], diver[10000], pradol, pinit, bunphase[70000]	
		// /* was [7][10000] */, dgamhist[10000], ffield[70000]	/* 
		// was [7][10000] */, pradoln[7], pgainhist[70000]	/* was [7][
		// 10000] */, pmidhist[70000]	/* was [7][10000] */;
		f2c_doublereal *error, *gain, *phimid, *whalf, *logp, powmid, *xrms, *yrms, *xpos, 
			*ypos, *pmodhist, *gamhist, *diver, pradol, pinit, 
			*bunphase, *dgamhist, *ffield, pradoln[7], 
			*pgainhist, *pmidhist;
		f2c_integer ihist;
	} diagcom_;

	//extern struct {
	//	//f2c_doublecomplex crwork3[116964], cpart1[70001], cpart2[70001];
	//	//f2c_doublereal k2gg[70001], k2pp[70001], k3gg[70001], k3pp[70001], p1[70001], p2[70001];
	//	f2c_doublecomplex *crwork3, *cpart1, *cpart2;
	//	f2c_doublereal *k2gg, *k2pp, *k3gg, *k3pp, *p1, *p2;
	//} workspace_;
	extern struct {
		//doublecomplex crwork3[1907388], cpart1[7000007], cpart2[1000001], cpart3[1000001];
		//doublereal k2gg[1000001], k2pp[1000001], k3gg[1000001], k3pp[1000001], p1[1000001], p2[1000001];
		//integer iwork[1000001];
		f2c_doublecomplex *crwork3, *cpart1, *cpart2, *cpart3;
		f2c_doublereal *k2gg, *k2pp, *k3gg, *k3pp, *p1, *p2;
		f2c_integer *iwork;
	} workspace_;

	//extern struct {
	//	//f2c_doublecomplex crtime[15960700];
	//	f2c_doublecomplex *crtime;
	//	f2c_integer ntmp;
	//} tslipcom_;
	extern struct {
		//doublecomplex crtime[79803500];
		f2c_doublecomplex *crtime;
		f2c_integer ntmp;
	} tslipcom_;

	extern struct {
		f2c_integer mpi_id__, mpi_err__, mpi_size__, mpi_loop__, nfldmpi, nparmpi, nfldhmpi[6];
	} mpicom_;
}

#define inputcom_1 inputcom_
#define iocom_1 iocom_
#define tbunchcom_1 tbunchcom_
#define cartcom_1 cartcom_
#define simcom_1 simcom_
#define beamcom_1 beamcom_
#define wigcom_1 wigcom_
#define diagcom_1 diagcom_
#define workspace_1 workspace_
#define tslipcom_1 tslipcom_
#define mpicom_1 mpicom_

//#define outputcom_1 outputcom_
//#define timecom_1 timecom_
//#define loadcom_1 loadcom_
//#define radiatcom_1 radiatcom_
//#define scancom_1 scancom_
//#define cartcom_1 cartcom_

extern srTYield srYield;
extern srTIntVect gVectWarnNos;

//*************************************************************************

double srTSASE::EEV = 510998.902; //511003.4E0;//510999.06E0; //Energy units (mc^2) in eV
double srTSASE::VACIMP = 376.7303135; //Vacuum impedance in Ohms

double srTSASE::PlankConstReduced = 1.0545887E-34; //Reduced Plank constant in J*s
double srTSASE::SpeedOfLight = 2.99792458E08; // in m/s

double srTSASE::PI = 3.14159265358979;
double srTSASE::TWOPI = 6.28318530717958;//2.*(srTSASE::PI);// = 2.E0*PI

//modification of these parameters doesn't ensure that everything 
//would work consistently in the f2c-ed part
long srTSASE::NPMAX = 300001; //70001; //65536; //# of particles
long srTSASE::NDMAX = 1250000; //70001; //65536; //# of particles
long srTSASE::NHMAX = 7; //max. number of harmonics
long srTSASE::NCMAX = 201; //171 //513; //# of gridpoints of cartesian mesh
long srTSASE::NRGRID = 400; //100; //Number of radial points for s. c.
long srTSASE::NZMAX = 10000; //# of integration steps
long srTSASE::NSMAX = 50000; //# of slices

//*************************************************************************

double srTSASE::TDElecFieldConvConstSRW2GENESIS()
{//This gives a constant to be used for multiplication of SRW TD Electric Field, assumed to be in units of sqrt(W/mm^2),
	//in order to obtain the Electric Field in internal GENESIS's units.
	double resConst = (1E+03)*cartcom_.xks*sqrt(VACIMP)/(simcom_.xkper0*simcom_.xkper0*EEV);
	//double resConst2 = (1E+03)*sqrt(SeedRad.xStep*SeedRad.zStep)*cartcom_.xks*sqrt(VACIMP)/(cartcom_.dxy*simcom_.xkper0*EEV);
	return resConst;
}

//*************************************************************************

double srTSASE::RadFieldMultip()
{// after multiplication by this value, GENESIS's radiation field should be in units of 
 // sqrt(Photons/s/0.1%bw/mm^2)
 // ATTENTION: should be used only after GENESIS's common blocks are set up !!!

	//double RadWavenum = TWOPI/radiatcom_1.xlamds;
	//double UndWavenum = TWOPI/wigcom_1.xlamd;

	double RadWavenum = TWOPI/inputcom_1.xlamds;
	double UndWavenum = TWOPI/inputcom_1.xlamd;
	double AuxMult = EEV*UndWavenum*UndWavenum/RadWavenum;
	double PowerMultSI = AuxMult*AuxMult/(RadWavenum*SpeedOfLight*PlankConstReduced*VACIMP);
	double PowerMultPract = PowerMultSI*(1.E-9); // (0.1%bw), per mm^2
	double FieldMult = sqrt(PowerMultPract);
	return FieldMult;
}

//*************************************************************************

void srTSASE::ZeroPtrs_tbunchcom()
{
	tbunchcom_.tgam0 = NULL;
	tbunchcom_.tdgam = NULL;
	tbunchcom_.temitx = NULL;
	tbunchcom_.temity = NULL;
	tbunchcom_.txrms = NULL;
	tbunchcom_.tyrms = NULL;
	tbunchcom_.txpos = NULL;
	tbunchcom_.typos = NULL;
	tbunchcom_.tpxpos = NULL;
	tbunchcom_.tpypos = NULL;
	tbunchcom_.talphx = NULL;
	tbunchcom_.talphy = NULL;
	tbunchcom_.tcurrent = NULL;
	tbunchcom_.tpos = NULL;
	tbunchcom_.tloss = NULL;
	tbunchcom_.tradpos = NULL;
	tbunchcom_.tzrayl = NULL;
	tbunchcom_.tzwaist = NULL;
	tbunchcom_.tprad0 = NULL;
	tbunchcom_.tradphase = NULL;
	tbunchcom_.distgam = NULL;
	tbunchcom_.distx = NULL;
	tbunchcom_.disty = NULL;
	tbunchcom_.distpx = NULL;
	tbunchcom_.distpy = NULL;
	tbunchcom_.distt = NULL;
}

void srTSASE::ZeroPtrs_cartcom()
{
	cartcom_.crfield = NULL;
	cartcom_.crsource = NULL;
	cartcom_.crhm = NULL;

	cartcom_.crmatc = NULL;
	cartcom_.cwet = NULL;
	cartcom_.cbet = NULL;
}

void srTSASE::ZeroPtrs_beamcom()
{
	beamcom_.xpart = NULL;
	beamcom_.ypart = NULL;
	beamcom_.px = NULL;
	beamcom_.py = NULL;
	beamcom_.gamma = NULL;
	beamcom_.theta = NULL;
	beamcom_.xporb = NULL;
	beamcom_.yporb = NULL;
	beamcom_.btpar = NULL;
	beamcom_.btper = NULL;
	beamcom_.ez = NULL;
	beamcom_.wx = NULL;
	beamcom_.wy = NULL;
	beamcom_.lostid = NULL;
	beamcom_.ipos = NULL;
}

void srTSASE::ZeroPtrs_wigcom()
{
	wigcom_.awz = NULL;
	wigcom_.awdz = NULL;
	wigcom_.solz = NULL;
	wigcom_.awerx = NULL;
	wigcom_.awery = NULL;
	wigcom_.qfld = NULL;
	wigcom_.dqfx = NULL;
	wigcom_.dqfy = NULL;
	wigcom_.awdx = NULL;
	wigcom_.awdy = NULL;
}

void srTSASE::ZeroPtrs_diagcom()
{
	diagcom_.error = NULL;
	diagcom_.gain = NULL;
	diagcom_.phimid = NULL;
	diagcom_.whalf = NULL;
	diagcom_.logp = NULL;
	//diagcom_.powmid = NULL;
	diagcom_.xrms = NULL;
	diagcom_.yrms = NULL;
	diagcom_.xpos = NULL;
	diagcom_.ypos = NULL;
	diagcom_.pmodhist = NULL;
	diagcom_.gamhist = NULL;
	diagcom_.diver = NULL;
	diagcom_.bunphase = NULL;
	diagcom_.dgamhist = NULL;
	diagcom_.ffield = NULL;

	diagcom_.pgainhist = NULL;
	diagcom_.pmidhist = NULL;
}

void srTSASE::ZeroPtrs_workspace()
{
	workspace_.crwork3 = NULL;
	workspace_.cpart1 = NULL;
	workspace_.cpart2 = NULL;
	workspace_.cpart3 = NULL;

	workspace_.k2gg = NULL;
	workspace_.k2pp = NULL;
	workspace_.k3gg = NULL;
	workspace_.k3pp = NULL;
	workspace_.p1 = NULL;
	workspace_.p2 = NULL;

	workspace_.iwork = NULL;
}

void srTSASE::ZeroPtrs_tslipcom()
{
	tslipcom_.crtime = NULL;
}

//*************************************************************************

int srTSASE::Alloc_tbunchcom()
{
	//int n1 = 30000, n2 = 250000;
	//int n1 = 50000, n2 = 1250000;
	int n1 = NSMAX, n2 = NDMAX;
	if((tbunchcom_.tgam0 = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tgam0, n1);
	if((tbunchcom_.tdgam = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tdgam, n1);
	if((tbunchcom_.temitx = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.temitx, n1);
	if((tbunchcom_.temity = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.temity, n1);
	if((tbunchcom_.txrms = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.txrms, n1);
	if((tbunchcom_.tyrms = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tyrms, n1);
	if((tbunchcom_.txpos = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.txpos, n1);
	if((tbunchcom_.typos = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.typos, n1);
	if((tbunchcom_.tpxpos = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tpxpos, n1);
	if((tbunchcom_.tpypos = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tpypos, n1);
	if((tbunchcom_.talphx = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.talphx, n1);
	if((tbunchcom_.talphy = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.talphy, n1);
	if((tbunchcom_.tcurrent = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tcurrent, n1);
	if((tbunchcom_.tpos = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tpos, n1);
	if((tbunchcom_.tloss = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tloss, n1);
	if((tbunchcom_.tradpos = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tradpos, n1);
	if((tbunchcom_.tzrayl = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tzrayl, n1);
	if((tbunchcom_.tzwaist = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tzwaist, n1);
	if((tbunchcom_.tprad0 = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tprad0, n1);
	if((tbunchcom_.tradphase = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.tradphase, n1);

	if((tbunchcom_.distgam = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.distgam, n2);
	if((tbunchcom_.distx = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.distx, n2);
	if((tbunchcom_.disty = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.disty, n2);
	if((tbunchcom_.distpx = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.distpx, n2);
	if((tbunchcom_.distpy = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.distpy, n2);
	if((tbunchcom_.distt = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tbunchcom_.distt, n2);
	return 0;
}

//int srTSASE::Alloc_cartcom(int numHarm)
int srTSASE::Alloc_cartcom()
{
	int n0 = NCMAX*NHMAX;
	int n1 = NCMAX*NCMAX*NHMAX; //476847;
	//int n0 = NCMAX*numHarm; //NCMAX*NHMAX;
	//int n1 = NCMAX*NCMAX*numHarm; //NCMAX*NCMAX*NHMAX; //476847;
	int n2 = NCMAX*NCMAX; //68121;

	if((cartcom_.crfield = new f2c_doublecomplex[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(cartcom_.crfield, n1);
	if((cartcom_.crsource = new f2c_doublecomplex[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(cartcom_.crsource, n2);
	if((cartcom_.crhm = new f2c_doublecomplex[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(cartcom_.crhm, n2);
	
	if((cartcom_.crmatc = new f2c_doublecomplex[n0]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(cartcom_.crmatc, n0);
	if((cartcom_.cwet = new f2c_doublecomplex[n0]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(cartcom_.cwet, n0);
	if((cartcom_.cbet = new f2c_doublecomplex[n0]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(cartcom_.cbet, n0);

	return 0;
}

int srTSASE::Alloc_beamcom()
{
	//int n1 = 70001, n2 = 280004;
	//int n1 = 1000001, n2 = 4000004;

	if((beamcom_.xpart = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.xpart, NPMAX);
	if((beamcom_.ypart = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.ypart, NPMAX);
	if((beamcom_.px = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.px, NPMAX);
	if((beamcom_.py = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.py, NPMAX);
	if((beamcom_.gamma = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.gamma, NPMAX);
	if((beamcom_.theta = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.theta, NPMAX);
	if((beamcom_.xporb = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.xporb, NPMAX);
	if((beamcom_.yporb = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.yporb, NPMAX);
	if((beamcom_.btpar = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.btpar, NPMAX);
	if((beamcom_.btper = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.btper, NPMAX);
	if((beamcom_.ez = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.ez, NPMAX);
	if((beamcom_.wx = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.wx, NPMAX);
	if((beamcom_.wy = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.wy, NPMAX);
	if((beamcom_.lostid = new f2c_integer[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.lostid, NPMAX);
	if((beamcom_.ipos = new f2c_integer[4*NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(beamcom_.ipos, 4*NPMAX);
	
	//beamcom_.pTest = new double[n1];
	return 0;
}

int srTSASE::Alloc_wigcom()
{
	//int n1 = 10001, n2 = 10000;
	int n1 = NZMAX + 1, n2 = NZMAX;
	if((wigcom_.awz = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(wigcom_.awz, n1);
	if((wigcom_.awdz = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(wigcom_.awdz, n2);
	if((wigcom_.solz = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(wigcom_.solz, n2);
	if((wigcom_.awerx = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(wigcom_.awerx, n2);
	if((wigcom_.awery = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(wigcom_.awery, n2);
	if((wigcom_.qfld = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(wigcom_.qfld, n1);
	if((wigcom_.dqfx = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(wigcom_.dqfx, n1);
	if((wigcom_.dqfy = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(wigcom_.dqfy, n1);
	if((wigcom_.awdx = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(wigcom_.awdx, n1);
	if((wigcom_.awdy = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(wigcom_.awdy, n1);
	return 0;
}

//int srTSASE::Alloc_diagcom(int numHarm)
int srTSASE::Alloc_diagcom()
{
	int n1 = NZMAX, n2 = NHMAX*NZMAX; //numHarm*NZMAX;

	if((diagcom_.error = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.error, n1);
	if((diagcom_.gain = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.gain, n1);
	if((diagcom_.phimid = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.phimid, n2);
	if((diagcom_.whalf = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.whalf, n1);
	if((diagcom_.logp = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.logp, n1);
	//if((diagcom_.powmid = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.powmid, n1);
	if((diagcom_.xrms = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.xrms, n1);
	if((diagcom_.yrms = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.yrms, n1);
	if((diagcom_.xpos = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.xpos, n1);
	if((diagcom_.ypos = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.ypos, n1);
	if((diagcom_.pmodhist = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.pmodhist, n2);
	if((diagcom_.gamhist = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.gamhist, n1);
	if((diagcom_.diver = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.diver, n1);
	if((diagcom_.bunphase = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.bunphase, n2);
	if((diagcom_.dgamhist = new f2c_doublereal[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.dgamhist, n1);
	if((diagcom_.ffield = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.ffield, n2);

	if((diagcom_.pgainhist = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.pgainhist, n2);
	if((diagcom_.pmidhist = new f2c_doublereal[n2]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(diagcom_.pmidhist, n2);
	return 0;
}

int srTSASE::Alloc_workspace()
{
	//int n1 = 116964, n2 = 70001;
	//long int n1 = 1907388; //, n2 = 1000001; //, n3 = 7000007;

	if((workspace_.crwork3 = new f2c_doublecomplex[4*NCMAX*NCMAX*NHMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(workspace_.crwork3, 4*NCMAX*NCMAX*NHMAX);
	if((workspace_.cpart1 = new f2c_doublecomplex[NPMAX*NHMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(workspace_.cpart1, NPMAX*NHMAX);
	if((workspace_.cpart2 = new f2c_doublecomplex[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(workspace_.cpart2, NPMAX);
	if((workspace_.cpart3 = new f2c_doublecomplex[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(workspace_.cpart3, NPMAX);

	if((workspace_.k2gg = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(workspace_.k2gg, NPMAX);
	if((workspace_.k2pp = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(workspace_.k2pp, NPMAX);
	if((workspace_.k3gg = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(workspace_.k3gg, NPMAX);
	if((workspace_.k3pp = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(workspace_.k3pp, NPMAX);
	if((workspace_.p1 = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(workspace_.p1, NPMAX);
	if((workspace_.p2 = new f2c_doublereal[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(workspace_.p2, NPMAX);

	if((workspace_.iwork = new f2c_integer[NPMAX]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(workspace_.iwork, NPMAX);
	return 0;
}

//int srTSASE::Alloc_tslipcom()
int srTSASE::Alloc_tslipcom(int nWfr, int nSlip)
{
	//int n1 = 15960700;
	//if((tslipcom_.crtime = new f2c_doublecomplex[n1]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; ZeroArr(tslipcom_.crtime, n1);
	
	//long totSize = nWfr*nWfr*(nSlip + 3);
	//long totSize = nWfr*nWfr*(nSlip + 3)*NHMAX;
	long long totSize = ((long long)nWfr)*((long long)nWfr)*(nSlip + 3)*NHMAX;
	if((tslipcom_.crtime = new f2c_doublecomplex[totSize]) == NULL) return MEMORY_ALLOCATION_FAILURE_SASE; 
	ZeroArr(tslipcom_.crtime, totSize);
	return 0;
}

//*************************************************************************

void srTSASE::Free_tbunchcom()
{
	if(tbunchcom_.tgam0 != NULL) { delete[] tbunchcom_.tgam0; tbunchcom_.tgam0 = NULL;}
	if(tbunchcom_.tdgam != NULL) { delete[] tbunchcom_.tdgam; tbunchcom_.tdgam = NULL;}
	if(tbunchcom_.temitx != NULL) { delete[] tbunchcom_.temitx; tbunchcom_.temitx = NULL;}
	if(tbunchcom_.temity != NULL) { delete[] tbunchcom_.temity; tbunchcom_.temity = NULL;}
	if(tbunchcom_.txrms != NULL) { delete[] tbunchcom_.txrms; tbunchcom_.txrms = NULL;}
	if(tbunchcom_.tyrms != NULL) { delete[] tbunchcom_.tyrms; tbunchcom_.tyrms = NULL;}
	if(tbunchcom_.txpos != NULL) { delete[] tbunchcom_.txpos; tbunchcom_.txpos = NULL;}
	if(tbunchcom_.typos != NULL) { delete[] tbunchcom_.typos; tbunchcom_.typos = NULL;}
	if(tbunchcom_.tpxpos != NULL) { delete[] tbunchcom_.tpxpos; tbunchcom_.tpxpos = NULL;}
	if(tbunchcom_.tpypos != NULL) { delete[] tbunchcom_.tpypos; tbunchcom_.tpypos = NULL;}
	if(tbunchcom_.talphx != NULL) { delete[] tbunchcom_.talphx; tbunchcom_.talphx = NULL;}
	if(tbunchcom_.talphy != NULL) { delete[] tbunchcom_.talphy; tbunchcom_.talphy = NULL;}
	if(tbunchcom_.tcurrent != NULL) { delete[] tbunchcom_.tcurrent; tbunchcom_.tcurrent = NULL;}
	if(tbunchcom_.tpos != NULL) { delete[] tbunchcom_.tpos; tbunchcom_.tpos = NULL;}
	if(tbunchcom_.tloss != NULL) { delete[] tbunchcom_.tloss; tbunchcom_.tloss = NULL;}
	if(tbunchcom_.tradpos != NULL) { delete[] tbunchcom_.tradpos; tbunchcom_.tradpos = NULL;}
	if(tbunchcom_.tzrayl != NULL) { delete[] tbunchcom_.tzrayl; tbunchcom_.tzrayl = NULL;}
	if(tbunchcom_.tzwaist != NULL) { delete[] tbunchcom_.tzwaist; tbunchcom_.tzwaist = NULL;}
	if(tbunchcom_.tprad0 != NULL) { delete[] tbunchcom_.tprad0; tbunchcom_.tprad0 = NULL;}
	if(tbunchcom_.tradphase != NULL) { delete[] tbunchcom_.tradphase; tbunchcom_.tradphase = NULL;}

	if(tbunchcom_.distgam != NULL) { delete[] tbunchcom_.distgam; tbunchcom_.distgam = NULL;}
	if(tbunchcom_.distx != NULL) { delete[] tbunchcom_.distx; tbunchcom_.distx = NULL;}
	if(tbunchcom_.disty != NULL) { delete[] tbunchcom_.disty; tbunchcom_.disty = NULL;}
	if(tbunchcom_.distpx != NULL) { delete[] tbunchcom_.distpx; tbunchcom_.distpx = NULL;}
	if(tbunchcom_.distpy != NULL) { delete[] tbunchcom_.distpy; tbunchcom_.distpy = NULL;}
	if(tbunchcom_.distt != NULL) { delete[] tbunchcom_.distt; tbunchcom_.distt = NULL;}
}

void srTSASE::Free_cartcom()
{
	//if(cartcom_.wx != NULL) delete[] cartcom_.wx; cartcom_.wx = NULL;
	//if(cartcom_.wy != NULL) delete[] cartcom_.wy; cartcom_.wy = NULL;
	//if(cartcom_.ipos != NULL) delete[] cartcom_.ipos; cartcom_.ipos = NULL;

	if(cartcom_.crfield != NULL) { delete[] cartcom_.crfield; cartcom_.crfield = NULL;}
	if(cartcom_.crsource != NULL) { delete[] cartcom_.crsource; cartcom_.crsource = NULL;}
	if(cartcom_.crhm != NULL) { delete[] cartcom_.crhm; cartcom_.crhm = NULL;}

	if(cartcom_.crmatc != NULL) { delete[] cartcom_.crmatc; cartcom_.crmatc = NULL;}
	if(cartcom_.cwet != NULL) { delete[] cartcom_.cwet; cartcom_.cwet = NULL;}
	if(cartcom_.cbet != NULL) { delete[] cartcom_.cbet; cartcom_.cbet = NULL;}
}

void srTSASE::Free_beamcom()
{
	if(beamcom_.xpart != NULL) { delete[] beamcom_.xpart; beamcom_.xpart = NULL;}
	if(beamcom_.ypart != NULL) { delete[] beamcom_.ypart; beamcom_.ypart = NULL;}
	if(beamcom_.px != NULL) { delete[] beamcom_.px; beamcom_.px = NULL;}
	if(beamcom_.py != NULL) { delete[] beamcom_.py; beamcom_.py = NULL;}
	if(beamcom_.gamma != NULL) { delete[] beamcom_.gamma; beamcom_.gamma = NULL;}
	if(beamcom_.theta != NULL) { delete[] beamcom_.theta; beamcom_.theta = NULL;}
	if(beamcom_.xporb != NULL) { delete[] beamcom_.xporb; beamcom_.xporb = NULL;}
	if(beamcom_.yporb != NULL) { delete[] beamcom_.yporb; beamcom_.yporb = NULL;}
	if(beamcom_.btpar != NULL) { delete[] beamcom_.btpar; beamcom_.btpar = NULL;}
	if(beamcom_.btper != NULL) { delete[] beamcom_.btper; beamcom_.btper = NULL;}
	if(beamcom_.ez != NULL) { delete[] beamcom_.ez; beamcom_.ez = NULL;}
	if(beamcom_.wx != NULL) { delete[] beamcom_.wx; beamcom_.wx = NULL;}
	if(beamcom_.wy != NULL) { delete[] beamcom_.wy; beamcom_.wy = NULL;}
	if(beamcom_.lostid != NULL) { delete[] beamcom_.lostid; beamcom_.lostid = NULL;}
	if(beamcom_.ipos != NULL) { delete[] beamcom_.ipos; beamcom_.ipos = NULL;}
}

void srTSASE::Free_wigcom()
{
	if(wigcom_.awz != NULL) { delete[] wigcom_.awz; wigcom_.awz = NULL;}
	if(wigcom_.awdz != NULL) { delete[] wigcom_.awdz; wigcom_.awdz = NULL;}
	if(wigcom_.solz != NULL) { delete[] wigcom_.solz; wigcom_.solz = NULL;}
	if(wigcom_.awerx != NULL) { delete[] wigcom_.awerx; wigcom_.awerx = NULL;}
	if(wigcom_.awery != NULL) { delete[] wigcom_.awery; wigcom_.awery = NULL;}
	if(wigcom_.qfld != NULL) { delete[] wigcom_.qfld; wigcom_.qfld = NULL;}
	if(wigcom_.dqfx != NULL) { delete[] wigcom_.dqfx; wigcom_.dqfx = NULL;}
	if(wigcom_.dqfy != NULL) { delete[] wigcom_.dqfy; wigcom_.dqfy = NULL;}
	if(wigcom_.awdx != NULL) { delete[] wigcom_.awdx; wigcom_.awdx = NULL;}
	if(wigcom_.awdy != NULL) { delete[] wigcom_.awdy; wigcom_.awdy = NULL;}
}

void srTSASE::Free_diagcom()
{
	if(diagcom_.error != NULL) { delete[] diagcom_.error; diagcom_.error = NULL;}
	if(diagcom_.gain != NULL) { delete[] diagcom_.gain; diagcom_.gain = NULL;}
	if(diagcom_.phimid != NULL) { delete[] diagcom_.phimid; diagcom_.phimid = NULL;}
	if(diagcom_.whalf != NULL) { delete[] diagcom_.whalf; diagcom_.whalf = NULL;}
	if(diagcom_.logp != NULL) { delete[] diagcom_.logp; diagcom_.logp = NULL;}
	//if(diagcom_.powmid != NULL) { delete[] diagcom_.powmid; diagcom_.powmid = NULL;}
	if(diagcom_.xrms != NULL) { delete[] diagcom_.xrms; diagcom_.xrms = NULL;}
	if(diagcom_.yrms != NULL) { delete[] diagcom_.yrms; diagcom_.yrms = NULL;}
	if(diagcom_.xpos != NULL) { delete[] diagcom_.xpos; diagcom_.xpos = NULL;}
	if(diagcom_.ypos != NULL) { delete[] diagcom_.ypos; diagcom_.ypos = NULL;}
	if(diagcom_.pmodhist != NULL) { delete[] diagcom_.pmodhist; diagcom_.pmodhist = NULL;}
	if(diagcom_.gamhist != NULL) { delete[] diagcom_.gamhist; diagcom_.gamhist = NULL;}
	if(diagcom_.diver != NULL) { delete[] diagcom_.diver; diagcom_.diver = NULL;}
	if(diagcom_.bunphase != NULL) { delete[] diagcom_.bunphase; diagcom_.bunphase = NULL;}
	if(diagcom_.dgamhist != NULL) { delete[] diagcom_.dgamhist; diagcom_.dgamhist = NULL;}
	if(diagcom_.ffield != NULL) { delete[] diagcom_.ffield; diagcom_.ffield = NULL;}

	if(diagcom_.pgainhist != NULL) { delete[] diagcom_.pgainhist; diagcom_.pgainhist = NULL;}
	if(diagcom_.pmidhist != NULL) { delete[] diagcom_.pmidhist; diagcom_.pmidhist = NULL;}
}

void srTSASE::Free_workspace()
{
	if(workspace_.crwork3 != NULL) { delete[] workspace_.crwork3; workspace_.crwork3 = NULL;}
	if(workspace_.cpart1 != NULL) { delete[] workspace_.cpart1; workspace_.cpart1 = NULL;}
	if(workspace_.cpart2 != NULL) { delete[] workspace_.cpart2; workspace_.cpart2 = NULL;}
	if(workspace_.cpart3 != NULL) { delete[] workspace_.cpart3; workspace_.cpart3 = NULL;}

	if(workspace_.k2gg != NULL) { delete[] workspace_.k2gg; workspace_.k2gg = NULL;}
	if(workspace_.k2pp != NULL) { delete[] workspace_.k2pp; workspace_.k2pp = NULL;}
	if(workspace_.k3gg != NULL) { delete[] workspace_.k3gg; workspace_.k3gg = NULL;}
	if(workspace_.k3pp != NULL) { delete[] workspace_.k3pp; workspace_.k3pp = NULL;}
	if(workspace_.p1 != NULL) { delete[] workspace_.p1; workspace_.p1 = NULL;}
	if(workspace_.p2 != NULL) { delete[] workspace_.p2; workspace_.p2 = NULL;}

	if(workspace_.iwork != NULL) { delete[] workspace_.iwork; workspace_.iwork = NULL;}
}

void srTSASE::Free_tslipcom()
{
	if(tslipcom_.crtime != NULL) { delete[] tslipcom_.crtime; tslipcom_.crtime = NULL;}
}

//*************************************************************************

int srTSASE::InitMainGenesisStructs()
{
	int result = 0;
	if(Alloc_tbunchcom() != 0) return MEMORY_ALLOCATION_FAILURE_SASE;
	if(Alloc_cartcom() != 0) return MEMORY_ALLOCATION_FAILURE_SASE;
	//if(Alloc_cartcom(numHarm) != 0) return MEMORY_ALLOCATION_FAILURE_SASE;
	if(Alloc_beamcom() != 0) return MEMORY_ALLOCATION_FAILURE_SASE;
	if(Alloc_wigcom() != 0) return MEMORY_ALLOCATION_FAILURE_SASE;
	if(Alloc_diagcom() != 0) return MEMORY_ALLOCATION_FAILURE_SASE;
	if(Alloc_workspace() != 0) return MEMORY_ALLOCATION_FAILURE_SASE;
	//if(Alloc_tslipcom() != 0) return MEMORY_ALLOCATION_FAILURE_SASE; //to be allocated separately

	return result;
}

//*************************************************************************

void srTSASE::ReleaseGenesisStructs()
{
    Free_tbunchcom();
	Free_cartcom();
	Free_beamcom();
	Free_wigcom();
	Free_diagcom();
	Free_workspace();
	Free_tslipcom();
}

//*************************************************************************

int srTSASE::CheckInputConsistency()
{// re-written chk_bnd__()
	//int result = 0;

	//if(iocom_1.npin <= 0) 
	//{
	//	inputcom_1.convharm = 1;
	//}
	if((!m_ElecDistribShouldBeUsed) && (inputcom_1.convharm > 1))
	{
		inputcom_1.convharm = 1;
		CErrWarn::AddWarningMessage(&gVectWarnNos, GENESIS_RAD_HARM_CALC_NEEDS_ELEC_DISTRIB);
	}

    if(inputcom_1.nslice <= 0) 
    {
		if (inputcom_1.curlen < 0.) 
		{
			/* step profile */
			inputcom_1.ntail = 0;
			inputcom_1.nslice = (int)(fabs(inputcom_1.curlen) / inputcom_1.xlamds / inputcom_1.zsep);
		} 
		else 
		{
			/* gaussian */
			inputcom_1.ntail = -((int)(inputcom_1.curlen * 3. / inputcom_1.xlamds / inputcom_1.zsep));
			inputcom_1.nslice = (int)(inputcom_1.curlen * 6. / inputcom_1.xlamds / inputcom_1.zsep);
		}
    }

	if(inputcom_1.aw0 <= 0.) return GENESIS_NO_WIGGLER_FIELD;
	if((inputcom_1.iwityp > 1) || (inputcom_1.iwityp < 0)) return GENESIS_INVALID_WIGGLER;
	if(inputcom_1.xlamd <= 0) return GENESIS_WIGGLER_PERIOD_MUST_BE_POSITIVE;
	//if((inputcom_1.gamma0 - 4*inputcom_1.delgam) < 1) return GENESIS_ELECTRON_BEAM_ENERGY_TOO_SMALL;
	if(((inputcom_1.gamma0 - 4*inputcom_1.delgam) < 1) && (EbmDat.pElecDistr == 0)) return GENESIS_ELECTRON_BEAM_ENERGY_TOO_SMALL;

	if(((inputcom_1.npart >> 3) << 3) != inputcom_1.npart) return GENESIS_NUMBER_OF_PARTICLES_IS_NOT_MULTIPLE_OF_8;
	if(inputcom_1.xlamds <= 0.) return GENESIS_RADIATION_WAVELENGTH_MUST_BE_POSITIVE;

	//if(inputcom_1.prad0 <= 0.) return GENESIS_RADIATION_POWER_MUST_BE_POSITIVE;
	if((inputcom_1.prad0 <= 0.) && (!SeedRad.ElecFieldIsDefined())) return GENESIS_RADIATION_POWER_MUST_BE_POSITIVE;

	if((inputcom_1.delaw != 0.) && (::fabs((double)(inputcom_1.delz - 0.5)) > 1.E-9) && (::fabs((double)inputcom_1.iertyp) != 0)) return GENESIS_DELZ_CONSTRAINT;
	//if(inputcom_1.iscan > 8) return GENESIS_ISCAN_CONSTRAINT;
	if(inputcom_1.iscan > 25) return GENESIS_ISCAN_CONSTRAINT;

	if(inputcom_1.iscan > 22) 
	{
		//if (tbunchcom_1.ndata <= 0) {
		//    i__ = printerr_(&c_n8, "BEAMFILE for scan not defined", (ftnlen)
		//	    29);
		//    last_();
		//}
		inputcom_1.nslice = tbunchcom_1.ndata;
		inputcom_1.nscan = tbunchcom_1.ndata;
	}

	if(inputcom_1.nwig <= 0) return GENESIS_NWIG_MUST_BE_POSITIVE;
	if(inputcom_1.delz <= 0) return 1; //"DELZ must be positive and non-zero" //GENESIS_DELZ_MUST_BE_POSITIVE;
	if(inputcom_1.zsep <= 0) return 1; //"ZSEP must be al least 1";
	
	//if(inputcom_1.npart > 70001) return 1; //"NPART > NPMAX - setting NPART=NPMAX";
	if(inputcom_1.npart > NPMAX) 
	{
		CErrWarn::AddWarningMessage(&gVectWarnNos, GENESIS_TOO_MANY_PARTICLES_REQUESTED);
		inputcom_1.npart = NPMAX;
	}
	
	if(inputcom_1.nbins < 4) inputcom_1.nbins = 4;
	
	if(inputcom_1.npart % (inputcom_1.nbins << 2) != 0) return 1; //"NPART not a multiple of 4*NBINS"
	//to programm returning error number

    //for(int i1 = inputcom_1.nharm + 1; i1 <= 5; ++i1) 
  //TEST  for(int i1 = inputcom_1.nharm + 1; i1 <= NHMAX; ++i1) 
  //  {
		//if(inputcom_1.lout[i1 + 13] != 0) 
		//{
		//	//CErrWarn::AddWarningMessage(&gVectWarnNos, GENESIS_TOO_MANY_PARTICLES_REQUESTED);
		//	//i__ = printerr_(&c_n9, "harmonic output not calculated", (ftnlen)30);
		//}
  //  }
    
	if(inputcom_1.idmppar > inputcom_1.nharm) 
	{
		//i__ = printerr_(&c_n9, "No dump possible (IDMPPAR > IHARM)", (ftnlen)34);
		inputcom_1.idmppar = 0;
    }

	//int IBAS[] = {inputcom_1.ildpsi, inputcom_1.ildx, inputcom_1.ildy, inputcom_1.ildpx, inputcom_1.ildpy, inputcom_1.ildgam, inputcom_1.ildgam + 1};
	//for(int i1=0; i1<7; i1++)
	//TEST for(int i1=0; i1<NHMAX; i1++)
	//{
	//	//for(int i2=i1+1; i2<7; i2++)
	//	for(int i2=i1+1; i2<NHMAX; i2++)
	//	{
	//		if(IBAS[i1] == IBAS[i2]) CErrWarn::AddWarningMessage(&gVectWarnNos, GENESIS_SAME_BASE_FOR_DIFFERENT_LOADING);
	//	}
	//}

	if((inputcom_1.iscan > 0) && (inputcom_1.nscan <= 1 || fabs(inputcom_1.svar) < 1.e-7)) return GENESIS_PARAMETER_SCAN_NOT_DEFINED;
	
	//if(inputcom_1.itdp != 0)
	//{//time-dependent computation
		//double test_zsep = inputcom_1.zsep, test_delz = inputcom_1.delz;
		//double testVal = inputcom_1.zsep / inputcom_1.delz - (f2c_integer) (inputcom_1.zsep / inputcom_1.delz) * (float)1.;
	if(inputcom_1.zsep / inputcom_1.delz - (f2c_integer) (inputcom_1.zsep / inputcom_1.delz) * (float)1. > (float)1e-8) return GENESIS_SLICE_SEPAR_NOT_INTEGER_OF_STEPSIZE;
	//}
	
	if(inputcom_1.nslice > 50000) return GENESIS_TOO_MANY_SLICES;
	if(inputcom_1.nslice <= 0) return 1; //"NSLICE < 1";
	
	if((inputcom_1.zrayl <= (float)0.) && (!SeedRad.ElecFieldIsDefined())) return 1; //"ZRAYL must be larger than 0"

	if(inputcom_1.ncar % 2 == 0) return GENESIS_GRID_POINTS_OF_MESH_NOT_ODD;
	if(inputcom_1.ncar > NCMAX) 
	{
		CErrWarn::AddWarningMessage(&gVectWarnNos, GENESIS_TOO_MANY_GRID_POINTS);
		inputcom_1.ncar = NCMAX;
	}

	if(inputcom_1.nptr > NRGRID - 1) 
	{
		CErrWarn::AddWarningMessage(&gVectWarnNos, GENESIS_TOO_MANY_RADIAL_GRID_POINTS);
		inputcom_1.nptr = NRGRID - 1;
	}
	if(inputcom_1.nptr < 2) 
	{
		CErrWarn::AddWarningMessage(&gVectWarnNos, GENESIS_NOT_ENOUGH_RADIAL_GRID_POINTS);
		//inputcom_1.nptr = NRGRID - 1;
		inputcom_1.nscz = 0;
		inputcom_1.nscr = 0;
	}
	
	if(inputcom_1.nscz >= inputcom_1.nbins / 2 + 1) 
	{
		/* somehow empirical boundary */
		//i__ = printerr_(&c_n9, "NSCZ too large - setting NSCZ=2", (ftnlen)31);
		inputcom_1.nscz = 2;
    }
    if(inputcom_1.nharm > inputcom_1.nbins / 2 + 1) 
    {
		//CErrWarn::AddWarningMessage(&gVectWarnNos, GENESIS_NOT_ENOUGH_RADIAL_GRID_POINTS);
		//i__ = printerr_(&c_n9, "Higher harmonics are inaccurate (NHARM)", (ftnlen)39);
    }

    if(CheckIfElecDistrShouldBeUsed() && (inputcom_1.nslice*inputcom_1.npart > EbmDat.nTotMacroPart)) 
    {
		return GENESIS_INCONSISTENT_NSLICE_NPART;
	}

	//if((timecom_1.itdp != 0) && (scancom_1.iscan > 0)) return GENESIS_TIME_DEPENDENCY_AND_SCAN_ARE_SELECTED_SIMULT;
	//if((timecom_1.itdp != 0) && (simcom_1.delz != 0.) && ((f2c_integer)(tbunchcom_1.zsep / simcom_1.delz) <= 0)) return GENESIS_STEP_SIZE_GREATER_THAN_SLICE_SEPARATION;

	//if(wigcom_1.fl / simcom_1.delz - (f2c_integer) (wigcom_1.fl / simcom_1.delz) * (float)1. > (float)1e-8) return GENESIS_FOC_LENGTH_NOT_INTEGER_OF_STEPSIZE;
	//if(wigcom_1.dl / simcom_1.delz - (f2c_integer) (wigcom_1.dl / simcom_1.delz) * (float)1. > (float)1e-8) return GENESIS_DEFOC_LENGTH_NOT_INTEGER_OF_STEPSIZE;
	//if(wigcom_1.drl / simcom_1.delz - (f2c_integer) (wigcom_1.drl / simcom_1.delz) * (float)1. > (float)1e-8) return GENESIS_DRIFT_LENGTH_NOT_INTEGER_OF_STEPSIZE;
	//if(wigcom_1.f1st / simcom_1.delz - (f2c_integer) (wigcom_1.f1st / simcom_1.delz) * (float)1. > (float)1e-8) return GENESIS_START_OF_FODO_NOT_INTEGER_OF_STEPSIZE;
	//
	//if((f2c_real) (simcom_1.nwig * simcom_1.nsec) / simcom_1.delz > 2501.) return GENESIS_TOO_MANY_INTEGRATION_STEPS;
	//
	//if(radiatcom_1.nscz > 3) 
	//{
	//	//pSend->AddWarningMessage(&gVectWarnNos, GENESIS_TOO_MANY_LONGITUDINAL_MODES);
	//	CErrWarn::AddWarningMessage(&gVectWarnNos, GENESIS_TOO_MANY_LONGITUDINAL_MODES);
	//	radiatcom_1.nscz = 3;
	//}

	return 0;
}

//*************************************************************************

int srTSASE::ConvertInputDataToGenesisFormat(int numHarm) //OC191108
{
	DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat = PrecDat.AllowAutoChoiceOfNxNzForPropagat;
	DistrInfoDat.NxNzOversamplingParam = PrecDat.NxNzOversamplingParam;

	CGenMathFFT2D FFT;
	FFT.NextCorrectNumberForFFT(PrecDat.ncar);
	PrecDat.ncar++; // since Genesis requires odd ncar, and in SRW it is even

	srTWigComSASE WigCom;
	MagDat.SetupWigSASE(WigCom);

    iocom_1.nprobe = 30; // filenumber for filetype probin 
    iocom_1.nlog = 6;

			//test !!!
			//PrecDat.ncar = 151;
			//end test
	inputcom_1.ncar = PrecDat.ncar; //"NCAR"
	inputcom_1.nsec = WigCom.simcom_nsec; //"NSEC"
	inputcom_1.iorb = PrecDat.iorb; // "IORB" - flag for orbit correction
	inputcom_1.delz = PrecDat.delz; // "DELZ"
	inputcom_1.qfdx = WigCom.qfdx; // "QFDX" - min. transverse misplacement of quads in x direction
	inputcom_1.qfdy = WigCom.qfdy; // "QFDY"
	inputcom_1.nwig = WigCom.simcom_nwig; // "NWIG"
	inputcom_1.nscr = PrecDat.nscr; // "NSCR"
	inputcom_1.itdp = PrecDat.itdp; // "ITDP" - A non-zero value enables time-dependent simulation. Time-dependence is not allowed if the scan-feature is enabled.

	inputcom_1.nscz = PrecDat.nscz; // "NSCZ"
	inputcom_1.zsep = PrecDat.zsep; // "ZSEP"
	inputcom_1.nptr = PrecDat.nptr; // "NPTR"
	//inputcom_1.lout // "LOUT" - defines which parameters to include into input file, to modify ?
	
	inputcom_1.convharm = EstimHarmNumFromPhotEn(WigCom); //1; // "CONVHARM" - harmonic number - to modify for HGHG !!!???
	inputcom_1.multconv = 0; // "MULTCONV" - affects accuracy of computation at higher harmonics
	//If an imported particle distribution from a PARTFILE is up-converted to a higher harmonics
	//the dault behavior is that the number of slices is preserved. This requires that ZSEP is
	//adjusted together with XLAMDS. However if frequency resolution is a problem then a particle
	//distribution can be converted and used multiple times to keep ZSEP constant. The disadvantage is that the CPU execution time is increased as well.

	inputcom_1.ibfield = WigCom.chic_bfield; //dispersion section parameters for HGHG // "IBFIELD" - When the PARTFILE features is used the imported particle distribution can be tracked through a generic 4 magnet chicane before running the Genesis simulation. The chicane consists out of 4 bending magnets of the field strength IBFIELD and length IMAGL separated by 5 drifts of length IDRIL. If the field strength of the magnet is set to zero the feature of a chicane is disabled (default behavior).
	inputcom_1.idril = WigCom.chic_dril; // "IDRIL" - the length of the 5 drift lengths of the magnetic chicane (three between the magnets and one before and after the magnets) 
	inputcom_1.imagl = WigCom.chic_magl; // "IMAGL" - length of bending magnet of chicane in [m]

	//inputcom_1.filetype // "FILETYPE"
	inputcom_1.prad0 = InRad.Power; // "PRAD0"
	inputcom_1.rmax0 = PrecDat.rmax0; // "RMAX0"

    inputcom_1.delaw = WigCom.delaw; // "DELAW"
	inputcom_1.xbeam = EbmDat.x0; // "XBEAM" check long. pos.
	inputcom_1.ybeam = EbmDat.z0; // "YBEAM" check long. pos.
	inputcom_1.bunch = 0; // "BUNCH" initial value for bunching factor
	inputcom_1.quadf = WigCom.quadf; // "QUADF"
	inputcom_1.quadd = WigCom.quadd; // "QUADD"
	
	inputcom_1.xlamd = WigCom.xlamd; // "XLAMD"

	inputcom_1.emitx = EbmDat.NormalizedEmittanceX(); // "EMITX" check this function !!!
	inputcom_1.emity = EbmDat.NormalizedEmittanceZ(); // "EMITY" check this function !!!
	inputcom_1.npart = PrecDat.npart; // "NPART"
	inputcom_1.ntail = PrecDat.ntail; // "NTAIL" - (-253 – integer – unitless) Position of the first simulated slice in measures of ZSEP*XLAMDS. GENESIS 1.3 starts with the tail side of the time window, progressing towards the head. Thus a negative or positive value shifts the slices towards the tail or head region of the beam, respectively. For a constant profile (CURLEN < 0) NTAIL has no impact.

	inputcom_1.gamma0 = EbmDat.Gamma; // "GAMMA0"
	//double W0 = InRad.WaistDiam; // or 0.5* ?
	//inputcom_1.zrayl = PI*W0*W0/inputcom_1.xlamds; // "ZRAYL"

	inputcom_1.zstop = PrecDat.zstop; // "ZSTOP"
	inputcom_1.fbess0 = WigCom.fbess0; // "FBESS0"
	inputcom_1.shotnoise = EbmDat.ShotNoiseFactor; // "SHOTNOISE" (1.0 – float – unitless) GENESIS 1.3 applies a random o_set to each macro particle phase to generate the correct statistic for the bunching factor. Each o_set is scaled prior by SHOTNOISE, thus SHOTNOISE can be set to zero to disable shot noise.
	inputcom_1.dl = PrecDat.delz*RoundDoubleToInt(WigCom.dl/PrecDat.delz); // "DL" re-calc to integer of step size
	inputcom_1.fl = PrecDat.delz*RoundDoubleToInt(WigCom.fl/PrecDat.delz); // "FL"
	inputcom_1.delgam = EbmDat.SigmaRelE*EbmDat.Gamma; // "DELGAM"
	
	//inputcom_1.pxbeam = EbmDat.dxds0; // "PXBEAM"
	//inputcom_1.pybeam = EbmDat.dzds0; // "PYBEAM" check long. pos.
	inputcom_1.pxbeam = EbmDat.Gamma*EbmDat.dxds0; // "PXBEAM" //OC051108 corrected to *gamma
	inputcom_1.pybeam = EbmDat.Gamma*EbmDat.dzds0; // "PYBEAM" check long. pos.

	inputcom_1.alphax = EbmDat.AlphaX(); // "ALPHAX" check this !!! 
	inputcom_1.alphay = EbmDat.AlphaZ(); // "ALPHAY" check this !!!
	inputcom_1.rxbeam = sqrt(::fabs(EbmDat.Mxx)); // "RXBEAM"
	inputcom_1.rybeam = sqrt(::fabs(EbmDat.Mzz)); // "RYBEAM"

	//inputcom_1.xlamds = 1.239842E-06/DistrInfoDat.LambStart; // "XLAMDS" - wavelength in m assuming input in eV ! The resonant radiation wavelength. Due to the bandwidth of time-dependent simulation SASE simulations do not require a very precise value for XLAMDS
	inputcom_1.xlamds = 1.239842E-06/PrecDat.photEn_xlamds; // "XLAMDS" - wavelength in m assuming input in eV ! The resonant radiation wavelength. Due to the bandwidth of time-dependent simulation SASE simulations do not require a very precise value for XLAMDS
	
	double W0 = InRad.WaistDiam; // or 0.5* ?
	inputcom_1.zrayl = PI*W0*W0/inputcom_1.xlamds; // "ZRAYL"
	
	// "CURLEN" - Bunch length of the current profile:
	if(EbmDat.TypeDistrLongitudinal == 2) inputcom_1.curlen = sqrt(::fabs(EbmDat.Mss)); // gaussian
	else inputcom_1.curlen = -1; //uniform current
	inputcom_1.nslice = PrecDat.nslice; // "NSLICE" - Total number of simulated slices. It defines the time window of the simulation with NSLICE * ZSEP * XLAMDS/c. Note that the output does not start with the first slice unless the parameter IOTAIL is set. If NSLICE set to zero it automatically adjust NSLICE and NTAIL to the time-window given by the external input files (BEAMFILE or DISTFILE). It assumes 6 standard deviation for a Gaussian distribution or the absolute value of CURLEN for a step profile.
	inputcom_1.itgaus = 3 - EbmDat.TypeDistrTransverse; // "ITGAUS" - Defines distribution profile in the transverse variables. The available distributions are: - Gaussian (ITGAUS=1), - Uniform (ITGAUS=2), - Parabolic (otherwise)
	
	inputcom_1.wcoefz[0] = 0; // "WCOEFZ" [0]- Start of undulator tapering [m]. Note that tapering is applied, even the magnetic lattice is defined by an external file.
	inputcom_1.wcoefz[1] = 0; // [1] - The relative change of the undulator field over the entire taper length (AW(exit) = (1 -WCOEFZ(2)) AW(entrance)). In the case of a multi section undulator GENESIS 1.3 tapers the magnetic field over the gaps as well, resulting in a jump of the magnetic field AW(z) between two modules.
	inputcom_1.wcoefz[2] = 0; // [2] - The taper model: = 1 for linear taper, = 2 for quadratic taper, or no taper otherwise.
	inputcom_1.sl = 0; // "SL" - solenoid length
	inputcom_1.solen = 0; // "SOLEN" - strength of solenoid field
	
	inputcom_1.dgrid = 0; // "DGRID" - grid size(-dgrid to dgrid) if dgrid > 0 - overrides scaling from rmax0
	
	inputcom_1.igamgaus = 1; // "IGAMGAUS" mesh discretization: gaussian (<>0) or uniform (=0) enegy load
	inputcom_1.iall = 0; //"IALL" - A non-zero value of IALL enforces that all slices are starting with the same element of the Hammersley sequences. It is recommend for time-dependent simulation excluding shot noise. IALL is set automatically for scans

	inputcom_1.nbins = 16; //12; //4; // "NBINS" number of bins in the particle phase, should be > 2*n + 2, where n is harmonic number
	
	//to modify!!!
	inputcom_1.iallharm = 1; //!<>0 => all higher harmonics are calculated up to nharm 
	inputcom_1.nharm = numHarm; //7; //5; // "NHARM" - number of harmonics in the bunching factor to be calculated for output (max. 7)
	
	inputcom_1.emodphase = 0; // "EMODPHASE" - radiation phase for energy modulation
	inputcom_1.emod = 0; //"EMOD" - initial energy modulation
	inputcom_1.ildpx = 3; // "ILDPX" - Hammerseley base for x momentum loading parameter
	inputcom_1.ildpy = 4; // "ILDPY" - Hammerseley base for x momentum loading parameter
	inputcom_1.ildx = 1; // "ILDX" - Hammersley base for loading distribution in x
	inputcom_1.ildy = 2; // "ILDY" - Hammersley base for loading distribution in y
	inputcom_1.ipseed = -1; // "IPSEED" - initial seed for the random number generator used for particle phase fluctuation (shot noise). GENESIS 1.3 requires a re-initialization of the random number generator to guarantee the same loading whether magnetic field errors are used or not.
	inputcom_1.iseed = -1; // "ISEED" - initial seed for wiggler error generation

	inputcom_1.ildgam = 5; // "ILDGAM" - Hammerseley base for energy distribution
	inputcom_1.bunchphase = 0.; // "BUNCHPHASE" - phase of initial bunching, when quite loading is used
	inputcom_1.isravg = 0; // "ISRAVG" - if set to a non-zero value the energy loss due to spontaneous synchrotron radiation is included in the calculation.
	inputcom_1.isrsig = 0; // "ISRSIG" - if set to a non-zero value the increase of the energy spread due to the quantum fluctuation of the spontaneous synchrotron radiation is included in the calculation.

	// "OFFSETRADF" - If the automatic alignment of the radiaiton field is disabled by setting ALIGNRADF 
	//to a nonzero value, the default alignment is that the first slice of the radiaiton field overlaps 
	//with the first slice of the electron beam. However the relative position can be controlled by OFFSETRADF. 
	//The value of OFFSETRADF defines the number of slice to skip before filling it for the first electorn 
	//slice. E.g. a value of 4 will uses the 5th slice for the simulation of the first slice. slices one 
	//to 4 will be used to fill up the slippage field. If Genesis 1.3. require to fill a slice which is not 
	//defined by the FIELDFILE then it uses the internal method of a fundamental Gauss-Hermite mode.
	inputcom_1.alignradf = PrecDat.alignradf; // "ALIGNRADF" <>0 imported radfile is aligned to electron beam
	inputcom_1.offsetradf = PrecDat.offsetradf; 	

	inputcom_1.cuttail = -1; // "CUTTAIL" - Cut in the transverse phase space in measures of the rms size to collimate transverse beam tails/halos. The cut is applied after the loading and beam current is set accordingly to the reduced number of macro particles. It is disabled if the value is negative or the electron beam is imported from an external file.
	inputcom_1.conditx = 0; // "CONDITX" - [1/m] correlation strength between the amplitude of the eletron’s betatron oscillation and its energy. If the condition is applied correctly any emittance e_ects are removed from the FEL amplification and the focal strength can be increased. However it requires a conditioning beamline prior to the FEL to apply the correlation. Refer to the paper of Sessler (A.N. Sessler, et al., Phys. Rev. Lett 68 (1992) 309) for more information.
	inputcom_1.condity = 0; // "CONDITY" - [1/m] same as CONDITX but for the correlation in the y-plane
	inputcom_1.eloss = 0; // "ELOSS" - energy loss per meter - Externally applied energy loss of the electron beam.
	inputcom_1.ndcut = -1; // "NDCUT" =<0 self optimized binning of ext. dist.
	//When loading a slice, only particles of the external distribution are used, which falls within a small time-window centered around the current position of the slice. If NDCUT has a value larger than zero the width is calculated by (tmax-tmin)/NDCUT, where tmax and tmin are determined, while scanning through the particle distribution. If NDCUT is zero, the timewindow is adjusted, so that in average NPART/NBINS particles fall in each slice.
	inputcom_1.awx = 0; // "AWX" - max offset in x for undulator misalignment
	inputcom_1.awy = 0; // "AWY" - max offset in y for undulator misalignment
	inputcom_1.isntyp = 0; // "ISNTYP" - For VERSION below 1.0 only the Pennman algorighm is used for the shot noise, which is only correct for the fundamental wavelength, but for a higher version the default is the Fawley algorithm which applies also the correct shotnoise to all higher harmonics. If the user wants to use the Pennman algorith at a higher version number (which is not recommended), the value of ISNTYP has to be set to a non-zero value.
	inputcom_1.ffspec = 0; // "FFSPEC" - to calculate spectrum at post-processing; <0 => on-axis intensity in far field, =0 => on-axis intensity in near field, >0 => total power in near field
	//To calculate a spectrum a post-processing step requires amplitude and phase, which are writen to the output file, defined by LOUT of column 3 and 4. The values depend on the
	//choice of FFSPEC. If the value is equal the near field on-axis intensity and phase is written, while for a negative value it is the same but in the far field zone. For a positive value the
	//near field intensity is replaced by the total radiation power, assuming transverse coherence.

	inputcom_1.ildpsi = 7; // "ILDPSI" - index of the Hammersley sequence bases for loading the particle phase
	inputcom_1.iotail = 1; //0; // "IOTAIL" - If set to a non-zero value the output time window is the same as the simulated time window. Otherwise the output for the first slices covered by the slippage length is suppressed. Needed for bunches which are completely covered by the time-window
	inputcom_1.iscan = 0; // "ISCAN" >0 -> iscanth parameter selected
	//(0 – integer – unitless) Selects the parameter for a scan over a certain range of its value:
		//1. GAMMA0
		//2. DELGAM
		//3. CURPEAK
		//4. XLAMDS
		//5. AW0
		//6. ISEED
		//7. PXBEAM
		//8. PYBEAM
		//9. XBEAM
		//10. YBEAM
		//11. RXBEAM
		//12. RYBEAM
		//13. XLAMD
		//14. DELAW
		//15. ALPHAX
		//16. ALPHAY
		//17. EMITX
		//18. EMITY
		//19. PRAD0
		//20. ZRAYL
		//21. ZWAIST
		//22. AWD
		//23. BEAMFILE
		//24. BEAMOPT
		//25. BEAMGAM
	//The last three enable a user-defined scan, where a steady-state run is performed for each entry of the supplied BEAMFILE. BEAMOPT di_ers from BEAMFILE that it automatically
	//adjust XLAMDS according to the resonance condition and the beam energy, defined in the BEAMFILE. The BEAMGAM option is similar to BEAMFILE, but overwrites the value of energy with GAMMA0 from the main input file.
	inputcom_1.scan[0] = '\0'; //"SCAN" - By supplying the parameter name to scan over it overrules the setting of ISCAN.
	inputcom_1.nscan = 3; // "NSCAN" number of steps per scan (is taken into account only if ISCAN > 0)
	inputcom_1.svar = 0.01; // "SVAR" (0.01 – float – unitless) - Defines the scan range of the selected scan parameter. The parameter is varied between (1-SVAR) and (1+SVAR) of its initial value. One exception is the scan over ISEED where the random number generator is not reinitialized.
	inputcom_1.version = 0.1; // "VERSION" - Used for backward compatibility of the input decks. Some parameters might change their behavior for older versions of GENESIS 1.3. The current version is 1.0.

	inputcom_1.fieldfile[0] = '\0'; // "FIELDFILE"
	inputcom_1.beamfile[0] = 0; //"BEAMFILE"
	inputcom_1.maginfile[0] = '\0'; // "MAGINFILE"
	inputcom_1.magoutfile[0] = '\0'; // "MAGOUTFILE" - Defines the file to which the magnetic field lattice is written to, bypassing the interactive request of MAGOUT
	inputcom_1.outputfile[0] = '\0'; // "OUTPUTFILE" - ??? The name of the main output file. The prompt for the output filename at the beginning of a GENESIS 1.3 run will be suppressed.
	inputcom_1.radfile[0] = '\0'; // "RADFILE" - Specifying a file containing a lookup table for the seeding radiation pulse at different position along the bunch.
	//BESSY simulations use such file
	inputcom_1.partfile[0] = '\0'; // "PARTFILE"
	inputcom_1.distfile[0] = '\0'; // "DISTFILE"

	inputcom_1.ispart = 0; // "ISPART" - write the particle distribution to file for every ISPART slice.
	inputcom_1.ippart = 0; // "IPPART" - write the particle distribution to file at each IPPARTth integration step. To disable output, set IPPART to zero. The filename is the same of the main outputfile + the extension ’.par’.
	inputcom_1.ipradi = 0; // "IPRADI" - write the radiation field into file at each step
	inputcom_1.iphsty = 1; // "IPHSTY" - Generate output in the main output file at each IPHSTYth integration step. To disable output set IPHSTY to zero.
	inputcom_1.ishsty = 1; // "ISHSTY" - Generate output in the main output file for each ISHSTYth slice.
    inputcom_1.magin = 0; // "MAGIN" - (0 – integer – unitless) - If set to a non-zero value the user is prompted to type in the file name containing a explicit description of the magnetic field.
	inputcom_1.magout = 0; // "MAGOUT" - Similar to MAGIN to write out the magnetic field lattice used for the simulation.
	inputcom_1.idump = 0; // "IDUMP" - If set to a non-zero value the complete particle and field distribution is dumped at the undulator exit into two outputfiles. The filenames are the filename of the main output file plus the extension ’.dpa’ and ’.dfl’, respectively.
	inputcom_1.idmpfld = 0; // "IDMPFLD" - similar to IDUMP but only for the field distribution
	inputcom_1.idmppar = 0; // "IDMPPAR" - similar to IDUMP but only for the particle distribution
	inputcom_1.isradi = 0; // "ISRADI" - write the field distribution to file for every ISRADI slice, if !=0
	inputcom_1.ilog = 0; // "ILOG" - log file (not used)

	inputcom_1.aw0 = WigCom.aw0; // "AW0" - The normalized, dimensionless rms undulator parameter, defined by AW0 = (e/mc)(Bu/ku), where e is the electron charge, m is electron mass, c is speed of light, ku=2_/_u is the undulator wave number, _u is the undulator period. Bu is the rms undulator field with Bu = Bp/2 for a planar undulator and Bu = Bp for a helical undulator, where Bp is the on-axis peak field.
	inputcom_1.iertyp = WigCom.iertyp; // "IERTYP" - Type of undulator field errors. Either a uniform (+/-1) or Gaussian (+/- 2) distribution can be chosen. If IERTYP is negative the errors are correlated to minimize the first and second field integral. IERTYP =0 disables field errors. Field errors requires a integration step size of half an undulator period (DELZ = 0.5). Note that field errors are applied even if the magnetic lattice is defined by an external file.
	
	//inputcom_1.zwaist = InRad.WaistLongPos; // "ZWAIST"
    double TotStructLength = WigCom.xlamd*WigCom.simcom_nsec*WigCom.simcom_nwig + (WigCom.simcom_nsec - 1)*WigCom.GapLen;
    inputcom_1.zwaist = InRad.WaistLongPos + TotStructLength; 
    //assuming 0 position at the end of the und. structure in the input

	inputcom_1.iwityp = WigCom.iwityp; // "IWITYP" - Flag indicating the undulator type. A value of zero indicates a planar undulator, any other value a helical one.
	inputcom_1.lbc = PrecDat.lbc; // "LBC" - Flag to set the boundary condition for the field solver. The options are Direchlet boundary condition (LBC = 0) and Neumann boundary condition (otherwise). Anyhow the grid should be chosen large enough so that the choice of LBC should not change the results.
	inputcom_1.awd = WigCom.awd; // "AWD" - A virtual undulator parameter for the gap between undulator modules. The only purpose of this parameter is to delay the longitudinal motion of the electrons in the same manner as AW0 does within the undulator modules. It is used to keep the electron and radiation phases synchronize up to the point when the interaction at the next undulator module starts again. AWD has typically the same value as AW0, but might vary for optimum matching between the modules.
	
	inputcom_1.drl = PrecDat.delz*RoundDoubleToInt(WigCom.drl/PrecDat.delz); // "DRL"
	
	inputcom_1.curpeak = EbmDat.Current; // "CURPEAK" - Peak current of the electron beam. Time-independent simulations enforces a constant current.
	if(EbmDat.CurrentPeak > 0) inputcom_1.curpeak = EbmDat.CurrentPeak;
	
	inputcom_1.xkx = WigCom.xkx; // "XKX" - Normalized natural focusing of the undulator in x. Common values are XKX = 0.0, XKY = 1.0 for a planar undulator or XKX, XKY = 0.5 for a helical undulator, but might vary if focusing by curved pole faces is simulated. The values should fulfill the constraint XKX + XKY = 1.0.
	inputcom_1.xky = WigCom.xky; // "XKY" - Normalized natural focusing of the undulator in y. 
		
	inputcom_1.f1st = PrecDat.delz*RoundDoubleToInt(WigCom.f1st/PrecDat.delz); // "F1ST" - Position within a FODO cell, where GENESIS 1.3 starts the FODO cell lattice. To start with a full F-quadrupole set F1ST to zero while a value of FL/2 begins the cell in the middle of the focusing quadrupole.
		
		//debug
		//inputcom_1.f1st = inputcom_1.drl; //0; //inputcom_1.drl; //test
		//end debug
		
		//inputcom_1.fl = inputcom_1.drl/3.;
		//inputcom_1.dl = inputcom_1.drl/3.;
		//inputcom_1.drl = inputcom_1.fl;
		//
		//inputcom_1.fl = PrecDat.delz*RoundDoubleToInt(inputcom_1.fl/PrecDat.delz);
		//inputcom_1.dl = PrecDat.delz*RoundDoubleToInt(inputcom_1.dl/PrecDat.delz);
		//inputcom_1.drl = PrecDat.delz*RoundDoubleToInt(inputcom_1.drl/PrecDat.delz);
		
		//debug
		//inputcom_1.f1st = 0.5*inputcom_1.dl; //inputcom_1.drl; //test
		//end debug
		
	if(inputcom_1.quadf == 0) inputcom_1.quadf = 0.0001; //to ensure that FODO exits (necessary to respect eventual gap bw undulator sections)
	if(inputcom_1.quadd == 0) inputcom_1.quadd = 0.0001;
		//inputcom_1.quadf = 0.01;
		//inputcom_1.quadd = 0.01;

//	outputcom_.ifldseed = 0; // to ensure that there is no i/o in fieldseed.f, loadrad.f
//	simcom_1.nwig0 = WigCom.simcom_nwig;

	if((InRad.Power <= 0) && SeedRad.ElecFieldIsDefined()) 
	{//seeding radiation will be read from SRW struct
	 //in that case, dgrid is set from seed radiation grid

		iocom_1.nfin = 1; // default behavior, i.e. opening/closing radiation file should never be used

		//inputcom_1.xlamds = 1.239842E-06/DistrInfoDat.LambStart;
		double cenPhotEnFEL_eV = PrecDat.photEn_xlamds; //DistrInfoDat.LambStart;
		const double redPlanckConst = 6.5821188926E-16; //eV*s
		m_dwSeedRad = (SeedRad.avgPhotEn == cenPhotEnFEL_eV)? 0. : (SeedRad.avgPhotEn - cenPhotEnFEL_eV)/redPlanckConst;
	
		double xRangeSeed = (SeedRad.nx - 1)*SeedRad.xStep;
		double zRangeSeed = (SeedRad.nz - 1)*SeedRad.zStep;
		double maxRange = (xRangeSeed > zRangeSeed)? xRangeSeed : zRangeSeed;
		// "DGRID" - grid size(-dgrid to dgrid) if dgrid > 0 - overrides scaling from rmax0
		inputcom_1.dgrid = 0.5*maxRange; 

		//calculate basic params from tabulated electric field of seed
		//if(PrecDat.alignradf != 0)
		//{
		double sigRe2 = 0;
		SeedRad.CalcBasicRadParamsFromTab(inputcom_1.prad0, inputcom_1.zwaist, sigRe2);
		inputcom_1.zrayl = PI*sigRe2/inputcom_1.xlamds; 
		//double W0 = InRad.WaistDiam; // or 0.5* ?
		//inputcom_1.zrayl = PI*W0*W0/inputcom_1.xlamds; // "ZRAYL"
		//}
	}
	else iocom_1.nfin = 0;

//fill-in/add more params here

	return 0;
}

//*************************************************************************

int srTSASE::AuxvalGenesisPlus()
{
	int result = 0;
////OC-TMP	
//	if(result = auxval_()) return result; /* Auxilliary values (Normalization, mesh, etc...) */
//	//fix mesh size here ??

	return result;
}

//*************************************************************************

int srTSASE::EstimHarmNumFromPhotEn(srTWigComSASE& WigCom)
{
	double K = WigCom.aw0*sqrt(2.);
	double FundPhotEn_eV = 9.5*EbmDat.Energy*EbmDat.Energy/((1 + 0.5*K*K)*WigCom.xlamd);
	double PhotEnObs_eV = PrecDat.photEn_xlamds; //DistrInfoDat.LambStart;
	//Wavelength in [m]: 1.239842E-06/DistrInfoDat.LambStart

	double dHarmNum = PhotEnObs_eV/FundPhotEn_eV;
	int iHarmNum = (int)dHarmNum;
	double deltaHarmNum = dHarmNum - iHarmNum;
	if(deltaHarmNum >= 0.5) iHarmNum++;
	return iHarmNum;
}

//*************************************************************************

int srTSASE::SetupOutputControlStruct(int numHarm)
{
	int result = 0;

	//ControlSASE.ns = timecom_1.nslp + 1;
	//ControlSASE.ns = tbunchcom_1.nslp + 1;
	ControlSASE.ns = wigcom_1.nstepz + 1;

	//ControlSASE.sStep = inputcom_1.xlamd*inputcom_1.delz/TWOPI;
	ControlSASE.sStep = inputcom_1.xlamd*inputcom_1.delz;
	ControlSASE.sStart = -ControlSASE.sStep*(ControlSASE.ns - 1);

	ControlSASE.nt = inputcom_1.nslice;
	ControlSASE.tStep = inputcom_1.xlamds*inputcom_1.zsep/SpeedOfLight;
	ControlSASE.tStart = inputcom_1.xlamds*inputcom_1.zsep*inputcom_1.ntail/SpeedOfLight;

	if(result = pSend->SetupSASEControlStruct(ControlSASE, numHarm)) return result;
	ControlSASE.AllocLocalCont(numHarm);
	//DLL_IMPLEMENT

	double UpdateTime_s = 0.5;
	long TotNumOfCyclesOfStepZ = inputcom_1.nslice*tbunchcom_1.nslp*tbunchcom_1.nsep;
	CompProgressInd.InitializeIndicator(TotNumOfCyclesOfStepZ, UpdateTime_s);

	return result;
}

//*************************************************************************

int srTSASE::DiagnoGenesisPlus(long PassCount, long StepNoVsS, long islice, int numHarm)
{//PassCount = slice number
	int result = 0;
	//const double RadSizeMultip = 1./sqrt(2.);
	const double RadSizeMultip = 1.;

//OC-TMP
	f2c_integer cur_istepz = (f2c_integer)StepNoVsS;
	if(result = diagno_(&cur_istepz)) return result;

	//long StepNoVsS = simcom_1.istepz;
	
	for(int iHarm = 0; iHarm < numHarm; iHarm++)
	{
		if(ControlSASE.ar_pBasePower_vs_s[iHarm] != 0)
		{
			*(ControlSASE.ar_pBasePower_vs_s[iHarm] + StepNoVsS) = (float)(*(diagcom_1.pgainhist + NHMAX*StepNoVsS + iHarm));
		}
		if(ControlSASE.ar_pBasePower_vs_t[iHarm] != 0)
		{
			*(ControlSASE.ar_pBasePower_vs_t[iHarm] + (StepNoVsS*ControlSASE.nt) + (islice - 1)) = (float)(*(diagcom_1.pgainhist + NHMAX*StepNoVsS + iHarm));
		}
		if(ControlSASE.ar_pBaseBunchFact_vs_s[iHarm] != 0)
		{
			*(ControlSASE.ar_pBaseBunchFact_vs_s[iHarm] + StepNoVsS) = (float)(*(diagcom_1.pmodhist + NHMAX*StepNoVsS + iHarm));
		}

		if(StepNoVsS == (ControlSASE.ns - 1))
		{
			if((ControlSASE.ar_pBasePeakPower_vs_s[iHarm] != 0) && (ControlSASE.ar_pBasePower_vs_t[iHarm] != 0))
			{
				for(int is=0; is<ControlSASE.ns; is++)
				{
					float peakPowerForThisLongLos=0;
					FindMaxArrElemInd(ControlSASE.ar_pBasePower_vs_t[iHarm] + (is*ControlSASE.nt), islice, peakPowerForThisLongLos);
					*(ControlSASE.ar_pBasePeakPower_vs_s[iHarm] + is) = peakPowerForThisLongLos;
				}
			}

			if(ControlSASE.ar_pBaseEnergy_vs_s[iHarm] != 0)
			{
				float *t_BaseEnergy_vs_s = ControlSASE.ar_pBaseEnergy_vs_s[iHarm];
				float *t_BasePower_vs_s = ControlSASE.ar_pBasePower_vs_s[iHarm];
				float *t_LocEnergy_vs_s = ControlSASE.ar_pLocEnergy_vs_s[iHarm];
				if(islice <= 1)
				{
					for(int is=0; is<ControlSASE.ns; is++)
					{
						*(t_BaseEnergy_vs_s++) = 0;
						*(t_LocEnergy_vs_s++) = *(t_BasePower_vs_s++);
					}
				}
				else
				{
					double multPowInt = 0.5*ControlSASE.tStep;
					for(int is=0; is<ControlSASE.ns; is++)
					{
						double auxPower = *(t_BasePower_vs_s++);
						*(t_BaseEnergy_vs_s++) += (float)(multPowInt*(*t_LocEnergy_vs_s + auxPower));
						*(t_LocEnergy_vs_s++) = (float)auxPower;
					}
				}
			}
		}
	}

	//if(ControlSASE.pBasePower_vs_s != 0)
	//{
	//	//*(ControlSASE.pBasePower_vs_s + StepNoVsS) = (float)(*(diagcom_1.gain + StepNoVsS));
	//	*(ControlSASE.pBasePower_vs_s + StepNoVsS) = (float)(*(diagcom_1.pgainhist + StepNoVsS*7 + 0)); //fundamental
	//	//diagcom_1.pgainhist[cartcom_1.hloop[n - 1] + diagcom_1.ihist * 7 - 8] = pradn;
	//}
	if(ControlSASE.pBaseRadPhase_vs_s != 0)
	{
		*(ControlSASE.pBaseRadPhase_vs_s + StepNoVsS) = (float)(*(diagcom_1.phimid + StepNoVsS));
	}
	if(ControlSASE.pBaseRadSize_vs_s != 0)
	{
		//*(ControlSASE.pBaseRadSize_vs_s + StepNoVsS) = *(diagcom_1.whalf + StepNoVsS);
		//since Genesis computes sqrt(sigma_x^2 + sigma_z^2)
		*(ControlSASE.pBaseRadSize_vs_s + StepNoVsS) = (float)((*(diagcom_1.whalf + StepNoVsS))*RadSizeMultip);
	}
	if(ControlSASE.pBaseBmSizeX_vs_s != 0)
	{
		*(ControlSASE.pBaseBmSizeX_vs_s + StepNoVsS) = (float)(*(diagcom_1.xrms + StepNoVsS));
	}
	if(ControlSASE.pBaseBmSizeZ_vs_s != 0)
	{
		*(ControlSASE.pBaseBmSizeZ_vs_s + StepNoVsS) = (float)(*(diagcom_1.yrms + StepNoVsS));
	}
	//if(ControlSASE.pBasePower_vs_t != 0)
	//{
	//	*(ControlSASE.pBasePower_vs_t + (StepNoVsS*ControlSASE.nt) + (islice - 1)) = (float)(*(diagcom_1.gain + StepNoVsS));
	//}

	//if(StepNoVsS == (ControlSASE.ns - 1))
	//{
	//	if((ControlSASE.pBasePeakPower_vs_s != 0) && (ControlSASE.pBasePower_vs_t != 0))
	//	{
	//		for(int is=0; is<ControlSASE.ns; is++)
	//		{
	//			float peakPowerForThisLongLos=0;
	//			FindMaxArrElemInd(ControlSASE.pBasePower_vs_t + (is*ControlSASE.nt), islice, peakPowerForThisLongLos);
	//			*(ControlSASE.pBasePeakPower_vs_s + is) = peakPowerForThisLongLos;
	//		}
	//	}
	//	if(ControlSASE.pBaseEnergy_vs_s != 0)
	//	{
	//		float *t_BaseEnergy_vs_s = ControlSASE.pBaseEnergy_vs_s;
	//		float *t_BasePower_vs_s = ControlSASE.pBasePower_vs_s;
	//		float *t_LocEnergy_vs_s = ControlSASE.pLocEnergy_vs_s;
	//		if(islice <= 1)
	//		{
	//			for(int is=0; is<ControlSASE.ns; is++)
	//			{
	//				*(t_BaseEnergy_vs_s++) = 0;
	//				*(t_LocEnergy_vs_s++) = *(t_BasePower_vs_s++);
	//			}
	//		}
	//		else
	//		{
	//			double multPowInt = 0.5*ControlSASE.tStep;
	//			for(int is=0; is<ControlSASE.ns; is++)
	//			{
	//				double auxPower = *(t_BasePower_vs_s++);
	//				*(t_BaseEnergy_vs_s++) += multPowInt*(*t_LocEnergy_vs_s + auxPower);
	//				*(t_LocEnergy_vs_s++) = auxPower;
	//			}
	//		}
	//	}
	//}

	//int maxNumHarm = (inputcom_1.nharm < 5)? inputcom_1.nharm : 5;
	////int maxNumHarm = (inputcom_1.nharm < 7)? inputcom_1.nharm : 7;
	//float *arBaseBunchFact[] = {
	//	ControlSASE.pBaseBunchFact1_vs_s,
	//	ControlSASE.pBaseBunchFact2_vs_s,
	//	ControlSASE.pBaseBunchFact3_vs_s,
	//	ControlSASE.pBaseBunchFact4_vs_s,
	//	ControlSASE.pBaseBunchFact5_vs_s
	//};
	//for(int iHarm = 0; iHarm < maxNumHarm; iHarm++)
	//{
	//	if(arBaseBunchFact[iHarm] != 0)
	//	{
	//		//*(arBaseBunchFact[iHarm] + StepNoVsS) = (float)(*(diagcom_1.pmodhist + 5*StepNoVsS + iHarm));
	//		*(arBaseBunchFact[iHarm] + StepNoVsS) = (float)(*(diagcom_1.pmodhist + 7*StepNoVsS + iHarm));
	//		//diagcom_1.pmodhist[i__ + diagcom_1.ihist * 7 - 8] = sqrt(d__1 * d__1 + d__2 * d__2);
	//	}
	//}

	if(result = UpdateInterface(PassCount)) return result;

	return 0;
}

//*************************************************************************

double srTSASE::RadMeshRange()
{
	double Dmax = 0;
	if(inputcom_1.dgrid != 0.) Dmax = 2*inputcom_1.dgrid;
	else
	{
		double W0 = 0.5*InRad.WaistDiam;
		double Zwaist_d_Zrayl = inputcom_1.zwaist/inputcom_1.zrayl;
		double We2 = W0*W0*(1 + Zwaist_d_Zrayl*Zwaist_d_Zrayl);
		//	double Rmax = (cartcom_1.rmax0)*sqrt((loadcom_1.rxbeam)*(loadcom_1.rxbeam) + (loadcom_1.rybeam)*(loadcom_1.rybeam) + We2);
		Dmax = (inputcom_1.rmax0)*(sqrt((inputcom_1.rxbeam)*(inputcom_1.rxbeam) + (inputcom_1.rybeam)*(inputcom_1.rybeam)) + sqrt(We2));
		//maybe wrong (?), but in agreement with GENESIS
	}
	return Dmax;

/**
		RWAIST = dSQRT(2.d0*ZRAYL/XKS)

		RW0 = RWAIST * dSQRT(1.0d0+(ZWAIST/ZRAYL)**2)
c
c		Maximum radius for the radial mesh (normalized)
		XYMAX = RMAX0 * XKPER0 * (RW0+dSQRT(RXBEAM**2+RYBEAM**2))/2.d0
		DXY=2.D0*XYMAX/float(NCAR-1)            !grid seperation
**/
}

//*************************************************************************

int srTSASE::UpdateInterface(long PassCount)
{
	int result = 0;

	if(result = pSend->FinishWorkingWithControlSASEStruct(ControlSASE)) return result; //to notify wave modifications
//	if(result = pSend->SetupSASEControlStruct(ControlSASE)) return result;
//OC: SetupSASEControlStruct doesn't seem necessary for the undate - to check

	pSend->UpdateInterface();
	//DLL_IMPLEMENT

	//to update control graphs on line

	if(result = CompProgressInd.UpdateIndicator(PassCount)) return result;
	if(result = srYield.Check()) return result;

	return result;
}

//*************************************************************************
//Fill-in output data structures, in order them to be accessed from interface (e.g. Igor Pro)
//int srTSASE::OutResDistrib(f2c_integer _istepz, f2c_integer _islice, f2c_doublereal _xkw0)
int srTSASE::OutResDistrib(long _istepz, long _islice, double _xkw0)
{
	int result = 0;
	if(result = OutElecDistrib(_istepz, _islice)) return result;

	return result;
}

//*************************************************************************

int srTSASE::OutElecDistrib(long _istepz, long _islice)
{
	int result = 0;
	if((EbmDat.pElecDistr == 0) || (EbmDat.nTotMacroPart == 0)) return result;

    if(inputcom_1.ippart <= 0 || inputcom_1.ispart <= 0) return 0;
    if(_istepz % inputcom_1.ippart != 0) return 0;
    if(_islice % inputcom_1.ispart != 0) return 0;

    if(_istepz == 0) rpos_(0, beamcom_1.xpart, beamcom_1.ypart);
	getpsi_(workspace_1.p1);

	if(inputcom_1.npart < simcom_1.npart0) //check for particle loss
	{
		for(long iz = inputcom_1.npart + 1; iz <= simcom_1.npart0; ++iz) 
		{
			beamcom_1.gamma[iz - 1] = (float)-1.; //indicate lost particles with neg.
		}
	}

	float *tElecDistr = EbmDat.pElecDistr + ((_islice - 1)*6*PrecDat.npart);
	for(long iPart=0; iPart<PrecDat.npart; iPart++)
	{
		*(tElecDistr++) = (float)beamcom_1.gamma[iPart];
		*(tElecDistr++) = (float)workspace_1.p1[iPart];
		*(tElecDistr++) = (float)(beamcom_1.xpart[iPart] / simcom_1.xkper0);
		*(tElecDistr++) = (float)(beamcom_1.ypart[iPart] / simcom_1.xkper0);
		*(tElecDistr++) = (float)beamcom_1.px[iPart];
		*(tElecDistr++) = (float)beamcom_1.py[iPart];
	}

//    io___62.ciunit = iocom_1.npar;
//    io___62.cirec = iocom_1.irecpar;
//    s_wdue(&io___62);
//    i__1 = simcom_1.npart0;
//    for (iz = 1; iz <= i__1; ++iz) {
//	do_uio(&c__1, (char *)&beamcom_1.gamma[iz - 1], (ftnlen)sizeof(
//		doublereal));
//    }
//    e_wdue();
//    io___63.ciunit = iocom_1.npar;
//    io___63.cirec = iocom_1.irecpar + 1;
//    s_wdue(&io___63);
//    i__1 = simcom_1.npart0;
//    for (iz = 1; iz <= i__1; ++iz) {
//	do_uio(&c__1, (char *)&workspace_1.p1[iz - 1], (ftnlen)sizeof(
//		doublereal));
//    }
//    e_wdue();
//    io___64.ciunit = iocom_1.npar;
//    io___64.cirec = iocom_1.irecpar + 2;
//    s_wdue(&io___64);
//    i__1 = simcom_1.npart0;
//    for (iz = 1; iz <= i__1; ++iz) {
//	d__1 = beamcom_1.xpart[iz - 1] / simcom_1.xkper0;
//	do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
//    }
//    e_wdue();
//    io___65.ciunit = iocom_1.npar;
//    io___65.cirec = iocom_1.irecpar + 3;
//    s_wdue(&io___65);
//    i__1 = simcom_1.npart0;
//    for (iz = 1; iz <= i__1; ++iz) {
//	d__1 = beamcom_1.ypart[iz - 1] / simcom_1.xkper0;
//	do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
//    }
//    e_wdue();
//    io___66.ciunit = iocom_1.npar;
//    io___66.cirec = iocom_1.irecpar + 4;
//    s_wdue(&io___66);
//    i__1 = simcom_1.npart0;
//    for (iz = 1; iz <= i__1; ++iz) {
//	do_uio(&c__1, (char *)&beamcom_1.px[iz - 1], (ftnlen)sizeof(
//		doublereal));
//    }
//    e_wdue();
//    io___67.ciunit = iocom_1.npar;
//    io___67.cirec = iocom_1.irecpar + 5;
//    s_wdue(&io___67);
//    i__1 = simcom_1.npart0;
//    for (iz = 1; iz <= i__1; ++iz) {
//	do_uio(&c__1, (char *)&beamcom_1.py[iz - 1], (ftnlen)sizeof(
//		doublereal));
//    }
//    e_wdue();
//    iocom_1.irecpar += 6;

/**
      if ((ippart.le.0).or.(ispart.le.0)) return   !no output at all
      if (mod(istepz,ippart).ne.0) return          !output ippartth step
      if (mod(islice,ispart).ne.0) return          !output ispartth slice
c
      if (istepz.eq.0) call rpos(0,xpart,ypart)
      call getpsi(p1)
c
      if (npart.lt.npart0) then     ! check for particle loss
         do iz=npart+1,npart0       ! indicate lost particles with neg. energy
            gamma(iz)=-1.
         enddo
      endif
c
      write(npar,rec=irecpar) (gamma(iz),iz=1,npart0)
      write(npar,rec=irecpar+1) (p1(iz),iz=1,npart0)
      write(npar,rec=irecpar+2) (xpart(iz)/xkper0,iz=1,npart0)
      write(npar,rec=irecpar+3) (ypart(iz)/xkper0,iz=1,npart0)
      write(npar,rec=irecpar+4) (px(iz),iz=1,npart0)
      write(npar,rec=irecpar+5) (py(iz),iz=1,npart0)
      irecpar=irecpar+6
**/

	return result;
}

//*************************************************************************

int srTSASE::OutDumpResDistrib(int _islice)
{
	int result = 0;
	if((EbmDat.pElecDistr == 0) || (EbmDat.nTotMacroPart == 0)) return result;

	if(_islice <= iocom_1.firstout) return 0;
	//particle distribution

	//suppress for IOTAIL=0
    //if(inputcom_1.idmppar != 0) 
	//{
	if(inputcom_1.npart < simcom_1.npart0) 
	{
		//i__1 = simcom_1.npart0;
		for(long i__ = inputcom_1.npart + 1; i__ <= simcom_1.npart0; ++i__) 
		{/* check for particle losses */
			beamcom_1.gamma[i__ - 1] = (float)-1.; /* indicate lost particles wit */
		}
	}

	float *tElecDistr = EbmDat.pElecDistr + ((_islice - 1)*6*PrecDat.npart);
	for(long iPart=0; iPart<PrecDat.npart; iPart++)
	{
		*(tElecDistr++) = (float)beamcom_1.gamma[iPart];
		*(tElecDistr++) = (float)beamcom_1.theta[iPart]; //this phase is correct for future runs
		*(tElecDistr++) = (float)(beamcom_1.xpart[iPart] / simcom_1.xkper0);
		*(tElecDistr++) = (float)(beamcom_1.ypart[iPart] / simcom_1.xkper0);
		*(tElecDistr++) = (float)beamcom_1.px[iPart];
		*(tElecDistr++) = (float)beamcom_1.py[iPart];
	}

///*       writing the record */
//
//	j = (*islice - iocom_1.firstout - 1) * 6 + 1;
//	io___77.ciunit = iocom_1.ndmp2;
//	io___77.cirec = j;
//	s_wdue(&io___77);
//	i__1 = simcom_1.npart0;
//	for (i__ = 1; i__ <= i__1; ++i__) {
//	    do_uio(&c__1, (char *)&beamcom_1.gamma[i__ - 1], (ftnlen)sizeof(
//		    doublereal));
//	}
//	e_wdue();
//	io___78.ciunit = iocom_1.ndmp2;
//	io___78.cirec = j + 1;
//	s_wdue(&io___78);
//	i__1 = simcom_1.npart0;
//	for (i__ = 1; i__ <= i__1; ++i__) {
//	    do_uio(&c__1, (char *)&beamcom_1.theta[i__ - 1], (ftnlen)sizeof(
//		    doublereal));
//	}
//	e_wdue();
//	io___79.ciunit = iocom_1.ndmp2;
//	io___79.cirec = j + 2;
//	s_wdue(&io___79);
//	i__1 = simcom_1.npart0;
//	for (i__ = 1; i__ <= i__1; ++i__) {
//	    d__1 = beamcom_1.xpart[i__ - 1] / simcom_1.xkper0;
//	    do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
//	}
//	e_wdue();
//	io___80.ciunit = iocom_1.ndmp2;
//	io___80.cirec = j + 3;
//	s_wdue(&io___80);
//	i__1 = simcom_1.npart0;
//	for (i__ = 1; i__ <= i__1; ++i__) {
//	    d__1 = beamcom_1.ypart[i__ - 1] / simcom_1.xkper0;
//	    do_uio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
//	}
//	e_wdue();
//	io___81.ciunit = iocom_1.ndmp2;
//	io___81.cirec = j + 4;
//	s_wdue(&io___81);
//	i__1 = simcom_1.npart0;
//	for (i__ = 1; i__ <= i__1; ++i__) {
//	    do_uio(&c__1, (char *)&beamcom_1.px[i__ - 1], (ftnlen)sizeof(
//		    doublereal));
//	}
//	e_wdue();
//	io___82.ciunit = iocom_1.ndmp2;
//	io___82.cirec = j + 5;
//	s_wdue(&io___82);
//	i__1 = simcom_1.npart0;
//	for (i__ = 1; i__ <= i__1; ++i__) {
//	    do_uio(&c__1, (char *)&beamcom_1.py[i__ - 1], (ftnlen)sizeof(
//		    doublereal));
//	}
//	e_wdue();
//    }
//
///*     field distribution */
//
//    if (inputcom_1.idmpfld == 0) {
//	return 0;
//    }
//
///*     problems arise if the dump is used for another run */
///*     any change in undulator period. the field has twice a scaling */
///*     with xkper0 - 1. eikonal equation + 2. normalization of dxy */
//
//    scltmp = cartcom_1.dxy * 510999.06 * simcom_1.xkper0 / cartcom_1.xks / 
//	    sqrt(376.73);
//
//
//    io___84.ciunit = iocom_1.ndump;
//    io___84.cirec = *islice - iocom_1.firstout;
//    s_wdue(&io___84);
//    i__1 = inputcom_1.ncar * inputcom_1.ncar;
//    for (i__ = 1; i__ <= i__1; ++i__) {
//	i__2 = i__ - 1;
//	z__2.r = scltmp * cartcom_1.crfield[i__2].r, z__2.i = scltmp * 
//		cartcom_1.crfield[i__2].i;
//	z__1.r = z__2.r, z__1.i = z__2.i;
//	do_uio(&c__2, (char *)&z__1, (ftnlen)sizeof(doublereal));
//    }
//    e_wdue();
//
//    if (inputcom_1.itdp == 0 || *islice < inputcom_1.nslice) {
//	return 0;
//    }
//
//    for (j = tbunchcom_1.nslp - 1; j >= 1; --j) {
///* dump field , escaping beam */
//	pulltimerec_(workspace_1.crwork3, &inputcom_1.ncar, &j);
//	i0 = inputcom_1.nslice + tbunchcom_1.nslp - j - iocom_1.firstout;
//	io___86.ciunit = iocom_1.ndump;
//	io___86.cirec = i0;
//	s_wdue(&io___86);
//	i__2 = inputcom_1.ncar * inputcom_1.ncar;
//	for (i__ = 1; i__ <= i__2; ++i__) {
//	    i__1 = i__ - 1;
//	    z__2.r = scltmp * workspace_1.crwork3[i__1].r, z__2.i = scltmp * 
//		    workspace_1.crwork3[i__1].i;
//	    z__1.r = z__2.r, z__1.i = z__2.i;
//	    do_uio(&c__2, (char *)&z__1, (ftnlen)sizeof(doublereal));
//	}
//	e_wdue();
//    }

/**
      if (islice.le.firstout) return       ! suppress for IOTAIL=0
c
c     particle distribution
c
      if (idmppar.ne.0) then
        if (npart.lt.npart0) then
           do i=npart+1,npart0             ! check for particle losses
              gamma(i)=-1.                 ! indicate lost particles with neg. energy
           enddo
        endif
c                              
c       writing the record
c
          j=6*(islice-firstout-1)+1
          write(ndmp2,rec=j)   (gamma(i),i=1,npart0)
          write(ndmp2,rec=j+1) (theta(i),i=1,npart0)
          write(ndmp2,rec=j+2) (xpart(i)/xkper0,i=1,npart0)
          write(ndmp2,rec=j+3) (ypart(i)/xkper0,i=1,npart0)
          write(ndmp2,rec=j+4) (px(i),i=1,npart0)
          write(ndmp2,rec=j+5) (py(i),i=1,npart0)
      endif
c
c     field distribution
c
      if (idmpfld.eq.0) return     
c
c     problems arise if the dump is used for another run
c     any change in undulator period. the field has twice a scaling
c     with xkper0 - 1. eikonal equation + 2. normalization of dxy
c
      scltmp=dxy*eev*xkper0/xks/dsqrt(vacimp)   !
c
      write(ndump,rec=islice-firstout) (crfield(i)*scltmp,i=1,ncar*ncar)
c
      if ((itdp.eq.0).or.(islice.lt.nslice)) return 
c
      do j=nslp-1,1,-1                     ! dump field , escaping beam
        call pulltimerec(crwork3,ncar,j)
        i0=nslice+nslp-j-firstout
        write(ndump,rec=i0) (crwork3(i)*scltmp,i=1,ncar*ncar)
      enddo
**/
	return result;
}

//*************************************************************************
//Hacked version of the GENESIS's loadbeam: to allow electron distribution
//input from SRW / Igor
//*************************************************************************
int srTSASE::loadbeam_srw(f2c_integer *islice, f2c_doublereal *xkper0)
{
    /* System generated locals */
	f2c_integer i__1, i__2;
	f2c_doublereal d__1, d__2, d__3;

    /* Builtin functions */
    //double sin(), sqrt();

    /* Local variables */
    //extern integer readpart_();
    //extern /* Subroutine */ int loaddist_();
    //extern integer printerr_();
    static f2c_integer i__;
    //extern /* Subroutine */ int shotnoise_penman__();
    //extern doublereal hammv_();
    static f2c_integer mpart;
    //extern /* Subroutine */ int shotnoise_fawley__(), loadquiet_();
    static f2c_integer ip;

/*     =================================================================== */
/*     this routine fills the phase space for one slice of phase space */
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
/*     ------------------------------------------------------------------ */
/*     input/output control */
/*     initialize particle loss */

    beamcom_1.lost = 0;
    i__1 = inputcom_1.npart;
    for (i__ = 1; i__ <= i__1; ++i__) {
	beamcom_1.lostid[i__ - 1] = 0;
    }

    if (inputcom_1.iall != 0) {
	inputcom_1.ildpsi = -abs(inputcom_1.ildpsi);
	//reinitialize all hammersley sequences
    }
/*     fill phase */
    mpart = inputcom_1.npart / inputcom_1.nbins;

/* particles per bin */
    i__1 = mpart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.theta[ip - 1] = hammv_(&inputcom_1.ildpsi) * (float)2. * 3.14159265358979 / (f2c_doublereal) inputcom_1.nbins - 3.14159265358979;
/* load in first bin */
    }

/*     branches for different loading methods */
    //if (iocom_1.npin > 0) {
    if (m_ElecDistribShouldBeUsed) { //OC port
	//i__ = readpart_(islice);
	i__ = readpart_srw(islice); //OC hacked version; to allow electron distribution input from SRW

/* if reading from file */
	if (inputcom_1.convharm > 1 && inputcom_1.multconv != 0) {
	    shotnoise_fawley__();
	}
	return 0;
/* skip the rest (loading is done) */
    }

/*     btpar is calculated in readpart */
    if (iocom_1.ndis > 0) {
	loaddist_(&mpart, islice);
/* load from distribution file */
    } else {
	loadquiet_(&mpart);
/* internal load (quiet start) */
    }

/*     normalized transverse position */
    i__1 = mpart;
    for (ip = 1; ip <= i__1; ++ip) {
	beamcom_1.xpart[ip - 1] *= *xkper0;
	beamcom_1.ypart[ip - 1] *= *xkper0;
    }

/*     mirror particle in remaining bins */
    i__1 = inputcom_1.nbins - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = mpart;
	for (ip = 1; ip <= i__2; ++ip) {
	    beamcom_1.xpart[ip + i__ * mpart - 1] = beamcom_1.xpart[ip - 1];
	    beamcom_1.ypart[ip + i__ * mpart - 1] = beamcom_1.ypart[ip - 1];
	    beamcom_1.px[ip + i__ * mpart - 1] = beamcom_1.px[ip - 1];
	    beamcom_1.py[ip + i__ * mpart - 1] = beamcom_1.py[ip - 1];
	    beamcom_1.gamma[ip + i__ * mpart - 1] = beamcom_1.gamma[ip - 1];
	    beamcom_1.theta[ip + i__ * mpart - 1] = beamcom_1.theta[ip - 1] + 
		    (f2c_doublereal) i__ * (float)2. * 3.14159265358979 / (f2c_doublereal) inputcom_1.nbins;
	    beamcom_1.lostid[ip + i__ * mpart - 1] = beamcom_1.lostid[ip - 1];
	}
    }

/*     add shotnoise */
    if (beamcom_1.xcuren < 0.) { //OC to implement
	//i__ = printerr_(&c_n23, "xcuren <=0 in loadbeam", (ftnlen)22);
	beamcom_1.xcuren = 0.;
    } else {
	if (inputcom_1.isntyp == 0) {
	    shotnoise_fawley__();
	} else {
	    shotnoise_penman__();
	}
    }

/*     add longitudinal correlations (energy modulation, prebunching) */
    if (inputcom_1.bunch != (float)0. || inputcom_1.emod != (float)0.) {
	i__1 = inputcom_1.npart;
	for (ip = 1; ip <= i__1; ++ip) {
	    beamcom_1.gamma[ip - 1] -= inputcom_1.emod * sin(beamcom_1.theta[
		    ip - 1] - inputcom_1.emodphase);
	    beamcom_1.theta[ip - 1] -= inputcom_1.bunch * (float)2. * sin(beamcom_1.theta[ip - 1] - inputcom_1.bunchphase);
/*              bunching is J_1(2.*bunch), approximately = bunch */
	}
    }

/*     calculate init. parallel velocity (needed in first call of track) */
    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
/* Computing 2nd power */
	d__1 = beamcom_1.px[ip - 1];
/* Computing 2nd power */
	d__2 = beamcom_1.py[ip - 1];
/* Computing 2nd power */
	d__3 = beamcom_1.gamma[ip - 1];
	beamcom_1.btpar[ip - 1] = sqrt(1. - (d__1 * d__1 + d__2 * d__2 + (float)1.) / (d__3 * d__3));
/* parallel velocity */
    }
	return 0;
}

//int srTSASE::loadbeam_srw(f2c_integer *islice, f2c_doublereal *xkper0)
//{
//	// System generated locals 
//	f2c_integer i__1, i__2;
//	f2c_doublereal d__1, d__2, d__3;
//
//	// Builtin functions
//	//double sin(), sqrt();
//
//    // Local variables
//	//extern f2c_integer readpart_();
//	//extern loaddist_();
//	//extern loaddist_(f2c_integer*, f2c_integer*);
//
//	//extern f2c_integer printerr_();
//	static f2c_integer i__;
//	//extern int shotnoise_penman__();
//	//extern f2c_doublereal hammv_();
//	//extern f2c_doublereal hammv_(f2c_integer*);
//	static f2c_integer mpart;
//	//extern int shotnoise_fawley__(), loadquiet_();
//	//extern int shotnoise_fawley__();
//	//extern int loadquiet_(f2c_integer*);
//
//	static f2c_integer ip;
//	static f2c_integer c_n23 = -23;
//
//	//     ===================================================================
//	//     this routine fills the phase space for one slice of phase space
//	//     ------------------------------------------------------------------
//
//	//initialize particle loss
//
//    beamcom_1.lost = 0;
//    i__1 = inputcom_1.npart;
//    for (i__ = 1; i__ <= i__1; ++i__) {
//	beamcom_1.lostid[i__ - 1] = 0;
//    }
//
//    if (inputcom_1.iall != 0) {
//	inputcom_1.ildpsi = -abs(inputcom_1.ildpsi);
//	//reinitialize all hammersley sequences
//    }
//
//	//fill phase
//
//    mpart = inputcom_1.npart / inputcom_1.nbins;
//
//	//particles per bin
//    i__1 = mpart;
//    for (ip = 1; ip <= i__1; ++ip) {
//	beamcom_1.theta[ip - 1] = hammv_(&inputcom_1.ildpsi) * (float)2. * 
//		3.14159265358979 / (f2c_doublereal) inputcom_1.nbins - 3.14159265358979;
//	//load in first bin
//    }
//
//	//branches for different loading methods
//
//    //if (iocom_1.npin > 0) {
//    if(m_ElecDistribShouldBeUsed) 
//	{
//		//i__ = readpart_(islice);
//		i__ = readpart_srw(islice); //OC hacked version; to allow electron distribution input from SRW
//		//if reading from file
//		if (inputcom_1.convharm > 1 && inputcom_1.multconv != 0) {
//			shotnoise_fawley__();
//		}
//		return 0;
//		//skip the rest (loading is done)
//    }
//
//	//btpar is calculated in readpart
//
//    if (iocom_1.ndis > 0) {
//	loaddist_(&mpart, islice);
//	//load from distribution file
//    } else {
//	loadquiet_(&mpart); //internal load (quiet start)
//    }
//
//    //normalized transverse position
//
//    i__1 = mpart;
//    for (ip = 1; ip <= i__1; ++ip) {
//	beamcom_1.xpart[ip - 1] *= *xkper0;
//	beamcom_1.ypart[ip - 1] *= *xkper0;
//    }
//
//    //mirror particle in remaining bins
//
//    i__1 = inputcom_1.nbins - 1;
//    for (i__ = 1; i__ <= i__1; ++i__) {
//	i__2 = mpart;
//	for (ip = 1; ip <= i__2; ++ip) {
//	    beamcom_1.xpart[ip + i__ * mpart - 1] = beamcom_1.xpart[ip - 1];
//	    beamcom_1.ypart[ip + i__ * mpart - 1] = beamcom_1.ypart[ip - 1];
//	    beamcom_1.px[ip + i__ * mpart - 1] = beamcom_1.px[ip - 1];
//	    beamcom_1.py[ip + i__ * mpart - 1] = beamcom_1.py[ip - 1];
//	    beamcom_1.gamma[ip + i__ * mpart - 1] = beamcom_1.gamma[ip - 1];
//	    beamcom_1.theta[ip + i__ * mpart - 1] = beamcom_1.theta[ip - 1] + 
//		    (f2c_doublereal) i__ * (float)2. * 3.14159265358979 / (f2c_doublereal) inputcom_1.nbins;
//	    beamcom_1.lostid[ip + i__ * mpart - 1] = beamcom_1.lostid[ip - 1];
//	}
//    }
//
//    //add shotnoise
//
//    if (beamcom_1.xcuren < 0.) {
//		//i__ = printerr_(&c_n23, "xcuren <=0 in loadbeam", (f2c_ftnlen)22);
//		beamcom_1.xcuren = 0.;
//    } else {
//	if (inputcom_1.isntyp == 0) {
//	    shotnoise_fawley__();
//	} else {
//	    shotnoise_penman__();
//	}
//    }
//
//    //add longitudinal correlations (energy modulation, prebunching)
//
//    if (inputcom_1.bunch != (float)0. || inputcom_1.emod != (float)0.) {
//	i__1 = inputcom_1.npart;
//	for (ip = 1; ip <= i__1; ++ip) {
//	    beamcom_1.gamma[ip - 1] -= inputcom_1.emod * sin(beamcom_1.theta[
//		    ip - 1] - inputcom_1.emodphase);
//	    beamcom_1.theta[ip - 1] -= inputcom_1.bunch * (float)2. * sin(
//		    beamcom_1.theta[ip - 1] - inputcom_1.bunchphase);
//			//bunching is J_1(2.*bunch), approximately = bunch
//	}
//    }
//
//	//calculate init. parallel velocity (needed in first call of track)
//
//    i__1 = inputcom_1.npart;
//    for (ip = 1; ip <= i__1; ++ip) {
//	//Computing 2nd power
//	d__1 = beamcom_1.px[ip - 1];
//	//Computing 2nd power
//	d__2 = beamcom_1.py[ip - 1];
//	//Computing 2nd power
//	d__3 = beamcom_1.gamma[ip - 1];
//	beamcom_1.btpar[ip - 1] = sqrt(1. - (d__1 * d__1 + d__2 * d__2 + (float)1.) / (d__3 * d__3));
//    //parallel velocity
//    }
//    return 0;
//}//loadbeam_srw

//*************************************************************************
//Hacked version of the GENESIS's loadbeam: to allow electron distribution
//input from SRW / Igor
//*************************************************************************
int srTSASE::readpart_srw(f2c_integer *islice)
{
	//System generated locals
	f2c_integer ret_val; //, i__1, i__2;
	//f2c_doublereal d__1, d__2, d__3;

    //Builtin functions
	//f2c_integer s_rdue(), do_uio(), e_rdue();

    //Local variables

	//static f2c_integer c__1 = 1;
	static f2c_integer idel;
    extern int last_();
    //extern f2c_integer printerr_();
    //static f2c_integer i__, j;
    //extern int importdispersion_();

	//Fortran I/O blocks
	//static cilist io___58 = { 1, 0, 0, 0, 0 };
	//static cilist io___60 = { 1, 0, 0, 0, 0 };
	//static cilist io___61 = { 1, 0, 0, 0, 0 };
	//static cilist io___62 = { 1, 0, 0, 0, 0 };
	//static cilist io___63 = { 1, 0, 0, 0, 0 };
	//static cilist io___64 = { 1, 0, 0, 0, 0 };

	//     =================================================================
	//     load complete set of particle from file
	//     -----------------------------------------------------------------

	ret_val = 0;
	long offsetElecDistr = 0;

	//OC (the number of slices and number of particles should be the same as in input particle distribution)
	long indPart = (*islice - 1)*simcom_1.npart0;
	if(indPart > EbmDat.nTotMacroPart) return ret_val; 

	if (inputcom_1.multconv == 0) { 
		//harmonics can be calculated accuretely even with this option
		//only the accuracy is affected

		//j = (*islice - 1) * 6 + 1;
		offsetElecDistr = (*islice - 1)*6*simcom_1.npart0;
		//reading every slice
	} else {
		//j = (*islice - 1) / inputcom_1.convharm * 6 + 1;
		offsetElecDistr = (long)((*islice - 1)/inputcom_1.convharm*6*simcom_1.npart0);
		//reading every (1/convharm)th slice
    }

	float *tElecDistr = EbmDat.pElecDistr + offsetElecDistr;
	for(long i=0; i<simcom_1.npart0; i++)
	{
		beamcom_1.gamma[i] = *(tElecDistr++);
		beamcom_1.theta[i] = *(tElecDistr++);
		beamcom_1.xpart[i] = *(tElecDistr++);
		beamcom_1.ypart[i] = *(tElecDistr++);
		beamcom_1.px[i] = *(tElecDistr++);
		beamcom_1.py[i] = *(tElecDistr++);
	}

	//apply transfermatrix to particle distribution
    importtransfer_();

	//apply dispersive section to it
	//OC: to check how dispersion section is defined (IDRIL)
	importdispersion_();

	//convert to higher harmonic
	if (inputcom_1.convharm > 1) {
		for (long i__ = 0; i__ < inputcom_1.npart; i__++) {
			beamcom_1.theta[i__] = (f2c_real) inputcom_1.convharm * beamcom_1.theta[i__];
		}
    }

	//calculate init. perpendicular velocity (needed in first call of track)
	for (long i__ = 0; i__ < simcom_1.npart0; i__++) {
		beamcom_1.xpart[i__] *= simcom_1.xkper0;
		beamcom_1.ypart[i__] *= simcom_1.xkper0;
		
        double d__1 = beamcom_1.px[i__]; //computing 2nd power
		double d__2 = beamcom_1.py[i__]; //computing 2nd power
		double d__3 = beamcom_1.gamma[i__]; //computing 2nd power
		beamcom_1.btpar[i__] = sqrt(1. - (d__1*d__1 + d__2*d__2 + (double)1.)/(d__3*d__3));
        //parallel velocity
	}

	//check for particle losses from previous run
	idel = 0;
	for (long i__ = 0; i__ < inputcom_1.npart; i__++) {
		if (beamcom_1.gamma[i__] > (float)0.) {
			beamcom_1.gamma[i__ - idel] = beamcom_1.gamma[i__];
			beamcom_1.theta[i__ - idel] = beamcom_1.theta[i__];
			beamcom_1.xpart[i__ - idel] = beamcom_1.xpart[i__];
			beamcom_1.ypart[i__ - idel] = beamcom_1.ypart[i__];
			beamcom_1.px[i__ - idel] = beamcom_1.px[i__];
			beamcom_1.py[i__ - idel] = beamcom_1.py[i__];
		} else {
			++idel;
		}
	}
	inputcom_1.npart = simcom_1.npart0 - idel;
	beamcom_1.xcuren = beamcom_1.xcuren*((f2c_real)inputcom_1.npart)/((f2c_real)simcom_1.npart0);

//L100:
//    ret_val = printerr_(&c_n2, inputcom_1.partfile, (ftnlen)30);
//    last_();
	return ret_val;
}

//*************************************************************************

int srTSASE::initrun_srw()
{
	static f2c_integer c__1 = 1;

    /* System generated locals */
	f2c_integer i__1;
	f2c_doublereal d__1, d__2;

    /* Local variables */
    //extern /* Subroutine */ int magfield_(), last_(), scaninit_();
    //extern integer printerr_();
	static f2c_integer ip;
	static f2c_doublereal xi;
    //extern doublereal ran1_();
    //extern /* Subroutine */ int loadslpfld_(), getdiag_();

/*     ================================================================== */
/*     initialize the run by setting up */
/*     the precalculated matrices for the field solver and */
/*     claculating/normalizing some auxiliary variables */
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
/*     simulation control and normalisation parameter */
/*     -------------------------------------------------------------------- */
/*     cartesian mesh */
/*     -------------------------------------------------------------------- */
/*     time dependency of bunch + slippage field */
/*     ------------------------------------------------------------------ */
/*     wiggler parameters */
/*     seeding of random number generator */

    xi = ran1_(&inputcom_1.iseed);

/*      initiate loop for higher harmonics */

/* init ran1 */
    beamcom_1.dedz = inputcom_1.eloss * inputcom_1.delz * inputcom_1.xlamd / 510999.06;
    beamcom_1.xcuren = inputcom_1.curpeak;
    simcom_1.npart0 = inputcom_1.npart;

/*     normalizations */
/*     ------------------------------------------------------------------ */

/* save total number or part */
    simcom_1.xkw0 = 6.28318530717958 / inputcom_1.xlamd;
/* wiggler wavenumber */
    simcom_1.xkper0 = simcom_1.xkw0;

/*     magnetic field */
/*     ------------------------------------------------------------------ */

/* transverse normalisation */
    magfield_(&simcom_1.xkper0, &c__1);

/* magnetic field descript */
    if (simcom_1.inorun != 0) {//OC: to re-implement?
	//ip = printerr_(&c_n24, "Termination enforced by user", (ftnlen)28);
	//last_();
    }

/*     slipping length */

    tbunchcom_1.nsep = (f2c_integer) (inputcom_1.zsep / inputcom_1.delz);
/* steps between field slip */
    tbunchcom_1.nslp = wigcom_1.nstepz / tbunchcom_1.nsep;
/* total slippage steps */
    if (wigcom_1.nstepz % tbunchcom_1.nsep != 0) {
	++tbunchcom_1.nslp;
    }
/* would be shorter */

/*     contruct grid properties (grid spacing, precalculated matrices) */

/* if not added the effecti */
    cartcom_1.dxy = simcom_1.xkw0 * 2. * inputcom_1.dgrid / (f2c_real) (inputcom_1.ncar - 1);
    cartcom_1.xks = 6.28318530717958 / inputcom_1.xlamds;

	int result = 0;
	if(inputcom_1.itdp != 0) //OC191108
	{//to check!
		if(result = Alloc_tslipcom(PrecDat.ncar, tbunchcom_1.nslp)) return result; //OC port: separate allocation of GENESIS crtime buffer
	}

/*     time dependencies */
	//loadslpfld_(&tbunchcom_1.nslp);
	loadslpfld_srw(&tbunchcom_1.nslp); //OC port

/*     scanning */
/* input field for first slice and seed */
    scaninit_();

/*     matrix initialization */
/* initialize scanning */
    d__1 = inputcom_1.delz * inputcom_1.xlamd;
    d__2 = cartcom_1.dxy / simcom_1.xkper0;
    getdiag_(&d__1, &d__2, &cartcom_1.xks);

/*     clear space charge field for case that space charge is disabled */

    i__1 = inputcom_1.npart;
    for (ip = 1; ip <= i__1; ++ip) {
/* clear space charge term */
	beamcom_1.ez[ip - 1] = 0.;
    }
/* ip */
    return 0;
}

//int srTSASE::initrun_srw()
//{
//	static f2c_integer c__1 = 1;
//	static f2c_integer c_n24 = -24;
//
//	/* System generated locals */
//	int i__1;
//	double d__1, d__2;
//
//	/* Local variables */
//	//extern /* Subroutine */ int magfield_(f2c_doublereal*, f2c_integer*), scaninit_(); //, last_();
//	//extern integer printerr_();
//	//extern f2c_doublereal bessj0_(f2c_doublereal*), bessj1_(f2c_doublereal*);
//	static f2c_integer ip;
//	static f2c_doublereal xi;
//	//extern /* Subroutine */ int loadslpfld_(f2c_integer*), getdiag_(f2c_doublereal*, f2c_doublereal*, f2c_doublereal*);
//	//extern f2c_doublereal ran1_(f2c_integer*);
//
//	/*     ================================================================== */
//	/*     initialize the run by setting up */
//	/*     the precalculated matrices for the field solver and */
//	/*     claculating/normalizing some auxiliary variables */
//	/*     ------------------------------------------------------------------ */
//	/*     error codes */
//	/* genesis version */
//	/* platform */
//	/* indicator for original fil */
//	/* indicator for sdds filetyp */
//	/* # of particles */
//	/* # of integration steps */
//	/* # of slices */
//	/* maximum of harmonics */
//	/* maximum of particle i */
//	/* <> 0 keeps distribution in */
//	/* energy units (mc^2) in ev */
//	/* vacuum impedence in ohms */
//	/* speed of light * electron */
//	/* pi */
//	/* pi/2 */
//	/* 2*pi */
//	/* check i for precission */
//	/* check ii for precission */
//	/* number of radial points fo */
//	/* # of gridpoints of cartesi */
//	/*     function prototypes */
//	/*     ------------------------------------------------------ */
//	/*     all input variables */
//	/*     wiggler */
//	/*     electron beam */
//	/*     radiation */
//	/*     grid-quantities */
//	/*     control */
//	/*     strong focusing */
//	/*     loading */
//	/*     output */
//	/*     external files */
//	/*     time-dependency */
//	/*     scan */
//	/*     extension */
//	/*     ------------------------------------------------------------------ */
//	/*     electron beam */
//	/*     simulation control and normalisation parameter */
//	/*     -------------------------------------------------------------------- */
//	/*     cartesian mesh */
//	/*     -------------------------------------------------------------------- */
//	/*     time dependency of bunch + slippage field */
//	/*     ------------------------------------------------------------------ */
//	/*     wiggler parameters */
//	/*     seeding of random number generator */
//
//	xi = ran1_(&inputcom_1.iseed);
//
//	/*     calculate coupling factor if necessary */
//	/* init ran1 */
//	if (inputcom_1.fbess0 == (float)0.) {
//		if (inputcom_1.iwityp == 0) {
//			/* Computing 2nd power */
//			d__1 = inputcom_1.aw0;
//			/* Computing 2nd power */
//			d__2 = inputcom_1.aw0;
//			xi = d__1 * d__1 / (d__2 * d__2 + 1.) / 2;
//			inputcom_1.fbess0 = bessj0_(&xi) - bessj1_(&xi);
//		} else {
//			inputcom_1.fbess0 = (float)1.;
//		}
//	}
//
//	beamcom_1.dedz = inputcom_1.eloss * inputcom_1.delz * inputcom_1.xlamd / 510999.06;
//	beamcom_1.xcuren = inputcom_1.curpeak;
//	simcom_1.npart0 = inputcom_1.npart;
//
//	/*     normalizations */
//	/*     ------------------------------------------------------------------ */
//	/* save total number or part */
//	simcom_1.xkw0 = 6.28318530717958 / inputcom_1.xlamd;
//	/* wiggler wavenumber */
//	simcom_1.xkper0 = simcom_1.xkw0;
//
//	/*     magnetic field */
//	/*     ------------------------------------------------------------------ */
//	/* transverse normalisation */
//	magfield_(&simcom_1.xkper0, &c__1);
//
//	/* magnetic field descript */
//	if (simcom_1.inorun != 0) {//OC: to re-implement?
//		//ip = printerr_(&c_n24, "Termination enforced by user", (ftnlen)28);
//		//last_();
//	}
//
//	/*     slipping length */
//	tbunchcom_1.nsep = (f2c_integer) (inputcom_1.zsep / inputcom_1.delz);
//	/* steps between field slip */
//	tbunchcom_1.nslp = wigcom_1.nstepz / tbunchcom_1.nsep;
//	/* total slippage steps */
//	if (wigcom_1.nstepz % tbunchcom_1.nsep != 0) {
//		++tbunchcom_1.nslp;
//	}
//	/* would be shorter */
//	/*     contruct grid properties (grid spacing, precalculated matrices) */
//	/* if not added the effecti */
//
//	cartcom_1.dxy = simcom_1.xkw0 * 2. * inputcom_1.dgrid / (f2c_real) (inputcom_1.ncar - 1);
//	cartcom_1.xks = 6.28318530717958 / inputcom_1.xlamds;
//
//	int result = 0;
//	if(result = Alloc_tslipcom(PrecDat.ncar, tbunchcom_1.nslp)) return result; //OC: separate allocation of GENESIS crtime buffer
//
//	/*     time dependencies */
//	//loadslpfld_(&tbunchcom_1.nslp);
//	loadslpfld_srw(&tbunchcom_1.nslp);
//
//	/*     scanning */
//	/* input field for first slice and seed */
//	scaninit_();
//
//	/*     matrix initialization */
//	/* initialize scanning */
//	d__1 = inputcom_1.delz * inputcom_1.xlamd;
//	d__2 = cartcom_1.dxy / simcom_1.xkper0;
//	getdiag_(&d__1, &d__2, &cartcom_1.xks);
//
//	/*     clear space charge field for case that space charge is disabled */
//	i__1 = inputcom_1.npart;
//	for (ip = 1; ip <= i__1; ++ip) {
//		/* clear space charge term */
//		beamcom_1.ez[ip - 1] = 0.;
//	}
//	/* ip */
//	return 0;
//} /* initrun_srw */

//*************************************************************************

int srTSASE::loadslpfld_srw(f2c_integer* nslp)
{
	static f2c_integer c__1 = 1;

    /* System generated locals */
    f2c_integer i__1, i__2, i__3;
    f2c_doublereal d__1, d__2;

    /* Local variables */
	static f2c_integer irec, islp, ierr;
	//extern integer printerr_();
	//extern /* Subroutine */ int last_();
	static f2c_integer i__;
	//extern integer readfield_();
	//extern /* Subroutine */ int gauss_hermite__(), dotimerad_();
	static f2c_integer ix;
	//extern /* Subroutine */ int pushtimerec_();

/*     ========================================================= */
/*     fills the array crtime with a seeding field for the first */
/*     slice. */
/*     --------------------------------------------------------- */
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
/*     ------------------------------------------------------------------ */
/*     input/output control */
/*     file:   timerec.com */
/*     time dependent common block - huge!!!! */
/*     ------------------------------------------------------------------ */
/*     -------------------------------------------------------------------- */
/*     slippage field */
/* must be > than ncar*ncar*nslp(98) */
/* <>0 -> slippage is stored on disk */
/*     ------------------------------------------------------------------ */
/*     diagnostic */

	if (inputcom_1.itdp == 0) {
	return 0;
	}

/*     check for limitation in timerecord */
/*     --------------------------------------------------------------- */
    if (inputcom_1.nslice < *nslp * (1 - inputcom_1.iotail)) { //OC: to implement
	//ierr = printerr_(&c_n20, "no output - ntail too small", (ftnlen)27);
    }
    //if (TRUE_) { //OC: to implement
	if (*nslp * inputcom_1.ncar * inputcom_1.ncar * cartcom_1.nhloop > 79803500) {
	    //ierr = printerr_(&c_n22, " ", (ftnlen)1);
	    //last_();
	}
    //}

/*     load initial slippage field from file or internal generation */
/*     ---------------------------------------------------------------- */
    i__1 = *nslp - 1;
    for (islp = 1; islp <= i__1; ++islp) {

	i__2 = inputcom_1.ncar * inputcom_1.ncar * cartcom_1.nhloop;
	for (ix = 1; ix <= i__2; ++ix) {
	    i__3 = ix - 1;
	    workspace_1.crwork3[i__3].r = 0., workspace_1.crwork3[i__3].i = 0.; /* initialize the radiation field */
	}

	if (iocom_1.nfin > 0) {
	    irec = islp;
	    if (inputcom_1.alignradf != 0) {
		irec = irec - *nslp + 1 + inputcom_1.offsetradf;
	    }
	    if (irec > 0) {
		//ierr = readfield_(workspace_1.crwork3, &irec); /* get field from file (record */
		ierr = readfield_srw(workspace_1.crwork3, &irec); //OC port

		if (ierr < 0) { //OC port
		    //last_();
		}
/* record nslp is loaded with */
	    } else {
		diagcom_1.pradoln[0] = inputcom_1.prad0;
		for (i__ = 2; i__ <= 7; ++i__) {
		    diagcom_1.pradoln[i__ - 1] = (float)0.;
		    if (i__ == inputcom_1.nharm) {
			diagcom_1.pradoln[i__ - 1] = inputcom_1.pradh0;
		    }
		}
		gauss_hermite__(workspace_1.crwork3, &inputcom_1.prad0, &inputcom_1.zrayl, &inputcom_1.zwaist, 
			&cartcom_1.xks, &cartcom_1.radphase, &c__1);
		if (inputcom_1.nharm > 1 && inputcom_1.pradh0 > 0.) {
		    d__1 = inputcom_1.zrayl * (f2c_doublereal) inputcom_1.nharm;
		    d__2 = cartcom_1.xks * (f2c_doublereal) inputcom_1.nharm;
		    gauss_hermite__(workspace_1.crwork3, &inputcom_1.pradh0, &d__1, &inputcom_1.zwaist, &d__2, 
				&cartcom_1.radphase, &inputcom_1.nharm);
/* lo */
		}
	    }
	} else {
	    i__2 = islp + 1 - *nslp;
	    dotimerad_(&i__2);
/* get time-dependence of sli */
	    gauss_hermite__(workspace_1.crwork3, &inputcom_1.prad0, &inputcom_1.zrayl, &inputcom_1.zwaist, 
			&cartcom_1.xks, &cartcom_1.radphase, &c__1);
	    if (inputcom_1.nharm > 1 && inputcom_1.pradh0 > 0.) {
		d__1 = inputcom_1.zrayl * (f2c_doublereal) inputcom_1.nharm;
		d__2 = cartcom_1.xks * (f2c_doublereal) inputcom_1.nharm;
		gauss_hermite__(workspace_1.crwork3, &inputcom_1.pradh0, &d__1, &inputcom_1.zwaist, &d__2, 
			&cartcom_1.radphase, &inputcom_1.nharm);
/* load */
	    }
	    diagcom_1.pradoln[0] = inputcom_1.prad0;
	}
	i__2 = *nslp - islp;
	pushtimerec_(workspace_1.crwork3, &inputcom_1.ncar, &i__2);
    }
    return 0;
}

//int srTSASE::loadslpfld_srw(f2c_integer* nslp)
//{
//    /* System generated locals */
//    f2c_integer i__1, i__2;
//
//    /* Local variables */
//	static f2c_integer irec, islp, ierr;
//    //extern /* Subroutine */ int last_();
//    //extern integer printerr_(), readfield_();
//    //extern /* Subroutine */ int dotimerad_(), gauss_hermite__(), pushtimerec_();
//
///*     ========================================================= */
///*     fills the array crtime with a seeding field for the first */
///*     slice. */
///*     --------------------------------------------------------- */
///*     error codes */
///* genesis version */
///* platform */
///* indicator for original fil */
///* indicator for sdds filetyp */
///* # of particles */
///* # of integration steps */
///* # of slices */
///* maximum of harmonics */
///* maximum of particle i */
///* <> 0 keeps distribution in */
///* energy units (mc^2) in ev */
///* vacuum impedence in ohms */
///* speed of light * electron */
///* pi */
///* pi/2 */
///* 2*pi */
///* check i for precission */
///* check ii for precission */
///* number of radial points fo */
///* # of gridpoints of cartesi */
///*     function prototypes */
///*     ------------------------------------------------------ */
///*     all input variables */
///*     wiggler */
///*     electron beam */
///*     radiation */
///*     grid-quantities */
///*     control */
///*     strong focusing */
///*     loading */
///*     output */
///*     external files */
///*     time-dependency */
///*     scan */
///*     extension */
///*     -------------------------------------------------------------------- */
///*     cartesian mesh */
///*     --------------------------------------------------------------------- */
///*     working-space (no cross-use of array values in call. subroutines) */
///*     ------------------------------------------------------------------ */
///*     input/output control */
///*     file:   timerec.com */
///*     time dependent common block - huge!!!! */
///*     ------------------------------------------------------------------ */
///*     -------------------------------------------------------------------- */
///*     slippage field */
///* must be > than ncar*ncar*nslp */
///* <>0 -> slippage is stored on disk */
//
//    if (inputcom_1.itdp == 0) {
//		return 0;
//    }
//
///*     check for limitation in timerecord */
///*     --------------------------------------------------------------- */
//	if (inputcom_1.nslice < *nslp * (1 - inputcom_1.iotail)) {
//	//	ierr = printerr_(&c_n20, "no output - ntail too small", (ftnlen)27);
//	}
//	//if (TRUE_) {
//	//if (*nslp * inputcom_1.ncar * inputcom_1.ncar > 15960700) {
//	//		ierr = printerr_(&c_n22, " ", (ftnlen)1);
//	//		last_();
//	//}
//	//}
//
///*     seding of the random number generator */
///*     ---------------------------------------------------------------- */
//    i__1 = *nslp - 1;
//	for (islp = 1; islp <= i__1; ++islp) {
//		if (iocom_1.nfin > 0) {
//			irec = islp;
//			if (inputcom_1.alignradf != 0) {
//				irec = irec - *nslp + 1 + inputcom_1.offsetradf;
//			}
//			if (irec > 0) {
//				//ierr = readfield_(workspace_1.crwork3, &irec);
//				ierr = readfield_srw(workspace_1.crwork3, &irec);
//				/* get field from file (record */
//				if (ierr < 0) {
//					//last_();
//				}
//				/* record nslp is loaded with */
//			} else {
//				gauss_hermite__(workspace_1.crwork3, &inputcom_1.prad0, &
//					inputcom_1.zrayl, &inputcom_1.zwaist, &cartcom_1.xks, 
//					&cartcom_1.radphase);
//			}
//		} else {
//			i__2 = 1 - islp;
//			dotimerad_(&i__2);
//			/* get time-dependence of slippage */
//			gauss_hermite__(workspace_1.crwork3, &inputcom_1.prad0, &
//				inputcom_1.zrayl, &inputcom_1.zwaist, &cartcom_1.xks, &
//				cartcom_1.radphase);
//		}
//		i__2 = *nslp - islp;
//		pushtimerec_(workspace_1.crwork3, &inputcom_1.ncar, &i__2);
//	}
//	return 0;
//} /* loadslpfld_srw */

//*************************************************************************

int srTSASE::loadrad_srw(f2c_integer* islice)
{
	static f2c_integer c__1 = 1;

    /* System generated locals */
    f2c_integer i__1, i__2;
    f2c_doublereal d__1, d__2;
    f2c_doublecomplex z__1, z__2;

    /* Builtin functions */
    //void d_cnjg();

    /* Local variables */
    static f2c_integer irec, ierr;
    //extern /* Subroutine */ int last_();
    static f2c_integer n;
    //extern integer readfield_();
    //extern /* Subroutine */ int gauss_hermite__();
    static f2c_integer ix;
    static f2c_doublereal pradin;

/*     ========================================================= */
/*     fills the array crfield with initial field */
/*     for start up from noise this should be small */
/*     --------------------------------------------------------- */
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
/*     --------------------------------------------------------------------- */
/*     working-space (no cross-use of array values in call. subroutines) */
/*     ------------------------------------------------------------------ */
/*     input/output control */
/*     ------------------------------------------------------------------ */
/*     diagnostic */
/*     -------------------------------------------------------------------- */
/*     time dependency of bunch + slippage field */
/*     simulation control and normalisation parameter */
/*     initialize field */
/*     --------------------------- */
    i__1 = inputcom_1.ncar * inputcom_1.ncar * cartcom_1.nhloop;
    for (ix = 1; ix <= i__1; ++ix) {
	i__2 = ix - 1;
	cartcom_1.crfield[i__2].r = 0., cartcom_1.crfield[i__2].i = 0.;
    }

    diagcom_1.pradoln[0] = inputcom_1.prad0;
/* first halfstep no gain (see diagno) */
    for (n = 2; n <= 7; ++n) {
	diagcom_1.pradoln[n - 1] = 0.;
/* kg */
	if (n == inputcom_1.nharm) {
	    diagcom_1.pradoln[n - 1] = inputcom_1.pradh0;
	}
    }

/*     load initial field */
/*     ------------------------------------------------------------------ */
    if (iocom_1.nfin <= 0) {
	gauss_hermite__(cartcom_1.crfield, &inputcom_1.prad0, &inputcom_1.zrayl, &inputcom_1.zwaist, 
		&cartcom_1.xks, &cartcom_1.radphase, &c__1);
/* load gauss-hermite mode for all harmonic */
	if (inputcom_1.nharm > 1 && inputcom_1.pradh0 > 0.) {
	    d__1 = inputcom_1.zrayl * (f2c_doublereal) inputcom_1.nharm;
	    d__2 = cartcom_1.xks * (f2c_doublereal) inputcom_1.nharm;
	    gauss_hermite__(cartcom_1.crfield, &inputcom_1.pradh0, &d__1, &inputcom_1.zwaist, 
			&d__2, &cartcom_1.radphase, &inputcom_1.nharm);
/* load gauss-hermite mod */
	}
    } else {
	irec = tbunchcom_1.nslp - 1 + *islice;
/* get record number */
	if (inputcom_1.alignradf != 0) {
	    irec = inputcom_1.offsetradf + *islice;
	}
/* add offset, when se */
	if (inputcom_1.itdp == 0) {
	    irec = 1;
	}
/* scan+ss -> use 1sr record */
	if (irec > 0) {
/* physical record? */
	    //ierr = readfield_(cartcom_1.crfield, &irec);
	    ierr = readfield_srw(cartcom_1.crfield, &irec); //OC port
/* get fundamental field from */
	    if (ierr < 0) { //OC to implement
		//last_();
	    }
/* stop if error occured */
	} else {
	    d__1 = inputcom_1.zrayl * (f2c_doublereal) inputcom_1.nharm;
	    d__2 = cartcom_1.xks * (f2c_doublereal) inputcom_1.nharm;
	    gauss_hermite__(cartcom_1.crfield, &inputcom_1.prad0, &d__1, &inputcom_1.zwaist, &d__2, 
			&cartcom_1.radphase, &c__1);
	    if (inputcom_1.nharm > 1 && inputcom_1.pradh0 > 0.) {
		d__1 = inputcom_1.zrayl * (f2c_doublereal) inputcom_1.nharm;
		d__2 = cartcom_1.xks * (f2c_doublereal) inputcom_1.nharm;
		gauss_hermite__(cartcom_1.crfield, &inputcom_1.pradh0, &d__1, &inputcom_1.zwaist, &d__2, 
			&cartcom_1.radphase, &inputcom_1.nharm);
/* load gauss-hermite m */
	    }
	}
	pradin = 0.;
	i__1 = inputcom_1.ncar * inputcom_1.ncar;
	for (ix = 1; ix <= i__1; ++ix) {
/* copy to crfield */
	    i__2 = ix - 1;
	    d_cnjg(&z__2, &cartcom_1.crfield[ix - 1]);
	    z__1.r = cartcom_1.crfield[i__2].r * z__2.r - cartcom_1.crfield[i__2].i * z__2.i;//, 
		z__1.i = cartcom_1.crfield[i__2].r * z__2.i + cartcom_1.crfield[i__2].i * z__2.r;
	    pradin += z__1.r;
	}
/* Computing 2nd power */
	d__1 = cartcom_1.dxy * 510999.06 * simcom_1.xkper0 / cartcom_1.xks;
	inputcom_1.prad0 = pradin * (d__1 * d__1) / 376.73;
	diagcom_1.pradoln[0] = inputcom_1.prad0;
    }
	return 0;
} /* loadrad_srw */

//int srTSASE::loadrad_srw(f2c_integer* islice)
//{
//	/* System generated locals */
//	f2c_integer i__1, i__2;
//	f2c_doublereal d__1;
//	f2c_doublecomplex z__1, z__2;
//
//	/* Builtin functions */
//	//void d_cnjg();
//
//	/* Local variables */
//	static f2c_integer irec, ierr;
//	//extern /* Subroutine */ int last_();
//	//extern integer readfield_();
//	//extern /* Subroutine */ int gauss_hermite__();
//	static f2c_integer ix;
//	static f2c_doublereal pradin;
//
//	/*     ========================================================= */
//	/*     fills the array crfield with initial field */
//	/*     for start up from noise this should be small */
//	/*     --------------------------------------------------------- */
//	/*     error codes */
//	/* genesis version */
//	/* platform */
//	/* indicator for original fil */
//	/* indicator for sdds filetyp */
//	/* # of particles */
//	/* # of integration steps */
//	/* # of slices */
//	/* maximum of harmonics */
//	/* maximum of particle i */
//	/* <> 0 keeps distribution in */
//	/* energy units (mc^2) in ev */
//	/* vacuum impedence in ohms */
//	/* speed of light * electron */
//	/* pi */
//	/* pi/2 */
//	/* 2*pi */
//	/* check i for precission */
//	/* check ii for precission */
//	/* number of radial points fo */
//	/* # of gridpoints of cartesi */
//	/*     function prototypes */
//	/*     ------------------------------------------------------ */
//	/*     cartesian mesh */
//	/*     ------------------------------------------------------ */
//	/*     all input variables */
//	/*     wiggler */
//	/*     electron beam */
//	/*     radiation */
//	/*     grid-quantities */
//	/*     control */
//	/*     strong focusing */
//	/*     loading */
//	/*     output */
//	/*     external files */
//	/*     time-dependency */
//	/*     scan */
//	/*     extension */
//	/*     ------------------------------------------------------ */
//	/*     working-space (no cross-use of array values in call. subroutines) */
//	/*     ------------------------------------------------------ */
//	/*     input/output control */
//	/*     ------------------------------------------------------ */
//	/*     diagnostic */
//	/*     ------------------------------------------------------ */
//	/*     time dependency of bunch + slippage field */
//	/*     simulation control and normalisation parameter */
//
//	diagcom_1.pradol = inputcom_1.prad0;
//
//	/*     load initial field */
//	/*     ------------------------------------------------------------------ */
//	/* first halfstep no gain (see diagno) */
//	if (iocom_1.nfin <= 0) {
//		gauss_hermite__(cartcom_1.crfield, &inputcom_1.prad0, &
//			inputcom_1.zrayl, &inputcom_1.zwaist, &cartcom_1.xks, &
//			cartcom_1.radphase);
//		/* loa */
//	} else {
//		irec = tbunchcom_1.nslp - 1 + *islice;
//		/* get record number */
//		if (inputcom_1.alignradf != 0) {
//			irec = inputcom_1.offsetradf + *islice;
//		}
//		/* add offset, when se */
//		if (inputcom_1.itdp == 0) {
//			irec = 1;
//		}
//		/* scan+ss -> use 1sr record */
//		if (irec > 0) {
//			/* physical record? */
//			//ierr = readfield_(cartcom_1.crfield, &irec);
//			ierr = readfield_srw(cartcom_1.crfield, &irec);
//			/* get field from file */
//			if (ierr < 0) {
//				//last_();
//			}
//			/* stop if error occured */
//		} else {
//			gauss_hermite__(cartcom_1.crfield, &inputcom_1.prad0, &
//				inputcom_1.zrayl, &inputcom_1.zwaist, &cartcom_1.xks, &
//				cartcom_1.radphase);
//		}
//		pradin = (float)0.;
//		i__1 = inputcom_1.ncar * inputcom_1.ncar;
//		for (ix = 1; ix <= i__1; ++ix) {
//			/* copy to crfield */
//			i__2 = ix - 1;
//			d_cnjg(&z__2, &cartcom_1.crfield[ix - 1]);
//			z__1.r = cartcom_1.crfield[i__2].r * z__2.r - cartcom_1.crfield[i__2].i * z__2.i;
//			z__1.i = cartcom_1.crfield[i__2].r * z__2.i + cartcom_1.crfield[i__2].i * z__2.r;
//			pradin += z__1.r;
//		}
//		/* Computing 2nd power */
//		d__1 = cartcom_1.dxy * 510999.06 * simcom_1.xkper0 / cartcom_1.xks;
//		inputcom_1.prad0 = pradin * (d__1 * d__1) / 376.73;
//		diagcom_1.pradol = inputcom_1.prad0;
//	}
//	return 0;
//} /* loadrad_srw */

//*************************************************************************

int srTSASE::readfield_srw(f2c_doublecomplex* cin, f2c_integer* irec)
{//to implement: read electric field data for a slice from SeedRad (3D srw rad. struct)
//the first slice is at *irec == 1

	if(*irec <= 0) return 0;
	if((SeedRad.PresT <= 0) || (!SeedRad.ElecFieldIsDefined())) return 0;

	double tStart = inputcom_1.xlamds*inputcom_1.zsep*inputcom_1.ntail/SpeedOfLight;
	double tStep = inputcom_1.xlamds*inputcom_1.zsep/SpeedOfLight;

	bool OutOfRangeT = false;

	double tMoment = tStart + (*irec - 1)*tStep;
	int it0_srw = (int)((tMoment - SeedRad.eStart)/SeedRad.eStep + 0.000001);
	int nt_srw_mi_1 = SeedRad.ne - 1;
	//if((it0_srw < 0) || (it0_srw >= nt_srw_mi_1)) return 0; //don't fill-in anything if out of range (- to check if it's OK)
	if((it0_srw < 0) || (it0_srw >= nt_srw_mi_1)) OutOfRangeT = true; //don't fill-in anything if out of range (- to check if it's OK)
	double rt = (tMoment - (SeedRad.eStart + it0_srw*SeedRad.eStep))/SeedRad.eStep; //0 <= rt <= 1
	int it1_srw = it0_srw + 1;

	double cos_phaseShiftDueToFreqDif = 1, sin_phaseShiftDueToFreqDif = 0;
	if(m_dwSeedRad != 0)
	{//add eventual linear phase shift due to difference of the two frequences
		double phaseShiftDueToFreqDif = m_dwSeedRad*tMoment;
		cos_phaseShiftDueToFreqDif = cos(phaseShiftDueToFreqDif);
		sin_phaseShiftDueToFreqDif = sin(phaseShiftDueToFreqDif);
	}

	double ElecFieldConvConst = TDElecFieldConvConstSRW2GENESIS(); //to check!!!
	const double inv_sqrt2 = 1./sqrt(2.);

	int Ncar = inputcom_1.ncar;
	bool IsPlanar = (inputcom_1.iwityp == 0);
	bool IsHelical = (inputcom_1.iwityp == 1);

	double Dmax = RadMeshRange();
	double xzStep = Dmax/double(Ncar - 1);
	double xzStart = -0.5*Dmax;
	int nz_srw_mi_1 = SeedRad.nz - 1, nx_srw_mi_1 = SeedRad.nx - 1;

	//long perT = 2;
	//long perX = perT*SeedRad.ne;
	//long perZ = perX*SeedRad.nx;
	long long perT = 2;
	long long perX = perT*SeedRad.ne;
	long long perZ = perX*SeedRad.nx;
	float *pEX0 = SeedRad.pBaseRadX;
	float *pEZ0 = SeedRad.pBaseRadZ;
	double ReEX_srw_interp = 0, ImEX_srw_interp = 0, ReEZ_srw_interp = 0, ImEZ_srw_interp = 0;

	//long perT_it0_srw = perT*it0_srw;
	//long perT_it1_srw = perT*it1_srw;
	long long perT_it0_srw = perT*it0_srw;
	long long perT_it1_srw = perT*it1_srw;

	double z_genesis = xzStart - xzStep;
	for(int iz=0; iz<inputcom_1.ncar; iz++)
	{
		z_genesis += xzStep;
		int iz0_srw = (int)((z_genesis - SeedRad.zStart)/SeedRad.zStep + 0.000001);
		if((iz0_srw < 0) || (iz0_srw >= nz_srw_mi_1)) continue; //don't fill-in anything if out of range (- to check if it's OK)
		double rz = (z_genesis - (SeedRad.zStart + iz0_srw*SeedRad.zStep))/SeedRad.zStep; //0 <= rz <= 1
		int iz1_srw = iz0_srw + 1;

		//long perZ_iz0_srw = perZ*iz0_srw;
		//long perZ_iz1_srw = perZ*iz1_srw;
		long long perZ_iz0_srw = perZ*iz0_srw;
		long long perZ_iz1_srw = perZ*iz1_srw;

		//f2c_doublecomplex *p_cin_p_1_p_Ncar_iz = cin + (1 + Ncar*iz);
		f2c_doublecomplex *p_cin_p_Ncar_iz = cin + (Ncar*iz);

		double x_genesis = xzStart - xzStep;
		for(int ix=0; ix<inputcom_1.ncar; ix++)
		{
			//f2c_doublecomplex *p_cin = p_cin_p_1_p_Ncar_iz + ix; //assuming cin to be 1-based flat array
			f2c_doublecomplex *p_cin = p_cin_p_Ncar_iz + ix; //assuming cin to be 1-based flat array

			if(OutOfRangeT)
			{
				p_cin->r = p_cin->i = 0;
			}
			else
			{
				x_genesis += xzStep;
				int ix0_srw = (int)((x_genesis - SeedRad.xStart)/SeedRad.xStep + 0.000001);
				if((ix0_srw < 0) || (ix0_srw >= nx_srw_mi_1)) continue; //don't fill-in anything if out of range (- to check if it's OK)
				double rx = (x_genesis - (SeedRad.xStart + ix0_srw*SeedRad.xStep))/SeedRad.xStep; //0 <= rz <= 1
				int ix1_srw = ix0_srw + 1;

				//long perX_ix0_srw = perX*ix0_srw;
				//long perX_ix1_srw = perX*ix1_srw;
				long long perX_ix0_srw = perX*ix0_srw;
				long long perX_ix1_srw = perX*ix1_srw;

				//long of000 = perZ_iz0_srw + perX_ix0_srw + perT_it0_srw;
				//long of100 = perZ_iz0_srw + perX_ix0_srw + perT_it1_srw;
				//long of010 = perZ_iz0_srw + perX_ix1_srw + perT_it0_srw;
				//long of001 = perZ_iz1_srw + perX_ix0_srw + perT_it0_srw;
				//long of110 = perZ_iz0_srw + perX_ix1_srw + perT_it1_srw;
				//long of101 = perZ_iz1_srw + perX_ix0_srw + perT_it1_srw;
				//long of011 = perZ_iz1_srw + perX_ix1_srw + perT_it0_srw;
				//long of111 = perZ_iz1_srw + perX_ix1_srw + perT_it1_srw;
				//long of000p = of000 + 1;
				//long of100p = of100 + 1;
				//long of010p = of010 + 1;
				//long of001p = of001 + 1;
				//long of110p = of110 + 1;
				//long of101p = of101 + 1;
				//long of011p = of011 + 1;
				//long of111p = of111 + 1;
				long long of000 = perZ_iz0_srw + perX_ix0_srw + perT_it0_srw;
				long long of100 = perZ_iz0_srw + perX_ix0_srw + perT_it1_srw;
				long long of010 = perZ_iz0_srw + perX_ix1_srw + perT_it0_srw;
				long long of001 = perZ_iz1_srw + perX_ix0_srw + perT_it0_srw;
				long long of110 = perZ_iz0_srw + perX_ix1_srw + perT_it1_srw;
				long long of101 = perZ_iz1_srw + perX_ix0_srw + perT_it1_srw;
				long long of011 = perZ_iz1_srw + perX_ix1_srw + perT_it0_srw;
				long long of111 = perZ_iz1_srw + perX_ix1_srw + perT_it1_srw;
				long long of000p = of000 + 1;
				long long of100p = of100 + 1;
				long long of010p = of010 + 1;
				long long of001p = of001 + 1;
				long long of110p = of110 + 1;
				long long of101p = of101 + 1;
				long long of011p = of011 + 1;
				long long of111p = of111 + 1;

				if(pEX0 != 0)
				{
					ReEX_srw_interp = ElecFieldConvConst*InterpBilin3D(rt, rx, rz, 
						*(pEX0 + of000), *(pEX0 + of100), *(pEX0 + of010), *(pEX0 + of001),
						*(pEX0 + of110), *(pEX0 + of101), *(pEX0 + of011), *(pEX0 + of111));
					ImEX_srw_interp = ElecFieldConvConst*InterpBilin3D(rt, rx, rz, 
						*(pEX0 + of000p), *(pEX0 + of100p), *(pEX0 + of010p), *(pEX0 + of001p),
						*(pEX0 + of110p), *(pEX0 + of101p), *(pEX0 + of011p), *(pEX0 + of111p));
					if(m_dwSeedRad != 0) 
					{//add eventual linear phase shift due to difference of the two frequences
						double newReEX = ReEX_srw_interp*cos_phaseShiftDueToFreqDif - ImEX_srw_interp*sin_phaseShiftDueToFreqDif;
						ImEX_srw_interp = ImEX_srw_interp*cos_phaseShiftDueToFreqDif + ReEX_srw_interp*sin_phaseShiftDueToFreqDif;
						ReEX_srw_interp = newReEX;
					}
				}
				if(pEZ0 != 0)
				{
					ReEZ_srw_interp = ElecFieldConvConst*InterpBilin3D(rt, rx, rz, 
						*(pEZ0 + of000), *(pEZ0 + of100), *(pEZ0 + of010), *(pEZ0 + of001),
						*(pEZ0 + of110), *(pEZ0 + of101), *(pEZ0 + of011), *(pEZ0 + of111));
					ImEZ_srw_interp = ElecFieldConvConst*InterpBilin3D(rt, rx, rz, 
						*(pEZ0 + of000p), *(pEZ0 + of100p), *(pEZ0 + of010p), *(pEZ0 + of001p),
						*(pEZ0 + of110p), *(pEZ0 + of101p), *(pEZ0 + of011p), *(pEZ0 + of111p));
					if(m_dwSeedRad != 0) 
					{//add eventual linear phase shift due to difference of the two frequences
						double newReEZ = ReEZ_srw_interp*cos_phaseShiftDueToFreqDif - ImEZ_srw_interp*sin_phaseShiftDueToFreqDif;
						ImEZ_srw_interp = ImEZ_srw_interp*cos_phaseShiftDueToFreqDif + ReEZ_srw_interp*sin_phaseShiftDueToFreqDif;
						ReEZ_srw_interp = newReEZ;
					}
				}

				if(IsPlanar)
				{
					if(pEX0 != 0)
					{
						p_cin->r = ReEX_srw_interp;
						p_cin->i = ImEX_srw_interp;
					}
					else if(pEZ0 != 0)
					{
						p_cin->r = ReEZ_srw_interp;
						p_cin->i = ImEZ_srw_interp;
					}
				}
				else if(IsHelical)
				{
					p_cin->r = inv_sqrt2*(ReEX_srw_interp - ImEZ_srw_interp); //or inv_sqrt2*(ReEX_srw_interp + ImEZ_srw_interp);
					p_cin->i = inv_sqrt2*(ImEX_srw_interp + ReEZ_srw_interp); //or inv_sqrt2*(ImEX_srw_interp - ReEZ_srw_interp);
				}
				else
				{
					p_cin->r = p_cin->i = 0.;
				}
			}
		}
	}
	return 0;

//    /* System generated locals */
//    integer ret_val, i__1, i__2, i__3;
//    doublecomplex z__1;
//
//    /* Builtin functions */
//    double sqrt();
//    integer s_rdue(), do_uio(), e_rdue();
//
//    /* Local variables */
//    extern integer printerr_();
//    static integer ix;
//    static doublereal scltmp;
//
//    /* Fortran I/O blocks */
//    static cilist io___45 = { 1, 0, 0, 0, 0 };
//
///*     ============================================== */
///*     read field from input file */
///*     ---------------------------------------------- */
///*     error codes */
///* genesis version */
///* platform */
///* indicator for original fil */
///* indicator for sdds filetyp */
///* # of particles */
///* # of integration steps */
///* # of slices */
///* maximum of harmonics */
///* maximum of particle i */
///* <> 0 keeps distribution in */
///* energy units (mc^2) in ev */
///* vacuum impedence in ohms */
///* speed of light * electron */
///* pi */
///* pi/2 */
///* 2*pi */
///* check i for precission */
///* check ii for precission */
///* number of radial points fo */
///* # of gridpoints of cartesi */
///*     function prototypes */
///*     ------------------------------------------------------------------ */
///*     input/output control */
///*     ------------------------------------------------------ */
///*     all input variables */
///*     wiggler */
///*     electron beam */
///*     radiation */
///*     grid-quantities */
///*     control */
///*     strong focusing */
///*     loading */
///*     output */
///*     external files */
///*     time-dependency */
///*     scan */
///*     extension */
///*     -------------------------------------------------------------------- */
///*     cartesian mesh */
///*     simulation control and normalisation parameter */
//    /* Parameter adjustments */
//    --cin;
//
//    /* Function Body */
//    scltmp = cartcom_1.xks * sqrt(376.73) / (cartcom_1.dxy * 510999.06 * 
//	    simcom_1.xkper0);
//
//    ret_val = 0;
//    io___45.ciunit = iocom_1.nfin;
//    io___45.cirec = *irec;
//    i__1 = s_rdue(&io___45);
//    if (i__1 != 0) {
//	goto L1;
//    }
//    i__2 = inputcom_1.ncar * inputcom_1.ncar;
//    for (ix = 1; ix <= i__2; ++ix) {
//	i__1 = do_uio(&c__2, (char *)&cin[ix], (ftnlen)sizeof(doublereal));
//	if (i__1 != 0) {
//	    goto L1;
//	}
//    }
//    i__1 = e_rdue();
//    if (i__1 != 0) {
//	goto L1;
//    }
//    i__1 = inputcom_1.ncar * inputcom_1.ncar;
//    for (ix = 1; ix <= i__1; ++ix) {
//	i__2 = ix;
//	i__3 = ix;
//	z__1.r = scltmp * cin[i__3].r, z__1.i = scltmp * cin[i__3].i;
//	cin[i__2].r = z__1.r, cin[i__2].i = z__1.i;
//    }
//
//    return ret_val;
//L1:
//    ret_val = printerr_(&c_n19, inputcom_1.fieldfile, (ftnlen)30);
//    return ret_val;
} /* readfield_ */

//*************************************************************************

int srTSASE::EstimateGenesisNSLP()
{
	int tbunchcom_1_nsep = (int) (inputcom_1.zsep / inputcom_1.delz);
	int tbunchcom_1_nslp = wigcom_1.nstepz / tbunchcom_1_nsep;
    if(wigcom_1.nstepz % tbunchcom_1.nsep != 0) ++tbunchcom_1_nslp;
	return tbunchcom_1_nslp;
}

//*************************************************************************

int srTSASE::GenesisCompDrive(srTSRWRadStructAccessData *arRadAccessData, int numHarm)
{
	int result = 0;

    /* System generated locals */
    f2c_integer i__1, i__2, i__2n, i__3;

    /* Builtin functions */
    /* Subroutine */ //int s_stop();

    /* Local variables */
    static f2c_integer islp, isep;
    //extern /* Subroutine */ int loadbeam_(), chk_loss__(), last_(), stepz_(), swapfield_(), doscan_();
    static f2c_integer islice, lstepz, istepz;
    //extern /* Subroutine */ int initio_(), dotime_(), output_(), loadrad_(), outglob_(), initrun_(), outhist_(), outdump_();
	static f2c_integer c__1 = 1;

	m_ElecDistribShouldBeUsed = CheckIfElecDistrShouldBeUsed(); //should be executed before anything

	//emulating MPI (wrappers)
    mpi_init__(&mpicom_1.mpi_err__);
    mpi_comm_rank__(&c__1, &mpicom_1.mpi_id__, &mpicom_1.mpi_err__);
    mpi_comm_size__(&c__1, &mpicom_1.mpi_size__, &mpicom_1.mpi_err__);

    //if(result = initio_()) return result; /* open files, read in input */
	if(result = readin_()) return result; //OC_port: readin_() was modified (reading from files removed) - to re-implement reading from files !
	//input data from files is not done !!!
	if(result = ConvertInputDataToGenesisFormat(numHarm)) return result;
	//if(result = Alloc_tslipcom(PrecDat.ncar, PrecDat.nslice)) return result; //separate allocation of GENESIS crtime buffer

	if(result = chk_input__()) return result;
	
	if(result = CheckInputConsistency()) return result; /* Check boundaries of input variables */
	//if(result = chk_bnd__()) return result;
	
	if(result = initrun_srw()) return result; /* initialize simulation */
	if(result = outglob_()) return result;
	
	//--------------SRW--------------
	//if(result = PrepSRWRadStructTD(RadAccessData)) return result;
	if(result = PrepSRWRadStructTD(arRadAccessData, numHarm)) return result;

	if(result = SetupOutputControlStruct(numHarm)) return result;
	long PassCount = 0; //for Progress Indicator
	//--------------SRW--------------

	//start loop for each slice
	//output of global parameter (t-independent)
	i__1 = inputcom_1.nslice;
    i__2n = mpicom_1.mpi_size__;

    //for(islice = 1; islice <= i__1; ++islice) 
    //{
    for (islice = mpicom_1.mpi_id__ + 1; i__2n < 0 ? islice >= i__1 : islice <= i__1; islice += i__2n) 
	{
		/* loop for each slice */
		mpicom_1.mpi_loop__ = mpicom_1.mpi_size__;
		if(islice - mpicom_1.mpi_id__ >= inputcom_1.nslice - mpicom_1.mpi_size__ + 1) 
		{
			mpicom_1.mpi_loop__ = inputcom_1.nslice % mpicom_1.mpi_size__;
		}
		if(mpicom_1.mpi_loop__ == 0) 
		{
			mpicom_1.mpi_loop__ = mpicom_1.mpi_size__;
		}

		//initial loading
		//loop for each slice
		istepz = 0;
	
		if(result = doscan_(&islice)) return result; //update scan value
		if(result = dotime_(&islice)) return result; //calculate time-dependent parameters

		if(result = loadrad_srw(&islice)) return result; //radiation field loading
		//modified to allow using radiation field from SRW

		if(result = loadbeam_srw(&islice, &simcom_1.xkw0)) return result; //particle loading
		//modified to allow using electron distribution from SRW

		if(result = chk_loss__()) return result; //remove cut particle
		//--------------SRW--------------
		if(result = DiagnoGenesisPlus(PassCount++, istepz, islice, numHarm)) return result;
		//--------------SRW--------------

		//if(result = openoutputbinmpi_(&islice)) return result; //open binary outputfile for mpi fe
		//if(result = output_(&istepz, &islice, &simcom_1.xkw0)) return result;

		//propagate beam for nu wiggler periods
		i__2 = tbunchcom_1.nslp;

		for(islp = 1; islp <= i__2; ++islp) //loop over slippage (advance field)
		{
			lstepz = tbunchcom_1.nsep;
			if(islp == tbunchcom_1.nslp) 
			{
				lstepz = wigcom_1.nstepz - (islp - 1) * tbunchcom_1.nsep; //correct for last section
			}
			
			i__3 = lstepz;

			for(isep = 1; isep <= i__3; ++isep) //loop 'steady state' simulation
			{
				++istepz;
				if(result = stepz_(&istepz, &simcom_1.xkw0)) return result; //advance one step in z
				
				//--------------SRW--------------
				if(result = DiagnoGenesisPlus(PassCount++, istepz, islice, numHarm)) return result; //debug
				//--------------SRW--------------

				//if(result = output_(&istepz, &islice, &simcom_1.xkw0)) return result;
				//if(result = OutResDistrib((long)istepz, (long)islice, (double)simcom_1.xkw0)) return result;
			}

					//DEBUG
					//char ErrorMesTitle[] = "SRW Debug";
					//char ErrorStr[100];
					//int j = sprintf(ErrorStr, "inputcom_1.itdp: %d", inputcom_1.itdp);
					////j += sprintf(ErrorStr + j, "          Nln: %d", Nln);
					//UINT DlgStyle = MB_OK | MB_ICONSTOP | MB_DEFBUTTON1 | MB_SYSTEMMODAL;
					//int MesBoxInf = MessageBox(NULL, ErrorStr, ErrorMesTitle, DlgStyle); 
					//END DEBUG

			if((inputcom_1.itdp != 0) && (islp < tbunchcom_1.nslp)) 
			{
				if(result = swapfield_(&islp)) return result; //advance field in time dep. simulation
			}
		}

		//outhist_(&islice);
		//outdump_(&islice); //dump rad + part if needed
		//closeoutputbinmpi_(); //close binary files for m
		if(PrecDat.CreateElecDistr)
		{
			if(result = OutDumpResDistrib((int)islice)) return result;
		}
		if(result = CopyRadSliceToSRWRadStructTD((int)islice, arRadAccessData, numHarm)) return result;
	}

	//merge has to be done first so that the field dump is already
	//created and the records from the slippage field are written at the right position
	mpi_merge__(); 
	//outdumpslippage_(); //merge single outputfiles into one

	//last_();
	///* end run */
	//s_stop("", (ftnlen)0);

////OC-TMP	
//	if(result = preset_()) return result; /* Set default values to the input parameters */
//
//	if(result = ConvertInputDataToGenesisFormat()) return result;
//	//if(result = readin_()) return result; /* I/O action */
//
//	if(result = CheckInputConsistency()) return result; /* Check boundaries of input variables */
//	//if(result = chk_bnd__()) return result; /* Check boundaries of input variables */
//
//	//if(result = auxval_()) return result; /* Auxilliary values (Normalization, mesh, etc...) */
//	if(result = AuxvalGenesisPlus()) return result; /* Auxilliary values (Normalization, mesh, etc...): fixing probable bugs */
//
//	//if(result = outglob_()) return result; /* output of global parameter (t-independent) */
//
// 	//--------------SRW--------------
//	if(result = SetupOutputControlStruct()) return result;
//	long PassCount = 0; //for Progress Indicator
//	//--------------SRW--------------
//
//	long i__1 = timecom_1.nslice;
//    for(timecom_1.islice = 1; timecom_1.islice <= i__1; ++timecom_1.islice) // Loop for each slice
//	{
//		//Initial loading
////OC-TMP		
//		doscan_();
////OC-TMP		
//		loading_(); // Load particles and beam
//		
//		//diagno_(); // Diagnostics
//		//--------------SRW--------------
//		if(result = DiagnoGenesisPlus(PassCount++)) return result;
//		//--------------SRW--------------
//
//		//outfield_(); // output of field // change it: write to a wave, if necessary
//		//outpart_(); // output of particles // change it: write to a wave, if necessary
//		//status_(); // Tell user how much of the program has run // ??
//
//		// Propagate beam for NU wiggler periods
//		long i__2 = timecom_1.nslp - 1;
//		for(long islp = 0; islp <= i__2; ++islp) // Loop over slippage (advance field)
//		{
//			long i__3 = tbunchcom_1.nsep;
//			for(long isep = 1; isep <= i__3; ++isep) // Loop 'steady state' simulation 
//			{
//				simcom_1.istepz = isep + islp * tbunchcom_1.nsep;
////OC-TMP				
//				stepz_(); // Advance one step in Z
//		
//				//diagno_(); // Diagnostic of Beam and Radiation
//				//--------------SRW--------------
//				if(result = DiagnoGenesisPlus(PassCount++)) return result;
//				//--------------SRW--------------
//
//				//outfield_(); // Dump field when selected! // change it: write to a wave, if necessary
//				//outpart_(); // dump particles when selected // change it: write to a wave, if necessary
//				//status_(); // Output how much is done // ??
//			}
//
//			// Advance field after NSEP periods if time depending mode is selected
//			if(timecom_1.itdp != 0) 
//			{
//				i__3 = cartcom_1.ncar * cartcom_1.ncar;
//				for(long it = 1; it <= i__3; ++it) 
//				{
//					long i__4 = it - 1;
//
//					f2c_doublecomplex ctmp;
//					ctmp.r = cartcom_1.crfield[i__4].r; 
//					ctmp.i = cartcom_1.crfield[i__4].i;
//		    
//					i__4 = it - 1;
//					long i__5 = islp * cartcom_1.ncar * cartcom_1.ncar + it - 1;
//		    
//					f2c_doublecomplex z__1;
//					z__1.r = timecom_1.crtime[i__5].r; z__1.i = timecom_1.crtime[i__5].i;
//					cartcom_1.crfield[i__4].r = z__1.r; cartcom_1.crfield[i__4].i = z__1.i;
//
//					i__4 = islp * cartcom_1.ncar * cartcom_1.ncar + it - 1;
//					timecom_1.crtime[i__4].r = (f2c_real)ctmp.r; timecom_1.crtime[i__4].i = (f2c_real)ctmp.i;
//				}
//			}
//		}// ISLP
//		//outhist_(); // Final output // change it: write to a wave, if necessary
//    }// ISLICE
//
//	//last_(); // Closes all open files // change it: close waves, if necessary

	return 0;
}

//*************************************************************************

int srTSASE::FillInSRWRadStruct(srTSRWRadStructAccessData& RadAccessData)
{// This does not take the outer right and top points, since Genesis required odd nx and nz
	int result = 0;
	const double InvSqRt2 = 1./sqrt(2.);

	long Ncar = inputcom_1.ncar;

	bool IsPlanar = (inputcom_1.iwityp == 0);
	bool IsHelical = (inputcom_1.iwityp == 1);
	
	RadAccessData.ne = 1;
	RadAccessData.nx = RadAccessData.nz = Ncar - 1;

	double Dmax = RadMeshRange();
	double MeshStep = Dmax/double(Ncar - 1);
	RadAccessData.xStart = -0.5*Dmax; // check what offset in transverse pos. does here
	RadAccessData.zStart = -0.5*Dmax;

	RadAccessData.xStep = RadAccessData.zStep = MeshStep;

	//if(result = pSend->ModifyRadNeNxNz(RadAccessData)) return result;
	if(result = RadAccessData.ModifyWfrNeNxNz()) return result;

	//double GenesisSCLTMP = cartcom_1.dxy*EEV*simcom_1.xkper0/simcom_1.xks/sqrt(VACIMP);
	//double MultForE = GenesisSCLTMP; // this is in compliance with GENESIS
	double MultForE = RadFieldMultip(); // this is in compliance with SRW

	float *pEx0 = RadAccessData.pBaseRadX, *pEz0 = RadAccessData.pBaseRadZ;
	//long PerX = 2;
	//long PerZ = PerX*RadAccessData.nx;
	//long PerZ_CR = Ncar;
	long long PerX = 2;
	long long PerZ = PerX*RadAccessData.nx;
	long long PerZ_CR = Ncar;

	f2c_doublecomplex *pCRFIELD_0 = cartcom_1.crfield;

	for(long iz=0; iz<(Ncar - 1); iz++)
	{
		//long izPerZ = iz*PerZ;
		//long izPerZ_CR = iz*PerZ_CR;
		long long izPerZ = iz*PerZ;
		long long izPerZ_CR = iz*PerZ_CR;

		for(long ix=0; ix<(Ncar - 1); ix++)
		{//skip internal loop over e since only 1 energy slice is assumed here
			//long ixPerX = ix*PerX;
			long long ixPerX = ix*PerX;

			//long Offset = izPerZ + ixPerX;
			long long Offset = izPerZ + ixPerX;
			float *pEx = pEx0 + Offset, *pEz = pEz0 + Offset;

			//long OffsetCR = izPerZ_CR + ix;
			long long OffsetCR = izPerZ_CR + ix;
			f2c_doublecomplex *pCRFIELD = pCRFIELD_0 + OffsetCR;

			if(IsPlanar) // fill horizontal
			{
				*pEx = (float)(MultForE*(pCRFIELD->r)); *(pEx + 1) = (float)(MultForE*(pCRFIELD->i));
				*pEz = *(pEz + 1) = 0.;
			}
			else if(IsHelical) // fill circular-right component into Ex, Ez
			{
				double Ecrr = MultForE*(pCRFIELD->r), Ecri = MultForE*(pCRFIELD->i);
				*pEx = (float)(InvSqRt2*Ecrr); *(pEx + 1) = (float)(InvSqRt2*Ecri); //??
				*pEz = (float)(InvSqRt2*Ecri); *(pEz + 1) = (float)(-InvSqRt2*Ecrr); //??
			}
			else
			{
				*pEx = *(pEx + 1) = 0.;
				*pEz = *(pEz + 1) = 0.;
			}
		}
	}
	return result;
}

//*************************************************************************

int srTSASE::CopyRadSliceToSRWRadStructTD(int iSlice, srTSRWRadStructAccessData* arRadAccessData, int numHarm)
//int srTSASE::CopyRadSliceToSRWRadStructTD(int iSlice, srTSRWRadStructAccessData& RadAccessData)
{//This does not take the outer right and top points, since Genesis required odd nx and nz
 //iSlice is assumed 1-based !
	int result = 0;
	const double InvSqRt2 = 1./sqrt(2.);

	long Ncar = inputcom_1.ncar;
	bool IsPlanar = (inputcom_1.iwityp == 0);
	bool IsHelical = (inputcom_1.iwityp == 1);

	//double GenesisSCLTMP = cartcom_1.dxy*EEV*simcom_1.xkper0/simcom_1.xks/sqrt(VACIMP);
	//double MultForE = GenesisSCLTMP; // this is in compliance with GENESIS
	//To modify !!!
	//double MultForE = RadFieldMultip(); // this is in compliance with SRW
	double MultForE = 1./TDElecFieldConvConstSRW2GENESIS(); // this will give TD Electric field in sqrt(W/mm^2) in SRW

	//long PerZ_CR = Ncar;
	//long PerHarm_CR = PerZ_CR*Ncar;
	long long PerZ_CR = Ncar;
	long long PerHarm_CR = PerZ_CR*Ncar;

	//long PerT = 2;
	//long PerX = PerT*(arRadAccessData->ne); //same for all harmonics
	//long PerZ = PerX*(arRadAccessData->nx);
	long long PerT = 2;
	long long PerX = PerT*(arRadAccessData->ne); //same for all harmonics
	long long PerZ = PerX*(arRadAccessData->nx);
	int it = iSlice - 1;
	//long itPerT = it*PerT;
	long long itPerT = it*PerT;

	f2c_doublecomplex *pCRFIELD_0 = cartcom_1.crfield;

	srTSRWRadStructAccessData *t_arRad = arRadAccessData;

	for(int iHarm=0; iHarm<numHarm; iHarm++)
	{
		float *pEx0 = t_arRad->pBaseRadX, *pEz0 = t_arRad->pBaseRadZ;
		//long iHarmPerHarm_CR = iHarm*PerHarm_CR;
		long long iHarmPerHarm_CR = iHarm*PerHarm_CR;

		for(long iz=0; iz<(Ncar - 1); iz++)
		{
			//long izPerZ_p_itPerT = iz*PerZ + itPerT;
			//long izPerZ_CR = iz*PerZ_CR;
			//long iHarmPerHarm_p_izPerZ_CR = iHarmPerHarm_CR + izPerZ_CR;
			long long izPerZ_p_itPerT = iz*PerZ + itPerT;
			long long izPerZ_CR = iz*PerZ_CR;
			long long iHarmPerHarm_p_izPerZ_CR = iHarmPerHarm_CR + izPerZ_CR;

			for(long ix=0; ix<(Ncar - 1); ix++)
			{//skip internal loop over e since only 1 energy slice is assumed here
				//long ixPerX = ix*PerX;
				//long Offset = izPerZ_p_itPerT + ixPerX;
				long long ixPerX = ix*PerX;
				long long Offset = izPerZ_p_itPerT + ixPerX;
				float *pEx = pEx0 + Offset, *pEz = pEz0 + Offset;

				//long OffsetCR = izPerZ_CR + ix;
				//long OffsetCR = iHarmPerHarm_p_izPerZ_CR + ix;
				long long OffsetCR = iHarmPerHarm_p_izPerZ_CR + ix;
				f2c_doublecomplex *pCRFIELD = pCRFIELD_0 + OffsetCR;

				if(IsPlanar) // fill horizontal
				{
					*pEx = (float)(MultForE*(pCRFIELD->r)); *(pEx + 1) = (float)(MultForE*(pCRFIELD->i));
					*pEz = *(pEz + 1) = 0.;
				}
				else if(IsHelical) // fill circular-right component into Ex, Ez
				{
					double Ecrr = MultForE*(pCRFIELD->r), Ecri = MultForE*(pCRFIELD->i);
					*pEx = (float)(InvSqRt2*Ecrr); *(pEx + 1) = (float)(InvSqRt2*Ecri); //??
					*pEz = (float)(InvSqRt2*Ecri); *(pEz + 1) = (float)(-InvSqRt2*Ecrr); //??
				}
				else
				{
					*pEx = *(pEx + 1) = 0.;
					*pEz = *(pEz + 1) = 0.;
				}
			}
		}
		t_arRad++;
	}
	return result;
}

//*************************************************************************

int srTSASE::PrepSRWRadStructTD(srTSRWRadStructAccessData* arRadAccessData, int numHarm)
//int srTSASE::PrepSRWRadStructTD(srTSRWRadStructAccessData& RadAccessData)
{// This does not take the outer right and top points, since Genesis required odd nx and nz
	int result = 0;
	//const double InvSqRt2 = 1./sqrt(2.);
			//ControlSASE.nt = inputcom_1.nslice;
			//ControlSASE.tStep = inputcom_1.xlamds*inputcom_1.zsep/SpeedOfLight;

	long Ncar = inputcom_1.ncar;

	//bool IsPlanar = (inputcom_1.iwityp == 0);
	//bool IsHelical = (inputcom_1.iwityp == 1);

	double Dmax = RadMeshRange();
	double MeshStep = Dmax/double(Ncar - 1);
	double tStep = inputcom_1.xlamds*inputcom_1.zsep/SpeedOfLight;
	//double tRange = (inputcom_1.nslice - 1)*tStep;
	double tStart = inputcom_1.xlamds*inputcom_1.zsep*inputcom_1.ntail/SpeedOfLight;
	
	srTSRWRadStructAccessData *t_arRad = arRadAccessData;
	for(int iHarm=0; iHarm<numHarm; iHarm++)
	{
		t_arRad->ne = inputcom_1.nslice; //actually, it's nt
		if(t_arRad->ne < 1) t_arRad->ne = 1;
		t_arRad->nx = t_arRad->nz = Ncar - 1;

		t_arRad->xStart = -0.5*Dmax; // check what offset in transverse pos. does here
		t_arRad->zStart = -0.5*Dmax;
		t_arRad->xStep = t_arRad->zStep = MeshStep;

		if((inputcom_1.nslice > 1) && (inputcom_1.itdp != 0))
		{
			t_arRad->eStep = tStep;
			t_arRad->eStart = tStart;
			t_arRad->PresT = 1; //to ensure time-domain pres.
			t_arRad->avgPhotEn = PrecDat.photEn_xlamds*(iHarm + 1); //DistrInfoDat.LambStart;
			//inputcom_1.xlamds1.239842E-06/ = DistrInfoDat.LambStart;
		}
		else
		{
			t_arRad->eStart = PrecDat.photEn_xlamds*(iHarm + 1); //check units [eV]
			t_arRad->eStep = 0;
			t_arRad->PresT = 0; //to ensure frequency-domain pres.
		}

		if(result = t_arRad->ModifyWfrNeNxNz()) return result;

		//double GenesisSCLTMP = cartcom_1.dxy*EEV*simcom_1.xkper0/simcom_1.xks/sqrt(VACIMP);
		//double MultForE = GenesisSCLTMP; // this is in compliance with GENESIS
		//double MultForE = RadFieldMultip(); // this is in compliance with SRW

		t_arRad->SetNonZeroWavefrontLimitsToFullRange();
		t_arRad++;
	}
	return result;
}

//*************************************************************************

int srTSASE::InitAuxWfrParams(srTSRWRadStructAccessData *arSRWRadAccessData, int numHarm)
{
	int result = 0;
	if((arSRWRadAccessData == 0) || (numHarm <= 0)) return result;

	double OneUndSegmLength = inputcom_1.nwig*inputcom_1.xlamd;
	double Rinit = 0.2*OneUndSegmLength; //assuming saturation happens there
	double RinitErr = 1.0*Rinit; //to steer

	srTMomentsPtrs MomPtrsX, MomPtrsZ;
	srTSRWRadStructAccessData *t_arSRWRadAccessData = arSRWRadAccessData;
	for(int i=0; i<numHarm; i++)
	{
		t_arSRWRadAccessData->InitialSetupOf4x4PropMatr(Rinit);
		t_arSRWRadAccessData->RobsX = t_arSRWRadAccessData->RobsZ = Rinit;
		t_arSRWRadAccessData->RobsXAbsErr = t_arSRWRadAccessData->RobsZAbsErr = RinitErr; // To steer

		t_arSRWRadAccessData->ComputeRadMoments();
		t_arSRWRadAccessData->SetupRadMomentsPtrs(MomPtrsX, MomPtrsZ);
		t_arSRWRadAccessData->xc = *(MomPtrsX.pX);
		t_arSRWRadAccessData->zc = *(MomPtrsX.pZ);
		t_arSRWRadAccessData->UnderSamplingX = t_arSRWRadAccessData->UnderSamplingZ = 1;
		t_arSRWRadAccessData->SetupNonZeroWavefrontLimitsAtCreation();

		t_arSRWRadAccessData++;
	}

	return result;
}

//*************************************************************************

int srTSASE::PropagateWavefrontToObservationPlane(srTSRWRadStructAccessData *arSRWRadAccessData, int numHarm)
//int srTSASE::PropagateWavefrontToObservationPlane(srTSRWRadStructAccessData& RadAccessData)
{
	int result = 0;

	//to implement eventual change of presentation and propagation for Steady State and TD modes !!!
	//
//If some "Observation" parameters are defined (in frequency domain !), then
//- in case of steady-state simulations: straightforward propagation to the observation plane, and then resizing are performed;
//- in case of TD simulations: transformation to frequency domain is performed, and then, for each frequency, propagation to the observation plane, and resizing (if necessary);

/**
	double OneUndSegmLength = inputcom_1.nwig*inputcom_1.xlamd;
	double Rinit = 0.25*OneUndSegmLength; //assuming saturation happens there
	RadAccessData.InitialSetupOf4x4PropMatr(Rinit);
	RadAccessData.RobsX = RadAccessData.RobsZ = Rinit;
	RadAccessData.RobsXAbsErr = RadAccessData.RobsZAbsErr = 0.01*Rinit; // To steer

	//long MaxStepNo = simcom_1.istepz;

	//RadAccessData.xc = diagcom_1.xpos[MaxStepNo];
	//RadAccessData.zc = diagcom_1.ypos[MaxStepNo];
	//RadAccessData.UnderSamplingX = RadAccessData.UnderSamplingZ = 1;
	//RadAccessData.SetupNonZeroWavefrontLimitsAtCreation();

	//srTGenOptElem AuxOptElem;
	//if(result = AuxOptElem.ComputeRadMoments(&RadAccessData)) return result;
	//RadAccessData.RobsX = RadAccessData.RobsZ = 0.;

	//double EndLongPosOfTheUndulator = 0.;
	//double DriftLength = DistrInfoDat.yStart - EndLongPosOfTheUndulator;
	//bool ThereIsDrift = (::fabs(DriftLength) > 1.E-06);
	//if(!ThereIsDrift) // no propagation needed
	//{
	//	if(DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat) 
	//	{
	//		if(result = ResizeForOtimalSamplingInSpot(RadAccessData)) return result;
	//	}
	//	else
	//	{
	//		if(result = ResizeForSpecifiedSampling(RadAccessData)) return result;
	//	}
	//}
	//else
	//{
	//	if(result = ResizeAndPropagateFromWaist(RadAccessData)) return result;

	//	//dddddddddddddddddd
	//	//propagate using a special method
	//}
	////DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat = 0;
**/
	return result;
}

//*************************************************************************

int srTSASE::ResizeForOtimalSamplingInSpot(srTSRWRadStructAccessData& RadAccessData)
{
	int result = 0;

	double xRangeCur = RadAccessData.xStep*RadAccessData.nx;
	double xRangeFin = DistrInfoDat.xEnd - DistrInfoDat.xStart;
	double Pxm = xRangeFin/xRangeCur;
	double zRangeCur = RadAccessData.zStep*RadAccessData.nz;
	double zRangeFin = DistrInfoDat.zEnd - DistrInfoDat.zStart;
	double Pzm = zRangeFin/zRangeCur;

	const int NominalNumberOfPointsInSpot = 13; //to steer
	const double AcceptableOffsetFromNominalResolution = 0.2; //to steer
	const double RelZeroTolForIntens = 0.1;
	//the function sets a "nominal" resolution in the spot:
	//- if original resolution is better or worse than this for more than AcceptableOffsetFromNominalResolution, it is modified;

	srTGenOptElem AuxOptElem;
	srTRadSect1D Sect1D[2];
	if(result = AuxOptElem.SetupCharacteristicSections1D(&RadAccessData, Sect1D)) return result;

	srTMomentsPtrs MomX(RadAccessData.pMomX), MomZ(RadAccessData.pMomZ);
	char TreatExOrEz = (*(MomX.pTotPhot) > *(MomZ.pTotPhot))? 'x' : 'z';

	double AuxPd[] = {1., 1.};
	for(int i=0; i<2; i++)
	{
		long iFirst = -1, iLast = -1;
		AuxOptElem.FindIntensityBorders1D(Sect1D[i], TreatExOrEz, RelZeroTolForIntens, iFirst, iLast);
		if((iFirst > 0) && (iLast < Sect1D[i].np - 1) && (iFirst < iLast))
		{
			double dNpInSpot = double(iLast - iFirst);
			double TestPd = double(NominalNumberOfPointsInSpot)/dNpInSpot;
			if((TestPd < (1 - AcceptableOffsetFromNominalResolution)) ||
			   (TestPd > (1 + AcceptableOffsetFromNominalResolution)))
			{
				AuxPd[i] = TestPd;
			}
		}
	}

	srTRadResize RadRes;
	RadRes.pxd = AuxPd[0];
	RadRes.pzd = AuxPd[1];
	//RadRes.UseOtherSideFFT = 1;
	RadRes.useOtherSideFFT(1); //OC090311
	//RadRes.DoNotTreatSpherTerm = 1;
	RadRes.doNotTreatSpherTerm(1); //OC090311
	if(result = AuxOptElem.RadResizeGen(RadAccessData, RadRes)) return result;

	RadRes.pxd = RadRes.pzd = 1.;
	RadRes.pxm = Pxm;
	RadRes.pzm = Pzm;
	//RadRes.UseOtherSideFFT = 0;
	RadRes.useOtherSideFFT(0); //OC090311
	if(result = AuxOptElem.RadResizeGen(RadAccessData, RadRes)) return result;

	return result;
}

//*************************************************************************

int srTSASE::ResizeForSpecifiedSampling(srTSRWRadStructAccessData& RadAccessData)
{
	int result = 0;
	srTGenOptElem AuxOptElem;
	srTRadResize RadRes;

	double xRangeCur = RadAccessData.xStep*RadAccessData.nx;
	double xRangeFin = DistrInfoDat.xEnd - DistrInfoDat.xStart;
	double Pxm = xRangeFin/xRangeCur;
	double zRangeCur = RadAccessData.zStep*RadAccessData.nz;
	double zRangeFin = DistrInfoDat.zEnd - DistrInfoDat.zStart;
	double Pzm = zRangeFin/zRangeCur;

	RadRes.pxd = double(DistrInfoDat.nx)/double(RadAccessData.nx)/Pxm;
	RadRes.pzd = double(DistrInfoDat.nz)/double(RadAccessData.nz)/Pzm;
	//RadRes.UseOtherSideFFT = 1;
	RadRes.useOtherSideFFT(1); //OC090311
	//RadRes.DoNotTreatSpherTerm = 1;
	RadRes.doNotTreatSpherTerm(1); //OC090311
	if(result = AuxOptElem.RadResizeGen(RadAccessData, RadRes)) return result;

	RadRes.pxd = RadRes.pzd = 1.;

	RadRes.pxm = Pxm;
	RadRes.pzm = Pzm;
	//RadRes.UseOtherSideFFT = 0;
	RadRes.useOtherSideFFT(0); //OC090311
	if(result = AuxOptElem.RadResizeGen(RadAccessData, RadRes)) return result;
	return result;
}

//*************************************************************************

int srTSASE::ResizeAndPropagateFromWaist(srTSRWRadStructAccessData& RadAccessData)
{//Determines resize coefficients to have necessary range after propagation from waist
	int result = 0;
	const double DontResizeThresh = 0.15;
  	//const int SmallestN = 8;
	//const int LargestN = 1024; // to steer

	double pxd, pzd, pxm, pzm;

	double EndLongPosOfTheUndulator = 0.;
	double DriftLength = DistrInfoDat.yStart - EndLongPosOfTheUndulator;
	//double Wavelength_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 1.239854E-06/DistrInfoDat.LambStart : 1.E-06*DistrInfoDat.LambEnd;
	double Wavelength_m = (DistrInfoDat.TreatLambdaAsEnergyIn_eV)? 1.239842E-06/DistrInfoDat.LambStart : 1.E-06*DistrInfoDat.LambEnd;

//Resizing param. for current resolution -> future range
	double xRangeReq = DistrInfoDat.xEnd - DistrInfoDat.xStart;
	double xStepReq = DriftLength*Wavelength_m/xRangeReq;
	pxd = RadAccessData.xStep/xStepReq;
	if(::fabs(pxd - 1.) <= DontResizeThresh) pxd = 1.;

	double zRangeReq = DistrInfoDat.zEnd - DistrInfoDat.zStart;
	double zStepReq = DriftLength*Wavelength_m/zRangeReq;
	pzd = RadAccessData.zStep/zStepReq;
	if(::fabs(pzd - 1.) <= DontResizeThresh) pzd = 1.;

//Resizing param. for current range -> future resolution
	//double UnderSamplingX = 1., UnderSamplingZ = 1.;
	double HalfLambR = 0.5*Wavelength_m*DriftLength;
	double dx, dz;
	if(DistrInfoDat.AllowAutoChoiceOfNxNzForPropagat)
	{
		double MultFactor = DistrInfoDat.NxNzOversamplingParam*1.2; // to steer
		
		double xStartRel = DistrInfoDat.xStart - RadAccessData.xc;
		double xEndRel = DistrInfoDat.xEnd - RadAccessData.xc;
		double dxStart = ::fabs(HalfLambR/xStartRel);
		double dxEnd = ::fabs(HalfLambR/xEndRel);
		dx = ((dxStart < dxEnd)? dxStart : dxEnd)/MultFactor; // req. step
		
		double zStartRel = DistrInfoDat.zStart - RadAccessData.zc;
		double zEndRel = DistrInfoDat.zEnd - RadAccessData.zc;
		double dzStart = ::fabs(HalfLambR/zStartRel);
		double dzEnd = ::fabs(HalfLambR/zEndRel);
		dz = ((dzStart < dzEnd)? dzStart : dzEnd)/MultFactor; // req. step
	}
	else
	{
		dx = (DistrInfoDat.xEnd - DistrInfoDat.xStart)/double(DistrInfoDat.nx);
		dz = (DistrInfoDat.zEnd - DistrInfoDat.zStart)/double(DistrInfoDat.nz);
	}
	xRangeReq = 2.*HalfLambR/dx;
	double xRangeCur = RadAccessData.xStep*RadAccessData.nx;
	pxm = xRangeReq/xRangeCur;
	if(::fabs(pxm - 1.) <= DontResizeThresh) pxm = 1.;
	
	zRangeReq = 2.*HalfLambR/dz;
	double zRangeCur = RadAccessData.zStep*RadAccessData.nz;
	pzm = zRangeReq/zRangeCur;
	if(::fabs(pzm - 1.) <= DontResizeThresh) pzm = 1.;

	srTRadResize RadResBefore_Resol, RadResBefore_Range, RadResAfter;
	//RadResBefore_Resol.DoNotTreatSpherTerm = 1;
	RadResBefore_Resol.doNotTreatSpherTerm(1); //OC090311
	//RadResBefore_Resol.UseOtherSideFFT = 1;
	RadResBefore_Resol.useOtherSideFFT(1); //OC090311
	//RadResBefore_Range.DoNotTreatSpherTerm = 1;
	RadResBefore_Range.doNotTreatSpherTerm(1); //OC090311

	if(pxd < 1.) RadResAfter.pxm = pxd; // no mistake here
	else if(pxd > 1.) RadResBefore_Resol.pxd = pxd;
	if(pzd < 1.) RadResAfter.pzm = pzd; // no mistake here
	else if(pzd > 1.) RadResBefore_Resol.pzd = pzd;

	if(pxm < 1.) RadResAfter.pxd = pxm; // no mistake here
	else if(pxm > 1.) RadResBefore_Range.pxm = pxm;
	if(pzm < 1.) RadResAfter.pzd = pzm; // no mistake here
	else if(pzm > 1.) RadResBefore_Range.pzm = pzm;

	if(result = CheckCorrectParamAndResize(RadAccessData, RadResBefore_Resol, 0)) return result;
	if(result = CheckCorrectParamAndResize(RadAccessData, RadResBefore_Range, 1)) return result;

	srTDriftSpace AuxDrift(DriftLength);
	srTMomentsRatios MomRatios;
	if(result = AuxDrift.PropagateRadMoments(&RadAccessData, &MomRatios)) return result;

	if(result = AuxDrift.PropagateRadiationSimple_PropFromWaist(&RadAccessData)) return result;
	
	if(result = AuxDrift.PropagateWaveFrontRadius(&RadAccessData)) return result;
	RadAccessData.InitialSetupOf4x4PropMatr(DriftLength);

	if(result = CheckCorrectParamAndResize(RadAccessData, RadResAfter, 0)) return result;
	if(result = AuxDrift.RecomputeRadMomentsIfPossible(&RadAccessData)) return result; // Make it at a condition

	return result;
}

//*************************************************************************

int srTSASE::CheckCorrectParamAndResize(srTSRWRadStructAccessData& RadAccessData, srTRadResize& RadResize, char AllowUnderSampling)
{
	int result = 0;
  	const long SmallestN = 8;
	const double DontResizeThresh = 0.15;

	double UnderSamplingX, UnderSamplingZ;
	UnderSamplingX = UnderSamplingZ = 1.;

	double px = RadResize.pxd*RadResize.pxm;
	double pz = RadResize.pzd*RadResize.pzm;
	if(px > 1.0001)
	{
		double ResNp = RadAccessData.nx*px;
		if(ResNp < SmallestN) 
		{
			px = double(SmallestN)/double(RadAccessData.nx);
			if(RadResize.pxd <= RadResize.pxm) RadResize.pxd = px/RadResize.pxm;
			else RadResize.pxm = px/RadResize.pxd;
		}
	}
	if(pz > 1.0001)
	{
		double ResNp = RadAccessData.nz*pz;
		if(ResNp < SmallestN) 
		{
			pz = double(SmallestN)/double(RadAccessData.nz);
			if(RadResize.pzd <= RadResize.pzm) RadResize.pzd = pz/RadResize.pzm;
			else RadResize.pzm = pz/RadResize.pzd;
		}
	}

	srTGenOptElem AuxOptElem;
	if(AllowUnderSampling)
	{
		if(result = AuxOptElem.ReduceBiggerResizeParamUntilFitMemory(RadAccessData, RadResize, UnderSamplingX, UnderSamplingZ)) return result;
	}

	if(result = AuxOptElem.RadResizeGen(RadAccessData, RadResize)) return result;
	if(AllowUnderSampling)
	{
		if(UnderSamplingX > 1. + DontResizeThresh) RadAccessData.UnderSamplingX *= UnderSamplingX;
		if(UnderSamplingZ > 1. + DontResizeThresh) RadAccessData.UnderSamplingZ *= UnderSamplingZ;
	}
	return result;
}

//*************************************************************************

int srTSASE::CreateWavefrontElField(srTSRWRadStructAccessData& RadAccessData)
{//If some "Observation" parameters are defined (in frequency domain !), then
 //- in case of steady-state simulations: straightforward propagation to the observation plane, and then resizing are performed;
 //- in case of TD simulations: transformation to frequency domain is performed, and then, for each frequency, propagation to the observation plane, and resizing (if necessary);
	int result = 0;
	try
	{
		if(result = InitMainGenesisStructs()) throw result; /* Allocate memory; Set default values to the input parameters */
		
		//srTSRWRadStructAccessData arSRWRadAccessData[] = {RadAccessData};
		
		//if(result = GenesisCompDrive(RadAccessData)) throw result;
		if(result = GenesisCompDrive(&RadAccessData, 1)) throw result;

		InitAuxWfrParams(&RadAccessData, 1); //OC041108

		ReleaseGenesisStructs(); /* De-allocate memory; Set default values to the input parameters */
		ControlSASE.DeallocLocalCont();
		SeedRad.Initialize();

		//RadAccessData.SetNonZeroWavefrontLimitsToFullRange();
		//To modify!!!
		//RadAccessData.NormalizeElFieldToArbUnits();

		if(DistrInfoDat.IsDefined())
		{
			//if(result = PropagateWavefrontToObservationPlane(RadAccessData)) throw result;
			if(result = PropagateWavefrontToObservationPlane(&RadAccessData, 1)) throw result;
		}
	}
	//catch(int ErrNo)
	catch(...)
	{
		ReleaseGenesisStructs();
		ControlSASE.DeallocLocalCont();
		SeedRad.Initialize();
		return result;
	}
	return result;
}

//*************************************************************************

int srTSASE::CreateWavefrontElField(srTSRWRadStructAccessData* arSRWRadAccessData, int numHarm)
{//If some "Observation" parameters are defined (in frequency domain !), then
 //- in case of steady-state simulations: straightforward propagation to the observation plane, and then resizing are performed;
 //- in case of TD simulations: transformation to frequency domain is performed, and then, for each frequency, propagation to the observation plane, and resizing (if necessary);
	int result = 0;
	try
	{
		if(result = InitMainGenesisStructs()) throw result; //Allocate memory; Set default values to the input parameters
		//if(result = InitMainGenesisStructs(numHarm)) throw result; //Allocate memory; Set default values to the input parameters

		if(result = GenesisCompDrive(arSRWRadAccessData, numHarm)) throw result;

		InitAuxWfrParams(arSRWRadAccessData, numHarm); //OC041108

		ReleaseGenesisStructs(); //De-allocate memory; Set default values to the input parameters
		ControlSASE.DeallocLocalCont();
		SeedRad.Initialize();

		//RadAccessData.SetNonZeroWavefrontLimitsToFullRange();
		//To modify!!!
		//RadAccessData.NormalizeElFieldToArbUnits();

		if(DistrInfoDat.IsDefined())
		{
			//if(result = PropagateWavefrontToObservationPlane(RadAccessData)) throw result;
			if(result = PropagateWavefrontToObservationPlane(arSRWRadAccessData, numHarm)) throw result;
		}
	}
	//catch(int ErrNo)
	catch(...)
	{
		ReleaseGenesisStructs();
		ControlSASE.DeallocLocalCont();
		SeedRad.Initialize();
		return result;
	}
	return result;
}

//*************************************************************************
