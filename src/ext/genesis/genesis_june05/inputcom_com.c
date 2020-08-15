
#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

struct {
	doublereal aw0,	xkx, xky, wcoefz[3], xlamd,	fbess0,	delaw, awd,	awx, awy, 
		gamma0,	delgam,	rxbeam,	rybeam,	alphax,	alphay,	emitx, emity, 
		xbeam, ybeam, pxbeam, pybeam, cuttail, curpeak,	conditx, condity, 
		bunch, bunchphase, emod, emodphase,	xlamds,	prad0, zrayl, zwaist, 
		rmax0, zsep, delz, zstop, quadf, quadd,	fl,	dl,	drl, f1st, qfdx, 
		qfdy, sl, solen, curlen, shotnoise,	svar, dgrid, eloss,	version, 
		ibfield, imagl,	idril;
	integer	iseed, nwig, nsec, npart, ncar,	lbc, nscr, nscz, nptr, ildgam, 
		ildpsi,	ildx, ildy,	ildpx, ildpy, itgaus, nbins, iphsty, ishsty, 
		ippart,	ispart,	ipradi,	isradi,	iertyp,	iwityp,	idump, iotail, 
		nharm, magin, magout, lout[35],	ffspec,	ntail, nslice, iall, itdp,
		ipseed,	iscan,	nscan, isntyp, isravg, isrsig, iorb, ndcut,	
		idmppar, idmpfld, ilog,	igamgaus, convharm,	alignradf, offsetradf,
		multconv;
	char beamfile[30], fieldfile[30], maginfile[30], magoutfile[30], 
		outputfile[30], inputfile[30], scan[30], distfile[30], partfile[30], filetype[30], radfile[30];
} inputcom_;

#ifdef __cplusplus
}
#endif
