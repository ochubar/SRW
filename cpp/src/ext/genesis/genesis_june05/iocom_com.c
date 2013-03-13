
#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

struct {
	doublereal distversion, distrev;
	integer iout[24], nout, nfld, npar, ndump, firstout, nfin, irecpar, 
		irecfld, kout, ndmp2, npin, nlog, ndis, ncoldis, iconv2t, iconv2g,
		icolpar[10], ndistsize, iconv2px, iconv2py, nprobe, ftype, 
		ftdist, ftpart, ftfield;
} iocom_;

#ifdef __cplusplus
}
#endif
