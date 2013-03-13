
#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

struct {
	doublereal *error, *gain, *phimid, *whalf, *logp, *powmid, *xrms, *yrms, *xpos, 
		*ypos, *pmodhist, *gamhist, *diver, pradol, pinit, 
		*bunphase, *dgamhist, *ffield;
	integer ihist;
} diagcom_;

#ifdef __cplusplus
}
#endif
