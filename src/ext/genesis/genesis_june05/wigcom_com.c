
#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

struct {
	doublereal *awz, *awdz, *solz, *awerx, *awery, *qfld, 
		fbess, magversion, unitlength, *dqfx, *dqfy, *awdx, *awdy;
	integer iran, nstepz, itap;
} wigcom_;

#ifdef __cplusplus
}
#endif
