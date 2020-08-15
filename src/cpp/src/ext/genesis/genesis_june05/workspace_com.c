
#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

struct {
	doublecomplex *crwork3, *cpart1, *cpart2;
	doublereal *k2gg, *k2pp, *k3gg, *k3pp, *p1, *p2;
} workspace_;

#ifdef __cplusplus
}
#endif
