
#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

struct {
	doublecomplex *crfield, *crsource, crmatc[513], cstep, *crhm, cwet[513], cbet[513];
	doublereal dxy, xks, radphase;
} cartcom_;

#ifdef __cplusplus
}
#endif
