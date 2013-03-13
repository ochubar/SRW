
#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

struct {
	doublereal *xpart, *ypart, *px, *py, *gamma,
		*theta, *xporb, *yporb, *btpar, *btper, *ez, *wx, *wy, 
		xcuren, dedz, tdmin, tdmax, delcharge, dtd, charge;
	integer *lostid, lost, losttot, *ipos;
} beamcom_;

#ifdef __cplusplus
}
#endif
