
#include "f2c.h"

#ifdef __cplusplus
extern "C" {
#endif

struct {
	doublereal *tgam0, *tdgam, *temitx, *temity, *txrms, *tyrms, *txpos, 
		*typos, *tpxpos, *tpypos, *talphx, *talphy, *tcurrent, *tpos, *tloss, 
		*distgam, *distx, *disty, *distpx, *distpy, *distt,
		*tradpos, *tzrayl, *tzwaist, *tprad0, *tradphase;
	integer ndata, nsep, nslp, ndist, nraddata;
} tbunchcom_;

#ifdef __cplusplus
}
#endif
