/* timerec.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

Extern struct {
    //doublecomplex crtime[15960700];
	doublecomplex *crtime;
    integer ntmp;
} tslipcom_;

#define tslipcom_1 tslipcom_

/* Table of constant values */

static integer c_n1 = -1;
static integer c__2 = 2;
static integer c_n2 = -2;

/* Subroutine */ int opentimerec_(n)
integer *n;
{
    /* System generated locals */
    integer i__1;
    olist o__1;

    /* Builtin functions */
    integer f_open();

    /* Local variables */
    extern /* Subroutine */ int last_();
    extern integer printerr_();

/*     ================================================================== */
/*     opens a scratch-file for the time record if necessary */
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
/* maximum of particle i */
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


/*     file:   timerec.com */

/*     time dependent common block - huge!!!! */
/*     ------------------------------------------------------------------ */




/*     -------------------------------------------------------------------- */
/*     slippage field */

/* must be > than ncar*ncar*nslp */
/* <>0 -> slippage is stored on disk */




    if (FALSE_) {
/* write crtime to disk? */
	tslipcom_1.ntmp = 13;
	o__1.oerr = 1;
	o__1.ounit = tslipcom_1.ntmp;
	o__1.ofnm = 0;
	o__1.orl = *n * *n << 4;
	o__1.osta = "scratch";
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
    tslipcom_1.ntmp = printerr_(&c_n1, "CRTIME-scratch", (ftnlen)14);
    last_();
    return 0;
} /* opentimerec_ */

/* Subroutine */ int pushtimerec_(cpush, n, irec)
doublecomplex *cpush;
integer *n, *irec;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_wdue(), do_uio(), e_wdue();

    /* Local variables */
    extern /* Subroutine */ int last_();
    extern integer printerr_();
    static integer it;

    /* Fortran I/O blocks */
    static cilist io___2 = { 1, 0, 0, 0, 0 };


/*     ================================================================== */
/*     copies cpush into crtime array/file */
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
/* maximum of particle i */
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

/*     file:   timerec.com */

/*     time dependent common block - huge!!!! */
/*     ------------------------------------------------------------------ */




/*     -------------------------------------------------------------------- */
/*     slippage field */

/* must be > than ncar*ncar*nslp */
/* <>0 -> slippage is stored on disk */



    /* Parameter adjustments */
    --cpush;

    /* Function Body */
    if (TRUE_) {
	i__1 = *n * *n;
	for (it = 1; it <= i__1; ++it) {
	    i__2 = (*irec - 1) * *n * *n + it - 1;
	    i__3 = it;
	    tslipcom_1.crtime[i__2].r = cpush[i__3].r, tslipcom_1.crtime[i__2]
		    .i = cpush[i__3].i;
	}
    } else {
	io___2.ciunit = tslipcom_1.ntmp;
	io___2.cirec = *irec;
	i__1 = s_wdue(&io___2);
	if (i__1 != 0) {
	    goto L10;
	}
	i__2 = *n * *n;
	for (it = 1; it <= i__2; ++it) {
	    i__1 = do_uio(&c__2, (char *)&cpush[it], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L10;
	    }
	}
	i__1 = e_wdue();
	if (i__1 != 0) {
	    goto L10;
	}
    }
    return 0;
L10:
    it = printerr_(&c_n2, "CRTIME-scratch", (ftnlen)14);
    last_();
    return 0;
} /* pushtimerec_ */



/* of pushtimerec */
/* Subroutine */ int pulltimerec_(cpull, n, irec)
doublecomplex *cpull;
integer *n, *irec;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_rdue(), do_uio(), e_rdue();

    /* Local variables */
    extern /* Subroutine */ int last_();
    extern integer printerr_();
    static integer it;

    /* Fortran I/O blocks */
    static cilist io___4 = { 1, 0, 0, 0, 0 };


/*     ================================================================== */
/*     copies crtime array/file into cpush */
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
/* maximum of particle i */
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

/*     file:   timerec.com */

/*     time dependent common block - huge!!!! */
/*     ------------------------------------------------------------------ */




/*     -------------------------------------------------------------------- */
/*     slippage field */

/* must be > than ncar*ncar*nslp */
/* <>0 -> slippage is stored on disk */



    /* Parameter adjustments */
    --cpull;

    /* Function Body */
    if (TRUE_) {
	i__1 = *n * *n;
	for (it = 1; it <= i__1; ++it) {
	    i__2 = it;
	    i__3 = (*irec - 1) * *n * *n + it - 1;
	    cpull[i__2].r = tslipcom_1.crtime[i__3].r, cpull[i__2].i = 
		    tslipcom_1.crtime[i__3].i;
	}
    } else {
	io___4.ciunit = tslipcom_1.ntmp;
	io___4.cirec = *irec;
	i__1 = s_rdue(&io___4);
	if (i__1 != 0) {
	    goto L10;
	}
	i__2 = *n * *n;
	for (it = 1; it <= i__2; ++it) {
	    i__1 = do_uio(&c__2, (char *)&cpull[it], (ftnlen)sizeof(
		    doublereal));
	    if (i__1 != 0) {
		goto L10;
	    }
	}
	i__1 = e_rdue();
	if (i__1 != 0) {
	    goto L10;
	}
    }
    return 0;
L10:
    it = printerr_(&c_n2, "CRTIME-scratch", (ftnlen)14);
    last_();
    return 0;
} /* pulltimerec_ */

/* of pulltimerec */
/* Subroutine */ int closetimerec_()
{
    /* System generated locals */
    cllist cl__1;

    /* Builtin functions */
    integer f_clos();

/*     ==================================================================== */
/*     close the scratch file if necessary */
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
/* maximum of particle i */
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


/*     file:   timerec.com */

/*     time dependent common block - huge!!!! */
/*     ------------------------------------------------------------------ */




/*     -------------------------------------------------------------------- */
/*     slippage field */

/* must be > than ncar*ncar*nslp */
/* <>0 -> slippage is stored on disk */



    if (FALSE_) {
	cl__1.cerr = 0;
	cl__1.cunit = tslipcom_1.ntmp;
	cl__1.csta = 0;
	f_clos(&cl__1);
    }
    return 0;
} /* closetimerec_ */

