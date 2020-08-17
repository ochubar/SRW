/* mpi.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

Extern struct {
    integer mpi_id__, mpi_err__, mpi_size__, mpi_loop__, nfldmpi, nparmpi, nfldhmpi[6];
} mpicom_;

#define mpicom_1 mpicom_

/*     wrapper for stand-alone version  of genesis */
/*     most of the MPI routines are empty except for the initialization */
/*     where the size=1 and id=0 is assigned */

/*     to compile for mpi, replce this file with mpi.f from the mpi */
/*     subdirectory and remove the file mpif.h */

/*     ---------------------------------- */
/* Subroutine */ int mpi_init__(i1)
integer *i1;
{





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


/*     ------------------------------------------------------------------ */
/*     common block for mpi implementation for genesis */





    mpicom_1.mpi_size__ = 1;
    mpicom_1.mpi_id__ = 0;
    return 0;
} /* mpi_init__ */

/*     --------------------------------- */
/* Subroutine */ int mpi_bcast__(c__, i1, i2, i3, i4, i5, c_len)
char *c__;
integer *i1, *i2, *i3, *i4, *i5;
ftnlen c_len;
{


    return 0;
} /* mpi_bcast__ */

/*     --------------------------------- */
/* Subroutine */ int mpi_send__(c__, i1, i2, i3, i4, i5, i6)
doublecomplex *c__;
integer *i1, *i2, *i3, *i4, *i5, *i6;
{


    return 0;
} /* mpi_send__ */

/*     ---------------------------------- */
/* Subroutine */ int mpi_recv__(c__, i1, i2, i3, i4, i5, i6, i7)
doublecomplex *c__;
integer *i1, *i2, *i3, *i4, *i5, *i6, *i7;
{


    /* Parameter adjustments */
    --i6;

    /* Function Body */
    return 0;
} /* mpi_recv__ */

/*     ---------------------------------- */
/* Subroutine */ int mpi_comm_rank__(i1, i2, i3)
integer *i1, *i2, *i3;
{


    return 0;
} /* mpi_comm_rank__ */

/*     --------------------------------- */
/* Subroutine */ int mpi_comm_size__(i1, i2, i3)
integer *i1, *i2, *i3;
{


    return 0;
} /* mpi_comm_size__ */

/*     ---------------------------------- */
/* Subroutine */ int mpi_finalize__(ierr)
integer *ierr;
{
    return 0;
} /* mpi_finalize__ */

/*     ----------------------------------- */
/* Subroutine */ int mpi_merge__()
{
    return 0;
} /* mpi_merge__ */

