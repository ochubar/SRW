/* string.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;

integer extractnumber_(line, val, line_len)
char *line;
doublereal *val;
ftnlen line_len;
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer i_len(), i_indx();
    /* Subroutine */ int s_copy();
    integer s_rsli(), do_lio(), e_rsli();

    /* Local variables */
    extern /* Subroutine */ int getfirstchar_();
    static integer i__, j, n;
    static char cline[255];

    /* Fortran I/O blocks */
    static icilist io___5 = { 0, cline, 0, 0, 255, 1 };


/*     ====================================================================== */
/*     extract float from line, ignoring preceeding '='-signs */
/*     ====================================================================== */


    ret_val = 0;
    n = i_len(line, line_len);
    i__ = i_indx(line, "=", line_len, (ftnlen)1) + 1;
/* find equal sign */
    s_copy(cline, line + (i__ - 1), (ftnlen)255, n - (i__ - 1));
/* cut string */
    getfirstchar_(cline, &j, (ftnlen)255);
/* check if string is empty */
    if (j == 0) {
	ret_val = -1;
	return ret_val;
    }
    s_rsli(&io___5);
    do_lio(&c__5, &c__1, (char *)&(*val), (ftnlen)sizeof(doublereal));
    e_rsli();
    return ret_val;
} /* extractnumber_ */


integer extractval_(line, values, nval, line_len)
char *line;
doublereal *values;
integer *nval;
ftnlen line_len;
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    icilist ici__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy();
    integer s_rsli(), do_lio(), e_rsli();

    /* Local variables */
    extern /* Subroutine */ int getfirstchar_();
    static integer i__, j;
    static char cline[255];
    static integer ix1, ix2;

/*     ====================================================================== */
/*     extract nval data out of line */
/*     ====================================================================== */


    /* Parameter adjustments */
    --values;

    /* Function Body */
    ret_val = 0;
    s_copy(cline, line, (ftnlen)255, line_len);
    i__1 = *nval;
    for (i__ = 1; i__ <= i__1; ++i__) {
	getfirstchar_(cline, &ix1, (ftnlen)255);
/* check for characters */
	if (ix1 == 0) {
/* empty string */
	    ret_val = -1;
	    return ret_val;
	}
	ix2 = 255;
/* search backwards */
	i__2 = ix1 + 1;
	for (j = 255; j >= i__2; --j) {
/* for first space after ix1 */
	    if (*(unsigned char *)&cline[j - 1] == ' ') {
		ix2 = j;
	    }
	}
	s_copy(line, cline + (ix1 - 1), line_len, ix2 - (ix1 - 1));
/* copy word */
	ici__1.icierr = 0;
	ici__1.iciend = 0;
	ici__1.icirnum = 1;
	ici__1.icirlen = line_len;
	ici__1.iciunit = line;
	ici__1.icifmt = 0;
	s_rsli(&ici__1);
	do_lio(&c__5, &c__1, (char *)&values[i__], (ftnlen)sizeof(doublereal))
		;
	e_rsli();
/* get value */
	s_copy(line, cline + (ix2 - 1), line_len, 255 - (ix2 - 1));
/* copy remaining part of the line */
	s_copy(cline, line, (ftnlen)255, line_len);
/* to cline */
    }
    return ret_val;
} /* extractval_ */


/* Subroutine */ int getfirstchar_(line, idx, line_len)
char *line;
integer *idx;
ftnlen line_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len();

    /* Local variables */
    static integer i__;

/*     ====================================================================== */
/*     get the index of the first non space character */
/*     ====================================================================== */


    *idx = 0;
    i__1 = i_len(line, line_len);
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*(unsigned char *)&line[i__ - 1] > ' ' && *idx < 1) {
	    *idx = i__;
	}
    }
    return 0;
} /* getfirstchar_ */


/* Subroutine */ int touppercase_(c__, c_len)
char *c__;
ftnlen c_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer i_len();

    /* Local variables */
    static integer i__, ic;

/*     ====================================================================== */
/*     convert string to upper case letters */
/*     ====================================================================== */

    i__1 = i_len(c__, c_len);
    for (i__ = 1; i__ <= i__1; ++i__) {
	ic = *(unsigned char *)&c__[i__ - 1];
	if (ic > 96 && ic < 123) {
	    ic += -32;
	    *(unsigned char *)&c__[i__ - 1] = (char) ic;
/* replace with uppercase character */
	}
    }
    return 0;
} /* touppercase_ */


integer strlen_(line, line_len)
char *line;
ftnlen line_len;
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer i_len(), i_indx();

    /* Local variables */
    static integer nchar, nchar1;

/*     =================================================================== */
/*     check length of given string */
/*     =================================================================== */

    ret_val = 1;
    nchar1 = i_len(line, line_len);
L999:
    nchar = i_indx(line + (nchar1 - 1), " ", (ftnlen)1, (ftnlen)1);
    if (nchar != 0) {
	if (nchar1 > 1) {
	    --nchar1;
	    goto L999;
	}
    }
    ret_val = nchar1;

    return ret_val;
} /* strlen_ */

