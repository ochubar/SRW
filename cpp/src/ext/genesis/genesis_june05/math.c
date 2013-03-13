/* math.f -- translated by f2c (version 20000118).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b10 = 1.;

doublereal gasham_(j)
integer *j;
{
    /* Initialized data */

    static integer icall = 0;
    static integer iset[26] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	    0,0,0 };
    static doublereal gset[26] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	    0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double log(), sqrt();

    /* Local variables */
    static doublereal r__;
    extern doublereal hammv_();
    static doublereal v1, v2;
    static integer jd;
    static doublereal fac;

/*     ================================================================== */
/*     gaussian hammersley sequence */
/*     reinitializable */
/*     ------------------------------------------------------------------ */


    if (icall == 0 || *j < 0) {
	for (jd = 1; jd <= 26; ++jd) {
	    iset[jd - 1] = 0;
	}
	icall = 1;
    }
    if (iset[*j - 1] == 0) {
L1:
	v1 = hammv_(j) * 2. - 1.;
	i__1 = *j + 1;
	v2 = hammv_(&i__1) * 2. - 1.;
/* Computing 2nd power */
	d__1 = v1;
/* Computing 2nd power */
	d__2 = v2;
	r__ = d__1 * d__1 + d__2 * d__2;
	if (r__ >= (float)1. || r__ == (float)0.) {
	    goto L1;
	}
	fac = sqrt(log(1. - r__) * -2. / r__);
	gset[*j - 1] = v1 * fac;
	ret_val = v2 * fac;
	iset[*j - 1] = 1;
    } else {
	ret_val = gset[*j - 1];
	iset[*j - 1] = 0;
    }

    return ret_val;
} /* gasham_ */




/* gasham */
doublereal hammv_(j)
integer *j;
{
    /* Initialized data */

    static integer icall = 0;
    static integer i__[26] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	    0,0,0 };
    static integer nbase[26] = { 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,
	    59,61,67,71,73,79,83,89,97,101 };

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer i1[26], i2[26], jd;
    static doublereal xs[26], xsi[26];

/*     ================================================================== */
/*     uniform hammersley sequence */
/*     reinitializable */
/*     ------------------------------------------------------------------ */


    if (icall == 0 || *j < 0) {
	for (jd = 1; jd <= 26; ++jd) {
	    i__[jd - 1] = 0;
	}
	icall = 1;
	*j = abs(*j);
    }
    xs[*j - 1] = 0.;
    xsi[*j - 1] = 1.;
    ++i__[*j - 1];
    i2[*j - 1] = i__[*j - 1];
L10:
    xsi[*j - 1] /= (real) nbase[*j - 1];
    i1[*j - 1] = i2[*j - 1] / nbase[*j - 1];
    xs[*j - 1] += (i2[*j - 1] - nbase[*j - 1] * i1[*j - 1]) * xsi[*j - 1];
    i2[*j - 1] = i1[*j - 1];
    if (i2[*j - 1] > 0) {
	goto L10;
    }
    ret_val = xs[*j - 1];

    return ret_val;
} /* hammv_ */




/* hammv */
doublereal ran1_(idum)
integer *idum;
{
    /* Initialized data */

    static integer idum2 = 123456789;
    static integer iv[32] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	    0,0,0,0,0,0,0,0 };
    static integer iy = 0;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer j, k;

/*     ================================================================== */
/*     random number generator from numerical recipes (p. 272f). */
/*     ------------------------------------------------------------------ */


    if (*idum <= 0) {
/* Computing MAX */
	i__1 = -(*idum);
	*idum = max(i__1,1);
	idum2 = *idum;
	for (j = 40; j >= 1; --j) {
	    k = *idum / 53668;
	    *idum = (*idum - k * 53668) * 40014 - k * 12211;
	    if (*idum < 0) {
		*idum += 2147483563;
	    }
	    if (j <= 32) {
		iv[j - 1] = *idum;
	    }
	}
	iy = iv[0];
    }
    k = *idum / 53668;
    *idum = (*idum - k * 53668) * 40014 - k * 12211;
    if (*idum < 0) {
	*idum += 2147483563;
    }
    k = idum2 / 52774;
    idum2 = (idum2 - k * 52774) * 40692 - k * 3791;
    if (idum2 < 0) {
	idum2 += 2147483399;
    }
    j = iy / 67108862 + 1;
    iy = iv[j - 1] - idum2;
    iv[j - 1] = *idum;
    if (iy < 1) {
	iy += 2147483562;
    }
/* Computing MIN */
    d__1 = iy * 4.6566130573917691e-10;
    ret_val = min(d__1,1.);
    return ret_val;
} /* ran1_ */



doublereal gasran_(idum)
integer *idum;
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double log(), sqrt();

    /* Local variables */
    static doublereal gset;
    static integer iset;
    static doublereal r__, v1, v2, fac;
    extern doublereal ran1_();

/*     ================================================================== */
/*     random number generator from numerical recipes (p. 272f). */
/*     ------------------------------------------------------------------ */


    if (iset == 0) {
L1:
	v1 = ran1_(idum) * 2. - 1.;
	v2 = ran1_(idum) * 2. - 1.;
/* Computing 2nd power */
	d__1 = v1;
/* Computing 2nd power */
	d__2 = v2;
	r__ = d__1 * d__1 + d__2 * d__2;
	if (r__ >= (float)1. || r__ == (float)0.) {
	    goto L1;
	}
	fac = sqrt(log(1. - r__) * -2. / r__);
	gset = v1 * fac;
	ret_val = v2 * fac;
	iset = 1;
    } else {
	ret_val = gset;
	iset = 0;
    }

    return ret_val;
} /* gasran_ */




doublereal bessj0_(x)
doublereal *x;
{
    /* Initialized data */

    static doublereal p1 = 1.;
    static doublereal p2 = -.001098628627;
    static doublereal p3 = 2.734510407e-5;
    static doublereal p4 = -2.073370639e-6;
    static doublereal p5 = 2.093887211e-7;
    static doublereal q1 = -.01562499995;
    static doublereal q2 = 1.430488765e-4;
    static doublereal q3 = -6.911147651e-6;
    static doublereal q4 = 7.621095161e-7;
    static doublereal q5 = -9.34945152e-8;
    static doublereal r1 = 57568490574.;
    static doublereal r2 = -13362590354.;
    static doublereal r3 = 651619640.7;
    static doublereal r4 = -11214424.18;
    static doublereal r5 = 77392.33017;
    static doublereal r6 = -184.9052456;
    static doublereal s1 = 57568490411.;
    static doublereal s2 = 1029532985.;
    static doublereal s3 = 9494680.718;
    static doublereal s4 = 59272.64853;
    static doublereal s5 = 267.8532712;
    static doublereal s6 = 1.;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double cos(), sin(), sqrt();

    /* Local variables */
    static doublereal y, z__, ax, xx;

/*     ================================================================== */
/*     bessel function j0 - numerical rec. */
/*     ------------------------------------------------------------------ */

    if (abs(*x) < (float)8.) {
/* Computing 2nd power */
	d__1 = *x;
	y = d__1 * d__1;
	ret_val = (r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * r6))))) / 
		(s1 + y * (s2 + y * (s3 + y * (s4 + y * (s5 + y * s6)))));
    } else {
	ax = abs(*x);
	z__ = (float)8. / ax;
/* Computing 2nd power */
	d__1 = z__;
	y = d__1 * d__1;
	xx = ax - (float).785398164;
	ret_val = sqrt((float).636619772 / ax) * (cos(xx) * (p1 + y * (p2 + y 
		* (p3 + y * (p4 + y * p5)))) - z__ * sin(xx) * (q1 + y * (q2 
		+ y * (q3 + y * (q4 + y * q5)))));
    }

    return ret_val;
} /* bessj0_ */




doublereal bessj1_(x)
doublereal *x;
{
    /* Initialized data */

    static doublereal r1 = 72362614232.;
    static doublereal r2 = -7895059235.;
    static doublereal r3 = 242396853.1;
    static doublereal r4 = -2972611.439;
    static doublereal r5 = 15704.4826;
    static doublereal r6 = -30.16036606;
    static doublereal s1 = 144725228442.;
    static doublereal s2 = 2300535178.;
    static doublereal s3 = 18583304.74;
    static doublereal s4 = 99447.43394;
    static doublereal s5 = 376.9991397;
    static doublereal s6 = 1.;
    static doublereal p1 = 1.;
    static doublereal p2 = .00183105;
    static doublereal p3 = -3.516396496e-5;
    static doublereal p4 = 2.457520174e-6;
    static doublereal p5 = -2.40337019e-7;
    static doublereal q1 = .04687499995;
    static doublereal q2 = -2.002690873e-4;
    static doublereal q3 = 8.449199096e-6;
    static doublereal q4 = -8.8228987e-7;
    static doublereal q5 = 1.05787412e-7;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double cos(), sin(), sqrt(), d_sign();

    /* Local variables */
    static doublereal y, z__, ax, xx;

/*     ================================================================== */
/*     bessel function j1 - numerical rec. */
/*     ------------------------------------------------------------------ */

    if (abs(*x) < (float)8.) {
/* Computing 2nd power */
	d__1 = *x;
	y = d__1 * d__1;
	ret_val = *x * (r1 + y * (r2 + y * (r3 + y * (r4 + y * (r5 + y * r6)))
		)) / (s1 + y * (s2 + y * (s3 + y * (s4 + y * (s5 + y * s6)))))
		;
    } else {
	ax = abs(*x);
	z__ = (float)8. / ax;
/* Computing 2nd power */
	d__1 = z__;
	y = d__1 * d__1;
	xx = ax - (float)2.356194491;
	ret_val = sqrt((float).636619772 / ax) * (cos(xx) * (p1 + y * (p2 + y 
		* (p3 + y * (p4 + y * p5)))) - z__ * sin(xx) * (q1 + y * (q2 
		+ y * (q3 + y * (q4 + y * q5))))) * d_sign(&c_b10, x);
    }

    return ret_val;
} /* bessj1_ */




integer luf_(x, table, n)
doublereal *x, *table;
integer *n;
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer jl, jm, ju;

/*     ================================================================== */
/*     luf is a table lookup function that locates a value x between */
/*     elements of an increasing table of size n. */
/*     luf is the value of the index after the table location which */
/*     x corresponds to. */
/*     the routine uses a bisection methode (numerical rec.) */
/*     the array table must be monotonic */
/*     luf=1 or luf=n+2 is returned to indicate out of range */
/*     ------------------------------------------------------------------ */


    /* Parameter adjustments */
    --table;

    /* Function Body */
    jl = 0;
/* lower limit */
    ju = *n + 1;
/* upper limit */
L10:
    if (ju - jl > 1) {
	jm = (ju + jl) / 2;
/* midpoint */
	if (table[*n] > table[1] == *x > table[jm]) {
	    jl = jm;
	} else {
	    ju = jm;
	}
	goto L10;
    }
    ret_val = jl + 1;

    return ret_val;
} /* luf_ */




/* luf */
/* Subroutine */ int fourn_(data, nn, ndim, isign)
doublereal *data;
integer *nn, *ndim, *isign;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double sin();

    /* Local variables */
    static integer idim, ibit, nrem, ntot, i2rev, i3rev, n;
    static doublereal theta, tempi, tempr;
    static integer i1, i2, i3, k1, k2, nprev;
    static doublereal wtemp, wi, wr;
    static integer ip1, ip2, ip3;
    static doublereal wpi, wpr;
    static integer ifp1, ifp2;

/*     ================================================================= */
/*     multidimensional fft of complex values (num. rec.) */
/*     number of elements of data must be a power of 2! */
/*     nn = number of dimension */
/*     ndim(nn) = elements per dimension */
/*     isign = 1 fft, = -1 inverse fft */
/*     ----------------------------------------------------------------- */

    /* Parameter adjustments */
    --data;
    --nn;

    /* Function Body */
    ntot = 1;
    i__1 = *ndim;
    for (idim = 1; idim <= i__1; ++idim) {
	ntot *= nn[idim];
/* L11: */
    }
    nprev = 1;
    i__1 = *ndim;
    for (idim = 1; idim <= i__1; ++idim) {
	n = nn[idim];
	nrem = ntot / (n * nprev);
	ip1 = nprev << 1;
	ip2 = ip1 * n;
	ip3 = ip2 * nrem;
	i2rev = 1;
	i__2 = ip2;
	i__3 = ip1;
	for (i2 = 1; i__3 < 0 ? i2 >= i__2 : i2 <= i__2; i2 += i__3) {
	    if (i2 < i2rev) {
		i__4 = i2 + ip1 - 2;
		for (i1 = i2; i1 <= i__4; i1 += 2) {
		    i__5 = ip3;
		    i__6 = ip2;
		    for (i3 = i1; i__6 < 0 ? i3 >= i__5 : i3 <= i__5; i3 += 
			    i__6) {
			i3rev = i2rev + i3 - i2;
			tempr = data[i3];
			tempi = data[i3 + 1];
			data[i3] = data[i3rev];
			data[i3 + 1] = data[i3rev + 1];
			data[i3rev] = tempr;
			data[i3rev + 1] = tempi;
/* L12: */
		    }
/* L13: */
		}
	    }
	    ibit = ip2 / 2;
L1:
	    if (ibit >= ip1 && i2rev > ibit) {
		i2rev -= ibit;
		ibit /= 2;
		goto L1;
	    }
	    i2rev += ibit;
/* L14: */
	}
	ifp1 = ip1;
L2:
	if (ifp1 < ip2) {
	    ifp2 = ifp1 << 1;
	    theta = *isign * 6.28318530717959 / (ifp2 / ip1);
/* Computing 2nd power */
	    d__1 = sin(theta * .5);
	    wpr = d__1 * d__1 * -2.;
	    wpi = sin(theta);
	    wr = 1.;
	    wi = 0.;
	    i__3 = ifp1;
	    i__2 = ip1;
	    for (i3 = 1; i__2 < 0 ? i3 >= i__3 : i3 <= i__3; i3 += i__2) {
		i__4 = i3 + ip1 - 2;
		for (i1 = i3; i1 <= i__4; i1 += 2) {
		    i__6 = ip3;
		    i__5 = ifp2;
		    for (i2 = i1; i__5 < 0 ? i2 >= i__6 : i2 <= i__6; i2 += 
			    i__5) {
			k1 = i2;
			k2 = k1 + ifp1;
			tempr = (real) wr * data[k2] - (real) wi * data[k2 + 
				1];
			tempi = (real) wr * data[k2 + 1] + (real) wi * data[
				k2];
			data[k2] = data[k1] - tempr;
			data[k2 + 1] = data[k1 + 1] - tempi;
			data[k1] += tempr;
			data[k1 + 1] += tempi;
/* L15: */
		    }
/* L16: */
		}
		wtemp = wr;
		wr = wr * wpr - wi * wpi + wr;
		wi = wi * wpr + wtemp * wpi + wi;
/* L17: */
	    }
	    ifp1 = ifp2;
	    goto L2;
	}
	nprev = n * nprev;
/* L18: */
    }
    return 0;
} /* fourn_ */

