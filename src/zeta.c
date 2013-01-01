/* zeta.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Subroutine */ int zetawr_(doublereal *x, doublereal *ans, integer *deriv, 
	integer *nn)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal b2[12];
    extern doublereal zeta_(doublereal *, doublereal *), dzeta_(doublereal *, 
	    doublereal *);
    extern /* Subroutine */ int becoef_(doublereal *);
    extern doublereal ddzeta_(doublereal *, doublereal *);

    /* Parameter adjustments */
    --ans;
    --x;

    /* Function Body */
    becoef_(b2);
    if (! (*deriv == 0)) {
	goto L23000;
    }
    i__1 = *nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ans[i__] = zeta_(&x[i__], b2);
/* L23002: */
    }
L23000:
    if (! (*deriv == 1)) {
	goto L23004;
    }
    i__1 = *nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ans[i__] = dzeta_(&x[i__], b2);
/* L23006: */
    }
L23004:
    if (! (*deriv == 2)) {
	goto L23008;
    }
    i__1 = *nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ans[i__] = ddzeta_(&x[i__], b2);
/* L23010: */
    }
L23008:
    return 0;
} /* zetawr_ */

doublereal zeta_(doublereal *s, doublereal *b2)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer a, k, m, n;
    static doublereal p, a2;
    static integer m2;
    static doublereal sum, fred;

    /* Parameter adjustments */
    --b2;

    /* Function Body */
    a = 12;
    k = 8;
    a2 = (doublereal) (a * a);
    p = *s / 2. / a2;
    sum = 1. / (*s - 1.) + .5 / a + b2[1] * p;
    i__1 = k;
    for (m = 2; m <= i__1; ++m) {
	m2 = m + m;
	p = p * (*s + m2 - 3.) * (*s + m2 - 2.) / (m2 - 1.) / m2 / a2;
	sum += p * b2[m];
/* L23012: */
    }
    fred = exp((*s - 1.) * log(a * 1.));
    sum = sum / fred + 1.;
    i__1 = a - 1;
    for (n = 2; n <= i__1; ++n) {
	sum += exp(-(*s) * log(n * 1.)) * 1.;
/* L23014: */
    }
    ret_val = sum;
    return ret_val;
} /* zeta_ */

doublereal dzeta_(doublereal *s, doublereal *b2)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer a, k, m, n;
    static doublereal p, q, a2;
    static integer m2;
    static doublereal sum, fred, loga, logn;

    /* Parameter adjustments */
    --b2;

    /* Function Body */
    a = 12;
    k = 8;
    loga = log(a * 1.);
    a2 = (doublereal) (a * a);
    p = *s / 2. / a2;
    q = 1. / *s - loga;
    sum = b2[1] * p * q;
    i__1 = k;
    for (m = 2; m <= i__1; ++m) {
	m2 = m + m;
	p = p * (*s + m2 - 3.) * (*s + m2 - 2.) / (m2 - 1.) / m2 / a2;
	q = q + 1. / (*s + m2 - 3.) + 1. / (*s + m2 - 2.);
	sum += b2[m] * p * q;
/* L23016: */
    }
    fred = exp((1. - *s) * loga);
/* Computing 2nd power */
    d__1 = *s - 1.;
    sum = sum - 1. / (d__1 * d__1) - loga * (1. / (*s - 1.) + .5 / a);
    sum *= fred;
    i__1 = a - 1;
    for (n = 2; n <= i__1; ++n) {
	logn = log(n * 1.);
	sum -= logn / exp(logn * *s);
/* L23018: */
    }
    ret_val = sum;
    return ret_val;
} /* dzeta_ */

doublereal upsilon_(doublereal *s, doublereal *b2)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    extern doublereal zeta_(doublereal *, doublereal *), dzeta_(doublereal *, 
	    doublereal *);

    /* Parameter adjustments */
    --b2;

    /* Function Body */
    ret_val = -dzeta_(s, &b2[1]) / zeta_(s, &b2[1]);
    return ret_val;
} /* upsilon_ */

doublereal ddzeta_(doublereal *s, doublereal *b2)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer a, k, m, n;
    static doublereal p, q, r__, a2;
    static integer m2;
    static doublereal sum, fred, loga, logn, fred2;

    /* Parameter adjustments */
    --b2;

    /* Function Body */
    a = 12;
    k = 8;
    loga = log(a * 1.);
    a2 = (doublereal) (a * a);
    p = *s / 2. / a2;
    q = 1. / *s - loga;
    r__ = 1. / *s / *s;
    sum = b2[1] * p * (q * q - r__);
    i__1 = k;
    for (m = 2; m <= i__1; ++m) {
	m2 = m + m;
	p = p * (*s + m2 - 3.) * (*s + m2 - 2.) / (m2 - 1.) / m2 / a2;
	q = q + 1. / (*s + m2 - 3.) + 1. / (*s + m2 - 2.);
/* Computing 2nd power */
	d__1 = *s + m2 - 3.;
/* Computing 2nd power */
	d__2 = *s + m2 - 2.;
	r__ = r__ + 1. / (d__1 * d__1) + 1. / (d__2 * d__2);
	sum += b2[m] * p * (q * q - r__);
/* L23020: */
    }
    fred = exp((1. - *s) * loga);
/* Computing 2nd power */
    d__1 = loga;
    fred2 = d__1 * d__1 * (1. / (*s - 1.) + .5 / a);
/* Computing 3rd power */
    d__1 = *s - 1.;
/* Computing 2nd power */
    d__2 = *s - 1.;
    sum = sum + 2. / (d__1 * (d__1 * d__1)) + loga * 2. / (d__2 * d__2) + 
	    fred2;
    sum *= fred;
    i__1 = a - 1;
    for (n = 2; n <= i__1; ++n) {
	logn = log(n * 1.);
/* Computing 2nd power */
	d__1 = logn;
	sum += d__1 * d__1 / exp(logn * *s);
/* L23022: */
    }
    ret_val = sum;
    return ret_val;
} /* ddzeta_ */

doublereal duds_(doublereal *s, doublereal *b2)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal zs;
    extern doublereal zeta_(doublereal *, doublereal *), dzeta_(doublereal *, 
	    doublereal *), ddzeta_(doublereal *, doublereal *);

    /* Parameter adjustments */
    --b2;

    /* Function Body */
    zs = zeta_(s, &b2[1]);
/* Computing 2nd power */
    d__1 = dzeta_(s, &b2[1]) / zs;
    ret_val = d__1 * d__1 - ddzeta_(s, &b2[1]) / zs;
    return ret_val;
} /* duds_ */

/* Subroutine */ int becoef_(doublereal *b2)
{
    /* Parameter adjustments */
    --b2;

    /* Function Body */
    b2[1] = .16666666666666666;
    b2[2] = -.033333333333333333;
    b2[3] = .023809523809523808;
    b2[4] = -.033333333333333333;
    b2[5] = .07575757575757576;
    b2[6] = -.2531135531135531;
    b2[7] = 1.1666666666666667;
    b2[8] = -7.0921568627450977;
    b2[9] = 54.971177944862156;
    b2[10] = -529.12424242424242;
    b2[11] = 6192.123188405797;
    b2[12] = -86580.253113553117;
    return 0;
} /* becoef_ */

