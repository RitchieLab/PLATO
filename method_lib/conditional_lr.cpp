//conditional_lr.cpp

/// Functions converted from fortran
#include <iostream>
using namespace std;
/* likelihood.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/* ********************************************************************* */
/* This file was edited by bill bardsley on 18/08/2002 to correct four */
/* serious errors. To locate these, search for bill.bardsley@man.ac.uk. */
/* Unless the errors are corrected this code returns incorrect results. */
/* ********************************************************************* */

#include "conditional_lr.h"

namespace Methods{
/* Subroutine */ int logcch_(integer *ns, integer *nca, integer *nct, integer 
	*nimax, integer *nmax, integer *nmax1, integer *nvmax, integer *
	nvmax1, integer *nv, doublereal *z__, integer *ivar, doublereal *covi, doublereal *cntr,
	 doublereal *w, doublereal *wb, doublereal *wdb, doublereal *wd2b, doublereal *u, integer *ins, 
	doublereal *db, doublereal *d2b, doublereal *dl, doublereal *b, doublereal *cov, doublereal *chi2, doublereal *
	st, integer *ifault)
{
    /* Initialized data */

    static integer maxit = 50;
    static doublereal eps = (doublereal)1e-5;
    static doublereal zero = (doublereal)0.;
    static doublereal one = (doublereal)1.;
    static doublereal two = (doublereal)2.;

    /* System generated locals */
    integer z_dim1, z_offset, cntr_dim1, cntr_offset, wdb_dim1, wdb_offset, 
	    wd2b_dim1, wd2b_offset, i__1, i__2, i__3;
    doublereal r__1;

    /* Local variables */
    static integer i__, j, k, l, m, n;
    static doublereal t, c1;
    static integer i1, j1, k1, i2;
    static doublereal fm;
    static integer jj, kk, ir, is, nid;
    static doublereal bmn;
    static integer its;
    static doublereal ttt, rlik, rlikp, rliks, const__;

    static integer nullty;


/*        ALGORITHM AS 196 APPL. STATIST. (1984) VOL.33, NO.1 */

/*        LOGISTIC ANALYSIS OF CASE-CONTROL STUDIES */

/*     *** WARNING  This file has been input using a scanner and may */
/*                  contain errors. */

// output the parameters
// cout << "PARAMETERS IN logcch_" << endl;
// cout << "logcch_ ns=" << *ns << endl;
// cout << "logcch_ nimax=" << *nimax << endl;
// cout << "logcch_ nmax=" << *nmax << endl;
// cout << "logcch_ nmax1=" << *nmax1 << endl;
// cout << "logcch_ nvmax=" << *nvmax << endl;
// cout << "logcch_ nvmax1=" << *nvmax1 << endl;
// cout << "logcch_ nv=" << *nv << endl;
// // for(int i=0; i<*ns; i++)
// //   cout << i << " nca[i] " << nca[i] << " nct[i] " << nct[i] << endl;
// for(int i=0; i<*nvmax; i++)
//   cout << i << " ivar[i] " << ivar[i] << endl;
// for(int i=0; i<*nvmax; i++)
//   cout << i << " b[i] " << b[i] << endl;
// cout << endl;
*chi2=0.0;


    /* Parameter adjustments */
    --ins;
    --nct;
    --nca;
    --u;
    --wb;
    --b;
    --dl;
    --db;
    wdb_dim1 = *nvmax;
    wdb_offset = 1 + wdb_dim1;
    wdb -= wdb_offset;
    --w;
    cntr_dim1 = *nvmax;
    cntr_offset = 1 + cntr_dim1;
    cntr -= cntr_offset;
    --ivar;
    z_dim1 = *nvmax;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --cov;
    --d2b;
    wd2b_dim1 = *nvmax1;
    wd2b_offset = 1 + wd2b_dim1;
    wd2b -= wd2b_offset;
    --covi;
    
// cout << "z_offset=" << z_offset << endl;
// for(int i=1; i<= *nv; ++i){
//   cout << "logcch_ = b[" << i << "]=" << b[i] << endl;
//   cout << "logcch_ = ivar[" << i << "]=" << ivar[i] << endl;
// }

    /* Function Body */

/*        INITIAL SETTINGS */

    rlikp = one;
    *ifault = 0;
    its = 0;
    if (*nv > *nvmax) {
	goto L21;
    }
    if (*nvmax * (*nvmax + 1) / 2 > *nvmax1) {
	goto L21;
    }
    if (*nmax + 1 > *nmax1) {
	goto L21;
    }
    ins[1] = 0;
    if (nca[1] + nct[1] > *nmax || nca[1] + nct[1] > *nimax) {
	goto L21;
    }
    if (*ns == 1) {
	goto L2;
    }
    i__1 = *ns;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (nca[i__] + nct[i__] > *nmax) {
	    goto L21;
	}
	i1 = i__ - 1;
	ins[i__] = nca[i1] + nct[i1] + ins[i1];
/* L1: */
    }
    if (ins[*ns] + nca[*ns] + nct[*ns] > *nimax) {
	goto L21;
    }

/*        CENTRE THE INDEPENDENT VARIABLES ABOUT THE MEAN OF THE */
/*        COVARIATES FOR THE CASES (C.F. S.HOWARDS COMMENT TO COX (1972)) */

L2:
    i__1 = *ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nca[i__] * nct[i__] == 0) {
	    goto L6;
	}
	m = nca[i__];
	n = m + nct[i__];
	fm = (doublereal)1. / (doublereal) m;
	i1 = ins[i__];
	i__2 = *nv;
	for (k = 1; k <= i__2; ++k) { // k counts the number of covariates
	    cntr[k + i__ * cntr_dim1] = zero;
	    k1 = ivar[k]; // k1 is simply the index of the first covariate included in analysis
	    j1 = i1; // initially j1 equals the current stratus value of ins
//  cout << "i1=" << i1 << endl;
	    i__3 = m; // m is number of affected individuals in the stratum
	    for (j = 1; j <= i__3; ++j) {
		++j1;
//  cout << "j1=" << j1 << " k1=" << k1 << " z is " << z__[k1+j1*z_dim1] << endl;
		cntr[k + i__ * cntr_dim1] += z__[k1 + j1 * z_dim1];
/* L3: */
	    }
	    cntr[k + i__ * cntr_dim1] *= fm;
	    j1 = i1;
	    i__3 = n;
	    for (j = 1; j <= i__3; ++j) {
		++j1;
		z__[k1 + j1 * z_dim1] -= cntr[k + i__ * cntr_dim1];
//  cout << "z value when index= " << k1 + j1 * z_dim1 << " = " << z__[k1 + j1 * z_dim1] << endl;
/* L4: */
	    }
/* L5: */
	}
L6:
	;
    }
L7:
    ++its;
// cout << "its=" << its << endl;
    if (its > maxit) {
// cout << "its=" << its << endl;
	goto L23;
    }
    rlik = zero;
    k = 0;
    i__1 = *nv;
    for (j = 1; j <= i__1; ++j) {
	dl[j] = zero;
	i__2 = j;
	for (jj = 1; jj <= i__2; ++jj) {
	    ++k;
	    covi[k] = zero;
/* L8: */
	}
    }

/*        LOOP THROUGH STRATA */

    i__2 = *ns;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (nca[i__] * nct[i__] == 0) {
	    goto L14;
	}
	m = nca[i__];
	n = m + nct[i__];
	nid = ins[i__];

/*        FIND U(J)=EXP(Z(J)*BETA) FOR EACH INDIVIDUAL J IN THE STRATA. */
/*        ALSO, FIHD TTT, THE TOTAL OF THE (Z(J)*BETA)S. */

	ttt = zero;
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    j1 = j + nid;
	    t = zero;
// cout << "t should be zero " << endl;
	    i__3 = *nv;
// cout << "i__3=" << i__3 << endl;
	    for (k = 1; k <= i__3; ++k) {
		l = ivar[k];
// cout << "k=" << k << " l=" << ivar[k] << endl;
// cout << "b[k]=" << b[k] << endl;
// cout << "z__="<< z__[l + j1 * z_dim1] << endl;
		t += b[k] * z__[l + j1 * z_dim1];
/* L9: */
	    }
// cout << "t=" << t << endl;
	    ttt += t;
	    u[j] = exp(t);
/* L10: */
	}

/*        CALC LATE THE CONSTANT CONST= */
/*        ((N(C)M) * EXP(M*BETA*XBAR))**(1/M) */
/*        WHERE XBAR IS THE MEAN OF THE COVARIATES OVER THE CASES AND */
/*        CONTROLS, AHD DIVIDE EACH U(J) BY THIS CONSTANT. */
/*        THIS KEEPS THE SUMS CALCULATED BY SUBROUTINE HOWARD FROM */
/*        BECOMING TOO LARGE IN ABSOLUTE VALUE. */
/*        NOTE - ALGFAC(X)=LN((X)FACTORIAL) */

	i__1 = n - m;
// cout << "c1 before algfac=" << c1 << endl;
// cout << "ttt=" << ttt << " algfac(&n) =" << algfac_(&n) << endl;
// cout << "m=" << m << " algfac(&m) = " << algfac_(&m) << endl;
// cout << "i__1=" << i__1 << " algfac(&i__1) = " << algfac_(&i__1) << endl;
	c1 = algfac_(&n) - (algfac_(&m) + algfac_(&i__1)) + ttt * (doublereal) m / (
		doublereal) n;
// cout << "c1 after algfac=" << c1 << " n=" << n << endl;
	const__ = exp(-c1 / (doublereal) m);
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
/* L11: */
	    u[j] *= const__;
	}

/*        CALL TO HOWARD TO CALCULATE SUM(EXP(S(L)*BETA)) OVER ALL */
/*        COMBINATIONS OF N LABELS TAKEN M AT A TIME, AND ITS FIRST AND */
/*        SECOND DERIVATIVES WITH RESPECT TO BETA. */
/*        NOTE - THE VALUE FOR THE SUM AND ITS DERIVATIVES RETURNED BY */
/*        HOWARD ARE THE TRUE VALUES DIVIDED BY THE CONSTANT CONST**M. */
/*        THEREFORE RLIK IS CORRECTED FOR THIS CONSTANT SO THAT THE */
/*        LIKELIHOOD RATIO TEST STATISTIC WILL BE CORRECT. */
/*        NOTE - SC=BETA*(SUM(Z(J)) OVER THE CASES) = 0 BY DEFINITION. */
// cout << "bmn=" << bmn << " c1=" << c1 << endl;
	howard_(&m, &n, &u[1], &z__[z_offset], &nid, &ivar[1], nvmax, nimax, 
		nmax, nvmax1, nv, nmax1, &wb[1], &wdb[wdb_offset], &wd2b[
		wd2b_offset], &db[1], &d2b[1], &bmn);
// cout << "rlik before resetting with info from call to howard = " <<rlik << endl;
// cout << "bmn=" << bmn << " c1=" << c1 << endl;
	rlik = rlik - log(bmn) - c1;
// cout << "rlik after resetting with info from call to howard = " << rlik << endl;
	l = 0;
	ir = 1;
	is = 0;

/*        CALCULATE THE CUMULATIVE SCORE UP TO THIS STRATUM. */

	i__1 = *nv;
	for (k = 1; k <= i__1; ++k) {
	    dl[k] -= db[k] / bmn;
	    i__3 = k;
	    for (kk = 1; kk <= i__3; ++kk) {
		++l;
		++is;
		if (is <= ir) {
		    goto L12;
		}
		++ir;
		is = 1;

/*       CALCULATE THE CUMULATIVE INFORMATION UP TO THIS STRATUM. */

/* error 1 : bill.bardsley@man.ac.uk corrected the next line 18/08/2002 */
/*  12 COVI(L) = COVI(L) + D2B(L) / BMN - DB(IR) */
L12:
// cout << "l=" << l << " covi[l]=" << covi[l] << endl;
		covi[l] = covi[l] + d2b[l] / bmn - db[ir] * db[is] / (bmn * 
			bmn);
// cout << "l=" << l << " covi[l]=" << covi[l] << endl;
/* L13: */
	    }
	}
L14:
	;
    }
    if (its == 1) {
	rliks = rlik;
    }

/*        CALCULATE THE INVERSE OF THE INFORMATION MATRIX */

/* error 2: bill.bardsley@man.ac.uk corrected the next line 18/08/2002 */
/*      CALL SYMINV(COVI, NV, COV, W, NULLTY, IFAULT, NVMAX1) */
// cout << "cov[1]=" << cov[1] << endl;
    syminv_(&covi[1], nv, nvmax1, &cov[1], &w[1], &nullty, ifault);
// cout << "cov[1]=" << cov[1] << endl;
    if (*ifault != 0) {
	goto L22;
    }

/*        CALCULATE NEW PARAMETER ESTIMATES */

    i__2 = *nv;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w[i__] = zero;
	i2 = i__ * (i__ - 1) / 2;
	i__3 = i__;
	for (j = 1; j <= i__3; ++j) {
	    k = i2 + j;
// cout << "k=" << k << " cov[k]=" << cov[k] << " j=" << j <<  " dl[j]=" << dl[j] << endl;
	    w[i__] += dl[j] * cov[k];
/* L15: */
	}
	i1 = i__ + 1;
	if (i1 > *nv) {
	    goto L17;
	}
	i__3 = *nv;
	for (k = i1; k <= i__3; ++k) {
	    j = k * (k - 1) / 2 + i__;
	    w[i__] += dl[k] * cov[j];
/* L16: */
	}
L17:
	;
    }
    i__2 = *nv;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L18: */
// cout << "b=>" << b[i__] << " w[" << i__ <<"]=" << w[i__] << endl;
	b[i__] += w[i__];
    }
    if (its != 1) {
	goto L20;
    }

/*        CALCULATE THE TEST SCORE */
    *st = zero;
    i__2 = *nv;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L19: */
	*st += w[i__] * dl[i__];
    }

/*        TEST FOR CONVERGENCE */

L20:
// cout << "rlik=" << rlik << " rliks=" << rliks << endl;
    rlik -= rliks;
//     if ((r__1 = rlikp - rlik, dabs(r__1)) <= eps) {
    r__1 = rlikp - rlik;
// cout << "rlikp=" << rlikp << " rlik=" << rlik << " r__1=" << r__1 << endl;
    if(dabs(r__1) <= eps){
    
	goto L24;
    }
    rlikp = rlik;
    goto L7;
L21:
    *ifault = 1;
    goto L25;
L22:
    *ifault = 2;
    goto L25;
L23:
    *ifault = 3;
    goto L25;
L24:
    *chi2 = two * rlik;

/*        RETURN MATRIX Z TO THE FORM IT WAS IN WHEN LOGCCH WAS CALLED */

L25:
    i__2 = *ns;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (nca[i__] * nct[i__] == 0) {
	    goto L27;
	}
	n = nca[i__] + nct[i__];
	i1 = ins[i__];
	i__3 = *nv;
	for (k = 1; k <= i__3; ++k) {
	    k1 = ivar[k];
	    j1 = i1;
	    i__1 = n;
	    for (j = 1; j <= i__1; ++j) {
		++j1;
		z__[k1 + j1 * z_dim1] += cntr[k + i__ * cntr_dim1];
/* L26: */
	    }
	}
L27:
	;
    }
    return 0;
} /* logcch_ */


/* Subroutine */ int howard_(integer *m, integer *n, doublereal *u, doublereal *z__, 
	integer *nid, integer *ivar, integer *nvmax, integer *nimax, integer *
	nmax, integer *nvmax1, integer *nv, integer *nmax1, doublereal *wb, doublereal *
	wdb, doublereal *wd2b, doublereal *db, doublereal *d2b, doublereal *bmn)
{
    /* Initialized data */

    static doublereal zero = (doublereal)0.;
    static doublereal one = (doublereal)1.;

    /* System generated locals */
    integer z_dim1, z_offset, wdb_dim1, wdb_offset, wd2b_dim1, wd2b_offset, 
	    i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, i1, j1, l1, m1, im, ll, ir, is, ir1, is1, nv1, 
	    nmmp1;


/*        ALGORITHM AS 122.1 APPL. STATIST. (1984) VOL.33, NO.1 */



    /* Parameter adjustments */
    --db;
    --ivar;
    z_dim1 = *nvmax;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --u;
    --d2b;
    wd2b_dim1 = *nvmax1;
    wd2b_offset = 1 + wd2b_dim1;
    wd2b -= wd2b_offset;
    wdb_dim1 = *nvmax;
    wdb_offset = 1 + wdb_dim1;
    wdb -= wdb_offset;
    --wb;

    /* Function Body */

    nv1 = *nv * (*nv + 1) / 2;

/*        INITIALISE THE REQUIRED MATRICES AND VECTORS FOR THE RECURSION */

    wb[1] = one;
    i__ = 0;
    i__1 = *nv;
    for (l = 1; l <= i__1; ++l) {
	wdb[l + wdb_dim1] = zero;
	i__2 = l;
	for (ll = 1; ll <= i__2; ++ll) {
	    ++i__;
	    wd2b[i__ + wd2b_dim1] = zero;
/* L1: */
	}
    }
    m1 = *m + 1;
    if (*m == 0) {
	goto L9;
    }
    i__2 = m1;
    for (j = 2; j <= i__2; ++j) {
	wb[j] = zero;
	i__ = 0;
	i__1 = *nv;
	for (l = 1; l <= i__1; ++l) {
	    wdb[l + j * wdb_dim1] = zero;
	    i__3 = l;
	    for (ll = 1; ll <= i__3; ++ll) {
		++i__;
		wd2b[i__ + j * wd2b_dim1] = zero;
/* L2: */
	    }
	}
/* L3: */
    }
    nmmp1 = *n - *m + 1;
    i__2 = nmmp1;
    for (im = 1; im <= i__2; ++im) {
	i__3 = *m;
	for (j = 1; j <= i__3; ++j) {
	    j1 = j + 1;
	    i__ = j + im - 1;
	    i1 = i__ + *nid;
	    ir = 1;
	    is = 0;
	    i__1 = nv1;
	    for (l = 1; l <= i__1; ++l) {
		++is;
		is1 = ivar[is];
/* error 3: bill.bardsley@man.ac.uk corrected the next line 18/08/2002 */
/*     IR1 = IVAR(IS) */
		ir1 = ivar[ir];
		if (is <= ir) {
		    goto L4;
		}
		++ir;
		is = 1;
		is1 = ivar[is];
		ir1 = ivar[ir];

/*        CALCULATE SUM(EXP(S(L)*BETA)) OVER ALL COMBINATIONS OF N */
/*        LABELS TAKEN M AT A TIME?  AND ITS FIRST AND SECOND DERIVATIVES */
/*        WITH RESPECT TO BETA. */
/*        NOTE - THE VALUES RETURNED BT THIS ROUTINE ARE THE TRUE VALUES */
/*        OF THE SUMS DIVIDED BY TME CONSTANT CONST**M, WHERE CONST */
/*        IS AS CALCULATED IN 2UBROUTINE LOGCCH. */

/* error 4: bill.bardsley@man.ac.uk corrected the next line 18/08/2002 */
/*   4 WD2B(L, I1) = WD2B(L, J1) + U(I) * WD2B(L, J) + Z(IR1, I1) * */
L4:
		wd2b[l + j1 * wd2b_dim1] = wd2b[l + j1 * wd2b_dim1] + u[i__] *
			 wd2b[l + j * wd2b_dim1] + z__[ir1 + i1 * z_dim1] * 
			z__[is1 + i1 * z_dim1] * u[i__] * wb[j] + z__[ir1 + 
			i1 * z_dim1] * u[i__] * wdb[is + j * wdb_dim1] + z__[
			is1 + i1 * z_dim1] * u[i__] * wdb[ir + j * wdb_dim1];
/* L5: */
	    }
	    i__1 = *nv;
	    for (l = 1; l <= i__1; ++l) {
		l1 = ivar[l];
		wdb[l + j1 * wdb_dim1] = wdb[l + j1 * wdb_dim1] + u[i__] * 
			wdb[l + j * wdb_dim1] + z__[l1 + i1 * z_dim1] * u[i__]
			 * wb[j];
/* L6: */
	    }
	    wb[j1] += u[i__] * wb[j];
/* L7: */
	}
/* L8: */
    }
L9:
    *bmn = wb[m1];
    i__ = 0;
    i__2 = *nv;
    for (l = 1; l <= i__2; ++l) {
	db[l] = wdb[l + m1 * wdb_dim1];
	i__3 = l;
	for (ll = 1; ll <= i__3; ++ll) {
	    ++i__;
	    d2b[i__] = wd2b[i__ + m1 * wd2b_dim1];
/* L10: */
	}
    }
    return 0;
} /* howard_ */

double algfac_(integer *i__)
{
    /* Initialized data */

    static doublereal half = (doublereal).5;
    static doublereal const__ = (doublereal).9189385333;
    static doublereal one = (doublereal)1.;
    static doublereal twelve = (doublereal)12.;
    static doublereal tree60 = (doublereal)360.;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static doublereal b;
    static integer k;
    static doublereal aa, ak;


/*        ALGORITHM AS 196.2 APPL. STATIST. (1984) VOL.33, NO.1 */

/*        EVALUATES THE NATURAL LOGARITHM OF GAMMA(I+1)=LN(I-FACTORIAL) */
/*        (FORTRAN VERSION OF ACM291 - M.C.PIKE AND I.D.HILL,CACM,9,1966) */



    if (*i__ > 7) {
	goto L3;
    }
    b = one;
    if (*i__ <= 1) {
	goto L2;
    }
    i__1 = *i__;
    for (k = 2; k <= i__1; ++k) {
/* L1: */
	b *= (doublereal) k;
    }
L2:
    ret_val = log(b);
    return ret_val;
L3:
    aa = (doublereal) (*i__);
    ak = one / aa;
    ret_val = (aa + half) * log(aa) + ak * (one / twelve - ak * ak / tree60) 
	    - aa + const__;
    return ret_val;
} /* algfac_ */



// chol_(&a[1], &nrow, nvmax1, &c__[1], nullty, ifault);

/* Subroutine */ int chol_(doublereal *a, integer *n, integer *nvmax1, doublereal *u, 
	integer *nullty, integer *ifault)
{
    /* Initialized data */

    static real eta = 1e-9f;

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real w;
    static integer icol, irow;


/*       ALGORITHM AS 6 J,R.ST,ATIST.SOC. */

/*       GIVEN A SYMMETRIC MATRIX ORDER N */
/*       CALCULATES AN UPPER T P I A N G L E , U ( */
/*       U( ) MAY C O I N C I D E WITH A( ) . A ( */
/*       D- -F.F -T N. -I T F- .. */
/*       ETA I S SET TO M U L T I P L Y I N G FACTOR */
/*       ZERO FOR P I V O T . */
/*       I ( 1 9 6 8 ) V O L . 1 7 , N O t 2 . */
/*       AS LOWER TRIANGLE I N A( */
/*       ) , SUCH THAT UPRIME*U = A. */
/*       ) MUST RE P O S I T I V E S E M I - */
/*       DETERMINING E F F E C T I V E */
/*       NULLTY I S RETURNED AS NO. OF F F F E C T I V E ZERO P I V O T S . */
/*       I F A U L T I S RETURNED AS 1 I F N,I.E.O. 2 I F A( ) I S NOT P O S I T I V E */
/*       S E M I - D E F I N I T E W I T H I N THE. TOLERANCE D E F I N E D BY ETA, OTHERWISE */
/*       ZERO. */


/*      DIMENSION A(1),U(1) */
/*      Changed according to Remark AS R12 in applied statistics */
/* **** */
    /* Parameter adjustments */
    --u;
    --a;

    /* Function Body */

/*      THE VALUE OF ETA W I L L DEPEND ON THE WORD LENGTH OF THE */
/*       COMPUTER BEING USED */


/* **** */
    *ifault = 1;
    if (*n <= 0) {
	goto L100;
    }
    *ifault = 2;
    *nullty = 0;
    j = 1;
    k = 0;
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	l = 0;
	i__2 = icol;
	for (irow = 1; irow <= i__2; ++irow) {
	    ++k;
	    w = a[k];
	    m = j;
	    i__3 = irow;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		l += 1;
		if (i__ == irow) {
		    goto L13;
		}
		w -= u[l] * u[m];
		++m;
/* L12: */
	    }
L13:
	    if (irow == icol) {
		goto L14;
	    }
	    if (u[l] == 0.f) {
		goto L21;
	    }
	    u[k] = w / u[l];
	    goto L11;
L21:
	    u[k] = 0.f;
	    if (fabs(w) > (r__1 = eta * a[k], fabs(r__1))) {
		goto L100;
	    }
L11:
	    ;
	}
L14:
	if (fabs(w) < (r__1 = eta * a[k], fabs(r__1))) {
	    goto L20;
	}
	if (w < 0.f) {
	    goto L100;
	}
	u[k] = sqrt(w);
	goto L15;
L20:
	u[k] = 0.f;
	++(*nullty);
L15:
	j += icol;
/* L10: */
    }
    *ifault = 0;
L100:
    return 0;
} /* chol_ */

//    syminv_(&covi[1], nv, nvmax1, &cov[1], &w[1], &nullty, ifault);
/* Subroutine */ int syminv_(doublereal *a, integer *n, integer *nvmax1, doublereal *c__, 
	doublereal *w, integer *nullty, integer *ifault)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, l;
    static real x;
    static integer nn;
//    extern /* Subroutine */ int chol_(real *, integer *, integer *, real *, 
//	    integer *, integer *);
    static integer icol, jcol, irow, nrow, mdiag, ndiag;


/*       ALGORITHM AS 7 J.R,STATIST.SOC. C I (1968) VOL.17, N0.2. */

/*       FORMS IN C( AS LOWER TRIANGLE, A GENERALISED INVERSE */
/*       OF THE POSITIVE SEMI-DEFINITE SYMMETRIC MATRIX A( ) */
/*       ORDER N> STORED AS LOWER TRIANGLE. */
/*       C( ) MAY COINCIDE WITH A( ). NULLTY IS RETURNED AS THE NULLITY */
/*       OF A( ) , IFAULT IS RETURNED AS 1 IF NQLT.1, OTHERWISE ZERO. */
/*       W( ) IS A WORK ARRAY OF LENGTH AT LEAST N THAT IS ALLOCATED B Y */
/*       THE CALLING ROUTINE. */

/*      Below is original declaration */
/*      DIMENSION A(1),C(1),W(1) */
/*      Changed according to Remark AS R12 in applied statistics */

    /* Parameter adjustments */
    --w;
    --c__;
    --a;

    /* Function Body */
    nrow = *n;
    *ifault = 1;
    if ((real) nrow <=0) {
	goto L100;
    }
    *ifault = 0;
// cout << "a[1] " << a[1] << endl;
// cout << "calling chol_" << endl;
// cout << "c__[1]=" << c__[1] << endl;
    chol_(&a[1], &nrow, nvmax1, &c__[1], nullty, ifault);
// cout << "c__[1]=" << c__[1] << endl;
    if (*ifault != 0) {
	goto L100;
    }
    nn = nrow * (nrow + 1) / 2;
    irow = nrow;
    ndiag = nn;
L16:
    if (c__[ndiag] == 0.f) {
	goto L11;
    }
    l = ndiag;
    i__1 = nrow;
    for (i__ = irow; i__ <= i__1; ++i__) {
	w[i__] = c__[l];
	l += i__;
/* L10: */
    }
    icol = nrow;
    jcol = nn;
    mdiag = nn;
L15:
    l = jcol;
    x = 0.f;
    if (icol == irow) {
	x = 1.f / w[irow];
    }
    k = nrow;
L13:
    if (k == irow) {
	goto L12;
    }
    x -= w[k] * c__[l];
    --k;
    --l;
    if (l > mdiag) {
	l = l - k + 1;
    }
    goto L13;
L12:
    c__[l] = x / w[irow];
    if (icol == irow) {
	goto L14;
    }
    mdiag -= icol;
    --icol;
    --jcol;
    goto L15;
L11:
    l = ndiag;
    i__1 = nrow;
    for (j = irow; j <= i__1; ++j) {
	c__[l] = 0.f;
	l += j;
/* L17: */
    }
L14:
    ndiag -= irow;
    --irow;
    if (irow != 0) {
	goto L16;
    }
L100:
    return 0;
} /* syminv_ */

}
