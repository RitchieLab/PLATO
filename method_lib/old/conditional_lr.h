//conditional_lr.h
namespace Methods{
typedef long integer;
typedef unsigned long int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;

#include <math.h>

#define dabs(x) (double)fabs(x)

int logcch_(integer *ns, integer *nca, integer *nct, integer *nimax, integer *nmax, integer *nmax1, integer *nvmax, integer *nvmax1, integer *nv, doublereal *z__, integer *ivar, doublereal *covi, doublereal *cntr, doublereal *w, doublereal *wb, doublereal *wdb, doublereal *wd2b, doublereal *u, integer *ins, doublereal *db, doublereal *d2b, doublereal *dl, doublereal *b, doublereal *cov, doublereal *chi2, doublereal *st, integer *ifault);
int howard_(integer *m, integer *n, doublereal *u, doublereal *z__, integer *nid, integer *ivar, integer *nvmax, integer *nimax, integer *nmax, integer *nvmax1, integer *nv, integer *nmax1, doublereal *wb, doublereal *wdb, doublereal *wd2b, doublereal *db, doublereal *d2b, doublereal *bmn);
double algfac_(integer *i__);
int chol_(doublereal *a, integer *n, integer *nvmax1, doublereal *u, integer *nullty, integer *ifault);
int syminv_(doublereal *a, integer *n, integer *nvmax1, doublereal *c__, doublereal *w, integer *nullty, integer *ifault);
};
