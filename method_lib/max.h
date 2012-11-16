/* NEWLIP */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include <string>

namespace Methods{
#define MAXINDFAM   300
#define MAXFAM	    3001/*200*/		
#define MAXALL	    35
#define MAXLEN	    250
#define MAXNUC	    50
#define YES	    1  
#define NO	    0 
#define BOTH	    0
#define MALE	    1
#define FEMALE	    2

typedef struct {
	string id;//char	   *id;
	int        alle1; 
	int        alle2; 
	} Genotype;

typedef struct {
	int	   famnum; /* family ID number */
	int	   type; /* 1=single affected, 2=two affecteds */
	Genotype   P[3]; /* parents */
	Genotype   C[3]; /* children */
	} FFamily;

typedef struct {
	int        sing;  /* number of single affected families */
	int        mult;  /* number of multiply affected families */
	vector<FFamily> fam;//FFamily     *fam[MAXFAM];
	} Familyfile;

typedef struct {
	string disease;
	//	char	disease[15];	/* name of disease locus */
	int	liabclass;	/* is there liability class */
	string marker;
	//	char	marker[20];	/* name of the marker locus */
	int	nall;		/* num alleles at marker locus */
	int	BN;		/* band names available Y/N? */
	vector<string> band;//char	*band[MAXALL];	/* band names, if available */
	} Loci;

typedef struct {
	int	numalle; /* number of alleles in this marker */
	
/* for triads: */
	/* # of parents with genotype MiMj */
	int 	hy[MAXALL][MAXALL];
	/* sum of all heterozygous parents of triads with allele i */
	int	hyt[MAXALL];
	/* # of parents who transmitted alle i and not j to aff kid */
	int	syij[MAXALL][MAXALL];

/* for sibpairs */
	/* total number of heter parents of sibpairs */
	int	hxall;
	/* # of parents with genotype MiMj */
	int	hx[MAXALL][MAXALL];
	/* total of heter paretns of sibpairs (sum over i != j) */
	int	hxt[MAXALL];
	/* # of parents with genotype MiMj who transmit i to both*/
	int	sxstar;
	/* # parents who transmit alle i and not j to both aff kids */
	int	sxij[MAXALL][MAXALL];

	double PVAL;
	double CHI;
	} Totals;

typedef struct {
	int fam,id,rel[7],aff,st,alle1,alle2;
	string idnum;//char idnum[8];
	} Person;
/* rel[0]=fa  1=ma  2=offs  3=nfa  4=nma  5=sex  6=pro */

/*struct nlist {
     char  *cind;
     char  *cid;
     struct nlist *next;
     } *np;
*/
typedef struct {
     char  *cind;
     char  *cid;
     void *next;
     } Np;//*np;

typedef struct {
	vector<string> ind;//char *ind[MAXINDFAM];
	} Nuclear;

//static char *SEXLIST[] = {"both parents","fathers only","mothers only"};

//char lineno[256];
};
#include "proto.h"
