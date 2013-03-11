// earth.h: externs for earth.c
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// A copy of the GNU General Public License is available at
// http://www.r-project.org/Licenses

#ifndef EARTH_H
#define EARTH_H

#include <string.h>
#include "DataSet.h"
#include "StepOptions.h"

using namespace std;

namespace Methods{
static const double *pxGlobal;

class earth{

private:
	DataSet* data_set;

	string EARTH_VERSION;// change if you modify this file!
	double BX_TOL;
	double QR_TOL;
	double MIN_GRSQ;
	double ALMOST_ZERO;
	int    ONE;        // parameter for BLAS routines
	double POS_INF;
	int    MAX_DEGREE;


	// Global copies of some input parameters.  These stay constant for the entire MARS fit.
	int nTraceGlobal;        // copy of nTrace parameter
	int nMinSpanGlobal;      // copy of nMinSpan parameter


	//-----------------------------------------------------------------------------
	// These are malloced blocks.  They unfortunately have to be declared globally so
	// under R if the user interrupts we can free them using on.exit(.C("FreeR"))

	int    *xOrder;          // local to FindTerm
	bool   *WorkingSet;      // local to FindTerm and EvalSubsets
	double *xbx;             // local to FindTerm
	double *CovSx;           // local to FindTerm
	double *CovCol;          // local to FindTerm
	double *CovSy;           // local to FindTerm
	double *ycboSum;         // local to FindTerm
	double *bxOrth;          // local to ForwardPass
	double *yMean;           // local to ForwardPass
	double *Weights;         // local to ForwardPass and EvalSubsetsUsingXtx

	// Transposed and mean centered copy of bxOrth, for fast update in FindKnot.
	// It's faster because there is better data locality as iTerm increases, so
	// better L1 cache use.  This is used only if USE_BLAS is true.

	double *bxOrthCenteredT; // local to ForwardPass

	double *bxOrthMean;      // local to ForwardPass
	int  *nFactorsInTerm;    // local to Earth or ForwardPassR
	int  *nUses;             // local to Earth or ForwardPassR


	//-----------------------------------------------------------------------------
	// The BetaCache is used when searching for a new term pair, via FindTerm.
	// Most of the calculation for the orthogonal regression betas is repeated
	// with the same data, and thus we can save time by caching betas.
	//
	// iParent    is the term that forms the base for the new term
	// iPred      is the predictor for the new term
	// iOrthCol   is the column index in the bxOrth matrix

	double *BetaCacheGlobal; // [iOrthCol,iParent,iPred]
	                                // dim nPreds x nMaxTerms x nMaxTerms

	//FAST_MAX
	typedef struct tQueue {
	    int     iParent;            // parent term
	    double  RssDelta;
	    int     nTermsForRssDelta;  // number of terms when RssDelta was calculated
	    double  AgedRank;
	} tQueue;

	tQueue *Q;           // indexed on iTerm (this Q is used for queue updates)
	tQueue *SortedQ;     // indexed on iParent rank (this Q is used to get next iParent)
	int    nQMax;        // number of elements in Q

	StepOptions* options;

	bool checkForMissing(Sample* samp, vector<int> loci, vector<int> covs);

	double grsq;
	double rsq;
	string best_model;


public:
	earth();
	earth(DataSet*);
	~earth(){};

	double getGeneralRsq(){return grsq;}
	double getRsq(){return rsq;}
	string getBestModel(){return best_model;}
	void set_parameters(StepOptions* opts){options = opts;}
	void resetDataSet(DataSet* ds){data_set = ds;}
	void resetData(DataSet* ds){data_set = ds;}
	int runMe(void);

	void* malloc1(size_t size);

	void* calloc1(size_t num, size_t size);

	static int Compare(const void *p1, const void *p2);  // for qsort

	void Order(int sorted[],                   // out: vector with nx elements
	                  const double x[], const int nx);   // in: x is a vector with nx elems

	int* OrderArray(const double x[], const int nRows, const int nCols);

	int GetNbrUsedCols(const bool UsedCols[], const int nLen);

	int CopyUsedCols(double **pxUsed,                // out
	                    const double x[],                   // in: nCases x nCols
	                    const int nCases, const int nCols,  // in
	                    const bool UsedCols[]);              // in

	void PrintSummary(
	    const int    nMaxTerms,         // in
	    const int    nTerms,            // in: number of cols in bx, some may be unused
	    const int    nPreds,            // in: number of predictors
	    const int    nResp,             // in: number of cols in y
	    const bool   UsedCols[],        // in: specifies used colums in bx
	    const int    Dirs[],            // in
	    const double Cuts[],            // in
	    const double Betas[],           // in: if NULL will print zeroes
	    const int    nFactorsInTerm[]);  // in: number of hinge funcs in basis term

	void CalcDiags(
	    double Diags[],     // out: nCols x 1
	    const double R[],   // in: nCases x nCols, QR from prev call to dqrsl
	    const int nCases,   // in
	    const int nCols);    // in

	void Regress(
	    double       Betas[],       // out: nUsedCols * nResp, can be NULL
	    double       Residuals[],   // out: nCases * nResp, can be NULL
	    double       *pRss,         // out: RSS, summed over all nResp, can be NULL
	    double       Diags[],       // out: diags of inv(transpose(x) * x), can be NULL
	    int          *pnRank,       // out: nbr of indep cols in x
	    int          iPivots[],     // out: nCols, can be NULL
	    const double x[],           // in: nCases x nCols, must include intercept
	    const double y[],           // in: nCases x nResp
	    const double weights[],     // in: nCases x 1, can be NULL
	    const int    nCases,        // in: number of rows in x and in y
	    const int    nResp,         // in: number of cols in y
	    int          nCols,         // in: number of columns in x, some may not be used
	    const bool   UsedCols[]);    // in: specifies used columns in x

	void RegressAndFix(
	    double Betas[],         // out: nMaxTerms x nResp, can be NULL
	    double Residuals[],     // out: nCases x nResp, can be NULL
	    double Diags[],         // out: if !NULL set to diags of inv(transpose(bx) * bx)
	    bool   UsedCols[],      // io:  will remove cols if necessary, nMaxTerms x 1
	    const  double bx[],     // in:  nCases x nMaxTerms
	    const double y[],       // in:  nCases x nResp
	    const double weights[], // in: nCases x 1, can be NULL
	    const int nCases,       // in
	    const int nResp,        // in: number of cols in y
	    const int nTerms);       // in: number of cols in bx, some may not be used

	double Mean(const double x[], int n);

	double SumOfSquares(const double x[], const double mean, int n);

	double GetGcv(const int nTerms, // nbr basis terms including intercept
	                const int nCases, double Rss, const double Penalty);

	int GetMinSpan(int nCases, int nPreds, const double *bx,
	                             const int iTerm);


	int GetEndSpan(const int nCases, const int nPreds);

	bool GetNewFormFlag(const int iPred, const int iTerm,
	                        const int Dirs[], const bool UsedCols[],
	                        const int nTerms, const int nPreds, const int nMaxTerms);

	double GetCut(int iCase, const int iPred, const int nCases,
	                        const double x[], const int xOrder[]);

	void InitBetaCache(const bool UseBetaCache,
	                          const int nMaxTerms, const int nPreds);

	void FreeBetaCache(void);

	void OrthogResiduals(
	    double bxOrthCol[],     // out: nCases x 1      { bxOrth[,nTerms] }
	    const double y[],       // in:  nCases x nResp  { bx[,nTerms], xbx }
	    const double bxOrth[],  // in:  nTerms x nPreds { bxOrth }
	    const int nCases,       // in
	    const int nTerms,       // in: nTerms in model, i.e. number of used cols in bxOrth
	    const bool UsedTerms[], // in: UsedTerms[i] is true if col is used, unused cols ignored
	                            //     Following parameters are only for the beta cache
	    const int iParent,      // in: if >= 0, use BetaCacheGlobal {FindTerm iTerm, addTermP -1}
	    const int iPred,        // in: predictor index i.e. col index in input matrix x
	    const int nMaxTerms);    // in:

	void InitBxOrthCol(
	    double bxOrth[],         // io: col nTerms is changed, other cols not touched
	    double bxOrthCenteredT[],// io: kept in sync with bxOrth
	    double bxOrthMean[],     // io: element at nTerms is updated
	    bool   *pGoodCol,        // io: set to false if col sum-of-squares is under BX_TOL
	    const double *y,         // in: { AddCandLinTerm xbx, addTermPair bx[,nTerms] }
	    const int nTerms,        // in: column goes in at index nTerms, 0 is the intercept
	    const bool WorkingSet[], // in
	    const int nCases,        // in
	    const int nMaxTerms,     // in
	    const int iCacheTerm,    // in: if >= 0, use BetaCacheGlobal {FindTerm iTerm, AddTermP -1}
	                             //     if < 0 then recalc Betas from scratch
	    const int iPred,         // in: predictor index i.e. col index in input matrix x
	    const double weights[]);  // in:

	void AddTermPair(
	    int    Dirs[],              // io
	    double Cuts[],              // io
	    double bx[],                // io: MARS basis matrix
	    double bxOrth[],            // io
	    double bxOrthCenteredT[],   // io
	    double bxOrthMean[],        // io
	    bool   FullSet[],           // io
	    int    nFactorsInTerm[],    // io
	    int    nUses[],             // io: nbr of times each predictor is used in the model
	    const int nTerms,           // in: new term pair goes in at index nTerms and nTerms1
	    const int iBestParent,      // in: parent term
	    const int iBestCase,        // in
	    const int iBestPred,        // in
	    const int nPreds,           // in
	    const int nCases,           // in
	    const int nMaxTerms,        // in
	    const bool IsNewForm,       // in
	    const bool IsLinPred,       // in: pred was discovered by search to be linear
	    const int LinPreds[],       // in: user specified preds which must enter linearly
	    const double x[],           // in
	    const int xOrder[],         // in
	    const double weights[]);     // in:
	void FindKnot(
	    int    *piBestCase,         // out: possibly updated, index into the ORDERED x's
	    double *pRssDeltaForParentPredPair, // io: updated if knot is better
	    double CovCol[],            // scratch buffer, overwritten, nTerms x 1
	    double ycboSum[],           // scratch buffer, overwritten, nMaxTerms x nResp
	    double CovSx[],             // scratch buffer, overwritten, nTerms x 1
	    double *ybxSum,             // scratch buffer, overwritten, nResp x 1
	    const int nTerms,           // in
	    const int iParent,          // in: parent term
	    const int iPred,            // in: predictor index
	    const int nCases,           // in
	    const int nResp,            // in: number of cols in y
	    const int nMaxTerms,        // in
	    const double RssDeltaLin,   // in: change in RSS if predictor iPred enters linearly
	    const double MaxAllowedRssDelta, // in: FindKnot rejects any changes in Rss greater than this
	    const double bx[],          // in: MARS basis matrix
	    const double bxOrth[],      // in
	    const double bxOrthCenteredT[], // in
	    const double bxOrthMean[],  // in
	    const double x[],           // in: nCases x nPreds
	    const double y[],           // in: nCases x nResp
	    const double Weights[],     // in: nCases x 1, must not be NULL
	    const int xOrder[],         // in
	    const double yMean[],       // in: vector nResp x 1
	    const int nMinSpan,
	    const int nEndSpan,
	    const double NewVarAdjust);  // in: 1 if not a new var, 1+NewVarPenalty if new var

	void AddCandidateLinearTerm(
	    double *pRssDeltaLin,       // out: change to RSS caused by adding new term
	    bool   *pIsNewForm,         // io:
	    double xbx[],               // out: nCases x 1
	    double CovCol[],            // out: nMaxTerms x 1
	    double ycboSum[],           // io: nMaxTerms x nResp
	    double bxOrth[],            // io
	    double bxOrthCenteredT[],   // io
	    double bxOrthMean[],        // io
	    const int iPred,            // in
	    const int iParent,          // in
	    const double x[],           // in: nCases x nPreds
	    const double y[],           // in: nCases x nResp
	    const double Weights[],     // in: nCases x 1, must not be NULL
	    const int nCases,           // in
	    const int nResp,            // in: number of cols in y
	    const int nTerms,           // in
	    const int nMaxTerms,        // in
	    const double yMean[],       // in: vector nResp x 1
	    const double bx[],          // in: MARS basis matrix
	    const bool FullSet[],       // in
	    const double NewVarAdjust);  // in

	void FindPred(
	    int    *piBestCase,         // out: return -1 if no new term available
	                                //      else return an index into the ORDERED x's
	    int    *piBestPred,         // out
	    int    *piBestParent,       // out: existing term on which we are basing the new term
	    double *pBestRssDeltaForTerm,   // io: updated if new predictor is better
	    double *pBestRssDeltaForParent, // io: used only by FAST_MARS
	    bool   *pIsNewForm,         // out
	    bool   *pIsLinPred,         // out: true if knot is at min x val so x enters linearly
	    double MaxRssPerPred[],     // io: nPreds x 1, max RSS for each predictor over all parents
	    double xbx[],               // io: nCases x 1
	    double CovSx[],             // io
	    double CovCol[],            // io
	    double ycboSum[],           // io: nMaxTerms x nResp
	    double bxOrth[],            // io
	    double bxOrthCenteredT[],   // io
	    double bxOrthMean[],        // io
	    const int iBestPred,        // in: if -1 then search for best predictor, else use this predictor
	    const int iParent,          // in
	    const double x[],           // in: nCases x nPreds
	    const double y[],           // in: nCases x nResp
	    const double Weights[],     // in: nCases x 1
	    const int nCases,           // in
	    const int nResp,            // in: number of cols in y
	    const int nPreds,           // in
	    const int nTerms,           // in
	    const int nMaxTerms,        // in
	    const double yMean[],       // in: vector nResp x 1
	    const double MaxAllowedRssDelta, // in: FindKnot rejects any changes in Rss greater than this
	    const double bx[],          // in: MARS basis matrix
	    const bool FullSet[],       // in
	    const int xOrder[],         // in
	    const int nUses[],          // in: nbr of times each predictor is used in the model
	    const int Dirs[],           // in
	    const double NewVarPenalty, // in: penalty for adding a new variable (default is 0)
	    const int LinPreds[]);       // in: nPreds x 1, 1 if predictor must enter linearly

	void FindTerm(
	    int    *piBestCase,         // out: return -1 if no new term available
	                                //      else return an index into the ORDERED x's
	    int    *piBestPred,         // out:
	    int    *piBestParent,       // out: existing term on which we are basing the new term
	    double *pBestRssDeltaForTerm, // out: adding new term reduces RSS this much
	                                  //      will be set to 0 if no possible new term
	    bool   *pIsNewForm,         // out
	    bool   *pIsLinPred,         // out: true if knot is at min x val so x enters linearly
	    double MaxRssPerPred[],     // io: nPreds x 1, max RSS for each predictor over all parents
	    double bxOrth[],            // io: column nTerms overwritten
	    double bxOrthCenteredT[],   // io: kept in sync with bxOrth
	    double bxOrthMean[],        // io: element nTerms overwritten
	    const int iBestPred,        // in: if -1 then search for best predictor, else use this predictor
	    const double x[],           // in: nCases x nPreds
	    const double y[],           // in: nCases x nResp
	    const double Weights[],     // in: nCases x 1
	    const int nCases,           // in:
	    const int nResp,            // in: number of cols in y
	    const int nPreds,           // in:
	    const int nTerms,           // in:
	    const int nMaxDegree,       // in:
	    const int nMaxTerms,        // in:
	    const double yMean[],       // in: vector nResp x 1
	    const double MaxAllowedRssDelta, // in: FindKnot rejects any changes in Rss greater than this
	    const double bx[],          // in: MARS basis matrix
	    const bool FullSet[],       // in:
	    const int xOrder[],         // in:
	    const int nFactorsInTerm[], // in:
	    const int nUses[],          // in: nbr of times each predictor is used in the model
	    const int Dirs[],           // in:
	    const int nFastK,           // in: Fast MARS K
	    const double NewVarPenalty, // in: penalty for adding a new variable (default is 0)
	    const int LinPreds[]);       // in: nPreds x 1, 1 if predictor must enter linearly

	void PrintForwardProlog(const int nCases, const int nPreds,
	        const char *sPredNames[]);   // in: predictor names, can be NULL

	void PrintForwardStep(
	        const int nTerms,
	        const int nUsedTerms,
	        const int iBestCase,
	        const int iBestPred,
	        const double RSq,
	        const double RSqDelta,
	        const double Gcv,
	        const double GcvNull,
	        const int nCases,
	        const int xOrder[],
	        const double x[],
	        const bool IsLinPred,
	        const bool IsNewForm,
	        const char *sPredNames[]);   // in: predictor names, can be NULL

	void PrintForwardEpilog(
	            const int nTerms, const int nMaxTerms,
	            const double Thresh,
	            const double RSq, const double RSqDelta,
	            const double Gcv, const double GcvNull,
	            const int iBestCase,
	            const bool FullSet[]);

	void CheckVec(const double x[], int nCases, int nCols, const char sVecName[]);

	void PrintEstimatedMemoryUse(int nMaxTerms, int nCases, int nPreds, int nResp);

	void CheckRssNull(double RssNull, const double y[], int iResp, int nCases);

	double* pInitWeights(const double weightsArg[], int nCases);

	void ForwardPass(
	    int    *pnTerms,            // out: highest used term number in full model
	    bool   FullSet[],           // out: 1 * nMaxTerms, indices of lin indep cols of bx
	    double bx[],                // out: MARS basis matrix, nCases * nMaxTerms
	    int    Dirs[],              // out: nMaxTerms * nPreds, -1,0,1,2 for iTerm, iPred
	    double Cuts[],              // out: nMaxTerms * nPreds, cut for iTerm, iPred
	    int    nFactorsInTerm[],    // out: number of hockey stick funcs in each MARS term
	    int    nUses[],             // out: nbr of times each predictor is used in the model
	    const double x[],           // in: nCases x nPreds
	    const double y[],           // in: nCases x nResp
	    const double weightsArg[],  // in: nCases x 1, can be NULL, currently ignored
	    const int nCases,           // in: number of rows in x and elements in y
	    const int nResp,            // in: number of cols in y
	    const int nPreds,           // in:
	    const int nMaxDegree,       // in:
	    const int nMaxTerms,        // in:
	    const double Penalty,       // in: GCV penalty per knot
	    double Thresh,              // in: forward step threshold
	    int nFastK,                 // in: Fast MARS K
	    const double FastBeta,      // in: Fast MARS ageing coef
	    const double NewVarPenalty, // in: penalty for adding a new variable
	    const int LinPreds[],       // in: 1 x nPreds, 1 if predictor must enter linearly
	    const bool UseBetaCache,    // in: 1 to use the beta cache, for speed
	    const char *sPredNames[]);   // in: predictor names, can be NULL

	void EvalSubsetsUsingXtx(
	    bool   PruneTerms[],    // out: nMaxTerms x nMaxTerms
	    double RssVec[],        // out: nMaxTerms x 1, RSS of each subset
	    const int    nCases,    // in
	    const int    nResp,     // in: number of cols in y
	    const int    nMaxTerms, // in: number of MARS terms in full model
	    const double bx[],      // in: nCases x nMaxTerms, all cols must be indep
	    const double y[],       // in: nCases * nResp
	    const double weightsArg[]); // in: nCases x 1, can be NULL

	void BackwardPass(
	    double *pBestGcv,       // out: GCV of the best model i.e. BestSet columns of bx
	    bool   BestSet[],       // out: nMaxTerms x 1, indices of best set of cols of bx
	    double Residuals[],     // out: nCases x nResp
	    double Betas[],         // out: nMaxTerms x nResp
	    const double bx[],      // in: nCases x nMaxTerms
	    const double y[],       // in: nCases x nResp
	    const double weightsArg[], // in: nCases x 1, can be NILL
	    const int nCases,       // in: number of rows in bx and elements in y
	    const int nResp,        // in: number of cols in y
	    const int nMaxTerms,    // in: number of cols in bx
	    const double Penalty);   // in: GCV penalty per knot

	int DiscardUnusedTerms(
	    double bx[],             // io: nCases x nMaxTerms
	    int    Dirs[],           // io: nMaxTerms x nPreds
	    double Cuts[],           // io: nMaxTerms x nPreds
	    bool   WhichSet[],       // io: tells us which terms to discard
	    int    nFactorsInTerm[], // io
	    const int nMaxTerms,
	    const int nPreds,
	    const int nCases);

	int GetMaxKnotsPerTerm(
	    const bool   UsedCols[],    // in
	    const int    Dirs[],        // in
	    const int    nPreds,        // in
	    const int    nTerms,        // in
	    const int    nMaxTerms);     // in

	void FormatOneResponse(
	    const bool   UsedCols[],// in: nMaxTerms x 1, indices of best set of cols of bx
	    const int    Dirs[],    // in: nMaxTerms x nPreds, -1,0,1,2 for iTerm, iPred
	    const double Cuts[],    // in: nMaxTerms x nPreds, cut for iTerm, iPred
	    const double Betas[],   // in: nMaxTerms x nResp
	    const int    nPreds,
	    const int    iResp,
	    const int    nTerms,
	    const int    nMaxTerms,
	    const int    nDigits,   // number of significant digits to print
	    const double MinBeta);   // terms with fabs(betas) less than this are not printed, 0 for all


	double PredictOneResponse(
	    const double x[],        // in: vector nPreds x 1 of input values
	    const bool   UsedCols[], // in: nMaxTerms x 1, indices of best set of cols of bx
	    const int    Dirs[],     // in: nMaxTerms x nPreds, -1,0,1,2 for iTerm, iPred
	    const double Cuts[],     // in: nMaxTerms x nPreds, cut for iTerm, iPred
	    const double Betas[],    // in: nMaxTerms x 1
	    const int    nPreds,     // in: number of cols in x
	    const int    nTerms,
	    const int    nMaxTerms);


void Earth(
    double *pBestGcv,       // out: GCV of the best model i.e. BestSet columns of bx
    int    *pnTerms,        // out: max term nbr in final model, after removing lin dep terms
    bool   BestSet[],       // out: nMaxTerms x 1, indices of best set of cols of bx
    double bx[],            // out: nCases x nMaxTerms
    int    Dirs[],          // out: nMaxTerms x nPreds, 1,0,-1 for term iTerm, predictor iPred
    double Cuts[],          // out: nMaxTerms x nPreds, cut for term iTerm, predictor iPred
    double Residuals[],     // out: nCases x nResp
    double Betas[],         // out: nMaxTerms x nResp
    const double x[],       // in: nCases x nPreds
    const double y[],       // in: nCases x nResp
    const double weightsArg[], // in: nCases, can be NULL
    const int nCases,       // in: number of rows in x and elements in y
    const int nResp,        // in: number of cols in y
    const int nPreds,       // in: number of cols in x
    const int nMaxDegree,   // in: Friedman's mi
    const int nMaxTerms,    // in: includes the intercept term
    const double Penalty,   // in: GCV penalty per knot
    double Thresh,          // in: forward step threshold
    const int nMinSpan,     // in: set to non zero to override internal calculation
    const bool Prune,       // in: do backward pass
    const int nFastK,       // in: Fast MARS K
    const double FastBeta,  // in: Fast MARS ageing coef
    const double NewVarPenalty, // in: penalty for adding a new variable
    const int LinPreds[],       // in: 1 x nPreds, 1 if predictor must enter linearly
    const bool UseBetaCache,    // in: 1 to use the beta cache, for speed
    const int nTrace,           // in: 0 none 1 overview 2 forward 3 pruning 4 more pruning
    const char *sPredNames[]);  // in: predictor names in trace printfs, can be NULL

void FormatEarth(
    const bool   UsedCols[],// in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],    // in: nMaxTerms x nPreds, 1,0,-1 for term iTerm, predictor iPred
    const double Cuts[],    // in: nMaxTerms x nPreds, cut for term iTerm, predictor iPred
    const double Betas[],   // in: nMaxTerms x nResp
    const int    nPreds,
    const int    nResp,     // in: number of cols in y
    const int    nTerms,
    const int    nMaxTerms,
    const int    nDigits,   // number of significant digits to print
    const double MinBeta);  // terms with abs(beta) less than this are not printed, 0 for all

void PredictEarth(
    double       y[],           // out: vector nResp
    const double x[],           // in: vector nPreds x 1 of input values
    const bool   UsedCols[],    // in: nMaxTerms x 1, indices of best set of cols of bx
    const int    Dirs[],        // in: nMaxTerms x nPreds, 1,0,-1 for iTerm iPred
    const double Cuts[],        // in: nMaxTerms x nPreds, cut for term iTerm predictor iPred
    const double Betas[],       // in: nMaxTerms x nResp
    const int    nPreds,        // in: number of cols in x
    const int    nResp,         // in: number of cols in y
    const int    nTerms,
    const int    nMaxTerms);


    void error(const char *args, ...);

    void calculate(vector<int>, vector<int>);
    void calculate(vector<int> v){vector<int> temp; calculate(v, temp);}

void FreeQ(void);
void InitQ(const int nMaxTerms);
void PrintSortedQ(int nFastK);     // for debugging
static int CompareQ(const void *p1, const void *p2);     // for qsort
static int CompareAgedQ(const void *p1, const void *p2); // for qsort
void AddTermToQ(
    const int iTerm,        // in
    const int nTerms,       // in
    const double RssDelta,  // in
    const bool Sort,        // in
    const int nMaxTerms,    // in
    const double FastBeta);  // in: ageing Coef, 0 is no ageing, FastMARS recommends 1
void UpdateRssDeltaInQ(const int iParent, const int nTermsForRssDelta,
                              const double RssDelta);

int GetNextParent(   // returns -1 if no more parents
    const bool InitFlag,    // use true to init, thereafter false
    const int  nFastK);

};
};

#endif // EARTH_H
