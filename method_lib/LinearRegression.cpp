//Liberally pulled from Plink and adapted to work in the Plato environment

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2006 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <boost/math/distributions/fisher_f.hpp>
#include "LinearRegression.h"
#include "Helpers.h"
#include "Options.h"

using namespace boost::math;

///#include "stats.h"
namespace Methods {
void LinearRegression::display(vector<vector<double> > & m) {
	cout << "\n";
	for (int i = 0; i < (int) m.size(); i++) {
		cout << i << ")\t";
		for (int j = 0; j < (int) m[i].size(); j++)
			cout << m[i][j] << " ";
		cout << "\n";
	}
	cout << "\n";
}

void LinearRegression::display(vector<double> & m) {
	cout << "\n";
	for (int i = 0; i < (int) m.size(); i++)
		cout << i << ")\t" << m[i] << "\n";
	cout << "\n";
	cout << "\n";
}

void LinearRegression::display(vector<int> & m) {
	cout << "\n";
	for (int i = 0; i < (int) m.size(); i++)
		cout << i << ")\t" << m[i] << "\n";
	cout << "\n";
	cout << "\n";
}

LinearRegression::LinearRegression()//Plink * p_)
{
	//  P = p_;

	nc = 0;
	cluster = false;
	RSS = -1;
}

void LinearRegression::setDependent() {
	// Set dependent variable and intercept
	Y.clear();

	for (int i = 0; i < data_set->num_inds(); i++)
		//P->n; i++)
		if (!miss[i]) {
			//    	if(options.getTraits().size() > 0){
			if (options.getUsePheno())
			{
				int index = options.getPhenoLoc();

				if (options.getPhenoName() != "") {
					index = data_set->get_trait_index(options.getPhenoName());
				}
				if (index < 0) {
					throw MethodException(
							"Internal Error: Trait/Phenotype index value < 0 in Linear Regression!");
				}
				Y.push_back(data_set->get_sample(i)->getPheno(index));//data_set->get_trait_index(options.getTraits().at(0))));//P->sample[i]->pperson->phenotype ) ;
			}
			else
			{
				Y.push_back(data_set->get_sample(i)->getPheno());//P->sample[i]->pperson->phenotype ) ;
			}
		}
}

void LinearRegression::pruneY() {

	//////////////////////////////////
	// Prune out rows that are missing

	vector<double> Y2;
	for (int i = 0; i < (int) Y.size(); i++)
		if (!miss[i])
			Y2.push_back(Y[i]);
	Y = Y2;
}

void covsrt(vector<vector<double> > & covar, vector<bool> &ia, const int mfit) {
	int i, j, k;

	int ma = ia.size();
	for (i = mfit; i < ma; i++)
		for (j = 0; j < i + 1; j++)
			covar[i][j] = covar[j][i] = 0.0;
	k = mfit - 1;
	for (j = ma - 1; j >= 0; j--) {
		if (ia[j]) {
			for (i = 0; i < ma; i++)
				SWAP(covar[i][k], covar[i][j]);
			for (i = 0; i < ma; i++)
				SWAP(covar[k][i], covar[j][i]);
			k--;
		}
	}
}

void gaussj(vector<vector<double> > & a, vector<vector<double> > & b) {
	int i, icol = 0, irow = 0, j, k, l, ll;
	double big, dum, pivinv;

	int n = a.size();
	int m = b[0].size();
	vector<double> indxc(n), indxr(n), ipiv(n);
	for (j = 0; j < n; j++)
		ipiv[j] = 0;
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++)
			if (ipiv[j] != 1)
				for (k = 0; k < n; k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l = 0; l < n; l++)
				SWAP(a[irow][l], a[icol][l]);
			for (l = 0; l < m; l++)
				SWAP(b[irow][l], b[icol][l]);
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0)
			throw MethodException("gaussj: Singular Matrix");
		pivinv = 1.0 / a[icol][icol];
		a[icol][icol] = 1.0;
		for (l = 0; l < n; l++)
			a[icol][l] *= pivinv;
		for (l = 0; l < m; l++)
			b[icol][l] *= pivinv;
		for (ll = 0; ll < n; ll++)
			if (ll != icol) {
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l = 0; l < n; l++)
					a[ll][l] -= a[icol][l] * dum;
				for (l = 0; l < m; l++)
					b[ll][l] -= b[icol][l] * dum;
			}
	}
	for (l = n - 1; l >= 0; l--) {
		if (indxr[l] != indxc[l])
			for (k = 0; k < n; k++)
				SWAP(a[k][(int) indxr[l]], a[k][(int) indxc[l]]);
	}
}

void lfit(vector<double> &x, vector<double> &y, vector<double> &sig, vector<
		double> &a, vector<bool> &ia, vector<vector<double> > &covar,
		double &chisq, vector<vector<double> > & X) {
	int i, j, k, l, m, mfit = 0;
	double ym, wt, sum, sig2i;

	int ndat = x.size();
	int ma = a.size();
	vector<double> afunc(ma);
	vector<vector<double> > beta;
	Helpers::sizeMatrix(beta, ma, 1);
	for (j = 0; j < ma; j++)
		if (ia[j])
			mfit++;
	if (mfit == 0)
		throw MethodException("lfit: no parameters to be fitted");
	for (j = 0; j < mfit; j++) {
		for (k = 0; k < mfit; k++)
			covar[j][k] = 0.0;
		beta[j][0] = 0.0;
	}
	for (i = 0; i < ndat; i++) {
		afunc = X[i];

		ym = y[i];
		if (mfit < ma) {
			for (j = 0; j < ma; j++)
				if (!ia[j])
					ym -= a[j] * afunc[j];
		}
		sig2i = 1.0 / SQR(sig[i]);
		for (j = 0, l = 0; l < ma; l++) {
			if (ia[l]) {
				wt = afunc[l] * sig2i;
				for (k = 0, m = 0; m <= l; m++)
					if (ia[m])
						covar[j][k++] += wt * afunc[m];
				beta[j++][0] += ym * wt;
			}
		}
	}
	for (j = 1; j < mfit; j++)
		for (k = 0; k < j; k++)
			covar[k][j] = covar[j][k];
	vector<vector<double> > temp;
	Helpers::sizeMatrix(temp, mfit, mfit);
	for (j = 0; j < mfit; j++)
		for (k = 0; k < mfit; k++)
			temp[j][k] = covar[j][k];
	gaussj(temp, beta);
	for (j = 0; j < mfit; j++)
		for (k = 0; k < mfit; k++)
			covar[j][k] = temp[j][k];
	for (j = 0, l = 0; l < ma; l++)
		if (ia[l])
			a[l] = beta[j++][0];
	chisq = 0.0;
	for (i = 0; i < ndat; i++) {
		afunc = X[i];
		sum = 0.0;
		for (j = 0; j < ma; j++)
			sum += a[j] * afunc[j];
		chisq += SQR((y[i] - sum) / sig[i]);
	}
	covsrt(covar, ia, mfit);
}

void LinearRegression::setCluster(vector<int> & cl) {
	cluster = true;
	clst = cl;
	nc = 0;
	if ((int) cl.size() != nind)
		throw MethodException("Setting LM cluster to wrong N -- internal error");
	map<int, int> novelc;
	int k = 0;
	for (int i = 0; i < nind; i++)
		if (novelc.find(clst[i]) == novelc.end())
			novelc.insert(make_pair(clst[i], k++));
	nc = novelc.size();
	if (nc < 2)
		throw MethodException("Must have at least two clusters");
	for (int i = 0; i < nind; i++)
		clst[i] = novelc.find(cl[i])->second;

}

void LinearRegression::noCluster() {
	cluster = false;
	clst.clear();
	nc = 0;
}

void LinearRegression::setVariance() {
	varY = 0;
	double mean = 0;
	int actualN = 0;

	for (int i = 0; i < nind; i++) {
		actualN++;
		mean += Y[i];
	}
	if (actualN == 0) {
		varY = 0;
		return;
	}
	mean /= (double) actualN;

	for (int i = 0; i < nind; i++) {
		varY += (Y[i] - mean) * (Y[i] - mean);
	}
	varY /= (double) (actualN - 1);

	if (actualN != nind)
		throw MethodException("actualN <> nind...");
}

void LinearRegression::fitLM() {
	bool verbose = false;
	if (verbose) {
		cout << "fitLM NP = " << np << endl;
		for (int i = 0; i < nind; i++) {
			cout << "VO " << i << "\t" << Y[i] << "\t";
			for (int j = 0; j < np; j++)
				cout << X[i][j] << "\t";
			cout << "\n";
		}
	}

	//cout << "Made it here A\n";
	//cout << "LM VIEW\n";
	//    display(Y);
	//    display(X);
	//    cout << "---\n";

	coef.resize(np);
	Helpers::sizeMatrix(S, np, np);

	//cout << "Made it here B\n";

	//cout << "np: " << getString<int>(np) << ", nind: " << getString<int>(nind) << ", all_valid: " << getString<bool>(all_valid) << "\n";
	if (np == 0 || nind == 0 || !all_valid) {
		return;
	}

	//cout << "Made it here C\n";

	setVariance();
	sig.resize(nind, sqrt(1.0 / sqrt((double) nind)));

	w.resize(np);
	Helpers::sizeMatrix(u, nind, np);
	Helpers::sizeMatrix(v, np, np);

	//cout << "Made it here D\n";

	//  Perform "svdfit(C,Y,sig,b,u,v,w,chisq,function)"

	int i, j;
	const double TOL = 1.0e-13;
	double wmax, tmp, thresh, sum;

	vector<double> b(nind), afunc(np);
	for (i = 0; i < nind; i++) {
		afunc = X[i];
		tmp = 1.0 / sig[i];
		for (j = 0; j < np; j++)
			u[i][j] = afunc[j] * tmp;
		b[i] = Y[i] * tmp;
	}

	Helpers::svdcmp(u, w, v);

	wmax = 0.0;
	for (j = 0; j < np; j++)
		if (w[j] > wmax)
			wmax = w[j];
	thresh = TOL * wmax;
	for (j = 0; j < np; j++)
		if (w[j] < thresh)
			w[j] = 0.0;

	Helpers::svbksb(u, w, v, b, coef);

	chisq = 0.0;
	for (i = 0; i < nind; i++) {
		afunc = X[i];
		sum = 0.0;
		for (j = 0; j < np; j++)
			sum += coef[j] * afunc[j];
		chisq += (tmp = (Y[i] - sum) / sig[i], tmp * tmp);
	}

	//cout << "Made it here E\n";

	/////////////////////////////////////////
	// Obtain covariance matrix of estimates

	// OLS variance estimator = s^2 * ( X'X )^-1
	// where s^2 = (1/(N-k)) \sum_i=1^N e_i^2

	// Robust cluster variance estimator
	// V_cluster = (X'X)^-1 * \sum_{j=1}^{n_C} u_{j}' * u_j * (X'X)^-1
	// where u_j = \sum_j cluster e_i * x_i

	// Above, e_i is the residual for the ith observation and x_i is a
	// row vector of predictors including the constant.

	// For simplicity, I omitted the multipliers (which are close to 1)
	// from the formulas for Vrob and Vclusters.

	// The formula for the clustered estimator is simply that of the
	// robust (unclustered) estimator with the individual ei*s replaced
	// by their sums over each cluster. xi

	// Williams, R. L. 2000.  A note on robust variance estimation for
	// cluster-correlated data. Biometrics 56: 64


	//  t ( y - yhat X  ) %*%  ( y - yhat)  / nind - np

	// = variance of residuals

	// Variance of residuals

	// j <- ( t( y- m %*% t(b) ) %*% ( y - m %*% t(b) ) ) / ( N - p )
	// print( sqrt(kronecker( solve( t(m) %*% m ) , j )  ))


	// Calcuate S = (XtX)^-1


	vector<vector<double> > Xt;
	Helpers::sizeMatrix(Xt, np, nind);
	for (int i = 0; i < nind; i++)
		for (int j = 0; j < np; j++)
			Xt[j][i] = X[i][j];

	vector<vector<double> > S0;
	Helpers::multMatrix(Xt, X, S0);

	/// if (par::verbose)
	///    {
	//      cout << "beta...\n";
	//      display(coef);
	//      cout << "Sigma(S0a)\n";
	//      display(S0);
	//      cout << "\n";
	///    }

	//cout << "Made it here F\n";

	S0 = Helpers::svd_inverse(S0);

	///  if (par::verbose)
	///    {
//	      cout << "beta...\n";
//	      display(coef);
//	      cout << "Sigma(S0b)\n";
//	      display(S0);
//	      cout << "\n";
	///    }


	////////////////////
	// Standard OLS s^2

	if (!cluster) {

		double sigma = 0.0;
		for (int i = 0; i < nind; i++) {
			double partial = 0.0;
			for (int j = 0; j < np; j++)
				partial += coef[j] * X[i][j];
			partial -= Y[i];
			sigma += partial * partial;
		}
		sigma /= nind - np;

		for (int i = 0; i < np; i++)
			for (int j = 0; j < np; j++)
				S[i][j] = S0[i][j] * sigma;
	}

	///////////////////////////
	// Robust-cluster variance

	if (cluster) {

		vector<vector<double> > sc(nc);
		for (int i = 0; i < nc; i++)
			sc[i].resize(np, 0);

		for (int i = 0; i < nind; i++) {
			double partial = 0.0;
			for (int j = 0; j < np; j++)
				partial += coef[j] * X[i][j];
			partial -= Y[i];

			for (int j = 0; j < np; j++)
				sc[clst[i]][j] += partial * X[i][j];
		}

		vector<vector<double> > meat;
		Helpers::sizeMatrix(meat, np, np);
		for (int k = 0; k < nc; k++) {

			for (int i = 0; i < np; i++)
				for (int j = 0; j < np; j++)
					meat[i][j] += sc[k][i] * sc[k][j];

		}

		vector<vector<double> > tmp1;
		Helpers::multMatrix(S0, meat, tmp1);
		Helpers::multMatrix(tmp1, S0, S);

	}

	///  if (par::verbose)
	///    {
//	      cout << "beta...\n";
//	      display(coef);
//	     cout << "Sigma\n";
//	      display(S);
//	      cout << "\n";
	///    }

}

void LinearRegression::fitUnivariateLM() {

	///  if (par::verbose)
	///    {
	//      cout << "LM VIEW\n";
	//     display(Y);
	//     display(X);
	//     cout << "---\n";
	///    }

	// Speed-up version for univariate case Has set set coef and S

	if (np != 2 || nind == 0)
		return;

	coef.resize(2);
	Helpers::sizeMatrix(S, 2, 2);

	double x_mean = 0, x_var = 0;
	double y_mean = 0, y_var = 0;
	double y_x_covar = 0;

	/////////////////////////////
	// Iterate over individuals

	// X and Y

	for (int i = 0; i < nind; i++) {
		y_mean += Y[i];
		x_mean += X[i][1];
	}

	x_mean /= (double) nind;
	y_mean /= (double) nind;

	for (int i = 0; i < nind; i++) {
		double ty = Y[i] - y_mean;
		double tx = X[i][1] - x_mean;
		y_var += ty * ty;
		x_var += tx * tx;
		y_x_covar += tx * ty;
	}

	y_var /= (double) nind - 1;
	x_var /= (double) nind - 1;
	y_x_covar /= (double) nind - 1;

	// Do not set intercept; only the univariate coefficient
	coef[1] = y_x_covar / x_var;
	S[1][1] = (y_var / x_var - (y_x_covar * y_x_covar) / (x_var * x_var))
			/ (nind - 2);

}

double LinearRegression::getFStat(){
	double r2 = calculateRSquared();
	double result = (r2 * (nind - (np - 1) - 1)) / ((1 - r2) * (np - 1));
	return result;
}

double LinearRegression::findF(){
	double result = 0;

	//TODO: TESTING 01-02-2011
	if (! nind == 0)
	{
		double r2 = calculateRSquared();
		int tempnp = np - 1;
		result = (r2 * (nind - tempnp - 1)) / ((1 - r2) * (tempnp));


		fisher_f dist(tempnp, (nind - tempnp - 1));

		result = cdf(complement(dist, result));
	}
	return result;
}

vector_t LinearRegression::getCoefs() {
	return coef;
}

vector<double> LinearRegression::getVar() {

	double multiplier = 1;

	if (cluster)
		multiplier = (double) (nc) / ((double) (nc - np));

	vector_t var(np);
	for (int i = 0; i < np; i++)
		var[i] = multiplier * S[i][i];

	return var;
}

void LinearRegression::reset() {
	np = 0;
	nind = 0;
	additive.clear();
	dominance.clear();
	haplotype.clear();
	covariate.clear();
	interaction.clear();
	xchr.clear();

	haploid.clear();
	valid.clear();
	coef.clear();
	w.clear();
	S.clear();
	Y.clear();
	X.clear();
	miss.clear();
	label.clear();

	clst.clear();
	C.clear();

	se.clear();
	type.clear();

	sig.clear();
	u.clear();
	v.clear();

	ZS.clear();
	pvalues.clear();
	haploid.resize(0);
	xchr.resize(0);
	order.clear();
	sex_effect = false;
	all_valid = true;
	has_snps = true;
	testParameter = 1; // Permutation test parameter
	RSS = -1;

	// Automatically add intercept now
	label.push_back("M"); // Intercept
	type.push_back(INTERCEPT);
	order.push_back(0);
}

//Is this calculate method only used by the epistasis module?
void LinearRegression::calculate(vector<Marker*> model)
{
	if(model.size() != 2)
	{
		throw MethodException("Linear Regression interaction model requires 2 loci.\n");
	}
	setMissing();

	string mainEffect = "EPI";
	bool genotypic = false;

	// Main effect of SNP 1

	addAdditiveSNP(model[0]);
	label.push_back("ADD1");

	// Main effect of SNP 2

	addAdditiveSNP(model[1]);
	label.push_back("ADD2");

	// Epistasis

	addInteraction(1, 2);
	label.push_back("EPI");


	//COVARIATE ADDITIONS
	if (data_set->num_covariates() > 0)
	{
		for (int i = 0; i < data_set->num_covariates(); i++)
		{
			//cout << "Adding Covariate to Linear Regression: " << getString<int>(data_set->num_covariates()) << " \n";
			addCovariate(i);
			label.push_back(data_set->get_covariate_name(i));
		}
	}
	if (options.getLinRInteraction() && !options.getLinRNoMainSnp())
	{
		int cindex = 2;
		if (genotypic)
		{
			cindex = 3;
		}

		//covariates
		for (int c = 0; c < data_set->num_covariates(); c++)
		{
			//cout << "Adding Interaction to Linear Regression" << "\n";
			addInteraction(1, cindex);
			label.push_back(mainEffect + "x" + data_set->get_covariate_name(c));
			if (genotypic)
			{
				//cout << "Adding genotypic Interaction to Linear Regression" << "\n";
				addInteraction(2, cindex);
				//				if(par::twoDFmodel_hethom){
				//					label.push_back("HETx" + clistname[c]);
				//				}
				//				else{
				label.push_back("DOMDEVx" + data_set->get_covariate_name(c));
				//				}
			}
			cindex++;
		}//end for loop for adding interactions to linear model....
	}
	//END COVARIATE ADDITIONS

	// Build design matrix

	//cout << "Builing Design Matrix " << "\n";
	buildDesignMatrix();

	// Prune out any remaining missing individuals
	// No longer needed

	//         pruneY();


	// Fit linear model
	//cout << "Fitting Linear Model" << "\n";
	fitLM();

	// Did model fit okay?

	//cout << "Valid Parameters? " << "\n";
	validParameters();

	// Obtain estimates and statistic
	//TODO:  should the following code be in here???
	vector<double> var;
	if (all_valid)
	{
		var = getVar();
	}
	else
	{
		var.clear();
		var.resize(np, 0);
	}

	ZS.push_back(0);
	pvalues.push_back(0);
	for (int p = 1; p < np; p++)
	{
		bool okay = var[p] < 1e-20 || !Helpers::realnum(var[p]) ? false : all_valid;
		double se = 0;
		double Z = 0;
		double pvalue = 1;
		if (okay)
		{
			se = sqrt(var[p]);
			Z = coef[p] / se;
			ZS.push_back(Z);
			pvalue = Helpers::pT(Z, Y.size() - np);///pT(Z,Y.size()-np);
			pvalues.push_back(pvalue);
		}
	}
	//end TODO...
	testParameter = 3; // interaction
}

void LinearRegression::calculate(Marker* l) {
	reset();
	Marker* mark = l;//data_set->get_locus(l);

	bool X = false;
	bool automaticSex = false;
	bool variationInSex = data_set->num_males() > 0 && data_set->num_females()
			> 0;

	if (!options.getLinRNoMainSnp())//par::assoc_glm_without_main_smp)
	{
		if (options.get_xchr_model() == 0) {
			//if(par::chr_sex[locus[l]->chr] || par::chr_haploid[locus[l]->chr])
			//	continue;
			if (opts::_CHRX_ == mark->getChrom())
				return;
		} else {
			if (mark->getChrom() == opts::_CHRX_) {
				X = true;
				//cout << "X is true for " << mark->toString() << endl;
			}
		}
	}

	this->setMissing();

	if (options.getLinRModelType() == "DOMINANT") {
		setDominant();
	} else if (options.getLinRModelType() == "RECESSIVE") {
		setRecessive();
	}
	//	else if(options.getLinRModelType() == "ADDITIVE"){

	//	}

	string mainEffect = "";
	bool genotypic = false;

	if (!options.getLinRNoMainSnp()) {//par::assoc_glm_without_main_snp){
		//genotypic = par::chr_haploid[locus[l]->chr] ? false : par::twoDFmodel;

		if (options.getLinRModelType() == "RECESSIVE") {
			mainEffect = "REC";
		} else if (options.getLinRModelType() == "DOMINANT") {
			mainEffect = "DOM";
		}
		//else if(par::twoDFmodel_hethom){
		//	mainEffect = "HOM";
		//}
		else {
			mainEffect = "ADD";
		}
		addAdditiveSNP(l);
		label.push_back(mainEffect);

		if (genotypic) {
			addDominanceSNP(l);

			//			if(par::twoDFmodel_hethom)
			//				label.push_back("HET");
			//			else
			//				label.push_back("DOMDEV");
		}
	}

	//	if(par::chap_test){
	//		for(int h=1; h < whap->current->group.size(); h++){
	//			lm->addHaplotypeDosage(whap->current->group[h]);
	//			lm->label.push_back("WHAP" + int2str(h+1));
	//		}
	//	}

	//	if(par::proxy_glm){
	//		set<int> t1 = haplo->makeSetFromMap(haplo->testSet);
	//		ml->addHaplotypeDosage(t1);
	//		lm->label.push_back("PROXY");
	//
	//	}

	if (options.getLinRCondition()) {
		//		if(par::chap_test){
		//			for(int c = 0; c < conditioner.size(); c++){
		//				if(whap->current->masked_conditioning_snps[c]){
		//					lm->addAdditiveSNP(conditioner[c]);
		//					lm->label.push_back(locus[conditioner[c]]->name);
		//				}
		//			}
		//		}
		//		else{
		vector<int> list = options.getLinRConditionList();
		for (unsigned int i = 0; i < list.size(); i++) {
			addAdditiveSNP(list[i]);
			label.push_back(data_set->get_locus(list[i])->getRSID());
		}
		//		}
	}

	//	if(X){
	//		cout << "I HAVE AN X!!!" << endl;
	//	}
	//	if(!options.getAutoSexEffect()){
	//		cout << "I DONT HAVE AUTO SEX EFFECT!!!" << endl;
	//	}

	if ((options.getSexEffect() || (X && !options.getAutoSexEffect()))
			&& variationInSex) {
		automaticSex = true;
		addSexEffect();
		//cout << "ADDING SEX EFFECT!!!!!\n";
		//cout << mA << "\t" << mB << endl;
		label.push_back("SEX");
	}
	if (data_set->num_covariates() > 0) {
		for (int i = 0; i < data_set->num_covariates(); i++) {
			addCovariate(i);
			label.push_back(data_set->get_covariate_name(i));
		}
	}

	//interactions
	if (options.getLinRInteraction() && !options.getLinRNoMainSnp()) {
		int cindex = 2;
		if (genotypic) {
			cindex = 3;
		}
		vector<int> list = options.getLinRConditionList();
		for (int c = 0; c < (int) list.size(); c++) {
			addInteraction(1, cindex);
			label.push_back(mainEffect + "xCSNP" + getString<int> (c + 1));

			if (genotypic) {
				addInteraction(2, cindex);
				//				if(par::twoDFmodel_hethom){
				//					label.push_back("HETxCSNP"+getString<int>(c+1));
				//				}
				//				else{
				label.push_back("DOMDEVxCSNP" + getString<int> (c + 1));
				//				}
			}
			cindex++;
		}

		if (automaticSex) {
			addInteraction(1, cindex);
			label.push_back(mainEffect + "xSEX");
			if (genotypic) {
				addInteraction(2, cindex);
				//				if(par::twoDFmodel_hethom)
				//					label.push_back("HETxSEX");
				//				else
				label.push_back("DOMDEVxSEX");
			}
			cindex++;
		}
		//covariates
		for (int c = 0; c < data_set->num_covariates(); c++) {
			addInteraction(1, cindex);
			label.push_back(mainEffect + "x" + data_set->get_covariate_name(c));
			if (genotypic) {
				addInteraction(2, cindex);
				//				if(par::twoDFmodel_hethom){
				//					label.push_back("HETx" + clistname[c]);
				//				}
				//				else{
				label.push_back("DOMDEVx" + data_set->get_covariate_name(c));
				//				}
			}
			cindex++;
		}
	}
	//fancy X chrom models
	if (X && automaticSex && options.get_xchr_model() > 2) {
		//interaction between allelic term and sex
		int sindex = 2;
		if (genotypic) {
			sindex++;
		}
		sindex += options.getLinRConditionList().size();
		addInteraction(2, sindex);
		label.push_back("XxSEX");

		//xchr model3 : test ADD + XxSEX
		//xchr model4 : test ADD + DOM + XxSEX
	}

	buildDesignMatrix();
	fitLM();
	validParameters();
	vector<double> var;
	if (all_valid) {
		var = getVar();
	} else {
		var.clear();
		var.resize(np, 0);
	}

	ZS.push_back(0);
	pvalues.push_back(0);
	for (int p = 1; p < np; p++)
	{
		bool okay = var[p] < 1e-20 || !Helpers::realnum(var[p]) ? false : all_valid;
		double se = 0;
		double Z = 0;
		double pvalue = 1;
		if (okay)
		{
			se = sqrt(var[p]);
			Z = coef[p] / se;
			ZS.push_back(Z);
			pvalue = Helpers::pT(Z, Y.size() - np);///pT(Z,Y.size()-np);
			pvalues.push_back(pvalue);
		}
	}
}

void LinearRegression::calculate(vector<int> m)
{
	  vector<Marker*> snps;
	  for(unsigned int i = 0; i < m.size(); i++)
	  {
		  snps.push_back(data_set->get_locus(m[i]));
	  }

	  calculate(snps);
}

void LinearRegression::displayResults(ofstream & OUT, Marker * loc) {

	vector<double> var;

	if (all_valid)
		var = getVar();
	else {
		var.clear();
		var.resize(np, 0);
	}

	for (int p = 1; p < np; p++) // skip intercept
	{

		bool okay = var[p] < 1e-20 || !Helpers::realnum(var[p]) ? false : all_valid;

		double se = 0;
		double Z = 0;
		double pvalue = 1;

		if (okay) {
			se = sqrt(var[p]);
			Z = coef[p] / se;
			pvalue = Helpers::chiprobP(Z, Y.size() - np);///pT(Z,Y.size()-np);
		}

		// If filtering p-values
		///      if ( (!par::pfilter) || pvalue <= par::pfvalue )
		///	{

		OUT << loc->getChrom() << " " << loc->getRSID() << " "
				<< loc->getBPLOC() << " " << loc->getAllele1() << " "
				<< label[p] << " " << Y.size() << " ";

		///	  if (okay)
		///	    {
		OUT << coef[p] << " ";

		///	      if (par::display_ci)

		double zt = Helpers::ltqnorm(1 - (1 - options.getCI()) / 2);

		OUT << se << " " << coef[p] - zt * se << " " << coef[p] + zt * se
				<< " ";

		OUT << Z << " " << pvalue;
		///	    }
		///	  else
		///	    {
		///	      OUT << setw(10) << "NA" << " ";

		///	      if (par::display_ci)
		///		OUT << setw(8) << "NA" << " "
		///		    << setw(8) << "NA" << " "
		///		    << setw(8) << "NA" << " ";

		///	      OUT << setw(12) << "NA" << " "
		///		  << setw(12) << "NA";
		///	    }

		OUT << "\n";
		///	}

	}

}

double LinearRegression::calculateRSS()
{

	// Might already be calculated?

	if (RSS >= 0)
	{
		return RSS;
	}

	// Calculate residual sum of squares (RSS)

	RSS = 0;

	for (int i = 0; i < nind; i++)
	{
		double t = Y[i];

		for (int p = 0; p < np; p++)
		{
			t -= coef[p] * X[i][p];
		}
		t *= t;
		RSS += t;
	}

	return RSS;
}

double LinearRegression::calculateRSquared()
{
	// Return coefficient of determination. Ifnot already calculated,
	// get residual sum of squares first (set to -1)

	if (RSS < 0)
	{
		RSS = calculateRSS();
	}

	double SSy = varY * (nind - 1);

	double r = (SSy - RSS) / SSy;

	return r > 0 ? (r > 1 ? 1 : r) : 0;

}

double LinearRegression::calculateAdjustedRSquared() {

	double ra = 1 - ((double) (nind - 1) / (double) (nind - np - 1)) * (1
			- calculateRSquared());

	return ra > 0 ? (ra > 1 ? 1 : ra) : 0;

}

double LinearRegression::calculateMallowC(LinearRegression * submodel) {

	// Mallow's C = RSSm / S^2 + 2(m+1)-n
	// where S^2 = RSSk / (n-k-1);

	double Sk = calculateRSS() / (nind - np - 1);
	return (submodel->calculateRSS() / Sk) + 2 * (submodel->np + 1) - nind;
}

double LinearRegression::calculateFTest(LinearRegression * submodel) {

	double RSSk = calculateRSS();
	double RSSm = submodel->calculateRSS();

	return ((RSSm - RSSk) / (double) (np - submodel->np)) / (RSSk
			/ (double) (nind - np - 1));
}

double LinearRegression::getZ() {
	vector<double> var = getVar();
	bool okay =
			var[testParameter] < 1e-20 || !Helpers::realnum(var[testParameter]) ? false
					: all_valid;

	if (okay) {
		double se = sqrt(var[testParameter]);
		double Z = coef[testParameter] / se;
		return Z;
	} else
		return 0;
}

double LinearRegression::getTestCoef() {
	return coef[testParameter];
}

int LinearRegression::getCalcMissing() {
	return (int) Y.size();
}

double LinearRegression::getPValue() {
	vector<double> var = getVar();
	bool okay =
			var[testParameter] < 1e-20 || !Helpers::realnum(var[testParameter]) ? false
					: all_valid;

	if (okay) {
		double se = sqrt(var[testParameter]);
		double Z = coef[testParameter] / se;
		//cout << "coef: " << coef[testParameter] << endl;
		//cout << "se: " << se << endl;
		//cout << "z: " << Z << endl;
		//cout << "testParameter: " << testParameter << endl;
		//cout << "Y: " << Y.size() << endl;
		//cout << "np: " << np << endl;
		return Helpers::pT(Z, Y.size() - np);
	} else
		return 1;
}

vector<double> LinearRegression::getPVals() {
	int tmp = testParameter;
	vector_t res;
	for (testParameter = 1; testParameter < np; testParameter++)
		res.push_back(getPValue());
	testParameter = tmp;
	return res;
}

}
