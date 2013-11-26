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

namespace Methods {
void LinearRegression::display(vector<vector<double> > & m) {
	cout << "\n";
	for (int i = 0; i < (int) m.size(); i++) {
		cout << i << ")\t";
		for (int j = 0; j < (int) m.at(i).size(); j++)
			cout << m.at(i).at(j) << " ";
		cout << "\n";
	}
	cout << "\n";
}

void LinearRegression::display(vector<double> & m) {
	cout << "\n";
	for (int i = 0; i < (int) m.size(); i++)
		cout << i << ")\t" << m.at(i) << "\n";
	cout << "\n";
	cout << "\n";
}

void LinearRegression::display(vector<int> & m) {
	cout << "\n";
	for (int i = 0; i < (int) m.size(); i++)
		cout << i << ")\t" << m.at(i) << "\n";
	cout << "\n";
	cout << "\n";
}

LinearRegression::LinearRegression()
{
	nc = 0;
	cluster = false;
	RSS = -1;
}

void LinearRegression::setDependent() {
	// Set dependent variable and intercept
	Y.clear();

	for (int i = 0; i < data_set->num_inds(); i++)
		if (!miss.at(i)) {
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
				Y.push_back(data_set->get_sample(i)->getPheno(index));
			}
			else
			{
				Y.push_back(data_set->get_sample(i)->getPheno());
			}
		}
}

void LinearRegression::pruneY() {

	//////////////////////////////////
	// Prune out rows that are missing

	vector<double> Y2;
	for (int i = 0; i < (int) Y.size(); i++)
		if (!miss.at(i))
			Y2.push_back(Y.at(i));
	Y = Y2;
}

void covsrt(vector<vector<double> > & covar, vector<bool> &ia, const int mfit) {
	int i, j, k;

	int ma = ia.size();
	for (i = mfit; i < ma; i++)
		for (j = 0; j < i + 1; j++)
			covar.at(i).at(j) = covar.at(j).at(i) = 0.0;
	k = mfit - 1;
	for (j = ma - 1; j >= 0; j--) {
		if (ia.at(j)) {
			for (i = 0; i < ma; i++)
				SWAP(covar.at(i).at(k), covar.at(i).at(j));
			for (i = 0; i < ma; i++)
				SWAP(covar.at(k).at(i), covar.at(j).at(i));
			k--;
		}
	}
}

void gaussj(vector<vector<double> > & a, vector<vector<double> > & b) {
	int i, icol = 0, irow = 0, j, k, l, ll;
	double big, dum, pivinv;

	int n = a.size();
	int m = b.at(0).size();
	vector<double> indxc(n), indxr(n), ipiv(n);
	for (j = 0; j < n; j++)
		ipiv.at(j) = 0;
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++)
			if (ipiv.at(j) != 1)
				for (k = 0; k < n; k++) {
					if (ipiv.at(k) == 0) {
						if (fabs(a.at(j).at(k)) >= big) {
							big = fabs(a.at(j).at(k));
							irow = j;
							icol = k;
						}
					}
				}
		++(ipiv.at(icol));
		if (irow != icol) {
			for (l = 0; l < n; l++)
				SWAP(a.at(irow).at(l), a.at(icol).at(l));
			for (l = 0; l < m; l++)
				SWAP(b.at(irow).at(l), b.at(icol).at(l));
		}
		indxr.at(i) = irow;
		indxc.at(i) = icol;
		if (a.at(icol).at(icol) == 0.0)
			throw MethodException("gaussj: Singular Matrix");
		pivinv = 1.0 / a.at(icol).at(icol);
		a.at(icol).at(icol) = 1.0;
		for (l = 0; l < n; l++)
			a.at(icol).at(l) *= pivinv;
		for (l = 0; l < m; l++)
			b.at(icol).at(l) *= pivinv;
		for (ll = 0; ll < n; ll++)
			if (ll != icol) {
				dum = a.at(ll).at(icol);
				a.at(ll).at(icol) = 0.0;
				for (l = 0; l < n; l++)
					a.at(ll).at(l) -= a.at(icol).at(l) * dum;
				for (l = 0; l < m; l++)
					b.at(ll).at(l) -= b.at(icol).at(l) * dum;
			}
	}
	for (l = n - 1; l >= 0; l--) {
		if (indxr.at(l) != indxc.at(l))
			for (k = 0; k < n; k++)
				SWAP(a.at(k).at((int) indxr.at(l)), a.at(k).at((int) indxc.at(l)));
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
		if (ia.at(j))
			mfit++;
	if (mfit == 0)
		throw MethodException("lfit: no parameters to be fitted");
	for (j = 0; j < mfit; j++) {
		for (k = 0; k < mfit; k++)
			covar.at(j).at(k) = 0.0;
		beta.at(j).at(0) = 0.0;
	}
	for (i = 0; i < ndat; i++) {
		afunc = X.at(i);

		ym = y.at(i);
		if (mfit < ma) {
			for (j = 0; j < ma; j++)
				if (!ia.at(j))
					ym -= a.at(j) * afunc.at(j);
		}
		sig2i = 1.0 / SQR(sig.at(i));
		for (j = 0, l = 0; l < ma; l++) {
			if (ia.at(l)) {
				wt = afunc.at(l) * sig2i;
				for (k = 0, m = 0; m <= l; m++)
					if (ia.at(m))
						covar.at(j).at(k++) += wt * afunc.at(m);
				beta.at(j++).at(0) += ym * wt;
			}
		}
	}
	for (j = 1; j < mfit; j++)
		for (k = 0; k < j; k++)
			covar.at(k).at(j) = covar.at(j).at(k);
	vector<vector<double> > temp;
	Helpers::sizeMatrix(temp, mfit, mfit);
	for (j = 0; j < mfit; j++)
		for (k = 0; k < mfit; k++)
			temp.at(j).at(k) = covar.at(j).at(k);
	gaussj(temp, beta);
	for (j = 0; j < mfit; j++)
		for (k = 0; k < mfit; k++)
			covar.at(j).at(k) = temp.at(j).at(k);
	for (j = 0, l = 0; l < ma; l++)
		if (ia.at(l))
			a.at(l) = beta.at(j++).at(0);
	chisq = 0.0;
	for (i = 0; i < ndat; i++) {
		afunc = X.at(i);
		sum = 0.0;
		for (j = 0; j < ma; j++)
			sum += a.at(j) * afunc.at(j);
		chisq += SQR((y.at(i) - sum) / sig.at(i));
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
		if (novelc.find(clst.at(i)) == novelc.end())
			novelc.insert(make_pair(clst.at(i), k++));
	nc = novelc.size();
	if (nc < 2)
		throw MethodException("Must have at least two clusters");
	for (int i = 0; i < nind; i++)
		clst.at(i) = novelc.find(cl.at(i))->second;

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
		mean += Y.at(i);
	}
	if (actualN == 0) {
		varY = 0;
		return;
	}
	mean /= (double) actualN;

	for (int i = 0; i < nind; i++) {
		varY += (Y.at(i) - mean) * (Y.at(i) - mean);
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
			cout << "VO " << i << "\t" << Y.at(i) << "\t";
			for (int j = 0; j < np; j++)
				cout << X.at(i).at(j) << "\t";
			cout << "\n";
		}
	}
	coef.resize(np);
	Helpers::sizeMatrix(S, np, np);
	if (np == 0 || nind == 0 || !all_valid) {
		return;
	}

	setVariance();
	sig.resize(nind, sqrt(1.0 / sqrt((double) nind)));

	w.resize(np);
	Helpers::sizeMatrix(u, nind, np);
	Helpers::sizeMatrix(v, np, np);
	int i, j;
	const double TOL = 1.0e-13;
	double wmax, tmp, thresh, sum;

	vector<double> b(nind), afunc(np);
	for (i = 0; i < nind; i++) {
		afunc = X.at(i);
		tmp = 1.0 / sig.at(i);
		for (j = 0; j < np; j++)
			u.at(i).at(j) = afunc.at(j) * tmp;
		b.at(i) = Y.at(i) * tmp;
	}

	Helpers::svdcmp(u, w, v);

	wmax = 0.0;
	for (j = 0; j < np; j++)
		if (w.at(j) > wmax)
			wmax = w.at(j);
	thresh = TOL * wmax;
	for (j = 0; j < np; j++)
		if (w.at(j) < thresh)
			w.at(j) = 0.0;

	Helpers::svbksb(u, w, v, b, coef);

	chisq = 0.0;
	for (i = 0; i < nind; i++) {
		afunc = X.at(i);
		sum = 0.0;
		for (j = 0; j < np; j++)
			sum += coef.at(j) * afunc.at(j);
		chisq += (tmp = (Y.at(i) - sum) / sig.at(i), tmp * tmp);
	}

	vector<vector<double> > Xt;
	Helpers::sizeMatrix(Xt, np, nind);
	for (int i = 0; i < nind; i++)
		for (int j = 0; j < np; j++)
			Xt.at(j).at(i) = X.at(i).at(j);

	vector<vector<double> > S0;
	Helpers::multMatrix(Xt, X, S0);

	S0 = Helpers::svd_inverse(S0);


	////////////////////
	// Standard OLS s^2

	if (!cluster) {

		double sigma = 0.0;
		for (int i = 0; i < nind; i++) {
			double partial = 0.0;
			for (int j = 0; j < np; j++)
				partial += coef.at(j) * X.at(i).at(j);
			partial -= Y.at(i);
			sigma += partial * partial;
		}
		sigma /= nind - np;

		for (int i = 0; i < np; i++)
			for (int j = 0; j < np; j++)
				S.at(i).at(j) = S0.at(i).at(j) * sigma;
	}

	///////////////////////////
	// Robust-cluster variance

	if (cluster) {

		vector<vector<double> > sc(nc);
		for (int i = 0; i < nc; i++)
			sc.at(i).resize(np, 0);

		for (int i = 0; i < nind; i++) {
			double partial = 0.0;
			for (int j = 0; j < np; j++)
				partial += coef.at(j) * X.at(i).at(j);
			partial -= Y.at(i);

			for (int j = 0; j < np; j++)
				sc.at(clst.at(i)).at(j) += partial * X.at(i).at(j);
		}

		vector<vector<double> > meat;
		Helpers::sizeMatrix(meat, np, np);
		for (int k = 0; k < nc; k++) {

			for (int i = 0; i < np; i++)
				for (int j = 0; j < np; j++)
					meat.at(i).at(j) += sc.at(k).at(i) * sc.at(k).at(j);

		}

		vector<vector<double> > tmp1;
		Helpers::multMatrix(S0, meat, tmp1);
		Helpers::multMatrix(tmp1, S0, S);

	}
}

void LinearRegression::fitUnivariateLM() {

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
		y_mean += Y.at(i);
		x_mean += X.at(i).at(1);
	}

	x_mean /= (double) nind;
	y_mean /= (double) nind;

	for (int i = 0; i < nind; i++) {
		double ty = Y.at(i) - y_mean;
		double tx = X.at(i).at(1) - x_mean;
		y_var += ty * ty;
		x_var += tx * tx;
		y_x_covar += tx * ty;
	}

	y_var /= (double) nind - 1;
	x_var /= (double) nind - 1;
	y_x_covar /= (double) nind - 1;

	// Do not set intercept; only the univariate coefficient
	coef.at(1) = y_x_covar / x_var;
	S.at(1).at(1) = (y_var / x_var - (y_x_covar * y_x_covar) / (x_var * x_var))
			/ (nind - 2);

}

double LinearRegression::getFStat(){
	double r2 = calculateRSquared();
	double result = (r2 * (nind - (np - 1) - 1)) / ((1 - r2) * (np - 1));
	return result;
}

double LinearRegression::findF(){
	double result = 0;

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
		var.at(i) = multiplier * S.at(i).at(i);

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

	addAdditiveSNP(model.at(0));
	label.push_back("ADD1");

	// Main effect of SNP 2

	addAdditiveSNP(model.at(1));
	label.push_back("ADD2");

	// Epistasis

	addInteraction(1, 2);
	label.push_back("EPI");


	//COVARIATE ADDITIONS
	if (data_set->num_covariates() > 0)
	{
		vector<string> covsToUse = options.getCovars();
		opts::printLog(getString<int>(covsToUse.size()) + " Covariates being used in Linear Model.\n");

		if(options.doCovars())
		{
			if(options.doCovarsName())
			{
				//only add covariates if name matches
				for(int i = 0; i < (int)covsToUse.size(); i++)
				{
					addCovariate(data_set->get_covariate_index(covsToUse.at(i)));
					label.push_back(covsToUse.at(i));
				}
			}
			else if(options.doCovarsNumber())
			{
				//only add covariates if number matches (column number)
				for(int i = 0; i < (int)covsToUse.size(); i++)
				{
					addCovariate(atoi(covsToUse.at(i).c_str()));
					label.push_back(data_set->get_covariate_name(atoi(covsToUse.at(i).c_str())));
				}
			}
		}
		else
		{
			for(int i=0; i<data_set->num_covariates(); i++)
			{
				//user did not specify which covariates to use, use all covariates found in covariate input file
				addCovariate(i);
				label.push_back(data_set->get_covariate_name(i));
			}
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
		if(options.doCovars())
		{
			vector<string> covsToUse = options.getCovars();

			for(int c = 0; c < (int)covsToUse.size(); c++)
			{
				if(options.doCovarsName())
				{
					addInteraction(1, cindex);
					label.push_back(mainEffect + "x" + covsToUse.at(c));
					if(genotypic)
					{
						addInteraction(2, cindex);
						label.push_back("DOMDEVx" + covsToUse.at(c));
					}
				}
				else if (options.doCovarsNumber())
				{
					addInteraction(1, cindex);
					label.push_back(mainEffect + "x" + data_set->get_covariate_name(atoi(covsToUse.at(c).c_str())));
					if(genotypic)
					{
						addInteraction(2, cindex);
						label.push_back("DOMDEVx" + data_set->get_covariate_name(atoi(covsToUse.at(c).c_str())));
					}
				}
				cindex++;
			}
		}
		else
		{ //user did not specify which covariates to use, so use all covariates contained in -covar-file
			for (int c = 0; c < data_set->num_covariates(); c++)
			{
				addInteraction(1, cindex);
				label.push_back(mainEffect + "x" + data_set->get_covariate_name(c));
				if (genotypic)
				{
					addInteraction(2, cindex);
					label.push_back("DOMDEVx" + data_set->get_covariate_name(c));
				}
				cindex++;
			}//end for loop for adding interactions to linear model....
		}
	}
	//END COVARIATE ADDITIONS

	// Build design matrix

	buildDesignMatrix();

	// Fit linear model
	fitLM();

	// Did model fit okay?
	validParameters();

	// Obtain estimates and statistic
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
		bool okay = var.at(p) < 1e-20 || !Helpers::realnum(var.at(p)) ? false : all_valid;
		double se = 0;
		double Z = 0;
		double pvalue = 1;
		if (okay)
		{
			se = sqrt(var.at(p));
			Z = coef.at(p) / se;
			ZS.push_back(Z);
			pvalue = Helpers::pT(Z, Y.size() - np);
			pvalues.push_back(pvalue);
		}
	}
	testParameter = 3; // interaction
}

void LinearRegression::calculate(Marker* l)
{
	reset();
	Marker* mark = l;

	bool X = false;
	bool automaticSex = false;
	bool variationInSex = data_set->num_males() > 0 && data_set->num_females() > 0;

	if (!options.getLinRNoMainSnp())
	{
		if (options.get_xchr_model() == 0)
		{
			if (opts::_CHRX_ == mark->getChrom())
				return;
		} else {
			if (mark->getChrom() == opts::_CHRX_)
			{
				X = true;
			}
		}
	}

	this->setMissing();

	if (options.getLinRModelType() == "DOMINANT")
	{
		setDominant();
	}
	else if (options.getLinRModelType() == "RECESSIVE")
	{
		setRecessive();
	}


	string mainEffect = "";
	bool genotypic = false;

	if (!options.getLinRNoMainSnp()) {

		if (options.getLinRModelType() == "RECESSIVE")
		{
			mainEffect = "REC";
		}
		else if (options.getLinRModelType() == "DOMINANT")
		{
			mainEffect = "DOM";
		}
		else {
			mainEffect = "ADD";
		}
		addAdditiveSNP(l);
		label.push_back(mainEffect);

		if (genotypic)
		{
			addDominanceSNP(l);
		}
	}

	if (options.getLinRCondition())
	{
		vector<int> list = options.getLinRConditionList();
		for (unsigned int i = 0; i < list.size(); i++)
		{
			addAdditiveSNP(list.at(i));
			label.push_back(data_set->get_locus(list.at(i))->getRSID());
		}
	}

	if ((options.getSexEffect() || (X && !options.getAutoSexEffect())) && variationInSex)
	{
		automaticSex = true;
		addSexEffect();
		label.push_back("SEX");
	}

	if (data_set->num_covariates() > 0)
	{
		vector<string> covsToUse = options.getCovars();

		if(options.doCovars())
		{
			if(options.doCovarsName())
			{
				//only add covariates if name matches
				for(int i = 0; i < (int)covsToUse.size(); i++)
				{
					addCovariate(data_set->get_covariate_index(covsToUse.at(i)));
					label.push_back(covsToUse.at(i));
				}
			}
			else if(options.doCovarsNumber())
			{
				//only add covariates if number matches (column number)
				for(int i = 0; i < (int)covsToUse.size(); i++)
				{
					addCovariate(atoi(covsToUse.at(i).c_str()));
					label.push_back(data_set->get_covariate_name(atoi(covsToUse.at(i).c_str())));
				}
			}
		}
		else
		{
			for(int i=0; i<data_set->num_covariates(); i++)
			{
				//user did not specify which covariates to use, use all covariates found in covariate input file
				addCovariate(i);
				label.push_back(data_set->get_covariate_name(i));
			}
		}

	}

	//interactions
	if (options.getLinRInteraction() && !options.getLinRNoMainSnp())
	{
		int cindex = 2;
		if (genotypic)
		{
			cindex = 3;
		}
		vector<int> list = options.getLinRConditionList();
		for (int c = 0; c < (int) list.size(); c++)
		{
			addInteraction(1, cindex);
			label.push_back(mainEffect + "xCSNP" + getString<int> (c + 1));

			if (genotypic)
			{
				addInteraction(2, cindex);
				label.push_back("DOMDEVxCSNP" + getString<int> (c + 1));
			}
			cindex++;
		}

		if (automaticSex)
		{
			addInteraction(1, cindex);
			label.push_back(mainEffect + "xSEX");
			if (genotypic)
			{
				addInteraction(2, cindex);
				label.push_back("DOMDEVxSEX");
			}
			cindex++;
		}


		//covariates
		//Fixed 02-28-2011, allowed for the use of -covars-name, -covars-number for
		//specifying which covariates to use...
		if(options.doCovars())
		{
			vector<string> covsToUse = options.getCovars();

			for(int c = 0; c < (int)covsToUse.size(); c++)
			{
				if(options.doCovarsName())
				{
					addInteraction(1, cindex);
					label.push_back(mainEffect + "x" + covsToUse.at(c));
					if(genotypic)
					{
						addInteraction(2, cindex);
						label.push_back("DOMDEVx" + covsToUse.at(c));
					}
				}
				else if (options.doCovarsNumber())
				{
					addInteraction(1, cindex);
					label.push_back(mainEffect + "x" + data_set->get_covariate_name(atoi(covsToUse.at(c).c_str())));
					if(genotypic)
					{
						addInteraction(2, cindex);
						label.push_back("DOMDEVx" + data_set->get_covariate_name(atoi(covsToUse.at(c).c_str())));
					}
				}
				cindex++;
			}
		}
		else
		{
			//user did not specify which covariates to use, so use all covariates contained in -covar-file
			for (int c = 0; c < data_set->num_covariates(); c++)
			{
				addInteraction(1, cindex);
				label.push_back(mainEffect + "x" + data_set->get_covariate_name(c));
				if (genotypic)
				{
					addInteraction(2, cindex);
					label.push_back("DOMDEVx" + data_set->get_covariate_name(c));
				}
				cindex++;
			}//end for loop for adding interactions to linear model....
		}
	}//end add interactions...

	//fancy X chrom models
	if (X && automaticSex && options.get_xchr_model() > 2)
	{
		//interaction between allelic term and sex
		int sindex = 2;
		if (genotypic)
		{
			sindex++;
		}
		sindex += options.getLinRConditionList().size();
		addInteraction(2, sindex);
		label.push_back("XxSEX");
	}

	buildDesignMatrix();
	fitLM();
	validParameters();
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
		bool okay = var.at(p) < 1e-20 || !Helpers::realnum(var.at(p)) ? false : all_valid;
		double se = 0;
		double Z = 0;
		double pvalue = 1;
		if (okay)
		{
			se = sqrt(var.at(p));
			Z = coef.at(p) / se;
			ZS.push_back(Z);
			pvalue = Helpers::pT(Z, Y.size() - np);
			pvalues.push_back(pvalue);
		}
	}
}

void LinearRegression::calculate(vector<int> m)
{
	  vector<Marker*> snps;
	  for(unsigned int i = 0; i < m.size(); i++)
	  {
		  snps.push_back(data_set->get_locus(m.at(i)));
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

		bool okay = var.at(p) < 1e-20 || !Helpers::realnum(var.at(p)) ? false : all_valid;

		double se = 0;
		double Z = 0;
		double pvalue = 1;

		if (okay) {
			se = sqrt(var.at(p));
			Z = coef.at(p) / se;
			pvalue = Helpers::chiprobP(Z, Y.size() - np);
		}
		OUT << loc->getChrom() << " " << loc->getRSID() << " "
				<< loc->getBPLOC() << " " << loc->getAllele1() << " "
				<< label.at(p) << " " << Y.size() << " ";

		OUT << coef.at(p) << " ";

		double zt = Helpers::ltqnorm(1 - (1 - options.getCI()) / 2);

		OUT << se << " " << coef.at(p) - zt * se << " " << coef.at(p) + zt * se
				<< " ";

		OUT << Z << " " << pvalue;

		OUT << "\n";
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
		double t = Y.at(i);

		for (int p = 0; p < np; p++)
		{
			t -= coef.at(p) * X.at(i).at(p);
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
			var.at(testParameter) < 1e-20 || !Helpers::realnum(var.at(testParameter)) ? false
					: all_valid;

	if (okay) {
		double se = sqrt(var.at(testParameter));
		double Z = coef.at(testParameter) / se;
		return Z;
	} else
		return 0;
}

double LinearRegression::getTestCoef() {
	return coef.at(testParameter);
}

int LinearRegression::getCalcMissing() {
	return (int) Y.size();
}

double LinearRegression::getPValue() {
	vector<double> var = getVar();
	bool okay =
			var.at(testParameter) < 1e-20 || !Helpers::realnum(var.at(testParameter)) ? false
					: all_valid;

	if (okay) {
		double se = sqrt(var.at(testParameter));
		double Z = coef.at(testParameter) / se;
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
