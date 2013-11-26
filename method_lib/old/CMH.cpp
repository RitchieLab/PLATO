///Liberally pulled from Plink and adapted to work in the Plato/wasp evironment


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
#include <fstream>
#include <cmath>
#include "General.h"
#include "Helpers.h"
#include "CMH.h"
namespace Methods{

/*
 * change data sets to use for calculation
 * need to reflag samples to use
 */
void CMH::resetDataSet(DataSet* ds){
	data_set = ds;

	flagSamples();
}

/*
 * marks samples to use for calculation
 */
void CMH::flagSamples(){
	map<Sample*, int> samp_clusters = options.getSampleClusters();
	map<Sample*, int>::iterator siter;
	samp_flags.resize(data_set->num_inds(), false);
	for(int s = 0; s < data_set->num_inds(); s++){
		siter = samp_clusters.find(data_set->get_sample(s));
		if(siter == samp_clusters.end()){
			samp_flags.at(s) = true;
		}
	}

}

/*
 * main calculation method
 * calculate CMH for one marker
 */
void CMH::calculate(Marker* mark){

	if(!options.doClusterFile()){
		throw MethodException("-cluster-file is a required argument for CMH.  Please specify a cluster file.\n");
	}
	if(options.getClusters().size() == 0){
		throw MethodException("CMH -- number of clusters = 0!\n");
	}

	if(options.getCMH2x2xK()){
		calcMantelHaenszel_2x2xK(mark);
	}
	else{
		calcCMH(mark);
	}
}


/*
 * standard CMH calculation
 */
void CMH::calcCMH(Marker* mark) {

	///////////////////////////////////
	// Basic 2 x 2 x K CMH test
	// i.e. Disease x allele x strata
	// is taken care of in assoc.cpp
	// (i.e. allows for permutation, sets, etc)

	//////////////////////////////////
	// Any individual not assigned to a cluster,
	// making missing phenotype

	map<Sample*, int> samp_clusters = options.getSampleClusters();

	///////////////////////////////////
	// Generalized I x J x K CMH test

	// Either ordinal or normal
	// i.e. test strata X SNP controlling for disease

	if (options.doCMHIxJxK() || options.doCMHOrdinal())
	{

		if (options.doCMHOrdinal() && !data_set->get_binary_trait())
			throw MethodException("-cmh-ordinal specified but the phenotype is only binary: use -mh\n");

		if (options.getClusters().size() < 2)
			throw MethodException("No groups defined for -cmh-ijk test, i.e. K=1\n");

			/////////////////////////
			// Autosomal or haploid?

			bool Xchr = false, haploid = false;
			if (opts::_CHRX_ == mark->getChrom())
				Xchr = true;
			else if (opts::_CHRXY_ == mark->getChrom())
				haploid = true;

			if (haploid || Xchr)
				throw MethodException("-cmh-ijk / -cmh-ordinal cannot handle X/Y markers currently...\n");

			vector<int> X(0); // SNP
			vector<int> Y(0); // Cluster
			vector<int> Z(0); // Phenotype

			int l = mark->getLoc();
			for(int s = 0; s < data_set->num_inds(); s++){
				Sample* person = data_set->get_sample(s);

				if (!person->isEnabled() || samp_flags.at(s)){
					// Next person
					continue;
				}

				// Only consider individuals who have been assigned to a cluster
				if (!person->getAone(l) && !person->getAtwo(l)){
					X.push_back(1);
					X.push_back(1);
				} else if (!person->getAone(l) && person->getAtwo(l)){
					X.push_back(1);
					X.push_back(2);
				} else if (person->getAone(l) && person->getAtwo(l) && !person->getAmissing(l)){
					X.push_back(2);
					X.push_back(2);
				} else {
					continue;
				}
				Y.push_back(samp_clusters[person]); //currently groups are strings
				Y.push_back(samp_clusters[person]);
				if (options.doCMHOrdinal())
					Z.push_back((int) person->getPheno());
				else {
					if (person->getPheno() == 2) {
						Z.push_back(2);
						Z.push_back(2);
					} else {
						Z.push_back(1);
						Z.push_back(1);
					}
				}

			}

			vector<double> res;

			if (options.doCMHOrdinal())
				res = calcMantelHaenszel_ORD(X, Z, Y);
			else
				res = calcMantelHaenszel_IxJxK(X, Y, Z);

			chisq = res.at(0);
			pval = Helpers::p_from_chi(res.at(0), res.at(1));

	}

}

/*
 * CMH using 2x2xK model
 *
 */
void CMH::calcMantelHaenszel_2x2xK(Marker* mark) {
	map<Sample*, int> samp_clusters = options.getSampleClusters();
	// Should we perform BD test (K>1)
	// nk = number of clusters
	if (options.getClusters().size() < 2)
		options.setBreslowDay(false);

		// Warnings,
		if (options.doBreslowDay() && options.getClusters().size() > 10)
			opts::printLog("** Warning ** Breslow-Day statistics require large N per cluster ** \n");

	double zt = Helpers::ltqnorm(1 - (1 - options.getCI()) / 2);

	// Cochran-Mantel-Haenszel 2x2xK test
	// Disease X allele X strata

	// Calculate mean of 11 cell for each strata
	vector<double> mean_11(options.getClusters().size(), 0);
	vector<double> var_11(options.getClusters().size(), 0);

	// Calculate statistic
	vector<double> n_11(options.getClusters().size(), 0);
	vector<double> n_12(options.getClusters().size(), 0);
	vector<double> n_21(options.getClusters().size(), 0);
	vector<double> n_22(options.getClusters().size(), 0);

	// Disease marginals
	vector<double> n_1X(options.getClusters().size(), 0); // disease
	vector<double> n_2X(options.getClusters().size(), 0); // no disease

	vector<double> n_X1(options.getClusters().size(), 0); // F allele
	vector<double> n_X2(options.getClusters().size(), 0); // T allele

	vector<double> n_TT(options.getClusters().size(), 0); // Total allele count


	/////////////////////////
	// Autosomal or haploid?

	bool X = false, haploid = false;
	if (opts::_CHRX_ == mark->getChrom())
		X = true;
	else if (opts::_CHRXY_ == mark->getChrom())
		haploid = true;

	////////////////////////
	// Consider each person

	int l = mark->getLoc();
	for(int s = 0; s < data_set->num_inds(); s++){
		Sample* pperson = data_set->get_sample(s);

		bool s1 = pperson->getAone(l);
		bool s2 = pperson->getAtwo(l);
		bool s3 = pperson->getAmissing(l);

		// Affected individuals
		if (pperson->getAffected() && pperson->isEnabled() && !samp_flags.at(s)){

			// Haploid?
			if (haploid || (X && pperson->getSex())) {

				// Allelic marginal
				if (!s1) {
					// FF hom
					n_11.at(samp_clusters[pperson])++;
					n_X1.at(samp_clusters[pperson])++;
				} else {
					if (s2 && s3) // FT
					{
						continue; // skip missing genotypes
					} else if(s2 && !s3)// TT
					{
						n_12.at(samp_clusters[pperson])++;
						n_X2.at(samp_clusters[pperson])++;
					}
				}

				// Disease marginal
				n_1X.at(samp_clusters[pperson])++;
				n_TT.at(samp_clusters[pperson])++;

			} else // autosomal
			{

				// Allelic marginal
				if (!s1) {
					if (!s2) // FF hom
					{
						n_11.at(samp_clusters[pperson]) += 2;
						n_X1.at(samp_clusters[pperson]) += 2;
					} else {
						n_11.at(samp_clusters[pperson])++; // FT het
						n_12.at(samp_clusters[pperson])++;
						n_X1.at(samp_clusters[pperson])++;
						n_X2.at(samp_clusters[pperson])++;
					}
				} else {
					if (s2 && s3) // FT
					{
						continue; // skip missing genotypes
					} else if(s2 && !s3) // TT
					{
						n_12.at(samp_clusters[pperson]) += 2;
						n_X2.at(samp_clusters[pperson]) += 2;
					}
				}

				// Disease marginal
				n_1X.at(samp_clusters[pperson]) += 2;
				n_TT.at(samp_clusters[pperson]) += 2;

			} // end autosomal

		} else if (pperson->isEnabled() && !samp_flags[s])
		{

			// Haploid?
			if (haploid || (X && pperson->getSex())) {

				// Allelic marginal
				if (!s1) {
					// FF hom
					n_21.at(samp_clusters[pperson])++;
					n_X1.at(samp_clusters[pperson])++;
				} else {
					if (s2 && s3) // FT
					{
						continue; // skip missing genotypes
					} else if(s2 && !s3) // TT
					{
						n_22.at(samp_clusters[pperson])++;
						n_X2.at(samp_clusters[pperson])++;
					}
				}

				// Disease marginal
				n_2X.at(samp_clusters[pperson])++;
				n_TT.at(samp_clusters[pperson])++;

			} else // autosomal
			{
				// Allelic marginal
				if (!s1) {
					if (!s2) // FF
					{
						n_X1.at(samp_clusters[pperson]) += 2;
						n_21.at(samp_clusters[pperson]) += 2;
					} else {
						n_X1.at(samp_clusters[pperson])++;
						n_X2.at(samp_clusters[pperson])++;
						n_21.at(samp_clusters[pperson])++;
						n_22.at(samp_clusters[pperson])++;
					}
				} else {
					if (s2 && s3) // FT
					{
						continue; // skip missing genotypes
					} else if(s2 && !s3) // TT
					{
						n_X2.at(samp_clusters[pperson]) += 2;
						n_22.at(samp_clusters[pperson]) += 2;
					}
				}

				// disease marginal
				n_2X.at(samp_clusters[pperson]) += 2;
				n_TT.at(samp_clusters[pperson]) += 2;

			} // end autosomal
		} // end unaffected


	} // count next individual


	// Finished iterating over individuals: cluster needs at least 2
	// nonmissing individuals

	vector<bool> validK(options.getClusters().size(), false);
	for (unsigned int k = 0; k < options.getClusters().size(); k++)
		if (n_TT.at(k) >= 2)
			validK.at(k) = true;

	for (unsigned int k = 0; k < options.getClusters().size(); k++) {
		if (validK.at(k)) {
			mean_11.at(k) = (n_X1.at(k) * n_1X.at(k)) / n_TT.at(k);
			var_11.at(k) = (n_X1.at(k) * n_X2.at(k) * n_1X.at(k) * n_2X.at(k)) / (n_TT.at(k)
					* n_TT.at(k) * (n_TT.at(k) - 1));
		}
	}

	double CMH = 0;
	double denom = 0;
	for (unsigned int k = 0; k < options.getClusters().size(); k++) {
		if (validK.at(k)) {
			CMH += n_11.at(k) - mean_11.at(k);
			denom += var_11.at(k);
		}
	}

	CMH *= CMH;
	CMH /= denom;

	// MH Odds ratio & CI
	double R = 0;
	double S = 0;
	vector<double> r2(options.getClusters().size());
	vector<double> s2(options.getClusters().size());

	for (unsigned int k = 0; k < options.getClusters().size(); k++) {
		if (validK.at(k)) {
			r2.at(k) = (n_11.at(k) * n_22.at(k)) / n_TT.at(k);
			s2.at(k) = (n_12.at(k) * n_21.at(k)) / n_TT.at(k);
			R += r2.at(k);
			S += s2.at(k);
		}
	}
	OR = R / S;

	double v1 = 0, v2 = 0, v3 = 0;
	for (unsigned int k = 0; k < options.getClusters().size(); k++) {
		if (validK.at(k)) {
			v1 += (1 / n_TT.at(k)) * (n_11.at(k) + n_22.at(k)) * r2.at(k);
			v2 += (1 / n_TT.at(k)) * (n_12.at(k) + n_21.at(k)) * s2.at(k);
			v3 += (1 / n_TT.at(k)) * ((n_11.at(k) + n_22.at(k)) * s2.at(k) + (n_12.at(k)
					+ n_21.at(k)) * r2.at(k));
		}
	}

	double varOR = (1 / (2 * R * R)) * v1 + (1 / (2 * S * S)) * v2 + (1
			/ (2 * R * S)) * v3;

	OR_lower = exp(log(OR) - zt * sqrt(varOR));
	OR_upper = exp(log(OR) + zt * sqrt(varOR));

	chisq = CMH;
	pval = Helpers::chiprobP(chisq, 1);

		// Optional Breslow-Day test of homogeneity of odds ratios
		if (options.doBreslowDay()){

			double amax;
			double bb;
			double determ;
			double as_plus;
			double as_minus;
			double Astar;
			double Bstar;
			double Cstar;
			double Dstar;
			double Var;
			double BDX2 = 0;
			int df = 0;
			for (unsigned int k = 0; k < options.getClusters().size(); k++) {
				if (validK.at(k)) {
					df++;
					amax = (n_1X.at(k) < n_X1.at(k)) ? n_1X.at(k) : n_X1.at(k);
					bb = n_2X.at(k) + n_1X.at(k) * OR - n_X1.at(k) * (1 - OR);
					determ = sqrt(bb * bb + 4 * (1 - OR) * OR * n_1X.at(k)
							* n_X1.at(k));
					as_plus = (-bb + determ) / (2 - 2 * OR);
					as_minus = (-bb - determ) / (2 - 2 * OR);
					Astar = as_minus <= amax && as_minus >= 0 ? as_minus
							: as_plus;
					Bstar = n_1X.at(k) - Astar;
					Cstar = n_X1.at(k) - Astar;
					Dstar = n_2X.at(k) - n_X1.at(k) + Astar;
					Var = 1 / (1 / Astar + 1 / Bstar + 1 / Cstar + 1
							/ Dstar);
					BDX2 += ((n_11.at(k) - Astar) * (n_11.at(k) - Astar)) / Var;
				}
			}

			double BDp = Helpers::chiprobP(BDX2, df - 1);

			chisq_bd = BDX2;
			pval_bd = BDp;
		}

}

/*
 * calculate CMH using IxJxK model
 */
vector<double> CMH::calcMantelHaenszel_IxJxK(vector<int> & X, vector<int> & Y,
		vector<int> & Z) {

	if (X.size() != Y.size() || Y.size() != Z.size() || X.size() != Z.size())
		throw MethodException("Internal problem:\n  problem in calcMantelHaenszel_IxJxK()...uneven input columns");

	// Determine unique elements
	int nx = 0, ny = 0, nz = 0;
	map<int, int> mx;
	map<int, int> my;
	map<int, int> mz;

	for (unsigned int i = 0; i < X.size(); i++) {
		if (mx.find(X.at(i)) == mx.end())
			mx.insert(make_pair(X.at(i), nx++));

		if (my.find(Y.at(i)) == my.end())
			my.insert(make_pair(Y.at(i), ny++));

		if (mz.find(Z.at(i)) == mz.end())
			mz.insert(make_pair(Z.at(i), nz++));
	}

	// Generic function to calculate generalized IxJxK CMH
	// Assumes no missing data

	vector<vector<double> > N(nz); // observed counts
	vector<vector<double> > U(nz); // expected
	vector<vector<vector<double> > > V(nz); // variance matrix

	vector<vector<double> > Tx(nz); // marginal totals
	vector<vector<double> > Ty(nz); // ..
	vector<double> T(nz); // totals (per K)

	for (int k = 0; k < nz; k++) {
		Tx.at(k).resize(nx);
		Ty.at(k).resize(ny);

		N.at(k).resize((nx - 1) * (ny - 1));
		U.at(k).resize((nx - 1) * (ny - 1));
		V.at(k).resize((nx - 1) * (ny - 1));
		for (int k2 = 0; k2 < (nx - 1) * (ny - 1); k2++) {
			N.at(k).at(k2) = U.at(k).at(k2) = 0;
			V.at(k).at(k2).resize((nx - 1) * (ny - 1));
			for (int k3 = 0; k3 < (nx - 1) * (ny - 1); k3++)
				V.at(k).at(k2).at(k3) = 0;
		}
	}

	// Consider each observation
	for (unsigned int i = 0; i < X.size(); i++) {
		int vx = mx.find(X.at(i))->second;
		int vy = my.find(Y.at(i))->second;
		int vz = mz.find(Z.at(i))->second;

		// exclude nx + ny (upper limits)
		if (vx < nx - 1 && vy < ny - 1)
			N.at(vz).at(vx + vy * (nx - 1))++;

		Tx.at(vz).at(vx)++;
		Ty.at(vz).at(vy)++;
		T.at(vz)++;
	}

	// Determine valid clusters (at least 2 people)
	vector<bool> validK(options.getClusters().size(), false);
	for (unsigned int k = 0; k < options.getClusters().size(); k++)
		if (T.at(k) >= 2)
			validK.at(k) = true;

	// Calculate expecteds
	for (int k = 0; k < nz; k++) {
		if (validK.at(k)) {
			for (int ix = 0; ix < nx - 1; ix++)
				for (int iy = 0; iy < ny - 1; iy++) {
					U.at(k).at(ix + iy * (nx - 1)) = (Tx.at(k).at(ix) * Ty.at(k).at(iy)) / T.at(k);

					for (int ix2 = 0; ix2 < nx - 1; ix2++)
						for (int iy2 = 0; iy2 < ny - 1; iy2++) {
							int dx = 0;
							int dy = 0;
							if (ix == ix2)
								dx = 1;
							if (iy == iy2)
								dy = 1;
							V.at(k).at(ix + iy * (nx - 1)).at(ix2 + iy2 * (nx - 1))
									= ((Tx.at(k).at(ix) * (dx * T.at(k) - Tx.at(k).at(ix2))
											* Ty.at(k).at(iy) * (dy * T.at(k)
											- Ty.at(k).at(iy2))) / (T.at(k) * T.at(k)
											* (T.at(k) - 1)));
							if (ix == ix2 && iy == iy2)
								V.at(k).at(ix + iy * (nx - 1)).at(ix2 + iy2 * (nx - 1))
										= abs(V.at(k).at(ix + iy * (nx - 1)).at(ix2 + iy2 * (nx - 1)));
						}
				}
		}
	}

	vector<vector<double> > V0((nx - 1) * (ny - 1));
	for (int k2 = 0; k2 < (nx - 1) * (ny - 1); k2++)
		V0.at(k2).resize((nx - 1) * (ny - 1));
	vector<double> N0((nx - 1) * (ny - 1));
	vector<double> U0((nx - 1) * (ny - 1));

	// Sum N, U and V over K
	for (int k = 0; k < nz; k++) {
		if (validK.at(k)) {
			for (int i = 0; i < (nx - 1) * (ny - 1); i++) {
				N0.at(k) += N.at(k).at(i);
				U0.at(i) += U.at(k).at(i);

				for (int i2 = 0; i2 < (nx - 1) * (ny - 1); i2++)
					V0.at(k).at(i2) += V.at(k).at(i).at(i2);
			}
		}
	}

	vector<double> tmp1((nx - 1) * (ny - 1), 0);
	vector<double> tmp2((nx - 1) * (ny - 1), 0);
	V0 = Helpers::svd_inverse(V0);
	for (int i = 0; i < (nx - 1) * (ny - 1); i++)
		tmp1.at(i) = N0.at(i) - U0.at(i);

	// Matrix mult -- rows by columns

	for (int i = 0; i < (nx - 1) * (ny - 1); i++)
		for (int j = 0; j < (nx - 1) * (ny - 1); j++)
			tmp2.at(j) += tmp1.at(i) * V0.at(i).at(j);

	vector<double> result(2);

	// CMH Chi-square
	result.at(0) = 0;
	for (int i = 0; i < (nx - 1) * (ny - 1); i++)
		result.at(0) += tmp2.at(i) * tmp1.at(i);

	// DF
	result.at(1) = (nx - 1) * (ny - 1);
	return result;

}



vector<double> CMH::calcMantelHaenszel_ORD(vector<int> & X, vector<int> & Y,
		vector<int> & Z) {
vector<double> temp;
	// X is SNP coding
	// Y is phenotype (assumed to be ordinal, integers)
	// Z is cluster

	if (X.size() != Y.size() || Y.size() != Z.size() || X.size() != Z.size())
		throw MethodException("Internal problem:\n  problem in calcMantelHaenszel_ORD()...uneven input columns");

	// Determine unique elements
	int nx = 0, ny = 0, nz = 0;
	map<int, int> mx;
	map<int, int> my;
	map<int, int> mz;

	for (unsigned int i = 0; i < X.size(); i++) {
		if (mx.find(X.at(i)) == mx.end())
			mx.insert(make_pair(X.at(i), nx++));

		if (my.find(Y.at(i)) == my.end())
			my.insert(make_pair(Y.at(i), ny++));

		if (mz.find(Z.at(i)) == mz.end())
			mz.insert(make_pair(Z.at(i), nz++));
	}

	// Generic function to calculate generalized ordinal IxJxK CMH
	// Assumes no missing data

	vector<vector<double> > N(nz); // observed counts
	vector<vector<double> > U(nz); // expected
	vector<vector<vector<double> > > V(nz); // variance matrix

	vector<vector<double> > Tx(nz); // marginal totals
	vector<vector<double> > Ty(nz); // ..
	vector<double> T(nz); // totals (per K)

	for (int k = 0; k < nz; k++) {
		Tx.at(k).resize(nx);
		Ty.at(k).resize(ny);

		N.at(k).resize((nx - 1) * (ny - 1));
		U.at(k).resize((nx - 1) * (ny - 1));
		V.at(k).resize((nx - 1) * (ny - 1));
		for (int k2 = 0; k2 < (nx - 1) * (ny - 1); k2++) {
			N.at(k).at(k2) = U.at(k).at(k2) = 0;
			V.at(k).at(k2).resize((nx - 1) * (ny - 1));
			for (int k3 = 0; k3 < (nx - 1) * (ny - 1); k3++)
				V.at(k).at(k2).at(k3) = 0;
		}
	}

	// Create counts
	// Consider each observation
	for (unsigned int i = 0; i < X.size(); i++) {
		int vx = mx.find(X.at(i))->second;
		int vy = my.find(Y.at(i))->second;
		int vz = mz.find(Z.at(i))->second;

		// exclude nx + ny (upper limits)
		if (vx < nx - 1 && vy < ny - 1)
			N.at(vz).at(vx + vy * (nx - 1))++;

		Tx.at(vz).at(vx)++;
		Ty.at(vz).at(vy)++;
		T.at(vz)++;
	}

	// Determine valid clusters (at least 2 people)
	vector<bool> validK(options.getClusters().size(), false);
	for (unsigned int k = 0; k < options.getClusters().size(); k++)
		if (T.at(k) >= 2)
			validK.at(k) = true;

	// Calculate expecteds
	for (int k = 0; k < nz; k++) {
		if (validK.at(k)) {
			for (int ix = 0; ix < nx - 1; ix++)
				for (int iy = 0; iy < ny - 1; iy++) {
					U.at(k).at(ix + iy * (nx - 1)) = (Tx.at(k).at(ix) * Ty.at(k).at(iy)) / T.at(k);
				}
		}
	}

	return temp;
}
}
