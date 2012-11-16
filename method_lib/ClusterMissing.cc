/**********************************************************************************
*                       Marker Genotype Efficiency Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all genotypes and generates a genotype efficiency for all markers.
*
*
* Files generated:
*	percent_breakdown_by_marker.txt
*	percent_breakdown_by_chrom.txt
*       post_marker_geno_eff_filter_summary.txt
*
*File: MarkerGenoEff.cc
**********************************************************************************/


#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include <set>
#include "ClusterMissing.h"
#include "Options.h"
#include "General.h"
#include "Helpers.h"
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
namespace Methods{
string ClusterMissing::stepname = "cluster-missing";

void ClusterMissing::FilterSummary(){
}

void ClusterMissing::PrintSummary(){

}

void ClusterMissing::filterOne(int m){
}

void ClusterMissing::filter(){
}

int ClusterMissing::calcPairLocus(int s1, int s2, int m){
	Marker* mark = data_set->get_locus(m);
	Sample* samp1 = data_set->get_sample(s1);
	Sample* samp2 = data_set->get_sample(s2);
	int count = 0;

	if(mark != NULL && samp1 != NULL && samp2 != NULL){
		int mloc = mark->getLoc();
		if(mark->isEnabled() && samp1->isEnabled() && samp2->isEnabled() && !samp1->isExcluded() && !samp2->isExcluded()){
			if((samp1->getAone(mloc) && samp1->getAtwo(mloc) && samp1->getAmissing(mloc))
					|| (samp2->getAone(mloc) && samp2->getAtwo(mloc) && samp2->getAmissing(mloc))){
				return -1;
			}
			if(samp1->getAone(mloc) == samp2->getAone(mloc)){
				count++;
			}
			if(samp1->getAtwo(mloc) == samp2->getAtwo(mloc)){
				count++;
			}
		}
	}
	else{
		return -1;
	}
	return count;
}

double ClusterMissing::calcPairAverage(int s1, int s2){
	comparisons = 0;
//	int prev_base = 0;
//	int prev_chrom = -1;
	Sample* samp1 = data_set->get_sample(s1);
	Sample* samp2 = data_set->get_sample(s2);

	int sum = 0;
	int num_loci = 0;

	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	int msize = good_markers.size();//data_set->num_loci();
	for(int m = 0; m < msize; m++){
		Marker* mark = good_markers[m];//data_set->get_locus(m);
		if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
			int mloc = mark->getLoc();
			if((samp1->getAone(mloc) && samp1->getAtwo(mloc) && samp1->getAmissing(mloc))
					|| (samp2->getAone(mloc) && samp2->getAtwo(mloc) && samp2->getAmissing(mloc))){
				continue;
			}
			num_loci++;
			if(samp1->getAone(mloc) == samp2->getAone(mloc)){
				sum++;
			}
			if(samp1->getAtwo(mloc) == samp2->getAtwo(mloc)){
				sum++;
			}

		}
	}
	double avg = 0.0f;
	if(num_loci > 0){
		avg = (double) sum / (double) num_loci;
	}
	comparisons = num_loci;
	return avg;
}

void ClusterMissing::calcOne(int m){
//	Marker* mark = data_set->get_locus(m);
//	int mloc = mark->getLoc();


}

void ClusterMissing::calculate(string file_prefix, string file_suffix){
	vector<vector<double> > mdist;

	vector<vector<int > > cl;
	int ninds = data_set->num_inds();
	int np = (int)((double)(ninds*(ninds-1))/(double)2);

	//if no file defining pairs
	for(int i = 0; i < ninds; i++){
		vector<int> t(1);
		t[0] = i;
		cl.push_back(t);
	}

	vector<vector<bool> > pairable(ninds);
	for(int i = 0; i < ninds; i++){
		vector<bool> tmp(ninds, true);
		pairable[i] = tmp;
	}

	opts::printLog("Clustering individuals based on genome-wide IBM\n");

	set<int> selcon;

	mdist.resize(ninds);
	for(int j = 0; j < ninds; j++){
		mdist[j].resize(j);
	}

	vector<double> prop_sig_diff(ninds);
	double pv = 0;
	double merge_p = 0;
	if(true){ //calculate
		int c = 0;
		int c2 = 0;

		for(int i1 = 0; i1 < ninds - 1; i1++){
			for(int i2 = i1 + 1; i2 < ninds; i2++){
				if(c == c2 || c == np){
					cout << "IBM calculation: " << c++ << " of " << np << "             \r";
					cout.flush();
					c2+= 100;
				}
				else{
					++c;
				}

				double dst = calcGenomeIBM(data_set->get_sample(i1), data_set->get_sample(i2));
				mdist[i2][i1] = dst;
//				cout << i2 << ":" << i1 << ":" << dst << endl;
				if(pv < merge_p && Helpers::realnum(pv)){
					pairable[i1][i2] = pairable[i2][i1] = false;
					prop_sig_diff[i1]++;
					prop_sig_diff[i2]++;
				}

			}
		}
		cout << "\n";
	}

	//Display matrix of IBS distance
	string f = file_prefix + "mdist.missing" + file_suffix;
	ofstream MAT(f.c_str(), ios::out);
	MAT.clear();
	opts::printLog("Writing IBM distance matrix.\n");
	for(int i = 0; i < (int)mdist.size(); i++){
		for(int j = 0; j < (int)mdist.size(); j++){
			//if distance matrix
			//else
			if(i > j){
				MAT << mdist[i][j] << " ";
			}
			else if(i == j){
				MAT << 1 << " ";
			}
			else{
				MAT << mdist[j][i] << " ";
			}
		}
		MAT << "\n";
	}
	MAT.close();

	int c = 1;

	bool done = false;

	// Matrix of solutions
	vector<vector<int> > sol(ninds);
	for(int i =0; i < ninds; i++){
		sol[i].resize(ninds);
	}

	vector<double> hist(1);

	// Build solution
//cl?????
	for(int i = 0; i < (int)cl.size(); i++){
		for(int j = 0; j < (int)cl[i].size(); j++){
			sol[cl[i][j]][0] = i;
		}
	}

	opts::printLog("Writing cluster progress.\n");
	f = file_prefix + "cluster0" + file_suffix;
	ofstream CLST(f.c_str(), ios::out);
	CLST.clear();

	while(!done){
		double dmin = -999;

		int imin = -1;
		int jmin = -1;

		// 1. Find min/max distance between pairable clusters
		for(int i = 0; i < (int)cl.size()-1; i++){
			for(int j = i + 1; j < (int)cl.size(); j++){
				// Cluster on IBS
				double d = cldist(mdist, cl[i],cl[j]);

				// Are these individuals/clusters more similar AND pairable?
				if(d > dmin && pairable_cluster(pairable, cl[i], cl[j])){
					//And will the max cluster size requirement be fulfilled?
					if(options.getMaxClusterSize() == 0 || ((int)(cl[i].size()+cl[j].size()) <= options.getMaxClusterSize())){ //0
						//And will the basic phenotypic matching requirement be fulfilled?
						if(!options.getClusterOnPheno()){// || (!homogeneous_clusters((*this), cl[i], cl[j])))
							//What about the --mcc clustering
							if(!options.getClusterOnMcc()){// || spec_clusters((*this), cl[i], cl[j])) //true/false
								//And what about pick1 constraints? (this must be final constraint)
								if(!options.getClusterSelcon()){ //||selcon_inds((*this), cl[i], cl[j], selcon)){ //true/false
									imin = i;
									jmin = j;
									dmin = d;
								}
							}
						}
					}
				}
			}
		}
		if(imin != -1){
			hist.push_back(dmin);

			for(int j = 0; j < (int)cl[jmin].size(); j++){
				cl[imin].push_back(cl[jmin][j]);
			}
			cl.erase(cl.begin()+jmin);
			if(cl.size() == 1 || (int)cl.size() == options.getMaxClusterN()){ //-1
				done = true;
			}

			CLST << "Merge step " << c << "\t" << hist[c];

			// Build solution
			for(int i = 0; i < (int)cl.size(); i++){
				for(int j = 0; j < (int)cl[i].size(); j++){
					sol[cl[i][j]][c] = i;
				}
			}

			//Calculate average within/between cluster distances
			double between = 0, within = 0;
			int withinN = 0, betweenN = 0;

			for(int j1 = 0; j1 < (int)sol.size(); j1++){
				for(int j2 = 0; j2 < (int)sol.size(); j2++){
					if(j1 < j2){
						if(sol[j1][c] == sol[j2][c]){
							within += mdist[j2][j1];
							withinN++;
						}
						else{
							between += mdist[j2][j1];
							betweenN++;
						}
					}
				}
			}
			CLST << "\t" << between/(double)betweenN
				<< "\t" << within/(double)withinN
				<< "\t" << (between/(double)betweenN) / (within/(double)withinN)
				<< "\n";

			//Next merge
			c++;

		}
		//Did we get a merge?
		if(imin == -1){
			done = true;
			//goto done_making_clusters;
		}
	}

	CLST.close();
	//Best solution is final solution
////	int best = hist.size() - 1;

	opts::printLog("Writing cluster solution (3)\n");
	f = file_prefix + "cluster3.missing" + file_suffix;
	CLST.open(f.c_str(), ios::out);
	CLST.clear();

	for(int j = 0; j < (int)sol.size(); j++){
		//Display...
		CLST << data_set->get_sample(j)->getFamID() << " "
			<< data_set->get_sample(j)->getInd() << "\t";
		for(int i = 0; i < (int)sol[0].size(); i++){
			CLST << sol[j][i] << " ";
		}
		CLST << "\n";
	}
	CLST << "\n";
	CLST.close();

}

double ClusterMissing::cldist(vector<vector<double> > & d,
          vector<int> & a,
          vector<int> & b)
{
  // Compare based on first metric, but also return paired second
  double l;
  l = a[0]>b[0] ? d[a[0]][b[0]] : d[b[0]][a[0]];

  for (int i=0; i< (int)a.size(); i++)
    for (int j=0; j< (int)b.size(); j++)
      {

    if ( a[i] > b[j] )
      {
        if ( d[a[i]][b[j]] < l ) l = d[a[i]][b[j]];
      }
    else
      {
        if ( d[b[j]][a[i]] < l ) l = d[b[j]][a[i]];
      }

      }
  return l;
}

bool ClusterMissing::pairable_cluster(vector<vector<bool> > & pairable, vector<int> & a, vector<int> & b)
{
  for (int i=0; i< (int)a.size(); i++)
    for (int j=0; j< (int)b.size(); j++)
       if (!pairable[a[i]][b[j]]) return false;
  return true;
}


double ClusterMissing::calcGenomeIBM(Sample* one, Sample* two){
	int cnt = 0;

	int msize = data_set->num_loci();

	int prev_bploc = 0;
	int prev_chrom = -1;

	int goodsnps = 0;
	for(int i = 0; i < msize; i++){
		Marker* mark = data_set->get_locus(i);
		if(mark->isEnabled() && Helpers::isValidMarker(mark, &options, prev_chrom, prev_bploc)){
			int mloc = mark->getLoc();

			if(one->getAone(mloc) && one->getAtwo(mloc) && one->getAmissing(mloc)){
				if(! (two->getAone(mloc) && two->getAtwo(mloc) && two->getAmissing(mloc))){
					cnt++;
				}
			}
			else if(two->getAone(mloc) && two->getAtwo(mloc) && two->getAmissing(mloc)){
				cnt++;
			}
			goodsnps++;
		}
	}

	return (1.0 - ( (double)cnt / (double) goodsnps ));
}

void ClusterMissing::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

}

}
