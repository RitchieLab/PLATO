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
	Sample* samp1 = data_set->get_sample(s1);
	Sample* samp2 = data_set->get_sample(s2);

	int sum = 0;
	int num_loci = 0;

	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	int msize = good_markers.size();
	for(int m = 0; m < msize; m++){
		Marker* mark = good_markers.at(m);
		if(mark->isEnabled()){
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
}

void ClusterMissing::calculate(string file_prefix, string file_suffix){
	vector<vector<double> > mdist;

	vector<vector<int > > cl;
	int ninds = data_set->num_inds();
	int np = (int)((double)(ninds*(ninds-1))/(double)2);

	//if no file defining pairs
	for(int i = 0; i < ninds; i++){
		vector<int> t(1);
		t.at(0) = i;
		cl.push_back(t);
	}

	vector<vector<bool> > pairable(ninds);
	for(int i = 0; i < ninds; i++){
		vector<bool> tmp(ninds, true);
		pairable.at(i) = tmp;
	}

	opts::printLog("Clustering individuals based on genome-wide IBM\n");

	set<int> selcon;

	mdist.resize(ninds);
	for(int j = 0; j < ninds; j++){
		mdist.at(j).resize(j);
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
				mdist.at(i2).at(i1) = dst;
				if(pv < merge_p && Helpers::realnum(pv)){
					pairable.at(i1).at(i2) = pairable.at(i2).at(i1) = false;
					prop_sig_diff.at(i1)++;
					prop_sig_diff.at(i2)++;
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
				MAT << mdist.at(i).at(j) << " ";
			}
			else if(i == j){
				MAT << 1 << " ";
			}
			else{
				MAT << mdist.at(j).at(i) << " ";
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
		sol.at(i).resize(ninds);
	}

	vector<double> hist(1);

	// Build solution
	for(int i = 0; i < (int)cl.size(); i++){
		for(int j = 0; j < (int)cl.at(i).size(); j++){
			sol.at(cl.at(i).at(j)).at(0) = i;
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
				double d = cldist(mdist, cl.at(i),cl.at(j));

				// Are these individuals/clusters more similar AND pairable?
				if(d > dmin && pairable_cluster(pairable, cl.at(i), cl.at(j))){
					//And will the max cluster size requirement be fulfilled?
					if(options.getMaxClusterSize() == 0 || ((int)(cl.at(i).size()+cl.at(j).size()) <= options.getMaxClusterSize())){
						//And will the basic phenotypic matching requirement be fulfilled?
						if(!options.getClusterOnPheno()){
							//What about the --mcc clustering
							if(!options.getClusterOnMcc()){
								//And what about pick1 constraints? (this must be final constraint)
								if(!options.getClusterSelcon()){
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

			for(int j = 0; j < (int)cl.at(jmin).size(); j++){
				cl.at(imin).push_back(cl.at(jmin).at(j));
			}
			cl.erase(cl.begin()+jmin);
			if(cl.size() == 1 || (int)cl.size() == options.getMaxClusterN()){ //-1
				done = true;
			}

			CLST << "Merge step " << c << "\t" << hist.at(c);

			// Build solution
			for(int i = 0; i < (int)cl.size(); i++){
				for(int j = 0; j < (int)cl.at(i).size(); j++){
					sol.at(cl.at(i).at(j)).at(c) = i;
				}
			}

			//Calculate average within/between cluster distances
			double between = 0, within = 0;
			int withinN = 0, betweenN = 0;

			for(int j1 = 0; j1 < (int)sol.size(); j1++){
				for(int j2 = 0; j2 < (int)sol.size(); j2++){
					if(j1 < j2){
						if(sol.at(j1).at(c) == sol.at(j2).at(c)){
							within += mdist.at(j2).at(j1);
							withinN++;
						}
						else{
							between += mdist.at(j2).at(j1);
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

	opts::printLog("Writing cluster solution (3)\n");
	f = file_prefix + "cluster3.missing" + file_suffix;
	CLST.open(f.c_str(), ios::out);
	CLST.clear();

	for(int j = 0; j < (int)sol.size(); j++){
		//Display...
		CLST << data_set->get_sample(j)->getFamID() << " "
			<< data_set->get_sample(j)->getInd() << "\t";
		for(int i = 0; i < (int)sol.at(0).size(); i++){
			CLST << sol.at(j).at(i) << " ";
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
  l = a.at(0)>b.at(0) ? d.at(a.at(0)).at(b.at(0)) : d.at(b.at(0)).at(a.at(0));

  for (int i=0; i< (int)a.size(); i++)
    for (int j=0; j< (int)b.size(); j++)
      {

    if ( a.at(i) > b.at(j) )
      {
        if ( d.at(a.at(i)).at(b.at(j)) < l ) l = d.at(a.at(i)).at(b.at(j));
      }
    else
      {
        if ( d.at(b.at(j)).at(a.at(i)) < l ) l = d.at(b.at(j)).at(a.at(i));
      }

      }
  return l;
}

bool ClusterMissing::pairable_cluster(vector<vector<bool> > & pairable, vector<int> & a, vector<int> & b)
{
  for (int i=0; i< (int)a.size(); i++)
    for (int j=0; j< (int)b.size(); j++)
       if (!pairable.at(a.at(i)).at(b.at(j))) return false;
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
