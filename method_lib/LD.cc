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
*File: LD.cc
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
#include "LD.h"
#include "Options.h"
#include "General.h"
#include "Helpers.h"

namespace Methods{
void LD::FilterSummary(){

	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void LD::PrintSummary(){
	int msize = markers->size();
	for(int m = 0; m < msize; m++){
		(*markers).at(m)->setFlag(false);
	}

}

void LD::filter(){
}

//courtesy of plink
void LD::calcSetMeanVariance(vector<Marker*> pSNP, vector<double> & mean, vector<vector<double> > & variance){
	int flagged_n = 0;

	int nss = pSNP.size();

	if (nss<1)
	   return;

	mean.resize(nss,0);

	vector<int> cnt(nss,0);
	variance.resize(nss);

	for (int j=0; j<nss; j++)
		variance.at(j).resize(nss,0);

	////////////////////////////
	//Iterate over SNPs in SET
	//
	int msize = pSNP.size();
	int fsize = families->size();
	for(int m = 0; m < msize; m++){
		Marker* mark = pSNP.at(m);
		if(mark->isEnabled() && !mark->isFlagged()){
			int mloc = mark->getLoc();
			for(int f = 0; f < fsize; f++){
				Family* fam = (*families).at(f);
				if(fam->isEnabled()){
					vector<Sample*>* founders = fam->getFounders();
					for(int s = 0; s < (int)founders->size(); s++){
						Sample* samp = (*founders).at(s);
						if(samp->isEnabled()){
							//increase general, flagged sample size (only count once)
							if(m == 0){
								flagged_n++;
							}
							if(!mark->isMicroSat()){
								if(samp->getAone(mloc)){
									if(samp->getAtwo(mloc) && !samp->getAmissing(mloc)){ //11 homozyg
										mean.at(m)++;
										cnt.at(m)++;
									}
								}
								else{
									cnt.at(m)++;
									if(!samp->getAtwo(mloc)){ //00 homozyg
										mean.at(m)--;
									}
								}
							}
							else{
								if(samp->getAbone(mloc) == samp->getAbtwo(mloc) && samp->getAbone(mloc) != -1){
									mean.at(m)++;
									cnt.at(m)++;
								}
								else{
									cnt.at(m)++;
								}
							}
						}
					}
				}
			}
		}
	}

	  // Having iterated over all individuals, we can now calculate the mean
	  //   // values, perform mean-substitution of missing data, and calculate the
	  //     // second order terms
	  //

	for(int j = 0; j < nss; j++){
		mean.at(j) /= (double)cnt.at(j);
	}


	  /////////////////////////////////////
	  //  // Iterate over pairs of SNPs in SET
	  //
	  //    // First SNP
	  //
	for(int m1 = 0; m1 < (int)pSNP.size(); m1++){
		Marker* mark1 = pSNP.at(m1);
		if(mark1->isEnabled() && !mark1->isFlagged()){
			int mloc1 = mark1->getLoc();
			for(int m2 = m1; m2 < (int)pSNP.size(); m2++){
				Marker* mark2 = pSNP.at(m2);
				if(mark2->isEnabled() && !mark2->isFlagged()){
					int mloc2 = mark2->getLoc();

					for(int f = 0; f < fsize; f++){
						Family* fam = (*families).at(f);
						if(fam->isEnabled()){
							vector<Sample*>* founders = fam->getFounders();
							for(int s = 0; s < (int)founders->size(); s++){
								Sample* samp = (*founders).at(s);
								if(samp->isEnabled()){
									double v1 = mean.at(m1), v2=mean.at(m2);
									//first snp
									if(!mark1->isMicroSat()){
										if(samp->getAone(mloc1)){
											if(samp->getAtwo(mloc1) && !samp->getAmissing(mloc1)){ //11 homozyg
												v1 = 1;
											}
										}
										else{
											if(!samp->getAtwo(mloc1)){ //00 homozyg
												v1 = -1;
											}
											else{
												v1 = 0; //01 het
											}
										}
									}
									else{
										if(samp->getAbone(mloc1) == samp->getAbtwo(mloc1) && samp->getAbone(mloc1) != -1){
											v1 = 1;
										}
										else if(samp->getAbone(mloc1) != -1){
											v1 = 0;
										}
									}


									//second snp
									if(!mark2->isMicroSat()){
										if(samp->getAone(mloc2)){
											if(samp->getAtwo(mloc2) && !samp->getAmissing(mloc2)){ //11 homozyg
												v2 = 1;
											}
										}
										else{
											if(!samp->getAtwo(mloc2)){ //00 homozyg
												v2 = -1;
											}
											else{
												v2 = 0; //01 het
											}
										}
									}
									else{
										if(samp->getAbone(mloc2) == samp->getAbtwo(mloc2) && samp->getAbone(mloc2) != -1){
											v2 = 1;
										}
										else if(samp->getAbone(mloc2) != -1){
											v2 = 0;
										}
									}
									//contribution to covariance term
									//mean substitution
									variance.at(m1).at(m2) += (v1 - mean.at(m1)) * (v2 - mean.at(m2));
								}
							}
						}
					}
				}
			}
		}
	}

	//make symmetric covariance matrix
	for(int i = 0; i < nss; i++){
		for(int j = i; j < nss; j++){
			variance.at(i).at(j) /= (double)(flagged_n);
			variance.at(j).at(i) = variance.at(i).at(j);
		}
	}

	return;
}

vector<bool> LD::vif_prune(vector<vector<double> > m, double threshold){
	//number of variables
	int p = m.size();
	//cout << "cur size = " << p << endl;
	vector<bool> cur(p, true);

	//this on ly is needed if we have 2+ snps
	if(p < 2){
		return cur;
	}

	vector<vector<double> > r = m;

	// Make 'm' a correlation matrix
	for (int i=0; i<p; i++)
		for (int j=0; j<p; j++)
	        r.at(i).at(j) = m.at(i).at(j) / sqrt(m.at(i).at(i) * m.at(j).at(j));

	// Number of excluded items
	int it = 0;

	// Any SNPs with zero variance should be automatically excluded
	for (int i=0; i<p; i++)
	    if ( r.at(i).at(i) == 0 || !Helpers::realnum(r.at(i).at(i)) )
	    {
	      cur.at(i) = false;
	      it++;
	    }


	// For any pair of perfectly correlated SNPs, exclude 1
	while(1)
	{
		bool done = true;
		for (int i=0; i<p-1; i++)
		{
			if (cur.at(i)){
			    for (int j=i+1;j<p; j++)
		        {
			        if (cur.at(j)) {
				        if ( fabs(r.at(i).at(j)) > options.getLDPWthreshold() )
				        {
					        cur.at(i) = false;
				            it++;
				            done = false;
				            break;
				        }
				    }
				}
	        }
    	}
     	if (done) break;
    }

	//skip vif caculation?
	if(options.doLDpairwise()){
		return cur;
	}

  // Calculate r^2 for each element versus all others
  // considering only the current non-pruned elements
  //
    while (1)
    {
      // Build correlation matrix all included items
		vector<vector<double> > u;

		for (int i=0;i<p;i++)
		{
			if ( cur.at(i) )
			{
				vector<double> mt;
		        for (int j=0;j<p;j++)
		 	       if ( cur.at(j) )
				       mt.push_back(r.at(i).at(j));
				u.push_back(mt);
			}
		}

		// Check enough markers left
		if (u.size()<2) break;

		// Get inverse
		u = Helpers::svd_inverse(u);//////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		//Calculate VIFs
        double maxVIF = 0;
        int maxI = 0;
        int cnt=0;
        for (int i=0;i<p;i++)
		    if ( cur.at(i) )
		    {
		        double vif = u.at(cnt).at(cnt);
		        if ( Helpers::dLess(maxVIF, vif) )
                {
		            maxVIF = vif;
		            maxI = i;
		        }
		        cnt++;
            }

		// How are we doing?
		if ( Helpers::dGreater(maxVIF,threshold) )
		{
			// exclude this item
			cur.at(maxI) = false;
		}
		else{
			break;
		}

		// Increase count of removed items
		it++;
		// Down to a single item or worse?
		if (it==p-1) break;
	}
	return cur;

}

vector<int> LD::getChromosomeMarkerRange(vector<Marker*>* marks, int chrom){
	vector<int> m(2);
	m.at(0) = -1;
	m.at(1) = -1;
	for(int i = 0; i < (int)marks->size(); i++){
		Marker* mark = (*marks).at(i);
		if(mark->isEnabled() && mark->getChrom() == chrom){
			if(i <= m.at(0) || m.at(0) == -1) m.at(0) = i;
			if(i >= m.at(1) || m.at(1) == -1) m.at(1) = i;
		}
	}
	return m;
}

void LD::setMarkerRange(){
		vector<int> m = getChromosomeMarkerRange(markers, run_chr);
		if(m.at(0) == -1 || m.at(1) == -1){
			cerr << "blah blah blha error getchrommarkrange\n";
			throw MethodException("blah blah blha error getchrommarkrange\n");
		}
		run_start = m.at(0);
		run_end = m.at(1);
}

void LD::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

	if(options.doLDCalcOnly()){
		calcLDStatistics();
		return;
	}

	map<int, int> chroms;
	int prev_base = 0;
	int prev_chrom = -1;
	for(int m = 0; m < (int)markers->size(); m++){
		if((*markers).at(m)->getChrom() < opts::_CHRX_ && Helpers::isValidMarker((*markers).at(m), &options, prev_base, prev_chrom)){
			if((*markers).at(m)->isEnabled() && !(*markers).at(m)->isFlagged()){
				chroms[(*markers).at(m)->getChrom()] = 1;
			}
		}
	}

	if(options.doLDpairwise() || options.doLDvif()){
		string filename = opts::_OUTPREFIX_ + "ld_keep" + options.getOut() + ".txt";
		if(!overwrite){
			filename += "." + getString<int>(order);
		}
		string filename2 = opts::_OUTPREFIX_ + "ld_removed_" + options.getOut() + ".txt";
		if(!overwrite){
			filename2 += "." + getString<int>(order);
		}
		ofstream LDIN(filename.c_str(), ios::out);
		ofstream LDOUT(filename2.c_str(), ios::out);
		if(!LDIN){
			opts::printLog("Error opening " + filename + ". Exiting!\n");
			throw MethodException("Error opening " + filename + ". Exiting!\n");
		}
		if(!LDOUT){
			opts::printLog("Error opening " + filename2 + ".  Exiting!\n");
			throw MethodException("Error opening " + filename2 + ". Exiting!\n");
		}

		vector<bool> include(markers->size(), true);
		map<int, int>::iterator chrom;
		opts::printLog("Total number of markers enabled at beginning of LD calculations: " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
		opts::printLog("Looking at ALL markers per chromosome, skipping those that are currently removed or disabled...\n");
		for(chrom = chroms.begin(); chrom != chroms.end(); chrom++){
			run_chr = chrom->first;

			setMarkerRange();
			int s1 = run_start;
			int s2 = run_start + options.getLDWin() - 1;
			while(s2 <= run_end){
				vector<Marker*> pSNP(0);
				for(int l = s1; l <= s2; l++){
					if(include.at(l) && (*markers).at(l)->isEnabled() && !(*markers).at(l)->isFlagged()){
						pSNP.push_back((*markers).at(l));
					}
				}

				if(pSNP.size() < 2){
					if(s2 == run_end){
						break;
					}
					s1 += options.getLDstep();
					s2 += options.getLDstep();

					if(s2 > run_end){
						s2 = run_end;
					}

					if(s2 - s1 < 1){
						break;
					}

					continue;
				}

				vector<double> mean;
				vector<vector<double> > variance;

				cout << "Pruning SNPs " << s1-run_start+1 << " to " << s2-run_start+1 << " of " << run_end - run_start+1 << "       \r";
				cout.flush();

				calcSetMeanVariance(pSNP,mean,variance);

				vector<bool> cur = vif_prune(variance, options.getLDthreshold());

				int k=0;
				for(int l=s1; l<=s2; l++){
					//update main list bu do not get back already excluded snps
					if(include.at(l) && !cur.at(k++)){
						include.at(l) = false;
					}
				}

				//advance window
				if(s2 == run_end){
					break;
				}

				s1 += options.getLDstep();
				s2 += options.getLDstep();

				if(s2 > run_end){
					s2 = run_end;
				}
				if(s2 - s1 < 1){
					break;
				}
			}// next window

			cout << "\n";

			// record what is in, what is out
			int cnt_in = 0, cnt_out = 0, cnt_dis = 0, cnt_din = 0, cnt_dout = 0;
			for(int l =run_start; l<= run_end; l++){
				if(!(*markers).at(l)->isEnabled()){
					cnt_dis++;
				}
			}
			for(int l = run_start; l<= run_end; l++){
				if(include.at(l)){
					//output to IN file
					LDIN << (*markers).at(l)->toString() << endl;
					cnt_in++;
					if((*markers).at(l)->isEnabled() && !(*markers).at(l)->isFlagged()){
						cnt_din++;
					}
				}
				else{
					//output to OUT file
					LDOUT << (*markers).at(l)->toString() << endl;
					cnt_out++;
					if((*markers).at(l)->isEnabled() && !(*markers).at(l)->isFlagged()){
						cnt_dout++;
					}
				}
			}

			opts::printLog("Overall results for chrom " + getString<int>(run_chr) + ", " + getString<int>(cnt_out) + " SNPs slated for removal, " + getString<int>(cnt_in) + " remaining\n");
			opts::printLog("For chrom " + getString<int>(run_chr) + ", " + getString<int>(cnt_dis) + " SNPs already removed or disabled, " + getString<int>(cnt_dout) + " more SNPs slated to be chopped out, " + getString<int>(cnt_din) + " would remain\n");
		}// next chrom

		if(options.doLDchop()){
			for(int i = 0; i < (int)include.size(); i++){
				if(!include.at(i) && (*markers).at(i)->isEnabled() && !(*markers).at(i)->isFlagged()){
					(*markers).at(i)->setEnabled(false);
					orig_num_markers++;
				}
			}
		}

		if(LDIN.is_open()){
			LDIN.close();
		}
		if(LDOUT.is_open()){
			LDOUT.close();
		}
	}

}

/*Borrowed and adapted from Plink*/
vector<double> LD::correlation2SNP(Marker* pri, Marker* sec){
	vector<double> results(2,-1);
	if(pri->isMicroSat() || sec->isMicroSat()){
		return results;
	}

	double X = 0; //sum count of primary minor allele
	double X2 = 0; // sum of pri squares
	double Y = 0; // sum count of secondary minor allele
	double Y2 = 0; // sec. MA squared sum
	double XY = 0; // multiplicative sum of minor allele counts
	double count = 0;

	int ploc = pri->getLoc();
	int sloc = sec->getLoc();

	for(int i = 0; i < (int)samples->size(); i++){
		Sample* samp = (*samples).at(i);
		if(!samp->isFounder()){
			continue;
		}

		bool p1 = samp->getAone(ploc);
		bool p2 = samp->getAtwo(ploc);
		bool pM = samp->getAmissing(ploc);
		if(p1 && p2 && pM) continue;

		bool q1 = samp->getAone(sloc);
		bool q2 = samp->getAtwo(sloc);
		bool qM = samp->getAmissing(sloc);
		if(q1 && q2 && qM) continue;

		count++;

		int sp = 0, ss = 0;
		if(pri->getChrom() == opts::_CHRX_ || sec->getChrom() == opts::_CHRX_){
			if(pri->getChrom() == opts::_CHRX_ && samp->getSex()){
				if(!p1){
					sp = 1;
				}
			}
			else{
				if(!p1){
					if(!p2){
						sp = 2;
					}
					else{
						sp = 1;
					}
				}
			}
			if(sec->getChrom() == opts::_CHRX_ && samp->getSex()){
				if(!q1){
					ss = 1;
				}
			}
			else{
				if(!q1){
					if(!q2){
						ss = 2;
					}
					else{
						ss = 1;
					}
				}
			}
		}
		else{
			if(!p1){
				if(!p2){
					sp = 2;
				}
				else{
					sp = 1;
				}
			}
			if(!q1){
				if(!q2){
					ss = 2;
				}
				else{
					ss = 1;
				}
			}
		}
		X += sp;
		Y += ss;
		XY += sp*ss;
		sp *= sp;
		ss *= ss;
		X2 += sp;
		Y2 += ss;
	}

	X /= count;
	X2 /= count;
	Y /= count;
	Y2 /= count;
	XY /= count;

	double var1 = X2 - X*X;
	double var2 = Y2 - Y*Y;
	double D = XY - X*Y;
	results.at(0) = (D*D)/(var1*var2);

	double dmax1 = (1-X)*Y;
	double dmax2 = (X)*(1-Y);

	double dmax = dmax1 < dmax2 ? dmax1 : dmax2;
	if(dmax == 0){
		results.at(1) = -1;
	}
	else{
		results.at(1) = D/dmax;
	}

	return results;
}

void LD::calcLDStatistics(){
	map<int, vector<int> > chroms;
	int prev_base = 0;
	int prev_chrom = -1;
	for(int m = 0; m < (int)markers->size(); m++){
		if((*markers).at(m)->getChrom() <= opts::_CHRX_ && Helpers::isValidMarker((*markers).at(m), &options, prev_base, prev_chrom)){
			if((*markers).at(m)->isEnabled() && !(*markers).at(m)->isFlagged()){
				chroms[(*markers).at(m)->getChrom()].push_back(m);
			}
		}
	}

	ofstream LD;
    string f = opts::_OUTPREFIX_ + "ld_calc" + options.getOut() + ".txt";
	if(!overwrite){
		f += "." + getString<int>(order);
	}
	LD.open(f.c_str(),ios::out);
	LD << "CHR1\tSNP1\tCHR2\tSNP2\tr^2\tD'\n";

	map<int, vector<int> >::iterator citer;
	cout << "Chrom size = " << chroms.size() << endl;
	for(citer = chroms.begin(); citer != chroms.end(); citer++){
		vector<int> snps = citer->second;
		for(int m = 0; m < (int)snps.size(); m++){
			Marker* primary = (*markers)[snps.at(m)];
			int s1 = m + 1;
			int s2 = m + options.getLDWin() - 1;
			if(s2 >= (int)snps.size()){
				s2 = snps.size() - 1;
			}

			while(s1 <= s2){
				Marker* secondary = (*markers).at(snps.at(s1));
				if((secondary->getBPLOC() - primary->getBPLOC()) > (options.getLDWinKB() * 1000)){
					break;
				}
				vector<double> results = correlation2SNP(primary, secondary);
				LD << primary->getChrom() << "\t" << primary->getProbeID() << "\t" << secondary->getChrom() << "\t" << secondary->getProbeID() << "\t" << results.at(0) << "\t" << results.at(1) << endl;
				s1++;
			}
		}
	}
	if(LD.is_open()){
		LD.close();
	}
}

}
