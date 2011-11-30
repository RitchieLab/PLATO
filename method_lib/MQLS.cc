#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <unistd.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#include <fenv.h>
#include <algorithm>
#include "MQLS.h"
#include "General.h"
#include "Helper.h"
#include "cdflib.h"

using namespace std;
namespace Methods{
string MQLS::stepname = "mqls";

void MQLS::PrintSummary(){
	string fname1 = opts::_OUTPREFIX_ + "tdt" + options.getOut() + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	ofstream output (fname1.c_str(), ios::out);
	if(!output){
		opts::printLog("Error opening " + fname1 + ". Exiting!\n");
		//exit(1);
		throw MethodException("Error opening " + fname1 + ". Exiting!\n");
	}
	opts::addFile("Marker", stepname, fname1);
	output.precision(4);
	output << "Chrom\trsID\tProbeID\tbploc";
	if((*markers).at(0)->getDetailHeaders().size() > 0){
		output << "\t" << (*markers).at(0)->getDetailHeaders();
	}
	output << "\tChi_square\tTDT_pvalue\tTDT_neglog(pvalue)\tNum_Fams\tT:U\tA1:A2\tOR\tL" + getString<double>(options.getCI()*100) + "\t" + "U" + getString<double>(options.getCI()*100) << endl;
	opts::addHeader(fname1, "Chi_square");
	opts::addHeader(fname1, "TDT_pvalue");
	opts::addHeader(fname1, "TDT_neglog(pvalue)");
	opts::addHeader(fname1, "Num_Fams");
	opts::addHeader(fname1, "T:U");
	opts::addHeader(fname1, "A1:A2");
	opts::addHeader(fname1, "OR");
	opts::addHeader(fname1, "L" + getString<double>(options.getCI()*100));
	opts::addHeader(fname1, "U" + getString<double>(options.getCI()*100));

	int msize = markers->size();

	double zt = ltqnorm(1 - (1 - options.getCI()) / 2);

	for(int i = 0; i < msize; i++){
		if((*markers).at(i)->isEnabled() && !(*markers).at(i)->isFlagged()){
            if(options.doChrom()){
                if(!options.checkChrom((*markers).at(i)->getChrom())){
                    continue;
                }
                if(!options.checkBp((*markers).at(i)->getBPLOC())){
                    continue;
                }
            }
			output << (*markers).at(i)->toString() << "\t";
			output << chi.at(i) << "\t"
				<< pval.at(i) << "\t"
				<< (double)abs(log10(pval.at(i))) << "\t";
			output << fams_used.at(i) << "\t"
				<< trans.at(i) << ":" << untrans.at(i) << "\t";
			if(!(*markers).at(i)->isMicroSat()){
				output << (*markers).at(i)->getAllele1() << ":" << (*markers).at(i)->getAllele2();
			}
			else{
				output << "NA";
			}
			double OR = (double)trans.at(i) / (double)untrans.at(i);
			output << "\t" << OR;
			double OR_lower = exp(log(OR) - zt * sqrt(1/trans.at(i) + 1/untrans.at(i)));
			double OR_upper = exp(log(OR) + zt * sqrt(1/trans.at(i) + 1/untrans.at(i)));
			output << "\t" << OR_lower << "\t" << OR_upper;

				output << endl;
		}
		(*markers).at(i)->setFlag(false);
	}

	if(output.is_open()){
		output.close();
	}
}

void MQLS::filter(){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int size = markers->size();
		for(int i = 0; i < size; i++){
			if((*markers).at(i)->isEnabled() && !(*markers).at(i)->isFlagged()){
	            if(options.doChrom()){
	                if(!options.checkChrom((*markers).at(i)->getChrom())){
	                    continue;
		            }
		            if(!options.checkBp((*markers).at(i)->getBPLOC())){
			            continue;
		            }
		        }

				bool inc = false;
				if(options.doThreshMarkersHigh() && dGreater(pval.at(i), options.getThreshMarkersHigh())){
					(*markers).at(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersLow() && dLess(pval.at(i),options.getThreshMarkersLow())){
					(*markers).at(i)->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}
		}
	}
}

void MQLS::setThreshold(string thresh){
	options.setUp(thresh);
}

void MQLS::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
	    getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
    opts::_MARKERS_WORKING_ -= orig_num_markers;
}

void MQLS::calculate(Marker* mark){
	vector<vector<double> >







	chi_one = 0;
	pval_one = 0;
	fams_used_one = 0;
	trans_one = 0;
	untrans_one = 0;
	ci_l_one = 0;
	ci_u_one = 0;

	if(!mark->isMicroSat()){
		int mloc = mark->getLoc();

		int ssize = samples->size();
		double zt = ltqnorm(1 - (1 - options.getCI()) / 2);
		//start new
		int fams = 0;
		double trans1 = 0;
		double trans2 = 0;
		for(int i = 0; i < ssize; i++){
			Sample* samp = (*samples).at(i);
			if(samp != NULL && samp->getDad() != NULL && samp->getMom() != NULL && samp->getDadID() != "0" && samp->getMomID() != "0" && samp->isEnabled() && samp->getDad()->isEnabled() && samp->getMom()->isEnabled()){
				Sample* dad = samp->getDad();
				Sample* mom = samp->getMom();

				int trA = 0; //transmitted allele from first het parent
				int unA = 0; //untransmitted allele from first het parent
				int trB = 0; //transmitted allele from second het parent
				int unB = 0; //untransmitted allele from second het parent

				if(!mark->isMicroSat()){
					bool pat1 = dad->getAone(mloc);
					bool pat2 = dad->getAtwo(mloc);
					bool pat3 = dad->getAmissing(mloc);
					bool mat1 = mom->getAone(mloc);
					bool mat2 = mom->getAtwo(mloc);
					bool mat3 = mom->getAmissing(mloc);

					if(pat1 == pat2 && mat1 == mat2){ //mono alleleic
						continue;
					}
					if((pat1 && pat2 && pat3) || (mat1 && mat2 && mat3)){ //no genotype
						continue;
					}

					if(samp->getPheno() != 2){//make sure child is affected
						continue;
					}

					bool kid1 = samp->getAone(mloc);
					bool kid2 = samp->getAtwo(mloc);
					bool kid3 = samp->getAmissing(mloc);

					//kid is 0/0 call
					if(kid1 && kid2 && kid3){
						continue;
					}

					//00 - homozygous allele 1
					if((!kid1) && (!kid2)){
						if(((!pat1) && pat2) && ((!mat1) && mat2)){
							trA = 1;
							unA = 2;
							trB = 1;
							unB = 2;
						}
						else{
							trA = 1;
							unA = 2;
						}
					}
					else if((!kid1) && kid2){ //01 - heterozygous
						//het dad
						if(pat1 != pat2){
							//het mom
							if(mat1 != mat2){
								trA = 1;
								trB = 2;
								unA = 2;
								unB = 1;
							}
							else if(!mat1){
								trA = 2;
								unA = 1;
							}
							else{
								trA = 1;
								unA = 2;
							}
						}
						else if(!pat1){
							trA = 2;
							unA = 1;
						}
						else{
							trA = 1;
							unA = 2;
						}
					}
					else{ //11 - homozygous allele 2
						if(((!pat1) && pat2) && ((!mat1) && mat2)){ //dad het & mom het
							trA = 2;
							unA = 1;
							trB = 2;
							unB = 1;
						}
						else{
							trA = 2;
							unA = 1;
						}
					}
					//increment transmission counts
					if(trA == 1) trans1++;
					if(trB == 1) trans1++;
					if(trA == 2) trans2++;
					if(trB == 2) trans2++;
					fams++;
				}//end !micro-sat
			}
		}//end sample iteration

		//skipped discordant parent

		double tdt_chisq, par_chisq, com_chisq;
		tdt_chisq = par_chisq = com_chisq = -1;

		if(trans1 + trans2 > 0){
			tdt_chisq = ((trans1 - trans2)*(trans1 - trans2))/(trans1 + trans2);
		}
		double pvalue, p, bound, df = 1;
		pvalue = -1;
		int code = 1, status;
		if(tdt_chisq > -1){
			pvalue = p_from_chi(tdt_chisq, df);
		}

		pval_one = pvalue;
		chi_one = tdt_chisq;
		fams_used_one = fams;
		trans_one = trans1;
		untrans_one = trans2;
		odds_ratio_one = (double)trans_one / (double)untrans_one;
		ci_l_one = exp(log(odds_ratio_one) - zt * sqrt(1/trans_one + 1/untrans_one));
		ci_u_one = exp(log(odds_ratio_one) + zt * sqrt(1/trans_one + 1/untrans_one));
	}
	else{
		//calculate pval for micro satellites
		//compares each allele against all others
		//2x2  A = allele, Oth = other alleles
		//-------------------
		//|        |        |
		//| A/A    |  A/Oth |
		//|        |        |
		//-------------------
		//|        |        |
		//| Oth/A  | Oth/Oth|
		//|        |        |
		//-------------------
		//
		//Finds Allele/Other combination with lowest pvalue and chooses that one to output.
		int numalleles = mark->getNumAlleles();
		vector<double> results;
		results.resize(numalleles, 0);
		for(int n = 0; n < numalleles; n++){
		}

	}
}

void MQLS::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;


	int msize = markers->size();
	chi.resize(msize);
	pval.resize(msize);
	fams_used.resize(msize);
	maf.resize(msize);
	trans.resize(msize);
	untrans.resize(msize);
	int ssize = samples->size();

	map<int, string> people;
	int prev_base = 0;
	int prev_chrom = -1;
	for(int m = 0; m < msize; m++){
		if((*markers).at(m)->isEnabled()){
            if(options.doChrom()){
                if(!options.checkChrom((*markers).at(m)->getChrom())){
                    continue;
                }
                if(!options.checkBp((*markers).at(m)->getBPLOC())){
                    continue;
                }
            }
            if(options.doBpSpace()){
	            if(prev_base == 0){
	                prev_base = (*markers).at(m)->getBPLOC();
                    prev_chrom = (*markers).at(m)->getChrom();
                }
                else{
 	               if((*markers).at(m)->getChrom() == prev_chrom && (((*markers).at(m)->getBPLOC() - prev_base) < options.getBpSpace())){
	               	   (*markers).at(m)->setFlag(true);
				   	   continue;
                   }
                   prev_base = (*markers).at(m)->getBPLOC();
                   prev_chrom = (*markers).at(m)->getChrom();
                }
            }
			if(!(*markers).at(m)->isMicroSat()){
				int mloc = (*markers).at(m)->getLoc();


				//start new
				int fams = 0;
				double trans1 = 0;
				double trans2 = 0;
				for(int i = 0; i < ssize; i++){
					Sample* samp = (*samples).at(i);
					if(samp != NULL && samp->getDad() != NULL && samp->getMom() != NULL && samp->getDadID() != "0" && samp->getMomID() != "0" && samp->isEnabled() && samp->getDad()->isEnabled() && samp->getMom()->isEnabled()){
						Sample* dad = samp->getDad();
						Sample* mom = samp->getMom();

						int trA = 0; //transmitted allele from first het parent
						int unA = 0; //untransmitted allele from first het parent
						int trB = 0; //transmitted allele from second het parent
						int unB = 0; //untransmitted allele from second het parent

						if(!(*markers).at(m)->isMicroSat()){
							bool pat1 = dad->getAone(mloc);
							bool pat2 = dad->getAtwo(mloc);
							bool mat1 = mom->getAone(mloc);
							bool mat2 = mom->getAtwo(mloc);

							if(pat1 == pat2 && mat1 == mat2){ //mono alleleic
								continue;
							}
							if((pat1 && !pat2) || (mat1 && !mat2)){ //no genotype
								continue;
							}

							if(samp->getPheno() != 2){//make sure child is affected
								continue;
							}

							bool kid1 = samp->getAone(mloc);
							bool kid2 = samp->getAtwo(mloc);

							//kid is 0/0 call
							if(kid1 && !kid2){
								continue;
							}

							//00 - homozygous allele 1
							if((!kid1) && (!kid2)){
								if(((!pat1) && pat2) && ((!mat1) && mat2)){
									trA = 1;
									unA = 2;
									trB = 1;
									unB = 2;
								}
								else{
									trA = 1;
									unA = 2;
								}
							}
							else if((!kid1) && kid2){ //01 - heterozygous
								//het dad
								if(pat1 != pat2){
									//het mom
									if(mat1 != mat2){
										trA = 1;
										trB = 2;
										unA = 2;
										unB = 1;
									}
									else if(!mat1){
										trA = 2;
										unA = 1;
									}
									else{
										trA = 1;
										unA = 2;
									}
								}
								else if(!pat1){
									trA = 2;
									unA = 1;
								}
								else{
									trA = 1;
									unA = 2;
								}
							}
							else{ //11 - homozygous allele 2
								if(((!pat1) && pat2) && ((!mat1) && mat2)){ //dad het & mom het
									trA = 2;
									unA = 1;
									trB = 2;
									unB = 1;
								}
								else{
									trA = 2;
									unA = 1;
								}
							}
							//increment transmission counts
							if(trA == 1) trans1++;
							if(trB == 1) trans1++;
							if(trA == 2) trans2++;
							if(trB == 2) trans2++;
							fams++;
						}//end !micro-sat
					}
				}//end sample iteration

				//skipped discordant parent

				double tdt_chisq, par_chisq, com_chisq;
				tdt_chisq = par_chisq = com_chisq = -1;

				if(trans1 + trans2 > 0){
					tdt_chisq = ((trans1 - trans2)*(trans1 - trans2))/(trans1 + trans2);
				}
				double pvalue, p, bound, df = 1;
				pvalue = -1;
				int code = 1, status;
				if(tdt_chisq > -1){
					pvalue = p_from_chi(tdt_chisq, df);
				}

				pval.at(m) = pvalue;
				chi.at(m) = tdt_chisq;
				fams_used.at(m) = fams;
				maf.at(m) = -1;
				trans.at(m) = trans1;
				untrans.at(m) = trans2;
			}
			else{
				//calculate pval for micro satellites
				//compares each allele against all others
				//2x2  A = allele, Oth = other alleles
				//-------------------
				//|        |        |
				//| A/A    |  A/Oth |
				//|        |        |
				//-------------------
				//|        |        |
				//| Oth/A  | Oth/Oth|
				//|        |        |
				//-------------------
				//
				//Finds Allele/Other combination with lowest pvalue and chooses that one to output.
				Marker* mark = (*markers).at(m);
				int numalleles = mark->getNumAlleles();
				vector<double> results;
				results.resize(numalleles, 0);
				for(int n = 0; n < numalleles; n++){
				}

			}
		}
	}
}
}
