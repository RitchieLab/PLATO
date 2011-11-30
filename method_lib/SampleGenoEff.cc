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
#include <vector>
#include <list>
#include <map>
#include "SampleGenoEff.h"
#include "Sample.h"
#include "Family.h"
#include "Options.h"
#include "General.h"
#include "Helpers.h"
namespace Methods{
string SampleGenoEff::stepname = "sample-geno-eff";

void SampleGenoEff::calculate(Sample* samp){
	zeros_one = 0;
	total_one = 0;
		int msize = markers->size();
		if(samp->isEnabled()){
			int prev_base = 0;
			int prev_chrom = -1;
			for(int k = 0; k < msize; k++){
				if((*markers)[k]->isEnabled()){
		            if(options.doChrom()){
        		        if(!options.checkChrom((*markers)[k]->getChrom())){
	            	        continue;
		            	}
                		if(!options.checkBp((*markers)[k]->getBPLOC())){
	                    	continue;
		                }
            		}
                    if(options.doBpSpace()){
	                    if(prev_base == 0){
	                        prev_base = (*markers)[k]->getBPLOC();
	                        prev_chrom = (*markers)[k]->getChrom();
                        }
                        else{
	                        if((*markers)[k]->getChrom() == prev_chrom && (((*markers)[k]->getBPLOC() - prev_base) < options.getBpSpace())){
								(*markers)[k]->setFlag(true);
	                            continue;
                            }
                            prev_base = (*markers)[k]->getBPLOC();
                            prev_chrom = (*markers)[k]->getChrom();
                        }
                    }

					int loc = (*markers)[k]->getLoc();
					if(!(*markers)[k]->isMicroSat()){
						if(samp->getAone(loc) && samp->getAtwo(loc) && samp->getAmissing(loc)){
							zeros_one++;
						}
					}
					else{
						if(samp->getAbone(loc) == -1){
							zeros_one++;
						}
					}
					total_one++;
				}
			}
		}

}

void SampleGenoEff::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

	int msize = markers->size();
	int ssize = samples->size();
	zeros.resize(ssize);
	total.resize(ssize);
	if(opts::_ENZYMES_){
		enzyme_zeros.resize(ssize);
		enzyme_total.resize(ssize);
	}
	orig_num_samples = 0;

	for(int i = 0; i < ssize; i++){
		if((*samples)[i]->isEnabled()){
			int prev_base = 0;
			int prev_chrom = -1;
			for(int k = 0; k < msize; k++){
				if((*markers)[k]->isEnabled()){
		            if(options.doChrom()){
        		        if(!options.checkChrom((*markers)[k]->getChrom())){
	            	        continue;
		            	}
                		if(!options.checkBp((*markers)[k]->getBPLOC())){
	                    	continue;
		                }
            		}
                    if(options.doBpSpace()){
	                    if(prev_base == 0){
	                        prev_base = (*markers)[k]->getBPLOC();
	                        prev_chrom = (*markers)[k]->getChrom();
                        }
                        else{
	                        if((*markers)[k]->getChrom() == prev_chrom && (((*markers)[k]->getBPLOC() - prev_base) < options.getBpSpace())){
								(*markers)[k]->setFlag(true);
	                            continue;
                            }
                            prev_base = (*markers)[k]->getBPLOC();
                            prev_chrom = (*markers)[k]->getChrom();
                        }
                    }

					int loc = (*markers)[k]->getLoc();
					if(!(*markers)[k]->isMicroSat()){
						if((*samples)[i]->getAone(loc) && !(*samples)[i]->getAtwo(loc)){
							zeros[i]++;
							if(opts::_ENZYMES_){
								enzyme_zeros[i][(*markers)[k]->getEnzyme()]++;
							}
						}
					}
					else{
						if((*samples)[i]->getAbone(loc) == -1){
							zeros[i]++;
							if(opts::_ENZYMES_){
								enzyme_zeros[i][(*markers)[k]->getEnzyme()]++;
							}
						}
					}
					total[i]++;
					if(opts::_ENZYMES_){
						enzyme_total[i][(*markers)[k]->getEnzyme()]++;
					}
				}
			}
		}
	}

}

void SampleGenoEff::PrintSummary(){

	string fname1 = opts::_OUTPREFIX_ + "sample_geno_eff" + options.getOut() + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	ofstream indeff (fname1.c_str());
	if(!indeff.is_open()){
		opts::printLog("Unable to open " + fname1 + "\n");
		throw MethodException("Unable to open " + fname1 + "\n");
	}
	opts::addFile("Sample",stepname, fname1);
	string sdetails = "";
	if(opts::_SAMPDESC_.length() > 0){
		sdetails = (*samples)[0]->getDetailHeaders();
	}
	opts::addHeader(fname1, "%GenoEff_All");
	indeff << "FamID\t"
		   << "IndID\tCenter\tSex\t"
		   << "Affection_Status\t"
		   << "Plate\tWell\t%GenoEff_All";
	vector<string> enzymes;
	if(opts::_ENZYMES_){
		int msize = markers->size();
		for(int i = 0; i < msize; i++){
			if((*markers)[i]->isEnabled() && !(*markers)[i]->isFlagged()){
	            if(options.doChrom()){
  		            if(!options.checkChrom((*markers)[i]->getChrom())){
			            continue;
		            }
		            if(!options.checkBp((*markers)[i]->getBPLOC())){
			            continue;
		            }
		        }

				vector<string>::iterator e_iter = find(enzymes.begin(), enzymes.end(), (*markers)[i]->getEnzyme());
				if(e_iter == enzymes.end()){
					enzymes.push_back((*markers)[i]->getEnzyme());
				}
			}
		}
		if(enzymes.size() > 1){
			sort(enzymes.begin(), enzymes.end());
		}
		for(int i = 0; i < (int)enzymes.size(); i++){
			indeff << "\tGenoEff_" << enzymes[i];
			opts::addHeader(fname1, "GenoEff_" + enzymes[i]);
		}
	}
	if(opts::_SAMPDESC_.length() > 0){
		indeff << "\t" << sdetails;
	}
	indeff << endl;

	indeff.precision(4);

	int ssize = samples->size();

	for(int i = 0; i < ssize; i++){
		if((*samples)[i]->isEnabled()){
			float percent = 0.0f;

			if(total[i] > 0){
				percent = (1.0f - ((float)zeros[i]/(float)total[i])) * 100.0f;
			}

			indeff << (*samples)[i]->getFamID() << "\t"
				   << (*samples)[i]->getInd() << "\t"
				   << (*samples)[i]->getFamily()->getCenter() << "\t";
			if((*samples)[i]->getSex()){
				indeff << "M\t";
			}
			else{
				indeff << "F\t";
			}
			indeff << (*samples)[i]->getPheno() << "\t";
			indeff << (*samples)[i]->getPlate() << "\t"
				   << (*samples)[i]->getWell() << "\t"
				   << percent;
			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < (int)enzymes.size(); e++){
					float epercent = 0.0f;
					if(enzyme_total[i][enzymes[e]] > 0){
						epercent = (1.0f - ((float)enzyme_zeros[i][enzymes[e]] / (float)enzyme_total[i][enzymes[e]])) * 100.0f;
					}
					indeff << "\t" << epercent;
				}
			}
			if(opts::_SAMPDESC_.length() > 0){
				indeff << "\t" << (*samples)[i]->getDetails();
			}
			indeff << endl;
		}
	}

	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}

	if(indeff.is_open()){
		indeff.close();
	}

}

void SampleGenoEff::filter(){
	if(options.doThreshSamplesLow() || options.doThreshSamplesHigh()){
		int ssize = samples->size();

		for(int i = 0; i < ssize; i++){
			if((*samples)[i]->isEnabled()){
				float percent = 0.0f;
				bool inc = false;
				if(total[i] > 0){
					percent = (1.0f - ((float)zeros[i]/(float)total[i])) * 100.0f;
				}

				if(options.doThreshSamplesLow() && Helpers::dLess(percent, options.getThreshSamplesLow())){
					(*samples)[i]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshSamplesHigh() && Helpers::dGreater(percent, options.getThreshSamplesHigh())){
					(*samples)[i]->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_samples++;
				}
			}
		}
	}
}

void SampleGenoEff::FilterSummary(){
	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Samples Passed:\t" + getString<int>(opts::_SAMPLES_WORKING_ - orig_num_samples) + " (" +
		getString<float>(((float) (opts::_SAMPLES_WORKING_ - orig_num_samples) / (float) opts::_SAMPLES_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_SAMPLES_WORKING_) + "\n");
	opts::_SAMPLES_WORKING_ -= orig_num_samples;

}

}
