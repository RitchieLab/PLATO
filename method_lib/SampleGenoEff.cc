#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
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
				if((*markers).at(k)->isEnabled()){
		            if(options.doChrom()){
        		        if(!options.checkChrom((*markers).at(k)->getChrom())){
	            	        continue;
		            	}
                		if(!options.checkBp((*markers).at(k)->getBPLOC())){
	                    	continue;
		                }
            		}
                    if(options.doBpSpace()){
	                    if(prev_base == 0){
	                        prev_base = (*markers).at(k)->getBPLOC();
	                        prev_chrom = (*markers).at(k)->getChrom();
                        }
                        else{
	                        if((*markers).at(k)->getChrom() == prev_chrom && (((*markers).at(k)->getBPLOC() - prev_base) < options.getBpSpace())){
								(*markers).at(k)->setFlag(true);
	                            continue;
                            }
                            prev_base = (*markers).at(k)->getBPLOC();
                            prev_chrom = (*markers).at(k)->getChrom();
                        }
                    }

					int loc = (*markers).at(k)->getLoc();
					if(!(*markers).at(k)->isMicroSat()){
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
		if((*samples).at(i)->isEnabled()){
			int prev_base = 0;
			int prev_chrom = -1;
			for(int k = 0; k < msize; k++){
				if((*markers).at(k)->isEnabled()){
		            if(options.doChrom()){
        		        if(!options.checkChrom((*markers).at(k)->getChrom())){
	            	        continue;
		            	}
                		if(!options.checkBp((*markers).at(k)->getBPLOC())){
	                    	continue;
		                }
            		}
                    if(options.doBpSpace()){
	                    if(prev_base == 0){
	                        prev_base = (*markers).at(k)->getBPLOC();
	                        prev_chrom = (*markers).at(k)->getChrom();
                        }
                        else{
	                        if((*markers).at(k)->getChrom() == prev_chrom && (((*markers).at(k)->getBPLOC() - prev_base) < options.getBpSpace())){
								(*markers).at(k)->setFlag(true);
	                            continue;
                            }
                            prev_base = (*markers).at(k)->getBPLOC();
                            prev_chrom = (*markers).at(k)->getChrom();
                        }
                    }

					int loc = (*markers).at(k)->getLoc();
					if(!(*markers).at(k)->isMicroSat()){
						if((*samples).at(i)->getAone(loc) && !(*samples).at(i)->getAtwo(loc)){
							zeros.at(i)++;
							if(opts::_ENZYMES_){
								enzyme_zeros.at(i)[(*markers).at(k)->getEnzyme()]++;
							}
						}
					}
					else{
						if((*samples).at(i)->getAbone(loc) == -1){
							zeros.at(i)++;
							if(opts::_ENZYMES_){
								enzyme_zeros.at(i)[(*markers).at(k)->getEnzyme()]++;
							}
						}
					}
					total.at(i)++;
					if(opts::_ENZYMES_){
						enzyme_total.at(i)[(*markers).at(k)->getEnzyme()]++;
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
		sdetails = (*samples).at(0)->getDetailHeaders();
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
			if((*markers).at(i)->isEnabled() && !(*markers).at(i)->isFlagged()){
	            if(options.doChrom()){
  		            if(!options.checkChrom((*markers).at(i)->getChrom())){
			            continue;
		            }
		            if(!options.checkBp((*markers).at(i)->getBPLOC())){
			            continue;
		            }
		        }

				vector<string>::iterator e_iter = find(enzymes.begin(), enzymes.end(), (*markers).at(i)->getEnzyme());
				if(e_iter == enzymes.end()){
					enzymes.push_back((*markers).at(i)->getEnzyme());
				}
			}
		}
		if(enzymes.size() > 1){
			sort(enzymes.begin(), enzymes.end());
		}
		for(int i = 0; i < (int)enzymes.size(); i++){
			indeff << "\tGenoEff_" << enzymes.at(i);
			opts::addHeader(fname1, "GenoEff_" + enzymes.at(i));
		}
	}
	if(opts::_SAMPDESC_.length() > 0){
		indeff << "\t" << sdetails;
	}
	indeff << endl;

	indeff.precision(4);

	int ssize = samples->size();

	for(int i = 0; i < ssize; i++){
		if((*samples).at(i)->isEnabled()){
			float percent = 0.0f;

			if(total.at(i) > 0){
				percent = (1.0f - ((float)zeros.at(i)/(float)total.at(i))) * 100.0f;
			}

			indeff << (*samples).at(i)->getFamID() << "\t"
				   << (*samples).at(i)->getInd() << "\t"
				   << (*samples).at(i)->getFamily()->getCenter() << "\t";
			if((*samples).at(i)->getSex()){
				indeff << "M\t";
			}
			else{
				indeff << "F\t";
			}
			indeff << (*samples).at(i)->getPheno() << "\t";
			indeff << (*samples).at(i)->getPlate() << "\t"
				   << (*samples).at(i)->getWell() << "\t"
				   << percent;
			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < (int)enzymes.size(); e++){
					float epercent = 0.0f;
					if(enzyme_total.at(i).at(enzymes.at(e)) > 0){
						epercent = (1.0f - ((float)enzyme_zeros.at(i).at(enzymes.at(e)) / (float)enzyme_total.at(i).at(enzymes.at(e)))) * 100.0f;
					}
					indeff << "\t" << epercent;
				}
			}
			if(opts::_SAMPDESC_.length() > 0){
				indeff << "\t" << (*samples).at(i)->getDetails();
			}
			indeff << endl;
		}
	}

	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		(*markers).at(i)->setFlag(false);
	}

	if(indeff.is_open()){
		indeff.close();
	}

}

void SampleGenoEff::filter(){
	if(options.doThreshSamplesLow() || options.doThreshSamplesHigh()){
		int ssize = samples->size();

		for(int i = 0; i < ssize; i++){
			if((*samples).at(i)->isEnabled()){
				float percent = 0.0f;
				bool inc = false;
				if(total.at(i) > 0){
					percent = (1.0f - ((float)zeros.at(i)/(float)total.at(i))) * 100.0f;
				}

				if(options.doThreshSamplesLow() && Helpers::dLess(percent, options.getThreshSamplesLow())){
					(*samples).at(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshSamplesHigh() && Helpers::dGreater(percent, options.getThreshSamplesHigh())){
					(*samples).at(i)->setEnabled(false);
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
