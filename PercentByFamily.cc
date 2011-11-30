/**********************************************************************************
*                       Family Genotype Efficiency Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all genotypes and generates a genotype efficiency for all families and markers
*
*
*
*File: PercentByFamily.cc
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
#include "PercentByFamily.h"
#include "Chrom.h"
#include <General.h>
#include <Helpers.h>

string PercentByFamily::stepname = "family-geno-eff";

/*
 *Function: FilterSummary
 *Description:
 *Outputs total number of markers and families remaining
 *
 */
void PercentByFamily::FilterSummary(){
	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

/*
 * Function: PrintSummary
 * Description:
 * Outputs results.
 */
void PercentByFamily::PrintSummary(){
	int msize = markers->size();
	int fsize = families->size();

	string fname1 = opts::_OUTPREFIX_ + "family_geno_eff" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	string fname2 = opts::_OUTPREFIX_ + "family_geno_eff_marker" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
		fname2 += "." + getString<int>(order);
	}
	ofstream byfamily (fname1.c_str());
	ofstream bymarker (fname2.c_str());
	if(!bymarker.is_open()){
		opts::printLog("Unable to open " + fname1 + "\n");
		exit(1);
	}
	opts::addFile("Family",stepname,fname1);
	if(!byfamily.is_open()){
		opts::printLog("Unable to open " + fname2 + "\n");
		exit(1);
	}
	opts::addFile("Marker",stepname,fname2);
	bymarker.precision(4);
	bymarker << "Chrom\trsID\tProbeID\tbploc";
	if((*markers)[0]->getDetailHeaders().size() > 0){
		bymarker << "\t" << (*markers)[0]->getDetailHeaders();
	}
	bymarker << "\t%GenoEff_Fam_All\tFamily_Zero_Count\tTotal_Families_Used" << endl;
	opts::addHeader(fname2, "%GenoEff_Fam_All");
	opts::addHeader(fname2, "Family_Zero_Count");
	opts::addHeader(fname2, "Total_Families_Used");

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
			float percent = 0.0f;
			if(mtotal[i] > 0){
				percent = (1.0f - ((float)mzeros[i]/(float)mtotal[i]));// * 100.0f;
			}
			bymarker << (*markers)[i]->getChrom() << "\t"
					 << (*markers)[i]->getRSID() << "\t"
					 << (*markers)[i]->getProbeID() << "\t"
					 << (*markers)[i]->getBPLOC() << "\t"
					 << percent << "\t"
					 << mzeros[i] << "\t"
					 << mtotal[i]
					 << endl;
		}
	}

	if(!byfamily.is_open()){
		opts::printLog("Unable to open efficiency_family.txt\n");
		exit(0);
	}
	byfamily.precision(4);
	byfamily << "FamID\tNumInds\tCenter\t%GenoEff_All";
	opts::addHeader(fname1, "%GenoEff_All");

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
            byfamily << "\tGenoEff_" << enzymes[i];
			opts::addHeader(fname1, enzymes[i]);
        }
    }
	byfamily << endl;

	for(int i = 0; i < fsize; i++){
		if((*families)[i]->isEnabled()){
			float percent = 0.0f;
			if(ftotal[i] > 0){
				percent = (1.0f - ((float)fzeros[i]/(float)ftotal[i]));// * 100.0f;
			}

			byfamily << (*families)[i]->getFamID() << "\t"
					 << (*families)[i]->getSamples()->size() << "\t"
					 << (*families)[i]->getCenter() << "\t"
					 << percent;
					 //<< fzeros[i] << "\t"
					 //<< ftotal[i];
			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < (int)enzymes.size(); e++){
					float eper = 0.0f;
					if(enzyme_zeros[i][enzymes[e]] > 0){
						eper = (1.0f - ((float)enzyme_zeros[i][enzymes[e]]/(float)enzyme_total[i][enzymes[e]]));// * 100.0f;
					}
					byfamily << "\t" << eper;
				}
			}
			byfamily << endl;
		}
	}

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}

}


/*
 * Function: filter
 * Description:
 * Performs filtering of markers based on threshold
 */
void PercentByFamily::filter(){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
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

				bool inc = false;
				float percent = 0.0f;
				if(mtotal[i] > 0){
					percent = (1.0f - ((float)mzeros[i]/(float)mtotal[i]));// * 100.0f;
				}
				if(options.doThreshMarkersLow() && Helpers::dLess(percent,options.getThreshMarkersLow())){
					(*markers)[i]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && Helpers::dGreater(percent, options.getThreshMarkersHigh())){
					(*markers)[i]->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}
		}
	}
}

/*
 * Function: process
 * Description:
 * Main function for processing
 */
void PercentByFamily::process(DataSet* ds){
	data_set = ds;
	markers = data_set->get_markers();
	families = data_set->get_families();
	samples = data_set->get_samples();
	marker_map =data_set->get_marker_map();

	int msize = markers->size();
/*	for(int i = 0; i < msize; i++){
		if((*markers)[i]->isEnabled()){
			orig_num_markers++;
		}
	}
*/
	int fsize = families->size();
/*	for(int i = 0; i < fsize; i++){
		if((*families)[i]->isEnabled()){
			orig_num_families++;
		}
	}
*/
	fzeros.resize(fsize);
	ftotal.resize(fsize);
	mzeros.resize(msize);
	mtotal.resize(msize);
	if(opts::_ENZYMES_){
		enzyme_zeros.resize(fsize);
		enzyme_total.resize(fsize);
	}
	for(int i = 0; i < fsize; i++){
		int fsamps = 0;
		vector<Sample*>* samps = (*families)[i]->getSamples();
		for(int f = 0; f < (int)samps->size(); f++){
			Sample* samp = (*samps)[f];
			if(samp->isEnabled()){
				fsamps++;
			}
		}
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
				int tempzeros = 0;
				int inds = 0;
				int numinds = 0;
				for(int j = 0; j< fsamps; j++){
					if((*samps)[j]->isEnabled()){
						if(!(*markers)[k]->isMicroSat()){
							if((*samps)[j]->getAone(loc) && (*samps)[j]->getAtwo(loc) && (*samps)[j]->getAmissing(loc)){
								tempzeros++;
							}
							else{
								inds++;
							}
						}
						else{
							if((*samps)[j]->getAbone(loc) == -1 && (*samps)[j]->getAmissing(loc)){
								tempzeros++;
							}
							else{
								inds++;
							}
						}
						numinds++;
					}
				}
				if(inds != fsamps || tempzeros > 0){
					fzeros[i]++;
					if(opts::_ENZYMES_ && (*markers)[k]->getEnzyme().length() > 0){
						enzyme_zeros[i][(*markers)[k]->getEnzyme()]++;
					}
					//if(inds != 3){
						mzeros[k]++;
					//}
				}
				ftotal[i]++;
				if(opts::_ENZYMES_ && (*markers)[k]->getEnzyme().length() > 0){
					enzyme_total[i][(*markers)[k]->getEnzyme()]++;
				}
				//if(numinds == 3){
					mtotal[k]++;
				//}
				tempzeros = 0;
				inds = 0;
				numinds = 0;
			}
		}
	}
}

