/**********************************************************************************
*                       BEAGLE Output Module
*
* Written by: Justin Giles
*             Vanderbilt University
*             Center for Human Genetics Research
*
* Outputs BEAGLE input files.
*
*
*File: BEAGLEOutput.cc
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
#include <iomanip>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include "BEAGLEOutput.h"
#include "General.h"
#include "Helper.h"

namespace Methods{
/*
 *Function: FilterSummary
 *Description:
 *Not used.
 */
void BEAGLEOutput::FilterSummary(){
}

/*
 * Function: PrintSummary
 * Description:
 * Resets marker flags
 */
void BEAGLEOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}

}

/*
 * Function: filter
 * Description:
 * not used
 */
void BEAGLEOutput::filter(){
}

/*
 * Function: process
 * Description:
 * Main function for producing output files.
 * Creates output files, one per chromosome.
 */
void BEAGLEOutput::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

   	int ssize = samples->size();
	int msize = markers->size();
	int fsize = families->size();

	vector<int> chrom_counts;
	chrom_counts.resize(24,0);
	if(opts::_DOG_){
		chrom_counts.resize(40, 0);
	}
	int prev_base = 0;
	int prev_chrom = -1;
	for(int i = 0; i < msize; i++){
		Marker* mark = (*markers)[i];
		if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
/*			if(options.doChrom()){
				if(!options.checkChrom(mark->getChrom())){
					mark->setFlag(true);
				    continue;
			    }
			    if(!options.checkBp(mark->getBPLOC())){
					mark->setFlag(true);
				    continue;
			    }
			}
            if(options.doBpSpace()){
	            if(prev_base == 0){
	                prev_base = mark->getBPLOC();
                    prev_chrom = mark->getChrom();
                }
                else{
	                if(mark->getChrom() == prev_chrom && ((mark->getBPLOC() - prev_base) < options.getBpSpace())){
		    			mark->setFlag(true);
						continue;
	                }
	                prev_base = mark->getBPLOC();
	                prev_chrom = mark->getChrom();
                }
			}
*/
			chrom_counts[mark->getChrom() - 1]++;
		}
	}


	int parents = 0;
	int stotal = 0;
	for(int i = 0; i < ssize; i++){
		Sample* samp = (*samples)[i];
		if((samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())) && !samp->isFlagged()){
			if(samp->getDad() != NULL && samp->getMom() != NULL && (samp->getDad()->isEnabled() || (samp->getDad()->isExcluded() && options.doIncExcludedSamples()) || (!samp->getDad()->isEnabled() && options.doIncDisabledSamples())) && (samp->getMom()->isEnabled() || (samp->getMom()->isExcluded() && options.doIncExcludedSamples()) || (!samp->getMom()->isEnabled() && options.doIncDisabledSamples())) && !samp->getDad()->isFlagged() && !samp->getMom()->isFlagged() && samp->getFamily()->getTotalInds() == 3){
				samp->getDad()->setFlag(true);
				samp->getMom()->setFlag(true);
				parents += 2;
			}
		}
		if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
			stotal++;
		}
	}

	string marker_types = "";

	for(int c = 0; c < (int)chrom_counts.size(); c++){
		if(chrom_counts[c] > 0){
			if(options.getChrom() == -1 || (options.getChrom() == (c + 1))){
				string fname1 = opts::_OUTPREFIX_ + "input_beagle_chr" + getString<int>(c + 1) + options.getOut() + ".txt";//"_" + getString<int>(order) + ".txt";
				if(options.getOverrideOut().size() > 0){
					fname1 = options.getOverrideOut() + "_chr" + getString<int>(c + 1) + ".txt";
				}
				if(!overwrite){
					fname1 += "." + getString<int>(order);
				}
				filenames.push_back(fname1);
				ofstream str (fname1.c_str());
				if(!str.is_open()){
					opts::printLog("Unable to open " + fname1 + " for output!\n");
					exit(1);
				}
				string fname2 = opts::_OUTPREFIX_ + "input_beagle_chr" + getString<int>(c + 1) + "_trait" + options.getOut() + ".txt";//_" + getString<int>(order) + ".txt";
				if(options.getOverrideOut().size() > 0){
					fname2 = options.getOverrideOut() + "_chr" + getString<int>(c + 1) + "_trait.txt";
				}
				if(!overwrite){
					fname2 += "." + getString<int>(order);
				}
				filenames.push_back(fname2);
				ofstream trait (fname2.c_str());
				if(!trait.is_open()){
					opts::printLog("Unable to open " + fname2 + " for output!\n");
					exit(1);
				}

				string disease = "A ";
				if(options.getDisease() == ""){
					disease += "model";
				}
				else{
					disease += options.getDisease();
				}

				str << "# sampleID";
				trait << "# sampleID";
				for(int f = 0; f < fsize; f++){
					Family* fam = (*families)[f];
					vector<Sample*>* fsamps = fam->getSamples();
					int fssize = fsamps->size();
					for(int s = 0; s < fssize; s++){
						Sample* samp = (*fsamps)[s];
						if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
							str << " " << samp->getFamID() << "_" << samp->getInd();
							str << " " << samp->getFamID() << "_" << samp->getInd();
							trait << " " << samp->getFamID() << "_" << samp->getInd();
							trait << " " << samp->getFamID() << "_" << samp->getInd();
							if(samp->getAffected()){
								disease += " 2 2";
							}
							else{
								disease += " 1 1";
							}
						}
					}
				}
				str << endl;
				str << disease << endl;
				trait << endl;
				trait << disease << endl;
				trait.close();

				for(int m = 0; m < msize; m++){
					Marker* mark = (*markers)[m];
					if(mark->getChrom() == (c + 1) && mark->isEnabled() && !mark->isFlagged()){
						str << "M " << mark->getRSID();
						int mloc = mark->getLoc();
						for(int f = 0; f < fsize; f++){
							Family* fam = (*families)[f];
							vector<Sample*>* fsamps = fam->getSamples();
							int fssize = fsamps->size();
							for(int s = 0; s < fssize; s++){
								Sample* samp = (*fsamps)[s];
								if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
									if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doZeroDisabled())){
										str << " 0 0";
										continue;
									}
									if(!mark->isMicroSat()){
										if(!samp->getAone(mloc) && !samp->getAtwo(mloc)){
											str << " " << map_allele(mark, mark->getAllele1(), &options);;
											str << " " << map_allele(mark, mark->getAllele1(), &options);
										}
										else if(!samp->getAone(mloc) && samp->getAtwo(mloc)){
											str << " " << map_allele(mark, mark->getAllele1(), &options);
											str << " " << map_allele(mark, mark->getAllele2(), &options);
										}
										else if(samp->getAone(mloc) && samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
											str << " " << map_allele(mark, mark->getAllele2(), &options);
											str << " " << map_allele(mark, mark->getAllele2(), &options);
										}
										else{
											str << " 0 0";
										}
									}//end !microsat
									else{
										if(samp->getAbone(mloc) == -1){
											str << " 0 0";
										}
										else{
											str << " " << map_allele(mark, mark->getAllele(samp->getAbone(mloc)), &options);
											str << " " << map_allele(mark, mark->getAllele(samp->getAbtwo(mloc)), &options);
										}
									}//end microsat
								}
							}
						}
						str << endl;
					}
				}
				str.close();
			}
		}
	}

}


int BEAGLEOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

}
