/**********************************************************************************
*                        GRR Output Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates GRR input files.
*
*
*File: GRROutput.cc
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
#include "GRROutput.h"
#include "General.h"
#include "Helpers.h"

namespace Methods{
void GRROutput::FilterSummary(){
}

void GRROutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}
}

void GRROutput::filter(){
}

void GRROutput::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
	orig_num_markers = markers->size();
	orig_num_families = families->size();
 	orig_num_individuals = samples->size();

	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	vector<vector<Marker*> > marker_sets;
	if(options.doRandomMarkers() || options.doRandomRepeat()){
		vector<Marker*> used_markers;
		for(int ss = 0; ss < options.getSets(); ss++){
			vector<Marker*> temp = Helpers::findRandomMarkers(good_markers, &used_markers, &options);
			if(temp.size() > 0){
				marker_sets.push_back(temp);
			}
		}
	}
	if(marker_sets.size() == 0){
		marker_sets.push_back(good_markers);
	}

   	int ssize = samples->size();
//	int msize = good_markers.size();

	for(int k = 0; k < (int)marker_sets.size(); k++){
		vector<Marker*> marks = marker_sets[k];
		int msize = marks.size();

		string fname1 = opts::_OUTPREFIX_ + "input_grr_" + getString<int>(k + 1) + options.getOut() + ".txt";//getString<int>(order) + ".txt";
        if(options.getOverrideOut().size() > 0){
			fname1 = options.getOverrideOut() + "_" + getString<int>(k + 1) + ".txt";
		}

		if(!overwrite){
			fname1 += "." + getString<int>(order);
		}
		filenames.push_back(fname1);
		ofstream grr (fname1.c_str());
		if(!grr.is_open()){
			opts::printLog("Unable to open " +fname1+" for output!\n");
	   		//exit(1);
			throw MethodException("Unable to open " +fname1+" for output!\n");
		}
		string fname2 = opts::_OUTPREFIX_ + "input_grr_" + getString<int>(k + 1) + options.getOut() + ".map";//getString<int>(order) + ".map";
        if(options.getOverrideOut().size() > 0){
			fname2 = options.getOverrideOut() + "_" + getString<int>(k + 1) + ".map";
		}
		if(!overwrite){
			fname2 += "." + getString<int>(order);
		}
		filenames.push_back(fname2);
		ofstream grr_map (fname2.c_str());
		if(!grr_map.is_open()){
			opts::printLog("Unable to open "+fname2+" for output!\n");
			//exit(1);
			throw MethodException("Unable to open "+fname2+" for output!\n");
		}
		grr_map << "Chrom\tBPLOC\tRsid\tA1/A2" << endl;
		bool map_done = false;

		for(int i = 0; i < ssize; i++){
			Sample* samp = (*samples)[i];
			if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
////				int prev_base = 0;
////				int prev_chrom = -1;
				grr << samp->getFamID() << "\t" << samp->getInd() << "\t" << samp->getDadID() << "\t" << samp->getMomID() << "\t";
				if(samp->getSex()){
					grr << "1";
				}
				else{
					grr << "2";
				}
				for(int m = 0; m < msize; m++){
					Marker* mark = (marks)[m];
					if(mark->isEnabled()){
//						if(options.doChrom()){
//							if(!options.checkChrom(mark->getChrom())){
//							    continue;
//						    }
//						    if(!options.checkBp(mark->getBPLOC())){
//							    continue;
//						    }
//						}

//						if(options.doBpSpace()){
//							if(prev_base == 0){
//								prev_base = mark->getBPLOC();
//								prev_chrom = mark->getChrom();
//							}
//							else{
//								if(mark->getChrom() == prev_chrom && ((mark->getBPLOC() - prev_base) < options.getBpSpace())){
//									mark->setFlag(true);
//									continue;
//								}
//								prev_base = mark->getBPLOC();
//								prev_chrom = mark->getChrom();
//							}
//						}
						int m_loc = mark->getLoc();
						if(!map_done){
							grr_map << mark->getChrom() << "\t" << mark->getBPLOC() << "\t" << mark->getRSID() << "\t" << Helpers::map_allele(mark, mark->getAllele1(), &options) << "/" << Helpers::map_allele(mark, mark->getAllele2(), &options) << endl;
						}
						if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doZeroDisabled())){
							grr << "\t0/0";
							continue;
						}
						if(!mark->isMicroSat()){
							if(!samp->getAone(m_loc) && !samp->getAtwo(m_loc)){
								grr << "\t" << Helpers::map_allele(mark, mark->getAllele1(), &options) << "/" << Helpers::map_allele(mark, mark->getAllele1(), &options);
							}
							else if(!samp->getAone(m_loc) && samp->getAtwo(m_loc)){
								if(mark->getAllele1() > mark->getAllele2()){
									grr << "\t" << Helpers::map_allele(mark, mark->getAllele2(), &options) << "/" << Helpers::map_allele(mark, mark->getAllele1(), &options);
								}
								else{
									grr << "\t" << Helpers::map_allele(mark, mark->getAllele1(), &options) << "/" << Helpers::map_allele(mark, mark->getAllele2(), &options);
								}
							}
							else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && !samp->getAmissing(m_loc)){
								grr << "\t" << Helpers::map_allele(mark, mark->getAllele2(), &options) << "/" << Helpers::map_allele(mark, mark->getAllele2(), &options);
							}
							else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && samp->getAmissing(m_loc)){
								grr << "\t0/0";
							}
						}
						else{
							if(samp->getAbone(m_loc) != -1){
								if(mark->getAllele(samp->getAbone(m_loc)) < mark->getAllele(samp->getAbtwo(m_loc))){
									grr << "\t" << Helpers::map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options) << "/" << Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options);
								}
								else{
									grr << "\t" << Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options) << "/" << Helpers::map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options);
								}
							}
							else{
								grr << "\t0/0";
							}
						}
					}
				}
				map_done = true;
				grr << endl;
			}
		}
		grr.close();
		grr_map.close();
	}
	good_markers.clear();
}


int GRROutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

}
