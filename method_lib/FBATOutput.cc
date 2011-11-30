/**********************************************************************************
*                       FBAT Output Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates FBAT input files.
*
*
*
*File: FBATOutput.cc
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
#include "FBATOutput.h"
#include "General.h"
#include "Helpers.h"

namespace Methods{
void FBATOutput::FilterSummary(){
}

void FBATOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}
}

void FBATOutput::filter(){
}

void FBATOutput::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
	orig_num_markers = markers->size();
	orig_num_families = families->size();
 	orig_num_individuals = samples->size();

   	int ssize = samples->size();
	int msize = markers->size();
	string fname1 = opts::_OUTPREFIX_ + "input_fbat" + options.getOut() + ".txt";
    if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".txt";
	}

	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	filenames.push_back(fname1);
	ofstream fbat (fname1.c_str());
	if(!fbat.is_open()){
		opts::printLog("Unable to open " + fname1 + " for output!\n");
		throw MethodException("Unable to open " + fname1 + " for output!\n");
	}
	string fname2 = opts::_OUTPREFIX_ + "input_fbat" + options.getOut() + ".map";
    if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".map";
	}
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	filenames.push_back(fname2);
	ofstream fbat_map (fname2.c_str());
	if(!fbat_map.is_open()){
		opts::printLog("Unable to open " + fname2 +" for output!\n");
		throw MethodException("Unable to open " + fname2 +" for output!\n");
	}
	bool map_done = false;

	int prev_base = 0;
	int prev_chrom = -1;
	prev_base = 0;
	prev_chrom = -1;
	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	msize = good_markers.size();
	for(int m = 0; m < msize; m++){
		Marker* mark = good_markers[m];
		if(mark->isEnabled()){
			fbat << mark->getRSID() << " ";
		}
	}
	fbat << endl;
	prev_base = 0;
	prev_chrom = -1;
	for(int i = 0; i < ssize; i++){
		Sample* samp = (*samples)[i];
		if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
			fbat << samp->getFamID() << "\t" << samp->getInd() << "\t" << samp->getDadID() << "\t" << samp->getMomID() << "\t";
			if(samp->getSex()){
				fbat << "1";
			}
			else{
				fbat << "2";
			}
			if(samp->getPheno() > 2){
				fbat << " 0";
			}
			else{
				fbat << " " << samp->getPheno();
			}
			for(int m = 0; m < msize; m++){
				Marker* mark = good_markers[m];
				if(mark->isEnabled()){
					int m_loc = mark->getLoc();
					if(!map_done){
						double cm = (double)mark->getBPLOC() / (double)1000000;
						fbat_map << mark->getRSID() << "\t" << mark->getChrom() << "\t" << cm << "\t" << mark->getBPLOC();
						if(mark->getChrom() == opts::_CHRX_){
							fbat_map << "\t1" << endl;
						}
						else{
							fbat_map << "\t0" << endl;
						}
					}
					if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doZeroDisabled())){
						fbat << " 0 0";
						continue;
					}
					if(!mark->isMicroSat()){
						if(!samp->getAone(m_loc) && !samp->getAtwo(m_loc)){
							fbat << " " << Helpers::map_allele(mark, mark->getAllele1(), &options) << " " << Helpers::map_allele(mark, mark->getAllele1(), &options);
						}
						else if(!samp->getAone(m_loc) && samp->getAtwo(m_loc)){
							if(mark->getAllele1() > mark->getAllele2()){
								fbat << " " << Helpers::map_allele(mark, mark->getAllele2(), &options) << " " << Helpers::map_allele(mark, mark->getAllele1(), &options);
							}
							else{
								fbat << " " << Helpers::map_allele(mark, mark->getAllele1(), &options) << " " << Helpers::map_allele(mark, mark->getAllele2(), &options);
							}
						}
						else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && !samp->getAmissing(m_loc)){
							fbat << " " << Helpers::map_allele(mark, mark->getAllele2(), &options) << " " << Helpers::map_allele(mark, mark->getAllele2(), &options);
						}
						else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && samp->getAmissing(m_loc)){
							fbat << " 0 0";
						}
					}
					else{
						if(samp->getAbone(m_loc) != -1){
							if(mark->getAllele(samp->getAbone(m_loc)) < mark->getAllele(samp->getAbtwo(m_loc))){
								fbat << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options) << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options);
							}
							else{
								fbat << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options) << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options);
							}
						}
						else{
							fbat << " 0 0";
						}
					}
				}
			}
			map_done = true;
			fbat << endl;
		}
	}
	fbat.close();
	fbat_map.close();
}


int FBATOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

}
