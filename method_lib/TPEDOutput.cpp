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
*File: PEDOutput.cc
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
#include "TPEDOutput.h"
#include "General.h"
#include "Helper.h"

namespace Methods{
Sample* TPEDOutput::find_sample(string i, string f){
	int ssize = samples->size();
	for(int s = 0; s < ssize; s++){
		if((*samples)[s]->getInd() == i && (*samples)[s]->getFamID() == f){
			return ((*samples)[s]);
		}
	}
	return NULL;
}

bool TPEDOutput::find_marker(string p){
	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		if((*markers)[i]->getProbeID() == p){
			return true;
		}
	}
	return false;
}

int TPEDOutput::get_marker_loc(int i){
	int msize = markers->size();
	for(int m = 0; m < msize; m++){
		if((*markers)[m]->getLoc() == i){
			return (*markers)[m]->getLoc();
		}
	}
	return -1;
}

void TPEDOutput::process(vector<Sample*>* ss, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = ss;
	marker_map = mm;

	int ssize = samples->size();
	int msize = markers->size();
	string fname1 = opts::_OUTPREFIX_ + "input_tped" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".ped";
	}
	string fname2 = opts::_OUTPREFIX_ + "input_tped_map" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".map";
	}
	if(!overwrite){
		fname1 += "." + getString<int>(order);
		fname2 += "." + getString<int>(order);
	}
	filenames.push_back(fname1);
	filenames.push_back(fname2);
	ofstream pout(fname1.c_str());
	ofstream mout (fname2.c_str());
	if(!pout){
		opts::printLog("Error opening " + fname1 + ". Exiting!\n");
		throw MethodException("Error opening " + fname1 + ".\n");
	}
	if(!mout){
		opts::printLog("Error opening " + fname2 + ". Exiting!\n");
		throw MethodException("Error opening " + fname2 + ".\n");
	}
	map<string, string> dummy;

	string male = "\t0\t0\t1\t0\t";
	string female = "\t0\t0\t1\t0\t";

	for(int s = 0; s < ssize; s++){
		Sample* samp = (*samples)[s];
		if(!samp->isEnabled() && !samp->isExcluded() && !options.doIncDisabledSamples()){
			continue;
		}
		int valid_parents = 0;
		if(samp->getDad() != NULL){
			valid_parents++;
		}
		if(samp->getMom() != NULL){
			valid_parents++;
		}
		mout << samp->getFamID() << "\t" << samp->getInd() << "\t";
		if((samp->getDad() == NULL || !samp->getDad()->isEnabled()) && samp->getDadID() != "0" && options.doRemMissingParents()){
			mout << "0\t";
		}
		else if(samp->getDad() == NULL && samp->getDadID() != "0" && options.doDummyMissingParents()){
			mout << samp->getDadID() << "\t";
			dummy[samp->getFamID() + "\t" + samp->getDadID()] =  "\t0\t0\t1\t0\t";
		}
		else if(options.doZeroIncompleteTrioIds() && valid_parents != 2){
			mout << "0\t";
		}
		else{
			mout << samp->getDadID() << "\t";
		}
		if((samp->getMom() == NULL || !samp->getMom()->isEnabled()) && samp->getMomID() != "0" && options.doRemMissingParents()){
			mout << "0\t";
		}
		else if(samp->getMom() == NULL && samp->getMomID() != "0" && options.doDummyMissingParents()){
			mout << samp->getMomID() << "\t";
			dummy[samp->getFamID() + "\t" + samp->getMomID()] = "\t0\t0\t2\t0\t";
		}
		else if(options.doZeroIncompleteTrioIds() && valid_parents != 2){
			mout << "0\t";
		}
		else{
		    mout << samp->getMomID() << "\t";
		}
		if(samp->getSex()){
			mout << "1\t";
		}
		else{
			mout << "2\t";
		}
		if(options.getUsePheno()){
			mout << samp->getPheno(options.getPhenoLoc()) << "\n";
		}
		else{
			mout << samp->getPheno() << "\n";
		}
	}//end fam map
	if(options.doDummyMissingParents()){
		map<string, string>::iterator iter;
		for(iter = dummy.begin(); iter != dummy.end(); iter++){
			mout << iter->first << iter->second << endl;
		}
	}
	mout.close();

	int prev_base = 0;
	int prev_chrom = -1;
	for(int i = 0; i < msize; i++){
		//Marker* mark = (*markers)[mloc];
		Marker* mark = (*markers)[i];
		if(mark == NULL){
			//cout << "Marker not found: " << i << endl;
			continue;
		}
		if(!mark->isEnabled()){
			continue;
		}
		if(options.doChrom()){
			if(!options.checkChrom(mark->getChrom())){
				continue;
			}
			if(!options.checkBp(mark->getBPLOC())){
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
					continue;
           		}
           		prev_base = mark->getBPLOC();
           		prev_chrom = mark->getChrom();
          	}
    	}
		pout << mark->getChrom() << "\t" << mark->getRSID() << "\t0\t" << mark->getBPLOC();
		if(options.getMapContainsReferent()){
			pout << "\t" << mark->getReferent();
		}
		pout << "\t";

		bool first_samp = true;
		int loc = mark->getLoc();
		for(unsigned int s = 0; s < samples->size(); s++){
			Sample* samp = (*samples)[s];
			if((samp->isExcluded() && options.doZeroExcluded()) || (options.doZeroDisabled() && !samp->isEnabled())){
				if(first_samp){
					pout << "0 0";
					first_samp = false;
				}
				else{
					pout << " 0 0";
				}
				continue;
			}
			else if(!samp->isEnabled()){
				continue;
			}
			if(!first_samp){
				pout << " ";
			}
			if(!mark->isMicroSat()){
				first_samp = false;
				if(!samp->getAone(loc) && !samp->getAtwo(loc)){
					pout << map_allele(mark, mark->getAllele1(),&options) << " " << map_allele(mark, mark->getAllele1(),&options);
				}
				else if(!samp->getAone(loc) && samp->getAtwo(loc)){
					if(mark->getAllele1() < mark->getAllele2()){
						pout << map_allele(mark, mark->getAllele1(),&options) << " " << map_allele(mark, mark->getAllele2(),&options);
					}
					else{
						pout << map_allele(mark, mark->getAllele2(),&options) << " " << map_allele(mark, mark->getAllele1(),&options);
					}
				}
				else if(samp->getAone(loc) && samp->getAtwo(loc) && !samp->getAmissing(loc)){
					pout << map_allele(mark, mark->getAllele2(),&options) << " " << map_allele(mark, mark->getAllele2(),&options);
				}
				else if(samp->getAone(loc) && samp->getAtwo(loc) && samp->getAmissing(loc)){
					pout << "0 0";
				}
			}
			else{
				first_samp = false;
				if(samp->getAbone(loc) != -1){
					if(mark->getAllele(samp->getAbone(loc)) < mark->getAllele(samp->getAbtwo(loc))){
						pout << map_allele(mark, mark->getAllele(samp->getAbone(loc)),&options) << " " << map_allele(mark, mark->getAllele(samp->getAbtwo(loc)),&options);
					}
					else{
						pout << map_allele(mark, mark->getAllele(samp->getAbtwo(loc)),&options) << " " << map_allele(mark, mark->getAllele(samp->getAbone(loc)),&options);
					}
				}
				else{
					pout << "0 0";
				}
			}
			if(options.doDummyMissingParents()){
				map<string, string>::iterator iter;
				for(iter = dummy.begin(); iter != dummy.end(); iter++){
					pout << " 0 0";
				}
			}
		}
		pout << "\n";
	}

	pout.close();
}


int TPEDOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}
}
