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
*File: PowerMarkerOutput.cc
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
#include "PowerMarkerOutput.h"
#include "General.h"
#include "Helpers.h"
namespace Methods{
string PowerMarkerOutput::stepname = "output-powermarker";


void PowerMarkerOutput::FilterSummary(){
}

void PowerMarkerOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}

}

void PowerMarkerOutput::filter(){
}

Sample* PowerMarkerOutput::find_sample(string i, string f){
	int ssize = samples->size();
	for(int s = 0; s < ssize; s++){
		if((*samples)[s]->getInd() == i && (*samples)[s]->getFamID() == f){
			return ((*samples)[s]);
		}
	}
	return NULL;
}

bool PowerMarkerOutput::find_marker(string p){
	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		if((*markers)[i]->getProbeID() == p){
			return true;
		}
	}
	return false;
}

int PowerMarkerOutput::get_marker_loc(int i){
	int msize = markers->size();
	for(int m = 0; m < msize; m++){
		if((*markers)[m]->getLoc() == i){
			return (*markers)[m]->getLoc();
		}
	}
	return -1;
}

void PowerMarkerOutput::process(vector<Sample*>* ss, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = ss;
	marker_map = mm;

	int ssize = samples->size();
	int msize = markers->size();

	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	msize = good_markers.size();

	string fname1 = opts::_OUTPREFIX_ + "input_powermarker" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".txt";
	}
	string fname2 = opts::_OUTPREFIX_ + "input_powermarker_map" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
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
		//exit(1);
		throw MethodException("Error opening " + fname1 + ". Exiting!\n");
	}
	if(!mout){
		opts::printLog("Error opening " + fname2 + ". Exiting!\n");
		//exit(1);
		throw MethodException("Error opening " + fname2 + ". Exiting!\n");
	}
	bool first = true;
	pout << "Center\tFamID\tIndID\tDadID\tMomID\tGender\tAffection_Status";
	int prev_base = 0;
	int prev_chrom = -1;
	for(int i = 0; i < msize; i++){
		//Marker* mark = (*markers)[mloc];
		Marker* mark = good_markers[i];//(*markers)[i];
		if(mark == NULL){
			//cout << "Marker not found: " << i << endl;
			continue;
		}
		if(!mark->isEnabled()){
			continue;
		}
/*		if(options.doChrom()){
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
           			mark->setFlag(true);
					continue;
           		}
           		prev_base = mark->getBPLOC();
           		prev_chrom = mark->getChrom();
           	}
        }
*/
		pout << "\t" << mark->getRSID();
	}
	pout << endl;
	for(int s = 0; s < ssize; s++){
		Sample* samp = (*samples)[s];
		if(!samp->isEnabled() && !samp->isExcluded() && !options.doIncDisabledSamples()){
			continue;
		}
		if(options.haveCenterCodes()){
			string cent = options.findCenterCode(samp->getFamID() + " " + samp->getInd());
			if(cent != ""){
				pout << cent << "\t";
			}
			else{
				pout << "CEN\t";
			}
		}
		else{
			pout << "CEN\t";
		}
		pout << samp->getFamID() << "\t" << samp->getInd() << "\t" << samp->getDadID() << "\t" << samp->getMomID() << "\t";
		if(samp->getSex()){
			pout << "M\t";
		}
		else{
			pout << "F\t";
		}
		if(samp->getPheno() < 1){
			pout << "U";
		}
		else if(samp->getPheno() == 1){
			pout << "N";
		}
		else if(samp->getPheno() > 1){
			pout << "A";
		}
		//pout << "\t";
		//if(samp->getAffected()){
		//	pout << "2\t";
		//}
		//else{
		//	pout << "1\t";
		//}
		prev_base = 0;
		prev_chrom = -1;
		for(int i = 0; i < msize; i++){
			//Marker* mark = (*markers)[mloc];
			Marker* mark = good_markers[i];//(*markers)[i];
			if(mark == NULL){
				//cout << "Marker not found: " << i << endl;
				continue;
			}
			if(!mark->isEnabled()){
				continue;
			}
/*			if(options.doChrom()){
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
            			mark->setFlag(true);
						continue;
            		}
            		prev_base = mark->getBPLOC();
            		prev_chrom = mark->getChrom();
            	}
            }
*/

			if(first){
				mout << mark->getChrom() << "\t" << mark->getProbeID() << "\t" << mark->getBPLOC() << endl;
		}
			int loc = mark->getLoc();
			if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doZeroDisabled())){
				pout << "\t0/0";
				continue;
			}
			if(!mark->isMicroSat()){
				if(!samp->getAone(loc) && !samp->getAtwo(loc)){
					pout << "\t" << Helpers::map_allele(mark, mark->getAllele1(), &options) << "/" << Helpers::map_allele(mark, mark->getAllele1(), &options);
				}
				else if(!samp->getAone(loc) && samp->getAtwo(loc)){
					if(mark->getAllele1() < mark->getAllele2()){
						pout << "\t" << Helpers::map_allele(mark, mark->getAllele1(), &options) << "/" << Helpers::map_allele(mark, mark->getAllele2(), &options);
					}
					else{
						pout << "\t" << Helpers::map_allele(mark, mark->getAllele2(), &options) << "/" << Helpers::map_allele(mark, mark->getAllele1(), &options);
					}
				}
				else if(samp->getAone(loc) && samp->getAtwo(loc) && !samp->getAmissing(loc)){
					pout << "\t" << Helpers::map_allele(mark, mark->getAllele2(), &options) << "/" << Helpers::map_allele(mark, mark->getAllele2(), &options);
				}
				else if(samp->getAone(loc) && samp->getAtwo(loc) && samp->getAmissing(loc)){
					pout << "\t0/0";
				}
			}
			else{
				if(samp->getAbone(loc) != -1){
					if(mark->getAllele(samp->getAbone(loc)) < mark->getAllele(samp->getAbtwo(loc))){
						pout << "\t" << Helpers::map_allele(mark, mark->getAllele(samp->getAbone(loc)), &options) << "/" << Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(loc)), &options);
					}
					else{
						pout << "\t" << Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(loc)), &options) << "/" << Helpers::map_allele(mark, mark->getAllele(samp->getAbone(loc)), &options);
					}
				}
				else{
					pout << "\t0/0";
				}
			}
		}
		first = false;
		pout << endl;
	}

	pout.close();
	mout.close();
}


int PowerMarkerOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

}
