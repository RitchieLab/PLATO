/**********************************************************************************
*                       Eigenstrat Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
*Outputs Eigenstrat formatted input files
*
*File: EigenstratOutput.cc
**********************************************************************************/


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
#include <list>
#include <algorithm>
#include <map>
#include "EigenstratOutput.h"
#include "General.h"
#include "Helpers.h"

namespace Methods{
//DEPRECATED
void EigenstratOutput::FilterSummary(){
}

//DEPRECATED
void EigenstratOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers).at(i)->setFlag(false);
	}

}

//DEPRECATED
void EigenstratOutput::filter(){
}

//DEPRECATED
Sample* EigenstratOutput::find_sample(string i, string f){
	int ssize = samples->size();
	for(int s = 0; s < ssize; s++){
		if((*samples).at(s)->getInd() == i && (*samples).at(s)->getFamID() == f){
			return ((*samples).at(s));
		}
	}
	return NULL;
}

//DEPRECATED
bool EigenstratOutput::find_marker(string p){
	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		if((*markers).at(i)->getProbeID() == p){
			return true;
		}
	}
	return false;
}

//DEPRECATED
int EigenstratOutput::get_marker_loc(int i){
	int msize = markers->size();
	for(int m = 0; m < msize; m++){
		if((*markers).at(m)->getLoc() == i){
			return (*markers).at(m)->getLoc();
		}
	}
	return -1;
}

//Main method to generate files
//called by this->calculate();
void EigenstratOutput::process(vector<Sample*>* ss, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = ss;
	marker_map = mm;

	int ssize = samples->size();
	int msize = markers->size();
	string fname1 = opts::_OUTPREFIX_ + "eigenstrat" + options.getOut() + ".geno";
	if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".geno";
	}
	string fname2 = opts::_OUTPREFIX_ + "eigenstrat" + options.getOut() + ".snp";
	if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".snp";
	}
	string fname3 = opts::_OUTPREFIX_ + "eigenstrat" + options.getOut() + ".ind";
	if(options.getOverrideOut().size() > 0){
		fname3 = options.getOverrideOut() + ".ind";
	}
	string fname4 = opts::_OUTPREFIX_ + "eigenstrat" + options.getOut() + ".phenoQTL";
	if(options.getOverrideOut().size() > 0){
		fname4 = options.getOverrideOut() + ".phenoQTL";
	}
	string fname5 = opts::_OUTPREFIX_ + "eigenstrat" + options.getOut() + ".ancestrymapgeno";
	if(options.getOverrideOut().size() > 0){
		fname5 = options.getOverrideOut() + ".ancestrymapgeno";
	}
	if(!overwrite){
		fname1 += "." + getString<int>(order);
		fname2 += "." + getString<int>(order);
		fname3 += "." + getString<int>(order);
		fname4 += "." + getString<int>(order);
		fname5 += "." + getString<int>(order);
	}
	filenames.push_back(fname1);
	filenames.push_back(fname2);
	filenames.push_back(fname3);
	if(options.doQTL()){
		filenames.push_back(fname4);
	}
	if(options.doAncestry()){
		filenames.push_back(fname5);
	}

	ofstream genoout(fname1.c_str());
	ofstream snpout (fname2.c_str());
	ofstream indout (fname3.c_str());
	ofstream phenoout;
	ofstream ancestryout;

	if(!genoout){
		opts::printLog("Error opening " + fname1 + ". Exiting!\n");
		throw MethodException("Error opening " + fname1 + ". Exiting!\n");
	}
	if(!snpout){
		opts::printLog("Error opening " + fname2 + ". Exiting!\n");
		throw MethodException("Error opening " + fname2 + ". Exiting!\n");
	}
	if(!indout){
		opts::printLog("Error opening " + fname3 + ". Exiting!\n");
		throw MethodException("Error opening " + fname3 + ". Exiting!\n");
	}
	if(options.doQTL()){
		phenoout.open(fname4.c_str());
		if(!phenoout){
			opts::printLog("Error opening " + fname4 + ". Exiting!\n");
			throw MethodException("Error opening " + fname4 + ". Exiting!\n");
		}
	}
	if(options.doAncestry()){
		ancestryout.open(fname5.c_str());
		if(!ancestryout){
			opts::printLog("Error opening " + fname5 + ". Exiting!\n");
			throw MethodException("Error opening " + fname5 + ". Exiting!\n");
		}
	}
	map<string, string> dummy;

	string male = "\t0\t0\t1\t0\t";
	string female = "\t0\t0\t1\t0\t";

	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	msize = good_markers.size();
	bool inddone = false;
	for(int i = 0; i < msize; i++){
		Marker* mark = good_markers.at(i);
		if(mark == NULL){
			continue;
		}
		if(mark->isEnabled()){
			double morgans = (double)mark->getBPLOC() / 1000000.0f;
			snpout << mark->getProbeID() << "\t" << mark->getChrom() << "\t" << morgans << "\t" << mark->getBPLOC() << endl;
			int mloc = mark->getLoc();
			bool first_samp = true;
			for(int s = 0; s < ssize; s++){
				Sample* samp = data_set->get_sample(s);
				if(samp->isEnabled() && !samp->isExcluded() && !options.doIncDisabledSamples()){

					if(!inddone){
						if(samp->getDad() == NULL && samp->getDadID() != "0" && options.doDummyMissingParents()){
							dummy[samp->getFamID() + "-" + samp->getDadID()] =  "\tM\tIgnore";
						}
						if(samp->getMom() == NULL && samp->getMomID() != "0" && options.doDummyMissingParents()){
							dummy[samp->getFamID() + "-" + samp->getMomID()] =  "\tF\tIgnore";
						}
						indout << samp->getFamID() << "-" << samp->getInd() << "\t";
						if(samp->getSex()){
							indout << "M";
						}
						else{
							indout << "F";
						}
						if(samp->getPheno() == 2){
							indout << "\tCase";
						}
						else if(samp->getPheno() == 1){
							indout << "\tControl";
						}
						else{
							indout << "\tIgnore";
						}
						indout << endl;

						if(options.doQTL()){
						  if (options.getUsePheno()) {
								int index = options.getPhenoLoc();

								if (options.getPhenoName() != "") {
									index = data_set->get_trait_index(options.getPhenoName());
								}
								if (index < 0) {
									throw MethodException(
											"Internal Error: Trait/Phenotype index value < 0 in Eigenstrat Output!");
								}
								if(data_set->get_sample(s)->getPheno(index) == options.getPhenoMissing()){
									phenoout << "-100.0\n";
								}
								else{
									phenoout << data_set->get_sample(s)->getPheno(index) << "\n";
								}
							} else {
								if(data_set->get_sample(s)->getPheno() == options.getPhenoMissing()){
									phenoout << "-100.0\n";
								}
								else{
									phenoout << data_set->get_sample(s)->getPheno() << endl;
								}
							}
						}
					}
					if(options.doAncestry()){
						ancestryout << mark->getRSID() << " " << samp->getFamID() << "-" << samp->getInd() << " ";
					}
					if(!samp->getAone(mloc) && !samp->getAtwo(mloc)){
						if(first_samp){
							genoout << "0";
							if(options.doAncestry()){
								ancestryout << "2";
							}
							first_samp = false;
						}
						else{
							genoout << "0";
							if(options.doAncestry()){
								ancestryout << "2";
							}
						}
					}
					else if(!samp->getAone(mloc) && samp->getAtwo(mloc)){
						if(first_samp){
							genoout << "1";
							if(options.doAncestry()){
								ancestryout << "1";
							}
							first_samp = false;
						}
						else{
							genoout << "1";
							if(options.doAncestry()){
								ancestryout << "1";
							}
						}
					}
					else if(samp->getAone(mloc) && samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
						if(first_samp){
							genoout << "2";
							if(options.doAncestry()){
								ancestryout << "0";
							}
							first_samp = false;
						}
						else{
							genoout << "2";
							if(options.doAncestry()){
								ancestryout << "0";
							}
						}
					}
					else if(samp->getAone(mloc) && samp->getAtwo(mloc) && samp->getAmissing(mloc)){
						if(first_samp){
							genoout << "9";
							if(options.doAncestry()){
								ancestryout << "-1";
							}
							first_samp = false;
						}
						else{
							genoout << "9";
							if(options.doAncestry()){
								ancestryout << "-1";
							}
						}
					}
					if(options.doAncestry()){
						ancestryout << endl;
					}
				}
			}
			if(options.doDummyMissingParents()){
				map<string, string>::iterator iter;
				for(iter = dummy.begin(); iter != dummy.end(); iter++){
					if(!inddone){
						indout << iter->first << iter->second << endl;
					}
					genoout << "9";
				}
			}
			genoout << endl;

			inddone = true;
		}
	}



	genoout.close();
	indout.close();
	snpout.close();
	if(phenoout){
		phenoout.close();
	}
}


int EigenstratOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

}
