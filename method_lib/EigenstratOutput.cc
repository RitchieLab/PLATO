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
*File: EigenstratOutput.cc
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
#include "EigenstratOutput.h"
#include "General.h"
#include "Helper.h"

namespace Methods{
void EigenstratOutput::FilterSummary(){
}

void EigenstratOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}

}

void EigenstratOutput::filter(){
}

Sample* EigenstratOutput::find_sample(string i, string f){
	int ssize = samples->size();
	for(int s = 0; s < ssize; s++){
		if((*samples)[s]->getInd() == i && (*samples)[s]->getFamID() == f){
			return ((*samples)[s]);
		}
	}
	return NULL;
}

bool EigenstratOutput::find_marker(string p){
	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		if((*markers)[i]->getProbeID() == p){
			return true;
		}
	}
	return false;
}

int EigenstratOutput::get_marker_loc(int i){
	int msize = markers->size();
	for(int m = 0; m < msize; m++){
		if((*markers)[m]->getLoc() == i){
			return (*markers)[m]->getLoc();
		}
	}
	return -1;
}

void EigenstratOutput::process(vector<Sample*>* ss, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = ss;
	marker_map = mm;

	int ssize = samples->size();
	int msize = markers->size();
	string fname1 = opts::_OUTPREFIX_ + "eigenstrat" + options.getOut() + ".geno";//getString<int>(order) + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".geno";
	}
	string fname2 = opts::_OUTPREFIX_ + "eigenstrat" + options.getOut() + ".snp";//getString<int>(order) + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".snp";
	}
	string fname3 = opts::_OUTPREFIX_ + "eigenstrat" + options.getOut() + ".ind";//getString<int>(order) + ".description";
	if(options.getOverrideOut().size() > 0){
		fname3 = options.getOverrideOut() + ".ind";
	}
	if(!overwrite){
		fname1 += "." + getString<int>(order);
		fname2 += "." + getString<int>(order);
		fname3 += "." + getString<int>(order);
	}
	filenames.push_back(fname1);
	filenames.push_back(fname2);
	filenames.push_back(fname3);
	ofstream genoout(fname1.c_str());
	ofstream snpout (fname2.c_str());
	ofstream indout (fname3.c_str());
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
////	bool first = true;
	map<string, string> dummy;

	string male = "\t0\t0\t1\t0\t";
	string female = "\t0\t0\t1\t0\t";

	int prev_base = 0;
	int prev_chrom = -1;
	bool inddone = false;
	for(int i = 0; i < msize; i++){
		//Marker* mark = (*markers)[mloc];
		Marker* mark = data_set->get_locus(i);
		if(mark == NULL){
			continue;
		}
		if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
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
					}
					if(!samp->getAone(mloc) && !samp->getAtwo(mloc)){
						if(first_samp){
							genoout << "0";
							first_samp = false;
						}
						else{
							genoout << "0";
						}
					}
					else if(!samp->getAone(mloc) && samp->getAtwo(mloc)){
						if(first_samp){
							genoout << "1";
							first_samp = false;
						}
						else{
							genoout << "1";
						}
					}
					else if(samp->getAone(mloc) && samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
						if(first_samp){
							genoout << "2";
							first_samp = false;
						}
						else{
							genoout << "2";
						}
					}
					else if(samp->getAone(mloc) && samp->getAtwo(mloc) && samp->getAmissing(mloc)){
						if(first_samp){
							genoout << "9";
							first_samp = false;
						}
						else{
							genoout << "9";
						}
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
}


int EigenstratOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

/*string EigenstratOutput::map_allele(string a){
	if(options.doAllele1234()){
		if(a == "A"){
			return "1";
		}
		else if(a == "C"){
			return "2";
		}
		else if(a == "G"){
			return "3";
		}
		else if(a == "T"){
			return "4";
		}
		else if(a == opts::_NOCALL_){
			return "0";
		}
		else{
			return "0";
		}
	}
	return a;
}*/

}
