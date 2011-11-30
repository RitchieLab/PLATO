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
*File: ImputeOutput.cc
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
#include "SampleGenoEff.h"
#include "ImputeOutput.h"
#include "General.h"
#include "Helpers.h"

namespace Methods{
void ImputeOutput::FilterSummary(){
}

void ImputeOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers).at(i)->setFlag(false);
	}

}

void ImputeOutput::filter(){
}

Sample* ImputeOutput::find_sample(string i, string f){
	int ssize = samples->size();
	for(int s = 0; s < ssize; s++){
		if((*samples).at(s)->getInd() == i && (*samples).at(s)->getFamID() == f){
			return ((*samples).at(s));
		}
	}
	return NULL;
}

bool ImputeOutput::find_marker(string p){
	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		if((*markers).at(i)->getProbeID() == p){
			return true;
		}
	}
	return false;
}

int ImputeOutput::get_marker_loc(int i){
	int msize = markers->size();
	for(int m = 0; m < msize; m++){
		if((*markers).at(m)->getLoc() == i){
			return (*markers).at(m)->getLoc();
		}
	}
	return -1;
}

void ImputeOutput::process(vector<Sample*>* ss, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = ss;
	marker_map = mm;

	SampleGenoEff sge(data_set);
	sge.setOptions(options);

	int ssize = samples->size();
	int msize = markers->size();
	string fname1 = opts::_OUTPREFIX_ + "impute_format" + options.getOut() + ".gen";
	if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".gen";
	}
	string fname2 = opts::_OUTPREFIX_ + "impute_format" + options.getOut() + ".map";
	if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".map";
	}
	string fname3 = opts::_OUTPREFIX_ + "impute_format" + options.getOut() + ".sample";
	if(options.getOverrideOut().size() > 0){
		fname3 = options.getOverrideOut() + ".sample";
	}
	if(!overwrite){
		fname1 += "." + getString<int>(order);
		fname2 += "." + getString<int>(order);
		fname3 += "." + getString<int>(order);
	}
	filenames.push_back(fname1);
	filenames.push_back(fname2);
	filenames.push_back(fname3);
	ofstream genout(fname1.c_str());
	ofstream mapout (fname2.c_str());
	ofstream sampout (fname3.c_str());
	if(!genout){
		opts::printLog("Error opening " + fname1 + ". Exiting!\n");
		throw MethodException("Error opening " + fname1 + ". Exiting!\n");
	}
	if(!mapout){
		opts::printLog("Error opening " + fname2 + ". Exiting!\n");
		throw MethodException("Error opening " + fname2 + ". Exiting!\n");
	}
	if(!sampout){
		opts::printLog("Error opening " + fname3 + ". Exiting!\n");
		throw MethodException("Error opening " + fname3 + ". Exiting!\n");
	}
	bool first = true;
	map<string, string> dummy;

	sampout << "ID_1 ID_2 missing sex";
	//initialize with 1 for zero proportion
	string dummy_covs = "1 -9";
	vector<string>* covs = data_set->get_covariates();
	for(int c = 0; c < (int)covs->size(); c++){
		dummy_covs += " -9";
		sampout << " " << covs->at(c);
	}
	//add dummy phenotype
	dummy_covs += " -9";
	sampout << " phenotype" << endl;
	sampout << "0 0 0 2";
	for(int c = 0; c < (int)covs->size(); c++){
		sampout << " " << 3;
	}
	sampout << " P" << endl;

	for(int s = 0; s < ssize; s++){
		Sample* samp = (*samples).at(s);
		if(!samp->isEnabled() && !samp->isExcluded() && !options.doIncDisabledSamples()){
			continue;
		}

		sge.calculate(samp);

		int valid_parents = 0;
		if(samp->getDad() != NULL){
			valid_parents++;
		}
		if(samp->getMom() != NULL){
			valid_parents++;
		}

		sampout << samp->getFamID() << " " << samp->getInd() << " " << (1.0f - sge.getPercent());

		if(samp->getSex()){
			sampout << " 1";
		}
		else{
			sampout << " 2";
		}
		for(int c = 0; c < (int)covs->size(); c++){
			double cov = samp->getCovariate(c);
			if(getString<double>(cov) == opts::_COVAR_MISSING_){
				cov = -9;
			}
			sampout << " " << cov;
		}

		if(samp->getDad() == NULL && samp->getDadID() != "0" && options.doDummyMissingParents()){
			dummy[samp->getFamID() + " " + samp->getDadID()] =  dummy_covs;
		}
		else if(samp->getDad() == NULL && samp->getMom() != NULL && options.doDummyIncompleteParentIds()){
			string did = "DUM999999";
			if(samp->getDadID() != "0"){
				did = samp->getDadID();
			}
			dummy[samp->getFamID() + " " + did] = dummy_covs;
		}
		if(samp->getMom() == NULL && samp->getMomID() != "0" && options.doDummyMissingParents()){
			dummy[samp->getFamID() + " " + samp->getMomID()] = dummy_covs;
		}
		else if(samp->getMom() == NULL && samp->getDad() != NULL && options.doDummyIncompleteParentIds()){
			string mid = "DUM999998";
			if(samp->getMomID() != "0"){
				mid = samp->getMomID();
			}
			dummy[samp->getFamID() + " " + mid] = dummy_covs;
		}
		if(options.getUsePheno()){
			if(options.getPhenoName() != ""){
				float pheno = samp->getPheno(data_set->get_trait_index(options.getPhenoName()));
				if(pheno <= 0){
					pheno = -9;
				}
				sampout << " " << pheno;
			}
			else{
				float pheno = samp->getPheno(options.getPhenoLoc());
				if(pheno <= 0){
					pheno = -9;
				}
				sampout << " " << pheno;
			}
		}
		else{
			if(opts::_BINTRAIT_ && samp->getPheno() == 0){
				sampout << " " << "-9";
			}
			else{
				sampout <<  " " << samp->getPheno();
			}
		}
		sampout << endl;
	}
	if(options.doDummyMissingParents() || options.doDummyIncompleteParentIds()){

		map<string, string>::iterator iter;
		for(iter = dummy.begin(); iter != dummy.end(); iter++){
			sampout << iter->first << " " << iter->second << endl;
		}
	}
		vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
		msize = good_markers.size();
		for(int i = 0; i < msize; i++){
			Marker* mark = good_markers.at(i);
			Marker* mark2 = NULL;
			if((i + 1) < msize){
				mark2 = good_markers.at(i + 1);
			}
			if(mark == NULL){
				continue;
			}
			if(!mark->isEnabled()){
				continue;
			}
            genout << mark->getChrom() << " " << mark->getProbeID() << " " << mark->getBPLOC() << " "
					<< mark->getAllele1() << " " << mark->getAllele2();

			if(first){
				mapout << mark->getBPLOC();
				if(mark2 == NULL || mark->getChrom() != mark2->getChrom()){
					mapout << " 0 0" << endl;
				}
				else{
					mapout << " " << (((double)mark2->getBPLOC() - (double)mark->getBPLOC()) / 1000000.0f) << " " << ((double)mark2->getBPLOC() / 1000000.0f) << endl;
				}
			}
			int loc = mark->getLoc();

			for(int s = 0; s < ssize; s++){
				Sample* samp = data_set->get_sample(s);
				if(!samp->isEnabled()){
					continue;
				}
				bool a1 = samp->getAone(loc);
				bool a2 = samp->getAtwo(loc);

				string value = "0 0 0";
				if(!a1 && !a2){
					value = "1 0 0";
				}
				else if(!a1 && a2){
					value = "0 1 0";
				}
				else if(a1 && a2){
					value = "0 0 1";
				}
				genout << " " << value;
			}

			if(options.doDummyMissingParents() || options.doDummyIncompleteParentIds()){

				map<string, string>::iterator iter;
				for(iter = dummy.begin(); iter != dummy.end(); iter++){
					genout << " 0 0 0";

				}
			}


			genout << endl;

		}

	sampout.close();
	genout.close();
	mapout.close();
}


int ImputeOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

}
