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
#include "PEDOutput.h"
#include "General.h"
#include "Helper.h"

namespace Methods{
void PEDOutput::FilterSummary(){
}

void PEDOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}

}

void PEDOutput::filter(){
}

Sample* PEDOutput::find_sample(string i, string f){
	int ssize = samples->size();
	for(int s = 0; s < ssize; s++){
		if((*samples)[s]->getInd() == i && (*samples)[s]->getFamID() == f){
			return ((*samples)[s]);
		}
	}
	return NULL;
}

bool PEDOutput::find_marker(string p){
	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		if((*markers)[i]->getProbeID() == p){
			return true;
		}
	}
	return false;
}

int PEDOutput::get_marker_loc(int i){
	int msize = markers->size();
	for(int m = 0; m < msize; m++){
		if((*markers)[m]->getLoc() == i){
			return (*markers)[m]->getLoc();
		}
	}
	return -1;
}

void PEDOutput::process(vector<Sample*>* ss, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = ss;
	marker_map = mm;

	int ssize = samples->size();
	int msize = markers->size();
	string fname1 = opts::_OUTPREFIX_ + "input_ped" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".txt";
	}
	string fname2 = opts::_OUTPREFIX_ + "input_ped_map" + options.getOut() + ".map";//getString<int>(order) + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".map";
	}
	string fname3 = opts::_OUTPREFIX_ + "input_ped_map" + options.getOut() + ".description";//getString<int>(order) + ".description";
	if(options.getOverrideOut().size() > 0){
		fname3 = options.getOverrideOut() + ".description";
	}
	if(!overwrite){
		fname1 += "." + getString<int>(order);
		fname2 += "." + getString<int>(order);
		fname3 += "." + getString<int>(order);
	}
	filenames.push_back(fname1);
	filenames.push_back(fname2);
	filenames.push_back(fname3);
	ofstream pout(fname1.c_str());
	ofstream mout (fname2.c_str());
	ofstream mdout (fname3.c_str());
	if(!pout){
		opts::printLog("Error opening " + fname1 + ". Exiting!\n");
		exit(1);
	}
	if(!mout){
		opts::printLog("Error opening " + fname2 + ". Exiting!\n");
		exit(1);
	}
	if(!mdout){
		opts::printLog("Error opening " + fname3 + ". Exiting!\n");
		exit(1);
	}
	bool first = true;
	map<string, string> dummy;

	string male = "\t0\t0\t1\t0\t";
	string female = "\t0\t0\t2\t0\t";

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
		pout << samp->getFamID() << "\t" << samp->getInd() << "\t";
		if((samp->getDad() == NULL || !samp->getDad()->isEnabled()) && samp->getDadID() != "0" && options.doRemMissingParents()){
			pout << "0\t";
		}
		else if(samp->getDad() == NULL && samp->getDadID() != "0" && options.doDummyMissingParents()){
			pout << samp->getDadID() << "\t";
			dummy[samp->getFamID() + "\t" + samp->getDadID()] =  "\t0\t0\t1\t0\t";
		}
		else if(samp->getDad() == NULL && samp->getMom() != NULL && options.doDummyIncompleteParentIds()){
			string did = "DUM999999";
			if(samp->getDadID() != "0"){
				did = samp->getDadID();
			}
			pout << did << "\t";
			dummy[samp->getFamID() + "\t" + did] = "\t0\t0\t1\t0\t";
		}
		else if(options.doZeroIncompleteTrioIds() && valid_parents != 2){
			pout << "0\t";
		}
		else{
			pout << samp->getDadID() << "\t";
		}
		if((samp->getMom() == NULL || !samp->getMom()->isEnabled()) && samp->getMomID() != "0" && options.doRemMissingParents()){
			pout << "0\t";
		}
		else if(samp->getMom() == NULL && samp->getMomID() != "0" && options.doDummyMissingParents()){
			pout << samp->getMomID() << "\t";
			dummy[samp->getFamID() + "\t" + samp->getMomID()] = "\t0\t0\t2\t0\t";
		}
		else if(samp->getMom() == NULL && samp->getDad() != NULL && options.doDummyIncompleteParentIds()){
			string mid = "DUM999998";
			if(samp->getMomID() != "0"){
				mid = samp->getMomID();
			}
			pout << mid << "\t";
			dummy[samp->getFamID() + "\t" + mid] = "\t0\t0\t2\t0\t";
		}
		else if(options.doZeroIncompleteTrioIds() && valid_parents != 2){
			pout << "0\t";
		}
		else{
		    pout << samp->getMomID() << "\t";
		}
		if(samp->getSex()){
			pout << "1\t";
		}
		else{
			pout << "2\t";
		}
		if(options.getUsePheno()){
			pout << samp->getPheno(options.getPhenoLoc()) << "\t";
		}
		else{
			pout << samp->getPheno() << "\t";
		}
		//if(samp->getAffected()){
		//	pout << "2\t";
		//}
		//else{
		//	pout << "1\t";
		//}
		int prev_base = 0;
		int prev_chrom = -1;
		bool first_marker = true;
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
            			mark->setFlag(true);
						continue;
            		}
            		prev_base = mark->getBPLOC();
            		prev_chrom = mark->getChrom();
            	}
            }


			if(first){
				mout << mark->getChrom() << "\t" << mark->getProbeID() << "\t" << mark->getBPLOC() << endl;
				mdout << mark->getChrom() << "\t" << mark->getProbeID() << "\t" << mark->getBPLOC();
			   	if(mark->getDetails() != ""){
					mdout << mark->getDetails();
				}
				mdout << "\t" << mark->getAllele1() << "\t" << mark->getAllele2() << endl;
			}
			int loc = mark->getLoc();
			if(!first_marker){
				pout << " ";
			}
			if((samp->isExcluded() && options.doZeroExcluded()) || (options.doZeroDisabled() && !samp->isEnabled())){
				first_marker = false;
				if(!options.get_allele_as_snp()){
					pout << "0 0";
				}
				else{
					pout << mark->getProbeID() << "_0 " << mark->getProbeID() << "_0";
				}
				continue;
			}
			if(!mark->isMicroSat()){
				first_marker = false;
				if(!samp->getAone(loc) && !samp->getAtwo(loc)){
					if(!options.get_allele_as_snp()){
						pout << map_allele(mark, mark->getAllele1(),&options) << " " << map_allele(mark, mark->getAllele1(),&options);
					}
					else{
						pout << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele1(), &options) << " " << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele1(), &options);
					}
				}
				else if(!samp->getAone(loc) && samp->getAtwo(loc)){
					if(mark->getAllele1() < mark->getAllele2()){
						if(!options.get_allele_as_snp()){
							pout << map_allele(mark, mark->getAllele1(),&options) << " " << map_allele(mark, mark->getAllele2(),&options);
						}
						else{
							pout << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele1(), &options) << " " << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele2(), &options);
						}
					}
					else{
						if(!options.get_allele_as_snp()){
							pout << map_allele(mark, mark->getAllele2(),&options) << " " << map_allele(mark, mark->getAllele1(),&options);
						}
						else{
							pout << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele2(), &options) << " " << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele1(), &options);
						}
					}
				}
				else if(samp->getAone(loc) && samp->getAtwo(loc) && !samp->getAmissing(loc)){
					if(!options.get_allele_as_snp()){
						pout << map_allele(mark, mark->getAllele2(),&options) << " " << map_allele(mark, mark->getAllele2(),&options);
					}
					else{
						pout << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele2(), &options) << " " << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele2(), &options);
					}
				}
				else if(samp->getAone(loc) && samp->getAtwo(loc) && samp->getAmissing(loc)){
					if(!options.get_allele_as_snp()){
						pout << "0 0";
					}
					else{
						pout << mark->getProbeID() << "_0" << " " << mark->getProbeID() << "_0";
					}
				}
			}
			else{
				first_marker = false;
				if(samp->getAbone(loc) != -1){
					if(mark->getAllele(samp->getAbone(loc)) < mark->getAllele(samp->getAbtwo(loc))){
						if(!options.get_allele_as_snp()){
							pout << map_allele(mark, mark->getAllele(samp->getAbone(loc)),&options) << " " << map_allele(mark, mark->getAllele(samp->getAbtwo(loc)),&options);
						}
						else{
							pout << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele(samp->getAbone(loc)), &options) << " " << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele(samp->getAbtwo(loc)), &options);
						}
					}
					else{
						if(!options.get_allele_as_snp()){
							pout << map_allele(mark, mark->getAllele(samp->getAbtwo(loc)),&options) << " " << map_allele(mark, mark->getAllele(samp->getAbone(loc)),&options);
						}
						else{
							pout << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele(samp->getAbtwo(loc)), &options) << " " << mark->getProbeID() << "_" << map_allele(mark, mark->getAllele(samp->getAbone(loc)), &options);
						}
					}
				}
				else{
					pout << "0 0";
				}
			}
		}
		first = false;
		pout << "\n";
	}
	if(options.doDummyMissingParents() || options.doDummyIncompleteParentIds()){

		map<string, string>::iterator iter;
		for(iter = dummy.begin(); iter != dummy.end(); iter++){
			pout << iter->first << iter->second;

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
            				mark->setFlag(true);
							continue;
            			}
            			prev_base = mark->getBPLOC();
            			prev_chrom = mark->getChrom();
            		}
            	}
				pout << "0 0 ";
			}
			pout << "\n";
		}
	}


	pout.close();
	mout.close();
	mdout.close();
}


int PEDOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

/*string PEDOutput::map_allele(string a){
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
