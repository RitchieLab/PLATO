/**********************************************************************************
*                       QTDT Input Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
*
*
* Files generated:
*
*File: QTDTOutput.cc
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
#include "QTDTOutput.h"
#include "General.h"
#include "Helper.h"

namespace Methods{
void QTDTOutput::FilterSummary(){
}

void QTDTOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}
}

void QTDTOutput::filter(){
}

void QTDTOutput::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
	orig_num_markers = markers->size();
	orig_num_families = families->size();
 	orig_num_individuals = samples->size();

	if(options.doCovarsFile()){
		options.readCovariates(samples);
	}
	if(options.doTraitsFile()){
		options.readTraits(samples);
	}

   	int ssize = samples->size();
	int msize = markers->size();
	string fname1 = opts::_OUTPREFIX_ + "input_qtdt" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".txt";
	}
	string fname2 = opts::_OUTPREFIX_ + "input_qtdt" + options.getOut() + ".map";//getString<int>(order) + ".map";
	if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".map";
	}
	string fname3 = opts::_OUTPREFIX_ + "input_qtdt" + options.getOut() + ".dat";//getString<int>(order) + ".dat";
	if(options.getOverrideOut().size() > 0){
		fname3 = options.getOverrideOut() + ".dat";
	}
	if(!overwrite){
		fname1 += "." + getString<int>(order);
		fname2 += "." + getString<int>(order);
		fname3 += "." + getString<int>(order);
	}
	filenames.push_back(fname1);
	filenames.push_back(fname2);
	filenames.push_back(fname3);
	ofstream qtdt (fname1.c_str());
	if(!qtdt.is_open()){
		opts::printLog("Unable to open "+fname1+" for output!\n");
	   	exit(1);
	}
	ofstream qtdt_map (fname2.c_str());
	if(!qtdt_map.is_open()){
		opts::printLog("Unable to open "+fname2+" for output!\n");
		exit(1);
	}
	ofstream qtdt_dat (fname3.c_str());
	if(!qtdt_dat.is_open()){
		opts::printLog("Unable to open "+fname3+" for output!\n");
		exit(1);
	}
////	bool map_done = false;

	int prev_base = 0;
	int prev_chrom = -1;
	qtdt_dat << "A AFF" << endl;
	//if file specified do that
	if(options.doCovarsFile()){
		vector<string> cov_map = options.getCovarMap();
		for(int i = 0; i < (int)cov_map.size(); i++){
			qtdt_dat << "C " << cov_map[i] << endl;
		}
	}
	//otherwise use global traits
	else if(options.doCovars()){
		for(int i = 0; i < (int)opts::cov_loc.size(); i++){
			qtdt_dat << "C " << opts::cov_loc[i] << endl;
		}
	}
	//if file specified do that
	if(options.doTraitsFile()){
		vector<string> trait_map = options.getTraitMap();
		for(int i = 0; i < (int)trait_map.size(); i++){
			qtdt_dat << "T " << trait_map[i] << endl;
		}
	}
	//otherwise use global traits
	else if(options.doTraits()){
		for(int i = 0; i < (int)opts::trait_loc.size(); i++){
			qtdt_dat << "T " << opts::trait_loc[i] << endl;
		}
	}
	for(int m = 0; m < msize; m++){
		Marker* mark = (*markers)[m];
		if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
////			int m_loc = mark->getLoc();
			qtdt_dat << "M" << " " << mark->getProbeID() << endl;
			qtdt_map << mark->getChrom() << "\t" << mark->getProbeID() << "\t" << mark->getBPLOC() << endl;
		}
	}

	map<string, string> dummy;
	for(int i = 0; i < ssize; i++){
		Sample* samp = (*samples)[i];
		if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
			int valid_parents = 0;
			if(samp->getDad() != NULL){
				valid_parents++;
			}
			if(samp->getMom() != NULL){
				valid_parents++;
			}
			int prev_base = 0;
			int prev_chrom = -1;
			qtdt << samp->getFamID() << "\t" << samp->getInd() << "\t";
			//<< samp->getDadID() << "\t" << samp->getMomID() << "\t";



			if((samp->getDad() == NULL || !samp->getDad()->isEnabled()) && samp->getDadID() != "0" && options.doRemMissingParents()){
				qtdt << "0\t";
			}
			else if(samp->getDad() == NULL && samp->getDadID() != "0" && options.doDummyMissingParents()){
				qtdt << samp->getDadID() << "\t";
				dummy[samp->getFamID() + "\t" + samp->getDadID()] =  "\t0\t0\t1\t0";
			}
			else if(samp->getDad() == NULL && samp->getMom() != NULL && options.doDummyIncompleteParentIds()){
				string did = "DUM999999";
				if(samp->getDadID() != "0"){
					did = samp->getDadID();
				}
				qtdt << did << "\t";
				dummy[samp->getFamID() + "\t" + did] = "\t0\t0\t1\t0";
			}
			else if(options.doZeroIncompleteTrioIds() && valid_parents != 2){
				qtdt << "0\t";
			}
			else{
				qtdt << samp->getDadID() << "\t";
			}
			if((samp->getMom() == NULL || !samp->getMom()->isEnabled()) && samp->getMomID() != "0" && options.doRemMissingParents()){
				qtdt << "0\t";
			}
			else if(samp->getMom() == NULL && samp->getMomID() != "0" && options.doDummyMissingParents()){
				qtdt << samp->getMomID() << "\t";
				dummy[samp->getFamID() + "\t" + samp->getMomID()] = "\t0\t0\t2\t0";
			}
			else if(samp->getMom() == NULL && samp->getDad() != NULL && options.doDummyIncompleteParentIds()){
				string mid = "DUM999998";
				if(samp->getMomID() != "0"){
					mid = samp->getMomID();
				}
				qtdt << mid << "\t";
				dummy[samp->getFamID() + "\t" + mid] = "\t0\t0\t2\t0";
			}
			else if(options.doZeroIncompleteTrioIds() && valid_parents != 2){
				qtdt << "0\t";
			}
			else{
			    qtdt << samp->getMomID() << "\t";
			}



			if(samp->getSex()){
				qtdt << "1";
			}
			else{
				qtdt << "2";
			}
			if(samp->getPheno() > 2){
				qtdt << " 0";
			}
			else{
				qtdt << " " << samp->getPheno();
			}
			//if file specified do that
			map<string, vector<double> >::iterator ctiter;
			if(options.doCovarsFile()){
				map<string, vector<double> > covs = options.getCovarData();
				ctiter = covs.find(samp->getFamID() +"#"+samp->getInd());
				if(ctiter != covs.end()){
					vector<double> data = ctiter->second;
					for(int c = 0; c < (int)data.size(); c++){
						if(data[c] == options.getDefaultCovarMissing()){
							qtdt << "\tx";
						}
						else{
							qtdt << "\t" << data[c];
						}
					}
				}
				else{
					for(int c = 0; c < (int)options.getCovarMap().size(); c++){
						qtdt << "\tx";
					}
				}
			}
			else if(options.doCovarsName() || options.doCovarsNumber()){
			}
			if(options.doTraitsFile()){
				map<string, vector<double> > traits = options.getTraitData();
				ctiter = traits.find(samp->getFamID() +"#"+samp->getInd());
				if(ctiter != traits.end()){
					vector<double> data = ctiter->second;
					for(int c = 0; c < (int)data.size(); c++){
						if(data[c] == options.getDefaultTraitMissing()){
							qtdt << "\tx";
						}
						else{
							qtdt << "\t" << data[c];
						}
					}
				}
				else{
					for(int c = 0; c < (int)options.getTraitMap().size(); c++){
						qtdt << "\tx";
					}
				}
			}
			else if(options.doTraitsName() || options.doTraitsNumber()){
			}

			prev_base = 0;
			prev_chrom = -1;
			for(int m = 0; m < msize; m++){
				Marker* mark = (*markers)[m];
				if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
					int m_loc = mark->getLoc();
					if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doZeroDisabled())){
						qtdt << "\t0/0";
						continue;
					}
					if(!mark->isMicroSat()){
						if(!samp->getAone(m_loc) && !samp->getAtwo(m_loc)){
							qtdt << "\t" << map_allele(mark, mark->getAllele1(), &options) << "/" << map_allele(mark, mark->getAllele1(), &options);
						}
						else if(!samp->getAone(m_loc) && samp->getAtwo(m_loc)){
							if(mark->getAllele1() > mark->getAllele2()){
								qtdt << "\t" << map_allele(mark, mark->getAllele2(), &options) << "/" << map_allele(mark, mark->getAllele1(), &options);
							}
							else{
								qtdt << "\t" << map_allele(mark, mark->getAllele1(), &options) << "/" << map_allele(mark, mark->getAllele2(), &options);
							}
						}
						else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && !samp->getAmissing(m_loc)){
							qtdt << "\t" << map_allele(mark, mark->getAllele2(), &options) << "/" << map_allele(mark, mark->getAllele2(), &options);
						}
						else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && samp->getAmissing(m_loc)){
							qtdt << "\t0/0";
						}
					}
					else{
						if(samp->getAbone(m_loc) != -1){
							if(mark->getAllele(samp->getAbone(m_loc)) < mark->getAllele(samp->getAbtwo(m_loc))){
								qtdt << "\t" << map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options) << "/" << map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options);
							}
							else{
								qtdt << "\t" << map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options) << "/" << map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options);
							}
						}
						else{
							qtdt << "\t0/0";
						}
					}
				}
			}
			qtdt << endl;
		}
	}

	if(options.doDummyMissingParents() || options.doDummyIncompleteParentIds()){

		map<string, string>::iterator iter;
		for(iter = dummy.begin(); iter != dummy.end(); iter++){
			qtdt << iter->first << iter->second;

			if(options.doCovarsFile()){
					for(int c = 0; c < (int)options.getCovarMap().size(); c++){
						qtdt << "\tx";
					}
			}
			if(options.doTraitsFile()){
					for(int c = 0; c < (int)options.getTraitMap().size(); c++){
						qtdt << "\tx";
					}
			}

			int prev_base = 0;
			int prev_chrom = -1;
			for(int i = 0; i < msize; i++){
			//Marker* mark = (*markers)[mloc];
				Marker* mark = (*markers)[i];
				if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
					qtdt << "\t0/0";
				}
			}
			qtdt << "\n";
		}
	}



	qtdt.close();
	qtdt_map.close();
	qtdt_dat.close();
}

int QTDTOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

}
