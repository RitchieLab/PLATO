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
*File: PDT2Output.cc
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
#include "AlleleFrequency.h"
#include "PDT2Output.h"
#include "General.h"
#include "Helper.h"

namespace Methods{
void PDT2Output::FilterSummary(){
}

void PDT2Output::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}
}

void PDT2Output::filter(){
}

void PDT2Output::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
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

////   	int ssize = samples->size();
	int fsize = families->size();
	int msize = markers->size();
	string fname1 = opts::_OUTPREFIX_ + "input_pdt2" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".txt";
	}
	string fname2 = opts::_OUTPREFIX_ + "input_pdt2" + options.getOut() + ".map";
	if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".map";
	}
	string fname3 = opts::_OUTPREFIX_ + "input_pdt2" + options.getOut() + ".dat";//getString<int>(order) + ".dat";
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
	ofstream pdt (fname1.c_str());
	if(!pdt.is_open()){
		opts::printLog("Unable to open "+fname1+" for output!\n");
	   	exit(1);
	}
	ofstream pdt_map(fname2.c_str());
	if(!pdt_map.is_open()){
		opts::printLog("Unable to open "+fname2+" for output!\n");
		exit(1);
	}
	ofstream pdt_dat (fname3.c_str());
	if(!pdt_dat.is_open()){
		opts::printLog("Unable to open "+fname3+" for output!\n");
		exit(1);
	}
////	bool map_done = false;

	bool alldigit = true;

    for(int f = 0; f < fsize; f++){
        Family* fam = (*families)[f];
        if(isAlphaNum(fam->getFamID())){
            alldigit = false;
        }
        else{
            fam->setAlphanumeric(false);
        }
	}

	if(!alldigit){
		opts::printLog("Remapping families to digit format.\n");
	    remapFamsToDigit(families);
	    printFamsToDigit(families, "input_pdt2", options);
	}


	int prev_base = 0;
	int prev_chrom = -1;
	//if file specified do that
//	if(options.doCovarsFile()){
//		vector<string> cov_map = options.getCovarMap();
//		for(int i = 0; i < cov_map.size(); i++){
//			pdt_dat << "C " << cov_map[i] << endl;
//		}
//	}
	//otherwise use global traits
//	else if(options.doCovars()){
//		for(int i = 0; i < opts::cov_loc.size(); i++){
//			pdt_dat << "C " << opts::cov_loc[i] << endl;
//		}
//	}
	//if file specified do that
//	if(options.doTraitsFile()){
//		vector<string> trait_map = options.getTraitMap();
//		for(int i = 0; i < trait_map.size(); i++){
//			pdt_dat << "T " << trait_map[i] << endl;
//		}
//	}
	//otherwise use global traits
//	else if(options.doTraits()){
//		for(int i = 0; i < opts::trait_loc.size(); i++){
//			pdt_dat << "T " << opts::trait_loc[i] << endl;
//		}
//	}
	AlleleFrequency* af = new AlleleFrequency(samples, families);
	af->setOptions(options);
	af->flagSamples();
	for(int m = 0; m < msize; m++){
		Marker* mark = (*markers)[m];
		if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
////			int m_loc = mark->getLoc();
			af->calcOne(mark);
			pdt_dat << "M" << " " << mark->getProbeID() << endl;
			pdt_map << mark->getChrom() << " " << mark->getProbeID() << " " << mark->getBPLOC() << endl;
			if(mark->getNumAlleles() < 3){
				pdt_dat << "F " << af->getAone_freq() << endl;
				pdt_dat << "F " << af->getAtwo_freq() << endl;
			}
			else{
				for(int a = 0; a < mark->getNumAlleles(); a++){
					pdt_dat << "F " << af->getMicroFreq(a) << endl;
				}
			}
//			pdt_map << mark->getChrom() << "\t" << mark->getProbeID() << "\t" << mark->getBPLOC() << endl;
		}
	}
	for(int f = 0; f < fsize; f++){
		Family* fam = (*families)[f];
		if(fam->isEnabled()){
			vector<Sample*>* famsamps = fam->getSamples();
			for(int i = 0; i < (int)famsamps->size(); i++){
				Sample* samp = (*famsamps)[i];
		if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
			int prev_base = 0;
			int prev_chrom = -1;
			pdt << fam->getFamID_digit() << " " << samp->getInd_digit() << " " << samp->getDadID_digit() << " " << samp->getMomID_digit() << " ";
			if(samp->getSex()){
				pdt << "1";
			}
			else{
				pdt << "2";
			}
			if(samp->getPheno() > 2){
				pdt << " 0";
			}
			else{
				pdt << " " << samp->getPheno();
			}
			//if file specified do that
/*			map<string, vector<double> >::iterator ctiter;
			if(options.doCovarsFile()){
				map<string, vector<double> > covs = options.getCovarData();
				ctiter = covs.find(samp->getFamID() +"#"+samp->getInd());
				if(ctiter != covs.end()){
					vector<double> data = ctiter->second;
					for(int c = 0; c < data.size(); c++){
						if(data[c] == options.getDefaultCovarMissing()){
							pdt << "\tx";
						}
						else{
							pdt << "\t" << data[c];
						}
					}
				}
				else{
					for(int c = 0; c < options.getCovarMap().size(); c++){
						pdt << "\tx";
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
					for(int c = 0; c < data.size(); c++){
						if(data[c] == options.getDefaultTraitMissing()){
							pdt << "\tx";
						}
						else{
							pdt << "\t" << data[c];
						}
					}
				}
				else{
					for(int c = 0; c < options.getTraitMap().size(); c++){
						pdt << "\tx";
					}
				}
			}
			else if(options.doTraitsName() || options.doTraitsNumber()){
			}
*/
			prev_base = 0;
			prev_chrom = -1;
			for(int m = 0; m < msize; m++){
				Marker* mark = (*markers)[m];
				if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
					int m_loc = mark->getLoc();
					if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doZeroDisabled())){
						pdt << " 0 0";
						continue;
					}
					if(!mark->isMicroSat()){
						if(!samp->getAone(m_loc) && !samp->getAtwo(m_loc)){
							pdt << " " << map_allele(mark, mark->getAllele1(), &options) << " " << map_allele(mark, mark->getAllele1(), &options);
						}
						else if(!samp->getAone(m_loc) && samp->getAtwo(m_loc)){
							if(mark->getAllele1() > mark->getAllele2()){
								pdt << " " << map_allele(mark, mark->getAllele2(), &options) << " " << map_allele(mark, mark->getAllele1(), &options);
							}
							else{
								pdt << " " << map_allele(mark, mark->getAllele1(), &options) << " " << map_allele(mark, mark->getAllele2(), &options);
							}
						}
						else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && !samp->getAmissing(m_loc)){
							pdt << " " << map_allele(mark, mark->getAllele2(), &options) << " " << map_allele(mark, mark->getAllele2(), &options);
						}
						else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && samp->getAmissing(m_loc)){
							pdt << " 0 0";
						}
					}
					else{
						if(samp->getAbone(m_loc) != -1){
							if(mark->getAllele(samp->getAbone(m_loc)) < mark->getAllele(samp->getAbtwo(m_loc))){
								pdt << " " << map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options) << " " << map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options);
							}
							else{
								pdt << " " << map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options) << " " << map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options);
							}
						}
						else{
							pdt << " 0 0";
						}
					}
				}
			}
			pdt << endl;
		}
	}
		}
	}
	pdt.close();
	pdt_dat.close();
	pdt_map.close();
}

int PDT2Output::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}
}
