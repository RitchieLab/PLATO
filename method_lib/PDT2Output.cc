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
#include "AlleleFrequency.h"
#include "PDT2Output.h"
#include "General.h"
#include "Helpers.h"

namespace Methods{
void PDT2Output::FilterSummary(){
}

void PDT2Output::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers).at(i)->setFlag(false);
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

	int fsize = families->size();
	int msize = markers->size();

	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	msize = good_markers.size();

	string fname1 = opts::_OUTPREFIX_ + "input_pdt2" + options.getOut() + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".txt";
	}
	string fname2 = opts::_OUTPREFIX_ + "input_pdt2" + options.getOut() + ".map";
	if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".map";
	}
	string fname3 = opts::_OUTPREFIX_ + "input_pdt2" + options.getOut() + ".dat";
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
		throw MethodException("Unable to open "+fname1+" for output!\n");
	}
	ofstream pdt_map(fname2.c_str());
	if(!pdt_map.is_open()){
		opts::printLog("Unable to open "+fname2+" for output!\n");
		throw MethodException("Unable to open "+fname2+" for output!\n");
	}
	ofstream pdt_dat (fname3.c_str());
	if(!pdt_dat.is_open()){
		opts::printLog("Unable to open "+fname3+" for output!\n");
		throw MethodException("Unable to open "+fname3+" for output!\n");
	}

	bool alldigit = true;

    for(int f = 0; f < fsize; f++){
        Family* fam = (*families).at(f);
        if(Helpers::isAlphaNum(fam->getFamID())){
            alldigit = false;
        }
        else{
            fam->setAlphanumeric(false);
        }
	}

	if(!alldigit){
		opts::printLog("Remapping families to digit format.\n");
	    Helpers::remapFamsToDigit(families);
	    Helpers::printFamsToDigit(families, "input_pdt2", options);
	}

	AlleleFrequency* af = new AlleleFrequency(samples, families);
	af->setOptions(options);
	af->flagSamples();
	for(int m = 0; m < msize; m++){
		Marker* mark = good_markers.at(m);
		if(mark->isEnabled()){
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
		}
	}
	for(int f = 0; f < fsize; f++){
		Family* fam = (*families).at(f);
		if(fam->isEnabled()){
			vector<Sample*>* famsamps = fam->getSamples();
			for(int i = 0; i < (int)famsamps->size(); i++){
				Sample* samp = (*famsamps).at(i);
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
			prev_base = 0;
			prev_chrom = -1;
			for(int m = 0; m < msize; m++){
				Marker* mark = good_markers.at(m);
				if(mark->isEnabled()){
					int m_loc = mark->getLoc();
					if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doZeroDisabled())){
						pdt << " 0 0";
						continue;
					}
					if(!mark->isMicroSat()){
						if(!samp->getAone(m_loc) && !samp->getAtwo(m_loc)){
							pdt << " " << Helpers::map_allele(mark, mark->getAllele1(), &options) << " " << Helpers::map_allele(mark, mark->getAllele1(), &options);
						}
						else if(!samp->getAone(m_loc) && samp->getAtwo(m_loc)){
							if(mark->getAllele1() > mark->getAllele2()){
								pdt << " " << Helpers::map_allele(mark, mark->getAllele2(), &options) << " " << Helpers::map_allele(mark, mark->getAllele1(), &options);
							}
							else{
								pdt << " " << Helpers::map_allele(mark, mark->getAllele1(), &options) << " " << Helpers::map_allele(mark, mark->getAllele2(), &options);
							}
						}
						else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && !samp->getAmissing(m_loc)){
							pdt << " " << Helpers::map_allele(mark, mark->getAllele2(), &options) << " " << Helpers::map_allele(mark, mark->getAllele2(), &options);
						}
						else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && samp->getAmissing(m_loc)){
							pdt << " 0 0";
						}
					}
					else{
						if(samp->getAbone(m_loc) != -1){
							if(mark->getAllele(samp->getAbone(m_loc)) < mark->getAllele(samp->getAbtwo(m_loc))){
								pdt << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options) << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options);
							}
							else{
								pdt << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options) << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options);
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
