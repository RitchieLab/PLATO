/**********************************************************************************
*                       Structure Output Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates STRUCTURE input files.
*
*
* Files generated:
*
*File: STRUCTOutput.cc
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
#include <iomanip>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include "STRUCTOutput.h"
#include "General.h"
#include "Helper.h"

namespace Methods{
void STRUCTOutput::FilterSummary(){
}

void STRUCTOutput::PrintSummary(){
	int msize = markers->size();

	for(int m = 0; m < msize; m++){
		(*markers)[m]->setFlag(false);
	}

}

void STRUCTOutput::filter(){
}

void STRUCTOutput::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
	orig_num_markers = markers->size();
	orig_num_families = families->size();
 	orig_num_individuals = samples->size();

   	int ssize = samples->size();
//	int msize = markers->size();
	options.setDoAllele1234(true);
	vector<Marker*> good_markers = findValidMarkers(markers, options);
	vector<vector<Marker*> > marker_sets;
	if(options.doRandomMarkers() || options.doRandomRepeat()){
		vector<Marker*> used_markers;
		for(int ss = 0; ss < options.getSets(); ss++){
			vector<Marker*> temp = findRandomMarkers(good_markers, &used_markers, &options);
			if(temp.size() > 0){
				marker_sets.push_back(temp);
			}
		}
	}
	if(marker_sets.size() == 0){
		marker_sets.push_back(good_markers);
	}

	map<string, int> pop;
	/*if(options.getStratFile().length() > 0){
		ifstream pop_f;
		pop_f.open(options.getStratFile().c_str(), ios::in);
		if(!pop_f){
			opts::printLog("Error opening structure statification file: " + options.getStratFile() + "\n");
			throw MethodException("Error opening structure stratification file: " + options.getStratFile() + "\n");
		}
		int line_num = 0;
		while(!pop_f.eof()){
			char buf[256];
			pop_f.getline(buf, 256, '\n');
			line_num++;
			string line = buf;
			if(line == ""){
				continue;
			}

			vector<string> tokens = General::ParseDelimitedLine(line);

			if(tokens.size() != 3){
				opts::printLog("Line " + getString<int>(line_num) + " has incorrect number of elements.  Required number is 3 (ie: Family Individual 1)\n");
				throw MethodException("Line " + getString<int>(line_num) + " has incorrect number of elements.  Required number is 3 (ie: Family Individual 1)\n");
			}
			pop[tokens[0] + " " + tokens[1]] = atoi(tokens[2].c_str());
		}
	}
	*/
	pop = options.getStratificationMap();
	if(options.getStratFile().length() == 0 || pop.size() == 0){
		opts::printLog("Stratification file is required for Structure output creation.  Use option -strat <filename>.\n");
		throw MethodException("Stratification file is required for Structure output creation.  Use option -strat <filename>.\n");
	}

	if(!options.doParentsOnly()){
		opts::printLog("Using ALL samples in STRUCTURE input file.  WARNING: STRUCTURE may produce strange results due to unexpected relationships since all samples are being used as input.  If you experience strange results, please rerun this step to generate a new STRUCTURE input file by using the -parents-only flag in order to use only unrelated parents.\n");
	}
	for(int k = 0; k < (int)marker_sets.size(); k++){
		vector<Marker*> marks = marker_sets[k];
		int msize = marks.size();

		string fname1 = opts::_OUTPREFIX_ + "input_struct_" + getString<int>(k + 1) + options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if(options.getOverrideOut().size() > 0){
			fname1 = options.getOverrideOut() + "_" + getString<int>(k + 1) + ".txt";
		}
		string fname2 = opts::_OUTPREFIX_ + "input_struct_" + getString<int>(k + 1) + options.getOut() + ".map";//getString<int>(order) + ".map";
		if(options.getOverrideOut().size() > 0){
			fname2 = options.getOverrideOut() + "_" + getString<int>(k + 1) + ".map";
		}
		if(!overwrite){
			fname1 += "." + getString<int>(order);
			fname2 += "." + getString<int>(order);
		}
		filenames.push_back(fname1);
		filenames.push_back(fname2);
		ofstream str (fname1.c_str());
		if(!str.is_open()){
			opts::printLog("Unable to open "+fname1+" for output!\n");
		   	throw MethodException("Unable to open " + fname1+ " for output!");
		}
		ofstream str_map (fname2.c_str());
		if(!str_map.is_open()){
			opts::printLog("Unable to open "+fname2+" for output!\n");
		   	throw MethodException("Unable to open " + fname2+ " for output!");
		}
		str_map << "Chrom\tBPLOC\tRsid\tA1/A2" << endl;
		bool map_done = false;
		for(int i = 0; i < ssize; i++){
			Sample* samp = (*samples)[i];
			if((samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())) && ((samp->getDadID() == "0" && samp->getMomID() == "0" && options.doParentsOnly()) || (!options.doParentsOnly()))){
				str << samp->getFamID() << "_" << samp->getInd() << "\t" << pop[samp->getFamID() + " " + samp->getInd()];
//				int prev_base = 0;
//				int prev_chrom = -1;
				for(int m = 0; m < msize; m++){
					Marker* mark = (marks)[m];
					if(mark->isEnabled()){
//			            if(options.doChrom()){
//			                if(!options.checkChrom(mark->getChrom())){
//			                    continue;
//			                }
 //			                if(!options.checkBp(mark->getBPLOC())){
//			                    continue;
//			                }
//			            }
//						if(options.doBpSpace()){
//							if(prev_base == 0){
//								prev_base = mark->getBPLOC();
//								prev_chrom = mark->getChrom();
//							}
//							else{
//								if(mark->getChrom() == prev_chrom && ((mark->getBPLOC() - prev_base) < options.getBpSpace())){
//									mark->setFlag(true);
//									continue;
//								}
//								prev_base = mark->getBPLOC();
//								prev_chrom = mark->getChrom();
//							}
//						}
						int m_loc = mark->getLoc();
						if(!map_done){
							str_map << mark->getChrom() << "\t" << mark->getBPLOC() << "\t" << mark->getRSID() << "\t" << map_allele(mark, mark->getAllele1(), &options) << "/" << map_allele(mark, mark->getAllele2(), &options) << endl;
						}
						if((samp->isExcluded() && options.doIncExcludedSamples() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doZeroDisabled())){
							str << "\t0 0";
							continue;
						}
						if(!mark->isMicroSat()){
							if(!samp->getAone(m_loc) && !samp->getAtwo(m_loc)){
								str << "\t" << map_allele(mark, mark->getAllele1(), &options) << " " << map_allele(mark, mark->getAllele1(), &options);
							}
							else if(!samp->getAone(m_loc) && samp->getAtwo(m_loc)){
								if(mark->getAllele1() > mark->getAllele2()){
									str << "\t" << map_allele(mark, mark->getAllele2(), &options) << " " << map_allele(mark, mark->getAllele1(), &options);
								}
								else{
									str << "\t" << map_allele(mark, mark->getAllele1(), &options) << " " << map_allele(mark, mark->getAllele2(), &options);
								}
							}
							else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && !samp->getAmissing(m_loc)){
								str << "\t" << map_allele(mark, mark->getAllele2(), &options) << " " << map_allele(mark, mark->getAllele2(), &options);
							}
							else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && samp->getAmissing(m_loc)){
								str << "\t0 0";
							}
						}
						else{
							if(samp->getAbone(m_loc) != -1){
								if(mark->getAllele(samp->getAbone(m_loc)) < mark->getAllele(samp->getAbtwo(m_loc))){
									str << "\t" << map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options) << " " << map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options);
								}
								else{
									str << "\t" << map_allele(mark, mark->getAllele(samp->getAbtwo(m_loc)), &options) << " " << map_allele(mark, mark->getAllele(samp->getAbone(m_loc)), &options);
								}
							}
							else{
								str << "\t0 0";
							}
						}
					}
				}
				map_done = true;
				str << endl;
			}
		}
		str.close();
		str_map.close();
	}
}


int STRUCTOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

}
