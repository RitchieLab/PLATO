/**********************************************************************************
*                       Lapis Input Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates Lapis input files
*
*
*
*File: LAPISOutput.cc
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
#include <time.h>
#include "LAPISOutput.h"
#include "General.h"
#include "Helpers.h"
namespace Methods{
string LAPISOutput::stepname = "output-lapis";

void LAPISOutput::FilterSummary(){
}

void LAPISOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}

}

void LAPISOutput::filter(){
}

void LAPISOutput::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

  //// 	int ssize = samples->size();
	int msize = markers->size();
	int fsize = families->size();

////	int prev_base = 0;
////	int prev_chrom = -1;
	int numgoodmarkers = 0;
	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
//	for(int i = 0; i < msize; i++){
//		Marker* mark = (*markers)[i];
//		if(mark->isEnabled()){
//			if(options.doChrom()){
//				if(!options.checkChrom(mark->getChrom())){
//					mark->setFlag(true);
//				    continue;
//			    }
//			    if(!options.checkBp(mark->getBPLOC())){
//					mark->setFlag(true);
//				    continue;
//			    }
//			}
//            if(options.doBpSpace()){
//	            if(prev_base == 0){
//	                prev_base = mark->getBPLOC();
  //                  prev_chrom = mark->getChrom();
//                }
//                else{
//	                if(mark->getChrom() == prev_chrom && ((mark->getBPLOC() - prev_base) < options.getBpSpace())){
//		    			mark->setFlag(true);
//						continue;
//	                }
//	                prev_base = mark->getBPLOC();
//	                prev_chrom = mark->getChrom();
 //               }
//			}
//			numgoodmarkers++;
//		}
//	}
	numgoodmarkers = good_markers.size();
	msize = good_markers.size();
	int numfams = 0;
	bool alldigit = true;
	for(int f = 0; f < fsize; f++){
		Family* fam = (*families)[f];
		if(Helpers::isAlphaNum(fam->getFamID())){
			alldigit = false;
		}
		else{
			fam->setAlphanumeric(false);
		}
		vector<Sample*>* fsamps = fam->getSamples();
		int fssize = fsamps->size();
		for(int s = 0; s < fssize; s++){
			Sample* samp = (*fsamps)[s];
			if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
				numfams++;
				break;
			}
		}
	}

	if(!alldigit){
		opts::printLog("Remapping families to digit format.\n");
		Helpers::remapFamsToDigit(families);
		Helpers::printFamsToDigit(families, "input_lapis", options);
	}

////	int parents = 0;
////	int stotal = 0;

	string fname1 = opts::_OUTPREFIX_ + "input_lapis" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
    if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".txt";
	}
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	filenames.push_back(fname1);
	ofstream str (fname1.c_str());
	if(!str.is_open()){
		opts::printLog("Unable to open " + fname1 + " for output!\n");
		//exit(1);
		throw MethodException("Unable to open " + fname1 + " for output!\n");
	}
    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
	string date = asctime(timeinfo);
	date = date.erase((date.size() - 1), 1);
	str << (numgoodmarkers + 1) << " " << numfams << " " << date << " " << fname1 << " " << opts::_PEDFILE_ << " " << opts::_BINPREFIX_ << endl << endl;
	//DEFAULT DISEASE MARKER
	str << " 1   2   TESTDX\n";
	str << "0.9990 0.0010\n";
	str << "  1\n";
	str << "0.0000 0.0001 0.0001\n";
	str << "A\n\n";

	AlleleFrequency* af = new AlleleFrequency(samples, families);
	af->setOptions(options);
	af->flagSamples();
	for(int m = 0; m < msize; m++){
		Marker* mark = good_markers[m];//(*markers)[m];
		if(mark->isEnabled()){// && !mark->isFlagged()){
			str << "3 " << mark->getNumAlleles() << " " << mark->getChrom() << " " << mark->getProbeID() << endl;
			af->calcOne(mark);
			if(mark->getNumAlleles() < 3){
				str << af->getAone_freq() << " " << af->getAtwo_freq() << endl;
				str << Helpers::map_allele(mark, mark->getAllele1(), &options) << " " << Helpers::map_allele(mark, mark->getAllele2(), &options) << endl;
			}
			else{
				str << af->getMicroFreq(0);
				for(int a = 1; a < mark->getNumAlleles(); a++){
					str << " " << af->getMicroFreq(a);
				}
				str << endl;
				str << Helpers::map_allele(mark, mark->getAllele(0), &options);
				for(int a = 1; a < mark->getNumAlleles(); a++){
					str << " " << Helpers::map_allele(mark, mark->getAllele(a), &options);
				}
				str << endl;
			}
			str << endl;
		}
	}
	delete(af);

	for(int f = 0; f < fsize; f++){
		Family* fam = (*families)[f];
		if(fam->isEnabled() || (fam->isExcluded() && options.doIncExcludedSamples()) || (!fam->isEnabled() && options.doIncDisabledSamples())){
			int toout = 0;
			int loopcount = 0;
			vector<Sample*>* fsamps = fam->getSamples();
			string type = "GROUP";
			if(fsamps->size() > 1){
				type = "FAMILY";
			}
			if(fsamps->size() == 0){
				continue;
			}
			int fssize = fsamps->size();
			for(int s = 0; s < fssize; s++){
				Sample* samp = (*fsamps)[s];
				if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
					toout++;
				}
			}
			str << toout << " " << loopcount << " " << type << " " << fam->getFamID_digit() << " ";
		   	if(options.haveCenterCodes()){
				Sample* samp = (*fsamps)[0];
				string cent = options.findCenterCode(samp->getFamID() + " " + samp->getInd());
				if(cent != ""){
					str << cent;
				}
				else{
					str << "CEN";
				}
			}
			else{
				str << "CEN";
			}
			str << endl;//fam->getCenter() << endl;
			for(int s = 0; s < fssize; s++){
				Sample* samp = (*fsamps)[s];
				if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
					str << samp->getInd_digit() << " " << samp->getDadID_digit() << " " << samp->getMomID_digit() << " ";
					if(samp->getSex()){
						str << "M";
					}
					else{
						str << "F";
					}
					if(samp->getPheno() == 2){
						str << " " << "A";
					}
					else if(samp->getPheno() == 1){
						str << " " << "N";
					}
					else{
						str << " " << "U";
					}
//					str << " " << "U";

					for(int m = 0; m < msize; m++){
						Marker* mark = good_markers[m];//(*markers)[m];
						if(mark->isEnabled()){// && !mark->isFlagged()){
							int mloc = mark->getLoc();
							if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doZeroDisabled())){
								str << " 0 0";
								continue;
							}
							if(!mark->isMicroSat()){
								if(!samp->getAone(mloc) && !samp->getAtwo(mloc)){
									str << " " << Helpers::map_allele(mark, mark->getAllele1(), &options);
									str << " " << Helpers::map_allele(mark, mark->getAllele1(), &options);
								}
								else if(!samp->getAone(mloc) && samp->getAtwo(mloc)){
									str << " " << Helpers::map_allele(mark, mark->getAllele1(), &options);
									str << " " << Helpers::map_allele(mark, mark->getAllele2(), &options);
								}
								else if(samp->getAone(mloc) && samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
									str << " " << Helpers::map_allele(mark, mark->getAllele2(), &options);
									str << " " << Helpers::map_allele(mark, mark->getAllele2(), &options);
								}
								else{
									str << " 0 0";
								}
							}//end !microsat
							else{
								if(samp->getAbone(mloc) == -1){
									str << " 0 0";
								}
								else{
									str << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbone(mloc)), &options);
									str << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(mloc)), &options);
								}
							}//end microsat
						}
					}
					str << endl;
				}
			}
			str << "7000" << endl;
		}
	}
	str.close();

}


int LAPISOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

}
