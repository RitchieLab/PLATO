/**********************************************************************************
*                       Superlink Input Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates Superlink input files
*
*
*
*File: SuperlinkOutput.cc
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
#include "SuperlinkOutput.h"
#include "General.h"
#include "Helper.h"
namespace Methods{
string SuperlinkOutput::stepname = "output-superlink";

void SuperlinkOutput::FilterSummary(){
}

void SuperlinkOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}

}

void SuperlinkOutput::filter(){
}

void SuperlinkOutput::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

////   	int ssize = samples->size();
	int msize = markers->size();
	int fsize = families->size();

////	int prev_base = 0;
////	int prev_chrom = -1;
	int numgoodmarkers = 0;
	vector<Marker*> good_markers = findValidMarkers(markers, options);
	numgoodmarkers = good_markers.size();
	msize = good_markers.size();
	int numfams = 0;
	bool alldigit = true;
	for(int f = 0; f < fsize; f++){
		Family* fam = (*families)[f];
		if(isAlphaNum(fam->getFamID())){
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

//	if(!alldigit){
		opts::printLog("Remapping families to digit format.\n");
		remapFamsToDigit(families);
		printFamsToDigit(families, "input_superlink", options);
//	}

////	int parents = 0;
////	int stotal = 0;

	string fname1 = opts::_OUTPREFIX_ + "input_superlink_locus" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".locus";
	}
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	filenames.push_back(fname1);
	ofstream locus (fname1.c_str());
	if(!locus.is_open()){
		opts::printLog("Unable to open " + fname1 + " for output!\n");
		throw MethodException("Unable to open " + fname1 + " for output!\n");
	}
	string fname2 = opts::_OUTPREFIX_ + "input_superlink_ped" + options.getOut() + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".ped";
	}
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	filenames.push_back(fname2);
	ofstream ped (fname2.c_str());

    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
	string date = asctime(timeinfo);
	date = date.erase((date.size() - 1), 1);
	//#loci, 0, sexlinked, #4=superlinkmap;5=supermlink,#complex affection locs (0/2)
	locus << (numgoodmarkers + 1) << " 0 0 5 0" << endl;
	//filler
	locus << "0 0.0 0.0 0" << endl;

	//print order of markers
	for(int i = 0; i < (numgoodmarkers + 1); i++){
		locus << (i+1) << " ";
	}
	locus << endl;

	locus << "###########\n#Insert affection status locus information here" << endl << "#Including penetrance information\n############" << endl;


	AlleleFrequency* af = new AlleleFrequency(samples, families);
	af->setOptions(options);
	af->flagSamples();
	for(int m = 0; m < msize; m++){
		Marker* mark = good_markers[m];//(*markers)[m];
		if(mark->isEnabled()){// && !mark->isFlagged()){
			locus << "3 " << mark->getNumAlleles() << " #" << mark->getProbeID() << "#" << endl;
			af->calcOne(mark);
			if(mark->getNumAlleles() < 3){
				locus << af->getAone_freq() << " " << af->getAtwo_freq() << endl;
			}
			else{
				locus << af->getMicroFreq(0);
				for(int a = 1; a < mark->getNumAlleles(); a++){
					locus << " " << af->getMicroFreq(a);
				}
				locus << endl;
			}
		}
	}
	delete(af);
	locus.close();

	for(int f = 0; f < fsize; f++){
		Family* fam = (*families)[f];
		if(fam->isEnabled() || (fam->isExcluded() && options.doIncExcludedSamples()) || (!fam->isEnabled() && options.doIncDisabledSamples())){
////			int toout = 0;
////			int loopcount = 0;
			vector<Sample*>* fsamps = fam->getSamples();
			if(fsamps->size() == 0){
				continue;
			}
			int fssize = fsamps->size();
			for(int s = 0; s < fssize; s++){
				Sample* samp = (*fsamps)[s];
				if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
					ped << fam->getFamID_digit() << "\t" << samp->getInd_digit() << " " << samp->getDadID_digit() << " " << samp->getMomID_digit() << " ";
					vector<Sample*>* children = samp->getChildren();
					if(children->size() == 0){
						ped << "0 ";
					}
					else{
						ped << (*children)[0]->getInd_digit() << " ";
					}
					//paternal sib
					Sample* sib = samp->getPatSib();
					if(sib == NULL){
						ped << "0 ";
					}
					else{
						ped << sib->getInd_digit() << " ";
					}
					//maternal sib
					sib = samp->getMatSib();
					if(sib == NULL){
						ped << "0 ";
					}
					else{
						ped << sib->getInd_digit() << " ";
					}

					//gender
					if(samp->getSex()){
						ped << "1";
					}
					else{
						ped << "2";
					}
					//filler
					ped << " 0";
					//Disease status
					if(samp->getPheno() == 2){
						ped << " " << "2";
					}
					else if(samp->getPheno() == 1){
						ped << " " << "1";
					}
					else{
						ped << " " << "0";
					}

					//penetrance info
					ped << " " << options.findPenetranceCode(samp->getFamID() + " " + samp->getInd());

					for(int m = 0; m < msize; m++){
						Marker* mark = good_markers[m];//(*markers)[m];
						if(mark->isEnabled()){// && !mark->isFlagged()){
							int mloc = mark->getLoc();
							if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doZeroDisabled())){
								ped << " 0 0";
								continue;
							}
							if(!mark->isMicroSat()){
								if(!samp->getAone(mloc) && !samp->getAtwo(mloc)){
									ped << " " << map_allele(mark, mark->getAllele1(), &options);
									ped << " " << map_allele(mark, mark->getAllele1(), &options);
								}
								else if(!samp->getAone(mloc) && samp->getAtwo(mloc)){
									ped << " " << map_allele(mark, mark->getAllele1(), &options);
									ped << " " << map_allele(mark, mark->getAllele2(), &options);
								}
								else if(samp->getAone(mloc) && samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
									ped << " " << map_allele(mark, mark->getAllele2(), &options);
									ped << " " << map_allele(mark, mark->getAllele2(), &options);
								}
								else{
									ped << " 0 0";
								}
							}//end !microsat
							else{
								if(samp->getAbone(mloc) == -1){
									ped << " 0 0";
								}
								else{
									ped << " " << map_allele(mark, mark->getAllele(samp->getAbone(mloc)), &options);
									ped << " " << map_allele(mark, mark->getAllele(samp->getAbtwo(mloc)), &options);
								}
							}//end microsat
						}
					}
					ped << endl;
				}
			}
		}
	}
	ped.close();


}


int SuperlinkOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}
}
