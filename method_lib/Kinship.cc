#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <unistd.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#include <fenv.h>
#include <algorithm>
#include "Kinship.h"
#include "General.h"
#include "Helpers.h"
#include "cdflib.h"

using namespace std;
namespace Methods{
string Kinship::stepname = "Kinship";

void Kinship::PrintSummary(){

}

void Kinship::filter(){

}

void Kinship::setThreshold(string thresh){
	options.setUp(thresh);
}

void Kinship::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
	    getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
    opts::_MARKERS_WORKING_ -= orig_num_markers;
}




void Kinship::check(Sample* dad, Sample* mom, Family* fam){
	if(dad->getGeneration() >= mom->getGeneration()){
		return;
	}
	else{

	}

}

double Kinship::phi2(Sample* dad, Sample* mom, Family* fam){
	if((dad == NULL || dad->getInd() == "0") || (mom == NULL || mom->getInd() == "0")){
		return 0.0;
	}
	else if(dad->getInd() == mom->getInd()){
		return ((1 + phi2(dad->getDad(), dad->getMom(), fam)) / (double) 2);
	}
	else{
		check(dad, mom, fam);
		return (((phi2(dad->getDad(), mom, fam)) + phi2(dad->getMom(), mom, fam)) / (double) 2);
	}
}

void Kinship::create_generation(Family* fam){
    vector<Sample*>* samples = fam->getSamples();

    for(unsigned int i = 0; i < samples->size(); i++){
    	Sample* samp = (*samples).at(i);
    	Sample* dad = samp->getDad();
    	Sample* mom = samp->getMom();

    	if(samp->getGeneration() == 0){
    		samp->setGeneration(make_generation(samp, dad, mom, fam));
    	}
    }

}

int Kinship::make_generation(Sample* samp, Sample* dad, Sample* mom, Family* fam){
	if(dad != NULL && dad->getInd() != "0" && dad->getGeneration() == 0){
		dad->setGeneration(make_generation(dad, dad->getDad(), dad->getMom(), fam));
	}
	if(mom != NULL && mom->getInd() != "0" && mom->getGeneration() == 0){
		mom->setGeneration(make_generation(mom, mom->getDad(), mom->getMom(), fam));
	}

	if((dad == NULL || dad->getInd() == "0") && (mom == NULL || mom->getInd() == "0")){
		return 1;
	}
	else{
		if(dad == NULL && mom != NULL && mom->getInd() != "0"){
			return (mom->getGeneration() + 1);
		}
		if(mom == NULL && dad != NULL && dad->getInd() != "0"){
			return (dad->getGeneration() + 1);
		}
		if(dad->getGeneration() >= mom->getGeneration()){
			return (dad->getGeneration() + 1);
		}
		else{
			return (mom->getGeneration() + 1);
		}
	}
}

void Kinship::calculate(Family* family){
	coefficients.clear();

	create_generation(family);

	vector<Sample*>* samples = family->getSamples();
	for(unsigned int i = 0; i < samples->size(); i++){
		coefficients[getString<int>(i) + " " + getString<int>(i)] = phi2((*samples).at(i)->getDad(), (*samples).at(i)->getMom(), family);
		for(unsigned int j = i + 1; j < samples->size(); j++){
			coefficients[getString<int>(i) + " " + getString<int>(j)] = phi2((*samples).at(i), (*samples).at(j), family);
		}
	}

}

void Kinship::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

	int msize = markers->size();
	chi.resize(msize);
	pval.resize(msize);
	fams_used.resize(msize);
	maf.resize(msize);
	trans.resize(msize);
	untrans.resize(msize);
	int ssize = samples->size();

	map<int, string> people;
	int prev_base = 0;
	int prev_chrom = -1;
	for(int m = 0; m < msize; m++){
		if((*markers).at(m)->isEnabled()){
            if(options.doChrom()){
                if(!options.checkChrom((*markers).at(m)->getChrom())){
                    continue;
                }
                if(!options.checkBp((*markers).at(m)->getBPLOC())){
                    continue;
                }
            }
            if(options.doBpSpace()){
	            if(prev_base == 0){
	                prev_base = (*markers).at(m)->getBPLOC();
                    prev_chrom = (*markers).at(m)->getChrom();
                }
                else{
 	               if((*markers).at(m)->getChrom() == prev_chrom && (((*markers).at(m)->getBPLOC() - prev_base) < options.getBpSpace())){
	               	   (*markers).at(m)->setFlag(true);
				   	   continue;
                   }
                   prev_base = (*markers).at(m)->getBPLOC();
                   prev_chrom = (*markers).at(m)->getChrom();
                }
            }
			if(!(*markers).at(m)->isMicroSat()){
				int mloc = (*markers).at(m)->getLoc();


				//start new
				int fams = 0;
				double trans1 = 0;
				double trans2 = 0;
				for(int i = 0; i < ssize; i++){
					Sample* samp = (*samples).at(i);
					if(samp != NULL && samp->getDad() != NULL && samp->getMom() != NULL && samp->getDadID() != "0" && samp->getMomID() != "0" && samp->isEnabled() && samp->getDad()->isEnabled() && samp->getMom()->isEnabled()){
						Sample* dad = samp->getDad();
						Sample* mom = samp->getMom();

						int trA = 0; //transmitted allele from first het parent
						int unA = 0; //untransmitted allele from first het parent
						int trB = 0; //transmitted allele from second het parent
						int unB = 0; //untransmitted allele from second het parent

						if(!(*markers).at(m)->isMicroSat()){
							bool pat1 = dad->getAone(mloc);
							bool pat2 = dad->getAtwo(mloc);
							bool mat1 = mom->getAone(mloc);
							bool mat2 = mom->getAtwo(mloc);

							if(pat1 == pat2 && mat1 == mat2){ //mono alleleic
								continue;
							}
							if((pat1 && !pat2) || (mat1 && !mat2)){ //no genotype
								continue;
							}

							if(samp->getPheno() != 2){//make sure child is affected
								continue;
							}

							bool kid1 = samp->getAone(mloc);
							bool kid2 = samp->getAtwo(mloc);

							//kid is 0/0 call
							if(kid1 && !kid2){
								continue;
							}

							//00 - homozygous allele 1
							if((!kid1) && (!kid2)){
								if(((!pat1) && pat2) && ((!mat1) && mat2)){
									trA = 1;
									unA = 2;
									trB = 1;
									unB = 2;
								}
								else{
									trA = 1;
									unA = 2;
								}
							}
							else if((!kid1) && kid2){ //01 - heterozygous
								//het dad
								if(pat1 != pat2){
									//het mom
									if(mat1 != mat2){
										trA = 1;
										trB = 2;
										unA = 2;
										unB = 1;
									}
									else if(!mat1){
										trA = 2;
										unA = 1;
									}
									else{
										trA = 1;
										unA = 2;
									}
								}
								else if(!pat1){
									trA = 2;
									unA = 1;
								}
								else{
									trA = 1;
									unA = 2;
								}
							}
							else{ //11 - homozygous allele 2
								if(((!pat1) && pat2) && ((!mat1) && mat2)){ //dad het & mom het
									trA = 2;
									unA = 1;
									trB = 2;
									unB = 1;
								}
								else{
									trA = 2;
									unA = 1;
								}
							}
							//increment transmission counts
							if(trA == 1) trans1++;
							if(trB == 1) trans1++;
							if(trA == 2) trans2++;
							if(trB == 2) trans2++;
							fams++;
						}//end !micro-sat
					}
				}//end sample iteration

				//skipped discordant parent

				double tdt_chisq, par_chisq, com_chisq;
				tdt_chisq = par_chisq = com_chisq = -1;

				if(trans1 + trans2 > 0){
					tdt_chisq = ((trans1 - trans2)*(trans1 - trans2))/(trans1 + trans2);
				}
				double pvalue, df = 1;
				pvalue = -1;
				if(tdt_chisq > -1){
					pvalue = Helpers::p_from_chi(tdt_chisq, df);
				}

				pval.at(m) = pvalue;
				chi.at(m) = tdt_chisq;
				fams_used.at(m) = fams;
				maf.at(m) = -1;
				trans.at(m) = trans1;
				untrans.at(m) = trans2;
			}
			else{
				//calculate pval for micro satellites
				//compares each allele against all others
				//2x2  A = allele, Oth = other alleles
				//-------------------
				//|        |        |
				//| A/A    |  A/Oth |
				//|        |        |
				//-------------------
				//|        |        |
				//| Oth/A  | Oth/Oth|
				//|        |        |
				//-------------------
				//
				//Finds Allele/Other combination with lowest pvalue and chooses that one to output.
				Marker* mark = (*markers).at(m);
				int numalleles = mark->getNumAlleles();
				vector<double> results;
				results.resize(numalleles, 0);
				for(int n = 0; n < numalleles; n++){
				}

			}
		}
	}
}
}
