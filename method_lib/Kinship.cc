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
//#include "Chrom.h"
#include "General.h"
//#include "ChiSquare.h"
#include "Helper.h"
#include "cdflib.h"

using namespace std;
namespace Methods{
string Kinship::stepname = "Kinship";

void Kinship::PrintSummary(){
	string fname1 = opts::_OUTPREFIX_ + "tdt" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	ofstream output (fname1.c_str(), ios::out);
	if(!output){
		opts::printLog("Error opening " + fname1 + ". Exiting!\n");
		exit(1);
	}
	opts::addFile("Marker", stepname, fname1);
	output.precision(4);
	output << "Chrom\trsID\tProbeID\tbploc";
	if((*markers)[0]->getDetailHeaders().size() > 0){
		output << "\t" << (*markers)[0]->getDetailHeaders();
	}
	output << "\tChi_square\tTDT_pvalue\tTDT_neglog(pvalue)\tNum_Fams\tT:U\tA1:A2\tOR\tL" + getString<double>(options.getCI()*100) + "\t" + "U" + getString<double>(options.getCI()*100) << endl;
	opts::addHeader(fname1, "Chi_square");
	opts::addHeader(fname1, "TDT_pvalue");
	opts::addHeader(fname1, "TDT_neglog(pvalue)");
	opts::addHeader(fname1, "Num_Fams");
	opts::addHeader(fname1, "T:U");
	opts::addHeader(fname1, "A1:A2");
	opts::addHeader(fname1, "OR");
	opts::addHeader(fname1, "L" + getString<double>(options.getCI()*100));
	opts::addHeader(fname1, "U" + getString<double>(options.getCI()*100));

	int msize = markers->size();

	double zt = ltqnorm(1 - (1 - options.getCI()) / 2);

	for(int i = 0; i < msize; i++){
		if((*markers)[i]->isEnabled() && !(*markers)[i]->isFlagged()){
            if(options.doChrom()){
                if(!options.checkChrom((*markers)[i]->getChrom())){
                    continue;
                }
                if(!options.checkBp((*markers)[i]->getBPLOC())){
                    continue;
                }
            }
			output << (*markers)[i]->toString() << "\t";
				//<< ((float)(*markers)[i]->getBPLOC()/1000000.0f) << "\t"
			output << chi[i] << "\t"
				<< pval[i] << "\t"
				<< (double)abs(log10(pval[i])) << "\t";
			//output.precision(4);
			output << fams_used[i] << "\t"
				<< trans[i] << ":" << untrans[i] << "\t";
			if(!(*markers)[i]->isMicroSat()){
				output << (*markers)[i]->getAllele1() << ":" << (*markers)[i]->getAllele2();
			}
			else{
				output << "NA";
			}
			double OR = (double)trans[i] / (double)untrans[i];
			output << "\t" << OR;
			double OR_lower = exp(log(OR) - zt * sqrt(1/trans[i] + 1/untrans[i]));
			double OR_upper = exp(log(OR) + zt * sqrt(1/trans[i] + 1/untrans[i]));
			output << "\t" << OR_lower << "\t" << OR_upper;

				output << endl;
			//output.precision(100);
		}
		(*markers)[i]->setFlag(false);
	}

	if(output.is_open()){
		output.close();
	}
}

void Kinship::filter(){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int size = markers->size();
		for(int i = 0; i < size; i++){
			if((*markers)[i]->isEnabled() && !(*markers)[i]->isFlagged()){
	            if(options.doChrom()){
	                if(!options.checkChrom((*markers)[i]->getChrom())){
	                    continue;
		            }
		            if(!options.checkBp((*markers)[i]->getBPLOC())){
			            continue;
		            }
		        }

				bool inc = false;
				if(options.doThreshMarkersHigh() && dGreater(pval[i], options.getThreshMarkersHigh())){
					(*markers)[i]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersLow() && dLess(pval[i],options.getThreshMarkersLow())){
					(*markers)[i]->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}
		}
	}
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




void Kinship::check(Sample* dad, Sample* mom, Family* fam){//int *x, int *y, int fam){
	if(dad->getGeneration() >= mom->getGeneration()){
		return;
	}
	else{

	}

//	int z;

//    if (Global[fam].ped[*x].generation>=Global[fam].ped[*y].generation)
//    	return;
//    else{
//        z = *x;
//        *x = *y;
//        *y = z;
//    }

}

double Kinship::phi2(Sample* dad, Sample* mom, Family* fam){//int i, int j, int fam){
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

//	if (i * j == 0)
//	  {         return 0.0;
//	  }
//	  else if (i == j)
//	  {    return ((1 + phi2(Global[fam].ped[i].father, Global[fam].ped[i].mother,fam)) / 2);

//	  }
//	  else {
//	      check(&i, &j,fam);
//	    return ((phi2(Global[fam].ped[i].father, j,fam) + phi2(Global[fam].ped[i].mother, j,fam)) / 2);

//	  }

}

void Kinship::create_generation(Family* fam){//int fam, struct info *Global){
//    int i=0,per=0,mer=0,generation=0,pas=0,ind=0;
//    map<int,parents>::iterator iter,iterend;

    vector<Sample*>* samples = fam->getSamples();

//    iter=Global[fam].ped.begin();
//    iterend=Global[fam].ped.end();

//    while (iter!=iterend)
    for(unsigned int i = 0; i < samples->size(); i++){
    	Sample* samp = (*samples)[i];
    	Sample* dad = samp->getDad();
    	Sample* mom = samp->getMom();

    	if(samp->getGeneration() == 0){
    		samp->setGeneration(make_generation(samp, dad, mom, fam));
    	}

//    ind=iter->first;
 //   per=iter->second.father;
 //   mer=iter->second.mother;

//    if (iter->second.generation==0)
//    {
//        iter->second.generation=make_generation(ind,per,mer,Global[fam]);
//    }
//    ++iter;
    }

}

int Kinship::make_generation(Sample* samp, Sample* dad, Sample* mom, Family* fam){//int ind, int per, int mer, struct info G){
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


	////////////
//	if (per!=0 && G.ped[per].generation==0)
//    G.ped[per].generation=make_generation(per,G.ped[per].father,G.ped[per].mother,G);
//    if (mer!=0 && G.ped[mer].generation==0)
//    G.ped[mer].generation=make_generation(mer,G.ped[mer].father,G.ped[mer].mother,G);

//    if (per==0 && mer==0)
//        return(1);

//    else {
//    if (G.ped[per].generation>=G.ped[mer].generation)
//        return(G.ped[per].generation+1);

//    else return(G.ped[mer].generation+1);
//    }

}

void Kinship::calculate(Family* family){
	coefficients.clear();

	create_generation(family);

	vector<Sample*>* samples = family->getSamples();
	for(unsigned int i = 0; i < samples->size(); i++){
		coefficients[getString<int>(i) + " " + getString<int>(i)] = phi2((*samples)[i]->getDad(), (*samples)[i]->getMom(), family);
		for(unsigned int j = i + 1; j < samples->size(); j++){
			coefficients[getString<int>(i) + " " + getString<int>(j)] = phi2((*samples)[i], (*samples)[j], family);
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
		if((*markers)[m]->isEnabled()){
            if(options.doChrom()){
                if(!options.checkChrom((*markers)[m]->getChrom())){
                    continue;
                }
                if(!options.checkBp((*markers)[m]->getBPLOC())){
                    continue;
                }
            }
            if(options.doBpSpace()){
	            if(prev_base == 0){
	                prev_base = (*markers)[m]->getBPLOC();
                    prev_chrom = (*markers)[m]->getChrom();
                }
                else{
 	               if((*markers)[m]->getChrom() == prev_chrom && (((*markers)[m]->getBPLOC() - prev_base) < options.getBpSpace())){
	               	   (*markers)[m]->setFlag(true);
				   	   continue;
                   }
                   prev_base = (*markers)[m]->getBPLOC();
                   prev_chrom = (*markers)[m]->getChrom();
                }
            }
			if(!(*markers)[m]->isMicroSat()){
				int mloc = (*markers)[m]->getLoc();


				//start new
				int fams = 0;
				double trans1 = 0;
				double trans2 = 0;
				for(int i = 0; i < ssize; i++){
					Sample* samp = (*samples)[i];
					if(samp != NULL && samp->getDad() != NULL && samp->getMom() != NULL && samp->getDadID() != "0" && samp->getMomID() != "0" && samp->isEnabled() && samp->getDad()->isEnabled() && samp->getMom()->isEnabled()){
						Sample* dad = samp->getDad();
						Sample* mom = samp->getMom();

						int trA = 0; //transmitted allele from first het parent
						int unA = 0; //untransmitted allele from first het parent
						int trB = 0; //transmitted allele from second het parent
						int unB = 0; //untransmitted allele from second het parent

						if(!(*markers)[m]->isMicroSat()){
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
				double pvalue, p, bound, df = 1;
				pvalue = -1;
				int code = 1, status;
				if(tdt_chisq > -1){
					cdfchi(&code, &p, &pvalue, &tdt_chisq, &df, &status, &bound);
				}

				pval[m] = pvalue;
				chi[m] = tdt_chisq;
				fams_used[m] = fams;
				maf[m] = -1;
				trans[m] = trans1;
				untrans[m] = trans2;
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
				Marker* mark = (*markers)[m];
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
