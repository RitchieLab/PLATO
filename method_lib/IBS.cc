/**********************************************************************************
*         IBS calculation module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Performs a variety of IBS calculations compairing samples or trios or transmissions
*
*
*
*File: IBS.cc
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
#include "IBS.h"
#include "Options.h"
#include "General.h"
#include "Helpers.h"
namespace Methods{
string IBS::stepname = "ibs";

//DEPRECATED
void IBS::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

//DEPRECATED
void IBS::PrintSummary(){


}

//DEPRECATED
void IBS::filterOne(int m){

}

//DEPRECATED
void IBS::filter(){

}

//calculates transmissions by trio.
//Special coding used developed by Will Bush
//-1 and 1 for heterozygous paternal and maternal transmissions to children
vector<int> IBS::calcTrioTransmission(int f1, Marker* m){
	Marker* mark = m;
	int mloc = mark->getLoc();
	Family* fam = data_set->get_pedigree(f1);
	vector<int> count (2, -9);
	vector<Sample*>* fsamps = fam->getSamples();
	if(mark != NULL && fam != NULL && fsamps->size() == 3){
		Sample* mom = (*fsamps)[0];
		Sample* dad = (*fsamps)[1];
		Sample* child = (*fsamps)[2];
		if(child->getDad() == NULL || child->getMom() == NULL){
			if(mom->getDad() != NULL && mom->getMom() != NULL && mom->getChildren()->size() == 0){
				Sample* temp = mom;
				mom = child;
				child = temp;
			}
			else if(dad->getDad() != NULL && dad->getMom() != NULL && dad->getChildren()->size() == 0){
				Sample* temp = dad;
				dad = child;
				child = temp;
			}
		}
		if(dad->getSex() == 2 && mom->getSex() == 1){
			Sample* temp = dad;
			dad = mom;
			mom = temp;
		}

		bool f1pa1 = dad->getAone(mloc);
		bool f1pa2 = dad->getAtwo(mloc);
		bool f1pa3 = dad->getAmissing(mloc);
		bool f1ma1 = mom->getAone(mloc);
		bool f1ma2 = mom->getAtwo(mloc);
		bool f1ma3 = mom->getAmissing(mloc);
		bool f1kid1 = child->getAone(mloc);
		bool f1kid2 = child->getAtwo(mloc);
		bool f1kid3 = child->getAmissing(mloc);

		//zero geno for fam1?
		if((f1pa1 && f1pa2 && f1pa3) || (f1ma1 && f1ma2 && f1ma3) || (f1kid1 && f1kid2 && f1kid3)){
			return count;
		}

		if(f1pa1 != f1pa2){
			if(f1kid1 != f1kid2 && f1ma1 && f1ma2){
				count[1] = 1;
			}
			else if(f1kid1 && f1kid2 && f1ma1 && f1ma2){
				count[1] = -1;
			}
			else if(f1kid1 && f1kid2 && f1ma1 != f1ma2){
				count[1] = -1;
			}
			else if(!f1kid1 && !f1kid2 && f1ma1 != f1ma2){
				count[1] = 1;
			}
			else if(f1kid1 != f1kid2 && !f1ma1 && !f1ma2){
				count[1] = -1;
			}
			else if(!f1kid1 && !f1kid2 && !f1ma1 && !f1ma2){
				count[1] = 1;
			}
		}
		if(f1ma1 != f1ma2){
			if(f1kid1 != f1kid2 && f1pa1 && f1pa2){
				count[0] = 1;
			}
			else if(f1kid1 && f1kid2 && f1pa1 && f1pa2){
				count[0] = -1;
			}
			else if(f1kid1 && f1kid2 && f1pa1 != f1pa2){
				count[0] = -1;
			}
			else if(!f1kid1 && !f1kid2 && f1pa1 != f1pa2){
				count[0] = 1;
			}
			else if(f1kid1 != f1kid2 && !f1pa1 && !f1pa2){
				count[0] = -1;
			}
			else if(!f1kid1 && !f1kid2 && !f1pa1 && !f1pa2){
				count[0] = 1;
			}
		}
	}
	return count;
}

//calculates IBS based on sample/sample pair and a locus
int IBS::calcPairLocus(int s1, int s2, Marker* m){
	Marker* mark = m;
	Sample* samp1 = data_set->get_sample(s1);
	Sample* samp2 = data_set->get_sample(s2);
	int count = 0;

	if(mark != NULL && samp1 != NULL && samp2 != NULL){
		int mloc = mark->getLoc();
		if(mark->isEnabled() && samp1->isEnabled() && samp2->isEnabled() && !samp1->isExcluded() && !samp2->isExcluded()){
			if((samp1->getAone(mloc) && samp1->getAtwo(mloc) && samp1->getAmissing(mloc))
					|| (samp2->getAone(mloc) && samp2->getAtwo(mloc) && samp2->getAmissing(mloc))){
				return -1;
			}
			if(samp1->getAone(mloc) == samp2->getAone(mloc)){
				count++;
			}
			if(samp1->getAtwo(mloc) == samp2->getAtwo(mloc)){
				count++;
			}
		}
	}
	else{
		return -1;
	}
	return count;
}

//Calculates IBS based on two trios and a locus
vector<double> IBS::calcTriosLocus(int f1, int f2, Marker* m){
	Marker* mark = m;
	Family* fam1 = data_set->get_pedigree(f1);
	Family* fam2 = data_set->get_pedigree(f2);

	vector<double> counts(2, 0);

	if(mark != NULL && fam1 != NULL && fam2 != NULL && fam1->getSamples()->size() == 3 && fam2->getSamples()->size() == 3){
		//confirm only 1 child
		if(fam1->isEnabled() && fam1->getNonFounders()->size() == 1 && fam2->isEnabled() && fam2->getNonFounders()->size() == 1){
			Sample* child1 = fam1->getNonFounders()->at(0);
			Sample* child2 = fam2->getNonFounders()->at(0);
			Sample* dad1 = child1->getDad();
			Sample* dad2 = child2->getDad();
			Sample* mom1 = child1->getMom();
			Sample* mom2 = child2->getMom();

			if(dad1 != NULL && dad2 != NULL && mom1 != NULL && mom2 != NULL){
				if(child1->isEnabled() && child2->isEnabled() && dad1->isEnabled() && dad2->isEnabled() && mom1->isEnabled() && mom2->isEnabled()){
					int mloc = mark->getLoc();
					//family 1
					bool f1pa1 = dad1->getAone(mloc);
					bool f1pa2 = dad1->getAtwo(mloc);
					bool f1pa3 = dad1->getAmissing(mloc);
					bool f1ma1 = mom1->getAone(mloc);
					bool f1ma2 = mom1->getAtwo(mloc);
					bool f1ma3 = mom1->getAmissing(mloc);
					bool f1kid1 = child1->getAone(mloc);
					bool f1kid2 = child1->getAtwo(mloc);
					bool f1kid3 = child1->getAmissing(mloc);
					int f1ptr = 0;
					int f1pun = 0;
					int f1mtr = 0;
					int f1mun = 0;
					bool f1_allhet = false;

					//zero geno for fam1?
					if((f1pa1 && f1pa2 && f1pa3) || (f1ma1 && f1ma2 && f1ma3) || (f1kid1 && f1kid2 && f1kid3)){
						counts[0] = -1;
						counts[1] = -1;
						return counts;
					}

					if(!f1kid1 && !f1kid2){
						if(!f1pa1){
							f1ptr = 1;
							f1pun = 2;
						}
						if(!f1ma1){
							f1mtr = 1;
							f1mun = 2;
						}
					}
					else if(!f1kid1 && f1kid2){
						if(f1pa1 != f1pa2 && f1ma1 != f1ma2){
							f1_allhet = true;
						}
						else{
							if(f1pa1 != f1pa2){
								if(!f1ma1){
									f1ptr = 2;
									f1pun = 1;
								}
								else{
									f1ptr = 1;
									f1pun = 2;
								}
							}
							else{
								if(!f1pa1){
									f1mtr = 2;
									f1mun = 1;
								}
								else{
									f1mtr = 1;
									f1mun = 2;
								}
							}
						}
					}
					else{ //kid is 1/1
						if(!f1pa1 && f1pa2){
							f1ptr = 2;
							f1pun = 1;
						}
						else if(f1pa1 && f1pa2){
							f1ptr = 1;
							f1pun = 2;
						}
						if(!f1ma1 && f1ma2){
							f1mtr = 2;
							f1mun = 1;
						}
						else if(f1ma1 && f1ma2){
							f1mtr = 1;
							f1mun = 2;
						}
					}

					//family 2
					bool f2pa1 = dad2->getAone(mloc);
					bool f2pa2 = dad2->getAtwo(mloc);
					bool f2pa3 = dad2->getAmissing(mloc);
					bool f2ma1 = mom2->getAone(mloc);
					bool f2ma2 = mom2->getAtwo(mloc);
					bool f2ma3 = mom2->getAmissing(mloc);
					bool f2kid1 = child2->getAone(mloc);
					bool f2kid2 = child2->getAtwo(mloc);
					bool f2kid3 = child2->getAmissing(mloc);
					int f2ptr = 0;
					int f2pun = 0;
					int f2mtr = 0;
					int f2mun = 0;

					bool f2_allhet = false;

					//zero geno for fam2?
					if((f2pa1 && f2pa2 && f2pa3) || (f2ma1 && f2ma2 && f2ma3) || (f2kid1 && f2kid2 && f2kid3)){
						counts[0] = -1;
						counts[1] = -1;
						return counts;
					}
					if(!f2kid1 && !f2kid2){
						if(!f2pa1){
							f2ptr = 1;
							f2pun = 2;
						}
						if(!f2ma1){
							f2mtr = 1;
							f2mun = 2;
						}
					}
					else if(!f2kid1 && f2kid2){
						if(f2pa1 != f2pa2 && f2ma1 != f2ma2){
							f2_allhet = true;
						}
						else{
							if(f2pa1 != f2pa2){
								if(!f1ma1){
									f2ptr = 2;
									f2pun = 1;
								}
								else{
									f2ptr = 1;
									f2pun = 2;
								}
							}
							else{
								if(!f2pa1){
									f2mtr = 2;
									f2mun = 1;
								}
								else{
									f2mtr = 1;
									f2mun = 2;
								}
							}
						}
					}
					else{ //kid is 1/1
						if(!f2pa1 && f2pa2){
							f2ptr = 2;
							f2pun = 1;
						}
						else if(f2pa1 && f2pa2){
							f2ptr = 1;
							f2pun = 2;
						}
						if(!f2ma1 && f2ma2){
							f2mtr = 2;
							f2mun = 1;
						}
						else if(f2ma1 && f2ma2){
							f2mtr = 1;
							f2mun = 2;
						}
					}


					//compare pedigrees transmissions
					if(f1_allhet && f2_allhet){
						counts[0] += 0.25;
						counts[1] += 0.25;
					}
					else if(f1_allhet){
						if(f2ptr != f2mtr){
							counts[0] += 0.5;
							counts[1] += 0.5;
						}
						else{
							counts[0] += 0.25;
							counts[1] += 0.25;
						}
					}
					else if(f2_allhet){
						if(f1ptr != f2mtr){
							counts[0] += 0.5;
							counts[1] += 0.5;
						}
						else{
							counts[0] += 0.25;
							counts[1] += 0.25;
						}
					}
					else{
						if(f1ptr == f2ptr){
							counts[0]++;
						}
						if(f1mtr == f2mtr){
							counts[1]++;
						}
					}

				}
			}
		}
	}
	return counts;
}


//calculates average IBS over a pair of samples
double IBS::calcPairAverage(int s1, int s2){
	comparisons = 0;
	Sample* samp1 = data_set->get_sample(s1);
	Sample* samp2 = data_set->get_sample(s2);

	int sum = 0;
	int num_loci = 0;
	//pair of inds, snp, and 0, 1 or 2
	int msize = data_set->num_loci();
	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	msize = good_markers.size();

	for(int m = 0; m < msize; m++){
		Marker* mark = good_markers[m];
		if(mark->isEnabled()){
			int mloc = mark->getLoc();
			if((samp1->getAone(mloc) && samp1->getAtwo(mloc) && samp1->getAmissing(mloc))
					|| (samp2->getAone(mloc) && samp2->getAtwo(mloc) && samp2->getAmissing(mloc))){
				continue;
			}
			num_loci++;
			if(samp1->getAone(mloc) == samp2->getAone(mloc)){
				sum++;
			}
			if(samp1->getAtwo(mloc) == samp2->getAtwo(mloc)){
				sum++;
			}

		}
	}
	double avg = 0.0f;
	if(num_loci > 0){
		avg = (double) sum / (double) num_loci;
	}
	comparisons = num_loci;
	return avg;
}

//DEPRECATED
void IBS::calcOne(int m){
}

//NOT USED
void IBS::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

}

}
