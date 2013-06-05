/**********************************************************************************
*			Allele Frequency Module
*
* Written by: Justin Giles
*	          Vanderbilt University
*	          Center for Human Genetics Research
*
* Iterates over all genotypes and generates a Major/Minor allele count including
* frequencies as well as genotype frequencies.
*
*
*File: AlleleFrequency.cc
**********************************************************************************/

#include <unistd.h>
#include <sstream>
#include <stdio.h>
#include <iostream>
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
#include <map>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "AlleleFrequency.h"
#include "General.h"
#include "Helpers.h"
namespace Methods{
string AlleleFrequency::stepname = "allele-freq";


/*DEPRECATED
 *
 *Function: FilterSummary
 *Description:
 *Outputs total markers remaining after filtering
 *
 */
void AlleleFrequency::FilterSummary(){

	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;
}

/*DEPRECATED
 *
 *Function: PrintSummary
 *Description:
 *Resets marker flags in preparation for next step
 */
void AlleleFrequency::PrintSummary(){
	int msize = markers->size();

	for(int m = 0; m < msize; m++){
		(*markers).at(m)->setFlag(false);
	}
	return;
}

/*DEPRECATED
 *
 *Function: filter
 *Description:
 *Not used.
 */
void AlleleFrequency::filter(){
	return;
}


/*
 *Function: filterOne
 *Description:
 *Filters markers based on overall minor allele frequency
 */
void AlleleFrequency::filterOne(Marker* mark){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh() || options.doRmMono() || options.doRmHetOnly()){
		int notfoundcount = 0;
		double majfreq = 0;
		double minfreq = 0;

		if(!options.doFilterOverall() && !options.doFilterFile()){
			//overall is stored in a1/a2_countP
			if(!mark->isMicroSat()){
				//if frequencies are flipped?
				if(a1_countP > a2_countP){
					majfreq = ((double)a1_countP/(double)(a1_countP + a2_countP));
					minfreq = 1.0f - majfreq;
				}
				else{
					majfreq = ((double)a2_countP/(double)(a1_countP + a2_countP));
					minfreq = 1.0f - majfreq;
				}
				bool inc = false;

				//remove monozygous snps
				if((Helpers::dEquals(minfreq, 0) || Helpers::dEquals(majfreq, 1)) && options.doRmMono()){
					mark->setEnabled(false);
					orig_num_markers++;
					inc = true;
				}

				//remove heterozygous only snps
				if(a1_homo_countP == 0 && a2_homo_countP == 0 && options.doRmHetOnly()){
					mark->setEnabled(false);
					if(inc == false){
						orig_num_markers++;
						inc = true;
					}
				}

				if(Helpers::dLess(minfreq, options.getThreshMarkersLow()) && options.doThreshMarkersLow()){
					mark->setEnabled(false);
					if(inc == false){
						orig_num_markers++;
						inc = true;
					}
				}

				if(Helpers::dGreater(minfreq, options.getThreshMarkersHigh()) && options.doThreshMarkersHigh()){
					mark->setEnabled(false);
					if(inc == false){
						orig_num_markers++;
						inc = true;
					}
				}

			}
		}
		else if(options.doFilterOverall() && !options.doFilterFile()){
			//if filter overall enabled, overall stored in a1/a2_count
			if(!mark->isMicroSat()){
				//if frequencies are flipped?
				if(a1_count > a2_count){
					majfreq = ((double)a1_count/(double)(a1_count + a2_count));
					minfreq = 1.0f - majfreq;
				}
				else{
					majfreq = ((double)a2_count/(double)(a1_count + a2_count));
					minfreq = 1.0f - majfreq;
				}
				bool inc = false;

				//remove monozygous snps
				if((Helpers::dEquals(minfreq, 0) || Helpers::dEquals(majfreq, 1)) && options.doRmMono()){
					mark->setEnabled(false);
					orig_num_markers++;
					inc = true;
				}
				//remove het only snps
				if(a1_homo_count == 0 && a2_homo_count == 0 && options.doRmHetOnly()){
					mark->setEnabled(false);
					if(inc == false){
						orig_num_markers++;
						inc = true;
					}
				}

				if(Helpers::dLess(minfreq, options.getThreshMarkersLow()) && options.doThreshMarkersLow()){
					mark->setEnabled(false);
					if(inc == false){
						orig_num_markers++;
						inc = true;
					}
				}

				if(Helpers::dGreater(minfreq, options.getThreshMarkersHigh()) && options.doThreshMarkersHigh()){
					mark->setEnabled(false);
					if(inc == false){
						orig_num_markers++;
						inc = true;
					}
				}
			}
		}
		//if we are filtering by filter file
		else if(options.doFilterFile() && !options.doFilterOverall()){
			if(mark->hasMAF()){
				double majfreq = 0.0;
				double minfreq = 0.0;
				//still look at calculated values in order to weed out monozygous and het only
				if(a1_countP > a2_countP){
					majfreq = ((double)a1_countP/(double)(a1_countP + a2_countP));
					minfreq = 1.0f - majfreq;
				}
				else{
					majfreq = ((double)a2_countP/(double)(a1_countP + a2_countP));
					minfreq = 1.0f - majfreq;
				}
				bool inc = false;

				//remove mono only
				if((Helpers::dEquals(minfreq, 0) || Helpers::dEquals(majfreq, 1)) && options.doRmMono()){
					mark->setEnabled(false);
					orig_num_markers++;
					inc = true;
				}
				//remove het only
				if(a1_homo_countP == 0 && a2_homo_countP == 0 && options.doRmHetOnly()){
					mark->setEnabled(false);
					if(inc == false){
						orig_num_markers++;
						inc = true;
					}
				}

				//filter by filter file data
				if(Helpers::fLess(mark->getMAF(), options.getThreshMarkersLow()) && options.doThreshMarkersLow()){
					mark->setEnabled(false);
					if(inc == false){
						orig_num_markers++;
						inc = true;
					}
				}
				if(Helpers::fGreater(mark->getMAF(), options.getThreshMarkersHigh()) && options.doThreshMarkersHigh()){
					mark->setEnabled(false);
					if(inc == false){
						orig_num_markers++;
						inc = true;
					}
				}
			}
			else{
				notfoundcount++;
			}
		}
		if(options.doFilterFile()){
			if(notfoundcount > 0){
				opts::printLog("Markers skipped due to missing predefined MAF: " + getString<int>(notfoundcount) + "\n");
			}
		}
	}
}

/*
 *Function: flagSamples
 *Description:
 *Finds the samples that will be counted in the overall frequency calculation
 *
 */
void AlleleFrequency::flagSamples(){
	useoverall = false;
	if(options.doAll() || options.doFilterOverall()){
		useoverall = true;
	}

	if(sample_flags.size() == 0){
		sample_flags.resize(samples->size(), false);
	}
	int random;
	int ssize = samples->size();
	int fsize = families->size();
	for(int s = 0; s < ssize; s++){
		Sample* samp = (*samples).at(s);
		//set enabled initially

		//set everyone to usuable
		sample_flags.at(s) = true;

		if(options.doAll()){
			continue;
		}
		else if(options.doAllChildren()){
			//this is a founder, disable them
			if(samp->getDad() == NULL && samp->getMom() == NULL){
				sample_flags.at(s) = false;
			}
		}
		else if(options.doFoundersOnly()){
			//this is a child, disable it
			if(!samp->isFounder()){
				sample_flags.at(s) = false;
			}
		}
		else if(options.doUnaffSpousesOnly()){
			Sample* dad = samp->getDad();
			Sample* mom = samp->getMom();
			if(dad && dad->isEnabled() && dad->getDad() == NULL && dad->getMom() == NULL){
				if(dad->getPheno() != 1){
					sample_flags.at(dad->getLoc()) = false;
				}
			}
			if(mom && mom->isEnabled() && mom->getDad() == NULL && mom->getMom() == NULL){
				if(mom->getPheno() != 1){
					sample_flags.at(mom->getLoc()) = false;
				}
			}
			sample_flags.at(s) = false;
		}
		else if(options.doUnknownSpouses()){
			Sample* dad = samp->getDad();
			Sample* mom = samp->getMom();
			if(dad && dad->isEnabled() && dad->getDad() == NULL && dad->getMom() == NULL){
				if(dad->getPheno() == 0 || dad->getPheno() == 1){
					sample_flags.at(dad->getLoc()) = true;
				}
			}
			if(mom && mom->isEnabled() && mom->getDad() == NULL && mom->getMom() == NULL){
				if(mom->getPheno() == 0 || mom->getPheno() == 1){
					sample_flags.at(mom->getLoc()) = true;
				}
			}
			sample_flags.at(s) = false;
		}
	}

	//foreach family find a random child
	if(options.doRandomChild() && !options.doAll() && !options.doAllChildren() && !options.doFoundersOnly() && !options.doUnknownSpouses() && !options.doUnaffSpousesOnly()){
		for(int f = 0; f < fsize; f++){
			Family* fam = (*families).at(f);
			if(fam->isEnabled()){
				vector<Sample*>* nonfound = fam->getNonFounders();
				vector<Sample*> nf_enabled;
				for(int i = 0; i < (int)nonfound->size(); i++){
					Sample* samp = (*nonfound).at(i);
					if(samp->isEnabled()){
						nf_enabled.push_back(samp);
					}
				}
				if(nf_enabled.size() == 0){
					continue;
				}
				//find random child
				random = int(nf_enabled.size() * rand() / (RAND_MAX + 1.0));
				Sample* samp = nf_enabled.at(random);
				sample_flags[samp->getLoc()] = false;
				vector<Sample*>* samples = fam->getSamples();
				for(int i = 0; i < (int)samples->size(); i++){
					Sample* temp = (*samples).at(i);
					if(sample_flags.at(temp->getLoc())){
						sample_flags.at(temp->getLoc()) = false;
					}
					else{
						sample_flags.at(temp->getLoc()) = true;
					}
				}
			}
		}
	}
}

/*
 *Function: initializeCounts
 *Description:
 *Sets counts to the specified value (usually 0)
 *
 */
void AlleleFrequency::initializeCounts(int v){
	 a1_count = a2_count = a1_homo_count = a2_homo_count = a12_count = v;
	 a1_countM = a2_countM = a1_homo_countM = a2_homo_countM = a12_countM = v;
	 a1_countF = a2_countF = a1_homo_countF = a2_homo_countF = a12_countF = v;
	 a1_countP = a2_countP = a1_homo_countP = a2_homo_countP = a12_countP = v;
	 a1_countPM = a2_countPM = a1_homo_countPM = a2_homo_countPM = a12_countPM = v;
	 a1_countPF = a2_countPF = a1_homo_countPF = a2_homo_countPF = a12_countPF = v;
	 a1_countC = a2_countC = a1_homo_countC = a2_homo_countC = a12_countC = v;
	 a1_countCM = a2_countCM = a1_homo_countCM = a2_homo_countCM = a12_countCM = v;
	 a1_countCF = a2_countCF = a1_homo_countCF = a2_homo_countCF = a12_countCF = v;
	 a1_countCa = a2_countCa = a1_homo_countCa = a2_homo_countCa = a12_countCa = v;
	 a1_countCaM = a2_countCaM = a1_homo_countCaM = a2_homo_countCaM = a12_countCaM = v;
	 a1_countCaF = a2_countCaF = a1_homo_countCaF = a2_homo_countCaF = a12_countCaF = v;
	 a1_countCon = a2_countCon = a1_homo_countCon = a2_homo_countCon = a12_countCon = v;
	 a1_countConM = a2_countConM = a1_homo_countConM = a2_homo_countConM = a12_countConM = v;
	 a1_countConF = a2_countConF = a1_homo_countConF = a2_homo_countConF = a12_countConF = v;

	 ga1_count.clear();
	 ga2_count.clear();
	 ga1_homo_count.clear();
	 ga2_homo_count.clear();
	 ga12_count.clear();

	 gm_allele_counts_o.clear();
	 gm_geno_counts_o.clear();

	 m_allele_counts_o.clear();
	 m_allele_counts_om.clear();
	 m_allele_counts_of.clear();
	 m_geno_counts_o.clear();
	 m_geno_counts_om.clear();
	 m_geno_counts_of.clear();
	 m_allele_counts_p.clear();
	 m_allele_counts_pm.clear();
	 m_allele_counts_pf.clear();
	 m_geno_counts_p.clear();
	 m_geno_counts_pm.clear();
	 m_geno_counts_pf.clear();
	 m_allele_counts_c.clear();
	 m_allele_counts_cm.clear();
	 m_allele_counts_cf.clear();
	 m_geno_counts_c.clear();
	 m_geno_counts_cm.clear();
	 m_geno_counts_cf.clear();
	 m_allele_counts_ca.clear();
	 m_allele_counts_cam.clear();
	 m_allele_counts_caf.clear();
	 m_geno_counts_ca.clear();
	 m_geno_counts_cam.clear();
	 m_geno_counts_caf.clear();
	 m_allele_counts_con.clear();
	 m_allele_counts_conm.clear();
	 m_allele_counts_conf.clear();
	 m_geno_counts_con.clear();
	 m_geno_counts_conm.clear();
	 m_geno_counts_conf.clear();

}

/*
 *Function: findRandomSample
 *Description:
 *Returns a random sample, null if no samples are available
 *Used when no founders have genotypes typically.
 */
Sample* AlleleFrequency::findRandomSample(Sample* tsamp, vector<int>& rsamps, Marker* mark){
	Family* fam = tsamp->getFamily();
	vector<Sample*>* samps = fam->getSamples();
	int random;
	random = int(rsamps.size() * rand() / (RAND_MAX + 1.0));
	Sample* samp = (*samps)[rsamps.at(random)];
	vector<int>::iterator found = find(rsamps.begin(), rsamps.end(), rsamps.at(random));
	rsamps.erase(found);
	if(sample_flags[samp->getLoc()]){
		return NULL;
	}
	if(!mark->isMicroSat()){
		if(samp->getAone(mark->getLoc()) && samp->getAtwo(mark->getLoc()) && samp->getAmissing(mark->getLoc())){
			return NULL;
		}
	}
	else{
		if(samp->getAbone(mark->getLoc()) == -1 && samp->getAbtwo(mark->getLoc()) == -1){
			return NULL;
		}
	}
	return samp;

}

/*
 *Function: calcOneGroups
 *Description:
 *Performs calculations for a single marker for sample groups.
 */
void AlleleFrequency::calcOneGroups(Marker* mark){
//	int ssize = samples->size();
	map<string, vector<Sample*> > groups = options.getGroups();
	map<string, vector<Sample*> >::iterator giter;
	for(giter = groups.begin(); giter != groups.end(); giter++){
		string mygroup = giter->first;
		vector<Sample*> fsamps = giter->second;
		int ssize = fsamps.size();
		for(int i = 0; i < ssize; i++){
			if(fsamps.at(i)->isEnabled()){
				int loc = mark->getLoc();
				bool isX = false;
				//if X chrom, then count males once
				if(mark->getChrom() == opts::_CHRX_){
					isX = true;
				}
				if(!mark->isMicroSat()){
					if(isX && fsamps.at(i)->getSex()){
						if(!fsamps.at(i)->getAone(loc) && !fsamps.at(i)->getAtwo(loc)){
						}
						else if(fsamps.at(i)->getAone(loc) && fsamps.at(i)->getAtwo(loc) && !fsamps.at(i)->getAmissing(loc)){
						}
					}
					else{
						if(!fsamps.at(i)->getAone(loc)){
							ga1_count[mygroup]++;
						}
						else if(fsamps.at(i)->getAone(loc) && fsamps.at(i)->getAtwo(loc) && !fsamps.at(i)->getAmissing(loc)){
							ga2_count[mygroup]++;
						}

						if(!fsamps.at(i)->getAtwo(loc) && !fsamps.at(i)->getAone(loc)){
							ga1_count[mygroup]++;
						}
						else if(fsamps.at(i)->getAtwo(loc) && !fsamps.at(i)->getAmissing(loc)){
							ga2_count[mygroup]++;
						}
					}

					if(isX && fsamps.at(i)->getSex()){
						continue;
					}

					if(!fsamps.at(i)->getAone(loc) && !fsamps.at(i)->getAtwo(loc)){
						ga1_homo_count[mygroup]++;
					}
					if(!fsamps.at(i)->getAone(loc) && fsamps.at(i)->getAtwo(loc)){
						ga12_count[mygroup]++;
					}
					if(fsamps.at(i)->getAone(loc) && fsamps.at(i)->getAtwo(loc) && !fsamps.at(i)->getAmissing(loc)){
						ga2_homo_count[mygroup]++;
					}
				}//end !microsat
				else
				{//is microsat
					if(gm_allele_counts_o[mygroup].size() == 0){
						int numalleles = mark->getNumAlleles();
						gm_allele_counts_o[mygroup].resize(numalleles, 0);
					}
					if(fsamps.at(i)->getAbone(loc) != -1){
						int a1 = fsamps.at(i)->getAbone(loc);
						int a2 = fsamps.at(i)->getAbtwo(loc);
						if(isX && fsamps.at(i)->getSex()){
							if(a1 == a2){
							}
							else{
							}
						}
						else{
							gm_allele_counts_o[mygroup][a1]++;
							gm_allele_counts_o[mygroup][a2]++;
						}

						if(isX && fsamps.at(i)->getSex()){
							continue;
						}
						//genotype counts
						gm_geno_counts_o[mygroup][getString<int>(a1) + "_" + getString<int>(a2)]++;
					}
				}
			}
		}//end sample iter
	}//end family iter
}

/*
 *Function: calcOne
 *Description:
 *Performs calculations for a single marker.
 */
void AlleleFrequency::calcOne(Marker* mark){
	initializeCounts(0);
	int fsize = families->size();
	for(int f = 0; f < fsize; f++){
		Family* fam = (*families).at(f);
		if(!fam->isEnabled()){
			continue;
		}
		vector<Sample*>* fsamps = fam->getSamples();
		int ssize = fsamps->size();
		bool famdone = false;

		if(options.doFoundersOnly()){
			vector<Sample*>* founders = fam->getFounders();
			for(int i = 0; i < (int)founders->size(); i++){
				int loc = mark->getLoc();
				if(!mark->isMicroSat()){
					if((*founders).at(i)->getAone(loc) && (*founders).at(i)->getAtwo(loc) && (*founders).at(i)->getAmissing(loc)){
						famdone = false;
					}
					else{
						famdone = true;
					}
				}
				else{
					if((*founders).at(i)->getAbone(loc) == -1 && (*founders).at(i)->getAbtwo(loc) == -1){
						famdone = false;
					}
					else{
						famdone = true;
					}
				}
				if(famdone){
					break;
				}
			}
		}

		//iterate over samples
		for(int i = 0; i < ssize; i++){
			if((*fsamps).at(i)->isEnabled()){
				int loc = mark->getLoc();
				bool isX = false;
				//if X chrom, then count males once
				if(mark->getChrom() == opts::_CHRX_){
					isX = true;
				}
				if(!mark->isMicroSat()){
					if((*fsamps).at(i)->getAone(loc) && (*fsamps).at(i)->getAtwo(loc) && (*fsamps).at(i)->getAmissing(loc)){
						if((*fsamps).at(i)->isFounder() && options.doFoundersOnly() && !famdone){
							Sample* randsamp = NULL;
							int indcount = 0;
							int famsamps = (*fsamps).at(i)->getFamily()->getSamples()->size();
							vector<int> rsamps;
							for(int r = 0; r < famsamps; r++){
								rsamps.push_back(r);
							}
							while(randsamp == NULL && rsamps.size() > 0){
								randsamp = findRandomSample((*fsamps).at(i), rsamps, mark);
								indcount++;
							}
							if(randsamp == NULL){
								continue;
							}
							else{
								famdone = true;
								if(isX && randsamp->getSex()){
									if(!randsamp->getAone(loc) && !randsamp->getAtwo(loc)){
									}
									else if(randsamp->getAone(loc) && randsamp->getAtwo(loc) && !randsamp->getAmissing(loc)){
									}
								}
								else{
									if(!randsamp->getAone(loc) && !randsamp->getAtwo(loc)){
										a1_count++;
										a1_count++;
										a1_homo_count++;
										if(randsamp->getSex()){
											a1_countM++;
											a1_countM++;
											a1_homo_countM++;
										}
										else{
											a1_countF++;
											a1_countF++;
											a1_homo_countF++;
										}
									}
									else if(!randsamp->getAone(loc) && randsamp->getAtwo(loc)){
										a1_count++;
										a2_count++;
										a12_count++;
										if(randsamp->getSex()){
											a1_countM++;
											a2_countM++;
											a12_countM++;
										}
										else{
											a1_countF++;
											a2_countF++;
											a12_countF++;
										}
									}
									else if(randsamp->getAone(loc) && randsamp->getAtwo(loc) && !randsamp->getAmissing(loc)){
										a2_count++;
										a2_count++;
										a2_homo_count++;
										if(randsamp->getSex()){
											a2_countM++;
											a2_countM++;
											a2_homo_countM++;
										}
										else{
											a2_countF++;
											a2_countF++;
											a2_homo_countF++;
										}
									}
								}
							}
						}
							continue;
					}

					//look at a family only once
					if(sample_flags[(*fsamps).at(i)->getLoc()] &&  (*fsamps).at(i)->isFounder() && options.doFoundersOnly()){
						famdone = true;
					}

					if(isX && (*fsamps).at(i)->getSex()){
						if(!(*fsamps).at(i)->getAone(loc) && !(*fsamps).at(i)->getAtwo(loc)){
							if(sample_flags[(*fsamps).at(i)->getLoc()]){
							}
							if((*fsamps).at(i)->getPheno() == 2){
							}
							else if((*fsamps).at(i)->getPheno() == 1){
							}
							if((*fsamps).at(i)->isFounder()){
							}
							else{
							}
						}
						else if((*fsamps).at(i)->getAone(loc) && (*fsamps).at(i)->getAtwo(loc) && (*fsamps).at(i)->getAmissing(loc)){
							if(sample_flags[(*fsamps).at(i)->getLoc()]){
							}
							if((*fsamps).at(i)->getPheno() == 2){
							}
							else if((*fsamps).at(i)->getPheno() == 1){
							}
							if((*fsamps).at(i)->isFounder()){
							}
							else{
							}
						}
					}
					else{
						if(!(*fsamps).at(i)->getAone(loc)){
							if(sample_flags[(*fsamps).at(i)->getLoc()]){
								a1_count++;
								if((*fsamps).at(i)->getSex()){
									a1_countM++;
								}
								else{
									a1_countF++;
								}
							}
							if((*fsamps).at(i)->getPheno() == 2){
								a1_countCa++;
								if((*fsamps).at(i)->getSex()){
									a1_countCaM++;
								}
								else{
									a1_countCaF++;
								}
							}
							else if((*fsamps).at(i)->getPheno() == 1){
								a1_countCon++;
								if((*fsamps).at(i)->getSex()){
									a1_countConM++;
								}
								else{
									a1_countConF++;
								}
							}
							if((*fsamps).at(i)->isFounder()){
								a1_countP++;
								if((*fsamps).at(i)->getSex()){
									a1_countPM++;
								}
								else{
									a1_countPF++;
								}
							}
							else if((*fsamps).at(i)->getSex()){
								a1_countC++;
								a1_countCM++;
							}
							else{
								a1_countC++;
								a1_countCF++;
							}
						}
						else if((*fsamps).at(i)->getAone(loc) && (*fsamps).at(i)->getAtwo(loc) && !(*fsamps).at(i)->getAmissing(loc)){
							if(sample_flags[(*fsamps).at(i)->getLoc()]){
								a2_count++;
								if((*fsamps).at(i)->getSex()){
									a2_countM++;
								}
								else{
									a2_countF++;
								}
							}
							if((*fsamps).at(i)->getPheno() == 2){
								a2_countCa++;
								if((*fsamps).at(i)->getSex()){
									a2_countCaM++;
								}
								else{
									a2_countCaF++;
								}
							}
							else if((*fsamps).at(i)->getPheno() == 1){
								a2_countCon++;
								if((*fsamps).at(i)->getSex()){
									a2_countConM++;
								}
								else{
									a2_countConF++;
								}
							}

							if((*fsamps).at(i)->isFounder()){
								a2_countP++;
								if((*fsamps).at(i)->getSex()){
									a2_countPM++;
								}
								else{
									a2_countPF++;
								}
							}
							else if((*fsamps).at(i)->getSex()){
								a2_countC++;
								a2_countCM++;
							}
							else{
								a2_countC++;
								a2_countCF++;
							}
						}

						if(!(*fsamps).at(i)->getAtwo(loc) && !(*fsamps).at(i)->getAone(loc)){
							if(sample_flags[(*fsamps).at(i)->getLoc()]){
								a1_count++;
								if((*fsamps).at(i)->getSex()){
									a1_countM++;
								}
								else{
									a1_countF++;
								}
							}
							if((*fsamps).at(i)->getPheno() == 2){
								a1_countCa++;
								if((*fsamps).at(i)->getSex()){
									a1_countCaM++;
								}
								else{
									a1_countCaF++;
								}
							}
							else if((*fsamps).at(i)->getPheno() == 1){
								a1_countCon++;
								if((*fsamps).at(i)->getSex()){
									a1_countConM++;
								}
								else{
									a1_countConF++;
								}
							}
							if((*fsamps).at(i)->isFounder()){
								a1_countP++;
								if((*fsamps).at(i)->getSex()){
									a1_countPM++;
								}
								else{
									a1_countPF++;
								}
							}
							else if((*fsamps).at(i)->getSex()){
								a1_countC++;
								a1_countCM++;
							}
							else{
								a1_countC++;
								a1_countCF++;
							}
						}
						else if((*fsamps).at(i)->getAtwo(loc)){
							if(sample_flags[(*fsamps).at(i)->getLoc()]){
								a2_count++;
								if((*fsamps).at(i)->getSex()){
									a2_countM++;
								}
								else{
									a2_countF++;
								}
							}
							if((*fsamps).at(i)->getPheno() == 2){
								a2_countCa++;
								if((*fsamps).at(i)->getSex()){
									a2_countCaM++;
								}
								else{
									a2_countCaF++;
								}
							}
							else if((*fsamps).at(i)->getPheno() == 1){
								a2_countCon++;
								if((*fsamps).at(i)->getSex()){
									a2_countConM++;
								}
								else{
									a2_countConF++;
								}
							}
							if((*fsamps).at(i)->isFounder()){
								a2_countP++;
								if((*fsamps).at(i)->getSex()){
									a2_countPM++;
								}
								else{
									a2_countPF++;
								}
							}
							else if((*fsamps).at(i)->getSex()){
								a2_countC++;
								a2_countCM++;
							}
							else{
								a2_countC++;
								a2_countCF++;
							}
						}
					}

					if(isX && (*fsamps).at(i)->getSex()){
						continue;
					}
					if(!(*fsamps).at(i)->getAone(loc) && !(*fsamps).at(i)->getAtwo(loc)){
						if(sample_flags[(*fsamps).at(i)->getLoc()]){
							a1_homo_count++;
							if((*fsamps).at(i)->getSex()){
								a1_homo_countM++;
							}
							else{
								a1_homo_countF++;
							}
						}
						if((*fsamps).at(i)->getPheno() == 2){
							a1_homo_countCa++;
							if((*fsamps).at(i)->getSex()){
								a1_homo_countCaM++;
							}
							else{
								a1_homo_countCaF++;
							}
						}
						else if((*fsamps).at(i)->getPheno() == 1){
							a1_homo_countCon++;
							if((*fsamps).at(i)->getSex()){
								a1_homo_countConM++;
							}
							else{
								a1_homo_countConF++;
							}
						}
						if((*fsamps).at(i)->isFounder()){
							a1_homo_countP++;
							if((*fsamps).at(i)->getSex()){
								a1_homo_countPM++;
							}
							else{
								a1_homo_countPF++;
							}
						}
						else if((*fsamps).at(i)->getSex()){
							a1_homo_countC++;
							a1_homo_countCM++;
						}
						else{
							a1_homo_countC++;
							a1_homo_countCF++;
						}
					}
					if(!(*fsamps).at(i)->getAone(loc) && (*fsamps).at(i)->getAtwo(loc)){
						if(sample_flags[(*fsamps).at(i)->getLoc()]){
							a12_count++;
							if((*fsamps).at(i)->getSex()){
								a12_countM++;
							}
							else{
								a12_countF++;
							}
						}
						if((*fsamps).at(i)->getPheno() == 2){
							a12_countCa++;
							if((*fsamps).at(i)->getSex()){
								a12_countCaM++;
							}
							else{
								a12_countCaF++;
							}
						}
						else if((*fsamps).at(i)->getPheno() == 1){
							a12_countCon++;
							if((*fsamps).at(i)->getSex()){
								a12_countConM++;
							}
							else{
								a12_countConF++;
							}
						}
						if((*fsamps).at(i)->isFounder()){
							a12_countP++;
							if((*fsamps).at(i)->getSex()){
								a12_countPM++;
							}
							else{
								a12_countPF++;
							}
						}
						else if((*fsamps).at(i)->getSex()){
							a12_countC++;
							a12_countCM++;
						}
						else{
							a12_countC++;
							a12_countCF++;
						}
					}
					if((*fsamps).at(i)->getAone(loc) && (*fsamps).at(i)->getAtwo(loc) && !(*fsamps).at(i)->getAmissing(loc)){
						if(sample_flags[(*fsamps).at(i)->getLoc()]){
							a2_homo_count++;
							if((*fsamps).at(i)->getSex()){
								a2_homo_countM++;
							}
							else{
								a2_homo_countF++;
							}
						}
						if((*fsamps).at(i)->getPheno() == 2){
							a2_homo_countCa++;
							if((*fsamps).at(i)->getSex()){
								a2_homo_countCaM++;
							}
							else{
								a2_homo_countCaF++;
							}
						}
						else if((*fsamps).at(i)->getPheno() == 1){
							a2_homo_countCon++;
							if((*fsamps).at(i)->getSex()){
								a2_homo_countConM++;
							}
							else{
								a2_homo_countConF++;
							}
						}
						if((*fsamps).at(i)->isFounder()){
							a2_homo_countP++;
							if((*fsamps).at(i)->getSex()){
								a2_homo_countPM++;
							}
							else{
								a2_homo_countPF++;
							}
						}
						else if((*fsamps).at(i)->getSex()){
							a2_homo_countC++;
							a2_homo_countCM++;
						}
						else{
							a2_homo_countC++;
							a2_homo_countCF++;
						}
					}
				}//end !microsat
				else{//is microsat
					if(m_allele_counts_o.size() == 0){
						int numalleles = mark->getNumAlleles();
						m_allele_counts_o.resize(numalleles, 0);
						m_allele_counts_om.resize(numalleles, 0);
						m_allele_counts_of.resize(numalleles, 0);
						m_allele_counts_p.resize(numalleles, 0);
						m_allele_counts_pm.resize(numalleles, 0);
						m_allele_counts_pf.resize(numalleles, 0);
						m_allele_counts_c.resize(numalleles, 0);
						m_allele_counts_cm.resize(numalleles, 0);
						m_allele_counts_cf.resize(numalleles, 0);
						m_allele_counts_ca.resize(numalleles, 0);
						m_allele_counts_cam.resize(numalleles, 0);
						m_allele_counts_caf.resize(numalleles, 0);
						m_allele_counts_con.resize(numalleles, 0);
						m_allele_counts_conm.resize(numalleles, 0);
						m_allele_counts_conf.resize(numalleles, 0);
					}
					if((*fsamps).at(i)->getAbone(loc) == -1){
						if((*fsamps).at(i)->isFounder() && options.doFoundersOnly() && !famdone){
							Sample* randsamp = NULL;
							int indcount = 0;
							int famsamps = (*fsamps).at(i)->getFamily()->getSamples()->size();
							vector<int> rsamps;
							for(int r = 0; r < famsamps; r++){
								rsamps.push_back(r);
							}
							while(randsamp == NULL && rsamps.size() > 0){
								randsamp = findRandomSample((*fsamps).at(i), rsamps, mark);
								indcount++;
							}
							if(randsamp == NULL){
								continue;
							}
							else{
								famdone = true;
								int a1 = randsamp->getAbone(loc);
								int a2 = randsamp->getAbtwo(loc);
								if(isX && randsamp->getSex()){
									if(a1 == a2){
									}
									else{
									}
								}
								else{
									if(a1 == a2){
										m_allele_counts_o.at(a1)++;
										m_allele_counts_o.at(a1)++;
										m_geno_counts_o[getString<int>(a1) + "_" + getString<int>(a1)]++;
										if(randsamp->getSex()){
											m_allele_counts_om.at(a1)++;
											m_allele_counts_om.at(a1)++;
											m_geno_counts_om[getString<int>(a1) + "_" + getString<int>(a1)]++;
										}
										else{
											m_allele_counts_of.at(a1)++;
											m_allele_counts_of.at(a1)++;
											m_geno_counts_of[getString<int>(a1) + "_" + getString<int>(a1)]++;
										}
									}
									else{
										m_allele_counts_o.at(a1)++;
										m_allele_counts_o.at(a2)++;
										m_geno_counts_o[getString<int>(a1) + "_" + getString<int>(a2)]++;
										if(randsamp->getSex()){
											m_allele_counts_om.at(a1)++;
											m_allele_counts_om.at(a2)++;
											m_geno_counts_om[getString<int>(a1) + "_" + getString<int>(a2)]++;
										}
										else{
											m_allele_counts_of.at(a1)++;
											m_allele_counts_of.at(a2)++;
											m_geno_counts_of[getString<int>(a1) + "_" + getString<int>(a2)]++;
										}
									}
								}
							}
						}
						else{
							continue;
						}
					}
					if((*fsamps).at(i)->getAbone(loc) != -1){
						int a1 = (*fsamps).at(i)->getAbone(loc);
						int a2 = (*fsamps).at(i)->getAbtwo(loc);
						if(isX && (*fsamps).at(i)->getSex()){
							if(a1 == a2){
								if(sample_flags[(*fsamps).at(i)->getLoc()]){
								}
								if((*fsamps).at(i)->getPheno() == 2){
								}
								else if((*fsamps).at(i)->getPheno() == 1){
								}
								if((*fsamps).at(i)->isFounder()){
								}
								else{
								}
							}
							else{
								if(sample_flags[(*fsamps).at(i)->getLoc()]){
								}
								if((*fsamps).at(i)->getPheno() == 2){
								}
								else if((*fsamps).at(i)->getPheno() == 1){
								}
								if((*fsamps).at(i)->isFounder()){
								}
								else{
								}
							}
						}
						else{
							if(sample_flags[(*fsamps).at(i)->getLoc()]){
								m_allele_counts_o.at(a1)++;
								m_allele_counts_o.at(a2)++;
							}
							if((*fsamps).at(i)->getSex()){
								if(sample_flags[(*fsamps).at(i)->getLoc()]){
									m_allele_counts_om.at(a1)++;
									m_allele_counts_om.at(a2)++;
								}
								if((*fsamps).at(i)->getPheno() == 2){
									m_allele_counts_ca.at(a1)++;
									m_allele_counts_ca.at(a2)++;
									m_allele_counts_cam.at(a1)++;
									m_allele_counts_cam.at(a2)++;
								}
								else if((*fsamps).at(i)->getPheno() == 1){
									m_allele_counts_con.at(a1)++;
									m_allele_counts_con.at(a2)++;
									m_allele_counts_conm.at(a1)++;
									m_allele_counts_conm.at(a2)++;
								}
								if((*fsamps).at(i)->isFounder()){
									m_allele_counts_p.at(a1)++;
									m_allele_counts_p.at(a2)++;
									m_allele_counts_pm.at(a1)++;
									m_allele_counts_pm.at(a2)++;
								}
								else{
									m_allele_counts_c.at(a1)++;
									m_allele_counts_c.at(a2)++;
									m_allele_counts_cm.at(a1)++;
									m_allele_counts_cm.at(a2)++;
								}

							}
							else{
								if(sample_flags[(*fsamps).at(i)->getLoc()]){
									m_allele_counts_of.at(a1)++;
									m_allele_counts_of.at(a2)++;
								}
								if((*fsamps).at(i)->getPheno() == 2){
									m_allele_counts_ca.at(a1)++;
									m_allele_counts_ca.at(a2)++;
									m_allele_counts_caf.at(a1)++;
									m_allele_counts_caf.at(a2)++;
								}
								else if((*fsamps).at(i)->getPheno() == 1){
									m_allele_counts_con.at(a1)++;
									m_allele_counts_con.at(a2)++;
									m_allele_counts_conf.at(a1)++;
									m_allele_counts_conf.at(a2)++;
								}
								if((*fsamps).at(i)->isFounder()){
									m_allele_counts_p.at(a1)++;
									m_allele_counts_p.at(a2)++;
									m_allele_counts_pf.at(a1)++;
									m_allele_counts_pf.at(a2)++;
								}
								else{
									m_allele_counts_c.at(a1)++;
									m_allele_counts_c.at(a2)++;
									m_allele_counts_cf.at(a1)++;
									m_allele_counts_cf.at(a2)++;
								}

							}
						}
						if(isX && (*fsamps).at(i)->getSex()){
							continue;
						}
						//genotype counts
						if(sample_flags[(*fsamps).at(i)->getLoc()]){
							m_geno_counts_o[getString<int>(a1) + "_" + getString<int>(a2)]++;
						}
						if((*fsamps).at(i)->getSex()){
							if(sample_flags[(*fsamps).at(i)->getLoc()]){
								m_geno_counts_om[getString<int>(a1) + "_" + getString<int>(a2)]++;
							}
							if((*fsamps).at(i)->getPheno() == 2){
								m_geno_counts_ca[getString<int>(a1) + "_" + getString<int>(a2)]++;
								m_geno_counts_cam[getString<int>(a1) + "_" + getString<int>(a2)]++;
							}
							else if((*fsamps).at(i)->getPheno() == 1){
								m_geno_counts_con[getString<int>(a1) + "_" + getString<int>(a2)]++;
								m_geno_counts_conm[getString<int>(a1) + "_" + getString<int>(a2)]++;
							}
							if((*fsamps).at(i)->isFounder()){
								m_geno_counts_p[getString<int>(a1) + "_" + getString<int>(a2)]++;
								m_geno_counts_pm[getString<int>(a1) + "_" + getString<int>(a2)]++;
							}
							else{
								m_geno_counts_c[getString<int>(a1) + "_" + getString<int>(a2)]++;
								m_geno_counts_cm[getString<int>(a1) + "_" + getString<int>(a2)]++;
							}

						}
						else{
							if(sample_flags[(*fsamps).at(i)->getLoc()]){
								m_geno_counts_of[getString<int>(a1) + "_" + getString<int>(a2)]++;
							}
							if((*fsamps).at(i)->getPheno() == 2){
								m_geno_counts_ca[getString<int>(a1) + "_" + getString<int>(a2)]++;
								m_geno_counts_caf[getString<int>(a1) + "_" + getString<int>(a2)]++;
							}
							else if((*fsamps).at(i)->getPheno() == 1){
								m_geno_counts_con[getString<int>(a1) + "_" + getString<int>(a2)]++;
								m_geno_counts_conf[getString<int>(a1) + "_" + getString<int>(a2)]++;
							}
							if((*fsamps).at(i)->isFounder()){
								m_geno_counts_p[getString<int>(a1) + "_" + getString<int>(a2)]++;
								m_geno_counts_pf[getString<int>(a1) + "_" + getString<int>(a2)]++;
							}
							else{
								m_geno_counts_c[getString<int>(a1) + "_" + getString<int>(a2)]++;
								m_geno_counts_cf[getString<int>(a1) + "_" + getString<int>(a2)]++;
							}

						}

					}
				}
			}
		}//end sample iter
	}//end family iter
}


/*
 *Function: processtest
 *Description:
 *Main function to perform allele frequency test
 *
 */
void AlleleFrequency::processtest(){
	string afname = opts::_OUTPREFIX_ + "allele_freq" + options.getOut() + ".txt";
	if(!overwrite){
		afname += "." + getString<int>(order);
	}
	string gfname = opts::_OUTPREFIX_ + "allele_freq_genotype" + options.getOut() + ".txt";
	if(!overwrite){
		gfname += "." + getString<int>(order);
	}
	string gafname;
	string ggfname;
	if(options.doGroupFile()){
		gafname = opts::_OUTPREFIX_ + "allele_freq_group" + options.getOut() + ".txt";
		if(!overwrite){
			gafname += "." + getString<int>(order);
		}
		ggfname = opts::_OUTPREFIX_ + "allele_freq_genotype_group" + options.getOut() + ".txt";
		if(!overwrite){
			ggfname += "." + getString<int>(order);
		}
	}

    string pfname = opts::_OUTPREFIX_ + "allele_freq_parental" + options.getOut() + ".txt";
	if(!overwrite){
		pfname += "." + getString<int>(order);
	}
    string pgfname = opts::_OUTPREFIX_ + "allele_freq_parental_genotype" + options.getOut() + ".txt";
	if(!overwrite){
		pgfname += "." + getString<int>(order);
	}
    string gendfname = opts::_OUTPREFIX_ + "allele_freq_gender" + options.getOut() + ".txt";
	if(!overwrite){
		gendfname += "." + getString<int>(order);
	}
    string gendgfname = opts::_OUTPREFIX_ + "allele_freq_gender_genotype" + options.getOut() + ".txt";
	if(!overwrite){
		gendgfname += "." + getString<int>(order);
	}
    string ccfname = opts::_OUTPREFIX_ + "allele_freq_casecontrol" + options.getOut() + ".txt";
	if(!overwrite){
		ccfname += "." + getString<int>(order);
	}
    string ccgfname = opts::_OUTPREFIX_ + "allele_freq_casecontrol_genotype" + options.getOut() + ".txt";
	if(!overwrite){
		ccgfname += "." + getString<int>(order);
	}
    ofstream paren;
    ofstream pareng;
    ofstream gend;
    ofstream gendg;
    ofstream cc;
    ofstream ccg;

	ofstream myoutput (afname.c_str());
	ofstream mygeno (gfname.c_str());
	opts::addFile("Marker", stepname, afname);
	opts::addFile("Marker", stepname, gfname);
    if(!myoutput){
        opts::printLog("Error opening: " + afname + ".  Exiting!\n");
        throw MethodException("Error opening: " + afname + ".  Exiting!\n");
    }
    if(!mygeno){
        opts::printLog("Error opening: " + gfname + ".  Exiting!\n");
        throw MethodException("Error opening: " + gfname + ".  Exiting!\n");
    }

	ofstream gmyoutput; //group file output
	ofstream gmygeno; //group file output
	if(options.doGroupFile()){
		gmyoutput.open(gafname.c_str());
		gmygeno.open(ggfname.c_str());
		opts::addFile("Marker", stepname, gafname);
		opts::addFile("Marker", stepname, ggfname);
    	if(!gmyoutput){
	        opts::printLog("Error opening: " + gafname + ".  Exiting!\n");
	        throw MethodException("Error opening: " + gafname + ".  Exiting!\n");
    	}
    	if(!gmygeno){
        	opts::printLog("Error opening: " + ggfname + ".  Exiting!\n");
        	throw MethodException("Error opening: " + ggfname + ".  Exiting!\n");
    	}
		gmyoutput.precision(4);
		gmygeno.precision(4);
	}

	int msize = markers->size();

    int maxalleles = 0;
    for(int i = 0; i < (int)markers->size(); i++){
        if((*markers).at(i)->isEnabled()){
            if((*markers).at(i)->getNumAlleles() > maxalleles){
                maxalleles = (*markers).at(i)->getNumAlleles();
            }
        }
    }

	myoutput.precision(4);
	mygeno.precision(4);
	if((*markers).at(0)->getDetailHeaders().size() > 0){
    	myoutput << "Chrom\trsID\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders();
    	mygeno << "Chrom\trsID\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders() << "\tGenotype11\tGenotype12\tGenotype22\tOverall_Freq_Genotype11\tOverall_Freq_Genotype12\tOverall_Freq_Genotype22\tOverall_Count_Genotype11\tOverall_Count_Genotype12\tOverall_Count_Genotype22\tCase_Freq_Genotype11\tCase_Freq_Genotype12\tCase_Freq_Genotype22\tCase_Count_Genotype11\tCase_Count_Genotype12\tCase_Count_Genotype22\tControl_Freq_Genotype11\tControl_Freq_Genotype12\tControl_Freq_Genotype22\tControl_Count_Genotype11\tControl_Count_Genotype12\tControl_Count_Genotype22";
		if(options.doGroupFile()){
			gmyoutput << "Chrom\trsID\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders();
			gmygeno << "Chrom\trsID\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders() << "\tGenotype11\tGenotype12\tGenotype22";
		}
	}
	else{
    	myoutput << "Chrom\trsID\tProbeID\tbploc";
    	mygeno << "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22\tOverall_Freq_Genotype11\tOverall_Freq_Genotype12\tOverall_Freq_Genotype22\tOverall_Count_Genotype11\tOverall_Count_Genotype12\tOverall_Count_Genotype22\tCase_Freq_Genotype11\tCase_Freq_Genotype12\tCase_Freq_Genotype22\tCase_Count_Genotype11\tCase_Count_Genotype12\tCase_Count_Genotype22\tControl_Freq_Genotype11\tControl_Freq_Genotype12\tControl_Freq_Genotype22\tControl_Count_Genotype11\tControl_Count_Genotype12\tControl_Count_Genotype22";
		if(options.doGroupFile()){
			gmyoutput << "Chrom\trsID\tProbeID\tbploc";
			gmygeno << "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22";
		}
	}
	if(options.doGroupFile()){
		opts::addHeader(ggfname, "Genotype11");
		opts::addHeader(ggfname, "Genotype12");
		opts::addHeader(ggfname, "Genotype22");
	}

	opts::addHeader(gfname, "Genotype11");
	opts::addHeader(gfname, "Genotype12");
	opts::addHeader(gfname, "Genotype22");
	opts::addHeader(gfname, "Overall_Freq_Genotype11");
	opts::addHeader(gfname, "Overall_Freq_Genotype12");
	opts::addHeader(gfname, "Overall_Freq_Genotype22");
	opts::addHeader(gfname, "Overall_Count_Genotype11");
	opts::addHeader(gfname, "Overall_Count_Genotype12");
	opts::addHeader(gfname, "Overall_Count_Genotype22");
	opts::addHeader(gfname, "Case_Freq_Genotype11");
	opts::addHeader(gfname, "Case_Freq_Genotype12");
	opts::addHeader(gfname, "Case_Freq_Genotype22");
	opts::addHeader(gfname, "Case_Count_Genotype11");
	opts::addHeader(gfname, "Case_Count_Genotype12");
	opts::addHeader(gfname, "Case_Count_Genotype22");
	opts::addHeader(gfname, "Control_Freq_Genotype11");
	opts::addHeader(gfname, "Control_Freq_Genotype12");
	opts::addHeader(gfname, "Control_Freq_Genotype22");
	opts::addHeader(gfname, "Control_Count_Genotype11");
	opts::addHeader(gfname, "Control_Count_Genotype12");
	opts::addHeader(gfname, "Control_Count_Genotype22");

    for(int i = 0; i < maxalleles; i++){
		myoutput << "\tAllele" << (i+1);
		opts::addHeader(afname, "Allele" + getString<int>(i+1));
		if(options.doGroupFile()){
			gmyoutput << "\tAllele" << (i+1);
			opts::addHeader(gafname, "Allele" +getString<int>(i+1));
		}
    }
    for(int i = 0; i < maxalleles; i++){
        myoutput << "\t" << "Overall_Allele" << (i+1) << "_freq";
		opts::addHeader(afname, "Overall_Allele" + getString<int>(i+1) + "_freq");
    }
    for(int i = 0; i < maxalleles; i++){
        myoutput << "\t" << "Overall_Allele" << (i+1) << "_count";
		opts::addHeader(afname, "Overall_Allele" + getString<int>(i+1) + "_count");
    }
    for(int i = 0; i < maxalleles; i++){
        myoutput << "\t" << "Case_Allele" << (i+1) << "_freq";
		opts::addHeader(afname, "Case_Allele" + getString<int>(i+1) + "_freq");
    }
    for(int i = 0; i < maxalleles; i++){
        myoutput << "\t" << "Case_Allele" << (i+1) << "_count";
		opts::addHeader(afname, "Case_Allele" + getString<int>(i+1) + "_count");
    }
    for(int i = 0; i < maxalleles; i++){
        myoutput << "\t" << "Control_Allele" << (i+1) << "_freq";
		opts::addHeader(afname, "Control_Allele" + getString<int>(i+1) + "_freq");
    }
    for(int i = 0; i < maxalleles; i++){
        myoutput << "\t" << "Control_Allele" << (i+1) << "_count";
		opts::addHeader(afname, "Control_Allele" + getString<int>(i+1) + "_count");
    }

	if(options.doGroupFile()){
		map<string, vector<Sample*> >::iterator giter;
		map<string, vector<Sample*> > groups = options.getGroups();
		for(giter = groups.begin(); giter != groups.end(); giter++){
			string mygroup = giter->first;
    		for(int i = 0; i < maxalleles; i++){
        		gmyoutput << "\t" << mygroup << "_Allele" << (i+1) << "_freq";
				opts::addHeader(gafname, mygroup + "_Allele" + getString<int>(i+1) + "_freq");
    		}
    		for(int i = 0; i < maxalleles; i++){
        		gmyoutput << "\t" << mygroup << "_Allele" << (i+1) << "_count";
				opts::addHeader(gafname, mygroup + "_Allele" + getString<int>(i+1) + "_count");
    		}
			gmygeno << "\t" << mygroup << "_Freq_Genotype11\t" << mygroup << "_Freq_Genotype12\t" << mygroup << "_Freq_Genotype22\t" << mygroup << "_Count_Genotype11\t" << mygroup << "_Count_Genotype12\t" << mygroup << "_Count_Genotype22";
			opts::addHeader(ggfname, mygroup + "_Freq_Genotype11");
			opts::addHeader(ggfname, mygroup + "_Freq_Genotype12");
			opts::addHeader(ggfname, mygroup + "_Freq_Genotype22");
			opts::addHeader(ggfname, mygroup + "_Count_Genotype11");
			opts::addHeader(ggfname, mygroup + "_Count_Genotype12");
			opts::addHeader(ggfname, mygroup + "_Count_Genotype22");
		}
		gmygeno << endl;
		gmyoutput << endl;
	}
	mygeno << endl;
    myoutput << endl;


    if(options.doParental()){
        paren.open(pfname.c_str(), ios::out);
        pareng.open(pgfname.c_str(), ios::out);
        if(!paren){
	        opts::printLog("Error opening: " + pfname + ". Exiting!\n");
            //exit(1);
	        throw MethodException("Error opening: " + pfname + ".  Exiting!\n");
        }
        if(!pareng){
	        opts::printLog("Error opening: " + pgfname + ". Exiting!\n");
	        //exit(1);
	        throw MethodException("Error opening: " + pgfname + ".  Exiting!\n");
	    }
		opts::addFile("Marker", stepname, pfname);
		opts::addFile("Marker", stepname, pgfname);
	    paren.precision(4);
	    pareng.precision(4);
		if((*markers).at(0)->getDetailHeaders().size() > 0){
        	paren << "Chrom\trsID\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders();
        	pareng << "Chrom\trsID\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders() << "\tGenotype11\tGenotype12\tGenotype22\tParent_Male_Freq_Genotype11\tParent_Male_Freq_Genotype12\tParent_Male_Freq_Genotype22\tParent_Male_Count_Genotype11\tParent_Male_Count_Genotype12\tParent_Male_Count_Genotype22\tParent_Female_Freq_Genotype11\tParent_Female_Freq_Genotype12\tParent_Female_Freq_Genotype22\tParent_Female_Count_Genotype11\tParent_Female_Count_Genotype12\tParent_Female_Count_Genotype22\n";
		}
		else{
			paren << "Chrom\trsID\tProbeID\tbploc";
        	pareng << "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22\tParent_Male_Freq_Genotype11\tParent_Male_Freq_Genotype12\tParent_Male_Freq_Genotype22\tParent_Male_Count_Genotype11\tParent_Male_Count_Genotype12\tParent_Male_Count_Genotype22\tParent_Female_Freq_Genotype11\tParent_Female_Freq_Genotype12\tParent_Female_Freq_Genotype22\tParent_Female_Count_Genotype11\tParent_Female_Count_Genotype12\tParent_Female_Count_Genotype22\n";
		}
		opts::addHeader(pgfname, "Parent_Male_Freq_Genotype11");
		opts::addHeader(pgfname, "Parent_Male_Freq_Genotype12");
		opts::addHeader(pgfname, "Parent_Male_Freq_Genotype22");
		opts::addHeader(pgfname, "Parent_Male_Count_Genotype11");
		opts::addHeader(pgfname, "Parent_Male_Count_Genotype12");
		opts::addHeader(pgfname, "Parent_Male_Count_Genotype22");
		opts::addHeader(pgfname, "Parent_Female_Freq_Genotype11");
		opts::addHeader(pgfname, "Parent_Female_Freq_Genotype12");
		opts::addHeader(pgfname, "Parent_Female_Freq_Genotype22");
		opts::addHeader(pgfname, "Parent_Female_Count_Genotype11");
		opts::addHeader(pgfname, "Parent_Female_Count_Genotype12");
		opts::addHeader(pgfname, "Parent_Female_Count_Genotype22");


        for(int i = 0; i < maxalleles; i++){
	       	paren << "\t" << "Allele" << (i+1);
	    }
	    for(int i = 0; i < maxalleles; i++){
		    paren << "\t" << "Parent_Male_Allele" << (i+1) << "_freq";
			opts::addHeader(pfname, "Parent_Male_Allele" + getString<int>(i+1) + "_freq");
		}
		for(int i = 0; i < maxalleles; i++){
			paren << "\t" << "Parent_Male_Allele" << (i+1) << "_count";
			opts::addHeader(pfname, "Parent_Male_Allele" + getString<int>(i+1) + "_count");
		}
		for(int i = 0; i < maxalleles; i++){
            paren << "\t" << "Parent_Female_Allele" << (i+1) << "_freq";
			opts::addHeader(pfname, "Parent_Female_Allele" + getString<int>(i+1) + "_freq");
        }
        for(int i = 0; i < maxalleles; i++){
            paren << "\t" << "Parent_Female_Allele" << (i+1) << "_count";
			opts::addHeader(pfname, "Parent_Female_Allele" + getString<int>(i+1) + "_count");
        }
        paren << endl;
    }
    if(options.doGender()){
        gend.open(gendfname.c_str(), ios::out);
        gendg.open(gendgfname.c_str(), ios::out);
		opts::addFile("Marker", stepname, gendfname);
		opts::addFile("Marker", stepname, gendgfname);
        if(!gend){
	        opts::printLog("Error opening: " + gendfname + ". Exiting!\n");
	        //exit(1);
	        throw MethodException("Error opening: " + gendfname + ".  Exiting!\n");
	    }
	    if(!gendg){
		    opts::printLog("Error opening: " + gendgfname + ". Exiting!\n");
		    //exit(1);
		    throw MethodException("Error opening: " + gendgfname + ".  Exiting!\n");
		}
		gend.precision(4);
		gendg.precision(4);
		if((*markers).at(0)->getDetailHeaders().size() > 0){
        	gend << "Chrom\trsID\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders();
        	gendg << "Chrom\trsID\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders() << "\tGenotype11\tGenotype12\tGenotype22\tOverall_Male_Freq_Genotype11\tOverall_Male_Freq_Genotype12\tOverall_Male_Freq_Genotype22\tOverall_Male_Count_Genotype11\tOverall_Male_Count_Genotype12\tOverall_Male_Count_Genotype22\tOverall_Female_Freq_Genotype11\tOverall_Female_Freq_Genotype12\tOverall_Female_Freq_Genotype22\tOverall_Female_Count_Genotype11\tOverall_Female_Count_Genotype12\tOverall_Female_Count_Genotype22\n";
		}
		else{
        	gend << "Chrom\trsID\tProbeID\tbploc";
        	gendg << "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22\tOverall_Male_Freq_Genotype11\tOverall_Male_Freq_Genotype12\tOverall_Male_Freq_Genotype22\tOverall_Male_Count_Genotype11\tOverall_Male_Count_Genotype12\tOverall_Male_Count_Genotype22\tOverall_Female_Freq_Genotype11\tOverall_Female_Freq_Genotype12\tOverall_Female_Freq_Genotype22\tOverall_Female_Count_Genotype11\tOverall_Female_Count_Genotype12\tOverall_Female_Count_Genotype22\n";
		}
		opts::addHeader(gendgfname, "Overall_Male_Freq_Genotype11");
		opts::addHeader(gendgfname, "Overall_Male_Freq_Genotype12");
		opts::addHeader(gendgfname, "Overall_Male_Freq_Genotype22");
		opts::addHeader(gendgfname, "Overall_Male_Count_Genotype11");
		opts::addHeader(gendgfname, "Overall_Male_Count_Genotype12");
		opts::addHeader(gendgfname, "Overall_Male_Count_Genotype22");
		opts::addHeader(gendgfname, "Overall_Female_Freq_Genotype11");
		opts::addHeader(gendgfname, "Overall_Female_Freq_Genotype12");
		opts::addHeader(gendgfname, "Overall_Female_Freq_Genotype22");
		opts::addHeader(gendgfname, "Overall_Female_Count_Genotype11");
		opts::addHeader(gendgfname, "Overall_Female_Count_Genotype12");
		opts::addHeader(gendgfname, "Overall_Female_Count_Genotype22");

        for(int i = 0; i < maxalleles; i++){
	        gend << "\t" << "Allele" << (i+1);
	    }
        for(int i = 0; i < maxalleles; i++){
            gend << "\t" << "Overall_Male_Allele" << (i+1) << "_freq";
			opts::addHeader(gendfname, "Overall_Male_Allele" + getString<int>(i+1) + "_freq");
        }
        for(int i = 0; i < maxalleles; i++){
            gend << "\t" << "Overall_Male_Allele" << (i+1) << "_count";
			opts::addHeader(gendfname, "Overall_Male_Allele" + getString<int>(i+1) + "_count");
        }
        for(int i = 0; i < maxalleles; i++){
            gend << "\t" << "Overall_Female_Allele" << (i+1) << "_freq";
			opts::addHeader(gendfname, "Overall_Female_Allele" + getString<int>(i+1) + "_freq");
        }
        for(int i = 0; i < maxalleles; i++){
            gend << "\t" << "Overall_Female_Allele" << (i+1) << "_count";
			opts::addHeader(gendfname, "Overall_Female_Allele" + getString<int>(i+1) + "_count");
        }
        gend << endl;
    }
    if(options.doCaseControl()){
        cc.open(ccfname.c_str(), ios::out);
        ccg.open(ccgfname.c_str(), ios::out);
		opts::addFile("Marker", stepname, ccfname);
		opts::addFile("Marker", stepname, ccgfname);
        if(!cc){
            opts::printLog("Error opening: " + ccfname + ". Exiting!\n");
            throw MethodException("Error opening: " + ccfname + ".  Exiting!\n");
        }
        if(!ccg){
            opts::printLog("Error opening: " + ccgfname + ". Exiting!\n");
            throw MethodException("Error opening: " + ccgfname + ".  Exiting!\n");
        }
        cc.precision(4);
        ccg.precision(4);

		if((*markers).at(0)->getDetailHeaders().size() > 0){
        	cc << "Chrom\trsID\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders();
        	ccg << "Chrom\trsID\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders() << "\tGenotype11\tGenotype12\tGenotype22\tCase_Male_Freq_Genotype11\tCase_Male_Freq_Genotype12\tCase_Male_Freq_Genotype22\tCase_Male_Count_Genotype11\tCase_Male_Count_Genotype12\tCase_Male_Count_Genotype22\tCase_Female_Freq_Genotype11\tCase_Female_Freq_Genotype12\tCase_Female_Freq_Genotype22\tCase_Female_Count_Genotype11\tCase_Female_Count_Genotype12\tCase_Female_Count_Genotype22\tControl_Male_Freq_Genotype11\tControl_Male_Freq_Genotype12\tControl_Male_Freq_Genotype22\tControl_Male_Count_Genotype11\tControl_Male_Count_Genotype12\tControl_Male_Count_Genotype22\tControl_Female_Freq_Genotype11\tControl_Female_Freq_Genotype12\tControl_Female_Freq_Genotype22\tControl_Female_Count_Genotype11\tControl_Female_Count_Genotype12\tControl_Female_Count_Genotype22\n";
		}
		else{
        	cc << "Chrom\trsID\tProbeID\tbploc";
        	ccg << "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22\tCase_Male_Freq_Genotype11\tCase_Male_Freq_Genotype12\tCase_Male_Freq_Genotype22\tCase_Male_Count_Genotype11\tCase_Male_Count_Genotype12\tCase_Male_Count_Genotype22\tCase_Female_Freq_Genotype11\tCase_Female_Freq_Genotype12\tCase_Female_Freq_Genotype22\tCase_Female_Count_Genotype11\tCase_Female_Count_Genotype12\tCase_Female_Count_Genotype22\tControl_Male_Freq_Genotype11\tControl_Male_Freq_Genotype12\tControl_Male_Freq_Genotype22\tControl_Male_Count_Genotype11\tControl_Male_Count_Genotype12\tControl_Male_Count_Genotype22\tControl_Female_Freq_Genotype11\tControl_Female_Freq_Genotype12\tControl_Female_Freq_Genotype22\tControl_Female_Count_Genotype11\tControl_Female_Count_Genotype12\tControl_Female_Count_Genotype22\n";
		}
		opts::addHeader(ccgfname, "Case_Male_Freq_Genotype11");
		opts::addHeader(ccgfname, "Case_Male_Freq_Genotype12");
		opts::addHeader(ccgfname, "Case_Male_Freq_Genotype22");
		opts::addHeader(ccgfname, "Case_Male_Count_Genotype11");
		opts::addHeader(ccgfname, "Case_Male_Count_Genotype12");
		opts::addHeader(ccgfname, "Case_Male_Count_Genotype22");
		opts::addHeader(ccgfname, "Case_Female_Freq_Genotype11");
		opts::addHeader(ccgfname, "Case_Female_Freq_Genotype12");
		opts::addHeader(ccgfname, "Case_Female_Freq_Genotype22");
		opts::addHeader(ccgfname, "Case_Female_Count_Genotype11");
		opts::addHeader(ccgfname, "Case_Female_Count_Genotype12");
		opts::addHeader(ccgfname, "Case_Female_Count_Genotype22");
		opts::addHeader(ccgfname, "Control_Male_Freq_Genotype11");
		opts::addHeader(ccgfname, "Control_Male_Freq_Genotype12");
		opts::addHeader(ccgfname, "Control_Male_Freq_Genotype22");
		opts::addHeader(ccgfname, "Control_Male_Count_Genotype11");
		opts::addHeader(ccgfname, "Control_Male_Count_Genotype12");
		opts::addHeader(ccgfname, "Control_Male_Count_Genotype22");
		opts::addHeader(ccgfname, "Control_Female_Freq_Genotype11");
		opts::addHeader(ccgfname, "Control_Female_Freq_Genotype12");
		opts::addHeader(ccgfname, "Control_Female_Freq_Genotype22");
		opts::addHeader(ccgfname, "Control_Female_Count_Genotype11");
		opts::addHeader(ccgfname, "Control_Female_Count_Genotype12");
		opts::addHeader(ccgfname, "Control_Female_Count_Genotype22");

        for(int i = 0; i < maxalleles; i++){
            cc << "\t" << "Allele" << (i+1);
        }
        for(int i = 0; i < maxalleles; i++){
            cc << "\t" << "Case_Male_Allele" << (i+1) << "_freq";
			opts::addHeader(ccfname, "Case_Male_Allele" + getString<int>(i+1) + "_freq");
        }
        for(int i = 0; i < maxalleles; i++){
            cc << "\t" << "Case_Male_Allele" << (i+1) << "_count";
			opts::addHeader(ccfname, "Case_Male_Allele" + getString<int>(i+1) + "_count");
        }
        for(int i = 0; i < maxalleles; i++){
            cc << "\t" << "Case_Female_Allele" << (i+1) << "_freq";
			opts::addHeader(ccfname, "Case_Female_Allele" + getString<int>(i+1) + "_freq");
        }
        for(int i = 0; i < maxalleles; i++){
            cc << "\t" << "Case_Female_Allele" << (i+1) << "_count";
			opts::addHeader(ccfname, "Case_Female_Allele" + getString<int>(i+1) + "_count");
        }
        for(int i = 0; i < maxalleles; i++){
            cc << "\t" << "Control_Male_Allele" << (i+1) << "_freq";
			opts::addHeader(ccfname, "Control_Male_Allele" + getString<int>(i+1) + "_freq");
        }
        for(int i = 0; i < maxalleles; i++){
            cc << "\t" << "Control_Male_Allele" << (i+1) << "_count";
			opts::addHeader(ccfname, "Control_Male_Allele" + getString<int>(i+1) + "_count");
        }
        for(int i = 0; i < maxalleles; i++){
            cc << "\t" << "Control_Female_Allele" << (i+1) << "_freq";
			opts::addHeader(ccfname, "Control_Female_Allele" + getString<int>(i+1) + "_freq");
        }
	    for(int i = 0; i < maxalleles; i++){
	        cc << "\t" << "Control_Female_Allele" << (i+1) << "_count";
			opts::addHeader(ccfname, "Control_Female_Allele" + getString<int>(i+1) + "_count");
	    }
	    cc << endl;
	}

	//begin processing
	int prev_base = 0;
	int prev_chrom = -1;
	for(int k = 0; k < msize; k++){
		if((*markers).at(k)->isEnabled() && Helpers::isValidMarker((*markers).at(k), &options, prev_base, prev_chrom)){
			//flag markers to be ignored

			//perform calculations
			calcOne((*markers).at(k));
			if(options.doGroupFile()){
				calcOneGroups((*markers).at(k));
			}

			if((*markers).at(k)->isMicroSat()){
				myoutput << (*markers).at(k)->toString();
				if(options.doGroupFile()){
					gmyoutput << (*markers).at(k)->toString();
				}
				if(options.doParental()){
					paren << (*markers).at(k)->toString();
				}
				if(options.doGender()){
					gend << (*markers).at(k)->toString();
				}
				if(options.doCaseControl()){
					cc << (*markers).at(k)->toString();
				}
                int total_o = 0;
                int total_ca = 0;
                int total_con = 0;
                int numalleles = (*markers).at(k)->getNumAlleles();
                for(int a = 0; a < numalleles; a++){
					if(useoverall){
						total_o += m_allele_counts_o.at(a);
					}
					else{
	                	total_o += m_allele_counts_p.at(a);
					}
	                total_ca += m_allele_counts_ca.at(a);
	                total_con += m_allele_counts_con.at(a);
	            }
	            for(int a = 0; a < numalleles; a++){
		            myoutput << "\t" << (*markers).at(k)->getAllele(a);
					if(options.doGroupFile()){
						gmyoutput << "\t"  << (*markers).at(k)->getAllele(a);
					}
		            if(options.doParental()){
			           	paren << "\t" << (*markers).at(k)->getAllele(a);
			        }
			        if(options.doGender()){
				       	gend << "\t" << (*markers).at(k)->getAllele(a);
				    }
				    if(options.doCaseControl()){
				       	cc << "\t" << (*markers).at(k)->getAllele(a);
				    }
                }
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						myoutput << "\tNA";
						if(options.doGroupFile()){
							gmyoutput << "\tNA";
						}
						if(options.doParental()){
							paren << "\tNA";
						}
						if(options.doGender()){
							gend << "\tNA";
						}
						if(options.doCaseControl()){
							cc << "\tNA";
						}
					}
				}
                //overall
                for(int a = 0; a < numalleles; a++){
					float freq = 0.0f;
					if(useoverall){
						freq = ((float) m_allele_counts_o.at(a) / (float) total_o);
					}
					else{
						freq = ((float) m_allele_counts_p.at(a) / (float)total_o);
					}
					myoutput << "\t" << freq;
				}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						myoutput << "\tNA";
					}
				}
				for(int a = 0; a < numalleles; a++){
                    if(useoverall){
						myoutput << "\t" << m_allele_counts_o.at(a);
					}
					else{
						myoutput << "\t" << m_allele_counts_p.at(a);
					}
				}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						myoutput << "\tNA";
					}
				}
				//case
				for(int a = 0; a < numalleles; a++){
                    float freq = ((float) m_allele_counts_ca.at(a) / (float)total_ca);
                    myoutput << "\t" << freq;
                }
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						myoutput << "\tNA";
					}
				}
                for(int a = 0; a < numalleles; a++){
                    myoutput << "\t" << m_allele_counts_ca.at(a);
                }
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						myoutput << "\tNA";
					}
				}
				//control
				for(int a = 0; a < numalleles; a++){
                    float freq = ((float) m_allele_counts_con.at(a) / (float)total_con);
                    myoutput << "\t" << freq;
                }
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						myoutput << "\tNA";
					}
				}
                for(int a = 0; a < numalleles; a++){
                    myoutput << "\t" << m_allele_counts_con.at(a);
                }
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						myoutput << "\tNA";
					}
				}
				//groups?
				mygeno << (*markers).at(k)->toString();
				if(options.doGroupFile()){
					gmygeno << (*markers).at(k)->toString();
					int gm_total = 0;
					map<string, vector<Sample*> > groups = options.getGroups();
					map<string, vector<Sample*> >::iterator giter;
					for(giter = groups.begin(); giter != groups.end(); giter++){
						string mygroup = giter->first;
						for(int a = 0; a < numalleles; a++){
							gm_total += gm_allele_counts_o[mygroup][a];
						}
						for(int a = 0; a < numalleles; a++){
							float freq = ((float) gm_allele_counts_o[mygroup][a] / (float)gm_total);
							gmyoutput << "\t" << freq;
						}
						if(maxalleles > numalleles){
							for(int b = 0; b < (maxalleles - numalleles); b++){
								gmyoutput << "\tNA";
							}
						}
						for(int a = 0; a < numalleles; a++){
							gmyoutput << "\t" << gm_allele_counts_o[mygroup][a];
						}
						if(maxalleles > numalleles){
							for(int b = 0; b < (maxalleles - numalleles); b++){
								gmyoutput << "\tNA";
							}
						}
						gmygeno << "\tNA\tNA\tNA\tNA\tNA\tNA";
					}
					gmyoutput << endl;
				}
				myoutput << endl;

				for(int l = 0; l < 21; l++){
					mygeno << "\tNA";
					if(options.doGroupFile()){
						gmygeno << "\tNA";
					}
				}

				mygeno << endl;
				if(options.doGroupFile()){
					gmygeno << endl;
				}

                if(options.doParental()){
					pareng << (*markers).at(k)->toString();
                    int total_pm = 0;
                    int total_pf = 0;
                    for(int a = 0; a < numalleles; a++){
	                    total_pm += m_allele_counts_pm.at(a);
                        total_pf += m_allele_counts_pf.at(a);
					}
					//Male
					for(int a = 0; a < numalleles; a++){
						float freq = ((float) m_allele_counts_pm.at(a) / (float)total_pm);
						paren << "\t" << freq;
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						paren << "\tNA";
					}
				}
					for(int a = 0; a < numalleles; a++){
						paren << "\t" << m_allele_counts_pm.at(a);
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						paren << "\tNA";
					}
				}
					//Female
					for(int a = 0; a < numalleles; a++){
						float freq = ((float) m_allele_counts_pf.at(a) / (float)total_pf);
						paren << "\t" << freq;
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						paren << "\tNA";
					}
				}
					for(int a = 0; a < numalleles; a++){
						paren << "\t" << m_allele_counts_pf.at(a);
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						paren << "\tNA";
					}
				}
					paren << endl;

					for(int l = 0; l < 15; l++){
						pareng << "\tNA";
					}
					pareng << endl;

				}
				if(options.doGender()){
					gendg << (*markers).at(k)->toString();
					int total_pm = 0;
					int total_pf = 0;
					for(int a = 0; a < numalleles; a++){
						if(useoverall){
							total_pm += m_allele_counts_om.at(a);
							total_pf += m_allele_counts_of.at(a);
						}
						else{
							total_pm += m_allele_counts_pm.at(a);
							total_pf += m_allele_counts_pf.at(a);
						}
					}
					//Male
					for(int a = 0; a < numalleles; a++){
						float freq = 0.0f;
						if(useoverall){
							freq = ((float) m_allele_counts_om.at(a) / (float)total_pm);
						}
						else{
							freq = ((float) m_allele_counts_pm.at(a) / (float)total_pm);
						}
						gend << "\t" << freq;
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						gend << "\tNA";
					}
				}
					for(int a = 0; a < numalleles; a++){
						if(useoverall){
							gend << "\t" << m_allele_counts_om.at(a);
						}
						else{
							gend << "\t" << m_allele_counts_pm.at(a);
						}
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						gend << "\tNA";
					}
				}
					//Female
					for(int a = 0; a < numalleles; a++){
						float freq = 0.0f;
						if(useoverall){
							freq = ((float) m_allele_counts_of.at(a) / (float)total_pf);
						}
						else{
							freq = ((float) m_allele_counts_pf.at(a) / (float)total_pf);
						}
						gend << "\t" << freq;
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						gend << "\tNA";
					}
				}
					for(int a = 0; a < numalleles; a++){
						if(useoverall){
							gend << "\t" << m_allele_counts_of.at(a);
						}
						else{
							gend << "\t" << m_allele_counts_pf.at(a);
						}
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						gend << "\tNA";
					}
				}
					gend << endl;

					for(int l = 0; l < 15; l++){
						gendg << "\tNA";
					}
					gendg << endl;

				}
				if(options.doCaseControl()){
					ccg << (*markers).at(k)->toString();
					int total_cam = 0;
					int total_caf = 0;
					int total_conm = 0;
					int total_conf = 0;
					for(int a = 0; a < numalleles; a++){
						total_cam += m_allele_counts_cam.at(a);
						total_caf += m_allele_counts_caf.at(a);
						total_conm += m_allele_counts_conm.at(a);
						total_conf += m_allele_counts_conf.at(a);
					}
					//Case Male
					for(int a = 0; a < numalleles; a++){
						float freq = ((float) m_allele_counts_cam.at(a) / (float)total_cam);
						cc << "\t" << freq;
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						cc << "\tNA";
					}
				}
					for(int a = 0; a < numalleles; a++){
						cc << "\t" << m_allele_counts_cam.at(a);
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						cc << "\tNA";
					}
				}
					//Case Female
					for(int a = 0; a < numalleles; a++){
						float freq = ((float) m_allele_counts_caf.at(a) / (float)total_caf);
						cc << "\t" << freq;
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						cc << "\tNA";
					}
				}
					for(int a = 0; a < numalleles; a++){
						cc << "\t" << m_allele_counts_caf.at(a);
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						cc << "\tNA";
					}
				}
					//Control Male
					for(int a = 0; a < numalleles; a++){
						float freq = ((float) m_allele_counts_conm.at(a) / (float)total_conm);
						cc << "\t" << freq;
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						cc << "\tNA";
					}
				}
					for(int a = 0; a < numalleles; a++){
						cc << "\t" << m_allele_counts_conm.at(a);
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						cc << "\tNA";
					}
				}
					//Control Female
					for(int a = 0; a < numalleles; a++){
						float freq = ((float) m_allele_counts_conf.at(a) / (float)total_conf);
						cc << "\t" << freq;
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						cc << "\tNA";
					}
				}
					for(int a = 0; a < numalleles; a++){
						cc << "\t" << m_allele_counts_conf.at(a);
					}
				if(maxalleles > numalleles){
					for(int b = 0; b < (maxalleles - numalleles); b++){
						cc << "\tNA";
					}
				}
					cc << endl;

					for(int l = 0; l < 27; l++){
						ccg << "\tNA";
					}
					ccg << endl;

				}
			}
			else{
            	myoutput << (*markers).at(k)->toString() << "\t" << (*markers).at(k)->getAllele1() << "\t" << (*markers).at(k)->getAllele2();
				if(options.doGroupFile()){
            		gmyoutput << (*markers).at(k)->toString() << "\t" << (*markers).at(k)->getAllele1() << "\t" << (*markers).at(k)->getAllele2();
				}
				for(int l = 2; l < maxalleles; l++){
					myoutput << "\tNA";
					if(options.doGroupFile()){
						gmyoutput << "\tNA";
					}
				}

            	//overall
            	float majfreq = 0.0f;
				if(useoverall){
					majfreq = ((float) a1_count / (float)(a1_count + a2_count));
				}
				else{
					majfreq = ((float)a1_countP/(float)(a1_countP + a2_countP));
				}
            	float minfreq = 1.0f - majfreq;
            	myoutput << "\t" << majfreq << "\t" << minfreq;
				for(int l = 2; l < maxalleles; l++){
					myoutput << "\tNA";
				}
				if(useoverall){
					myoutput << "\t" << a1_count << "\t" << a2_count;
				}
				else{
					myoutput << "\t" << a1_countP << "\t" << a2_countP;
				}
				for(int t = 2; t < maxalleles; t++){
					myoutput << "\tNA";
				}
				//case
				majfreq = ((float)a1_countCa/(float)(a1_countCa + a2_countCa));
				minfreq = 1.0f - majfreq;
				myoutput << "\t" << majfreq << "\t" << minfreq;
				for(int l = 2; l < maxalleles; l++){
					myoutput << "\tNA";
				}
				myoutput << "\t" << a1_countCa << "\t" << a2_countCa;
				for(int t = 2; t < maxalleles; t++){
	 				myoutput << "\tNA";
				}
				//control
				majfreq = ((float)a1_countCon/(float)(a1_countCon + a2_countCon));
				minfreq = 1.0f - majfreq;
				myoutput << "\t" << majfreq << "\t" << minfreq;
				for(int t = 2; t < maxalleles; t++){
					 myoutput << "\tNA";
				}
				myoutput << "\t" << a1_countCon << "\t" << a2_countCon;
				for(int t = 2; t < maxalleles; t++){
					myoutput << "\tNA";
				}

				//groups?
				if(options.doGroupFile()){
					int gm_total = 0;
					map<string, vector<Sample*> > groups = options.getGroups();
					map<string, vector<Sample*> >::iterator giter;
					for(giter = groups.begin(); giter != groups.end(); giter++){
						string mygroup = giter->first;
						gm_total = ga1_count[mygroup] + ga2_count[mygroup];
						float freq = ((float) ga1_count[mygroup] / (float)gm_total);
						float freq2 = 1.0f - freq;
						gmyoutput << "\t" << freq << "\t" << freq2;
						for(int b = 2; b < maxalleles; b++){
							gmyoutput << "\tNA";
						}
						gmyoutput << "\t" << ga1_count[mygroup] << "\t" << ga2_count[mygroup];
						for(int b = 2; b < maxalleles; b++){
							gmyoutput << "\tNA";
						}
					}
					gmyoutput << endl;
				}

				myoutput << endl;

				mygeno << (*markers).at(k)->toString() << "\t" << (*markers).at(k)->getAllele1() << "_" << (*markers).at(k)->getAllele1() << "\t" << (*markers).at(k)->getAllele1() << "_" << (*markers).at(k)->getAllele2() << "\t" << (*markers).at(k)->getAllele2() << "_" << (*markers).at(k)->getAllele2();
				if(options.doGroupFile()){
					gmygeno << (*markers).at(k)->toString() << "\t" << (*markers).at(k)->getAllele1() << "_" << (*markers).at(k)->getAllele1() << "\t" << (*markers).at(k)->getAllele1() << "_" << (*markers).at(k)->getAllele2() << "\t" << (*markers).at(k)->getAllele2() << "_" << (*markers).at(k)->getAllele2();
				}
				//overall
				float freq1 = 0.0f;
				float freq2 = 0.0f;
				float freq3 = 0.0f;
				if(useoverall){
					freq1 = ((float) a1_homo_count)/(a1_homo_count + a12_count + a2_homo_count);
					freq2 = ((float) a12_count)/(a1_homo_count + a12_count + a2_homo_count);
					freq3 = ((float) a2_homo_count)/(a1_homo_count + a12_count + a2_homo_count);
					mygeno << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_count << "\t"
						<< a12_count << "\t"
						<< a2_homo_count << "\t";
				}
				else{
					freq1 = ((float) a1_homo_countP)/(a1_homo_countP + a12_countP + a2_homo_countP);
					freq2 = ((float) a12_countP)/(a1_homo_countP + a12_countP + a2_homo_countP);
					freq3 = ((float) a2_homo_countP)/(a1_homo_countP + a12_countP + a2_homo_countP);
					mygeno << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countP << "\t"
						<< a12_countP << "\t"
						<< a2_homo_countP << "\t";
				}

				//case overall
				freq1 = ((float) a1_homo_countCa)/(a1_homo_countCa + a12_countCa + a2_homo_countCa);
				freq2 = ((float) a12_countCa)/(a1_homo_countCa + a12_countCa + a2_homo_countCa);
				freq3 = ((float) a2_homo_countCa)/(a1_homo_countCa + a12_countCa + a2_homo_countCa);
				mygeno << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countCa << "\t"
					<< a12_countCa << "\t"
					<< a2_homo_countCa << "\t";
				//control overall
				freq1 = ((float) a1_homo_countCon)/(a1_homo_countCon + a12_countCon + a2_homo_countCon);
				freq2 = ((float) a12_countCon)/(a1_homo_countCon + a12_countCon + a2_homo_countCon);
				freq3 = ((float) a2_homo_countCon)/(a1_homo_countCon + a12_countCon + a2_homo_countCon);
				mygeno << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countCon << "\t"
					<< a12_countCon << "\t"
					<< a2_homo_countCon;

				//groups?
				if(options.doGroupFile()){
					map<string, vector<Sample*> > groups = options.getGroups();
					map<string, vector<Sample*> >::iterator giter;
					for(giter = groups.begin(); giter != groups.end(); giter++){
						string mygroup = giter->first;
						float genotot = ga1_homo_count[mygroup] + ga12_count[mygroup] + ga2_homo_count[mygroup];
						freq1 = ((float) ga1_homo_count[mygroup]/genotot);
						freq2 = ((float) ga12_count[mygroup]/genotot);
						freq3 = ((float) ga2_homo_count[mygroup]/genotot);
						gmygeno << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << ga1_homo_count[mygroup] << "\t" << ga12_count[mygroup] << "\t" << ga2_homo_count[mygroup];
					}
					gmygeno << endl;
				}
				mygeno << endl;


				if(options.doParental()){
					paren << (*markers).at(k)->toString() << "\t" << (*markers).at(k)->getAllele1() << "\t" << (*markers).at(k)->getAllele2();
					for(int l = 2; l < maxalleles; l++){
						paren << "\tNA";
					}
					//parent male
					float majfreq = ((float)a1_countPM/(float)(a1_countPM + a2_countPM));
					float minfreq = 1.0f - majfreq;
					paren << "\t" << majfreq << "\t" << minfreq;
					for(int t = 2; t < maxalleles; t++){
						paren << "\tNA";
					}
					paren << "\t" << a1_countPM << "\t" << a2_countPM;
					for(int t = 2; t < maxalleles; t++){
						paren << "\tNA";
					}
					//parent female
					majfreq = ((float)a1_countPF/(float)(a1_countPF + a2_countPF));
					minfreq = 1.0f - majfreq;
					paren << "\t" << majfreq << "\t" << minfreq;
					for(int t = 2; t < maxalleles; t++){
						paren << "\tNA";
					}
					paren << "\t" << a1_countPF << "\t" << a2_countPF;
					for(int t = 2; t < maxalleles; t++){
						paren << "\tNA";
					}
					paren << endl;

					pareng << (*markers).at(k)->toString() << "\t" << (*markers).at(k)->getAllele1() << "_" << (*markers).at(k)->getAllele1() << "\t" << (*markers).at(k)->getAllele1() << "_" << (*markers).at(k)->getAllele2() << "\t" << (*markers).at(k)->getAllele2() << "_" << (*markers).at(k)->getAllele2();
					float freq1 = ((float) a1_homo_countPM)/(a1_homo_countPM + a12_countPM + a1_homo_countPM);
					float freq2 = ((float) a12_countPM)/(a1_homo_countPM + a12_countPM + a2_homo_countPM);
					float freq3 = ((float) a2_homo_countPM)/(a1_homo_countPM + a12_countPM + a2_homo_countPM);
					pareng << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countPM << "\t"
						<< a12_countPM << "\t"
						<< a2_homo_countPM;
					freq1 = ((float) a1_homo_countPF)/(a1_homo_countPF + a12_countPF + a2_homo_countPF);
					freq2 = ((float) a12_countPF)/(a1_homo_countPF + a12_countPF + a2_homo_countPF);
					freq3 = ((float) a2_homo_countPF)/(a1_homo_countPF + a12_countPF + a2_homo_countPF);
					pareng << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countPF << "\t"
						<< a12_countPF << "\t"
						<< a2_homo_countPF;
					pareng << endl;
				}
				if(options.doGender()){
					//overall male
					gend << (*markers).at(k)->toString() << "\t" << (*markers).at(k)->getAllele1() << "\t" << (*markers).at(k)->getAllele2();
					for(int l = 2; l < maxalleles; l++){
						gend << "\tNA";
					}
					float majfreq = 0.0f;
					if(useoverall){
						majfreq = ((float)a1_countM/(float)(a1_countM + a2_countM));
					}
					else{
						majfreq = ((float)a1_countPM/(float)(a1_countPM + a2_countPM));
					}
					float minfreq = 1.0f - majfreq;
					gend << "\t" << majfreq << "\t" << minfreq;
					for(int t = 2; t < maxalleles; t++){
						gend << "\tNA";
					}
					if(useoverall){
						gend << "\t" << a1_countM << "\t" << a2_countM;
					}
					else{
						gend << "\t" << a1_countPM << "\t" << a2_countPM;
					}
					for(int t = 2; t < maxalleles; t++){
						gend << "\tNA";
					}
					//overall female
					if(useoverall){
						majfreq = ((float)a1_countF/(float)(a1_countF + a2_countF));
					}
					else{
						majfreq = ((float)a1_countPF/(float)(a1_countPF + a2_countPF));
					}
					minfreq = 1.0f - majfreq;
					gend << "\t" << majfreq << "\t" << minfreq;
					for(int t = 2; t < maxalleles; t++){
						gend << "\tNA";
					}
					if(useoverall){
						gend << "\t" << a1_countF << "\t" << a2_countF;
					}
					else{
						gend << "\t" << a1_countPF << "\t" << a2_countPF;
					}
					for(int t = 2; t < maxalleles; t++){
						gend << "\tNA";
					}
					gend << endl;

					gendg << (*markers).at(k)->toString() << "\t" << (*markers).at(k)->getAllele1() << "_" << (*markers).at(k)->getAllele1() << "\t" << (*markers).at(k)->getAllele1() << "_" << (*markers).at(k)->getAllele2() << "\t" << (*markers).at(k)->getAllele2() << "_" << (*markers).at(k)->getAllele2();
					if(useoverall){
						float freq1 = ((float) a1_homo_countM)/(a1_homo_countM + a12_countM + a2_homo_countM);
						float freq2 = ((float) a12_countM)/(a1_homo_countM + a12_countM + a2_homo_countM);
						float freq3 = ((float) a2_homo_countM)/(a1_homo_countM + a12_countM + a2_homo_countM);
						gendg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countM << "\t"
							<< a12_countM << "\t"
							<< a2_homo_countM;
						freq1 = ((float) a1_homo_countF)/(a1_homo_countF + a12_countF + a2_homo_countF);
						freq2 = ((float) a12_countF)/(a1_homo_countF + a12_countF + a2_homo_countF);
						freq3 = ((float) a2_homo_countF)/(a1_homo_countF + a12_countF + a2_homo_countF);
						gendg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countF << "\t"
							<< a12_countF << "\t"
							<< a2_homo_countF;
					}
					else{
						float freq1 = ((float) a1_homo_countPM)/(a1_homo_countPM + a12_countPM + a2_homo_countPM);
						float freq2 = ((float) a12_countPM)/(a1_homo_countPM + a12_countPM + a2_homo_countPM);
						float freq3 = ((float) a2_homo_countPM)/(a1_homo_countPM + a12_countPM + a2_homo_countPM);
						gendg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countPM << "\t"
							<< a12_countPM << "\t"
							<< a2_homo_countPM;
						freq1 = ((float) a1_homo_countPF)/(a1_homo_countPF + a12_countPF + a2_homo_countPF);
						freq2 = ((float) a12_countPF)/(a1_homo_countPF + a12_countPF + a2_homo_countPF);
						freq3 = ((float) a2_homo_countPF)/(a1_homo_countPF + a12_countPF + a2_homo_countPF);
						gendg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countPF << "\t"
							<< a12_countPF << "\t"
							<< a2_homo_countPF;
					}
					gendg << endl;

				}
				if(options.doCaseControl()){
					//case male
					cc << (*markers).at(k)->toString() << "\t" << (*markers).at(k)->getAllele1() << "\t" << (*markers).at(k)->getAllele2();
					for(int l = 2; l < maxalleles; l++){
						cc << "\tNA";
					}
					float majfreq = ((float)a1_countCaM/(float)(a1_countCaM + a2_countCaM));
					float minfreq = 1.0f - majfreq;
					cc << "\t" << majfreq << "\t" << minfreq;
					for(int t = 2; t < maxalleles; t++){
						cc << "\tNA";
					}
					cc << "\t" << a1_countCaM << "\t" << a2_countCaM;
					for(int t = 2; t < maxalleles; t++){
						cc << "\tNA";
					}
					//case female
					majfreq = ((float)a1_countCaF/(float)(a1_countCaF + a2_countCaF));
					minfreq = 1.0f - majfreq;
					cc << "\t" << majfreq << "\t" << minfreq;
					for(int t = 2; t < maxalleles; t++){
						cc << "\tNA";
					}
					cc << "\t" << a1_countCaF << "\t" << a2_countCaF;
					for(int t = 2; t < maxalleles; t++){
						cc << "\tNA";
					}
					//control male
					majfreq = ((float)a1_countConM/(float)(a1_countConM + a2_countConM));
					minfreq = 1.0f - majfreq;
					cc << "\t" << majfreq << "\t" << minfreq;
					for(int t = 2; t < maxalleles; t++){
						cc << "\tNA";
					}
					cc << "\t" << a1_countConM << "\t" << a2_countConM;
					for(int t = 2; t < maxalleles; t++){
						cc << "\tNA";
					}
					//control female
					majfreq = ((float)a1_countConF/(float)(a1_countConF + a2_countConF));
					minfreq = 1.0f - majfreq;
					cc << "\t" << majfreq << "\t" << minfreq;
					for(int t = 2; t < maxalleles; t++){
						cc << "\tNA";
					}
					cc << "\t" << a1_countConF << "\t" << a2_countConF;
					for(int t = 2; t < maxalleles; t++){
						cc << "\tNA";
					}
					cc << endl;

					ccg << (*markers).at(k)->toString() << "\t" << (*markers).at(k)->getAllele1() << "_" << (*markers).at(k)->getAllele1() << "\t" << (*markers).at(k)->getAllele1() << "_" << (*markers).at(k)->getAllele2() << "\t" << (*markers).at(k)->getAllele2() << "_" << (*markers).at(k)->getAllele2();
					float freq1 = ((float) a1_homo_countCaM)/(a1_homo_countCaM + a12_countCaM + a2_homo_countCaM);
					float freq2 = ((float) a12_countCaM)/(a1_homo_countCaM + a12_countCaM + a2_homo_countCaM);
					float freq3 = ((float) a2_homo_countCaM)/(a1_homo_countCaM + a12_countCaM + a2_homo_countCaM);
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countCaM << "\t"
						<< a12_countCaM << "\t"
						<< a2_homo_countCaM;
					freq1 = ((float) a1_homo_countCaF)/(a1_homo_countCaF + a12_countCaF + a2_homo_countCaF);
					freq2 = ((float) a12_countCaF)/(a1_homo_countCaF + a12_countCaF + a2_homo_countCaF);
					freq3 = ((float) a2_homo_countCaF)/(a1_homo_countCaF + a12_countCaF + a2_homo_countCaF);
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countCaF << "\t"
						<< a12_countCaF << "\t"
						<< a2_homo_countCaF;
					freq1 = ((float) a1_homo_countConM)/(a1_homo_countConM + a12_countConM + a2_homo_countConM);
					freq2 = ((float) a12_countConM)/(a1_homo_countConM + a12_countConM + a2_homo_countConM);
					freq3 = ((float) a2_homo_countConM)/(a1_homo_countConM + a12_countConM + a2_homo_countConM);
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countConM << "\t"
						<< a12_countConM << "\t"
						<< a2_homo_countConM;
					freq1 = ((float) a1_homo_countConF)/(a1_homo_countConF + a12_countConF + a2_homo_countConF);
					freq2 = ((float) a12_countConF)/(a1_homo_countConF + a12_countConF + a2_homo_countConF);
					freq3 = ((float) a2_homo_countConF)/(a1_homo_countConF + a12_countConF + a2_homo_countConF);
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countConF << "\t"
						<< a12_countConF << "\t"
						<< a2_homo_countConF;
					 ccg << endl;

				}
			}

			//filter Markers
			filterOne((*markers).at(k));
		}
	}

}


/*
 *Function: process
 *Description:
 *Main method to begin the whole process.  Flags samples then diverts work to processtest
 *
 *
 */
void AlleleFrequency::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

	sample_flags.resize(samples->size(), false);

	if(options.doGroupFile()){
		options.readGroups(samples);
	}

	if(options.doRandomChild() || options.doAll() || options.doAllChildren() || options.doUnaffSpousesOnly() || options.doUnknownSpouses()){
		flagSamples();
	}
	else{
		options.setFoundersOnly();
		flagSamples();
	}
	processtest();
	return;

}

//Overall
float AlleleFrequency::getAonehomo_exp(){
	float p2 = getAone_freq() * getAone_freq();
	int pop = getPop();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomo_exp(){
	float p2 = getAtwo_freq() * getAtwo_freq();
	int pop = getPop();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHet_exp(){
	float pq = 2 * getAone_freq() * getAtwo_freq();
	int pop = getPop();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAone_freq(){
	float percent = 0.0f;
	int denom = a1_count + a2_count;
	if(denom > 0){
		percent = ((float)a1_count/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwo_freq(){
	float percent = 0.0f;
	int denom = a1_count + a2_count;
	if(denom > 0){
		percent = ((float)a2_count/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getMicroFreq(int l){
	float percent = 0.0f;
	int denom = getMicroDenom();
	if(denom > 0){
		percent = ((float)m_allele_counts_o.at(l)/(float)denom);
	}
	return percent;
}
int AlleleFrequency::getMicroDenom(){
	int total = 0;
	for(int i = 0; i < (int)m_allele_counts_o.size(); i++){
		total += m_allele_counts_o.at(i);
	}
	return total;
}
//Overall Male
float AlleleFrequency::getAonehomoM_exp(){
	float p2 = getAoneM_freq() * getAoneM_freq();
	int pop = getPopM();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoM_exp(){
	float p2 = getAtwoM_freq() * getAtwoM_freq();
	int pop = getPopM();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetM_exp(){
	float pq = 2 * getAoneM_freq() * getAtwoM_freq();
	int pop = getPopM();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneM_freq(){
	float percent = 0.0f;
	int denom = a1_countM + a2_countM;
	if(denom > 0){
		percent = ((float)a1_countM/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoM_freq(){
	float percent = 0.0f;
	int denom = a1_countM + a2_countM;
	if(denom > 0){
		percent = ((float)a2_countM/(float)denom);
	}
	return percent;
}
//Overall Female
float AlleleFrequency::getAonehomoF_exp(){
	float p2 = getAoneF_freq() * getAoneF_freq();
	int pop = getPopF();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoF_exp(){
	float p2 = getAtwoF_freq() * getAtwoF_freq();
	int pop = getPopF();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetF_exp(){
	float pq = 2 * getAoneF_freq() * getAtwoF_freq();
	int pop = getPopF();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneF_freq(){
	float percent = 0.0f;
	int denom = a1_countF + a2_countF;
	if(denom > 0){
		percent = ((float)a1_countF/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoF_freq(){
	float percent = 0.0f;
	int denom = a1_countF + a2_countF;
	if(denom > 0){
		percent = ((float)a2_countF/(float)denom);
	}
	return percent;
}
//Parent
float AlleleFrequency::getAonehomoP_exp(){
	float p2 = getAoneP_freq() * getAoneP_freq();
	int pop = getPopP();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoP_exp(){
	float p2 = getAtwoP_freq() * getAtwoP_freq();
	int pop = getPopP();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetP_exp(){
	float pq = 2 * getAoneP_freq() * getAtwoP_freq();
	int pop = getPopP();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneP_freq(){
	float percent = 0.0f;
	int denom = a1_countP + a2_countP;
	if(denom > 0){
		percent = ((float)a1_countP/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoP_freq(){
	float percent = 0.0f;
	int denom = a1_countP + a2_countP;
	if(denom > 0){
		percent = ((float)a2_countP/(float)denom);
	}
	return percent;
}
//Parent Mom
float AlleleFrequency::getAonehomoPF_exp(){
	float p2 = getAonePF_freq() * getAonePF_freq();
	int pop = getPopPF();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoPF_exp(){
	float p2 = getAtwoPF_freq() * getAtwoPF_freq();
	int pop = getPopPF();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetPF_exp(){
	float pq = 2 * getAonePF_freq() * getAtwoPF_freq();
	int pop = getPopPF();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAonePF_freq(){
	float percent = 0.0f;
	int countone = a1_countPF;
	int counttwo = a2_countPF;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)countone/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoPF_freq(){
	float percent = 0.0f;
	int countone = a1_countPF;
	int counttwo = a2_countPF;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)counttwo/(float)denom);
	}
	return percent;
}

//Parent Dad
float AlleleFrequency::getAonehomoPM_exp(){
	float p2 = getAonePM_freq() * getAonePM_freq();
	int pop = getPopPM();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoPM_exp(){
	float p2 = getAtwoPM_freq() * getAtwoPM_freq();
	int pop = getPopPM();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetPM_exp(){
	float pq = 2 * getAonePM_freq() * getAtwoPM_freq();
	int pop = getPopPM();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAonePM_freq(){
	float percent = 0.0f;
	int denom = a1_countPM + a2_countPM;
	if(denom > 0){
		percent = ((float)a1_countPM/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoPM_freq(){
	float percent = 0.0f;
	int denom = a1_countPM + a2_countPM;
	if(denom > 0){
		percent = ((float)a2_countPM/(float)denom);
	}
	return percent;
}

//Child
float AlleleFrequency::getAonehomoC_exp(){
	float p2 = getAoneC_freq() * getAoneC_freq();
	int pop = getPopC();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoC_exp(){
	float p2 = getAtwoC_freq() * getAtwoC_freq();
	int pop = getPopC();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetC_exp(){
	float pq = 2 * getAoneC_freq() * getAtwoC_freq();
	int pop = getPopC();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneC_freq(){
	float percent = 0.0f;
	int countone = a1_countC;
	int counttwo = a2_countC;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)countone/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoC_freq(){
	float percent = 0.0f;
	int countone = a1_countC;
	int counttwo = a2_countC;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)counttwo/(float)denom);
	}
	return percent;
}

//Child Male
float AlleleFrequency::getAonehomoCM_exp(){
	float p2 = getAoneCM_freq() * getAoneCM_freq();
	int pop = getPopCM();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoCM_exp(){
	float p2 = getAtwoCM_freq() * getAtwoCM_freq();
	int pop = getPopCM();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetCM_exp(){
	float pq = 2 * getAoneCM_freq() * getAtwoCM_freq();
	int pop = getPopCM();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneCM_freq(){
	float percent = 0.0f;
	int denom = a1_countCM + a2_countCM;
	if(denom > 0){
		percent = ((float)a1_countCM/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoCM_freq(){
	float percent = 0.0f;
	int denom = a1_countCM + a2_countCM;
	if(denom > 0){
		percent = ((float)a2_countCM/(float)denom);
	}
	return percent;
}

//Child Female
float AlleleFrequency::getAonehomoCF_exp(){
	float p2 = (float)pow(getAoneCF_freq(), 2);
	int pop = getPopCF();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoCF_exp(){
	float p2 = getAtwoCF_freq() * getAtwoCF_freq();
	int pop = getPopCF();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetCF_exp(){
	float pq = 2 * getAoneCF_freq() * getAtwoCF_freq();
	int pop = getPopCF();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneCF_freq(){
	float percent = 0.0f;
	int countone = a1_countCF;
	int counttwo = a2_countCF;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)countone/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoCF_freq(){
	float percent = 0.0f;
	int countone = a1_countCF;
	int counttwo = a2_countCF;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)counttwo/(float)denom);
	}
	return percent;
}

//Case
float AlleleFrequency::getAonehomoCa_exp(){
	float p2 = getAoneCa_freq() * getAoneCa_freq();
	int pop = getPopCa();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoCa_exp(){
	float p2 = getAtwoCa_freq() * getAtwoCa_freq();
	int pop = getPopCa();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetCa_exp(){
	float pq = 2 * getAoneCa_freq() * getAtwoCa_freq();
	int pop = getPopCa();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneCa_freq(){
	float percent = 0.0f;
	int denom = a1_countCa + a2_countCa;
	if(denom > 0){
		percent = ((float)a1_countCa/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoCa_freq(){
	float percent = 0.0f;
	int denom = a1_countCa + a2_countCa;
	if(denom > 0){
		percent = ((float)a2_countCa/(float)denom);
	}
	return percent;
}

//Case Female
float AlleleFrequency::getAonehomoCaF_exp(){
	float p2 = getAoneCaF_freq() * getAoneCaF_freq();
	int pop = getPopCaF();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoCaF_exp(){
	float p2 = getAtwoCaF_freq() * getAtwoCaF_freq();
	int pop = getPopCaF();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetCaF_exp(){
	float pq = 2 * getAoneCaF_freq() * getAtwoCaF_freq();
	int pop = getPopCaF();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneCaF_freq(){
	float percent = 0.0f;
	int countone = a1_countCaF;
	int counttwo = a2_countCaF;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)countone/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoCaF_freq(){
	float percent = 0.0f;
	int countone = a1_countCaF;
	int counttwo = a2_countCaF;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)counttwo/(float)denom);
	}
	return percent;
}

//Case Male
float AlleleFrequency::getAonehomoCaM_exp(){
	float p2 = getAoneCaM_freq() * getAoneCaM_freq();
	int pop = getPopCaM();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoCaM_exp(){
	float p2 = getAtwoCaM_freq() * getAtwoCaM_freq();
	int pop = getPopCaM();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetCaM_exp(){
	float pq = 2 * getAoneCaM_freq() * getAtwoCaM_freq();
	int pop = getPopCaM();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneCaM_freq(){
	float percent = 0.0f;
	int denom = a1_countCaM + a2_countCaM;
	if(denom > 0){
		percent = ((float)a1_countCaM/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoCaM_freq(){
	float percent = 0.0f;
	int denom = a1_countCaM + a2_countCaM;
	if(denom > 0){
		percent = ((float)a2_countCaM/(float)denom);
	}
	return percent;
}

//Control
float AlleleFrequency::getAonehomoCon_exp(){
	float p2 = getAoneCon_freq() * getAoneCon_freq();
	int pop = getPopCon();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoCon_exp(){
	float p2 = getAtwoCon_freq() * getAtwoCon_freq();
	int pop = getPopCon();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetCon_exp(){
	float pq = 2 * getAoneCon_freq() * getAtwoCon_freq();
	int pop = getPopCon();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneCon_freq(){
	float percent = 0.0f;
	int countone = a1_countCon;
	int counttwo = a2_countCon;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)countone/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoCon_freq(){
	float percent = 0.0f;
	int countone = a1_countCon;
	int counttwo = a2_countCon;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)counttwo/(float)denom);
	}
	return percent;
}

//Control Female
float AlleleFrequency::getAonehomoConF_exp(){
	float p2 = getAoneConF_freq() * getAoneConF_freq();
	int pop = getPopConF();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoConF_exp(){
	float p2 = getAtwoConF_freq() * getAtwoConF_freq();
	int pop = getPopConF();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetConF_exp(){
	float pq = 2 * getAoneConF_freq() * getAtwoConF_freq();
	int pop = getPopConF();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneConF_freq(){
	float percent = 0.0f;
	int countone = a1_countConF;
	int counttwo = a2_countConF;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)countone/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoConF_freq(){
	float percent = 0.0f;
	int countone = a1_countConF;
	int counttwo = a2_countConF;
	int denom = countone + counttwo;
	if(denom > 0){
		percent = ((float)counttwo/(float)denom);
	}
	return percent;
}

//Control Male
float AlleleFrequency::getAonehomoConM_exp(){
	float p2 = getAoneConM_freq() * getAoneConM_freq();
	int pop = getPopConM();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getAtwohomoConM_exp(){
	float p2 = getAtwoConM_freq() * getAtwoConM_freq();
	int pop = getPopConM();
	float exp = p2 * pop;
	return exp;
}
float AlleleFrequency::getHetConM_exp(){
	float pq = 2 * getAoneConM_freq() * getAtwoConM_freq();
	int pop = getPopConM();
	float exp = pq * pop;
	return exp;
}
float AlleleFrequency::getAoneConM_freq(){
	float percent = 0.0f;
	int denom = a1_countConM + a2_countConM;
	if(denom > 0){
		percent = ((float)a1_countConM/(float)denom);
	}
	return percent;
}
float AlleleFrequency::getAtwoConM_freq(){
	float percent = 0.0f;
	int denom = a1_countConM + a2_countConM;
	if(denom > 0){
		percent = ((float)a2_countConM/(float)denom);
	}
	return percent;
}
}
