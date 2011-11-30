/**********************************************************************************
*                       Deletion detection
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Performs deletion detection based on mendelian error as described in
* Conrad et al.
*
*
*File: Deletions.cc
**********************************************************************************/


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include "Globals.h"
#include <string>
#include <list>
#include <map>
#include <algorithm>
#include <bitset>
#include "Deletions.h"
#include "General.h"
#include "Helpers.h"

using namespace std;
namespace Methods{
string Deletions::stepname="deletions";

void Deletions::setThreshold(string thresh){
	options->setUp(thresh);
}

/*
 * Function: FilterSummary
 * Description:
 * Outputs remaining markers and families
 */
void Deletions::FilterSummary(){
	opts::printLog("Options:\t" + options->toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::printLog("Families Passed:\t" + getString<int>(opts::_FAMILIES_WORKING_ - orig_num_families) + " (" +
        getString<float>(((float)(opts::_FAMILIES_WORKING_ - orig_num_families) / (float)opts::_FAMILIES_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_FAMILIES_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;
	opts::_FAMILIES_WORKING_ -= orig_num_families;

}

void Deletions::calcThreshold(){
}

/*
 * Function: PrintSummary
 * Description:
 * Resets marker flags
 */

/*
 * Function: process
 * Description:
 * Main lead into processing, diverts to perform_evaliation
 */
void Deletions::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

	calcThreshold();
cout << "HERE!!!\n";
	perform_evaluation(false);
}

void Deletions::zeroErrors(){

	int esize = error_map.size();
	for(int s = 0; s < esize; s++){
		if(error_map.at(s).size() > 0){
			for(int m = 0; m < (int)error_map.at(s).size(); m++){
				int aloc = error_map.at(s).at(m);
				int mloc = (*markers).at(aloc)->getLoc();
				(*samples).at(s)->addAone(mloc, true);
				(*samples).at(s)->addAtwo(mloc, false);
				if((*markers).at(aloc)->isMicroSat()){
					(*samples).at(s)->addAbone(mloc, -1);
					(*samples).at(s)->addAbtwo(mloc, -1);
				}
			}
		}
	}
}

vector<int> Deletions::thisfindValidMarkersIndexes(vector<Methods::Marker*>* marks, StepOptions* options){
	int prev_base = 0;
	int prev_chrom = -1;
	vector<int> good_markers;
	int msize = marks->size();
	for(int m = 0; m < msize; m++){
		Methods::Marker* mark = (*marks).at(m);
		if(mark->isEnabled()){
			if(options->doChrom()){
				if(!options->checkChrom(mark->getChrom())){
					continue;
				}
				if(!options->checkBp(mark->getBPLOC())){
					continue;
				}
			}
			if(options->doBpSpace()){
				if(prev_base == 0){
					prev_base = mark->getBPLOC();
					prev_chrom = mark->getChrom();
				}
				else{
					if(mark->getChrom() == prev_chrom && ((mark->getBPLOC() - prev_base) < options->getBpSpace())){
						continue;
					}
					prev_base = mark->getBPLOC();
					prev_chrom = mark->getChrom();
				}
			}
			good_markers.push_back(m);
		}
	}

	return good_markers;
}

/*
 * Function: perform_evaluation
 * Description:
 * Does the main processing of deletion detection
 */
void Deletions::perform_evaluation(bool output){
	ofstream myoutput;
	ofstream myoutputd;
	string filename = opts::_OUTPREFIX_ + "deletion" + options->getOut() + ".txt";
	if(!overwrite){
		filename += "." + getString<int>(order);
	}

	default_filename = filename;
	myoutput.open(filename.c_str(), ios::out);
	if(!myoutput){
		opts::printLog("Error opening " + filename + "\n");
		throw MethodException("Error opening " + filename + "\n");
	}
	myoutput << "FamID\tIndID\tChrom\tStart_rsID\tEnd_rsID\tStart_bploc\tEnd_bploc\tBploc_range\tNum_snps\tType\n";
	string filename2 = opts::_OUTPREFIX_ + "deletion_density" + options->getOut() + ".txt";
	if(!overwrite){
		filename2 += "." + getString<int>(order);
	}
	density_filename = filename2;
	myoutputd.open(filename2.c_str(), ios::out);
	opts::addFile("Marker", stepname, filename2);
	if(!myoutputd){
		opts::printLog("Error opening " + filename2 + "\n");
		throw MethodException("Error opening " + filename2 + "\n");
	}
	myoutputd << "Chrom\trsID\tProbeID\tbploc\t";
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		myoutputd << (*markers).at(0)->getDetailHeaders() << "\t";
	}
	myoutputd << "Ind_Count" << endl;
	opts::addHeader(filename2, "Ind_Count");

	int msize = data_set->num_loci();
	int ssize = data_set->num_inds();
	vector<int> marker_counts;
	marker_counts.resize(0,0);
	good_markers.resize(0,0);

	marker_counts.resize(msize, 0);

	vector<bool> marker_flags;
	marker_flags.resize(msize, false);

	for(int s = 0; s < ssize; s++){
		Sample* samp = data_set->get_sample(s);
		if(samp->isEnabled() && samp->getDad() != NULL && samp->getMom() != NULL){
			Sample* mom = samp->getMom();
			Sample* dad = samp->getDad();
			Sample* child = samp;
			vector<char> cats(msize);

			int prev_base = 0;
			int prev_chrom = -1;

			for(int m = 0; m < msize; m++){
					Methods::Marker* mark = (*markers).at(m);
					if(mark->isEnabled()){
						if(options->doChrom()){
							if(!options->checkChrom(mark->getChrom())){
								marker_flags.at(m) = true;
								continue;
							}
							if(!options->checkBp(mark->getBPLOC())){
								marker_flags.at(m) = true;
								continue;
							}
						}
						if(options->doBpSpace()){
							if(prev_base == 0){
								prev_base = mark->getBPLOC();
								prev_chrom = mark->getChrom();
							}
							else{
								if(mark->getChrom() == prev_chrom && ((mark->getBPLOC() - prev_base) < options->getBpSpace())){
									marker_flags.at(m) = true;
									continue;
								}
								prev_base = mark->getBPLOC();
								prev_chrom = mark->getChrom();
							}
						}
					}

				if(data_set->get_locus(m)->isEnabled()){

					int loc = data_set->get_locus(m)->getLoc();

					if(mom && dad && mom->isEnabled() && dad->isEnabled()){
						Marker* mark = data_set->get_locus(m);
						if(!mark->isMicroSat()){
							if(dad->getAone(loc) == dad->getAtwo(loc) && !dad->getAmissing(loc) && dad->getAone(loc) == child->getAone(loc) && child->getAone(loc) == child->getAtwo(loc) && !child->getAmissing(loc) && mom->getAone(loc) == mom->getAtwo(loc) && !mom->getAmissing(loc) && mom->getAone(loc) != child->getAone(loc)){
								cats.at(loc) =  'A';
							}
							else if(dad->getAone(loc) && dad->getAtwo(loc) && dad->getAmissing(loc) && child->getAone(loc) == child->getAtwo(loc) && !child->getAmissing(loc) && mom->getAone(loc) == mom->getAtwo(loc) && !mom->getAmissing(loc) && mom->getAone(loc) != child->getAone(loc)){
								cats.at(loc) = 'A';
							}
							else if(mom->getAone(loc) == mom->getAtwo(loc) && !mom->getAmissing(loc) && mom->getAone(loc) == child->getAone(loc) && child->getAone(loc) == child->getAtwo(loc) && !child->getAmissing(loc) && dad->getAone(loc) == dad->getAtwo(loc) && !dad->getAmissing(loc) && dad->getAone(loc) != child->getAone(loc)){
								cats.at(loc) = 'B';
							}
							else if(mom->getAone(loc) && mom->getAtwo(loc) && mom->getAmissing(loc) && child->getAone(loc) == child->getAtwo(loc) && !child->getAmissing(loc) && dad->getAone(loc) == dad->getAtwo(loc) && !dad->getAone(loc) && dad->getAone(loc) != child->getAtwo(loc)){
								cats.at(loc)='B';
							}
							else if(dad->getAone(loc) == dad->getAtwo(loc) && !dad->getAmissing(loc) && mom->getAone(loc) == mom->getAtwo(loc) && !mom->getAmissing(loc) && dad->getAone(loc) == mom->getAone(loc) && ((!child->getAone(loc) && child->getAtwo(loc)) || (child->getAone(loc) == child->getAtwo(loc) && !child->getAmissing(loc) && dad->getAone(loc) != child->getAone(loc)))){
								cats.at(loc)= 'C';
							}
							else if(dad->getAone(loc) == dad->getAtwo(loc) && !dad->getAmissing(loc) && mom->getAone(loc) && mom->getAtwo(loc) && mom->getAmissing(loc) && dad->getAone(loc) != child->getAone(loc) && child->getAone(loc) == child->getAtwo(loc) && !child->getAmissing(loc)){
								cats.at(loc) ='C';
							}
							else if(mom->getAone(loc) == mom->getAtwo(loc) && !mom->getAmissing(loc) && dad->getAone(loc) && dad->getAtwo(loc) && dad->getAmissing(loc) && mom->getAone(loc) != child->getAone(loc) && child->getAone(loc) == child->getAtwo(loc) && !child->getAmissing(loc)){
								cats.at(loc)= 'C';
							}
							else if(((child->getAone(loc) == child->getAtwo(loc) && !child->getAmissing(loc)) || (child->getAone(loc) && child->getAtwo(loc) && child->getAmissing(loc))) && ((dad->getAone(loc) == dad->getAtwo(loc) && !dad->getAmissing(loc)) || (dad->getAone(loc) && dad->getAtwo(loc) && dad->getAmissing(loc))) && ((mom->getAone(loc) == mom->getAtwo(loc) && !mom->getAmissing(loc)) || (mom->getAone(loc) && mom->getAtwo(loc) && mom->getAmissing(loc)))){
								cats.at(loc)= 'D';
							}
							else if(((child->getAone(loc) == child->getAtwo(loc) && !child->getAmissing(loc)) || (child->getAone(loc) && child->getAtwo(loc) && child->getAmissing(loc))) && !dad->getAone(loc) && dad->getAtwo(loc) && ((mom->getAone(loc) == mom->getAtwo(loc) && !mom->getAmissing(loc)) || (mom->getAone(loc) && mom->getAtwo(loc) && mom->getAmissing(loc)))){
								cats.at(loc)= 'E';
							}
							else if(((child->getAone(loc) == child->getAtwo(loc) && !child->getAmissing(loc)) || (child->getAone(loc) && child->getAtwo(loc) && child->getAmissing(loc))) && !mom->getAone(loc) && mom->getAtwo(loc) && ((dad->getAone(loc) == dad->getAtwo(loc) && !dad->getAmissing(loc)) || (dad->getAone(loc) && dad->getAtwo(loc) && dad->getAmissing(loc)))){
								cats.at(loc)= 'F';
							}
							else if((!child->getAone(loc) && child->getAtwo(loc)) || (!mom->getAone(loc) && mom->getAtwo(loc) && !dad->getAone(loc) && dad->getAtwo(loc))){
								cats.at(loc)= 'G';
							}
						}
						else{
							if(dad->getAbone(loc) == dad->getAbtwo(loc) && dad->getAbone(loc) == child->getAbone(loc) && child->getAbone(loc) == child->getAbtwo(loc) && mom->getAbone(loc) == mom->getAbtwo(loc) && mom->getAbone(loc) != child->getAbone(loc)){
								cats.at(loc)= 'A';
							}
							else if(dad->getAbone(loc) == -1 && dad->getAbtwo(loc) == -1 && child->getAbone(loc) == child->getAbtwo(loc) && mom->getAbone(loc) == mom->getAbtwo(loc) && mom->getAbone(loc) != child->getAbone(loc)){
								cats.at(loc)= 'A';
							}
							else if(mom->getAbone(loc) == mom->getAbtwo(loc) && mom->getAbone(loc) == child->getAbone(loc) && child->getAbone(loc) == child->getAbtwo(loc) && dad->getAbone(loc) == dad->getAbtwo(loc) && dad->getAbone(loc) != child->getAbone(loc)){
								cats.at(loc)= 'B';
							}
							else if(mom->getAbone(loc) == -1 && mom->getAbtwo(loc) == -1 && child->getAbone(loc) == child->getAbtwo(loc) && dad->getAbone(loc) == dad->getAbtwo(loc) && dad->getAbone(loc) != child->getAbtwo(loc)){
								cats.at(loc)= 'B';
							}
							else if(dad->getAbone(loc) == dad->getAbtwo(loc) && mom->getAbone(loc) == mom->getAbtwo(loc) && dad->getAbone(loc) == mom->getAbone(loc) && ((child->getAbone(loc) != child->getAbtwo(loc)) || (child->getAbone(loc) == child->getAbtwo(loc) && dad->getAbone(loc) != child->getAbone(loc)))){
								cats.at(loc)= 'C';
							}
							else if(dad->getAbone(loc) == dad->getAbtwo(loc) && mom->getAbone(loc) == -1 && mom->getAbtwo(loc) == -1 && dad->getAbone(loc) != child->getAbone(loc) && child->getAbone(loc) == child->getAbtwo(loc)){
								cats.at(loc)= 'C';
							}
							else if(mom->getAbone(loc) == mom->getAbtwo(loc) && dad->getAbone(loc) == -1 && dad->getAbtwo(loc) == -1 && mom->getAbone(loc) != child->getAbone(loc) && child->getAbone(loc) == child->getAbtwo(loc)){
								cats.at(loc)= 'C';
							}
							else if(((child->getAbone(loc) == child->getAbtwo(loc)) || (child->getAbone(loc) == -1 && child->getAbtwo(loc) == -1)) && ((dad->getAbone(loc) == dad->getAbtwo(loc)) || (dad->getAbone(loc) == -1 && dad->getAbtwo(loc) == -1)) && ((mom->getAbone(loc) == mom->getAbtwo(loc)) || (mom->getAbone(loc) == -1 && mom->getAbtwo(loc) == -1))){
								cats.at(loc)= 'D';
							}
							else if(((child->getAbone(loc) == child->getAbtwo(loc)) || (child->getAbone(loc) == -1 && child->getAbtwo(loc) == -1)) && dad->getAbone(loc) != dad->getAbtwo(loc) && ((mom->getAbone(loc) == mom->getAbtwo(loc)) || (mom->getAbone(loc) == -1 && mom->getAbtwo(loc) == -1))){
								cats.at(loc)= 'E';
							}
							else if(((child->getAbone(loc) == child->getAbtwo(loc)) || (child->getAbone(loc) == -1 && child->getAbtwo(loc) == -1)) && mom->getAbone(loc) != mom->getAbtwo(loc) && ((dad->getAbone(loc) == dad->getAbtwo(loc)) || (dad->getAbone(loc) == -1 && dad->getAbtwo(loc) == -1))){
								cats.at(loc)= 'F';
							}
							else if((child->getAbone(loc) != child->getAbtwo(loc)) || (mom->getAbone(loc) != mom->getAbtwo(loc) && dad->getAbone(loc) != dad->getAbtwo(loc))){
								cats.at(loc)= 'G';
							}
						}
					}
				}//end else
			}//end foreach marker

			prev_chrom = 0;
			prev_base = 0;
			Marker* markstartmat = NULL;
			Marker* markstartpat = NULL;
			Marker* markendmat = NULL;
			Marker* markendpat = NULL;
			int markstartmat_loc = -1;
			int markstartpat_loc = -1;
			int markendmat_loc = -1;
			int markendpat_loc = -1;
			int matspan = 0;
			int patspan = 0;
			int currchrom = -1;

			for(int m = 0; m < msize; m++){
				Marker* mark = data_set->get_locus(m);
				if(mark->isEnabled() && !marker_flags.at(m)){
					int loc = mark->getLoc();
					if(currchrom == -1){
						currchrom = mark->getChrom();
					}

					if(currchrom != mark->getChrom()){
						//output span
						if(matspan >= options->getDeletionSpan()){
							for(int l = markstartmat_loc; l <= markendmat_loc; l++){
								Marker* temp = data_set->get_locus(l);
								if(temp->isEnabled() && !marker_flags.at(l)){
									marker_counts.at(l)++;
								}
							}
							myoutput << samp->getFamID() << "\t" << samp->getInd() << "\t" << currchrom << "\t" << markstartmat->getProbeID() << "\t" << markendmat->getProbeID() << "\t" << markstartmat->getBPLOC() << "\t" << markendmat->getBPLOC() << "\t" << (markendmat->getBPLOC() - markstartmat->getBPLOC()) << "\t" << matspan << "\tM" << endl;
						}
						if(patspan >= options->getDeletionSpan()){
							for(int l = markstartpat_loc; l <= markendpat_loc; l++){
								Marker* temp = data_set->get_locus(l);
								if(temp->isEnabled() && !marker_flags.at(l)){
									marker_counts.at(l)++;
								}
							}
							myoutput << samp->getFamID() << "\t" << samp->getInd() << "\t" << currchrom << "\t" << markstartpat->getProbeID() << "\t" << markendpat->getProbeID() << "\t" << markstartpat->getBPLOC() << "\t" << markendpat->getBPLOC() << "\t" << (markendpat->getBPLOC() - markstartpat->getBPLOC()) << "\t" << patspan << "\tP" << endl;
						}
						markstartmat = NULL;
						markstartpat = NULL;
						markendmat = NULL;
						markendpat = NULL;
						markstartmat_loc = -1;
						markstartpat_loc = -1;
						markendmat_loc = -1;
						markendpat_loc = -1;
						matspan = 0;
						patspan = 0;
						currchrom = mark->getChrom();
					}

					if(cats.at(loc) == 'A' || cats.at(loc) == 'D' || (cats.at(loc) == 'E' && markstartpat == NULL)){
						if(markstartmat == NULL){
							markstartmat = mark;
							markstartmat_loc = m;
						}
						markendmat = mark;
						markendmat_loc = m;
						matspan++;

						if(patspan >= options->getDeletionSpan()){
							for(int l = markstartpat_loc; l <= markendpat_loc; l++){
								Marker* temp = data_set->get_locus(l);
								if(temp->isEnabled() && !marker_flags.at(l)){
									marker_counts.at(l)++;
								}
							}
							myoutput << samp->getFamID() << "\t" << samp->getInd() << "\t" << currchrom << "\t" << markstartpat->getProbeID() << "\t" << markendpat->getProbeID() << "\t" << markstartpat->getBPLOC() << "\t" << markendpat->getBPLOC() << "\t" << (markendpat->getBPLOC() - markstartpat->getBPLOC()) << "\t" << patspan << "\tP" << endl;
						}
						patspan = 0;
						markstartpat = NULL;
						markendpat = NULL;
						markstartpat_loc = -1;
						markendpat_loc = -1;
					}
					else if(cats.at(loc) == 'B' || cats.at(loc) == 'D' || (cats.at(loc) == 'F' && markstartmat == NULL)){
						if(markstartpat == NULL){
							markstartpat = mark;
							markstartpat_loc = m;
						}
						markendpat = mark;
						markendpat_loc = m;
						patspan++;

						if(matspan >= options->getDeletionSpan()){
							for(int l = markstartmat_loc; l <= markendmat_loc; l++){
								Marker* temp = data_set->get_locus(l);
								if(temp->isEnabled() && !marker_flags.at(l)){
									marker_counts.at(l)++;
								}
							}
							myoutput << samp->getFamID() << "\t" << samp->getInd() << "\t" << currchrom << "\t" << markstartmat->getProbeID() << "\t" << markendmat->getProbeID() << "\t" << markstartmat->getBPLOC() << "\t" << markendmat->getBPLOC() << "\t" << (markendmat->getBPLOC() - markstartmat->getBPLOC()) << "\t" << matspan << "\tM" << endl;
						}
						matspan = 0;
						markstartmat = NULL;
						markendmat = NULL;
						markstartmat_loc = -1;
						markendmat_loc = -1;
					}
					else{
						if(matspan >= options->getDeletionSpan()){
							for(int l = markstartmat_loc; l <= markendmat_loc; l++){
								Marker* temp = data_set->get_locus(l);
								if(temp->isEnabled() && !marker_flags.at(l)){
									marker_counts.at(l)++;
								}
							}
							myoutput << samp->getFamID() << "\t" << samp->getInd() << "\t" << currchrom << "\t" << markstartmat->getProbeID() << "\t" << markendmat->getProbeID() << "\t" << markstartmat->getBPLOC() << "\t" << markendmat->getBPLOC() << "\t" << (markendmat->getBPLOC() - markstartmat->getBPLOC()) << "\t" << matspan << "\tM" << endl;
						}
						if(patspan >= options->getDeletionSpan()){
							for(int l = markstartpat_loc; l <= markendpat_loc; l++){
								Marker* temp = data_set->get_locus(l);
								if(temp->isEnabled() && !marker_flags.at(l)){
									marker_counts.at(l)++;
								}
							}
							myoutput << samp->getFamID() << "\t" << samp->getInd() << "\t" << currchrom << "\t" << markstartpat->getProbeID() << "\t" << markendpat->getProbeID() << "\t" << markstartpat->getBPLOC() << "\t" << markendpat->getBPLOC() << "\t" << (markendpat->getBPLOC() - markstartpat->getBPLOC()) << "\t" << patspan << "\tP" << endl;
						}
						markstartmat = NULL;
						markendmat = NULL;
						markstartpat = NULL;
						markendpat = NULL;
						markstartmat_loc = -1;
						markstartpat_loc = -1;
						markendmat_loc = -1;
						markendpat_loc = -1;
						matspan = 0;
						patspan = 0;
					}

				}
			}

			if(patspan >= options->getDeletionSpan()){
				for(int l = markstartpat_loc; l <= markendpat_loc; l++){
					Marker* temp = data_set->get_locus(l);
					if(temp->isEnabled() && !marker_flags.at(l)){
						marker_counts.at(l)++;
					}
				}
				myoutput << samp->getFamID() << "\t" << samp->getInd() << "\t" << currchrom << "\t" << markstartpat->getProbeID() << "\t" << markendpat->getProbeID() << "\t" << markstartpat->getBPLOC() << "\t" << markendpat->getBPLOC() << "\t" << (markendpat->getBPLOC() - markstartpat->getBPLOC()) << "\t" << patspan << "\tP" << endl;
			}
			if(matspan >= options->getDeletionSpan()){
				for(int l = markstartmat_loc; l <= markendmat_loc; l++){
					Marker* temp = data_set->get_locus(l);
					if(temp->isEnabled() && !marker_flags.at(l)){
						marker_counts[l]++;
					}
				}
				myoutput << samp->getFamID() << "\t" << samp->getInd() << "\t" << currchrom << "\t" << markstartmat->getProbeID() << "\t" << markendmat->getProbeID() << "\t" << markstartmat->getBPLOC() << "\t" << markendmat->getBPLOC() << "\t" << (markendmat->getBPLOC() - markstartmat->getBPLOC()) << "\t" << matspan << "\tM" << endl;
			}
		}
	}
	for(int i = 0; i < (int)marker_counts.size(); i++){
		Marker* mark = data_set->get_locus(i);
		if(mark != NULL){
		if(mark->isEnabled() && !marker_flags.at(i)){
			myoutputd << mark->toString() << "\t" << marker_counts.at(i) << endl;
		}
		}
	}

	if(myoutput && myoutput.is_open()){
		myoutput.close();
	}
	if(myoutputd && myoutputd.is_open()){
		myoutputd.close();
	}
	good_markers.clear();
}

void Deletions::filter_markers(){
	return;
}

void Deletions::filter(){
	return;
}
}
