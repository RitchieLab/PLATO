/**********************************************************************************
*                       Gender Check Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all Male genotypes on the X Chromosome and generates a count of
* all genotyping errors (not homozygous).  Outputs include counts by individual
* and counts by marker.
*
* Performs an initial scan.  Filters out any Markers with counts exceeding
* threshold.  Performs a secondary scan with bad markers removed.
*
*File: GenderCheck.cc
**********************************************************************************/


#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include "GenderCheck.h"
#include "General.h"
#include "Helpers.h"
namespace Methods{
string GenderCheck::stepname = "gender-error";

void GenderCheck::calcThreshold(){
}

void GenderCheck::setThreshold(string thresh){
	options.setUp(thresh);
}

/*
 * Function: PrintSummary
 * Description:
 * Outputs results of calculations
 */
void GenderCheck::PrintSummary(){
	int msize = markers->size();
	int ssize = samples->size();

	string fname1 = opts::_OUTPREFIX_ + "gender_check_marker" + options.getOut() + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	ofstream myoutput (fname1.c_str());
	if(!myoutput){
		opts::printLog("Error opening " + fname1 + "!  Exiting!\n");
		//exit(1);
		throw MethodException("Error opening " + fname1 + "!  Exiting!\n");
	}
	opts::addFile("Marker", stepname, fname1);
	myoutput.precision(4);
	myoutput << "Chrom\trsID\tProbeID\tbploc\t";
	if((*markers)[0]->getDetailHeaders().size() > 0){
		myoutput << (*markers)[0]->getDetailHeaders() << "\t";
	}
	myoutput << "Count_Geno_all_Males\tCount_Geno_HET_Males\t%HET_Males" << endl;
	opts::addHeader(fname1, "Count_Geno_all_Males");
	opts::addHeader(fname1, "Count_Geno_HET_Males");
	opts::addHeader(fname1, "%HET_Males");

	msize = good_markers.size();
	for(int i = 0; i < msize; i++){
		if(good_markers[i]->isEnabled() && (good_markers[i]->getChrom() == opts::_CHRX_ || options.doAll())){

			myoutput << good_markers[i]->toString() << "\t"
				<< mtotal[i] << "\t"
				<< merrors[i] << "\t";
			if(mtotal[i] > 0){
				myoutput << (((float)merrors[i]/(float)mtotal[i]) * 100.0f);
			}
			else{
				myoutput << "0";
			}
			myoutput << endl;
		}
	}
	string fname2 = opts::_OUTPREFIX_ + "gender_check_ind" + options.getOut() + ".txt";
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	ofstream hetout (fname2.c_str());
	if(!hetout){
		opts::printLog("Error opening " + fname2 + "!  Exiting!!\n");
		//exit(1);
		throw MethodException("Error opening " + fname2 + "!  Exiting!!\n");
	}
	opts::addFile("Sample", stepname, fname2);
	hetout.precision(4);
	hetout << "FamID\t";
	hetout << "IndID\tCenter\tSex\t";
	hetout << "Affection_Status\t";
	hetout << "Plate\tWell\tCount_Geno_total_all\tCount_Geno_HET_all\t%HET_all";
	opts::addHeader(fname2, "Count_Geno_total_all");
	opts::addHeader(fname2, "Count_Geno_HET_all");
	opts::addHeader(fname2, "%HET_all");

	vector<string> enzymes;
    if(opts::_ENZYMES_){
        int msize = good_markers.size();
        for(int i = 0; i < msize; i++){
            if(good_markers[i]->isEnabled()){
				vector<string>::iterator e_iter = find(enzymes.begin(), enzymes.end(), good_markers[i]->getEnzyme());
                if(e_iter == enzymes.end()){
                    enzymes.push_back(good_markers[i]->getEnzyme());
                }
            }
        }
        if(enzymes.size() > 1){
            sort(enzymes.begin(), enzymes.end());
        }
        for(int i = 0; i < (int)enzymes.size(); i++){
			hetout << "\tCount_Geno_total_" << enzymes[i] << "\tCount_Geno_HET_" << enzymes[i] << "\t%HET_" << enzymes[i];
        }
    }
    string sdetails = "";
    if(opts::_SAMPDESC_.length() > 0){
        sdetails = (*samples)[0]->getDetailHeaders();
    	hetout << "\t" << sdetails;
	}

    hetout << endl;

	for(int i = 0; i < ssize; i++){
		if((*samples)[i]->isEnabled()){
			hetout << (*samples)[i]->toString() << "\t"
				<< (*samples)[i]->getFamily()->getCenter() << "\t";
			if((*samples)[i]->getSex()){
				hetout << "M\t";
			}
			else{
				hetout << "F\t";
			}
			if((*samples)[i]->getPheno() == 2){
				hetout << "Y\t";
			}
			else if((*samples)[i]->getPheno() == 1){
				hetout << "N\t";
			}
			else{
				hetout << "U\t";
			}

			hetout << (*samples)[i]->getPlate() << "\t"
				<< (*samples)[i]->getWell() << "\t"
				<< stotal[i] << "\t"
				<< shets[i] << "\t";
			if(stotal[i] > 0){
				hetout << (((float)shets[i]/(float)stotal[i]) * 100.0f);
			}
			else{
				hetout << "0";
			}

			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < (int)enzymes.size(); e++){
					hetout << "\t" << senzyme_tot[i][enzymes[e]]
						<< "\t" << senzyme_hets[i][enzymes[e]];
					if(senzyme_tot[i][enzymes[e]] > 0){
						hetout << "\t" << ((1.0f - ((float)senzyme_hets[i][enzymes[e]]/(float)senzyme_tot[i][enzymes[e]])) * 100.0f);
					}
					else{
						hetout << "\t" << "0";
					}
				}
			}
			if(opts::_SAMPDESC_.length() > 0){
				hetout << "\t" << (*samples)[i]->getDetails();
			}
			hetout << endl;
		}
	}

	if(hetout.is_open()){
		hetout.close();
	}
	if(myoutput.is_open()){
		myoutput.close();
	}

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}
}

/*
 * Function: FilterSummary
 * Description:
 * Outputs remaining marker and sample counts
 */
void GenderCheck::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::printLog("Individuals Passed:\t" + getString<int>(opts::_SAMPLES_WORKING_ - orig_num_samples) + " (" +
	    getString<float>(((float) (opts::_SAMPLES_WORKING_ - orig_num_samples) / (float) opts::_SAMPLES_WORKING_) * 100.0) +
	    "%) of " + getString<int>(opts::_SAMPLES_WORKING_) + "\n");

	opts::_MARKERS_WORKING_ -= orig_num_markers;
	opts::_SAMPLES_WORKING_ -= orig_num_samples;

}

/*
 * Function: filter
 * Description:
 * Performs filtering based on marker and sample error counts
 */
void GenderCheck::filter(){

	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int msize = markers->size();
		for(int m = 0; m < msize; m++){
			if((*markers)[m]->isEnabled() && !(*markers)[m]->isFlagged()){

				bool inc = false;
				if(options.doThreshMarkersLow() && merrors[m] < options.getThreshMarkersLow()){
					(*markers)[m]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && merrors[m] > options.getThreshMarkersHigh()){
					(*markers)[m]->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}
		}
	}

	if(options.doThreshSamplesLow() || options.doThreshSamplesHigh()){
		int ssize = samples->size();
		for(int i = 0; i < ssize; i++){
			if((*samples)[i]->isEnabled() && (*samples)[i]->getSex()){
				bool inc = false;
				if(options.doThreshSamplesLow() && shets[i] < options.getThreshSamplesLow()){
					(*samples)[i]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshSamplesHigh() && shets[i] > options.getThreshSamplesHigh()){
					(*samples)[i]->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_samples++;
				}
			}
		}
	}
}

void GenderCheck::filter_markers(){
}

/*
 * Function: process
 * Description:
 * Main entry into step.  Diverts to perform_evaluation
 */
void GenderCheck::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

	calcThreshold();

	perform_evaluation(false);
}

/*
 * Function: perform_evaluation
 * Description:
 * Checks X chromosome for heterozygosity and homozygosity
 */
void GenderCheck::perform_evaluation(bool dofams){
	int msize = markers->size();
	int ssize = samples->size();

	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, options);
	msize = good_markers.size();

	merrors.resize(msize);
	shets.resize(ssize);
	mtotal.resize(msize);
	stotal.resize(ssize);
	senzyme_hets.resize(ssize);
	senzyme_tot.resize(ssize);

	for(int m = 0; m < msize; m++){
		if(good_markers[m]->isEnabled()){
			if(good_markers[m]->getChrom() != opts::_CHRX_ && !options.doAll()){
				continue;
			}

			int mloc = good_markers[m]->getLoc();
			bool micro = good_markers[m]->isMicroSat();
			for(int i = 0; i < ssize; i++){
				if((*samples)[i]->isEnabled()){
					if(!micro){
						if((*samples)[i]->getAone(mloc) && (*samples)[i]->getAtwo(mloc) && (*samples)[i]->getAmissing(mloc)){
						continue;
						}
						if(!(*samples)[i]->getAone(mloc) && (*samples)[i]->getAtwo(mloc)){
							shets[i]++;
							if(opts::_ENZYMES_){
								senzyme_hets[i][(*markers)[m]->getEnzyme()]++;
							}
							if((*samples)[i]->getSex()){
								merrors[m]++;
							}
						}
					}
					else{
						if((*samples)[i]->getAbone(mloc) == -1){
							continue;
						}
						if((*samples)[i]->getAbone(mloc) != (*samples)[i]->getAbtwo(mloc)){
							shets[i]++;
							if(opts::_ENZYMES_){
								senzyme_hets[i][(*markers)[m]->getEnzyme()]++;
							}
							if((*samples)[i]->getSex()){
								merrors[m]++;
							}
						}
					}
					if((*samples)[i]->getSex()){
						mtotal[m]++;
					}
					stotal[i]++;
					if(opts::_ENZYMES_){
						senzyme_tot[i][(*markers)[m]->getEnzyme()]++;
					}
				}
			}
		}
	}

}

}
