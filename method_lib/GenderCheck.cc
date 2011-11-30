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
	if((*markers).at(0)->getDetailHeaders().size() > 0){
		myoutput << (*markers).at(0)->getDetailHeaders() << "\t";
	}
	myoutput << "Count_Geno_all_Males\tCount_Geno_HET_Males\t%HET_Males" << endl;
	opts::addHeader(fname1, "Count_Geno_all_Males");
	opts::addHeader(fname1, "Count_Geno_HET_Males");
	opts::addHeader(fname1, "%HET_Males");

	msize = good_markers.size();
	for(int i = 0; i < msize; i++){
		if(good_markers.at(i)->isEnabled() && (good_markers.at(i)->getChrom() == opts::_CHRX_ || options.doAll())){

			myoutput << good_markers.at(i)->toString() << "\t"
				<< mtotal.at(i) << "\t"
				<< merrors.at(i) << "\t";
			if(mtotal.at(i) > 0){
				myoutput << (((float)merrors.at(i)/(float)mtotal.at(i)) * 100.0f);
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
            if(good_markers.at(i)->isEnabled()){
				vector<string>::iterator e_iter = find(enzymes.begin(), enzymes.end(), good_markers.at(i)->getEnzyme());
                if(e_iter == enzymes.end()){
                    enzymes.push_back(good_markers.at(i)->getEnzyme());
                }
            }
        }
        if(enzymes.size() > 1){
            sort(enzymes.begin(), enzymes.end());
        }
        for(int i = 0; i < (int)enzymes.size(); i++){
			hetout << "\tCount_Geno_total_" << enzymes.at(i) << "\tCount_Geno_HET_" << enzymes.at(i) << "\t%HET_" << enzymes.at(i);
        }
    }
    string sdetails = "";
    if(opts::_SAMPDESC_.length() > 0){
        sdetails = (*samples).at(0)->getDetailHeaders();
    	hetout << "\t" << sdetails;
	}

    hetout << endl;

	for(int i = 0; i < ssize; i++){
		if((*samples).at(i)->isEnabled()){
			hetout << (*samples).at(i)->toString() << "\t"
				<< (*samples).at(i)->getFamily()->getCenter() << "\t";
			if((*samples).at(i)->getSex()){
				hetout << "M\t";
			}
			else{
				hetout << "F\t";
			}
			if((*samples).at(i)->getPheno() == 2){
				hetout << "Y\t";
			}
			else if((*samples).at(i)->getPheno() == 1){
				hetout << "N\t";
			}
			else{
				hetout << "U\t";
			}

			hetout << (*samples).at(i)->getPlate() << "\t"
				<< (*samples).at(i)->getWell() << "\t"
				<< stotal.at(i) << "\t"
				<< shets.at(i) << "\t";
			if(stotal.at(i) > 0){
				hetout << (((float)shets.at(i)/(float)stotal.at(i)) * 100.0f);
			}
			else{
				hetout << "0";
			}

			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < (int)enzymes.size(); e++){
					hetout << "\t" << senzyme_tot.at(i)[enzymes.at(e)]
						<< "\t" << senzyme_hets.at(i)[enzymes.at(e)];
					if(senzyme_tot.at(i)[enzymes.at(e)] > 0){
						hetout << "\t" << ((1.0f - ((float)senzyme_hets.at(i)[enzymes.at(e)]/(float)senzyme_tot.at(i)[enzymes.at(e)])) * 100.0f);
					}
					else{
						hetout << "\t" << "0";
					}
				}
			}
			if(opts::_SAMPDESC_.length() > 0){
				hetout << "\t" << (*samples).at(i)->getDetails();
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
		(*markers).at(i)->setFlag(false);
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
			if((*markers).at(m)->isEnabled() && !(*markers).at(m)->isFlagged()){

				bool inc = false;
				if(options.doThreshMarkersLow() && merrors.at(m) < options.getThreshMarkersLow()){
					(*markers).at(m)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && merrors.at(m) > options.getThreshMarkersHigh()){
					(*markers).at(m)->setEnabled(false);
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
			if((*samples).at(i)->isEnabled() && (*samples).at(i)->getSex()){
				bool inc = false;
				if(options.doThreshSamplesLow() && shets.at(i) < options.getThreshSamplesLow()){
					(*samples).at(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshSamplesHigh() && shets.at(i) > options.getThreshSamplesHigh()){
					(*samples).at(i)->setEnabled(false);
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
		if(good_markers.at(m)->isEnabled()){
			if(good_markers.at(m)->getChrom() != opts::_CHRX_ && !options.doAll()){
				continue;
			}

			int mloc = good_markers.at(m)->getLoc();
			bool micro = good_markers.at(m)->isMicroSat();
			for(int i = 0; i < ssize; i++){
				if((*samples).at(i)->isEnabled()){
					if(!micro){
						if((*samples).at(i)->getAone(mloc) && (*samples).at(i)->getAtwo(mloc) && (*samples).at(i)->getAmissing(mloc)){
						continue;
						}
						if(!(*samples).at(i)->getAone(mloc) && (*samples).at(i)->getAtwo(mloc)){
							shets.at(i)++;
							if(opts::_ENZYMES_){
								senzyme_hets.at(i)[(*markers)[m]->getEnzyme()]++;
							}
							if((*samples).at(i)->getSex()){
								merrors.at(m)++;
							}
						}
					}
					else{
						if((*samples).at(i)->getAbone(mloc) == -1){
							continue;
						}
						if((*samples).at(i)->getAbone(mloc) != (*samples).at(i)->getAbtwo(mloc)){
							shets.at(i)++;
							if(opts::_ENZYMES_){
								senzyme_hets.at(i)[(*markers).at(m)->getEnzyme()]++;
							}
							if((*samples).at(i)->getSex()){
								merrors.at(m)++;
							}
						}
					}
					if((*samples).at(i)->getSex()){
						mtotal.at(m)++;
					}
					stotal.at(i)++;
					if(opts::_ENZYMES_){
						senzyme_tot.at(i)[(*markers).at(m)->getEnzyme()]++;
					}
				}
			}
		}
	}

}

}
