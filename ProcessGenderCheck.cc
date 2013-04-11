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


#include "ProcessGenderCheck.h"

#include <iostream>

#include <Options.h>
#include <GenderCheck.h>
#include <Helpers.h>
#include <MethodException.h>

using std::string;
using std::ofstream;
using std::vector;
using std::map;
using std::getString;
using Methods::opts;
using Methods::GenderCheck;
using Methods::Helpers;
using Methods::MethodException;
using Methods::DataSet;

const string ProcessGenderCheck::stepname = ProcessGenderCheck::doRegister("gender-error");

/*
 * Function: PrintSummary
 * Description:
 * Outputs results of calculations
 */
void ProcessGenderCheck::PrintSummary(){
	int msize = good_markers.size();//data_set->num_loci();
	int ssize = data_set->num_inds();

	string fname1 = opts::_OUTPREFIX_ + "gender_check_marker" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	ofstream myoutput (fname1.c_str());
	if(!myoutput){
		opts::printLog("Error opening " + fname1 + "!  Exiting!\n");
		throw MethodException("");
	}
	opts::addFile("Marker", stepname, fname1);
	myoutput.precision(4);
	myoutput << "Chrom\trsID\tProbeID\tbploc\t";
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		myoutput << data_set->get_locus(0)->getDetailHeaders() << "\t";
	}
	myoutput << "Count_Geno_all_Males\tCount_Geno_HET_Males\t%HET_Males" << endl;
	opts::addHeader(fname1, "Count_Geno_all_Males");
	opts::addHeader(fname1, "Count_Geno_HET_Males");
	opts::addHeader(fname1, "%HET_Males");

	for(int i = 0; i < msize; i++){
		if(good_markers[i]->isEnabled() && (good_markers[i]->getChrom() == opts::_CHRX_ || options.doAll())){// && !data_set->get_locus(i)->isFlagged()){

			myoutput << good_markers[i]->toString() << "\t"
				<< mtotal[i] << "\t"
				<< merrors[i] << "\t";
			if(mtotal[i] > 0){
				myoutput << (((float)merrors[i]/(float)mtotal[i]));// * 100.0f);
			}
			else{
				myoutput << "0";
			}
			myoutput << endl;
		}
	}
	string fname2 = opts::_OUTPREFIX_ + "gender_check_ind" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	ofstream hetout (fname2.c_str());
	if(!hetout){
		opts::printLog("Error opening " + fname2 + "!  Exiting!!\n");
		throw MethodException("");
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
        int msize = data_set->num_loci();
        for(int i = 0; i < msize; i++){
            if(good_markers[i]->isEnabled()){// && !data_set->get_locus(i)->isFlagged()){
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
        sdetails = data_set->get_sample(0)->getDetailHeaders();
    	hetout << "\t" << sdetails;
	}

    hetout << endl;

	for(int i = 0; i < ssize; i++){
		if(data_set->get_sample(i)->isEnabled()){
			hetout << data_set->get_sample(i)->toString() << "\t"
				<< data_set->get_sample(i)->getFamily()->getCenter() << "\t";
			if(data_set->get_sample(i)->getSex()){
				hetout << "M\t";
			}
			else{
				hetout << "F\t";
			}
			if(data_set->get_sample(i)->getPheno() == 2){
				hetout << "Y\t";
			}
			else if(data_set->get_sample(i)->getPheno() == 1){
				hetout << "N\t";
			}
			else{
				hetout << "U\t";
			}

			hetout << data_set->get_sample(i)->getPlate() << "\t"
				<< data_set->get_sample(i)->getWell() << "\t"
				<< stotal[i] << "\t"
				<< shets[i] << "\t";
			if(stotal[i] > 0){
				hetout << (((float)shets[i]/(float)stotal[i]));// * 100.0f);
			}
			else{
				hetout << "0";
			}

			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < (int)enzymes.size(); e++){
					hetout << "\t" << senzyme_tot[i][enzymes[e]]
						<< "\t" << senzyme_hets[i][enzymes[e]];
					if(senzyme_tot[i][enzymes[e]] > 0){
						hetout << "\t" << ((1.0f - ((float)senzyme_hets[i][enzymes[e]]/(float)senzyme_tot[i][enzymes[e]])));// * 100.0f);
					}
					else{
						hetout << "\t" << "0";
					}
				}
			}
			if(opts::_SAMPDESC_.length() > 0){
				hetout << "\t" << data_set->get_sample(0)->getDetails();
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
		data_set->get_locus(i)->setFlag(false);
	}
}

/*
 * Function: filter
 * Description:
 * Performs filtering based on marker and sample error counts
 */
void ProcessGenderCheck::filter(){

	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int msize = good_markers.size();//data_set->num_loci();
		for(int m = 0; m < msize; m++){
			if(good_markers[m]->isEnabled()){// && !data_set->get_locus(m)->isFlagged()){

				bool inc = false;
				if(options.doThreshMarkersLow() && merrors[m] < options.getThreshMarkersLow()){
					good_markers[m]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && merrors[m] > options.getThreshMarkersHigh()){
					good_markers[m]->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}
		}
	}

	if(options.doThreshSamplesLow() || options.doThreshSamplesHigh()){
		int ssize = data_set->num_inds();
		for(int i = 0; i < ssize; i++){
			if(data_set->get_sample(i)->isEnabled() && data_set->get_sample(i)->getSex()){
				bool inc = false;
				if(options.doThreshSamplesLow() && shets[i] < options.getThreshSamplesLow()){
					data_set->get_sample(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshSamplesHigh() && shets[i] > options.getThreshSamplesHigh()){
					data_set->get_sample(i)->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_samples++;
				}
			}
		}
	}
}

/*
 * Function: process
 * Description:
 * Main entry into step.  Diverts to perform_evaluation
 */
void ProcessGenderCheck::process(DataSet* ds){
	data_set = ds;

	good_markers = Helpers::findValidMarkers(data_set->get_markers(), &options);

	GenderCheck gc(data_set);

	gc.setOptions(options);
	gc.calculate();

	merrors = gc.getMarkerErrors();
	shets = gc.getSampleHets();
	mtotal = gc.getMarkerTotalSamplesUsed();
	stotal = gc.getSampleTotalMarkersUsed();
	senzyme_hets = gc.getEnzymeHets();
	senzyme_tot = gc.getEnzymeTotals();
}

