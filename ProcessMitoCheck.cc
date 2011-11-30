/**********************************************************************************
*                       Mito Check Module
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
* Files generated:
*	gender_errors_ind.txt
*	gender_errors_markers.txt
*       post_gender_error_filter_summary.txt
*       post_gender_error_filter_summary_chrom.txt
*
*File: MitoCheck.cc
**********************************************************************************/


#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include "ProcessMitoCheck.h"
#include "Chrom.h"
#include <General.h>
#include <Helper.h>
using namespace Methods;
string ProcessMitoCheck::stepname = "mito-check";

void ProcessMitoCheck::Tokenize(const string& str, vector<string>& tokens, const string& delimiter){
    string::size_type lastPos = str.find_first_not_of(delimiter, 0);
    string::size_type pos = str.find_first_of(delimiter, lastPos);

    while(string::npos != pos || string::npos != lastPos){
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiter, pos);
        pos = str.find_first_of(delimiter, lastPos);
    }
}

void ProcessMitoCheck::setThreshold(string thresh){
	options.setUp(thresh);
}

void ProcessMitoCheck::PrintSummary(){
	int msize = data_set->num_loci();
	int ssize = data_set->num_inds();

	string fname1 = opts::_OUTPREFIX_ + "mito_check_marker" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	ofstream myoutput (fname1.c_str());
	if(!myoutput){
		opts::printLog("Error opening " + fname1 + "!  Exiting!\n");
		throw MethodException("");
	}
	opts::addFile("Marker",stepname,fname1);
	myoutput.precision(4);
	myoutput << "Chrom\trsID\tProbeID\tbploc";
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		myoutput << "\t" << data_set->get_locus(0)->getDetailHeaders();
	}
	myoutput << "\tNum_errors" << endl;
	opts::addHeader(fname1, "Num_errors");

	for(int i = 0; i < msize; i++){
		if(data_set->get_locus(i)->isEnabled() && data_set->get_locus(i)->getChrom() == opts::_MITO_ && !data_set->get_locus(i)->isFlagged()){
			if(options.doChrom()){
				if(!options.checkChrom(data_set->get_locus(i)->getChrom())){
					continue;
				}
				if(!options.checkBp(data_set->get_locus(i)->getBPLOC())){
					continue;
				}
			}

			myoutput << data_set->get_locus(i)->toString() << "\t"
				<< merrors[i] << endl;
		}
	}
	string fname2 = opts::_OUTPREFIX_ + "mito_check_ind" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	ofstream hetout (fname2.c_str());
	if(!hetout){
		opts::printLog("Error opening " + fname2 + "!  Exiting!!\n");
		throw MethodException("");
	}
	opts::addFile("Sample",stepname,fname2);
	hetout.precision(4);
	hetout << "FamID\t";
	hetout << "IndID\tCenter\tSex\t";
	hetout << "Affection_Status\t";
	hetout << "Plate\tWell\tNum_errors";
	opts::addHeader(fname2, "Num_errors");

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
				<< serrors[i];
			if(opts::_SAMPDESC_.length() > 0){
				hetout << "\t" << data_set->get_sample(i)->getDetails();
			}
			hetout << endl;
		}
	}
	string fname3 = opts::_OUTPREFIX_ + "mito_check_errors_" + getString<int>(order) + ".txt";
	ofstream errorout (fname3.c_str());
	if(!errorout){
		opts::printLog("Error opening " + fname3 + "!  Exiting!!\n");
		exit(1);
	}
	errorout << "Chrom\trsID\tProbeID\tbploc\tFamID\tMom\tMom_Genotype\tChild\tChild_Genotype\n";

	for(int i = 0; i < (int)error_map.size(); i++){
		Sample* child = data_set->get_sample(i);
		Sample* mom = child->getMom();
		for(int j = 0; j < (int)error_map[i].size(); j++){
			Marker* mark = data_set->get_locus(error_map[i][j]);
			int loc = mark->getLoc();
			string ma1;
			string ma2;
			string ca1;
			string ca2;
			if(!mark->isMicroSat()){
				if(child->getAone(loc)){
					if(child->getAtwo(loc) && !child->getAmissing(loc)){
						ca1 = mark->getAllele2();
						ca2 = mark->getAllele2();
					}
					else{
						ca1 = "0";
						ca2 = "0";
					}
				}
				else{
					if(!child->getAtwo(loc)){
						ca1 = mark->getAllele1();
						ca2 = mark->getAllele1();
					}
					else{
						ca1 = mark->getAllele1();
						ca2 = mark->getAllele2();
					}
				}
				if(mom->getAone(loc)){
					if(mom->getAtwo(loc) && !mom->getAmissing(loc)){
						ma1 = mark->getAllele2();
						ma2 = mark->getAllele2();
					}
					else{
						ca1 = "0";
						ca2 = "0";
					}
				}
				else{
					if(!mom->getAtwo(loc)){
						ma1 = mark->getAllele1();
						ma2 = mark->getAllele1();
					}
					else{
						ma1 = mark->getAllele1();
						ma2 = mark->getAllele2();
					}
				}
			}
			else{
				ca1 = mark->getAllele(child->getAbone(loc));
				ca2 = mark->getAllele(child->getAbtwo(loc));
				ma1 = mark->getAllele(child->getAbone(loc));
				ma2 = mark->getAllele(child->getAbtwo(loc));
			}
			errorout << mark->toString() << "\t" << child->getFamID() << "\t" << mom->getInd() << "\t" << ma1 << "_" << ma2 << "\t" << child->getInd() << "\t" << ca1 << "_" << ca2 << endl;
		}

	}

	if(hetout.is_open()){
		hetout.close();
	}
	if(myoutput.is_open()){
		myoutput.close();
	}

	if(errorout.is_open()){
		errorout.close();
	}
	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessMitoCheck::FilterSummary(){
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

void ProcessMitoCheck::filter(){

	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int msize = data_set->num_loci();
		for(int m = 0; m < msize; m++){
			if(data_set->get_locus(m)->isEnabled() && data_set->get_locus(m)->getChrom() == 25 && !data_set->get_locus(m)->isFlagged()){
				if(options.doChrom()){
					if(!options.checkChrom(data_set->get_locus(m)->getChrom())){
					    continue;
				    }
				    if(!options.checkBp(data_set->get_locus(m)->getBPLOC())){
					    continue;
				    }
				}

				bool inc = false;
				if(options.doThreshMarkersLow() && merrors[m] < options.getThreshMarkersLow()){
					data_set->get_locus(m)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && merrors[m] > options.getThreshMarkersHigh()){
					data_set->get_locus(m)->setEnabled(false);
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
			if(data_set->get_sample(i)->isEnabled()){
				bool inc = false;
				if(options.doThreshSamplesLow() && serrors[i] < options.getThreshSamplesLow()){
					data_set->get_sample(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshSamplesHigh() && serrors[i] > options.getThreshSamplesHigh()){
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

void ProcessMitoCheck::filter_markers(){


/*	MKR* mymarkers = markers->getList();
	MKR::iterator m_iter;

ofstream myoutput ("gender_errors_markers.txt", ios::out | ios::app);
	for(m_iter = mymarkers->begin(); m_iter != mymarkers->end();){
		myoutput << m_iter->second.getChrom() << "\t" <<  m_iter->second.getRSID() << "\t" << m_iter->second.getBPLOC() << "\t" << m_iter->second.getMitoErrors() << endl;
		if(m_iter->second.getMitoErrors() > marker_thresh){
			mymarkers->erase(m_iter++);
		}
		else{
			++m_iter;
		}
	}
	myoutput.close();
*/
}

void ProcessMitoCheck::process(DataSet* ds){
	data_set = ds;
	MitoCheck mc;
	mc.setOptions(options);
	mc.calculate();
	merrors = mc.getNumMarkerErrors();
	serrors = mc.getNumSampleErrors();
	error_map = mc.getErrorMap();
}
