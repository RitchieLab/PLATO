/**********************************************************************************
*                       Mendelian Error Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all family genotypes and calculates mendelian errors.  Generates
* counts by family and marker.  Also generates a list of deletion candidates.
* Errors are checked from Father to Child and Mother to child for a total of up
* to 3 mendelian errors per trio (or more if larger family structure).
* Performs initial scan and removes bad markers based on threshold.  A secondary
* scan is then performed.
*
*File: MendelianErrors.cc
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
#include "ProcessMendelianErrors.h"
#include "Chrom.h"
#include <General.h>
#include <Helpers.h>

using namespace std;
using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif

string ProcessMendelianErrors::stepname = "mendelian-error";

void ProcessMendelianErrors::setThreshold(string thresh){
	options.setUp(thresh);
}

void ProcessMendelianErrors::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::printLog("Families Passed:\t" + getString<int>(opts::_FAMILIES_WORKING_ - orig_num_families) + " (" +
        getString<float>(((float)(opts::_FAMILIES_WORKING_ - orig_num_families) / (float)opts::_FAMILIES_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_FAMILIES_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;
	opts::_FAMILIES_WORKING_ -= orig_num_families;

}

void ProcessMendelianErrors::PrintSummary(){
	string fname1 = opts::_OUTPREFIX_ + "mendelian_error_family" + options.getOut() + ".txt";//+ getString<int>(order) + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	string fname2 = opts::_OUTPREFIX_ + "mendelian_error_individual" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	string fname3 = opts::_OUTPREFIX_ + "mendelian_error_marker" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname3 += "." + getString<int>(order);
	}
	ofstream myoutputf (fname1.c_str());
	ofstream myoutputi (fname2.c_str());
	ofstream myoutputm (fname3.c_str());
	if(!myoutputf){
		opts::printLog("Error opening " + fname1 + ".  Exiting!\n");
		throw MethodException("");
	}
	opts::addFile("Family",stepname, fname1);
	if(!myoutputi){
		opts::printLog("Error opening " + fname2 + ".  Exiting!\n");
		throw MethodException("");
	}
	opts::addFile("Sample",stepname,fname2);
	if(!myoutputm){
		opts::printLog("Error opening " + fname3 + ".  Exiting!\n");
		throw MethodException("");
	}
	opts::addFile("Marker",stepname,fname3);
	myoutputi << "FamID\tIndID\tCenter\tSex\tAffection Satus\tPlate\tWell\tME_count_All";
	opts::addHeader(fname2, "ME_count_All");

	myoutputf << "FamID\tNumInds\tCenter\tME_count_All";
	opts::addHeader(fname1, "ME_count_All");

	myoutputm << "Chrom\trsID\tProbeID\tbploc\t";
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		myoutputm << data_set->get_locus(0)->getDetailHeaders() << "\t";
	}
	myoutputm << "ME_count_All" << endl;
	opts::addHeader(fname3, "ME_count_All");

	int msize = good_markers.size();//data_set->num_loci();
	int fsize = data_set->num_pedigrees();
	int ssize = data_set->num_inds();

	vector<string> enzymes;
    if(opts::_ENZYMES_){
        for(int i = 0; i < msize; i++){
            if(good_markers[i]->isEnabled()){//data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged()){
/*				if(options.doChrom()){
					if(!options.checkChrom(data_set->get_locus(i)->getChrom())){
					    continue;
				    }
				    if(!options.checkBp(data_set->get_locus(i)->getBPLOC())){
					    continue;
				    }
				}
*/
                vector<string>::iterator e_iter = find(enzymes.begin(), enzymes.end(), good_markers[i]->getEnzyme());//data_set->get_locus(i)->getEnzyme());
                if(e_iter == enzymes.end()){
                    enzymes.push_back(good_markers[i]->getEnzyme());//data_set->get_locus(i)->getEnzyme());
                }
            }
        }
        if(enzymes.size() > 1){
            sort(enzymes.begin(), enzymes.end());
        }
        for(int i = 0; i < (int)enzymes.size(); i++){
            myoutputi << "\tME_count_" << enzymes[i];
			myoutputf << "\tME_count_" << enzymes[i];
			opts::addHeader(fname1, "ME_count_" + enzymes[i]);
			opts::addHeader(fname2, "ME_count_" + enzymes[i]);
        }
    }
	string sdetails = "";
	if(opts::_SAMPDESC_.length() > 0){
		sdetails = data_set->get_sample(0)->getDetailHeaders();
	}
	myoutputi << "\t" << sdetails;
	myoutputi << endl;
	myoutputf << endl;

	for(int i = 0; i < msize; i++){
		if(good_markers[i]->isEnabled()){//data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged()){
/*			if(options.doChrom()){
				if(!options.checkChrom(data_set->get_locus(i)->getChrom())){
				    continue;
			    }
			    if(!options.checkBp(data_set->get_locus(i)->getBPLOC())){
				    continue;
			    }
			}
*/
			myoutputm << good_markers[i]->toString() << "\t"//data_set->get_locus(i)->toString() << "\t"
					  << merrors[i] << endl;
		}
	}

	for(int i = 0; i < fsize; i++){
		if(data_set->get_pedigree(i)->isEnabled()){
			myoutputf << data_set->get_pedigree(i)->getFamID() << "\t"
				<< data_set->get_pedigree(i)->getSamples()->size() << "\t"
				<< data_set->get_pedigree(i)->getCenter() << "\t"
				<< ferrors[i];
			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < (int)enzymes.size(); e++){
					myoutputf << "\t" << fenzyme[i][enzymes[e]];
				}
			}
			myoutputf << endl;
		}
	}

	for(int i = 0; i < ssize; i++){
		if(data_set->get_sample(i)->isEnabled()){
			myoutputi << data_set->get_sample(i)->getFamID() << "\t"
				<< data_set->get_sample(i)->getInd() << "\t"
				<< data_set->get_sample(i)->getFamily()->getCenter() << "\t";
			if(data_set->get_sample(i)->getSex()){
				myoutputi << "M\t";
			}
			else{
				myoutputi << "F\t";
			}
			if(data_set->get_sample(i)->getPheno() == 2){//getAffected()){
				myoutputi << "Y\t";
			}
			else if(data_set->get_sample(i)->getPheno() == 1){
				myoutputi << "N\t";
			}
			else{
				myoutputi << "U\t";
			}
			myoutputi << data_set->get_sample(i)->getPlate() << "\t"
				<< data_set->get_sample(i)->getWell() << "\t"
				<< serrors[i];
			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < (int)enzymes.size(); e++){
					myoutputi << "\t" << senzyme[i][enzymes[e]];
				}
			}
			if(opts::_SAMPDESC_.length() > 0){
				myoutputi << "\t" << data_set->get_sample(i)->getDetails();
			}
			myoutputi << endl;
		}
	}

	if(myoutputm.is_open()){
		myoutputm.close();
	}
	if(myoutputi.is_open()){
		myoutputi.close();
	}
	if(myoutputf.is_open()){
		myoutputf.close();
	}

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessMendelianErrors::process(DataSet* ds){
	data_set = ds;

	good_markers = Helpers::findValidMarkers(data_set->get_markers(), &options);

	MendelianErrors me(data_set);
	me.setOrder(this->order);
	me.setOverwrite(this->overwrite);
	me.setOptions(options);
	me.calculate();
	error_map = me.getErrorMap();
	merrors = me.getNumMarkerErrors();
	ferrors = me.getNumFamilyErrors();
	serrors = me.getNumSampleErrors();
	if(options.zeroGenos()){
		zeroErrors();
	}
	filter_markers();
	//resetCounts();
	//perform_evaluation(false);
}

void ProcessMendelianErrors::zeroErrors(){
	int esize = error_map.size();
	for(int s = 0; s < esize; s++){
		if(error_map[s].size() > 0){
			for(int m = 0; m < (int)error_map[s].size(); m++){
				Marker* aloc = error_map[s][m];
				int mloc = aloc->getLoc();//data_set->get_locus(aloc)->getLoc();
				data_set->get_sample(s)->addAone(mloc, true);
				data_set->get_sample(s)->addAtwo(mloc, true);
				data_set->get_sample(s)->addAmissing(mloc, true);
				if(aloc->isMicroSat()){//data_set->get_locus(aloc)->isMicroSat()){
					data_set->get_sample(s)->addAbone(mloc, -1);
					data_set->get_sample(s)->addAbtwo(mloc, -1);
				}
			}
		}
	}
}


void ProcessMendelianErrors::filter_markers(){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int msize = good_markers.size();//data_set->num_loci();

		for(int i = 0; i < msize; i++){
			if(good_markers[i]->isEnabled()){//data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged()){
/*				if(options.doChrom()){
					if(!options.checkChrom(data_set->get_locus(i)->getChrom())){
				    	continue;
			    	}
			    	if(!options.checkBp(data_set->get_locus(i)->getBPLOC())){
				    	continue;
			    	}
				}
*/
				bool inc = false;
				if(options.doThreshMarkersHigh() && merrors[i] > options.getThreshMarkersHigh()){
					good_markers[i]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersLow() && merrors[i] < options.getThreshMarkersLow()){
					good_markers[i]->setEnabled(false);
					inc = true;
				}

				if(inc){
					orig_num_markers++;
				}
			}
		}
	}
}

void ProcessMendelianErrors::filter(){
	if(options.doThreshFamiliesLow() || options.doThreshFamiliesHigh()){
		int fsize = data_set->num_pedigrees();

		for(int i = 0; i < fsize; i++){
			if(data_set->get_pedigree(i)->isEnabled()){
				bool inc = false;
				if(options.doThreshFamiliesLow() && ferrors[i] < options.getThreshFamiliesLow()){
					data_set->get_pedigree(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshFamiliesHigh() && ferrors[i] > options.getThreshFamiliesHigh()){
					data_set->get_pedigree(i)->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_families++;
				}
			}
		}
	}
}
#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
