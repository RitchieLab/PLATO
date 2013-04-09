/**********************************************************************************
*                       Marker Genotype Efficiency Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all genotypes and generates a genotype efficiency for all markers.
*
*
* Files generated:
*	percent_breakdown_by_marker.txt
*	percent_breakdown_by_chrom.txt
*       post_marker_geno_eff_filter_summary.txt
*
*File: LD.cc
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
#include "ProcessMDRPDT.h"
#include <Options.h>
#include <General.h>
#include <Helpers.h>
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
using namespace Methods;

string ProcessMDRPDT::stepname = ProcessMDRPDT::doRegister("mdrpdt");

void ProcessMDRPDT::FilterSummary(){

	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessMDRPDT::PrintSummary(){
	int msize = data_set->num_loci();
	for(int m = 0; m < msize; m++){
		data_set->get_locus(m)->setFlag(false);
	}

}

void ProcessMDRPDT::filter(){
}

void ProcessMDRPDT::doFilter(Methods::Marker* mark, double value){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		if(mark->isEnabled() && !mark->isFlagged()){
			bool inc = false;
			if(options.doThreshMarkersLow() && Helpers::dLess(value, options.getThreshMarkersLow())){
				mark->setEnabled(false);
				inc = true;
			}
			if(options.doThreshMarkersHigh() && Helpers::dGreater(value, options.getThreshMarkersHigh())){
				mark->setEnabled(false);
				inc = true;
			}
			if(inc){
				orig_num_markers++;
			}
		}
	}

}


void ProcessMDRPDT::process(DataSet* ds){
	data_set = ds;
	vector<int> good_markers = Helpers::findValidMarkersIndexes(data_set->get_markers(), &options);

	//check if new covariate file is listed...or covariate name.
	//create vector of covariate indexes to use if specified.

    string fname = opts::_OUTPREFIX_ + "mdrpdt" + options.getOut() + ".txt";
    if(!overwrite){
        fname += "." + getString<int>(order);
    }
    ofstream mdrout (fname.c_str());
    if(!mdrout){
        opts::printLog("Error opening " + fname + "!  Exiting!\n");
        throw MethodException("");
    }
    mdrout.precision(4);


    mdrout << "Chrom\trsID\tProbeID\tBPLOC\tFold\tTstat\tMatched_OR\tPvalue\n";

    MDRPDT mdrpdt;
    mdrpdt.set_parameters(&options);
    mdrpdt.resetDataSet(data_set);

//   int prev_base = 0;
//    int prev_chrom = -1;
    int msize = good_markers.size();//data_set->num_loci();
    for(int m = 0; m < msize; m++){
    	Marker* mark = data_set->get_locus(good_markers[m]);
    	if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_chrom, prev_base)){
    		mdrpdt.calculate(good_markers[m]);

    		if(options.getMDRDPTNumCrossvals() != (int)mdrpdt.getMatchedOddsRatio().size()){
    			throw MethodException("Num crossvals != matched odds ratio size?");
    		}
    		double tstat = 0.0;
    		for(int f = 0; f < options.getMDRDPTNumCrossvals(); f++){
    			mdrout << mark->toString() << "\t"
					<< (f + 1) << "\t";
    			if(options.getMDRDPTNumCrossvals() == 1){
    				mdrout << mdrpdt.getTstatTraining()[f] << "\t";
    				tstat = mdrpdt.getTstatTraining()[f];
    			}
    			else{
    				mdrout << mdrpdt.getTstatTesting()[f] << "\t";
    				tstat = mdrpdt.getTstatTesting()[f];
    			}
    			mdrout << mdrpdt.getMatchedOddsRatio()[f] << "\t";
    			if(options.getMDRPDTNumPTests() > 0){
					mdrout << mdrpdt.getPvalue(mdrpdt.getAvgMatchedOddsRatio());
    			}
    			else{
    				mdrout << 0;
    			}
    			mdrout << endl;
    		}
    		doFilter(mark, tstat);

    	}
    }


    if(mdrout.is_open()){
    	mdrout.close();
    }
}

