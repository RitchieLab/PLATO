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


#include "ProcessMDRPDT.h"
#include <MDRPDT.h>

#include <iostream>
#include <vector>
#include <map>

#include <Helpers.h>
#include <Options.h>
#include <MethodException.h>

using std::string;
using std::ofstream;
using std::vector;

using Methods::DataSet;
using Methods::MDRPDT;
using Methods::opts;
using Methods::MethodException;
using Methods::Marker;
using Methods::Helpers;

const string ProcessMDRPDT::stepname = ProcessMDRPDT::doRegister("mdrpdt");

void ProcessMDRPDT::PrintSummary(){
	int msize = data_set->num_loci();
	for(int m = 0; m < msize; m++){
		data_set->get_locus(m)->setFlag(false);
	}

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

