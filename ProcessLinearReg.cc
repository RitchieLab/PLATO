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
#include <LinearRegression.h>
#include <MultComparison.h>
#include "ProcessLinearReg.h"
#include <Options.h>
#include <General.h>
#include <Helper.h>
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
using namespace Methods;

void ProcessLinearReg::FilterSummary(){

	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessLinearReg::PrintSummary(){
	int msize = data_set->num_loci();
	for(int m = 0; m < msize; m++){
		data_set->get_locus(m)->setFlag(false);
	}

}

void ProcessLinearReg::filter(){
}

void ProcessLinearReg::doFilter(Methods::Marker* mark, double value){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		if(mark->isEnabled() && !mark->isFlagged()){
			bool inc = false;
			if(options.doThreshMarkersLow() && dLess(value, options.getThreshMarkersLow())){
				mark->setEnabled(false);
				inc = true;
			}
			if(options.doThreshMarkersHigh() && dGreater(value, options.getThreshMarkersHigh())){
				mark->setEnabled(false);
				inc = true;
			}
			if(inc){
				orig_num_markers++;
			}
		}
	}
}


void ProcessLinearReg::process(DataSet* ds){
	data_set = ds;

	//check if new covariate file is listed...or covariate name.
	//create vector of covariate indexes to use if specified.

    string fname = opts::_OUTPREFIX_ + "linearreg" + options.getOut() + ".txt";
    if(!overwrite){
        fname += "." + getString<int>(order);
    }
    ofstream lrout (fname.c_str());
    if(!lrout){
        opts::printLog("Error opening " + fname + "!  Exiting!\n");
        throw MethodException("");
    }
    lrout.precision(4);

    lrout << "Chrom\trsID\tProbeID\tBPLOC\tReference_Allele\tTest\tNMISS\tBeta\tSTAT\tPvalue\n";

	LinearRegression lr;
	lr.setOptions(options);
	lr.resetDataSet(ds);

	vector<double> chis(ds->num_loci(), 0);
	vector<double> main_pvals(ds->num_loci(), 0);

	int prev_base = 0;
	int prev_chrom = -1;

	for(int i = 0; i < (int)ds->num_loci(); i++){
		Marker* mark = ds->get_locus(i);
		if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
			lr.calculate(i);
			vector<double>pvals = lr.getPvalues();
			vector<double>coefs = lr.getCoefs();
			vector<string>labels = lr.getLabels();
			vector<double> zs = lr.getZs();
			for(int l = 1; l < (int)labels.size(); l++){
				lrout << mark->getChrom() << "\t" << mark->getRSID() << "\t" << mark->getProbeID() << "\t";
				lrout << mark->getBPLOC() << "\t" << mark->getReferent() << "\t";
				lrout << labels[l] << "\t" << lr.getCalcMissing() << "\t" << coefs[l] << "\t" << zs[l] << "\t" << pvals[l] << endl;
			}

			chis[i] = lr.getStatistic();
			main_pvals[i] = pvals[1];
		}
//		lrout << options.getLinRModelType() << "\t";
//		lrout << lr.getCalcMissing() << "\t" << lr.getTestCoef() << "\t" << lr.getZ() << "\t" << lr.getPValue() << endl;
	}

	if(options.doMultCompare()){
		string fcomp = opts::_OUTPREFIX_ + "linearreg_comparisons" + options.getOut() + ".txt";
		if (!overwrite) {
			fcomp += "." + getString<int>(order);
		}
		ofstream COMP;

		COMP.open(fcomp.c_str(), ios::out);
		if(!COMP){
			throw MethodException("Could not open " + fcomp + " for output.\n");
		}

		   COMP << "Chrom"
				  << "\trsID"
				  << "\tProbeID"
				  << "\tbploc";
			if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
				COMP << "\t" << data_set->get_locus(0)->getDetailHeaders();
			}

			COMP  << "\tCALC"
				  << "\tOriginal_Pval"
				  << "\tGC"
				  << "\tBONF"
				  << "\tHOLM"
				  << "\tSIDAK_SS"
				  << "\tSIDAK_SD"
				  << "\tFDR_BH"
				  << "\tFDR_BY"
				  << endl;
			opts::addHeader(fcomp, "CALC");
			opts::addHeader(fcomp, "Original_Pval");
			opts::addHeader(fcomp, "GC");
			opts::addHeader(fcomp, "BONF");
			opts::addHeader(fcomp, "HOLM");
			opts::addHeader(fcomp, "SIDAK_SS");
			opts::addHeader(fcomp, "SIDAK_SD");
			opts::addHeader(fcomp, "FDR_BH");
			opts::addHeader(fcomp, "FDR_BY");


		MultComparison mc(options);
		vector<int> tcnt;
		mc.calculate(chis, tcnt);

		prev_base = 0;
		prev_chrom = -1;
		for(unsigned int m = 0; m < data_set->num_loci(); m++){
			Marker* mark = data_set->get_locus(m);

			if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
				COMP << mark->toString() << "\tLINREG\t";
				COMP << main_pvals[m] << "\t"
				<< mc.get_genomic_control(m) << "\t"
				<< mc.get_bonferroni(m) << "\t"
				<< mc.get_holm(m) << "\t"
				<< mc.get_sidak_single_step(m) << "\t"
				<< mc.get_sidak_step_down(m) << "\t"
				<< mc.get_fdr_bh(m) << "\t"
				<< mc.get_fdr_by(m)
				<< endl;

			}
		}
		if(COMP.is_open()){
			COMP.close();
		}


	}

	prev_base = 0;
	prev_chrom = -1;
	for(unsigned int m = 0; m < data_set->num_loci(); m++){
		Marker* mark = data_set->get_locus(m);

		if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
			doFilter(mark, main_pvals[m]);
		}
	}

	if(lrout.is_open()){
		lrout.close();
	}
}
