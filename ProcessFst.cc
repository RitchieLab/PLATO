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
#include <MultComparison.h>
#include "ProcessFst.h"
#include <Options.h>
#include <General.h>
#include <Helper.h>
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
using namespace Methods;

string ProcessFst::stepname = "fst";


void ProcessFst::FilterSummary() {

	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int> (
			opts::_MARKERS_WORKING_ - orig_num_markers) + " (" + getString<
			float> (((float) (opts::_MARKERS_WORKING_ - orig_num_markers)
			/ (float) opts::_MARKERS_WORKING_) * 100.0) + "%) of " + getString<
			int> (opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessFst::PrintSummary() {
	int msize = data_set->num_loci();
	for (int m = 0; m < msize; m++) {
		data_set->get_locus(m)->setFlag(false);
	}

}

void ProcessFst::filter() {
}

void ProcessFst::doFilter(Methods::Marker* mark, double value) {
	if (options.doThreshMarkersLow() || options.doThreshMarkersHigh()) {
		if (mark->isEnabled()){// && !mark->isFlagged()) {
			bool inc = false;
			if (options.doThreshMarkersLow() && dLess(value,
					options.getThreshMarkersLow())) {
				mark->setEnabled(false);
				inc = true;
			}
			if (options.doThreshMarkersHigh() && dGreater(value,
					options.getThreshMarkersHigh())) {
				mark->setEnabled(false);
				inc = true;
			}
			if (inc) {
				orig_num_markers++;
			}
		}
	}
}

void ProcessFst::process(DataSet* ds) {
	data_set = ds;

	//check if new covariate file is listed...or covariate name.
	//create vector of covariate indexes to use if specified.

	string fname = opts::_OUTPREFIX_ + "fst" + options.getOut() + ".txt";
	if (!overwrite) {
		fname += "." + getString<int> (order);
	}
	ofstream eout(fname.c_str());
	if (!eout) {
		opts::printLog("Error opening " + fname + "!  Exiting!\n");
		throw MethodException("");
	}
	eout.precision(4);

	//    DataSet* trimmed_data = new DataSet();
	//    trimmed_data->set_markers(ds->get_markers());
	//    trimmed_data->set_covariates(ds->get_covariates());
	//    trimmed_data->set_traits(ds->get_traits());
	//	for(int s = 0; s < data_set->num_inds(); s++){
	//		Sample* samp = data_set->get_sample(s);
	//		if(samp->isEnabled() && (samp->getPheno() == 1 || samp->getPheno() == 2)){
	//			trimmed_data->add_ind(samp);
	//		}
	//	}
	//	trimmed_data->recreate_family_vector();
	//	trimmed_data->set_affection_vectors();


	Fst fst;
	fst.set_parameters(&options);
	fst.resetDataSet(ds);
	int prev_base = 0;
	int prev_chrom = -1;

	vector<Marker*> good_markers = findValidMarkers(ds->get_markers(), &options);
	int msize = good_markers.size();

	eout << "Chrom\trsID\tProbeID\tbploc\tFSTWC\tFSTRH\n";//\tFSTHM\n";
		opts::addFile("Marker", stepname, fname);

		opts::addHeader(fname, "FSTWC");
		opts::addHeader(fname, "FSTRH");
//		opts::addHeader(fname, "FSTHM");


		for (int m = 0; m < (int) msize; m++){//ds->num_loci(); m++) {
			Marker* mark = good_markers[m];//ds->get_locus(m);
			if (mark->isEnabled()){// && isValidMarker(mark, &options, prev_base,prev_chrom)) {
				fst.calculate(m);
				eout << mark->toString() << "\t" << fst.getFst() << "\t" << fst.getFstRH()// << "\t" << fst.getFstHM()
						<< endl;

			}
		}

	if (options.doMultCompare()) {
		/*		string fcomp = opts::_OUTPREFIX_ + "Fst_comparisons" + options.getOut() + ".txt";
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
		 if(ds->get_locus(0)->getDetailHeaders().size() > 0){
		 COMP << "\t" << ds->get_locus(0)->getDetailHeaders();
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
		 for(unsigned int m = 0; m < ds->num_loci(); m++){
		 Marker* mark = ds->get_locus(m);

		 if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
		 COMP << mark->toString() << "\tFst\t";
		 COMP << pvals[m] << "\t"
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

		 */
	}
	/*
	 prev_base = 0;
	 prev_chrom = -1;
	 for(unsigned int m = 0; m < ds->num_loci(); m++){
	 Marker* mark = ds->get_locus(m);

	 if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
	 doFilter(mark, pvals[m]);
	 }
	 }
	 */
	eout.close();
}
