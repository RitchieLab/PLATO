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
#include "ProcessKinship.h"
#include <Options.h>
#include <General.h>
#include <Helpers.h>
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif

void ProcessKinship::FilterSummary() {

	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int> (
			opts::_MARKERS_WORKING_ - orig_num_markers) + " (" + getString<
			float> (((float) (opts::_MARKERS_WORKING_ - orig_num_markers)
			/ (float) opts::_MARKERS_WORKING_) * 100.0) + "%) of " + getString<
			int> (opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessKinship::PrintSummary() {
	int msize = data_set->num_loci();
	for (int m = 0; m < msize; m++) {
		data_set->get_locus(m)->setFlag(false);
	}

}

void ProcessKinship::filter() {
}

void ProcessKinship::doFilter(Methods::Marker* mark, double value) {
	if (options.doThreshMarkersLow() || options.doThreshMarkersHigh()) {
		if (mark->isEnabled() && !mark->isFlagged()) {
			bool inc = false;
			if (options.doThreshMarkersLow() && Helpers::dLess(value,
					options.getThreshMarkersLow())) {
				mark->setEnabled(false);
				inc = true;
			}
			if (options.doThreshMarkersHigh() && Helpers::dGreater(value,
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

void ProcessKinship::process(DataSet* ds) {
	data_set = ds;

	//check if new covariate file is listed...or covariate name.
	//create vector of covariate indexes to use if specified.

	string fname = opts::_OUTPREFIX_ + "Kinship" + options.getOut() + ".txt";
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


	Kinship kin;
//	kin.set_parameters(&options);
	kin.resetDataSet(ds);
//	int prev_base = 0;
	//int prev_chrom = -1;

		eout << "FamID\tIndID\tCMP_IndID\tKinshipCoef\n";

		for(unsigned int i = 0; i < data_set->num_pedigrees(); i++){
			kin.calculate(i);
			map<string, double> coefs = kin.getCoefficients();
			map<string, double>::iterator iter;
			for(iter = coefs.begin(); iter != coefs.end(); iter++){
				string keys = iter->first;
				vector<string> vals;
				General::Tokenize(keys, vals, " ");

				int a = atoi(vals[0].c_str());
				int b = atoi(vals[1].c_str());

				Sample* one = data_set->get_sample(a);
				Sample* two = data_set->get_sample(b);

				eout << one->toString() << "\t" << two->getInd() << "\t" << iter->second << endl;
			}
		}
	eout.close();
}
#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
