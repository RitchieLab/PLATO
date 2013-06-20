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
 *File: PDT2Output.cc
 **********************************************************************************/

#include "ProcessPDT2Output.h"
#include <PDT2Output.h>

#include <vector>

#include <Marker.h>
#include <Sample.h>
#include <Helpers.h>

using std::string;
using std::vector;
using Methods::Helpers;
using Methods::DataSet;
using Methods::Marker;
using Methods::Sample;

using Methods::PDT2Output;

const string ProcessPDT2Output::stepname = ProcessPDT2Output::doRegister("output-pdt2");

void ProcessPDT2Output::PrintSummary() {
	int msize = data_set->num_loci();

	for (int i = 0; i < msize; i++) {
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessPDT2Output::process(DataSet* ds) {
	data_set = ds;

	PDT2Output pdt2;
	if (options.getRandSamps() > 0 || options.getSetsSamps() > 0) {
		vector<vector<Sample*> > sample_sets = Helpers::generateSampleSets(
				data_set, &options);
		for (int i = 0; i < (int) sample_sets.size(); i++) {
			//			cout << "Sample vect size: " << sample_sets[i].size() << endl;
			DataSet ds;
			ds.set_samples(&sample_sets[i]);
			ds.set_markers(data_set->get_markers());
			ds.set_affection_vectors();
			ds.set_missing_value(data_set->get_missing_value());
			ds.set_covariates(data_set->get_covariates());
			ds.set_traits(data_set->get_traits());
			ds.recreate_family_vector();
			string tempout = options.getOut();
			options.setOut("_random_set_" + getString<int> (i + 1) + tempout);
			pdt2.setOptions(options);
			pdt2.calculate(&ds);
			options.setOut(tempout);
			ds.clear_all();
		}
	} else {
		pdt2.setOptions(options);
		pdt2.calculate(data_set);
	}
}//end method process(DataSet*)


