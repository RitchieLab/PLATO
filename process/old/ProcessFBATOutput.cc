/**********************************************************************************
 *                       FBAT Output Module
 *
 *Written by: Justin Giles
 *            Vanderbilt University
 *            Center for Human Genetics Research
 *
 * Generates FBAT input files.
 *
 *
 *
 *File: FBATOutput.cc
 **********************************************************************************/

#include "ProcessFBATOutput.h"
#include <FBATOutput.h>
#include <Helpers.h>

using std::string;
using std::getString;
using Methods::DataSet;
using Methods::Sample;
using Methods::FBATOutput;
using Methods::Helpers;

const string ProcessFBATOutput::stepname = ProcessFBATOutput::doRegister("output-fbat");

void ProcessFBATOutput::PrintSummary() {
	int msize = data_set->num_loci();

	for (int i = 0; i < msize; i++) {
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessFBATOutput::process(DataSet* ds) {
	data_set = ds;

	FBATOutput fbat;
	fbat.setOrder(this->order);
	fbat.setOverwrite(this->overwrite);
	if (options.getRandSamps() > 0 || options.getSetsSamps() > 0) {
		vector<vector<Sample*> > sample_sets = Helpers::generateSampleSets(
				data_set, &options);
		for (int i = 0; i < (int) sample_sets.size(); i++) {
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
			fbat.setOptions(options);
			fbat.calculate(&ds);
			options.setOut(tempout);
			ds.clear_all();
		}
	} else {
		fbat.setOptions(options);
		fbat.calculate(data_set);
	}
}

