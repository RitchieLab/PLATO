/**********************************************************************************
*                        GRR Output Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates GRR input files.
*
*
*File: GRROutput.cc
**********************************************************************************/

#include "ProcessGRROutput.h"
#include <GRROutput.h>
#include <Helpers.h>

using Methods::GRROutput;
using Methods::DataSet;
using std::string;
using std::vector;
using Methods::Sample;
using Methods::Helpers;

const string ProcessGRROutput::stepname = ProcessGRROutput::doRegister("output-grr");

void ProcessGRROutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessGRROutput::process(DataSet* ds){
	data_set = ds;

	GRROutput grr;
	grr.setOverwrite(this->overwrite);
	grr.setOrder(this->order);
	if(options.getRandSamps() > 0 || options.getSetsSamps() > 0){
		vector<vector<Sample*> > sample_sets = Helpers::generateSampleSets(data_set, &options);
		for(int i = 0; i < (int)sample_sets.size(); i++){
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
			options.setOut("_random_set_" + getString<int>(i + 1) + tempout);
			grr.setOptions(options);
			grr.calculate(&ds);
			options.setOut(tempout);
			ds.clear_all();
		}
	}
	else
	{
		grr.setOptions(options);
		grr.calculate(data_set);
	}
}//end method process(DataSet*)

