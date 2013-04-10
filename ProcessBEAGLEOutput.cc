/**********************************************************************************
*                       ProcessBEAGLE Output Module
*
* Written by: Justin Giles
*             Vanderbilt University
*             Center for Human Genetics Research
*
* Outputs ProcessBEAGLE input files.
*
*
*File: ProcessBEAGLEOutput.cc
**********************************************************************************/

#include "ProcessBEAGLEOutput.h"
#include <BEAGLEOutput.h>
#include <General.h>
#include <Helpers.h>

using std::string;
using std::vector;

using Methods::DataSet;
using Methods::BEAGLEOutput;
using Methods::Helpers;
using Methods::Sample;


const string ProcessBEAGLEOutput::stepname = ProcessBEAGLEOutput::doRegister("output-beagle");

/*
 * Function: PrintSummary
 * Description:
 * Resets marker flags
 */
void ProcessBEAGLEOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}

}

/*
 * Function: process
 * Description:
 * Main function for producing output files.
 * Creates output files, one per chromosome.
 */
void ProcessBEAGLEOutput::process(DataSet* ds){
	data_set = ds;

	BEAGLEOutput bo;
	bo.setOrder(this->order);
	bo.setOverwrite(this->overwrite);
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
			bo.setOptions(options);
			bo.calculate(&ds);
			options.setOut(tempout);
			ds.clear_all();
		}
	}
	else{


	bo.setOptions(options);

	bo.calculate(data_set);
	}
}


