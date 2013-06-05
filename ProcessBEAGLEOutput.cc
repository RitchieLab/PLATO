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


#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <iomanip>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include "ProcessBEAGLEOutput.h"
#include <General.h>
#include <Helpers.h>

using namespace Methods;

//define the PlatoLib namespace for use with Plato as a library
#ifdef PLATOLIB
namespace PlatoLib
{
#endif
#ifdef PLATOLIB
//Constructor to allow for use as Library with Plato-viewer
ProcessBEAGLEOutput::ProcessBEAGLEOutput(string bn, int pos, Database* pdb, string projPath)
{
	name = "Output Beagle";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
	projectPath = projPath;
}
#endif

/*
 *Function: FilterSummary
 *Description:
 *Not used.
 */
void ProcessBEAGLEOutput::FilterSummary(){
}

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
 * Function: filter
 * Description:
 * not used
 */
void ProcessBEAGLEOutput::filter(){
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

#ifdef PLATOLIB
void ProcessBEAGLEOutput::dump2db(){};
/*
 * Function: create_tables
 * Description: create db tables for use in process
 * not used
 */
void ProcessBEAGLEOutput::create_tables(){};
/*
 * Method: ProcessBEAGLEOutput::run
 * used by Plato-viewer
 */
void ProcessBEAGLEOutput::run(DataSetObject* ds)
{
	data_set = ds;

	Methods::BEAGLEOutput bo;
#ifdef WIN
	//use windows version of project path + batchname
	options.setOverrideOut(projectPath + "\\" + options.convertString(batchname + "_" + name + "_" + getString<int>(position)));
#else
	//use non-windows version of project path + batchname
	options.setOverrideOut(projectPath + "/" + options.convertString(batchname + "_" + name + "_" + getString<int>(position)));
#endif
	bo.setOrder(position);
	bo.setOverwrite(true);
	bo.setOptions(options);
	bo.calculate(data_set);

	filenames = bo.get_filenames();
}
#endif
#ifdef PLATOLIB
}
#endif

