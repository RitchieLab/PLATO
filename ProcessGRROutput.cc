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
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include "ProcessGRROutput.h"
#include <General.h>
#include <Helpers.h>
using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif
#ifdef PLATOLIB
ProcessGRROutput::ProcessGRROutput(string bn, int pos, Database* pdb, string projPath)
{
	name = "Output GRR";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
	projectPath = projPath;
}
#endif

void ProcessGRROutput::FilterSummary(){}

void ProcessGRROutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessGRROutput::filter(){
}

void ProcessGRROutput::process(DataSet* ds){
	data_set = ds;

	GRROutput grr;
#ifdef PLATOLIB
	grr.setOverwrite(true);
#else
	grr.setOverwrite(this->overwrite);
#endif
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
#ifdef PLATOLIB
			FixOutputName(i, tempout);
#else
			options.setOut("_random_set_" + getString<int>(i + 1) + tempout);
#endif
			grr.setOptions(options);
			grr.calculate(&ds);
			options.setOut(tempout);
#ifdef PLATOLIB
			filenames = grr.get_filenames();
#endif
			ds.clear_all();
		}
	}
	else
	{
#ifdef PLATOLIB
		FixOutputName(1, options.getOut());
#endif
		grr.setOptions(options);
		grr.calculate(data_set);
	}
}//end method process(DataSet*)
#ifdef PLATOLIB
void ProcessGRROutput::dump2db(){}

void ProcessGRROutput::create_tables(){}

void ProcessGRROutput::run(DataSetObject* ds)
{
	process(ds);
}

void ProcessGRROutput::FixOutputName(int i, string tempout)
{
	#ifdef WIN
		options.setOverrideOut(projectPath + "\\" + "_random_set_" + getString<int>(i + 1) + tempout);
	#else
		options.setOverrideOut(projectPath + "/" + "_random-set_" + getString<int>(i + 1) + tempout);
	#endif
}
#endif

#ifdef PLATOLIB
};//end namespace PlatoLib
#endif
