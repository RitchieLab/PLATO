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
#include "ProcessFBATOutput.h"
#include <General.h>
#include <Helpers.h>
using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

ProcessFBATOutput::ProcessFBATOutput(string bn, int pos, Database* pdb, string projPath)
{
	name = "Output FBAT";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
	projectPath = projPath;
}

void ProcessFBATOutput::FilterSummary(){}

void ProcessFBATOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessFBATOutput::filter(){}

void ProcessFBATOutput::process(DataSet* ds){
	data_set = ds;

	FBATOutput fbat;
	fbat.setOrder(this->order);
	#ifdef PLATOLIB
		fbat.setOverwrite(true);
	#else
		fbat.setOverwrite(this->overwrite);
	#endif
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
				FixOutputName();
			#endif
			options.setOut("_random_set_" + getString<int>(i + 1) + tempout);
			fbat.setOptions(options);
			fbat.calculate(&ds);
			options.setOut(tempout);
			ds.clear_all();
		}
	}
	else{
	#ifdef PLATOLIB
		FixOutputName();
	#endif
	fbat.setOptions(options);
	fbat.calculate(data_set);
	}
}

#ifdef PLATOLIB
void ProcessFBATOutput::dump2db(){}

void ProcessFBATOutput::create_tables(){}

void ProcessFBATOutput::run(DataSetObject* ds)
{
	process(ds);
}

void ProcessFBATOutput::FixOutputName()
{

	#ifdef WIN
			//use windows version of project path + batchname
			options.setOverrideOut(projectPath + "\\" + batchname + "_" + name +  "_" + getString<int>(position));
	#else
			//use non-windows version of project path + batchname
			options.setOverrideOut(projectPath + "/" + batchname + "_" + name +  "_" + getString<int>(position));
	#endif

}
#endif
#ifdef PLATOLIB
}//end namespace PlatoLib
#endif

