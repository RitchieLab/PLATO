/**********************************************************************************
*                       Lapis Input Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates Lapis input files
*
*
*
*File: LAPISOutput.cc
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
#include <iomanip>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include <time.h>
#include "ProcessLAPISOutput.h"
#include <General.h>
#include <Helpers.h>
using namespace Methods;

string ProcessLAPISOutput::stepname = ProcessLAPISOutput::doRegister("output-lapis");

void ProcessLAPISOutput::FilterSummary(){
}

void ProcessLAPISOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}

}

void ProcessLAPISOutput::filter(){}

void ProcessLAPISOutput::process(DataSet* ds){
	data_set = ds;

	LAPISOutput lapis;
	lapis.setOrder(this->order);
	lapis.setOverwrite(this->overwrite);
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
			lapis.setOptions(options);
			lapis.calculate(&ds);
			options.setOut(tempout);
			ds.clear_all();
		}
	}
	else{
	lapis.setOptions(options);
	lapis.calculate(data_set);
	}
}

#ifdef PLATOLIB
void ProcessLAPISOutput::dump2db(){}
void ProcessLAPISOutput::create_tables(){}
void ProcessLAPISOutput::run(DataSetObject* ds)
{
	#ifdef WIN
		options.setOverrideOut(projectPath + "\\" + "input_lapis");
	#else
		options.setOverrideOut(projectPath + "/" + "input_lapis");
	#endif
	process(ds);
}
#endif

#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
