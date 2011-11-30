/**********************************************************************************
*                       QTDT Input Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
*
*
* Files generated:
*
*File: QTDTOutput.cc
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
#include "ProcessQTDTOutput.h"
#include <General.h>
#include <Helpers.h>
using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif

void ProcessQTDTOutput::FilterSummary(){
}

void ProcessQTDTOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessQTDTOutput::filter(){
}

void ProcessQTDTOutput::process(DataSet* ds){
	data_set = ds;

	QTDTOutput qtdt;
	qtdt.setOverwrite(this->overwrite);
	qtdt.setOrder(this->order);
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
			qtdt.setOptions(options);
			qtdt.calculate(&ds);
			options.setOut(tempout);
			ds.clear_all();
		}
	}
	else{
	qtdt.setOptions(options);
	qtdt.calculate(data_set);
	}
}
#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
