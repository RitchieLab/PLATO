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
*File: PartialOutput.cc
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
#include "ProcessBINOutput.h"
#include "Chrom.h"
#include <General.h>
#include <Helpers.h>
using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif
#ifdef PLATOLIB
ProcessBINOutput::ProcessBINOutput(string bn, int pos, Database* pdb, string projPath)
{
	name = "Output Bin";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
	projectPath = projPath;
}
#endif

ProcessBINOutput::~ProcessBINOutput(){}
void ProcessBINOutput::filter(){}
void ProcessBINOutput::FilterSummary(){}

void ProcessBINOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessBINOutput::process(DataSet* ds){
	data_set = ds;

	BINOutput BIN;
	BIN.setOrder(this->order);
	BIN.setOverwrite(this->overwrite);
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
			BIN.setOptions(options);
			BIN.calculate(&ds);
			options.setOut(tempout);
			ds.clear_all();
		}
	}
	else{
	BIN.setOptions(options);
	BIN.calculate(data_set);
	}
}
#ifdef PLATOLIB
void ProcessBINOutput::dump2db(){}

void ProcessBINOutput::create_tables(){}

void ProcessBINOutput::run(DataSetObject* ds)
{
	data_set = ds;

	Methods::BINOutput Bin;
#ifdef WIN
	//use windows version of project path + batchname
	options.setOverrideOut(projectPath + "\\" + options.convertString(batchname + "_" + name + "_" + getString<int>(position)));
#else
	//use non-windows version of project path + batchname
	options.setOverrideOut(projectPath + "/" + options.convertString(batchname + "_" + name + "_" + getString<int>(position)));
#endif
	Bin.setOrder(position);
	Bin.setOverwrite(true);
	Bin.setOptions(options);
	Bin.calculate(data_set);

	filenames = Bin.get_filenames();
}
#endif

#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
