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
#include "ProcessEigenstratOutput.h"
#include "Chrom.h"
#include <General.h>
#include <Helpers.h>
using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

string ProcessEigenstratOutput::stepname = "output-eigenstrat";
#ifdef PLATOLIB
ProcessEigenstratOutput::ProcessEigenstratOutput(string bn, int pos, Database* pdb, string projPath)
{
	name = "Output Eigenstrature";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
	projectPath = projPath;
}
#endif

void ProcessEigenstratOutput::FilterSummary(){}

void ProcessEigenstratOutput::filter(){}

void ProcessEigenstratOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessEigenstratOutput::process(DataSet* ds){
	data_set = ds;

	EigenstratOutput ped;
	ped.setOrder(this->order);
	ped.setOverwrite(this->overwrite);
	string tempout = options.getOut();

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
			#ifdef PLATOLIB
				FixOutputName(i, tempout);
			#endif
			ped.setOptions(options);
			ped.calculate(&ds);
			#ifdef PLATOLIB
				filenames = ped.get_filenames();
			#endif
			options.setOut(tempout);
			ds.clear_all();
		}//end for
	}//end if
	else{
	#ifdef PLATOLIB
		FixOutputName(1, tempout);
	#endif
		ped.setOptions(options);
		ped.calculate(data_set);
	#ifdef PLATOLIB
		filenames = ped.get_filenames();
	#endif
	}
}//end method process(DataSet*)

#ifdef PLATOLIB
void ProcessEigenstratOutput::dump2db(){}

void ProcessEigenstratOutput::create_tables(){}

void ProcessEigenstratOutput::run(DataSetObject* ds)
{
	process(ds);
}

void ProcessEigenstratOutput::FixOutputName(int i, string tempout)
{
	#ifdef WIN
	//use windows version of project path + batchname
	options.setOverrideOut(projectPath + "\\" + "_random_set_" + getString<int>(i + 1) + tempout);
	#else
	//use non-windows version of project path + batchname
	options.setOverrideOut(projectPath + "/" + "_random_set_" + getString<int>(i + 1) + tempout);
	#endif
}
#endif

#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
