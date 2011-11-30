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
*File: PHASEOutput.cc
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
#include "ProcessPHASEOutput.h"
#include <General.h>
#include <Helpers.h>
using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif

ProcessPHASEOutput::ProcessPHASEOutput(string bn, int pos, Database* pdb, string projPath)
{
    name = "Output Phase";
    batchname = bn;
    position = pos;
    hasresults = false;
    db = pdb;
    projectPath = projPath;
}

void ProcessPHASEOutput::FilterSummary(){}

void ProcessPHASEOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessPHASEOutput::filter(){}

void ProcessPHASEOutput::process(DataSet* ds)
{
	data_set = ds;

	PHASEOutput phase;
	phase.setOrder(this->order);
	#ifdef PLATOLIB
		phase.setOverwrite(true);
	#else
		phase.setOverwrite(this->overwrite);
	#endif
	if(options.getRandSamps() > 0 || options.getSetsSamps() > 0)
	{
		vector<vector<Sample*> > sample_sets = Helpers::generateSampleSets(data_set, &options);
		for(int i = 0; i < (int)sample_sets.size(); i++)
		{
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
			phase.setOptions(options);
			phase.calculate(&ds);
			options.setOut(tempout);
			ds.clear_all();
		}
	}
	else
	{
		phase.setOptions(options);
		phase.calculate(data_set);
	}
	#ifdef PLATOLIB
		filenames = phase.get_filenames();
	#endif
}//end method Process(DataSet* ds)

#ifdef PLATOLIB
void ProcessPHASEOutput::create_tables(){}

void ProcessPHASEOutput::dump2db(){}

void ProcessPHASEOutput::run(DataSetObject* ds)
{
	#ifdef WIN
		options.setOverrideOut(projectPath + "\\" + batchname + "_" + name + "_" + getString<int>(position));
	#else
		options.setOverrideOut(projectPath + "/" + batchname + "_" + name + "_" + getString<int>(position));
	#endif
		process(ds);
}
#endif

#ifdef PLATOLIB
}//end namespace PlatoLib
#endif


