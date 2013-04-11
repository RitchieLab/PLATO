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

#include "ProcessEigenstratOutput.h"
#include <EigenstratOutput.h>
#include <DataSet.h>
#include <Helpers.h>

using Methods::DataSet;
using Methods::Sample;
using Methods::Helpers;
using Methods::EigenstratOutput;
using std::string;

const string ProcessEigenstratOutput::stepname = ProcessEigenstratOutput::doRegister("output-eigenstrat");

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
			ped.setOptions(options);
			ped.calculate(&ds);
			options.setOut(tempout);
			ds.clear_all();
		}//end for
	}//end if
	else{
		ped.setOptions(options);
		ped.calculate(data_set);
	}
}//end method process(DataSet*)
