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
#include "ProcessEpistasis.h"
#include <Epistasis.h>

using Methods::DataSet;
using Methods::Epistasis;
using std::string;

const string ProcessEpistasis::stepname = ProcessEpistasis::doRegister("epistasis");

void ProcessEpistasis::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessEpistasis::process(DataSet* ds){
	data_set = ds;

	Epistasis epi;
	epi.setOrder(this->order);
	epi.setOverwrite(this->overwrite);
	epi.setOptions(options);
	epi.calculate(data_set);
}

