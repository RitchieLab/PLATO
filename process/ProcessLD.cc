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
*File: LD.cc
**********************************************************************************/


#include "ProcessLD.h"
#include <LD.h>

using std::string;
using Methods::LD;
using Methods::DataSet;

const string ProcessLD::stepname = ProcessLD::doRegister("ld");

void ProcessLD::PrintSummary(){
	int msize = data_set->num_loci();
	for(int m = 0; m < msize; m++){
		data_set->get_locus(m)->setFlag(false);
	}
}

void ProcessLD::process(DataSet* ds){
	data_set = ds;

	LD ld;
	ld.setOptions(options);
	ld.setOrder(this->order);
	ld.setOverwrite(this->overwrite);
	ld.calculate(data_set);

	orig_num_markers = ld.getOrigNumMarkers();
}

