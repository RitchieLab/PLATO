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
#include "ProcessClusterMissing.h"
#include "Chrom.h"
#include <General.h>
using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

void ProcessClusterMissing::FilterSummary(){
}

void ProcessClusterMissing::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}

}

void ProcessClusterMissing::filter(){}

void ProcessClusterMissing::process(DataSet* ds){
	data_set = ds;

	ClusterMissing miss(data_set);
	miss.setOrder(this->order);
	miss.setOverwrite(this->overwrite);
	miss.setOptions(options);
	miss.calculate(opts::_OUTPREFIX_, (options.getOut() + ".txt"));

}
#ifdef PLATOLIB
void ProcessClusterMissing::dump2db(){}
void ProcessClusterMissing::create_tables(){}
void ProcessClusterMissing::run(DataSetObject* ds){}
#endif

#ifdef PLATOLIB
};//end namespace PlatoLib
#endif
