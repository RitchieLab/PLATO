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
#include "ProcessFilterProcess.h"
#include "Chrom.h"
#include <General.h>
#include <Helpers.h>
using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

void ProcessFilterProcess::FilterSummary(){
}

void ProcessFilterProcess::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}

}

void ProcessFilterProcess::filter(){
}


void ProcessFilterProcess::process(DataSet* ds){
	data_set = ds;

	string fcomp = opts::_OUTPREFIX_ + "plato_filters" + options.getOut() + ".txt";

	if(!options.hasFilterProcessConfig()){
		opts::printLog("No configuration file specified for filter-process!  Please use -config option.");
		throw MethodException("No configuration file specified for filter-process!  Please use -config option.");
	}
	FilterProcess filter;
	filter.setOutputName(fcomp);
	filter.calculate(data_set, options.getFilterProcessConfig());


}
#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
