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
*File: MDROutput.cc
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
#include "ProcessMDROutput.h"
#include <General.h>
#include <Helper.h>
using namespace Methods;

void ProcessMDROutput::FilterSummary(){
}

void ProcessMDROutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessMDROutput::filter(){
}

void ProcessMDROutput::process(DataSet* ds){
	data_set = ds;
	MDROutput mdr;
	mdr.setOrder(this->order);
	mdr.setOverwrite(this->overwrite);
	mdr.setOptions(options);
	mdr.calculate(data_set);

}
