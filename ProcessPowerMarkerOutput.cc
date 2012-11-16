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
*File: PowerMarkerOutput.cc
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
#include "ProcessPowerMarkerOutput.h"
#include "Chrom.h"
#include <General.h>
#include <Helper.h>
using namespace Methods;
string ProcessPowerMarkerOutput::stepname = "output-powermarker";


void ProcessPowerMarkerOutput::FilterSummary(){
}

void ProcessPowerMarkerOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}

}

void ProcessPowerMarkerOutput::filter(){
}


void ProcessPowerMarkerOutput::process(DataSet* ds){
	data_set = ds;

	PowerMarkerOutput pmo;
	pmo.setOrder(this->order);
	pmo.setOverwrite(this->overwrite);
	pmo.setOptions(options);
	pmo.calculate(data_set);

}
