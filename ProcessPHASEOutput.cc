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
#include <Helper.h>


void ProcessPHASEOutput::FilterSummary(){
}

void ProcessPHASEOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}

}

void ProcessPHASEOutput::filter(){
}

void ProcessPHASEOutput::process(DataSet* ds){
	data_set = ds;

	PHASEOutput phase;
	phase.setOrder(this->order);
	phase.setOverwrite(this->overwrite);
	phase.setOptions(options);
	phase.calculate(data_set);
}



