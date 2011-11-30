/**********************************************************************************
*                        GRR Output Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates GRR input files.
*
*
*File: GRROutput.cc
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
#include "ProcessGRROutput.h"
#include <General.h>
#include <Helper.h>


void ProcessGRROutput::FilterSummary(){
}

void ProcessGRROutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessGRROutput::filter(){
}

void ProcessGRROutput::process(DataSet* ds){
	data_set = ds;

	GRROutput grr;
	grr.setOverwrite(this->overwrite);
	grr.setOrder(this->order);
	grr.setOptions(options);
	grr.calculate(data_set);
}
