/**********************************************************************************
*                       Structure Output Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates STRUCTURE input files.
*
*
* Files generated:
*
*File: STRUCTOutput.cc
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
#include "ProcessSTRUCTOutput.h"
#include <General.h>
#include <Helper.h>
using namespace Methods;

void ProcessSTRUCTOutput::FilterSummary(){
}

void ProcessSTRUCTOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int m = 0; m < msize; m++){
		data_set->get_locus(m)->setFlag(false);
	}

}

void ProcessSTRUCTOutput::filter(){
}

void ProcessSTRUCTOutput::process(DataSet* ds){
	data_set = ds;

	STRUCTOutput str;
	str.setOrder(this->order);
	str.setOverwrite(this->overwrite);
	str.setOptions(options);
	str.calculate(data_set);
}

