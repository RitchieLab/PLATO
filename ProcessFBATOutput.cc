/**********************************************************************************
*                       FBAT Output Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates FBAT input files.
*
*
*
*File: FBATOutput.cc
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
#include "ProcessFBATOutput.h"
#include <General.h>
#include <Helper.h>
using namespace Methods;

void ProcessFBATOutput::FilterSummary(){
}

void ProcessFBATOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessFBATOutput::filter(){
}

void ProcessFBATOutput::process(DataSet* ds){
	data_set = ds;

	FBATOutput fbat;
	fbat.setOrder(this->order);
	fbat.setOverwrite(this->overwrite);
	fbat.setOptions(options);
	fbat.calculate(data_set);
}


