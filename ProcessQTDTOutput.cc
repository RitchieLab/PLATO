/**********************************************************************************
*                       QTDT Input Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
*
*
* Files generated:
*
*File: QTDTOutput.cc
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
#include "ProcessQTDTOutput.h"
#include <General.h>
#include <Helper.h>
using namespace Methods;

void ProcessQTDTOutput::FilterSummary(){
}

void ProcessQTDTOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessQTDTOutput::filter(){
}

void ProcessQTDTOutput::process(DataSet* ds){
	data_set = ds;

	QTDTOutput qtdt;
	qtdt.setOverwrite(this->overwrite);
	qtdt.setOrder(this->order);
	qtdt.setOptions(options);
	qtdt.calculate(data_set);
}

