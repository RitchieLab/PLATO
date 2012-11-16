/**********************************************************************************
*                       Lapis Input Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates Lapis input files
*
*
*
*File: LAPISOutput.cc
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
#include <time.h>
#include "ProcessLAPISOutput.h"
#include <General.h>
#include <Helper.h>
using namespace Methods;
string ProcessLAPISOutput::stepname = "output-lapis";

void ProcessLAPISOutput::FilterSummary(){
}

void ProcessLAPISOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}

}

void ProcessLAPISOutput::filter(){
}

void ProcessLAPISOutput::process(DataSet* ds){
	data_set = ds;

	LAPISOutput lapis;
	lapis.setOrder(this->order);
	lapis.setOverwrite(this->overwrite);
	lapis.setOptions(options);
	lapis.calculate(data_set);


}

