/**********************************************************************************
*                       Concordance Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Performs concordance check with another set of ped files.
*
*
* Files generated:
*
*File: Concordance.cc
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
#include "ProcessConcordance.h"
#include <Options.h>
#include <General.h>
#include <Helpers.h>
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
using namespace Methods;
string ProcessConcordance::stepname = "concordance";

void ProcessConcordance::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessConcordance::PrintSummary(){

}

void ProcessConcordance::filter(){
}

void ProcessConcordance::process(DataSet* ds){
	data_set = ds;
	Concordance con(data_set);
	con.setOrder(this->order);
	con.setOverwrite(this->overwrite);
	con.setOptions(options);
	con.calculate();

}

