/**********************************************************************************
*                       Marker Genotype Efficiency Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Example modlue framework
*
* Files generated:
*
*File: ExampleModule.cc
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
#include "ExampleModule.h" //////CHANGE TO REAL MODULE NAME
#include <Options.h>
#include <General.h>
#include <Helper.h>
#include <cdflib.h>
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"

//FilterSummary()
//used to output the number of elements remaining after process
void ExampleModule::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;
}

//PrintSummary()
//used to output results after processing the data
void ExampleModule::PrintSummary(){
	//output for process goes here
	//
	cout << "Hi, I'm in the ExampleModule PrintSummary method\n";
}

//filter()
//used to filter markers, samples, families after processing the data
void ExampleModule::filter(){
	//if markers, samples, or families need to be filtered, do it here (ie: use "setEnabled(flase);" method in classes)
	//
	cout << "Hi, I'm in the ExampleModule filter method\n";
}

//process()
//main method to get the process going and to the work
void ExampleModule::process(DataSet* ds){
	data_set = ds;
	//main processing of step goes here.
	cout << "Hi, I'm in the ExampleModule process method\n";
	cout << "I'm going to process " << data_set->num_loci() << " markers, " << data_set->num_pedigrees() << " families, and " << data_set->num_inds() << " samples!\n";

}

