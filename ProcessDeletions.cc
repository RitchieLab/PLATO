/**********************************************************************************
*                       Deletion detection
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Performs deletion detection based on mendelian error as described in
* Conrad et al.
*
*
*File: Deletions.cc
**********************************************************************************/


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include "Globals.h"
#include <string>
#include <list>
#include <map>
#include <algorithm>
#include <bitset>
#include "ProcessDeletions.h"
#include "Chrom.h"
#include <General.h>
#include <Helper.h>

using namespace std;
using namespace Methods;
string ProcessDeletions::stepname="deletions";

void ProcessDeletions::setThreshold(string thresh){
	options.setUp(thresh);
}

/*
 * Function: FilterSummary
 * Description:
 * Outputs remaining markers and families
 */
void ProcessDeletions::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::printLog("Families Passed:\t" + getString<int>(opts::_FAMILIES_WORKING_ - orig_num_families) + " (" +
        getString<float>(((float)(opts::_FAMILIES_WORKING_ - orig_num_families) / (float)opts::_FAMILIES_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_FAMILIES_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;
	opts::_FAMILIES_WORKING_ -= orig_num_families;

}


/*
 * Function: PrintSummary
 * Description:
 * Resets marker flags
 */
void ProcessDeletions::PrintSummary(){
	cout << "In print summary\n";
	int msize = data_set->num_loci();
	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

/*
 * Function: process
 * Description:
 * Main lead into processing, diverts to perform_evaliation
 */
void ProcessDeletions::process(DataSet* ds){
	data_set = ds;

	Deletions dels(data_set);
	dels.setOptions(&options);
	dels.setOrder(this->order);
	dels.setOverwrite(this->overwrite);

	dels.calculate();
	cout << "AFter calculate\n";

}


void ProcessDeletions::filter_markers(){
	return;
}

void ProcessDeletions::filter(){
	return;
}

