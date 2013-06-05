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
#include "config.h"
#ifdef HAVE_MALLOC_H
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
#include <Helpers.h>

using namespace std;
using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

string ProcessDeletions::stepname="deletions";
#ifdef PLATOLIB
ProcessDeletions::ProcessDeletions(string bn, int pos, Database* pdb, string projPath)
{
	name = "Deletion";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
	projectPath = projPath;
}
#endif

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
#ifdef PLATOLIB
	options.setOverrideOut(projectPath + "\\");
	dels.setOverwrite(true);
#else
	dels.setOverwrite(this->overwrite);
#endif

	dels.calculate();

#ifdef PLATOLIB
	filenames.push_back(dels.get_default_filename());
	filenames.push_back(dels.get_density_filename());
#endif

}


void ProcessDeletions::filter_markers(){
	return;
}

void ProcessDeletions::filter(){
	return;
}


#ifdef PLATOLIB
void ProcessDeletions::create_tables(){}

void ProcessDeletions::dump2db(){}

void ProcessDeletions::resize(int i){}

void ProcessDeletions::run(DataSetObject* ds)
{
	process(ds);
}
#endif

#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
