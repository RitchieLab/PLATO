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

#ifdef PLATOLIB
namespace PlatoLib
{
#endif
string ProcessConcordance::stepname = "concordance";
#ifdef PLATOLIB
ProcessConcordance::ProcessConcordance(string bn, int pos, Database* pdb, string projPath)
{
	name = "Concordance";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
	projectPath = projPath;
}
#endif

void ProcessConcordance::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessConcordance::PrintSummary(){}

void ProcessConcordance::filter(){}

void ProcessConcordance::process(DataSet* ds){
	data_set = ds;
	Concordance con(data_set);
	con.setOrder(this->order);
#ifdef PLATOLIB
	options.setOverrideOut(projectPath + "\\");
	con.setOverwrite(true);
#else
	con.setOverwrite(this->overwrite);
#endif
	con.setOptions(options);
	con.calculate();

#ifdef PLATOLIB
	filenames.push_back(con.get_sample_error_file());
	filenames.push_back(con.get_main_file());
	filenames.push_back(con.get_error_file());
#endif
}//end method process(DataSet* ds)

#ifdef PLATOLIB
void ProcessConcordance::create_tables(){}

void ProcessConcordance::dump2db(){}

void ProcessConcordance::resize(int i){}

void ProcessConcordance::run(DataSetObject* ds)
{
	process(ds);
}//end method run(DataSetObject* ds)

}//end namespace PlatoLib
#endif
