/**********************************************************************************
*                       ProcessBEAGLE Output Module
*
* Written by: Justin Giles
*             Vanderbilt University
*             Center for Human Genetics Research
*
* Outputs ProcessBEAGLE input files.
*
*
*File: ProcessBEAGLEOutput.cc
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
#include "ProcessBEAGLEOutput.h"
#include <General.h>
#include <Helper.h>

using namespace Methods;
/*
 *Function: FilterSummary
 *Description:
 *Not used.
 */
void ProcessBEAGLEOutput::FilterSummary(){
}

/*
 * Function: PrintSummary
 * Description:
 * Resets marker flags
 */
void ProcessBEAGLEOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}

}

/*
 * Function: filter
 * Description:
 * not used
 */
void ProcessBEAGLEOutput::filter(){
}

/*
 * Function: process
 * Description:
 * Main function for producing output files.
 * Creates output files, one per chromosome.
 */
void ProcessBEAGLEOutput::process(DataSet* ds){
	data_set = ds;

	BEAGLEOutput bo;
	bo.setOrder(this->order);
	bo.setOverwrite(this->overwrite);
	bo.setOptions(options);

	bo.calculate(data_set);

}


