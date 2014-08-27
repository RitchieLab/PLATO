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

#include "ProcessDeletions.h"
#include <iostream>
#include <Deletions.h>

using std::cout;
using std::string;
using Methods::DataSet;
using Methods::Deletions;

const string ProcessDeletions::stepname=ProcessDeletions::doRegister("deletions");

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
}




