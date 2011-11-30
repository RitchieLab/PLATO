/**********************************************************************************
*                       Superlink Input Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates Superlink input files
*
*
*
*File: SuperlinkOutput.cc
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
#include "ProcessSuperlinkOutput.h"
#include <General.h>
#include <Helper.h>

string ProcessSuperlinkOutput::stepname = "output-superlink";

void ProcessSuperlinkOutput::FilterSummary(){
}

void ProcessSuperlinkOutput::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}

}

void ProcessSuperlinkOutput::filter(){
}

void ProcessSuperlinkOutput::process(DataSet* ds){
	data_set = ds;

	SuperlinkOutput so;
	so.setOrder(this->order);
	so.setOverwrite(this->overwrite);
	so.setOptions(options);
	so.calculate(data_set);

}

