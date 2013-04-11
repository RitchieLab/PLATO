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
#include "ProcessInteraction.h"
#include <General.h>
using namespace Methods;


string ProcessInteraction::stepname = ProcessInteraction::doRegister("interaction");

void ProcessInteraction::FilterSummary(){
}

void ProcessInteraction::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessInteraction::filter(){}

void ProcessInteraction::process(DataSet* ds){

	data_set = ds;

	Interactions inter;
	inter.setOrder(this->order);
	inter.setOverwrite(this->overwrite);
	inter.setOptions(options);
	inter.calculate(data_set);
}

