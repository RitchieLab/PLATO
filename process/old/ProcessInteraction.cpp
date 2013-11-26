#include "ProcessInteraction.h"
#include <Interactions.h>

using std::string;
using Methods::DataSet;
using Methods::Interactions;

const string ProcessInteraction::stepname = ProcessInteraction::doRegister("interaction");

void ProcessInteraction::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}
}

void ProcessInteraction::process(DataSet* ds){

	data_set = ds;

	Interactions inter;
	inter.setOrder(this->order);
	inter.setOverwrite(this->overwrite);
	inter.setOptions(options);
	inter.calculate(data_set);
}

