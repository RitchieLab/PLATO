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

#include "ProcessConcordance.h"
#include <Concordance.h>

using Methods::Concordance;
using Methods::DataSet;

using std::string;

const string ProcessConcordance::stepname = ProcessConcordance::doRegister("concordance");

void ProcessConcordance::process(DataSet* ds){
	data_set = ds;
	Concordance con(data_set);
	con.setOrder(this->order);
	con.setOverwrite(this->overwrite);
	con.setOptions(options);
	con.calculate();

}//end method process(DataSet* ds)

