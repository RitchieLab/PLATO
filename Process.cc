#include "Process.h"
#include <Helpers.h>

using std::getString;
using Methods::opts;

void Process::FilterSummary(){
	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;
}

void Process::run(Methods::DataSet* ds){
	this->process(ds);
	this->PrintSummary();
	this->filter();
	this->FilterSummary();
}


