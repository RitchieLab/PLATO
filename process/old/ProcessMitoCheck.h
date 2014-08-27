#ifndef PROCESSMITOCHECK_H
#define PROCESSMITOCHECK_H

#include "Process.h"
#include <vector>
#include <Marker.h>

class ProcessMitoCheck : public ProcessImpl<ProcessMitoCheck>{
private:
	static const std::string stepname;

	std::vector<int> merrors;
	std::vector<int> serrors;
	std::vector<std::vector<Methods::Marker*> > error_map;
	std::vector<Methods::Marker*> good_markers;

	int orig_num_samples;

public:
	ProcessMitoCheck() : orig_num_samples(0) {name="Mitochondrial Error Checking (Chrom 26)";}
	virtual ~ProcessMitoCheck(){};

protected:
	virtual void PrintSummary();
	virtual void filter();
	virtual void process(Methods::DataSet*);
	virtual void FilterSummary();
};
#endif
