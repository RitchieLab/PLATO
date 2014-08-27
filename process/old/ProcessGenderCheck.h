#ifndef PROCESSGENDERCHECK_H
#define PROCESSGENDERCHECK_H

#include "Process.h"
#include <vector>
#include <map>

class ProcessGenderCheck : public ProcessImpl<ProcessGenderCheck>{
private:
	static const std::string stepname;

	int orig_num_samples;

	std::vector<int> merrors;
	std::vector<int> shets;
	std::vector<int> mtotal;
	std::vector<int> stotal;
	std::vector<std::map<std::string, int> > senzyme_hets;
	std::vector<std::map<std::string, int> > senzyme_tot;

	std::vector<Methods::Marker*> good_markers;

public:
	ProcessGenderCheck() : orig_num_samples(0) {name="Gender Correctness (using X-chromosome markers)";}
	virtual ~ProcessGenderCheck(){};

protected:
	virtual void PrintSummary();
	virtual void filter();
	virtual void process(Methods::DataSet*);

};
#endif
