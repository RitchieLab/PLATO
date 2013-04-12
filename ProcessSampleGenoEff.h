#ifndef PROCESSSAMPLEGENOEFF_H
#define PROCESSSAMPLEGENOEFF_H
#include "Process.h"

#include <vector>
#include <map>

class ProcessSampleGenoEff : public ProcessImpl<ProcessSampleGenoEff>{
private:
	static const std::string stepname;

	std::vector<int> zeros;
	std::vector<int> total;
	std::vector< std::map<std::string, int> > enzyme_zeros;
	std::vector< std::map<std::string, int> > enzyme_total;

	int orig_num_samples;

public:
	ProcessSampleGenoEff() : orig_num_samples(0) {name="Individual Sample Genotyping Efficiency";}
	virtual ~ProcessSampleGenoEff(){}

protected:

	virtual void PrintSummary();
	virtual void filter();
	virtual void process(Methods::DataSet*);
	virtual void FilterSummary();
};
#endif
