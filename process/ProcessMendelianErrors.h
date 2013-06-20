#ifndef PROCESSMENDELIANERRORS_H
#define PROCESSMENDELIANERRORS_H


#include "Process.h"
#include <vector>
#include <map>

class ProcessMendelianErrors : public ProcessImpl<ProcessMendelianErrors>{
private:
	static const std::string stepname;

private:

	int orig_num_families;

	std::vector<int> ferrors;
	std::vector<std::map<std::string, int> > fenzyme;
	std::vector<int> merrors;
	std::vector<int> serrors;
	std::vector<std::map<std::string, int> > senzyme;
	std::vector<std::vector<Methods::Marker*> > error_map;
	std::vector<Methods::Marker*> good_markers;

public:
	ProcessMendelianErrors() : orig_num_families(0) {name="Mendelian Errors";}
	virtual ~ProcessMendelianErrors(){};

protected:
	virtual void FilterSummary();
	virtual void PrintSummary();
	virtual void filter();
	virtual void process(Methods::DataSet*);

private:
	void filter_markers();
	void zeroErrors();
};
#endif
