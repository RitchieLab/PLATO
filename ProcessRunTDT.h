#ifndef PROCESSRUNTDT_H
#define PROCESSRUNTDT_H

#include "Process.h"
#include <vector>
#include <Marker.h>

class ProcessRunTDT : public ProcessImpl<ProcessRunTDT>{

private:
	static const std::string stepname;
	std::vector<double> chi;
	std::vector<long double> pval;
	std::vector<int> fams_used;
	std::vector<float> maf;
	std::vector<double> trans;
	std::vector<double> untrans;

	//by group variables
	std::vector<std::vector<double> > gchi;
	std::vector<std::vector<long double> > gpval;
	std::vector<std::vector<int> > gfams_used;
	std::vector<std::vector<float> > gmaf;
	std::vector<std::vector<double> > gtrans;
	std::vector<std::vector<double> > guntrans;

	std::vector<Methods::Marker*> good_markers;

public:
	ProcessRunTDT(){name="TDT";};
	virtual ~ProcessRunTDT(){};


protected:
	virtual void PrintSummary();
	virtual void filter();
	virtual void process(Methods::DataSet*);
};
#endif
