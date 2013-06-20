#ifndef PROCESSFST_H
#define PROCESSFST_H

#include "Process.h"
#include <cdflib.h>

class ProcessFst : public ProcessImpl<ProcessFst>{
private:
	static const std::string stepname;

public:
	ProcessFst(){name="FST";}
	virtual ~ProcessFst(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);

private:
	void doFilter(Methods::Marker*, double);
};
#endif
