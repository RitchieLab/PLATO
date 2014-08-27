#ifndef PROCESSLOGREG_H
#define PROCESSLOGREG_H

#include "Process.h"

class ProcessLogReg : public ProcessImpl<ProcessLogReg>{
private:
	static const string stepname;

public:
	ProcessLogReg(){name="Logistic Regression";}
	virtual ~ProcessLogReg(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);

private:
	void doFilter(Methods::Marker*, double);
};
#endif
