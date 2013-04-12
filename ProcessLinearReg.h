#ifndef PROCESSLINEARREG_H
#define PROCESSLINEARREG_H

#include "Process.h"

class ProcessLinearReg : public ProcessImpl<ProcessLinearReg>{
private:
	static const string stepname;

public:
	ProcessLinearReg(){name="Linear Regression";}
	virtual ~ProcessLinearReg(){};

protected:
	void PrintSummary();
	void process(Methods::DataSet*);

private:
	void doFilter(Methods::Marker*, double);
};
#endif
