#ifndef PROCESSLD_H
#define PROCESSLD_H


#include "Process.h"

class ProcessLD : public ProcessImpl<ProcessLD>{

private:
	static const std::string stepname;

public:
	ProcessLD(){name="LD calculations";}
	virtual ~ProcessLD(){}

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
