#ifndef PROCESSKINSHIP_H
#define PROCESSKINSHIP_H

#include "Process.h"

class ProcessKinship : public ProcessImpl<ProcessKinship>{

private:
	static const std::string stepname;

public:
	ProcessKinship(){name="Kinship";}
	virtual ~ProcessKinship(){}
protected:
	void PrintSummary();
	void process(Methods::DataSet*);
};
#endif
