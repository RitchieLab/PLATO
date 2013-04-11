#ifndef PROCESSIBS_H
#define PROCESSIBS_H

#include "Process.h"

class ProcessIBS : public ProcessImpl<ProcessIBS>{
private:
	static const std::string stepname;

public:
	ProcessIBS(){name="IBS";}
	virtual ~ProcessIBS(){};

	virtual void setThreshold(std::string s);

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
