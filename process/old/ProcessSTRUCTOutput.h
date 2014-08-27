#ifndef PROCESSSTRUCTOUTPUT_H
#define PROCESSSTRUCTOUTPUT_H

#include "Process.h"

class ProcessSTRUCTOutput : public ProcessImpl<ProcessSTRUCTOutput>{
private:
	static const std::string stepname;

public:
	ProcessSTRUCTOutput(){name="Create Structure input file";}
	virtual ~ProcessSTRUCTOutput(){};

	virtual void setThreshold(string s);

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
