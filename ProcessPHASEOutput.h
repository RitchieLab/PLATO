#ifndef PROCESSPHASEOUTPUT_H
#define PROCESSPHASEOUTPUT_H


#include "Process.h"

class ProcessPHASEOutput : public ProcessImpl<ProcessPHASEOutput>{
private:
	static const std::string stepname;

public:
	ProcessPHASEOutput(){name="Create Phase input file";}
	virtual ~ProcessPHASEOutput(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
