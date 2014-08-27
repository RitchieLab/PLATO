#ifndef PROCESSFBATOUTPUT_H
#define PROCESSFBATOUTPUT_H

#include "Process.h"

class ProcessFBATOutput : public ProcessImpl<ProcessFBATOutput>{
private:
	static const std::string stepname;

public:
	ProcessFBATOutput(){name="Create FBAT input file";}
	virtual ~ProcessFBATOutput(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
