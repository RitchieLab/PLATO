#ifndef PROCESSBINOUTPUT_H
#define PROCESSBINOUTPUT_H

#include "Process.h"

class ProcessBINOutput: public ProcessImpl<ProcessBINOutput> {

private:
	const static std::string stepname;

public:
	ProcessBINOutput() {name="Binary (.bed, .fam, .bim) file output";}
	virtual ~ProcessBINOutput(){}

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
