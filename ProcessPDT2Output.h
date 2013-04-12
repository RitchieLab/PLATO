#ifndef PROCESSPDT2OUTPUT_H
#define PROCESSPDT2OUTPUT_H


#include "Process.h"

class ProcessPDT2Output : public ProcessImpl<ProcessPDT2Output>{
private:
	static const std::string stepname;

public:
	ProcessPDT2Output(){name="Create PDT2 input file";}
	virtual ~ProcessPDT2Output(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
