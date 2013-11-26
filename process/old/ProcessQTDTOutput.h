#ifndef PROCESSQTDTOUTPUT_H
#define PROCESSQTDTOUTPUT_H

#include "Process.h"

class ProcessQTDTOutput : public ProcessImpl<ProcessQTDTOutput>{

private:
	static const std::string stepname;

public:
	ProcessQTDTOutput(){name="Create QTDT input file";}
	virtual ~ProcessQTDTOutput(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
