#ifndef PROCESSDELETIONS_H
#define PROCESSDELETIONS_H

#include "Process.h"

class ProcessDeletions : public ProcessImpl<ProcessDeletions>{
private:
	static const std::string stepname;

public:
	ProcessDeletions(){name="Deletion Detection";}
	virtual ~ProcessDeletions(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
