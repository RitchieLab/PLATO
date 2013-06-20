#ifndef PROCESSTPEDOUTPUT_H
#define PROCESSTPEDOUTPUT_H

#include "Process.h"

class ProcessTPEDOutput : public ProcessImpl<ProcessTPEDOutput>{
private:
	static const std::string stepname;

public:
	ProcessTPEDOutput(){name="Create TPED file";}
	virtual ~ProcessTPEDOutput(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
