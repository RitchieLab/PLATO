#ifndef PROCESSSUPERLINKOUTPUT_H
#define PROCESSSUPERLINKOUTPUT_H

#include "Process.h"

class ProcessSuperlinkOutput : public ProcessImpl<ProcessSuperlinkOutput>{
private:
	static const std::string stepname;

public:
	ProcessSuperlinkOutput(){name="Create Superlink input file";}
	virtual ~ProcessSuperlinkOutput(){};

	virtual void setThreshold(string s);

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
