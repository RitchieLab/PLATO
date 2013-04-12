#ifndef PROCESSPOWERMARKEROUTPUT_H
#define PROCESSPOWERMARKEROUTPUT_H

#include "Process.h"

class ProcessPowerMarkerOutput : public ProcessImpl<ProcessPowerMarkerOutput>{
private:
	static const string stepname;

public:
	ProcessPowerMarkerOutput(){name="Create PowerMarker input file";}
	virtual ~ProcessPowerMarkerOutput(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
