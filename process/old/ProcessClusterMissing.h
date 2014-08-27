#ifndef PROCESSCLUSTERMISSING_H
#define PROCESSCLUSTERMISSING_H

#include "Process.h"
#include <string>

class ProcessClusterMissing : public ProcessImpl<ProcessClusterMissing>{

private:
	static const std::string stepname;

public:
	ProcessClusterMissing(){name="Cluster Missing";}
	virtual ~ProcessClusterMissing(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
