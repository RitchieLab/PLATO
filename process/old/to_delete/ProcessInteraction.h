#ifndef PROCESSINTERACTION_H
#define PROCESSINTERACTION_H

#include "Process.h"

class ProcessInteraction : public ProcessImpl<ProcessInteraction>{
private:
	static const std::string stepname;

public:
	ProcessInteraction(){name="Interaction";}
	virtual ~ProcessInteraction(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
