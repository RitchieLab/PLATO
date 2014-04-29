#ifndef PROCESSCONCORDANCE_H
#define PROCESSCONCORDANCE_H

#include "Process.h"
#include <Options.h>
#include <General.h>
#include <StepOptions.h>
#include <DataSet.h>
#include <MethodException.h>

class ProcessConcordance : public ProcessImpl<ProcessConcordance>{

private:
	static const std::string stepname;

public:
	ProcessConcordance(){name="Concordance Check";}
	virtual ~ProcessConcordance(){}

protected:
	virtual void process(Methods::DataSet*);
    };
#endif
