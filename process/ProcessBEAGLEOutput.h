#ifndef PROCESSBEAGLEOUTPUT_H
#define PROCESSBEAGLEOUTPUT_H


#include "Process.h"
#include <Options.h>
#include <StepOptions.h>
#include <DataSet.h>
#include <MethodException.h>
#include <BEAGLEOutput.h>

//define the PlatoLib namespace for use with Plato as a library
class ProcessBEAGLEOutput : public ProcessImpl<ProcessBEAGLEOutput>{
private:
	static const std::string stepname;

public:
	ProcessBEAGLEOutput(){name="Create BEAGLE input file";}
	virtual ~ProcessBEAGLEOutput(){};

protected:
	virtual void process(Methods::DataSet*);
	virtual void PrintSummary();
};
#endif
