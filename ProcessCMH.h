#ifndef PROCESSCMH_H
#define PROCESSCMH_H

#include "Process.h"
#include <Options.h>
#include <General.h>
#include <StepOptions.h>
#include <DataSet.h>
#include <CMH.h>

using namespace std;
using namespace Methods;


class ProcessCMH : public ProcessImpl<ProcessCMH>{

private:
	static const string stepname;

	int run_chr;
	int run_start;
	int run_end;
	string defaultinsert;


public:
	ProcessCMH(){name="Cochran-Mantel-Haenszel test";}
	virtual ~ProcessCMH(){};

protected:
	virtual void PrintSummary();
	virtual void process(DataSet*);

private:
	void doFilter(Methods::Marker*, double);
};
#endif
