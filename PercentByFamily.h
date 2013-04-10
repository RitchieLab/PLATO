#ifndef PERCENTBYFAMILY_H
#define PERCENTBYFAMILY_H

#include "Process.h"

class PercentByFamily : public ProcessImpl<PercentByFamily>{
private:
	const static std::string stepname;

public:
	PercentByFamily() {name = "Family Genotyping Efficiency";};
	virtual ~PercentByFamily(){};

protected:
	virtual void PrintSummary();
	virtual void filter();
	virtual void process(Methods::DataSet*);

};

#endif
