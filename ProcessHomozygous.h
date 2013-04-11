#ifndef PROCESSHOMOZYGOUS_H
#define PROCESSHOMOZYGOUS_H

#include "Process.h"
#include <vector>

class ProcessHomozygous : public ProcessImpl<ProcessHomozygous>{
private:
	static const std::string stepname;

	std::vector<int> homoaffcount;
	std::vector<int> homounaffcount;
	std::vector<int> homoallcount;
	std::vector<int> homominallcount;
	std::vector<int> homomajallcount;

	std::vector<Methods::Marker*> good_markers;


public:
	ProcessHomozygous(){name="Homozygous spans";}
	virtual ~ProcessHomozygous(){};

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
};
#endif
