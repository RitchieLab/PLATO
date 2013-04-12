#ifndef PROCESSMARKERGENOEFF_H
#define PROCESSMARKERGENOEFF_H

#include "Process.h"
#include <vector>
#include <map>

class ProcessMarkerGenoEff : public ProcessImpl<ProcessMarkerGenoEff>{
private:
	static const std::string stepname;

	std::vector<int> zeros;
	std::vector<int> total;
	std::vector<int> casezeros;
	std::vector<int> casetotal;
	std::vector<int> controlzeros;
	std::vector<int> controltotal;
	std::map<std::string, std::vector<int> > groupzeros;
	std::map<std::string, std::vector<int> > grouptotal;

	std::vector<Methods::Marker*> good_markers;

public:
	ProcessMarkerGenoEff(){name="Marker Genotyping Efficiency";}
	virtual ~ProcessMarkerGenoEff(){};

protected:
	virtual void PrintSummary();
	virtual void filter();
	virtual void process(Methods::DataSet*);};
#endif
