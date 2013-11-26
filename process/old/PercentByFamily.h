#ifndef PERCENTBYFAMILY_H
#define PERCENTBYFAMILY_H

#include "Process.h"

#include <Sample.h>
#include <Marker.h>
#include <Family.h>

#include <vector>
#include <map>

class PercentByFamily : public ProcessImpl<PercentByFamily>{
private:
	const static std::string stepname;

	std::vector<Methods::Sample*>* samples;
	std::vector<Methods::Marker*>* markers;
	std::vector<Methods::Family*>* families;
	std::vector<int>* marker_map;

	std::vector<int> fzeros;
	std::vector<int> ftotal;
	std::vector<int> mzeros;
	std::vector<int> mtotal;
	std::vector< std::map<std::string, int> > enzyme_zeros;
	std::vector< std::map<std::string, int> > enzyme_total;

public:
	PercentByFamily() {name = "Family Genotyping Efficiency";};
	virtual ~PercentByFamily(){};

protected:
	virtual void PrintSummary();
	virtual void filter();
	virtual void process(Methods::DataSet*);

};

#endif
