#ifndef PROCESS_H
#define PROCESS_H

#include <vector>
#include <string>
#include <map>

#include <Sample.h>
#include <Family.h>
#include <Marker.h>
#include <DataSet.h>
#include <StepOptions.h>
#include <Options.h>

#include "ProcessFactory.h"

class Process{

public:
	Process(){options.setCovarMissing(Methods::opts::_COVAR_MISSING_); options.setTraitMissing(Methods::opts::_TRAIT_MISSING_);}
	virtual ~Process(){}

	void run(Methods::DataSet*);

	bool hasIncExc(){return options.doIncExcludedSamples();}
	void setOrder(int o){order = o;}
	void setOverwrite(bool b){overwrite = b;}
	virtual void setThreshold(std::string s){options.setUp(s);}

	const std::string& getName(){return name;}
	void setMarkerList(){_MARKERLIST_ = true;};
	void setStratify(){_STRATIFY_ = true;};

	Methods::StepOptions* getOptions(){return &options;}
	Methods::StepOptions get_options(){return options;}
	void set_options(Methods::StepOptions* opts){options = *opts;}

protected:
	virtual void process(Methods::DataSet*) = 0;
	virtual void PrintSummary(){};
	virtual void filter(){};

	virtual void FilterSummary();

protected:
	Methods::StepOptions options;

	std::string name;
	int position;
	float threshold;
	int orig_num_markers;
	int orig_num_families;
	int orig_num_samples;
	bool overwrite;
	int order;

	bool _MARKERLIST_;
	bool _STRATIFY_;

	Methods::DataSet* data_set;
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
};

template <class T>
class ProcessImpl : public Process {
public:
	static Process* create(){return new T();}

protected:
	static const std::string& doRegister(const std::string& key_in);
};

template<typename T>
const std::string& ProcessImpl<T>::doRegister(const std::string& key_in){
	return ProcessFactory::getFactory().RegisterProcess(key_in, &T::create);
}

#endif
