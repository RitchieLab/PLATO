#ifndef PROCESS_H
#define PROCESS_H

#include <string>

#include <StepOptions.h>
#include <Options.h>

#include "ProcessFactory.h"

class Process{

public:
	Process();
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
	Methods::DataSet* data_set;

	std::string name;
	bool overwrite;
	int order;

	int orig_num_markers;

	bool _MARKERLIST_;
	bool _STRATIFY_;
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
