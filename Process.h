#ifndef PROCESS_H
#define PROCESS_H

#include <string>

#include <boost/program_options.hpp>

#include "ProcessFactory.h"

namespace Methods{
	class DataSet;
}

class Process{

public:
	Process(const std::string& name_in) : name(name_in){}
	virtual ~Process(){}

	void run(Methods::DataSet*);
	boost::program_options::options_description& addOptions(boost::program_options::options_description& opts);

	virtual void parseOptions(const boost::program_options::variables_map& vm) = 0;

	void printHelp(std::ostream& o){o << *opt_ptr;}
	const std::string& getName(){return name;}

protected:
	virtual void process(Methods::DataSet*) = 0;
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts) = 0;

	virtual void PrintSummary(){};

protected:
	Methods::DataSet* data_set;
	std::string name;

private:

	boost::program_options::options_description* opt_ptr;
};

template <class T>
class ProcessImpl : public Process {
public:
	static Process* create(){return new T(T::stepname);}

protected:
	static const std::string& doRegister(const std::string& key_in);
};

template<typename T>
const std::string& ProcessImpl<T>::doRegister(const std::string& key_in){
	return ProcessFactory::getFactory().RegisterProcess(key_in, &T::create);
}

#endif
