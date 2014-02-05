#ifndef PROCESS_H
#define PROCESS_H

#include <string>

#include <boost/program_options.hpp>

#include "ProcessFactory.h"

namespace PLATO{

namespace Data{
	class DataSet;
}

class Process{

public:
	Process(const std::string& name, const std::string& desc) : _name(name), _desc(desc){}
	virtual ~Process(){}

	void run(PLATO::Data::DataSet&);
	boost::program_options::options_description& addOptions(boost::program_options::options_description& opts);

	virtual void parseOptions(const boost::program_options::variables_map& vm) = 0;

	void printHelp(std::ostream& o){if(opt_ptr){o << *opt_ptr;}}
	const std::string& getName() const {return _name;}
	const std::string& getDesc() const {return _desc;}

protected:
	virtual void process(PLATO::Data::DataSet&) = 0;
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts) = 0;

	virtual void PrintSummary(){};

protected:
	//PLATO::Data::DataSet* data_set;
	std::string _name;
	std::string _desc;

private:

	boost::program_options::options_description* opt_ptr;
};

template <class T>
class ProcessImpl : public Process {
public:
	ProcessImpl(const std::string& n, const std::string& d) : Process(n, d) {}

public:
	static Process* create(){return new T();}

protected:
	static const std::string& doRegister(const std::string& key_in);
};

template<typename T>
const std::string& ProcessImpl<T>::doRegister(const std::string& key_in){
	return ProcessFactory::getFactory().RegisterProcess(key_in, &T::create);
}

}

#endif
