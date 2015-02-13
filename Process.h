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
	Process() {}
	virtual ~Process(){}

	void run(PLATO::Data::DataSet&);
	boost::program_options::options_description& addOptions(boost::program_options::options_description& opts);

	virtual void parseOptions(const boost::program_options::variables_map& vm) = 0;

	void printHelp(std::ostream& o){if(opt_ptr){o << *opt_ptr;}}
	virtual const std::string& getName() const = 0;
	virtual const std::string& getDesc() const = 0;

protected:
	virtual void process(PLATO::Data::DataSet&) = 0;
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts) = 0;

	virtual void PrintSummary(){};

private:

	boost::program_options::options_description* opt_ptr;
};

template <class T>
class ProcessImpl : public virtual Process {
public:
	ProcessImpl(const std::string& n, const std::string& d) : Process(), _name(n), _desc(d) {}

public:
	static Process* create(){return new T();}
	virtual const std::string& getName() const {return _name;}
	virtual const std::string& getDesc() const {return _desc;}

protected:
	static const std::string& doRegister(const std::string& key_in);

protected:
	std::string _name;
	std::string _desc;
};

template<typename T>
const std::string& ProcessImpl<T>::doRegister(const std::string& key_in){
	return ProcessFactory::getFactory().RegisterProcess(key_in, &T::create);
}

}

#endif
