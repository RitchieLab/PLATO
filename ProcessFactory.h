#ifndef PROCESS_FACTORY_H
#define PROCESS_FACTORY_H

#include <map>
#include <string>

namespace PLATO{

class Process;

typedef Process* (createFunc)();

class ProcessFactory{
public:
	typedef std::map<const std::string, createFunc*>::const_iterator const_iterator;

private:
	ProcessFactory(){}
	ProcessFactory(const ProcessFactory&);
	ProcessFactory& operator=(const ProcessFactory&);

public:
	const std::string& RegisterProcess(const std::string& key, createFunc* ptr);
	Process* Create(const std::string& key);

	const_iterator begin() const{return creation_map.begin();}
	const_iterator end() const{return creation_map.end();}

	static ProcessFactory& getFactory(){static ProcessFactory f; return f;}

private:
	std::map<const std::string, createFunc*> creation_map;
};

}

#endif
