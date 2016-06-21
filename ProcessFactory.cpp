#include "ProcessFactory.h"
#include "Process.h"

using std::string;
using std::map;

namespace PLATO{

const string& ProcessFactory::RegisterProcess(const string& key, createFunc* ptr){
	creation_map[key] = ptr;
	return key;
}

Process* ProcessFactory::Create(const string& key){
#ifdef HAVE_OSX
	map<const std::string, createFunc*>::iterator it=creation_map.find(key);
#else
	map<std::string, createFunc*>::const_iterator it=creation_map.find(key);
#endif
	if(it != creation_map.end()){
		return (*it).second();
	}else{
		return NULL;
	}
}

}
