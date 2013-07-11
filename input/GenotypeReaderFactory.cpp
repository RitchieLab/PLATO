#include "GenotypeReaderFactory.h"
#include "GenotypeReader.h"

#include <boost/unordered_set.hpp>

namespace po=boost::program_options;

using std::string;
using std::map;
using boost::unordered_set;
using po::value;

void GenotypeReaderFactory::RegisterReader(const string& key, createGenoReader* ptr){
	creation_map[key] = ptr;
}

GenotypeReader* GenotypeReaderFactory::Create(const po::variables_map& vm){
	createGenoReader* toCreate = 0;

	for(const_iterator itr=creation_map.begin(); itr != creation_map.end(); itr++){
		if(vm.count((*itr).first)){
			// Uh-oh! we have a problem! If this is true, we're trying to load
			// from two different genotype sources!
			if(toCreate != 0 && toCreate != (*itr).second){
				toCreate = 0;
				// TODO: probably want to throw an exception here
				break;
			}else if(toCreate == 0){
				toCreate = (*itr).second;
			}
		}
	}

	if(toCreate){
		return (*toCreate)(vm);
	}else{
		return NULL;
	}
}

po::options_description& GenotypeReaderFactory::addOptions(po::options_description& opts){

	unordered_set<createGenoReader*> uniq_ctors;
	for(const_iterator itr=creation_map.begin(); itr != creation_map.end(); itr++){
		uniq_ctors.insert((*itr).second);
	}

	po::variables_map empty_vm;
	for(unordered_set<createGenoReader*> c_itr = uniq_ctors.begin(); c_itr != uniq_ctors.end(); c_itr++){
		GenotypeReader* gr = (*c_itr)(empty_vm);
		opts = gr->addOptions(opts);
		delete gr;
	}

	opts = GenotypeReader::addOptions(opts);


	return opts;
}
