/*
 * MPIProcessFactory.cpp
 *
 *  Created on: Jul 7, 2014
 *      Author: jrw32
 */

#include "MPIProcessFactory.h"

#include <algorithm>

using std::map;
using std::pair;
using std::string;

namespace PLATO{

const string& MPIProcessFactory::RegisterMPIProcess(const string& key, calcFunc* ptr){
	calc_map[key] = ptr;
	return key;
}

pair<unsigned int, const char*> MPIProcessFactory::calculate(unsigned int key_pos, unsigned int bufsz, const char* buf){
	map<string, calcFunc*>::const_iterator fn_itr = calc_map.begin();
	unsigned int i=0;
	while(fn_itr != calc_map.end() && ++i < key_pos){
		++fn_itr;
	}
	return (fn_itr == calc_map.end()) ? pair<unsigned int, const char*>(0, 0) : (*fn_itr).second(bufsz, buf);
}

unsigned int MPIProcessFactory::getKeyPos(const string& key){
	map<const string, calcFunc*>::iterator fn_itr = calc_map.find(key);
	return (fn_itr == calc_map.end()) ? 0 : (std::distance(calc_map.begin(), fn_itr) + 1);
}

}
