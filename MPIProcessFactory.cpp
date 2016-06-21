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

void MPIProcessFactory::calculate(unsigned int key_pos,
		unsigned int bufsz, const char* buf,
		std::deque<std::pair<unsigned int, const char*> >& resp_queue, boost::mutex& resp_mutex){
#ifdef HAVE_OSX
map<const string, calcFunc*>::const_iterator fn_itr = calc_map.begin();
#else
	map<string, calcFunc*>::const_iterator fn_itr = calc_map.begin();
#endif
	unsigned int i=0;
	while(fn_itr != calc_map.end() && ++i < key_pos){
		++fn_itr;
	}
	if(fn_itr != calc_map.end()){
		(*fn_itr).second(bufsz, buf, resp_queue, resp_mutex, cv);
	}
}

unsigned int MPIProcessFactory::getKeyPos(const string& key){
	map<const string, calcFunc*>::iterator fn_itr = calc_map.find(key);
	return (fn_itr == calc_map.end()) ? 0 : (std::distance(calc_map.begin(), fn_itr) + 1);
}

}
