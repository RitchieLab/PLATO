/*
 * MPIProcessFactory.h
 *
 *  Created on: Jul 7, 2014
 *      Author: jrw32
 */

#ifndef MPIPROCESSFACTORY_H
#define MPIPROCESSFACTORY_H

#include <utility>
#include <map>
#include <string>

typedef std::pair<unsigned int, const char*> (calcFunc)(unsigned int, const char*);

namespace PLATO{

class MPIProcessFactory{

private:
	MPIProcessFactory(){}
	MPIProcessFactory(const MPIProcessFactory&);
	MPIProcessFactory& operator=(const MPIProcessFactory&);

public:
	const std::string& RegisterMPIProcess(const std::string& key, calcFunc* ptr);

	std::pair<unsigned int, const char*> calculate(unsigned int key_pos, unsigned int bufsz, const char* buf);
	unsigned int getKeyPos(const std::string& key);

	static MPIProcessFactory& getFactory(){static MPIProcessFactory f; return f;}

private:
	std::map<const std::string, calcFunc*> calc_map;
};





}
#endif /* MPIPROCESSFACTORY_H */
