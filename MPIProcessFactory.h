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
#include <deque>
#include <utility>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>

typedef void (calcFunc)(unsigned int, const char*, std::deque<std::pair<unsigned int, const char*> >&, boost::mutex&);

namespace PLATO{

class MPIProcessFactory{

private:
	MPIProcessFactory(){}
	MPIProcessFactory(const MPIProcessFactory&);
	MPIProcessFactory& operator=(const MPIProcessFactory&);

public:
	const std::string& RegisterMPIProcess(const std::string& key, calcFunc* ptr);

	void calculate(unsigned int key_pos, unsigned int bufsz, const char* buf, std::deque<std::pair<unsigned int, const char*> >&, boost::mutex&);
	unsigned int getKeyPos(const std::string& key);
	boost::condition_variable& getConditionVar() { return cv; }

	static MPIProcessFactory& getFactory(){static MPIProcessFactory f; return f;}

private:
	std::map<const std::string, calcFunc*> calc_map;
	//! This condition variable is set when a calculate finishes or when a receive is ready
	boost::condition_variable cv;

};





}
#endif /* MPIPROCESSFACTORY_H */
