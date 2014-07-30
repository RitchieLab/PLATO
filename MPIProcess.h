/*
 * MPIProcess.h
 *
 *  Created on: Jul 7, 2014
 *      Author: jrw32
 */

#ifndef MPIPROCESS_H
#define MPIPROCESS_H

#include <string>
#include <utility>
#include <deque>

#include "MPIProcessFactory.h"

namespace PLATO{

class MPIProcess {
protected:
	MPIProcess();

protected:
	virtual void processResponse(unsigned int bufsz, const char* buf) = 0;
	virtual std::pair<unsigned int, const char*> nextQuery() = 0;

protected:
	virtual const std::string& getMPIName() const = 0;
	void processMPI(unsigned int threads);
	// waits for all currently running models to return
	void collect();
	void sendAll(unsigned int, const char*) const;

private:
	void sendMPI(const std::pair<unsigned int, const char*>& query);
	std::deque<int> _idle_queue;
	int n_procs;
	//! Number of cores per node.
	unsigned int _mpi_threads;

protected:
	int _tag;
};

template <class T>
class MPIProcessImpl : virtual public MPIProcess{

protected:
	MPIProcessImpl(const std::string& name) : _mpi_name(name) {
		_tag = MPIProcessFactory::getFactory().getKeyPos(_mpi_name);
	}

protected:
	static const std::string& registerMPI(const std::string& key_in);
	virtual const std::string& getMPIName() const { return _mpi_name;}

private:
	std::string _mpi_name;
};

template<class T>
const std::string& MPIProcessImpl<T>::registerMPI(const std::string& key_in){
	return MPIProcessFactory::getFactory().RegisterMPIProcess(key_in, &T::calculate_MPI);
}


}
#endif /* MPIPROCESS_H_ */
