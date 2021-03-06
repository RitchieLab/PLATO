/*
 * MPIProcess.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: jrw32
 */

#include "MPIProcess.h"

#include "util/Logger.h"
#include "util/MPIUtils.h"

#include "config.h"

#ifdef HAVE_CXX_MPI
#include <mpi.h>
#endif

using std::pair;

namespace PLATO{

MPIProcess::MPIProcess() : n_procs(1), _mpi_threads(1), _tag(0) {
#ifdef HAVE_CXX_MPI
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
#endif
}

void MPIProcess::sendMPI(const pair<unsigned int, const char*>& nextval){
	// Note: if we don't have MPI, you REALLY shouldn't be here!
#ifdef HAVE_CXX_MPI
	if(nextval.second != 0){
		int recv = _idle_queue.front();
		_idle_queue.pop_front();
		//std::cout << "Sending " << nextval.first << " bytes to " << recv << std::endl;
		MPI_Ssend(const_cast<char*>(nextval.second), nextval.first, MPI_CHAR, recv, _tag, MPI_COMM_WORLD);
	}
#endif
}

void MPIProcess::sendAll(unsigned int bufsz, const char* buf) const{
#ifdef HAVE_CXX_MPI
	MPI_Request req_array[n_procs - 1];
	MPI_Status stat_array[n_procs - 1];

	for(int i=1; i<n_procs; i++){
		MPI_Isend(const_cast<char*>(buf), bufsz, MPI_CHAR, i, _tag, MPI_COMM_WORLD, &req_array[i-1]);
	}
	
	MPI_Waitall(n_procs - 1, req_array, stat_array);
#endif
}

void MPIProcess::collect(){

#ifdef HAVE_CXX_MPI
	MPI_Status m_stat;
	char* buf;
	int bufsz;

	while(_idle_queue.size() < static_cast<unsigned int>(n_procs - 1)*_mpi_threads){
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &m_stat);
		MPI_Get_count(&m_stat, MPI_CHAR, &bufsz);
		buf = new char[bufsz];

		MPI_Recv(buf, bufsz, MPI_CHAR, m_stat.MPI_SOURCE, _tag, MPI_COMM_WORLD, &m_stat);
		_idle_queue.push_back(m_stat.MPI_SOURCE);
		processResponse(bufsz, buf);
		delete[] buf;
	}
#endif
}

void MPIProcess::processMPI(unsigned int threads){

	// make sure to have at least 1 thread per core!
	_mpi_threads = threads < 1 ? 1 : threads;
	// However, if we aren't threadsafe, force to 1 thread!
	if(!Utility::MPIUtils::threadsafe_mpi){
		_mpi_threads = 1;
	}

	pair<unsigned int, const char*> nextval = nextQuery();

#ifdef HAVE_CXX_MPI
	int n_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	// set up the list of processors currently idle
	for(unsigned int i=0; i < (n_procs - 1) * _mpi_threads; i++){
		// we're going to initially try to round-robin everything
		// i.e. 3 threads on 2 nodes should have the initial queue: [1,2,1,2,1,2]
		_idle_queue.push_back(1 + (i % (n_procs - 1)));
	}

	char * buf;
	int bufsz;
	MPI_Status m_stat;
	while(nextval.first != 0){
		// send all the messages I can:
		while(nextval.second != 0 && !_idle_queue.empty()){
			sendMPI(nextval);
			delete[] nextval.second;
			nextval = nextQuery();
		}

		// at this point, I'm in need of a response...
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &m_stat);
		MPI_Get_count(&m_stat, MPI_CHAR, &bufsz);
		buf = new char[bufsz];

		MPI_Recv(buf, bufsz, MPI_CHAR, m_stat.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &m_stat);
		_idle_queue.push_back(m_stat.MPI_SOURCE);

		processResponse(bufsz, buf);
		delete[] buf;

		// if I was waiting on a response, see if I have a new query to send
		if(nextval.second == 0){
			nextval = nextQuery();
		}

	}

	// If we're here, we have nothing more to process, so please wait for all
	// outstanding responses
	collect();


#else
	Utility::Logger::log_err("WARNING: MPI requested, but no MPI found during compilation!");
	pair<unsigned int, const char*> calc_val;
	std::deque<pair<unsigned int, const char*> > resp_queue;
	boost::mutex resp_mutex;

	while(nextval.first != 0){
		if(nextval.second != 0){
			MPIProcessFactory::getFactory().calculate(_tag, nextval.first, nextval.second, resp_queue, resp_mutex);
			delete[] nextval.second;
			while(resp_queue.size() > 0){
				calc_val = resp_queue.front();
				resp_queue.pop_front();
				processResponse(calc_val.first, calc_val.second);
				delete[] calc_val.second;
			}

		}
		nextval = nextQuery();
	}

#endif
}

}
