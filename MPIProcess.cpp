/*
 * MPIProcess.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: jrw32
 */

#include "MPIProcess.h"
#include "util/Logger.h"

#include "config.h"

#ifdef HAVE_CXX_MPI
#include <mpi.h>
#endif

using std::pair;

namespace PLATO{

MPIProcess::MPIProcess() :_tag(0), n_procs(1) {
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
		MPI_Send(const_cast<char*>(nextval.second), nextval.first, MPI_CHAR, recv, _tag, MPI_COMM_WORLD);
	}
#endif
}

void MPIProcess::sendAll(unsigned int bufsz, const char* buf) const{
#ifdef HAVE_CXX_MPI
	for(unsigned int i=1; i<n_procs; i++){
		MPI_Request req;
		MPI_Isend(const_cast<char*>(buf), bufsz, MPI_CHAR, i, _tag, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);
	}
#endif
}

void MPIProcess::collect(){

#ifdef HAVE_CXX_MPI
	MPI_Status m_stat;
	char* buf;
	int bufsz;

	for(unsigned int j=1; j<n_procs - _idle_queue.size(); j++){
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &m_stat);
		MPI_Get_count(&m_stat, MPI_CHAR, &bufsz);
		buf = new char[bufsz];

		MPI_Recv(buf, bufsz, MPI_CHAR, m_stat.MPI_SOURCE, _tag, MPI_COMM_WORLD, &m_stat);
		processResponse(bufsz, buf);
		delete[] buf;
	}
#endif
}

void MPIProcess::processMPI(){

	pair<unsigned int, const char*> nextval = nextQuery();

#ifdef HAVE_CXX_MPI
	int n_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	// set up the list of processors currently idle
	for(int i=0; ++i < n_procs && nextval.first != 0; ){
		_idle_queue.push_back(i);
	}

	// at this point, everyone is working, so we'll wait for a response to send
	// another message
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

	while(nextval.first != 0){
		if(nextval.second != 0){
			calc_val = MPIProcessFactory::getFactory().calculate(_tag, nextval.first, nextval.second);
			delete[] nextval.second;
			processResponse(calc_val.first, calc_val.second);
			delete[] calc_val.second;
		}
		nextval = nextQuery();
	}

#endif
}

}
