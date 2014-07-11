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
#incldue <mpi.h>
#endif

using std::pair;

namespace PLATO{

void MPIProcess::sendMPI(const pair<unsigned int, const char*>& nextval){
	// Note: if we don't have MPI, you REALLY shouldn't be here!
#ifdef HAVE_CXX_MPI
	if(nextval.second == 0){
		MPI_Send(nextval.second, nextval.first, MPI_CHAR, _idle_queue.pop_front(), tag, MPI_COMM_WORLD);
	}
#endif
}

void MPIProcess::processMPI(){

	pair<unsigned int, const char*> nextval = nextQuery();
	unsigned int tag = MPIProcessFactory::getFactory().getKeyPos(getMPIName());

#ifdef HAVE_CXX_MPI
	int n_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

	int i=0;
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
		MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &m_stat);
		MPI_Get_count(&m_stat, MPI_CHAR, &bufsz);
		buf = new char[bufsz];
		MPI_Recv(buf, bufsz, MPI_CHAR, m_stat.MPI_SOURCE, tag, MPI_COMM_WORLD, &m_stat);
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
	for(int j=1; j<n_procs - _idle_queue.size(); j++){
		MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &m_stat);
		MPI_Get_count(&m_stat, MPI_CHAR, &bufsz);
		buf = new char[bufsz];
		MPI_Recv(buf, bufsz, MPI_CHAR, m_stat.MPI_SOURCE, tag, MPI_COMM_WORLD, &m_stat);

		processResponse(bufsz, buf);
		delete[] buf;
	}

	// Now we need to send the special "please die now" message
	// (in our case, this is a 0-length message with tag 0
	for(int j=1; j<n_procs; j++){
		MPI_Send(0, 0, MPI_CHAR, j, 0, MPI_COMM_WORLD);
	}

#else
	Utility::Logger::log_err("WARNING: MPI requested, but no MPI found during compilation!");
	pair<unsigned int, const char*> calc_val;

	while(nextval.first != 0){
		if(nextval.second != 0){
			calc_val = MPIProcessFactory::getFactory().calculate(tag, nextval.first, nextval.second);
			delete[] nextval.second;
			processResponse(calc_val.first, calc_val.second);
			delete[] calc_val.second;
		}
		nextval = nextQuery();
	}

#endif
}

}
