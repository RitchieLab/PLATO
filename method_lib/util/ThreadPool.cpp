/*
 * ThreadPool.cpp
 *
 *  Created on: Jul 29, 2014
 *      Author: jrw32
 */

#include "ThreadPool.h"

#include <boost/bind.hpp>

namespace PLATO{
namespace Utility{

void ThreadPool::run(boost::function<void()>& f){
	// if everything is running, create a new thread
	if(available == 0 && tg.size() < max_threads){
		createThread();
	}

	// Add to the list of available threads
	pool_mutex.lock();
	tasks.push_back(f);
	if(available){
		notifier.notify_one();
	}
	pool_mutex.unlock();
}

ThreadPool::~ThreadPool(){
	pool_mutex.lock();
	running = false;
	pool_mutex.unlock();

	tg.join_all();

}

void ThreadPool::createThread(){
	pool_mutex.lock();
	tg.create_thread( boost::bind( &ThreadPool::pool_main, this) );
	pool_mutex.unlock();
}

void ThreadPool::join_all(){

	boost::unique_lock<boost::mutex> lock(pool_mutex);

	// juist wait until we have all threads available for processing
	while (available < tg.size()) {
		notifier.wait(lock);
	}
}


void ThreadPool::pool_main(){

	pool_mutex.lock();
	++available;
	pool_mutex.unlock();

	while (running && !(exiting || tasks.empty())) {
		// Wait on condition variable while the task is empty and the pool is
		// still running.
		boost::unique_lock<boost::mutex> lock(pool_mutex);
		while (tasks.empty() && running) {
			notifier.wait(lock);
		}
		// If pool is no longer running, break out.
		if (running){
			// Copy task locally and remove from the queue.  This is done within
			// its own scope so that the task object is destructed immediately
			// after running the task.  This is useful in the event that the
			// function contains shared_ptr arguments bound via bind.
			boost::function<void()> task = tasks.front();
			tasks.pop_front();
			--available;
			lock.unlock();

			task();

			// Task has finished, so increment count of available threads.
			lock.lock();
			++available;
			lock.unlock();
		}
	} // while running

	// thread is exiting, so task is unavailable
	pool_mutex.lock();
	--available;
	pool_mutex.unlock();
}

}
}

