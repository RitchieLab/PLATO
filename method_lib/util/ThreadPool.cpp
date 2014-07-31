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
	pool_mutex.lock();
	tasks.push(f);

	if(tasks.size() - available == 1 && tg.size() < max_threads){
		createThread();
		pool_mutex.unlock();

		// give the newly created thread a chance to run
		boost::this_thread::yield();

		pool_mutex.lock();
	}

	// Add to the list of available threads
	if(available){
		notifier.notify_one();
	}
	pool_mutex.unlock();
}

ThreadPool::~ThreadPool(){
	pool_mutex.lock();
	running = false;
	notifier.notify_all();
	pool_mutex.unlock();

	tg.join_all();

}

void ThreadPool::createThread(){
	// no need to lock here, this is only called in the run() method
	// and is within a lock there.
	++available;
	tg.create_thread( boost::bind( &ThreadPool::pool_main, this) );
}

void ThreadPool::join_all(){

	boost::unique_lock<boost::mutex> lock(pool_mutex);

	// juist wait until we have all threads available for processing
	while (available < tg.size()) {
		notifier.wait(lock);
	}
}


void ThreadPool::pool_main(){

//	pool_mutex.lock();
//	++available;
//	pool_mutex.unlock();

	while (running) {
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
			tasks.pop();
			--available;
			lock.unlock();

			task();

			// Task has finished, so increment count of available threads.
			lock.lock();
			++available;
			// only send a notification if there's nothing to do
			if(tasks.empty()){
				notifier.notify_all();
			}
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

