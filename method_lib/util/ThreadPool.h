/*
 * ThreadPool.h
 *
 *  Created on: Jul 29, 2014
 *      Author: jrw32
 */

#ifndef UTILITY_THREADPOOL_H
#define UTILITY_THREADPOOL_H

#include <limits>
#include <deque>
#include <boost/thread.hpp>
#include <boost/function.hpp>

namespace PLATO{
namespace Utility{


class ThreadPool{
public:
	ThreadPool(unsigned int max=std::numeric_limits<unsigned int>::max())
		: max_threads(max < 1 ? 1 : max), available(0), running(true), exiting(false) {}
	~ThreadPool();

	/*!
	 * Wait for all threads in the pool to stop (stopping the ASIO service as well)
	 */
	void join_all();

	/*!
	 * Run a given function inside a thread
	 */
	void run(boost::function<void ()>& f);


private:

	void pool_main();

	/*!
	 * Create and start a new thread
	 */
	void createThread();

	unsigned int max_threads;

	boost::mutex pool_mutex;
	unsigned int available;
	std::deque<boost::function<void()> > tasks;
	boost::condition_variable notifier;
	boost::thread_group tg;
	bool running;
	bool exiting;


};

}
}


#endif /* THREADPOOL_H_ */
