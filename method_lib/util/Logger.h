#ifndef UTIL_LOGGER_H
#define UTIL_LOGGER_H

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include <boost/thread.hpp>

namespace PLATO{
namespace Utility{

class Logger{
private:
	Logger();
	~Logger();
	Logger(const Logger&);
	Logger& operator=(const Logger&);

public:
	static Logger& getLogger();
	static void setLogFile(const std::string& fn);

	static void log(const std::string& msg);
	template <class E>
	static void log_err(const std::string& msg, bool fatal = false);
	static void log_err(const std::string& msg, bool fatal = false){log_err<std::logic_error>(msg, fatal);}

private:
	void print(const std::string& msg, std::ostream& out);
	void resetFile(const std::string& fn);

	static std::string logfn;
	static Logger* _log;

	// We need to make the log thread-safe
	boost::mutex _log_mutex;
	std::ofstream logstream;
};

template <class E>
void Logger::log_err(const std::string& msg, bool fatal){
	getLogger().print(msg, std::cerr);
	if(fatal){
		throw E(msg);
	}
}

}
}

#endif
