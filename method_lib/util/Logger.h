#ifndef UTIL_LOGGER_H
#define UTIL_LOGGER_H

#include <string>
#include <iostream>

namespace Utility{

class Logger{
private:
	Logger();
	~Logger();
	Logger(const Logger&);
	Logger& operator=(const Logger&);

public:
	static Logger& getLogger(){return (logger == 0) ? logger = new Logger() : *logger;}
	static void setLogFile(const std::string& fn);

	void log(const std::string& msg){print(msg, std::cout);}
	void log_err(const std::string& msg){print(msg, std::cerr);}

private:
	void print(const std::string& msg, std::ostream& out);

	static std::string logfn;
	static Logger* logger;

	std::ofstream logstream;
};
}
