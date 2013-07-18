#include "Logger.h"

using std::string;
using std::ofstream;
using std::ostream;
using std::endl;

namespace Utility{

string Logger::logfn = "plato.log";
Logger* Logger::logger = 0;

Logger::Logger() : logstream(logfn.c_str()){}

Logger::~Logger(){
	logstream.close();
}

void Logger::setLogFile(const string& fn){
	if(fn != logfn){
		logfn = fn;
		if(logger){
			delete logger;
			logger = new Logger();
		}
	}
}

void Logger::print(const string& msg, ostream& out){
	logstream << msg << endl;

	out << msg << endl;

}

}
