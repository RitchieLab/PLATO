#include "Logger.h"

using std::string;
using std::ofstream;
using std::ostream;
using std::endl;

namespace Utility{

string Logger::logfn = "plato.log";
Logger* Logger::_log = 0;

Logger::Logger() : logstream(logfn.c_str()){}

Logger::~Logger(){
	logstream.close();
}

void Logger::setLogFile(const string& fn){
	if(fn != logfn){
		logfn = fn;
		if(_log){
			delete _log;
			_log = new Logger();
		}
	}
}

void Logger::print(const string& msg, ostream& out){
	logstream << msg << endl;

	out << msg << endl;

}

Logger& Logger::getLogger(){return *((_log == 0) ? _log = new Logger() : _log);}

void Logger::log(const std::string& msg){getLogger().print(msg, std::cout);}

}
