#include "Logger.h"

using std::string;
using std::ofstream;
using std::ostream;
using std::endl;

namespace PLATO{
namespace Utility{

string Logger::logfn = "plato.log";
Logger* Logger::_log = 0;

Logger::Logger() : logstream(logfn.c_str()){}

Logger::~Logger(){
	logstream.close();
}

void Logger::setLogFile(const string& fn){
	if(_log){
		_log->resetFile(fn);
	} else {
		logfn = fn;
	}
}

void Logger::resetFile(const string& fn){
	_log_mutex.lock();
	if(fn != logfn){
		logfn = fn;
		logstream.close();
		logstream.open(logfn.c_str());
	}
	_log_mutex.unlock();
}

void Logger::print(const string& msg, ostream& out){
	_log_mutex.lock();
	logstream << msg << endl;
	out << msg << endl;
	_log_mutex.unlock();

}

Logger& Logger::getLogger(){return *((_log == 0) ? _log = new Logger() : _log);}

void Logger::log(const std::string& msg){getLogger().print(msg, std::cout);}

}
}
