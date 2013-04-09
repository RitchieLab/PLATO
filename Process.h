#ifndef PROCESS_H
#define PROCESS_H

#include <vector>
#include <string>
#include <map>

#include <Sample.h>
#include <Family.h>
#include <Marker.h>
#include <DataSet.h>
#include <StepOptions.h>

#include "ProcessFactory.h"

using namespace Methods;

class Process{

	public:
		Process(){options.setCovarMissing(opts::_COVAR_MISSING_); options.setTraitMissing(opts::_TRAIT_MISSING_);}
		virtual ~Process(){}

		virtual void process(DataSet*) = 0;
		virtual void PrintSummary() = 0;
		virtual void filter() = 0;
		virtual void setThreshold(std::string) = 0;
		virtual void FilterSummary() = 0;
		virtual void setRank(int) = 0;
		virtual void setDBOUT() = 0;
		virtual void setMarkerList() = 0;
		virtual void setStratify() = 0;
		virtual void setOrder(int) = 0;
		virtual void setOverwrite(bool) = 0;
		virtual bool hasIncExc() = 0;

		virtual const std::string& getName(){return name;}

		StepOptions* getOptions(){return &options;}
		StepOptions get_options(){return options;}
		void set_options(StepOptions* opts){options = *opts;}
		bool has_results(){return hasresults;}
		std::vector<std::string> get_tablename(){return tablename;}
		std::map<std::string, std::vector<std::string> > get_headers(){return headers;}
		std::vector<std::string> get_headers(std::string table){return headers[table];}
		std::vector<std::string> get_primary_table(std::string table){return primary_table[table];}
		std::string get_table_nickname(int i){return tablenicknames[i];}
		std::string get_name(){return name;}
		int get_position(){return position;}
		std::vector<std::string> get_filenames(){return filenames;}

	protected:
		StepOptions options;

		std::string name;
		std::string batchname;
		int position;
		bool hasresults;
		std::vector<std::string> tablename;
		std::map<std::string, std::vector<std::string> > headers;
		std::map<std::string, std::vector<std::string> > primary_table;
		std::vector<std::string> tablenicknames;
		std::vector<std::string> filenames;
};

template <class T>
class ProcessImpl : public Process {
public:
	static Process* create(){return new T();}

protected:
	static const std::string& doRegister(const std::string& key_in);
};

template<typename T>
const std::string& ProcessImpl<T>::doRegister(const std::string& key_in){
	return ProcessFactory::getFactory().RegisterProcess(key_in, &T::create);
}

#endif
