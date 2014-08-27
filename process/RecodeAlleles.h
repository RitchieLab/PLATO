#ifndef PROCESSLIB_RECODEALLELES_H   //change to actual process name
#define PROCESSLIB_RECODEALLELES_H   //change to actual process name

#include <string>

#include "Process.h"

namespace PLATO{
namespace ProcessLib{

class RecodeAlleles : public ProcessImpl<RecodeAlleles>{
private:
	const static std::string stepname;

public:
	RecodeAlleles() : ProcessImpl<RecodeAlleles>(stepname, "Recode the referent / alternate alleles in markers") {};
	virtual ~RecodeAlleles(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);


protected:
	virtual void PrintSummary();
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:
	std::string recode_fn;
	std::string out_fn;
	bool auto_recode;
	bool map_input;
	bool map3;
	bool map_output;

};

}
}

#endif
