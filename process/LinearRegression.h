#ifndef PROCESSLIB_LOGISTICREGRESSION_H
#define PROCESSLIB_LOGISTICREGRESSION_H

#include "Process.h"
#include "analysis/Regression.h"

namespace PLATO{
namespace ProcessLib{

class LinearRegression : public ProcessImpl<LinearRegression>, public Analysis::Regression {
private:
	const static std::string stepname;

public:
	LinearRegression() : ProcessImpl<LinearRegression>(stepname, "Run Linear Regression") {};
	virtual ~LinearRegression(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);
	virtual Result* calculate(double* data, unsigned int n_cols, unsigned int n_rows, const Result* null_result);

};

}
}

#endif
