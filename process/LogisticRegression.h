#ifndef PROCESSLIB_LOGISTICREGRESSION_H
#define PROCESSLIB_LOGISTICREGRESSION_H

#include "Process.h"
#include "analysis/Regression.h"

namespace PLATO{
namespace ProcessLib{

class LogisticRegression : public ProcessImpl<LogisticRegression>, public Analysis::Regression {
private:
	const static std::string stepname;

public:
	LogisticRegression() : ProcessImpl<LogisticRegression>(stepname, "Run Logistic Regression") {};
	virtual ~LogisticRegression(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);
	virtual Result* calculate(const double* Y, const double* data,
			unsigned int n_cols, unsigned int n_rows, unsigned int offset,
			unsigned int n_covars);

	virtual void printVarHeader(const std::string& var_name);
	virtual bool initData(const Data::DataSet& ds);

	virtual void printExtraHeader();
	virtual void printExtraResults(const Result& r);

private:
	bool show_odds;
	unsigned int maxIterations;

};

}
}

#endif
