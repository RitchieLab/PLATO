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
	virtual Analysis::Regression::Result* calculate(double* data, unsigned int n_cols, unsigned int n_rows);

	virtual void printVarHeader(const std::string& var_name);
	virtual void initData(const std::string& model_str, const Data::DataSet& ds);

private:
	bool show_odds;
	unsigned int maxIterations;

};

}
}

#endif
