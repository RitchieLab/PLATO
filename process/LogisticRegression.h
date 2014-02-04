#ifndef PROCESSLIB_LOGISTICREGRESSION_H
#define PROCESSLIB_LOGISTICREGRESSION_H

#include "Process.h"
#include "analysis/Regression.h"

namespace ProcessLib{

class LogisticRegression : public ProcessImpl<LogisticRegression>, public Methods::Analysis::Regression {
private:
	const static std::string stepname;

public:
	LogisticRegression() : ProcessImpl<LogisticRegression>(stepname, "Run Logistic Regression") {};
	virtual ~LogisticRegression(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Methods::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);
	virtual Methods::Analysis::Regression::Result* calculate(double* data, unsigned int n_cols, unsigned int n_rows);

	virtual void printVarHeader(const std::string& var_name);
	virtual void initData(const std::string& model_str, const Methods::DataSet& ds);

private:
	bool show_odds;
	unsigned int maxIterations;

};

}

#endif
