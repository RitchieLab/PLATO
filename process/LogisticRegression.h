#ifndef PROCESSLIB_LOGISTICREGRESSION_H
#define PROCESSLIB_LOGISTICREGRESSION_H

#include "Process.h"
#include "analysis/Regression.h"

#include <vector>

namespace ProcessLib{

class LogisticRegression : public ProcessImpl<LogisticRegression>, public Methods::Analysis::Regression {
private:
	const static std::string stepname;

public:
	LogisticRegression() : ProcessImpl<LogisticRegression>(stepname, "Run Logistic Regression") {};
	virtual ~LogisticRegression(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);
	virtual void initData(const std::string& model_str, const Methods::DataSet& ds);
	virtual void printResults();

protected:
	virtual void process(Methods::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);
	virtual Methods::Analysis::Regression::Result* calculate(float* data, unsigned int n_cols, unsigned int n_rows);

private:
	bool show_odds;
	unsigned int maxIterations;

};

}

#endif
