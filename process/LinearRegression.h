#ifndef PROCESSLIB_LOGISTICREGRESSION_H
#define PROCESSLIB_LOGISTICREGRESSION_H

#include "Process.h"
#include "analysis/Regression.h"
#include "MPIProcess.h"

namespace PLATO{
namespace ProcessLib{

class LinearRegression : public ProcessImpl<LinearRegression>,
	public Analysis::Regression, public MPIProcessImpl<LinearRegression> {
private:
	const static std::string stepname;
	const static std::string MPIName;

public:
	LinearRegression() : ProcessImpl<LinearRegression>(stepname, "Run Linear Regression"),
		MPIProcessImpl<LinearRegression>(stepname) {};
	virtual ~LinearRegression(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Data::DataSet&);
	virtual bool initData();
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

	virtual calc_fn& getCalcFn() const;

	static Result* calculate(const double* result, const double* data,
			unsigned int n_cols, unsigned int n_rows, unsigned int offset,
			unsigned int n_covars, bool run_null, const Regression::ExtraData* extra_data);

private:
	static std::pair<float, float> calcPVal(Result* r, Result* submodel, double chisq, double tss,
			                         unsigned int n_rows, unsigned int n_indep);
public:
	static std::pair<unsigned int, const char *> calculate_MPI(unsigned int bufsz, const char* buf);


};

}
}

#endif
