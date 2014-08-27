#ifndef PROCESSLIB_LINEARREGRESSION_H
#define PROCESSLIB_LINEARREGRESSION_H

#include "Process.h"
#include "analysis/Regression.h"
#include "MPIProcess.h"

namespace PLATO{
namespace ProcessLib{

class LinearRegression : public ProcessImpl<LinearRegression>,
	public virtual Analysis::Regression, public MPIProcessImpl<LinearRegression> {

	friend class AutoRegression;

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
	virtual boost::program_options::options_description getExtraOptions();

	virtual calc_fn& getCalcFn() const;

	static Result* calculate(const double* result, const double* data,
			unsigned int n_cols, unsigned int n_rows, unsigned int offset,
			unsigned int n_covars, bool run_null, const Regression::ExtraData* extra_data);

private:
	static std::pair<float, float> calcPVal(Result* r, Result* submodel, double chisq, double tss,
			                         unsigned int n_rows, unsigned int n_indep);
public:
	static void calculate_MPI(unsigned int bufsz, const char* buf, 
		std::deque<std::pair<unsigned int, const char*> >& result_queue, boost::mutex& result_mutex, boost::condition_variable& cv);


};

}
}

#endif
