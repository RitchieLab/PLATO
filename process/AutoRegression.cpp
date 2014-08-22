#include "AutoRegression.h"

#include <limits>
#include <utility>
#include <cmath>
#include <cstring>
#include <numeric>

#include <boost/lexical_cast.hpp>
#include <boost/serialization/export.hpp>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include "util/Logger.h"
#include "util/GSLUtils.h"

using PLATO::Data::DataSet;
using PLATO::Analysis::Regression;
using PLATO::Analysis::Encoding;
using PLATO::Utility::Logger;

using std::vector;
using std::string;
using std::set;
using std::numeric_limits;
using std::fabs;
using std::log;
using std::exp;
using std::pow;
using std::pair;
using std::map;
using std::deque;

using boost::array;

namespace po=boost::program_options;


BOOST_CLASS_EXPORT(PLATO::ProcessLib::AutoRegression::ExtraData)

namespace PLATO{
namespace ProcessLib{

const std::string AutoRegression::stepname = AutoRegression::doRegister("regress-auto");
const std::string AutoRegression::MPIname = AutoRegression::registerMPI("regress-auto");

po::options_description& AutoRegression::appendOptions(po::options_description& opts){
	Regression::addOptions(opts);

	po::options_description autoreg_opts("Auto Regression Options");

	opts.add(logreg_opts);

	return opts;
}

void AutoRegression::parseOptions(const po::variables_map& vm){
	Regression::parseOptions(vm);
}

void AutoRegression::printVarHeader(const std::string& var_name, std::ofstream& of) const{
	LogisticRegression::printVarHeader(var_name, of);
}

bool AutoRegression::initData(){

	// make sure the class_data is initialized, please!
	getExtraData();

	bool good_pheno = true;
	//transform the phenotype from [0,1]
	set<float> uniq_pheno;
	for(unsigned int i=0; i<_pheno.size(); i++){
		if(!std::isnan(_pheno[i])){
			uniq_pheno.insert(_pheno[i]);
		}
	}

	if(uniq_pheno.size() > 2){
		class_data->analysis_type = 0;
	} else if(uniq_pheno.size() < 2){
		Utility::Logger::log_err("ERROR: Desired phenotype has only " +
				boost::lexical_cast<string>(uniq_pheno.size()) +
				" unique value(s); Regression will almost certainly fail!", outcome_names.size() <= 1);
		good_pheno = false;
	} else {
		LogisticRegression::initData();
		class_data->analysis_type = 1;
	}

	return good_pheno;
}

Regression::calc_fn& AutoRegression::getCalcFn() const{
	return (AutoRegression::calculate);
}

const Regression::ExtraData* AutoRegression::getExtraData() const{
	if(!class_data){
		class_data = new ExtraData(*LogisticRegression::getExtraData());
		class_data->analysis_type = 0;
	}

	return class_data;
}

Regression::Result* AutoRegression::calculate(
		const double* Y, const double* data,
		unsigned int n_cols, unsigned int n_rows, unsigned int offset,
		unsigned int n_covars, bool run_null, const Regression::ExtraData* other_data){

	const ExtraData* extra_data = dynamic_cast<const ExtraData*>(other_data);
	if(!extra_data){
		// something went VERY wrong here!
		return 0;
	}

	Result* r = 0;
	string analysis;
	// use the approprate calculate method
	if(extra_data->analysis_type == 0){
		analysis="linear";
		r = LinearRegression::calculate(Y, data, n_cols, n_rows, offset, n_covars, run_null, other_data);
	} else if (extra_data->analysis_type == 1){
		analysis="logistic";
		r = LogisticRegression::calculate(Y, data, n_cols, n_rows, offset, n_covars, run_null, other_data);
	}

	// add the analysis type to the result
	if(r){
		r->prefix += analysis + extra_data->sep;
	}

	return r;

}

void AutoRegression::process(DataSet& ds){
	runRegression(ds);
}

void AutoRegression::printExtraHeader(std::ofstream& of){
	of << "Analysis_Type" << sep;
	LinearRegression::printExtraHeader(of);
}

string AutoRegression::printExtraResults(const Result& r){
	return boost::lexical_cast<string>(r.converged) + sep;
}

void AutoRegression::calculate_MPI(unsigned int bufsz, const char* buf,
	deque<pair<unsigned int, const char*> >& result_queue, boost::mutex& result_mutex, boost::condition_variable& cv){
	Regression::calculate_MPI(bufsz, buf, result_queue, result_mutex, cv, AutoRegression::calculate);
}

}
}
