#include "AutoRegression.h"

#include <limits>
#include <utility>
#include <cmath>
#include <cstring>
#include <numeric>
#include <algorithm>

#include <boost/lexical_cast.hpp>
#include <boost/serialization/export.hpp>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include "util/Logger.h"
#include "util/GSLUtils.h"
#include "util/InputManager.h"

using PLATO::Data::DataSet;
using PLATO::Analysis::Regression;
using PLATO::Analysis::Encoding;
using PLATO::Utility::Logger;
using PLATO::Utility::InputManager;

using std::vector;
using std::string;
using std::set;
using std::deque;
using std::pair;

using boost::array;

namespace po=boost::program_options;


BOOST_CLASS_EXPORT(PLATO::ProcessLib::AutoRegression::ExtraData)

namespace PLATO{
namespace ProcessLib{

const std::string AutoRegression::stepname = AutoRegression::doRegister("regress-auto");
const std::string AutoRegression::MPIname = AutoRegression::registerMPI("regress-auto");

po::options_description& AutoRegression::appendOptions(po::options_description& opts){
	Regression::addOptions(opts);

	opts.add(LogisticRegression::getExtraOptions());
	opts.add(LinearRegression::getExtraOptions());
	opts.add(getExtraOptions());

	return opts;
}

po::options_description AutoRegression::getExtraOptions(){
	po::options_description autoreg_opts("Auto Regression Options");

	autoreg_opts.add_options()
		("linear", po::value<vector<string> >(&linear_traits)->composing(), "Outcomes to force linear regression")
		("logistic", po::value<vector<string> >(&logistic_traits)->composing(), "Outcomes to force logistic regression")
		;

	return autoreg_opts;
}

void AutoRegression::parseOptions(const po::variables_map& vm){
	Regression::parseOptions(vm);

	if(vm.count("linear")){
		InputManager::parseInput(vm["linear"].as<vector<string> >(), linear_set);
	}

	if(vm.count("logistic")){
		InputManager::parseInput(vm["logistic"].as<vector<string> >(), logistic_set);
	}

	// make sure the 2 sets are mutually exclusive
	set<string> tmp_set;
	std::set_intersection(linear_set.begin(), linear_set.end(),
			logistic_set.begin(), logistic_set.end(),
			std::inserter(tmp_set, tmp_set.begin()));

	if(tmp_set.size() > 0){
		Logger::log_err("ERROR: --linear and --logistic must be mututally exclusive", true);
	}

	std::set_union(linear_set.begin(), linear_set.end(),
				logistic_set.begin(), logistic_set.end(),
				std::inserter(tmp_set, tmp_set.begin()));

	set<string> tmp_set2;
	std::set_intersection(tmp_set.begin(), tmp_set.end(),
			outcome_names.begin(), outcome_names.end(),
			std::inserter(tmp_set2, tmp_set2.begin()));

	if(tmp_set.size() != tmp_set2.size() && !_phewas){
		Logger::log_err("ERROR: Outcomes given in --linear and --logistic must also be given in --outcomes (or use --phewas)", true);
	}

}

void AutoRegression::printVarHeader(const std::string& var_name, std::ofstream& of) const{
	LogisticRegression::printVarHeader(var_name, of);
}

bool AutoRegression::initData(){

	// make sure the class_data is initialized, please!
	getExtraData();

	bool good_pheno = true;

	// first, check to see if this is in either linear or logistic
	if(linear_set.count(*output_itr)){
		analysis_type = class_data->analysis_type = LINEAR;
	} else if(logistic_set.count(*output_itr)){
		analysis_type = class_data->analysis_type = LOGISTIC;
	} else {
		set<float> uniq_pheno;
		for(unsigned int i=0; i<_pheno.size() && uniq_pheno.size() <= 2; i++){
			if(!std::isnan(_pheno[i])){
				uniq_pheno.insert(_pheno[i]);
			}
		}

		if(uniq_pheno.size() > 2){
			analysis_type = class_data->analysis_type = LINEAR;
		} else if(uniq_pheno.size() == 2){
			analysis_type = class_data->analysis_type = LOGISTIC;
		} else {
			Utility::Logger::log_err("ERROR: Desired phenotype has only " +
					boost::lexical_cast<string>(uniq_pheno.size()) +
					" unique value(s); Regression will almost certainly fail!", outcome_names.size() <= 1);
			good_pheno = false;
		}
	}

	//transform the phenotype from [0,1]
	if(good_pheno && analysis_type == LOGISTIC){
		good_pheno = LogisticRegression::initData();
	} else if(good_pheno && analysis_type == LINEAR){
		good_pheno = LinearRegression::initData();
	}

	return good_pheno;
}

Regression::calc_fn& AutoRegression::getCalcFn() const{
	return AutoRegression::calculate;
}

const Regression::ExtraData* AutoRegression::getExtraData() const{
	if(!class_data){
		class_data = new ExtraData(*LogisticRegression::getExtraData());
		class_data->analysis_type = LINEAR;
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
		r->prefix += analysis + extra_data->sep + r->prefix;
	}

	return r;

}

void AutoRegression::process(DataSet& ds){
	Regression::runRegression(ds);
}

void AutoRegression::printExtraHeader(std::ofstream& of){
	of << "Analysis_Type" << sep;
	LogisticRegression::printExtraHeader(of);
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
