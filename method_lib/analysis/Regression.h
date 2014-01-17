/*
 * Regression.h
 *
 *  Created on: Jan 8, 2014
 *      Author: jrw32
 */

#ifndef METHODS_ANALYSIS_REGRESSION_H
#define METHODS_ANALYSIS_REGRESSION_H

#include <set>
#include <vector>
#include <string>

#include <boost/program_options.hpp>

#include "Correction.h"

//#include "method_lib/util/Container.h"

namespace Methods{

namespace Analysis{

/*!
 * \brief A base class for any regression analysis - really any statistical test
 */
class Regression {

protected:
	/*!
	 * \brief A class that holds all data needed to process a single regression
	 */
	class Model{

	};

	/*!
	 * \brief A class that holds the result of the regression and can be used
	 * to print the results.
	 */
	class Result{

	};

public:

	Regression();
	virtual ~Regression();

	boost::program_options::options_description& addOptions(boost::program_options::options_description& opts);
	void parseOptions(const boost::program_options::variables_map& vm);

protected:

	//! a set of correction methods to apply to the p-values
	std::set<CorrectionModel> corr_methods;
	//! covariates to use in the regression
	std::set<std::string> covar_names;
	//! A string of the outcome name
	std::string outcome_name;
	//! a file of models to use
	std::vector<std::string> model_files;


	//! raw p-value cutoff for displaying models
	float cutoff_p;
	//! exclude markers? (only use covariates?)
	bool exclude_markers;
	//! build models with traits as well?
	bool include_traits;
	//! include interactions?
	bool interactions;
	//! autogenerate pairwise models
	bool pairwise;


};

}

}


#endif /* REGRESSION_H_ */
