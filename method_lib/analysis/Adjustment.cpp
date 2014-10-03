#include "Adjustment.h"

#include "util/Logger.h"

#include <algorithm>
#include <string>
#include <sstream>

#include <boost/program_options.hpp>

#include <gsl/gsl_cdf.h>

using std::map;
using std::vector;
using std::sort;
using std::stringstream;

using PLATO::Utility::Logger;

using boost::program_options::validation_error;

namespace PLATO{
namespace Analysis{

AdjustmentModel::AdjustmentModel(const std::string& s){
	std::istringstream ss(s);
	ss >> (*this);
}

AdjustmentModel::AdjustmentModel(const char * s){
	std::istringstream ss(s);
	ss >> (*this);
}

Adjustment* Adjustment::getAdjustmentMethod(const AdjustmentModel& c){
	static Adjustment* reg = new RegressionAdjustment();
	static Adjustment* gc = new GCAdjustment();

	switch((int) c){
	case REGRESSION:
		return reg;
	case GC:
		return gc;
	default:
		return 0;
	}
}

std::string Adjustment::listAdjustmentMethods(){

	stringstream ss;
	ss << AdjustmentModel(REGRESSION) << "," << AdjustmentModel(GC);

	return ss.str();
}

double RegressionAdjustment::adjust(Wrapper cont) const{
	double num = 0;
	double den = 0;
	unsigned int sz = cont.size();
	// I'll need this for the "expected" p-value
	double yi = 0.5 / static_cast<double>(sz);
	double y_increment = 1 / static_cast<double>(sz);

	for(unsigned int i=0; i<sz; i++){
		double xi = cont[i];
		num += xi*yi;
		den += xi*xi;
		yi += y_increment;
	}

	return num / den;
}

double RegressionAdjustment::correct(double pval, double gif_recip) const{
	return std::min(pval / gif_recip, 1.0);
}

double GCAdjustment::adjust(Wrapper cont) const{
	// first, get the two flanking p-values that represent the median
	// (note that they may be identical)
	unsigned int sz = cont.size();
	double pv1 = cont[sz/2];
	double pv2 = cont[sz/2 - (sz + 1) % 2];

	// convert pv1 and pv2 into chi-squared values
	pv1 = gsl_cdf_chisq_Qinv(pv1, 1);
	pv2 = gsl_cdf_chisq_Qinv(pv2, 1);

	// average them
	double median = (pv1+pv2)/2;

	// and compare with the theoretic median
	return median / gsl_cdf_chisq_Qinv(0.5,1);
}

double GCAdjustment::correct(double pval, double gif_recip) const {
	return gsl_cdf_chisq_Q(gif_recip * gsl_cdf_chisq_Qinv(pval, 1), 1);
}

}
}

namespace std{

istream& operator>>(istream& in,
		PLATO::Analysis::AdjustmentModel& model_out) {
	string token;
	in >> token;
	if (token.size() > 0) {
		char s = token[0];
		if (s == 'r' || s == 'R') {
			model_out = PLATO::Analysis::Adjustment::REGRESSION;
		} else if (s == 'g' || s == 'G') {
			model_out = PLATO::Analysis::Adjustment::GC;
		} else {
			throw validation_error(validation_error::invalid_option_value);
		}
	} else {
		throw validation_error(validation_error::invalid_option_value);
	}
	//    else throw boost::program_options::validation_error("Invalid unit");
	return in;
}
ostream& operator<<(ostream& o, const PLATO::Analysis::AdjustmentModel& m){
	switch(m){
	case PLATO::Analysis::Adjustment::REGRESSION:
		return o << "Regression";
	case PLATO::Analysis::Adjustment::GC:
		return o << "GC";
	default:
		return o << "unknown";
	}
}

}
