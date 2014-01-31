#include "Encoding.h"

#include <boost/program_options.hpp>
#include <sstream>

using boost::program_options::validation_error;
using std::istream;
using std::ostream;
using std::string;

namespace Methods{
namespace Analysis{

EncodingModel::EncodingModel(const std::string& s){
	std::istringstream ss(s);
	ss >> (*this);
}

EncodingModel::EncodingModel(const char * s){
	std::istringstream ss(s);
	ss >> (*this);
}

int EncodingModel::operator()(unsigned int add) const{
	switch(_data){
	case DOMINANT:
		return add > 0;
	case RECESSIVE:
		return add == 2;
	default:
		return add;
	}
}

istream& operator>>(istream& in,
		Methods::Analysis::EncodingModel& model_out) {
	string token;
	in >> token;
	if (token.size() > 0) {
		char s = token[0];
		if (s == 'a' || s == 'A') {
			model_out = Methods::Analysis::ADDITIVE;
		} else if (s == 'd' || s == 'D') {
			model_out = Methods::Analysis::DOMINANT;
		} else if (s == 'r' || s == 'R') {
			model_out = Methods::Analysis::RECESSIVE;
		} else if (s == 'c' || s == 'C') {
			model_out = Methods::Analysis::CATEGORICAL;
		} else {
			throw validation_error(validation_error::invalid_option_value);
		}
	} else {
		throw validation_error(validation_error::invalid_option_value);
	}
	//    else throw boost::program_options::validation_error("Invalid unit");
	return in;
}

ostream& operator<<(ostream& o, const Methods::Analysis::EncodingModel& m){
	switch(m){
	case Methods::Analysis::ADDITIVE:
		return o << "additive";
	case Methods::Analysis::DOMINANT:
		return o << "dominant";
	case Methods::Analysis::RECESSIVE:
		return o << "recessive";
	case Methods::Analysis::CATEGORICAL:
		return o << "categorical";
	default:
		return o << "unknown";
	}
}

}
}
