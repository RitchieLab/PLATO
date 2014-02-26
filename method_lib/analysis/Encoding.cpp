#include "Encoding.h"

#include <boost/program_options.hpp>
#include <sstream>

using boost::program_options::validation_error;
using std::istream;
using std::ostream;
using std::string;

namespace PLATO{
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
	case Encoding::DOMINANT:
		return add > 0;
	case Encoding::RECESSIVE:
		return add == 2;
	default:
		return add;
	}
}

istream& operator>>(istream& in,
		PLATO::Analysis::EncodingModel& model_out) {
	string token;
	in >> token;
	if (token.size() > 0) {
		char s = token[0];
		if (s == 'a' || s == 'A') {
			model_out = PLATO::Analysis::Encoding::ADDITIVE;
		} else if (s == 'd' || s == 'D') {
			model_out = PLATO::Analysis::Encoding::DOMINANT;
		} else if (s == 'r' || s == 'R') {
			model_out = PLATO::Analysis::Encoding::RECESSIVE;
		} else if (s == 'w' || s == 'W') {
			model_out = PLATO::Analysis::Encoding::WEIGHTED;
		} else if (s == 'c' || s == 'C') {
			model_out = PLATO::Analysis::Encoding::CODOMINANT;
		} else {
			throw validation_error(validation_error::invalid_option_value);
		}
	} else {
		throw validation_error(validation_error::invalid_option_value);
	}
	//    else throw boost::program_options::validation_error("Invalid unit");
	return in;
}

ostream& operator<<(ostream& o, const PLATO::Analysis::EncodingModel& m){
	switch(m){
	case PLATO::Analysis::Encoding::ADDITIVE:
		return o << "additive";
	case PLATO::Analysis::Encoding::DOMINANT:
		return o << "dominant";
	case PLATO::Analysis::Encoding::RECESSIVE:
		return o << "recessive";
	case PLATO::Analysis::Encoding::WEIGHTED:
		return o << "weighted";
	case PLATO::Analysis::Encoding::CODOMINANT:
		return o << "codominant";
	default:
		return o << "unknown";
	}
}

}
}
