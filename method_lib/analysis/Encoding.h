#ifndef METHODS_ANALYSIS_ENCODING_H
#define METHODS_ANALYSIS_ENCODING_H

#include <iostream>
#include <string>

namespace Methods{
namespace Analysis{

enum encoding_ENUM { ADDITIVE, DOMINANT, RECESSIVE, CATEGORICAL};

class EncodingModel{
public:
	EncodingModel(const char* s);
	EncodingModel(const std::string& s);
	EncodingModel() : _data(ADDITIVE) {}
	EncodingModel(encoding_ENUM c) : _data(c) {}

	operator int() const{return _data;}

	int operator()(unsigned int add) const;

private:
	encoding_ENUM _data;
};

std::istream& operator>>(std::istream& in, Methods::Analysis::EncodingModel& model_out);
std::ostream& operator<<(std::ostream& o, const Methods::Analysis::EncodingModel& m);

}
}

#endif
