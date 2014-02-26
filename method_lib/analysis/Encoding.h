#ifndef ANALYSIS_ENCODING_H
#define ANALYSIS_ENCODING_H

#include <iostream>
#include <string>

namespace PLATO{
namespace Analysis{

class Encoding{
public:
	enum encoding_ENUM { ADDITIVE, DOMINANT, RECESSIVE, WEIGHTED, CODOMINANT};
};

class EncodingModel{
public:
	EncodingModel(const char* s);
	EncodingModel(const std::string& s);
	EncodingModel() : _data(Encoding::ADDITIVE) {}
	EncodingModel(Encoding::encoding_ENUM c) : _data(c) {}

	operator int() const{return _data;}

	int operator()(unsigned int add) const;

private:
	Encoding::encoding_ENUM _data;
};

std::istream& operator>>(std::istream& in, PLATO::Analysis::EncodingModel& model_out);
std::ostream& operator<<(std::ostream& o, const PLATO::Analysis::EncodingModel& m);

}
}

#endif
