#ifndef ANALYSIS_ENCODING_H
#define ANALYSIS_ENCODING_H

#include <iostream>
#include <string>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace PLATO{
namespace Analysis{

class Encoding{
public:
	enum encoding_ENUM { ADDITIVE, DOMINANT, RECESSIVE, WEIGHTED, CODOMINANT};
};

class EncodingModel{
public:
	friend class boost::serialization::access;

	EncodingModel(const char* s);
	EncodingModel(const std::string& s);
	EncodingModel() : _data(Encoding::ADDITIVE) {}
	EncodingModel(Encoding::encoding_ENUM c) : _data(c) {}

	template <class Archive>
	void serialize(Archive& ar, const unsigned int){
		ar & _data;
	}

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
