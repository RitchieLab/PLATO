#ifndef UTILITY_MISSING_H
#define UTILITY_MISSING_H

#include <string>

namespace PLATO{

namespace Data{
class DataSet;
class Marker;
class Sample;
}


namespace Utility{
/*!
 * \brief A class to calculate the missingness by sample or by marker
 */
class Missing{
public:
	static double markerMissing(const Data::DataSet& ds, const Data::Marker& m);
	static double sampleMissing(const Data::DataSet& ds, const Data::Sample& s);
	static double traitMissing(const Data::DataSet& ds, const std::string& t);
};

}
}

#endif
