/*
 * OutputPLINK.h
 *
 *  Created on: Dec 3, 2013
 *      Author: jrw32
 */

#ifndef UTILITY_OUTPUTPLINK_H
#define UTILITY_OUTPUTPLINK_H

#include <ostream>

namespace PLATO {

namespace Data{
class Marker;
class Sample;
}

namespace Utility{
class OutputPLINK {
public:
	OutputPLINK(){}
	~OutputPLINK(){}

protected:

	void printPEDHeader(std::ostream&, const PLATO::Data::Sample&) const;
	void printMAPInfo(std::ostream&, const PLATO::Data::Marker&, bool print_alleles=false) const;


};

}
}

#endif /* OUTPUTPLINK_H_ */
