/*
 * OutputPLINK.h
 *
 *  Created on: Dec 3, 2013
 *      Author: jrw32
 */

#ifndef METHODS_OUTPUTPLINK_H
#define METHODS_OUTPUTPLINK_H

#include <ostream>

namespace Methods {

class Marker;
class Sample;

class OutputPLINK {
public:
	OutputPLINK(){}
	~OutputPLINK(){}

protected:

	void printPEDHeader(std::ostream&, const Sample&) const;
	void printMAPInfo(std::ostream&, const Marker&, bool print_alleles=false) const;


};

}

#endif /* OUTPUTPLINK_H_ */
