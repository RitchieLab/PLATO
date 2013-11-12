/*
 * PhasedBiallelicSample.h
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#ifndef METHODS_PHASED_BIALLELICSAMPLE_H
#define METHODS_PHASED_BIALLELICSAMPLE_H

#include "Sample.h"

#include <boost/dynamic_bitset.hpp>

namespace Methods {

class PhasedBiallelicSample : public Sample {

public:
	PhasedBiallelicSample(const std::string& famid, const std::string& id, unsigned int n_genos);

	virtual void appendGenotype(unsigned char geno1, unsigned char geno2);
	virtual bool isMissing(unsigned int pos) const;
	virtual std::pair<unsigned char, unsigned char> getGeno(unsigned int pos) const;

	virtual ~PhasedBiallelicSample(){}

private:

	boost::dynamic_bitset _genotype;
	boost::dynamic_bitset _missing;

};

}

#endif /* PHASEDBIALLELICSAMPLE_H_ */
