/*
 * BiallelicSample.h
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#ifndef METHODS_BIALLELIC_SAMPLE_H
#define METHODS_BIALLELIC_SAMPLE_H

#include "Sample.h"

#include <boost/dynamic_bitset.hpp>


namespace Methods{

class BiallelicSample : public Sample {

public:
	BiallelicSample(const std::string& famid, const std::string& id, unsigned int n_genos);

	virtual void appendGenotype(unsigned char geno1, unsigned char geno2);
	virtual bool isMissing(unsigned int pos) const;
	virtual std::pair<unsigned char, unsigned char> getGeno(unsigned int pos) const;

	virtual ~BiallelicSample(){}

private:

	boost::dynamic_bitset _genotype;

};

}

#endif /* BIALLELICSAMPLE_H_ */
