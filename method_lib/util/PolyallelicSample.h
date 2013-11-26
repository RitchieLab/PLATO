/*
 * PolyallelicSample.h
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#ifndef METHODS_POLYALLELIC_SAMPLE_H
#define METHODS_POLYALLELIC_SAMPLE_H

#include "Sample.h"

namespace Methods {

class PolyallelicSample : public Sample {
public:
	PolyallelicSample(const std::string& famid, const std::string& id, unsigned int n_genos);

	virtual void appendGenotype(unsigned char geno1, unsigned char geno2);
	virtual void appendMissingGenotype();
	virtual bool isMissing(unsigned int pos) const;
	virtual std::pair<unsigned char, unsigned char> getGeno(unsigned int pos) const;

	virtual ~PolyallelicSample(){}

private:

	std::deque<unsigned char> _genotype;
	static unsigned char _missing_val;

};

}

#endif /* POLYALLELICSAMPLE_H_ */
