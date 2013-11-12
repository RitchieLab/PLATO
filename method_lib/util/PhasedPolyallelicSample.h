/*
 * PhasedPolyallelicSample.h
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#ifndef METHDOS_PHASED_POLYALLELIC_SAMPLE_H
#define METHDOS_PHASED_POLYALLELIC_SAMPLE_H

#include <deque>

namespace Methods {

class PhasedPolyallelicSample {
public:
	PhasedPolyallelicSample(const std::string& famid, const std::string& id, unsigned int n_genos);

	virtual void appendGenotype(unsigned char geno1, unsigned char geno2);
	virtual bool isMissing(unsigned int pos) const;
	virtual std::pair<unsigned char, unsigned char> getGeno(unsigned int pos) const;

	virtual ~PhasedPolyallelicSample(){}

private:

	std::deque<unsigned char> _genotype;
	static unsigned char _missing_val = static_cast<unsigned char>(-1);

};

}

#endif /* PHASEDPOLYALLELICSAMPLE_H_ */
