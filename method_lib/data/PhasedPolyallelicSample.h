/*
 * PhasedPolyallelicSample.h
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#ifndef DATA_PHASED_POLYALLELIC_SAMPLE_H
#define DATA_PHASED_POLYALLELIC_SAMPLE_H

#include "Sample.h"

#include <deque>

namespace PLATO {
namespace Data{

class PhasedPolyallelicSample : public Sample {
public:
	PhasedPolyallelicSample(const std::string& famid, const std::string& id, unsigned int n_genos);

	virtual void appendGenotype(unsigned char geno1, unsigned char geno2);
	virtual void appendMissingGenotype();
	virtual void setGenotype(const Marker& m, unsigned char geno1, unsigned char geno2);
	virtual void setMissingGenotype(const Marker& m);
	virtual bool isMissing(const Marker&) const;
	virtual std::pair<unsigned char, unsigned char> getGeno(const Marker&) const;

	virtual ~PhasedPolyallelicSample(){}

private:

	std::deque<unsigned char> _genotype;
};

}
}

#endif /* PHASEDPOLYALLELICSAMPLE_H_ */
