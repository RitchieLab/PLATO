/*
 * PolyallelicSample.h
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#ifndef DATA_POLYALLELIC_SAMPLE_H
#define DATA_POLYALLELIC_SAMPLE_H

#include "Sample.h"

namespace PLATO{
namespace Data{

class PolyallelicSample : public Sample {
public:
	PolyallelicSample(const std::string& famid, const std::string& id, unsigned int n_genos);

	virtual void appendGenotype(unsigned char geno1, unsigned char geno2);
	virtual void appendMissingGenotype();
	virtual bool isMissing(const Marker&) const;
	virtual std::pair<unsigned char, unsigned char> getGeno(const Marker&) const;

	virtual ~PolyallelicSample(){}

private:

	std::deque<unsigned char> _genotype;

};

}
}

#endif /* POLYALLELICSAMPLE_H_ */
