/*
 * BiallelicSample.h
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#ifndef DATA_BIALLELIC_SAMPLE_H
#define DATA_BIALLELIC_SAMPLE_H

#include "Sample.h"

#include <boost/dynamic_bitset.hpp>


namespace PLATO{
namespace Data{

class BiallelicSample : public Sample {

public:
	BiallelicSample(const std::string& famid, const std::string& id, unsigned int n_genos);

	virtual void appendGenotype(unsigned char geno1, unsigned char geno2);
	virtual void appendMissingGenotype();
	virtual bool isMissing(const Marker&) const;
	virtual std::pair<unsigned char, unsigned char> getGeno(const Marker&) const;

	virtual ~BiallelicSample(){}

private:

	boost::dynamic_bitset<> _genotype;
	unsigned int _max_pos;

};

}
}

#endif /* BIALLELICSAMPLE_H_ */
