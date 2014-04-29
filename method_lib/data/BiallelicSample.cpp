/*
 * BiallelicSample.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "BiallelicSample.h"

#include "Marker.h"

#include <algorithm>

using std::string;
using std::pair;
using std::make_pair;
using boost::dynamic_bitset;

namespace PLATO {
namespace Data{

BiallelicSample::BiallelicSample(const string& famid, const string& id, unsigned int n_genos) :
	Sample(famid, id) {
	//_genotype.resize(2*n_genos);
	// We want to initialize to all missing, so we'll precompute our mask, which
	// should be 10101010...
	dynamic_bitset<>::block_type mask = 0;
	for(unsigned int i=0; i<dynamic_bitset<>::bits_per_block; i+=2){
		mask <<= 2;
		mask |= 2;
	}

	// now, just tack on the masks an appropriate number of times
	for(unsigned int i=0; i<2*n_genos; i+=dynamic_bitset<>::bits_per_block){
		_genotype.append(mask);
	}
	_genotype.resize(2*n_genos);

}

void BiallelicSample::appendGenotype(unsigned char geno1, unsigned char geno2){

	// If either is "missing", then both are missing!
	if(geno1 == missing_allele || geno2 == missing_allele){
		appendMissingGenotype();
	}else{
		if(geno1 && !geno2){
			std::swap(geno1, geno2);
		}
		// Note: this can still be "missing" if we give [1,0], so this will take
		// Binary PED without modification.
		_genotype.push_back(geno1);
		_genotype.push_back(geno2);

	}
}

void BiallelicSample::appendMissingGenotype(){
	_genotype.push_back(1);
	_genotype.push_back(0);
}

void BiallelicSample::setGenotype(const Marker& m, unsigned char geno1, unsigned char geno2){
	unsigned int pos = m.getIndex();

	// If either is "missing", then both are missing!
	if(geno1 == missing_allele || geno2 == missing_allele){
		_genotype.set(2*pos);
		_genotype.reset(2*pos + 1);
	}else{
		if(geno1 && !geno2){
			std::swap(geno1, geno2);
		}
		// Note: this can still be "missing" if we give [1,0], so this will take
		// Binary PED without modification.
		_genotype.set(2*pos, geno1);
		_genotype.set(2*pos + 1, geno2);

	}
}

void BiallelicSample::setMissingGenotype(const Marker& m){
	unsigned int pos = m.getIndex();
	_genotype.set(2*pos);
	_genotype.reset(2*pos + 1);
}

bool BiallelicSample::isMissing(const Marker& m) const{
	unsigned int pos = m.getIndex();
	return _genotype[2*pos] && (!_genotype[2*pos + 1]);
}

pair<unsigned char, unsigned char> BiallelicSample::getGeno(const Marker& m) const{
	unsigned int pos = m.getIndex();
	if(isMissing(m)){
		return missing_geno;
	}else{
		return make_pair(_genotype[2*pos] ? m.getAltIdx() : m.getRefIdx(),
				_genotype[2*pos + 1] ? m.getAltIdx() : m.getRefIdx());
	}
}

}
}
