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

namespace PLATO {
namespace Data{

BiallelicSample::BiallelicSample(const string& famid, const string& id, unsigned int n_genos) :
	Sample(famid, id) {
	//_genotype.resize(2*n_genos);

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
