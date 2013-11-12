/*
 * BiallelicSample.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "BiallelicSample.h"

using std::string;
using std::pair;
using std::make_pair;

namespace Methods {

BiallelicSample::BiallelicSample(string famid, string id, unsigned int n_genos) :
	Sample(famid, id) {

	_genotype.resize(2*n_genos);

}

void BiallelicSample::appendGenotype(unsigned char geno1, unsigned char geno2){

	// If either is "missing", then both are missing!
	if(geno1 == static_cast<unsigned char>(-1) || geno2 == static_cast<unsigned char>(-1)){
		_genotype.push_back(1);
		_genotype.push_back(0);
	}else{
		// Note: this can still be "missing" if we give [1,0], so this will take
		// Binary PED without modification.
		_genotype.push_back(geno1);
		_genotype.push_back(geno2);
	}
}

bool BiallelicSample::isMissing(unsigned int pos) const{
	return _genotype[2*pos] && (!_genotype[2*pos + 1]);
}

pair<unsigned char, unsigned char> BiallelicSample::getGeno(unsigned int pos) const{
	if(isMissing(pos)){
		return missing_geno;
	}else{
		return make_pair(_genotype[2*pos], _genotype[2*pos + 1]);
	}
}

}
