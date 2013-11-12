/*
 * PolyallelicSample.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "PolyallelicSample.h"

using std::string;
using std::pair;
using std::make_pair;

namespace Methods{

PolyallelicSample::PolyallelicSample(string famid, string id, unsigned int n_genos) :
	Sample(famid, id) {
}

void PolyallelicSample::appendGenotype(unsigned char geno1, unsigned char geno2){

	// If either is "missing", then both are missing!
	if(geno1 == _missing_val || geno2 == _missing_val){
		_genotype.push_back(_missing_val);
		_genotype.push_back(_missing_val);
	}else{
		_genotype.push_back(geno1);
		_genotype.push_back(geno2);
	}
}

bool PolyallelicSample::isMissing(unsigned int pos) const{
	return _genotype[2*pos] == _missing_val || _genotype[2*pos + 1] == _missing_val;
}

pair<unsigned char, unsigned char> PolyallelicSample::getGeno(unsigned int pos) const{
	return make_pair(_genotype[2*pos], _genotype[2*pos + 1]);
}

}
