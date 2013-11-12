/*
 * PhasedPolyallelicSample.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "PhasedPolyallelicSample.h"

namespace Methods{

PhasedPolyallelicSample::PhasedPolyallelicSample(string famid, string id, unsigned int n_genos) :
	Sample(famid, id) {
}

void PhasedPolyallelicSample::appendGenotype(unsigned char geno1, unsigned char geno2){
	_genotype.push_back(geno1);
	_genotype.push_back(geno2);
}

bool PhasedPolyallelicSample::isMissing(unsigned int pos) const{
	return _genotype[2*pos] == _missing_val || _genotype[2*pos + 1] == _missing_val;
}

pair<unsigned char, unsigned char> PhasedPolyallelicSample::getGeno(unsigned int pos) const{
	return make_pair(_genotype[2*pos], _genotype[2*pos + 1]);
}

}
