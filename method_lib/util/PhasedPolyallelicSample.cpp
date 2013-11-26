/*
 * PhasedPolyallelicSample.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "PhasedPolyallelicSample.h"

using std::string;
using std::pair;

namespace Methods{

unsigned char PhasedPolyallelicSample::_missing_val = static_cast<unsigned char>(-1);

PhasedPolyallelicSample::PhasedPolyallelicSample(const string& famid, const string& id, unsigned int n_genos) :
	Sample(famid, id) {
}

void PhasedPolyallelicSample::appendGenotype(unsigned char geno1, unsigned char geno2){
	_genotype.push_back(geno1);
	_genotype.push_back(geno2);
}

void PhasedPolyallelicSample::appendMissingGenotype(){
	_genotype.push_back(_missing_val);
	_genotype.push_back(_missing_val);
}

bool PhasedPolyallelicSample::isMissing(unsigned int pos) const{
	return _genotype[2*pos] == _missing_val || _genotype[2*pos + 1] == _missing_val;
}

pair<unsigned char, unsigned char> PhasedPolyallelicSample::getGeno(unsigned int pos) const{
	return std::make_pair(_genotype[2*pos], _genotype[2*pos + 1]);
}

}
