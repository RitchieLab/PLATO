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

PhasedPolyallelicSample::PhasedPolyallelicSample(const string& famid, const string& id, unsigned int n_genos) :
	Sample(famid, id) {
}

void PhasedPolyallelicSample::appendGenotype(unsigned char geno1, unsigned char geno2){
	_genotype.push_back(geno1);
	_genotype.push_back(geno2);
}

void PhasedPolyallelicSample::appendMissingGenotype(){
	_genotype.push_back(missing_allele);
	_genotype.push_back(missing_allele);
}

bool PhasedPolyallelicSample::isMissing(unsigned int pos) const{
	return _genotype[2*pos] == missing_allele || _genotype[2*pos + 1] == missing_allele;
}

pair<unsigned char, unsigned char> PhasedPolyallelicSample::getGeno(unsigned int pos) const{
	return std::make_pair(_genotype[2*pos], _genotype[2*pos + 1]);
}

}
