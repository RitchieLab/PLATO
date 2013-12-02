/*
 * PhasedBiallelicSample.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "PhasedBiallelicSample.h"

using std::string;
using std::pair;
using std::make_pair;

namespace Methods{

PhasedBiallelicSample::PhasedBiallelicSample(const string& famid, const string& id, unsigned int n_genos) :
	Sample(famid, id) {

	//_genotype.resize(2*n_genos);
	//_missing.resize(2*n_genos);

}

void PhasedBiallelicSample::appendGenotype(unsigned char geno1, unsigned char geno2){

	_missing.push_back(geno1 == missing_allele);
	_missing.push_back(geno2 == missing_allele);

	_genotype.push_back(geno1);
	_genotype.push_back(geno2);
}

void PhasedBiallelicSample::appendMissingGenotype(){
	// If either is "missing", then both are missing!
	_missing.push_back(true);
	_missing.push_back(true);

	_genotype.push_back(0);
	_genotype.push_back(0);
}

bool PhasedBiallelicSample::isMissing(unsigned int pos) const{
	return _missing[2*pos] || _missing[2*pos + 1];
}

pair<unsigned char, unsigned char> PhasedBiallelicSample::getGeno(unsigned int pos) const{

	return make_pair(_missing[2*pos] ? missing_allele : _genotype[2*pos],
			_missing[2*pos+1] ? missing_allele : _genotype[2*pos+1]);
}

}
