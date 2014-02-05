/*
 * PhasedBiallelicSample.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "PhasedBiallelicSample.h"

#include "Marker.h"

using std::string;
using std::pair;
using std::make_pair;

namespace PLATO{
namespace Data{

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

bool PhasedBiallelicSample::isMissing(const Marker& m) const{
	unsigned int pos = m.getIndex();
	return _missing[2*pos] || _missing[2*pos + 1];
}

pair<unsigned char, unsigned char> PhasedBiallelicSample::getGeno(const Marker& m) const{
	unsigned int pos = m.getIndex();
	unsigned char ref_idx = m.getRefIdx();
	unsigned char alt_idx = m.getAltIdx();

	return make_pair(_missing[2*pos] ? missing_allele : (_genotype[2*pos] ? alt_idx : ref_idx),
			_missing[2*pos+1] ? missing_allele : (_genotype[2*pos+1] ? alt_idx : ref_idx));
}

}
}
