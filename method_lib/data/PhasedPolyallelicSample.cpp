/*
 * PhasedPolyallelicSample.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "PhasedPolyallelicSample.h"

#include "Marker.h"

using std::string;
using std::pair;

namespace PLATO{
namespace Data{

PhasedPolyallelicSample::PhasedPolyallelicSample(const string& famid, const string& id, unsigned int n_genos) :
	Sample(famid, id) {
	for(unsigned int i=0; i<2*n_genos; i++){
		_genotype.push_back(missing_allele);
	}
}

void PhasedPolyallelicSample::appendGenotype(unsigned char geno1, unsigned char geno2){
	_genotype.push_back(geno1);
	_genotype.push_back(geno2);
}

void PhasedPolyallelicSample::appendMissingGenotype(){
	_genotype.push_back(missing_allele);
	_genotype.push_back(missing_allele);
}

void PhasedPolyallelicSample::setGenotype(const Marker& m, unsigned char geno1, unsigned char geno2){
	unsigned int pos = m.getIndex();
	_genotype[2*pos] = geno1;
	_genotype[2*pos + 1] = geno2;
}

void PhasedPolyallelicSample::setMissingGenotype(const Marker& m){
	unsigned int pos = m.getIndex();
	_genotype[2*pos] = missing_allele;
	_genotype[2*pos + 1] = missing_allele;
}

bool PhasedPolyallelicSample::isMissing(const Marker& m) const{
	unsigned int pos = m.getIndex();
	return _genotype[2*pos] == missing_allele || _genotype[2*pos + 1] == missing_allele;
}

pair<unsigned char, unsigned char> PhasedPolyallelicSample::getGeno(const Marker& m) const{
	unsigned int pos = m.getIndex();
	return std::make_pair(_genotype[2*pos], _genotype[2*pos + 1]);
}

}
}
