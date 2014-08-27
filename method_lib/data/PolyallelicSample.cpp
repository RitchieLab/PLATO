/*
 * PolyallelicSample.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "PolyallelicSample.h"

#include "Marker.h"

using std::string;
using std::pair;
using std::make_pair;

namespace PLATO{
namespace Data{

PolyallelicSample::PolyallelicSample(const string& famid, const string& id, unsigned int n_genos) :
	Sample(famid, id) {
	for(unsigned int i=0; i<2*n_genos; i++){
		_genotype.push_back(missing_allele);
	}
}

void PolyallelicSample::appendGenotype(unsigned char geno1, unsigned char geno2){

	// If either is "missing", then both are missing!
	if(geno1 == missing_allele || geno2 == missing_allele){
		appendMissingGenotype();
	}else{
		_genotype.push_back(geno1);
		_genotype.push_back(geno2);
	}
}

void PolyallelicSample::appendMissingGenotype(){
	_genotype.push_back(missing_allele);
	_genotype.push_back(missing_allele);
}

void PolyallelicSample::setGenotype(const Marker& m, unsigned char geno1, unsigned char geno2){

	// If either is "missing", then both are missing!
	if(geno1 == missing_allele || geno2 == missing_allele){
		setMissingGenotype(m);
	}else{
		unsigned int pos = m.getIndex();
		_genotype[2*pos] = geno1;
		_genotype[2*pos + 1] = geno2;
	}
}

void PolyallelicSample::setMissingGenotype(const Marker& m){
	unsigned int pos = m.getIndex();
	_genotype[2*pos] = missing_allele;
	_genotype[2*pos + 1] = missing_allele;
}

bool PolyallelicSample::isMissing(const Marker& m) const{
	unsigned int pos = m.getIndex();
	return _genotype[2*pos] == missing_allele || _genotype[2*pos + 1] == missing_allele;
}

pair<unsigned char, unsigned char> PolyallelicSample::getGeno(const Marker& m) const{
	unsigned int pos = m.getIndex();
	return make_pair(_genotype[2*pos], _genotype[2*pos + 1]);
}

}
}
