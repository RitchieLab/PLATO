#include "Sample.h"

using std::string;
using std::set;

namespace Methods{

Sample::Sample(const string& famid, const string& id, unsigned int n_genos) :
	_famid(famid), _id(id), _mom(NULL), _dad(NULL){

	_genotype.resize(N_FLAGS+2*n_genos);

}

Sample::Sample(const string& id, unsigned int n_genos) :
	_famid(id), _id(id), _mom(NULL), _dad(NULL){

	_genotype.resize(N_FLAGS+2*n_genos);

}

void Sample::appendGenotype(bool geno1, bool geno2){
	_genotype.push_back(geno1);
	_genotype.push_back(geno2);
}

}
