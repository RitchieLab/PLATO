#include "Sample.h"

#include "PhasedBiallelicSample.h"
#include "BiallelicSample.h"
#include "PolyallelicSample.h"
#include "PhasedPolyallelicSample.h"

using std::string;
using std::set;

namespace Methods{

bool Sample::_biallelic = true;
bool Sample::_phased = false;
std::pair<unsigned char, unsigned char> Sample::missing_geno = std::make_pair(-1, -1);

Sample::Sample(const string& famid, const string& id) :
	_famid(famid), _id(id), _mom(NULL), _dad(NULL){
}

Sample* Sample::create(const string& famid, const string& id, unsigned int n_genos){
	if(_phased){
		if(_biallelic){
			return new PhasedBiallelicSample(famid, id, n_genos);
		}else{
			return new PhasedPolyallelicSample(famid, id, n_genos);
		}
	}else{
		if(_biallelic){
			return new BiallelicSample(famid, id, n_genos);
		}else{
			return new PolyallelicSample(famid, id, n_genos);
		}
	}
}

}
