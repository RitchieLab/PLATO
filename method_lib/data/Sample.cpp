#include "Sample.h"

#include <limits>

#include "Marker.h"

#include "PhasedBiallelicSample.h"
#include "BiallelicSample.h"
#include "PolyallelicSample.h"
#include "PhasedPolyallelicSample.h"

using std::string;
using std::set;

namespace PLATO{
namespace Data{

bool Sample::_biallelic = true;
bool Sample::_phased = false;
const unsigned char Sample::missing_allele = static_cast<unsigned char>(-1);
const std::pair<unsigned char, unsigned char> Sample::missing_geno = std::make_pair(Sample::missing_allele, Sample::missing_allele);

Sample::Sample(const string& famid, const string& id) :
	_famid(famid), _id(id), _mom(NULL), _dad(NULL),
	_pheno(std::numeric_limits<float>::quiet_NaN()),
	_sex_known(false), _affected_known(false), _founder(true), _enabled(true){
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

unsigned char Sample::getAdditiveGeno(const Marker& m) const{
	unsigned char geno = missing_allele;

	std::pair<unsigned char, unsigned char> allele_pair = getGeno(m);

	if(allele_pair.first != missing_allele && allele_pair.second != missing_allele){
		geno = (allele_pair.first != m.getRefIdx()) + (allele_pair.second != m.getRefIdx());
	}

	return geno;
}

}
}
