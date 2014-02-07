#include "Marker.h"

#include "DataSet.h"
#include "Sample.h"

#include <algorithm>

using std::string;
using std::find;
using std::vector;

using PLATO::Utility::InputManager;

namespace PLATO{
namespace Data{

string Marker::_missing_allele = "0";

Marker::Marker(const string& chr, unsigned int loc, const string& id, unsigned int idx) :
		_chr(InputManager::chrStringToInt(chr) | ENABLED_MASK), _loc(loc), _idx(idx), _id(id){
	_ref_idx = _alt_idx = static_cast<unsigned char>(-1);
}

unsigned int Marker::addAllele(const string& allele){
	vector<string>::const_iterator allele_itr = find(_alleles.begin(), _alleles.end(), allele);
	unsigned int idx = allele_itr - _alleles.begin();
	if(allele_itr == _alleles.end()){
		_alleles.push_back(allele);
		if(_alleles.size() == 1){
			_ref_idx = 0;
		}else if(_alleles.size() == 2){
			_alt_idx = 1;
		}
	}
	return idx;
}

bool Marker::setRefAllele(const string& allele){
	bool to_ret = false;
	vector<string>::const_iterator itr = find(_alleles.begin(), _alleles.end(), allele);
	if(itr != _alleles.end()){
		unsigned int old_ref = _ref_idx;
		_ref_idx = itr - _alleles.begin();
		if(_ref_idx == _alt_idx){
			_alt_idx = old_ref;
		}
		to_ret = true;
	}
	return to_ret;
}

bool Marker::setAltAllele(const string& allele){
	bool to_ret = false;
	vector<string>::const_iterator itr = find(_alleles.begin(), _alleles.end(), allele);
	if(_alleles.size() > 1 && itr != _alleles.end()){
		unsigned int old_alt = _alt_idx;
		_alt_idx = itr - _alleles.begin();
		if(_ref_idx == _alt_idx){
			_ref_idx = old_alt;
		}
		to_ret = true;
	}
	return to_ret;
}

float Marker::calcMAF(const DataSet& ds) const{
	DataSet::const_sample_iterator si = ds.beginSample();
	int n_sample = 0;
	int n_allele = 0;
	while(si != ds.endSample()){
		n_allele += (*si)->getAdditiveGeno(*this);
		++n_sample;
		++si;
	}

	float maf = n_allele / static_cast<float>(2 * n_sample);

	return std::min(maf, 1-maf);
}

}
}
