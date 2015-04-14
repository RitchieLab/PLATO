#include "DataSet.h"

#include <algorithm>
#include <limits>
#include <vector>
#include <sstream>

#include <boost/algorithm/string.hpp>

#include "util/InputManager.h"

#include "Marker.h"
#include "Sample.h"
#include "Family.h"

using std::vector;
using std::deque;
using std::string;
using std::map;
using std::pair;

using PLATO::Utility::InputManager;

namespace PLATO{
namespace Data{

const unsigned char DataSet::MISSING_CATEGORICAL = static_cast<unsigned char>(-1);

DataSet::~DataSet(){
	vector<Marker*>::iterator m_it = _markers.begin();
	while(m_it != _markers.end()){
		delete *m_it;
		++m_it;
	}
	_markers.clear();

	vector<Sample*>::iterator s_it = _samples.begin();
	while(s_it != _samples.end()){
		delete *s_it;
		++s_it;
	}
	_samples.clear();

	vector<Family*>::iterator f_it = _families.begin();
	while(f_it != _families.end()){
		delete *f_it;
		++f_it;
	}
	_families.clear();

}

Marker* DataSet::addMarker(const std::string& chrom, unsigned int loc, const std::string& id){
	// maybe do some sanity checks here...

	Marker* new_marker = new Marker(chrom, loc, id, _marker_idx++);

	_markers.push_back(new_marker);

	if(id != "." && id != ""){
		_marker_map[id] = new_marker;
	}

	_marker_pos_map[std::make_pair(InputManager::chrStringToInt(chrom), loc)] = new_marker;
	return new_marker;
}

Sample* DataSet::addSample(const std::string& famid, const std::string& id, unsigned int n_genos){
	Sample* new_samp = Sample::create(*this,famid,id,n_genos);

	_sample_idx_map[new_samp] = _samples.size();
	_samples.push_back(new_samp);
	_sample_map[std::make_pair(famid, id)] = new_samp;

	// Now, add a NaN for everything in the trait map
	map<string, pair<bool, deque<float> > >::iterator t_itr = _trait_map.begin();
	while(t_itr != _trait_map.end()){
		(*t_itr).second.second.push_back(std::numeric_limits<float>::quiet_NaN());
		++t_itr;
	}


	return new_samp;
}

Sample* DataSet::addSample(const std::string& id, unsigned int n_genos){
	return addSample(id, id, n_genos);
}

Family* DataSet::addFamily(const std::string& id){
	Family* new_fam = new Family(id);

	_families.push_back(new_fam);
	_family_map[id] = new_fam;

	return new_fam;
}

Sample* DataSet::getSample(const string& id){
	map<pair<string, string>,Sample*>::const_iterator s_itr = _sample_map.find(std::make_pair(id, id));
	return s_itr == _sample_map.end() ? 0 : (*s_itr).second;
}

const Sample* DataSet::getSample(const string& id) const{
	map<pair<string, string>,Sample*>::const_iterator s_itr = _sample_map.find(std::make_pair(id, id));
	return s_itr == _sample_map.end() ? 0 : (*s_itr).second;
}

Sample* DataSet::getSample(const string& fid, const string& iid){
	map<pair<string, string>,Sample*>::const_iterator s_itr = _sample_map.find(std::make_pair(fid, iid));
	return s_itr == _sample_map.end() ? 0 : (*s_itr).second;
}

const Sample* DataSet::getSample(const string& fid, const string& iid) const{
	map<pair<string, string>,Sample*>::const_iterator s_itr = _sample_map.find(std::make_pair(fid, iid));
	return s_itr == _sample_map.end() ? 0 : (*s_itr).second;
}

Marker* DataSet::getMarker(const std::string& id){
	MarkerIDMap::const_iterator m_itr = _marker_map.find(id);
	Marker* m = 0;
	if(m_itr != _marker_map.end()){
		m = (*m_itr).second;
	}else{
		// try to find with chr:pos instead!
		vector<string> chr_pos;
		boost::algorithm::split(chr_pos, id, boost::is_any_of(":"), boost::token_compress_off);
		if(chr_pos.size() == 2){
			unsigned int pos;
			std::stringstream pos_ss(chr_pos[1]);
			pos_ss >> pos;
			if(pos_ss.eof()){
				m = getMarker(chr_pos[0], pos);
			}
		}
	}
	return m;
}

const Marker* DataSet::getMarker(const std::string& id) const{
	map<string,Marker*>::const_iterator m_itr = _marker_map.find(id);
	const Marker* m = 0;
	if(m_itr != _marker_map.end()){
		m = (*m_itr).second;
	}else{
		// try to find with chr:pos instead!
		vector<string> chr_pos;
		boost::algorithm::split(chr_pos, id, boost::is_any_of(":"), boost::token_compress_off);
		if(chr_pos.size() == 2){
			unsigned int pos;
			std::stringstream pos_ss(chr_pos[1]);
			pos_ss >> pos;
			if(pos_ss.eof()){
				m = getMarker(chr_pos[0], pos);
			}
		}
	}
	return m;
}

Family* DataSet::getFamily(const std::string& id){
	map<string,Family*>::const_iterator f_itr = _family_map.find(id);
	return f_itr == _family_map.end() ? 0 : (*f_itr).second;
}

const Family* DataSet::getFamily(const std::string& id) const{
	map<string,Family*>::const_iterator f_itr = _family_map.find(id);
	return f_itr == _family_map.end() ? 0 : (*f_itr).second;
}

Marker* DataSet::getMarker(const std::string& chrom, unsigned int loc){
	map<pair<unsigned short, unsigned int>,Marker*>::const_iterator m_itr =
			_marker_pos_map.find(std::make_pair(InputManager::chrStringToInt(chrom), loc));
	return m_itr == _marker_pos_map.end() ? 0 : (*m_itr).second;
}

const Marker* DataSet::getMarker(const std::string& chrom, unsigned int loc) const{
	map<pair<unsigned short, unsigned int>,Marker*>::const_iterator m_itr =
			_marker_pos_map.find(std::make_pair(InputManager::chrStringToInt(chrom), loc));
	return m_itr == _marker_pos_map.end() ? 0 : (*m_itr).second;
}

bool DataSet::addTrait(const std::string& trait, const Sample* samp, float val){
	map<const Sample*, unsigned int>::const_iterator s_itr = _sample_idx_map.find(samp);

	// Could not find the position of the given sample - something is bad here!
	if(s_itr == _sample_idx_map.end()){
		return false;
	}

	map<string, pair<bool, deque<float> > >::iterator itr = _trait_map.find(trait);

	//If the mapping isn't found for the given trait, add a completely NaN entry
	if(itr == _trait_map.end()){
		itr = _trait_map.insert(_trait_map.end(), make_pair(trait, make_pair(true,
				deque<float> (_samples.size(),
						std::numeric_limits<float>::quiet_NaN()))));
	}

	(*itr).second.second[(*s_itr).second] = val;
	return true;
}

float DataSet::getTrait(const std::string& trait, const Sample* samp) const{
	map<const Sample*, unsigned int>::const_iterator s_itr = _sample_idx_map.find(samp);
	map<string, pair<bool, deque<float> > >::const_iterator itr = _trait_map.find(trait);

	if(s_itr == _sample_idx_map.end() || itr == _trait_map.end()){
		return std::numeric_limits<float>::quiet_NaN();
	}

	return (*itr).second.second[(*s_itr).second];
}

bool DataSet::setTraitEnabled(const std::string& trait, bool isEnabled){
	map<string, pair<bool, deque<float> > >::iterator m_itr = _trait_map.find(trait);

	if(m_itr != _trait_map.end()){
		(*m_itr).second.first = isEnabled;
		return true;
	}
	// return false if the given trait is not found
	return false;
}


bool DataSet::addCategorical(const std::string& trait, const Sample* samp, unsigned char val){
	map<const Sample*, unsigned int>::const_iterator s_itr = _sample_idx_map.find(samp);

	// Could not find the position of the given sample - something is bad here!
	if(s_itr == _sample_idx_map.end()){
		return false;
	}

	map<string, pair<bool, deque<unsigned char> > >::iterator itr = _categorical_map.find(trait);

	if(val != MISSING_CATEGORICAL){
		// make sure to invalidate the cache if needed
		map<string, unsigned char>::iterator cache_itr = _categorical_sz.find(trait);
		if(cache_itr == _categorical_sz.end()){
			_categorical_sz[trait] = val;
		} else if(val > (*cache_itr).second){
			(*cache_itr).second = val;
		}
	}

	//If the mapping isn't found for the given trait, add a completely NaN entry
	if(itr == _categorical_map.end()){
		itr = _categorical_map.insert(_categorical_map.end(), make_pair(trait, make_pair(true,
				deque<unsigned char> (_samples.size(), MISSING_CATEGORICAL))));
	}

	(*itr).second.second[(*s_itr).second] = val;
	return true;
}

unsigned char DataSet::getCategorical(const std::string& trait, const Sample* samp) const{
	map<const Sample*, unsigned int>::const_iterator s_itr = _sample_idx_map.find(samp);
	map<string, pair<bool, deque<unsigned char> > >::const_iterator itr = _categorical_map.find(trait);

	if(s_itr == _sample_idx_map.end() || itr == _categorical_map.end()){
		return MISSING_CATEGORICAL;
	}

	return (*itr).second.second[(*s_itr).second];
}

bool DataSet::setCategoricalEnabled(const std::string& trait, bool isEnabled){
	map<string, pair<bool, deque<unsigned char> > >::iterator m_itr = _categorical_map.find(trait);

	if(m_itr != _categorical_map.end()){
		(*m_itr).second.first = isEnabled;
		return true;
	}
	// return false if the given trait is not found
	return false;
}

unsigned char DataSet::numCategories(const std::string& trait) const{
	map<string, unsigned char>::iterator c_itr = _categorical_sz.find(trait);
	map<string, pair<bool, deque<unsigned char> > >::const_iterator m_itr = _categorical_map.find(trait);

	if(c_itr != _categorical_sz.end()){
		return (*c_itr).second;
	} else if(m_itr == _categorical_map.end()) {
		return MISSING_CATEGORICAL;
	} else {
		// WE SHOULD NOT BE HERE!!!

		unsigned char n_elements = 0;
		for(std::deque<unsigned char>::const_iterator d_itr = (*m_itr).second.second.begin();
		    d_itr != (*m_itr).second.second.end(); d_itr++){
			n_elements = std::max(n_elements,
					static_cast<unsigned char>((*d_itr) * (*d_itr != MISSING_CATEGORICAL)));
		}
		_categorical_sz[trait] = n_elements;
		return n_elements;
	}
}

void DataSet::sortMarkers() {
	 std::sort(_markers.begin(), _markers.end(), std::less<Marker*>());
}

void DataSet::sortSamples() {
	 std::sort(_samples.begin(), _samples.end(), std::less<Sample*>());
}

}
}
