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

DataSet::~DataSet(){
	deque<Marker*>::iterator m_it = _markers.begin();
	while(m_it != _markers.end()){
		delete *m_it;
		++m_it;
	}
	_markers.clear();

	deque<Sample*>::iterator s_it = _samples.begin();
	while(s_it != _samples.end()){
		delete *s_it;
		++s_it;
	}
	_samples.clear();

	deque<Family*>::iterator f_it = _families.begin();
	while(f_it != _families.end()){
		delete *f_it;
		++f_it;
	}
	_families.clear();

}

void DataSet::sortMarkers(){
	// Sort the Marker deque (note that operator< is defined for Marker* objects)
	std::sort(_markers.begin(), _markers.end());
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
	Sample* new_samp = Sample::create(famid,id,n_genos);

	_sample_idx_map[new_samp] = _samples.size();
	_samples.push_back(new_samp);
	_sample_map[std::make_pair(famid, id)] = new_samp;

	// Now, add a NaN for everything in the trait map
	map<string, deque<float> >::iterator t_itr = _trait_map.begin();
	while(t_itr != _trait_map.end()){
		(*t_itr).second.push_back(std::numeric_limits<float>::quiet_NaN());
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

Sample* const DataSet::getSample(const string& id) const{
	map<pair<string, string>,Sample*>::const_iterator s_itr = _sample_map.find(std::make_pair(id, id));
	return s_itr == _sample_map.end() ? 0 : (*s_itr).second;
}

Sample* const DataSet::getSample(const string& fid, const string& iid) const{
	map<pair<string, string>,Sample*>::const_iterator s_itr = _sample_map.find(std::make_pair(fid, iid));
	return s_itr == _sample_map.end() ? 0 : (*s_itr).second;
}

Marker* const DataSet::getMarker(const std::string& id) const{
	map<string,Marker*>::const_iterator m_itr = _marker_map.find(id);
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

Family* const DataSet::getFamily(const std::string& id) const{
	map<string,Family*>::const_iterator f_itr = _family_map.find(id);
	return f_itr == _family_map.end() ? 0 : (*f_itr).second;
}

Marker* const DataSet::getMarker(const std::string& chrom, unsigned int loc) const{
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

	map<string, deque<float> >::iterator itr = _trait_map.find(trait);

	//If the mapping isn't found for the given trait, add a completely NaN entry
	if(itr == _trait_map.end()){
		itr = _trait_map.insert(_trait_map.end(), make_pair(trait,
				deque<float> (_samples.size(),
						std::numeric_limits<float>::quiet_NaN())));
	}

	(*itr).second[(*s_itr).second] = val;
	return true;
}

float DataSet::getTrait(const std::string& trait, const Sample* samp) const{
	map<const Sample*, unsigned int>::const_iterator s_itr = _sample_idx_map.find(samp);
	map<string, deque<float> >::const_iterator itr = _trait_map.find(trait);

	if(s_itr == _sample_idx_map.end() || itr == _trait_map.end()){
		return std::numeric_limits<float>::quiet_NaN();
	}

	return (*itr).second[(*s_itr).second];
}

}
}
