#include "DataSet.h"

#include <algorithm>

#include "InputManager.h"

#include "Marker.h"
#include "Sample.h"
#include "Family.h"

using std::deque;
using std::string;
using std::map;
using std::pair;

namespace Methods{

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

	_samples.push_back(new_samp);
	_sample_map[std::make_pair(famid, id)] = new_samp;

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

Marker* const DataSet::getMaker(const std::string& id) const{
	map<string,Marker*>::const_iterator m_itr = _marker_map.find(id);
	return m_itr == _marker_map.end() ? 0 : (*m_itr).second;
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
}
