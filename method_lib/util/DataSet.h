#ifndef METHODS_DATASET_H
#define METHODS_DATASET_H

#include <vector>
#include <deque>
#include <map>
#include <string>

#include <boost/iterator/iterator_facade.hpp>

namespace Methods{

class Marker;
class Sample;
class Family;

class DataSet{

public:
	~DataSet();

	template <class T>
	class const_iterator : public boost::iterator_facade<const_iterator<T>, T* const, boost::forward_traversal_tag>{

	public:
		const_iterator(std::deque<T*>::const_iterator& itr, std::deque<T*>::const_iterator& end) :
			_itr(itr), _end(end) {}

	private:
		friend class boost::iterator_core_access;

		// Iterate only over enabled samples
		void increment() { while(_itr != _end && !((*(++_itr))->isEnabled()));}
		bool equal(const const_iterator& other) const { return _itr == other._itr;}
		T& dereference() const { return (**_itr);}

		std::deque<T*>::const_iterator _itr;
		const std::deque<T*>::const_iterator _end;
	};
	template <class T>
	class iterator : public boost::iterator_facade<iterator<T>, T*, boost::forward_traversal_tag>{

	public:
		iterator(std::deque<T*>::iterator& itr, std::deque<T*>::iterator& end) :
			_itr(itr), _end(end) {}

	private:
		friend class boost::iterator_core_access;

		// Iterate only over enabled samples
		void increment() { while(_itr != _end && !((*(++_itr))->isEnabled()));}
		bool equal(const const_iterator& other) const { return _itr == other._itr;}
		T& dereference() const { return (**_itr);}

		std::deque<T*>::iterator _itr;
		const std::deque<T*>::iterator _end;
	};

	typedef const_iterator<Sample> const_sample_iterator;
	typedef const_iterator<Marker> const_marker_iterator;
	typedef const_iterator<Family> const_family_iterator;

	typedef iterator<Sample> sample_iterator;
	typedef iterator<Marker> marker_iterator;
	typedef iterator<Family> family_iterator;

	void setBiallelic(bool enabled=true){Sample::_biallelic = enabled;}
	void setPhased(bool enabled=true){Sample::_phased = enabled;}

	const_sample_iterator beginSample() const{
		return const_sample_iterator(_samples.begin(), _samples.end());}
	const_marker_iterator beginMarker() const{
		return const_marker_iterator(_markers.begin(), _markers.end());}
	const_family_iterator beginFamily() const{
		return const_family_iterator(_families.begin(), _families.end());}

	const_sample_iterator endSample() const{
		return const_sample_iterator(_samples.end(), _samples.end());}
	const_marker_iterator endMarker() const{
		return const_marker_iterator(_markers.end(), _markers.end());}
	const_family_iterator endFamily() const{
		return const_family_iterator(_families.end(), _families.end());}

	sample_iterator beginSample() {
		return sample_iterator(_samples.begin(), _samples.end());}
	marker_iterator beginMarker() {
		return marker_iterator(_markers.begin(), _markers.end());}
	family_iterator beginFamily() {
		return family_iterator(_families.begin(), _families.end());}

	sample_iterator endSample() {
		return sample_iterator(_samples.end(), _samples.end());}
	marker_iterator endMarker() {
		return marker_iterator(_markers.end(), _markers.end());}
	family_iterator endFamily() {
		return family_iterator(_families.end(), _families.end());}

	Marker* addMarker(const std::string& chrom, unsigned int loc, const std::string& id);
	Sample* addSample(const std::string& famid, const std::string& id, unsigned int n_genos=0);
	Sample* addSample(const std::string& id, unsigned int n_genos=0);
	Family* addFamily(const std::string& id);

	void sortMarkers();

	Sample* const getSample(const std::string& id) const;
	Marker* const getMaker(const std::string& id) const;
	Family* const getFamily(const std::string& id) const;

	Marker* const getMarker(const std::string& chrom, unsigned int loc) const;

private:
	DataSet(const DataSet& other);
	DataSet& operator=(const DataSet& other);

private:
	std::deque<Marker*> _markers;
	std::deque<Sample*> _samples;
	std::deque<Family*> _families;
	std::deque<std::string> _trait_name;
	std::map<std::string, Marker*> _marker_map;
	std::map<std::pair<std::string, unsigned int>, Marker*> _marker_pos_map;
	std::map<std::string, Sample*> _sample_map;
	std::map<std::string, Family*> _family_map;

	unsigned int _marker_idx;

};

}

#endif
