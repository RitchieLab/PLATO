#ifndef DATA_DATASET_H
#define DATA_DATASET_H

#include <deque>
#include <map>
#include <string>
#include <vector>

#include <boost/iterator/iterator_facade.hpp>

#include "Sample.h"

namespace PLATO{

namespace Data{

class Marker;
class Family;

class DataSet{

public:
	DataSet() : _marker_idx(0) {}
	~DataSet();

	template <class T>
	class const_iterator : public boost::iterator_facade<const_iterator<T>, const T* const, boost::forward_traversal_tag>{

	public:
		const_iterator(
				const typename std::vector<T*>::const_iterator& itr,
				const typename std::vector<T*>::const_iterator& end) :
			_itr(itr), _end(end) {}

	private:
		friend class boost::iterator_core_access;

		// Iterate only over enabled samples
		void increment() { while(++_itr != _end && !((*(_itr))->isEnabled()));}
		bool equal(const const_iterator& other) const { return _itr == other._itr;}
		const T* const & dereference() const { return _tmp = *_itr;}

		mutable T const* _tmp;
		typename std::vector<T*>::const_iterator _itr;
		typename std::vector<T*>::const_iterator _end;
	};

	template <class T>
	class iterator : public boost::iterator_facade<iterator<T>, T*, boost::forward_traversal_tag>{

	public:
		iterator(const typename std::vector<T*>::iterator& itr,
				const typename std::vector<T*>::iterator& end) :
			_itr(itr), _end(end) {}

		// Make sure I can convert an iterator<T> to a const_iterator<T>
		operator const_iterator<T>() {
			return const_iterator<T>(_itr, _end);
		}

	private:
		friend class boost::iterator_core_access;

		// Iterate only over enabled samples
		void increment() { while(  !(++_itr == _end) && !((*(_itr))->isEnabled()));}
		bool equal(const iterator& other) const { return _itr == other._itr;}
		T* & dereference() const { return (*_itr);}

		typename std::vector<T*>::iterator _itr;
		typename std::vector<T*>::iterator _end;
	};

	class const_trait_iterator: public boost::iterator_facade<
			const_trait_iterator, const std::string&,
			boost::forward_traversal_tag> {

	public:
		const_trait_iterator(std::map<std::string, std::pair<bool, std::deque<float> > >::const_iterator itr,
				std::map<std::string, std::pair<bool, std::deque<float> > >::const_iterator end) :
			_itr(itr), _end(end) {}

	private:
		friend class boost::iterator_core_access;

		void increment() { while(_itr != _end && !(*_itr).second.first) {++_itr;}}
		bool equal(const const_trait_iterator& other) const { return _itr == other._itr; }
		const std::string& dereference() const { return (*_itr).first;}

		std::map<std::string, std::pair<bool, std::deque<float> > >::const_iterator _itr;
		std::map<std::string, std::pair<bool, std::deque<float> > >::const_iterator _end;
	};

	class trait_iterator: public boost::iterator_facade<
			trait_iterator, const std::string&,
			boost::forward_traversal_tag> {

	public:
		trait_iterator(std::map<std::string, std::pair<bool, std::deque<float> > >::iterator itr,
				std::map<std::string, std::pair<bool, std::deque<float> > >::iterator end) :
			_itr(itr), _end(end) {}

		operator const_trait_iterator() const{
			return const_trait_iterator(_itr, _end);
		}
	private:
		friend class boost::iterator_core_access;

		void increment() { while(_itr != _end && !(*_itr).second.first) {++_itr;}}
		bool equal(const trait_iterator& other) const { return _itr == other._itr; }
		const std::string& dereference() const { return (*_itr).first;}

		std::map<std::string, std::pair<bool, std::deque<float> > >::iterator _itr;
		std::map<std::string, std::pair<bool, std::deque<float> > >::iterator _end;
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
	const_trait_iterator beginTrait() const{
		return const_trait_iterator(_trait_map.begin(), _trait_map.end());}

	const_sample_iterator endSample() const{
		return const_sample_iterator(_samples.end(), _samples.end());}
	const_marker_iterator endMarker() const{
		return const_marker_iterator(_markers.end(), _markers.end());}
	const_family_iterator endFamily() const{
		return const_family_iterator(_families.end(), _families.end());}
	const_trait_iterator endTrait() const{
		return const_trait_iterator(_trait_map.end(), _trait_map.end());}


	sample_iterator beginSample() {
		return sample_iterator(_samples.begin(), _samples.end());}
	marker_iterator beginMarker() {
		return marker_iterator(_markers.begin(), _markers.end());}
	family_iterator beginFamily() {
		return family_iterator(_families.begin(), _families.end());}
	trait_iterator beginTrait() {
		return trait_iterator(_trait_map.begin(), _trait_map.end());}

	sample_iterator endSample() {
		return sample_iterator(_samples.end(), _samples.end());}
	marker_iterator endMarker() {
		return marker_iterator(_markers.end(), _markers.end());}
	family_iterator endFamily() {
		return family_iterator(_families.end(), _families.end());}
	trait_iterator endTrait() {
		return trait_iterator(_trait_map.end(), _trait_map.end());}

	Marker* addMarker(const std::string& chrom, unsigned int loc, const std::string& id);
	Sample* addSample(const std::string& famid, const std::string& id, unsigned int n_genos=0);
	Sample* addSample(const std::string& id, unsigned int n_genos=0);
	Family* addFamily(const std::string& id);
	bool addTrait(const std::string& trait, const Sample* samp, float val);
	float getTrait(const std::string& trait, const Sample* samp) const;
	bool setTraitEnabled(const std::string& trait, bool isEnabled=true);
	bool isTrait(const std::string& trait) const { return _trait_map.find(trait) != _trait_map.end();}

	void sortMarkers();

	const Sample* getSample(const std::string& id) const;
	const Sample* getSample(const std::string& fid, const std::string& id) const;
	const Marker* getMarker(const std::string& id) const;
	const Family* getFamily(const std::string& id) const;
	const Marker* getMarker(const std::string& chrom, unsigned int loc) const;

	Sample* getSample(const std::string& id);
	Sample* getSample(const std::string& fid, const std::string& id);
	Marker* getMarker(const std::string& id);
	Family* getFamily(const std::string& id);
	Marker* getMarker(const std::string& chrom, unsigned int loc);

	unsigned int num_loci() const {return _markers.size();}
	unsigned int num_pedigrees() const {return _families.size();}
	unsigned int num_inds() const {return _samples.size();}

private:
	DataSet(const DataSet& other);
	DataSet& operator=(const DataSet& other);

private:
	std::vector<Marker*> _markers;
	std::vector<Sample*> _samples;
	std::vector<Family*> _families;
	std::map<std::string, Marker*> _marker_map;
	std::map<std::pair<unsigned short, unsigned int>, Marker*> _marker_pos_map;
	std::map<std::pair<std::string, std::string>, Sample*> _sample_map;
	std::map<std::string, Family*> _family_map;

	std::map<const Sample*, unsigned int> _sample_idx_map;
	std::map<std::string, std::pair<bool, std::deque<float> > > _trait_map;

	unsigned int _marker_idx;

};

}
}
#endif
