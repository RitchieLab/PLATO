#ifndef DATA_DATASET_H
#define DATA_DATASET_H

#include <deque>
#include <map>
#include <string>
#include <vector>
#include <algorithm>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/algorithm/string.hpp>

#include "Sample.h"
#include "Marker.h"

namespace PLATO{

namespace Data{

class Family;

class DataSet{

public:
	DataSet() : _marker_idx(0), _biallelic(false), _phased(false) {}
	~DataSet();

	template <class T>
	class const_vec_iterator : public boost::iterator_facade<const_vec_iterator<T>, const T* const, boost::forward_traversal_tag>{

	public:
		const_vec_iterator(
				const typename std::vector<T*>::const_iterator& itr,
				const typename std::vector<T*>::const_iterator& end) :
			_itr(itr), _end(end) {}

	private:
		friend class boost::iterator_core_access;

		// Iterate only over enabled samples
		void increment() { while(++_itr != _end && !((*(_itr))->isEnabled()));}
		bool equal(const const_vec_iterator& other) const { return _itr == other._itr;}
		const T* const & dereference() const { return _tmp = *_itr;}

		mutable T const* _tmp;
		typename std::vector<T*>::const_iterator _itr;
		typename std::vector<T*>::const_iterator _end;
	};

	template <class T>
	class vec_iterator : public boost::iterator_facade<vec_iterator<T>, T*, boost::forward_traversal_tag>{

	public:
		vec_iterator(const typename std::vector<T*>::iterator& itr,
				const typename std::vector<T*>::iterator& end) :
			_itr(itr), _end(end) {}

		// Make sure I can convert an iterator<T> to a const_iterator<T>
		operator const_vec_iterator<T>() const {
			return const_vec_iterator<T>(_itr, _end);
		}

	private:
		friend class boost::iterator_core_access;

		// Iterate only over enabled samples
		void increment() { while(  !(++_itr == _end) && !((*(_itr))->isEnabled()));}
		bool equal(const vec_iterator& other) const { return _itr == other._itr;}
		T* & dereference() const { return (*_itr);}

		typename std::vector<T*>::iterator _itr;
		typename std::vector<T*>::iterator _end;
	};

	template <class T>
	class const_map_iterator: public boost::iterator_facade<
			const_map_iterator<T>, const std::string,
			boost::forward_traversal_tag> {

	public:
		const_map_iterator(const typename std::map<std::string, std::pair<bool, T> >::const_iterator& itr,
				const typename std::map<std::string, std::pair<bool, T> >::const_iterator& end) :
			_itr(itr), _end(end) {}

	private:
		friend class boost::iterator_core_access;

		void increment() { while(_itr != _end && ++_itr != _end && !(*_itr).second.first);}
		bool equal(const const_map_iterator& other) const { return _itr == other._itr; }
		const std::string& dereference() const { return (*_itr).first;}

		typename std::map<std::string, std::pair<bool, T> >::const_iterator _itr;
		typename std::map<std::string, std::pair<bool, T> >::const_iterator _end;
	};

	template <class T>
	class map_iterator: public boost::iterator_facade<
			map_iterator<T>, const std::string,
			boost::forward_traversal_tag> {

	public:
		map_iterator(const typename std::map<std::string, std::pair<bool, T> >::iterator& itr,
				const typename std::map<std::string, std::pair<bool, T> >::iterator& end) :
			_itr(itr), _end(end) {}

		operator const_map_iterator<T>() const{
			return const_map_iterator<T>(_itr, _end);
		}
	private:
		friend class boost::iterator_core_access;

		void increment() { while(_itr != _end && !(*_itr).second.first) {++_itr;}}
		bool equal(const map_iterator& other) const { return _itr == other._itr; }
		const std::string& dereference() const { return (*_itr).first;}

		typename std::map<std::string, std::pair<bool, T> >::iterator _itr;
		typename std::map<std::string, std::pair<bool, T> >::iterator _end;
	};

	typedef const_vec_iterator<Sample> const_sample_iterator;
	typedef const_vec_iterator<Marker> const_marker_iterator;
	typedef const_vec_iterator<Family> const_family_iterator;
	typedef const_map_iterator<std::deque<float> > const_trait_iterator;
	typedef const_map_iterator<std::deque<unsigned char> > const_categorical_iterator;

	typedef vec_iterator<Sample> sample_iterator;
	typedef vec_iterator<Marker> marker_iterator;
	typedef vec_iterator<Family> family_iterator;
	typedef map_iterator<std::deque<float> > trait_iterator;
	typedef map_iterator<std::deque<unsigned char> > categorical_iterator;

	void setBiallelic(bool enabled=true){_biallelic = enabled;}
	void setPhased(bool enabled=true){_phased = enabled;}
	bool isPhased() const {return _phased;}
	bool isBiallelic() const {return _biallelic;}


	const_sample_iterator beginSample() const{
		return const_sample_iterator(_samples.begin(), _samples.end());}
	const_marker_iterator beginMarker() const{
		return const_marker_iterator(_markers.begin(), _markers.end());}
	const_family_iterator beginFamily() const{
		return const_family_iterator(_families.begin(), _families.end());}
	const_trait_iterator beginTrait() const{
		return const_trait_iterator(_trait_map.begin(), _trait_map.end());}
	const_categorical_iterator beginCategorical() const{
		return const_categorical_iterator(_categorical_map.begin(), _categorical_map.end());}


	const_sample_iterator endSample() const{
		return const_sample_iterator(_samples.end(), _samples.end());}
	const_marker_iterator endMarker() const{
		return const_marker_iterator(_markers.end(), _markers.end());}
	const_family_iterator endFamily() const{
		return const_family_iterator(_families.end(), _families.end());}
	const_trait_iterator endTrait() const{
		return const_trait_iterator(_trait_map.end(), _trait_map.end());}
	const_categorical_iterator endCategorical() const{
		return const_categorical_iterator(_categorical_map.end(), _categorical_map.end());}


	sample_iterator beginSample() {
		return sample_iterator(_samples.begin(), _samples.end());}
	marker_iterator beginMarker() {
		return marker_iterator(_markers.begin(), _markers.end());}
	family_iterator beginFamily() {
		return family_iterator(_families.begin(), _families.end());}
	trait_iterator beginTrait() {
		return trait_iterator(_trait_map.begin(), _trait_map.end());}
	categorical_iterator beginCategorical() {
		return categorical_iterator(_categorical_map.begin(), _categorical_map.end());}

	sample_iterator endSample() {
		return sample_iterator(_samples.end(), _samples.end());}
	marker_iterator endMarker() {
		return marker_iterator(_markers.end(), _markers.end());}
	family_iterator endFamily() {
		return family_iterator(_families.end(), _families.end());}
	trait_iterator endTrait() {
		return trait_iterator(_trait_map.end(), _trait_map.end());}
	categorical_iterator endCategorical(){
		return categorical_iterator(_categorical_map.end(), _categorical_map.end());}

	Marker* addMarker(const std::string& chrom, unsigned int loc, const std::string& id);
	Sample* addSample(const std::string& famid, const std::string& id, unsigned int n_genos=0);
	Sample* addSample(const std::string& id, unsigned int n_genos=0);
	Family* addFamily(const std::string& id);
	bool addTrait(const std::string& trait, const Sample* samp, float val);
	float getTrait(const std::string& trait, const Sample* samp) const;
	bool setTraitEnabled(const std::string& trait, bool isEnabled=true);
	bool isTrait(const std::string& trait) const { return _trait_map.find(trait) != _trait_map.end();}

	bool addCategorical(const std::string& trait, const Sample* samp, unsigned char val);
	unsigned char getCategorical(const std::string& trait, const Sample* samp) const;
	bool setCategoricalEnabled(const std::string& trait, bool isEnabled=true);
	unsigned char numCategories(const std::string& trait) const;


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

	void sortMarkers();
	void sortSamples();

private:
	DataSet(const DataSet& other);
	DataSet& operator=(const DataSet& other);

private:
	template<typename T>
	struct ci_less:std::binary_function<T,T,bool>
	  { bool operator() (const T& s1,const T& s2) const { return boost::ilexicographical_compare(s1,s2); }};

	typedef std::map<std::string, Marker*, ci_less<std::string> > MarkerIDMap;

	std::vector<Marker*> _markers;
	std::vector<Sample*> _samples;
	std::vector<Family*> _families;
	MarkerIDMap _marker_map;
	std::map<std::pair<unsigned short, unsigned int>, Marker*> _marker_pos_map;
	std::map<std::pair<std::string, std::string>, Sample*> _sample_map;
	std::map<std::string, Family*> _family_map;

	std::map<const Sample*, unsigned int> _sample_idx_map;
	std::map<std::string, std::pair<bool, std::deque<float> > > _trait_map;
	std::map<std::string, std::pair<bool, std::deque<unsigned char> > > _categorical_map;
	mutable std::map<std::string, unsigned char> _categorical_sz;

	unsigned int _marker_idx;

	bool _biallelic;
	bool _phased;

public:
	static const unsigned char MISSING_CATEGORICAL;
};

}
}
#endif
