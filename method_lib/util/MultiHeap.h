/*
 * MultiHeap.h
 *
 *  Created on: Feb 13, 2015
 *      Author: jrw32
 */

#ifndef UTILITY_MULTIHEAP_H
#define UTILITY_MULTIHEAP_H

#include <vector>
#include <queue>
#include <utility>

namespace PLATO{
namespace Utility{

template <class T, class Cont, class Cmp>
class MultiHeap {
private:

	typedef std::pair<T, size_t> val_idx;

	struct pair_cmp{
		bool operator() (const std::pair<T, size_t>& x, const std::pair<T, size_t>& y){
			return c(x.first, y.first) || (!c(y.first, x.first) && x.second < y.second);
		}
		Cmp c;
	};

public:
	MultiHeap() {}
	~MultiHeap() {}

	void push(const T& val, size_t idx){
		while(!_superheap.empty()){
			_superheap.pop();
		}
		// If we don't have enough, construct empties until we do!
		if(idx >= _array_of_heaps.size()){
			_array_of_heaps.insert(_array_of_heaps.end(),
					idx + 1 - _array_of_heaps.size(), std::priority_queue<T, Cont, Cmp>());
		}
		_array_of_heaps[idx].push(val);
	}
	void pop(){
		if(_superheap.empty()){
			construct_superheap();
		}
		val_idx top_elt = _superheap.top();
		_superheap.pop();
		_array_of_heaps[top_elt.second].pop();
		if(!_array_of_heaps[top_elt.second].empty()){
			_superheap.push(std::make_pair(_array_of_heaps[top_elt.second].top(), top_elt.second));
		}
	}
	const T& top() const{
		if(_superheap.empty()){
			construct_superheap();
		}
		return _superheap.top().first;
	}
	size_t size() const{
		size_t retv = 0;
		for(size_t i=0; i<_array_of_heaps.size(); i++){
			retv += _array_of_heaps[i].size();
		}
		return retv;
	}

	// Added to get "rank" permuted value
	size_t numLess(size_t heapSz) const{
		size_t retv = 0;
		for(size_t i=0; i<_array_of_heaps.size(); i++){
			retv += (_array_of_heaps[i].size() < heapSz);
		}
		return retv;
	}

private:

	void construct_superheap() const{
		for(size_t i=0; i<_array_of_heaps.size(); i++){
			if(!_array_of_heaps[i].empty()){
				_superheap.push(std::make_pair(_array_of_heaps[i].top(), i));
			}
		}
	}

	mutable std::priority_queue<val_idx, std::vector<val_idx>, pair_cmp > _superheap;
	std::vector<std::priority_queue<T, Cont, Cmp> > _array_of_heaps;

};


}
}

#endif /* MULTIHEAP_H_ */
