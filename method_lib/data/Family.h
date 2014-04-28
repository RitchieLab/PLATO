#ifndef DATA_FAMILY_H
#define DATA_FAMILY_H

#include <deque>
#include <vector>
#include <string>
#include <algorithm>

#include <boost/iterator/iterator_facade.hpp>

namespace PLATO{
namespace Data{

class Sample;

class Family{
public:

	class const_member_iterator : public boost::iterator_facade<const_member_iterator, const Sample*, boost::forward_traversal_tag>{
	public:
		const_member_iterator(std::deque<Sample*>::const_iterator si,
				        const std::deque<Sample*>& samples,
				        const std::vector<bool>& enabled, bool all=true)
			: _si(si), _samples(samples), _enabled(enabled), _all(all) {
			while(_si != _samples.end() && (_all || !_enabled[_si - _samples.begin() + 1])){
				++_si;
			}

		}

		// we need this b/c of the reference items in the class object... boo!
		const_member_iterator& operator=(const_member_iterator other){
			std::swap(*this, other);
			return *this;
		}

	private:
		friend class boost::iterator_core_access;

		void increment(){
			// NOTE: the increment is in the
			while(_si != _samples.end() && ++_si != _samples.end()
				&& (_all || !_enabled[_si - _samples.begin() + 1]));
		}
		bool equal(const const_member_iterator& o) const {return _si == o._si;}

		const Sample*& dereference() const { return _tmp = (*_si);}

		mutable const Sample* _tmp;
		std::deque<Sample*>::const_iterator _si;
		const std::deque<Sample*>& _samples;
		const std::vector<bool>& _enabled;
		bool _all;
	};

public:
	Family(const std::string& fid) : _fam_id(fid) {_founder.push_back(true);}

private:
	Family(const Family& o);
	Family& operator=(const Family& o);

public:

	void addFounder(Sample* mem) {addMember(mem, true);}
	void addNonFounder(Sample* mem) {addMember(mem, false);};

	void setEnabled(bool enabled=true){_founder[0] = enabled;}

	bool isEnabled() const {return _founder[0];}
	const std::string& getFamID() const {return _fam_id;}

	const_member_iterator beginMembers(bool nonfounders=true) const{
		return const_member_iterator(_members.begin(), _members, _founder, nonfounders);
	}

	const_member_iterator endMembers() const{
		return const_member_iterator(_members.end(), _members, _founder);
	}



private:
	void addMember(Sample* mem, bool found){_members.push_back(mem); _founder.push_back(found);}

	std::deque<Sample*> _members;
	// 1st bit is if the family is enabled as a whole
	std::vector<bool> _founder;

	std::string _fam_id;

};

}
}
#endif
