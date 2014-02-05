#ifndef DATA_FAMILY_H
#define DATA_FAMILY_H

#include <deque>
#include <vector>
#include <string>

namespace PLATO{
namespace Data{

class Sample;

class Family{

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
