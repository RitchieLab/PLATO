#ifndef METHODS_FAMILY_H
#define METHODS_FAMILY_H

#include <deque>
#include <vector>

namespace Methods{

class Sample;

class Family{
	Family(){_founder.push_back(true);}// default ctor is fine for now

private:
	Family(const Family& o);
	Family& operator=(const Family& o);

public:
	void addFounder(Sample* mem) {addMember(mem, true);}
	void addNonFounder(Sample* mem) {addMember(mem, false);};

	void setEnabled(bool enabled=true){_founder[0] = enabled;}

	bool isEnabled() const {return _founder[0];}


private:
	void addMember(Sample* mem, bool found){_members.push_back(mem); _founder.push_back(found);}

	std::deque<Sample*> _members;
	// 1st bit is if the family is enabled as a whole
	std::vector<bool> _founder;
};

}

#endif
