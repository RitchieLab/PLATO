#ifndef DATA_SAMPLE_H
#define DATA_SAMPLE_H

#include <deque>
#include <string>
#include <set>

#include <utility>

namespace PLATO{
namespace Data{

class DataSet;
class Marker;

class Sample{

public:
	virtual ~Sample(){}

protected:

	Sample(const std::string& famid, const std::string& id);

public:
	static Sample* create(const std::string& famid, const std::string& id, unsigned int n_genos=0);
	//static Sample* create(const std::string& id, unsigned int n_genos=0){return create(id, id, n_genos);}

private:
	// No copying or assignment!!
	Sample(const Sample&);
	Sample& operator=(const Sample&);

public:
	bool operator<(const Sample& o) const{
		return std::make_pair(_id, _famid) < std::make_pair(o._id, o._famid);
	}

	virtual void appendGenotype(unsigned char geno1, unsigned char geno2) = 0;
	virtual void appendMissingGenotype() = 0;
	virtual bool isMissing(const Marker&) const = 0;
	virtual std::pair<unsigned char, unsigned char> getGeno(const Marker&) const = 0;

	unsigned char getAdditiveGeno(const Marker&) const;


	bool addMother(Sample* mom) {return mom == (_mom = (_mom == NULL) ? mom : _mom);}
	bool addFather(Sample* dad) {return dad == (_dad = (_dad == NULL) ? dad : _dad);}
	bool addChild(Sample* child) {return _children.insert(child).second;}

	void setFounder(bool founder){_founder = founder;}
	void setAffected(bool affected){_affected_known = true; _pheno = _affected = affected;}
	void setGender(bool is_male){_sex_known = true; _male = is_male;}
	void setPheno(float pheno){_pheno = pheno;}

	const std::string& getFID() const {return _famid;}
	const std::string& getID() const {return _id;}

	Sample* getFather() const {return _dad;}
	Sample* getMother() const {return _mom;}

	bool isFounder() const {return _founder;}
	bool isGenderKnown() const {return _sex_known;}
	bool isMale() const {return _sex_known && _male;}
	bool isFemale() const {return _sex_known && !_male;}
	bool isAffected() const {return _affected_known && _affected;}
	bool isAffectedKnown() const {return _affected_known;}

	float getPheno() const {return _pheno;}

	void setEnabled(bool enabled=true){_enabled = enabled;}
	bool isEnabled() const {return _enabled;}

	friend class DataSet;

private:

	//! Family ID
	std::string _famid;
	//! Individual ID
	std::string _id;

	Sample* _mom;
	Sample* _dad;
	std::set<Sample*> _children;

	float _pheno;

	bool _sex_known;
	bool _affected_known;
	bool _male;
	bool _affected;
	bool _founder;
	bool _enabled;

	static bool _biallelic;
	static bool _phased;

public:
	static const unsigned char missing_allele;
	static const std::pair<unsigned char, unsigned char> missing_geno;
};

}
}


#endif
