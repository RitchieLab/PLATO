#ifndef METHODS_SAMPLE_H
#define METHODS_SAMPLE_H

#include <string>
#include <set>

namespace Methods{

class DataSet;

class Sample{

public:
	virtual ~Sample();

protected:

	Sample(const std::string& famid, const std::string& id);

public:
	static Sample* create(const std::string& famid, const std::string& id, unsigned int n_genos=0);
	static Sample* create(const std::string& id, unsigned int n_genos=0){return create(id, id, n_genos);}

private:
	// No copying or assignment!!
	Sample(const Sample&);
	Sample& operator=(const Sample&);

public:
	virtual void appendGenotype(unsigned char geno1, unsigned char geno2) = 0;
	virtual void appendMissingGenotype() = 0;
	virtual bool isMissing(unsigned int pos) const = 0;
	virtual std::pair<unsigned char, unsigned char> getGeno(unsigned int pos) const = 0;

	bool addMother(const Sample* mom) {return mom == (_mom = (_mom == NULL) ? mom : _mom);}
	bool addFather(const Sample* dad) {return dad == (_dad = (_dad == NULL) ? dad : _dad);}
	bool addChild(const Sample* child) {return _children.insert(child).second;}

	void setFounder(bool founder){_founder = founder;}
	void setAffected(bool affected){_affected_known = true; _affected = affected;}
	void setGender(bool is_male){_sex_known = true; _male = is_male;}

	void addTrait(float trait){ _traits.push_back(trait);}
	void setPhenoPos(unsigned int pos){ _pheno_pos = static_cast<unsigned char>(pos);}

	float getTrait(unsigned int pos) const {return pos < _traits.size() ? _traits[pos] : missing_trait;}
	float getPhenotype() const {return getTrait(_pheno_pos);}

	bool isFounder() const {return _founder;}
	bool isGenderKnown() const {return _sex_known;}
	bool isMale() const {return _sex_known && _male;}
	bool isFemale() const {return _sex_known && !_male;}
	bool isAffected() const {return _affected_known && _affected;}
	bool isAffectedKnown() const {return _affected_known;}

	friend class DataSet;

private:

	//! Family ID
	std::string _famid;
	//! Individual ID
	std::string _id;

	Sample* _mom;
	Sample* _dad;
	std::set<Sample*> _children;


	bool _sex_known;
	bool _affected_known;
	bool _male;
	bool _affected;
	bool _founder;
	bool _enabled;

	unsigned char _pheno_pos;
	deque<float> _traits;

	static float missing_trait;
	static bool _biallelic = true;
	static bool _phased = false;

protected:
	static std::pair<unsigned char, unsigned char> missing_geno;

};

}

#endif
