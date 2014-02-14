#ifndef ANALYSIS_CORRECTION_H
#define ANALYSIS_CORRECTION_H

#include <vector>
#include <iostream>

namespace PLATO{
namespace Analysis{

class CorrectionModel;

class Correction {
public:
	enum correction_ENUM { BONFERRONI, FDR };

private:

	template<class T>
	class idx_sorter{
	public:
		idx_sorter(const std::vector<T>& v) : _v(v) {}

		bool operator() (size_t i, size_t j){
			return _v[i] < _v[j];
		}

	private:
		const std::vector<T>& _v;
	};

protected:
	Correction(){}

public:
	virtual ~Correction(){};
	virtual void correct(const std::vector<float>& pval_in, std::vector<float>& pval_out) = 0;

	static Correction* getCorrectionMethod(const CorrectionModel& c);

	static std::string listCorrectionMethods();

protected:
	void initOutVec(const std::vector<float>& pval_in, std::vector<float>& pval_out, std::vector<size_t>& idx_out);

};

class CorrectionModel{
public:
	CorrectionModel(const std::string& s);
	CorrectionModel(const char* s);
	CorrectionModel() : _data(Correction::BONFERRONI) {}
	CorrectionModel(Correction::correction_ENUM c) : _data(c) {}

	operator int() const{return _data;}

//	bool operator<(const CorrectionModel& o) {return static_cast<int>(_data) < static_cast<int>(o._data);}

private:
	Correction::correction_ENUM _data;
};


class BonferroniCorrection : public Correction {
// default ctor fine!
public:
	virtual ~BonferroniCorrection() {}

	virtual void correct(const std::vector<float>& pval_in, std::vector<float>& pval_out);

};

class FDRCorrection : public Correction {
// default ctor fine!
public:
	virtual ~FDRCorrection() {}

	virtual void correct(const std::vector<float>& pval_in, std::vector<float>& pval_out);

};

}
}

namespace std{
istream& operator>>(istream& in, PLATO::Analysis::CorrectionModel& model_out);
ostream& operator<<(ostream& o, const PLATO::Analysis::CorrectionModel& m);
}

#endif
