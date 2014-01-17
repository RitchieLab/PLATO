#ifndef METHODS_ANALYSIS_CORRECTION_H
#define METHODS_ANALYSIS_CORRECTION_H

#include <vector>
#include <iostream>

namespace Methods{
namespace Analysis{

enum correction_ENUM { BONFERRONI, FDR };

class CorrectionModel{
public:
	CorrectionModel() : _data(BONFERRONI) {}
	CorrectionModel(correction_ENUM c) : _data(c) {}

	operator int() const{return _data;}

//	bool operator<(const CorrectionModel& o) {return static_cast<int>(_data) < static_cast<int>(o._data);}

private:
	correction_ENUM _data;
};

class Correction {
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
	virtual ~Correction(){}
	virtual void correct(const std::vector<float>& pval_in, std::vector<float>& pval_out) = 0;
	static Correction* getCorrectionMethod(const CorrectionModel& c);

protected:
	void initOutVec(const std::vector<float>& pval_in, std::vector<float>& pval_out, std::vector<size_t>& idx_out);

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
istream& operator>>(istream& in, Methods::Analysis::CorrectionModel& model_out);
ostream& operator<<(ostream& o, const Methods::Analysis::CorrectionModel& m);
}

#endif
