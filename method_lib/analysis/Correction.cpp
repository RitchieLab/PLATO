#include "Correction.h"

#include <algorithm>
#include <string>

#include <boost/program_options.hpp>

using std::map;
using std::vector;
using std::sort;

using boost::program_options::validation_error;

namespace PLATO{
namespace Analysis{

CorrectionModel::CorrectionModel(const std::string& s){
	std::istringstream ss(s);
	ss >> (*this);
}

CorrectionModel::CorrectionModel(const char * s){
	std::istringstream ss(s);
	ss >> (*this);
}

Correction* Correction::getCorrectionMethod(const CorrectionModel& c){
	static Correction* bonf = new BonferroniCorrection();
	static Correction* fdr = new FDRCorrection();

	switch((int) c){
	case BONFERRONI:
		return bonf;
	case FDR:
		return fdr;
	default:
		return 0;
	}
}

void Correction::initOutVec(const vector<float>& pval_in, vector<float>& pval_out, vector<size_t>& idx_out){
	pval_out.clear();
	pval_out.reserve(pval_in.size());
	pval_out.insert(pval_out.end(), pval_in.begin(), pval_in.end());

	idx_out.clear();
	idx_out.reserve(pval_in.size());
	for(unsigned int i =0; i < pval_in.size(); i++){
		idx_out.push_back(i);
	}

	// Now, sort the index vector based on the values in the input vector
	sort(idx_out.begin(), idx_out.end(), idx_sorter<float>(pval_in));
}

void BonferroniCorrection::correct(const vector<float>& pval_in, vector<float>& pval_out){
	vector<size_t> idx;
	initOutVec(pval_in, pval_out, idx);
	for(unsigned int i=0; i < pval_out.size(); i++){
		pval_out[idx[i]] = std::min(pval_out[idx[i]] * pval_out.size(), static_cast<float>(1));
	}
}

void FDRCorrection::correct(const vector<float>& pval_in, vector<float>& pval_out){
	vector<size_t> idx;
	initOutVec(pval_in, pval_out, idx);
	for(int i=(pval_out.size()-1); i > 0; --i){
		pval_out[idx[i]] = std::min(pval_out[idx[i+1]], pval_out[idx[i]] * pval_out.size() / (i+1));
	}
}

}
}

namespace std{

istream& operator>>(istream& in,
		PLATO::Analysis::CorrectionModel& model_out) {
	string token;
	in >> token;
	if (token.size() > 0) {
		char s = token[0];
		if (s == 'b' || s == 'B') {
			model_out = PLATO::Analysis::Correction::BONFERRONI;
		} else if (s == 'f' || s == 'F') {
			model_out = PLATO::Analysis::Correction::FDR;
		} else {
			throw validation_error(validation_error::invalid_option_value);
		}
	} else {
		throw validation_error(validation_error::invalid_option_value);
	}
	//    else throw boost::program_options::validation_error("Invalid unit");
	return in;
}
ostream& operator<<(ostream& o, const PLATO::Analysis::CorrectionModel& m){
	switch(m){
	case PLATO::Analysis::Correction::BONFERRONI:
		return o << "Bonferroni";
	case PLATO::Analysis::Correction::FDR:
		return o << "FDR";
	default:
		return o << "unknown";
	}
}

}
