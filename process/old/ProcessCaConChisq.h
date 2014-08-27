#ifndef PROCESSCACONCHISQ_H
#define PROCESSCACONCHISQ_H

#include "Process.h"
#include <vector>
#include <Marker.h>

class ProcessCaConChisq: public ProcessImpl<ProcessCaConChisq> {

private:
	static std::string stepname;

	std::vector<double> chi_geno;
	std::vector<double> chi_allele;
	std::vector<double> chi_arm;
	std::vector<double> pval_arm;
	std::vector<double> pval_allele;
	std::vector<double> pval_geno;
	std::vector<double> pval_geno_exact;
	std::vector<double> pval_allele_exact;
	std::vector<double> odds_ratio;
	std::vector<double> ci_l;
	std::vector<double> ci_u;
	std::vector<int> geno_df;
	std::vector<int> allele_df;
	std::vector<int> arm_df;

	//groups
	double gchi_geno;
	double gchi_allele;
	double gchi_arm;
	double gpval_arm;
	double gpval_allele;
	double gpval_geno;
	double gpval_geno_exact;
	double gpval_allele_exact;
	double godds_ratio;

	std::vector<Methods::Marker*> good_markers;

public:
	ProcessCaConChisq() {name = "Chisquare test";}
	virtual ~ProcessCaConChisq() {}

protected:
	virtual void PrintSummary();
	virtual void filter();
	virtual void process(Methods::DataSet*);

};
#endif
