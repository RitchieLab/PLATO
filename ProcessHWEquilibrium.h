#ifndef PROCESSHWEQUILIBRIUM_H
#define PROCESSHWEQUILIBRIUM_H


#include "Process.h"
#include <HWEquilibrium.h>
#include <vector>

class ProcessHWEquilibrium : public ProcessImpl<ProcessHWEquilibrium>{
private:
	static const std::string stepname;

	std::string defaultinsert;
	std::string hweptinsert;
	std::string parentalinsert;
	std::string genderinsert;
	std::string ccinsert;

	float hw_O;
	float hw_OM;
	float hw_OF;
	float hw_P;
	float hw_PM;
	float hw_PF;
	float hw_C;
	float hw_CM;
	float hw_CF;
	float hw_Ca;
	float hw_CaF;
	float hw_CaM;
	float hw_Con;
	float hw_ConM;
	float hw_ConF;

	std::vector<float> hwO;
	std::vector<float> hwP;
	std::vector<float> hwPM;
	std::vector<float> hwPD;
	std::vector<float> hwC;
	std::vector<float> hwCM;
	std::vector<float> hwCF;
	std::vector<float> hwCa;
	std::vector<float> hwCaF;
	std::vector<float> hwCaM;
	std::vector<float> hwCon;
	std::vector<float> hwConF;
	std::vector<float> hwConM;

	std::vector<Methods::Marker*> good_markers;


public:
	ProcessHWEquilibrium(){name="Hardy-Weinberg Calculations";}
	virtual ~ProcessHWEquilibrium(){}

protected:
	virtual void PrintSummary();
	virtual void process(Methods::DataSet*);
	virtual void FilterSummary();

private:
	void doFilter(Methods::Marker*, Methods::HWEquilibrium*);

};
#endif
