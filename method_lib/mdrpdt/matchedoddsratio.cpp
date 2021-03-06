//
// C++ Implementation: matchedoddsratio
//
// Description: 
//
//
// Author: Eric Torstenson <torstenson@chgr.mc.vanderbilt.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "matchedoddsratio.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string.h>

namespace MdrPDT {
namespace Evaluation {

bool MatchedOddsRatio::Verbose = false;


MatchedOddsRatio::MatchedOddsRatio() {
	Reset();
}

MatchedOddsRatio::MatchedOddsRatio(const MatchedOddsRatio& other) {
	memcpy((void*)&exposures, (void*)&(other.exposures), 4);
}

MatchedOddsRatio::~MatchedOddsRatio() { }

void MatchedOddsRatio::Reset() {
	familyGroups = 0;
	for (uint i=0; i<4; i++) 
		exposures[i]=0;
}


//These are just for easing the amount of typing required

void MatchedOddsRatio::GenerateReport(ostream *os) {

	if (os == NULL)
		os = &cout;
	int cw = 13;
	uint ExpExp  = exposures[0];
	uint ExpUxp	 = exposures[1];
	uint UxpExp	 = exposures[2];
	uint UxpUxp	 = exposures[3];

	*os<<"\n\nMatched Odds Ratio: \n";
	if (Verbose)
		*os<<setupDetails;
	*os<<setw(cw)<<" "<<"          Unaffected"<<endl;
	*os<<setw(cw)<<"Affected"<<setw(cw)<<"Exposed"<<setw(cw)<<"Unexposed"<<setw(cw)<<"(Total)"<<endl;
	*os<<setw(cw)<<"Exposed"<<setw(cw)<<(ExpExp)<<setw(cw)<<(ExpUxp)<<setw(cw)<<(ExpExp+ExpUxp)<<endl;
	*os<<setw(cw)<<"Unexposed"<<setw(cw)<<UxpExp<<setw(cw)<<UxpUxp<<setw(cw)<<(UxpExp+UxpUxp)<<endl;
	*os<<setw(cw)<<"(Total)"<<setw(cw)<<(ExpExp+UxpExp)<<setw(cw)<<(ExpUxp+UxpUxp)<<setw(cw)<<(ExpExp+UxpExp+ExpUxp+UxpUxp)<<endl;
	*os<<"Matched Odds Ratio: "<<setw(cw)<<(float)ExpUxp/(float)UxpExp<<"\n";
	
}
MatchedOddsRatio& MatchedOddsRatio::operator+(const MatchedOddsRatio& other) {
	exposures[0]+=other.exposures[0];
	exposures[1]+=other.exposures[1];
	exposures[2]+=other.exposures[2];
	exposures[3]+=other.exposures[3];
	familyGroups+=other.familyGroups;	
	return *this;
}
float MatchedOddsRatio::GetRatio() {
	if (exposures[1] == 0 || exposures[2] == 0)
		return 0;
	return (float)exposures[1]/(float)exposures[2];
}

void MatchedOddsRatio::AddFamily(uint truePos, uint falsePos, uint falseNeg, uint trueNeg) {
	uint affecteds[2];
	uint unaffecteds[2];
	affecteds[0] 	= truePos;
	affecteds[1] 	= falsePos;
	unaffecteds[0]	= falseNeg;
	unaffecteds[1]  = trueNeg;
	int cw = 13;

	uint pairsToAssign = truePos+falseNeg;
	
	if (Verbose) {
		stringstream ss;
		if (familyGroups == 0)
			ss<<"Building Paired Associations"<<endl<<setw(cw)<<left<<" True "<<setw(cw)<<" False "<<setw(cw)<<"  False "<<setw(cw)<<"  True "<<endl<<setw(cw)<<"Positive"<<setw(cw)<<" Positive"<<setw(cw)<<"Negative"<<setw(cw)<<"Negative"<<endl;
		familyGroups++;
		ss<<"  "<<left<<setw(cw)<<truePos<<setw(cw)<<falsePos<<setw(cw)<<falseNeg<<setw(cw)<<trueNeg<<endl;
		setupDetails=ss.str();
	}
	//This should spread the pairings out equally, since we are talking about parental groups
	while (pairsToAssign > 0) {
		if (truePos > 0 && falsePos > 0) {
			truePos--;
			falsePos--;
			exposures[0]++;
			pairsToAssign--;
		}
		if (truePos > 0 && trueNeg > 0) {
			truePos--;
			trueNeg--;
			exposures[1]++;
			pairsToAssign--;
		}
 		if (falsePos > 0 && falseNeg > 0) {
			falsePos--;
			falseNeg--;
			exposures[2]++;
			pairsToAssign--;
		}
		if (falseNeg > 0 && trueNeg > 0) {
			falseNeg--;
			trueNeg--;
			exposures[3]++;
			pairsToAssign--;
		}
	}	
}

}


}
