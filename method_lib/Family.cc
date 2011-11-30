#include <stdio.h>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include "Helpers.h"
#include "Options.h"
#include "Family.h"
using namespace std;

namespace Methods{
string Family::getFamID(){
	if(opts::_TODIGIT_ && famid_digit != -1){
		return getString<int>(famid_digit);
	}
	return famid;
}

string Family::toString(){
	if(opts::_TODIGIT_ && famid_digit != -1){
		return getString<int>(famid_digit);
	}
	return (famid);
}

void Family::setEnabled(bool v){
	enabled = v;
    for(int i = 0; i < (int)samples.size(); i++){
	    Sample* samp = samples[i];
	    if(samp->isEnabled()){
		    samp->setEnabled(v);
		}
	}

}

void Family::addNonFounder(Sample* s){
	int fsize = nonfounders.size();
	for(int f = 0; f < fsize; f++){
		Sample* fsamp = nonfounders[f];
		if(fsamp->getInd() == s->getInd() && fsamp->getSex() == s->getSex()){
			return;
		}
	}
	nonfounders.push_back(s);
}

void Family::addFounder(Sample* s){
	int fsize = founders.size();
	for(int f = 0; f < fsize; f++){
		Sample* fsamp = founders[f];

		if(fsamp->getInd() == s->getInd() && fsamp->getSex() == s->getSex()){
			return;
		}
	}
	founders.push_back(s);
}

void Family::AddInd(Sample* newind){
	vector<Sample*>::iterator s_iter;
	s_iter = find(samples.begin(), samples.end(), newind);
	if(s_iter == samples.end()){
		samples.push_back(newind);
	}
	if((newind->getDad() != NULL || newind->getMom() != NULL)
			&& (newind->getDadID().length() != 0 || newind->getMomID().length() != 0)){
		has_children = true;
	}
	if(newind->getDad() == NULL && newind->getMom() == NULL && newind->getDadID().length() == 0 && newind->getMomID().length() == 0){
		has_parents = true;
	}
}

int Family::getTotalEnabledInds(){
	int total = 0;
	vector<Sample*>::iterator s_iter;
	for(s_iter = samples.begin(); s_iter != samples.end(); s_iter++){
		if((*s_iter)->isEnabled()){
			total++;
		}
	}
	return total;
}

}
