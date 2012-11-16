#include "Chrom.h"
using namespace Methods;
bool Chrom::haveMarker(int v){
	MKR::iterator m_iter;
	m_iter = markers.find(v);
	if(m_iter == markers.end()){
		return false;
	}
	return true;
}

Chrom::Chrom(string c){
	chrom = c;
	Percent* p = new Percent(100);
	percent_counts[100] = *p;
	delete p;
	p = new Percent(99.99);
	percent_counts[99.99] = *p;
	delete p;
	p = new Percent(99);
	percent_counts[99] = *p;
	delete p;
	p = new Percent(95);
	percent_counts[95] = *p;
	delete p;
	p = new Percent(90);
	percent_counts[90] = *p;
	delete p;
	p = new Percent(85);
	percent_counts[85] = *p;
	delete p;
	p = new Percent(80);
	percent_counts[80] = *p;
	delete p;
	/*
	p = new Percent(70);
	percent_counts[70] = *p;
	delete p;
	p = new Percent(65);
	percent_counts[65] = *p;
	delete p;
	p = new Percent(60);
	percent_counts[60] = *p;
	delete p;
	p = new Percent(55);
	percent_counts[55] = *p;
	delete p;
	p = new Percent(50);
	percent_counts[50] = *p;
	delete p;
	*/
	total = 0;
}

void Chrom::incPC(float val){
	PERCOUNT::iterator p_iter;
	for(p_iter = percent_counts.begin(); p_iter != percent_counts.end(); p_iter++){
		if(val == 100.0 && p_iter->first == 100){
			p_iter->second.incCount();
		}
		else if(p_iter->first == 99.99 && val < 100.0){
			p_iter->second.incCount();
		}
		else if(val < p_iter->second.getPercent() && p_iter->first != 100 && p_iter->first != 101){
			p_iter->second.incCount();
		}
	}
}
