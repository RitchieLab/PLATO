#ifndef CHROM_H
#define CHROM_H
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
#include <map>
#include "Percent.h"
#include "Marker.h"
#include "General.h"

using namespace Methods;
typedef map<float,Percent, Percent::mysort> PERCOUNT;
typedef map<int,Marker> MKR;
typedef vector<Marker> BPMKR;

class Chrom{
	private:
		string chrom;
		PERCOUNT percent_counts;		
		int total;
		MKR markers;
		BPMKR bp_markers;
		
	public:
		Chrom(){};
		Chrom(string);
		~Chrom(){};
		void incPC(float);
		string getChrom(){return chrom;};
		PERCOUNT getPC(){return percent_counts;};
		void incTotal(){total++;};
		int getTotal(){return total;};
		struct mysort{
			bool operator()(const string s1, const string s2) const{
				if(s1 == "X"){
					return false;
				}
				else if(s2 == "X"){
					return true;
				}
				else{
					return atoi(s1.c_str()) < atoi(s2.c_str());
				}
			}
		};
		void insert_marker(Marker m){
		};
		MKR* getMarkers(){
			return &markers;
		};
		BPMKR* getBPMarkers(){
			return &bp_markers;
		}
		bool haveMarker(int);
};
#endif
