#ifndef CHROM_H
#define CHROM_H

#include <string>
#include <vector>
#include <map>

#include "Percent.h"
#include "Marker.h"
#include "General.h"

//using namespace Methods;

typedef std::map<float,Percent, Percent::mysort> PERCOUNT;
typedef std::map<int,Methods::Marker> MKR;
typedef std::vector<Methods::Marker> BPMKR;

class Chrom{
	private:
		std::string chrom;
		PERCOUNT percent_counts;		
		int total;
		MKR markers;
		BPMKR bp_markers;
		
	public:
		Chrom(){};
		Chrom(string);
		~Chrom(){};
		void incPC(float);
		std::string getChrom(){return chrom;};
		PERCOUNT getPC(){return percent_counts;};
		void incTotal(){total++;};
		int getTotal(){return total;};
		struct mysort{
			bool operator()(const std::string s1, const std::string s2) const{
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
		void insert_marker(Methods::Marker m){
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
