#ifndef PERCENT_H
#define PERCENT_H

#include <stdio.h>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <map>

class Percent{
	private:
		float percent;
		int count;
	
	public:
		Percent(){};
		Percent(float p) : percent(p){count = 0;};
		~Percent(){};
		void incCount(){count++;};
		int getCount(){return count;};
		float getPercent(){return percent;};
		void incCount(int val){count += val;};
		struct mysort{
			bool operator() (const float s1, const float s2) const{
				return s2 < s1;
			}
		};
};
#endif
