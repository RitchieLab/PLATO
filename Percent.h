#ifndef PERCENT_H
#define PERCENT_H

class Percent{
	private:
		float percent;
		int count;
	
	public:
		Percent(){};
		Percent(float p) : percent(p){count = 0;};
		virtual ~Percent(){};
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
