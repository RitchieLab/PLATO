#ifndef COVTRAIT_H
#define COVTRAIT_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <strings.h>
#include <string>
#include "Globals.h"
//#include "Families.h"
//#include "Markers.h"
using namespace std;

namespace Methods{
class CovTrait{

protected:
	string name;
	bool enabled;
	double value;

public:
		CovTrait(){name = ""; enabled = true; value = NAN;}
		CovTrait(string nam){name = nam;}
		CovTrait(string nam, double val){name = nam; value = val;}
		CovTrait(string nam, double val, bool use){name = nam; value = val; enabled = use;}
		~CovTrait(){}

		string get_name(){return name;}
		void set_name(string nam){name = nam;}
		double get_value(){return value;}
		void set_value(double val){value = val;}
		bool isEnabled(){return enabled;}
		void set_enabled(bool v){enabled = v;}


};
};

#endif
