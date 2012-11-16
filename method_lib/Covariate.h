#ifndef COVARIATE_H
#define COVARIATE_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "CovTrait.h"
#include "Globals.h"
//#include "Families.h"
//#include "Markers.h"
using namespace std;

namespace Methods{
class Covariate : public CovTrait{
	public:
		Covariate(){}
		Covariate(string nam){name = nam;}
		Covariate(string nam, double val){name = nam; value = val;}
		Covariate(string nam, double val, bool use){name = nam; value = val; enabled = use;}
		~Covariate(){}
};
};

#endif
