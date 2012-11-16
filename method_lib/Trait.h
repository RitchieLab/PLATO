#ifndef TRAIT_H
#define TRAIT_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include "Globals.h"
#include "CovTrait.h"
//#include "Families.h"
//#include "Markers.h"
using namespace std;

namespace Methods{
class Trait : public CovTrait{
	public:
		Trait(){};
		Trait(string nam){name = nam;}
		Trait(string nam, double val){name = nam; value = val;}
		Trait(string nam, double val, bool use){name = nam; value = val; enabled = use;}
		~Trait(){}

};
};

#endif
