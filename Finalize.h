#ifndef FINALIZE_H
#define FINALIZE_H

#include <stdio.h>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
//#include <occi.h>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include <bitset>
#include "Globals.h"
#include "General.h"
#include "Marker.h"
#include "Sample.h"
#include "Family.h"
//#include "Helper.h"
//using namespace oracle::occi;
using namespace std;


class Finalize{
	private:
		
	public:
		Finalize(){
		};
		~Finalize(){};
		void finish(vector<Marker*>*, vector<Sample*>*, vector<Family*>*);
};
#endif
