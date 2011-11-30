#ifndef FINALIZE_H
#define FINALIZE_H

#include <stdio.h>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
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
using namespace std;
using namespace Methods;

class Finalize{
	private:
		
	public:
		Finalize(){
		};
		~Finalize(){};
		void finish(vector<Marker*>*, vector<Sample*>*, vector<Family*>*);
};
#endif
