/*
 * AnnealParameters.cpp
 *
 *  Created on: Nov 3, 2009
 *      Author: gilesjt
 */

#include "AnnealParameters.h"

namespace Methods {

AnnealParameters::AnnealParameters() {
	// TODO Auto-generated constructor stub
	start = 2;
	end = 1;
	iter = 50000;
	update = 1000;
}

AnnealParameters::AnnealParameters(int s, int e, int i, bool e, int u){
	start = s;
	end = e;
	iter = i;
	earlyout = e;
	update = u;
}

AnnealParameters::~AnnealParameters() {
	// TODO Auto-generated destructor stub
}

}
