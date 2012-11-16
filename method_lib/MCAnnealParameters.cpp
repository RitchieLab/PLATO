/*
 * MCAnnealParameters.cpp
 *
 *  Created on: Nov 3, 2009
 *      Author: gilesjt
 */

#include "MCAnnealParameters.h"

namespace Methods {

MCAnnealParameters::MCAnnealParameters() {
	// TODO Auto-generated constructor stub
	nburn = 1000;
	niter = 25000;
	hyperpars = 0;
	update = 0;
	output = 4;
}

MCAnnealParameters::MCAnnealParameters(int nb, int ni, int hp, int u, int o){
	nburn = nb;
	niter = ni;
	hyperpars = hp;
	update = u;
	output = o;
}


MCAnnealParameters::~MCAnnealParameters() {
	// TODO Auto-generated destructor stub
}

}
