/*
 * TreeAnnealParameters.cpp
 *
 *  Created on: Nov 3, 2009
 *      Author: gilesjt
 */

#include "TreeAnnealParameters.h"

namespace Methods {

TreeAnnealParameters::TreeAnnealParameters() {
	// TODO Auto-generated constructor stub
	treesize = 8;
	opers = 1;
	minmass = 0;

}

TreeAnnealParameters::TreeAnnealParameters(int t, int o, int mm, int n){
	treesize = t;
	opers = o;
	minmass = mm;
	n1 = n;
}


TreeAnnealParameters::~TreeAnnealParameters() {
	// TODO Auto-generated destructor stub
}

}
