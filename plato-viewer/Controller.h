/*
 * Controller.h
 *
 *  Created on: Jan 20, 2009
 *      Author: gilesjt
 */

#ifndef CONTROLLER_H_
#define CONTROLLER_H_

//#include <Helper.h>
#include <vector>
#include "DataSetObject.h"
#include "PlatoProject.h"

using namespace std;

class Controller {
public:
	Controller();
	virtual ~Controller();

	static void load_data(PlatoProject*);
};

#endif /* CONTROLLER_H_ */
