/*
 * Controller.h
 *
 *  Created on: Jun 24, 2010
 *      Author: cozartc
 */

#ifndef CONTROLLER_H_
#define CONTROLLER_H_

#include "Process.h"
//#include "PlatoProject.h"
#include "libsqlitewrapped.h"

using namespace std;

class Controller
{
public:
		Controller();
		virtual ~Controller();
		static void drop_table(Database*, string);
		static long long execute_sql(Query&, string);
};

#endif /* CONTROLLER_H_ */
