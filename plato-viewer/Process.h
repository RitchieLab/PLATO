/*
 * Process.h
 *
 *  Created on: Jan 23, 2009
 *      Author: gilesjt
 */

#ifndef PROCESS_H_
#define PROCESS_H_

#include <DataSet.h>
#include <StepOptions.h>

class Process {
public:
	Process();
	virtual ~Process();
	void set_options(StepOptions opts){options = opts;}
	StepOptions get_options(){return options;}
	string get_name(){return name;}

protected:
	DataSet* data_set;
	StepOptions options;
	string name;
};

#endif /* PROCESS_H_ */
