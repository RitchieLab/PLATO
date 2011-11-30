/*
 * Controller.cpp
 *
 *  Created on: Jan 20, 2009
 *      Author: gilesjt
 */

#include <Helper.h>
#include <MethodException.h>
#include "Controller.h"

Controller::Controller() {
	// TODO Auto-generated constructor stub

}

Controller::~Controller() {
	// TODO Auto-generated destructor stub
}

void Controller::load_data(PlatoProject* proj){
	vector<DataSetObject*> data_sets = proj->get_datasets();
cout << "Has data_sets\n";
cout << "Data set size: " << data_sets.size() << endl;
	for(int i = 0; i < data_sets.size(); i++){
		DataSetObject* ds = data_sets[i];
		if(data_sets[i]->get_type() == PEDMAP){
			try{
				readMapM(ds, ds->get_options());
				readPedM_3vec_set(ds, ds->get_options());
			}catch(MethodException ex){
				cout << ex.what() << endl;
			}
		}
	}
cout << "done iterating data sets\n";
}
