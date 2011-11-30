/*
 * OptionsList.h
 *
 *  Created on: Jan 12, 2009
 *      Author: gilesjt
 */

#ifndef OPTIONSLIST_H_
#define OPTIONSLIST_H_

#include <stdio.h>
#include <string.h>
#include <string>
#include <iostream>
#include <vector>
#include <gtkmm.h>
#include <libglademm/xml.h>
#include <General.h>
#include <MethodException.h>
#include "OptionData.h"
#include "Vars.h"

using namespace std;

class OptionsList {
public:
	OptionsList(string f);
	virtual ~OptionsList();
	int num_opts(){return (int)option_data.size();}
	OptionData* get_option_data(int i){
		if(!option_data[i]){
			return NULL;
		}
		else{
			return option_data[i];
		}
	}

protected:
	string design_file;
	vector<OptionData*> option_data;
	Glib::RefPtr<Gnome::Glade::Xml> optionsXml;
};

#endif /* OPTIONSLIST_H_ */
