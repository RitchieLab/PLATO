/*
 * OptionData.h
 *
 *  Created on: Jan 14, 2009
 *      Author: gilesjt
 */

#ifndef OPTIONDATA_H_
#define OPTIONDATA_H_

#include <vector>
#include <string.h>
#include <string>
#include <gtkmm.h>

using namespace std;

class OptionData {
public:
	OptionData();
	virtual ~OptionData();

	void add_required(string s){required.push_back(s);}
	void add_disabled(string s){disabled.push_back(s);}
	Gtk::Widget* get_layout(){
		if(!layout){
			return NULL;
		}
		return (Gtk::Widget*)layout;
	}
	void set_layout(Gtk::Layout* l){layout = l;}
	void set_name(string s){name = s;}

protected:
	vector<string> required;
	vector<string> disabled;

	string name;

	Gtk::Layout* layout;
};

#endif /* OPTIONDATA_H_ */
