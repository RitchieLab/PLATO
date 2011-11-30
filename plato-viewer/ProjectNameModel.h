/*
 * ProjectNameModel.h
 *
 *  Created on: Jan 16, 2009
 *      Author: gilesjt
 */

#ifndef PROJECTNAMEMODEL_H_
#define PROJECTNAMEMODEL_H_

#include <gtkmm.h>
#include <libglademm/xml.h>
#include <string>
#include "Vars.h"

using namespace std;

class ProjectNameModel {
public:
	ProjectNameModel();
	virtual ~ProjectNameModel();

	Gtk::Layout* get_layout(){return layout;}
	string get_name(){return name_entry->get_text();}

protected:
	Glib::RefPtr<Gnome::Glade::Xml> refXml;

	Gtk::Layout* layout;
	Gtk::Entry* name_entry;
};

#endif /* PROJECTNAMEMODEL_H_ */
