/*
 * OtherModel.h
 *
 *  Created on: Jan 16, 2009
 *      Author: gilesjt
 */

#ifndef OTHERMODEL_H_
#define OTHERMODEL_H_
#include <gtkmm.h>
#include <libglademm/xml.h>
#include <string>
#include "Vars.h"

using namespace std;

class OtherModel {
public:
	OtherModel();
	virtual ~OtherModel();

	Gtk::Alignment* get_layout(){return layout;}

protected:
	Glib::RefPtr<Gnome::Glade::Xml> refXml;

	Gtk::Alignment* layout;

};

#endif /* OTHERMODEL_H_ */
