/*
 * LocusFilterModel.h
 *
 *  Created on: Jan 16, 2009
 *      Author: gilesjt
 */

#ifndef LOCUSFILTERMODEL_H_
#define LOCUSFILTERMODEL_H_
#include <gtkmm.h>
#include <libglademm/xml.h>
#include <string>
#include "Vars.h"

using namespace std;

class LocusFilterModel {
public:
	LocusFilterModel();
	virtual ~LocusFilterModel();
	Gtk::VBox* get_layout(){return layout;}

protected:
	Glib::RefPtr<Gnome::Glade::Xml> refXml;

	Gtk::VBox* layout;

};

#endif /* LOCUSFILTERMODEL_H_ */
