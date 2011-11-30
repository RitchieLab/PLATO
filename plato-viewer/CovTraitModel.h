/*
 * CovTraitModel.h
 *
 *  Created on: Jan 16, 2009
 *      Author: gilesjt
 */

#ifndef COVTRAITMODEL_H_
#define COVTRAITMODEL_H_
#include <gtkmm.h>
#include <libglademm/xml.h>
#include <string>
#include "Vars.h"

using namespace std;

class CovTraitModel {
public:
	CovTraitModel();
	virtual ~CovTraitModel();

	Gtk::VBox* get_layout(){return layout;}

protected:
	Glib::RefPtr<Gnome::Glade::Xml> refXml;

	Gtk::VBox* layout;

};

#endif /* COVTRAITMODEL_H_ */
