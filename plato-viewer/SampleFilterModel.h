/*
 * SampleFilterModel.h
 *
 *  Created on: Jan 16, 2009
 *      Author: gilesjt
 */

#ifndef SAMPLEFILTERMODEL_H_
#define SAMPLEFILTERMODEL_H_
#include <gtkmm.h>
#include <libglademm/xml.h>
#include <string>
#include "Vars.h"

using namespace std;

class SampleFilterModel {
public:
	SampleFilterModel();
	virtual ~SampleFilterModel();

	Gtk::VBox* get_layout(){return layout;}

protected:
	Glib::RefPtr<Gnome::Glade::Xml> refXml;

	Gtk::VBox* layout;

};

#endif /* SAMPLEFILTERMODEL_H_ */
