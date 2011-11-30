/*
 * ThresholdsModel.h
 *
 *  Created on: Jan 23, 2009
 *      Author: gilesjt
 */

#ifndef THRESHOLDSMODEL_H_
#define THRESHOLDSMODEL_H_

#include <gtkmm.h>
#include <libglademm.h>
#include "Vars.h"

class ThresholdsModel {
public:
	ThresholdsModel();
	virtual ~ThresholdsModel();

	Gtk::VBox* get_layout(){return layout;}
	void set_fields(vector<int>);

protected:
	Glib::RefPtr<Gnome::Glade::Xml> refXml;

	Gtk::VBox* layout;
	Gtk::Entry* locus_min_entry;
	Gtk::Entry* locus_max_entry;
	Gtk::Entry* sample_min_entry;
	Gtk::Entry* sample_max_entry;
	Gtk::Entry* pedigree_min_entry;
	Gtk::Entry* pedigree_max_entry;

};

#endif /* THRESHOLDSMODEL_H_ */
