/*
 * ProjectNameModel.cpp
 *
 *  Created on: Jan 16, 2009
 *      Author: gilesjt
 */

#include "ProjectNameModel.h"

ProjectNameModel::ProjectNameModel() {
	try{
		refXml = Gnome::Glade::Xml::create(Vars::PLATO_MODELS_GLADE);
	}
	catch(const Gnome::Glade::XmlError& ex){
	}

	layout = 0;
	name_entry = 0;

	refXml->get_widget("project_name_layout", layout);
	refXml->get_widget("project_name_entry", name_entry);
	name_entry->set_text("My Project");
}

ProjectNameModel::~ProjectNameModel() {
	// TODO Auto-generated destructor stub
	delete name_entry;
	delete layout;
	//refXml->unreference();
}
