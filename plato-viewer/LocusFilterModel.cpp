/*
 * LocusFilterModel.cpp
 *
 *  Created on: Jan 16, 2009
 *      Author: gilesjt
 */

#include "LocusFilterModel.h"

LocusFilterModel::LocusFilterModel() {
	try{
		refXml = Gnome::Glade::Xml::create(Vars::PLATO_MODELS_GLADE);
	}
	catch(const Gnome::Glade::XmlError& ex){
	}

	layout = 0;

	refXml->get_widget("locus_filter_vbox", layout);


}

LocusFilterModel::~LocusFilterModel() {
	// TODO Auto-generated destructor stub

	delete layout;
//	refXml->unreference();
}
