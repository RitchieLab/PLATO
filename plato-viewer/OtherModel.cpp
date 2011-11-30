/*
 * OtherModel.cpp
 *
 *  Created on: Jan 16, 2009
 *      Author: gilesjt
 */

#include "OtherModel.h"

OtherModel::OtherModel() {
	try{
		refXml = Gnome::Glade::Xml::create(Vars::PLATO_MODELS_GLADE);
	}
	catch(const Gnome::Glade::XmlError& ex){
	}

	layout = 0;

	refXml->get_widget("other_alignment", layout);

}

OtherModel::~OtherModel() {
	// TODO Auto-generated destructor stub
	delete layout;
	refXml->unreference();
}
