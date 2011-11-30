/*
 * SampleFilterModel.cpp
 *
 *  Created on: Jan 16, 2009
 *      Author: gilesjt
 */

#include "SampleFilterModel.h"

SampleFilterModel::SampleFilterModel() {
	try{
		refXml = Gnome::Glade::Xml::create(Vars::PLATO_MODELS_GLADE);
	}
	catch(const Gnome::Glade::XmlError& ex){
	}

	layout = 0;

	refXml->get_widget("sample_filters_vbox", layout);


}

SampleFilterModel::~SampleFilterModel() {
	// TODO Auto-generated destructor stub
	delete layout;
//	refXml->unreference();
}
