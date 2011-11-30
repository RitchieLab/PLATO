/*
 * CovTraitModel.cpp
 *
 *  Created on: Jan 16, 2009
 *      Author: gilesjt
 */

#include "CovTraitModel.h"

CovTraitModel::CovTraitModel() {
	try{
		refXml = Gnome::Glade::Xml::create(Vars::PLATO_MODELS_GLADE);
	}
	catch(const Gnome::Glade::XmlError& ex){
	}

	layout = 0;

	refXml->get_widget("covtrait_vbox", layout);

}

CovTraitModel::~CovTraitModel() {
	// TODO Auto-generated destructor stub

	delete layout;
//	refXml->unreference();
}
