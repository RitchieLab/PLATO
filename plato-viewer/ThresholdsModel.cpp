/*
 * ThresholdsModel.cpp
 *
 *  Created on: Jan 23, 2009
 *      Author: gilesjt
 */

#include "ThresholdsModel.h"

ThresholdsModel::ThresholdsModel() {
	try{
		refXml = Gnome::Glade::Xml::create(Vars::PLATO_MODELS_GLADE);
	}
	catch(const Gnome::Glade::XmlError& ex){
	}

	layout = 0;

	refXml->get_widget("thresholds_vbox", layout);
	refXml->get_widget("locus_threshold_min_entry", locus_min_entry);
	refXml->get_widget("locus_threshold_max_entry", locus_max_entry);
	refXml->get_widget("sample_threshold_min_entry", sample_min_entry);
	refXml->get_widget("sample_threshold_max_entry", sample_max_entry);
	refXml->get_widget("pedigree_threshold_min_entry", pedigree_min_entry);
	refXml->get_widget("pedigree_threshold_max_entry", pedigree_max_entry);

}

ThresholdsModel::~ThresholdsModel() {
	// TODO Auto-generated destructor stub

	delete layout;
}

void ThresholdsModel::set_fields(vector<int> opts){
	for(int oi = 0; oi < opts.size(); oi++){
		switch(opts[oi]){
		case OptionsItem::oi_locus_threshold_min:
			locus_min_entry->set_sensitive(true);
			break;
		case OptionsItem::oi_locus_threshold_max:
			locus_max_entry->set_sensitive(true);
			break;
		case OptionsItem::oi_sample_threshold_min:
			sample_min_entry->set_sensitive(true);
			break;
		case OptionsItem::oi_sample_threshold_max:
			sample_max_entry->set_sensitive(true);
			break;
		case OptionsItem::oi_pedigree_threshold_min:
			pedigree_min_entry->set_sensitive(true);
			break;
		case OptionsItem::oi_pedigree_threshold_max:
			pedigree_max_entry->set_sensitive(true);
			break;
		default:
			break;
		}
	}

}
