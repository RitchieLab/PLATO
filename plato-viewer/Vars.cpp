/*
 * Vars.cpp
 *
 *  Created on: Sep 29, 2008
 *      Author: gilesjt
 */

#include <algorithm>
#include "Vars.h"


 GMutex* Vars::mutex_Project = NULL;
 string Vars::RESOURCES_DIR = "resources/";
 string Vars::PLATO_GUI_GLADE = Vars::RESOURCES_DIR + "plato-gui.glade";
 string Vars::PLATO_MODELS_GLADE = Vars::RESOURCES_DIR + "plato-models.glade";
 string Vars::PLATO_OPTIONS_GLADE = Vars::RESOURCES_DIR + "plato-options.glade";
 string Vars::LOCUS_FILTER_GUI = "locus_filter";
 double Vars::COVTRAITMISSINGDOUB = -99999;
 string Vars::REQUIRE_TOKEN = "*R*";
 string Vars::DISABLE_TOKEN = "*D*";
 map<string, int> Vars::processes;
 //Processes, OptionsGroup, OptionsItem
 map<int, map<int, vector<int> > > Vars::process_options_map;

 void Vars::initialize(){
//	 Vars::PLATO_DESIGN_FILES[Vars::LOCUS_FILTER_GUI] = Vars::RESOURCES_DIR + "locus-filter.txt";

	 {
		 Vars::processes["Sample Efficiency"] = Processes::p_sampleefficiency;
	 }

	 {
		 Vars::processes["Marker Efficiency"] = Processes::p_markerefficiency;
		 map<int, vector<int> > grouping;
		 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
		 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
		 Vars::process_options_map[Processes::p_markerefficiency] = grouping;
	 }
	 {
		 Vars::processes[""] = -1;
	 }
 }
