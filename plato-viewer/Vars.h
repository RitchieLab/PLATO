/*
 * Vars.h
 *
 *  Created on: Sep 29, 2008
 *      Author: gilesjt
 */

#ifndef VARS_H_
#define VARS_H_
#include <gtkmm.h>
#include <glibmm.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <map>
using namespace std;

//namespace Mappings{
typedef enum DataTypes{
	PEDMAP,
	TPEDFAM,
	BINARY
} DataTypes;

class Processes{
public:
	enum{
		p_markerefficiency,
		p_sampleefficiency
	};
};

class OptionsGroup{
public:
	enum {
	og_thresholds
	};
};

class OptionsItem{
public:
	enum{
	oi_locus_threshold_min,
	oi_locus_threshold_max,
	oi_sample_threshold_min,
	oi_sample_threshold_max,
	oi_pedigree_threshold_min,
	oi_pedigree_threshold_max
	};
};
//};

class Vars {
public:

	static string PLATO_GUI_GLADE;
	static string PLATO_OPTIONS_GLADE;
	static string PLATO_MODELS_GLADE;
//	static map<string, string> PLATO_DESIGN_FILES;
	static string LOCUS_FILTER_GUI;
	static string RESOURCES_DIR;
	static string REQUIRE_TOKEN;
	static string DISABLE_TOKEN;

	static GMutex *mutex_Project;
	static double COVTRAITMISSINGDOUB;

	static map<string, int> processes;
	static map<int, map<int, vector<int> > > Vars::process_options_map;

	static void initialize();
};

#endif /* VARS_H_ */
