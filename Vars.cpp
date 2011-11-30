/*
 * Vars.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: cozartc
 */
#include <algorithm>
#include <string>
#include "Vars.h"

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

// GMutex* Vars::mutex_Project = NULL;
 string Vars::BLANK = "BLANK";
 string Vars::BATCH = "BATCH";
 string Vars::RESOURCES_DIR = "resources/";
 string Vars::PLATO_GUI_GLADE = Vars::RESOURCES_DIR + "plato-gui.glade";
 string Vars::PLATO_MODELS_GLADE = Vars::RESOURCES_DIR + "plato-models.glade";
 string Vars::PLATO_OPTIONS_GLADE = Vars::RESOURCES_DIR + "plato-options.glade";
 string Vars::LOCUS_FILTER_GUI = "locus_filter";
 double Vars::COVTRAITMISSINGDOUB = -99999;
 string Vars::REQUIRE_TOKEN = "*R*";
 string Vars::DISABLE_TOKEN = "*D*";
 string Vars::PLOT_GENOME = "Genome";
 string Vars::PLOT_PEDIGREE = "Pedigree";
 string Vars::PLOT_SAMPLE = "Sample";
 string Vars::PLOT_IMPORT_FILE = "Imported File";
 string Vars::LOCUS_LOCATION = "bploc";
 string Vars::LOCUS_TABLE = "LOCI";
 string Vars::SAMPLE_TABLE = "SAMPLES";
 string Vars::PEDIGREE_TABLE = "PEDIGREES";
 int Vars::FILE_TYPE = 1;
 int Vars::TABLE_TYPE = 2;
 int Vars::IMPORT_TABLE_TYPE = 3;

 map<string, int> Vars::processes;
 //Processes, OptionsGroup, OptionsItem
 map<int, map<int, vector<int> > > Vars::process_options_map;

 vector<double> Vars::distributionNegLog(int count){
     vector<double> distro;

     distro.resize(count);

     for(int i = 1; i <= count; i++){
         distro[i] = (double)abs(log10((double((double)i / (double)count))));
     }

     sort(distro.begin(), distro.end());

     return distro;
 }

 void Vars::initialize()
 {
	 //Vars::PLATO_DESIGN_FILES[Vars::LOCUS_FILTER_GUI] = Vars::RESOURCES_DIR + "locus-filter.txt";
	 {
		 Vars::processes["Sample Efficiency"] = Processes::p_sampleefficiency;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_sample_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_sample_threshold_max);
                 Vars::process_options_map[Processes::p_sampleefficiency] = grouping;
	 }

         {
                 Vars::processes["Chisquare Tests"] = Processes::p_chisquare;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_all);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_all_children);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_unaffected_spouses);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_unknown_spouses);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_random_child);
                 Vars::process_options_map[Processes::p_chisquare] = grouping;
         }
         {
                 Vars::processes["Logistic Regression"] = Processes::p_logreg;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_logreg);
                 grouping[OptionsGroup::og_covariates].push_back(OptionsItem::oi_covariates);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_transform);
                 Vars::process_options_map[Processes::p_logreg] = grouping;
         }
         {
                 Vars::processes["Linear Regression"] = Processes::p_linreg;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_linreg);
                 grouping[OptionsGroup::og_covariates].push_back(OptionsItem::oi_covariates);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_xchr);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_transform);
                 Vars::process_options_map[Processes::p_linreg] = grouping;
         }

        {
                 Vars::processes["Gender Check"] = Processes::p_gendererror;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_sample_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_sample_threshold_max);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_gendererror] = grouping;
         }

        {
                 Vars::processes["CMH"] = Processes::p_cmh;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_clusters);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_cmh);
                 Vars::process_options_map[Processes::p_cmh] = grouping;
         }

        {
                 Vars::processes["HWE"] = Processes::p_hwe;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_all);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_all_children);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_unaffected_spouses);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_unknown_spouses);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_random_child);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_parental);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_casecontrol);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_gender);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_hwe);
                 grouping[OptionsGroup::og_filters].push_back(OptionsItem::oi_filter_overall);

                 Vars::process_options_map[Processes::p_hwe] = grouping;
         }
        {
                 Vars::processes["Mendelian Error"] = Processes::p_mendelianerror;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_pedigree_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_pedigree_threshold_max);
                 grouping[OptionsGroup::og_filters].push_back(OptionsItem::oi_zero);
                 grouping[OptionsGroup::og_filters].push_back(OptionsItem::oi_zero_l2);
                 grouping[OptionsGroup::og_filters].push_back(OptionsItem::oi_zero_l2_fams);
                 Vars::process_options_map[Processes::p_mendelianerror] = grouping;
         }

	 {
		 Vars::processes["Marker Efficiency"] = Processes::p_markerefficiency;
		 map<int, vector<int> > grouping;
		 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
		 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
		 Vars::process_options_map[Processes::p_markerefficiency] = grouping;
	 }
     {
                 Vars::processes["TDT"] = Processes::p_tdt;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_tdt] = grouping;
      }
     {
                 Vars::processes["Allele Frequency"] = Processes::p_allelefrequency;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_all);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_all_children);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_unaffected_spouses);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_unknown_spouses);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_random_child);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_parental);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_casecontrol);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_gender);
                 grouping[OptionsGroup::og_filters].push_back(OptionsItem::oi_filter_overall);
                 grouping[OptionsGroup::og_filters].push_back(OptionsItem::oi_freq_file);
                 grouping[OptionsGroup::og_filters].push_back(OptionsItem::oi_removehet);
                 grouping[OptionsGroup::og_filters].push_back(OptionsItem::oi_removemono);

                 Vars::process_options_map[Processes::p_allelefrequency] = grouping;
     }
     {
                 Vars::processes["Output Beagle"] = Processes::p_outputbeagle;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_disease_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);

//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputbeagle] = grouping;
     }
     {
                 Vars::processes["Output FBAT"] = Processes::p_outputfbat;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputfbat] = grouping;
     }
     {
                 Vars::processes["Output GRR"] = Processes::p_outputgrr;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_random_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);

//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputgrr] = grouping;
     }
     {
                 Vars::processes["Output Lapis"] = Processes::p_outputlapis;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_center_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputlapis] = grouping;
     }     {
                 Vars::processes["Output MDR"] = Processes::p_outputmdr;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputmdr] = grouping;
     }
     {
                 Vars::processes["Output PDT2"] = Processes::p_outputpdt2;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputpdt2] = grouping;
     }
     {
                 Vars::processes["Output PED"] = Processes::p_outputped;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputped] = grouping;
     }
     {
                 Vars::processes["Output BIN"] = Processes::p_outputbin;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputbin] = grouping;
     }
     {
                 Vars::processes["Output Phase"] = Processes::p_outputphase;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_trios);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputphase] = grouping;
     }
     {
                 Vars::processes["Output PowerMarker"] = Processes::p_outputpowermarker;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_center_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);

//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputpowermarker] = grouping;
     }
     {
                 Vars::processes["Output QTDT"] = Processes::p_outputqtdt;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputqtdt] = grouping;
     }
     {
                 Vars::processes["Output Structure"] = Processes::p_outputstructure;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_parents_only);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_random_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_stratification_options);

//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputstructure] = grouping;
     }
     {
                 Vars::processes["Output Superlink"] = Processes::p_outputsuperlink;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_penetrance_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputsuperlink] = grouping;
     }
     {
                 Vars::processes["Output Eigenstrat"] = Processes::p_outputeigenstrat;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_penetrance_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputeigenstrat] = grouping;
     }
     {
                 Vars::processes["Output TPED"] = Processes::p_outputtped;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_allele_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_sample_options);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_excluded);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_output_disabled);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_dummy_missing_parents);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_zero_incomplete_trio_ids);
                 grouping[OptionsGroup::og_outputs].push_back(OptionsItem::oi_rem_missing_parents);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
//                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_outputtped] = grouping;
     }
    {
                 Vars::processes["Concordance"] = Processes::p_concordance;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_inputs].push_back(OptionsItem::oi_pedinput);
                 grouping[OptionsGroup::og_inputs].push_back(OptionsItem::oi_tpedinput);
                 grouping[OptionsGroup::og_inputs].push_back(OptionsItem::oi_bininput);
                 grouping[OptionsGroup::og_inputs].push_back(OptionsItem::oi_mdrinput);
                 Vars::process_options_map[Processes::p_concordance] = grouping;
     }
    {
                 Vars::processes["Mitochondrial Error"] = Processes::p_mitoerror;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_sample_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_sample_threshold_max);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_min);
                 grouping[OptionsGroup::og_thresholds].push_back(OptionsItem::oi_locus_threshold_max);
                 Vars::process_options_map[Processes::p_mitoerror] = grouping;
     }
    {
                 Vars::processes["Plato Filter Processing"] = Processes::p_filterprocess;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_inputs].push_back(OptionsItem::oi_filterprocess_input);
                 Vars::process_options_map[Processes::p_filterprocess] = grouping;
     }
    {
                 Vars::processes["Deletion Detection"] = Processes::p_deletion;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_deletion);
                 Vars::process_options_map[Processes::p_deletion] = grouping;
     }
    {
                 Vars::processes["Homozygous"] = Processes::p_homozygous;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_homozygous);
                 Vars::process_options_map[Processes::p_homozygous] = grouping;
     }
    {
                 Vars::processes["LD"] = Processes::p_ld;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_ld);
                 Vars::process_options_map[Processes::p_ld] = grouping;
     }
    {
                 Vars::processes["IBS"] = Processes::p_ibs;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_ibs);
                 Vars::process_options_map[Processes::p_ibs] = grouping;
     }
    {
                 Vars::processes["Kinship"] = Processes::p_kinship;
                 map<int, vector<int> > grouping;
                 Vars::process_options_map[Processes::p_kinship] = grouping;
     }
    {
                 Vars::processes["FST"] = Processes::p_fst;
                 map<int, vector<int> > grouping;
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_fst);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_all);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_all_children);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_unaffected_spouses);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_unknown_spouses);
                 grouping[OptionsGroup::og_other].push_back(OptionsItem::oi_random_child);

                 Vars::process_options_map[Processes::p_fst] = grouping;
     }


    {
		 Vars::processes[""] = -1;
     }
 }
#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
