/*
 * Vars.h
 *
 *  Created on: Jun 25, 2010
 *      Author: cozartc
 */

#ifndef VARS_H_
#define VARS_H_
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <math.h>
#include <cmath>
//#include <sqlite3.h>
//#include <QMutex>
#include <map>
using namespace std;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

class Process;

typedef map<string, vector<Process*> > BATCHES;

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
		p_sampleefficiency,
		p_allelefrequency,
		p_concordance,
		p_familyefficiency,
		p_gendererror,
		p_hwe,
		p_mendelianerror,
		p_mitoerror,
		p_chisquare,
		p_cmh,
		p_deletion,
		p_homozygous,
		p_ld,
		p_logreg,
		p_linreg,
		p_tdt,
		p_mars,
		p_filterprocess,
		p_logicreg,
		p_mdr,
		p_mdrpdt,
		p_ibs,
		p_kinship,
		p_fst,
		p_clusermising,
		p_outputbeagle,
		p_outputeigenstrat,
		p_outputfbat,
		p_outputgrr,
		p_outputlapis,
		p_outputmdr,
		p_outputpdt2,
		p_outputped,
		p_outputphase,
		p_outputpowermarker,
		p_outputqtdt,
		p_outputstructure,
		p_outputsuperlink,
		p_outputtped,
		p_outputbin
	};
};

class OptionsGroup{
public:
	enum {
	og_thresholds,
	og_locusfilters,
	og_samplefilters,
	og_outputs,
	og_other,
	og_filters,
	og_covariates,
	og_inputs
	};
};

class OptionsItem{
public:
	enum{
	//og_thresholds
	oi_locus_threshold_min,
	oi_locus_threshold_max,
	oi_sample_threshold_min,
	oi_sample_threshold_max,
	oi_pedigree_threshold_min,
	oi_pedigree_threshold_max
	};

	enum{
        //og_locusfilters/og_filters
            oi_removemono,
            oi_removehet,
            oi_filter_overall,
            oi_freq_file,
            oi_zero,
            oi_zero_l2,
            oi_zero_l2_fams,
            oi_chrom,
            oi_bpmin,
            oi_bpmax,
            oi_bpspace
	};

	enum{
        //og_samplefilters/og_filters
	};

        enum{
            //og_inputs
            oi_pedinput,
            oi_tpedinput,
            oi_bininput,
            oi_mdrinput,
            oi_filterprocess_input
        };

        enum{
        //og_other
            oi_all,
            oi_all_children,
            oi_unaffected_spouses,
            oi_unknown_spouses,
            oi_random_child,
            oi_casecontrol,
            oi_gender,
            oi_parental,
            oi_deletion,
            oi_hwe,
            oi_homozygous,
            oi_ld,
            oi_logreg,
            oi_linreg,
            oi_xchr,
            oi_clusters,
            oi_cmh,
            oi_transform,
            oi_hwept,
            oi_prevalence,
            oi_disease,
            oi_ibs,
            oi_fst
        };

	enum{
        //og_outputs
            oi_allele_options,
            oi_allele_1234,
            oi_allele_acgt,
            oi_allele_custom,
            oi_sample_options,
            oi_output_disabled,
            oi_output_excluded,
            oi_zero_disabled,
            oi_zero_excluded,
            oi_parents_only,
            oi_trios,
            oi_stratification_options,
            oi_rem_missing_parents,
            oi_dummy_missing_parents,
            oi_zero_incomplete_trio_ids,
            oi_penetrance_options,
            oi_center_options,
            oi_disease_options,
            oi_random_options
	};

	enum{
	//og_covariates
            oi_covariates
	};

};

class Equality{
public:
    enum{
        greater_than,
        less_than,
        greater_equal_than,
        less_equal_than,
        equal_to
    };
};
//};

class Vars {
public:

	static string BLANK;
	static string BATCH;
	static string PLATO_GUI_GLADE;
	static string PLATO_OPTIONS_GLADE;
	static string PLATO_MODELS_GLADE;
//	static map<string, string> PLATO_DESIGN_FILES;
	static string LOCUS_FILTER_GUI;
	static string RESOURCES_DIR;
	static string REQUIRE_TOKEN;
	static string DISABLE_TOKEN;
	static string PLOT_GENOME;
	static string PLOT_SAMPLE;
	static string PLOT_PEDIGREE;
	static string PLOT_IMPORT_FILE;
	static string LOCUS_LOCATION;
	static string LOCUS_TABLE;
	static string SAMPLE_TABLE;
	static string PEDIGREE_TABLE;
	static int FILE_TYPE;
	static int TABLE_TYPE;
	static int IMPORT_TABLE_TYPE;
	//static GMutex *mutex_Project;
	//static QMutex mutex_Project;
	static double COVTRAITMISSINGDOUB;

	static map<string, int> processes;
	static map<int, map<int, vector<int> > > process_options_map;

	static void initialize();
	static vector<double> distributionNegLog(int count);
};
#ifdef PLATOLIB
};//end namespace PlatoLib
#endif

#endif /* VARS_H_ */
