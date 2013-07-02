#ifndef STEPOPTIONS_H
#define STEPOPTIONS_H

#include <string>
#include <vector>
#include <map>
#include "Marker.h"
#include "General.h"
#include "DataSet.h"
#include "Options.h"
#include "MethodException.h"


using namespace std;

namespace Methods{
class StepOptions {
	private:
		string options;
		double default_missing_value;
		bool dothresh_markers_min;
		bool dothresh_markers_max;
		double thresh_markers_min;
		double thresh_markers_max;
		bool dothresh_samples_min;
		bool dothresh_samples_max;
		double thresh_samples_min;
		double thresh_samples_max;
		bool dothresh_families_min;
		bool dothresh_families_max;
		double thresh_families_min;
		double thresh_families_max;
		bool dochrom;
		int chrom;
		bool trio;
		bool me_zero;
		bool me_zero_l2;
		bool me_zero_l2_fams;
		bool dobp_min;
		bool dobp_max;
		int bp_min;
		int bp_max;
		bool dobp_space;
		int bp_space;
		string disease;
		string struct_strat_file;
		bool parents_only;
		bool rem_missing_parents;
		bool dummy_missing_parents;
		bool zero_incomplete_trio_ids;
		bool dummy_incomplete_parent_ids;
		int homozyg_zeros;
		bool homozyg_raw;
		int homozyg_span;
		int homozyg_min_samp;
		bool homozyg_wgha;
		bool homozyg_seq_test;
		double homozyg_seq_val;
		int homozyg_permute_count;
		bool homozyg_permute;
		bool ld_pairwise;
		bool ld_vif;
		bool ld_regular;
		bool ld_chop;
		bool ld_calc_only;
		int ld_window;
		int ld_window_kb;
		int ld_step;
		double ld_threshold;
		double ld_pw_threshold;
		bool filter_overall;
		bool filter_file;
		bool center_codes;
		string center_file;
		map<string, string> center_map;
		bool stratification;
		map<string, int> strat_map;
		bool penetrance_codes;
		string penetrance_file;
		map<string, float> penetrance_map;
		bool unaff_spouses_only;
		bool unk_spouses;
		bool founders_only;
		bool random_child;
		bool all_children;
		bool all;
		bool deletion;
		int deletion_span;
		bool parental;
		bool gender;
		bool casecontrol;
		bool do_random_markers;
		bool do_random_repeat;
		bool do_sets;
		bool do_random_samples;
		bool do_sample_sets;
		bool rand_samps_repeat;
		vector<float> percent_samps;
		float percent_cases;
		float percent_controls;
		int rand_samps;
		int sets_samps;
		int random_markers;
		int sets;
		bool rm_mono;
		bool rm_het_only;
		bool inc_excluded_samples;
		bool inc_disabled_samples;
		bool zero_excluded;
		bool zero_disabled;
		string out;
		string override_out;

		//frequencies
		map<string, float> frequencies;
		map<string, float> group_frequencies;
		bool do_group_freq;
		string group_frequency_file;

		//cluster missing
		int max_cluster_size;
		int max_cluster_N;
		bool cluster_on_pheno;
		bool cluster_on_mcc;
		bool cluster_selcon;

		//plato filtering
		bool has_config;
		string config;

		//covariates & traits & loci (loci for modeling)
		bool do_covs;
		bool do_covs_number;
		bool do_covs_name;
		bool do_gxe, do_gxe_number, do_gxe_name;
		bool do_covars_only;
		string cov_file;
		vector<string> cov_map;
		vector<string> cov_use;
		vector<string> gxe_use;
		map<string, vector<double> > covs;
		bool do_traits;
		bool do_traits_number;
		bool do_traits_name;
		string trait_file;
		vector<string> trait_map;
		vector<string> trait_use;
		map<string, vector<double> > traits;
		bool do_covs_file;
		bool do_traits_file;
		string covar_missing;
		string trait_missing;
		bool alternate_pheno;
		int pheno_loc;
		int pheno_missing;
		string pheno_name;
		vector<string> loci_use;

		//groups
		bool do_group_file;
		string group_file;
		map<string, vector<Sample*> > groups;
		map<Sample*, string> sample_groups;
		map<string, vector<Family*> > group_families;

		bool do_cluster_file;
		string cluster_file;
		map<string, int> cluster_map;
		map<string, vector<Sample*> > clusters;
		map<Sample*, int> sample_clusters;
		map<string, string> sample_clusters_string;
		bool do_pedfile;
		bool do_mapfile;
		bool do_tpedfile;
		bool do_tfamfile;
		bool do_binfile;
		bool do_mdrfile;
		bool do_mdrpedfile;
		bool do_mdrmapfile;
		bool do_lgenfile;
		bool do_referencefile;
		string ped_file;
		string map_file;
		string tped_file;
		string tfam_file;
		string mdr_file;
		string mdr_map_file;
		string mdr_ped_file;
		string bin_prefix;
		string lgen_file;
		string reference_file;
		double confidence_interval;
		bool allele1234;
		bool alleleACGT;
		bool allele12;
		bool allelecustom;
		map<string, string> custom_alleles;
		bool hwept;
		double prevalence;
		int perms;

		//outputped
		bool allele_as_snp;

		//outputmdr
		bool no_rules;
		bool mdr_gui_output;
		bool mdr_pedigree_output;

		//Linear Regresssion method parameters
		string linr_modType;
		bool linr_condition;
		string linr_condition_string;
		string linr_condition_file;
		bool linr_interaction;
		bool linr_no_main_snp;
		vector<int> linr_condition_list;

		// Logistic Regression method parameters
    bool lr_fullInteraction;
    bool lr_includeInteractions;
    string lr_modType;
    unsigned int lr_maxIterations;
		// Conditional Logistic Regression method parameters
	bool cond_lr;
    bool cond_lr_fullInteraction;
    bool cond_lr_includeInteractions;
    string cond_lr_modType;
    unsigned int cond_lr_maxIterations;
    bool mdr_set_only_threshold;
    // Uncertainty Coefficient method parameters
    string uncert_coeff_total_type;
    // Likelihood Ratio method parameters
    string llr_total_type;
    // Odds Ratio method parameters
    string oddsratio_total_type;
    // NMI method parameters
    string nmi_total_type;
    bool nmi_transposed;

    // true if map file contains referent allele
    bool map_includes_ref_allele;

    //CMH
    bool cmh_ijk;
    bool cmh_ordinal;
    bool cmh_22k;
    bool breslowday;

    //MODEL REQUIREMENTS (PLINK)
	bool sex_effect;
	bool auto_sex_effect;
    int xchr_model;
    bool qfam_total;
    bool qfam_between;
    bool qfam_within1;
    bool qfam_within2;
    bool qfam_adaptive;
    double vif_threshold;

	//MULT COMPARISON
	double fixed_lambda;
	bool do_fixed_lambda;

    // MDRPDT
    int mdrpdt_ptests;
    int mdrpdt_seed;
    int mdrpdt_xv;
    int mdrpdt_mincombo;
    int mdrpdt_maxcombo;

    //transforms
    bool sqrt_transform;
    bool boxcox_transform;
    vector<double> boxcox_covar_lambdas;
    vector<double> boxcox_trait_lambdas;
    bool log_transform;
    vector<string> to_transform;

    //multiple comparisions
    bool mult_compare;

    //ibs
    bool do_ibs_pairs;
    bool do_ibs_trio_pairs;
    bool do_ibs_all_pairs;
    bool do_ibs_all_trio_pairs;
    bool do_ibs_trio_transmissions;
    bool do_ibs_trios_raw;
    map<string, vector<string> > ibs_pairs;
    map<string, vector<string> > ibs_trio_pairs;
    string ibs_file;
    string ibs_trios_file;

    //concordance
    bool unique_id;
	bool inc_missing;

    //MARS
    int mars_maxterms;
    int mars_k;
    int mars_degree;
    int mars_nprune;

    //sample bprange filter
    string sampbprange_file;
    vector<Sample*> sampbprange_samples;
    vector<vector<Marker*> > sampbprange_markers;

    //epistasis
    bool epi_quickscan;
    bool epi_fast;
    bool epi_caseonly;
    double epi_alpha1;
    double epi_alpha2;
    int epi_caseonly_kb_gap;
    map<string, vector<string> > epi_sets;
    bool epi_set_by_set;
    bool do_epi_sets;
    string epi_sets_filename;
    bool epi_filter;
    
    // regression interaction testing
    double lrt_pval;
    bool lrt_pval_filter;

    //athena/biofilter comparison input
    string bio_comparison_file, gxe_list_file;
    int bio_offset_begin;
    int bio_offset_end;
    bool bio_file_binary;
    map<string, vector<string> > bio_pairs, gxe_pairs;

    //eigenstrat options
    bool qtl;
    bool ancestry;

    //autosome only
    bool autosome_only;

    //synthesis view
    bool output_synthview;


		enum Argument{
			s_thresh_markers_max,
			s_thresh_markers_min,
			s_thresh_samples_max,
			s_thresh_samples_min,
			s_thresh_families_max,
			s_thresh_families_min,
			s_ped_file,
			s_map_file,
			s_tped_file,
			s_tfam_file,
			s_mdr_file,
			s_lgen_file,
			s_mdr_ped_file,
			s_mdr_map_file,
			s_mdr_gui_output,
			s_bin_prefix,
			s_hwept,
			s_chrom,
			s_trio,
			s_prevalence,
			s_me_zero,
			s_me_zero_l2,
			s_me_zero_l2_fams,
			s_bp_min,
			s_bp_max,
			s_bp_space,
			s_disease,
			s_struct_strat_file,
			s_parents_only,
			s_homozyg_raw,
			s_homozyg_wgha,
			s_homozyg_zeros,
			s_homozyg_span,
			s_homozyg_min_samp,
			s_homozyg_seq,
			s_homozyg_permute,
			s_ld_pairwise,
			s_ld_regular,
			s_ld_chop,
			s_ld_calc_only,
			s_ld_window,
			s_ld_window_kb,
			s_ld_step,
			s_filter_overall,
			s_filter_file,
			s_center_file,
			s_penetrance_file,
			s_rm_mono,
			s_rm_het_only,
			s_unaff_spouses_only,
			s_unk_spouses,
			s_founders_only,
			s_random_child,
			s_all_children,
			s_all,
			s_deletion,
			s_deletion_span,
			s_parental,
			s_gender,
			s_casecontrol,
			s_random_markers,
			s_random_repeat,
			s_sets,
			s_rand_samps,
			s_rand_samps_repeat,
			s_sets_samps,
			s_percent_samps,
			s_percent_cases,
			s_percent_controls,
			s_out,
			s_override_out,
			s_inc_excluded_samples,
			s_inc_disabled_samples,
			s_zero_excluded,
			s_zero_disabled,
			//covariates & traits
			s_covar_file,
			s_covars_name,
			s_covars_number,
			s_trait_file,
			s_traits_name,
			s_traits_number,
			s_covar_missing,
			s_trait_missing,
			s_pheno,
			s_pheno_index,
			s_pheno_missing,
			s_loci,
			s_gxe_number,
			s_gxe_name,
			s_run_covars_only,

			//insert cluster missing here
			s_max_cluster_size,
			s_max_cluster_N,
			s_cluster_on_pheno,
			s_cluster_on_mcc,
			s_cluster_selcon,

			//plato filtering
			s_config,

			s_group_file,
			s_cluster_file,
			s_ci,
			s_allele1234,
			s_alleleACGT,
			s_allele12,
			s_allelecustom,
			s_perms,
			//linear regression
			s_linr_mod_type,
			s_linr_condition,
			s_linr_condition_file,
			s_linr_interaction,
			s_linr_no_main_snp,
			//logistic regression
			s_lr_full_interaction,
			s_lr_include_interaction,
			s_lr_no_interaction,
			s_lr_reduced_interaction,
			s_lr_mod_type,
			s_lr_max_iter,
			s_cond_lr_full_interaction,
			s_cond_lr_include_interaction,
			s_cond_lr_no_interaction,
			s_cond_lr_reduced_interaction,
			s_cond_lr_mod_type,
			s_cond_lr_max_iter,
			//mdr?
			s_mdr_only_set_thresh,
			s_uncert_coeff_total,
			s_llr_total,
			s_or_total,
			s_nmi_total,
			s_nmi_transposed,
			s_rem_missing_parents,
			s_dummy_missing_parents,
			s_zero_incomplete_trio_ids,
			s_dummy_incomplete_parent_ids,
			s_ref_allele_included,
			//CMH
			s_cmh_ijk,
			s_cmh_22k,
			s_cmh_ordinal,
			//MODEL REQS (PLINK)
			s_xchr_model,
			s_sex_effect,
			s_auto_sex_effect,
			s_qfam_total,
			s_qfam_between,
			s_qfam_within1,
			s_qfam_within2,
			s_qfam_adaptives,
			s_vif_threshold,
			//MULT COMPARISON
			s_lambda,
			s_mult_compare,
			//MDRPDT
			s_mdrpdt_ptests,
			s_mdrpdt_randseed,
			s_mdrpdt_xvcount,
			s_mdrpdt_mincombo,
			s_mdrpdt_maxcombo,
			//data transforms
			s_sqrt_transform,
			s_boxcox_transform,
			s_log_transform,
			//ibs
			s_ibs_pairs,
			s_ibs_trio_pairs,
			s_ibs_all_pairs,
			s_ibs_all_trio_pairs,
			s_ibs_trios_raw,
			s_ibs_trio_trans,
			//outputmdr
			s_no_rules,
			s_mdr_pedigree_output,
			//outputped
			s_allele_as_snp,
			//concordance
			s_unique_id,
			s_inc_missing,
			//MARS
			s_mars_maxterms,
			s_mars_degree,
			s_mars_k,
			s_mars_nprune,
			//samp bprange filter
			s_sampbprange_filter,

			//group frequencies
			s_group_freq,

			//epistasis
		    s_epi_quickscan,
		    s_epi_fast,
		    s_epi_caseonly,
		    s_epi_alpha1,
		    s_epi_alpha2,
		    s_epi_caseonly_kb_gap,
		    s_epi_sets,
		    s_epi_set_by_set,
		    s_do_epi_sets,
		    s_epi_sets_filename,
		    s_epi_filter,

		    //biofilter input
		    s_bio_comparison_file,
		    s_bio_offset_begin,
		    s_bio_offset_end,
		    s_bio_file_binary,

				s_gxe_file,

		    //autosome_only
		    s_autosome_only,

		    //eigenstrat
		    s_qtl,
		    s_ancestry,
		    
		    // regression interaction testing
		    s_lrt_pval,

		    //output synthesis view
		    s_output_synthview,

		};
		map<string, Argument> s_ArgVals;
		map<string, Argument> s_subsetVals;

	public:
		StepOptions(){
			dobp_space = false;
			default_missing_value = -99999;
			do_pedfile = false;
			do_mapfile = false;
			do_tpedfile = false;
			do_tfamfile = false;
			do_binfile = false;
			do_mdrfile = false;
			do_mdrpedfile = false;
			do_mdrmapfile = false;
			do_lgenfile = false;
			mdr_gui_output = false;
			allele1234 = false;
			alleleACGT = false;
			allele12 = false;
			allelecustom = false;
			prevalence = 0.25;
			ped_file = "";
			map_file = "";
			tped_file = "";
			tfam_file = "";
			mdr_file = "";
			mdr_ped_file = "";
			mdr_map_file = "";
			bin_prefix = "";
			do_covs = false;
			do_gxe = false;
			do_covars_only = false;
			do_gxe_name = false;
			do_gxe_number =false;
			do_traits = false;
			do_covs_file = false;
			do_traits_file = false;
			do_traits_name = false;
			do_covs_name = false;
			do_group_file = false;
			group_file = "";
			cov_file = "";
			trait_file = "";
			alternate_pheno = false;
			pheno_loc = -1;
			pheno_name = "";
			pheno_missing = opts::_PHENO_MISS_;

			hwept = false;
			covar_missing = "-99999";
			trait_missing = "-99999";
			dothresh_markers_min = false;
			dothresh_markers_max = false;
			dothresh_samples_min = false;
			dothresh_samples_max = false;
			dothresh_families_min = false;
			dothresh_families_max = false;
			center_codes = false;
			center_file = "";
			penetrance_codes = false;
			penetrance_file = "";
			options = "";
			rm_mono = false;
			rm_het_only = false;

			//group frequencies
			do_group_freq = false;
			group_frequency_file = "";

			//cluster missing
			max_cluster_size = 0;
			max_cluster_N = -1;
			cluster_on_mcc = false;
			cluster_on_pheno = false;
			cluster_selcon = false;

			//plato filtering
			has_config = false;
			config = "";

			//outputmdr
			no_rules = false;
			mdr_pedigree_output = false;
			//outputped
			allele_as_snp = false;
			//CMH
			cmh_ijk = false;
			cmh_22k = false;
			cmh_ordinal = false;
			breslowday = false;

			//MULT COMPARISON
			fixed_lambda = 0;
			do_fixed_lambda = false;
			mult_compare = false;

			//transforms
			sqrt_transform = false;
			boxcox_transform = false;
			log_transform = false;

			//ibs
			do_ibs_pairs = false;
			do_ibs_all_pairs = false;
			do_ibs_all_trio_pairs = false;
			do_ibs_trio_transmissions = false;
			do_ibs_trio_pairs = false;
			do_ibs_trios_raw = false;
			ibs_file = "";
			ibs_trios_file = "";

			//MARS
			mars_maxterms = 21;
			mars_degree = 1;
			mars_nprune = 2147483647;	//initialize to max int value to simulate infinity

			//concordance
			unique_id = false;
			inc_missing = false;

			//samp bprange filter
			sampbprange_file = "";
			sampbprange_samples.clear();
			sampbprange_markers.clear();

			//epistasis
			epi_quickscan = false;
			epi_fast = false;
			epi_caseonly = false;
			epi_alpha1 = 0.0001;
			epi_alpha2 = 0.01;
			epi_caseonly_kb_gap = 0;
			epi_set_by_set = false;
			do_epi_sets = false;
			epi_sets_filename = "";
			epi_filter = true;

			//biofilter input
			bio_comparison_file = "";
			bio_offset_begin = -1;
			bio_offset_end = -1;
			bio_file_binary = false;

			gxe_list_file = "";

			//autosome_only
			autosome_only = false;

      // regression interactions
      lrt_pval = 0.05;

			//eigenstrat
			qtl = false;
			ancestry = false;

			//output synthesis view
			output_synthview = false;

			//output synthesis view
			s_ArgVals["-output-synthview"] = s_output_synthview;

			//eigenstrat
			s_ArgVals["-qtl"] = s_qtl;
			s_ArgVals["-ancestry"] = s_ancestry;

			s_ArgVals["-thresh-markers-max"] = s_thresh_markers_max;
			s_ArgVals["-thresh-markers-min"] = s_thresh_markers_min;
			s_ArgVals["-thresh-samples-max"] = s_thresh_samples_max;
			s_ArgVals["-thresh-samples-min"] = s_thresh_samples_min;
			s_ArgVals["-thresh-families-max"] = s_thresh_families_max;
			s_ArgVals["-thresh-families-min"] = s_thresh_families_min;
			s_ArgVals["-allele1234"] = s_allele1234;
			s_ArgVals["-alleleACGT"] = s_alleleACGT;
			s_ArgVals["-allele12"] = s_allele12;
			s_ArgVals["-allele-custom"] = s_allelecustom;
			s_ArgVals["-rm-mono"] = s_rm_mono;
			s_ArgVals["-rm-het-only"] = s_rm_het_only;
			s_ArgVals["-trio"] = s_trio;
			s_ArgVals["-chrom"] = s_chrom;
			s_ArgVals["-zero"] = s_me_zero;
			s_ArgVals["-zero-l2"] = s_me_zero_l2;
			s_ArgVals["-zero-l2-fams"] = s_me_zero_l2_fams;
			s_ArgVals["-bp-min"] = s_bp_min;
			s_ArgVals["-bp-max"] = s_bp_max;
			s_ArgVals["-bp-space"] = s_bp_space;
			s_ArgVals["-hwept"] = s_hwept;
			s_ArgVals["-prevalence"] = s_prevalence;
			s_ArgVals["-disease"] = s_disease;
			s_ArgVals["-strat-file"] = s_struct_strat_file;
			s_ArgVals["-parents-only"] = s_parents_only;
			s_ArgVals["-homozyg-raw"] = s_homozyg_raw;
			s_ArgVals["-homozyg-zeros"] = s_homozyg_zeros;
			s_ArgVals["-homozyg-wgha"] = s_homozyg_wgha;
			s_ArgVals["-homozyg-span"] = s_homozyg_span;
			s_ArgVals["-homozyg-min-samp"] = s_homozyg_min_samp;
			s_ArgVals["-homozyg-seq-prob"] = s_homozyg_seq;
			s_ArgVals["-homozyg-permute"] = s_homozyg_permute;
			s_ArgVals["-ld-pairwise"] = s_ld_pairwise;
			s_ArgVals["-ld-vif"] = s_ld_regular;
			s_ArgVals["-ld-chop"] = s_ld_chop;
			s_ArgVals["-ld-calc-only"] = s_ld_calc_only;
			s_ArgVals["-ld-window"] = s_ld_window;
			s_ArgVals["-ld-window-kb"] = s_ld_window_kb;
			s_ArgVals["-ld-step"] = s_ld_step;
			s_ArgVals["-filter-overall"] = s_filter_overall;
			s_ArgVals["-filter-file"] = s_filter_file;
			s_ArgVals["-center-file"] = s_center_file;
			s_ArgVals["-penetrance-file"] = s_penetrance_file;
			s_ArgVals["-unaff-spouses-only"] = s_unaff_spouses_only;
			s_ArgVals["-unk-spouses"] = s_unk_spouses;
			s_ArgVals["-founders-only"] = s_founders_only;
			s_ArgVals["-all-children"] = s_all_children;
			s_ArgVals["-random-child"] = s_random_child;
			s_ArgVals["-all"] = s_all;
			s_ArgVals["-deletion"] = s_deletion;
			s_ArgVals["-deletion-span"] = s_deletion_span;
			s_ArgVals["-parental"] = s_parental;
			s_ArgVals["-gender"] = s_gender;
			s_ArgVals["-casecontrol"] = s_casecontrol;
			s_ArgVals["-random-markers"] = s_random_markers;
			s_ArgVals["-random-repeat"] = s_random_repeat;
			s_ArgVals["-sets"] = s_sets;
			s_ArgVals["-rand-samps"] = s_rand_samps;
			s_ArgVals["-rand-samps-repeat"] = s_rand_samps_repeat;
			s_ArgVals["-sets-samps"] = s_sets_samps;
			s_ArgVals["-percent-samps"] = s_percent_samps;
			s_ArgVals["-percent-cases"] = s_percent_cases;
			s_ArgVals["-percent-controls"] = s_percent_controls;
			s_ArgVals["-out"] = s_out;
			s_ArgVals["-override-out"] = s_override_out;
			s_ArgVals["-inc-excluded-samples"] = s_inc_excluded_samples;
			s_ArgVals["-inc-disabled-samples"] = s_inc_disabled_samples;
			s_ArgVals["-zero-excluded"] = s_zero_excluded;
			s_ArgVals["-zero-disabled"] = s_zero_disabled;
			s_ArgVals["-covar-file"] = s_covar_file;
			s_ArgVals["-covars-name"] = s_covars_name;
			s_ArgVals["-covars-number"] = s_covars_number;
			s_ArgVals["-gxe-name"] = s_gxe_name;
			s_ArgVals["-gxe-number"] = s_gxe_number;			
			s_ArgVals["-trait-file"] = s_trait_file;
			s_ArgVals["-traits-name"] = s_traits_name;
			s_ArgVals["-traits-number"] = s_traits_number;
			s_ArgVals["-covar-missing"] = s_covar_missing;
			s_ArgVals["-run-covars-only"] = s_run_covars_only;
			s_ArgVals["-trait-missing"] = s_trait_missing;
			s_ArgVals["-loci"] = s_loci;
			s_ArgVals["-pheno"] = s_pheno;
			s_ArgVals["-pheno-index"] = s_pheno_index;
			s_ArgVals["-pheno-missing"] = s_pheno_missing;

			s_ArgVals["-group-file"] = s_group_file;
			s_ArgVals["-cluster-file"] = s_cluster_file;
			s_ArgVals["-ped"] = s_ped_file;
			s_ArgVals["-map"] = s_map_file;
			s_ArgVals["-mdr"] = s_mdr_file;
			s_ArgVals["-mdrped"] = s_mdr_ped_file;
			s_ArgVals["-mdrmap"] = s_mdr_map_file;
			s_ArgVals["-mdr-gui-output"] = s_mdr_gui_output;
			s_ArgVals["-tped"] = s_tped_file;
			s_ArgVals["-tfam"] = s_tfam_file;
			s_ArgVals["-bin-input"] = s_bin_prefix;
			s_ArgVals["-perms"] = s_perms;
			s_ArgVals["-ci"] = s_ci;
			s_ArgVals["-rem-missing-parents"] = s_rem_missing_parents;
			s_ArgVals["-dummy-missing-parents"] = s_dummy_missing_parents;
			s_ArgVals["-zero-incomplete-trio-ids"] = s_zero_incomplete_trio_ids;
			s_ArgVals["-dummy-incomplete-parent-ids"] = s_dummy_incomplete_parent_ids;
			//concordance
			s_ArgVals["-unique-id"] = s_unique_id;
			s_ArgVals["-inc-missing"] = s_inc_missing;

			//group frequencies
			s_ArgVals["-group-freq"] = s_group_freq;

			//cluster missing
			s_ArgVals["-max-cluster-size"] = s_max_cluster_size; //undocumented
			s_ArgVals["-max-cluster-N"] = s_max_cluster_N;       //undocumented
			s_ArgVals["-cluster-on-pheno"] = s_cluster_on_pheno; //undocumented
			s_ArgVals["-cluster-on-mcc"] = s_cluster_on_mcc;     //undocumented
			s_ArgVals["-cluster-selcon"] = s_cluster_selcon;     //undocumented

			//plato filtering
			s_ArgVals["-config"] = s_config;

			//Linear regression options
			s_ArgVals["-linr-model-type"] = s_linr_mod_type;
			s_ArgVals["-linr-condition"] = s_linr_condition;
			s_ArgVals["-linr-contition-file"] = s_linr_condition_file;
			s_ArgVals["-linr-no-main-snp"] = s_linr_no_main_snp;
			s_ArgVals["-linr-interaction"] = s_linr_interaction;
			//outputmdr
			s_ArgVals["-no-rules"] = s_no_rules;
			s_ArgVals["-mdr-pedigree"] = s_mdr_pedigree_output;
			//outputped
			s_ArgVals["-allele-as-snp"] = s_allele_as_snp;     //undocumented
			// Logistic regression options
			s_ArgVals["-lr-full-interact"] = s_lr_full_interaction;
			s_ArgVals["-lr-inc-interact"] = s_lr_include_interaction;
			s_ArgVals["-lr-model-type"] = s_lr_mod_type;
			s_ArgVals["-lr-max-iter"] = s_lr_max_iter;
			s_ArgVals["-lr-no-interact"] = s_lr_no_interaction;
			s_ArgVals["-lr-reduced-interact"] = s_lr_reduced_interaction;
			// Conditional logistic regression options
			s_ArgVals["-cond-lr-full-interact"] = s_cond_lr_full_interaction;
			s_ArgVals["-cond-lr-inc-interact"] = s_cond_lr_include_interaction;
			s_ArgVals["-cond-lr-model-type"] = s_cond_lr_mod_type;
			s_ArgVals["-cond-lr-max-iter"] = s_cond_lr_max_iter;
			// MARS
			s_ArgVals["-mars-maxterms"] = s_mars_maxterms;
			s_ArgVals["-mars-degree"] = s_mars_degree;
			s_ArgVals["-mars-nprune"] = s_mars_nprune;

			// MDR option
			s_ArgVals["-mdr-set-only-thresh"] = s_mdr_only_set_thresh;
			// Uncertainty Coefficient
			s_ArgVals["-s-uncert-coeff-total"] = s_uncert_coeff_total; //undocumented
			// Likelihood Ratio
			s_ArgVals["-s-llr-total"] = s_llr_total;                   //undocumented
			// Odds Ratio
			s_ArgVals["-s-or-total"] = s_or_total;                     //undocumented
			// NMI
			s_ArgVals["-s-nmi-total"] = s_nmi_total;                   //undocumented
			s_ArgVals["-s-nmi-tranposed"] = s_nmi_transposed;          //undocumented
			// flag for including referent alleles in map in the last column of map file
			s_ArgVals["-map-includes-ref"] = s_ref_allele_included;
			//CMH
			s_ArgVals["-cmh-22k"] = s_cmh_22k;
			s_ArgVals["-cmh-ijk"] = s_cmh_ijk;
			s_ArgVals["-cmh-ordinal"] = s_cmh_ordinal;
			//MODEL REQS (PLINK)
			s_ArgVals["-xchr-model"] = s_xchr_model;                  //undocumented
			s_ArgVals["-sex-effect"] = s_sex_effect;
			s_ArgVals["-no-sex-effect"] = s_auto_sex_effect;
			s_ArgVals["-qfam-total"] = s_qfam_total;            //undocumented
			s_ArgVals["-qfam-between"] = s_qfam_between;        //undocumented
			s_ArgVals["-qfam-within1"] = s_qfam_within1;        //undocumented
			s_ArgVals["-qfam-within2"] = s_qfam_within2;        //undocumented
			s_ArgVals["-qfam-adaptive"] = s_qfam_adaptives;     //undocumented
			s_ArgVals["-vif-threshold"] = s_vif_threshold;

			//transforms
			s_ArgVals["-sqrt-transform"] = s_sqrt_transform;
			s_ArgVals["-boxcox-transform"] = s_boxcox_transform;
			s_ArgVals["-log-transform"] = s_log_transform;

			// MDRPDT
			s_ArgVals["-mdrpdt-ptests"] = s_mdrpdt_ptests;            //undocumented
			s_ArgVals["-mdrpdt-randseed"] = s_mdrpdt_randseed;        //undocumented
			s_ArgVals["-mdrpdt-xvcount"] = s_mdrpdt_xvcount;          //undocumented
			s_ArgVals["-mdrpdt-mincombo"] = s_mdrpdt_mincombo;        //undocumented
			s_ArgVals["-mdrpdt-maxcombo"] = s_mdrpdt_maxcombo;        //undocumented

			//MULT COMPARISON
			s_ArgVals["-lambda"] = s_lambda;
			s_ArgVals["-mult-compare"] = s_mult_compare;

			//ibs
			s_ArgVals["-ibs-pairs"] = s_ibs_pairs;
			s_ArgVals["-ibs-trio-pairs"] = s_ibs_trio_pairs;
			s_ArgVals["-ibs-all-pairs"] = s_ibs_all_pairs;
			s_ArgVals["-ibs-all-trio-pairs"] = s_ibs_all_trio_pairs;
			s_ArgVals["-ibs-trios-raw"] = s_ibs_trios_raw;
			s_ArgVals["-ibs-trio-trans"] = s_ibs_trio_trans;

			//samp bprange filter
			s_ArgVals["-sample-bprange-filter"] = s_sampbprange_filter;

			//epistasis
			s_ArgVals["-epi-alpha1"] = s_epi_alpha1;
			s_ArgVals["-epi-alpha2"] = s_epi_alpha2;
			s_ArgVals["-epi-sets"] = s_epi_sets;
			s_ArgVals["-epi-set-by-set"] = s_epi_set_by_set;
			s_ArgVals["-epi-fast"] = s_epi_fast;

			//biofilter input
			s_ArgVals["-bio-snp-file"] = s_bio_comparison_file;
			s_ArgVals["-bio-file-begin"] = s_bio_offset_begin;
			s_ArgVals["-bio-file-end"] = s_bio_offset_end;
			s_ArgVals["-bio-file-binary"] = s_bio_file_binary;

			s_ArgVals["-gxe-list-file"] = s_gxe_file;

			//autosome only
			s_ArgVals["-auto-only"] = s_autosome_only;

      // regression interaction testing
      s_ArgVals["-lrt-pval"] = s_lrt_pval;

			s_subsetVals["-parental"] = s_parental;
			s_subsetVals["-gender"] = s_gender;
			s_subsetVals["-casecontrol"] = s_casecontrol;

			rem_missing_parents = false;
			dummy_missing_parents = false;
			zero_incomplete_trio_ids = false;
			dummy_incomplete_parent_ids = false;
			confidence_interval = 0.95;
			chrom = -1;
			trio = false;
			me_zero = false;
			me_zero_l2 = false;
			me_zero_l2_fams = false;
			dobp_min = false;
			dobp_max = false;
			dochrom = false;
			inc_excluded_samples = false;
			inc_disabled_samples = false;
			zero_excluded = false;
			zero_disabled = false;
			bp_space = 0;
			disease = "";
			struct_strat_file = "";
			stratification = false;
			parents_only = false;
			homozyg_raw = false;
			homozyg_zeros = -1;
			homozyg_wgha = false;
			homozyg_span = 50;
			homozyg_min_samp = 10;
			homozyg_seq_test = false;
			homozyg_seq_val = 1;
			homozyg_permute = false;
			homozyg_permute_count = 0;
			ld_pairwise = false;
			ld_vif = false;
			ld_chop = false;
			ld_step = 0;
			ld_window = 0;
			ld_window_kb = 1000;
			ld_threshold = 2;
			ld_pw_threshold = 1 - 1e-6;
			ld_calc_only = false;
			filter_overall = false;
			filter_file = false;
			unaff_spouses_only = false;
			unk_spouses = false;
			founders_only = false;
			random_child = false;
			all_children = false;
			all = false;
			deletion = false;
			deletion_span = 25;
			parental = false;
			gender = false;
			casecontrol = false;
			do_random_markers = false;
			do_random_repeat = false;
			do_sets = false;
			random_markers = 0;
			sets = 1;
			//random samples
			rand_samps = 0;
			rand_samps_repeat = false;
			sets_samps = 0;
			percent_samps.clear();
			percent_cases = 0;
			percent_controls = 0;

			out = convertString(options);
			override_out = "";
			perms = 1000;

			//linear regression default values
			linr_modType = "";
			linr_condition = false;
			linr_interaction = false;
			linr_no_main_snp = false;
			linr_condition_string = "";
			linr_condition_file = "";

			// Logistic Regression default values
			lr_fullInteraction = false;
			lr_includeInteractions = false;
			lr_modType = "ADDITIVE";
			lr_maxIterations = 20;
			// Conditional logistic regression default values
			cond_lr = false;
			cond_lr_fullInteraction = true;
			cond_lr_includeInteractions = true;
			cond_lr_modType = "ADDITIVE";
			cond_lr_maxIterations = 20;
			// MDR default options
			mdr_set_only_threshold = false;
			// Uncertainty Coefficient default values
			uncert_coeff_total_type = "ALLELE";
			// Odds Ratio default value
			oddsratio_total_type = "ALLELE";
			// NMI
			nmi_total_type = "ALLELE";
			nmi_transposed = false;
			// MDRPDT
			mdrpdt_ptests = 0;
			mdrpdt_seed = 7;
			mdrpdt_xv = 1;
			mdrpdt_maxcombo = 2;
			mdrpdt_mincombo = 1;

			//MODEL REQS (PLINK)
			xchr_model = 1;
			sex_effect = false;
			auto_sex_effect = false;
			qfam_total = false;
			qfam_between = false;
			qfam_within1 = false;
			qfam_within2 = false;
			qfam_adaptive = false;
			vif_threshold = 50;

			map_includes_ref_allele = false;
		};
		~StepOptions(){};

		void printOptions();

		void setUp(string s);
		string toString(){return options;};

		///get/set output synthesis view
		void setOutputSynthView(bool b){output_synthview = b;}
		bool getOutputSynthView(){return output_synthview;}
		bool doOutputSynthView(){return output_synthview;}

		///get/set eigenstrat
		void setQTL(bool b){qtl = b;}
		bool getQTL(){return qtl;}
		bool doQTL(){return qtl;}
		void setAncestry(bool b){ancestry = b;}
		bool getAncestry(){return ancestry;}
		bool doAncestry(){return ancestry;}

		///get/set autosome only
		void setAutosomeOnly(bool b){autosome_only = b;}
		bool getAutosomeOnly(){return autosome_only;}

		///get/set biofilter input files
		void readBioTextFile(string file);
		map<string, vector<string> > getBioPairs(){return bio_pairs;}
		void readGXETextFile(string file);
		map<string, vector<string> > getGXEPairs(){return gxe_pairs;}
		void setBioSnpFile(string s){bio_comparison_file = s;}
		string getBioSnpFile(){return bio_comparison_file;}
		void setGXEFile(string s){gxe_list_file=s;}
		string getGXEFile(){return gxe_list_file;}
		void setBioOffsetBegin(int s){bio_offset_begin = s;}
		int getBioOffsetBegin(){return bio_offset_begin;}
		void setBioOffsetEnd(int s){bio_offset_end = s;}
		int getBioOffsetEnd(){return bio_offset_end;}
		void setBioFileBinary(bool b){bio_file_binary = b;}
		bool getBioFileBinary(){return bio_file_binary;}

		///get/set epistasis
		void setEpiQuickscan(bool b){epi_quickscan = b;}
		bool doEpiQuickscan(){return epi_quickscan;}
		bool getEpiQuickscan(){return epi_quickscan;}
		void setEpiFast(bool b){epi_fast = b;}
		bool doEpiFast(){return epi_fast;}
		bool getEpiFast(){return epi_fast;}
		void setEpiCaseOnly(bool b){epi_caseonly = b;}
		bool doEpiCaseOnly(){return epi_caseonly;}
		bool getEpiCaseOnly(){return epi_caseonly;}
		void set_epi_alpha1(double d){epi_alpha1 = d;}
		void set_epi_alpha2(double d){epi_alpha2 = d;}
		double get_epi_alpha1(){return epi_alpha1;}
		double get_epi_alpha2(){return epi_alpha2;}
		void setEpiCaseOnlyKbGap(int v){epi_caseonly_kb_gap = v;}
		int getEpiCaseOnlyKbGap(){return epi_caseonly_kb_gap;}
		map<string, vector<string> > getEpiSets(){return epi_sets;}
		void readEpiSets(string);
		bool doEpiSets(){return do_epi_sets;}
		void setDoEpiSets(bool b){do_epi_sets = b;}
		bool doEpiSetBySet(){return epi_set_by_set;}
		void setDoEpiSetBySet(bool b){epi_set_by_set = b;}
		void setEpiSetsFilename(string s){epi_sets_filename = s;}
		string getEpiSetsFilename(){return epi_sets_filename;}
		void setEpiFilter(bool b){epi_filter = b;}
		bool getEpiFilter(){return epi_filter;}
		bool doEpiFilter(){return epi_filter;}

    //get/set interaction regression
    void setLRTPval(double val){lrt_pval=val;}
    double getLRTPval(){return lrt_pval;}
		void setLRTFilter(bool b){lrt_pval_filter = b;}
		bool doLRTFilter(){return lrt_pval_filter;}

		///get/set group frequencies
		void setDoGroupFreq(bool b){do_group_freq = b;}
		bool getDoGroupFreq(){return do_group_freq;}
		map<string, float> getGroupFreq(){return group_frequencies;}
		void setGroupFreq(map<string, float> gf){group_frequencies = gf;}
		string getGroupFreqFile(){return group_frequency_file;}
		void setGroupFreqFile(string f){group_frequency_file = f;}
		void readGroupFrequencies();

		///get/set concorddance options
		void setUniqueId(bool b){unique_id = b;}
		bool getUniqueId(){return unique_id;}
		void setIncMissing(bool b){inc_missing = b;}
		bool getIncMissing(){return inc_missing;}

		///get/set mars options
		void setMarsMaxTerms(int i){mars_maxterms = i;}
		int getMarsMaxTerms(){return mars_maxterms;}
		void setMarsDegree(int i){mars_degree = i;}
		int getMarsDegree(){return mars_degree;}
		void setMarsNPrune(int i){mars_nprune = i;}
		int getMarsNPrune(){return mars_nprune;}

		///get/set plato filtering filter-process
		string getFilterProcessConfig(){return config;}
		void setFilterProcessConfig(string s){config = s;}
		bool hasFilterProcessConfig(){return has_config;}
		void hasFilterProcessConfig(bool b){has_config = b;}

		///get/set outputmdr
		void set_no_rules(bool v){no_rules = v;}
		bool get_no_rules(){return no_rules;}
		void setMDRPedigreeOutput(bool v){mdr_pedigree_output = v;}
		bool getMDRPedigreeOutput(){return mdr_pedigree_output;}

		///get/set outputped
		void set_allele_as_snp(bool v){allele_as_snp = v;}
		bool get_allele_as_snp(){return allele_as_snp;}

		///get/set transforms
		void set_sqrt_transform(bool v){sqrt_transform = v;}
		void set_boxcox_transform(bool v){boxcox_transform = v;}
		void set_log_transform(bool v){log_transform = v;}
		bool doTransform(){return (sqrt_transform || boxcox_transform || log_transform);}
		bool get_sqrt_transform(){return sqrt_transform;}
		bool get_boxcox_transform(){return boxcox_transform;}
		bool get_log_transform(){return log_transform;}
		void performTransforms(DataSet*);
		void undoTransforms(DataSet*);
		vector<string> getToTransform(){return to_transform;}
		void setToTransform(vector<string> s){to_transform = s;}

		///get/set mult comparison fixed lambda
		void set_fixed_lambda(double v){fixed_lambda = v; do_fixed_lambda = true;}
		double get_fixed_lambda(){return fixed_lambda;}
		bool doFixedLambda(){return do_fixed_lambda;}

		///get/set xchr model var (model req plink)
		int get_xchr_model(){return xchr_model;} //if dominant or recessive, set to 0.  Default = 1
		void set_xchr_model(int v){xchr_model = v;};
		bool getSexEffect(){return sex_effect;}
		void setSexEffect(bool b){sex_effect = b;}
		bool getAutoSexEffect(){return auto_sex_effect;}
		void setAutoSexEffect(bool b){auto_sex_effect = b;}
		bool getQFAM_total(){return qfam_total;}
		void setQFAM_total(bool v){qfam_total = v;}
		bool getQFAM_between(){return qfam_between;}
		void setQFAM_between(bool v){qfam_between = v;}
		bool getQFAM_within1(){return qfam_within1;}
		void setQFAM_within1(bool v){qfam_within1 = v;}
		bool getQFAM_within2(){return qfam_within2;}
		void setQFAM_within2(bool v){qfam_within2 = v;}
		bool getQFAM_adaptive(){return qfam_adaptive;}
		void setQFAM_adaptive(bool v){qfam_adaptive = v;}
		double get_vif_threshold(){return vif_threshold;}
		void set_vif_threshold(double v){vif_threshold = v;}

		///get/set boolean to perform marker thresholding low
		bool doThreshMarkersLow(){return dothresh_markers_min;};
		void setDoThreshMarkersLow(bool b){dothresh_markers_min = b;};

		///get/set boolean to perform marker thresholding high
		bool doThreshMarkersHigh(){return dothresh_markers_max;};
		void setDoThreshMarkersHigh(bool b){dothresh_markers_max = b;};

		///get/set marker low threshold
		double getThreshMarkersLow(){return thresh_markers_min;};
		void setThreshMarkersLow(double d){thresh_markers_min = d; setDoThreshMarkersLow(true);};

		//get/set marker high threshold
		double getThreshMarkersHigh(){return thresh_markers_max;};
		void setThreshMarkersHigh(double d){thresh_markers_max = d; setDoThreshMarkersHigh(true);};

		//get/set boolean to perform sample threshold low
		bool doThreshSamplesLow(){return dothresh_samples_min;};
		void setDoThreshSamplesLow(bool b){dothresh_samples_min = b;};

		//get/set boolean to perform sample threshold high
		bool doThreshSamplesHigh(){return dothresh_samples_max;};
		void setDoThreshSamplesHigh(bool b){dothresh_samples_max = b;};

		//get/set sample thresh low
		double getThreshSamplesLow(){return thresh_samples_min;};
		void setThreshSamplesLow(double d){thresh_samples_min = d; setDoThreshSamplesLow(true);};

		//get/set sample thresh high
		double getThreshSamplesHigh(){return thresh_samples_max;};
		void setThreshSamplesHigh(double d){thresh_samples_max = d; setDoThreshSamplesHigh(true);};

		//get/set boolean to perform family threshold low
		bool doThreshFamiliesLow(){return dothresh_families_min;};
		void setDoThreshFamiliesLow(bool b){dothresh_families_min = b;};

		//get/set boolean to perform family threshold high
		bool doThreshFamiliesHigh(){return dothresh_families_max;};
		void setDoThreshFamiliesHigh(bool b){dothresh_families_max = b;};

		//get/set family thresh low
		double getThreshFamiliesLow(){return thresh_families_min;};
		void setThreshFamiliesLow(double d){thresh_families_min = d; setDoThreshFamiliesLow(true);};
		//get/set family thresh high
		double getThreshFamiliesHigh(){return thresh_families_max;};
		void setThreshFamiliesHigh(double d){thresh_families_max = d; setDoThreshFamiliesHigh(true);};

		//get/set chrom restriction
		int getChrom(){return chrom;};
		void setChrom(int c){chrom = c; setDoChrom(true);};
		bool doChrom(){return dochrom;};
		void setDoChrom(bool b){dochrom = b;};

		//get/set trio
		bool doTriosOnly(){return trio;};
		void setTriosOnly(bool b){trio = b;};

		//get/set zero genotypes
		bool zeroGenos(){return me_zero;};
		void setZeroGenos(bool b){me_zero = b;};

		//get/set zero L2 genos
		bool zeroL2Genos(){return me_zero_l2;};
		void setZeroL2Genos(bool b){me_zero_l2 = b;};

		//get/set zero L2 family genos
		bool zeroL2FamGenos(){return me_zero_l2_fams;};
		void setZeroL2FamGenos(bool d){me_zero_l2_fams = d;};

		//get/set basepair low
		bool doBpLow(){return dobp_min;};
		void setDoBpLow(bool b){dobp_min = b;};
		int getBpLow(){return bp_min;};
		void setBpLow(int b){bp_min = b; setDoBpLow(true);};

		//get/set basepair high
		bool doBpHigh(){return dobp_max;};
		void setDoBpHigh(bool b){dobp_max = b;};
		int getBpHigh(){return bp_max;};
		void setBpHigh(int b){bp_max = b; setDoBpHigh(true);};

		//get/set basepair spacing
		bool doBpSpace(){return dobp_space;};
		void setDoBpSpace(bool b){dobp_space = b;};
		int getBpSpace(){return bp_space;};
		void setBpSpace(int b){bp_space = b; setDoBpSpace(true);};


		//checks to see if parameter is equal to the chromosome set for this option
		bool checkChrom(int);
		//checks to see if parameter is within basepair specifications based on options set
		bool checkBp(int);

		//get/set disease name
		string getDisease(){return disease;};
		void setDisease(string n){disease = n;};

		//get/set stratification file name
		string getStratFile(){return struct_strat_file;};
		void setStratFile(string n){struct_strat_file = n;};
		bool haveStratification(){return stratification;};
		void setHaveStratification(bool b){stratification = b;};
		void readStratificationFile(string f);
		map<string, int> getStratificationMap(){return strat_map;};
		void setStratification(map<string, int> m){strat_map = m; setHaveStratification(true);};
		//searches center_map for center code
		int findStratification(string s);

		//get/set parents only
		bool doParentsOnly(){return parents_only;};
		void setDoParentsOnly(bool b){parents_only = b;};

		//get/set perform homozygous permutations
		bool doHomozygPermute(){return homozyg_permute;};
		void setDoHomozygPermute(bool b){homozyg_permute = b;};

		//get/set homozyg permutations count
		int getHomozygPermutationCount(){return homozyg_permute_count;};
		void setHomozygPermutationCount(int i){homozyg_permute_count = i;};

		bool doHomozygRaw(){return homozyg_raw;};
		void setDoHomozygRaw(bool b){homozyg_raw = b;};

		//get/set number of zeros to allow in homozygous
		int getHomozygZeros(){return homozyg_zeros;};
		void setHomozygZeros(int b){homozyg_zeros = b;};

		//get/set number of markers in homozyg span
		int getHomozygSpan(){return homozyg_span;};
		void setHomozygSpan(int b){homozyg_span = b;};

		//get/set minimum # of samples for valid span
		int getHomozygMinSamp(){return homozyg_min_samp;};
		void setHomozygMinSamp(int b){homozyg_min_samp = b;};

		//get/set perform WGHA version of homozygous span detection
		bool doHomozygWGHA(){return homozyg_wgha;};
		void setDoHomozygWGHA(bool b){homozyg_wgha = b;};

		//get/set perform ld calculations only
		bool doLDCalcOnly(){return ld_calc_only;};
		void setDoLDCalcOnly(bool b){ld_calc_only = b;};

		//get/set perform ld pairwise pruning
		bool doLDpairwise(){return ld_pairwise;};
		void setDoLDpairwise(bool b){ld_pairwise = b;};

		//get/set perform ld vif pruning
		bool doLDvif(){return ld_vif;};
		void setDoLDvif(bool b){ld_vif = b;};

		//get/set prune snps
		bool doLDchop(){return ld_chop;};
		void setLDchop(bool b){ld_chop = b;};

		//get/set LD window (snp count)
		int getLDWin(){return ld_window;};
		void setLDWin(int w){ld_window = w;};

		//get/set LD window (kb spacing)
		int getLDWinKB(){return ld_window_kb;};
		void setLDWinKB(int k){ld_window_kb = k;};

		//get/set LD threshold
		double getLDthreshold(){return ld_threshold;};
		void setLDthreshold(double d){ld_threshold = d;};

		//get/set LD pairwise threshold
		double getLDPWthreshold(){return ld_pw_threshold;};
		void setLDPWthreshold(double d){ld_pw_threshold = d;};

		//get/set LD step
		int getLDstep(){return ld_step;};
		void setLDstep(int s){ld_step = s;};

		//get/set filter based on overall results (overall column in output)
		bool doFilterOverall(){return filter_overall;};
		void setDoFilterOverall(bool b){filter_overall = b;};

		//get/set filter based on file values (allele-freq MAF)
		bool doFilterFile(){return filter_file;};
		void setDoFilterFile(bool b){filter_file = b;};

		//CMH
		//get/set cmh 2x2xk
		bool doCMH2x2xK(){return cmh_22k;}
		bool getCMH2x2xK(){return cmh_22k;}
		void setCMH2x2xK(bool v){cmh_22k = v;}
		//get/set cmh ixjxk
		bool doCMHIxJxK(){return cmh_ijk;}
		bool getCMHIxJxK(){return cmh_ijk;}
		void setCMHIxJxK(bool v){cmh_ijk = v;}
		//get/set cmh ordinal
		bool doCMHOrdinal(){return cmh_ordinal;}
		bool getCMHOrdinal(){return cmh_ordinal;}
		void setCMHOrdinal(bool v){cmh_ordinal = v;}
		//get/set breslowday
		bool doBreslowDay(){return breslowday;}
		bool getBreslowDay(){return breslowday;}
		void setBreslowDay(bool v){breslowday = v;}


		//get/set homozyg probability threshold for sequence
		bool doHomozygSeqTest(){return homozyg_seq_test;};
		void setDoHomozygSeqTest(bool b){homozyg_seq_test = b;};
		double getHomozygSeqVal(){return homozyg_seq_val;};
		void setHomozygSeqVal(double d){homozyg_seq_val = d; setDoHomozygSeqTest(true);};

		//get/set center code info
		bool haveCenterCodes(){return center_codes;};
		void setHaveCenterCodes(bool b){center_codes = b;};
		string getCenterFile(){return center_file;};
		void setCenterFile(string s){center_file = s; setHaveCenterCodes(true);};
		void readCenterFile(string f);
		map<string, string> getCenterCodes(){return center_map;};
		void setCenterCodes(map<string, string> m){center_map = m; setHaveCenterCodes(true);};
		//searches center_map for center code
		string findCenterCode(string s);

		//get/set penetrance info
		bool havePenetrance(){return penetrance_codes;};
		void setHavePenetrance(bool b){penetrance_codes = b;};
		string getPenetranceFile(){return penetrance_file;};
		void setPenetranceFile(string s){penetrance_file = s; setHavePenetrance(true);};
		void readPenetranceFile(string f);
		//searches penetrance_map for penetrance code
		float findPenetranceCode(string s);
		map<string, float> getPenetranceCodes(){return penetrance_map;};
		void setPenetranceCodes(map<string, float> m){penetrance_map = m; setHavePenetrance(true);};

		//use unaffected spouses
		bool doUnaffSpousesOnly(){return unaff_spouses_only;};
		void setDoUnaffSpousesOnly(bool b){unaff_spouses_only = b;};

		//use unknown spouses
		bool doUnknownSpouses(){return unk_spouses;};
		void setDoUnknownSpouses(bool b){unk_spouses = b;};

		//use founders only
		bool doFoundersOnly(){return founders_only;};
		void setDoFoundersOnly(bool b){founders_only = b;};
		void setFoundersOnly(){setDoFoundersOnly(true);};

		//use random child
		bool doRandomChild(){return random_child;};
		void setDoRandomChild(bool b){random_child = b;};

		//use all children
		bool doAllChildren(){return all_children;};
		void setDoAllChildren(bool b){all_children = b;};

		//use all samples
		bool doAll(){return all;};
		void setDoAll(bool b){all = b;};

		//deprecated
		bool doDeletion(){return deletion;};
		void setDoDeletion(bool b){deletion = b;};

		//get/set deletion span
		int getDeletionSpan(){return deletion_span;};
		void setDeletionSpan(int s){deletion_span = s;};

		//get/set output parental info
		bool doParental(){return parental;};
		void setDoParental(bool b){parental = b;};

		//get/set output gender info
		bool doGender(){return gender;};
		void setDoGender(bool b){gender = b;};

		//get/set output case control info
		bool doCaseControl(){return casecontrol;};
		void setDoCaseControl(bool b){casecontrol = b;};

		//get/set remove mono-allelic snps
		bool doRmMono(){return rm_mono;};
		void setDoRmMono(bool b){rm_mono = b;};

		//get/set remove heterozygous only snps
		bool doRmHetOnly(){return rm_het_only;};
		void setDoRmHetOnly(bool b){rm_het_only = b;};

		//get/set include excluded samples
		bool doIncExcludedSamples(){return inc_excluded_samples;};
		void setDoIncExcludedSamples(bool b){inc_excluded_samples = b;};

		//get/set include disabled samples
		bool doIncDisabledSamples(){return inc_disabled_samples;};
		void setDoIncDisabledSamples(bool b){inc_disabled_samples = b;};

		//get/set zero excluded sample genotypes
		bool doZeroExcluded(){return zero_excluded;};
		void setDoZeroExcluded(bool b){zero_excluded = b;};

		//get/set zero disabled sample genotypes
		bool doZeroDisabled(){return zero_disabled;};
		void setDoZeroDisabled(bool b){zero_disabled = b;};

		//get/set number of random markers to use
		bool doRandomMarkers(){return do_random_markers;};
		void setDoRandomMarkers(bool b){do_random_markers = b;};
		int getRandomMarkers(){return random_markers;};
		void setRandomMarkers(int s){random_markers = s; setDoRandomMarkers(true);};
		void setRandomRepeat(bool b){do_random_repeat = b;}
		bool doRandomRepeat(){return do_random_repeat;}

		//get/set number of sets of snps
		bool doSets(){return do_sets;};
		void setDoSets(bool b){do_sets = b;};
		int getSets(){return sets;};
		void setSets(int s){sets = s; setDoSets(true);};

		//get/set number of random samples to use & sets
		void setRandSamps(int i){rand_samps = i;}
		int getRandSamps(){return rand_samps;}
		void setRandSampsRepeat(bool b){rand_samps_repeat = b;}
		bool getRandSampsRepeat(){return rand_samps_repeat;}
		void setSetsSamps(int i){sets_samps = i;}
		int getSetsSamps(){return sets_samps;}
		void setPercentSamps(vector<float> s){percent_samps = s;}
		vector<float> getPercentSamps(){return percent_samps;}
		void setPercentCases(float i){percent_cases = i;}
		float getPercentCases(){return percent_cases;}
		void setPercentControls(float i){percent_controls = i;}
		float getPercentControls(){return percent_controls;}

		//get/set output file suffix
		string getOut(){return out;};
		void setOut(string s){out = s;};
		void setOverrideOut(string s){override_out = s;}
		string getOverrideOut(){return convertString(override_out);}

		//parses string of options
		string convertString(string);

		//covariate and trait methods
		string getCovarMissing(){return covar_missing;};
		void setCovarMissing(string s){covar_missing = s;};
		string getTraitMissing(){return trait_missing;};
		void setTraitMissing(string s){trait_missing = s;};
		double getDefaultCovarMissing(){return default_missing_value;};//no set method
		double getDefaultTraitMissing(){return default_missing_value;};//no set method
		bool doCovars(){return do_covs;};
		void setDoCovars(bool b){do_covs = b;};
		//use covars by name
		bool doCovarsName(){return do_covs_name;};
		void setDoCovarsName(bool b){do_covs_name = b;};
		//use covars by number
		bool doCovarsNumber(){return do_covs_number;};
		void setDoCovarsNumber(bool b){do_covs_number = b;};
		bool doGXE(){return do_gxe;};
		void setDoGXE(bool b){do_gxe = b;};
		bool doGXEName(){return do_gxe_name;};
		void setDoGXEName(bool b){do_gxe_name = b;};
		bool doGXENumber(){return do_gxe_number;};
		void setDoGXENumber(bool b){do_gxe_number = b;};		
		//set/get covar file
		bool doCovarsFile(){return do_covs_file;};
		void setDoCovarsFile(bool b){do_covs_file = b; setDoCovars(true);};
		//set run only covars
		void setRunCovarsOnly(){do_covars_only=true;}
		bool runCovarsOnly(){return do_covars_only;}
		//set/get trait file
		bool doTraitsFile(){return do_traits_file;};
		void setDoTraitsFile(bool b){do_traits_file = b; setDoTraits(true);};
		bool doTraits(){return do_traits;};
		void setDoTraits(bool b){do_traits = b;};
		//use traits by name
		bool doTraitsName(){return do_traits_name;};
		void setDoTraitsName(bool b){do_traits_name = b;};
		//use traits by number
		bool doTraitsNumber(){return do_traits_number;};
		void setDoTraitsNumber(bool b){do_traits_number = b;};
		//get/set covariates to use
		vector<string> getCovars(){return cov_use;};
		void setCovars(vector<string> s){cov_use = s;};
		//get/set GXE variables to use (from covariates file)
		vector<string> getGXEcovars(){return gxe_use;}
		void setGXEcovars(vector<string> s){gxe_use = s;}
		//get/set covariate map
		vector<string> getCovarMap(){return cov_map;};
		void setCovarMap(vector<string> s){cov_map = s;};
		//get/set covariate data points
		map<string, vector<double> > getCovarData(){return covs;};
		void setCovarData(map<string, vector<double> > m){covs = m;};
		//get/set traits to use
		vector<string> getTraits(){return trait_use;};
		void setTraits(vector<string> s){trait_use = s;};
		//get/set trait map
		vector<string> getTraitMap(){return trait_map;};
		void setTraitMap(vector<string> s){trait_map = s;};
		//get/set trait data points
		map<string, vector<double> > getTraitData(){return traits;};
		void setTraitData(map<string, vector<double> > m){traits = m;};
		//read covariates from file specified in options
		void readCovariates(vector<Sample*>*);
		//read traits from file specified in options
		void readTraits(vector<Sample*>*);
		void setUsePheno(bool b){alternate_pheno = b;}
		bool getUsePheno(){return alternate_pheno;}
		int getPhenoLoc(){return pheno_loc;}
		void setPhenoLoc(int l){pheno_loc = l;}
		void setPhenoName(string s){pheno_name = s;}
		string getPhenoName(){return pheno_name;}
		vector<string> getUseLoci(){return loci_use;}
		void setPhenoMissing(int v){pheno_missing = v;}
		int getPhenoMissing(){return pheno_missing;}

		//get/set grouping information
		bool doGroupFile(){return do_group_file;};
		void setDoGroupFile(bool b){do_group_file = b;};
		string getGroupFile(){return group_file;};
		void setGroupFile(string s){group_file = s; setDoGroupFile(true);};
		map<string, vector<Sample*> > getGroups(){return groups;};
		map<Sample*, string> getSampleGroups(){return sample_groups;};
		void setGroups(map<string, vector<Sample*> > m){groups = m; setDoGroupFile(true);};
		void setSampleGroups(map<Sample*, string> m){sample_groups = m; setDoGroupFile(true);};
		map<string, vector<Family*> > getGroupFamilies(){return group_families;};
		void setGroupFamilies(map<string, vector<Family*> > m){group_families = m; setDoGroupFile(true);};
		//reads group file
		void readGroups(vector<Sample*>*);
		//Cluster info
		bool doClusterFile(){return do_cluster_file;};
		void setDoClusterFile(bool b){do_cluster_file = b;};
		string getClusterFile(){return cluster_file;};
		void setClusterFile(string f){cluster_file = f; setDoClusterFile(true);};
		map<string, vector<Sample*> > getClusters(){return clusters;};
		map<Sample*, int> getSampleClusters(){return sample_clusters;};
		void setSampleClustersString(map<string, string> v){sample_clusters_string = v;}
		map<string, string> getSampleClustersString(){return sample_clusters_string;}
		void readClusters(vector<Sample*>*);
		void readClustersFromString(vector<Sample*>*);

		//create a vector of Families based on sample vector
		vector<Family*> generateFamilySet(vector<Sample*>* samps);

		//get/set Ped file and map file
		bool doPedFile(){return do_pedfile;};
		void setDoPedFile(bool b){do_pedfile = b;};
		bool doMapFile(){return do_mapfile;};
		void setDoMapFile(bool b){do_mapfile = b;};
		bool doMDRFIle(){return do_mdrfile;};
		void setDoMDRFile(bool b){do_mdrfile = b;};
		bool doMdrPedFIle(){return do_mdrpedfile;};
		void setDoMdrPedFile(bool b){do_mdrpedfile = b;};
		bool doMdrMapFIle(){return do_mdrmapfile;};
		void setDoMdrMapFile(bool b){do_mdrmapfile = b;};
		bool getMDRGuiOutput(){return mdr_gui_output;};
		void setMDRGuiOutput(bool b){mdr_gui_output = b;};

		string getPedFile(){return ped_file;};
		void setPedFile(string s){ped_file = s; setDoPedFile(true);};
		string getMapFile(){return map_file;};
		void setMapFile(string s){map_file = s; setDoMapFile(true);};
		string getMDRFile(){return mdr_file;};
		void setMDRFile(string s){mdr_file = s; setDoMDRFile(true);};
		string getMdrPedFile(){return mdr_ped_file;};
		void setMdrPedFile(string s){mdr_ped_file = s; setDoMdrPedFile(true);};
		string getMdrMapFile(){return mdr_map_file;};
		void setMdrMapFile(string s){mdr_map_file = s; setDoMdrMapFile(true);};

		//get/set Transposed ped file and family map file
		bool doTPedFile(){return do_tpedfile;};
		void setDoTPedFile(bool b){do_tpedfile = b;};
		bool doTFamFile(){return do_tfamfile;};
		void setDoTFamFile(bool b){do_tfamfile = b;};
		string getTPedFile(){return tped_file;};
		void setTPedFile(string s){tped_file = s; setDoTPedFile(true);};
		string getTFamFile(){return tfam_file;};
		void setTFamFile(string s){tfam_file = s; setDoTFamFile(true);};

		//get/set binary input file prefix
		bool doBinInput(){return do_binfile;};
		void setDoBinInput(bool b){do_binfile = b;};
		string getBinInput(){return bin_prefix;};
		void setBinInput(string b){bin_prefix = b; setDoBinInput(true);};

		//get/set confidence interval
		double getCI(){return confidence_interval;};
		void setCI(double d){confidence_interval = d;};

		//get/set perform ACGT -> 1234 conversion
		bool doAllele1234(){return allele1234;};
		void setDoAllele1234(bool b){allele1234 = b;};

		//get/set perform 1234 -> ACGT conversion
		bool doAlleleACGT(){return alleleACGT;};
		void setDoAlleleACGT(bool b){alleleACGT = b;};

		//get/set perform ACGT, 1234 -> 1 2 mapping (not implemented yet)
		bool doAllele12(){return allele12;};
		void setDoAllele12(bool b){allele12 = b;};

		//get/set perform custom allele mapping
		bool doAlleleCustom(){return allelecustom;};
		void setDoAlleleCustom(bool b){allelecustom = b;};

		//converts string allele to mapping as specified by options
		string getAlleleMapping(string s);

		//get/set perform HWEPT hardy-weinberg
		bool doHWEPT(){return hwept;};
		void setDoHWEPT(bool b){hwept = b;};

		//get/set permutations
		int getPermutations(){return perms;};
		void setPermutations(int b){perms = b;};

		//get/set prevalance
		double getPrevalance(){return prevalence;};
		void setPrevalance(double d){prevalence = d;};

		//get/set Linear Regression options
		string getLinRModelType(){return linr_modType;}
		void setLinRModelType(string modelType){linr_modType = modelType;}
		bool getLinRInteraction(){return linr_interaction;}
		void setLinRInteraction(bool b){linr_interaction = b;}
		bool getLinRNoMainSnp(){return linr_no_main_snp;}
		void setLinRNoMainSnp(bool b){linr_no_main_snp = b;}
		bool getLinRCondition(){return linr_condition;}
		void setLinRCondition(bool v){linr_condition = v;}
		void setLinRConditionString(string v){linr_condition_string = v;}
		string getLinRConditionString(){return linr_condition_string;}
		void setLinRConditionFile(string v){linr_condition_file = v;}
		string getLinRConditionFile(){return linr_condition_file;}
		void addLinRConditionList(int l){linr_condition_list.push_back(l); setLinRCondition(true);}
		vector<int> getLinRConditionList(){return linr_condition_list;}
		void parseLinRConditionList(vector<Marker*>*);
		void readLinRConditionFile(vector<Marker*>*);


		//get/set Logistic Regression options
		void setLRFullInteraction(bool interact){lr_fullInteraction = interact;}
		void setLRIncludeInteractions(bool interact){lr_includeInteractions = interact;}
		void setLRModelType(string modelType){lr_modType = modelType;}
		void setLRMaximumIterations(unsigned int maxIter){lr_maxIterations = maxIter;}
		unsigned int getLRMaximumIterations(){return lr_maxIterations;}
		bool getLRFullInteraction(){return lr_fullInteraction;}
		bool getLRIncludeInteractions(){return lr_includeInteractions;}
		string getLRModelType(){return lr_modType;}

		// get/set Conditional Logistic Regression options
		void setDoCondLR(bool b){cond_lr = b;}
		bool doCondLR(){return cond_lr;}
		void setCondLRFullInteraction(bool interact){cond_lr_fullInteraction = interact;}
		void setCondLRIncludeInteractions(bool interact){cond_lr_includeInteractions = interact;}
		void setCondLRModelType(string modelType){cond_lr_modType = modelType;}
		void setCondLRMaximumIterations(unsigned int maxIter){cond_lr_maxIterations = maxIter;}
		unsigned int getCondLRMaximumIterations(){return cond_lr_maxIterations;}
		bool getCondLRFullInteraction(){return cond_lr_fullInteraction;}
		bool getCondLRIncludeInteractions(){return cond_lr_includeInteractions;}
		string getCondLRModelType(){return cond_lr_modType;}

		// get/set MDR options
		void setOnlySetThreshold(bool option){mdr_set_only_threshold = option;}
		bool getOnlySetThreshold(){return mdr_set_only_threshold;}

		// UncertaintyCoefficient options
		void setUncertaintyCoeffTotalType(string tot_type){uncert_coeff_total_type = tot_type;}
		string getUncertaintyCoeffTotalType(){return uncert_coeff_total_type;}

		// LikelihoodRatio options
		void setLikelihoodRatioTotalType(string tot_type){llr_total_type = tot_type;}
		string getLikelihoodRatioTotalType(){return llr_total_type;}

		//get/set remove missing parents
		bool doRemMissingParents(){return rem_missing_parents;};
		void setRemMissingParents(bool b){rem_missing_parents = b;};
		bool getRemMissingParents(){return rem_missing_parents;};

		//get/set dummy missing parents
		bool doDummyMissingParents(){return dummy_missing_parents;};
		void setDummyMissingParents(bool b){dummy_missing_parents = b;};
		bool getDummyMissingParents(){return dummy_missing_parents;};

		//get/set zero incomplete trio ids
		bool doZeroIncompleteTrioIds(){return zero_incomplete_trio_ids;};
		void setZeroIncompleteTrioIds(bool b){zero_incomplete_trio_ids = b;};
		bool getZeroIncompleteTrioIds(){return zero_incomplete_trio_ids;};

		//get/set dummy incomplete parent ids
		bool doDummyIncompleteParentIds(){return dummy_incomplete_parent_ids;};
		void setDummyIncompleteParentIds(bool b){dummy_incomplete_parent_ids = b;};
		bool getDummyIncompleteParentIds(){return dummy_incomplete_parent_ids;};

		// OddsRatio options
		void setOddsRatioTotalType(string tot_type){oddsratio_total_type = tot_type;}
		string getOddsRatioTotalType(){return oddsratio_total_type;}

		// NMI options
		void setNMITotalType(string tot_type){nmi_total_type = tot_type;}
		string getNMITotalType(){return nmi_total_type;}
		void setNMITransposed(bool trans){nmi_transposed = trans;}
		bool getNMITransposed(){return nmi_transposed;}

		// Referent allele in map file
		bool getMapContainsReferent(){return map_includes_ref_allele;}
		void setMapContainsReferent(bool b){map_includes_ref_allele = b;}

		// MDRPDT options
		void setMDRPDTNumPTests(int p){mdrpdt_ptests = p;}
		int getMDRPDTNumPTests(){return mdrpdt_ptests;}
		void setMDRPDTRandSeed(int rs){mdrpdt_seed = rs;}
		int getMDRPDTRandSeed(){return mdrpdt_seed;}
		void setMDRPDTNumCrossVals(int xv){mdrpdt_xv = xv;}
		int getMDRDPTNumCrossvals(){return mdrpdt_xv;}
		void setMDRPDTMaxCombo(int max){mdrpdt_maxcombo=max;}
		int getMDRPDTMaxCombo(){return mdrpdt_maxcombo;}
		void setMDRPDTMinCombo(int min){mdrpdt_mincombo=min;}
		int getMDRPDTMinCombo(){return mdrpdt_mincombo;}

		/// get/set multiple comparisons
		bool getMultCompare(){return mult_compare;}
		bool doMultCompare(){return mult_compare;}
		void setMultCompare(bool v){mult_compare = v;}

		/// get/set ibs options
		map<string, vector<string> > get_ibs_pairs(){return ibs_pairs;}
		map<string, vector<string> > get_ibs_trio_pairs(){return ibs_trio_pairs;}
		bool getDoIBSPairs(){return do_ibs_pairs;}
		bool getDoIBSTrioPairs(){return do_ibs_trio_pairs;}
		bool getDoIBSAllPairs(){return do_ibs_all_pairs;}
		bool getDoIBSAllTrioPairs(){return do_ibs_all_trio_pairs;}
		bool getDoIBSTrioTransmissions(){return do_ibs_trio_transmissions;}
		void setDoIBSTrioTransmissions(bool b){do_ibs_trio_transmissions = b;}
		void setDoIBSPairs(bool b){do_ibs_pairs = b;}
		void setDoIBSTrioPairs(bool b){do_ibs_trio_pairs = b;}
		void setDoIBSAllPairs(bool b){do_ibs_all_pairs = b;}
		void setDoIBSAllTrioPairs(bool b){do_ibs_all_trio_pairs = b;}
		void setIBSPairsFile(string s){ibs_file = s;}
		void setIBSTrioPairsFile(string s){ibs_trios_file = s;}
		void setIbsTriosRaw(bool b){do_ibs_trios_raw = b;}
		bool getIbsTriosRaw(){return do_ibs_trios_raw;}
		string getIBSPairsFile(){return ibs_file;}
		string getIBSTrioPairsFile(){return ibs_trios_file;}
		void set_ibs_pairs(map<string, vector<string> > v){ibs_pairs = v;}
		void set_ibs_trio_pairs(map<string, vector<string> > v){ibs_trio_pairs = v;}
		void readIBSPairsFile(string);
		void readIBSTrioPairsFile(string);

		/// cluster missing options
		int getMaxClusterSize(){return max_cluster_size;}
		int getMaxClusterN(){return max_cluster_N;}
		bool getClusterOnPheno(){return cluster_on_pheno;}
		bool getClusterSelcon(){return cluster_selcon;}
		bool getClusterOnMcc(){return cluster_on_mcc;}

		/// frequency file info
		map<string, float> getFrequencies(){return frequencies;}
		void setFrequencies(map<string, float> v){frequencies = v;}

		/// sample/bprange filtering
		void readSampleBprangeFile();
		void setSampleBprangeFile(string s){sampbprange_file = s;}
		vector<Sample*> getSampleBprangeSamples(){return sampbprange_samples;}
		vector<vector<Marker*> > getSampleBprangeMarkers(){return sampbprange_markers;}
	};
};
#endif
