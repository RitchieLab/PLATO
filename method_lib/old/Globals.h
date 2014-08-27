#ifndef GLOBALS_H
#define GLOBALS_H

namespace Methods{
/***IF USING MPI, UNCOMMENT BELOW***/
//#define _USEMPI_

/*#define SR_NUMPERNODE 10
#define SR_MARKERDATA 20
#define SR_ALLELEDATA 30
//#define SR_FAMILYDATA 40
static const int SR_FAMILYDATA = 40;
#define SR_INDDATA 50
#define SR_CHREFFDATA 60
#define DONE 999
#define RECEIVE 1
*/
#define NUMCHROMS 23
#define ENZ_FME 10
#define ENZ_FZERO 20
#define ENZ_FTOT 30
#define ENZ_IME 40
#define ENZ_IZERO 50
#define ENZ_ITOT 60
#define ENZ_IGHET 70
#define ENZ_IGTOT 80

	typedef struct{
		int famid;
		int ind;
		int gender_errors;
		int zero_calls;
		int total_calls;
		float geno_eff;
		int enabled;
		char sex[2];
		int mother;
		int father;
		int child;
		int mymother;
		int myfather;
		int mendelian_errors;
		int casecontrol;
		int is_case;
		int is_control;
		int affected;
		char locale[4];
		char plate[21];
		char well[6];
	} ind_data;

		typedef struct {
			int famid;
			int ind;
			int zero;
			int total;
			char chrom[3];
		} chr_eff_data;

		typedef struct {
			int famid;
			int ind;
			int type;
			int count;
			char enzyme[10];
		} enzyme_count;
		
		typedef struct {
			int famid;
			int mendelian_errors;
		} fam_data;

		typedef struct {
			int sysprobe;
			int mendelian_errors;
			int gender_errors;
			int zero_calls;
			int zero_calls_case;
			int zero_calls_control;
			int total_calls;
			int total_calls_case;
			int total_calls_control;
			float percent_eff;
			float percent_eff_case;
			float percent_eff_control;
			double pval;
			double chi;
			int bploc;
			float minor_allele_freq_tdt;
			float quality_score;
			int qual_score_total;
			int ge_male_genos;
			int qs_total_het[50];
			int qs_total_maj[50];
			int qs_total_min[50];
			int qs_delta_fm[50];
			int qs_delta_fc[50];
			int qs_delta_mc[50];
			char enzyme[4];
			char dbsnp_rsid[30];
			char chrom[3];
			char probe[30];
			int tdt_fams_used;
		} marker_data;

		typedef struct {
			int sysprobe;
			char allele1[2];
			char allele2[2];
			float parent_major_freq;
			float parent_minor_freq;        
			float father_major_freq;
			float mother_major_freq;
			float father_minor_freq;
			float mother_minor_freq;
			float child_major_freq;
			float child_minor_freq;
			float child_major_freq_f;
			float child_minor_freq_f;
			float child_major_freq_m;
			float child_minor_freq_m;
			//casecontrol
			float case_major_freq_f;
			float case_minor_freq_f;
			float case_major_freq_m;
			float case_minor_freq_m;
			float control_major_freq_f;
			float control_minor_freq_f;
			float control_major_freq_m;
			float control_minor_freq_m;

			int parent_major_count;
			int parent_minor_count;
			int father_major_count;
			int mother_major_count;
			int father_minor_count;
			int mother_minor_count;
			int child_major_count;
			int child_minor_count;
			int child_major_count_f;
			int child_minor_count_f;
			int child_major_count_m;
			int child_minor_count_m;
			//casecontrol
			int case_major_count_f;
			int case_minor_count_f;
			int case_major_count_m;
			int case_minor_count_m;
			int control_major_count_f;
			int control_minor_count_f;
			int control_major_count_m;
			int control_minor_count_m;

			int parent_het_count;
			int father_het_count;
			int mother_het_count;
			int child_het_count;
			int child_het_count_f;
			int child_het_count_m;
			//casecontrol
			int case_het_count_f;
			int case_het_count_m;
			int control_het_count_f;
			int control_het_count_m;

			int parent_homo1_count;
			int parent_homo2_count;
			int father_homo1_count;
			int mother_homo1_count;
			int father_homo2_count;
			int mother_homo2_count;
			int child_homo1_count;
			int child_homo2_count;
			int child_homo1_count_f;
			int child_homo2_count_f;
			int child_homo1_count_m;
			int child_homo2_count_m;
			//casecontrol
			int case_homo1_count_f;
			int case_homo2_count_f;
			int case_homo1_count_m;
			int case_homo2_count_m;
			int control_homo1_count_f;
			int control_homo2_count_f;
			int control_homo1_count_m;
			int control_homo2_count_m;

			float parent_hw;
			float father_hw;
			float mother_hw;
			float child_hw;
			float child_hw_f;
			float child_hw_m;
			//casecontrol
			float case_hw_f;
			float case_hw_m;
			float control_hw_f;
			float control_hw_m;

			int parents_used;
			int fathers_used;
			int mothers_used;
			int children_used;
			int children_used_f;
			int children_used_m;
			//casecontrol
			int cases_used_f;
			int cases_used_m;
			int controls_used_f;
			int controls_used_m;

		} allele_struct_data;
};												
#endif
