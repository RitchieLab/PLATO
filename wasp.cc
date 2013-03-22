/***************************************************************
*		wasp - Genome Association Study Pipeline
* Written by: Justin Giles
*	      Vanderbilt University
*	      Center for Human Genetics Research
*
*wasp.cc - Main class for wasp application
*
* See README file for complete overview.
***************************************************************/
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <string.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <bitset>
#include <ctime>
#include <Globals.h>
#include <vector>
#include <General.h>//"General.h"
#include "wasp.h"


using namespace std;
using namespace Methods;
#ifdef PLATOLIB
using namespace PlatoLib;
#endif

static string _WASPVER_ = "1.2.0";
//enum for batch file step switch
enum StepValue{
					   e_examplemodule,    //examplemodule
					   e_marker_geno_eff,  //marker-geno-eff
					   e_sample_geno_eff,  //sample-geno-eff
					   e_family_geno_eff,  //family-geno-eff
					   e_allele_freq,      //allele-freq
					   e_mendelian_error,  //mendelian-error
					   e_hw,               //hw
					   e_gender_error,     //gender-error
					   e_tdt,              //tdt
					   e_grr_output,       //grr-output
					   e_ped_output,       //ped-output
					   e_tped_output,		//tped-output
					   e_chisquare, //ca-con-chisquare
					   e_structure_output, //structure-output
					   e_phase_output,	   //phase-output
					   e_beagle_output,	   //beagle-output
					   e_lapis_output,	   //lapis-output
					   e_superlink_output, //superlink-output
					   e_homozygous, 	   //homozygous
					   e_ld, 	           //ld
					   e_powermarker_output,//powermarker-output
					   e_deletions,	       //deletions
					   e_fbat_output,	   //fbat-output
					   e_qtdt_output,	   //qtdt-output
					   e_mdr_output, 	   //mdr-output
					   e_pdt2_output,	   //output-pdt2
					   e_eigenstrat_output, //output-eigenstrat
					   e_concordance,	   //concordance
					   e_logreg,			//logreg
					   e_linreg,			//linear-reg
					   e_mitocheck,         //mito-check
					   e_ibs,				//IBS
					   e_cmh,				//CMH
					   e_mdr,				//MDR
					   e_mdrpdt,			//MDRPDT
					   e_cluster_missing,	//cluster-missing
					   e_mars,				//mars
					   e_filter_process,		//plato filter-process,
					   e_fst,				//fst
					   e_kinship,			//kinship
					   e_bin_output,			//binary output
					   e_epistasis,			//epistasis
					   e_interaction, // interaction
					   e_impute_output		//IMPUTE
					 };
enum cmdArgs{
	a_h,
	a_S,
	a_ped,
	a_map,
	//new 12-06-2010
	a_lgen_file,
	a_compound_genotypes,
	a_reference,
	a_bin_input,
	a_incmarkers,
	a_excmarkers,
	a_incsamples,
	a_excsamples,
	a_incfamilies,
	a_excfamilies,
	a_chrom,
	a_bp_min,
	a_bp_max,
	a_bp_span,
	a_pedinfo,
	a_freq_file,
	a_dog,
	a_noweb,
	a_out,
	a_no_summary,
	a_micro_sats,
	a_sampdesc,
	a_mapdesc,
	a_tfam,
	a_tped,
	a_mdrped,
	a_mdrmap,
	a_flip,
//	a_make_bin,
	a_print_families,
	a_zero_geno_file,
	a_remaining,
	a_inccenters,
	a_missing_geno,
	a_covar_missing,
	a_trait_missing,
	a_bp_space,
	a_todigit,
	a_exccovs,
	a_inccovs,
	a_exctraits,
	a_inctraits,
	a_exccovs_name,
	a_inccovs_name,
	a_exctraits_name,
	a_inctraits_name,
	a_covar_file,
	a_trait_file,
	a_pheno_missing,
	a_map_includes_ref,
	a_samplebprangefilter,
	a_threads,
	a_numthreads,
	a_autoonly,
	a_update_ids

};
static map<string, StepValue> s_mapStepValues;
static map<string, cmdArgs> s_mapcmdArgs;



/*
 *Function: Initialize()
 *Return: void
 *Description:
 *Initializes the batch file options enumeration.
 */
void Initialize(){
	srand((unsigned) time(0));

	s_mapStepValues["examplemodule"] = e_examplemodule;
	s_mapStepValues["marker-geno-eff"] = e_marker_geno_eff;
	s_mapStepValues["sample-geno-eff"] = e_sample_geno_eff;
	s_mapStepValues["family-geno-eff"] = e_family_geno_eff;
	s_mapStepValues["allele-freq"] = e_allele_freq;
	s_mapStepValues["mendelian-error"] = e_mendelian_error;
	s_mapStepValues["hw"] = e_hw;
	s_mapStepValues["gender-error"] = e_gender_error;
	s_mapStepValues["tdt"] = e_tdt;
	s_mapStepValues["output-grr"] = e_grr_output;
	s_mapStepValues["output-ped"] = e_ped_output;
	s_mapStepValues["output-eigenstrat"] = e_eigenstrat_output;
	s_mapStepValues["ibs"] = e_ibs;
	s_mapStepValues["output-tped"] = e_tped_output;
	s_mapStepValues["chisquare"] = e_chisquare;
	s_mapStepValues["output-structure"] = e_structure_output;
	s_mapStepValues["output-phase"] = e_phase_output;
	s_mapStepValues["output-beagle"] = e_beagle_output;
	s_mapStepValues["output-superlink"] = e_superlink_output;
	s_mapStepValues["homozygous"] = e_homozygous;
	s_mapStepValues["ld"] = e_ld;
	s_mapStepValues["output-powermarker"] = e_powermarker_output;
	s_mapStepValues["deletions"] = e_deletions;
	s_mapStepValues["concordance"] = e_concordance;
	s_mapStepValues["mito-check"] = e_mitocheck;
	s_mapStepValues["output-lapis"] = e_lapis_output;
	s_mapStepValues["output-fbat"] = e_fbat_output;
	s_mapStepValues["output-mdr"] = e_mdr_output;
	s_mapStepValues["output-qtdt"] = e_qtdt_output;
	s_mapStepValues["output-pdt2"] = e_pdt2_output;
	s_mapStepValues["logreg"] = e_logreg;
	s_mapStepValues["linear-reg"] = e_linreg;
	s_mapStepValues["cmh"] = e_cmh;
	s_mapStepValues["mdr"] = e_mdr;
	s_mapStepValues["mdrpdt"] = e_mdrpdt;
#ifdef HAVE_R
	s_mapStepValues["mars"] = e_mars;
#endif
	s_mapStepValues["cluster-missing"] = e_cluster_missing;
	s_mapStepValues["filter-process"] = e_filter_process;
	s_mapStepValues["fst"] = e_fst;
	s_mapStepValues["kinship"] = e_kinship;
	s_mapStepValues["output-bin"] = e_bin_output;
	s_mapStepValues["epistasis"] = e_epistasis;
	s_mapStepValues["interaction"] = e_interaction;
	s_mapStepValues["output-impute"] = e_impute_output;

	s_mapcmdArgs["-h"] = a_h;
	s_mapcmdArgs["-S"] = a_S;
	s_mapcmdArgs["-bin-input"] = a_bin_input;
	s_mapcmdArgs["-bp-max"] = a_bp_max;
	s_mapcmdArgs["-bp-min"] = a_bp_min;
	s_mapcmdArgs["-bp-span"] = a_bp_span;
	s_mapcmdArgs["-chrom"] = a_chrom;
	s_mapcmdArgs["-excfamilies"] = a_excfamilies;
	s_mapcmdArgs["-excmarkers"] = a_excmarkers;
	s_mapcmdArgs["-excsamples"] = a_excsamples;
	s_mapcmdArgs["-exccovs"] = a_exccovs;
	s_mapcmdArgs["-exccovs-name"] = a_exccovs_name;
	s_mapcmdArgs["-inccovs-name"] = a_inccovs_name;
	s_mapcmdArgs["-inccovs"] = a_inccovs;
	s_mapcmdArgs["-exctraits-name"] = a_exctraits_name;
	s_mapcmdArgs["-inctraits-name"] = a_inctraits_name;
	s_mapcmdArgs["-exctraits"] = a_exctraits;
	s_mapcmdArgs["-inctraits"] = a_inctraits;
	s_mapcmdArgs["-freq-file"] = a_freq_file;
	s_mapcmdArgs["-incfamilies"] = a_incfamilies;
	s_mapcmdArgs["-incmarkers"] = a_incmarkers;
	s_mapcmdArgs["-incsamples"] = a_incsamples;
	s_mapcmdArgs["-map"] = a_map;
	s_mapcmdArgs["-ped"] = a_ped;
	//new 12-06-2010
	s_mapcmdArgs["-lgen-file"] = a_lgen_file;
	s_mapcmdArgs["-compound-genotypes"] = a_compound_genotypes;
	s_mapcmdArgs["-reference"] = a_reference;
	s_mapcmdArgs["-pedinfo"] = a_pedinfo;
	s_mapcmdArgs["-dog"] = a_dog;
	s_mapcmdArgs["-flip"] = a_flip;
//	s_mapcmdArgs["-make-bin"] = a_make_bin;
	s_mapcmdArgs["-mapdesc"] = a_mapdesc;
	s_mapcmdArgs["-micro-sats"] = a_micro_sats;
	s_mapcmdArgs["-no-summary"] = a_no_summary;
	s_mapcmdArgs["-noweb"] = a_noweb;
	s_mapcmdArgs["-out"] = a_out;
	s_mapcmdArgs["-sampdesc"] = a_sampdesc;
	s_mapcmdArgs["-tped"] = a_tped;
	s_mapcmdArgs["-tfam"] = a_tfam;
	s_mapcmdArgs["-mdrped"] = a_mdrped;
	s_mapcmdArgs["-mdrmap"] = a_mdrmap;
	s_mapcmdArgs["-print-families"] = a_print_families;
	s_mapcmdArgs["-zero-geno-file"] = a_zero_geno_file;
	s_mapcmdArgs["-remaining"] = a_remaining;
	s_mapcmdArgs["-inccenters"] = a_inccenters;
	s_mapcmdArgs["-missing-geno"] = a_missing_geno;
	s_mapcmdArgs["-covar-missing"] = a_covar_missing;
	s_mapcmdArgs["-trait-missing"] = a_trait_missing;
	s_mapcmdArgs["-covar-file"] = a_covar_file;
	s_mapcmdArgs["-trait-file"] = a_trait_file;
	s_mapcmdArgs["-pheno-missing"] = a_pheno_missing;
	s_mapcmdArgs["-bp-space"] = a_bp_space;
	s_mapcmdArgs["-todigit"] = a_todigit;
	s_mapcmdArgs["-map-includes-ref"] = a_map_includes_ref;
	s_mapcmdArgs["-sample-bprange-filter"] = a_samplebprangefilter;
	s_mapcmdArgs["-threads"] = a_threads;
	s_mapcmdArgs["-numthreads"] = a_numthreads;
	s_mapcmdArgs["-auto-only"] = a_autoonly;
	//Added 02-23-2011 to support new feature -update-ids
	s_mapcmdArgs["-update-ids"] = a_update_ids;
}

int
main (int argc, char* argv[])
{
	try{
	Initialize();
	int myrank;
	myrank = 0;
	STEPS steps;
	ORDER proc_order;
	vector<string> usechroms;
	vector<long> usembrange;
	InputFilter cmd_filters;


	//Parse command line arguments
	//assume first argument is -h, -S or the batch file
	if(argc >= 2){
	   	string arg = argv[1];
	   	if(arg == "-h"){
			print_help();
			exit(1);
		}
		if(arg == "-S"){
			steps = initializeSteps();
			print_steps(steps);
			exit(1);
		}
		if(arg == "-O"){
			printOptions();
			exit(1);
		}
		if(arg == "-SO"){
			StepOptions printopts;
			printopts.printOptions();
			exit(1);
		}
		vector<string> arguments;

		for(int i = 2; i < argc; i++){
			string list = "";
			list.assign(argv[i]);
			arguments.push_back(list);
		}
		vector<string>::iterator niter = find(arguments.begin(), arguments.end(), "-noweb");
		map<string, vector<string> > batchargs = getBatchArgs(arg);
		if(niter == arguments.end()){
			webcheck(arguments, batchargs);
		}

		for(unsigned int j = 0; j < arguments.size(); j++){
			switch(s_mapcmdArgs[arguments[j]]){
				case a_map_includes_ref:
				{
					opts::_MAP_INCLUDES_REF_ = true;
					break;
			    }
				case a_ped:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-ped option requires specifying a filename. (Ex: -ped pedinput.ped).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-ped option requires specifying a filename. (Ex: -ped pedinput.ped).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-ped option missing file name?  Halting." << endl;
						exit(1);
					}

					opts::_PEDFILE_ = arguments[++j];
					break;
				}
				//new 12-06-2010
				case a_lgen_file:
				{
					if(j + 1 >= arguments.size())
					{
						cerr << "-lgen_file option requires specifying a filename. (Ex: -lgen-file lgeninput.lgen). Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() ==0)
					{
						cerr << "-lgen_file option requires specifying a filename. (Ex: -lgen-file lgeninput.lgen). Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-')
					{
						cerr << "-lgen option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_LGENFILE_ = arguments[++j];
					break;
				}
				//new 12-09-2010
				case a_compound_genotypes:
				{
					opts::_COMPOUND_GENOTYPES_ = true;
					break;
				}
				//new 12-10-2010
				case a_reference:
				{
					if(j + 1 >= arguments.size())
					{
						cerr << "-reference option requires specifying a filename.  (Ex: -reference references.txt). Halting!" << endl;
						exit(1);
					}
					string test = arguments[j + 1];
					if(test.length() == 0)
					{
						cerr << "-reference option requires specifying a filename.  (Ex: -reference references.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-')
					{
						cerr << "-reference option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_REFERENCE_FILE_ = arguments[++j];
					break;
				}
				case a_noweb:
					break;
				case a_autoonly:
					opts::_AUTOONLY_ = true;
					break;
				case a_missing_geno:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-missing-geno option requires specifying a missing code. (Ex: -missing-geno N). Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-missing-geno option requires specifying a missing code. (Ex: -missing-geno N). Halting!" << endl;
						exit(1);
					}
					//if(test[0] == '-'){
					//	cerr << "-missing-geno option missing code?  Halting." << endl;
					//	exit(1);
					//}
					opts::_NOCALL_ = arguments[++j];
					opts::_NOCALLMDR_ = opts::_NOCALL_;
					break;
				}
				case a_pheno_missing:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-pheno-missing option requires specifying a missing code. (Ex: -pheno-missing -9). Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-pheno-missing option requires specifying a missing code. (Ex: -pheno-missing -9). Halting!" << endl;
						exit(1);
					}
					opts::_PHENO_MISS_ = atoi(test.c_str());
					j++;
					break;
				}
				case a_covar_missing:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-covar-missing option requires specifying a missing code. (Ex: -covar-missing -9). Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-covar-missing option requires specifying a missing code. (Ex: -covar-missing -9). Halting!" << endl;
						exit(1);
					}
					opts::_COVAR_MISSING_ = test;
					j++;
					break;
				}
				case a_trait_missing:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-trait-missing option requires specifying a missing code. (Ex: -trait-missing -9). Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-trait-missing option requires specifying a missing code. (Ex: -trait-missing -9). Halting!" << endl;
						exit(1);
					}
					opts::_TRAIT_MISSING_ = test;
					j++;
					break;
				}
				case a_pedinfo:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-pedinfo option requires specifying a filename. (Ex: -pedinfo pedinfo.txt). Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-pedinfo option requires specifying a filename. (Ex: -pedinfo pedinfo.txt). Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-pedinfo option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_PEDINFO_ = arguments[++j];
					break;
				}
				case a_samplebprangefilter:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-sample-bprange-filter option requires specifying a filename. (Ex: -sample-bprange-filter filter.txt). Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-sample-bprange-filter option requires specifying a filename. (Ex: -sample-bprange-filter filter.txt). Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-sample-bprange-filter option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_SAMPLEBPRANGEFILTER_ = arguments[++j];
					break;
				}
				case a_tped:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-tped option requires specifying a filename. (Ex: -tped pedinput.ped).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-tped option requires specifying a filename. (Ex: -tped pedinput.ped).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-tped option missing file name?  Halting." << endl;
						exit(1);
					}

					opts::_TPEDFILE_ = arguments[++j];
					break;
				}
				case a_tfam:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-tfam option requires specifying a filename. (Ex: -tfam family.map).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-tfam option requires specifying a filename. (Ex: -tfam family.map).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-tfam option missing file name?  Halting." << endl;
						exit(1);
					}

					opts::_FAMFILE_ = arguments[++j];
					break;
				}
				case a_mdrped:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-mdrped option requires specifying a filename. (Ex: -mdrped pedinput.ped).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-mdrped option requires specifying a filename. (Ex: -mdrped pedinput.ped).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-mdrped option missing file name?  Halting." << endl;
						exit(1);
					}

					opts::_MDRPEDFILE_ = arguments[++j];
					break;
				}
				case a_mdrmap:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-mdrmap option requires specifying a filename. (Ex: -mdrmap mdrped.map).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-mdrmap option requires specifying a filename. (Ex: -mdrmap mdrped.map).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-mdrmap option missing file name?  Halting." << endl;
						exit(1);
					}

					opts::_MDRMAPFILE_ = arguments[++j];
					break;
				}
				case a_dog:
					opts::_DOG_ = true;
					opts::_CHRX_ = 39;
					opts::_CHRY_ = 40;
					opts::_CHRXY_ = 41;
					break;
				case a_map:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-map option requires specifying a filename. (Ex: -map pedmap.map).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-map option requires specifying a filename. (Ex: -map pedmap.map).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-map option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_MAPFILE_ = arguments[++j];
					break;
				}
				case a_zero_geno_file:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-zero-geno-file option requires specifying a filename. (Ex: -zero-geno-file myzeros.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-zero-geno-file option requires specifying a filename. (Ex: -zero-geno-file myzeros.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-zero-geno-file option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_ZEROGENOFILE_ = arguments[++j];
					break;
				}
				case a_mapdesc:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-mapdesc option requires specifying a filename. (Ex: -mapdesc pedmap.desc).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-mapdesc option requires specifying a filename. (Ex: -mapdesc pedmap.desc).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-mapdesc option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_MAPDESC_ = arguments[++j];
					break;
				}
				case a_sampdesc:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-sampdesc option requires specifying a filename. (Ex: -sampdesc sampleinfo.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-sampdesc option requires specifying a filename. (Ex: -sampdesc sampleinfo.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-sampdesc option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_SAMPDESC_ = arguments[++j];
					break;
				}
				case a_excsamples:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-excsamples option requires specifying a filename. (Ex: -excsamples excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-excsamples option requires specifying a filename. (Ex: -excsamples excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-excsamples option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_EXCSAMPS_ = arguments[++j];
					break;
				}
				case a_exccovs_name:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-exccovs-name option requires specifying one or more covariates.  (Ex: -exccovs-name cov1,cov2,cov3). Halting!" << endl;
						exit(0);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-exccovs-name option requires specifying one or more covariates.  (Ex: -exccovs-name cov1,cov2,cov3). Halting!" << endl;
						exit(0);
					}
					if(test[0] == '-' || test[0] == ','){
						cerr << "-exccovs-name option requires specifying one or more covariates.  (Ex: -exccovs-name cov1,cov2,cov3). Halting!" << endl;
						exit(0);
					}
					vector<string>* cov_use = new vector<string>;
					General::Tokenize(test, (*cov_use), ",");
					if(cov_use->size() == 0){
						cerr << "-exccovs-name option requires specifying one or more covariates.  (Ex: -exccovs-name cov1,cov2,cov3). Halting!" << endl;
						exit(0);
					}
					cmd_filters.add_covariate_list(cov_use);
					cmd_filters.add_covariate_filter(InputFilter::ExcludeCovariateFilter);
					j++;
					break;
				}
				case a_inccovs_name:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-inccovs-name option requires specifying one or more covariates.  (Ex: -inccovs-name cov1,cov2,cov3). Halting!" << endl;
						exit(0);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-inccovs-name option requires specifying one or more covariates.  (Ex: -inccovs-name cov1,cov2,cov3). Halting!" << endl;
						exit(0);
					}
					if(test[0] == '-' || test[0] == ','){
						cerr << "-inccovs-name option requires specifying one or more covariates.  (Ex: -inccovs-name cov1,cov2,cov3). Halting!" << endl;
						exit(0);
					}
					vector<string>* cov_use = new vector<string>;
					General::Tokenize(test, (*cov_use), ",");
					if(cov_use->size() == 0){
						cerr << "-inccovs-name option requires specifying one or more covariates.  (Ex: -inccovs-name cov1,cov2,cov3). Halting!" << endl;
						exit(0);
					}
					cmd_filters.add_covariate_list(cov_use);
					cmd_filters.add_covariate_filter(InputFilter::IncludeCovariateFilter);
					j++;
					break;
				}
				case a_exctraits_name:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-exctraits-name option requires specifying one or more traits.  (Ex: -exctraits-name trait1,trait2,trait3). Halting!" << endl;
						exit(0);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-exctraits-name option requires specifying one or more traits.  (Ex: -exctraits-name trait1,trait2,trait3). Halting!" << endl;
						exit(0);
					}
					if(test[0] == '-' || test[0] == ','){
						cerr << "-exctraits-name option requires specifying one or more traits.  (Ex: -exctraits-name trait1,trait2,trait3). Halting!" << endl;
						exit(0);
					}
					vector<string>* cov_use = new vector<string>;
					General::Tokenize(test, (*cov_use), ",");
					if(cov_use->size() == 0){
						cerr << "-exctraits-name option requires specifying one or more traits.  (Ex: -exctraits-name trait1,trait2,trait3). Halting!" << endl;
						exit(0);
					}
					cmd_filters.add_trait_list(cov_use);
					cmd_filters.add_trait_filter(InputFilter::ExcludeTraitFilter);
					j++;
					break;
				}
				case a_inctraits_name:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-inctraits-name option requires specifying one or more traits.  (Ex: -inctraits-name trait1,trait2,trait3). Halting!" << endl;
						exit(0);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-inctraits-name option requires specifying one or more traits.  (Ex: -inctraits-name trait1,trait2,trait3). Halting!" << endl;
						exit(0);
					}
					if(test[0] == '-' || test[0] == ','){
						cerr << "-inctraits-name option requires specifying one or more traits.  (Ex: -inctraits-name trait1,trait2,trait3). Halting!" << endl;
						exit(0);
					}
					vector<string>* cov_use = new vector<string>;
					General::Tokenize(test, (*cov_use), ",");
					if(cov_use->size() == 0){
						cerr << "-inctraits-name option requires specifying one or more traits.  (Ex: -inctraits-name trait1,trait2,trait3). Halting!" << endl;
						exit(0);
					}
					cmd_filters.add_trait_list(cov_use);
					cmd_filters.add_trait_filter(InputFilter::IncludeTraitFilter);
					j++;
					break;
				}
				case a_covar_file:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-covar-file option requires specifying a filename. (Ex: -covar-file covarfile.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-covar-file option requires specifying a filename. (Ex: -covar-file covarfile.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-covar-file option requires specifying a filename. (Ex: -covar-file covarfile.txt).  Halting!" << endl;
						exit(1);
					}
					opts::_COVFILE_ = arguments[++j];
					break;
				}
				case a_trait_file:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-trait-file option requires specifying a filename. (Ex: -trait-file traitfile.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-trait-file option requires specifying a filename. (Ex: -trait-file traitfile.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-trait-file option requires specifying a filename. (Ex: -trait-file traitfile.txt).  Halting!" << endl;
						exit(1);
					}
					opts::_TRAITFILE_ = arguments[++j];
					break;
				}
				case a_exccovs:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-exccovs option requires specifying a filename. (Ex: -exccovs excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-exccovs option requires specifying a filename. (Ex: -exccovs excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-exccovs option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_EXCCOVS_ = arguments[++j];
					break;
				}
				case a_inccovs:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-inccovs option requires specifying a filename. (Ex: -inccovs excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-inccovs option requires specifying a filename. (Ex: -inccovs excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-inccovs option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_INCCOVS_ = arguments[++j];
					break;
				}
				case a_exctraits:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-exctraits option requires specifying a filename. (Ex: -exctraits excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-exctraits option requires specifying a filename. (Ex: -exctraits excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-exctraits option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_EXCTRAITS_ = arguments[++j];
					break;
				}
				case a_inctraits:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-inctraits option requires specifying a filename. (Ex: -inctraits excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-inctraits option requires specifying a filename. (Ex: -inctraits excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-inctraits option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_INCTRAITS_ = arguments[++j];
					break;
				}
				case a_excmarkers:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-excmarkers option requires specifying a filename. (Ex: -excmarkers excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-excmarkers option requires specifying a filename. (Ex: -excmarkers excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-excmarkers option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_EXCMARKERS_ = arguments[++j];
					break;
				}
				case a_excfamilies:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-excfamilies option requires specifying a filename. (Ex: -excfamilies excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-excfamilies option requires specifying a filename. (Ex: -excfamilies excludelist.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-excfamilies option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_EXCFAMILIES_ = arguments[++j];
					break;
				}
				case a_incmarkers:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-incmarkers option requires specifying a filename. (Ex: -incmarkers includelist.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-incmarkers option requires specifying a filename. (Ex: -incmarkers includelist.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-incmarkers option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_INCMARKERS_ = arguments[++j];
					break;
				}
				case a_inccenters:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-map option requires specifying a filename. (Ex: -map pedmap.map).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-map option requires specifying a filename. (Ex: -map pedmap.map).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-map option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_INCCENTERS_ = arguments[++j];
					break;
				}
				case a_incsamples:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-incsamples option requires specifying a filename. (Ex: -incsamples includelist.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-incsamples option requires specifying a filename. (Ex: -incsamples includelist.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-incsamples option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_INCSAMPLES_ = arguments[++j];
					break;
				}
				case a_incfamilies:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-incfamilies option requires specifying a filename. (Ex: -incfamilies includelist.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-incfamilies option requires specifying a filename. (Ex: -incfamilies includelist.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-incfamilies option missing file name?  Halting." << endl;
						exit(1);
					}
					opts::_INCFAMILIES_ = arguments[++j];
					break;
				}
				case a_bin_input:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-bin-input option requires specifying a filename prefix. (Ex: -bin-input data).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-bin-input option requires specifying a filename prefix. (Ex: -bin-input data).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-bin-input option missing filename prefix?  Halting." << endl;
						exit(1);
					}
					opts::_BINPREFIX_ = arguments[++j];
					break;
				}
/*				case a_make_bin:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-make-bin option requires specifying a filename prefix. (Ex: -make-bin data).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-make-bin option requires specifying a filename prefix. (Ex: -make-bin data).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-make-bin option missing filename prefix?  Halting." << endl;
						exit(1);
					}
					opts::_BINPREFIX_ = arguments[++j];
					opts::_MAKEBIN_ = true;
					break;
				}
*/
				case a_out:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-out option requires specifying a filename prefix. (Ex: -out run1).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-out option requires specifying a filename prefix. (Ex: -out run1).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-out option missing filename prefix?  Halting." << endl;
						exit(1);
					}
					opts::_OUTPREFIX_ = arguments[++j];
					opts::_OUTPREFIX_ += ".";
					break;
				}
				case a_micro_sats:
					opts::_MICROSATS_ = true;
					break;
				case a_flip:
				{
					opts::_FLIPSTRAND_ = true;

					if(j + 1 >= arguments.size()){
						cerr << "-flip option requires specifying a filename. (Ex: -flip strand_info.txt).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-flip option requires specifying a filename. (Ex: -flip strand_info.txt).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-flip option missing filename?  Halting." << endl;
						exit(1);
					}
					opts::_FLIPFILE_ = arguments[++j];
					break;
				}
				case a_chrom:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-chrom option requires specifying a chromosome. (Ex: -chrom 12).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-chrom option requires specifying a chromosome. (Ex: -chrom 12).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-chrom option missing value?  Halting." << endl;
						exit(1);
					}
					opts::_CHROM_LIMIT_ = true;
					string val = arguments[++j];
					if(val != "X" && val != "x" && val != "Y" && val != "y"){
						opts::_CHROM_ = atoi(val.c_str());
					}
					else if(val == "X" || val == "x"){
						opts::_CHROM_ = opts::_CHRX_;
						vector<string>::iterator df = find(arguments.begin(), arguments.end(), "-dog");
						if(df != arguments.end()){
							opts::_CHROM_ = 39;
						}
					}
					else if(val == "Y" || val == "y"){
						opts::_CHROM_ = opts::_CHRY_;
						vector<string>::iterator df = find(arguments.begin(), arguments.end(), "-dog");
						if(df != arguments.end()){
							opts::_CHROM_ = 40;
						}
					}
					else if(val == "Xy" || val == "XY" || val == "xY" || val == "xy"){
						opts::_CHROM_ = opts::_CHRXY_;
						vector<string>::iterator df = find(arguments.begin(), arguments.end(), "-dog");
						if(df != arguments.end()){
							opts::_CHROM_ = 41;
						}
					}
					break;
				}
				case a_bp_min:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-bp-min option requires specifying a positive number. (Ex: -bp-min 1234567).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-bp-min option requires specifying a positive number. (Ex: -bp-min 1234567).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-bp-min option missing value?  Halting." << endl;
						exit(1);
					}
					opts::_BP_LOW_LIMIT_ = true;
					opts::_BP_LOW_ = atoi(arguments[++j].c_str());
					break;
				}
				case a_bp_max:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-bp-max option requires specifying a value. (Ex: -bp-max 987654321).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-bp-max option requires specifying a value. (Ex: -bp-max 987654321).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-bp-max option missing value?  Halting." << endl;
						exit(1);
					}
					opts::_BP_HIGH_LIMIT_ = true;
					opts::_BP_HIGH_ = atoi(arguments[++j].c_str());
					break;
				}
				case a_bp_space:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-bp-space option requires specifying a value. (Ex: -bp-space 1000000).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-bp-space option requires specifying a value. (Ex: -bp-space 1000000).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-bp-space option missing value?  Halting." << endl;
						exit(1);
					}
					opts::_BP_SPACE_LIMIT_ = true;
					opts::_BP_SPACE_ = atoi(arguments[++j].c_str());
					break;
				}
				case a_threads:
				{
					opts::_THREADS_ = true;
					break;
				}
				case a_numthreads:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-numthreads option requires specifying a value. (Ex: -numthreads 4).  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-numthreads option requires specifying a value. (Ex: -numthreads 4).  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-numthreads option missing value?  Halting." << endl;
						exit(1);
					}
					opts::_THREADS_ = true;
					opts::_NUMTHREADS_ = atoi(arguments[++j].c_str());
					break;
				}
				case a_freq_file:
				{
					if(j + 1 >= arguments.size()){
						cerr << "-freq-file option requires specifying a filename. (Ex: -freq-file mymaffile.txt.  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0){
						cerr << "-freq-file option requires specifying a filename. (Ex: -freq-file mymaffile.txt.  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-'){
						cerr << "-freq-file option requires specifying a filename. (Ex: -freq-file mymaffile.txt.  Halting!" << endl;
						exit(1);
					}
					opts::_FREQ_FILE_EXISTS_ = true;
					opts::_FREQ_FILE_ = arguments[++j];
					break;
				}
				//added 02-23-2011 to support new feature to update FamID and IndID based on file input
				case a_update_ids:
				{
					if(j + 1 >= arguments.size())
					{
						cerr << "-update-ids option requires specifying a filename. (Ex: -update-ids newIDsFile.txt.  Halting!" << endl;
						exit(1);
					}
					string test = arguments[j+1];
					if(test.length() == 0)
					{
						cerr << "-update-ids option requires specifying a filename. (Ex: -update-ids newIDsFile.txt.  Halting!" << endl;
						exit(1);
					}
					if(test[0] == '-')
					{
						cerr << "-update-ids option requires specifying a filename. (Ex: -update-ids newIDsFile.txt.  Halting!" << endl;
						exit(1);
					}
					opts::_ID_FILE_EXISTS_ = true;
					opts::_ID_FILE_ = arguments[++j];
					break;
				}
				case a_remaining:
					opts::_OUTPUTREMAINS_ = true;
					break;
				case a_print_families:
					opts::_PRINTFAMS_ = true;
					break;
				case a_no_summary:
					opts::_COMPILE_OUTPUTS_ = false;
					break;
				case a_todigit:
					opts::_TODIGIT_ = true;
					break;
				default:
					cerr << "Unknown argument: " << arguments[j] << endl;
					exit(1);
					break;
			}
		}
		string logoutput = opts::_OUTPREFIX_ + "plato.log";
//		ofstream log(logoutput.c_str());
		opts::_LOG_.open(logoutput.c_str());
		if(!opts::_LOG_){
			cerr << "Error opening plato.log.  Exiting!\n";
			exit(1);
		}
		string args = "";
		for(int j = 0; j < argc; j++){
			string t = argv[j];
			args += t + " ";
//			log << argv[j] << " ";
		}
		args += "\n";
		opts::printLog(args);
//		log << endl;
//		log.close();

		if(opts::_PEDFILE_.length() == 0 && opts::_BINPREFIX_.length() == 0 && opts::_TPEDFILE_.length() == 0
				&& opts::_MDRPEDFILE_.length() == 0){
			opts::printLog("No pedfile, binary, transposed, or mdr-style input files specified.  Exiting!\n");
			exit(0);
		}
		error_check();
		//steps = initializeSteps();
		if(!opts::_MAKEBIN_){
			proc_order = parseInput(arg);//, &steps);
		}

	}
	else{
		usage();
		exit(1);
	}
	time_t curr = time(0);
	string tdstamp = ctime(&curr);
	opts::printLog("\nPlato started: " + tdstamp + "\n");

		//Begin batch file steps and reading in data
		startProcess(&proc_order, NULL, myrank, &cmd_filters);
	}catch(MethodException & ex){
		opts::printLog(string(ex.what()) + "\n");
	}
//	catch(std::exception ex){
//		opts::printLog(string(ex.what()) + "\n");
//	}
	time_t curr = time(0);
	string tdstamp = ctime(&curr);
	opts::printLog("\nPlato finished: " + tdstamp + "\n");
	//	cerr << "Unknown exception caught..." << endl;
	//	exit(0);
	//}
    return 1;
}

/*
 *Function: error_check
 *Return: void
 *Description:
 *Performs basic error checking on completeness of command line arguments.
 *
 */
void error_check(){
	if((opts::_BP_LOW_LIMIT_ || opts::_BP_HIGH_LIMIT_) && !opts::_CHROM_LIMIT_){
		opts::printLog("-chrom is required when specifying -bp-min and/or -bp-max.  Exiting!\n");
		exit(0);
	}
}

/*
 *Function: print_help
 *Return: void
 *Description:
 *Outputs general help statement when -h option is used.
 *
 */
void print_help(){
	//int w = 5;
	cout << "plato - PLatform for the Analysis, Translation, and Organization of large-scale data (Version: " << _WASPVER_ << ")" << endl << endl;
	cout << "Center for Human Genetics Research" << endl;
	cout << "Vanderbilt University Medical Center" << endl << endl;
	cout << "!!!!  For most up-to-date options and details  !!!!" << endl;
	cout << "!!!!  visit http://chgr.mc.vanderbilt.edu/plato !!!!" << endl << endl;
/*
	cout << "usage: plato <batchfile> [<option1> <option2>...]" << endl << endl;
	cout << "   Special options: " << endl;
	cout << "     -h                 Prints this help" << endl;
	cout << "     -S                 Displays list of valid processing steps" << endl << endl;
	cout << "   Options: " << endl;
	cout << "     -ped {pedfile}     Specify pedfile <fam> <ind> <dad> <mom> <sex> <aff> <genotypes...>" << endl;
	cout << "     -map {mapfile}     Specify mapfile <chrom> <rsid> <bploc>" << endl << endl;
	cout << "     -tped {tpedfile}   Specify tpedfile <chr> <snp> <cm> <bploc> <genotypes...>" << endl;
	cout << "     -tfam {famfile}    Specify family map file for transposed set" << endl;
	cout << "              Format: <fam> <ind> <dad> <mom> <gend> <aff>" << endl;
	cout << "     -bin-input {file_prefix} Read inputs in binary format (based on -make-bin output) (optional)" << endl;
	cout << "     -micro-sats        Use if maker(s) have more than 2 alleles." << endl;
	cout << "     -mapdesc {mapdescriptive_file} Specify more details about map file (optional)" << endl;
	cout << "                                    Mapfile descriptor format: Chrom Marker Basepair_location Other_marker_id Enzyme" << endl;
	cout << "     -sampdesc {sample_Descriptive_file} Specify more details about samples (optional)" << endl;
	cout << "                                         Sample descriptor format: Family Ind Center Plate Well" << endl << endl;
	cout << "     -excsamples {samplefile} File containing samples to exclude (optional)" << endl;
	cout << "     -excmarkers {markerfile} File containing markers to exclude (optional)" << endl << endl;
	cout << "     -incsamples {samplefile} File containing samples to include (optional)" << endl;
	cout << "     -incmarkers {markerfile} File containing markers to include (optional)" << endl << endl;
	cout << "     -make-bin {name}         Convert Ped/Map files into binary input file format with {name} as the file prefix (optional)" << endl;
	cout << endl;
	cout << "     -out {name}         Prepends all output files with {name} (optional)" << endl;
	cout << endl;
	cout << "     -freq-file {frequencyfile} File containing predefined minor allele frequencies (optional)" << endl;
	cout << "                                Format: Marker_id Minor_allele_frequency" << endl;
	cout << "     -flip {filename}  Specifies a list of markers to flip to the opposite strand (optional)" << endl;
	cout << endl;
	cout << "     -chrom {chrom}  Limit markers read to specified chromosome (optional)" << endl;
	cout << "     -bp-min {number}  Specifies minimum chromosome specific base-pair location to load (optional)" << endl;
	cout << "     -bp-max {number} Specifies maximum chromosome specific base-pair location to load (optional)" << endl;
	cout << "     -bp-space {number} Specifies to load only markers spaced by {number} base-pair across all chromosomes. (optional)" << endl;
	cout << endl;
	cout << endl;
	cout << "   Step specific options (to be used in the batch file):" << endl;
	cout << "     -chrom {chrom}  Valid on all steps." << endl;
	cout << "     -bp-min {number}  Valid on all steps." << endl;
	cout << "     -bp-max {number}  Valid on all steps." << endl;
	cout << "     -bp-space {number}  Valid on all steps." << endl;
	cout << "     -thresh-markers-max {number}  Valid on marker-geno-eff, family-geno-eff, allele-freq, mendelian-error, hw, gender-error, tdt, chisquare.." << endl;
	cout << "     -thresh-markers-min {number}  Valid on marker-geno-eff, family-geno-eff, allele-freq, mendelian-error, hw, gender-error, tdt, chisquare." << endl;
	cout << "     -thresh-samples-max {chrom}  Valid on sample-geno-eff, mendelian-error." << endl;
	cout << "     -thresh-samples-min {chrom}  Valid on sample-geno-eff, mendelian-error." << endl;
	cout << "     -thresh-families-max {chrom}  Valid on mendelian-error." << endl;
	cout << "     -thresh-families-min {chrom}  Valid on mendelian-error." << endl;
	cout << "     -zero  Valid on mendelian-error." << endl;
	cout << "     -zero-l2  Valid on mendelian-error." << endl;
	cout << "     -zero-l2-fams Valid on mendelian-error." << endl;
	cout << "     -penetrance-file {filename} Valid on output-superlink." << endl;
	cout << "     -rm-mono  Valid on allele-freq." << endl;
	cout << "     -rm-het-only Valid on allele-freq." << endl;
	cout << "     -center-file {filename}  Valid on output-lapis." << endl;
	cout << "     -out {fileprefix}  Valid on all steps." << endl;
	cout << "     -covar-file {filename}" << endl;
	cout << "     -covars-name {list of headers}" << endl;
	cout << "     -covars-number {list of header numbers}" << endl;
	cout << "     -trait-file {filename}" << endl;
	cout << "     -traits-name {list of headers}" << endl;
	cout << "     -traits-number {list of header numbers}" << endl;
	cout << "     -covar-missing {missing value}" << endl;
	cout << "     -trait-missing {missing value}" << endl;
	cout << "     -group-file {filename}" << endl;
	cout << "     -ped {filename}  Valid on concordance" << endl;
	cout << "     -map {filename}  Valid on concordance" << endl;
	cout << "     -tped {filename} Valid on concordance" << endl;
	cout << "     -tfam {filename} Valid on concordance" << endl;
	cout << "     -bin-input {fileprefix} Valid on concordance" << endl;
	cout << "     -strat-file {filename}  Valid on output-structure." << endl;
	cout << "     -parents-only  Valid on output-structure." << endl;
	cout << "     -trio  Valid on output-phase." << endl;
	cout << "     -disease {name}  Valid on output-beagle." << endl;
	cout << "     -homozyg-zeros {number}  Valid on homozygous." << endl;
	cout << "     -homozyg-span {number}  Valid on homozygous." << endl;
	cout << "     -homozyg-seq-prob {number}  Valid on homozygous." << endl;
	cout << "     -homozyg-wgha  Alternate homozygous span algorithm.  Valid on homozygous." << endl;
	cout << "     -homozyg-min-samp {number}  Valid on homozygous." << endl;
	cout << "     -ld-pairwise {number} {number} {number} Valid on ld." << endl;
	cout << "     -ld-vif {number} {number} {number}  Valid on ld." << endl;
	cout << "     -ld-chop  Valid on ld." << endl;
	cout << "     -deletion-span {number}   Valid on deletions." << endl;
	cout << "     -deletion {number}  Valid on deletions." << endl;
	cout << "     -filter-overall  Valid on allele-freq." << endl;
	cout << "     -filter-file  Valid on allele-freq." << endl;
	cout << "     -parental    Provides further breakdown. Valid on allele-freq and hw." << endl;
	cout << "     -gender      Provides further breakdown. Valid on allele-freq and hw." << endl;
	cout << "     -casecontrol Provides further breakdown. Valid on allele-freq and hw." << endl;
	cout << "     -all  Override overall calc to include all samples.  Valid on allele-freq and hw." << endl;
	cout << "     -all-children  Override to include only all children.  Valid on allele-freq and hw." << endl;
	cout << "     -random-child  Override to include one random child from each family." << endl;
	cout << "                    Valid on allele-freq and hw." << endl;
	cout << "     -unk-spouses   Override to use only spouses with unknown phenotype." << endl;
	cout << "                    Valid on allele-freq and hw." << endl;
	cout << "     -unaff-spouses-only  Override to use only spouses with unaffected phenotype." << endl;
	cout << "                    Valid on allele-freq and hw." << endl;
	cout << "     -no-summary    Disables production of marker, sample, family summary files." << endl;
	cout << endl;
	cout << endl << endl;
	cout << "BATCH FILE FORMATTING:" << endl << endl;
	cout << "The batch file consists of a list of steps and options for the process" << endl;
	cout << "to follow in order of appearance in the file." << endl << endl;
	cout << "Example batch file:" << endl;
	cout << "marker-geno-eff -thresh-markers-min 90" << endl;
	cout << "allele-freq -thresh-markers-min 0.1" << endl;
	cout << endl;
*/
}

/*
 * printOptions
 *
 * outputs command line options
 *
 */
void printOptions(){
	map<string, cmdArgs>::iterator iter;

	cout << "Command line arguments:\n-------------------\n";
	for(iter = s_mapcmdArgs.begin(); iter != s_mapcmdArgs.end(); iter++){
		cout << iter->first << endl;
	}
}

/*
 *Function: print_steps
 *Return: void
 *Parameters: STEPS object
 *Description:
 *Outputs valid batch file steps and their descriptions when the -S command line argument is used
 */
void print_steps(STEPS s){
	STEPS::iterator s_iter;
	unsigned int field = 0;
	for(s_iter = s.begin(); s_iter != s.end(); s_iter++){
		Step mystep =  s_iter->second;
		if(s_iter->first.size() > field){
			field = s_iter->first.size();
		}
	}
	field+=5;
	cout << left << setw(field) << "Step:" << "Description:" << endl;
	cout << left << setw(field) << "----------" << "---------------" << endl;
	for(s_iter = s.begin(); s_iter != s.end(); s_iter++){
		Step mystep =  s_iter->second;
		cout << left << setw(field) << s_iter->first <<  mystep.getName() << endl;
		//cout << "\tThresh: " << mystep.getThreshold() << endl;
	}
}

/*
 *Function: usage
 *Description:
 *outputs general usage statement
 */
void usage(){
	cout << "plato - PLatform for the Analysis, Translation, and Organization of large-scale data (Version: " << _WASPVER_ << ")" << endl << endl;
	cout << "Center for Human Genetics Research" << endl;
	cout << "Vanderbilt University Medical Center" << endl << endl;
	cout << "!!!!  For most up-to-date options and details  !!!!" << endl;
	cout << "!!!!  visit http://chgr.mc.vanderbilt.edu/plato !!!!" << endl << endl;
	cout << "usage: plato <batchfile> [<option1> <option2>....]" << endl
		 << endl;
	cout << "For a list of valid steps to be inserted into the batch file:" << endl;
	cout << "\t\t> plato -S" << endl << endl;
	cout << "For help: " << endl;
	cout << "\t\t> plato -h" << endl << endl;
}

/*
 *Function: startProcess
 *Parameters:
 *ORDER object
 *database connection (null)
 *rank
 *
 *Description:
 *Starts wasp process
 */
void startProcess(ORDER* order, void* con, int myrank, InputFilter* filters){
	ORDER::iterator o_iter;
	vector<ORDER> optimized = optimize(order);

	DataSet data_set;
		//set up filters
//		InputFilter filters;
		if(opts::_EXCCOVS_.length() > 0){
			vector<string>* list = new vector<string>;
			Helpers::readCovTraitFile(opts::_EXCCOVS_, list);
			filters->add_covariate_list(list);
			filters->add_covariate_filter(InputFilter::ExcludeCovariateFilter);
		}
		if(opts::_INCCOVS_.length() > 0){
			vector<string>* list = new vector<string>;
			Helpers::readCovTraitFile(opts::_INCCOVS_, list);
			filters->add_covariate_list(list);
			filters->add_covariate_filter(InputFilter::IncludeCovariateFilter);
		}
		if(opts::_EXCTRAITS_.length() > 0){
			vector<string>* list = new vector<string>;
			Helpers::readCovTraitFile(opts::_EXCTRAITS_, list);
			filters->add_covariate_list(list);
			filters->add_covariate_filter(InputFilter::ExcludeTraitFilter);
		}
		if(opts::_INCTRAITS_.length() > 0){
			vector<string>* list = new vector<string>;
			Helpers::readCovTraitFile(opts::_INCTRAITS_, list);
			filters->add_covariate_list(list);
			filters->add_covariate_filter(InputFilter::IncludeTraitFilter);
		}
		if(opts::_EXCMARKERS_.length() > 0){
			vector<Marker*>* list = new vector<Marker*>;
			Helpers::readLocusFile(opts::_EXCMARKERS_, list);
			filters->add_locus_list(list);
			filters->add_locus_filter(InputFilter::ExcludeLocusFilter);
		}
		if(opts::_INCMARKERS_.length() > 0){
			vector<Marker*>* list = new vector<Marker*>;
			Helpers::readLocusFile(opts::_INCMARKERS_, list);
			filters->add_locus_list(list);
			filters->add_locus_filter(InputFilter::IncludeLocusFilter);
		}
		if(opts::_EXCSAMPS_.length() > 0){
			vector<Sample*>* list = new vector<Sample*>;
			Helpers::readSampleFile(opts::_EXCSAMPS_, list);
			filters->add_sample_list(list);
			filters->add_sample_filter(InputFilter::ExcludeSampleFilter);
		}
		if(opts::_INCSAMPLES_.length() > 0){
			vector<Sample*>* list = new vector<Sample*>;
			Helpers::readSampleFile(opts::_INCSAMPLES_, list);
			filters->add_sample_list(list);
			filters->add_sample_filter(InputFilter::IncludeSampleFilter);
		}
		if(opts::_EXCFAMILIES_.length() > 0){
			vector<Family*>* list = new vector<Family*>;
			Helpers::readFamilyFile(opts::_EXCFAMILIES_, list);
			filters->add_family_list(list);
			filters->add_family_filter(InputFilter::ExcludeFamilyFilter);
		}
		if(opts::_INCFAMILIES_.length() > 0){
			vector<Family*>* list = new vector<Family*>;
			Helpers::readFamilyFile(opts::_INCFAMILIES_, list);
			filters->add_family_list(list);
			filters->add_family_filter(InputFilter::IncludeFamilyFilter);
		}
		if(opts::_CHROM_LIMIT_){
			Marker* temp = new Marker();
			temp->setChrom(opts::_CHROM_);
			vector<Marker*>* list = new vector<Marker*>;
			list->push_back(temp);
			filters->add_locus_list(list);
			filters->add_locus_filter(InputFilter::LocusChromFilter);
		}
		if(opts::_BP_LOW_LIMIT_ || opts::_BP_HIGH_LIMIT_){
			vector<Marker*>* list = new vector<Marker*>;
			if(opts::_BP_LOW_LIMIT_){
				Marker* temp = new Marker();
				temp->setBPLOC(opts::_BP_LOW_);
				list->push_back(temp);
			}
			else{
				list->push_back(NULL);
			}
			if(opts::_BP_HIGH_LIMIT_){
				Marker* temp = new Marker();
				temp->setBPLOC(opts::_BP_HIGH_);
				list->push_back(temp);
			}
			else{
				list->push_back(NULL);
			}
			filters->add_locus_list(list);
			filters->add_locus_filter(InputFilter::LocusBplocRangeFilter);
		}
		//set options
		StepOptions options;
		options.setMapFile(opts::_MAPFILE_);
		options.setMapContainsReferent(opts::_MAP_INCLUDES_REF_);
		options.setPedFile(opts::_PEDFILE_);
		options.setBinInput(opts::_BINPREFIX_);
		options.setTFamFile(opts::_FAMFILE_);
		options.setTPedFile(opts::_TPEDFILE_);
		options.setMdrMapFile(opts::_MDRMAPFILE_);
		options.setMdrPedFile(opts::_MDRPEDFILE_);
		options.setCovarMissing(opts::_COVAR_MISSING_);
		options.setTraitMissing(opts::_TRAIT_MISSING_);

		if(opts::_SAMPLEBPRANGEFILTER_.length() > 0){
			options.setSampleBprangeFile(opts::_SAMPLEBPRANGEFILTER_);
			opts::printLog("Reading sample/bprange file: " + opts::_SAMPLEBPRANGEFILTER_ + "\n");
			options.readSampleBprangeFile();
		}

		if(opts::_FREQ_FILE_EXISTS_){
			opts::printLog("Reading marker frequencies from: " + opts::_FREQ_FILE_ + "\n");
			options.setFrequencies(Helpers::readFreqFile(opts::_FREQ_FILE_));
		}

		//read zerogenofile
		if(opts::_ZEROGENOFILE_.length() > 0){
			Helpers::readZeroGenoFile(opts::_ZEROGENOFILE_);
		}
		//new 12-06-2010
		if(opts::_LGENFILE_.length() > 0)
		{
			if(opts::_PEDFILE_.length() ==0 || opts::_MAPFILE_.length() ==0)
			{
				opts::printLog("You specified -lgen-file.  Corresponding ped and map files must be specified using -ped and -map.  Exiting!\n");
				exit(1);
			}
		}

		//read Binary input
		if(opts::_BINPREFIX_.length() > 0 && !opts::_MAKEBIN_){
			if(opts::_PEDINFO_.length() > 0){
				opts::printLog("Reading Pedigree information file: " + opts::_PEDINFO_ + "\n");
				Helpers::readPedInfo();
			}
			opts::printLog("Reading data using Binary inputs: prefix: " + opts::_BINPREFIX_ + "\n");
			if(opts::_MICROSATS_){
				opts::printLog("You specified microsatellite markers exist.  This option is not compatible with binary input files.  Please use the standard PED file format for your input when using microsatellite markers.  Exiting...\n");
				exit(1);
			}
			Helpers::readBinM(&data_set, options, filters);//(data_set.get_samples(), data_set.get_families(), data_set.get_markers(), data_set.get_marker_map());
			//assignLinks(data_set.get_families());
			//reorderAlleles(data_set.get_samples(), data_set.get_markers());
			if(opts::_FLIPSTRAND_){
				flipStrand(data_set.get_markers());
			}
		}
		//read ped & map file
		else if((opts::_PEDFILE_.length() > 0 || opts::_MAPFILE_.length() > 0) && opts::_TPEDFILE_.length() == 0
				&& opts::_MDRPEDFILE_.length() == 0){
			if(opts::_MAPFILE_.length() > 0){
				opts::printLog("Reading MAP file: " + opts::_MAPFILE_ + "\n");
				//readMap(data_set.get_markers(), data_set.get_marker_map());
				Helpers::readMapM(&data_set, options, filters);
			}
			else{
				opts::printLog("You specified -ped.  A corresponding map file needs to be specified using -map. Exiting!\n");
				exit(1);
			}
			if(opts::_PEDFILE_.length() > 0){
				if(opts::_PEDINFO_.length() > 0){
					opts::printLog("Reading Pedigree information file: " + opts::_PEDINFO_ + "\n");
					Helpers::readPedInfo();
				}
				opts::printLog("Reading PED file: " + opts::_PEDFILE_ + "\n");
				//readPed(data_set.get_samples(), data_set.get_families(), data_set.get_markers(), data_set.get_marker_map());
				//assignLinks(data_set.get_families());
				//reorderAlleles(data_set.get_samples(), data_set.get_markers());
				Helpers::readPedM_3vec_set(&data_set, options, filters);
				if(opts::_FLIPSTRAND_){
					flipStrand(data_set.get_markers());
				}
				//new 12-06-2010
				if(opts::_LGENFILE_.length() > 0)
				{
					opts::printLog("Reading LGEN file: " + opts::_LGENFILE_ + "\n");
					Helpers::readLgenFile(&data_set, options, filters);
					if(opts::_FLIPSTRAND_)
					{
						flipStrand(data_set.get_markers());
					}
				}
			}
			else{
				opts::printLog("You specified -map.  A corresponding PED file needs to be specified using -ped. Exiting!\n");
				exit(1);
			}
		}
		else if(opts::_PEDFILE_.length() == 0 && opts::_MDRPEDFILE_.length() == 0 &&
				(opts::_TPEDFILE_.length() > 0 || opts::_FAMFILE_.length() > 0)){
			if(opts::_FAMFILE_.length() > 0){
				if(opts::_PEDINFO_.length() > 0){
					opts::printLog("Reading Pedigree information file: " + opts::_PEDINFO_ + "\n");
					Helpers::readPedInfo();
				}
				opts::printLog("Reading Pedigree file: " + opts::_FAMFILE_ + "\n");
				Helpers::readTFamM(&data_set, options, filters);
				//assignLinks(data_set.get_families());
			}
			else{
				opts::printLog("You specified -tped.  A corresponding Pedigree Information file needs to be specified using -tfam. Exiting!\n");
				exit(1);
			}
			if(opts::_TPEDFILE_.length() > 0){
				opts::printLog("Reading TPED file: " + opts::_TPEDFILE_ + "\n");
				Helpers::readTPedM(&data_set, options, filters);//data_set.get_markers(), data_set.get_samples(), data_set.get_marker_map());
				//reorderAlleles(data_set.get_samples(), data_set.get_markers());
				if(opts::_FLIPSTRAND_){
					flipStrand(data_set.get_markers());
				}
			}
			else{
				opts::printLog("You specified -tfam.  A corresponding transposed PED file needs to be specified using -tped. Exiting!\n");
				exit(1);
			}
		}
		else if(opts::_MDRPEDFILE_.length() > 0 || opts::_MDRMAPFILE_.length() > 0){
			if(opts::_MDRMAPFILE_.length() > 0){
				opts::printLog("Reading MDRMAP file: " + opts::_MDRMAPFILE_ + "\n");
				//readMap(data_set.get_markers(), data_set.get_marker_map());
				Helpers::readMapMdr(&data_set, options, filters);

			}
			else{
				opts::printLog("You specified -mdrped.  A corresponding MAP file needs to be specified using -mdrmap. Exiting!\n");
				exit(1);
			}

			if(opts::_MDRPEDFILE_.length() > 0){
				if(opts::_PEDINFO_.length() > 0){
					opts::printLog("Reading Pedigree information file: " + opts::_PEDINFO_ + "\n");
					Helpers::readPedInfo();
				}
				opts::printLog("Reading MDRPED file: " + opts::_MDRPEDFILE_ + "\n");
				//readPed(data_set.get_samples(), data_set.get_families(), data_set.get_markers(), data_set.get_marker_map());
				//assignLinks(data_set.get_families());
				//reorderAlleles(data_set.get_samples(), data_set.get_markers());
				Helpers::readMdr(&data_set, options, filters);
				if(opts::_FLIPSTRAND_){
					flipStrand(data_set.get_markers());
				}

			}
			else{
				opts::printLog("You specified -mdrmap.  A corresponding MDR-PED file needs to be specified using -mdrped. Exiting!\n");
				exit(1);
			}
		}
		else if(!opts::_DBINPUT_){
			cerr << "No input method specified!" << endl;
			exit(1);
		}
		if(opts::_COVFILE_.length() > 0){
			opts::printLog("Reading covariate file: " + opts::_COVFILE_ + "\n");
			Helpers::readCovariateFile(opts::_COVFILE_, &data_set, options, filters);
		}
		if(opts::_TRAITFILE_.length() > 0){
			opts::printLog("Reading trait file: " + opts::_TRAITFILE_ + "\n");
			Helpers::readTraitFile(opts::_TRAITFILE_, &data_set, options, filters);
		}
		if(opts::_ID_FILE_EXISTS_)
		{
			opts::printLog("Reading ID file: " + opts::_ID_FILE_ + "\n");
			Helpers::readIDFile(opts::_ID_FILE_, &data_set);
		}

		if(opts::_SAMPLEBPRANGEFILTER_.length() > 0){
			opts::printLog("Running sample/bprange filter...\n");
			filters->run_sample_bprange_filter(data_set.get_samples(), data_set.get_markers(), (options.getSampleBprangeSamples()), (options.getSampleBprangeMarkers()));
		}

		//int goodsamps = 0;
		//int goodfams = 0;
		//int goodmarkers = 0;
		if(opts::_TODIGIT_){
			opts::printLog("Mapping Families/Samples to numeric format.");
			Helpers::remapFamsToDigit(data_set.get_families());
			Helpers::printFamsToDigit(data_set.get_families(), "global", options);
		}
		for(unsigned int i = 0; i < data_set.get_samples()->size(); i++){
			Sample* samp = data_set.get_sample(i);
			if(samp->isEnabled()){
				//goodsamps++;
				opts::_SAMPLES_WORKING_++;
			}
		}
		for(unsigned int i = 0; i < data_set.get_families()->size(); i++){
			Family* fam = data_set.get_pedigree(i);
			if(fam->isEnabled()){
				//goodfams++;
				opts::_FAMILIES_WORKING_++;
			}
		}
		for(unsigned int i = 0; i < data_set.get_markers()->size(); i++){
			Marker* mark = data_set.get_locus(i);
			if(mark->isEnabled()){
				//goodmarkers++;
				opts::_MARKERS_WORKING_++;
			}
		}
		opts::_COVS_WORKING_ = data_set.num_covariates();
		opts::_TRAITS_WORKING_ = data_set.num_traits();
		string text = "Markers found: ";
		text += getString<int>(opts::_MARKERS_FOUND_);
		text += ", ";
		text += getString<int>(opts::_MARKERS_WORKING_);
	   	text += " are enabled!\n";
		opts::printLog(text);
		text = "Samples found: ";
	   	text += getString<int>(opts::_SAMPLES_FOUND_);
	   	text += ", ";
	   	text += getString<int>(opts::_SAMPLES_WORKING_);
	   	text += " are enabled!\n";
		opts::printLog(text);
		text = "Families found: ";
	    text += getString<int>(opts::_FAMILIES_FOUND_);
	   	text += ", ";
	   	text += getString<int>(opts::_FAMILIES_WORKING_);
	    text += " are enabled!\n";
	    if(opts::_COVS_FOUND_ > 0){
	    	text += "Covariates found: " + getString<int>(opts::_COVS_FOUND_) + ", " + getString<int>(opts::_COVS_WORKING_) + " are enabled!\n";
	    }
	    if(opts::_TRAITS_FOUND_ > 0){
	    	text += "Traits found: " + getString<int>(opts::_TRAITS_FOUND_) + ", " + getString<int>(opts::_TRAITS_WORKING_) + " are enabled!\n";
	    }
		int founders = 0;
		int nonfounders = 0;
		for(unsigned int f = 0; f < data_set.get_families()->size(); f++){
			Family* fam = data_set.get_pedigree(f);
			if(fam->isEnabled()){
				vector<Sample*>* founder = fam->getFounders();
				vector<Sample*>* nonfound = fam->getNonFounders();
				for(unsigned int s = 0; s < founder->size(); s++){
					Sample* samp = (*founder)[s];
					if(samp->isEnabled()){
						founders++;
					}
				}
				for(unsigned int s = 0; s < nonfound->size(); s++){
					Sample* samp = (*nonfound)[s];
					if(samp->isEnabled()){
						nonfounders++;
					}
				}
			}
		}
		text += "Founders found: ";
		text += getString<int>(founders);
		text += "\n";
		opts::printLog(text);
		opts::printLog("Non-founders found: " + getString<int>(nonfounders) + "\n");

		if(opts::_MARKERS_WORKING_ == 0){
			opts::printLog("No enabled markers are found.  Please check your input files and arguments.\n");
			exit(1);
		}
		if(opts::_SAMPLES_WORKING_ == 0){
			opts::printLog("No enabled samples are found.  Please check your input files and arguments.\n");
			exit(1);
		}
		if(opts::_FAMILIES_WORKING_ == 0){
			opts::printLog("No enabled families are found.  Please check your input files and arguments.\n");
			exit(1);
		}

		if(opts::zerogenoinfo.size() > 0){
			opts::printLog("Zero-ing specified Sample/SNP pairs.\n");
			Helpers::zeroSingleGenos(data_set.get_markers(), data_set.get_samples());
		}

//		if(opts::_MAKEBIN_){
//			writeBit(data_set.get_samples(), data_set.get_families(), data_set.get_markers(), data_set.get_marker_map());
//		}
		if(opts::_PRINTFAMS_){
			opts::printLog("Printing family structure diagram.\n");
			printFamilies(data_set.get_families());
		}
		data_set.create_marker_name_map();
	//}
	//set_me_up->summary();

	//Start processing steps by slaves.  MASTER gathers all data from slaves and does summary work
	//threading here?

//		if(opts::_THREADS_){
//		}
		int count = 0;
		vector<boost::thread*> threads_hash;
		//vector<boost::thread*>::iterator thread_iter;
//		boost::thread *threads[opts::_NUMTHREADS_];

	map<int, int> thread_step_map;
	ORDER yesthread = optimized[0];
	ORDER nothread = optimized[1];

	for(unsigned int o = 0; o < yesthread.size(); o++){//o_iter = order->begin(); o_iter != order->end(); o_iter++){
		Step current_step = yesthread.at(o);//(Step)(*o_iter);//(*steps)[(*o_iter)];
//		runStep(current_step, &data_set);

		if(opts::_THREADS_){
//			thread_step_map[count] = o;
//			threads[count++] = new boost::thread(boost::bind(&runStep, current_step, &data_set));
			threads_hash.push_back(new boost::thread(boost::bind(&runStep, current_step, &data_set)));
//			cout << "thread count " << count << endl;
			count++;
			if(count >= opts::_NUMTHREADS_){
//				for(int i = 0; i < count; i++){
/////				for(thread_iter = threads_hash.begin(); thread_iter < threads_hash.end(); thread_iter++){
				while(count >= opts::_NUMTHREADS_){
				for(int i = 0; i < (int)threads_hash.size(); i++){
//					if(threads[i]->timed_join(boost::posix_time::milliseconds(500))){
					if(threads_hash[i]->timed_join(boost::posix_time::milliseconds(500))){
/////					finalizeStep(order->at(thread_step_map[i]));
						delete(threads_hash[i]);//thread_iter->first]);
//						thread_step_map.erase(count - 1);
						threads_hash.erase(threads_hash.begin() + i);//thread_iter->first);
						count--;
						i--;
					}
				}
				}
//				count = 0;
//				thread_step_map.clear();
			}
		}
		else{
			runStep(current_step, &data_set);
		}
	}

	if(opts::_THREADS_){
	//for(thread_iter = threads_hash.begin(); thread_iter < threads_hash.end(); thread_iter++){//int i = 0; i < count; i++){
	for(int i = 0; i < (int)threads_hash.size(); i++){
		threads_hash[i]->join();//thread_iter->first]->join();
			delete(threads_hash[i]);//thread_iter->first]);

	}
	}

	for(unsigned int o = 0; o < nothread.size(); o++){
		Step current_step = nothread.at(o);
		runStep(current_step, &data_set);
	}



	if(myrank == 0){
		if(opts::_OUTPUTREMAINS_){
			Finalize* f = new Finalize();
			f->finish(data_set.get_markers(), data_set.get_samples(), data_set.get_families());
			delete(f);
		}
		if(opts::_COMPILE_OUTPUTS_ && opts::filenames.size() > 0){
			compileOutputs(data_set.get_markers(), data_set.get_families(), data_set.get_samples());
		}
	}
//	int fsize = families.size();
//	int ssize = samples.size();
//	int msize = markers.size();
//	for(int i = 0; i < msize; i++){
//		delete(markers[i]);
//	}
//	for(int i = 0; i < ssize; i++){
//		delete(samples[i]);
//	}
//	for(int i = 0; i < fsize; i++){
//		delete(families[i]);
//	}

}

vector<ORDER> optimize(ORDER* order){
	ORDER nothread;
	ORDER yesthread;

	vector<ORDER> results;

	opts::printLog("Optimizing processing order, thresholding and batch file rules will be maintained...\n");
	bool thresh = false;
	for(unsigned int o = 0; o < order->size(); o++){
		StepOptions* options = order->at(o).getOptions();
		if(options->doTransform() ||
				options->doThreshMarkersLow() ||
				options->doThreshMarkersHigh() ||
				options->doThreshSamplesLow() ||
				options->doThreshSamplesHigh() ||
				options->doThreshFamiliesHigh() ||
				options->doThreshFamiliesLow() ||
				options->doBpSpace() ||
				options->doLDchop() ||
				options->doRmMono() ||
				options->doRmHetOnly() ||
				options->doIncDisabledSamples() ||
				options->doZeroDisabled() ||
				thresh
		){
			nothread.push_back(order->at(o));
			thresh = true;
		}
		else{
			yesthread.push_back(order->at(o));
		}
	}

	results.push_back(yesthread);
	results.push_back(nothread);


	cout << "Threadable: \n";
	for(int i = 0; i < (int)yesthread.size(); i++){
		cout << yesthread[i].getName() << "\n";
	}
	cout << "Un-Threadable: \n";
	for(int i = 0; i < (int)nothread.size(); i++){
		cout << nothread[i].getName() << "\n";
	}

	return results;
}

void runStep(Step current_step, DataSet* data_set){
	opts::printLog("Working on " + current_step.getName() + "\n");
	//current_step.process(con, families, markers);
	try{
		StepOptions* step_options = current_step.getOptions();

		if(step_options->doTransform()){
			step_options->performTransforms(data_set);
		}
		current_step.process(data_set);
		current_step.PrintSummary();
		current_step.filter();
		current_step.FilterSummary();
		//current_step.close();
		if(step_options->doTransform()){
			step_options->undoTransforms(data_set);
		}
		//TODO: make sure that the close() statement being below doesn't break anything...
		current_step.close();
	}catch(MethodException & ex){
		opts::printLog(ex.what());
		exit(0);
	}catch(std::exception & ex){
		opts::printLog(ex.what());
		exit(0);
	}catch(...){
		opts::printLog("Unknown exception caught!");
		exit(0);
	}

}

/*
 * Function: getBatchArgs
 * Return: map<string, vector<string> > (string = step, vector = all step options, even repeating)
 * Parameters: Batch file
 * Description:
 * Places all steps and step options in a map to check against the web error checker
 */
map<string, vector<string> > getBatchArgs(string f){
	map<string, vector<string> > batchargs;

	ifstream infile;
	infile.open(f.c_str(), ios::in);
	if(!infile){
		cerr << "Error opening batch file: " << f << endl;
		exit(1);
	}
	string line = "";
	while(getline(infile, line)){
		if(line == ""){
			continue;
		}
		while(line.length() > 0 && (line.at(line.length() - 1) == ' ' || line.at(line.length() - 1) == '\t')){
			line.erase(line.length() - 1);
		}
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() == 0){
			continue;
		}
		else if(tokens[0] == "#"){
			continue;
		}
		else{
			string step = tokens[0];
			batchargs[step].push_back("NA");
			for(int i = 1; i < (int)tokens.size(); i++){
				batchargs[step].push_back(tokens[i]);
			}
		}
	}
	infile.close();

	return batchargs;

}

/*
 *Function: parseInput
 *Return: ORDER object (batch file contents)
 *Parameters: batch file
 *Description:
 *Parses batch file and sets up steps for processing.
 */
ORDER parseInput(string file){
	ifstream inFile;
	ORDER process;

	const char* temp = file.c_str();
	inFile.open(temp, ios::in);
	if(!inFile){
		cerr << "Error opening batch file: " << file << endl;
		exit(1);
	}

//	string step;
//	string thresh;
	int count = 0;
	string line;
	map<string, int> allsteps;
	while(getline(inFile, line)){
		if(line == ""){
			continue;
		}
		count++;
		bool overwrite = true;
		//remove white space at end of line (command to process)
		while(line.length() > 0 && (line.at(line.length() - 1) == ' ' || line.at(line.length() - 1) == '\t')){
			line.erase(line.length() - 1);
		}
		map<string, int>::iterator found = allsteps.find(line);
		//if the command is a duplicate
		if(found != allsteps.end()){
			overwrite = false;
		}
		allsteps[line]++;
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() == 0){
			continue;
		}
		else if(tokens[0][0] == '#'){
			continue;
		}
		else{
			string step = tokens[0];
			string thresh = "";
			for(unsigned int i = 1; i < tokens.size(); i++){
				thresh += tokens[i] + " ";
			}
			Step s = initializeSteps(step);
			s.setOrder(count);
			s.setThreshold(thresh);
			s.setOverwrite(overwrite);
			if(s.hasIncExc()){
				opts::_KEEP_EXC_SAMPLES_ = true;
			}
			process.push_back(s);
		}
	}
	inFile.close();


	return process;
}

/*
 *Function: initializeSteps
 *Parameter: batch step string
 *Return: Step object
 *Description:
 *Creates a Step object based on each step listed in the batch file
 *
 */
Step initializeSteps(string i){
		Step *newstep = NULL;
		Process* tempproc = NULL;
		map<string, StepValue>::iterator found = s_mapStepValues.find(i);
		if(found == s_mapStepValues.end()){
			opts::printLog("Step: " + i + " not recognized.  Exiting.\n");
			exit(1);
		}
		switch (s_mapStepValues[i]){
			case e_examplemodule:
				newstep = new Step("Example Module", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ExampleModule();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_cluster_missing:
				newstep = new Step("Cluster Missing", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessClusterMissing();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_marker_geno_eff:
				newstep = new Step("Marker Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMarkerGenoEff();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_sample_geno_eff:
				newstep = new Step("Individual Sample Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessSampleGenoEff();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;

			case e_family_geno_eff:
				newstep = new Step("Family Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PercentByFamily();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_filter_process:
				newstep = new Step("Filter Process", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessFilterProcess();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_mdr:
				newstep = new Step("MDR Process", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMDR();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_mdrpdt:
				newstep = new Step("MDRPDT Process", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMDRPDT();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_concordance:
				newstep = new Step("Concordance Check", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessConcordance();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_kinship:
				newstep = new Step("Kinship", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessKinship();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_fst:
				newstep = new Step("FST", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessFst();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_ibs:
				newstep = new Step("IBS", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessIBS();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_cmh:
				newstep = new Step("Cochran-Mantel-Haenszel test", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessCMH();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_logreg:
				newstep = new Step("Logistic Regression", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessLogReg();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_linreg:
				newstep = new Step("Linear Regression", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessLinearReg();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;

			case e_allele_freq:
				newstep = new Step("Allele Frequencies", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessAlleleFrequency();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;

			case e_mendelian_error:
				newstep = new Step("Mendelian Errors", "", true);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMendelianErrors();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;

			case e_hw:
				newstep = new Step("Hardy-Weinberg Calculations", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessHWEquilibrium();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;

			case e_gender_error:
				newstep = new Step("Gender Correctness (using X-chromosome markers)", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessGenderCheck();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;

//			case 8:
//				newstep = new Step("Genotype Quality Score Analysis (Platform Specific)", "", false);
//				if(tempproc != NULL){
//					delete(tempproc);
//				}
//				tempproc = new QualityScore();
//				if(opts::_DBOUTPUT_){
//					tempproc->setDBOUT();
//				}
//				if(opts::_MARKERLIST_){
//					tempproc->setMarkerList();
//				}
//				if(opts::_STRATIFY_){
//					tempproc->setStratify();
//				}
//				//if(opts::_QSFILE_.length() <= 0){
//				//	cerr << "Quality Score file not specified!  Exiting...\n";
//				//	exit(1);
//				//}
//				newstep->setProcess(tempproc);
//				break;

			case e_tdt:
				newstep = new Step("TDT", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessRunTDT();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;

			case e_superlink_output:
				newstep = new Step("Superlink Output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessSuperlinkOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_grr_output:
				newstep = new Step("GRR Output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessGRROutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
//			case 11:
//				newstep = new Step("PED Output", "", false);
//				if(tempproc != NULL){
//					delete(tempproc);
//				}
//				tempproc = new PEDOutput();
//				if(opts::_DBOUTPUT_){
//					tempproc->setDBOUT();
//				}
//				if(opts::_MARKERLIST_){
//					tempproc->setMarkerList();
//				}
//				if(opts::_STRATIFY_){
//					tempproc->setStratify();
//				}
//				newstep->setProcess(tempproc);
//				break;
//			case 12:
//				newstep = new Step("Quality Score Output", "", false);
//				if(tempproc != NULL){
//					delete(tempproc);
//				}
//				tempproc = new QSOutput();
//				if(opts::_DBOUTPUT_){
//					tempproc->setDBOUT();
//				}
//				if(opts::_MARKERLIST_){
//					tempproc->setMarkerList();
//				}
//				if(opts::_STRATIFY_){
//					tempproc->setStratify();
//				}
//				newstep->setProcess(tempproc);
//				break;
			case e_ped_output:
				newstep = new Step("PED file output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessPEDOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_bin_output:
				newstep = new Step("Binary (.bim, .fam, .bed) file output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessBINOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_impute_output:
				newstep = new Step("Create inpute files for IMPUTE", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessImputeOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_epistasis:
				newstep = new Step("Epistasis", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessEpistasis();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_interaction:
				newstep = new Step("Interaction", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessInteraction();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_eigenstrat_output:
				newstep = new Step("Eigenstrat file output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessEigenstratOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_tped_output:
				newstep = new Step("TPED file output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessTPEDOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_chisquare:
				newstep = new Step("Chisquare test", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessCaConChisq();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_structure_output:
				newstep = new Step("Create Structure input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessSTRUCTOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_phase_output:
				newstep = new Step("Create Phase input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessPHASEOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_beagle_output:
				newstep = new Step("Create BEAGLE input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessBEAGLEOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_lapis_output:
				newstep = new Step("Create LAPIS input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessLAPISOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_qtdt_output:
				newstep = new Step("Create QTDT input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessQTDTOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_mdr_output:
				newstep = new Step("Create MDR input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMDROutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_fbat_output:
				newstep = new Step("Create FBAT input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessFBATOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_homozygous:
				newstep = new Step("Homozygous spans", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessHomozygous();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_ld:
				newstep = new Step("LD calculations", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessLD();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_powermarker_output:
				newstep = new Step("Create PowerMarker input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessPowerMarkerOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_deletions:
				newstep = new Step("Deletion Detection", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessDeletions();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_pdt2_output:
				newstep = new Step("Create PDT2 input file.", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessPDT2Output();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_mitocheck:
				newstep = new Step("Mitochondrial Error Checking (Chrom 26)", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMitoCheck();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			default:
				opts::printLog("Step case " + i + " not found!! Exiting...\n");
				exit(1);
				break;
		}
	return *newstep;
}

/*Function: initializeSteps
 *Return: STEPS object
 *
 *Description:
 *Returns a list of valid steps.  Used when -S command line argument is specified
 *
 */
STEPS initializeSteps(){
	STEPS steps;
	map<string, StepValue>::iterator s_iter;
    for(s_iter = s_mapStepValues.begin(); s_iter != s_mapStepValues.end(); s_iter++){
		Step *newstep = NULL;
		Process* tempproc = NULL;
		switch (s_mapStepValues[s_iter->first]){
			case e_examplemodule:
				newstep = new Step("Example Module", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ExampleModule();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_cluster_missing:
				newstep = new Step("Cluster Missing", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessClusterMissing();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_kinship:
				newstep = new Step("Kinship", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessKinship();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_fst:
				newstep = new Step("FST", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessFst();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_marker_geno_eff:
				newstep = new Step("Marker Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMarkerGenoEff();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_filter_process:
				newstep = new Step("Filter Process", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessFilterProcess();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_mdr:
				newstep = new Step("MDR Process", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMDR();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_mdrpdt:
				newstep = new Step("MDRPDT Process", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMDRPDT();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_sample_geno_eff:
				newstep = new Step("Individual Sample Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessSampleGenoEff();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;

			case e_family_geno_eff:
				newstep = new Step("Family Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PercentByFamily();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_concordance:
				newstep = new Step("Concordance Check", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessConcordance();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_ibs:
				newstep = new Step("IBS", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessIBS();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_cmh:
				newstep = new Step("Cochran-Mantel-Haenszel test", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessCMH();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_logreg:
				newstep = new Step("Logistic Regression", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessLogReg();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_linreg:
				newstep = new Step("Linear Regression", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessLinearReg();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_allele_freq:
				newstep = new Step("Allele Frequencies", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessAlleleFrequency();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;

			case e_mendelian_error:
				newstep = new Step("Mendelian Errors", "", true);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMendelianErrors();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;

			case e_hw:
				newstep = new Step("Hardy-Weinberg Calculations", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessHWEquilibrium();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;

			case e_gender_error:
				newstep = new Step("Gender Correctness (using X-chromosome markers)", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessGenderCheck();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
/*
			case 8:
				newstep = new Step("Genotype Quality Score Analysis (Platform Specific)", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new QualityScore();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				//if(opts::_QSFILE_.length() <= 0){
				//	cerr << "Quality Score file not specified!  Exiting...\n";
				//	exit(1);
				//}
				newstep->setProcess(tempproc);
				steps[i] = *newstep;
				delete newstep;
				break;
*/
			case e_tdt:
				newstep = new Step("TDT", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessRunTDT();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;

			case e_superlink_output:
				newstep = new Step("Create Superlink input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessSuperlinkOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;

				break;
			case e_grr_output:
				newstep = new Step("Create GRR input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessGRROutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;

				break;
/*			case 11:
				newstep = new Step("PED Output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PEDOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[i] = *newstep;
				delete newstep;
				break;
			case 12:
				newstep = new Step("Quality Score Output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new QSOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[i] = *newstep;
				delete newstep;
				break;
*/
			case e_ped_output:
				newstep = new Step("Create PED file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessPEDOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_bin_output:
				newstep = new Step("Binary (.bed, .fam, .bim) file output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessBINOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_impute_output:
				newstep = new Step("Create input files for IMPUTE", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessImputeOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_epistasis:
				newstep = new Step("Epistasis", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessEpistasis();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_interaction:
				newstep = new Step("Interaction", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessInteraction();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_eigenstrat_output:
				newstep = new Step("Create Eigenstrat file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessEigenstratOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_tped_output:
				newstep = new Step("Create TPED file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessTPEDOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_chisquare:
				newstep = new Step("Chisquare test", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessCaConChisq();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_structure_output:
				newstep = new Step("Create Structure input files", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessSTRUCTOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
            case e_phase_output:
                newstep = new Step("Create Phase input file", "", false);
                if(tempproc != NULL){
                    delete(tempproc);
                }
                tempproc = new ProcessPHASEOutput();
                if(opts::_DBOUTPUT_){
                    tempproc->setDBOUT();
                }
                if(opts::_MARKERLIST_){
                    tempproc->setMarkerList();
                }
                if(opts::_STRATIFY_){
                   tempproc->setStratify();
                }
                newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
               break;
			case e_beagle_output:
				newstep = new Step("Create BEAGLE input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessBEAGLEOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;
			case e_lapis_output:
				newstep = new Step("Create LAPIS input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessLAPISOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;
			case e_qtdt_output:
				newstep = new Step("Create QTDT input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessQTDTOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;
			case e_mdr_output:
				newstep = new Step("Create MDR input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMDROutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;
			case e_fbat_output:
				newstep = new Step("Create FBAT input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessFBATOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;
			case e_homozygous:
				newstep = new Step("Homozygous spans", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessHomozygous();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;
			case e_ld:
				newstep = new Step("LD calculations", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessLD();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;
			case e_powermarker_output:
				newstep = new Step("Create PowerMarker input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessPowerMarkerOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;
			case e_deletions:
				newstep = new Step("Deletion Detection", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessDeletions();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;
			case e_pdt2_output:
				newstep = new Step("Create PDT2 input file.", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessPDT2Output();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;
			case e_mitocheck:
				newstep = new Step("Mitochondrial Error Checking (Chrom 26)", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ProcessMitoCheck();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;
			default:
				cerr << "Step case " << s_iter->first<< " not found!! Exiting..." << endl;
				exit(1);
		}
	}
	return steps;
}

/*
 *Function: flipStrand
 *Description:
 *Flips the strand of the specified markers
 *
 */
void flipStrand(vector<Marker*>* markers){
	if(opts::_FLIPFILE_.length() > 0){
		map<string, int> exclude;
		opts::printLog("Flipping markers found in: " + opts::_FLIPFILE_ + "\n");
		ifstream einput;
		einput.open(opts::_FLIPFILE_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker strand flip file: " + opts::_FLIPFILE_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";

		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Strand flip markers column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			exclude[tokens[0]] = 1;
			//exclude.push_back(probe);
		}
		einput.close();

		int msize = markers->size();

		for(int m = 0; m < msize; m++){
			Marker* mark = (*markers)[m];
			if(!mark->isMicroSat()){
				map<string, int>::iterator found = exclude.find(mark->getProbeID());
    	        if(found != exclude.end()){
					string a1 = mark->getAllele1();
					string a2 = mark->getAllele2();

					if(a1 == "A"){
						mark->resetAllele1("T");
					}
					else if(a1 == "C"){
						mark->resetAllele1("G");
					}
					else if(a1 == "G"){
						mark->resetAllele1("C");
					}
					else if(a1 == "T"){
						mark->resetAllele1("A");
					}
					else if(a1 == "1"){
						mark->resetAllele1("4");
					}
					else if(a1 == "2"){
						mark->resetAllele1("3");
					}
					else if(a1 == "3"){
						mark->resetAllele1("2");
					}
					else if(a1 == "4"){
						mark->resetAllele1("1");
					}

					if(a2 == "A"){
						mark->resetAllele2("T");
					}
					else if(a2 == "C"){
						mark->resetAllele2("G");
					}
					else if(a2 == "G"){
						mark->resetAllele2("C");
					}
					else if(a2 == "T"){
						mark->resetAllele2("A");
					}
					else if(a2 == "1"){
						mark->resetAllele2("4");
					}
					else if(a2 == "2"){
						mark->resetAllele2("3");
					}
					else if(a2 == "3"){
						mark->resetAllele2("2");
					}
					else if(a2 == "4"){
						mark->resetAllele2("1");
					}
				}
			}
		}
	}

}


/*
 *Function: printFamilies
 *Description:
 *outputs the structure of all families
 *
 */
void printFamilies(vector<Family*>* families) {
	int fsize = families->size();

	string fname = opts::_OUTPREFIX_ + "family_structure.txt";
	ofstream fout (fname.c_str());
	if(!fout.is_open()) {
		opts::printLog("Unable to open " + fname + "\n");
		exit(1);
	}

	fname = opts::_OUTPREFIX_ + "family_data.txt";
	ofstream fdout (fname.c_str());
	if(!fdout.is_open()) {
		opts::printLog("Unable to open " + fname + "\n");
		exit(1);
	}

	for(int f = 0; f < fsize; f++) {
		Family* fam = (*families)[f];
		if(fam->isEnabled()) {
			fout << "Family: " << fam->getFamID() << endl;
			fout << "----------------------------------------------------------" << endl;
			vector<Sample*>* founders = fam->getFounders();
			int size = founders->size();
			if(size> 0) {
				for(int fs = 0; fs < size; fs++) {
					Sample* founder = (*founders)[fs];
					if(founder->isEnabled() || (founder->isExcluded() && opts::_KEEP_EXC_SAMPLES_)) {
						fout << "Founder: " << founder->getInd() << endl;
						fout << descendTree(founder, 0) << endl;
					}
				}
			}
			else {
				fout << "No Founders...\n";
				vector<Sample*>* samps = fam->getSamples();
				int ssize = samps->size();
				for(int ss = 0; ss < ssize; ss++) {
					Sample* samp = (*samps)[ss];
					if(samp->isEnabled() || (samp->isExcluded() && opts::_KEEP_EXC_SAMPLES_)) {
						fout << descendTree(samp, 0) << endl;
					}
				}
			}
			fout << endl << endl;
		}
	}

	fdout << "FamID\tMultigenerational?\tMutigen_with_affecteds\tTotal_Generations\tTotal_Num_Affected\tGeneration1-N_Affected Counts\n";
	//output summary levels & affecteds
	for(int f = 0; f < fsize; f++){
		Family* fam = (*families)[f];
		if(fam->isEnabled()){
			map<int, vector<Sample*> > levels;
			fdout << fam->getFamID();
			vector<Sample*>* founders = fam->getFounders();
			int size = founders->size();
			if(size > 0){
				//vector<int> levels(size);
				for(int fs = 0; fs < size; fs++){
					Sample* founder = (*founders)[fs];
					if(founder->isEnabled() || (founder->isExcluded() && opts::_KEEP_EXC_SAMPLES_)){
						//levels[fs] = descendTree3(founder, 0);
						map<int, vector<Sample*> > levelstemp = descendTree3(founder, 1);
						//if(founder->getPheno() == 2){// && !founder->isFlagged()){
						map<int, vector<Sample*> >::iterator titer;
						for(titer = levelstemp.begin(); titer != levelstemp.end(); titer++){
							vector<Sample*> mysamps = titer->second;
							for(unsigned int ms = 0; ms < mysamps.size(); ms++){
								levels[titer->first].push_back(mysamps[ms]);
							}
						}
						levels[1].push_back(founder);
//							founder->setFlag(true);
						//}
					}
				}
				//sort(levels.begin(), levels.end());
				//fdout << "\t" << levels[levels.size() - 1];
			}
			else{
				vector<Sample*>* samps = fam->getSamples();
				int ssize = samps->size();
				//vector<int> levels(ssize);
				//map<int, int> levels;
				for(int ss = 0; ss < ssize; ss++){
					Sample* samp = (*samps)[ss];
					if(samp->isEnabled() || (samp->isExcluded() && opts::_KEEP_EXC_SAMPLES_)){
						//levels[ss] = descendTree3(samp, 0);
						map<int, vector<Sample*> > levelstemp = descendTree3(samp, 1);
						//if(samp->getPheno() == 2){// && !samp->isFlagged()){
						map<int, vector<Sample*> >::iterator titer;
						for(titer = levelstemp.begin(); titer != levelstemp.end(); titer++){
							vector<Sample*> mysamps = titer->second;
							for(unsigned int ms = 0; ms < mysamps.size(); ms++){
								levels[titer->first].push_back(mysamps[ms]);
							}
						}
						levels[1].push_back(samp);
//							samp->setFlag(true);
						//}
					}
				}

				//sort(levels.begin(), levels.end());
				//fdout << "\t" << levels[levels.size() - 1];
			}
//			vector<Sample*>* fsamps = fam->getSamples();
			int affected = 0;
//			for(int s = 0; s < fsamps->size(); s++){
//				if((*fsamps)[s]->getPheno() == 2 && (*fsamps)[s]->isFlagged()){
					//affected++;
//					(*fsamps)[s]->setFlag(false);
//				}
//			}
			map<int, vector<Sample*> >::iterator iter = levels.end();
			if(levels.size() > 1){
				fdout << "\tY";
			}
			else{
				fdout << "\tN";
			}
			for(iter = levels.begin(); iter != levels.end(); iter++){
				vector<Sample*> mysamps = iter->second;
				map<Sample*, bool> newsamps;
				for(unsigned int ms = 0; ms < mysamps.size(); ms++){
					newsamps[mysamps[ms]] = mysamps[ms]->getAffected();
				}
				map<Sample*, bool>::iterator niter;
				for(niter = newsamps.begin(); niter != newsamps.end(); niter++){

					if(niter->second){
						bool useme = true;
						map<int, vector<Sample*> >::iterator tempiter;
						for(tempiter = levels.begin(); tempiter != levels.end(); tempiter++){
							if(tempiter->first > iter->first){
								vector<Sample*> findsamples = tempiter->second;
								vector<Sample*>::iterator findme = find(findsamples.begin(), findsamples.end(), niter->first);
								if(findme != findsamples.end()){
									useme = false;
									break;
								}
							}
						}
						if(useme)
							affected++;
					}
				}
			}
			string affcounts = "";
			int numlevelsofaff = 0;
			for(iter = levels.begin(); iter != levels.end(); iter++){
				int laffected = 0;
				vector<Sample*> mysamps = iter->second;
				map<Sample*, bool> newsamps;
				for(unsigned int ms = 0; ms < mysamps.size(); ms++){
					newsamps[mysamps[ms]] = mysamps[ms]->getAffected();
				}
				bool incme = false;
				map<Sample*, bool>::iterator niter;
				for(niter = newsamps.begin(); niter != newsamps.end(); niter++){
					if(niter->second){
						bool useme = true;
						map<int, vector<Sample*> >::iterator tempiter;
						for(tempiter = levels.begin(); tempiter != levels.end(); tempiter++){
							if(tempiter->first > iter->first){
								vector<Sample*> findsamples = tempiter->second;
								vector<Sample*>::iterator findme = find(findsamples.begin(), findsamples.end(), niter->first);
								if(findme != findsamples.end()){
									useme = false;
									break;
								}
							}
						}
						if(useme){
							laffected++;
							if(!incme){
								numlevelsofaff++;
								incme = true;
							}
						}
					}
				}
				affcounts += "\t" + getString<int>(laffected);
			}
			if(numlevelsofaff > 1){
				fdout << "\tY";
			}
			else{
				fdout << "\tN";
			}
			fdout << "\t" << levels.size();
			fdout << "\t" << affected;
			fdout << affcounts;
			fdout << endl;
		}
	}

	fout.close();
	fdout.close();

}

/*
 *Function: descendTree
 *Description:
 *Recursively moves through pedigree until the child leaves are found to find structure
 */
string descendTree(Sample* sample, int level){
	string infospace = "     ";
	for(int i = 0; i < (level * 2); i++){
		infospace += "     ";
	}
	vector<Sample*>* children = sample->getChildren();
	if(children->size() == 0){
		string v = infospace + "->Ind: " + sample->getInd() + "\n" + infospace + "->Mom: " + sample->getMomID() + "\n" + infospace + "->Dad: " + sample->getDadID() + "\n";
		if(sample->getSib() != NULL){
			v += infospace + "->Sibling: " + sample->getSib()->getInd() + "\n";
		}
	   	v += infospace + "->Children: 0\n";
		return v;
	}
	int csize = children->size();

	string val = "";
	string v = infospace + "->Ind: " + sample->getInd() + "\n" + infospace + "->Mom: " + sample->getMomID() + "\n" + infospace + "->Dad: " + sample->getDadID() + "\n";
	if(sample->getSib() != NULL){
		v += infospace + "->Sibling: " + sample->getSib()->getInd() + "\n";
	}
	v += infospace + "->Children: " + getString<int>(sample->getChildren()->size()) + "\n";
	val += v;
	for(int c = 0; c < csize; c++){
		string temp = descendTree((*children)[c], (level + 1));

		val += temp;
	}
	return val;
}


/*
 *Function: descendTree3
 *Description: returns # of levels
 *Recursively moves through pedigree until the child leaves are found to find structure
 */
map<int, vector<Sample*> > descendTree3(Sample* sample, int level){
	map<int, vector<Sample*> > values;

	vector<Sample*>* children = sample->getChildren();
	if(children->size() == 0){
		return values;
		//return level;
	}
	int csize = children->size();
	vector<int> levels(csize);
	for(int c = 0; c < csize; c++){
		map<int, vector<Sample*> > tempvalues = descendTree3((*children)[c], (level + 1));
		map<int, vector<Sample*> >::iterator iter;
		for(iter = tempvalues.begin(); iter != tempvalues.end(); iter++){
			vector<Sample*> samps = iter->second;
			for(unsigned int is = 0; is < samps.size(); is++){
				values[iter->first].push_back(samps[is]);
			}
		}
//		if((*children)[c]->getPheno() == 2){// && !(*children)[c]->isFlagged()){
			values[(level + 1)].push_back((*children)[c]);
			//(*children)[c]->setFlag(true);
//		}
	}
//	sort(levels.begin(), levels.end());
//	level = levels[levels.size() - 1];
	return values;
}





/*
 *Function: compileOutputs
 *Description:
 *Compiles all QC & analysis outputs into single marker, sample, and family files
 *Performed when all steps are complete.
 */
void compileOutputs(vector<Marker*>* markers, vector<Family*>* families, vector<Sample*>* samples){
	opts::printLog("Compiling QC and Analysis output files...\n");


	map<string, vector<string> > filenames = opts::getFilenames();
	if(filenames["Marker"].size() > 0){
		map<Marker*, vector<string> > marker_output;
		vector<string> all_columns;
		map<string, Marker*> snpmap;
		for(unsigned int i = 0; i < markers->size(); i++){
			Marker* mark = (*markers)[i];
			string chr = getString<int>(mark->getChrom());
			string bp = getString<int>(mark->getBPLOC());
			snpmap[chr + "#" + bp] = mark;
		}

		cout << "Markers initialized...\n";

		for(unsigned int i = 0; i < filenames["Marker"].size(); i++){
			string file = filenames["Marker"][i];
			string step = opts::filesteps[file];

			for(unsigned int j = 0; j < opts::fileheaders[file].size(); j++){
				cout << "Pushing " << file << " : " << j << endl;
				all_columns.push_back("(" + file + ")" + opts::fileheaders[file][j]);
			}
		}
		string outfile = opts::_OUTPREFIX_ + "marker_summary.txt";
		opts::printLog("Working on compiling SNP information...[" + outfile + "]\n");
		for(unsigned int i = 0; i < filenames["Marker"].size(); i++){
			string filename = filenames["Marker"][i];
			string step = opts::filesteps[filename];
			vector<string> filecols = opts::fileheaders[filename];

			opts::printLog("\tParsing " + filename + "\n");


			ifstream input;
			input.open(filename.c_str(), ios::in);

			if(!input){
				cerr << "Error opening file: " << filename << endl;
				exit(1);
			}

			string header = "";
		   	getline(input, header);

			vector<string> columns;
			General::Tokenize(header, columns, "\t");
			int chrloc = -1;
			int bploc = -1;
			for(unsigned int l = 0; l < columns.size(); l++){
				if(columns[l] == "Chrom"){
					chrloc = l;
				}
				if(columns[l] == "bploc"){
					bploc = l;
				}
				if(chrloc > -1 && bploc > -1){
					break;
				}
			}
			string line;
			while(getline(input, line)){
				int count = 0;
				vector<string> elems;
				string token;
				istringstream isstream (line);
				while(getline(isstream, token, '\t')){
					elems.push_back(token);
				}

				//General::Tokenize(line, elems, "\t");
				if(elems.size() == 0){
					continue;
				}
				if(elems.size() != columns.size()){
					cout << "Error on line: " << count << endl;
					exit(1);
				}
				int chrom = atoi(elems[chrloc].c_str());
				int bp = atoi(elems[bploc].c_str());
				string key = elems[chrloc] + "#" + elems[bploc];

				Marker* itermark = snpmap[key];//(*markers)[count];
///				while(itermark->getChrom() != chrom && itermark->getBPLOC() != bploc){
///					if(count + 1 >= markers->size()){
///						itermark = NULL;
///						break;
///					}
///					itermark = (*markers)[++count];
///				}
				//vector<Marker*>::iterator itermark = find_if(markers->begin(), markers->end(), FindMarkerByChromBploc(chrom, bp));
				if(itermark == NULL){
					cout << "Cannot find marker with chrom = " << chrom << " and bploc = " << bp << endl;
					exit(1);
				}
				Marker* mark = itermark;
				map<Marker*, vector<string> >::iterator data = marker_output.find(mark);
				if(data == marker_output.end()){
					vector<string> mydata;
					mydata.resize(all_columns.size(), "N/A");
					marker_output[mark] = mydata;
				}
				for(unsigned int j = 0; j < filecols.size(); j++){
					vector<string>::iterator realcol = find(columns.begin(), columns.end(), filecols[j]);
					vector<string>::iterator allloc = find(all_columns.begin(), all_columns.end(), "(" +filename+")"+filecols[j]);
					if(realcol != columns.end()){
						int myloc = realcol - columns.begin();
						int myallloc = allloc - all_columns.begin();
						marker_output[mark][myallloc] = elems[myloc];
					}
				}
				count++;
			}
			input.close();

		}
		ofstream out(outfile.c_str());
		out << "Chrom\trsID\tProbeID\tbploc";
		for(unsigned int i = 0; i < all_columns.size(); i++){
			out << "\t" << all_columns[i];
		}
		out << "\n";
		map<Marker*, vector<string> >::iterator iter;
		for(iter = marker_output.begin(); iter != marker_output.end(); iter++){
			Marker* mark = iter->first;
			vector<string> data = iter->second;
			out << mark->getChrom() << "\t"
				<< mark->getRSID() << "\t"
				<< mark->getProbeID() << "\t"
				<< mark->getBPLOC();
			for(unsigned int i = 0; i < data.size(); i++){
				out << "\t" << data[i];
			}
			out << "\n";
		}
		out.close();
	}
	if(filenames["Family"].size() > 0){
		map<Family*, vector<string> > family_output;
		vector<string> all_columns;
		for(unsigned int i = 0; i < filenames["Family"].size(); i++){
			string file = filenames["Family"][i];
			string step = opts::filesteps[file];

			for(unsigned int j = 0; j < opts::fileheaders[file].size(); j++){
				all_columns.push_back("(" + file + ")" + opts::fileheaders[file][j]);
			}
		}
		string outfile = opts::_OUTPREFIX_ + "family_summary.txt";
		opts::printLog("Working on compiling Family information...[" + outfile + "]\n");
		for(unsigned int i = 0; i < filenames["Family"].size(); i++){
			string filename = filenames["Family"][i];
			string step = opts::filesteps[filename];
			vector<string> filecols = opts::fileheaders[filename];

			opts::printLog("\tParsing " + filename + "\n");


			ifstream input;
			input.open(filename.c_str(), ios::in);

			if(!input){
				cerr << "Error opening file: " << filename << endl;
				exit(1);
			}

			string header = "";
		   	getline(input, header);

			vector<string> columns;
			General::Tokenize(header, columns, "\t");
			int famloc = -1;
			for(unsigned int l = 0; l < columns.size(); l++){
				if(columns[l] == "FamID"){
					famloc = l;
				}
				if(famloc > -1){
					break;
				}
			}
			if(famloc < 0){
				cout << "FamID column not found!\n";
				exit(1);
			}
			string line;
			int count = 1;
			while(getline(input, line)){
				count++;
				vector<string> elems;
				string token;
				istringstream isstream (line);
				while(getline(isstream, token, '\t')){
					elems.push_back(token);
				}
				//General::Tokenize(line, elems, "\t");
				if(elems.size() == 0){
					continue;
				}
				string famid = elems[famloc];
				vector<Family*>::iterator iterfam = find_if(families->begin(), families->end(), FindFamily(famid));
				if(iterfam == families->end()){
					cout << "Cannot find family with famid = " << famid << endl;
					exit(1);
				}
				Family* fam = *iterfam;
				map<Family*, vector<string> >::iterator data = family_output.find(fam);
				if(data == family_output.end()){
					vector<string> mydata;
					mydata.resize(all_columns.size(), "N/A");
					family_output[fam] = mydata;
				}
				for(unsigned int j = 0; j < filecols.size(); j++){
					vector<string>::iterator realcol = find(columns.begin(), columns.end(), filecols[j]);
					vector<string>::iterator allloc = find(all_columns.begin(), all_columns.end(), "(" +filename+")"+filecols[j]);
					if(realcol != columns.end()){
						int myloc = realcol - columns.begin();
						int myallloc = allloc - all_columns.begin();
						family_output[fam][myallloc] = elems[myloc];
					}
				}
			}
			input.close();

		}
		ofstream out(outfile.c_str());
		out << "FamID\tNumInds";
		for(unsigned int i = 0; i < all_columns.size(); i++){
			out << "\t" << all_columns[i];
		}
		out << "\n";
		map<Family*, vector<string> >::iterator iter;
		for(iter = family_output.begin(); iter != family_output.end(); iter++){
			Family* fam = iter->first;
			vector<string> data = iter->second;
			out << fam->getFamID() << "\t"
				<< fam->getSamples()->size();
			for(unsigned int i = 0; i < data.size(); i++){
				out << "\t" << data[i];
			}
			out << "\n";
		}
		out.close();
	}
	if(filenames["Sample"].size() > 0){
	map<Sample*, vector<string> > sample_output;
		vector<string> all_columns;
		for(unsigned int i = 0; i < filenames["Sample"].size(); i++){
			string file = filenames["Sample"][i];
			string step = opts::filesteps[file];

			for(unsigned int j = 0; j < opts::fileheaders[file].size(); j++){
				all_columns.push_back("(" + file + ")" + opts::fileheaders[file][j]);
			}
		}
		string outfile = opts::_OUTPREFIX_ + "sample_summary.txt";
		opts::printLog("Working on compiling Sample information...[" + outfile + "]\n");
		for(unsigned int i = 0; i < filenames["Sample"].size(); i++){
			string filename = filenames["Sample"][i];
			string step = opts::filesteps[filename];
			vector<string> filecols = opts::fileheaders[filename];

			opts::printLog("\tParsing " + filename + "\n");


			ifstream input;
			input.open(filename.c_str(), ios::in);

			if(!input){
				cerr << "Error opening file: " << filename << endl;
				exit(1);
			}

			string header = "";
		   	getline(input, header);

			vector<string> columns;
			General::Tokenize(header, columns, "\t");
			int famloc = -1;
			int indloc = -1;
			for(unsigned int l = 0; l < columns.size(); l++){
				if(columns[l] == "FamID"){
					famloc = l;
				}
				if(columns[l] == "IndID"){
					indloc = l;
				}
				if(famloc > -1 && indloc > -1){
					break;
				}
			}
			string line;
			int count = 1;
			while(getline(input, line)){
				count++;
				vector<string> elems;
				string token;
				istringstream isstream (line);
				while(getline(isstream, token, '\t')){
					elems.push_back(token);
				}
				//General::Tokenize(line, elems, "\t");
				if(elems.size() == 0){
					continue;
				}
				string famid = elems[famloc];
				string indid = elems[indloc];
				vector<Sample*>::iterator itersamp = find_if(samples->begin(), samples->end(), FindSampleByFamAndID(famid, indid));
				if(itersamp == samples->end()){
					cout << "Cannot find sample with famid = " << famid << " and indid = " << indid << endl;
					exit(1);
				}
				Sample* samp = *itersamp;
				map<Sample*, vector<string> >::iterator data = sample_output.find(samp);
				if(data == sample_output.end()){
					vector<string> mydata;
					mydata.resize(all_columns.size(), "N/A");
					sample_output[samp] = mydata;
				}
				for(unsigned int j = 0; j < filecols.size(); j++){
					vector<string>::iterator realcol = find(columns.begin(), columns.end(), filecols[j]);
					vector<string>::iterator allloc = find(all_columns.begin(), all_columns.end(), "(" +filename+")"+filecols[j]);
					if(realcol != columns.end()){
						int myloc = realcol - columns.begin();
						int myallloc = allloc - all_columns.begin();
						sample_output[samp][myallloc] = elems[myloc];
					}
				}
			}
			input.close();

		}
		ofstream out(outfile.c_str());
		out << "FamID\tIndID\tSex\tAffection_Status";
		for(unsigned int i = 0; i < all_columns.size(); i++){
			out << "\t" << all_columns[i];
		}
		out << "\n";
		map<Sample*, vector<string> >::iterator iter;
		for(iter = sample_output.begin(); iter != sample_output.end(); iter++){
			Sample* samp = iter->first;
			vector<string> data = iter->second;
			out << samp->getFamID() << "\t"
				<< samp->getInd() << "\t";
			if(samp->getSex()){
				out << "M\t";
			}
			else{
				out << "F\t";
			}
			out << samp->getPheno();

			for(unsigned int i = 0; i < data.size(); i++){
				out << "\t" << data[i];
			}
			out << "\n";
		}
		out.close();
	}
	if(filenames["Batch"].size() > 0){
		map<string, vector<string> > batch_output;
		vector<string> all_columns;
		for(unsigned int i = 0; i < filenames["Batch"].size(); i++){
			string file = filenames["Batch"][i];
			string step = opts::filesteps[file];

			for(unsigned int j = 0; j < opts::fileheaders[file].size(); j++){
				all_columns.push_back("(" + file + ")" + opts::fileheaders[file][j]);
			}
		}
		string outfile = opts::_OUTPREFIX_ + "batch_summary.txt";
		opts::printLog("Working on compiling Batch information...[" + outfile + "]\n");
		for(unsigned int i = 0; i < filenames["Batch"].size(); i++){
			string filename = filenames["Batch"][i];
			string step = opts::filesteps[filename];
			vector<string> filecols = opts::fileheaders[filename];

			opts::printLog("\tParsing " + filename + "\n");


			ifstream input;
			input.open(filename.c_str(), ios::in);

			if(!input){
				cerr << "Error opening file: " << filename << endl;
				exit(1);
			}

			string header = "";
		   	getline(input, header);

			vector<string> columns;
			General::Tokenize(header, columns, "\t");
			int famloc = -1;
			for(unsigned int l = 0; l < columns.size(); l++){
				if(columns[l] == "Batch"){
					famloc = l;
				}
				if(famloc > -1){
					break;
				}
			}
			if(famloc < 0){
				cout << "Batch column not found!\n";
				exit(1);
			}
			string line;
			int count = 1;
			while(getline(input, line)){
				count++;
				vector<string> elems;
				string token;
				istringstream isstream (line);
				while(getline(isstream, token, '\t')){
					elems.push_back(token);
				}
				//General::Tokenize(line, elems, "\t");
				if(elems.size() == 0){
					continue;
				}
				string batchid = elems[famloc];
//				vector<>::iterator iterfam = find_if(families->begin(), families->end(), FindFamily(famid));
//				if(iterfam == families->end()){
//					cout << "Cannot find family with famid = " << famid << endl;
//					exit(1);
//				}
//				Family* fam = *iterfam;
				map<string, vector<string> >::iterator data = batch_output.find(batchid);
				if(data == batch_output.end()){
					vector<string> mydata;
					mydata.resize(all_columns.size(), "N/A");
					batch_output[batchid] = mydata;
				}
				for(unsigned int j = 0; j < filecols.size(); j++){
					vector<string>::iterator realcol = find(columns.begin(), columns.end(), filecols[j]);
					vector<string>::iterator allloc = find(all_columns.begin(), all_columns.end(), "(" +filename+")"+filecols[j]);
					if(realcol != columns.end()){
						int myloc = realcol - columns.begin();
						int myallloc = allloc - all_columns.begin();
						batch_output[batchid][myallloc] = elems[myloc];
					}
				}
			}
			input.close();

		}
		ofstream out(outfile.c_str());
		out << "Batch";
		for(unsigned int i = 0; i < all_columns.size(); i++){
			out << "\t" << all_columns[i];
		}
		out << "\n";
		map<string, vector<string> >::iterator iter;
		for(iter = batch_output.begin(); iter != batch_output.end(); iter++){
			string batchid = iter->first;
			vector<string> data = iter->second;
			out << batchid;
			for(unsigned int i = 0; i < data.size(); i++){
				out << "\t" << data[i];
			}
			out << "\n";
		}
		out.close();
	}

}





//////////////////////////////////////////////////////////////////
// Borrowed from Plink v0.99s                                   //
//                                                              //
//////////////////////////////////////////////////////////////////


#define  PORT_NUM                80
#define  IP_ADDR    "160.129.37.40"
#define  GET_STRING "GET /plato/files/version.txt HTTP/1.1\nHost: chgr.mc.vanderbilt.edu\nConnection: close\n\n"


void webcheck(vector<string> a, map<string, vector<string> > b)
{

#ifdef SKIP
	opts::printLog("Web-check not implemented on this system...\n");
  return;
#else

  opts::printLog("Web-based version check ( -noweb to skip )\n");

  vector<string> tokens = socketConnection(
					    IP_ADDR,
					    PORT_NUM,
					    GET_STRING);

  bool print = false;
  bool print2 = false;
  bool version_okay = true;
  vector<string>::iterator fiter;
  map<string, vector<string> >::iterator siter;

  for (unsigned int i=0; i<tokens.size(); i++)
    {

      if (tokens[i]=="END") break;

      if (tokens[i]=="END-MESSAGE")
	{
	  print2=false;
	  continue;
	}

      if (tokens[i]=="WARN")
	{
	  if ( i < tokens.size()-1 )
	    {
	      i++;
		  fiter = find(a.begin(), a.end(), tokens[i]);
	      if ( fiter != a.end())
		{
		  opts::printLog("\n*** ALERT ***\n*** A warning flag has been set for: "+tokens[i]+
			   "\n*** See http://chgr.mc.vanderbilt.edu/plato/\n");
		}
	    }
	  continue;
	}

	  if(tokens[i] == "WARNSTEP"){
	  	if(i < tokens.size() - 1){
			i++;
			siter = b.find(tokens[i]);
			if(siter != b.end()){
				opts::printLog("\n*** ALERT ***\n*** A warning flag has been set for STEP: " + tokens[i] +
						"\n*** See http://chgr.mc.vanderbilt.edu/plato/\n");
			}
		}
		continue;
	  }

	  if(tokens[i] == "WARNSTEPOPTION"){
		  if(i < tokens.size() - 1){
			  i++;
			  for(siter = b.begin(); siter != b.end(); siter++){
				  vector<string> second = siter->second;
				  fiter = find(second.begin(), second.end(), tokens[i]);
				  if(fiter != second.end()){
						opts::printLog("\n*** ALERT ***\n*** A warning flag has been set for STEP OPTION: " + tokens[i] +
							"\n*** See http://chgr.mc.vanderbilt.edu/plato/\n");
				  	break;
				  }
			  }
		  }
		  continue;
	  }

      if (tokens[i]=="FATAL")
	{
	  if ( i < tokens.size()-1 )
	    {
	      i++;
		  fiter = find(a.begin(), a.end(), tokens[i]);
	      if ( fiter != a.end()){
		opts::printLog("A serious warning flag has been set for: "+tokens[i]+
		    "\nPlato has been instructed to stop"+
 	            "\nPlease see http://chgr.mc.vanderbilt.edu/plato/\n");
		  	exit(1);
		  }
	    }
	  continue;
	}

	  if(tokens[i] == "FATALSTEP"){
	  	if(i < tokens.size() - 1){
			i++;
			siter = b.find(tokens[i]);
			if(siter != b.end()){
		opts::printLog("A serious warning flag has been set for STEP: "+tokens[i]+
		    "\nPlato has been instructed to stop"+
 	            "\nPlease see http://chgr.mc.vanderbilt.edu/plato/\n");
		  	exit(1);
			}
		}
		continue;
	  }

	  if(tokens[i] == "FATALSTEPOPTION"){
		  if(i < tokens.size() - 1){
			  i++;
			  for(siter = b.begin(); siter != b.end(); siter++){
				  vector<string> second = siter->second;
				  fiter = find(second.begin(), second.end(), tokens[i]);
				  if(fiter != second.end()){
		opts::printLog("A serious warning flag has been set for STEP OPTION: "+tokens[i]+
		    "\nPlato has been instructed to stop"+
 	            "\nPlease see http://chgr.mc.vanderbilt.edu/plato/\n");
		  	exit(1);
				  }
			  }
		  }
		  continue;
	  }

      if (tokens[i]=="MESSAGE-ALL")
	{
	  print2=true;
	  continue;
	}

      // Display any other messages
      // Either conditional on old version (print)
      // or a broadcast to all users (print2)

      if ( ( print && !version_okay) || print2 )
	{
	  if (tokens[i]=="\\n")
	    opts::printLog("\n");
	  else
	    opts::printLog(tokens[i]+" ");
	}

      // Check version code
      if (tokens[i]=="WASPVER")
	{
	  print=true;
	  if ( i < tokens.size() - 1)
	    {
		  float currver = atof(_WASPVER_.c_str());
		  float webver = atof(tokens[i+1].c_str());
	      if (currver >= webver)//tokens[i+1] >= _WASPVER_)
		opts::printLog(" OK, v"+_WASPVER_+" is current\n");
	      else
		{
		  opts::printLog("\n\n          *** UPDATE REQUIRED ***\n\n");
		  opts::printLog("\tThis version        : "+_WASPVER_+"\n");
		  opts::printLog("\tMost recent version : "+tokens[i+1]+"\n\n");
		  opts::printLog("Please upgrade your version of PLATO as soon as possible!\n");
		  opts::printLog("  (visit the above website for free download)\n\n");
		  version_okay=false;
		}

	      // Skip the version number
	      i++;
	    }
	}

    }

  // did we get the information we needed?
  if (!print) opts::printLog(" problem connecting to web\n");

  opts::printLog("\n");
#endif
}

