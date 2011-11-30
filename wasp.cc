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
#include <Globals.h>
#include <vector>
#include <General.h>//"General.h"
#include "wasp.h"


using namespace std;

static string _WASPVER_ = "0.83";
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
					   e_concordance,	   //concordance
					   e_logreg,			//logreg
					   e_mitocheck,         //mito-check
					   e_cmh				//CMH
					 };
enum cmdArgs{
	a_h,
	a_S,
	a_ped,
	a_map,
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
	a_flip,
	a_make_bin,
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
	a_trait_file

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
	s_mapStepValues["cmh"] = e_cmh;

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
	s_mapcmdArgs["-pedinfo"] = a_pedinfo;
	s_mapcmdArgs["-dog"] = a_dog;
	s_mapcmdArgs["-flip"] = a_flip;
	s_mapcmdArgs["-make-bin"] = a_make_bin;
	s_mapcmdArgs["-mapdesc"] = a_mapdesc;
	s_mapcmdArgs["-micro-sats"] = a_micro_sats;
	s_mapcmdArgs["-no-summary"] = a_no_summary;
	s_mapcmdArgs["-noweb"] = a_noweb;
	s_mapcmdArgs["-out"] = a_out;
	s_mapcmdArgs["-sampdesc"] = a_sampdesc;
	s_mapcmdArgs["-tped"] = a_tped;
	s_mapcmdArgs["-tfam"] = a_tfam;
	s_mapcmdArgs["-print-families"] = a_print_families;
	s_mapcmdArgs["-zero-geno-file"] = a_zero_geno_file;
	s_mapcmdArgs["-remaining"] = a_remaining;
	s_mapcmdArgs["-inccenters"] = a_inccenters;
	s_mapcmdArgs["-missing-geno"] = a_missing_geno;
	s_mapcmdArgs["-covar-missing"] = a_covar_missing;
	s_mapcmdArgs["-trait-missing"] = a_trait_missing;
	s_mapcmdArgs["-covar-file"] = a_covar_file;
	s_mapcmdArgs["-trait-file"] = a_trait_file;
	s_mapcmdArgs["-bp-space"] = a_bp_space;
	s_mapcmdArgs["-todigit"] = a_todigit;
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
				case a_noweb:
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
				case a_dog:
					opts::_DOG_ = true;
					opts::_CHRX_ = 39;
					opts::_CHRY_ = 40;
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
				case a_make_bin:
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
		string logoutput = opts::_OUTPREFIX_ + "wasp.log";
//		ofstream log(logoutput.c_str());
		opts::_LOG_.open(logoutput.c_str());
		if(!opts::_LOG_){
			cerr << "Error opening wasp.log.  Exiting!\n";
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

		if(opts::_PEDFILE_.length() == 0 && opts::_BINPREFIX_.length() == 0 && opts::_TPEDFILE_.length() == 0){
			opts::printLog("No pedfile specified and no binary input files specified.  Exiting!\n");
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
		//Begin batch file steps and reading in data
		startProcess(&proc_order, NULL, myrank, &cmd_filters);
	}catch(MethodException ex){
		cerr << ex.what() << endl;
		exit(0);
	}//catch(...){
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

/*depricated*/
vector<string> readMarkers(string f){
	ifstream input;
	input.open(f.c_str(), ios::in);

	if(!input){
		cerr << "Error opening marker file: " << f << endl;
		exit(1);
	}

	string rsid = "";
	vector<string> mrks;
	while(input >> rsid){
		mrks.push_back(rsid);
	}
	input.close();

	return mrks;
}

/*depricated*/
vector<string> readCenters(string f){
	ifstream input;
	input.open(f.c_str(), ios::in);
	if(!input){
		cerr << "Error opening center file: " << f << endl;
		exit(1);
	}

	string center = "";
	vector<string> cntrs;
	while(input >> center){
		cntrs.push_back(center);
	}

	input.close();

	return cntrs;
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
	cout << "wasp - Whole-Genome Associstion Study Pipeline (Version: " << _WASPVER_ << ")" << endl << endl;
	cout << "Center for Human Genetics Research" << endl;
	cout << "Vanderbilt University Medical Center" << endl << endl;
	cout << "!!!!  For most up-to-date options and details  !!!!" << endl;
	cout << "!!!!  visit http://chgr.mc.vanderbilt.edu/wasp !!!!" << endl << endl;
	cout << "usage: wasp <batchfile> [<option1> <option2>...]" << endl << endl;
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
	cout << "wasp - Whole-Genome Associstion Study Pipeline (Version: " << _WASPVER_ << ")" << endl << endl;
	cout << "usage: wasp <batchfile> [<option1> <option2>....]" << endl
		 << endl;
	cout << "For a list of valid steps to be inserted into the batch file:" << endl;
	cout << "\t\t> wasp -S" << endl << endl;
	cout << "For help: " << endl;
	cout << "\t\t> wasp -h" << endl << endl;
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

	DataSet data_set;
		//set up filters
//		InputFilter filters;
		if(opts::_EXCCOVS_.length() > 0){
			vector<string>* list = new vector<string>;
			readCovTraitFile(opts::_EXCCOVS_, list);
			filters->add_covariate_list(list);
			filters->add_covariate_filter(InputFilter::ExcludeCovariateFilter);
		}
		if(opts::_INCCOVS_.length() > 0){
			vector<string>* list = new vector<string>;
			readCovTraitFile(opts::_INCCOVS_, list);
			filters->add_covariate_list(list);
			filters->add_covariate_filter(InputFilter::IncludeCovariateFilter);
		}
		if(opts::_EXCTRAITS_.length() > 0){
			vector<string>* list = new vector<string>;
			readCovTraitFile(opts::_EXCTRAITS_, list);
			filters->add_covariate_list(list);
			filters->add_covariate_filter(InputFilter::ExcludeTraitFilter);
		}
		if(opts::_INCTRAITS_.length() > 0){
			vector<string>* list = new vector<string>;
			readCovTraitFile(opts::_INCTRAITS_, list);
			filters->add_covariate_list(list);
			filters->add_covariate_filter(InputFilter::IncludeTraitFilter);
		}
		if(opts::_EXCMARKERS_.length() > 0){
			vector<Marker*>* list = new vector<Marker*>;
			readLocusFile(opts::_EXCMARKERS_, list);
			filters->add_locus_list(list);
			filters->add_locus_filter(InputFilter::ExcludeLocusFilter);
		}
		if(opts::_INCMARKERS_.length() > 0){
			vector<Marker*>* list = new vector<Marker*>;
			readLocusFile(opts::_INCMARKERS_, list);
			filters->add_locus_list(list);
			filters->add_locus_filter(InputFilter::IncludeLocusFilter);
		}
		if(opts::_EXCSAMPS_.length() > 0){
			vector<Sample*>* list = new vector<Sample*>;
			readSampleFile(opts::_EXCSAMPS_, list);
			filters->add_sample_list(list);
			filters->add_sample_filter(InputFilter::ExcludeSampleFilter);
		}
		if(opts::_INCSAMPLES_.length() > 0){
			vector<Sample*>* list = new vector<Sample*>;
			readSampleFile(opts::_INCSAMPLES_, list);
			filters->add_sample_list(list);
			filters->add_sample_filter(InputFilter::IncludeSampleFilter);
		}
		if(opts::_EXCFAMILIES_.length() > 0){
			vector<Family*>* list = new vector<Family*>;
			readFamilyFile(opts::_EXCFAMILIES_, list);
			filters->add_family_list(list);
			filters->add_family_filter(InputFilter::ExcludeFamilyFilter);
		}
		if(opts::_INCFAMILIES_.length() > 0){
			vector<Family*>* list = new vector<Family*>;
			readFamilyFile(opts::_INCFAMILIES_, list);
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
		options.setPedFile(opts::_PEDFILE_);
		options.setBinInput(opts::_BINPREFIX_);
		options.setTFamFile(opts::_FAMFILE_);
		options.setTPedFile(opts::_TPEDFILE_);
		options.setCovarMissing(opts::_COVAR_MISSING_);
		options.setTraitMissing(opts::_TRAIT_MISSING_);

		//read zerogenofile
		if(opts::_ZEROGENOFILE_.length() > 0){
			readZeroGenoFile();
		}
		//read Binary input
		if(opts::_BINPREFIX_.length() > 0 && !opts::_MAKEBIN_){
			if(opts::_PEDINFO_.length() > 0){
				opts::printLog("Reading Pedigree information file: " + opts::_PEDINFO_ + "\n");
				readPedInfo();
			}
			opts::printLog("Reading data using Binary inputs: prefix: " + opts::_BINPREFIX_ + "\n");
			if(opts::_MICROSATS_){
				opts::printLog("You specified microsatellite markers exist.  This option is not compatible with binary input files.  Please use the standard PED file format for your input when using microsatellite markers.  Exiting...\n");
				exit(1);
			}
			readBinM(&data_set, options, filters);//(data_set.get_samples(), data_set.get_families(), data_set.get_markers(), data_set.get_marker_map());
			//assignLinks(data_set.get_families());
			//reorderAlleles(data_set.get_samples(), data_set.get_markers());
			if(opts::_FLIPSTRAND_){
				flipStrand(data_set.get_markers());
			}
		}
		//read ped & map file
		else if(opts::_PEDFILE_.length() > 0 && opts::_MAPFILE_.length() > 0 && opts::_TPEDFILE_.length() == 0){
			if(opts::_MAPFILE_.length() > 0){
				opts::printLog("Reading MAP file: " + opts::_MAPFILE_ + "\n");
				//readMap(data_set.get_markers(), data_set.get_marker_map());
				readMapM(&data_set, options, filters);
			}
			if(opts::_PEDFILE_.length() > 0){
				if(opts::_PEDINFO_.length() > 0){
					opts::printLog("Reading Pedigree information file: " + opts::_PEDINFO_ + "\n");
					readPedInfo();
				}
				opts::printLog("Reading PED file: " + opts::_PEDFILE_ + "\n");
				//readPed(data_set.get_samples(), data_set.get_families(), data_set.get_markers(), data_set.get_marker_map());
				//assignLinks(data_set.get_families());
				//reorderAlleles(data_set.get_samples(), data_set.get_markers());
				readPedM_3vec_set(&data_set, options, filters);
				if(opts::_FLIPSTRAND_){
					flipStrand(data_set.get_markers());
				}
			}
		}
		else if(opts::_PEDFILE_.length() == 0 && opts::_TPEDFILE_.length() > 0 && opts::_FAMFILE_.length() > 0){
			if(opts::_FAMFILE_.length() > 0){
				if(opts::_PEDINFO_.length() > 0){
					opts::printLog("Reading Pedigree information file: " + opts::_PEDINFO_ + "\n");
					readPedInfo();
				}
				opts::printLog("Reading Pedigree file: " + opts::_FAMFILE_ + "\n");
				readTFamM(&data_set, options, filters);
				//assignLinks(data_set.get_families());
			}
			if(opts::_TPEDFILE_.length() > 0){
				opts::printLog("Reading TPED file: " + opts::_TPEDFILE_ + "\n");
				readTPedM(&data_set, options, filters);//data_set.get_markers(), data_set.get_samples(), data_set.get_marker_map());
				//reorderAlleles(data_set.get_samples(), data_set.get_markers());
				if(opts::_FLIPSTRAND_){
					flipStrand(data_set.get_markers());
				}
			}
		}
		else if(!opts::_DBINPUT_){
			cerr << "No input method specified!" << endl;
			exit(1);
		}
		if(opts::_COVFILE_.length() > 0){
			readCovariateFile(opts::_COVFILE_, &data_set, options, filters);
		}
		if(opts::_TRAITFILE_.length() > 0){
			readTraitFile(opts::_TRAITFILE_, &data_set, options, filters);
		}

		//int goodsamps = 0;
		//int goodfams = 0;
		//int goodmarkers = 0;
		if(opts::_TODIGIT_){
			opts::printLog("Mapping Families/Samples to numeric format.");
			remapFamsToDigit(data_set.get_families());
			printFamsToDigit(data_set.get_families(), "global", options);
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

		if(opts::zerogenoinfo.size() > 0){
			zeroSingleGenos(data_set.get_markers(), data_set.get_samples());
		}

		if(opts::_MAKEBIN_){
			writeBit(data_set.get_samples(), data_set.get_families(), data_set.get_markers(), data_set.get_marker_map());
			exit(1);
		}
		if(opts::_PRINTFAMS_){
			printFamilies(data_set.get_families());
			exit(1);
		}
	//}
	//set_me_up->summary();

	//Start processing steps by slaves.  MASTER gathers all data from slaves and does summary work
	for(o_iter = order->begin(); o_iter != order->end(); o_iter++){
		Step current_step = (Step)(*o_iter);//(*steps)[(*o_iter)];
		opts::printLog("Working on " + current_step.getName() + "\n");
		//current_step.process(con, families, markers);
		try{
			current_step.process(&data_set);
			current_step.PrintSummary();
			current_step.filter();
			current_step.FilterSummary();
			current_step.close();
		}catch(MethodException ex){
			opts::printLog(ex.what());
			exit(0);
		}catch(std::exception ex){
			opts::printLog(ex.what());
			exit(0);
		}catch(...){
			opts::printLog("Unknown exception caught!");
			exit(0);
		}
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

	//cout << myrank << " leaving process function" << endl;
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
			for(unsigned int i = 1; i < tokens.size(); i++){
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
		while(line.length() > 0 && (line.at(line.length() - 1) == ' ' || line.at(line.length() - 1) == '\t')){
			line.erase(line.length() - 1);
		}
		map<string, int>::iterator found = allsteps.find(line);
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
				newstep = new Step("Mitochondrial Error Checking (Chrom 25)", "", false);
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
				newstep = new Step("Mitochondrial Error Checking (Chrom 25)", "", false);
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
 *Function: reorderAlleles
 *Return: none
 *Parameters: vector of samples, vector of markers
 *Description:
 *Calculates the raw allele counts and remaps the sample bool vector 1 to be the
 *minor allele
 *
 */
/*void reorderAlleles(vector<Sample*>* samples, vector<Marker*>* markers){
	int msize = markers->size();
	int ssize = samples->size();

	for(int m = 0; m < msize; m++){
		int mloc = (*markers)[m]->getLoc();
		if((*markers)[m]->getNumAlleles() == 0){
			(*markers)[m]->addAllele("0");
			(*markers)[m]->addAllele("0");
		}
		if((*markers)[m]->getNumAlleles() == 1){
			(*markers)[m]->addAllele("0");
		}
		if((*markers)[m]->getNumAlleles() <= 2){
			int a1 = 0;
			int a2 = 0;
			for(int s = 0; s < ssize; s++){
				if((*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc)){
					a2+=2;
				}
				else if(!(*samples)[s]->getAone(mloc) && !(*samples)[s]->getAtwo(mloc)){
					a1+=2;
				}
				else if(!(*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc)){
					a1++;
					a2++;
				}
			}
			if(a2 < a1 && (*markers)[m]->getAllele2() != "0"){
				string temp = (*markers)[m]->getAllele1();
				//(*markers)[m]->setAllele1((*markers)[m]->getAllele2());
				//(*markers)[m]->setAllele2(temp);
				(*markers)[m]->resetAllele1((*markers)[m]->getAllele2());
				(*markers)[m]->resetAllele2(temp);
				for(int s = 0; s < ssize; s++){
					if((*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc)){
						(*samples)[s]->addAone(mloc, false);
						(*samples)[s]->addAtwo(mloc, false);
					}
					else if(!(*samples)[s]->getAone(mloc) && !(*samples)[s]->getAtwo(mloc)){
						(*samples)[s]->addAone(mloc, true);
						(*samples)[s]->addAtwo(mloc, true);
					}
				}
			}
		}
	}
}*/

/*
 *Function: assignLinks
 *Description:
 *Creates the family structure by assigning each sample their respective
 *parents, siblings and children
 *
 */
/*void assignLinks(vector<Family*>* families){
    vector<Family*>::iterator f_iter;

	int fsize = families->size();

	for(int f = 0; f < fsize; f++){
		Family* fam = (*families)[f];
		vector<Sample*>* samps = fam->getSamples();
		bool good = false;
		bool excluded = false;
		int ssize = samps->size();

		for(int s = 0; s < ssize; s++){
			Sample* samp = (*samps)[s];
			if(samp->isEnabled()){
				good = true;
			}
			else if(samp->isExcluded()){
				excluded = true;
			}
			if(samp->getDadID() == "0" && samp->getMomID() == "0"){
				fam->addFounder(samp);
				samp->setFounder(true);
			}
			else{
				fam->addNonFounder(samp);
				if(samp->getDadID() != "0"){
                	vector<Sample*>::iterator temp_iter = find_if(samps->begin(), samps->end(), FindSampleByID(samp->getDadID()));
                	if(temp_iter != samps->end()){
						Sample* dad = (*temp_iter);
	 					samp->setDad(dad);
						Sample* last_child = dad->getLastChild();
						if(last_child != NULL && last_child->getSib() == NULL){
							last_child->setSib(samp);
						}
						dad->addChild(samp);
                	}
				}
				if(samp->getMomID() != "0"){
                	vector<Sample*>::iterator temp_iter = find_if(samps->begin(), samps->end(), FindSampleByID(samp->getMomID()));
                	if(temp_iter != samps->end()){
						Sample* mom = (*temp_iter);
	 					samp->setMom(mom);
						Sample* last_child = mom->getLastChild();
						if(last_child != NULL && last_child->getSib() == NULL){
							last_child->setSib(samp);
						}
						mom->addChild(samp);
                	}
				}
			}
		}
		fam->setEnabled(good);
		if(!good && excluded){
			fam->setExcluded(excluded);
		}
	}
}*/


/*
 *Function: readMap
 *Description:
 *Reads the map file and performs the appropriate inclusion/exclusion of markers
 *
 */
void readMap(vector<Marker*>* markers, vector<int>* marker_map){
	map<string, vector<string> > descinfo;
	vector<string> descheaders;
	//vector<string> exclude;
	//vector<string> include;
	map<string,int> exclude;
	map<string,int> include;
	map<string,float> frequencies;

	//read in frequency file
	if(opts::_FREQ_FILE_EXISTS_){
		opts::printLog("Marker Minor Allele Frequencies being read from: " + opts::_FREQ_FILE_ + "\n");
		ifstream finput;
		finput.open(opts::_FREQ_FILE_.c_str(), ios::in);
		if(!finput){
			opts::printLog("Error opening marker frequency file: " + opts::_FREQ_FILE_ + "\n");
			exit(1);
		}
		string line = "";
		int freqline = 1;
		while(getline(finput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2){
				opts::printLog("Marker frequency file column size != 2 on line: " + line + " Skipping line!!\n");
				continue;
			}
			try{
				float freq = -1.0f;
				for(unsigned int c = 0; c < tokens[1].size(); c++){
					if(!isdigit(tokens[1][c]) && tokens[1][c] != '.'){
						throw "Number format exception";
					}
				}
				freq = std::atof(tokens[1].c_str());
				frequencies[tokens[0]] = freq;
			}catch(...){
				opts::printLog(tokens[1] + " is not a valid value for line: " + getString<int>(freqline) + "\n");
				exit(1);
			}
			freqline++;
		}
	}
	//read in map descriptive file
	if(opts::_MAPDESC_.length() > 0){
		opts::printLog("Marker description being used: " + opts::_MAPDESC_ + "\n");
		ifstream dinput;
		dinput.open(opts::_MAPDESC_.c_str(), ios::in);
		if(!dinput){
			opts::printLog("Error opening marker description file: " + opts::_MAPDESC_ + "\n");
			exit(1);
		}
		string head = "";
		getline(dinput, head, '\n');
		descheaders = General::ParseDelimitedLine(head);
		while(!dinput.eof()){
			string id = "";
			string line = "";
			getline(dinput, line, '\n');
			if(line == ""){
				continue;
			}
			if(line[0] == '#'){
				continue;
			}

			vector<string> tokens;
			tokens = General::ParseDelimitedLine(line);
			if(tokens.size() < 2){
				opts::printLog("Marker description column size is < 2: " + line + "\n");
				exit(1);
			}
			if(tokens.size() != descheaders.size()){
				opts::printLog("Line is not the same size as header line: " + line + "\n");
				exit(1);
			}

			id = tokens[0];
			descinfo[id] = tokens;
		}
		dinput.close();
	}
	//read in marker exclusion file
	if(opts::_EXCMARKERS_.length() > 0){
		opts::printLog("Excluding markers from file: " + opts::_EXCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker exclusion file: " + opts::_EXCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";

		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Excluding markers column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			exclude[tokens[0]] = 1;
			//exclude.push_back(probe);
		}
		einput.close();
	}
	//read in marker inclusion file
	if(opts::_INCMARKERS_.length() > 0){
		opts::printLog("Including markers from file: " + opts::_INCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker inclusion file: " + opts::_INCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";
		while(getline(einput,line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Including markers column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			include[tokens[0]] = 1;
			//include.push_back(probe);
		}
		einput.close();
	}


	ifstream input;
    input.open(opts::_MAPFILE_.c_str(), ios::in);

    if(!input){
		opts::printLog("Error opening map file: " + opts::_MAPFILE_ + "\n");
        exit(1);
    }
	int count = 0;
	int prev_chrom = -1;
	int prev_bploc = -1;
    while(!input.eof()){
        char buf[256];
        input.getline(buf, 256, '\n');

        string line = buf;

        if(line == ""){
            continue;
        }

        string temp;
        stringstream s(line);
        vector<string> elems;
        while(s >> temp){
            elems.push_back(temp);
        }

        if(elems.size() == 0){
            continue;
        }
        else if(elems.size() > 3){
			opts::printLog("Map file line has more than 3 elements: " + line + "\n");
            exit(1);
        }
		else if(elems.size() < 3 && elems.size() > 0){
			opts::printLog("Map file line has fewer than 3 elements: " + line + "\n");
			exit(1);
		}


        string chr = elems[0];
        string probe_id = elems[1];
        int bploc = atoi(elems[2].c_str());

		bool use = true;
		opts::_MARKERS_FOUND_++;
		if(exclude.size() > 0){
			//vector<string>::iterator found = find(exclude.begin(), exclude.end(), probe_id);
			map<string, int>::iterator found = exclude.find(probe_id);
			if(found != exclude.end()){
				use = false;
			}
		}
		if(include.size() > 0){
			//vector<string>::iterator found = find(include.begin(), include.end(), probe_id);
			map<string, int>::iterator found = include.find(probe_id);
			if(found == include.end()){
				use = false;
			}
		}

        Marker* m = new Marker(chr, probe_id, bploc);
		if(opts::_CHROM_LIMIT_){
			if(m->getChrom() != opts::_CHROM_){
				use = false;
			}
			else{
				if(opts::_BP_LOW_LIMIT_){
					if(bploc < opts::_BP_LOW_){
						use = false;
					}
				}
				if(opts::_BP_HIGH_LIMIT_){
					if(bploc > opts::_BP_HIGH_){
						use = false;
					}
				}
			}
		}
		if(opts::_BP_SPACE_LIMIT_){
			if(prev_bploc == -1){
				prev_bploc = bploc;
				prev_chrom = m->getChrom();
			}
			else{
				if(m->getChrom() == prev_chrom && ((bploc - prev_bploc) < opts::_BP_SPACE_)){
					use = false;
				}
				else{
					prev_bploc = bploc;
					prev_chrom = m->getChrom();
				}
			}
		}
		m->setEnabled(use);
		m->setLoc(count);
		m->setRSID(probe_id);


		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe_id);
			if(found != frequencies.end()){
				m->setMAF(frequencies[probe_id]);
				m->setFreqFlag(true);
			}
		}
		if(opts::_MAPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[probe_id];
			for(unsigned int i = 1; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					m->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					m->assignDetail(descheaders[i], "NA");
				}
			}
		}
        markers->push_back(m);
		count++;
    }

    input.clear();
    input.close();

	marker_map->resize(markers->size());

	//put markers in chrom/bploc order
	stable_sort(markers->begin(), markers->end(), less<Marker*>());

	for(unsigned int i =0; i < markers->size(); i++){
		(*marker_map)[(*markers)[i]->getLoc()] = i;//(*markers)[i]->getLoc();
	}

}

/*
 *Function: readString
 *Description:
 *Reads the next string from a file, space/tab delimited
 *
 *
 */
/*bool readString(FILE* fp, string* s){
    bool done = false;
    *s="";
    while (1)
    {
        char ch = fgetc(fp);
        if ( ch==' ' || ch == '\t' )
        {
            if (done)
                return true;
        }
        else if ( ch=='\n' || ch=='\r' || feof(fp) )
            return false;
        else
        {
            *s += ch;
            done = true;
        }
    }
}*/


/*
 *Function: writeBit
 *Parameters: sample vector, family vector, marker vector, marker map vector
 *Description:
 *Writes genotype (bed) file as binary genotype file based on Plink output.
 *Family output file (fam) is first 6 columns of original ped file
 *Map file (bim) is extended map file including alleles.
 *
 *
 */
void writeBit(vector<Sample*>* samples, vector<Family*>* families, vector<Marker*>* markers, vector<int>* marker_map){
	opts::printLog("Writing family information to " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".fam\n");
	ofstream BIT((opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".fam").c_str(), ios::out);
	if(!BIT){
		opts::printLog("Error opening " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".fam for family file creation.  Exiting!\n");
		exit(1);
	}

	int ssize = samples->size();

	for(int s = 0; s < ssize; s++){
		Sample* samp = (*samples)[s];

		BIT << samp->getFamID() << " "
			<< samp->getInd() << " "
			<< samp->getDadID() << " "
			<< samp->getMomID() << " ";
		if(samp->getSex()){
			BIT << "1 ";
		}
		else{
			BIT << "2 ";
		}

		BIT << samp->getPheno();
		//if(samp->getAffected()){
		//	BIT << "2";
		//}
		//else{
		//	BIT << "0";
		//}
		BIT << endl;
	}
	BIT.clear();
	BIT.close();


	opts::printLog("Writing map information to " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bim\n");
	BIT.open((opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bim").c_str(), ios::out);
	if(!BIT){
		opts::printLog("Error opening " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bim for map file creation.  Exiting!\n");
		exit(1);
	}

	int msize = markers->size();

	for(int m = 0; m < msize; m++){
		Marker* mark = (*markers)[m];

		BIT << mark->getChrom() << "\t";
		BIT << mark->getProbeID() << "\t"
			<< "0\t" //centimorgan column
			<< mark->getBPLOC() << "\t";
		if(mark->getAllele1().length() == 0){
			BIT << "0\t";
		}
		else{
			BIT << mark->getAllele1() << "\t";
		}
		if(mark->getAllele2().length() == 0){
			BIT << "0\t";
		}
		else{
			BIT << mark->getAllele2() << "\t";
		}
//		if(mark->getRSID() == ""){
//			BIT << ".\t";
//		}
//		else{
//			BIT << mark->getRSID() << "\t";
//		}
//		if(mark->getEnzyme() == ""){
//			BIT << ".";
//		}
//		else{
//			BIT << mark->getEnzyme();
//		}
		BIT << endl;
	}
	BIT.clear();
	BIT.close();

	opts::printLog("Writing genotype bitfile to " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bed\n");
	BIT.open((opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bed").c_str(), ios::out | ios::binary);
	if(!BIT){
		opts::printLog("Error opening " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bed for genotype file creation.  Exiting!\n");
		exit(1);
	}

	bitset<8> b;
	char ch[1];

  	b.reset();
    b.set(2);  b.set(3);  b.set(5);  b.set(6);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);

	b.reset();
	b.set(0);  b.set(1);  b.set(3);  b.set(4);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);

	b.reset();
    ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);

	for(int s = 0; s < ssize; s++){
		Sample* samp = (*samples)[s];

		for(int m = 0; m < msize;){
			bitset<8> b;
			b.reset();
			int c = 0;
			while(c < 8 && m < msize){
				Marker* mark = (*markers)[m];
				int loc = mark->getLoc();

				if(samp->getAone(loc))
					b.set(c);

				c++;
				if(samp->getAtwo(loc))
					b.set(c);
				c++;
				m++;
			}

			char ch[1];
			ch[0] = (char)b.to_ulong();
			BIT.write(ch, 1);
		}
	}
	BIT.close();
}

/*
 *Function: readBin
 *Paramters: sample vector, family vector, marker vector, marker_map vector
 *Description:
 *Reads set of binary input files (bim, bed, fam)
 *Based on Plink
 *
 */

void readBin(vector<Sample*>* samples, vector<Family*>* families, vector<Marker*>* markers, vector<int>* marker_map){
	//vector<string> exclude;
	map<string, vector<string> > mdescinfo;
	map<string,int> exclude;
	map<string,int> include;
	vector<string> mdescheaders;
	//vector<string> include;
	map<string,float> frequencies;
	if(opts::_FREQ_FILE_EXISTS_){
		opts::printLog("Marker Minor Allele Frequencies being read from: " + opts::_FREQ_FILE_ + "\n");
		ifstream finput;
		finput.open(opts::_FREQ_FILE_.c_str(), ios::in);
		if(!finput){
			opts::printLog("Error opening marker frequency file: " + opts::_FREQ_FILE_ + "\n");
			exit(1);
		}
		string line = "";
		int freqline = 1;
		while(getline(finput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2){
				opts::printLog("Marker frequency file column size != 2 on line: " + line + " Skipping line!!\n");
				continue;
			}
			try{
				float freq = -1.0f;
				for(unsigned int c = 0; c < tokens[1].size(); c++){
					if(!isdigit(tokens[1][c]) && tokens[1][c] != '.'){
						throw "Number format exception";
					}
				}
				freq = std::atof(tokens[1].c_str());
				frequencies[tokens[0]] = freq;
			}catch(...){
				opts::printLog(tokens[1] + " is not a valid value for line: " + getString<int>(freqline) + "\n");
				exit(1);
			}
			freqline++;
		}
	}
	if(opts::_MAPDESC_.length() > 0){
        opts::printLog("Marker description being used: " + opts::_MAPDESC_ + "\n");
        ifstream dinput;
        dinput.open(opts::_MAPDESC_.c_str(), ios::in);
        if(!dinput){
            opts::printLog("Error opening marker description file: " + opts::_MAPDESC_ + "\n");
            exit(1);
        }
        string head = "";
        getline(dinput, head, '\n');
        mdescheaders = General::ParseDelimitedLine(head);
        while(!dinput.eof()){
            string id = "";
            string line = "";
            getline(dinput, line, '\n');
            if(line == ""){
                continue;
            }
            if(line[0] == '#'){
                continue;
            }
            vector<string> tokens;
            tokens = General::ParseDelimitedLine(line);
            if(tokens.size() < 2){
                opts::printLog("Marker description column size is < 2: " + line + "\n");
                exit(1);
            }
            if(tokens.size() != mdescheaders.size()){
                opts::printLog("Line is not the same size as header line: " + line + "\n");
                exit(1);
			}
            id = tokens[0];
	        mdescinfo[id] = tokens;
	    }
	    dinput.close();
	}
	if(opts::_EXCMARKERS_.length() > 0){
		opts::printLog("Excluding markers from file: " + opts::_EXCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker exclusion file: " + opts::_EXCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Marker exclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			exclude[tokens[0]] = 1;
			//exclude.push_back(probe);
		}
		einput.close();
	}
	if(opts::_INCMARKERS_.length() > 0){
		opts::printLog("Including markers from file: " + opts::_INCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker exclusion file: " + opts::_INCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Marker inclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			include[tokens[0]] = 1;
			//include.push_back(probe);
		}
		einput.close();
	}

	opts::printLog("Reading map from " + opts::_BINPREFIX_ + ".bim\n");

	ifstream MAP((opts::_BINPREFIX_ + ".bim").c_str(), ios::in);
	if(!MAP){
		opts::printLog("Error opening map file: " + opts::_BINPREFIX_ + ".bim.  Exiting!\n");
		exit(1);
	}
	//MAP.clear();

	int count = 0;

		string chrom = "";
		string probe = "";
		int bploc = 0;
		string a1 = "";
		string a2 = "";
		string centi = "";
		string rsid = "";
		string enzyme = "";
		string line = "";
	int prev_bploc = -1;
	int prev_chrom = -1;
	while(getline(MAP, line)){
			//>> chrom
			//>> probe
			//>> bploc
			//>> a1
			//>> a2
			//>> rsid
			//>> enzyme){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 6){
			opts::printLog(".bim file column size != 6: " + line + " stopping!!\n");
			exit(0);
		}
		chrom = tokens[0];
		probe = tokens[1];
		centi = tokens[2];
		bploc = atoi(tokens[3].c_str());
		a1 = tokens[4];
		a2 = tokens[5];
		//rsid = tokens[5];
		//enzyme = tokens[6];
		//opts::_ENZYMES_ = true;
		if(rsid == "."){
			rsid = "";
		}
		if(enzyme == "."){
			enzyme = "";
			opts::_ENZYMES_ = false;
		}
		bool use = true;
		opts::_MARKERS_FOUND_++;

        if(exclude.size() > 0){
            //vector<string>::iterator found = find(exclude.begin(), exclude.end(), probe);
            map<string,int>::iterator found = exclude.find(probe);
            if(found != exclude.end()){
                use = false;
            }
			else{
				found = exclude.find(rsid);
				if(found != exclude.end()){
					use = false;
				}
			}
        }
        if(include.size() > 0){
            //vector<string>::iterator found = find(include.begin(), include.end(), probe);
            map<string,int>::iterator found = include.find(probe);
            if(found == include.end()){
                use = false;
				found = include.find(rsid);
				if(found == include.end()){
					use = false;
				}
				else{
					use = true;
				}
            }
        }

		Marker* m = new Marker(chrom, probe, bploc);
		if(opts::_CHROM_LIMIT_){
			if(m->getChrom() != opts::_CHROM_){
				use = false;
			}
			else{
				if(opts::_BP_LOW_LIMIT_){
					if(bploc < opts::_BP_LOW_){
						use = false;
					}
				}
				if(opts::_BP_HIGH_LIMIT_){
					if(bploc > opts::_BP_HIGH_){
						use = false;
					}
				}
			}
		}
        if(opts::_BP_SPACE_LIMIT_){
            if(prev_bploc == -1){
                prev_bploc = bploc;
                prev_chrom = m->getChrom();
            }
            else{
                if(m->getChrom() == prev_chrom && ((bploc - prev_bploc) < opts::_BP_SPACE_)){
                    use = false;
                }
                else{
                    prev_bploc = bploc;
                    prev_chrom = m->getChrom();
                }
            }
        }

		m->setEnabled(use);
		m->setLoc(count);
		m->setAllele1(a1);
		m->setAllele2(a2);
		m->setRSID(rsid);
		if(rsid == "." || rsid == ""){
			m->setRSID(probe);
		}
		m->setEnzyme(enzyme);
		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe);
			if(found == frequencies.end()){
				found = frequencies.find(m->getRSID());
				if(found != frequencies.end()){
					m->setMAF(frequencies[m->getRSID()]);
					m->setFreqFlag(true);
				}
			}
			else{
				m->setMAF(frequencies[m->getRSID()]);
				m->setFreqFlag(true);
			}
		}
        if(opts::_MAPDESC_.length() > 0 && mdescinfo.size() > 0){
			vector<string> tokens = mdescinfo[probe];
            for(unsigned int i = 1; i < mdescheaders.size(); i++){
	            if(tokens.size() == mdescheaders.size()){
		            m->assignDetail(mdescheaders[i], tokens[i]);
                }
                else{
	                m->assignDetail(mdescheaders[i], "NA");
                }
            }
        }

		markers->push_back(m);

		count++;
	}
	MAP.close();

	marker_map->resize(markers->size());
	stable_sort(markers->begin(), markers->end(), less<Marker*>());

	for(unsigned int i = 0; i < markers->size(); i++){
		(*marker_map)[(*markers)[i]->getLoc()] = i;
	}


	map<string, vector<string> > descinfo;
	vector<string> sexclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> finclude;
	vector<string> fexclude;
	vector<string> descheaders;
	if(opts::_SAMPDESC_.length() > 0){
		opts::printLog("Sample description being used: " + opts::_SAMPDESC_ + "\n");
		ifstream dinput;
		dinput.open(opts::_SAMPDESC_.c_str(), ios::in);
		if(!dinput){
			opts::printLog("Error opening sample description file: " + opts::_SAMPDESC_ + "\n");
			exit(1);
		}
        string head = "";
        getline(dinput, head, '\n');
        descheaders = General::ParseDelimitedLine(head);
		while(!dinput.eof()){
			string fam = "";
			string id = "";
			string center = "";
			string plate = "";
			string well = "";
//old        	char buf[256];
//old	        dinput.getline(buf, 256, '\n');
            string line = "";
			getline(dinput, line, '\n');

    	    //old string line = buf;

        	if(line == ""){
            	continue;
        	}
			vector<string> tokens = General::ParseDelimitedLine(line);
            if(tokens.size() < 2){
				opts::printLog("Sample description column size is < 2: " + line + " Skipping line!!\n");
	            continue;
	        }
	        if(tokens.size() != descheaders.size()){
				opts::printLog("Line is not the same size as header line: " + line + " Skipping line!!\n");
	            continue;
	        }

			fam = tokens[0];
			id = tokens[1];
			descinfo[fam + " " + id] = tokens;
//			center = tokens[2];
//			plate = tokens[3];
//			well = tokens[4];
//			descinfo[fam + " " + id] = center + " " + plate + " " + well;

//			descinfo.push_back(line);
		}
		dinput.close();
	}
	if(opts::_EXCSAMPS_.length() > 0){
		opts::printLog("Excluding samples from file: " + opts::_EXCSAMPS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCSAMPS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening sample exclusion file: " + opts::_EXCSAMPS_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2 && line != ""){
				opts::printLog("Sample exclusion column size != 2: " + line + " Skipping line!!\n");
				continue;
			}
			sexclude.push_back(tokens[0] + " " + tokens[1]);
		}
		einput.close();
	}
	if(opts::_EXCFAMILIES_.length() > 0){
		opts::printLog("Excluding families from file: " + opts::_EXCFAMILIES_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCFAMILIES_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening family exclusion file: " + opts::_EXCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Family exclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			fexclude.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCCENTERS_.length() > 0){
		opts::printLog("Including centers from file: " + opts::_INCCENTERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCCENTERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening center inclusion file: " + opts::_INCCENTERS_ + "\n");
			exit(1);
		}
		string center = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Center inclusion columns != 1: " + line + " Skipping line!!\n");
				continue;
			}
			inccenters.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCSAMPLES_.length() > 0){
		opts::printLog("Including samples from file: " + opts::_INCSAMPLES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCSAMPLES_.c_str(), ios::in);
		if(!iinput){
			opts::printLog("Error opening sample inclusion file: " + opts::_INCSAMPLES_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		while(getline(iinput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2 && line != ""){
				opts::printLog("Sample inclusion columns != 2: " + line + " Skipping line!!\n");
				continue;
			}
			sinclude.push_back(tokens[0] + " " + tokens[1]);
		}
		iinput.close();
	}
	if(opts::_INCFAMILIES_.length() > 0){
		opts::printLog("Including families from file: " + opts::_INCFAMILIES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCFAMILIES_.c_str(), ios::in);
		if(!iinput){
			opts::printLog("Error opening family inclusion file: " + opts::_INCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(iinput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Family inclusion columns != 1: " + line + " Skipping line!!\n");
				continue;
			}
			finclude.push_back(tokens[0]);
		}
		iinput.close();
	}

	opts::printLog("Reading family information from " + opts::_BINPREFIX_ + ".fam\n");

	ifstream PED;
	PED.open((opts::_BINPREFIX_ + ".fam").c_str());
	if(!PED){
		opts::printLog("Error opening family information file: " + opts::_BINPREFIX_ + ".fam.  Exiting!\n");
		exit(1);
	}
	PED.clear();
		string fam = "";
		string ind = "";
		string dad = "";
		string mom = "";
		string sex = "";
		string aff = "";
		line = "";
	while(getline(PED, line)){
	   	//	>> fam
		//	>> ind
		//	>> dad
		//	>> mom
		//	>> sex
		//	>> aff){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 6){
			opts::printLog("Family information file column size != 6: " + line + " Exitting!!\n");
			exit(0);
		}
		fam = tokens[0];
		ind = tokens[1];
		dad = tokens[2];
		mom = tokens[3];
		sex = tokens[4];
		aff = tokens[5];
		Sample* samp = new Sample();
		samp->setFamID(fam);
		samp->setInd(ind);
		samp->setEnabled(true);
		opts::_SAMPLES_FOUND_++;
		if(sexclude.size() > 0){
			vector<string>::iterator found = find(sexclude.begin(), sexclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found != sexclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				if(opts::_KEEP_EXC_SAMPLES_){
					samp->setExcluded(true);
				}
			}
		}
		if(sinclude.size() > 0){
			vector<string>::iterator found = find(sinclude.begin(), sinclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found == sinclude.end()){
				samp->setEnabled(false);
				if(opts::_KEEP_EXC_SAMPLES_){
					samp->setExcluded(true);
				}
			}
		}

		samp->setDadID(dad);
		samp->setMomID(mom);
		if(sex == "1"){
			samp->setSex(true);
		}
		else if(sex == "2"){
			samp->setSex(false);
		}

		if(aff == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atoi(aff.c_str()));
		if(opts::pedinfo.size() > 0){
			map<string, Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
			if(sfind != opts::pedinfo.end()){
				Sample* sfound = sfind->second;
				samp->setDadID(sfound->getDadID());
				samp->setMomID(sfound->getMomID());
				samp->setSex(sfound->getSex());
				samp->setPheno(sfound->getPheno());
			}
		}

		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(unsigned int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
				    samp->assignDetail(descheaders[i], "NA");
				}
			}

			//if(tokens.size() > 0){
			//	center = tokens[0];
			//}
			//if(opts::_INCCENTERS_.length() > 0 && inccenters.size() > 0){
			//	bool found = false;
			//	for(int cent = 0; cent < inccenters.size(); cent++){
			//		if(inccenters[cent] == center){
			//			found = true;
			//		}
			//	}
			//	if(!found){
			//		samp->setEnabled(false);
			//	}
//
			//}
			//samp->getFamily()->setCenter(tokens[2]);
			//if(tokens.size() > 0){
			//	samp->setPlate(tokens[1]);
			//	samp->setWell(tokens[2]);
			//}
		}
        vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Family* fam = new Family();
/*          fam->Setcenter()*/
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
			opts::_FAMILIES_FOUND_++;
        }
		if(fexclude.size() > 0){
			vector<string>::iterator found = find(fexclude.begin(), fexclude.end(), samp->getFamID());
			if(found != fexclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		if(finclude.size() > 0){
			vector<string>::iterator found = find(finclude.begin(), finclude.end(), samp->getFamID());
			if(found == finclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		samp->resizeAlleles(markers->size());
        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
	}

	PED.clear();
	PED.close();

	bool ind_major = false;
	bool snp_major = false;
	opts::printLog("Reading genotype bitfile from " + opts::_BINPREFIX_ + ".bed\n");

	ifstream BIT;
	BIT.open((opts::_BINPREFIX_+".bed").c_str(), ios::in | ios::binary);
	if(!BIT){
	opts::printLog("Error opening genotype bitfile: " + opts::_BINPREFIX_ + ".bed.  Exiting!\n");
		exit(1);
	}
	char temp[1];
	//header
	BIT.read(temp, 1);
	bitset<8> tb;
	tb = temp[0];
	if((tb[2] && tb[3] && tb[5] && tb[6]) && !(tb[0] || tb[1] || tb[4] || tb[7])){
		BIT.read(temp, 1);
		tb = temp[0];
		if((tb[0] && tb[1] && tb[3] && tb[4]) && !(tb[2] || tb[5] || tb[6] || tb[7])){
			BIT.read(temp, 1);
			tb = temp[0];
			if(!tb[0]){
			opts::printLog("IND major mode\n");
				ind_major = true;
			}
			else{
			opts::printLog("SNP major mode\n");
				snp_major = true;
			}
		}
		else{
			opts::printLog("Incorrect bit file version (2nd code)!\n");
			exit(1);
		}
	}
	else{
	opts::printLog("Incorrect bit file version (1st code)!\n");
		exit(1);
	}

	int ssize = samples->size();
	int msize = markers->size();

	if(ind_major){
		for(int s = 0; s < ssize; s++){
			Sample* samp = (*samples)[s];
			for(int m = 0; m < msize;){
				char ch[1];
				//cout << samples->size() << "\t" << samp->getFamID() << "\t" << samp->getInd() << "\t" << s << "\t" << (*markers)[m]->getProbeID() << "\t" << m << endl;
				BIT.read(ch, 1);
				if(!BIT){
				opts::printLog("Problem with the bed file.\n");
					exit(1);
				}
				bitset<8> b;
				b = ch[0];
				int c = 0;
				while(c < 7 && m < msize){
					Marker* mark = (*markers)[m];
					if(mark->isEnabled()){
						int mloc = (*markers)[m]->getLoc();
						samp->addAone(mloc, b[c++]);
						samp->addAtwo(mloc, b[c++]);
					}
					else{
//					samp->addAone(m, true);
//					samp->addAtwo(m, false);
						c+=2;
					}
					m++;
				}
			}
		}
	}
	else{//SNP_MAJOR
		for(int m = 0; m < msize; m++){
			Marker* mark = (*markers)[m];
			int mloc = mark->getLoc();
			for(int s = 0; s < ssize;){
				char ch[1];
				//cout << samples->size() << "\t" << samp->getFamID() << "\t" << samp->getInd() << "\t" << s << "\t" << (*markers)[m]->getProbeID() << "\t" << m << endl;
				BIT.read(ch, 1);
				if(!BIT){
					opts::printLog("Problem with the bed file.\n");
					exit(1);
				}
				bitset<8> b;
				b = ch[0];
				int c = 0;
				while(c < 7 && s < ssize){
					Sample* samp = (*samples)[s];
					if(samp->isEnabled() || (samp->isExcluded() && opts::_KEEP_EXC_SAMPLES_)){
						samp->addAone(mloc, b[c++]);
						samp->addAtwo(mloc, b[c++]);
					}
					else{
//					samp->addAone(m, true);
//					samp->addAtwo(m, false);
						c+=2;
					}
					s++;
				}
			}
		}
	}
	BIT.clear();
	BIT.close();
}

/*
 *Function: readPed
 *Description:
 *Reads ped file and stores data into sample, family, and marker vectors
 *
 */
void readPed(vector<Sample*>* samples, vector<Family*>* families, vector<Marker*>* markers, vector<int>* marker_map){
	map<string, vector<string> > descinfo;
	vector<string> exclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> fexclude;
	vector<string> finclude;
	vector<string> descheaders;
	if(opts::_SAMPDESC_.length() > 0){
		opts::printLog("Sample description being used: " + opts::_SAMPDESC_ + "\n");
		ifstream dinput;
		dinput.open(opts::_SAMPDESC_.c_str(), ios::in);
		if(!dinput){
		opts::printLog("Error opening sample description file: " + opts::_SAMPDESC_ + ".  Exiting!\n");
			exit(1);
		}
		string head = "";
		getline(dinput, head, '\n');
		descheaders = General::ParseDelimitedLine(head);
		while(!dinput.eof()){
			string fam = "";
			string id = "";
			string center = "";
			string plate = "";
			string well = "";
     //old   	char buf[256];
	 //old       dinput.getline(buf, 256, '\n');
	 		string line = "";
			getline(dinput, line, '\n');

//old    	    string line = buf;

        	if(line == ""){
            	continue;
        	}
//			vector<string> temp_tokens;
//			Tokenize(line, temp_tokens, " ");
			vector<string> tokens;
			tokens = General::ParseDelimitedLine(line);
			if(tokens.size() < 2){
				opts::printLog("Sample description column size is < 2: " + line + " Skipping line!!\n");
				continue;
			}
			if(tokens.size() != descheaders.size()){
			opts::printLog("Line is not the same size as header line: " + line + " Skipping line!!\n");
				continue;
			}

			fam = tokens[0];
			id = tokens[1];
			descinfo[fam + " " + id] = tokens;
//			for(int i = 0; i < temp_tokens.size(); i++){
//				Tokenize(temp_tokens[i], tokens, "\t");
//			}

/*old			if(tokens.size() != 5){
				cerr << "Sample description column size != 5: " << line << " Skipping line!!" << endl;
				continue;
			}
			fam = tokens[0];
			id = tokens[1];
			center = tokens[2];
			plate = tokens[3];
			well = tokens[4];
			descinfo[fam + " " + id] = center + " " + plate + " " + well;
*/

//			descinfo.push_back(line);
		}
		dinput.close();
	}
	if(opts::_EXCSAMPS_.length() > 0){
		opts::printLog("Excluding samples from file: " + opts::_EXCSAMPS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCSAMPS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening sample exclusion file: " + opts::_EXCSAMPS_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		//while(einput >> fam >> ind){
		while(getline(einput,line)){
//			vector<string> temp_tokens;
			vector<string> tokens;
			tokens = General::ParseDelimitedLine(line);
//			Tokenize(line, temp_tokens, " ");
//			for(int i = 0; i < temp_tokens.size(); i++){
//				Tokenize(temp_tokens[i], tokens, "\t");
//			}
			if(tokens.size() != 2 && line != ""){
			opts::printLog("Sample Exclusion column size != 2: " + line + " Skipping line!!\n");
				continue;
			}
			exclude.push_back(tokens[0] + " " + tokens[1]);
		}
		einput.close();
	}
	if(opts::_EXCFAMILIES_.length() > 0){
	opts::printLog("Excluding families from file: " + opts::_EXCFAMILIES_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCFAMILIES_.c_str(), ios::in);
		if(!einput){
		opts::printLog("Error opening family exclusion file: " + opts::_EXCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(einput, line)){
//			vector<string> temp_tokens;
			vector<string> tokens;
			tokens = General::ParseDelimitedLine(line);
//			Tokenize(line, temp_tokens, " ");
//			for(int i = 0; i < temp_tokens.size(); i++){
//				Tokenize(temp_tokens[i], tokens, "\t");
//			}
			if(tokens.size() != 1 && line != ""){
			opts::printLog("Family Exclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			fexclude.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCCENTERS_.length() > 0){
	opts::printLog("Including centers from file: " + opts::_INCCENTERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCCENTERS_.c_str(), ios::in);
		if(!einput){
		opts::printLog("Error opening center inclusion file: " + opts::_INCCENTERS_ + "\n");
			exit(1);
		}
		string center = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
			opts::printLog("Center Inclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			inccenters.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCSAMPLES_.length() > 0){
	opts::printLog("Including samples from file: " + opts::_INCSAMPLES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCSAMPLES_.c_str(), ios::in);
		if(!iinput){
		opts::printLog("Error opening sample inclusion file: " + opts::_INCSAMPLES_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		while(getline(iinput,line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2 && line != ""){
			opts::printLog("Sample inclusion column size != 2:" + line + " Skipping line!!\n");
				continue;
			}
			sinclude.push_back(tokens[0] + " " + tokens[1]);
		}
		iinput.close();
	}
	if(opts::_INCFAMILIES_.length() > 0){
	opts::printLog("Including families from file: " + opts::_INCFAMILIES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCFAMILIES_.c_str(), ios::in);
		if(!iinput){
		opts::printLog("Error opening family inclusion file: " + opts::_INCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(iinput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
			opts::printLog("Family inclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			finclude.push_back(tokens[0]);
		}
		iinput.close();
	}

    FILE* input;
    input = fopen(opts::_PEDFILE_.c_str(), "r");
	if(!input){
	opts::printLog("Error opening pedfile: " + opts::_PEDFILE_ + ".  Exiting!\n");
		exit(1);
	}
	int onind = -1;
    while(!feof(input)){
		onind++;
        Sample* samp = new Sample();
        int f = 0;
        string temp = "";
        if(readString(input, &temp)){
           samp->setFamID(temp);
              f++;
             temp = "";
         }

		string ftemp = samp->getFamID();
		if(samp->getFamID() == ""){
			delete(samp);
			continue;
		}
        if(ftemp.at(0) == '#'){
			delete(samp);
			while(fgetc(input) != '\n' && !feof(input)){}
            continue;
        }
		opts::_SAMPLES_FOUND_++;
        /*check for comments?*/
        string sex = "";
        string pheno = "";
        if(readString(input, &temp)){
            samp->setInd(temp);
            f++;
            temp = "";
        }
		samp->setEnabled(true);

		if(exclude.size() > 0){
			vector<string>::iterator found = find(exclude.begin(), exclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found != exclude.end()){
				if(!opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					continue;
				}
				else{
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
		}
		if(sinclude.size() > 0){
			vector<string>::iterator found = find(sinclude.begin(), sinclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found == sinclude.end()){
				if(!opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					continue;
				}
				else{
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
		}
        if(readString(input, &temp)){
            samp->setDadID(temp);
            f++;
            temp = "";
        }
        if(readString(input, &temp)){
            samp->setMomID(temp);
            f++;
            temp = "";
        }
        if(readString(input, &sex)) f++;
        if(readString(input, &pheno)) f++;

        if(sex == "1"){
            samp->setSex(true);
        }
        else if(sex == "2"){
            samp->setSex(false);
        }

		if(pheno == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atoi(pheno.c_str()));

		if(opts::pedinfo.size() > 0){
			map<string, Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
			if(sfind != opts::pedinfo.end()){
				Sample* sfound = sfind->second;
				samp->setDadID(sfound->getDadID());
				samp->setMomID(sfound->getMomID());
				samp->setSex(sfound->getSex());
				samp->setPheno(sfound->getPheno());
			}
		}

		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(unsigned int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					samp->assignDetail(descheaders[i], "NA");
				}
			}
			//if(tokens.size() > 0){
			//	center = tokens[0];
			//}
			//if(opts::_INCCENTERS_.length() > 0 && inccenters.size() > 0){
			//	bool found = false;
			//	for(int cent = 0; cent < inccenters.size(); cent++){
			//		if(inccenters[cent] == center){
			//			found = true;
			//		}
			//	}
			//	if(!found){
			//		delete(samp);
			//		while(fgetc(input) != '\n' && !feof(input)){}
			//		continue;
			//	}

			//}
			//samp->getFamily()->setCenter(tokens[2]);
			//if(tokens.size() > 0){
			//	samp->setPlate(tokens[1]);
			//	samp->setWell(tokens[2]);
			//}
		}
        samp->resizeAlleles(markers->size());

        unsigned int gn = 0;
        unsigned int i = 0;
        bool linedone = false;
        //bool fatal = false;

        string fmsg;
        while(!linedone){
            string one = "";
            string two = "";

            while(1){
                char ch = fgetc(input);

                if(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input)){
                    if(ch == '\n' || ch == '\r' || feof(input)){
                        linedone = true;
                    }

                    if(one.length() > 0){
                        gn++;
                        break;
                    }
                    if(ch == '\n' || ch == '\r' || feof(input)){
                        break;
                    }
                }
                else{
                    one += ch;
                }
            }
            if(!linedone){
                while(1){
                    char ch = fgetc(input);
                    if(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input)){
                        if(ch == '\n' || ch == '\r' || feof(input)){
                            linedone = true;
                        }
                        if(two.length() > 0){
                            gn++;
                            break;
                        }
                        if(ch == '\n' || ch == '\r' || feof(input)){
                            break;
                        }
                    }
                    else{
                        two += ch;
                    }
                }
                if(linedone && one.length() == 0 && two.length() == 0){
                    break;
                }

				//Marker* m = (*markers)[i];
				if(i > markers->size()){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>((2 * markers->size()) + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					opts::printLog(text);
					exit(0);
				}
				Marker* m = (*markers)[(*marker_map)[i]];
				if(m->isEnabled()){
					int oldallelecount = m->getNumAlleles();
	                if(one != opts::_NOCALL_){
						//new
						if(m->getAlleleLoc(one) < 0){
							m->addAllele(one);
						}
/*old
    	                if(one != m->getAllele1() && one != m->getAllele2()){
        	                if(m->getAllele1() == ""){
            	                m->setAllele1(one);
                	        }
                    	    else if(m->getAllele2() == ""){
                        	    m->setAllele2(one);
	                        }
    	                }
old*/
        	        }

            	    if(two != one){
                	    if(two != opts::_NOCALL_){
							//new
							if(m->getAlleleLoc(two) < 0){
								m->addAllele(two);
							}
/*old                    	    if(two != m->getAllele1() && two != m->getAllele2()){
                        	    if(m->getAllele1() == ""){
                            	    m->setAllele1(two);
	                            }
    	                        else if(m->getAllele2() == ""){
        	                        m->setAllele2(two);
            	                }
                	        }
old*/
                    	}
	                }


					if(m->getNumAlleles() <= 2){
   		 	            if(one == m->getAllele1() && two == m->getAllele1()){
   	    	     	        samp->addAone(i, false);
   	     		            samp->addAtwo(i, false);
   	     		            samp->addAmissing(i, false);
		                }
	    	            else if(one != opts::_NOCALL_ && two != opts::_NOCALL_ && one != two){
	        	            samp->addAone(i, false);
 		           	        samp->addAtwo(i, true);
 		           	        samp->addAmissing(i, false);
	                	}
		                else if(one == m->getAllele2() && two == m->getAllele2()){
	    	                samp->addAone(i, true);
	        	            samp->addAtwo(i, true);
	        	            samp->addAmissing(i, false);
	            	    }
	                	else if(one == opts::_NOCALL_ || two == opts::_NOCALL_){
		                    samp->addAone(i, true);
	    	                samp->addAtwo(i, true);
	    	                samp->addAmissing(i, true);
	                	}
					}
					else if(opts::_MICROSATS_){
						samp->addMicroSat(i);
						int loc1 = m->getAlleleLoc(one);
						int loc2 = m->getAlleleLoc(two);
						//samp->addAbone(loc1);
						//samp->addAbtwo(loc2);
						samp->addAbone(i, loc1);
						samp->addAbtwo(i, loc2);
						if(oldallelecount <= 2){
							remapSamples(samples, markers, marker_map, i);
						}
					}
					else if(m->getNumAlleles() > 2 && !opts::_MICROSATS_){
						string alleles = "";
						for(unsigned int aa = 0; aa < m->getAlleles().size(); aa++){
							alleles += " " + m->getAllele(aa);
						}
						opts::printLog("More than 2 unique alleles found for SNP: " + m->getProbeID() + ", PED File line: " + getString<int>(onind + 1) + ".  Microsatellites not specified. (Offending alleles: " + alleles + "\n");
						exit(1);
					}
				}
				else{
					samp->addAone(i,true);
					samp->addAtwo(i, true);
					samp->addAmissing(i, true);
				}
    	    	i++;
				if(i > markers->size()){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>((2 * markers->size()) + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					opts::printLog(text);

					exit(0);
				}
            }/*end !linedone*/
        }/*end while(1)*/
		if(gn != (2* markers->size())){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>(((2 * markers->size()) + 6));
					text += " columns but found ";
					text += getString<int>((f + gn));
					text += "\n";
					opts::printLog(text);
			exit(0);
		}

        vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Family* fam = new Family();
/*          fam->Setcenter()*/
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
			opts::_FAMILIES_FOUND_++;
        }
		if(fexclude.size() > 0){
			vector<string>::iterator found = find(fexclude.begin(), fexclude.end(), samp->getFamID());
			if(found != fexclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		if(finclude.size() > 0){
			vector<string>::iterator found = find(finclude.begin(), finclude.end(), samp->getFamID());
			if(found == finclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
    }/*end while(eof)*/

    fclose(input);
}

/*
 *Function: remapSamples
 *Description:
 *If -micro-sats is enabled, when a 3rd allele is found for a marker, this function
 *moves the Sample genotype storage for the specific marker from boolean storage to integer storage.
 *
 */
/*void remapSamples(vector<Sample*>* samples, vector<Marker*>* markers, vector<int>* marker_map, int loc){

	Marker* m = (*markers)[(*marker_map)[loc]];
	for(int i = 0; i < samples->size(); i++){
		Sample* samp = (*samples)[i];
		if(!samp->haveMicroSat(loc)){
			samp->addMicroSat(loc);
			if(!samp->getAone(loc) && !samp->getAtwo(loc)){
				samp->addAbone(loc, 0);
				samp->addAbtwo(loc, 0);
			}
			else if(samp->getAone(loc) && samp->getAtwo(loc)){
				samp->addAbone(loc, 1);
				samp->addAbtwo(loc, 1);
			}
			else if(!samp->getAone(loc) && samp->getAtwo(loc)){
				samp->addAbone(loc, 0);
				samp->addAbtwo(loc, 1);
			}
			else{
				cout << "Remapping sample: " << samp->toString() << "\t-1/-1\n";
				samp->addAbone(loc, -1);
				samp->addAbtwo(loc, -1);
			}
			samp->addAone(loc, true);
			samp->addAtwo(loc, false);
		}
	}
}*/

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
void printFamilies(vector<Family*>* families){
	int fsize = families->size();

	string fname = opts::_OUTPREFIX_ + "family_structure.txt";
	ofstream fout (fname.c_str());
	if(!fout.is_open()){
		opts::printLog("Unable to open " + fname + "\n");
		exit(1);
	}
	for(int f = 0; f < fsize; f++){
		Family* fam = (*families)[f];
		if(fam->isEnabled()){
			fout << "Family: " << fam->getFamID() << endl;
			fout << "----------------------------------------------------------" << endl;
			vector<Sample*>* founders = fam->getFounders();
			int size = founders->size();
			if(size > 0){
				for(int fs = 0; fs < size; fs++){
					Sample* founder = (*founders)[fs];
					if(founder->isEnabled() || (founder->isExcluded() && opts::_KEEP_EXC_SAMPLES_)){
						fout << "Founder: " << founder->getInd() << endl;
						fout << descendTree(founder, 0) << endl;
					}
				}
			}
			else{
				fout << "No Founders...\n";
				vector<Sample*>* samps = fam->getSamples();
				int ssize = samps->size();
				for(int ss = 0; ss < ssize; ss++){
					Sample* samp = (*samps)[ss];
					if(samp->isEnabled() || (samp->isExcluded() && opts::_KEEP_EXC_SAMPLES_)){
						fout << descendTree(samp, 0) << endl;
					}
				}
			}
			fout << endl << endl;
		}
	}

	fout.close();
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

		for(unsigned int i = 0; i < filenames["Marker"].size(); i++){
			string file = filenames["Marker"][i];
			string step = opts::filesteps[file];

			for(unsigned int j = 0; j < opts::fileheaders[file].size(); j++){
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
}

void readTPed(vector<Marker*>* markers, vector<Sample*>* samples, vector<int>* marker_map){
	map<string, vector<string> > descinfo;
	vector<string> descheaders;
	//vector<string> exclude;
	//vector<string> include;
	map<string,int> exclude;
	map<string,int> include;
	map<string,float> frequencies;
	vector<bool> includeme;

	//read in frequency file
	if(opts::_FREQ_FILE_EXISTS_){
		opts::printLog("Marker Minor Allele Frequencies being read from: " + opts::_FREQ_FILE_ + "\n");
		ifstream finput;
		finput.open(opts::_FREQ_FILE_.c_str(), ios::in);
		if(!finput){
			opts::printLog("Error opening marker frequency file: " + opts::_FREQ_FILE_ + "\n");
			exit(1);
		}
		string line = "";
		int freqline = 1;
		while(getline(finput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2){
				opts::printLog("Marker frequency file column size != 2 on line: " + line + " Skipping line!!\n");
				continue;
			}
			try{
				float freq = -1.0f;
				for(unsigned int c = 0; c < tokens[1].size(); c++){
					if(!isdigit(tokens[1][c]) && tokens[1][c] != '.'){
						throw "Number format exception";
					}
				}
				freq = std::atof(tokens[1].c_str());
				frequencies[tokens[0]] = freq;
			}catch(...){
				opts::printLog(tokens[1] + " is not a valid value for line: " + getString<int>(freqline) + "\n");
				exit(1);
			}
			freqline++;
		}
	}
	//read in map descriptive file
	if(opts::_MAPDESC_.length() > 0){
		opts::printLog("Marker description being used: " + opts::_MAPDESC_ + "\n");
		ifstream dinput;
		dinput.open(opts::_MAPDESC_.c_str(), ios::in);
		if(!dinput){
			opts::printLog("Error opening marker description file: " + opts::_MAPDESC_ + "\n");
			exit(1);
		}
		string head = "";
		getline(dinput, head, '\n');
		descheaders = General::ParseDelimitedLine(head);
		while(!dinput.eof()){
			string id = "";
			string line = "";
			getline(dinput, line, '\n');
			if(line == ""){
				continue;
			}
			if(line[0] == '#'){
				continue;
			}

			vector<string> tokens;
			tokens = General::ParseDelimitedLine(line);
			if(tokens.size() < 2){
				opts::printLog("Marker description column size is < 2: " + line + "\n");
				exit(1);
			}
			if(tokens.size() != descheaders.size()){
				opts::printLog("Line is not the same size as header line: " + line + "\n");
				exit(1);
			}

			id = tokens[0];
			descinfo[id] = tokens;
		}
		dinput.close();
	}
	//read in marker exclusion file
	if(opts::_EXCMARKERS_.length() > 0){
		opts::printLog("Excluding markers from file: " + opts::_EXCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker exclusion file: " + opts::_EXCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";

		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Excluding markers column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			exclude[tokens[0]] = 1;
			//exclude.push_back(probe);
		}
		einput.close();
	}
	//read in marker inclusion file
	if(opts::_INCMARKERS_.length() > 0){
		opts::printLog("Including markers from file: " + opts::_INCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker inclusion file: " + opts::_INCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";
		while(getline(einput,line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Including markers column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			include[tokens[0]] = 1;
			//include.push_back(probe);
		}
		einput.close();
	}


	FILE* input;
    input = fopen(opts::_TPEDFILE_.c_str(), "r");

    if(!input){
		opts::printLog("Error opening ped file: " + opts::_TPEDFILE_ + "\n");
        exit(1);
    }
	int count = 0;
	int prev_chrom = -1;
	int prev_bploc = -1;
    while(!feof(input)){
		int f = 0;
		string chr;
		string probe_id;
		string cm;
		string bp;
		if(readString(input, &chr)) f++;

		if(chr==""){
			continue;
		}

		if(readString(input, &probe_id)) f++;
		if(readString(input, &cm)) f++;
		if(readString(input, &bp)) f++;

		while(fgetc(input) != '\n' && !feof(input)){}

		int bploc = atoi(bp.c_str());

		bool use = true;
		opts::_MARKERS_FOUND_++;
		if(exclude.size() > 0){
			//vector<string>::iterator found = find(exclude.begin(), exclude.end(), probe_id);
			map<string, int>::iterator found = exclude.find(probe_id);
			if(found != exclude.end()){
				use = false;
			}
		}
		if(include.size() > 0){
			//vector<string>::iterator found = find(include.begin(), include.end(), probe_id);
			map<string, int>::iterator found = include.find(probe_id);
			if(found == include.end()){
				use = false;
			}
		}

        Marker* m = new Marker(chr, probe_id, bploc);
		if(opts::_CHROM_LIMIT_){
			if(m->getChrom() != opts::_CHROM_){
				use = false;
			}
			else{
				if(opts::_BP_LOW_LIMIT_){
					if(bploc < opts::_BP_LOW_){
						use = false;
					}
				}
				if(opts::_BP_HIGH_LIMIT_){
					if(bploc > opts::_BP_HIGH_){
						use = false;
					}
				}
			}
		}
		if(opts::_BP_SPACE_LIMIT_){
			if(prev_bploc == -1){
				prev_bploc = bploc;
				prev_chrom = m->getChrom();
			}
			else{
				if(m->getChrom() == prev_chrom && ((bploc - prev_bploc) < opts::_BP_SPACE_)){
					use = false;
				}
				else{
					prev_bploc = bploc;
					prev_chrom = m->getChrom();
				}
			}
		}
		includeme.push_back(use);
		if(!use){
	//		delete(m);
	//		continue;
		}
		m->setEnabled(use);
		m->setLoc(count);
		m->setRSID(probe_id);


		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe_id);
			if(found != frequencies.end()){
				m->setMAF(frequencies[probe_id]);
				m->setFreqFlag(true);
			}
		}
		if(opts::_MAPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[probe_id];
			for(unsigned int i = 1; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					m->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					m->assignDetail(descheaders[i], "NA");
				}
			}
		}
        markers->push_back(m);
		count++;
    }

    //input.clear();
    fclose(input);

	marker_map->resize(markers->size());

	//put markers in chrom/bploc order
	stable_sort(markers->begin(), markers->end(), less<Marker*>());

	for(unsigned int i =0; i < markers->size(); i++){
		(*marker_map)[(*markers)[i]->getLoc()] = i;//(*markers)[i]->getLoc();
	}

	count = 0;
	for(unsigned int i = 0; i < samples->size(); i++){
		(*samples)[i]->resizeAlleles(markers->size());
	}

	FILE* PED;
	PED = fopen(opts::_TPEDFILE_.c_str(), "r");
	int i =0;

	while(!feof(PED)){
		string dummy;
		string tempdummy;
		int f =0;
		if(readString(PED, &dummy)) f++;

		if(dummy==""){
			continue;
		}

		if(dummy.substr(0,1) == "#"){
			while(fgetc(PED) != '\n' && !feof(PED)){}
			continue;
		}

		if(!includeme[i]){
			while(fgetc(PED) != '\n' && !feof(PED)){}
			i++;
			continue;
		}

		if(readString(PED,&tempdummy)) f++; //probe
		if(readString(PED,&dummy)) f++; //cm
		if(readString(PED,&dummy)) f++; //bploc

		unsigned int gn = 0;
		unsigned int c=0; //ind count
		bool linedone = false;
		//bool fatal = false;
		string fmsg;
		while(!linedone){
			Sample* samp = (*samples)[c];
			string one = "";
			string two = "";
			while(1){
				char ch = fgetc(PED);
				if(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(PED)){
					if(ch == '\n' || ch == '\r' || feof(PED)){
						linedone = true;
					}
					if(one.length() > 0){
						gn++;
						break;
					}
					if(ch == '\n' || ch == '\r' || feof(PED)){
						break;
					}
				}
				else{
					one += ch;
				}
			}

			if(!linedone){
				while(1){
					char ch = fgetc(PED);
					if(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(PED)){
						if(ch == '\n' || ch == '\r' || feof(PED)){
							linedone = true;
						}
						if(two.length() > 0){
							gn++;
							break;
						}
						if(ch == '\n' || ch == '\r' || feof(PED)){
							break;
						}
					}
					else{
						two += ch;
					}
				}
			}

			if(linedone && one.length() == 0 && two.length() == 0){
				break;
			}

			if(includeme[i]){
				int k = (*marker_map)[i];
				Marker* mark = (*markers)[k];

				int oldallelecount = mark->getNumAlleles();
                if(one != opts::_NOCALL_){
					//new
					if(mark->getAlleleLoc(one) < 0){
						mark->addAllele(one);
					}
        	    }

            	if(two != one){
                    if(two != opts::_NOCALL_){
						//new
						if(mark->getAlleleLoc(two) < 0){
							mark->addAllele(two);
						}
	                }
				}
//cout << mark->toString() << "\t" << mark->getAllele1() << mark->getAllele2() << "\t" << samp->toString() << "\t" << one << two;
				if(mark->getNumAlleles() <= 2){
   		            if(one == mark->getAllele1() && two == mark->getAllele1()){
   	         	        samp->addAone(i, false);
   	    	            samp->addAtwo(i, false);
		            }
	                else if(one != opts::_NOCALL_ && two != opts::_NOCALL_ && one != two){
	       	            samp->addAone(i, false);
 		      	        samp->addAtwo(i, true);
	               	}
		            else if(one == mark->getAllele2() && two == mark->getAllele2()){
	                    samp->addAone(i, true);
	       	            samp->addAtwo(i, true);
	           	    }
	               	else if(one == opts::_NOCALL_ || two == opts::_NOCALL_){
		                samp->addAone(i, true);
	                    samp->addAtwo(i, false);
	       	        }
//cout << "\t" << samp->getAone(i) << samp->getAtwo(i) << endl;
				}
				else if(opts::_MICROSATS_){
					samp->addMicroSat(i);
					int loc1 = mark->getAlleleLoc(one);
					int loc2 = mark->getAlleleLoc(two);
					//samp->addAbone(loc1);
					//samp->addAbtwo(loc2);
					samp->addAbone(i, loc1);
					samp->addAbtwo(i, loc2);
					if(oldallelecount <= 2){
						remapSamples(samples, markers, marker_map, i);
					}
				}
				else if(mark->getNumAlleles() > 2 && !opts::_MICROSATS_){
					opts::printLog("More than 2 unique alleles found for map location: " + getString<int>(k) + ", line: " + getString<int>(c + 1) + " [" + samp->toString() + ":::" + mark->toString() + " " + one + "/" + two + ".  Microsatellites not specified.\n");
					exit(1);
				}
			}
			c++;
			if(c > samples->size()){
				opts::printLog("Problem with line " + getString<int>(i+1) + " in " + opts::_TPEDFILE_ + "\n");
				opts::printLog("Expecting 4 + 2 * " + getString<int>(samples->size()) + " = " +
						getString<int>(4+2*samples->size()) + " columns, but found more\n");
				exit(1);
			}
		}//line done? next snp

		if(gn != 2 * samples->size()){
			opts::printLog("Problem with line " + getString<int>(i+1) + " in " + opts::_TPEDFILE_ + "\n");
			opts::printLog("Expecting 4 + 2 * " + getString<int>(samples->size()) + " = " +
					getString<int>(4+2*samples->size()) + " columns, but found more\n");
			exit(1);
		}
		i++;
	}
	fclose(PED);
}

void readTFam(vector<Sample*>* samples, vector<Family*>* families){
	map<string, vector<string> > descinfo;
	vector<string> sexclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> finclude;
	vector<string> fexclude;
	vector<string> descheaders;
	if(opts::_SAMPDESC_.length() > 0){
		opts::printLog("Sample description being used: " + opts::_SAMPDESC_ + "\n");
		ifstream dinput;
		dinput.open(opts::_SAMPDESC_.c_str(), ios::in);
		if(!dinput){
			opts::printLog("Error opening sample description file: " + opts::_SAMPDESC_ + "\n");
			exit(1);
		}
        string head = "";
        getline(dinput, head, '\n');
        descheaders = General::ParseDelimitedLine(head);
		while(!dinput.eof()){
			string fam = "";
			string id = "";
			string center = "";
			string plate = "";
			string well = "";
//old        	char buf[256];
//old	        dinput.getline(buf, 256, '\n');
            string line = "";
			getline(dinput, line, '\n');

    	    //old string line = buf;

        	if(line == ""){
            	continue;
        	}
			vector<string> tokens = General::ParseDelimitedLine(line);
            if(tokens.size() < 2){
				opts::printLog("Sample description column size is < 2: " + line + " Skipping line!!\n");
	            continue;
	        }
	        if(tokens.size() != descheaders.size()){
				opts::printLog("Line is not the same size as header line: " + line + " Skipping line!!\n");
	            continue;
	        }

			fam = tokens[0];
			id = tokens[1];
			descinfo[fam + " " + id] = tokens;
//			center = tokens[2];
//			plate = tokens[3];
//			well = tokens[4];
//			descinfo[fam + " " + id] = center + " " + plate + " " + well;

//			descinfo.push_back(line);
		}
		dinput.close();
	}
	if(opts::_EXCSAMPS_.length() > 0){
		opts::printLog("Excluding samples from file: " + opts::_EXCSAMPS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCSAMPS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening sample exclusion file: " + opts::_EXCSAMPS_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2 && line != ""){
				opts::printLog("Sample exclusion column size != 2: " + line + " Skipping line!!\n");
				continue;
			}
			sexclude.push_back(tokens[0] + " " + tokens[1]);
		}
		einput.close();
	}
	if(opts::_EXCFAMILIES_.length() > 0){
		opts::printLog("Excluding families from file: " + opts::_EXCFAMILIES_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCFAMILIES_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening family exclusion file: " + opts::_EXCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Family exclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			fexclude.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCCENTERS_.length() > 0){
		opts::printLog("Including centers from file: " + opts::_INCCENTERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCCENTERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening center inclusion file: " + opts::_INCCENTERS_ + "\n");
			exit(1);
		}
		string center = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Center inclusion columns != 1: " + line + " Skipping line!!\n");
				continue;
			}
			inccenters.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCSAMPLES_.length() > 0){
		opts::printLog("Including samples from file: " + opts::_INCSAMPLES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCSAMPLES_.c_str(), ios::in);
		if(!iinput){
			opts::printLog("Error opening sample inclusion file: " + opts::_INCSAMPLES_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		while(getline(iinput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2 && line != ""){
				opts::printLog("Sample inclusion columns != 2: " + line + " Skipping line!!\n");
				continue;
			}
			sinclude.push_back(tokens[0] + " " + tokens[1]);
		}
		iinput.close();
	}
	if(opts::_INCFAMILIES_.length() > 0){
		opts::printLog("Including families from file: " + opts::_INCFAMILIES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCFAMILIES_.c_str(), ios::in);
		if(!iinput){
			opts::printLog("Error opening family inclusion file: " + opts::_INCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(iinput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Family inclusion columns != 1: " + line + " Skipping line!!\n");
				continue;
			}
			finclude.push_back(tokens[0]);
		}
		iinput.close();
	}

	opts::printLog("Reading family information from " + opts::_FAMFILE_ + "\n");

	ifstream PED;
	PED.open((opts::_FAMFILE_).c_str());
	if(!PED){
		opts::printLog("Error opening family information file: " + opts::_FAMFILE_ + ".  Exiting!\n");
		exit(1);
	}
	PED.clear();
		string fam = "";
		string ind = "";
		string dad = "";
		string mom = "";
		string sex = "";
		string aff = "";
		string line = "";
	while(getline(PED, line)){
	   	//	>> fam
		//	>> ind
		//	>> dad
		//	>> mom
		//	>> sex
		//	>> aff){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 6){
			opts::printLog("Family information file column size != 6: " + line + " Exitting!!\n");
			exit(0);
		}
		fam = tokens[0];
		ind = tokens[1];
		dad = tokens[2];
		mom = tokens[3];
		sex = tokens[4];
		aff = tokens[5];
		Sample* samp = new Sample();
		samp->setFamID(fam);
		samp->setInd(ind);
		samp->setEnabled(true);
		opts::_SAMPLES_FOUND_++;
		if(sexclude.size() > 0){
			vector<string>::iterator found = find(sexclude.begin(), sexclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found != sexclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				if(opts::_KEEP_EXC_SAMPLES_){
					samp->setExcluded(true);
				}
			}
		}
		if(sinclude.size() > 0){
			vector<string>::iterator found = find(sinclude.begin(), sinclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found == sinclude.end()){
				samp->setEnabled(false);
				if(opts::_KEEP_EXC_SAMPLES_){
					samp->setExcluded(true);
				}
			}
		}

		samp->setDadID(dad);
		samp->setMomID(mom);
		if(sex == "1"){
			samp->setSex(true);
		}
		else if(sex == "2"){
			samp->setSex(false);
		}

		if(aff == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atoi(aff.c_str()));
		if(opts::pedinfo.size() > 0){
			map<string, Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
			if(sfind != opts::pedinfo.end()){
				Sample* sfound = sfind->second;
				samp->setDadID(sfound->getDadID());
				samp->setMomID(sfound->getMomID());
				samp->setSex(sfound->getSex());
				samp->setPheno(sfound->getPheno());
			}
		}

		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(unsigned int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
				    samp->assignDetail(descheaders[i], "NA");
				}
			}

			//if(tokens.size() > 0){
			//	center = tokens[0];
			//}
			//if(opts::_INCCENTERS_.length() > 0 && inccenters.size() > 0){
			//	bool found = false;
			//	for(int cent = 0; cent < inccenters.size(); cent++){
			//		if(inccenters[cent] == center){
			//			found = true;
			//		}
			//	}
			//	if(!found){
			//		samp->setEnabled(false);
			//	}
//
			//}
			//samp->getFamily()->setCenter(tokens[2]);
			//if(tokens.size() > 0){
			//	samp->setPlate(tokens[1]);
			//	samp->setWell(tokens[2]);
			//}
		}
        vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Family* fam = new Family();
/*          fam->Setcenter()*/
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
			opts::_FAMILIES_FOUND_++;
        }
		if(fexclude.size() > 0){
			vector<string>::iterator found = find(fexclude.begin(), fexclude.end(), samp->getFamID());
			if(found != fexclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		if(finclude.size() > 0){
			vector<string>::iterator found = find(finclude.begin(), finclude.end(), samp->getFamID());
			if(found == finclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		//samp->resizeAlleles(markers->size());
        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
	}

	PED.clear();
	PED.close();
}




//////////////////////////////////////////////////////////////////
// Borrowed from Plink v0.99s                                   //
//                                                              //
//////////////////////////////////////////////////////////////////


#define  PORT_NUM                80
#define  IP_ADDR    "160.129.37.40"
#define  GET_STRING "GET /wasp/files/version.txt HTTP/1.1\nHost: chgr.mc.vanderbilt.edu\nConnection: close\n\n"


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
			   "\n*** See http://chgr.mc.vanderbilt.edu/wasp/\n");
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
						"\n*** See http://chgr.mc.vanderbilt.edu/wasp/\n");
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
							"\n*** See http://chgr.mc.vanderbilt.edu/wasp/\n");
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
		    "\nWasp has been instructed to stop"+
 	            "\nPlease see http://chgr.mc.vanderbilt.edu/wasp/\n");
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
		    "\nWasp has been instructed to stop"+
 	            "\nPlease see http://chgr.mc.vanderbilt.edu/wasp/\n");
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
		    "\nWasp has been instructed to stop"+
 	            "\nPlease see http://chgr.mc.vanderbilt.edu/wasp/\n");
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
	      if (tokens[i+1] == _WASPVER_)
		opts::printLog(" OK, v"+_WASPVER_+" is current\n");
	      else
		{
		  opts::printLog("\n\n          *** UPDATE REQUIRED ***\n\n");
		  opts::printLog("\tThis version        : "+_WASPVER_+"\n");
		  opts::printLog("\tMost recent version : "+tokens[i+1]+"\n\n");
		  opts::printLog("Please upgrade your version of Wasp as soon as possible!\n");
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

