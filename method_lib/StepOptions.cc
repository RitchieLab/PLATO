#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include "SquareRootTransform.h"
#include "BoxCoxTransform.h"
#include "LogTransform.h"
#include "Helpers.h"
#include "StepOptions.h"

using namespace std;
namespace Methods{

void StepOptions::printOptions(){
	map<string, Argument>::iterator iter;
	vector<string> alpha;
	cout << "Step Options:\n--------------------\n";
	for(iter = s_ArgVals.begin(); iter != s_ArgVals.end(); iter++){
		cout << iter->first << endl;
		alpha.push_back(iter->first);
	}


}

void StepOptions::setUp(string s){
	options = s;
	out = convertString(options);

	vector<string> tokens = General::ParseDelimitedLine(options);

	if(tokens.size() > 0){
		for(int i = 0; i < (int)tokens.size(); i++){
			map<string, Argument>::iterator found = s_ArgVals.find(tokens.at(i));
			if(found == s_ArgVals.end()){
				opts::printLog("Argument: " + tokens.at(i) + " is an invalid argument in " + options + "\n");
				throw MethodException("Argument: " + tokens.at(i) + " is an invalid argument in " + options + "\n");
			}
			switch(s_ArgVals[tokens.at(i)]){
				case s_lambda:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						double v = 0.0;
						istringstream test(val);
						if(!(test >> v)){
							throw MethodException(val + " is not a valid number!\n");
						}
						fixed_lambda = v;
						do_fixed_lambda = true;
					}
					else{
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;

				case s_thresh_markers_min:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							thresh_markers_min = std::strtod(val.c_str(), NULL);
							dothresh_markers_min = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_thresh_markers_max:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							thresh_markers_max = std::strtod(val.c_str(), NULL);
							dothresh_markers_max = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_thresh_samples_min:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							thresh_samples_min = std::strtod(val.c_str(), NULL);
							dothresh_samples_min = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_thresh_samples_max:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							thresh_samples_max = std::strtod(val.c_str(), NULL);
							dothresh_samples_max = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_thresh_families_min:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							thresh_families_min = std::strtod(val.c_str(), NULL);
							dothresh_families_min = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_thresh_families_max:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							thresh_families_max = std::strtod(val.c_str(), NULL);
							dothresh_families_max = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_ci:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							confidence_interval = std::strtod(val.c_str(), NULL);
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				//output synthesis view
				case s_output_synthview:
					setOutputSynthView(true);
					break;
				//eigenstrat
				case s_qtl:
					setQTL(true);
					break;
				case s_ancestry:
					setAncestry(true);
					break;
					//epistasis opts
				case s_epi_sets:
					if(i + 1 < (int) tokens.size()){
						string val = tokens.at(++i);
						setEpiSetsFilename(val);
						setDoEpiSets(true);
					}
					else{
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_epi_set_by_set:
					setDoEpiSetBySet(true);
					break;
				case s_epi_fast:
					setEpiFast(true);
					break;
				case s_epi_alpha1:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							epi_alpha1 = std::strtod(val.c_str(), NULL);
							if(Helpers::dGreater(epi_alpha1, 1.0f)){
								epi_filter = false;
							}
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_epi_alpha2:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							epi_alpha2 = std::strtod(val.c_str(), NULL);
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				//biofilter
				case s_bio_comparison_file:
					if(i + 1 < (int) tokens.size()){
						string val = tokens.at(++i);
						if(val.size() == 0){
							throw MethodException("File name required for: " + tokens.at(i-1) + "\n");
						}
						bio_comparison_file = val;
					}
					break;
				case s_bio_offset_begin:
					if(i + 1 < (int) tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c])){
									throw "oops!";
								}
							}
							bio_offset_begin = std::atoi(val.c_str());

						}
						catch(...){
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_bio_offset_end:
					if(i + 1 < (int) tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c])){
									throw "oops!";
								}
							}
							bio_offset_end = std::atoi(val.c_str());

						}
						catch(...){
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_bio_file_binary:
					bio_file_binary = true;
					break;
					//autosome only
				case s_autosome_only:
					autosome_only = true;
					break;
				case s_chrom:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							if(val == "X" || val == "x"){
								chrom = opts::_CHRX_;
							}
							else if(val == "Y" || val == "y"){
								chrom = opts::_CHRY_;
							}
							else if(val == "XY"){
								chrom = opts::_CHRXY_;
							}
							else{
								for(int c = 0; c < (int)val.size(); c++){
									if(!isdigit(val[c]) && val[c] != '.'){
										throw "oops!";
									}
								}
								chrom = std::atoi(val.c_str());
							}
							dochrom = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;

				case s_group_freq:
					if(i + 1 < (int) tokens.size()){
						string val = tokens.at(++i);
						setGroupFreqFile(val);
						setDoGroupFreq(true);
					}
					else{
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_override_out:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						override_out = val;
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_out:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						out = "_" + val;
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_sampbprange_filter:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						sampbprange_file = val;
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
                case s_rem_missing_parents:
                    rem_missing_parents = true;
                    break;
                case s_dummy_missing_parents:
                    dummy_missing_parents = true;
                    break;
                case s_zero_incomplete_trio_ids:
                    zero_incomplete_trio_ids = true;
                    break;
                case s_dummy_incomplete_parent_ids:
                	setDummyIncompleteParentIds(true);
                	break;
				case s_rm_mono:
					rm_mono = true;
					break;
				case s_rm_het_only:
					rm_het_only = true;
					break;
				case s_inc_excluded_samples:
					inc_excluded_samples = true;
					break;
				case s_zero_excluded:
					zero_excluded = true;
					break;
				case s_inc_disabled_samples:
					inc_disabled_samples = true;
					break;
				case s_zero_disabled:
					zero_disabled = true;
					break;
				case s_trio:
					trio = true;
					break;
				case s_me_zero:
					me_zero = true;
					break;
				case s_me_zero_l2:
					me_zero_l2 = true;
					break;
				case s_allele1234:
					allele1234 = true;
					break;
				case s_alleleACGT:
					alleleACGT = true;
					break;
				case s_allele12:
					allele12 = true;
					break;
				case s_config:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						setFilterProcessConfig(val);
						hasFilterProcessConfig(true);
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;

			  case s_ref_allele_included:
			    map_includes_ref_allele = true;
				break;
			  case s_hwept:
					hwept = true;
					break;
				case s_prevalence:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							prevalence = std::strtod(val.c_str(), NULL);
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				//transforms
				case s_sqrt_transform:
				{
					sqrt_transform = true;
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						vector<string> elems;
						General::Tokenize(val, to_transform, ",");
					}
					break;
				}
				case s_boxcox_transform:
				{
					boxcox_transform = true;
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						vector<string> elems;
						General::Tokenize(val, to_transform, ",");
					}
					break;
				}
				case s_log_transform:
				{
					log_transform = true;

					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						vector<string> elems;
						General::Tokenize(val, to_transform, ",");
					}
					break;
				}
				case s_perms:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							perms = std::atoi(val.c_str());
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_allelecustom:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							allelecustom = true;
							custom_alleles = Helpers::readCustomAlleles(val);
							map<string, string>::iterator c;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_me_zero_l2_fams:
					me_zero_l2_fams = true;
					break;
				case s_parents_only:
					parents_only = true;
					break;
				case s_bp_min:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							bp_min = std::atoi(val.c_str());
							dobp_min = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_bp_max:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							bp_max = std::atoi(val.c_str());
							dobp_max = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_bp_space:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							bp_space = std::atoi(val.c_str());
							dobp_space = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_disease:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							disease = val;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_ped_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							ped_file = val;
							do_pedfile = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_mdr_gui_output:
					mdr_gui_output = true;
					break;
				case s_mdr_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							mdr_file = val;
							do_mdrfile = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_mdr_ped_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							mdr_ped_file = val;
							do_mdrpedfile = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_mdr_map_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							mdr_map_file = val;
							do_mdrmapfile = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_map_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							map_file = val;
							do_mapfile = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_tped_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							tped_file = val;
							do_tpedfile = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_tfam_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							tfam_file = val;
							do_tfamfile = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_bin_prefix:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							bin_prefix = val;
							do_binfile = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_struct_strat_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							struct_strat_file = val;
							stratification = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_covar_missing:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							covar_missing = val;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_trait_missing:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							trait_missing = val;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_group_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							group_file = val;
							do_group_file = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_cluster_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							cluster_file = val;
							do_cluster_file = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_covar_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							cov_file = val;
							do_covs_file = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_covars_number:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							General::Tokenize(val, cov_use, ",");
							do_covs = true;
							do_covs_number = true;
							if(cov_use.size() == 0){
								throw("oops...");
							}

						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_covars_name:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							General::Tokenize(val, cov_use, ",");
							do_covs = true;
							do_covs_name = true;
							if(cov_use.size() == 0){
								throw("oops...");
							}

						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_trait_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							trait_file = val;
							do_traits_file = true;
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_traits_number:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							General::Tokenize(val, trait_use, ",");
							do_traits = true;
							do_traits_number = true;
							if(trait_use.size() == 0){
								throw("oops...");
							}

						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_traits_name:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							General::Tokenize(val, trait_use, ",");
							do_traits = true;
							do_traits_name = true;
							if(trait_use.size() == 0){
								throw("oops...");
							}

						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_loci:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							General::Tokenize(val, loci_use, ",");
							if(loci_use.size() == 0){
								throw("oops...");
							}

						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_pheno_missing:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(val[c] != '-'){
									if(!isdigit(val[c])){
										throw "oops!";
									}
								}
							}
							pheno_missing = atoi(val.c_str());

						}
						catch(...){
							opts::printLog(val + " is not a valid value for -pheno-missing\n");
							throw MethodException(val + " is not a valid value for -pheno-missing\n");
						}
					}
					break;
				case s_pheno:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						setPhenoName(val);
						setUsePheno(true);
					}
					break;

				case s_pheno_index:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
                        try{
                        	for(int c = 0; c < (int)val.size(); c++){
                        		if(!isdigit(val[c])){
                        			throw "oops!";
                        		}
                        	}
                        	pheno_loc = atoi(val.c_str());
                        	if(pheno_loc < 0){
                        		throw "oops!";
                        	}
                        	setUsePheno(true);
                        }catch(...){
                        	opts::printLog(val + " is not a valid value for -pheno-index\n");
                        	throw MethodException(val + " is not a valid value for -pheno-index\n");
                        }

					}
					break;
				case s_homozyg_wgha:
					homozyg_wgha = true;
					break;
                case s_homozyg_permute:
                    homozyg_permute = true;
                    if(i + 1 < (int)tokens.size()){
                        string val = tokens.at(++i);
                        try{
                            for(int c = 0; c < (int)val.size(); c++){
                                if(!isdigit(val[c]) && val[c] != '.'){
                                    throw "oops!";
                                }
                            }
                            homozyg_permute_count = std::atoi(val.c_str());
                        }catch(...){
                            opts::printLog(val + " is not a valid value for line: " + options + "\n");
                            throw MethodException(val + " is not a valid value for line: " + options + "\n");
                        }
                    }
                    else{
                        opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
                        throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
                    }
                    break;
				case s_homozyg_raw:
					homozyg_raw = true;
					break;
				case s_homozyg_zeros:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							homozyg_zeros = std::atoi(val.c_str());
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_homozyg_min_samp:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							homozyg_min_samp = std::atoi(val.c_str());
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_homozyg_span:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							homozyg_span = std::atoi(val.c_str());
						}catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_homozyg_seq:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							homozyg_seq_test = true;
							homozyg_seq_val = std::strtod(val.c_str(),NULL);
						}
						catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_ld_regular:
					if(i + 3 < (int)tokens.size()){
						string win = tokens.at(++i);
						string step = tokens.at(++i);
						string threshold = tokens.at(++i);
						try{
							for(int c = 0; c < (int)win.size(); c++){
								if(!isdigit(win[c]) && win[c] != '.'){
									throw "oops!";
								}
							}
							for(int c = 0; c < (int)step.size(); c++){
								if(!isdigit(step[c]) && step[c] != '.'){
									throw "oops!";
								}
							}
							for(int c = 0; c < (int)threshold.size(); c++){
								if(!isdigit(threshold[c]) && threshold[c] != '.'){
									throw "oops!";
								}
							}
							ld_vif = true;
							ld_window = std::atoi(win.c_str());
							ld_step = std::atoi(step.c_str());
							ld_threshold = std::strtod(threshold.c_str(), NULL);

						}catch(...){
							opts::printLog(win + ", or " + step + ", or " + threshold  + " is not a valid value for line: " + options + "\n");
							throw MethodException(win + ", or " + step + ", or " + threshold  + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_ld_calc_only:
					ld_calc_only = true;
					ld_window = 10;
					ld_step = 5;
					break;
				case s_ld_window:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							ld_window = std::atoi(val.c_str());
						}
						catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_ld_window_kb:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							ld_window_kb = std::atoi(val.c_str());
						}
						catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_ld_step:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							ld_step = std::atoi(val.c_str());
						}
						catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_ld_pairwise:
					if(i + 3 < (int)tokens.size()){
						string win = tokens.at(++i);
						string step = tokens.at(++i);
						string threshold = tokens.at(++i);
						try{
							for(int c = 0; c < (int)win.size(); c++){
								if(!isdigit(win[c]) && win[c] != '.'){
									throw "oops!";
								}
							}
							for(int c = 0; c < (int)step.size(); c++){
								if(!isdigit(step[c]) && step[c] != '.'){
									throw "oops!";
								}
							}
							for(int c = 0; c < (int)threshold.size(); c++){
								if(!isdigit(threshold[c]) && threshold[c] != '.'){
									throw "oops!";
								}
							}
							ld_pairwise = true;
							ld_window = std::atoi(win.c_str());
							ld_step = std::atoi(step.c_str());
							ld_pw_threshold = std::strtod(threshold.c_str(), NULL);
							ld_pw_threshold = (double)sqrt(ld_pw_threshold);

						}catch(...){
							opts::printLog(win + ", or " + step + ", or " + threshold  + " is not a valid value for line: " + options + "\n");
							throw MethodException(win + ", or " + step + ", or " + threshold  + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_ld_chop:
					if(tokens.size() > 0){
						vector<string>::iterator myfound = find(tokens.begin(), tokens.end(), "-ld-pairwise");
						if(myfound == tokens.end()){
							myfound = find(tokens.begin(), tokens.end(), "-ld-vif");
							if(myfound == tokens.end()){
								opts::printLog(tokens.at(i) + " requires either -ld-pairwise or -ld-vif on line: " + options + "\n");
								throw MethodException(tokens.at(i) + " requires either -ld-pairwise or -ld-vif on line: " + options + "\n");
							}
						}
					}
					ld_chop = true;

					break;
				case s_filter_overall:
					if(filter_file){
						opts::printLog(tokens.at(i) + " is not allowed when -filter-file is specified on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " is not allowed when -filter-file is specified on line: " + options + "\n");
					}
					filter_overall = true;
					break;
				case s_filter_file:
					if(filter_overall){
						opts::printLog(tokens.at(i) + " is not allowed when -filter-overall is specified on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " is not allowed when -filter-overall is specified on line: " + options + "\n");
					}
					if(!opts::_FREQ_FILE_EXISTS_){
						opts::printLog(tokens.at(i) + " is not allowed when -freq-file is not specified on the command line on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " is not allowed when -freq-file is not specified on the command line on line: " + options + "\n");
					}
					filter_file = true;
					break;
				case s_penetrance_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						penetrance_file = val;
						penetrance_codes = true;
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_center_file:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						center_file = val;
						center_codes = true;
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_unk_spouses:
					if(unaff_spouses_only || founders_only || all_children || all || random_child){
						opts::printLog(tokens.at(i) + " cannot be used with -all, -all-children, -founders-only, -unaff-spouses-only, or -random-child on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " cannot be used with -all, -all-children, -founders-only, -unaff-spouses-only, or -random-child on line: " + options + "\n");
					}
					unk_spouses = true;
					break;
				case s_unaff_spouses_only:
					if(unk_spouses || founders_only || all_children || all || random_child){
						opts::printLog(tokens.at(i) + " cannot be used with -all, -all-children, -founders-only, -unk-spouses, or -random-child on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " cannot be used with -all, -all-children, -founders-only, -unk-spouses, or -random-child on line: " + options + "\n");
					}
					unaff_spouses_only = true;
					break;
				case s_founders_only:
					if(unaff_spouses_only || unk_spouses || all_children || all || random_child){
						opts::printLog(tokens.at(i) + " cannot be used with -all, -all-children, -unk-spouses, -unaff-spouses-only or -random-child on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " cannot be used with -all, -all-children, -unk-spouses, -unaff-spouses-only or -random-child on line: " + options + "\n");
					}
					founders_only = true;
					break;
				case s_random_child:
					if(unk_spouses || unaff_spouses_only || founders_only || all_children || all){
						opts::printLog(tokens.at(i) + " cannot be used with -all, -all-children, -founders-only, -unk-spouses, -unaff-spouses-only on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " cannot be used with -all, -all-children, -founders-only, -unk-spouses, -unaff-spouses-only on line: " + options + "\n");
					}
					random_child = true;
					break;
				case s_all_children:
					if(unk_spouses || unaff_spouses_only || founders_only || all || random_child){
						opts::printLog(tokens.at(i) + " cannot be used with -all, -founders-only, -unk-spouses, -unaff-spouses-only or -random-child on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " cannot be used with -all, -founders-only, -unk-spouses, -unaff-spouses-only or -random-child on line: " + options + "\n");
					}
					all_children = true;
					break;
				case s_all:
					if(unk_spouses || unaff_spouses_only || founders_only || all_children || random_child){
						opts::printLog(tokens.at(i) + " cannot be used with -all-children, -founders-only, unk-spouses, -unaff-spouses-only or -random-child on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " cannot be used with -all-children, -founders-only, unk-spouses, -unaff-spouses-only or -random-child on line: " + options + "\n");
					}
					all = true;
					break;
				case s_deletion:
					deletion = true;
					break;
				case s_parental:
					parental = true;
					break;
				case s_gender:
					gender = true;
					break;
				case s_casecontrol:
					casecontrol = true;
					break;
				case s_deletion_span:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							deletion_span = std::atoi(val.c_str());
							deletion = true;
						}
						catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_random_repeat:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							random_markers = std::atoi(val.c_str());
							do_random_repeat = true;
						}
						catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_random_markers:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							random_markers = std::atoi(val.c_str());
							do_random_markers = true;
						}
						catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " is not a valid value for line: " + options + "\n");
					}
					break;
				case s_sets:
					if(i + 1 < (int)tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int)val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							sets = std::atoi(val.c_str());
							do_sets = true;
						}
						catch(...){
							opts::printLog(val + " is not a valid value for line: " + options + "\n");
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						opts::printLog(tokens.at(i) + " requires a value on line: " + options + "\n");
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;

					//random samples
				case s_rand_samps:
					if(i + 1 < (int) tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int) val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							rand_samps = std::atoi(val.c_str());
						}catch(...){
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_rand_samps_repeat:
					rand_samps_repeat = true;
					break;
				case s_sets_samps:
					if(i + 1 < (int) tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int) val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							sets_samps = std::atoi(val.c_str());
						}
						catch(...){
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_percent_samps:
					if(i + 1 < (int) tokens.size()){
						percent_samps.clear();
						string val = tokens.at(++i);
						vector<string> valtoks = General::ParseDelimitedLine(val, ",");
						for(int vt = 0; vt < (int)valtoks.size(); vt++){
							try{
								for(int c = 0; c < (int) valtoks.at(vt).size(); c++){
									if(!isdigit(valtoks.at(vt)[c])){
										throw "oops!";
									}
								}
								float percent = std::atof(valtoks.at(vt).c_str());
								if(percent > 100.0f || percent < 0.0f){
									throw MethodException("-percent-samps expects a values between 0-100\n");
								}

								percent_samps.push_back(percent);
							}
							catch(...){
								throw MethodException(valtoks.at(vt) + " is not a valid value for line: " + options + "\n");
							}
						}

					}
					else{
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
				case s_percent_cases:
					if(i + 1 < (int) tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int) val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							percent_cases = std::atof(val.c_str());
							if(percent_cases > 100.0f || percent_cases < 0.0f){
								throw MethodException("-percent-cases expects a value between 0-100\n");
							}
						}
						catch(...){
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}

					break;
				case s_percent_controls:
					if(i + 1 < (int) tokens.size()){
						string val = tokens.at(++i);
						try{
							for(int c = 0; c < (int) val.size(); c++){
								if(!isdigit(val[c]) && val[c] != '.'){
									throw "oops!";
								}
							}
							percent_controls = std::atof(val.c_str());
							if(percent_controls > 100.0f || percent_controls < 0.0f){
								throw MethodException("-percent-controls expects a value between 0-100\n");
							}
						}
						catch(...){
							throw MethodException(val + " is not a valid value for line: " + options + "\n");
						}
					}
					else{
						throw MethodException(tokens.at(i) + " requires a value on line: " + options + "\n");
					}
					break;
    // Linear Regression options
      case s_linr_mod_type:
		if(i + 1 < (int)tokens.size()){
			setLRModelType(tokens.at(++i));
		}
		else{
			opts::printLog(tokens.at(i) + " requires a parameter.\n");
			throw MethodException(tokens.at(i) + " requires a parameter.\n");
		}
        break;
      case s_linr_condition:
    	  if(i + 1 < (int)tokens.size()){
    		  setLinRCondition(true);
    		  setLinRConditionString(tokens.at(++i));
    	  }
    	  else{
    		  opts::printLog(tokens.at(i) + " requires a parameter.\n");
    		  throw MethodException(tokens.at(i) + " requires a parameter.\n");
    	  }
    	  break;
      case s_linr_condition_file:
    	  if(i + 1 < (int)tokens.size()){
    		  setLinRCondition(true);
    		  setLinRConditionFile(tokens.at(++i));
    	  }
    	  else{
    		  opts::printLog(tokens.at(i) + " requires a parameter.\n");
    		  throw MethodException(tokens.at(i) + " requires a parameter.\n");
    	  }
		  break;
      case s_linr_interaction:
    	  setLinRInteraction(true);
    	  break;
      case s_linr_no_main_snp:
    	  setLinRNoMainSnp(true);
    	  break;

    	  //concordance
      case s_unique_id:
    	  setUniqueId(true);
    	  break;
	  case s_inc_missing:
		  setIncMissing(true);
		  break;

    	  //ibs
      case s_ibs_pairs:
    	  if(i + 1 < (int)tokens.size()){
    		  setDoIBSPairs(true);
    		  setIBSPairsFile(tokens.at(++i));
    	  }
    	  else{
    		  opts::printLog(tokens.at(i) + " requires a parameter.\n");
    		  throw MethodException(tokens.at(i) + " requires a parameter.\n");
    	  }
		  break;
      case s_ibs_trio_pairs:
    	  if(i + 1 < (int)tokens.size()){
    		  setDoIBSTrioPairs(true);
    		  setIBSTrioPairsFile(tokens.at(++i));
    	  }
    	  else{
    		  opts::printLog(tokens.at(i) + " requires a parameter.\n");
    		  throw MethodException(tokens.at(i) + " requires a parameter.\n");
    	  }
		  break;
      case s_ibs_all_pairs:
    	  setDoIBSAllPairs(true);
    	  break;
      case s_ibs_all_trio_pairs:
    	  setDoIBSAllTrioPairs(true);
    	  break;
      case s_ibs_trios_raw:
    	  setIbsTriosRaw(true);
    	  break;
      case s_ibs_trio_trans:
    	  setDoIBSTrioTransmissions(true);
    	  break;

		  //outputmdr
      case s_no_rules:
    	  set_no_rules(true);
    	  break;
      case s_mdr_pedigree_output:
    	  setMDRPedigreeOutput(true);
    	  break;

      //outputped
      case s_allele_as_snp:
    	  set_allele_as_snp(true);
    	  break;

	// Logistic Regression options
      case s_lr_full_interaction:
        setLRFullInteraction(true);
        break;
      case s_lr_reduced_interaction:
        setLRFullInteraction(false);
        break;
      case s_lr_include_interaction:
        setLRIncludeInteractions(true);
        break;
      case s_lr_no_interaction:
        setLRIncludeInteractions(false);
        break;
      case s_lr_mod_type:
		if(i + 1 < (int)tokens.size()){
			setLRModelType(tokens.at(++i));
		}
		else{
			opts::printLog(tokens.at(i) + " requires a parameter.\n");
			throw MethodException(tokens.at(i) + " requires a parameter.\n");
		}
        break;
      case s_lr_max_iter:
		    if(i + 1 < (int)tokens.size()){
				  string val = tokens.at(++i);
					try{
					  for(int c = 0; c < (int)val.size(); c++){
						  if(!isdigit(val[c]) && val[c] != '.'){
							  throw "oops!";
							}
						}
						setLRMaximumIterations(std::atoi(val.c_str()));
					}catch(...){
						opts::printLog(val + " is not a valid value for line: " + options + "\n");
						throw MethodException(val + " is not a valid value for line: " + options + "\n");
					}
			  }
        break;
        //MARS
      case s_mars_maxterms:
		    if(i + 1 < (int)tokens.size()){
				  string val = tokens.at(++i);
					try{
					  for(int c = 0; c < (int)val.size(); c++){
						  if(!isdigit(val[c]) && val[c] != '.'){
							  throw "oops!";
							}
						}
						setMarsMaxTerms(std::atoi(val.c_str()));
					}catch(...){
						opts::printLog(val + " is not a valid value for line: " + options + "\n");
						throw MethodException(val + " is not a valid value for line: " + options + "\n");
					}
			  }
        break;
      case s_mars_degree:
    	  if(i + 1 < (int) tokens.size()){
    		  string val = tokens.at(++i);
    		  try{
    			  for(int c = 0; c < (int) val.size(); c++){
    				  if(!isdigit(val[c]) && val[c] != '.'){
    					  throw "oops!";
    				  }
    			  }
    			  setMarsDegree(std::atoi(val.c_str()));
    		  }catch(...){
    			  opts::printLog(val + " is not a valid value for line: " + options + "\n");
    			  throw MethodException(val + " is not a valid value for line: " + options + "\n");
    		  }
    	  }
    	  break;
      case s_mars_nprune:		//added 09-28-2010 mcc
			if(i + 1 < (int)tokens.size()){
				  string val = tokens.at(++i);
					try{
					  for(int c = 0; c < (int)val.size(); c++){
						  if(!isdigit(val[c]) && val[c] != '.'){
							  throw "oops!";
							}
						}
						setMarsNPrune(std::atoi(val.c_str()));
						cout << "setting MARS NPRUNE=" + val + "\n";
						opts::printLog("setting MARS NPRUNE=" + val + "\n");
					}catch(...){
						opts::printLog(val + " is not a valid value for line: " + options + "\n");
						throw MethodException(val + " is not a valid value for line: " + options + "\n");
					}
			  }
              break;


        // Conditional Logistic Regression options
      case s_cond_lr_full_interaction:
        setCondLRFullInteraction(true);
        break;
      case s_cond_lr_reduced_interaction:
        setCondLRFullInteraction(false);
        break;
      case s_cond_lr_include_interaction:
        setCondLRIncludeInteractions(true);
        break;
      case s_cond_lr_no_interaction:
        setCondLRIncludeInteractions(false);
        break;
      case s_cond_lr_mod_type:
        setCondLRModelType(tokens.at(++i));
        break;
      case s_cond_lr_max_iter:
		    if(i + 1 < (int)tokens.size()){
				  string val = tokens.at(++i);
					try{
					  for(int c = 0; c < (int)val.size(); c++){
						  if(!isdigit(val[c]) && val[c] != '.'){
							  throw "oops!";
							}
						}
						setCondLRMaximumIterations(std::atoi(val.c_str()));
					}catch(...){
						opts::printLog(val + " is not a valid value for line: " + options + "\n");
						throw MethodException(val + " is not a valid value for line: " + options + "\n");
					}
			  }
        break;
      // MDR options
      case s_mdr_only_set_thresh:
        setOnlySetThreshold(true);
      break;
      case s_uncert_coeff_total:
        setUncertaintyCoeffTotalType(tokens.at(++i));
      break;
      case s_llr_total:
        setLikelihoodRatioTotalType(tokens.at(++i));
      break;
      case s_or_total:
        setOddsRatioTotalType(tokens.at(++i));
      break;
      case s_nmi_total:
        setNMITotalType(tokens.at(++i));
      break;
      case s_nmi_transposed:
        setNMITransposed(true);
      break;

      //CMH
      case s_cmh_22k:
    	  setCMH2x2xK(true);
    	  break;
      case s_cmh_ijk:
    	  setCMHIxJxK(true);
    	  break;
      case s_cmh_ordinal:
    	  setCMHOrdinal(true);
    	  break;
      case s_mult_compare:
    	  setMultCompare(true);
    	  break;

    	// MDRPDT
    	case s_mdrpdt_ptests:
    	  if(i + 1 < (int)tokens.size()){
				  string val = tokens.at(++i);
					try{
					  for(int c = 0; c < (int)val.size(); c++){
						  if(!isdigit(val[c]) && val[c] != '.'){
							  throw "oops!";
							}
						}
						setMDRPDTNumPTests(std::atoi(val.c_str()));
					}catch(...){
						opts::printLog(val + " is not a valid value for line: " + options + "\n");
						throw MethodException(val + " is not a valid value for line: " + options + "\n");
					}
			  }
    	  break;
    	case s_mdrpdt_randseed:
    	  if(i + 1 < (int)tokens.size()){
				  string val = tokens.at(++i);
					try{
					  for(int c = 0; c < (int)val.size(); c++){
						  if(!isdigit(val[c]) && val[c] != '.'){
							  throw "oops!";
							}
						}
						setMDRPDTRandSeed(std::atoi(val.c_str()));
					}catch(...){
						opts::printLog(val + " is not a valid value for line: " + options + "\n");
						throw MethodException(val + " is not a valid value for line: " + options + "\n");
					}
			  }
    	  break;
    	case s_mdrpdt_xvcount:
    	  if(i + 1 < (int)tokens.size()){
				  string val = tokens.at(++i);
					try{
					  for(int c = 0; c < (int)val.size(); c++){
						  if(!isdigit(val[c]) && val[c] != '.'){
							  throw "oops!";
							}
						}
						setMDRPDTNumCrossVals(std::atoi(val.c_str()));
					}catch(...){
						opts::printLog(val + " is not a valid value for line: " + options + "\n");
						throw MethodException(val + " is not a valid value for line: " + options + "\n");
					}
			  }
      	break;
      case s_mdrpdt_mincombo:
    	  if(i + 1 < (int)tokens.size()){
				  string val = tokens.at(++i);
					try{
					  for(int c = 0; c < (int)val.size(); c++){
						  if(!isdigit(val[c]) && val[c] != '.'){
							  throw "oops!";
							}
						}
						setMDRPDTMinCombo(std::atoi(val.c_str()));
					}catch(...){
						opts::printLog(val + " is not a valid value for line: " + options + "\n");
						throw MethodException(val + " is not a valid value for line: " + options + "\n");
					}
			  }
        break;
      case s_mdrpdt_maxcombo:
    	  if(i + 1 < (int)tokens.size()){
				  string val = tokens.at(++i);
					try{
					  for(int c = 0; c < (int)val.size(); c++){
						  if(!isdigit(val[c]) && val[c] != '.'){
							  throw "oops!";
							}
						}
						setMDRPDTMaxCombo(std::atoi(val.c_str()));
					}catch(...){
						opts::printLog(val + " is not a valid value for line: " + options + "\n");
						throw MethodException(val + " is not a valid value for line: " + options + "\n");
					}
			  }
        break;
        break;
		  default:
					opts::printLog("Error reading step specific arguments, line: " + options + "\n");
					throw MethodException("Error reading step specific arguments, line: " + options + "\n");
			}
		}

		//error check pairs of options?
		string error = "";
		if(sets_samps > 0 && rand_samps == 0){
			error = "-sets-samps requires -rand-samps option as well.";
		}
		if(rand_samps > 0 && sets_samps == 0){
			sets_samps = 1;
		}
		if((int) percent_samps.size() > 0 && sets_samps > 0 && (int) percent_samps.size() != sets_samps){
			error = "-percent-samps needs same number of percentages as there are sets specified.";
		}
		if(error.size() > 0){
			throw MethodException(error);
		}
	}
}

int StepOptions::findStratification(string s){
	map<string, int>::iterator found = strat_map.find(s);
	if(found == strat_map.end()){
		return -1;
	}
	else{
		return found->second;
	}
}
string StepOptions::findCenterCode(string s){
	map<string, string>::iterator found = center_map.find(s);
	if(found == center_map.end()){
		return "";
	}
	else{
		return found->second;
	}
}

float StepOptions::findPenetranceCode(string s){
	map<string, float>:: iterator found = penetrance_map.find(s);
	if(found == penetrance_map.end()){
		return 0;
	}
	else{
		return found->second;
	}
}

void StepOptions::readEpiSets(string f){
	epi_sets.clear();
	ifstream input;
	input.open(f.c_str(), ios::in);
	if(!input){
		throw MethodException("Error opening epistasis sets file: " + f + "\n");
	}
	string line = "";
	int count = 1;
	while(getline(input, line)){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() < 2){
			throw MethodException("Epistasis sets file column size < 2 on line: " + getString<int>(count) + "\n");
		}
		string snp = tokens.at(0);
		for(int i = 1; i < (int)tokens.size(); i++){
			epi_sets[tokens.at(i)].push_back(snp);
		}
	}
	if(input.is_open()){
		input.close();
	}
}

void StepOptions::readStratificationFile(string f){
    strat_map.clear();
	ifstream input;
    input.open(f.c_str(), ios::in);
    if(!input){
        opts::printLog("Error opening penetrance file: " + f + "\n");
        throw MethodException("Error opening penetrance file: " + f + "\n");
    }
    string line = "";
	int count = 1;
    while(getline(input, line)){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 3){
			opts::printLog("Stratification file column size != 3 on line: " + line + " Exiting!\n");
			throw MethodException("Stratification file column size != 3 on line: " + line + " Exiting!\n");
		}
        try{
			strat_map[tokens.at(0) + " " + tokens.at(1)] = atoi(tokens.at(2).c_str());
        }catch(...){
           	opts::printLog("Column " + getString<int>(3) +" on line: " + getString<int>(count) + " is not a number!?\n");
           	throw MethodException("Column " + getString<int>(3) +" on line: " + getString<int>(count) + " is not a number!?\n");
        }
		count++;
    }
    if(input.is_open()){
        input.close();
    }
}

void StepOptions::readPenetranceFile(string f){
    penetrance_map.clear();
	ifstream input;
    input.open(f.c_str(), ios::in);
    if(!input){
        opts::printLog("Error opening penetrance file: " + f + "\n");
        throw MethodException("Error opening penetrance file: " + f + "\n");
    }
    string line = "";
	int count = 1;
    while(getline(input, line)){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 3){
			opts::printLog("Penetrance file column size != 3 one line: " + line + " Exiting!\n");
			throw MethodException("Penetrance file column size != 3 one line: " + line + " Exiting!\n");
		}
        try{
			penetrance_map[tokens.at(0) + " " + tokens.at(1)] = atof(tokens.at(2).c_str());
        }catch(...){
           	opts::printLog("Column " + getString<int>(3) +" on line: " + getString<int>(count) + " is not a number!?\n");
           	throw MethodException("Column " + getString<int>(3) +" on line: " + getString<int>(count) + " is not a number!?\n");
        }
		count++;
    }
    if(input.is_open()){
        input.close();
    }
}

void StepOptions::readCenterFile(string f){
    center_map.clear();
	ifstream input;
    input.open(f.c_str(), ios::in);
    if(!input){
        opts::printLog("Error opening center code file: " + f + "\n");
        throw MethodException("Error opening center code file: " + f + "\n");
    }
    string line = "";
    while(getline(input, line)){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 3){
			opts::printLog("Center code file column size != 3 one line: " + line + " Exiting!\n");
			throw MethodException("Center code file column size != 3 one line: " + line + " Exiting!\n");
		}
		center_map[tokens.at(0) + " " + tokens.at(1)] = tokens.at(2);
    }
    if(input.is_open()){
        input.close();
    }
}

bool StepOptions::checkChrom(int c){
	if(c == chrom){
		return true;
	}
	return false;
}

bool StepOptions::checkBp(int b){
	if(dobp_min){
		if(b < bp_min){
			return false;
		}
	}
	if(dobp_max){
		if(b > bp_max){
			return false;
		}
	}
	return true;
}

string StepOptions::convertString(string o){
	vector<string> tokens = General::ParseDelimitedLine(o);
	string temp = "";
	for(int i = 0; i < (int)tokens.size(); i++){
		map<string, Argument>::iterator found = s_subsetVals.find(tokens.at(i));
		if(found != s_subsetVals.end()){
			continue;
		}
		for(int j = 0; j < (int)tokens.at(i).size();){
			if(tokens.at(i).at(j) == '-' || tokens.at(i).at(j) == '/' || tokens.at(i).at(j) == '.' || tokens.at(i).at(j) == '~'){
				tokens.at(i).erase(j,1);
			}
			else{
				j++;
			}
		}
		temp += "_" + tokens.at(i);
	}
	return temp;
}

void StepOptions::readClusters(vector<Sample*>* samps){
	ifstream in(cluster_file.c_str(), ios::in);
	if(!in){
		opts::printLog("Error opening cluster file: " + cluster_file + "\n");
		throw MethodException("Error opening cluster file: " + cluster_file + "\n");
	}

	int count = 0;
	int cluster_count = 0;
	while(!in.eof()){
		count++;
		string line;
		getline(in, line);
		if(line == ""){
			continue;
		}
		if(line.at(0) == '#'){
			continue;
		}
		vector<string> elems;
		stringstream ss(line);
		string tok;
		while(ss >> tok){
			elems.push_back(tok);
		}
        if(elems.size() <= 2){
            opts::printLog("Line: " + getString<int>(count) + " has incorrect number of columns, in file " + cluster_file + "\n");
            in.close();

            throw MethodException("Line: " + getString<int>(count) + " has incorrect number of columns, in file " + cluster_file + "\n");
        }
		string fam = elems.at(0);
		string ind = elems.at(1);
		string cluster = elems.at(2);

		vector<Sample*>::iterator found = find_if(samps->begin(), samps->end(), FindSampleByFamAndID(fam, ind));
		if(found != samps->end()){
			Sample* mysamp = *found;
			map<string, vector<Sample*> >::iterator mfound = clusters.find(cluster);
			if(mfound != clusters.end()){
				vector<Sample*> gsamps = clusters[cluster];
				found = find(gsamps.begin(), gsamps.end(), mysamp);
				if(found == gsamps.end()){
					clusters[cluster].push_back(mysamp);
					sample_clusters[mysamp] = cluster_map[cluster];
				}
			}
			else{
				clusters[cluster].push_back(mysamp);
				cluster_map.insert(make_pair(cluster, cluster_count));
				sample_clusters[mysamp] = cluster_count;
				cluster_count++;
			}
		}
	}
}

void StepOptions::readClustersFromString(vector<Sample*>* samps){
	if(sample_clusters_string.size() == 0){
		opts::printLog("No clusters defined!\n");
		throw MethodException("No clusters defined!\n");
	}

	map<string, string>::iterator citer;
	int count = 0;
	int cluster_count = 0;
	for(citer = sample_clusters_string.begin(); citer != sample_clusters_string.end(); citer++){
		count++;
		string key = citer->first;
		string cluster = citer->second;

		vector<string> elems;
		General::Tokenize(key, elems, "#");

		string fam = elems.at(0);
		string ind = elems.at(1);

		vector<Sample*>::iterator found = find_if(samps->begin(), samps->end(), FindSampleByFamAndID(fam, ind));
		if(found != samps->end()){
			Sample* mysamp = *found;
			map<string, vector<Sample*> >::iterator mfound = clusters.find(cluster);
			if(mfound != clusters.end()){
				vector<Sample*> gsamps = clusters[cluster];
				found = find(gsamps.begin(), gsamps.end(), mysamp);
				if(found == gsamps.end()){
					clusters[cluster].push_back(mysamp);
					sample_clusters[mysamp] = cluster_map[cluster];
				}
			}
			else{
				clusters[cluster].push_back(mysamp);
				cluster_map.insert(make_pair(cluster, cluster_count));
				sample_clusters[mysamp] = cluster_count;
				cluster_count++;
			}
		}
	}
}


void StepOptions::readGroups(vector<Sample*>* samps){
	ifstream in(group_file.c_str(), ios::in);
	if(!in){
		opts::printLog("Error opening group file: " + group_file + "\n");
		throw MethodException("Error opening group file: " + group_file + "\n");
	}

	groups.clear();
	sample_groups.clear();
	group_families.clear();

	int count = 0;
	while(!in.eof()){
		count++;
		string line;
		getline(in, line);
		if(line == ""){
			continue;
		}
		if(line.at(0) == '#'){
			continue;
		}
		vector<string> elems;
		stringstream ss(line);
		string tok;
		while(ss >> tok){
			elems.push_back(tok);
		}
        if(elems.size() <= 2){
            opts::printLog("Line: " + getString<int>(count) + " has incorrect number of columns, in file " + group_file + "\n");
            in.close();

            throw MethodException("Line: " + getString<int>(count) + " has incorrect number of columns, in file " + group_file + "\n");
        }
		string fam = elems.at(0);
		string ind = elems.at(1);
		string group = elems.at(2);
		group = "GROUP_" + group;
		vector<Sample*>::iterator found = find_if(samps->begin(), samps->end(), FindSampleByFamAndID(fam, ind));
		if(found != samps->end()){
			Sample* mysamp = *found;
			map<string, vector<Sample*> >::iterator mfound = groups.find(group);
			if(mfound != groups.end()){
				vector<Sample*> gsamps = groups[group];
				found = find(gsamps.begin(), gsamps.end(), mysamp);
				if(found == gsamps.end()){
					groups[group].push_back(mysamp);
					sample_groups[mysamp] = group;
				}
			}
			else{
				groups[group].push_back(mysamp);
				sample_groups[mysamp] = group;
			}
		}
	}

	map<string, vector<Sample*> >::iterator giter;
	for(giter = groups.begin(); giter != groups.end(); giter++){
		string group = giter->first;
		group_families[group] = generateFamilySet(&(giter->second));
	}

}

void StepOptions::readSampleBprangeFile(){
	ifstream in(sampbprange_file.c_str(), ios::in);
	if(!in){
		throw MethodException("Error opening sample bprange filter file: " + sampbprange_file + "\n");
	}

	int count = 0;
	while(!in.eof()){
		count++;
		string line;
		getline(in, line);
		if(line == ""){
			continue;
		}
		if(line.at(0) == '#'){
			continue;
		}

		vector<string> elems;
		stringstream ss(line);
		string tok;
		while(ss >> tok){
			elems.push_back(tok);
		}
		if(elems.size() != 5){
			in.close();
			throw MethodException("In file: " + sampbprange_file + " expecting 5 columns on line " + getString<int>(count) + "\n");
		}
		string fam = elems.at(0);
		string ind = elems.at(1);
		string chrom_str = elems.at(2);
		string min_str = elems.at(3);
		string max_str = elems.at(4);

		int min = atoi(min_str.c_str());
		int max = atoi(max_str.c_str());

		vector<Marker*> snps;
		Marker* min_mark = new Marker(chrom_str, "min", min);
		Marker* max_mark = new Marker(chrom_str, "max", max);
		snps.push_back(min_mark);
		snps.push_back(max_mark);

		Sample* samp = new Sample();
		samp->setFamID(fam);
		samp->setInd(ind);

		sampbprange_samples.push_back(samp);
		sampbprange_markers.push_back(snps);
	}

	in.close();

}

void StepOptions::readCovariates(vector<Sample*>* samps){
	ifstream in(cov_file.c_str(), ios::in);
	if(!in){
		opts::printLog("Error opening covariate file: " + cov_file + "\n");
		//exit(1);
		throw MethodException("Error opening covariate file: " + cov_file + "\n");
	}
	int count = 0;
	int numcovs = 0;
	vector<int> goodlocs;
	while(!in.eof()){
        count++;
        string line;
        getline(in, line);
        if(line == ""){
            continue;
        }
        if(line.at(0) == '#'){
            continue;
        }

        vector<string> elems;
        stringstream ss(line);
        string tok;
        while(ss >> tok){
           elems.push_back(tok);
        }
        if(elems.size() <= 2){
            opts::printLog("Line: " + getString<int>(count) + " has incorrect number of columns!");
            in.close();
            throw MethodException("Line: " + getString<int>(count) + " has incorrect number of columns!");
        }
        string fam = elems.at(0);
        string ind = elems.at(1);
		bool haveheader = false;
		if(fam == "FamID" && ind == "IndID"){
			haveheader = true;
		}
		if(do_covs_name && !haveheader && count == 1){
			opts::printLog("No headers found.  Please use -covar-number.\n");
			in.close();
			throw MethodException("No headers found.  Please use -covar-number.\n");
		}
        if(numcovs == 0){
            numcovs = elems.size() - 2;
			bool dorange = false;
			string end = "";
            for(int c = 0; c < numcovs; c++){
				if(do_covs_name || do_covs_number){
					if(dorange){
						if(do_covs_number){
							while(end != getString<int>(c+1)){
								if(!haveheader){
									cov_map.push_back("COV"+getString<int>(c+1));
									goodlocs.push_back(c+2);
								}
								else{
									cov_map.push_back(elems.at(c+2));
									goodlocs.push_back(c+2);
								}
								c++;
							}
							if(!haveheader){
								cov_map.push_back("COV"+getString<int>(c+1));
								goodlocs.push_back(c+2);
							}
							else{
								cov_map.push_back(elems.at(c+2));
								goodlocs.push_back(c+2);
							}
						}
						else if(do_covs_name){
							while(end != elems.at(c+2)){
								cov_map.push_back(elems.at(c+2));
								goodlocs.push_back(c+2);
								c++;
							}
							cov_map.push_back(elems.at(c+2));
							goodlocs.push_back(c+2);
						}
						dorange = false;
						end = "";
						continue;
					}
					string fnd = "";
					if(do_covs_number){
						for(int v = 0; v < (int)cov_use.size(); v++){
							string test = getString<int>(c+1);
							if(test == cov_use.at(v)){
								fnd = cov_use.at(v);
								break;
							}
							else if(cov_use.at(v).find('-') != string::npos){
								vector<string> range;
								General::Tokenize(cov_use.at(v), range, "-");
								if(test == range.at(0)){
									fnd = cov_use.at(v);
									break;
								}
							}
						}
					}
					else if(do_covs_name){
						for(int v = 0; v < (int)cov_use.size(); v++){
							string test = elems.at(c+2);
							if(test == cov_use.at(v)){
								fnd = cov_use.at(v);
								break;
							}
							else if(cov_use.at(v).find('-') != string::npos){
								vector<string> range;
								General::Tokenize(cov_use.at(v), range, "-");
								if(test == range.at(0)){
									fnd = cov_use.at(v);
									break;
								}
							}
						}
					}
					if(fnd != ""){
						string cov = fnd;
						if(cov.find('-') != string::npos){
							dorange = true;
							vector<string> range;
							General::Tokenize(cov, range, "-");
							end = range.at(1);
							if(do_covs_number){
								if(!haveheader){
									cov_map.push_back("COV"+getString<int>(c+1));
									goodlocs.push_back(c+2);
								}
								else{
									cov_map.push_back(elems.at(c+2));
									goodlocs.push_back(c+2);
								}
							}
							else if(do_covs_name){
								cov_map.push_back(elems.at(c+2));
								goodlocs.push_back(c+2);
							}
						}
						else{
							if(do_covs_number){
								if(!haveheader){
									cov_map.push_back("COV"+getString<int>(c+1));
									goodlocs.push_back(c+2);
								}
								else{
									cov_map.push_back(elems.at(c+2));
									goodlocs.push_back(c+2);
								}
							}
							else if(do_covs_name){
								cov_map.push_back(elems.at(c+2));
								goodlocs.push_back(c+2);
							}
						}
					}
				}
				else{
					if(!haveheader){
						cov_map.push_back("COV"+getString<int>(c+1));
						goodlocs.push_back(c+2);
					}
					else{
						cov_map.push_back(elems.at(c+2));
						goodlocs.push_back(c+2);
					}
				}
            }
        }
        else if((int)elems.size() != numcovs + 2){
            opts::printLog("Line: " + getString<int>(count) + " has incorrect number of columns!");
            in.close();
            throw MethodException("Line: " + getString<int>(count) + " has incorrect number of columns!");
        }
		if(fam != "FamID" && ind != "IndID"){
            for(int c = 0; c < (int)goodlocs.size(); c++){
                try{
					if(elems.at(goodlocs.at(c)) == getCovarMissing()){
						covs[fam+"#"+ind].push_back(-999999);
					}
					else{
                    	covs[fam+"#"+ind].push_back(atof(elems.at(goodlocs.at(c)).c_str()));
					}
                }catch(...){
                    opts::printLog("Column " + getString<int>(c) +" on line: " + getString<int>(count) + " is not a number!?\n");
                    throw MethodException("Column " + getString<int>(c) +" on line: " + getString<int>(count) + " is not a number!?\n");
                }
            }
        }
    }
    in.close();

}

void StepOptions::readTraits(vector<Sample*>* samps){
	//filter out covariates we don't want
	vector<string> good;
	if(do_traits_name || do_traits_number){
		if(do_traits_name){
			for(int c = 0; c < (int)trait_use.size(); c++){
				string trait = trait_use.at(c);
				if(trait.find('-') != string::npos){
					vector<string> range;
					General::Tokenize(trait, range, "-");
					good.push_back(range.at(0));
					good.push_back(range.at(1));
				}
				else{
					good.push_back(trait);
				}
			}
		}
		else if(do_traits_number){
			for(int c = 0; c < (int)trait_use.size(); c++){
				string trait = trait_use.at(c);
				try{
					if(trait.find('-') != string::npos){
						vector<string> range;
						General::Tokenize(trait, range, "-");
						int loc1 = atoi(range.at(0).c_str());
						int loc2 = atoi(range.at(1).c_str());
						for(int l = loc1; l <= loc2; l++){
							good.push_back(getString<int>(l));
						}
					}
					else{
						good.push_back(trait);
					}
				}catch(...){
					opts::printLog("Error converting " + trait + " to number.\n");
					throw MethodException("Error converting " + trait + " to number.\n");
				}
			}
		}
	}
	ifstream in(trait_file.c_str(), ios::in);
	if(!in){
		opts::printLog("Error opening trait file: " + trait_file + "\n");
		throw MethodException("Error opening trait file: " + trait_file + "\n");
	}
	int count = 0;
	int numtraits = 0;
	vector<int> goodlocs;
	while(!in.eof()){
        count++;
        string line;
        getline(in, line);
        if(line == ""){
            continue;
        }
        if(line.at(0) == '#'){
            continue;
        }

        vector<string> elems;
        stringstream ss(line);
        string tok;
        while(ss >> tok){
           elems.push_back(tok);
        }
        if(elems.size() <= 2){
            opts::printLog("Line: " + getString<int>(count) + " has incorrect number of columns!");
            in.close();
            throw MethodException("Line: " + getString<int>(count) + " has incorrect number of columns!");
        }
        string fam = elems.at(0);
        string ind = elems.at(1);
		bool haveheader = false;
		if(fam == "FamID" && ind == "IndID"){
			haveheader = true;
		}
		if(do_traits_name && !haveheader && count == 1){
			opts::printLog("No headers found.  Please use -traits-number.\n");
			in.close();
			throw MethodException("No headers found.  Please use -traits-number.\n");
		}
        if(numtraits == 0){
            numtraits = elems.size() - 2;
			bool dorange = false;
			string end = "";
            for(int c = 0; c < numtraits; c++){
				if(do_traits_name || do_traits_number){
					if(dorange){
						if(do_traits_number){
							while(end != getString<int>(c+1)){
								if(!haveheader){
									trait_map.push_back("TRAIT"+getString<int>(c+1));
									goodlocs.push_back(c+2);
								}
								else{
									trait_map.push_back(elems.at(c+2));
									goodlocs.push_back(c+2);
								}
								c++;
							}
							if(!haveheader){
								trait_map.push_back("TRAIT"+getString<int>(c+1));
								goodlocs.push_back(c+2);
							}
							else{
								trait_map.push_back(elems.at(c+2));
								goodlocs.push_back(c+2);
							}
						}
						else if(do_traits_name){
							while(end != elems.at(c+2)){
								trait_map.push_back(elems.at(c+2));
								goodlocs.push_back(c+2);
								c++;
							}
							trait_map.push_back(elems.at(c+2));
							goodlocs.push_back(c+2);
						}
						dorange = false;
						end = "";
						continue;
					}
					string fnd = "";
					if(do_traits_number){
						for(int v = 0; v < (int)trait_use.size(); v++){
							string test = getString<int>(c+1);
							if(test == trait_use.at(v)){
								fnd = trait_use.at(v);
								break;
							}
							else if(trait_use.at(v).find('-') != string::npos){
								vector<string> range;
								General::Tokenize(trait_use.at(v), range, "-");
								if(test == range.at(0)){
									fnd = trait_use.at(v);
									break;
								}
							}
						}
					}
					else if(do_traits_name){
						for(int v = 0; v < (int)trait_use.size(); v++){
							string test = elems.at(c+2);
							if(test == trait_use.at(v)){
								fnd = trait_use.at(v);
								break;
							}
							else if(trait_use.at(v).find('-') != string::npos){
								vector<string> range;
								General::Tokenize(trait_use.at(v), range, "-");
								if(test == range.at(0)){
									fnd = trait_use.at(v);
									break;
								}
							}
						}
					}
					if(fnd != ""){
						string trait = fnd;
						if(trait.find('-') != string::npos){
							dorange = true;
							vector<string> range;
							General::Tokenize(trait, range, "-");
							end = range.at(1);
							if(do_traits_number){
								if(!haveheader){
									trait_map.push_back("TRAIT"+getString<int>(c+1));
									goodlocs.push_back(c+2);
								}
								else{
									trait_map.push_back(elems.at(c+2));
									goodlocs.push_back(c+2);
								}
							}
							else if(do_traits_name){
								trait_map.push_back(elems.at(c+2));
								goodlocs.push_back(c+2);
							}
						}
						else{
							if(do_traits_number){
								if(!haveheader){
									trait_map.push_back("TRAIT"+getString<int>(c+1));
									goodlocs.push_back(c+2);
								}
								else{
									trait_map.push_back(elems.at(c+2));
									goodlocs.push_back(c+2);
								}
							}
							else if(do_traits_name){
								trait_map.push_back(elems.at(c+2));
								goodlocs.push_back(c+2);
							}
						}
					}
				}
				else{
					if(!haveheader){
						trait_map.push_back("TRAIT"+getString<int>(c+1));
						goodlocs.push_back(c+2);
					}
					else{
						trait_map.push_back(elems.at(c+2));
						goodlocs.push_back(c+2);
					}
				}
            }
        }
        else if((int)elems.size() != numtraits + 2){
            opts::printLog("Line: " + getString<int>(count) + " has incorrect number of columns!");
            in.close();
            throw MethodException("Line: " + getString<int>(count) + " has incorrect number of columns!");
        }
		if(fam != "FamID" && ind != "IndID"){
            for(int c = 0; c < (int)goodlocs.size(); c++){
                try{
					if(elems.at(goodlocs.at(c)) == getTraitMissing()){
						traits[fam+"#"+ind].push_back(-999999);
					}
					else{
                    	traits[fam+"#"+ind].push_back(atof(elems.at(goodlocs.at(c)).c_str()));
					}
                }catch(...){
                    opts::printLog("Column " + getString<int>(c) +" on line: " + getString<int>(count) + " is not a number!?\n");
                    throw MethodException("Column " + getString<int>(c) +" on line: " + getString<int>(count) + " is not a number!?\n");
                }
            }
        }
    }
    in.close();
}

vector<Family*> StepOptions::generateFamilySet(vector<Sample*>* samps){
    vector<Family*> fams;

    vector<Family*>::iterator famiter;

    for(int i = 0; i < (int)samps->size(); i++){
        Sample* samp = (*samps).at(i);
        famiter = find_if(fams.begin(), fams.end(), FindFamily(samp->getFamID()));
        if(famiter != fams.end()){
            (*famiter)->AddInd(samp);
        }
        else{
            Family* fam = new Family();
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
            fam->setEnabled(true);
            fams.push_back(fam);
            fam->setLoc((fams.size() - 1));
        }
    }

	return fams;
}

string StepOptions::getAlleleMapping(string s){
	map<string, string>::iterator it = custom_alleles.find(s);
	if(it != custom_alleles.end()){
		string r = it->second;
		return r;
	}
	return s;
}

void StepOptions::performTransforms(DataSet* ds){
	vector<int> transform_index;
	for(int i = 0; i < (int)to_transform.size(); i++){
		if(ds->num_covariates() > 0){
		int index = ds->get_covariate_index(to_transform.at(i));
		if(index > 0){
			transform_index.push_back(index);
		}
		}
	}
	if(transform_index.size() > 0){
	if(sqrt_transform){
		SquareRootTransform srt;
		if(transform_index.size() == 0){
			srt.TransformData(ds);
		}
		else{
			for(int i = 0; i < (int)transform_index.size(); i++){
				srt.TransformCovar(ds, transform_index.at(i));
			}
		}
	}
	if(boxcox_transform){
		BoxCoxTransform bct;
		if(transform_index.size() == 0){
			bct.TransformData(ds);
		}
		else{
			for(int i = 0; i < (int)transform_index.size(); i++){
				bct.TransformCovar(ds, transform_index.at(i));
			}
		}
		boxcox_covar_lambdas = bct.GetCovarLambda();
		boxcox_trait_lambdas = bct.GetTraitLambda();
	}
	if(log_transform){
		LogTransform lt;
		if(transform_index.size() == 0){
			lt.TransformData(ds);
		}
		else{
			for(int i = 0; i < (int)transform_index.size(); i++){
				lt.TransformCovar(ds, transform_index.at(i));
			}
		}
	}
	}
}

void StepOptions::undoTransforms(DataSet* ds){
	vector<int> transform_index;
	for(int i = 0; i < (int)to_transform.size(); i++){
		if(ds->num_covariates() > 0){
		int index = ds->get_covariate_index(to_transform.at(i));
		if(index > 0){
			transform_index.push_back(index);
		}
		}
	}
	if(transform_index.size() > 0){
	if(log_transform){
		LogTransform lt;
		if(transform_index.size() == 0){
			lt.UndoTransform(ds);
		}
		else{
			for(int i = 0; i < (int)transform_index.size(); i++){
				lt.UndoCovariate(ds, transform_index.at(i));
			}
		}
	}
	if(boxcox_transform){
		BoxCoxTransform bct;
		bct.SetLambdas(boxcox_covar_lambdas, boxcox_trait_lambdas);
		if(transform_index.size() == 0){
			bct.UndoTransform(ds);
		}
		else{
			for(int i = 0; i < (int)transform_index.size(); i++){
				bct.UndoCovariate(ds, transform_index.at(i));
			}
		}
	}
	if(sqrt_transform){
		SquareRootTransform srt;
		if(transform_index.size() == 0){
			srt.UndoTransform(ds);
		}
		else{
			for(int i = 0; i < (int)transform_index.size(); i++){
				srt.UndoCovariate(ds, transform_index.at(i));
			}
		}
	}
	}
}

void StepOptions::parseLinRConditionList(vector<Marker*>* loci){
	vector<string> tokens;
	General::Tokenize(getLinRConditionString(), tokens, ",");
	for(int i = 0; i < (int)tokens.size(); i++){
		vector<Marker*>::iterator iter = find_if(loci->begin(), loci->end(), FindMarker(tokens.at(i)));
		if(iter != loci->end()){
			linr_condition_list.push_back(iter - loci->begin());
		}
	}
}

void StepOptions::readIBSPairsFile(string file){
	ifstream in(file.c_str(), ios::in);
	if(!in){
		opts::printLog("Error opening IBS Pairs file: " + file + "\n");
		throw MethodException("Error opening IBS Pairs file: " + file + "\n");
	}

	int count = 0;
	while(!in.eof()){
		count++;
		string line;
		getline(in, line);
		if(line == ""){
			continue;
		}
		if(line.at(0) == '#'){
			continue;
		}
		vector<string> elems;
		stringstream ss(line);
		string tok;
		while(ss >> tok){
			elems.push_back(tok);
		}
		if(elems.size() != 4){
			opts::printLog("Line:: " + getString<int>(count) + " has incorrect number of columns.  Expecting 4!");
			in.close();
			throw MethodException("Line:: " + getString<int>(count) + " has incorrect number of columns.  Expecting 4!");
		}
		string samp1 = elems.at(0) + "\t" + elems.at(1);
		string samp2 = elems.at(2) + "\t" + elems.at(3);
		ibs_pairs[samp1].push_back(samp2);
	}
	if(in.is_open()){
		in.close();
	}
}

void StepOptions::readIBSTrioPairsFile(string file){
	ifstream in(file.c_str(), ios::in);
	if(!in){
		opts::printLog("Error opening IBS Pairs file: " + file + "\n");
		throw MethodException("Error opening IBS Pairs file: " + file + "\n");
	}

	int count = 0;
	while(!in.eof()){
		count++;
		string line;
		getline(in, line);
		if(line == ""){
			continue;
		}
		if(line.at(0) == '#'){
			continue;
		}
		vector<string> elems;
		stringstream ss(line);
		string tok;
		while(ss >> tok){
			elems.push_back(tok);
		}
		if(elems.size() != 2){
			opts::printLog("Line:: " + getString<int>(count) + " has incorrect number of columns.  Expecting 2!");
			in.close();
			throw MethodException("Line:: " + getString<int>(count) + " has incorrect number of columns.  Expecting 2!");
		}
		string fam1 = elems.at(0);
		string fam2 = elems.at(1);
		ibs_trio_pairs[fam1].push_back(fam2);
	}
	if(in.is_open()){
		in.close();
	}
}

void StepOptions::readBioTextFile(string file){
	int line_length = 30;
	bio_pairs.clear();
	ifstream in(file.c_str(), ios::in);
	if(!in){
		opts::printLog("Error opening Biofilter file: " + file + "\n");
		throw MethodException("Error opening Biofilter file: " + file + "\n");
	}

	int count = 1;

	string header;
	getline(in, header);
	if(getBioOffsetBegin() > 0){
		count = getBioOffsetBegin();
		in.seekg(((getBioOffsetBegin() - 1) * line_length), ios::cur);
	}

	while(!in.eof()){
		string line;
		getline(in, line);
		if(line == ""){
			continue;
		}
		if(line.at(0) == '#'){
			continue;
		}
		vector<string> elems;
		stringstream ss(line);
		string tok;
		while(ss >> tok){
			elems.push_back(tok);
		}
		if(elems.size() < 2){
			opts::printLog("Line:: " + getString<int>(count) + " has incorrect number of columns.  Expecting at least 2!\n");
			in.close();
			throw MethodException("Line:: " + getString<int>(count) + " has incorrect number of columns.  Expecting at least 2!\n");
		}
		if(elems.at(0).find("rs", 0, 2) == string::npos){
			elems.at(0) = "rs" + elems.at(0);
		}
		if(elems.at(1).find("rs", 0, 2) == string::npos){
			elems.at(1) = "rs" + elems.at(1);
		}
		string snp1 = elems.at(0);
		string snp2 = elems.at(1);

		bio_pairs[snp1].push_back(snp2);

		if(getBioOffsetEnd() > 0){
			if(count >= getBioOffsetEnd()){
				break;
			}
		}
		count++;
	}
	if(in.is_open()){
		in.close();
	}

}

void StepOptions::readLinRConditionFile(vector<Marker*>* loci){
	ifstream in(linr_condition_file.c_str(), ios::in);
	if(!in){
		opts::printLog("Error opening condition file: " + linr_condition_file + "\n");
		throw MethodException("Error opening condition file: " + linr_condition_file + "\n");
	}
	int count = 0;
	vector<int> goodlocs;
	while(!in.eof()){
        count++;
        string line;
        getline(in, line);
        if(line == ""){
            continue;
        }
        if(line.at(0) == '#'){
            continue;
        }

        vector<string> elems;
        stringstream ss(line);
        string tok;
        while(ss >> tok){
           elems.push_back(tok);
        }
        if(elems.size() != 1){
            opts::printLog("Line: " + getString<int>(count) + " has incorrect number of columns!");
            in.close();
            throw MethodException("Line: " + getString<int>(count) + " has incorrect number of columns!");
        }

        string rsid = elems.at(0);
        vector<Marker*>::iterator iter = find_if(loci->begin(), loci->end(), FindMarker(rsid));
        if(iter != loci->end()){
        	linr_condition_list.push_back(iter - loci->begin());
        }
	}
	in.close();
}

/* file format is SNP GRP FREQ
 * option -group-freq filename
 *
 *
 */
void StepOptions::readGroupFrequencies(){
	opts::printLog("Reading marker/group frequency file: " + group_frequency_file + "\n");
	ifstream finput;
	finput.open(group_frequency_file.c_str(), ios::in);

	if(!finput){
		throw MethodException("Error opening marker/group frequency file: " + group_frequency_file + "\n");
	}

	string line = "";

	int freqline = 1;
	while(getline(finput, line)){
		if(line.size() == 0){
			continue;
		}
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 3){
			throw MethodException("Marker frequency file column size != 3 on line: " + line + "\n");
		}
		if(tokens.at(0).at(0) == '#'){
			continue;
		}

		float freq = -1.0f;
		istringstream b(tokens.at(2));
		if(!(b >> freq)){
			throw MethodException(tokens.at(1) + " is not a valid number on line: " + getString<int>(freqline));
		}
		group_frequencies[tokens.at(0) + " GROUP_" + tokens.at(1)] = freq;
		freqline++;
	}

	if(finput){
		finput.close();
	}

}

}
