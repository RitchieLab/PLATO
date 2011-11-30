/**********************************************************************************
*                       Marker Genotype Efficiency Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all genotypes and generates a genotype efficiency for all markers.
*
*
* Files generated:
*	percent_breakdown_by_marker.txt
*	percent_breakdown_by_chrom.txt
*       post_marker_geno_eff_filter_summary.txt
*
*File: MDROutput.cc
**********************************************************************************/


#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include "MDROutput.h"
#include "General.h"
#include "Helpers.h"

namespace Methods{
void MDROutput::FilterSummary(){
}

void MDROutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers).at(i)->setFlag(false);
	}
}

void MDROutput::filter(){
}

void MDROutput::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
	orig_num_markers = markers->size();
	orig_num_families = families->size();
 	orig_num_individuals = samples->size();

   	int ssize = samples->size();
	int msize = markers->size();
	string fname1 = opts::_OUTPREFIX_ + "input_mdr" + options.getOut() + ".txt";
    if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".txt";
	}

	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	filenames.push_back(fname1);
	ofstream mdr (fname1.c_str());
	if(!mdr.is_open()){
		opts::printLog("Unable to open "+fname1+" for output!\n");
		throw MethodException("Unable to open "+fname1+" for output!\n");
	}
	string fname2 = opts::_OUTPREFIX_ + "input_mdr" + options.getOut() + ".log";
    if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".log";
	}
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	filenames.push_back(fname2);
	ofstream mdr_log (fname2.c_str());
	if(!mdr_log.is_open()){
		opts::printLog("Unable to open "+fname2+" for output!\n");
		throw MethodException("Unable to open "+fname2+" for output!\n");
	}
	string fname3 = opts::_OUTPREFIX_ + "input_mdr" + options.getOut() + ".map";
    if(options.getOverrideOut().size() > 0){
		fname3 = options.getOverrideOut() + ".map";
	}
	if(!overwrite){
		fname3 += "." + getString<int>(order);
	}
	filenames.push_back(fname3);
	ofstream mdr_map (fname3.c_str());
	if(!mdr_map.is_open()){
		opts::printLog("Unable to open "+fname3+" for output!\n");
		throw MethodException("Unable to open "+fname3+" for output!\n");
	}
	int prev_base = 0;
	int prev_chrom = -1;
	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	msize = good_markers.size();
	mdr_log << "Excluding the following variations:" << endl;
	for(int m = 0; m < msize; m++){
		Marker* mark = good_markers.at(m);
		if(mark->isEnabled()){
			if(mark->isMicroSat()){
				mdr_log << mark->toString() << "\tMore than 2 alleles." << endl;
				continue;
			}
			mdr_map << mark->getChrom() << "\t" << mark->getRSID() << "\t" << mark->getBPLOC() << "\t" << mark->getAllele1() << "\t" << mark->getAllele2() << endl;
		}
		else{
			mdr_log << mark->toString() << "\tDisabled or doesn't meet chrom/bploc restrictions or is disabled" << endl;
		}
	}
	prev_base = 0;
	prev_chrom = -1;
	mdr_log << endl;
	mdr_log << "Excluding the following samples:" << endl;
	int cases = 0;
	int controls = 0;
	vector<bool> samp_flags(ssize, false);
	
  int phenoloc = -9;
	if(options.getUsePheno()){
		if(options.getPhenoName() == ""){
			phenoloc = options.getPhenoLoc();
		}
		else{
			phenoloc = data_set->get_trait_index(options.getPhenoName());
		}
	}
	
	if(!options.get_no_rules()){
	for(int i = 0; i < ssize; i++){
		Sample* samp = (*samples).at(i);
		if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples())){
		
		  double pheno;
		  if(options.getUsePheno()){
    		pheno = samp->getPheno(phenoloc);
  	  }
  		else{
  		  pheno = samp->getPheno();
  		}
		
			if(pheno != 2 && pheno != 1){
				mdr_log << samp->toString() << "\tNot affected or unaffected." << endl;
				samp_flags.at(samp->getLoc()) = true;
				continue;
			}
			else if(pheno == 2){
				cases++;
			}
			else if(pheno == 1){
				controls++;
			}
			for(int m = 0; m < msize; m++){
				Marker* mark = good_markers.at(m);
				if(mark->isEnabled()){
					if(mark->isMicroSat()){
						continue;
					}
					if(!mark->isMicroSat()){
					}
				}
			}
		}
	}
	}
	bool alert = false;
	if(cases < 100 && !options.get_no_rules()){
		opts::printLog("ALERT: Number of cases is < 100: " + getString<int>(cases) + " found.");
		alert = true;
	}
	if(controls < 100 && !options.get_no_rules()){
		opts::printLog("ALERT: Number of controls is < 100: " + getString<int>(controls) + " found.");
		alert = true;
	}
	if(alert && !options.get_no_rules()){
		opts::printLog("ALERT: Due to the low number of cases/controls, MDR may not work properly with the input files generated here.");
	}

	if(options.getMDRGuiOutput()){
		if(options.getMDRPedigreeOutput()){
			mdr << "FamID\tIndID\tFather\tMother\tGender\tAff\t";
		}
		for(int m = 0; m < msize; m++){
			Marker* mark = good_markers.at(m);
			if(mark->isEnabled()){
				if(mark->isMicroSat()){
					continue;
				}
				mdr << mark->getRSID() << "\t";
			}
		}
		if(!options.getMDRPedigreeOutput()){
			mdr << "AFF";
		}
		mdr << "\n";
	}


	for(int i = 0; i < ssize; i++){
		Sample* samp = (*samples).at(i);
		if((samp->isEnabled() ||(samp->isExcluded() && options.doIncExcludedSamples())) && !samp_flags.at(samp->getLoc())){
		
		  double pheno;
	    if(options.getUsePheno()){
    		pheno = (samp->getPheno(phenoloc) - 1);
  	  }
  		else{
  		  pheno = (samp->getPheno() - 1);
  		}
		
			if(options.getMDRPedigreeOutput()){
				mdr << samp->getFamID() << "\t";
				mdr << samp->getInd() << "\t";
				mdr << samp->getDadID() << "\t";
				mdr << samp->getMomID() << "\t";
				if(samp->getSex()){
					mdr << "1\t";
				}
				else{
					mdr << "2\t";
				}
				// add one because pedigree format uses 1 and 2 for affection status
				mdr << pheno+1 << "\t";
// 				if(options.getUsePheno()){
// 					mdr << samp->getPheno(phenoloc) << "\t";
// 				}
// 				else{
// 					mdr << samp->getPheno() << "\t";
// 				}
			}
			if(!options.getMDRGuiOutput() && !options.get_no_rules() && !options.getMDRPedigreeOutput()){
			  mdr << pheno << "\t";
// 			  if(options.getUsePheno()){
//           if(samp->getPheno(phenoloc)==comissing){
//             opts::printLog("ALERT Skipping individual with missing phenotype: " + samp->getFamIDOrig() + " " + samp->getIndOrig() +"\n");
//             continue;
//           }
//           else
//     				mdr << (samp->getPheno(phenoloc) - 1) << "\t";
//   		  }
//   		  else{
//   		    mdr << (samp->getPheno() - 1) << "\t";
//   		  }
			}
			else if(!options.getMDRGuiOutput() && options.get_no_rules() && !options.getMDRPedigreeOutput()){
			  mdr << pheno << "\t";
// 				if(options.getUsePheno()){
// 					mdr << samp->getPheno(phenoloc) << "\t";
// 				}
// 				else{
// 					mdr << (samp->getPheno()) << "\t";
// 				}
			}
			bool first = true;
			for(int m = 0; m < msize; m++){
				Marker* mark = good_markers.at(m);
				if(mark->isEnabled()){
					if(mark->isMicroSat()){
						continue;
					}
					int m_loc = mark->getLoc();
					if(!mark->isMicroSat()){
						if(!first){
							mdr << "\t";
						}
						if(!samp->getAone(m_loc) && !samp->getAtwo(m_loc)){
							mdr << "0";
						}
						else if(!samp->getAone(m_loc) && samp->getAtwo(m_loc)){
							mdr << "1";
						}
						else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && !samp->getAmissing(m_loc)){
							mdr << "2";
						}
						else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && samp->getAmissing(m_loc)){
							mdr << "-1";
						}
						first = false;
					}
				}
			}
			if(options.getMDRGuiOutput() && !options.getMDRPedigreeOutput()){
				mdr << "\t" << (samp->getPheno() - 1);
			}
			mdr << endl;
		}
		else if(samp->isFlagged()){
			samp->setFlag(false);
		}
	}
	mdr.close();
	mdr_log.close();
	mdr_map.close();
}


int MDROutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

}
