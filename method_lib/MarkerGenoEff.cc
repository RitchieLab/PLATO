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
*
*File: MarkerGenoEff.cc
**********************************************************************************/


#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include "MarkerGenoEff.h"
#include "Options.h"
#include "General.h"
#include "Helpers.h"
namespace Methods{
string MarkerGenoEff::stepname = "marker-geno-eff";

//DEPRECATED
void MarkerGenoEff::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

//DEPRECATED
void MarkerGenoEff::PrintSummary(){
	string filename = opts::_OUTPREFIX_ + "marker_geno_eff" + options.getOut() + ".txt";
	string filenameg = opts::_OUTPREFIX_ + "marker_geno_eff_groups" + options.getOut() + ".txt";
	string filenamegm = opts::_OUTPREFIX_ + "marker_geno_eff_groups_missing" + options.getOut() + ".txt";
	if(!overwrite){
		filename += "." + getString<int>(order);
		filenameg += "." + getString<int>(order);
	}

	ofstream bymarker (filename.c_str());
	if(!bymarker.is_open()){
		opts::printLog("Unable to open " + filename + "\n");
		throw MethodException("Unable to open " + filename + "\n");
	}

	ofstream bygroup;
	ofstream bygroupmissing;
	if(options.doGroupFile()){
		bygroup.open(filenameg.c_str());
		if(!bygroup.is_open()){
			throw MethodException("Unable to open " + filenameg + "\n");
		}
		bygroup.precision(4);
		opts::addFile("Marker", stepname, filenameg);

		bygroupmissing.open(filenamegm.c_str());
		if(!bygroupmissing.is_open()){
			throw MethodException("Unable to open " + filenamegm + "\n");
		}
		opts::addFile("Batch", stepname, filenamegm);
		bygroupmissing.precision(4);

		if((*markers).at(0)->getDetailHeaders().size() > 0){
			bygroup << "Chrom\trsid\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders();
		}
		else{
			bygroup << "Chrom\trsid\tProbeID\tbploc";
		}
		map<string, vector<Sample*> > groups = options.getGroups();
		map<string, vector<Sample*> >::iterator giter;
		for(giter = groups.begin(); giter != groups.end(); giter++){
			string group = giter->first;
			opts::addHeader(filenameg, "%" + group);
			opts::addHeader(filenameg, group + "_Zero_Count");
			opts::addHeader(filenameg, group + "_Total_Used");
			bygroup << "\t" << "%" << group;
			bygroup << "\t" << group << "_Zero_Count";
			bygroup << "\t" << group << "_Total_Used";
		}
		bygroup << endl;

		bygroupmissing << "Batch\tMiss_avg\tLog_miss\tn" << endl;
	}

	opts::addFile("Marker", stepname, filename);
	bymarker.precision(4);
	if((*markers).at(0)->getDetailHeaders().size() > 0){
		bymarker << "Chrom\trsid\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders() << "\t%GenoEff_Ind_All\tInd_Zero_Count\tTotal_Individuals_Used\t%GenoEff_Ind_Cases\tInd_Zero_Count_Cases\tTotal_Individuals_Used_Cases\t%GenoEff_Ind_Controls\tInd_Zero_Count_Controls\tTotal_Individuals_Used_Controls";
	}
	else{
		bymarker << "Chrom\trsid\tProbeID\tbploc\t%GenoEff_Ind_All\tInd_Zero_Count\tTotal_Individuals_Used\t%GenoEff_Ind_Cases\tInd_Zero_Count_Cases\tTotal_Individuals_Used_Cases\t%GenoEff_Ind_Controls\tInd_Zero_Count_Controls\tTotal_Individuals_Used_Controls";
	}
	opts::addHeader(filename, "%GenoEff_Ind_All");
	opts::addHeader(filename, "Ind_Zero_Count");
	opts::addHeader(filename, "Total_Individuals_Used");
	opts::addHeader(filename, "%GenoEff_Ind_Cases");
	opts::addHeader(filename, "Ind_Zero_Count_Cases");
	opts::addHeader(filename, "Total_Individuals_Used_Cases");
	opts::addHeader(filename, "%GenoEff_Ind_Controls");
	opts::addHeader(filename, "Ind_Zero_Count_Controls");
	opts::addHeader(filename, "Total_Individuals_Used_Controls");
	bymarker << endl;


	int size = markers->size();
	for(int i = 0; i < size; i++){
		if((*markers).at(i)->isEnabled() && !(*markers).at(i)->isFlagged()){

			float percent = 0.0f;
			float caseper = 0.0f;
			float contper = 0.0f;
			if(total.at(i) > 0){
				percent = (1.0f - ((float)zeros.at(i)/(float)total.at(i))) * 100.0f;
			}
			if(casetotal.at(i) > 0){
				caseper = (1.0f - ((float)casezeros.at(i)/(float)casetotal.at(i))) * 100.0f;
			}
			if(controltotal.at(i) > 0){
				contper = (1.0f - ((float)controlzeros.at(i)/(float)controltotal.at(i))) * 100.0f;
			}

			bymarker << (*markers).at(i)->toString() << "\t"
					 << percent << "\t"
					 << zeros.at(i) << "\t"
					 << total.at(i) << "\t"
					 << caseper << "\t"
					 << casezeros.at(i) << "\t"
					 << casetotal.at(i) << "\t"
					 << contper << "\t"
					 << controlzeros.at(i) << "\t"
					 << controltotal.at(i);
			bymarker << endl;
			if(options.doGroupFile()){
				bygroup << (*markers).at(i)->toString();
				map<string, vector<Sample*> > groups = options.getGroups();
				map<string, vector<Sample*> >::iterator giter;
				for(giter = groups.begin(); giter != groups.end(); giter++){
					string group = giter->first;
					int gzero = groupzeros.at(group).at(i);
					int gtotal = grouptotal.at(group).at(i);
					float gper = 0.0f;
					if(gtotal > 0){
						gper = (float)(1.0f - ((float)gzero/(float)gtotal)) * 100.0f;
					}
					bygroup << "\t" << gper;
					bygroup << "\t" << gzero;
					bygroup << "\t" << gtotal;
				}
				bygroup << endl;
			}
		}
		(*markers).at(i)->setFlag(false);
	}

	if(bymarker.is_open()){
		bymarker.close();
	}

}

//DEPRECATED
void MarkerGenoEff::filterOne(int m){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
			if((*markers).at(m)->isEnabled() && !(*markers).at(m)->isFlagged()){
				if(options.doChrom()){
					if(!options.checkChrom((*markers).at(m)->getChrom())){
				    	return;
					}
				    if(!options.checkBp((*markers).at(m)->getBPLOC())){
				    	return;
					}
				}

				double percent = 0.0f;
				double caseper = 0.0f;
				double contper = 0.0f;
				bool inc = false;
				if(total_one > 0){
					percent = getPercent();
				}
				if(casetotal_one > 0){
					caseper = getCasePercent();
				}
				if(controltotal_one > 0){
					contper = getControlPercent();
				}

				if(options.doThreshMarkersLow() && Helpers::dLess(percent, options.getThreshMarkersLow())){
					(*markers).at(m)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && Helpers::dGreater(percent, options.getThreshMarkersHigh())){
					(*markers).at(m)->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}

	}
}

//DEPRECATED
void MarkerGenoEff::filter(){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int size = markers->size();
		for(int i = 0; i < size; i++){
			if((*markers).at(i)->isEnabled() && !(*markers).at(i)->isFlagged()){
				if(options.doChrom()){
					if(!options.checkChrom((*markers).at(i)->getChrom())){
					    continue;
				    }
				    if(!options.checkBp((*markers).at(i)->getBPLOC())){
					    continue;
				    }
				}

				double percent = 0.0f;
				float caseper = 0.0f;
				float contper = 0.0f;
				bool inc = false;
				if(total.at(i) > 0){
					percent = (double)((1.0f - ((double)zeros.at(i)/(double)total.at(i))) * 100.0f);
				}
				if(casetotal.at(i) > 0){
					caseper = (1.0f - ((float)casezeros.at(i)/(float)casetotal.at(i))) * 100.0f;
				}
				if(controltotal.at(i) > 0){
					contper = (1.0f - ((float)controlzeros.at(i)/(float)controltotal.at(i))) * 100.0f;
				}

				if(options.doThreshMarkersLow() && Helpers::dLess(percent, options.getThreshMarkersLow())){
					(*markers).at(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && Helpers::dGreater(percent, options.getThreshMarkersHigh())){
					(*markers).at(i)->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}
		}
	}
}

//main calculation method.  depends on calling object to iterate over markers
void MarkerGenoEff::calcOne(Marker* mark){
 	zeros_one = 0;
	total_one = 0;
	casezeros_one = 0;
	casetotal_one = 0;
	controlzeros_one = 0;
	controltotal_one = 0;
	vector<Sample*>::iterator s_iter;

	if(options.doGroupFile()){
		map<string, vector<Sample*> >::iterator giter;
		map<string, vector<Sample*> > groups = options.getGroups();

		for(giter = groups.begin(); giter != groups.end(); giter++){
			string mygroup = giter->first;
			vector<Sample*> mysamps = giter->second;
			groupzeros_one[mygroup]=0;
			grouptotal_one[mygroup]=0;

			for(int s = 0; s < (int)mysamps.size(); s++){
				Sample* samp = mysamps.at(s);
				if(samp->isEnabled()){
					if(mark->isEnabled()){
						int loc = mark->getLoc();
						if(!mark->isMicroSat()){
							if(samp->getAone(loc) && samp->getAtwo(loc) && samp->getAmissing(loc)){
								groupzeros_one.at(mygroup)++;
							}
						}
						else{
							if(samp->getAbone(loc) == -1){
								groupzeros_one.at(mygroup)++;
							}
						}
						grouptotal_one.at(mygroup)++;
					}
				}
			}
		}

	}

	for(s_iter = samples->begin(); s_iter != samples->end(); s_iter++){
		if((*s_iter)->isEnabled()){
				int loc = mark->getLoc();
				if(mark->isEnabled()){
					if(!mark->isMicroSat()){
						if((*s_iter)->getAone(loc) && (*s_iter)->getAtwo(loc) && (*s_iter)->getAmissing(loc)){
							zeros_one++;
							if((*s_iter)->getPheno() == 2){
								casezeros_one++;
							}
							else if((*s_iter)->getPheno() == 1){
								controlzeros_one++;
							}
						}
					}
					else{
						if((*s_iter)->getAbone(loc) == -1){
							zeros_one++;
							if((*s_iter)->getPheno() == 2){
								casezeros_one++;
							}
							else if((*s_iter)->getPheno() == 1){
								controlzeros_one++;
							}
						}
					}
					total_one++;
					if((*s_iter)->getPheno() == 2){
						casetotal_one++;
					}
					else if((*s_iter)->getPheno() == 1){
						controltotal_one++;
					}
				}
		}
	}

}

//DEPRECATED
void MarkerGenoEff::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
 	zeros.resize(markers->size());
	total.resize(markers->size());
	casezeros.resize(markers->size());
	casetotal.resize(markers->size());
	controlzeros.resize(markers->size());
	controltotal.resize(markers->size());
	vector<Sample*>::iterator s_iter;

	if(options.doGroupFile()){
		options.readGroups(samples);
		map<string, vector<Sample*> >::iterator giter;
		map<string, vector<Sample*> > groups = options.getGroups();

		for(giter = groups.begin(); giter != groups.end(); giter++){
			string mygroup = giter->first;
			vector<Sample*> mysamps = giter->second;
			groupzeros.at(mygroup).resize(markers->size());
			grouptotal.at(mygroup).resize(markers->size());

			for(int s = 0; s < (int)mysamps.size(); s++){
				Sample* samp = mysamps.at(s);
				if(samp->isEnabled()){
					int prev_base = 0;
					int prev_chrom = -1;
					for(int m = 0; m < (int)markers->size(); m++){
						Marker* mark = (*markers).at(m);
						if(Helpers::isValidMarker(mark, &options, prev_base, prev_chrom)){
							int loc = mark->getLoc();
							if(!mark->isMicroSat()){
								if(samp->getAone(loc) && !samp->getAtwo(loc)){
									groupzeros.at(mygroup).at(m)++;
								}
							}
							else{
								if(samp->getAbone(loc) == -1){
									groupzeros.at(mygroup).at(m)++;
								}
							}
							grouptotal.at(mygroup).at(m)++;
						}
					}
				}
			}
		}

	}

	for(s_iter = samples->begin(); s_iter != samples->end(); s_iter++){
		if((*s_iter)->isEnabled()){
			int end = markers->size();
			int prev_base = 0;
			int prev_chrom = -1;
			for(int i = 0; i < end; i++){
				int loc = (*markers).at(i)->getLoc();
				if(Helpers::isValidMarker((*markers).at(i), &options, prev_base, prev_chrom)){

					if(!(*markers).at(i)->isMicroSat()){
						if((*s_iter)->getAone(loc) && !(*s_iter)->getAtwo(loc)){
							zeros.at(i)++;
							if((*s_iter)->getPheno() == 2){
								casezeros.at(i)++;
							}
							else if((*s_iter)->getPheno() == 1){
								controlzeros.at(i)++;
							}
						}
					}
					else{
						if((*s_iter)->getAbone(loc) == -1){
							zeros.at(i)++;
							if((*s_iter)->getPheno() == 2){
								casezeros.at(i)++;
							}
							else if((*s_iter)->getPheno() == 1){
								controlzeros.at(i)++;
							}
						}
					}
					total.at(i)++;
					if((*s_iter)->getPheno() == 2){
						casetotal.at(i)++;
					}
					else if((*s_iter)->getPheno() == 1){
						controltotal.at(i)++;
					}
				}
			}
		}
	}

}
}
