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
*File: ProcessMarkerGenoEff.cc
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
#include "ProcessMarkerGenoEff.h"
#include <Options.h>
#include <General.h>
#include <Helper.h>
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
using namespace Methods;
string ProcessMarkerGenoEff::stepname = "marker-geno-eff";

void ProcessMarkerGenoEff::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessMarkerGenoEff::PrintSummary(){
	string filename = opts::_OUTPREFIX_ + "marker_geno_eff" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	string filenameg = opts::_OUTPREFIX_ + "marker_geno_eff_groups" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	string filenamegm = opts::_OUTPREFIX_ + "marker_geno_eff_groups_missing" + options.getOut() + ".txt";
	if(!overwrite){
		filename += "." + getString<int>(order);
		filenameg += "." + getString<int>(order);
	}

	ofstream bymarker (filename.c_str());
	if(!bymarker.is_open()){
		opts::printLog("Unable to open " + filename + "\n");
		throw MethodException("");
	}

	ofstream bygroup;
	ofstream bygroupmissing;
	if(options.doGroupFile()){
		bygroup.open(filenameg.c_str());
		if(!bygroup.is_open()){
			opts::printLog("Unable to open " + filenameg + "\n");
			throw MethodException("");
		}
		opts::addFile("Marker", stepname, filenameg);
		bygroup.precision(4);

		bygroupmissing.open(filenamegm.c_str());
		if(!bygroupmissing.is_open()){
			throw MethodException("Unable to open " + filenamegm + "\n");
		}
		opts::addFile("Batch", stepname, filenamegm);
		bygroupmissing.precision(4);

		if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
			bygroup << "Chrom\trsid\tProbeID\tbploc\t" << data_set->get_locus(0)->getDetailHeaders();
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

		bygroupmissing << "Batch\tMiss_avg\tLog_miss\tNmiss" << endl;
		opts::addHeader(filenamegm, "Miss_avg");
		opts::addHeader(filenamegm, "Log_miss");
		opts::addHeader(filenamegm, "Nmiss");

	}

	opts::addFile("Marker", stepname, filename);
	bymarker.precision(4);
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		bymarker << "Chrom\trsid\tProbeID\tbploc\t" << data_set->get_locus(0)->getDetailHeaders() << "\t%GenoEff_Ind_All\tInd_Zero_Count\tTotal_Individuals_Used\t%GenoEff_Ind_Cases\tInd_Zero_Count_Cases\tTotal_Individuals_Used_Cases\t%GenoEff_Ind_Controls\tInd_Zero_Count_Controls\tTotal_Individuals_Used_Controls";
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


	int size = good_markers.size();//data_set->num_loci();
	int prev_base = 0;
	int prev_chrom = -1;
	map<string, double> group_avgs;
	int total_snps = 0;
	for(int i = 0; i < size; i++){
		if(good_markers[i]->isEnabled()){//isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom) && !data_set->get_locus(i)->isFlagged()){

			float percent = 0.0f;
			float caseper = 0.0f;
			float contper = 0.0f;
	//		int contzero = 0;
	//		int conttot = 0;
			if(total[i] > 0){
				percent = (1.0f - ((float)zeros[i]/(float)total[i]));// * 100.0f;
			}
			if(casetotal[i] > 0){
				caseper = (1.0f - ((float)casezeros[i]/(float)casetotal[i]));// * 100.0f;
			}
	//		contzero = zeros[i] - casezeros[i];
	//		conttot = total[i] - casetotal[i];
			if(controltotal[i] > 0){
				contper = (1.0f - ((float)controlzeros[i]/(float)controltotal[i]));// * 100.0f;
			}

			bymarker << good_markers[i]->toString() << "\t"
					 << percent << "\t"
					 << zeros[i] << "\t"
					 << total[i] << "\t"
					 << caseper << "\t"
					 << casezeros[i] << "\t"
					 << casetotal[i] << "\t"
					 << contper << "\t"
					 << controlzeros[i] << "\t"
					 << controltotal[i];
			bymarker << endl;
			if(options.doGroupFile()){
				bygroup << good_markers[i]->toString();
				map<string, vector<Sample*> > groups = options.getGroups();
				map<string, vector<Sample*> >::iterator giter;
				for(giter = groups.begin(); giter != groups.end(); giter++){
					string group = giter->first;
					int gzero = groupzeros[group][i];
					int gtotal = grouptotal[group][i];
					float gper = 0.0f;
					if(gtotal > 0){
						gper = (float)(1.0f - ((float)gzero/(float)gtotal));// * 100.0f;
					}
					bygroup << "\t" << gper;
					bygroup << "\t" << gzero;
					bygroup << "\t" << gtotal;

					group_avgs[group] += ((double) gzero / (double) gtotal);
				}
				bygroup << endl;
			}
			total_snps++;
		}
		good_markers[i]->setFlag(false);
	}

	if(options.doGroupFile()){
		map<string, vector<Sample*> > groups = options.getGroups();
		map<string, vector<Sample*> >::iterator giter;
		for(giter = groups.begin(); giter != groups.end(); giter++){
			string group = giter->first;

			int goodsamps = 0;
			for(int gs = 0; gs < (int)groups[group].size(); gs++){
				if(groups[group][gs]->isEnabled()){
					goodsamps++;
				}
			}

			double avg = group_avgs[group] / (double) total_snps;
			double avg_log = log10(avg);
			bygroupmissing << group << "\t";
			bygroupmissing << (group_avgs[group] / (double) total_snps) << "\t";
			bygroupmissing << avg_log << "\t";
			bygroupmissing << goodsamps;
			bygroupmissing << endl;
		}
 	}

	if(bymarker.is_open()){
		bymarker.close();
	}
	if(bygroup.is_open()){
		bygroup.close();
	}
	if(bygroupmissing.is_open()){
		bygroupmissing.close();
	}

}

void ProcessMarkerGenoEff::filter(){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int size = good_markers.size();
		int prev_base = 0;
		int prev_chrom = -1;
		for(int i = 0; i < size; i++){
			if(good_markers[i]->isEnabled()){//isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom) && !data_set->get_locus(i)->isFlagged()){
				double percent = 0.0f;
				float caseper = 0.0f;
				float contper = 0.0f;
				bool inc = false;
//		int contzero = 0;
//		int conttot = 0;
				if(total[i] > 0){
					percent = (double)((1.0f - ((double)zeros[i]/(double)total[i])));// * 100.0f);
				}
				if(casetotal[i] > 0){
					caseper = (1.0f - ((float)casezeros[i]/(float)casetotal[i]));// * 100.0f;
				}
//		contzero = zeros[i] - casezeros[i];
//			//		conttot = total[i] - casetotal[i];
				if(controltotal[i] > 0){
					contper = (1.0f - ((float)controlzeros[i]/(float)controltotal[i]));// * 100.0f;
				}

				if(options.doThreshMarkersLow() && dLess(percent, options.getThreshMarkersLow())){
					good_markers[i]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && dGreater(percent, options.getThreshMarkersHigh())){
					good_markers[i]->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}
		}
	}
}

void ProcessMarkerGenoEff::process(DataSet* ds){
	data_set = ds;
	good_markers = findValidMarkers(data_set->get_markers(), &options);
	int msize = good_markers.size();
	zeros.resize(msize);
	total.resize(msize);
	casezeros.resize(msize);
	casetotal.resize(msize);
	controlzeros.resize(msize);
	controltotal.resize(msize);

	if(options.doGroupFile()){
		options.readGroups(data_set->get_samples());
		map<string, vector<Sample*> >::iterator giter;
		map<string, vector<Sample*> > groups = options.getGroups();

		for(giter = groups.begin(); giter != groups.end(); giter++){
			string mygroup = giter->first;
			groupzeros[mygroup].resize(msize);
			grouptotal[mygroup].resize(msize);
		}
	}

	MarkerGenoEff mge(data_set);
//	mge.resetDataSet(data_set);
	mge.set_parameters(&options);

	for(int i = 0; i < msize; i++){//(int)data_set->num_loci(); i++){
		mge.calcOne(good_markers[i]);
		zeros[i] = mge.getZeros();
		total[i] = mge.getTotal();
		casezeros[i] = mge.getCaseZeros();
		casetotal[i] = mge.getCaseTotal();
		controlzeros[i] = mge.getControlZeros();
		controltotal[i] = mge.getControlTotal();
		if(options.doGroupFile()){
			map<string, int> gzeros = mge.getGroupZeros();
			map<string, int> gtotal = mge.getGroupTotal();
			map<string, int>::iterator iter;
			for(iter = gzeros.begin(); iter!= gzeros.end(); iter++){
				groupzeros[iter->first][i] = iter->second;
			}
			for(iter = gtotal.begin(); iter != gtotal.end(); iter++){
				grouptotal[iter->first][i] = iter->second;
			}
		}
	}

}

