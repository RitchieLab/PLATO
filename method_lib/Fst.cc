/**********************************************************************************
*                       Fst Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
*
*File: Fst.cc
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
#include "AlleleFrequency.h"
#include "HWEquilibrium.h"
#include "Fst.h"
#include "Options.h"
#include "General.h"
#include "Helpers.h"
namespace Methods{
string Fst::stepname = "fst";

void Fst::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void Fst::PrintSummary(){
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
	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	size = good_markers.size();
	for(int i = 0; i < size; i++){
		if(good_markers.at(i)->isEnabled()){

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

			bymarker << good_markers.at(i)->toString() << "\t"
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
				bygroup << good_markers.at(i)->toString();
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
	}

	if(bymarker.is_open()){
		bymarker.close();
	}

}

void Fst::filterOne(int m){
}

void Fst::filter(){
}

void Fst::calcOne(int m){
	Marker* mark = data_set->get_locus(m);
	af.initializeCounts(0);
	af.calculate(m);
	af.calculateGroups(m);

	map<string, vector<Sample*> > groups = options.getGroups();
	map<string, vector<Sample*> >::iterator giter;
	map<string, float> group_freq = options.getGroupFreq();


	if (groups.size() == samples->size())
	{
		opts::printLog("Total number of samples should be different from number of groups specified.\n");
		throw MethodException("Total number of samples should be different from number of groups specified.\n");
	}


	//FSTWC
	double nc = 0;
	double den = 0;

	for(giter = groups.begin(); giter != groups.end(); giter++)
	{

		nc = nc + (af.getGroupPop(giter->first));
		den = den + af.getGroupPop(giter->first);

	}

	nc = (den - (nc/den))/(double)(groups.size() - 1);


	double MSG = 0;
	double pmed = 0;

	for(giter = groups.begin(); giter != groups.end(); giter++)
	{
		if(!options.getDoGroupFreq())
		{
			MSG = MSG + (af.getGroupPop(giter->first) * af.getGroupAone_freq(giter->first) * af.getGroupAtwo_freq(giter->first));
			pmed = pmed + (af.getGroupPop(giter->first) * af.getGroupAone_freq(giter->first));
		}
		else
		{
			MSG = MSG + (af.getGroupPop(giter->first) * group_freq[mark->getRSID() + " " + giter->first] * (1.0f - group_freq[mark->getRSID() + " " + giter->first]));
			pmed = pmed + (af.getGroupPop(giter->first) * group_freq[mark->getRSID() + " " + giter->first]);


		}
	}


	MSG = MSG / (den - groups.size());
	pmed = pmed / den;


	double MSP = 0;

	for(giter = groups.begin(); giter != groups.end(); giter++)
	{
		if(!options.getDoGroupFreq())
		{
			MSP = MSP + (af.getGroupPop(giter->first) * ((af.getGroupAone_freq(giter->first) - pmed)*(af.getGroupAone_freq(giter->first) - pmed)));
		}
		else
		{
			MSP = MSP + (af.getGroupPop(giter->first) * ((group_freq[mark->getRSID() + " " + giter->first] - pmed) * (group_freq[mark->getRSID() + " " + giter->first] - pmed)));

		}
	}


	MSP = MSP / (double)(groups.size() - 1);


	double FstWC = 0;
	FstWC= (MSP - MSG) / (MSP + ((nc - 1) * MSG));


	if(FstWC < 0)
	{
		FstWC = 0;
	}

	fst = FstWC;
	//end FSTWC


	//FSTRH
	double Hs = 0;
	double pbar = 0;
	double Ht = 0;
	double FstRH = 0;
	if(groups.size() > 1){
		for(giter = groups.begin(); giter != groups.end(); giter++){
			float het_exp = hwe.calcHW(af.getGroupAone_freq(giter->first), af.getGroupAtwo_freq(giter->first),
					af.getGroupAonehomo(giter->first), af.getGroupHet(giter->first), af.getGroupAtwohomo(giter->first),
					af.getGroupPop(giter->first));

			pbar = (af.getGroupAone_freq(giter->first) * 2 * af.getGroupPop(giter->first)) + pbar;
			Hs = Hs + (het_exp * af.getGroupPop(giter->first));
		}
		pbar = pbar / (2 * den);
		double qbar = 1 - pbar;
		Ht = 1 - ((pbar * pbar) + (qbar * qbar));
		Hs = Hs / den;
	}
	if(Ht == 0){
		FstRH = -1;
	}
	else{
		FstRH = (Ht - Hs) / Ht;
		if(FstRH < 0){
			FstRH = 0;
		}
	}
	fstrh = FstRH;
	//end FSTWH

	double FstHM = 0;
	if(groups.size() > 1){
		double binom = 0;
		double numer = 0;
		double a = 0;
		double b = 0;
		double c = 0;
		for(giter = groups.begin(); giter != groups.end(); giter++){
			double pop = af.getGroupPop(giter->first);
			a = ((2 * pop) * (2 * (pop - 1))) / (double)2;
			b = 2 * ((2 * pop) / (2 * (pop - 1)));
			c = af.getGroupAone_freq(giter->first) * af.getGroupAtwo_freq(giter->first);

			numer = numer + (a * b * c);
			binom = binom + a;
		}

		numer = numer / binom;
		pmed = 0;
		double k = 0;
		for(giter = groups.begin(); giter != groups.end(); giter++){
			pmed = pmed + af.getGroupAonehomo(giter->first);
			k = k + (2* af.getGroupPop(giter->first));
		}

		pmed = pmed / k;

		den = 2 * den;

		if(FstHM != -1){
			a = 2 * (den / (den - 1));
			b = pmed * (1 - pmed);
			c = a * b;
			FstHM = 1 - (numer / c);
			if(FstHM < 0){
				FstHM = 0;
			}
		}
	}
	fsthm = FstHM;

}

void Fst::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
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

	if(options.doGroupFile())
	{
		options.readGroups(samples);
		map<string, vector<Sample*> >::iterator giter;
		map<string, vector<Sample*> > groups = options.getGroups();

		//check if the overall sample size is equal to the total number of groups.
		//If so, throw an error and alert the user that the total number of samples should differ from the number of groups specified
		if (groups.size() == samples->size())
		{
			opts::printLog("Total number of samples should be different from number of groups specified.\n");
			throw MethodException("Total number of samples should be different from number of groups specified.\n");
		}

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
