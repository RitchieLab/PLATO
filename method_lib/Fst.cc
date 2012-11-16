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
#include "Helper.h"
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
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
		exit(1);
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

		if((*markers)[0]->getDetailHeaders().size() > 0){
			bygroup << "Chrom\trsid\tProbeID\tbploc\t" << (*markers)[0]->getDetailHeaders();
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
	if((*markers)[0]->getDetailHeaders().size() > 0){
		bymarker << "Chrom\trsid\tProbeID\tbploc\t" << (*markers)[0]->getDetailHeaders() << "\t%GenoEff_Ind_All\tInd_Zero_Count\tTotal_Individuals_Used\t%GenoEff_Ind_Cases\tInd_Zero_Count_Cases\tTotal_Individuals_Used_Cases\t%GenoEff_Ind_Controls\tInd_Zero_Count_Controls\tTotal_Individuals_Used_Controls";
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
		if((*markers)[i]->isEnabled() && !(*markers)[i]->isFlagged()){

			float percent = 0.0f;
			float caseper = 0.0f;
			float contper = 0.0f;
	//		int contzero = 0;
	//		int conttot = 0;
			if(total[i] > 0){
				percent = (1.0f - ((float)zeros[i]/(float)total[i])) * 100.0f;
			}
			if(casetotal[i] > 0){
				caseper = (1.0f - ((float)casezeros[i]/(float)casetotal[i])) * 100.0f;
			}
	//		contzero = zeros[i] - casezeros[i];
	//		conttot = total[i] - casetotal[i];
			if(controltotal[i] > 0){
				contper = (1.0f - ((float)controlzeros[i]/(float)controltotal[i])) * 100.0f;
			}

			bymarker << (*markers)[i]->toString() << "\t"
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
				bygroup << (*markers)[i]->toString();
				map<string, vector<Sample*> > groups = options.getGroups();
				map<string, vector<Sample*> >::iterator giter;
				for(giter = groups.begin(); giter != groups.end(); giter++){
					string group = giter->first;
					int gzero = groupzeros[group][i];
					int gtotal = grouptotal[group][i];
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
		(*markers)[i]->setFlag(false);
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
	af.initializeCounts(0);
	af.calculate(m);
	af.calculateGroups(m);

	map<string, vector<Sample*> > groups = options.getGroups();
	map<string, vector<Sample*> >::iterator giter;






	//FSTWC
	double nc = 0;
	double den = 0;

	for(giter = groups.begin(); giter != groups.end(); giter++){

		nc = nc + (af.getGroupPop(giter->first));
		den = den + af.getGroupPop(giter->first);
//		cout << data_set->get_locus(m)->toString() << "\t" << af.getGroupAone_freq(giter->first) << endl;
	}

	nc = (den - (nc/den))/(double)(groups.size() - 1);


//	for(giter = groups.begin(); giter != groups.end(); giter++){
//		float het_exp = hwe.calcHW(af.getGroupAone_freq(giter->first), af.getGroupAtwo_freq(giter->first),
//				af.getGroupAonehomo(giter->first), af.getGroupHet(giter->first), af.getGroupAtwohomo(giter->first),
//				af.getGroupPop(giter->first));
//	}

	double MSG = 0;
	double pmed = 0;

	for(giter = groups.begin(); giter != groups.end(); giter++){
		MSG = MSG + (af.getGroupPop(giter->first) * af.getGroupAone_freq(giter->first) * af.getGroupAtwo_freq(giter->first));
		pmed = pmed + (af.getGroupPop(giter->first) * af.getGroupAone_freq(giter->first));
	}


	MSG = MSG / (den - groups.size());
	pmed = pmed / den;

	double MSP = 0;

	for(giter = groups.begin(); giter != groups.end(); giter++){
		MSP = MSP + (af.getGroupPop(giter->first) * ((af.getGroupAone_freq(giter->first) - pmed)*(af.getGroupAone_freq(giter->first) - pmed)));
	}


	MSP = MSP / (double)(groups.size() - 1);

	double FstWC = 0;
	FstWC= (MSP - MSG) / (MSP + ((nc - 1) * MSG));

	if(FstWC < 0){
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
//		FstWC = -1;
//		FstHM = -1;
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

	if(options.doGroupFile()){
		options.readGroups(samples);
		map<string, vector<Sample*> >::iterator giter;
		map<string, vector<Sample*> > groups = options.getGroups();

		for(giter = groups.begin(); giter != groups.end(); giter++){
			string mygroup = giter->first;
			vector<Sample*> mysamps = giter->second;
			groupzeros[mygroup].resize(markers->size());
			grouptotal[mygroup].resize(markers->size());

			for(int s = 0; s < (int)mysamps.size(); s++){
				Sample* samp = mysamps[s];
				if(samp->isEnabled()){
					int prev_base = 0;
					int prev_chrom = -1;
					for(int m = 0; m < (int)markers->size(); m++){
						Marker* mark = (*markers)[m];
						if(isValidMarker(mark, &options, prev_base, prev_chrom)){
							int loc = mark->getLoc();
							if(!mark->isMicroSat()){
								if(samp->getAone(loc) && !samp->getAtwo(loc)){
									groupzeros[mygroup][m]++;
								}
							}
							else{
								if(samp->getAbone(loc) == -1){
									groupzeros[mygroup][m]++;
								}
							}
							grouptotal[mygroup][m]++;
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
				//int loc = (*marker_map)[i];
				int loc = (*markers)[i]->getLoc();
				if(isValidMarker((*markers)[i], &options, prev_base, prev_chrom)){
//				if((*markers)[i]->isEnabled()){
//					if(options.doChrom()){
//						if(!options.checkChrom((*markers)[i]->getChrom())){
//							continue;
//						}
//						if(!options.checkBp((*markers)[i]->getBPLOC())){
//							continue;
//						}
//					}
 //                   if(options.doBpSpace()){
//						if(prev_base == 0){
//							prev_base = (*markers)[i]->getBPLOC();
//							prev_chrom = (*markers)[i]->getChrom();
//						}
//						else{
//							if((*markers)[i]->getChrom() == prev_chrom && (((*markers)[i]->getBPLOC() - prev_base) < options.getBpSpace())){
//								(*markers)[i]->setFlag(true);
//								continue;
//							}
//							prev_base = (*markers)[i]->getBPLOC();
//							prev_chrom = (*markers)[i]->getChrom();
//						}
//					}

					if(!(*markers)[i]->isMicroSat()){
						if((*s_iter)->getAone(loc) && !(*s_iter)->getAtwo(loc)){
							zeros[i]++;
							if((*s_iter)->getPheno() == 2){//getAffected()){
								casezeros[i]++;
							}
							else if((*s_iter)->getPheno() == 1){
								controlzeros[i]++;
							}
						}
					}
					else{
		//				int marrloc = (*s_iter)->getMicroSat(loc);
						if((*s_iter)->getAbone(loc) == -1){
							zeros[i]++;
							if((*s_iter)->getPheno() == 2){
								casezeros[i]++;
							}
							else if((*s_iter)->getPheno() == 1){
								controlzeros[i]++;
							}
						}
					}
					total[i]++;
					if((*s_iter)->getPheno() == 2){//getAffected()){
						casetotal[i]++;
					}
					else if((*s_iter)->getPheno() == 1){
						controltotal[i]++;
					}
				}
			}
		}
	}

}

}
