/**********************************************************************************
*                       Average Quality Score Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all genotypes and generates an average quality score for all markers. 
* 
*
* Files generated:
*
*File: QualityScore.cc
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
#include "QualityScore.h"
#include "Helper.h"
#include "General.h"
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"


void QualityScore::FilterSummary(){
/*	ofstream myoutput ("post_quality_score_filter_summary.txt");
	myoutput.precision(4);

	myoutput << "Threshold:\t" << threshold << endl;
	myoutput << "Markers Passed:\t" << markers->getSize() << " (" << 
		((float)markers->getSize() / (float)orig_num_markers) * 100.0 <<
		"%) of " << orig_num_markers << endl;
	myoutput << "Families Passed:\t" << families->getSize() << " (" <<
		((float)families->getSize() / (float)orig_num_families) * 100.0 <<
		"%) of " << orig_num_families << endl;
    myoutput << "Individuals Passed:\t" << families->getNumInds() << " (" <<
		        ((float) families->getNumInds() / (float) orig_num_individuals) * 100.0 <<
		        "%) of " << orig_num_individuals << endl;
	
	myoutput.close();
*/
}

void QualityScore::PrintSummary(){
	//map<string,Chrom, Chrom::mysort> chroms;
//	MKR* mymarkers = markers->getList();
//	MKR::iterator m_iter;
//	ofstream bymarker ("ave_quality_score_by_marker.txt");
//	bymarker.precision(6);
	string fname1 = opts::_OUTPREFIX_ + "quality_score_by_bin_" + getString<int>(order) + ".txt";
	ofstream bins (fname1.c_str());
	if(!bins){
		cerr << "Error opening quality_score_by_bin.txt. Exiting!" << endl;
		exit(1);
	}

//	if(!bymarker.is_open()){
//		cout << "Unable to open ave_quality_score_by_marker.txt" << endl;
//		exit(1);
//	}
//	else{
//		bymarker.precision(4);
//		bymarker << "Chrom\trsID\tProbeID\tbploc\tEnzyme\tAverage_Quality_Score" << endl;
//	}
	if(!bins.is_open()){
		cout << "Unable to open quality_score_by_bin.txt" << endl;
		exit(1);
	}
	else{
		bins << "Score_Range\tAll_genotypes\tHET_genotypes\tHomoMajorAllele_genotypes\tHomoMinorAllele_genotypes" << endl;//\tDelta_all\tDelta_FatherMother\tDelta_FatherChild\tDelta_MotherChild" << endl;
	}
	map<float,int>::iterator b_iter;
	for(b_iter = qs_total_het.begin(); b_iter != qs_total_het.end(); b_iter++){
		bins << b_iter->first << "\t";
		int all = b_iter->second + qs_total_maj[b_iter->first] + qs_total_min[b_iter->first];
		bins << all << "\t" << b_iter->second << "\t" << qs_total_maj[b_iter->first] << "\t" << qs_total_min[b_iter->first] << endl;
	}

	bins.close();

	string fname2 = opts::_OUTPREFIX_ + "ave_quality_score_" + getString<int>(order) + ".txt";
	ofstream aveout (fname2.c_str());

	if(!aveout.is_open()){
		cout << "Unable to open ave_quality_score.txt" << endl;
		exit(1);
	}

	aveout << "Chrom\trsID\tProbe ID\tbploc\tEnzyme\tAverage QS\tTotal Genotypes" << endl;
	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		Marker* mark = (*markers)[i];

		float ave = (float) qs_ave[i] / (float) qs_tot[i];
		aveout << mark->getChrom() << "\t" << mark->getRSID() << "\t" << mark->getProbeID() << "\t" << mark->getBPLOC() << "\t" << mark->getEnzyme() << "\t" << ave << "\t" << qs_tot[i] << endl;
	}	
	aveout.close();
/*	map<int, vector<int> > qs_bin;
	for(m_iter = mymarkers->begin(); m_iter != mymarkers->end(); m_iter++){
		Marker m = m_iter->second;
//		map<string,Chrom>::iterator c_iter;
//		m.calcPerEff();
//		c_iter = chroms.find(m.getChrom());
//		if(c_iter == chroms.end()){
//			Chrom* c = new Chrom(m.getChrom());
//			c->incPC(m.getPerEff());
//			c->incTotal();
//			c->insert_marker(m);
//			chroms[m.getChrom()] = *c;
//			delete c;
//		}
//		else{
//			c_iter->second.incPC(m.getPerEff());
//			c_iter->second.incTotal();
//			c_iter->second.insert_marker(m);
//		}
	
	//	vector<int> qs_tot = m.getQSAll();
		vector<int> qs_het = m.getQSHet();
		vector<int> qs_maj = m.getQSMaj();
		vector<int> qs_min = m.getQSMin();
		vector<int> qs_d_fm = m.getDeltaFM();
		vector<int> qs_d_fc = m.getDeltaFC();
		vector<int> qs_d_mc = m.getDeltaMC();
		for(int i = 0; i < 50; i++){
			if(qs_bin[i].size() == 0){
				qs_bin[i].push_back(0);
				qs_bin[i].push_back(0);
				qs_bin[i].push_back(0);
				qs_bin[i].push_back(0);
				qs_bin[i].push_back(0);
				qs_bin[i].push_back(0);
				qs_bin[i].push_back(0);
				qs_bin[i].push_back(0);
			}
			else{
				int het = 0;
			   	int	maj = 0;
				int min = 0;
				int fm = 0;
				int fc = 0;
				int mc = 0;
				if(qs_het.size() > i){
					het = qs_het[i];
				}
				if(qs_maj.size() > i){
					maj = qs_maj[i];
				}
				if(qs_min.size() > i){
					min = qs_min[i];
				}
				if(qs_d_fm.size() > i){
					fm = qs_d_fm[i];
				}
				if(qs_d_fc.size() > i){
					fc = qs_d_fc[i];
				}
				if(qs_d_mc.size() > i){
					mc = qs_d_mc[i];
				}
				int all = het + maj + min;
				int d_all = fm + fc + mc;
				qs_bin[i][0] += all;
				qs_bin[i][1] += het;
				qs_bin[i][2] += maj;
				qs_bin[i][3] += min;
				qs_bin[i][4] += d_all;
				qs_bin[i][5] += fm;
				qs_bin[i][6] += fc;
				qs_bin[i][7] += mc;
			}
		}

	}

	map<int, vector<int> >::iterator bin_iter;
	for(bin_iter = qs_bin.begin(); bin_iter != qs_bin.end(); bin_iter++){
		if(bin_iter->first < 10){
			bins << "0.0" << bin_iter->first;
		}
		else{
			bins << "0." << bin_iter->first;
		}
		for(int i = 0; i < 8; i++){
			bins << "\t" << bin_iter->second[i];
		}
		bins << endl;
	}
		
		
//	map<string,Chrom>::iterator c_iter;
//	for(c_iter = chroms.begin(); c_iter != chroms.end(); c_iter++){
			
//		BPMKR* cmarks = c_iter->second.getBPMarkers();
//		sort(cmarks->begin(), cmarks->end());
//		BPMKR::iterator bm_iter;
//		for(bm_iter = cmarks->begin(); bm_iter != cmarks->end(); bm_iter++){
//			bymarker << c_iter->first << "\t" << bm_iter->getRSID() << "\t" << bm_iter->getProbe() << "\t" << bm_iter->getBPLOC() << "\t" << bm_iter->getEnzyme() << "\t" << bm_iter->getQualScore() << endl;
//		}	
//	}
//	bymarker.close();
	bins.close();
//	chroms.clear();
*/
}

void QualityScore::filter(){
/*	MKR* mymarkers = markers->getList();
	MKR::iterator m_iter;
	orig_num_markers = markers->getSize();
	orig_num_families = families->getSize();
 	orig_num_individuals = families->getNumInds();	
	cout << "Doing filter" << endl;
	int total = 0;
	for(m_iter = mymarkers->begin(); m_iter != mymarkers->end();){
		Marker m = m_iter->second;
		if(m.getPerEff() > threshold){
			mymarkers->erase(m_iter++);
		}
		else{
			++m_iter;
		}
		total++;
	}
	cout << "total reviewed: " << total << endl;
*/
}

//void QualityScore::process(Families* f, Markers* m){
//}


int QualityScore::get_marker_loc(int s){
	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		if(s == (*markers)[i]->getLoc()){
			return i;
		}
	}
	return -1;
}

void QualityScore::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
//	if(markers == NULL){
//		cout << "Markers are empty...Obtaining Markers from DB..." << endl;
//		markers = new Markers();
//		markers->fillFromDB(con);
//		cout << "Markers obtained: " << markers->getSize() << endl;
//		*m = *markers;
//	}
//	if(families == NULL){
//		cout << "Families are empty...Obtaining Families from DB..." << endl;
//		families = new Families();
//		families->fillFromDB(con);
//		cout << "Families obtained: " << families->getSize() << endl;
//		*f = *families;
//	}

//	orig_num_markers = markers->getSize();
//	orig_num_families = families->getSize();
//	orig_num_individuals = families->getNumInds();	

cout << rank << " Doing AlleleFrequency" << endl; 
        AlleleFrequency* af = new AlleleFrequency();
		af->setRank(rank);
        af->process(s, f, m, mm);
cout << rank << " Done with AlleleFrequency" << endl;
//	delete(af);

	ifstream qs_input;
	qs_input.open(opts::_QSFILE_.c_str(), ios::in);
	float elems[6];
	int size = markers->size();
	qs_ave.resize(size);
	qs_tot.resize(size);
	int ind_loc = 0;
	while(!qs_input.eof()){
		qs_input >> elems[0]; //famid
		qs_input >> elems[1]; //indid
		qs_input >> elems[2]; //dad
		qs_input >> elems[3]; //mom
		qs_input >> elems[4]; //gender
		qs_input >> elems[5]; //aff

		Sample* samp = (*samples)[ind_loc];
		if(samp->isEnabled()){
			float qs = 0.0f;
			for(int i = 0; i < size; i++){
				qs_input >> qs;
				if(qs >= 0){
					for(float j = min_bin; j < max_bin; j+= bin_space){
						if(qs < j){
							j -= bin_space;
							int m_loc = get_marker_loc(i);
							if(m_loc == -1){
								cerr << "Cannot find marker!!! " << i << endl;
								exit(1);
							}
							qs_ave[m_loc]+= qs;
							qs_tot[m_loc]++;
							if(!samp->getAone(m_loc) && samp->getAtwo(m_loc)){
								qs_total_het[j]++;
							}
							else if(!samp->getAone(m_loc) && !samp->getAtwo(m_loc) && af->getAoneP_freq(m_loc) >= af->getAtwoP_freq(m_loc)){
								qs_total_maj[j]++;
							}
							else if(!samp->getAone(m_loc) && !samp->getAtwo(m_loc) && af->getAoneP_freq(m_loc) <= af->getAtwoP_freq(m_loc)){
								qs_total_min[j]++;
							}
							else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && af->getAoneP_freq(m_loc) <= af->getAtwoP_freq(m_loc)){
								qs_total_maj[j]++;
							}
							else if(samp->getAone(m_loc) && samp->getAtwo(m_loc) && af->getAoneP_freq(m_loc) >= af->getAtwoP_freq(m_loc)){
								qs_total_min[j]++;
							}
						
							break;
						}
					}
				}
			}	
		}
		else{
			string temp = "";
			getline(qs_input, temp);
		}
		ind_loc++;
	}
	qs_input.close();


//	FAM::iterator fam_iter;
//	FAM* myfams = families->getList();
//	MKR* mymarkers = markers->getList();
//	bool _CASECONTROL_ = families->getCaseControl();


/*	
	MKR::iterator m_iter;
	//set major allele and cleanup extraneous allele freq. data
	for(m_iter = mymarkers->begin(); m_iter != mymarkers->end(); m_iter++){
		//m_iter->second.calcFreqs();
		m_iter->second.setMajorAllele();
		//m_iter->second.removeAlleleCalcs();
	}
cout << rank << " Doing qualityscore stuff" << endl;
	Statement* s = NULL;
	string sql = "SELECT ri.sysprobe, ri.allele_1, ri.allele_2, ri.METRIC_VALUE, ind.ind, ind.fam FROM clinic_user.result_info ri, clinic_user.ind ind, clinic_user.samp_def sd ";
#ifdef _USEMPI_
//	ostringstream ostr;
//	ostr << rank;
//	sql += " , TODO_MARKERS_" + ostr.str() + " tm ";
#endif
	
	sql += " WHERE ri.syssamp = sd.syssamp AND sd.sysind = ind.sysind and ri.allele_1 != '0' and ri.allele_2 != '0' AND ri.status != 'B' AND ri.metric_value < 0.1 ";
	if(_CASECONTROL_){
		sql += " AND ind.fam is NULL ";
	}
	else{
		sql += " AND ind.fam is not NULL  ";
	}
#ifndef _USEMPI_
    if(_MARKERLIST_){
        sql += " AND (";
        MKR::iterator msql_iter;
        for(msql_iter = mymarkers->begin(); msql_iter != mymarkers->end(); msql_iter++){
            ostringstream ostr;
            ostr << msql_iter->first;
            sql += "ri.sysprobe = " + ostr.str() + " OR ";
        }
        sql.erase(sql.length() - 3);
        sql += ") ";
    }
#endif

#ifdef _USEMPI_
//	sql += " AND ri.sysprobe = tm.sysprobe ";
	 	sql += " AND (";	
	MKR::iterator msql_iter;
	for(msql_iter = mymarkers->begin(); msql_iter != mymarkers->end(); msql_iter++){
		ostringstream ostr;
		ostr << msql_iter->first;
		sql += "ri.sysprobe = " + ostr.str() + " OR ";
	}
	sql.erase(sql.length() - 3);
	sql += ") ";
#endif
	sql += "ORDER BY ri.sysprobe ASC, ind.fam asc, ind.ind asc";
	try{
		s = con->createStatement(sql.c_str());
	}
	catch(SQLException ea){
		cerr << "Error creating statement: " << ea.what() << endl;
		exit(1);
	}
	
//	for(fam_iter = myfams->begin(); fam_iter != myfams->end(); fam_iter++){
//		cout << "Looking at family: " << fam_iter->getFamID() << endl;
		try{	
//			s->setNumber(1, fam_iter->getFamID());
			s->setPrefetchRowCount(MAXROWS);
			ResultSet* r = s->executeQuery();

			Marker* temp_marker = NULL;
			int prev_fam = 0;
			int prev_sysprobe = 0;
			float f_qs = -1;
			float m_qs = -1;
			float c_qs = -1;
			while(r->next()){
				int sysprobe = (int) r->getNumber(1);
				//string chrom = r->getString(2);
				string a1 = r->getString(2);
				string a2 = r->getString(3);
				a1 = a1.substr(0,1);
				a2 = a2.substr(0,1);
				float score = (float) r->getNumber(4);
				int ind = (int) r->getNumber(5);
				int famid = (int) r->getNumber(6); //
			
				if(!prev_fam){
					prev_fam = famid;
				}
				if(!prev_sysprobe){
					prev_sysprobe = sysprobe;
				}

				if(prev_fam != famid){
					if(f_qs > -1 && m_qs > -1 && c_qs > -1){
						MKR::iterator mkr_iter;
						mkr_iter = mymarkers->find(prev_sysprobe);
						if(mkr_iter != mymarkers->end()){
							mkr_iter->second.incDeltaFM(fabs(f_qs - m_qs));
							mkr_iter->second.incDeltaFC(fabs(f_qs - c_qs));
							mkr_iter->second.incDeltaMC(fabs(m_qs - c_qs));
						}
					}
					f_qs = -1;
					m_qs = -1;
					c_qs = -1;
					prev_fam = famid;
				}
				if(prev_sysprobe != sysprobe){
					prev_sysprobe = sysprobe;
				}
				
				Family* temp_fam = families->getFamily(famid);
				if(temp_fam != NULL){
					Individual* temp_ind = temp_fam->getInd(ind);
					if(temp_ind != NULL){
						if(temp_ind->isEnabled()){
							MKR::iterator mkr_iter;
							mkr_iter = mymarkers->find(sysprobe);
							if(mkr_iter != mymarkers->end()){
								if(temp_marker == NULL){
									temp_marker = &(*mymarkers)[sysprobe];
								}
								if(temp_marker->getSysprobe() != sysprobe){
									(*mymarkers)[temp_marker->getSysprobe()] = *temp_marker;
									temp_marker = &(*mymarkers)[sysprobe];
								}
								temp_marker->setQualScore(score);
								if(a1 == a2 && a1 != "0"){
									if(a1 == temp_marker->getMajorAllele()){
										temp_marker->incQSMaj(score);
									}
									else{
										temp_marker->incQSMin(score);
									}
								}
								else if(a1 != a2){
									temp_marker->incQSHet(score);
								}
								if(temp_ind->getMother()){
									m_qs = score;
								}
								if(temp_ind->getFather()){
									f_qs = score;
								}
								if(temp_ind->getChild()){
									c_qs = score;
								}	
								//temp_marker->incQSAll(score);
							}
						}
					}
				}
			}	
					if(f_qs > -1 && m_qs > -1 && c_qs > -1){
						MKR::iterator mkr_iter;
						mkr_iter = mymarkers->find(prev_sysprobe);
						if(mkr_iter != mymarkers->end()){
							mkr_iter->second.incDeltaFM(fabs(f_qs - m_qs));
							mkr_iter->second.incDeltaFC(fabs(f_qs - c_qs));
							mkr_iter->second.incDeltaMC(fabs(m_qs - c_qs));
						}
					}
			s->closeResultSet(r);
		}
		catch(SQLException ea){
			cerr << "Error connecting: " << ea.what() << endl;
			exit(1);
		}
	//}

	try{
		con->terminateStatement(s);
	}
	catch(SQLException ea){
		cerr << "Error closing: " << ea.what() << endl;
		exit(1);
	}
cout << rank << " Done with all quality score stuff" << endl;	
*/
}

