/**********************************************************************************
*                       Run of Homozygosity Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Performs a run of homozygosity search based on Lencz et al (PNAS 0710021104)
* Performs a ROH based on straight forward common homozygous spans.
*
*
*File: Homozygous.cc
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
#include "Homozygous.h"
#include "Options.h"
#include "General.h"
#include "Helpers.h"
#include "AlleleFrequency.h"
#include "ChiSquareAllelic.h"
#include "cdflib.h"
namespace Methods{
string Homozygous::stepname = "homozygous";

void Homozygous::FilterSummary(){

	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;
}

void Homozygous::PrintSummary(){
	if(options.doHomozygPermute()){
		return;
	}
	if(!options.doHomozygWGHA()){
		string filename = opts::_OUTPREFIX_ + "homozygous_marker" + options.getOut() + ".txt";
		if(!overwrite){
			filename += "." + getString<int>(order);
		}
		ofstream homo (filename.c_str());
		if(!homo.is_open()){
			opts::printLog("Unable to open " + filename + "\n");
			throw MethodException("Unable to open " + filename + "\n");
		}
		opts::addFile("Marker", stepname, filename);
		homo << "Chrom\trsID\tProbeID\tbploc\tUnAff\tAff\tTotal_Maj\tTotal_Min\tTotal\tTotal%\n";
		int msize = markers->size();
		opts::addHeader(filename, "UnAff");
		opts::addHeader(filename, "Aff");
		opts::addHeader(filename, "Total_Maj");
		opts::addHeader(filename, "Total_Min");
		opts::addHeader(filename, "Total");
		opts::addHeader(filename, "Total%");

		for(int i = 0; i < msize; i++){
			Marker* m = (*markers)[i];
			float per = (((float)homoallcount[i] / (float)opts::_SAMPLES_WORKING_) * 100.0f);
			homo.precision(4);
			homo << m->getChrom() << "\t" << m->getRSID() << "\t" << m->getProbeID() << "\t" << m->getBPLOC() << "\t" << homounaffcount[i] << "\t" << homoaffcount[i] << "\t" << homomajallcount[i] << "\t" << homominallcount[i] << "\t" << (homounaffcount[i] + homoaffcount[i]) << "\t" << per << endl;
		}
		if(homo.is_open()){
			homo.close();
		}
	}
}

void Homozygous::filter(){
}

double Homozygous::checkSequence(int begin, int end, int &numfounders, int &numchildren, double &perinseq){
	double prob = -1.0f;

	vector<double> freq1;
	vector<double> freq2;
	for(int m = begin; m <= end; m++){
		Marker* mark = good_markers[m];
		if(mark->isEnabled()){
			int mloc = mark->getLoc();
			int a1 = 0;
			int a2 = 0;

			int fsize = families->size();
			int sampcount = 0;
			for(int f = 0; f < fsize; f++){
				Family* fam = (*families)[f];
				if(fam->isEnabled()){
					vector<Sample*>* founders = fam->getFounders();
					for(int s = 0; s < (int)founders->size(); s++){
						Sample* samp = (*founders)[s];
						if(samp->isEnabled()){
							if(samp->getAone(mloc)){
								if(samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
									a1++;
									a1++;
									sampcount++;
								}
							}
							else{
								if(!samp->getAtwo(mloc)){
									a2++;
									a2++;
									sampcount++;
								}
								else{
									a1++;
									a2++;
									sampcount++;
								}
							}
						}
					}
				}
			}
			if( (a1 + a2) > 0){
				double a1f = (double)((double) a1 / (double)(a1 + a2));
				double a2f = 1.0 - a1f;
				if(prob < 0){
					prob = (a1f * a1f) + (a2f * a2f);
				}
				else{
					prob *= ((a1f * a1f) + (a2f * a2f));
				}
			}
		}
	}

	int totalchildren = 0;
	for(int s = 0; s < (int)samples->size(); s++){
		Sample* samp = (*samples)[s];
		if(samp->isEnabled()){
			int mcount = 0;
			int hcount = 0;
			bool bad = false;
			for(int m = begin; m <= end; m++){
				Marker* mark = good_markers[m];
				if(mark->isEnabled()){
					int mloc = mark->getLoc();
					if(samp->getAone(mloc)){
						if(samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
							hcount++;
						}
						else{
							bad = true;
							break;
						}
					}
					else{
						if(!samp->getAtwo(mloc)){
							hcount++;
						}
					}

					mcount++;
				}
			}
			if(!bad){
				if(hcount == mcount){
					if(samp->isFounder()){
						numfounders++;
					}
					else{
						numchildren++;
					}
				}
				if(!samp->isFounder()){
					totalchildren++;
				}
			}
		}
	}
	perinseq = (double)((double)numchildren/(double)totalchildren);
	return prob;
}

vector<Marker*> Homozygous::enabledMarkers(){
	vector<Marker*> enabled;
	for(int i = 0; i < (int)good_markers.size(); i++){
		if(good_markers[i]->isEnabled()){
			enabled.push_back(good_markers[i]);
		}
	}
	return enabled;
}

void Homozygous::process_wgha(){
	int msize = good_markers.size();
	int ssize = samples->size();
	homoallcount.resize(msize, 0);

	vector<Marker*> enabledmarkers = enabledMarkers();
	for(int s = 0; s < ssize; s++){
		Sample* samp = (*samples)[s];
		if(samp->isEnabled()){
			int spancount = 0;
			int currchrom = -1;
			int startloc = -1;
			int endloc = -1;
			int markerloc = 0;
			Marker* startmark = NULL;
			Marker* endmark = NULL;
			samp->resizeRoh(enabledmarkers.size());
			for(int i = 0; i < msize; i++){
				int loc = good_markers[i]->getLoc();
				if(good_markers[i]->getLoc()){
					if(currchrom == -1){
						currchrom = good_markers[i]->getChrom();
					}

					if(currchrom != good_markers[i]->getChrom()){
						//output data
						if(spancount >= options.getHomozygSpan()){
							for(int r = startloc; r <= endloc; r++){
								samp->setRoh(r, true);
							}
						}
						currchrom = good_markers[i]->getChrom();
						startmark = NULL;
						endmark = NULL;
						startloc = -1;
						endloc = -1;
						spancount = 0;
					}
					if(samp->getSex() && good_markers[i]->getChrom() == opts::_CHRX_){
						startmark = NULL;
						endmark = NULL;
						startloc = -1;
						endloc = -1;
						spancount = 0;
						continue;
					}
					if(!good_markers[i]->isMicroSat()){
						if((samp->getAone(loc) && samp->getAtwo(loc) && !samp->getAmissing(loc)) || (!samp->getAone(loc) && !samp->getAtwo(loc)) || (samp->getAone(loc) && samp->getAtwo(loc) && samp->getAmissing(loc))){
							if(startmark == NULL){
								startmark = good_markers[i];
								endmark = good_markers[i];
								startloc = markerloc;
								endloc = markerloc;
							}
							endloc = markerloc;
							endmark = good_markers[i];
							spancount++;
						}
						else if((!samp->getAone(loc) && samp->getAtwo(loc))){
							if(spancount >= options.getHomozygSpan()){
								for(int r = startloc; r <= endloc; r++){
									samp->setRoh(r, true);
								}
							}
							spancount = 0;
							startmark = NULL;
							endmark = NULL;
							startloc = -1;
							endloc = -1;
						}
					}
					else{//is microsatellite
						if((samp->getAbone(loc) == samp->getAtwo(loc)) || (samp->getAbone(loc) == -1 && samp->getAbtwo(loc) == -1)){
							if(startmark == NULL){
								startmark = good_markers[i];
								endmark = good_markers[i];
								startloc = markerloc;
								endloc = markerloc;
							}
							endloc = markerloc;
							endmark = good_markers[i];
							spancount++;
						}
						else if((samp->getAbone(loc) != samp->getAbtwo(loc))){
							if(spancount >= options.getHomozygSpan()){
								for(int r = startloc; r <= endloc; r++){
									samp->setRoh(r, true);
								}
							}
							spancount = 0;
							startmark = NULL;
							endmark = NULL;
							startloc = -1;
							endloc = -1;
						}
					}
					markerloc++;
				}//end marker enabled
			}//end foreach marker
			if(spancount >= options.getHomozygSpan()){
				for(int r = startloc; r <= endloc; r++){
					samp->setRoh(r, true);
				}
			}
		}//end samp enabled
	}//end foreach sample


	//remove < 10
	for(int m = 0; m < (int)enabledmarkers.size(); m++){
		int count = 0;
		for(int s = 0; s < ssize; s++){
			Sample* samp = (*samples)[s];
			if(samp->isEnabled()){
				if(samp->getRohVal(m)){
					count++;
				}
			}//end samp enabled
		}//end foreach sample
		if(count < options.getHomozygMinSamp()){
			for(int s = 0; s < ssize; s++){
				Sample* samp = (*samples)[s];
				if(samp->isEnabled()){
					samp->setRoh(m, false);
				}
			}
		}
	}//end foreach enabled markers


	//construct common ROH
	vector<int> commonroh;
	int count = 0;
	int start = -1;
	int end = -1;
	int prev_chrom = -1;
	for(int m = 0; m < (int)enabledmarkers.size(); m++){
		Marker* mark = enabledmarkers[m];
		int indcount = 0;
		if(prev_chrom == -1){
			prev_chrom = mark->getChrom();
		}

		if(prev_chrom != mark->getChrom()){
			if(count > options.getHomozygSpan()){
				commonroh.push_back(start);
				commonroh.push_back(end);
			}
			start = -1;
			end = -1;
			count = 0;
		}
		for(int s = 0; s < (int)samples->size(); s++){
			Sample* samp = (*samples)[s];
			if(samp->isEnabled()){
				if(samp->getRohVal(m)){
					indcount++;
				}
			}
		}
		if(indcount >= options.getHomozygMinSamp()){
			if(start == -1){
				start = m;
			}
			end = m;
			count++;
		}
		else{
			if(count > options.getHomozygSpan()){
				commonroh.push_back(start);
				commonroh.push_back(end);
			}
			start = -1;
			end = -1;
			count = 0;
		}
		homoallcount[m]+=indcount;
		prev_chrom = mark->getChrom();
	}
	if(count > options.getHomozygSpan()){
		commonroh.push_back(start);
		commonroh.push_back(end);
	}


	string filename = opts::_OUTPREFIX_ + "homozygous_wgha" + options.getOut() + ".txt";
	if(!overwrite){
		filename += "." + getString<int>(order);
	}

	filenames.push_back(filename);

	ofstream output (filename.c_str());
	if(!output.is_open()){
		opts::printLog("Unable to open " + filename + "\n");
		throw MethodException("Unable to open " + filename + "\n");
	}
	output << "Chrom\tStart_rsID\tStart_bploc\tEnd_rsID\tEnd_bploc\tbploc_span\tNum_snps\tCaseFound:CaseNotFound (expected)\tControlFound:ControlNotFound (expected)\tCase_Minor_Count\tCase_Major_Count\tControl_Minor_Count\tControl_Major_Count\tCase_Strict_Hom\tControl_Strict_Hom\tChisq\tPval\n";
	string cfilename = opts::_OUTPREFIX_ + "homozygous_wgha_counts" + options.getOut() + ".txt";
	if(!overwrite){
		cfilename += "." + getString<int>(order);
	}

	filenames.push_back(cfilename);

	ofstream coutput (cfilename.c_str());
	if(!coutput.is_open()){
		opts::printLog("Unable to open " + cfilename + "\n");
		throw MethodException("Unable to open " + cfilename + "\n");
	}
	opts::addFile("Marker",stepname, cfilename);
	coutput << "Chrom\trsID\tProbeID\tbploc\tInd_Count\n";
	opts::addHeader(cfilename, "Ind_Count");

	for(int m = 0; m < msize; m++){
		Marker* mark = good_markers[m];
		if(mark->isEnabled()){
			coutput << mark->toString() << "\t" << homoallcount[m] << endl;
		}
	}
	string sfilename = opts::_OUTPREFIX_ + "homozygous_wgha_cases" + options.getOut() + ".txt";
	if(!overwrite){
		sfilename += "." + getString<int>(order);
	}

	filenames.push_back(sfilename);

	ofstream soutput (sfilename.c_str());
	if(!soutput.is_open()){
		opts::printLog("Unable to open " + sfilename + "\n");
		throw MethodException("Unable to open " + sfilename + "\n");
	}
	soutput << "Chrom\tStartProbe\tStartBPLOC\tEndProbe\tEndBPLOC\tLength\tNumSnps\n";
	for(int i = 0; i < (int)commonroh.size(); i++){
		Marker* startmark = enabledmarkers[commonroh[i++]];
		Marker* endmark = enabledmarkers[commonroh[i]];
		int nummarks = commonroh[i] - commonroh[i - 1];
		int difference = endmark->getBPLOC() - startmark->getBPLOC();
		int controlcount = 0;
		int casecount = 0;
		int casetotal = 0;
		int controltotal = 0;
		int casemajcount = 0;
		int casemincount = 0;
		int controlmajcount = 0;
		int controlmincount = 0;
		int casestrict = 0;
		int controlstrict = 0;

		vector<Sample*> caselist;
		for(int s = 0; s < (int)samples->size(); s++){
			Sample* samp = (*samples)[s];
			if(samp->isEnabled()){
				if(samp->getSex() && startmark->getChrom() == opts::_CHRX_){
					continue;
				}
				if(samp->getPheno() == 2){
					casetotal++;
				}
				else{
					controltotal++;
				}
				for(int m = commonroh[i - 1]; m < commonroh[i]; m++){
					if(samp->getRohVal(m)){
						if(samp->getPheno() == 2){
							casecount++;
							caselist.push_back(samp);
						}
						else{
							controlcount++;
						}
						break;
					}
				}
				int mincount = 0;
				int majcount = 0;
				int strict = 0;
				for(int m = commonroh[i-1]; m < commonroh[i]; m++){
					if(samp->getRohVal(m)){
						Marker* mymark = enabledmarkers[m];
						int myloc = mymark->getLoc();
						if(!samp->getAone(myloc) && !samp->getAtwo(myloc)){
							mincount++;
							strict++;
						}
						else if(samp->getAone(myloc) && samp->getAtwo(myloc) && !samp->getAmissing(myloc)){
							majcount++;
							strict++;
						}
						else if(samp->getAone(myloc) && samp->getAtwo(myloc) && samp->getAmissing(myloc)){
							strict++;
						}
					}
				}
				if(samp->getPheno() == 2){
					casemincount += mincount;
					casemajcount += majcount;
				}
				else{
					controlmincount += mincount;
					controlmajcount += majcount;
				}
				if(strict == nummarks){
					if(samp->getPheno() == 2){
						casestrict++;
					}
					else{
						controlstrict++;
					}
				}
			}
		}
		ChiSquareAllelic csa;
		vector<vector<int> > chitotals;
		chitotals.resize(2);
		chitotals[1].push_back(casetotal - casecount);
		chitotals[1].push_back(casecount);
		chitotals[0].push_back(controltotal - controlcount);
		chitotals[0].push_back(controlcount);
		double chi1 = csa.chisquare(chitotals);
		vector<vector<double> > expected = csa.expecteds(chitotals);

		double results, df = 1;
		results = Helpers::p_from_chi(chi1, df);

		output << startmark->getChrom() << "\t" << startmark->getProbeID() << "\t" << startmark->getBPLOC() << "\t" << endmark->getProbeID() << "\t" << endmark->getBPLOC() << "\t" << difference << "\t" << nummarks << "\t" << casecount << ":" << (casetotal - casecount) << " (" << expected[1][1] << ":" << expected[1][0] << ")\t" << controlcount << ":" << (controltotal - controlcount) << " (" << expected[0][1] << ":" << expected[0][0] << ")\t" << casemincount << "\t" << casemajcount << "\t" << controlmincount << "\t" << controlmajcount << "\t" << casestrict << "\t" << controlstrict << "\t" << chi1 << "\t" << results << endl;
		soutput << startmark->getChrom() << "\t" << startmark->getProbeID() << "\t" << startmark->getBPLOC() << "\t" << endmark->getProbeID() << "\t" << endmark->getBPLOC() << "\t" << difference << "\t" << nummarks << "\t" << casecount << ":" << (casetotal - casecount) << " (" << expected[1][1] << ":" << expected[1][0] << ")\t" << controlcount << ":" << (controltotal - controlcount) << " (" << expected[0][1] << ":" << expected[0][0] << ")\t" << casemincount << "\t" << casemajcount << "\t" << controlmincount << "\t" << controlmajcount << "\t" << casestrict << "\t" << controlstrict << "\t" << chi1 << "\t" << results << endl;
		soutput << "-------------------------------------------------------------------------\n";
		soutput << "FamID\tIndID\tGender\n";
		for(int l = 0; l < (int)caselist.size(); l++){
			Sample* samp = caselist[l];
			soutput << samp->getFamID() << "\t" << samp->getInd() << "\t";
			if(samp->getSex()){
				soutput << "M\n";
			}
			else{
				soutput << "F\n";
			}
		}
	}

	if(output.is_open()){
		output.close();
	}
	if(soutput.is_open()){
		soutput.close();
	}


	//exploratory analysis
	string filename2 = opts::_OUTPREFIX_ + "homozygous_wgha_exploratory" + options.getOut() + ".txt";
	if(!overwrite){
		filename2 += "." + getString<int>(order);
	}

	filenames.push_back(filename2);

	ofstream outputexp (filename2.c_str());
	if(!outputexp.is_open()){
		opts::printLog("Unable to open " + filename2 + "\n");
		//exit(1);
		throw MethodException("Unable to open " + filename2 + "\n");
	}
	outputexp << "Chrom\tProbe\tBPLOC\tpval\n";
	count = 0;
	start = -1;
	end = -1;
	prev_chrom = -1;
	vector<double> mpvals;
	mpvals.resize(enabledmarkers.size());
	for(int m = 0; m < (int)enabledmarkers.size(); m++){
		Marker* mark = enabledmarkers[m];
		if(prev_chrom == -1){
			prev_chrom = mark->getChrom();
		}

		if(prev_chrom != mark->getChrom()){
			if(count > (options.getHomozygSpan() / 2)){
				outputexp << "------------------------------------------BREAK-------------------------------" << endl;
				for(int k = start; k <= end; k++){
					Marker* temp = enabledmarkers[k];
					outputexp << temp->getChrom() << "\t" << temp->getProbeID() << "\t" << temp->getBPLOC() << "\t" << mpvals[k] << endl;
				}
			}
			start = -1;
			end = -1;
			count = 0;
		}
		prev_chrom = mark->getChrom();
		int controlcount = 0;
		int casecount = 0;
		int casetotal = 0;
		int controltotal = 0;
		for(int s = 0; s < (int)samples->size(); s++){
			Sample* samp = (*samples)[s];
			if(samp->isEnabled()){
				if(samp->getSex() && mark->getChrom() == opts::_CHRX_){
					continue;
				}
				if(samp->getPheno() == 2){
					casetotal++;
				}
				else{
					controltotal++;
				}
				if(samp->getRohVal(m)){
					if(samp->getPheno() == 2){
						casecount++;
					}
					else{
						controlcount++;
					}
				}
			}
		}
		ChiSquareAllelic csa;
		vector<vector<int> > chitotals;
		chitotals.resize(2);
		chitotals[1].push_back(casetotal - casecount);
		chitotals[1].push_back(casecount);
		chitotals[0].push_back(controltotal - controlcount);
		chitotals[0].push_back(controlcount);
		double chi1 = csa.chisquare(chitotals);
		vector<vector<double> > expected = csa.expecteds(chitotals);

		double results, df = 1;
		results = Helpers::p_from_chi(chi1, df);
		mpvals[m] = results;
		if(results < 0.01){
			if(start == -1){
				start = m;
			}
			end = m;
			count++;
		}
		else{
			if(count > (options.getHomozygSpan() / 2)){
				outputexp << "------------------------------------------BREAK-------------------------------" << endl;
				for(int k = start; k <= end; k++){
					Marker* temp = enabledmarkers[k];
					outputexp << temp->getChrom() << "\t" << temp->getProbeID() << "\t" << temp->getBPLOC() << "\t" << mpvals[k] << endl;
				}
			}
			count = 0;
			start = -1;
			end = -1;
		}

	}
	if(outputexp.is_open()){
		outputexp.close();
	}
	if(coutput.is_open()){
		coutput.close();
	}

	for(int s = 0; s < (int)samples->size(); s++){
		Sample* test = (*samples)[s];
		test->resizeRoh(0);
	}

}


Sample* Homozygous::findRandomSample(map<Sample*, bool> rsamps){
	Sample* samp = NULL;

	while(samp == NULL){
		int random;
		random = rand() % samples->size();
		samp = (*samples)[random];
		if(!samp->isEnabled()){
			samp = NULL;
			continue;
		}
		if(rsamps[samp]){
			samp = NULL;
			continue;
		}
	}
	return samp;

}

void Homozygous::perform_homozyg_permutations(){
	int msize = good_markers.size();
	map<Sample*, bool> change_samples;

	string filename = opts::_OUTPREFIX_ + "homozygous_permute" + options.getOut() + ".txt";
	if(!overwrite){
		filename += "." + getString<int>(order);
	}

	filenames.push_back(filename);

	ofstream perm (filename.c_str());
	if(!perm.is_open()){
		opts::printLog("Unable to open " + filename + "\n");
		//exit(1);
		throw MethodException("Unable to open " + filename + "\n");
	}
	opts::addFile("Marker", stepname, filename);
	perm << "Chrom\trsID\tProbeID\tbploc\tObserved_Test_Stat";

	for(int i = 0; i < options.getHomozygPermutationCount(); i++){
		perm << "\tIter" << (i + 1);
	}


	perm << "\tPercentile_in_permuted\n";
	opts::addHeader(filename, "Observed_non-chisq");
	opts::addHeader(filename, "Percentile_in_permuted");

	vector<int> orig_samp_pheno(samples->size());
	for(int s = 0; s < (int)samples->size(); s++){
		orig_samp_pheno[s] = (*samples)[s]->getPheno();
	}
	AlleleFrequency* af = new AlleleFrequency(samples, families);
	af->setRank(rank);

	for (int m = 0; m < msize; m++) {
		Marker* hmark = good_markers[m];
		if (hmark->isEnabled()){
			int mloc = hmark->getLoc();
			af->calcOne(hmark);
			double obsaone = af->getAoneCon_freq();
			int obstotalcases = af->getAonehomoCa() + af->getAtwohomoCa() + af->getHetCa();
			double obsnq2 = ((double)obstotalcases * (obsaone * obsaone));
			double obsresult = ((af->getAonehomoCa() - obsnq2) * (af->getAonehomoCa() - obsnq2)) / obsnq2;

			perm << hmark->toString() << "\t" << obsresult;
			vector<double> marker_results;
			marker_results.resize(options.getHomozygPermutationCount());

			vector<bool> orig_a1(samples->size());
			vector<bool> orig_a2(samples->size());
			vector<bool> orig_aM(samples->size());
			for(unsigned int s = 0; s < samples->size(); s++){
				orig_a1[s] = (*samples)[s]->getAone(mloc);
				orig_a2[s] = (*samples)[s]->getAtwo(mloc);
				orig_aM[s] = (*samples)[s]->getAmissing(mloc);
			}
			for (int i = 0; i < options.getHomozygPermutationCount(); i++) {
				change_samples.clear();
				for (int s = 0; s < (int)samples->size(); s++) {
					Sample* samp = (*samples)[s];
					if (samp->isEnabled()) {
						//generate random genotype between 0-2
						int geno = rand() % 3;
						if(geno == 0){
							samp->addAone(mloc, false);
							samp->addAtwo(mloc, false);
							samp->addAmissing(mloc, false);
						}
						else if(geno == 1){
							samp->addAone(mloc, false);
							samp->addAtwo(mloc, true);
							samp->addAmissing(mloc, false);
						}
						else if(geno == 2){
							samp->addAone(mloc, true);
							samp->addAtwo(mloc, true);
							samp->addAmissing(mloc, false);
						}
						else{
							cout << "ACK!!!\n";
						}
					}
				}
				AlleleFrequency* paf = new AlleleFrequency(samples, families);
				paf->setRank(0);
				paf->calcOne(hmark);
				//perform calculation for permutation
				double aone = paf->getAoneCon_freq();
				int totalcases = paf->getAonehomoCa() + paf->getAtwohomoCa() + paf->getHetCa();
				double nq2 = ((double)totalcases * (aone * aone));
				double result = ((paf->getAonehomoCa() - nq2) * (paf->getAonehomoCa() - nq2)) / nq2;
				marker_results[i] = result;

				perm << "\t" << result;

				delete(paf);
			}

			//reinstate original gneotype for this snp
			for(unsigned int s = 0; s < samples->size(); s++){
				(*samples)[s]->addAone(mloc, orig_a1[s]);
				(*samples)[s]->addAtwo(mloc, orig_a2[s]);
				(*samples)[s]->addAmissing(mloc, orig_aM[s]);
			}

			//calculate percentile
			map<double, bool> unique;
			for(int u = 0; u < (int)marker_results.size(); u++){
				unique[marker_results[u]] = true;
			}
			map<double, bool>::iterator uiter;
			vector<double> tempdouble;
			for(uiter = unique.begin(); uiter != unique.end(); uiter++){
				tempdouble.push_back(uiter->first);
			}
			sort(tempdouble.begin(), tempdouble.end());
			int countless = 0;
			for(int t = 0; t < (int)tempdouble.size(); t++){
				if(Helpers::dGreater(obsresult, tempdouble[t])){
					countless++;
					continue;
				}
				else{
					break;
				}
			}
			countless = tempdouble.size() - countless;
			double percentile = ((double)countless / (double)tempdouble.size());

			perm << "\t" << percentile << endl;
		}
	}
	delete(af);
}

void Homozygous::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
	filenames.clear();

	good_markers = Helpers::findValidMarkers(markers, &options);

	if(options.doHomozygPermute()){
		perform_homozyg_permutations();
		return;
	}
	if(options.doHomozygRaw()){
		int msize = markers->size();
		msize = good_markers.size();

		homoaffcount.resize(msize,0);
		homounaffcount.resize(msize,0);
		homoallcount.resize(msize,0);
		homominallcount.resize(msize, 0);
		homomajallcount.resize(msize, 0);

		for(int m = 0; m < msize; m++){
			Marker* hmark = good_markers[m];
			if(hmark->isEnabled()){
				int mloc = hmark->getLoc();
				for(int s = 0; s < (int)samples->size(); s++){
					Sample* samp = (*samples)[s];
					if(samp->isEnabled()){
						if(samp->getPheno() == 2){
							if(!samp->getAone(mloc) && !samp->getAtwo(mloc)){
								homoaffcount[m]++;
								homominallcount[m]++;
								homoallcount[m]++;
							}
							else if(samp->getAone(mloc) && samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
								homoaffcount[m]++;
								homomajallcount[m]++;
								homoallcount[m]++;
							}
						}
						else{
							if(!samp->getAone(mloc) && !samp->getAtwo(mloc)){
								homounaffcount[m]++;
								homominallcount[m]++;
								homoallcount[m]++;
							}
							else if(samp->getAone(mloc) && samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
								homounaffcount[m]++;
								homomajallcount[m]++;
								homoallcount[m]++;
							}
						}
					}
				}
			}
			if(options.doHomozygPermute()){
				perform_homozyg_permutations();
			}
		}
		return;
	}
	if(options.doHomozygWGHA()){
		process_wgha();
		return;
	}
	vector<Sample*>::iterator s_iter;

	int msize = good_markers.size();
	homoaffcount.resize(msize,0);
	homounaffcount.resize(msize,0);
	homoallcount.resize(msize,0);
	homominallcount.resize(msize, 0);
	homomajallcount.resize(msize, 0);

	string filename = opts::_OUTPREFIX_ + "homozygous_plot" + options.getOut() + ".txt";
	if(!overwrite){
		filename += "." + getString<int>(order);
	}

	filenames.push_back(filename);

	ofstream plot (filename.c_str());
	if(!plot.is_open()){
		opts::printLog("Unable to open " + filename + "\n");
		//exit(1);
		throw MethodException("Unable to open " + filename + "\n");
	}
	plot << "ID1\tID\tFamID\tIndID\tChrom\tbploc\n";
	filename = opts::_OUTPREFIX_ + "homozygous" + options.getOut() + ".txt";
	if(!overwrite){
		filename += "." + getString<int>(order);
	}
	filenames.push_back(filename);
	ofstream homo (filename.c_str());
	if(!homo.is_open()){
		opts::printLog("Unable to open " + filename + "\n");
		//exit(1);
		throw MethodException("Unable to open " + filename + "\n");
	}
	homo << "FamID\tIndID\tChrom\tStart_ProbeID\tEnd_ProbeID\tStart_bploc\tEnd_bploc\tBP_range\tNum_markers\tNum_founders\tNum_children\tPercent_Children\tProbability\n";
	filename = opts::_OUTPREFIX_ + "homozygous_zero" + options.getOut() + ".txt";
	if(!overwrite){
		filename += "." + getString<int>(order);
	}
	filenames.push_back(filename);
	ofstream zero (filename.c_str());
	if(!zero.is_open()){
		opts::printLog("Unable to open " + filename + "\n");
		//exit(1);
		throw MethodException("Unable to open " + filename + "\n");
	}
	zero << "FamID\tIndID\tChrom\tStart_ProbeID\tEnd_ProbeID\tStart_bploc\tEnd_bploc\tBP_range\tNum_markers\n";
	int scount = 0;
	int total = 0;
	for(s_iter = samples->begin(); s_iter != samples->end(); s_iter++){
		if((*s_iter)->isEnabled()){
			scount++;
			Marker* markstarthomo = NULL;
			Marker* markendhomo = NULL;
			Marker* markstartzero = NULL;
			Marker* markendzero = NULL;
			int locstarthomo = -1;
			int locendhomo = -1;
			int locstartzero = -1;
			int locendzero = -1;
			int homospan = 0;
			int zerospan = 0;
			int zerocount = 0;
			int currchrom = -1;
			vector<int> zerolocs;
			for(int i = 0; i < msize; i++){
				int loc = good_markers[i]->getLoc();
				if(good_markers[i]->isEnabled()){
					if(currchrom == -1){
						currchrom = good_markers[i]->getChrom();
					}

					if(currchrom != good_markers[i]->getChrom()){
						//output data
						if(markstarthomo != NULL && markendhomo != NULL && locstarthomo > -1 && locendhomo > -1 && homospan > options.getHomozygSpan()){
							double prob = 0.0;
							int numfounders = 0;
							int numchildren = 0;
							double perinseq = 0.0;
							prob = checkSequence(locstarthomo, locendhomo, numfounders, numchildren, perinseq);
							bool goahead = true;
							if(options.doHomozygSeqTest()){
								if(Helpers::dGreater(prob, options.getHomozygSeqVal())){
									goahead = false;
								}
							}
							if(goahead){
								homo << (*s_iter)->getFamID() << "\t" << (*s_iter)->getInd() << "\t" << currchrom << "\t" << markstarthomo->getProbeID() << "\t" << markendhomo->getProbeID() << "\t" << markstarthomo->getBPLOC() << "\t" << markendhomo->getBPLOC() << "\t" << (markendhomo->getBPLOC() - markstarthomo->getBPLOC()) << "\t" << homospan << "\t" << numfounders << "\t" << numchildren << "\t" << perinseq << "\t" << prob << endl;
								for(int k = locstarthomo; k <= locendhomo; k++){
									if(good_markers[k]->isEnabled()){
										vector<int>::iterator found = find(zerolocs.begin(), zerolocs.end(), k);
										if(found == zerolocs.end()){
											plot << ++total << "\t" << scount << "\t" << (*s_iter)->getFamID() << "\t" << (*s_iter)->getInd() << "\t" << currchrom << "\t" << good_markers[k]->getBPLOC() << endl;
											Marker* tempmarker = good_markers[k];
											if((*s_iter)->getAone(tempmarker->getLoc()) && (*s_iter)->getAtwo(tempmarker->getLoc()) && !(*s_iter)->getAmissing(tempmarker->getLoc())){
												homomajallcount[k]++;
											}
											else if(!(*s_iter)->getAone(tempmarker->getLoc()) && !(*s_iter)->getAtwo(tempmarker->getLoc())){
												homominallcount[k]++;
											}
											homoallcount[k]++;
											if((*s_iter)->getPheno() == 2){
												homoaffcount[k]++;
											}
											else{
												homounaffcount[k]++;
											}
										}
									}
								}
							}
						}
						currchrom = good_markers[i]->getChrom();
						homospan = 0;
						markstarthomo = NULL;
						markendhomo = NULL;
						locstarthomo = -1;
						locendhomo = -1;
						locstartzero = -1;
						locendzero = -1;
						markstartzero = NULL;
						markendzero = NULL;
						zerospan = 0;
						zerocount = 0;
						zerolocs.clear();
					}
					if(good_markers[i]->getChrom() == opts::_CHRX_ && (*s_iter)->getSex()){
						homospan = 0;
						markstarthomo = NULL;
						markendhomo = NULL;
						locstarthomo = -1;
						locendhomo = -1;
						locstartzero = -1;
						locendzero = -1;
						markstartzero = NULL;
						markendzero = NULL;
						zerospan = 0;
						zerocount = 0;
						zerolocs.clear();
						continue;
					}
						if((!good_markers[i]->isMicroSat() && (((*s_iter)->getAone(loc) && (*s_iter)->getAtwo(loc) && !(*s_iter)->getAmissing(loc)) || (!(*s_iter)->getAone(loc) && !(*s_iter)->getAtwo(loc)))) ||
								(good_markers[i]->isMicroSat() && (((*s_iter)->getAbone(loc) == (*s_iter)->getAbtwo(loc))))){
							if(markstarthomo == NULL){
								markstarthomo = (*markers)[i];
								locstarthomo = i;
							}
							homospan++;
							markendhomo = (*markers)[i];
							locendhomo = i;
							zerospan = 0;
							markstartzero = NULL;
							markendzero = NULL;
						}
						else if((!good_markers[i]->isMicroSat() && (*s_iter)->getAone(loc) && (*s_iter)->getAtwo(loc) && (*s_iter)->getAmissing(loc)) ||
								(good_markers[i]->isMicroSat() && (*s_iter)->getAbone(loc) == -1 && (*s_iter)->getAbtwo(loc) == -1)
								){
							zerocount++;
							zerolocs.push_back(i);
							if(options.getHomozygZeros() > -1 && zerocount >= options.getHomozygZeros()){
								if(markstarthomo != NULL && markendhomo != NULL && locstarthomo > -1 && locendhomo > -1 && homospan > options.getHomozygSpan()){
									double prob = 0.0;
							int numfounders = 0;
							int numchildren = 0;
							double perinseq = 0.0;
							prob = checkSequence(locstarthomo, locendhomo, numfounders, numchildren, perinseq);
									bool goahead = true;
									if(options.doHomozygSeqTest()){
										if(Helpers::dGreater(prob, options.getHomozygSeqVal())){
											goahead = false;
										}
									}
									if(goahead){
								homo << (*s_iter)->getFamID() << "\t" << (*s_iter)->getInd() << "\t" << currchrom << "\t" << markstarthomo->getProbeID() << "\t" << markendhomo->getProbeID() << "\t" << markstarthomo->getBPLOC() << "\t" << markendhomo->getBPLOC() << "\t" << (markendhomo->getBPLOC() - markstarthomo->getBPLOC()) << "\t" << homospan << "\t" << numfounders << "\t" << numchildren << "\t" << perinseq << "\t" << prob << endl;
										for(int k = locstarthomo; k <= locendhomo; k++){
											if((*markers)[k]->isEnabled()){
												vector<int>::iterator found = find(zerolocs.begin(), zerolocs.end(), k);
												if(found == zerolocs.end()){
													Marker* tempmarker = good_markers[k];
													if((*s_iter)->getAone(tempmarker->getLoc()) && (*s_iter)->getAtwo(tempmarker->getLoc()) && !(*s_iter)->getAmissing(tempmarker->getLoc())){
														homomajallcount[k]++;
													}
													else if(!(*s_iter)->getAone(tempmarker->getLoc()) && !(*s_iter)->getAtwo(tempmarker->getLoc())){
														homominallcount[k]++;
													}
													homoallcount[k]++;
													plot << ++total << "\t" << scount << "\t" << (*s_iter)->getFamID() << "\t" << (*s_iter)->getInd() << "\t" << currchrom << "\t" << good_markers[k]->getBPLOC() << endl;
													if((*s_iter)->getPheno() == 2){
														homoaffcount[k]++;
													}
													else{
														homounaffcount[k]++;
													}
												}
											}
										}
									}
									markstarthomo = NULL;
									markendhomo = NULL;
									locstarthomo = -1;
									locendhomo = -1;
									homospan = 0;
								 	zerolocs.clear();
								}
							}
							if(markstartzero == NULL){
								markstartzero = good_markers[i];
								locstartzero = i;
							}
							zerospan++;
							markendzero = good_markers[i];
							locendzero = i;
						}
						else{
							if(markstarthomo != NULL && markendhomo != NULL && locstarthomo > -1 && locendhomo > -1 && homospan > options.getHomozygSpan()){
								//output data
								double prob = 0.0;
							int numfounders = 0;
							int numchildren = 0;
							double perinseq = 0.0;
							prob = checkSequence(locstarthomo, locendhomo, numfounders, numchildren, perinseq);
								bool goahead = true;
								if(options.doHomozygSeqTest()){
									if(Helpers::dGreater(prob, options.getHomozygSeqVal())){
										goahead = false;
									}
								}
								if(goahead){
								homo << (*s_iter)->getFamID() << "\t" << (*s_iter)->getInd() << "\t" << currchrom << "\t" << markstarthomo->getProbeID() << "\t" << markendhomo->getProbeID() << "\t" << markstarthomo->getBPLOC() << "\t" << markendhomo->getBPLOC() << "\t" << (markendhomo->getBPLOC() - markstarthomo->getBPLOC()) << "\t" << homospan << "\t" << numfounders << "\t" << numchildren << "\t" << perinseq << "\t" << prob << endl;
									for(int k = locstarthomo; k <= locendhomo; k++){
										if(good_markers[k]->isEnabled()){
											vector<int>::iterator found = find(zerolocs.begin(), zerolocs.end(), k);
											if(found == zerolocs.end()){
												Marker* tempmarker = good_markers[k];
												if((*s_iter)->getAone(tempmarker->getLoc()) && (*s_iter)->getAtwo(tempmarker->getLoc()) && !(*s_iter)->getAmissing(tempmarker->getLoc())){
													homomajallcount[k]++;
												}
												else if(!(*s_iter)->getAone(tempmarker->getLoc()) && !(*s_iter)->getAtwo(tempmarker->getLoc())){
													homominallcount[k]++;
												}
												homoallcount[k]++;
												plot << ++total << "\t" << scount << "\t" << (*s_iter)->getFamID() << "\t" << (*s_iter)->getInd() << "\t" << currchrom << "\t" << good_markers[k]->getBPLOC() << endl;
												if((*s_iter)->getPheno() == 2){
													homoaffcount[k]++;
												}
												else{
													homounaffcount[k]++;
												}
											}
										}
									}
								}
							}
							markstarthomo = NULL;
							markendhomo = NULL;
							zerocount = 0;
							markstartzero = NULL;
							markendzero = NULL;
							locstarthomo = -1;
							locendhomo = -1;
							locstartzero = -1;
							locendzero = -1;
							zerospan = 0;
							homospan = 0;
							zerolocs.clear();
						}
				}
			}
			if(markstarthomo != NULL && markendhomo != NULL && locstarthomo > -1 && locendhomo > -1 && homospan > options.getHomozygSpan()){
				//output data
				if(markstarthomo->getChrom() == opts::_CHRX_ && (*s_iter)->getSex()){
					continue;
				}
				double prob = 0.0;
							int numfounders = 0;
							int numchildren = 0;
							double perinseq = 0.0;
							prob = checkSequence(locstarthomo, locendhomo, numfounders, numchildren, perinseq);
				bool goahead = true;
				if(options.doHomozygSeqTest()){
					if(Helpers::dGreater(prob, options.getHomozygSeqVal())){
						goahead = false;
					}
				}
				if(goahead){
								homo << (*s_iter)->getFamID() << "\t" << (*s_iter)->getInd() << "\t" << currchrom << "\t" << markstarthomo->getProbeID() << "\t" << markendhomo->getProbeID() << "\t" << markstarthomo->getBPLOC() << "\t" << markendhomo->getBPLOC() << "\t" << (markendhomo->getBPLOC() - markstarthomo->getBPLOC()) << "\t" << homospan << "\t" << numfounders << "\t" << numchildren << "\t" << perinseq << "\t" << prob << endl;
					for(int k = locstarthomo; k <= locendhomo; k++){
						if(good_markers[k]->isEnabled()){
							vector<int>::iterator found = find(zerolocs.begin(), zerolocs.end(), k);
							if(found == zerolocs.end()){
								plot << ++total << "\t" << scount << "\t" << (*s_iter)->getFamID() << "\t" << (*s_iter)->getInd() << "\t" << currchrom << "\t" << good_markers[k]->getBPLOC() << endl;
								Marker* tempmarker = good_markers[k];
								if((*s_iter)->getAone(tempmarker->getLoc()) && (*s_iter)->getAtwo(tempmarker->getLoc()) && !(*s_iter)->getAmissing(tempmarker->getLoc())){
									homomajallcount[k]++;
								}
								else if(!(*s_iter)->getAone(tempmarker->getLoc()) && !(*s_iter)->getAtwo(tempmarker->getLoc())){
									homominallcount[k]++;
								}
								homoallcount[k]++;
								if((*s_iter)->getPheno() == 2){
									homoaffcount[k]++;
								}
								else{
									homounaffcount[k]++;
								}
							}
						}
					}
				}
			}
		}
	}

	if(homo.is_open()){
		homo.close();
	}
	if(zero.is_open()){
		zero.close();
	}
	if(plot.is_open()){
		plot.close();
	}
}
}
