/**********************************************************************************
*                       Concordance Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Performs concordance check with another set of ped files.
*
*
* Files generated:
*
*File: Concordance.cc
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
#include "Concordance.h"
#include "Options.h"
#include "General.h"
#include "Helper.h"
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
namespace Methods{
string Concordance::stepname = "concordance";

void Concordance::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void Concordance::PrintSummary(){
//	string filename = opts::_OUTPREFIX_ + "marker_geno_eff" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
//	string filenameg = opts::_OUTPREFIX_ + "marker_geno_eff_groups" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
//	if(!overwrite){
//		filename += "." + getString<int>(order);
//		filenameg += "." + getString<int>(order);
//	}

//	ofstream bymarker (filename.c_str());
//	if(!bymarker.is_open()){
//		opts::printLog("Unable to open " + filename + "\n");
//		exit(1);
//	}

//	opts::addFile("Marker", stepname, filename);
//	bymarker.precision(4);
//	if((*markers)[0]->getDetailHeaders().size() > 0){
//		bymarker << "Chrom\trsid\tProbeID\tbploc\t" << (*markers)[0]->getDetailHeaders() << "\t%GenoEff_Ind_All\tInd_Zero_Count\tTotal_Individuals_Used\t%GenoEff_Ind_Cases\tInd_Zero_Count_Cases\tTotal_Individuals_Used_Cases\t%GenoEff_Ind_Controls\tInd_Zero_Count_Controls\tTotal_Individuals_Used_Controls";
//	}
//	else{
//		bymarker << "Chrom\trsid\tProbeID\tbploc\t%GenoEff_Ind_All\tInd_Zero_Count\tTotal_Individuals_Used\t%GenoEff_Ind_Cases\tInd_Zero_Count_Cases\tTotal_Individuals_Used_Cases\t%GenoEff_Ind_Controls\tInd_Zero_Count_Controls\tTotal_Individuals_Used_Controls";
//	}
//	opts::addHeader(filename, "%GenoEff_Ind_All");
//	opts::addHeader(filename, "Ind_Zero_Count");
//	opts::addHeader(filename, "Total_Individuals_Used");
//	opts::addHeader(filename, "%GenoEff_Ind_Cases");
//	opts::addHeader(filename, "Ind_Zero_Count_Cases");
//	opts::addHeader(filename, "Total_Individuals_Used_Cases");
//	opts::addHeader(filename, "%GenoEff_Ind_Controls");
//	opts::addHeader(filename, "Ind_Zero_Count_Controls");
//	opts::addHeader(filename, "Total_Individuals_Used_Controls");
//	bymarker << endl;

}

void Concordance::filter(){
}

void Concordance::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
	calculate();
}

void Concordance::calculate(){
	if(options.doPedFile() && options.doMapFile()){
		opts::printLog("Reading concordance MAP file: " + options.getMapFile() + "\n");
		readMapM(check_data_set, options);
		opts::printLog("Reading concordance PED file: " + options.getPedFile() + "\n");
		readPedM_3vec_set(check_data_set, options);
//		assignLinks(&check_families);
//		reorderAlleles(&check_samples, &check_markers);
	}
	else if(options.doTPedFile() && options.doTFamFile()){
		opts::printLog("Reading concordance family Map file: " + options.getTFamFile() + "\n");
		readTFamM(check_data_set, options);
		opts::printLog("Reading concordance TPED File: " + options.getTPedFile() + "\n");
		readTPedM(check_data_set, options);
//		assignLinks(&check_families);
//		reorderAlleles(&check_samples, &check_markers);
	}
	else if(options.doBinInput()){
		opts::printLog("Reading concordance data using Binary inputs: Prefix: " + options.getBinInput() + "\n");
		readBinM(check_data_set, options);
//		assignLinks(&check_families);
//		reorderAlleles(&check_samples, &check_markers);
	}
	else{
		opts::printLog("No input specified for " + stepname + "\n");
		exit(1);
	}

	string sfile = opts::_OUTPREFIX_ + "concordance_samples" + options.getOut() + ".txt";
	sample_error_file = sfile;
	string mfile = opts::_OUTPREFIX_ + "concordance" + options.getOut() + ".txt";
	main_file = mfile;
	string efile = opts::_OUTPREFIX_ + "concordance_errors" + options.getOut() + ".txt";
	error_file = efile;

	if(!overwrite){
		sfile += "." + getString<int>(order);
		mfile += "." + getString<int>(order);
		efile += "." + getString<int>(order);
	}

	ofstream samp_data (sfile.c_str());
	ofstream snp_data(mfile.c_str());
	ofstream errors(efile.c_str());

	if(!samp_data.is_open()){
		opts::printLog("Unable to open " + sfile + "\n");
		exit(1);
	}
	if(!snp_data.is_open()){
		opts::printLog("Unable to open " + mfile + "\n");
		exit(1);
	}
	if(!errors.is_open()){
		opts::printLog("Unable to open " + efile + "\n");
		exit(1);
	}

	errors << "FamID\tIndID\t";
	if(options.getUniqueId()){
		errors << "Unique_orig\tUnique_Test\t";
	}
	errors << "Chrom\tProbeID\tbploc\tOriginal_Genotype\tNew_Genotype" << endl;

	opts::addFile("Marker", stepname, mfile);
	opts::addFile("Sample", stepname, sfile);
	if((*markers)[0]->getDetailHeaders().size() > 0){
		snp_data << "Chrom\trsID\tProbeID\tbploc\t" << (*markers)[0]->getDetailHeaders() << "\tErrors\tTotal_Compared\t%Concordance";
	}
	else{
		snp_data << "Chrom\trsID\tProbeID\tbploc\tErrors\tTotal_Compared\t";
		if(options.getIncMissing()){
			snp_data << "Missing\t";
		}
		snp_data << "%Concordance";
	}
	opts::addHeader(mfile, "Errors");
	opts::addHeader(mfile, "Total_Compared");
	opts::addHeader(mfile, "%Concordance");
	snp_data << endl;

	samp_data << "FamID\tIndID\t";
	if(options.getUniqueId()){
		samp_data << "Unique_id_orig\tUnique_id_comp\t";
	}
	samp_data << "Center\tSex\tAffection_Status\tErrors\tTotal_Compared\t%Concordance";

	if(opts::_SAMPDESC_.length() > 0){
		samp_data << "\t" << (*samples)[0]->getDetailHeaders();
	}
	opts::addHeader(sfile, "Errors");
	opts::addHeader(sfile, "Total_Compared");
	opts::addHeader(sfile, "%Concordance");
	samp_data << endl;

	int commonsamps = 0;
	int commonsnps = 0;

	map<int,int> sorig_conc_map;
	map<int,int> morig_conc_map;
	map<int,int> snp_error_count;
	map<int,int> snp_zero_count;
	map<int,int> snp_total_count;

	for(int i = 0; i < check_data_set->num_inds(); i++){
		Sample* csamp = check_data_set->get_sample(i);
		vector<Sample*>::iterator itersamp = find_if(samples->begin(), samples->end(), FindSampleByFamAndID(csamp->getFamID(), csamp->getInd()));
		if(itersamp != samples->end()){
			int sloc = itersamp - samples->begin();
			sorig_conc_map[sloc] = i;
			commonsamps++;
		}

	}

	unsigned int orig_count = 0;
	for(unsigned int i = 0; i < check_data_set->num_loci(); i++){
		if(orig_count >= markers->size()){
			break;
		}
		Marker* mark = check_data_set->get_locus(i);
		Marker* orig_mark = (*markers)[orig_count];

		while(mark->getChrom() > orig_mark->getChrom() && orig_count < (markers->size() - 1)){
			orig_mark = (*markers)[++orig_count];
		}
		if(orig_count >=  markers->size()){
			break;
		}
		while(mark->getBPLOC() > orig_mark->getBPLOC() && orig_mark->getChrom() == mark->getChrom() && orig_count < (markers->size() - 1)){
			orig_mark = (*markers)[++orig_count];
		}

		if(orig_count >= markers->size()){
			break;
		}

		while(mark->getRSID() != orig_mark->getRSID() && mark->getChrom() == orig_mark->getChrom() && mark->getBPLOC() == orig_mark->getBPLOC() && orig_count < (markers->size() - 1)){
			orig_mark = (*markers)[++orig_count];
		}

//		vector<Marker*>::iterator itermark = find_if(markers->begin(), markers->end(), FindMarker(mark->getProbeID()));
//		if(itermark != markers->end()){
		if(mark->getChrom() == orig_mark->getChrom() && mark->getBPLOC() == orig_mark->getBPLOC()){
		//int mloc = itermark - markers->begin();
			Marker* omark = orig_mark;//(*markers)[mloc];
			int mloc = omark->getLoc();
			map<string, int> genomap;
			genomap[omark->getAllele1() + omark->getAllele2()] = 1;
			genomap[omark->getAllele2() + omark->getAllele1()] = 1;
			genomap[omark->getAllele1() + "0"] = 1;
			genomap["0"+omark->getAllele1()] = 1;
			genomap[omark->getAllele2() + "0"] = 1;
			genomap["0" + omark->getAllele2()] = 1;
			genomap["00"] = 1;
			map<string, int>::iterator giter = genomap.find(mark->getAllele1() + mark->getAllele2());
			if(giter == genomap.end()){
				opts::printLog("Strandedness issue found for snp: " + mark->toString() + "\n");
				opts::printLog("Original alleles found: " + omark->getAllele1() + "/" + omark->getAllele2() + "\n");
				opts::printLog("Checked alleles found: " + mark->getAllele1() + "/" + mark->getAllele2() + "\n");
				exit(1);
			}
			morig_conc_map[mloc] = i;
			snp_error_count[mloc] = 0;
			snp_zero_count[mloc] = 0;
			snp_total_count[mloc] = 0;
			commonsnps++;
		}
	}

	opts::printLog("Samples found: " + getString<int>(check_data_set->num_inds()) + ", " + getString<int>(commonsamps) + " in common with primary dataset.\n");
	opts::printLog("Markers found: " + getString<int>(check_data_set->num_loci()) + ", " + getString<int>(commonsnps) + " in common with primary dataset.\n");
	opts::printLog("Performing concordance check...");

	int totalchecked = 0;
	int numbad = 0;

	map<int, int>::iterator siter;
	map<int, int>::iterator miter;

	for(int i = 0; i < check_data_set->num_inds(); i++){//siter = sorig_conc_map.begin(); siter != sorig_conc_map.end(); siter++){
//		int orig = siter->first;
//		int second = siter->second;
		Sample* sconc = check_data_set->get_sample(i);//second);
		for(unsigned int j = 0; j < samples->size(); j++){
		Sample* sorig = (*samples)[j];//(*samples)[orig];
		if(sorig->toString() != sconc->toString()){//sorig->getFamID() != sconc->getFamID() && sorig->getInd() != sconc->getInd()){
			continue;
		}
		int sample_error_count = 0;
		int sample_total_count = 0;
		for(miter = morig_conc_map.begin(); miter != morig_conc_map.end(); miter++){
			int morig = miter->first;
			int msecond = miter->second;
			Marker* mmorig = (*markers)[morig];
			Marker* mconc = check_data_set->get_locus(msecond);
			int oloc = mmorig->getLoc();
			int cloc = mconc->getLoc();
			bool inc = false;
			bool both_zero = false;
			string oa1 = "";
			string oa2 = "";
			string ca1 = "";
			string ca2 = "";

			if(!options.getIncMissing()){
				if((sorig->getAone(oloc) && sorig->getAtwo(oloc) && sorig->getAmissing(oloc)) || (sconc->getAone(cloc) && sconc->getAtwo(cloc) && sconc->getAmissing(cloc))){
					continue;
				}
			}

			bool boa1 = sorig->getAone(oloc);
			bool boa2 = sorig->getAtwo(oloc);
			bool boa3 = sorig->getAmissing(oloc);
			if(boa1 && boa2 && boa3){
				oa1 = "0";
				oa2 = "0";
			}
			else{
				if(!sorig->getAone(oloc)){
					oa1 =  mmorig->getAllele1();
				}
				else{
					oa1 = mmorig->getAllele2();
				}
				if(!sorig->getAtwo(oloc)){
					oa2 = mmorig->getAllele1();
				}
				else{
					oa2 = mmorig->getAllele2();
				}
			}

			bool bca1 = sconc->getAone(cloc);
			bool bca2 = sconc->getAtwo(cloc);
			bool bca3 = sconc->getAmissing(cloc);

			if(bca1 && bca2 && bca3){
				ca1 = "0";
				ca2 = "0";
			}
			else{
				if(!sconc->getAone(cloc)){
					ca1 = mconc->getAllele1();
				}
				else{
					ca1 = mconc->getAllele2();
				}
				if(!sconc->getAtwo(cloc)){
					ca2 = mconc->getAllele1();
				}
				else{
					ca2 = mconc->getAllele2();
				}
			}
			if(ca1 == "0" && ca2 == "0" && oa1 == "0" && oa2 == "0"){
				both_zero = true;
				snp_zero_count[morig]++;
			}
			else if(options.getIncMissing() && ((ca1 == "0" && ca2 == "0") || (oa1 == "0" && oa2 == "0"))){
				inc = true;
			}
			else{
				if(ca1 == mmorig->getAllele1() && ca2 == mmorig->getAllele1()){ //c = 11
					if(! ((!boa1) && (!boa2))) inc = true;
				}
				else if(ca1 == mmorig->getAllele1() && ca2 == mmorig->getAllele2()){ // c = 12
					if(! ((!boa1) && boa2)) inc = true;
				}
				else if(ca1 == mmorig->getAllele2() && ca2 == mmorig->getAllele1()){ // c = 21
					if(! ((!boa1) && boa2)) inc = true;
				}
				else if(ca1 == mmorig->getAllele2() && ca2 == mmorig->getAllele2()){ // c = 22
					if(! ((boa1) && boa2)) inc = true;
				}
			}

			if(inc && !both_zero){
				numbad++;
				snp_error_count[morig]++;
				sample_error_count++;
				errors << sorig->toString() << "\t";
				if(options.getUniqueId()){
					errors << sorig->getDadID() << "\t" << sconc->getDadID() << "\t";
				}
				errors << mmorig->getChrom() << "\t" << mmorig->getProbeID() << "\t" << mmorig->getBPLOC() << "\t" << oa1 << "_" << oa2 << "\t" << ca1 << "_" << ca2 << endl;

			}
			snp_total_count[morig]++;
			sample_total_count++;
			if(!both_zero)
				totalchecked++;
		}
		float percent = (1 - (float)((float)sample_error_count/(float)sample_total_count)) * 100.0f;
		samp_data << sorig->toString() << "\t";
		if(options.getUniqueId()){
			samp_data << sorig->getDadID() << "\t" << sconc->getDadID() << "\t";
		}
		samp_data << sorig->getFamily()->getCenter() << "\t";
		if(sorig->getSex()){
			samp_data << "M\t";
		}
		else{
			samp_data << "F\t";
		}
		if(sorig->getPheno() == 2){
			samp_data << "Y\t";
		}
		else if(sorig->getPheno() == 1){
			samp_data << "N\t";
		}
		else{
			samp_data << "U\t";
		}
		samp_data << sample_error_count << "\t" << sample_total_count << "\t" << percent << endl;
		}
	}

	for(miter = morig_conc_map.begin(); miter != morig_conc_map.end(); miter++){
		int morig = miter->first;
		////int msecond = miter->second;
		Marker* mmorig = (*markers)[morig];
		int total = snp_total_count[morig];
		int err = snp_error_count[morig];
		int zero = snp_zero_count[morig];

		float percent = (1 - ((float)err/(float)total));//*100.0f;
		snp_data << mmorig->toString() << "\t" << err << "\t" << total << "\t";
		if(options.getIncMissing()){
			snp_data << zero << "\t";
			if(zero == total){
				snp_data << "NA";
			}
			else{
				snp_data << percent;
			}
		}
		else{
			snp_data << percent;
		}
	    snp_data << endl;
	}

	float percent = (1 - ((float)numbad / (float)totalchecked));// * 100.0f;

	opts::printLog("Total genotypes checked: " + getString<int>(totalchecked) + "\n");
	opts::printLog("Total mis-matched genotypes: " + getString<int>(numbad) + "\n");
	opts::printLog("Overall concordance: " + getString<float>(percent) + "\n");

	snp_data.close();
	samp_data.close();
	errors.close();
}
}
