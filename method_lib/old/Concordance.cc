/**********************************************************************************
*                       Concordance Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Performs concordance check with another set of genotype files.
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
#include "Concordance.h"
#include "Options.h"
#include "General.h"
#include "Helpers.h"
namespace Methods{
string Concordance::stepname = "concordance";

/*DEPRECATED
 *
 */
void Concordance::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

/*DEPRECATED
 *
 */
void Concordance::PrintSummary(){

}

/*DEPRECATED
 *
 */
void Concordance::filter(){
}

/*DEPRECATED
 *
 */
void Concordance::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
	calculate();
}

/*
 * Main method to get concordance going
 * Reads in files to check and checks for strand correctness
 * Performs check against matching snps and samples
 * By default, it ignores missing genotypes
 *
 */
void Concordance::calculate(){
	if(options.doPedFile() && options.doMapFile()){
		opts::printLog("Reading concordance MAP file: " + options.getMapFile() + "\n");
		Helpers::readMapM(check_data_set, options);
		opts::printLog("Reading concordance PED file: " + options.getPedFile() + "\n");
		Helpers::readPedM_3vec_set(check_data_set, this->options);
		Helpers::assignLinks(check_data_set->get_families());
		Helpers::reorderAlleles(check_data_set->get_samples(), check_data_set->get_markers());
	}
	else if(options.doTPedFile() && options.doTFamFile()){
		opts::printLog("Reading concordance family Map file: " + options.getTFamFile() + "\n");
		Helpers::readTFamM(check_data_set, options);
		opts::printLog("Reading concordance TPED File: " + options.getTPedFile() + "\n");
		Helpers::readTPedM(check_data_set, options);
		Helpers::assignLinks(check_data_set->get_families());
		Helpers::reorderAlleles(check_data_set->get_samples(), check_data_set->get_markers());
	}
	else if(options.doBinInput()){
		opts::printLog("Reading concordance data using Binary inputs: Prefix: " + options.getBinInput() + "\n");
		Helpers::readBinM(check_data_set, options);
		Helpers::assignLinks(check_data_set->get_families());
		Helpers::reorderAlleles(check_data_set->get_samples(), check_data_set->get_markers());
	}
	else{
		opts::printLog("No input specified for " + stepname + "\n");
		throw MethodException("No input specified for " + stepname + "\n");
	}

	string sfile = opts::_OUTPREFIX_ + "concordance_samples" + options.getOut() + ".txt";
	sample_error_file = sfile;
	string mfile = opts::_OUTPREFIX_ + "concordance" + options.getOut() + ".txt";
	main_file = mfile;
	string efile = opts::_OUTPREFIX_ + "concordance_errors" + options.getOut() + ".txt";
	error_file = efile;
	string mismatchFile = opts::_OUTPREFIX_ + "concordance_mismatched_snps" + options.getOut() + ".txt";
	mismatch_file = mismatchFile;

	if(!overwrite){
		sfile += "." + getString<int>(order);
		mfile += "." + getString<int>(order);
		efile += "." + getString<int>(order);
		mismatchFile += "." + getString<int>(order);
	}

	ofstream samp_data (sfile.c_str());
	ofstream snp_data(mfile.c_str());
	ofstream errors(efile.c_str());
	ofstream mismatches(mismatchFile.c_str());

	if(!samp_data.is_open()){
		opts::printLog("Unable to open " + sfile + "\n");
		throw MethodException("Unable to open " + sfile + "\n");
	}
	if(!snp_data.is_open()){
		opts::printLog("Unable to open " + mfile + "\n");
		throw MethodException("Unable to open " + mfile + "\n");
	}
	if(!errors.is_open()){
		opts::printLog("Unable to open " + efile + "\n");
		throw MethodException("Unable to open " + efile + "\n");
	}
	if(!mismatches.is_open())
	{
		opts::printLog("Unable to open " + mismatchFile + "\n");
		throw MethodException("Unable to open " + mismatchFile + "\n");
	}

	errors << "FamID\tIndID\t";
	if(options.getUniqueId()){
		errors << "Unique_orig\tUnique_Test\t";
	}
	errors << "Chrom\tProbeID\tbploc\tOriginal_Genotype\tNew_Genotype" << endl;

	mismatches << "Chrom\trsID\tProbeID\tbploc" << endl;

	opts::addFile("Marker", stepname, mfile);
	opts::addFile("Sample", stepname, sfile);
	if((*markers).at(0)->getDetailHeaders().size() > 0){
		snp_data << "Chrom\trsID\tProbeID\tbploc\t" << (*markers).at(0)->getDetailHeaders() << "\tErrors\tTotal_Compared\t%Concordance";
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
		samp_data << "\t" << (*samples).at(0)->getDetailHeaders();
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

	//check for common samples
	for(int i = 0; i < check_data_set->num_inds(); i++){
		Sample* csamp = check_data_set->get_sample(i);
		vector<Sample*>::iterator itersamp = find_if(samples->begin(), samples->end(), FindSampleByFamAndID(csamp->getFamID(), csamp->getInd()));
		if(itersamp != samples->end()){
			int sloc = itersamp - samples->begin();
			sorig_conc_map[sloc] = i;
			commonsamps++;
		}

	}

	//check for common snps and strandedness
	unsigned int orig_count = 0;
	for(unsigned int i = 0; i < check_data_set->num_loci(); i++){
		orig_count = 0;
		if(orig_count >= markers->size()){
			break;
		}
		Marker* mark = check_data_set->get_locus(i);

		Marker* orig_mark = (*markers)[orig_count];

		//new 12-13-2010  check if marker is enabled before performing analysis...
		if (!mark->isEnabled() || !orig_mark->isEnabled())
		{
			cout << "mark: " << mark->toString() << " is not enabled.  Skipping..." << endl;
			orig_count++;
			continue;
		}

		cout << "check: " << mark->toString() << " orig: " << orig_mark->toString() << endl;

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

		if(mark->getChrom() == orig_mark->getChrom() && mark->getBPLOC() == orig_mark->getBPLOC())
		{
			Marker* omark = orig_mark;
			map<string, int> genomap;
			genomap[omark->getAllele1() + omark->getAllele2()] = 1;
			genomap[omark->getAllele2() + omark->getAllele1()] = 1;
			genomap[omark->getAllele1() + "0"] = 1;
			genomap["0"+omark->getAllele1()] = 1;
			genomap[omark->getAllele2() + "0"] = 1;
			genomap["0" + omark->getAllele2()] = 1;
			genomap["00"] = 1;
			map<string, int>::iterator giter = genomap.find(mark->getAllele1() + mark->getAllele2());
			if(giter == genomap.end())
			{
				bool bad = true;
				//check if the problem is that one file has a monomorphic snp and the other doesn't...
				if(omark->getAllele1() == "0" || omark->getAllele2() == "0")
				{
					if(mark->getAllele1() == omark->getAllele1() || mark->getAllele1() == omark->getAllele2() ||
							mark->getAllele2() == omark->getAllele2() || mark->getAllele2() == omark->getAllele1())
					{
						bad = false;
					}
				}
				if(mark->isMicroSat() || omark->isMicroSat())
				{
					bad = false;
				}
				if(bad)
				{
					opts::printLog("Strandedness issue found for snp: " + mark->toString() + "\n");
					opts::printLog("Original alleles found: " + omark->getAllele1() + "/" + omark->getAllele2() + "\n");
					opts::printLog("Checked alleles found: " + mark->getAllele1() + "/" + mark->getAllele2() + "\n");
					throw MethodException("No input specified for " + stepname + "\n" + "Original alleles found: " + omark->getAllele1() + "/" + omark->getAllele2() + "\n" + "Checked alleles found: " + mark->getAllele1() + "/" + mark->getAllele2() + "\n");
				}
			}
			morig_conc_map[orig_count] = i;
			snp_error_count[orig_count] = 0;
			snp_zero_count[orig_count] = 0;
			snp_total_count[orig_count] = 0;
			commonsnps++;
		}
		else
		{
			//this mark didn't match any marks contained in the original dataset, place information in the mismatches file...
			mismatches << mark->getChrom() << "\t" << mark->getRSID() << "\t" << mark->getProbeID() << "\t" << mark->getBPLOC() << endl;
		}
		orig_count++;
	}

	opts::printLog("Samples found: " + getString<int>(check_data_set->num_inds()) + ", " + getString<int>(commonsamps) + " in common with primary dataset.\n");
	opts::printLog("Markers found: " + getString<int>(check_data_set->num_loci()) + ", " + getString<int>(commonsnps) + " in common with primary dataset.\n");
	opts::printLog("Performing concordance check...");

	int totalchecked = 0;
	int numbad = 0;

	map<int, int>::iterator siter;
	map<int, int>::iterator miter;

	//perform concordance check
	for(int i = 0; i < check_data_set->num_inds(); i++){
		Sample* sconc = check_data_set->get_sample(i);
		for(unsigned int j = 0; j < samples->size(); j++){
		Sample* sorig = (*samples).at(j);
		if(sorig->toString() != sconc->toString()){
			continue;
		}
		int sample_error_count = 0;
		int sample_total_count = 0;
		for(miter = morig_conc_map.begin(); miter != morig_conc_map.end(); miter++){
			int morig = miter->first;
			int msecond = miter->second;
			Marker* mmorig = (*markers).at(morig);
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
				if((!mmorig->isMicroSat() && sorig->getAone(oloc) && sorig->getAtwo(oloc) && sorig->getAmissing(oloc))
						|| (!mconc->isMicroSat() && sconc->getAone(cloc) && sconc->getAtwo(cloc) && sconc->getAmissing(cloc))){
					continue;
				}
			}

			//new 12-13-2010
			if(!options.doIncDisabledSamples())
			{
				if (!mmorig->isEnabled() || !mconc->isEnabled())
				{
					cout << "Mark: " << mmorig->toString() << " is not enabled.  Skipping." << endl;
					continue;
				}
			}

			bool boa1 = sorig->getAone(oloc);
			bool boa2 = sorig->getAtwo(oloc);
			bool boa3 = sorig->getAmissing(oloc);

			bool bca1 = sconc->getAone(cloc);
			bool bca2 = sconc->getAtwo(cloc);
			bool bca3 = sconc->getAmissing(cloc);

			//new 12-14-2010  changed to allow for checking of microsatellites.
			//take care of the original allele
			if(opts::_MICROSATS_ && (mmorig->isMicroSat() || mconc->isMicroSat()))
			{
				if (mmorig->isMicroSat() && (sorig->getAbone(oloc) >=0 && sorig->getAbtwo(oloc) >= 0))
				{
					//there are Abone and Abtwo alleles for the original sample...
					oa1 = mmorig->getAllele(sorig->getAbone(oloc));
					oa2 = mmorig->getAllele(sorig->getAbtwo(oloc));
				}
				else if (mmorig->isMicroSat() && (sorig->getAbone(oloc) < 0 && sorig->getAbtwo(oloc) < 0) && sorig->getAmissing(oloc))
				{
					//this is a microsatellite marker and this sample is "missing", i.e. 0 0
					oa1 = "0";
					oa2 = "0";
				}
				else
				{
					//there are no Abone or Abtwo alleles for the original sample, use the Aone and Atwo alleles
					if(boa1 && boa2 && boa3)
					{
						oa1 = "0";
						oa2 = "0";
					}
					else
					{
						if(!sorig->getAone(oloc))
						{
							oa1 =  mmorig->getAllele1();
						}
						else
						{
							oa1 = mmorig->getAllele2();
						}
						if(!sorig->getAtwo(oloc))
						{
							oa2 = mmorig->getAllele1();
						}
						else
						{
							oa2 = mmorig->getAllele2();
						}
					}
				}

				//take care of the concordant allele...
				if (mconc->isMicroSat() && (sconc->getAbone(cloc) >= 0 && sconc->getAbtwo(cloc) >= 0))
				{
					//there are Abone and Abtwo alleles for the concordant sample, use them
					ca1 = mconc->getAllele(sconc->getAbone(cloc));
					ca2 = mconc->getAllele(sconc->getAbtwo(cloc));
				}
				else if (mconc->isMicroSat() && (sconc->getAbone(cloc) < 0 && sconc->getAbtwo(cloc) < 0) && sconc->getAmissing(cloc))
				{
					//this is a microsatellite marker and this sample is "missing", i.e. 0 0
					ca1 = "0";
					ca2 = "0";
				}
				else
				{
					//there are no Abone and Abtwo alleles for the concordant sample, use the Aone and Atwo alleles
					if(bca1 && bca2 && bca3)
					{
						ca1 = "0";
						ca2 = "0";
					}
					else
					{
						if(!sconc->getAone(cloc))
						{
							ca1 = mconc->getAllele1();
						}
						else
						{
							ca1 = mconc->getAllele2();
						}
						if(!sconc->getAtwo(cloc))
						{
							ca2 = mconc->getAllele1();
						}
						else
						{
							ca2 = mconc->getAllele2();
						}
					}
				}

				if(ca1 == "0" && ca2 == "0" && oa1 == "0" && oa2 == "0")
				{
					both_zero = true;
					snp_zero_count[morig]++;
				}
				else if(options.getIncMissing() && ((ca1 == "0" && ca2 == "0") || (oa1 == "0" && oa2 == "0")))
				{
					inc = true;
				}
				else
				{
					if((oa1 == ca1) && (oa2 == ca2))
					{
						inc = false;
					}
					else if((oa1 == ca2) && (oa2 == ca1))
					{
						inc = false;
					}
					else
					{
						inc = true;
					}
				}
			}
			else
			{
				//not a microsatellite
				if(boa1 && boa2 && boa3){
					oa1 = "0";
					oa2 = "0";
				}
				else
				{
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



				if(bca1 && bca2 && bca3)
				{
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

				if(ca1 == "0" && ca2 == "0" && oa1 == "0" && oa2 == "0")
				{
					both_zero = true;
					snp_zero_count[morig]++;
				}
				else if(options.getIncMissing() && ((ca1 == "0" && ca2 == "0") || (oa1 == "0" && oa2 == "0")))
				{
					inc = true;
				}
				else
				{
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
					//new 12-10-2010 to deal with monomorphic gene
					else if(mmorig->getAllele1() == "0" && (ca1 == mmorig->getAllele2() || ca2 == mmorig->getAllele2()))//o=10
					{
						if(! ((boa1) && boa2)) inc = true;
					}
					else if(mmorig->getAllele2() == "0" && (ca1 == mmorig->getAllele1() || ca2 == mmorig->getAllele1()))//o=01
					{
						if(! ((boa1) && boa2)) inc = true;
					}
					else if(ca1 == "0" && (mmorig->getAllele1() == ca2 || mmorig->getAllele2() == ca2))//c=01
					{
						if(! ((boa1) && boa2)) inc = true;
					}
					else if(ca2 == "0" && (mmorig->getAllele1() == ca1 || mmorig->getAllele2() == ca1))//c=01
					{
						if(! ((boa1) && boa2)) inc = true;
					}
				}
			}
			if(inc && !both_zero)
			{
				numbad++;
				snp_error_count[morig]++;
				sample_error_count++;
				errors << sorig->toString() << "\t";
				if(options.getUniqueId())
				{
					errors << sorig->getDadID() << "\t" << sconc->getDadID() << "\t";
				}
				errors << mmorig->getChrom() << "\t" << mmorig->getProbeID() << "\t" << mmorig->getBPLOC() << "\t" << oa1 << "_" << oa2 << "\t" << ca1 << "_" << ca2 << endl;
				if(options.zeroGenos())
				{
					sorig->addAone(mmorig->getLoc(), true);
					sorig->addAtwo(mmorig->getLoc(), true);
					sorig->addAmissing(mmorig->getLoc(), true);
				}
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

	//output some results
	for(miter = morig_conc_map.begin(); miter != morig_conc_map.end(); miter++){
		int morig = miter->first;
		Marker* mmorig = (*markers).at(morig);
		int total = snp_total_count[morig];
		int err = snp_error_count[morig];
		int zero = snp_zero_count[morig];

		float percent = (1 - ((float)err/(float)total));
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

	float percent = (1 - ((float)numbad / (float)totalchecked));

	opts::printLog("Total genotypes checked: " + getString<int>(totalchecked) + "\n");
	opts::printLog("Total mis-matched genotypes: " + getString<int>(numbad) + "\n");
	opts::printLog("Overall concordance: " + getString<float>(percent) + "\n");

	snp_data.close();
	samp_data.close();
	errors.close();
}
}
