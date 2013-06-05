/**********************************************************************************
 *			Allele Frequency Module
 *
 * Written by: Justin Giles
 *	          Vanderbilt University
 *	          Center for Human Genetics Research
 *
 * Iterates over all genotypes and generates a Major/Minor allele count including
 * frequencies as well as genotype frequencies.
 *
 *
 *File: ProcessAlleleFrequency.cc
 **********************************************************************************/

#include <unistd.h>
#include <sstream>
#include <stdio.h>
#include <iostream>
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
#include <map>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include "ProcessAlleleFrequency.h"
#include "Chrom.h"
#include <General.h>
#include "Helpers.h"
#include <sqlite3.h>
#ifdef PLATOLIB
#include <libsqlitewrapped.h>
#include "Controller.h"
#endif

//TODO:  did not import the Vars.h or Vars.cpp, instead replaced
//			Vars::LOCUS_TABLE with "LOCI" in 10 places



//create a namespace to use when using Plato as a library
#ifdef PLATOLIB
namespace PlatoLib
{
#endif

using namespace Methods;

string ProcessAlleleFrequency::stepname = "allele-freq";

/*
 *Function: FilterSummary
 *Description:
 *Outputs total markers remaining after filtering
 *
 */
void ProcessAlleleFrequency::FilterSummary() {

	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int> (
			opts::_MARKERS_WORKING_ - orig_num_markers) + " (" + getString<
			float> (((float) (opts::_MARKERS_WORKING_ - orig_num_markers)
			/ (float) opts::_MARKERS_WORKING_) * 100.0) + "%) of " + getString<
			int> (opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;
}

/*
 *Function: PrintSummary
 *Description:
 *Resets marker flags in preparation for next step
 */
void ProcessAlleleFrequency::PrintSummary() {
	int msize = data_set->num_loci();
	int ssize = data_set->num_inds();
	for(int s = 0; s < ssize; s++){
		data_set->get_sample(s)->setFlag(false);
	}
	for (int m = 0; m < msize; m++) {
		data_set->get_locus(m)->setFlag(false);
	}
	return;
}

/*
 *Function: filter
 *Description:
 *Not used.
 */
void ProcessAlleleFrequency::filter() {
	return;
}

/*
 *Function: doFilter
 *Description:
 *Filters markers based on overall minor allele frequency
 */
void ProcessAlleleFrequency::doFilter(Marker* mark, AlleleFrequency* af) {
	if (options.doThreshMarkersLow() || options.doThreshMarkersHigh()
			|| options.doRmMono() || options.doRmHetOnly()) {
		int notfoundcount = 0;
		double majfreq = 0;
		double minfreq = 0;

		if (!options.doFilterOverall() && !options.doFilterFile()) {
			//overall is stored in a1/a2_countP
			if (!mark->isMicroSat()) {
				if (af->getAoneP_count() > af->getAtwoP_count()) {
					majfreq = ((double) af->getAoneP_count() / (double) (af->getAoneP_count()
							+ af->getAtwoP_count()));
					minfreq = 1.0f - majfreq;
				} else {
					majfreq = ((double) af->getAtwoP_count() / (double) (af->getAoneP_count()
							+ af->getAtwoP_count()));
					minfreq = 1.0f - majfreq;
				}
				bool inc = false;
				//if(!double_comp(minfreq, options.getThreshMarkersLow()) && minfreq < (double)options.getThreshMarkersLow() && options.doThreshMarkersLow()){
				if ((Helpers::dEquals(minfreq, 0) || Helpers::dEquals(majfreq, 1))
						&& options.doRmMono()) {
					mark->setEnabled(false);
					orig_num_markers++;
					inc = true;
				}

				if (af->getAonehomoP() == 0 && af->getAtwohomoP() == 0
						&& options.doRmHetOnly()) {
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}
				if (Helpers::dLess(minfreq, options.getThreshMarkersLow())
						&& options.doThreshMarkersLow()) {
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}
				if (Helpers::dGreater(minfreq, options.getThreshMarkersHigh())
						&& options.doThreshMarkersHigh()) {
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}

			}
		} else if (options.doFilterOverall() && !options.doFilterFile()) {
			//if filter overall enabled, overall stored in a1/a2_count
			if (!mark->isMicroSat()) {
				if (af->getAone_count() > af->getAtwo_count()) {
					majfreq = ((double) af->getAone_count() / (double) (af->getAone_count()
							+ af->getAtwo_count()));
					minfreq = 1.0f - majfreq;
				} else {
					majfreq = ((double) af->getAtwo_count() / (double) (af->getAone_count()
							+ af->getAtwo_count()));
					minfreq = 1.0f - majfreq;
				}
				bool inc = false;
				//if(minfreq < (double)options.getThreshMarkersLow() && options.doThreshMarkersLow()){
				//if(!double_comp(minfreq, options.getThreshMarkersLow()) && minfreq < options.getThreshMarkersLow() && options.doThreshMarkersLow()){
				if ((Helpers::dEquals(minfreq, 0) || Helpers::dEquals(majfreq, 1))
						&& options.doRmMono()) {
					mark->setEnabled(false);
					orig_num_markers++;
					inc = true;
				}
				else{
					//cout << "Keeping " << mark->toString() << " " << minfreq << " " << majfreq << endl;
				}
				if (af->getAonehomo() == 0 && af->getAtwohomo() == 0
						&& options.doRmHetOnly()) {
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}
				if (Helpers::dLess(minfreq, options.getThreshMarkersLow())
						&& options.doThreshMarkersLow()) {
					//cout << (*markers)[i]->getRSID() << "\t" << minfreq << endl;
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}
				if (Helpers::dGreater(minfreq, options.getThreshMarkersHigh())
						&& options.doThreshMarkersHigh()) {
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}
			}
		} else if (options.doFilterFile() && !options.doFilterOverall()) {
			if (mark->hasMAF()) {
				double majfreq = 0.0;
				double minfreq = 0.0;
				if (af->getAoneP_count() > af->getAtwoP_count()) {
					majfreq = ((double) af->getAoneP_count() / (double) (af->getAoneP_count()
							+ af->getAtwoP_count()));
					minfreq = 1.0f - majfreq;
				} else {
					majfreq = ((double) af->getAtwoP_count() / (double) (af->getAoneP_count()
							+ af->getAtwoP_count()));
					minfreq = 1.0f - majfreq;
				}
				bool inc = false;

				//if(!double_comp(minfreq, options.getThreshMarkersLow()) && minfreq < (double)options.getThreshMarkersLow() && options.doThreshMarkersLow()){
				if ((Helpers::dEquals(minfreq, 0) || Helpers::dEquals(majfreq, 1))
						&& options.doRmMono()) {
					mark->setEnabled(false);
					orig_num_markers++;
					inc = true;
				}
				if (af->getAonehomoP() == 0 && af->getAtwohomoP() == 0
						&& options.doRmHetOnly()) {
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}
				if (Helpers::fLess(mark->getMAF(), options.getThreshMarkersLow())
						&& options.doThreshMarkersLow()) {
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}
				if (Helpers::fGreater(mark->getMAF(), options.getThreshMarkersHigh())
						&& options.doThreshMarkersHigh()) {
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}
			} else {
				notfoundcount++;
			}
		}
		if (options.doFilterFile()) {
			if (notfoundcount > 0) {
				opts::printLog(
						"Markers skipped due to missing predefined MAF: "
								+ getString<int> (notfoundcount) + "\n");
			}
		}
	}
}




/*
 *Function: initializeCounts
 *Description:
 *Sets counts to the specified value (usually 0)
 *
 */
void ProcessAlleleFrequency::initializeCounts(int v) {
	a1_count = a2_count = a1_homo_count = a2_homo_count = a12_count = v;
	a1_countM = a2_countM = a1_homo_countM = a2_homo_countM = a12_countM = v;
	a1_countF = a2_countF = a1_homo_countF = a2_homo_countF = a12_countF = v;
	a1_countP = a2_countP = a1_homo_countP = a2_homo_countP = a12_countP = v;
	a1_countPM = a2_countPM = a1_homo_countPM = a2_homo_countPM = a12_countPM
			= v;
	a1_countPF = a2_countPF = a1_homo_countPF = a2_homo_countPF = a12_countPF
			= v;
	a1_countC = a2_countC = a1_homo_countC = a2_homo_countC = a12_countC = v;
	a1_countCM = a2_countCM = a1_homo_countCM = a2_homo_countCM = a12_countCM
			= v;
	a1_countCF = a2_countCF = a1_homo_countCF = a2_homo_countCF = a12_countCF
			= v;
	a1_countCa = a2_countCa = a1_homo_countCa = a2_homo_countCa = a12_countCa
			= v;
	a1_countCaM = a2_countCaM = a1_homo_countCaM = a2_homo_countCaM
			= a12_countCaM = v;
	a1_countCaF = a2_countCaF = a1_homo_countCaF = a2_homo_countCaF
			= a12_countCaF = v;
	a1_countCon = a2_countCon = a1_homo_countCon = a2_homo_countCon
			= a12_countCon = v;
	a1_countConM = a2_countConM = a1_homo_countConM = a2_homo_countConM
			= a12_countConM = v;
	a1_countConF = a2_countConF = a1_homo_countConF = a2_homo_countConF
			= a12_countConF = v;

	ga1_count.clear();
	ga2_count.clear();
	ga1_homo_count.clear();
	ga2_homo_count.clear();
	ga12_count.clear();

	gm_allele_counts_o.clear();
	gm_geno_counts_o.clear();

	m_allele_counts_o.clear();
	m_allele_counts_om.clear();
	m_allele_counts_of.clear();
	m_geno_counts_o.clear();
	m_geno_counts_om.clear();
	m_geno_counts_of.clear();
	m_allele_counts_p.clear();
	m_allele_counts_pm.clear();
	m_allele_counts_pf.clear();
	m_geno_counts_p.clear();
	m_geno_counts_pm.clear();
	m_geno_counts_pf.clear();
	m_allele_counts_c.clear();
	m_allele_counts_cm.clear();
	m_allele_counts_cf.clear();
	m_geno_counts_c.clear();
	m_geno_counts_cm.clear();
	m_geno_counts_cf.clear();
	m_allele_counts_ca.clear();
	m_allele_counts_cam.clear();
	m_allele_counts_caf.clear();
	m_geno_counts_ca.clear();
	m_geno_counts_cam.clear();
	m_geno_counts_caf.clear();
	m_allele_counts_con.clear();
	m_allele_counts_conm.clear();
	m_allele_counts_conf.clear();
	m_geno_counts_con.clear();
	m_geno_counts_conm.clear();
	m_geno_counts_conf.clear();

}

/*
 *Function: processtest
 *Description:
 *Main function to perform allele frequency test
 *
 */
void ProcessAlleleFrequency::processtest()
{
	int total_snps = 0;
	map<string, double> group_avg;
	string afname;

#ifdef PLATOLIB
	//create a Query object if set to use a database
	Query myQuery(*db);
	//TODO:  create the following method...
	create_tables();
	string insert, pinsert, pgeninsert, gendinsert, ginsert;
	string ccgeninsert, gendgeninsert, ccinsert, ggeninsert;

#else
	//this section sets up the output files, headers and filestreams
	//and is not included if USE_DB is defined
	afname = opts::_OUTPREFIX_ + "allele_freq" + options.getOut()
			+ ".txt";//+ getString<int>(order) + ".txt";
	if (!overwrite)
	{
		afname += "." + getString<int> (order);
	}
	string gfname = opts::_OUTPREFIX_ + "allele_freq_genotype"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite)
	{
		gfname += "." + getString<int> (order);
	}
	string gafname;
	string ggfname;
	string gafnameavg;
	if (options.doGroupFile())
	{
		gafname = opts::_OUTPREFIX_ + "allele_freq_group" + options.getOut()
				+ ".txt";//+ getString<int>(order) + ".txt";
		if (!overwrite)
		{
			gafname += "." + getString<int> (order);
		}
		ggfname = opts::_OUTPREFIX_ + "allele_freq_genotype_group"
				+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if (!overwrite)
		{
			ggfname += "." + getString<int> (order);
		}
		gafnameavg = opts::_OUTPREFIX_ + "allele_freq_group_avg"
				+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if (!overwrite)
		{
			gafnameavg += "." + getString<int> (order);
		}
	}

	string pfname = opts::_OUTPREFIX_ + "allele_freq_parental"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite)
	{
		pfname += "." + getString<int> (order);
	}
	string pgfname = opts::_OUTPREFIX_ + "allele_freq_parental_genotype"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite)
	{
		pgfname += "." + getString<int> (order);
	}
	string gendfname = opts::_OUTPREFIX_ + "allele_freq_gender"
			+ options.getOut() + ".txt"; //getString<int>(order) + ".txt";
	if (!overwrite)
	{
		gendfname += "." + getString<int> (order);
	}
	string gendgfname = opts::_OUTPREFIX_ + "allele_freq_gender_genotype"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite)
	{
		gendgfname += "." + getString<int> (order);
	}
	string ccfname = opts::_OUTPREFIX_ + "allele_freq_casecontrol"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite)
	{
		ccfname += "." + getString<int> (order);
	}
	string ccgfname = opts::_OUTPREFIX_ + "allele_freq_casecontrol_genotype"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite)
	{
		ccgfname += "." + getString<int> (order);
	}
	ofstream paren;
	ofstream pareng;
	ofstream gend;
	ofstream gendg;
	ofstream cc;
	ofstream ccg;

	ofstream myoutput(afname.c_str());
	ofstream mygeno(gfname.c_str());
	opts::addFile("Marker", stepname, afname);
	opts::addFile("Marker", stepname, gfname);
	if (!myoutput)
	{
		opts::printLog("Error opening: " + afname + ".  Exiting!\n");
		throw MethodException("");
	}
	if (!mygeno)
	{
		opts::printLog("Error opening: " + gfname + ".  Exiting!\n");
		throw MethodException("");
	}

	ofstream gmyoutput; //group file output
	ofstream gmygeno; //group file output
	ofstream gmyoutputavg; //avg groups
	if (options.doGroupFile())
	{
		gmyoutput.open(gafname.c_str());
		gmygeno.open(ggfname.c_str());
		opts::addFile("Marker", stepname, gafname);
		opts::addFile("Marker", stepname, ggfname);
		if (!gmyoutput)
		{
			opts::printLog("Error opening: " + gafname + ".  Exiting!\n");
			throw MethodException("");
		}
		if (!gmygeno)
		{
			opts::printLog("Error opening: " + ggfname + ".  Exiting!\n");
			throw MethodException("");
		}
		gmyoutput.precision(4);
		gmygeno.precision(4);

		gmyoutputavg.open(gafnameavg.c_str());
		opts::addFile("Batch", stepname, gafnameavg);
		if(!gmyoutputavg)
		{
			throw MethodException("Error opening: " + gafnameavg + ".\n");
		}
		gmyoutputavg << "Batch\tMAF\tN\n";
		opts::addHeader(gafnameavg, "MAF");
		opts::addHeader(gafnameavg, "N");
	}

#endif

	int msize = data_set->num_loci();

	int maxalleles = 0;
	for (int i = 0; i < msize; i++)
	{
		if (data_set->get_locus(i)->isEnabled())
		{
			if (data_set->get_locus(i)->getNumAlleles() > maxalleles)
			{
				maxalleles = data_set->get_locus(i)->getNumAlleles();
			}
		}
	}

#ifndef PLATOLIB
	//the following code is needed only for the file output
	myoutput.precision(4);
	mygeno.precision(4);
	if (data_set->get_locus(0)->getDetailHeaders().size() > 0)
	{
		myoutput << "Chrom\trsID\tProbeID\tbploc\t"
				<< data_set->get_locus(0)->getDetailHeaders();
		mygeno << "Chrom\trsID\tProbeID\tbploc\t"
				<< data_set->get_locus(0)->getDetailHeaders()
				<< "\tGenotype11\tGenotype12\tGenotype22\tOverall_Freq_Genotype11\tOverall_Freq_Genotype12\tOverall_Freq_Genotype22\tOverall_Count_Genotype11\tOverall_Count_Genotype12\tOverall_Count_Genotype22\tCase_Freq_Genotype11\tCase_Freq_Genotype12\tCase_Freq_Genotype22\tCase_Count_Genotype11\tCase_Count_Genotype12\tCase_Count_Genotype22\tControl_Freq_Genotype11\tControl_Freq_Genotype12\tControl_Freq_Genotype22\tControl_Count_Genotype11\tControl_Count_Genotype12\tControl_Count_Genotype22";
		if (options.doGroupFile())
		{
			gmyoutput << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders();
			gmygeno << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders()
					<< "\tGenotype11\tGenotype12\tGenotype22";
		}
	}
	else
	{
		myoutput << "Chrom\trsID\tProbeID\tbploc";
		mygeno
				<< "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22\tOverall_Freq_Genotype11\tOverall_Freq_Genotype12\tOverall_Freq_Genotype22\tOverall_Count_Genotype11\tOverall_Count_Genotype12\tOverall_Count_Genotype22\tCase_Freq_Genotype11\tCase_Freq_Genotype12\tCase_Freq_Genotype22\tCase_Count_Genotype11\tCase_Count_Genotype12\tCase_Count_Genotype22\tControl_Freq_Genotype11\tControl_Freq_Genotype12\tControl_Freq_Genotype22\tControl_Count_Genotype11\tControl_Count_Genotype12\tControl_Count_Genotype22";
		if (options.doGroupFile())
		{
			gmyoutput << "Chrom\trsID\tProbeID\tbploc";
			gmygeno
					<< "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22";
		}
	}
	if (options.doGroupFile())
	{
		opts::addHeader(ggfname, "Genotype11");
		opts::addHeader(ggfname, "Genotype12");
		opts::addHeader(ggfname, "Genotype22");
	}

	opts::addHeader(gfname, "Genotype11");
	opts::addHeader(gfname, "Genotype12");
	opts::addHeader(gfname, "Genotype22");
	opts::addHeader(gfname, "Overall_Freq_Genotype11");
	opts::addHeader(gfname, "Overall_Freq_Genotype12");
	opts::addHeader(gfname, "Overall_Freq_Genotype22");
	opts::addHeader(gfname, "Overall_Count_Genotype11");
	opts::addHeader(gfname, "Overall_Count_Genotype12");
	opts::addHeader(gfname, "Overall_Count_Genotype22");
	opts::addHeader(gfname, "Case_Freq_Genotype11");
	opts::addHeader(gfname, "Case_Freq_Genotype12");
	opts::addHeader(gfname, "Case_Freq_Genotype22");
	opts::addHeader(gfname, "Case_Count_Genotype11");
	opts::addHeader(gfname, "Case_Count_Genotype12");
	opts::addHeader(gfname, "Case_Count_Genotype22");
	opts::addHeader(gfname, "Control_Freq_Genotype11");
	opts::addHeader(gfname, "Control_Freq_Genotype12");
	opts::addHeader(gfname, "Control_Freq_Genotype22");
	opts::addHeader(gfname, "Control_Count_Genotype11");
	opts::addHeader(gfname, "Control_Count_Genotype12");
	opts::addHeader(gfname, "Control_Count_Genotype22");

	for (int i = 0; i < maxalleles; i++)
	{
		myoutput << "\tAllele" << (i + 1);
		opts::addHeader(afname, "Allele" + getString<int> (i + 1));
		if (options.doGroupFile())
		{
			gmyoutput << "\tAllele" << (i + 1);
			opts::addHeader(gafname, "Allele" + getString<int> (i + 1));
		}
	}
	for (int i = 0; i < maxalleles; i++)
	{
		myoutput << "\t" << "Overall_Allele" << (i + 1) << "_freq";
		opts::addHeader(afname, "Overall_Allele" + getString<int> (i + 1)
				+ "_freq");
	}
	for (int i = 0; i < maxalleles; i++)
	{
		myoutput << "\t" << "Overall_Allele" << (i + 1) << "_count";
		opts::addHeader(afname, "Overall_Allele" + getString<int> (i + 1)
				+ "_count");
	}
	for (int i = 0; i < maxalleles; i++)
	{
		myoutput << "\t" << "Case_Allele" << (i + 1) << "_freq";
		opts::addHeader(afname, "Case_Allele" + getString<int> (i + 1)
				+ "_freq");
	}
	for (int i = 0; i < maxalleles; i++)
	{
		myoutput << "\t" << "Case_Allele" << (i + 1) << "_count";
		opts::addHeader(afname, "Case_Allele" + getString<int> (i + 1)
				+ "_count");
	}
	for (int i = 0; i < maxalleles; i++)
	{
		myoutput << "\t" << "Control_Allele" << (i + 1) << "_freq";
		opts::addHeader(afname, "Control_Allele" + getString<int> (i + 1)
				+ "_freq");
	}
	for (int i = 0; i < maxalleles; i++)
	{
		myoutput << "\t" << "Control_Allele" << (i + 1) << "_count";
		opts::addHeader(afname, "Control_Allele" + getString<int> (i + 1)
				+ "_count");
	}

	if (options.doGroupFile())
	{
		map<string, vector<Sample*> >::iterator giter;
		map<string, vector<Sample*> > groups = options.getGroups();
		for (giter = groups.begin(); giter != groups.end(); giter++)
		{
			string mygroup = giter->first;
			for (int i = 0; i < maxalleles; i++)
			{
				gmyoutput << "\t" << mygroup << "_Allele" << (i + 1) << "_freq";
				opts::addHeader(gafname, mygroup + "_Allele" + getString<int> (
						i + 1) + "_freq");
			}
			for (int i = 0; i < maxalleles; i++)
			{
				gmyoutput << "\t" << mygroup << "_Allele" << (i + 1)
						<< "_count";
				opts::addHeader(gafname, mygroup + "_Allele" + getString<int> (
						i + 1) + "_count");
			}
			gmygeno << "\t" << mygroup << "_Freq_Genotype11\t" << mygroup
					<< "_Freq_Genotype12\t" << mygroup << "_Freq_Genotype22\t"
					<< mygroup << "_Count_Genotype11\t" << mygroup
					<< "_Count_Genotype12\t" << mygroup << "_Count_Genotype22";
			opts::addHeader(ggfname, mygroup + "_Freq_Genotype11");
			opts::addHeader(ggfname, mygroup + "_Freq_Genotype12");
			opts::addHeader(ggfname, mygroup + "_Freq_Genotype22");
			opts::addHeader(ggfname, mygroup + "_Count_Genotype11");
			opts::addHeader(ggfname, mygroup + "_Count_Genotype12");
			opts::addHeader(ggfname, mygroup + "_Count_Genotype22");
		}
		gmygeno << endl;
		gmyoutput << endl;
	}
	mygeno << endl;
	myoutput << endl;

	if (options.doParental())
	{
		paren.open(pfname.c_str(), ios::out);
		pareng.open(pgfname.c_str(), ios::out);
		if (!paren)
		{
			opts::printLog("Error opening: " + pfname + ". Exiting!\n");
			throw MethodException("");
		}
		if (!pareng)
		{
			opts::printLog("Error opening: " + pgfname + ". Exiting!\n");
			throw MethodException("");
		}
		opts::addFile("Marker", stepname, pfname);
		opts::addFile("Marker", stepname, pgfname);
		paren.precision(4);
		pareng.precision(4);
		if (data_set->get_locus(0)->getDetailHeaders().size() > 0)
		{
			paren << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders();
			pareng << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders()
					<< "\tGenotype11\tGenotype12\tGenotype22\tParent_Male_Freq_Genotype11\tParent_Male_Freq_Genotype12\tParent_Male_Freq_Genotype22\tParent_Male_Count_Genotype11\tParent_Male_Count_Genotype12\tParent_Male_Count_Genotype22\tParent_Female_Freq_Genotype11\tParent_Female_Freq_Genotype12\tParent_Female_Freq_Genotype22\tParent_Female_Count_Genotype11\tParent_Female_Count_Genotype12\tParent_Female_Count_Genotype22\n";
		}
		else
		{
			paren << "Chrom\trsID\tProbeID\tbploc";
			pareng << "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22\tParent_Male_Freq_Genotype11\tParent_Male_Freq_Genotype12\tParent_Male_Freq_Genotype22\tParent_Male_Count_Genotype11\tParent_Male_Count_Genotype12\tParent_Male_Count_Genotype22\tParent_Female_Freq_Genotype11\tParent_Female_Freq_Genotype12\tParent_Female_Freq_Genotype22\tParent_Female_Count_Genotype11\tParent_Female_Count_Genotype12\tParent_Female_Count_Genotype22\n";
		}
		opts::addHeader(pgfname, "Parent_Male_Freq_Genotype11");
		opts::addHeader(pgfname, "Parent_Male_Freq_Genotype12");
		opts::addHeader(pgfname, "Parent_Male_Freq_Genotype22");
		opts::addHeader(pgfname, "Parent_Male_Count_Genotype11");
		opts::addHeader(pgfname, "Parent_Male_Count_Genotype12");
		opts::addHeader(pgfname, "Parent_Male_Count_Genotype22");
		opts::addHeader(pgfname, "Parent_Female_Freq_Genotype11");
		opts::addHeader(pgfname, "Parent_Female_Freq_Genotype12");
		opts::addHeader(pgfname, "Parent_Female_Freq_Genotype22");
		opts::addHeader(pgfname, "Parent_Female_Count_Genotype11");
		opts::addHeader(pgfname, "Parent_Female_Count_Genotype12");
		opts::addHeader(pgfname, "Parent_Female_Count_Genotype22");

		for (int i = 0; i < maxalleles; i++)
		{
			paren << "\t" << "Allele" << (i + 1);
		}
		for (int i = 0; i < maxalleles; i++)
		{
			paren << "\t" << "Parent_Male_Allele" << (i + 1) << "_freq";
			opts::addHeader(pfname, "Parent_Male_Allele" + getString<int> (i
					+ 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			paren << "\t" << "Parent_Male_Allele" << (i + 1) << "_count";
			opts::addHeader(pfname, "Parent_Male_Allele" + getString<int> (i
					+ 1) + "_count");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			paren << "\t" << "Parent_Female_Allele" << (i + 1) << "_freq";
			opts::addHeader(pfname, "Parent_Female_Allele" + getString<int> (i
					+ 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			paren << "\t" << "Parent_Female_Allele" << (i + 1) << "_count";
			opts::addHeader(pfname, "Parent_Female_Allele" + getString<int> (i
					+ 1) + "_count");
		}
		paren << endl;
	}
	if (options.doGender())
	{
		gend.open(gendfname.c_str(), ios::out);
		gendg.open(gendgfname.c_str(), ios::out);
		opts::addFile("Marker", stepname, gendfname);
		opts::addFile("Marker", stepname, gendgfname);
		if (!gend)
		{
			opts::printLog("Error opening: " + gendfname + ". Exiting!\n");
			throw MethodException("");
		}
		if (!gendg)
		{
			opts::printLog("Error opening: " + gendgfname + ". Exiting!\n");
			throw MethodException("");
		}
		gend.precision(4);
		gendg.precision(4);
		if (data_set->get_locus(0)->getDetailHeaders().size() > 0)
		{
			gend << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders();
			gendg << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders()
					<< "\tGenotype11\tGenotype12\tGenotype22\tOverall_Male_Freq_Genotype11\tOverall_Male_Freq_Genotype12\tOverall_Male_Freq_Genotype22\tOverall_Male_Count_Genotype11\tOverall_Male_Count_Genotype12\tOverall_Male_Count_Genotype22\tOverall_Female_Freq_Genotype11\tOverall_Female_Freq_Genotype12\tOverall_Female_Freq_Genotype22\tOverall_Female_Count_Genotype11\tOverall_Female_Count_Genotype12\tOverall_Female_Count_Genotype22\n";
		}
		else
		{
			gend << "Chrom\trsID\tProbeID\tbploc";
			gendg << "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22\tOverall_Male_Freq_Genotype11\tOverall_Male_Freq_Genotype12\tOverall_Male_Freq_Genotype22\tOverall_Male_Count_Genotype11\tOverall_Male_Count_Genotype12\tOverall_Male_Count_Genotype22\tOverall_Female_Freq_Genotype11\tOverall_Female_Freq_Genotype12\tOverall_Female_Freq_Genotype22\tOverall_Female_Count_Genotype11\tOverall_Female_Count_Genotype12\tOverall_Female_Count_Genotype22\n";
		}
		opts::addHeader(gendgfname, "Overall_Male_Freq_Genotype11");
		opts::addHeader(gendgfname, "Overall_Male_Freq_Genotype12");
		opts::addHeader(gendgfname, "Overall_Male_Freq_Genotype22");
		opts::addHeader(gendgfname, "Overall_Male_Count_Genotype11");
		opts::addHeader(gendgfname, "Overall_Male_Count_Genotype12");
		opts::addHeader(gendgfname, "Overall_Male_Count_Genotype22");
		opts::addHeader(gendgfname, "Overall_Female_Freq_Genotype11");
		opts::addHeader(gendgfname, "Overall_Female_Freq_Genotype12");
		opts::addHeader(gendgfname, "Overall_Female_Freq_Genotype22");
		opts::addHeader(gendgfname, "Overall_Female_Count_Genotype11");
		opts::addHeader(gendgfname, "Overall_Female_Count_Genotype12");
		opts::addHeader(gendgfname, "Overall_Female_Count_Genotype22");

		for (int i = 0; i < maxalleles; i++)
		{
			gend << "\t" << "Allele" << (i + 1);
		}
		for (int i = 0; i < maxalleles; i++)
		{
			gend << "\t" << "Overall_Male_Allele" << (i + 1) << "_freq";
			opts::addHeader(gendfname, "Overall_Male_Allele" + getString<int> (
					i + 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			gend << "\t" << "Overall_Male_Allele" << (i + 1) << "_count";
			opts::addHeader(gendfname, "Overall_Male_Allele" + getString<int> (
					i + 1) + "_count");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			gend << "\t" << "Overall_Female_Allele" << (i + 1) << "_freq";
			opts::addHeader(gendfname, "Overall_Female_Allele"
					+ getString<int> (i + 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			gend << "\t" << "Overall_Female_Allele" << (i + 1) << "_count";
			opts::addHeader(gendfname, "Overall_Female_Allele"
					+ getString<int> (i + 1) + "_count");
		}
		gend << endl;
	}
	if (options.doCaseControl())
	{
		cc.open(ccfname.c_str(), ios::out);
		ccg.open(ccgfname.c_str(), ios::out);
		opts::addFile("Marker", stepname, ccfname);
		opts::addFile("Marker", stepname, ccgfname);
		if (!cc)
		{
			opts::printLog("Error opening: " + ccfname + ". Exiting!\n");
			throw MethodException("");
		}
		if (!ccg)
		{
			opts::printLog("Error opening: " + ccgfname + ". Exiting!\n");
			throw MethodException("");
		}
		cc.precision(4);
		ccg.precision(4);

		if (data_set->get_locus(0)->getDetailHeaders().size() > 0)
		{
			cc << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders();
			ccg << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders()
					<< "\tGenotype11\tGenotype12\tGenotype22\tCase_Male_Freq_Genotype11\tCase_Male_Freq_Genotype12\tCase_Male_Freq_Genotype22\tCase_Male_Count_Genotype11\tCase_Male_Count_Genotype12\tCase_Male_Count_Genotype22\tCase_Female_Freq_Genotype11\tCase_Female_Freq_Genotype12\tCase_Female_Freq_Genotype22\tCase_Female_Count_Genotype11\tCase_Female_Count_Genotype12\tCase_Female_Count_Genotype22\tControl_Male_Freq_Genotype11\tControl_Male_Freq_Genotype12\tControl_Male_Freq_Genotype22\tControl_Male_Count_Genotype11\tControl_Male_Count_Genotype12\tControl_Male_Count_Genotype22\tControl_Female_Freq_Genotype11\tControl_Female_Freq_Genotype12\tControl_Female_Freq_Genotype22\tControl_Female_Count_Genotype11\tControl_Female_Count_Genotype12\tControl_Female_Count_Genotype22\n";
		}
		else
		{
			cc << "Chrom\trsID\tProbeID\tbploc";
			ccg	<< "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22\tCase_Male_Freq_Genotype11\tCase_Male_Freq_Genotype12\tCase_Male_Freq_Genotype22\tCase_Male_Count_Genotype11\tCase_Male_Count_Genotype12\tCase_Male_Count_Genotype22\tCase_Female_Freq_Genotype11\tCase_Female_Freq_Genotype12\tCase_Female_Freq_Genotype22\tCase_Female_Count_Genotype11\tCase_Female_Count_Genotype12\tCase_Female_Count_Genotype22\tControl_Male_Freq_Genotype11\tControl_Male_Freq_Genotype12\tControl_Male_Freq_Genotype22\tControl_Male_Count_Genotype11\tControl_Male_Count_Genotype12\tControl_Male_Count_Genotype22\tControl_Female_Freq_Genotype11\tControl_Female_Freq_Genotype12\tControl_Female_Freq_Genotype22\tControl_Female_Count_Genotype11\tControl_Female_Count_Genotype12\tControl_Female_Count_Genotype22\n";
		}
		opts::addHeader(ccgfname, "Case_Male_Freq_Genotype11");
		opts::addHeader(ccgfname, "Case_Male_Freq_Genotype12");
		opts::addHeader(ccgfname, "Case_Male_Freq_Genotype22");
		opts::addHeader(ccgfname, "Case_Male_Count_Genotype11");
		opts::addHeader(ccgfname, "Case_Male_Count_Genotype12");
		opts::addHeader(ccgfname, "Case_Male_Count_Genotype22");
		opts::addHeader(ccgfname, "Case_Female_Freq_Genotype11");
		opts::addHeader(ccgfname, "Case_Female_Freq_Genotype12");
		opts::addHeader(ccgfname, "Case_Female_Freq_Genotype22");
		opts::addHeader(ccgfname, "Case_Female_Count_Genotype11");
		opts::addHeader(ccgfname, "Case_Female_Count_Genotype12");
		opts::addHeader(ccgfname, "Case_Female_Count_Genotype22");
		opts::addHeader(ccgfname, "Control_Male_Freq_Genotype11");
		opts::addHeader(ccgfname, "Control_Male_Freq_Genotype12");
		opts::addHeader(ccgfname, "Control_Male_Freq_Genotype22");
		opts::addHeader(ccgfname, "Control_Male_Count_Genotype11");
		opts::addHeader(ccgfname, "Control_Male_Count_Genotype12");
		opts::addHeader(ccgfname, "Control_Male_Count_Genotype22");
		opts::addHeader(ccgfname, "Control_Female_Freq_Genotype11");
		opts::addHeader(ccgfname, "Control_Female_Freq_Genotype12");
		opts::addHeader(ccgfname, "Control_Female_Freq_Genotype22");
		opts::addHeader(ccgfname, "Control_Female_Count_Genotype11");
		opts::addHeader(ccgfname, "Control_Female_Count_Genotype12");
		opts::addHeader(ccgfname, "Control_Female_Count_Genotype22");

		for (int i = 0; i < maxalleles; i++)
		{
			cc << "\t" << "Allele" << (i + 1);
		}
		for (int i = 0; i < maxalleles; i++)
		{
			cc << "\t" << "Case_Male_Allele" << (i + 1) << "_freq";
			opts::addHeader(ccfname, "Case_Male_Allele"
					+ getString<int> (i + 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			cc << "\t" << "Case_Male_Allele" << (i + 1) << "_count";
			opts::addHeader(ccfname, "Case_Male_Allele"
					+ getString<int> (i + 1) + "_count");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			cc << "\t" << "Case_Female_Allele" << (i + 1) << "_freq";
			opts::addHeader(ccfname, "Case_Female_Allele" + getString<int> (i
					+ 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			cc << "\t" << "Case_Female_Allele" << (i + 1) << "_count";
			opts::addHeader(ccfname, "Case_Female_Allele" + getString<int> (i
					+ 1) + "_count");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			cc << "\t" << "Control_Male_Allele" << (i + 1) << "_freq";
			opts::addHeader(ccfname, "Control_Male_Allele" + getString<int> (i
					+ 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			cc << "\t" << "Control_Male_Allele" << (i + 1) << "_count";
			opts::addHeader(ccfname, "Control_Male_Allele" + getString<int> (i
					+ 1) + "_count");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			cc << "\t" << "Control_Female_Allele" << (i + 1) << "_freq";
			opts::addHeader(ccfname, "Control_Female_Allele" + getString<int> (
					i + 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++)
		{
			cc << "\t" << "Control_Female_Allele" << (i + 1) << "_count";
			opts::addHeader(ccfname, "Control_Female_Allele" + getString<int> (
					i + 1) + "_count");
		}
		cc << endl;
	}
#endif
	//begin processing
#ifdef PLATOLIB
	//need these strings to hold DB insert statements if using DB
	insert = "";
	ginsert = "";
#endif
	AlleleFrequency af;
	af.resetDataSet(data_set);
	af.initializeCounts(0);

	if (options.doRandomChild() || options.doAll() || options.doAllChildren()
			|| options.doUnaffSpousesOnly() || options.doUnknownSpouses()) {
		af.setOptions(options);
	} else {
		options.setFoundersOnly();
		af.setOptions(options);
	}

	if(options.doAll() || options.doFilterOverall()){
		useoverall = true;
	}

#ifdef PLATOLIB
	//open a transaction for the following insert statements if using DB
	myQuery.transaction();
#endif
	vector<Marker*> good_markers = Helpers::findValidMarkers(data_set->get_markers(), &options);
	msize = good_markers.size();

	for (int k = 0; k < msize; k++) {
		Marker* mark = good_markers[k];
		if (mark->isEnabled()){// && isValidMarker(data_set->get_locus(k),&options, prev_base, prev_chrom)) {

			total_snps++;

			//perform calculations
			af.calcOne(mark);//data_set->get_locus(k));
			if (options.doGroupFile()) {
				af.calcOneGroups(mark);//data_set->get_locus(k));
			}
			doFilter(mark, &af);//data_set->get_locus(k), &af);

			if (mark->isMicroSat()){//data_set->get_locus(k)->isMicroSat()) {
#ifndef PLATOLIB
				myoutput << mark->toString();//data_set->get_locus(k)->toString();
#else
				string insert = defaultinsert;
				//insert += "," + mark->toString();
				insert += "," + getString<int>(mark->getBPLOC());
				string geninsert = defaultgenoinsert;
				//geninsert += "," + mark->toString();
				geninsert += "," + getString<int>(mark->getBPLOC());
#endif
				if (options.doGroupFile()) {
#ifndef PLATOLIB
					gmyoutput << mark->toString();//data_set->get_locus(k)->toString();
#else
					string ginsert = groupinsert;
					//ginsert += "," + mark->toString();
					ginsert += "," + getString<int>(mark->getBPLOC());
					string ggeninsert = groupgenoinsert;
					//ggeninsert += "," + mark->toString();
					ggeninsert += "," + getString<int>(mark->getBPLOC());
#endif
				}
				if (options.doParental()) {
#ifndef PLATOLIB
					paren << mark->toString();//data_set->get_locus(k)->toString();
#else
					string pinsert = parentalinsert;
					//pinsert += "," + mark->toString();
					pinsert += "," + getString<int>(mark->getBPLOC());
					string pgeninsert = parentalgenoinsert;
					//pgeninsert += "," + mark->toString();
					pgeninsert += "," + getString<int>(mark->getBPLOC());
#endif
				}
				if (options.doGender()) {
#ifndef PLATOLIB
					gend << mark->toString();//data_set->get_locus(k)->toString();
#else
					string gendinsert = genderinsert;
					//gendinsert += "," + mark->toString();
					gendinsert += "," + getString<int>(mark->getBPLOC());
					string gendgeninsert = gendergenoinsert;
					//gendgeninsert += "," + mark->toString();
					gendgeninsert += "," + getString<int>(mark->getBPLOC());
#endif
				}
				if (options.doCaseControl()) {
#ifndef PLATOLIB
					cc << mark->toString();//data_set->get_locus(k)->toString();
#else
					ccinsert = casecontrolinsert;
					//ccinsert += "," + mark->toString();
					ccinsert += "," + getString<int>(mark->getBPLOC());
					ccgeninsert = casecontrolgenoinsert;
					//ccgeninsert += "," + mark->toString();
					ccgeninsert += "," + getString<int>(mark->getBPLOC());
#endif
				}
				int total_o = 0;
				int total_ca = 0;
				int total_con = 0;
				int numalleles = mark->getNumAlleles();//data_set->get_locus(k)->getNumAlleles();
				for (int a = 0; a < numalleles; a++) {
					if (useoverall) {
						total_o += af.getMicroCount(a);
					} else {
						total_o += af.getMicroCountP(a);
					}
					total_ca += af.getMicroCountCa(a);
					total_con += af.getMicroCountCon(a);
				}
				for (int a = 0; a < numalleles; a++) {
#ifndef PLATOLIB
					myoutput << "\t" << mark->getAllele(a);//data_set->get_locus(k)->getAllele(a);
#else
					insert += ",'" + mark->getAllele(a) + "'";
#endif
					if (options.doGroupFile()) {
#ifndef PLATOLIB
						gmyoutput << "\t" << mark->getAllele(a);//data_set->get_locus(k)->getAllele(a);
#else
						ginsert += ",'" + mark->getAllele(a) + "'";
#endif
					}
					if (options.doParental()) {
#ifndef PLATOLIB
						paren << "\t" << mark->getAllele(a);//data_set->get_locus(k)->getAllele(a);
#else
						pinsert += ",'" + mark->getAllele(a) + "'";
#endif
					}
					if (options.doGender()) {
#ifndef PLATOLIB
						gend << "\t" << mark->getAllele(a);//data_set->get_locus(k)->getAllele(a);
#else
						gendinsert += ",'" + mark->getAllele(a) + "'";
#endif
					}
					if (options.doCaseControl()) {
#ifndef PLATOLIB
						cc << "\t" << mark->getAllele(a);//data_set->get_locus(k)->getAllele(a);
#else
						ccinsert += ",'" + mark->getAllele(a) = "'";
#endif
					}
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
						myoutput << "\tNA";
#else
						insert += ",NULL";
#endif
						if (options.doGroupFile()) {
#ifndef PLATOLIB
							gmyoutput << "\tNA";
#else
							ginsert += ",NULL";
#endif
						}
						if (options.doParental()) {
#ifndef PLATOLIB
							paren << "\tNA";
#else
							pinsert += ",NULL";
#endif
						}
						if (options.doGender()) {
#ifndef PLATOLIB
							gend << "\tNA";
#else
							geninsert += ",NULL";
#endif
						}
						if (options.doCaseControl()) {
#ifndef PLATOLIB
							cc << "\tNA";
#else
							ccinsert += ",NULL";
#endif
						}
					}
				}
				//overall
				for (int a = 0; a < numalleles; a++) {
					float freq = 0.0f;
					if (useoverall) {
						freq = ((float) af.getMicroCount(a) / (float) total_o);
					} else {
						freq = ((float) af.getMicroCountP(a) / (float) total_o);
					}
#ifndef PLATOLIB
					myoutput << "\t" << freq;
#else
					insert += "," + getString<float>(freq);
#endif
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
						myoutput << "\tNA";
#else
						insert += ",NULL";
#endif
					}
				}
				for (int a = 0; a < numalleles; a++) {
					if (useoverall) {
#ifndef PLATOLIB
						myoutput << "\t" << af.getMicroCount(a);
#else
						insert += "," + getString<int>(af.getMicroCount(a));
#endif
					} else {
#ifndef PLATOLIB
						myoutput << "\t" << af.getMicroCountP(a);
#else
						insert += "," + getString<int>(af.getMicroCountP(a));
#endif
					}
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
						myoutput << "\tNA";
#else
						insert += ",NULL";
#endif
					}
				}
				//case
				for (int a = 0; a < numalleles; a++) {
					float freq = ((float) af.getMicroCountCa(a)
							/ (float) total_ca);
#ifndef PLATOLIB
					myoutput << "\t" << freq;
#else
					insert += "," + getString<float>(freq);
#endif
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
						myoutput << "\tNA";
#else
						insert += ",NULL";
#endif
					}
				}
				for (int a = 0; a < numalleles; a++) {
#ifndef PLATOLIB
					myoutput << "\t" << af.getMicroCountCa(a);
#else
					insert += "," + getString<int>(af.getMicroCountCa(a));
#endif
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
						myoutput << "\tNA";
#else
						insert += ",NULL";
#endif
					}
				}
				//control
				for (int a = 0; a < numalleles; a++) {
					float freq = ((float) af.getMicroCountCon(a)
							/ (float) total_con);
#ifndef PLATOLIB
					myoutput << "\t" << freq;
#else
					insert += "," + getString<float>(freq);
#endif
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
						myoutput << "\tNA";
#else
						insert += ",NULL";
#endif
					}
				}
				for (int a = 0; a < numalleles; a++) {
#ifndef PLATOLIB
					myoutput << "\t" << af.getMicroCountCon(a);
#else
					insert += "," + getString<int>(af.getMicroCountCon(a));
#endif
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
						myoutput << "\tNA";
#else
						insert += ",NULL";
#endif
					}
				}
				//groups?
#ifndef PLATOLIB
				mygeno << mark->toString();//data_set->get_locus(k)->toString();
#endif
				if (options.doGroupFile()) {
#ifndef PLATOLIB
					gmygeno << mark->toString();//data_set->get_locus(k)->toString();
#endif
					int gm_total = 0;
					map<string, vector<Sample*> > groups = options.getGroups();
					map<string, vector<Sample*> >::iterator giter;
					for (giter = groups.begin(); giter != groups.end(); giter++) {
						string mygroup = giter->first;
						for (int a = 0; a < numalleles; a++) {
							gm_total += af.getGroupMicroCount(mygroup, a);
						}
						for (int a = 0; a < numalleles; a++) {
							float freq =
									((float) af.getGroupMicroCount(mygroup, a)
											/ (float) gm_total);
#ifndef PLATOLIB
							gmyoutput << "\t" << freq;
#else
							ginsert += "," + getString<float>(freq);
#endif
						}
						if (maxalleles > numalleles) {
							for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
								gmyoutput << "\tNA";
#else
								ginsert += ",NULL";
#endif
							}
						}
						for (int a = 0; a < numalleles; a++) {
#ifndef PLATOLIB
							gmyoutput << "\t" << af.getGroupMicroCount(mygroup, a);
#else
							ginsert += "," + getString<int>(af.getGroupMicroCount(mygroup, a));
#endif
						}
						if (maxalleles > numalleles) {
							for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
								gmyoutput << "\tNA";
#else
								ginsert += ",NULL";
#endif
							}
						}
#ifndef PLATOLIB
						gmygeno << "\tNA\tNA\tNA\tNA\tNA\tNA";
#else
						ginsert += ",NULL,NULL,NULL,NULL,NULL,NULL";
#endif
					}
#ifndef PLATOLIB
					gmyoutput << endl;
#else
					ginsert += ")";
					//TODO: add Controller to Plato library
					Controller::execute_sql(myQuery, ginsert);
#endif
				}
#ifndef PLATOLIB
				myoutput << endl;
#else
				insert += ")";
				//TODO: add Controller to Plato library
				Controller::execute_sql(myQuery, insert);
#endif

				for (int l = 0; l < 21; l++) {
#ifndef PLATOLIB
					mygeno << "\tNA";
#else
					geninsert += ",NULL";
#endif
					if (options.doGroupFile()) {
#ifndef PLATOLIB
						gmygeno << "\tNA";
#else
						ggeninsert += ",NULL";
#endif
					}
				}
#ifndef PLATOLIB
				mygeno << endl;
#else
				geninsert += ")";
				//TODO:  Controller
				Controller::execute_sql(myQuery, geninsert);
#endif
				if (options.doGroupFile()) {
#ifndef PLATOLIB
					gmygeno << endl;
#else
					ggeninsert += ")";
					//TODO: controller
					Controller::execute_sql(myQuery, ggeninsert);
#endif
				}

				if (options.doParental()) {
#ifndef PLATOLIB
					pareng << mark->toString();//data_set->get_locus(k)->toString();
#endif
					int total_pm = 0;
					int total_pf = 0;
					for (int a = 0; a < numalleles; a++) {
						total_pm += af.getMicroCountPM(a);
						total_pf += af.getMicroCountPF(a);
					}
					//Male
					for (int a = 0; a < numalleles; a++) {
						float freq = ((float) af.getMicroCountPM(a)
								/ (float) total_pm);
#ifndef PLATOLIB
						paren << "\t" << freq;
#else
						pinsert += "," + getString<float>(freq);
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							paren << "\tNA";
#else
							pinsert += ",NULL";
#endif
						}
					}
					for (int a = 0; a < numalleles; a++) {
#ifndef PLATOLIB
						paren << "\t" << af.getMicroCountPM(a);
#else
						pinsert += "," + af.getMicroCountPM(a);
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							paren << "\tNA";
#else
							pinsert += ",NULL";
#endif
						}
					}
					//Female
					for (int a = 0; a < numalleles; a++) {
						float freq = ((float) af.getMicroCountPF(a)
								/ (float) total_pf);
#ifndef PLATOLIB
						paren << "\t" << freq;
#else
						pinsert += "," + getString<float>(freq);
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							paren << "\tNA";
#else
							pinsert += ",NULL";
#endif
						}
					}
					for (int a = 0; a < numalleles; a++) {
#ifndef PLATOLIB
						paren << "\t" << af.getMicroCountPF(a);
#else
						pinsert += "," + getString<int>(af.getMicroCountPF(a));
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							paren << "\tNA";
#else
							pinsert += ",NULL";
#endif
						}
					}
#ifndef PLATOLIB
					paren << endl;
#else
					pinsert += ")";
					//TODO: controller
					Controller::execute_sql(myQuery, pinsert);
#endif

					for (int l = 0; l < 15; l++) {
#ifndef PLATOLIB
						pareng << "\tNA";
#else
						pgeninsert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					pareng << endl;
#else
					pgeninsert += ")";
					//TODO: Controller
					Controller::execute_sql(myQuery, pgeninsert);
#endif

				}
				if (options.doGender()) {
#ifndef PLATOLIB
					gendg << mark->toString();//data_set->get_locus(k)->toString();
#endif
					int total_pm = 0;
					int total_pf = 0;
					for (int a = 0; a < numalleles; a++) {
						if (useoverall) {
							total_pm += af.getMicroCountM(a);
							total_pf += af.getMicroCountF(a);
						} else {
							total_pm += af.getMicroCountPM(a);
							total_pf += af.getMicroCountPF(a);
						}
					}
					//Male
					for (int a = 0; a < numalleles; a++) {
						float freq = 0.0f;
						if (useoverall) {
							freq = ((float) af.getMicroCountM(a)
									/ (float) total_pm);
						} else {
							freq = ((float) af.getMicroCountPM(a)
									/ (float) total_pm);
						}
#ifndef PLATOLIB
						gend << "\t" << freq;
#else
						gendinsert += "," + getString<float>(freq);
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							gend << "\tNA";
#else
							gendinsert += ",NULL";
#endif
						}
					}
					for (int a = 0; a < numalleles; a++) {
						if (useoverall) {
#ifndef PLATOLIB
							gend << "\t" << af.getMicroCountM(a);
#else
							gendinsert += "," + getString<int>(af.getMicroCountM(a));
#endif
						} else {
#ifndef PLATOLIB
							gend << "\t" << af.getMicroCountPM(a);
#else
							gendinsert += "," + getString<int>(af.getMicroCountPM(a));
#endif
						}
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							gend << "\tNA";
#else
							gendinsert += ",NULL";
#endif
						}
					}
					//Female
					for (int a = 0; a < numalleles; a++) {
						float freq = 0.0f;
						if (useoverall) {
							freq = ((float) af.getMicroCountF(a)
									/ (float) total_pf);
						} else {
							freq = ((float) af.getMicroCountPF(a)
									/ (float) total_pf);
						}
#ifndef PLATOLIB
						gend << "\t" << freq;
#else
						gendinsert += "," + getString<float>(freq);
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							gend << "\tNA";
#else
							gendinsert += ",NULL";
#endif
						}
					}
					for (int a = 0; a < numalleles; a++) {
						if (useoverall) {
#ifndef PLATOLIB
							gend << "\t" << af.getMicroCountF(a);
#else
							gendinsert += "," + getString<int>(af.getMicroCountF(a));
#endif
						} else {
#ifndef PLATOLIB
							gend << "\t" << af.getMicroCountPF(a);
#else
							gendinsert += "," + getString<int>(af.getMicroCountPF(a));
#endif
						}
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							gend << "\tNA";
#else
							gendinsert += ",NULL";
#endif
						}
					}
#ifndef PLATOLIB
					gend << endl;
#else
					geninsert += ")";
					//TODO: Controller
					Controller::execute_sql(myQuery, gendinsert);
#endif

					for (int l = 0; l < 15; l++) {
#ifndef PLATOLIB
						gendg << "\tNA";
#else
						gendgeninsert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					gendg << endl;
#else
					gendgeninsert += ")";
					//TODO: Controller
					Controller::execute_sql(myQuery, gendgeninsert);
#endif

				}
				if (options.doCaseControl()) {
#ifndef PLATOLIB
					ccg << mark->toString();//data_set->get_locus(k)->toString();
#endif
					int total_cam = 0;
					int total_caf = 0;
					int total_conm = 0;
					int total_conf = 0;
					for (int a = 0; a < numalleles; a++) {
						total_cam += af.getMicroCountCaM(a);
						total_caf += af.getMicroCountCaF(a);
						total_conm += af.getMicroCountConM(a);
						total_conf += af.getMicroCountConF(a);
					}
					//Case Male
					for (int a = 0; a < numalleles; a++) {
						float freq = ((float) af.getMicroCountCaM(a)
								/ (float) total_cam);
#ifndef PLATOLIB
						cc << "\t" << freq;
#else
						ccinsert += "," + getString<float>(freq);
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							cc << "\tNA";
#else
							ccinsert += ",NULL";
#endif
						}
					}
					for (int a = 0; a < numalleles; a++) {
#ifndef PLATOLIB
						cc << "\t" << af.getMicroCountCaM(a);
#else
						ccinsert += "," + getString<float>(af.getMicroCountCaM(a));
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							cc << "\tNA";
#else
							ccinsert += ",NULL";
#endif
						}
					}
					//Case Female
					for (int a = 0; a < numalleles; a++) {
						float freq = ((float) af.getMicroCountCaF(a)
								/ (float) total_caf);
#ifndef PLATOLIB
						cc << "\t" << freq;
#else
						ccinsert += "," + getString<float>(freq);
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							cc << "\tNA";
#else
							ccinsert += ",NULL";
#endif
						}
					}
					for (int a = 0; a < numalleles; a++) {
#ifndef PLATOLIB
						cc << "\t" << af.getMicroCountCaF(a);
#else
						ccinsert += "," + getString<int>(af.getMicroCountCaF(a));
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							cc << "\tNA";
#else
							ccinsert += ",NULL";
#endif
						}
					}
					//Control Male
					for (int a = 0; a < numalleles; a++) {
						float freq = ((float) af.getMicroCountConM(a)
								/ (float) total_conm);
#ifndef PLATOLIB
						cc << "\t" << freq;
#else
						ccinsert += "," + getString<float>(freq);
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							cc << "\tNA";
#else
							ccinsert += ",NULL";
#endif
						}
					}
					for (int a = 0; a < numalleles; a++) {
#ifndef PLATOLIB
						cc << "\t" << af.getMicroCountConM(a);
#else
						ccinsert += "," + getString<int>(af.getMicroCountConM(a));
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							cc << "\tNA";
#else
							ccinsert += ",NULL";
#endif
						}
					}
					//Control Female
					for (int a = 0; a < numalleles; a++) {
						float freq = ((float) af.getMicroCountConF(a)
								/ (float) total_conf);
#ifndef PLATOLIB
						cc << "\t" << freq;
#else
						ccinsert += "," + getString<float>(freq);
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							cc << "\tNA";
#else
							ccinsert += ",NULL";
#endif
						}
					}
					for (int a = 0; a < numalleles; a++) {
#ifndef PLATOLIB
						cc << "\t" << af.getMicroCountConF(a);
#else
						ccinsert += "," + getString<int>(af.getMicroCountConF(a));
#endif
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
#ifndef PLATOLIB
							cc << "\tNA";
#else
							ccinsert += ",NULL";
#endif
						}
					}
#ifndef PLATOLIB
					cc << endl;
#else
					ccinsert += ")";
					//TODO: controller
					Controller::execute_sql(myQuery, ccinsert);
#endif

					for (int l = 0; l < 27; l++) {
#ifndef PLATOLIB
						ccg << "\tNA";
#else
						ccgeninsert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					ccg << endl;
#else
					ccgeninsert += ")";
					//TODO: Controller
					Controller::execute_sql(myQuery, ccgeninsert);
#endif

				}
			} else { //not microsats
#ifndef PLATOLIB
				myoutput << mark->toString() << "\t"
						<< mark->getAllele1() << "\t"
						<< mark->getAllele2();
#else
				insert = defaultinsert;
				//insert += "," + mark->toString();
				insert += "," + getString<int>(mark->getBPLOC());
				insert += ",'" + mark->getAllele1() + "'";
				insert += ",'" + mark->getAllele2() + "'";

#endif
				if (options.doGroupFile()) {
#ifndef PLATOLIB
					gmyoutput << mark->toString() << "\t"
							<< mark->getAllele1() << "\t"
							<< mark->getAllele2();
#else
					ginsert = groupinsert;
					//ginsert += "," + mark->toString();
					ginsert += "," + getString<int>(mark->getBPLOC());
					ginsert += ",'" + mark->getAllele1() + "'";
					ginsert += ",'" + mark->getAllele2() + "'";
#endif
				}
				for (int l = 2; l < maxalleles; l++) {
#ifndef PLATOLIB
					myoutput << "\tNA";
#else
					insert += ",NULL";
#endif
					if (options.doGroupFile()) {
#ifndef PLATOLIB
						gmyoutput << "\tNA";
#else
						ginsert += ",NULL";
#endif
					}
				}

				//overall
				float majfreq = 0.0f;
				if (useoverall) {
					majfreq = af.getAone_freq();
				} else {
					majfreq = af.getAoneP_freq();
				}
				float minfreq = 1.0f - majfreq;
#ifndef PLATOLIB
				myoutput << "\t" << majfreq << "\t" << minfreq;
#else
				insert += "," + getString<float>(majfreq);
				insert += "," + getString<float>(minfreq);
#endif
				for (int l = 2; l < maxalleles; l++) {
#ifndef PLATOLIB
					myoutput << "\tNA";
#else
					insert += ",NULL";
#endif
				}
				if (useoverall) {
#ifndef PLATOLIB
					myoutput << "\t" << af.getAone_count() << "\t" << af.getAtwo_count();
#else
					insert += "," + getString<int>(af.getAone_count());
					insert += "," + getString<int>(af.getAtwo_count());
#endif
				} else {
#ifndef PLATOLIB
					myoutput << "\t" << af.getAoneP_count() << "\t" << af.getAtwoP_count();
#else
					insert += "," + getString<int>(af.getAoneP_count());
					insert += "," + getString<int>(af.getAtwoP_count());
#endif
				}
				for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
					myoutput << "\tNA";
#else
					insert += ",NULL";
#endif
				}
				//case
				majfreq = af.getAoneCa_freq();
				minfreq = 1.0f - majfreq;
#ifndef PLATOLIB
				myoutput << "\t" << majfreq << "\t" << minfreq;
#else
				insert += "," + getString<float>(majfreq);
				insert += "," + getString<float>(minfreq);
#endif
				for (int l = 2; l < maxalleles; l++) {
#ifndef PLATOLIB
					myoutput << "\tNA";
#else
					insert += ",NULL";
#endif
				}
#ifndef PLATOLIB
				myoutput << "\t" << af.getAoneCa_count() << "\t" << af.getAtwoCa_count();
#else
				insert += "," + getString<int>(af.getAoneCa_count());
				insert += "," + getString<int>(af.getAtwoCa_count());
#endif
				for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
					myoutput << "\tNA";
#else
					insert += ",NULL";
#endif
				}
				//control
				majfreq = af.getAoneCon_freq();
				minfreq = 1.0f - majfreq;
#ifndef PLATOLIB
				myoutput << "\t" << majfreq << "\t" << minfreq;
#else
				insert += "," + getString<float>(majfreq);
				insert += "," + getString<float>(minfreq);
#endif
				for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
					myoutput << "\tNA";
#else
					insert += ",NULL";
#endif
				}
#ifndef PLATOLIB
				myoutput << "\t" << af.getAoneCon_count() << "\t" << af.getAtwoCon_count();
#else
				insert += "," + getString<int>(af.getAoneCon_count());
				insert += "," + getString<int>(af.getAtwoCon_count());
#endif
				for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
					myoutput << "\tNA";
#else
					insert += ",NULL";
#endif
				}

				//groups?
				if (options.doGroupFile()) {
					int gm_total = 0;
					map<string, vector<Sample*> > groups = options.getGroups();
					map<string, vector<Sample*> >::iterator giter;
					for (giter = groups.begin(); giter != groups.end(); giter++) {
						string mygroup = giter->first;
						gm_total = af.getGroupAone_count(mygroup) + af.getGroupAtwo_count(mygroup);
						float freq = ((float) af.getGroupAone_count(mygroup)
								/ (float) gm_total);
						float freq2 = 1.0f - freq;
#ifndef PLATOLIB
						gmyoutput << "\t" << freq << "\t" << freq2;
#else
						ginsert += "," + getString<float>(freq);
						ginsert += "," + getString<float>(freq2);
#endif
						for (int b = 2; b < maxalleles; b++) {
#ifndef PLATOLIB
							gmyoutput << "\tNA";
#else
							ginsert += ",NULL";
#endif
						}
#ifndef PLATOLIB
						gmyoutput << "\t" << af.getGroupAone_count(mygroup) << "\t"
								<< af.getGroupAtwo_count(mygroup);
#else
						ginsert += "," + getString<int>(af.getGroupAone_count(mygroup));
						ginsert += "," + getString<int>(af.getGroupAtwo_count(mygroup));
#endif
						for (int b = 2; b < maxalleles; b++) {
#ifndef PLATOLIB
							gmyoutput << "\tNA";
#else
							ginsert += ",NULL";
#endif
						}

						if(freq < freq2){
							group_avg[mygroup] += freq;
						}
						else{
							group_avg[mygroup] += freq2;
						}
					}
#ifndef PLATOLIB
					gmyoutput << endl;
#else
					ginsert += ")";
					//TODO: controller
					Controller::execute_sql(myQuery, ginsert);
#endif
				}
#ifndef PLATOLIB
				myoutput << endl;
#else
				insert += ")";
				//TODO: Controller
				Controller::execute_sql(myQuery, insert);
#endif
#ifndef PLATOLIB
				mygeno << mark->toString() << "\t"
						<< mark->getAllele1() << "_"
						<< mark->getAllele1() << "\t"
						<< mark->getAllele1() << "_"
						<< mark->getAllele2() << "\t"
						<< mark->getAllele2() << "_"
						<< mark->getAllele2();
#else
				insert = defaultgenoinsert;
				//insert += "," + mark->toString();
				insert += "," + getString<int>(mark->getBPLOC());
				insert += ",'" + mark->getAllele1() + "_" + mark->getAllele1() + "'";
				insert += ",'" + mark->getAllele1() + "_" + mark->getAllele2() + "'";
				insert += ",'" + mark->getAllele2() + "_" + mark->getAllele2() + "'";
#endif
				if (options.doGroupFile()) {
#ifndef PLATOLIB
					gmygeno << mark->toString() << "\t"
							<< mark->getAllele1() << "_"
							<< mark->getAllele1() << "\t"
							<< mark->getAllele1() << "_"
							<< mark->getAllele2() << "\t"
							<< mark->getAllele2() << "_"
							<< mark->getAllele2();
#else
					ginsert = groupgenoinsert;
					//ginsert += "," + mark->toString();
					ginsert += "," + getString<int>(mark->getBPLOC());
					ginsert += ",'" + mark->getAllele1() + "_" + mark->getAllele1() + "'";
					ginsert += ",'" + mark->getAllele1() + "_" + mark->getAllele2() + "'";
					ginsert += ",'" + mark->getAllele2() + "_" + mark->getAllele2() + "'";
#endif
				}
				//overall
				float freq1 = 0.0f;
				float freq2 = 0.0f;
				float freq3 = 0.0f;
				if (useoverall) {
					freq1 = ((float) af.getAonehomo()) / (af.getPop());
					freq2 = ((float) af.getHet()) / (af.getPop());
					freq3 = ((float) af.getAtwohomo()) / (af.getPop());
#ifndef PLATOLIB
					mygeno << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomo() << "\t" << af.getHet()
							<< "\t" << af.getAtwohomo() << "\t";
#else
					insert += ",";
					insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
					insert += ",";
					insert += (isnan(freq2)  || isinf(freq2)) ? "NULL" : getString<float>(freq2);
					insert += ",";
					insert += (isnan(freq3) || isinf(freq3)) ? "NULL" : getString<float>(freq3);
					insert += "," + getString<int>(af.getAonehomo());
					insert += "," + getString<int>(af.getHet());
					insert += "," + getString<int>(af.getAtwohomo());
#endif
				} else {
					freq1 = ((float) af.getAonehomoP()) / (af.getPopP());
					freq2 = ((float) af.getHetP()) / (af.getPopP());
					freq3 = ((float) af.getAtwohomoP()) / (af.getPopP());
#ifndef PLATOLIB
					mygeno << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoP() << "\t" << af.getHetP()
							<< "\t" << af.getAtwohomoP() << "\t";
#else
					insert += ",";
					insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
					insert += ",";
					insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
					insert += ",";
					insert += (isnan(freq3) || isinf(freq3)) ? "NULL" : getString<float>(freq3);
					insert += "," + getString<int>(af.getAonehomoP());
					insert += "," + getString<int>(af.getHetP());
					insert += "," + getString<int>(af.getAtwohomoP());
#endif
				}

				//			mygeno << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countP << "\t"
				//				<< a12_countP << "\t"
				//				<< a2_homo_countP << "\t";
				//case overall
				freq1 = ((float) af.getAonehomoCa()) / (af.getPopCa());
				freq2 = ((float) af.getHetCa()) / (af.getPopCa());
				freq3 = ((float) af.getAtwohomoCa()) / (af.getPopCa());
#ifndef PLATOLIB
				mygeno << freq1 << "\t" << freq2 << "\t" << freq3 << "\t"
						<< af.getAonehomoCa() << "\t" << af.getHetCa() << "\t"
						<< af.getAtwohomoCa() << "\t";
#else
				insert += ",";
				insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
				insert += ",";
				insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
				insert += ",";
				insert += (isnan(freq3) || isinf(freq3)) ? "NULL" : getString<float>(freq3);
				insert += "," + getString<int>(af.getAonehomoCa());
				insert += "," + getString<int>(af.getHetCa());
				insert += "," + getString<int>(af.getAtwohomoCa());
#endif
				//control overall
				freq1 = ((float) af.getAonehomoCon()) / (af.getPopCon());
				freq2 = ((float) af.getHetCon()) / (af.getPopCon());
				freq3 = ((float) af.getAtwohomoCon()) / (af.getPopCon());
#ifndef PLATOLIB
				mygeno << freq1 << "\t" << freq2 << "\t" << freq3 << "\t"
						<< af.getAonehomoCon() << "\t" << af.getHetCon() << "\t"
						<< af.getAtwohomoCon();
#else
				insert += ",";
				insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
				insert += ",";
				insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
				insert += ",";
				insert += (isnan(freq3) || isinf(freq3)) ? "NULL" : getString<float>(freq3);
				insert += "," + getString<int>(af.getAonehomoCon());
				insert += "," + getString<int>(af.getHetCon());
				insert += "," + getString<int>(af.getAtwohomoCon());
#endif

				//groups?
				if (options.doGroupFile()) {
////					int gm_total = 0;
					map<string, vector<Sample*> > groups = options.getGroups();
					map<string, vector<Sample*> >::iterator giter;
					for (giter = groups.begin(); giter != groups.end(); giter++) {
						string mygroup = giter->first;
						float genotot = af.getGroupPop(mygroup);
						freq1 = ((float) af.getGroupAonehomo(mygroup) / genotot);
						freq2 = ((float) af.getGroupHet(mygroup) / genotot);
						freq3 = ((float) af.getGroupAtwohomo(mygroup) / genotot);
#ifndef PLATOLIB
						gmygeno << "\t" << freq1 << "\t" << freq2 << "\t"
								<< freq3 << "\t" << af.getGroupAonehomo(mygroup)
								<< "\t" << af.getGroupHet(mygroup) << "\t"
								<< af.getGroupAtwohomo(mygroup);
#else
						ginsert += ",";
						ginsert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
						ginsert += ",";
						ginsert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
						ginsert += ",";
						ginsert += (isnan(freq3) || isinf(freq3)) ? "NULL" : getString<float>(freq3);
						ginsert += "," + getString<int>(af.getGroupAonehomo(mygroup));
						ginsert += "," + getString<int>(af.getGroupHet(mygroup));
						ginsert += "," + getString<int>(af.getGroupAtwohomo(mygroup));
#endif
					}
#ifndef PLATOLIB
					gmygeno << endl;
#else
					ginsert += ")";
					//TODO: Controller
					Controller::execute_sql(myQuery, ginsert);
#endif
				}
#ifndef PLATOLIB
				mygeno << endl;
#else
				insert += ")";
				//TODO: Controller
				Controller::execute_sql(myQuery, insert);
#endif

				if (options.doParental()) {
#ifndef PLATOLIB
					paren << mark->toString() << "\t"
							<< mark->getAllele1() << "\t"
							<< mark->getAllele2();
#else
					insert = parentalinsert;
					insert + "," + mark->toString();
					insert += ",'" + mark->getAllele1() + "'";
					insert += ",'" + mark->getAllele2() + "'";
#endif
					for (int l = 2; l < maxalleles; l++) {
#ifndef PLATOLIB
						paren << "\tNA";
#else
						insert += ",NULL";
#endif
					}
					//parent male
					float majfreq = af.getAonePM_freq();
					float minfreq = 1.0f - majfreq;
					if(af.getAonePM_count() == 0 && af.getAtwoPM_count() == 0){
						majfreq = 0;
						minfreq = 0;
					}
#ifndef PLATOLIB
					paren << "\t" << majfreq << "\t" << minfreq;
#else
					insert += "," + getString<float>(majfreq);
					insert += "," + getString<float>(minfreq);
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						paren << "\tNA";
#else
						insert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					paren << "\t" << af.getAonePM_count() << "\t" << af.getAtwoPM_count();
#else
					insert += "," + getString<int>(af.getAonePM_count());
					insert += "," + getString<int>(af.getAtwoPM_count());
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						paren << "\tNA";
#else
						insert += ",NULL";
#endif
					}
					//parent female
					majfreq = af.getAonePF_freq();
					minfreq = 1.0f - majfreq;
#ifndef PLATOLIB
					paren << "\t" << majfreq << "\t" << minfreq;
#else
					insert += "," + getString<float>(majfreq);
					insert += "," + getString<float>(minfreq);
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						paren << "\tNA";
#else
						insert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					paren << "\t" << af.getAonePF_count() << "\t" << af.getAtwoPF_count();
#else
					insert += "," + getString<int>(af.getAonePF_count());
					insert += "," + getString<int>(af.getAtwoPF_count());
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						paren << "\tNA";
#else
						insert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					paren << endl;

					pareng << mark->toString() << "\t"
					<< mark->getAllele1() << "_"
					<< mark->getAllele1() << "\t"
					<< mark->getAllele1() << "_"
					<< mark->getAllele2() << "\t"
					<< mark->getAllele2() << "_"
					<< mark->getAllele2();
#else
					insert += ")";
					//TODO: Controller
					Controller::execute_sql(myQuery, insert);

					insert = parentalgenoinsert;
					//insert += "," + mark->toString();
					insert += "," + getString<int>(mark->getBPLOC());
					insert += ",'" + mark->getAllele1() + "_" + mark->getAllele1() + "'";
					insert += ",'" + mark->getAllele1() + "_" + mark->getAllele2() + "'";
					insert += ",'" + mark->getAllele2() + "_" + mark->getAllele2() + "'";
#endif


					float freq1 = ((float) af.getAonehomoPM()) / (af.getPopPM());
					float freq2 = ((float) af.getHetPM()) / (af.getPopPM());
					float freq3 = ((float) af.getAtwohomoPM()) / (af.getPopPM());
#ifndef PLATOLIB
					pareng << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoPM() << "\t" << af.getHetPM()
							<< "\t" << af.getAtwohomoPM();
#else
					insert += ",";
					insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
					insert += ",";
					insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
					insert += ",";
					insert += (isnan(freq3) || isinf(freq3)) ? "NULL" : getString<float>(freq3);
					insert += "," + getString<int>(af.getAonehomoPM());
					insert += "," + getString<int>(af.getHetPM());
					insert += "," + getString<int>(af.getAtwohomoPM());
#endif
					freq1 = ((float) af.getAonehomoPF()) / (af.getPopPF());
					freq2 = ((float) af.getHetPF()) / (af.getPopPF());
					freq3 = ((float) af.getAtwohomoPF()) / (af.getPopPF());
#ifndef PLATOLIB
					pareng << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoPF() << "\t" << af.getHetPF()
							<< "\t" << af.getAtwohomoPF();
					pareng << endl;
#else
					insert += ",";
					insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
					insert += ",";
					insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
					insert += ",";
					insert += (isnan(freq3) || isinf(freq3)) ? "NULL" : getString<float>(freq3);
					insert += "," + getString<int>(af.getAonehomoPF());
					insert += "," + getString<int>(af.getHetPF());
					insert += "," + getString<int>(af.getAtwohomoPF());

					insert += ")";
					//TODO: Controller
					Controller::execute_sql(myQuery, insert);
#endif
				}
				if (options.doGender()) {
					//overall male
#ifndef PLATOLIB
					gend << mark->toString() << "\t"
							<< mark->getAllele1() << "\t"
							<< mark->getAllele2();
#else
					insert = genderinsert;
					//insert += "," + mark->toString();
					insert += "," + getString<int>(mark->getBPLOC());
					insert += ",'" + mark->getAllele1() + "'";
					insert += ",'" + mark->getAllele2() + "'";
#endif
					for (int l = 2; l < maxalleles; l++) {
#ifndef PLATOLIB
						gend << "\tNA";
#else
						insert += ",NULL";
#endif
					}
					float majfreq = 0.0f;
					if (useoverall) {
						majfreq = af.getAoneM_freq();
					} else {
						majfreq = af.getAonePM_freq();
					}
					float minfreq = 1.0f - majfreq;
					if(useoverall && af.getAoneM_count() == 0 && af.getAtwoM_count() == 0){
						majfreq = 0;
						minfreq = 0;
					}else if(!useoverall && af.getAonePM_count() == 0 && af.getAtwoPM_count() == 0){
						majfreq = 0;
						minfreq = 0;
					}
#ifndef PLATOLIB
					gend << "\t" << majfreq << "\t" << minfreq;
#else
					insert += "," + getString<float>(majfreq);
					insert += "," + getString<float>(minfreq);
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						gend << "\tNA";
#else
						insert += ",NULL";
#endif
					}
					if (useoverall) {
#ifndef PLATOLIB
						gend << "\t" << af.getAoneM_count() << "\t" << af.getAtwoM_count();
#else
						insert += "," + getString<int>(af.getAoneM_count());
						insert += "," + getString<int>(af.getAtwoM_count());
#endif
					} else {
#ifndef PLATOLIB
						gend << "\t" << af.getAonePM_count() << "\t" << af.getAtwoPM_count();
#else
						insert += "," + getString<int>(af.getAonePM_count());
						insert += "," + getString<int>(af.getAtwoPM_count());
#endif
					}
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						gend << "\tNA";
#else
						insert += ",NULL";
#endif
					}
					//overall female
					if (useoverall) {
						majfreq = af.getAoneF_freq();
					} else {
						majfreq = af.getAonePF_freq();
					}
					minfreq = 1.0f - majfreq;
#ifndef PLATOLIB
					gend << "\t" << majfreq << "\t" << minfreq;
#else
					insert += "," + getString<float>(majfreq);
					insert += "," + getString<float>(minfreq);
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						gend << "\tNA";
#else
						insert += ",NULL";
#endif
					}
					if (useoverall) {
#ifndef PLATOLIB
						gend << "\t" << af.getAoneF_count() << "\t" << af.getAtwoF_count();
#else
						insert += "," + getString<int>(af.getAoneF_count());
						insert += "," + getString<int>(af.getAtwoF_count());
#endif
					} else {
#ifndef PLATOLIB
						gend << "\t" << af.getAonePF_count() << "\t" << af.getAtwoPF_count();
#else
						insert += "," + getString<int>(af.getAonePF_count());
						insert += "," + getString<int>(af.getAtwoPF_count());
#endif
					}
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						gend << "\tNA";
#else
						insert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					gend << endl;

					gendg << mark->toString() << "\t"
							<< mark->getAllele1() << "_"
							<< mark->getAllele1() << "\t"
							<< mark->getAllele1() << "_"
							<< mark->getAllele2() << "\t"
							<< mark->getAllele2() << "_"
							<< mark->getAllele2();
#else
					insert += ")";
					//TODO: Controller
					Controller::execute_sql(myQuery, insert);

					insert = gendergenoinsert;
					//insert += "," + mark->toString();
					insert += "," + getString<int>(mark->getBPLOC());
					insert += ",'" + mark->getAllele1() + "_" + mark->getAllele1() + "'";
					insert += ",'" + mark->getAllele1() + "_" + mark->getAllele2() + "'";
					insert += ",'" + mark->getAllele2() + "_" + mark->getAllele2() + "'";
#endif
					if (useoverall) {
						float freq1 =
								((float) af.getAonehomoM()) / (af.getPopM());
						float freq2 = ((float) af.getHetM()) / (af.getPopM());
						float freq3 =
								((float) af.getAtwohomoM()) / (af.getPopM());
#ifndef PLATOLIB
						gendg << "\t" << freq1 << "\t" << freq2 << "\t"
								<< freq3 << "\t" << af.getAonehomoM() << "\t"
								<< af.getHetM() << "\t" << af.getAtwohomoM();
#else
						insert += ",";
						insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
						insert += ",";
						insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
						insert += ",";
						insert += (isnan(freq3) || isinf(freq3)) ? "NULL" : getString<float>(freq3);
						insert += "," + getString<int>(af.getAonehomoM());
						insert += "," + getString<int>(af.getHetM());
						insert += "," + getString<int>(af.getAtwohomoM());
#endif
						freq1 = ((float) af.getAonehomoF()) / (af.getPopF());
						freq2 = ((float) af.getHetF()) / (af.getPopF());
						freq3 = ((float) af.getAtwohomoF()) / (af.getPopF());
#ifndef PLATOLIB
						gendg << "\t" << freq1 << "\t" << freq2 << "\t"
								<< freq3 << "\t" << af.getAonehomoF() << "\t"
								<< af.getHetF() << "\t" << af.getAtwohomoF();
#else
						insert += ",";
						insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
						insert += ",";
						insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
						insert += ",";
						insert += (isnan(freq3) || isinf(freq3)) ? "NULL" : getString<float>(freq3);
						insert += "," + getString<int>(af.getAonehomoF());
						insert += "," + getString<int>(af.getHetF());
						insert += "," + getString<int>(af.getAtwohomoF());
#endif
					} else {
						float freq1 = ((float) af.getAonehomoPM())
								/ (af.getPopPM());
						float freq2 = ((float) af.getHetPM()) / (af.getPopPM());
						float freq3 = ((float) af.getAtwohomoPM())
								/ (af.getPopPM());
#ifndef PLATOLIB
						gendg << "\t" << freq1 << "\t" << freq2 << "\t"
								<< freq3 << "\t" << af.getAonehomoPM() << "\t"
								<< af.getHetPM() << "\t" << af.getAtwohomoPM();
#else
						insert += ",";
						insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
						insert += ",";
						insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
						insert += ",";
						insert += (isnan(freq3) || isinf(freq3)) ? "NULL" : getString<float>(freq3);
						insert += "," + getString<int>(af.getAonehomoPM());
						insert += "," + getString<int>(af.getHetPM());
						insert += "," + getString<int>(af.getAtwohomoPM());
#endif
						freq1 = ((float) af.getAonehomoPF()) / (af.getPopPF());
						freq2 = ((float) af.getHetPF()) / (af.getPopPF());
						freq3 = ((float) af.getAtwohomoPF()) / (af.getPopPF());
#ifndef PLATOLIB
						gendg << "\t" << freq1 << "\t" << freq2 << "\t"
								<< freq3 << "\t" << af.getAonehomoPF() << "\t"
								<< af.getHetPF() << "\t" << af.getAtwohomoPF();
#else
						insert += ",";
						insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
						insert += ",";
						insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
						insert += ",";
						insert += (isnan(freq3) || isinf(freq3)) ? "NULL" : getString<float>(freq3);
						insert += "," + getString<int>(af.getAonehomoPF());
						insert += "," + getString<int>(af.getHetPF());
						insert += "," + getString<int>(af.getAtwohomoPF());
#endif
					}
#ifndef PLATOLIB
					gendg << endl;
#else
					insert += ")";
					//TODO: Controller
					Controller::execute_sql(myQuery, insert);
#endif

				}
				if (options.doCaseControl()) {
					//case male
#ifndef PLATOLIB
					cc << mark->toString() << "\t"
							<< mark->getAllele1() << "\t"
							<< mark->getAllele2();
#else
					insert = casecontrolinsert;
					//insert += "," + mark->toString();
					insert += "," + getString<int>(mark->getBPLOC());
					insert += ",'" + mark->getAllele1() + "'";
					insert += ",'" + mark->getAllele2() + "'";
#endif
					for (int l = 2; l < maxalleles; l++) {
#ifndef PLATOLIB
						cc << "\tNA";
#else
						insert += ",NULL";
#endif
					}
					float majfreq = af.getAoneCaM_freq();
					float minfreq = 1.0f - majfreq;
					if(af.getAoneCaM_count() == 0 && af.getAtwoCaM_count() == 0){
						majfreq = 0;
						minfreq = 0;
					}
#ifndef PLATOLIB
					cc << "\t" << majfreq << "\t" << minfreq;
#else
					insert += "," + getString<float>(majfreq);
					insert += "," + getString<float>(minfreq);
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						cc << "\tNA";
#else
						insert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					cc << "\t" << af.getAoneCaM_count() << "\t" << af.getAtwoCaM_count();
#else
					insert += "," + getString<int>(af.getAoneCaM_count());
					insert += "," + getString<int>(af.getAtwoCaM_count());
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						cc << "\tNA";
#else
						insert += ",NULL";
#endif
					}
					//case female
					majfreq = af.getAoneCaF_freq();
					minfreq = 1.0f - majfreq;
#ifndef PLATOLIB
					cc << "\t" << majfreq << "\t" << minfreq;
#else
					insert += "," + getString<float>(majfreq);
					insert += "," + getString<float>(majfreq);
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						cc << "\tNA";
#else
						insert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					cc << "\t" << af.getAoneCaF_count() << "\t" << af.getAtwoCaF_count();
#else
					insert += "," + getString<int>(af.getAoneCaF_count());
					insert += "," + getString<int>(af.getAtwoCaF_count());
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						cc << "\tNA";
#else
						insert += ",NULL";
#endif
					}
					//control male
					majfreq = af.getAoneConM_freq();
					minfreq = 1.0f - majfreq;
					if(af.getAoneConM_count() == 0 && af.getAtwoConM_count() == 0){
						majfreq = 0;
						minfreq = 0;
					}
#ifndef PLATOLIB
					cc << "\t" << majfreq << "\t" << minfreq;
#else
					insert + "," + getString<float>(majfreq);
					insert += "," + getString<float>(minfreq);
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						cc << "\tNA";
#else
						insert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					cc << "\t" << af.getAoneConM_count() << "\t" << af.getAtwoConM_count();
#else
					insert += "," + getString<int>(af.getAoneConM_count());
					insert += "," + getString<int>(af.getAtwoConM_count());
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						cc << "\tNA";
#else
						insert += ",NULL";
#endif
					}
					//control female
					majfreq = af.getAoneConF_freq();
					minfreq = 1.0f - majfreq;
#ifndef PLATOLIB
					cc << "\t" << majfreq << "\t" << minfreq;
#else
					insert += "," + getString<float>(majfreq);
					insert += "," + getString<float>(minfreq);
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						cc << "\tNA";
#else
						insert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					cc << "\t" << af.getAoneConF_count() << "\t" << af.getAtwoConF_count();
#else
					insert += "," + getString<int>(af.getAoneConF_count());
					insert += "," + getString<int>(af.getAtwoConF_count());
#endif
					for (int t = 2; t < maxalleles; t++) {
#ifndef PLATOLIB
						cc << "\tNA";
#else
						insert += ",NULL";
#endif
					}
#ifndef PLATOLIB
					cc << endl;

					ccg << mark->toString() << "\t"
							<< mark->getAllele1() << "_"
							<< mark->getAllele1() << "\t"
							<< mark->getAllele1() << "_"
							<< mark->getAllele2() << "\t"
							<< mark->getAllele2() << "_"
							<< mark->getAllele2();
#else
					insert += ")";
					//TODO: controller
					Controller::execute_sql(myQuery, insert);

					insert = casecontrolgenoinsert;
					//insert += "," + mark->toString();
					insert += "," + getString<int>(mark->getBPLOC());
					insert += ",'" + mark->getAllele1() + "_" + mark->getAllele1() + "'";
					insert += ",'" + mark->getAllele1() + "_" + mark->getAllele2() + "'";
					insert += ",'" + mark->getAllele2() + "_" + mark->getAllele2() + "'";
#endif
					float freq1 = ((float) af.getAonehomoCaM())
							/ (af.getPopCaM());
					float freq2 = ((float) af.getHetCaM()) / (af.getPopCaM());
					float freq3 = ((float) af.getAtwohomoCaM())
							/ (af.getPopCaM());
#ifndef PLATOLIB
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoCaM() << "\t" << af.getHetCaM()
							<< "\t" << af.getAtwohomoCaM();
#else
					insert += ",";
					insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
					insert += ",";
					insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
					insert += ",";
					insert += (isnan(freq3)|| isinf(freq3)) ? "NULL" : getString<float>(freq3);
					insert += "," + getString<int>(af.getAonehomoCaM());
					insert += "," + getString<int>(af.getHetCaM());
					insert += "," + getString<int>(af.getAtwohomoCaM());
#endif
					freq1 = ((float) af.getAonehomoCaF()) / (af.getPopCaF());
					freq2 = ((float) af.getHetCaF()) / (af.getPopCaF());
					freq3 = ((float) af.getAtwohomoCaF()) / (af.getPopCaF());
#ifndef PLATOLIB
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoCaF() << "\t" << af.getHetCaF()
							<< "\t" << af.getAtwohomoCaF();
#else
					insert += ",";
					insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
					insert += ",";
					insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
					insert += ",";
					insert += (isnan(freq3)|| isinf(freq3)) ? "NULL" : getString<float>(freq3);
					insert += "," + getString<int>(af.getAonehomoCaF());
					insert += "," + getString<int>(af.getHetCaF());
					insert += "," + getString<int>(af.getAtwohomoCaF());
#endif
					freq1 = ((float) af.getAonehomoConM()) / (af.getPopConM());
					freq2 = ((float) af.getHetConM()) / (af.getPopConM());
					freq3 = ((float) af.getAtwohomoConM()) / (af.getPopConM());
#ifndef PLATOLIB
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoConM() << "\t"
							<< af.getHetConM() << "\t" << af.getAtwohomoConM();
#else
					insert += ",";
					insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
					insert += ",";
					insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
					insert += ",";
					insert += (isnan(freq3)|| isinf(freq3)) ? "NULL" : getString<float>(freq3);
					insert += "," + getString<int>(af.getAonehomoConM());
					insert += "," + getString<int>(af.getHetConM());
					insert += "," + getString<int>(af.getAtwohomoConM());
#endif
					freq1 = ((float) af.getAonehomoConF()) / (af.getPopConF());
					freq2 = ((float) af.getHetConF()) / (af.getPopConF());
					freq3 = ((float) af.getAtwohomoConF()) / (af.getPopConF());
#ifndef PLATOLIB
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoConF() << "\t"
							<< af.getHetConF() << "\t" << af.getAtwohomoConF();
					ccg << endl;
#else
					insert += ",";
					insert += (isnan(freq1) || isinf(freq1)) ? "NULL" : getString<float>(freq1);
					insert += ",";
					insert += (isnan(freq2) || isinf(freq2)) ? "NULL" : getString<float>(freq2);
					insert += ",";
					insert += (isnan(freq3)|| isinf(freq3)) ? "NULL" : getString<float>(freq3);
					insert += "," + getString<int>(af.getAonehomoConF());
					insert += "," + getString<int>(af.getHetConF());
					insert += "," + getString<int>(af.getAtwohomoConF());

					insert += ")";
					//TODO: Controller
					Controller::execute_sql(myQuery, insert);
#endif

				}
			}

			//filter Markers
			//doFilter(mark, &af);
		}
	}
//TODO:  not sure if this block should have a database version or not...
#ifndef PLATOLIB
	if(options.doGroupFile()){
		map<string, vector<Sample*> > groups = options.getGroups();
		map<string, vector<Sample*> >::iterator giter;
		for (giter = groups.begin(); giter != groups.end(); giter++) {
			string mygroup = giter->first;
			double val = group_avg[mygroup];
			val = val / (double) total_snps;
			int goodsamps = 0;
			for(int gs = 0; gs < (int)groups[mygroup].size(); gs++){
				if(groups[mygroup][gs]->isEnabled()){
					goodsamps++;
				}
			}
			gmyoutputavg << mygroup << "\t" << val << "\t" << goodsamps << endl;
		}
	}

#endif
//need to commit the above DB statements to the database
#ifdef PLATOLIB
	myQuery.commit();
#endif
}

/*
 *Function: process
 *Description:
 *Main method to begin the whole process.  Flags samples then diverts work to processtest
 *
 *
 */
void ProcessAlleleFrequency::process(DataSet* ds)
{
	data_set = ds;
	if (options.doGroupFile()) {
		options.readGroups(data_set->get_samples());
	}

	cout << "Calling processtest()\n";
	processtest();
	return;

}

#ifdef PLATOLIB
void ProcessAlleleFrequency::dump2db(){}

void ProcessAlleleFrequency::create_tables()
{
    Query myQuery(*db);
    int msize = data_set->num_loci();

    int maxalleles = 0;
    for (int i = 0; i < msize; i++) {
        if (data_set->get_locus(i)->isEnabled()) {
                if (data_set->get_locus(i)->getNumAlleles() > maxalleles) {
                        maxalleles = data_set->get_locus(i)->getNumAlleles();
                }
        }
    }

    for(int i = 0; i < (int)tablename.size(); i++){
        Controller::drop_table(db, tablename[i]);
    }
    headers.clear();
    tablename.clear();
    primary_table.clear();

    string tempbatch = batchname;
    for(int i = 0; i < (int)tempbatch.size(); i++){
        if(tempbatch[i] == ' '){
            tempbatch[i] = '_';
        }
    }
    string mytablename = tempbatch + "_";
    tempbatch = name;
    for(int i = 0; i < (int)tempbatch.size(); i++){
        if(tempbatch[i] == ' '){
            tempbatch[i] = '_';
        }
    }
    string base = mytablename + tempbatch;

    mytablename = base + "_" + getString<int>(position);
    tablename.push_back(mytablename);
    tablenicknames.push_back("");
    primary_table[mytablename].push_back("LOCI");
    string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
    defaultinsert = "INSERT INTO " + mytablename + "(id, fkey";
    sql += "fkey integer not null,";
    for(int i = 0; i < maxalleles; i++){
        sql += "Allele" + getString<int>(i + 1) + " varchar(10),";
        defaultinsert += ",Allele" + getString<int>(i + 1);
    }
    for(int i = 0; i < maxalleles; i++){
        sql += "Overall_Allele" + getString<int>(i + 1) + "_freq REAL,";
        headers[mytablename].push_back("Overall_Allele" + getString<int>(i + 1) + "_freq");
        defaultinsert += ",Overall_Allele" + getString<int>(i + 1) + "_freq";
    }
    for(int i = 0; i < maxalleles; i++){
        sql += "Overall_Allele" + getString<int>(i + 1) + "_count integer,";
        headers[mytablename].push_back("Overall_Allele" + getString<int>(i + 1) + "_count");
        defaultinsert += ",Overall_Allele" + getString<int>(i + 1) + "_count";
    }
    for(int i = 0; i < maxalleles; i++){
        sql += "Case_Allele" + getString<int>(i + 1) + "_freq REAL,";
        headers[mytablename].push_back("Case_Allele" + getString<int>(i + 1) + "_freq");
        defaultinsert += ",Case_Allele" + getString<int>(i + 1) + "_freq";
    }
    for(int i = 0; i < maxalleles; i++){
        sql += "Case_Allele" + getString<int>(i + 1) + "_count integer,";
        headers[mytablename].push_back("Case_Allele" + getString<int>(i + 1) + "_count");
        defaultinsert += ",Case_Allele" + getString<int>(i + 1) + "_count";
    }
    for(int i = 0; i < maxalleles; i++){
        sql += "Control_Allele" + getString<int>(i + 1) + "_freq REAL,";
        headers[mytablename].push_back("Control_Allele" + getString<int>(i + 1) + "_freq");
        defaultinsert += ",Control_Allele" + getString<int>(i + 1) + "_freq";
    }
    for(int i = 0; i < maxalleles; i++){
        sql += "Control_Allele" + getString<int>(i + 1) + "_count integer,";
        headers[mytablename].push_back("Control_Allele" + getString<int>(i + 1) + "_count");
        defaultinsert += ",Control_Allele" + getString<int>(i + 1) + "_count";
    }
    defaultinsert += ") VALUES (NULL";
    sql = sql.replace(sql.size() - 1, 1, ")");

    //Controller::execute_sql(db, sql);
    myQuery.transaction();
    Controller::execute_sql(myQuery, sql);
    myQuery.commit();

    mytablename = base + "_genotype_" + getString<int>(position);
    tablename.push_back(mytablename);
    tablenicknames.push_back("Genotype");
    primary_table[mytablename].push_back("LOCI");
    sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
    sql += "fkey integer not null,";
    sql += "Genotype11 varchar(20),";
    sql += "Genotype12 varchar(20),";
    sql += "Genotype22 varchar(20),";
    defaultgenoinsert = "INSERT INTO " + mytablename + " (id, fkey, Genotype11, Genotype12, Genotype22";
    sql += "Overall_Freq_Genotype11 REAL,";
    headers[mytablename].push_back("Overall_Freq_Genotype11");
    defaultgenoinsert += ",Overall_Freq_Genotype11";
    sql += "Overall_Freq_Genotype12 REAL,";
    headers[mytablename].push_back("Overall_Freq_Genotype12");
    defaultgenoinsert += ",Overall_Freq_Genotype12";
    sql += "Overall_Freq_Genotype22 REAL,";
    headers[mytablename].push_back("Overall_Freq_Genotype22");
    defaultgenoinsert += ",Overall_Freq_Genotype22";
    sql += "Overall_Count_Genotype11 integer,";
    headers[mytablename].push_back("Overall_Count_Genotype11");
    defaultgenoinsert += ",Overall_Count_Genotype11";
    sql += "Overall_Count_Genotype12 integer,";
    headers[mytablename].push_back("Overall_Count_Genotype12");
    defaultgenoinsert += ",Overall_Count_Genotype12";
    sql += "Overall_Count_Genotype22 integer,";
    headers[mytablename].push_back("Overall_Count_Genotype22");
    defaultgenoinsert += ",Overall_Count_Genotype22";
    sql += "Case_Freq_Genotype11 REAL,";
    headers[mytablename].push_back("Case_Freq_Genotype11");
    defaultgenoinsert += ",Case_Freq_Genotype11";
    sql += "Case_Freq_Genotype12 REAL,";
    headers[mytablename].push_back("Case_Freq_Genotype12");
    defaultgenoinsert += ",Case_Freq_Genotype12";
    sql += "Case_Freq_Genotype22 REAL,";
    headers[mytablename].push_back("Case_Freq_Genotype22");
    defaultgenoinsert += ",Case_Freq_Genotype22";
    sql += "Case_Count_Genotype11 integer,";
    headers[mytablename].push_back("Case_Count_Genotype11");
    defaultgenoinsert += ",Case_Count_Genotype11";
    sql += "Case_Count_Genotype12 integer,";
    headers[mytablename].push_back("Case_Count_Genotype12");
    defaultgenoinsert += ",Case_Count_Genotype12";
    sql += "Case_Count_Genotype22 integer,";
    headers[mytablename].push_back("Case_Count_Genotype22");
    defaultgenoinsert += ",Case_Count_Genotype22";
    sql += "Control_Freq_Genotype11 REAL,";
    headers[mytablename].push_back("Control_Freq_Genotype11");
    defaultgenoinsert += ",Control_Freq_Genotype11";
    sql += "Control_Freq_Genotype12 REAL,";
    headers[mytablename].push_back("Control_Freq_Genotype12");
    defaultgenoinsert += ",Control_Freq_Genotype12";
    sql += "Control_Freq_Genotype22 REAL,";
    headers[mytablename].push_back("Control_Freq_Genotype22");
    defaultgenoinsert += ",Control_Freq_Genotype22";
    sql += "Control_Count_Genotype11 integer,";
    headers[mytablename].push_back("Control_Count_Genotype11");
    defaultgenoinsert += ",Control_Count_Genotype11";
    sql += "Control_Count_Genotype12 integer,";
    headers[mytablename].push_back("Control_Count_Genotype12");
    defaultgenoinsert += ",Control_Count_Genotype12";
    sql += "Control_Count_Genotype22 integer)";
    headers[mytablename].push_back("Control_Count_Genotype22");
    defaultgenoinsert += ",Control_Count_Genotype22";

    defaultgenoinsert += ") VALUES (NULL";
    //Controller::execute_sql(db, sql);
    myQuery.transaction();
    Controller::execute_sql(myQuery, sql);
    myQuery.commit();

    if(options.doCaseControl()){
        mytablename = base + "_casecontrol_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("Case/Control");
        primary_table[mytablename].push_back("LOCI");
        sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        casecontrolinsert = "INSERT INTO " + mytablename + " (id, fkey";
        for(int i = 0; i < maxalleles; i++){
            sql += "Allele" + getString<int>(i + 1) + " varchar(10),";
            casecontrolinsert += ",Allele" + getString<int>(i + 1);
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Case_Male_Allele" + getString<int>(i + 1) + "_freq REAL,";
            headers[mytablename].push_back("Case_Male_Allele" + getString<int>(i + 1) + "_freq");
            casecontrolinsert += ",Case_Male_Allele" + getString<int>(i + 1) + "_freq";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Case_Male_Allele" + getString<int>(i + 1) + "_count integer,";
            headers[mytablename].push_back("Case_Male_Allele" + getString<int>(i + 1) + "_count");
            casecontrolinsert += ",Case_Male_Allele" + getString<int>(i + 1) + "_count";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Case_Female_Allele" + getString<int>(i + 1) + "_freq REAL,";
            headers[mytablename].push_back("Case_Female_Allele" + getString<int>(i + 1) + "_freq");
            casecontrolinsert += ",Case_Female_Allele" + getString<int>(i + 1) + "_freq";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Case_Female_Allele" + getString<int>(i + 1) + "_count integer,";
            headers[mytablename].push_back("Case_Female_Allele" + getString<int>(i + 1) + "_count");
            casecontrolinsert += ",Case_Female_Allele" + getString<int>(i + 1) + "_count";
        }

        for(int i = 0; i < maxalleles; i++){
            sql += "Control_Male_Allele" + getString<int>(i + 1) + "_freq REAL,";
            headers[mytablename].push_back("Control_Male_Allele" + getString<int>(i + 1) + "_freq");
            casecontrolinsert += ",Control_Male_Allele" + getString<int>(i + 1) + "_freq";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Control_Male_Allele" + getString<int>(i + 1) + "_count integer,";
            headers[mytablename].push_back("Control_Male_Allele" + getString<int>(i + 1) + "_count");
            casecontrolinsert += ",Control_Male_Allele" + getString<int>(i + 1) + "_count";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Control_Female_Allele" + getString<int>(i + 1) + "_freq REAL,";
            headers[mytablename].push_back("Control_Female_Allele" + getString<int>(i + 1) + "_freq");
            casecontrolinsert += ",Control_Female_Allele" + getString<int>(i + 1) + "_freq";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Control_Female_Allele" + getString<int>(i + 1) + "_count integer,";
            headers[mytablename].push_back("Control_Female_Allele" + getString<int>(i + 1) + "_count");
            casecontrolinsert += ",Control_Female_Allele" + getString<int>(i + 1) + "_count";
        }
        casecontrolinsert += ") VALUES (NULL";
        sql = sql.replace(sql.size() - 1, 1, ")");

        //Controller::execute_sql(db, sql);
        myQuery.transaction();
        Controller::execute_sql(myQuery, sql);
        myQuery.commit();

        mytablename = base + "_genotype_casecontrol_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("Case/Control Genotype");
        primary_table[mytablename].push_back("LOCI");
        sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        sql += "Genotype11 varchar(20),";
    sql += "Genotype12 varchar(20),";
    sql += "Genotype22 varchar(20),";
        casecontrolgenoinsert = "INSERT INTO " + mytablename + "(id, fkey, Genotype11, Genotype12, Genotype22";
    sql += "Case_Male_Freq_Genotype11 REAL,";
    headers[mytablename].push_back("Case_Male_Genotype11");
    casecontrolgenoinsert += ",Case_Male_Freq_Genotype11";
    sql += "Case_Male_Freq_Genotype12 REAL,";
    headers[mytablename].push_back("Case_Male_Freq_Genotype12");
    casecontrolgenoinsert += ",Case_Male_Freq_Genotype12";
    sql += "Case_Male_Freq_Genotype22 REAL,";
    headers[mytablename].push_back("Case_Male_Freq_Genotype22");
    casecontrolgenoinsert += ",Case_Male_Freq_Genotype22";
    sql += "Case_Male_Count_Genotype11 integer,";
    headers[mytablename].push_back("Case_Male_Count_Genotype11");
    casecontrolgenoinsert += ",Case_Male_Count_Genotype11";
    sql += "Case_Male_Count_Genotype12 integer,";
    headers[mytablename].push_back("Case_Male_Count_Genotype12");
    casecontrolgenoinsert += ",Case_Male_Count_Genotype12";
    sql += "Case_Male_Count_Genotype22 integer,";
    headers[mytablename].push_back("Case_Male_Count_Genotype22");
    casecontrolgenoinsert += ",Case_Male_Count_Genotype22";
    sql += "Case_Female_Freq_Genotype11 REAL,";
    headers[mytablename].push_back("Case_Female_Freq_Genotype11");
    casecontrolgenoinsert += ",Case_Female_Freq_Genotype11";
    sql += "Case_Female_Freq_Genotype12 REAL,";
    headers[mytablename].push_back("Case_Female_Freq_Genotype12");
    casecontrolgenoinsert += ",Case_Female_Freq_Genotype12";
    sql += "Case_Female_Freq_Genotype22 REAL,";
    headers[mytablename].push_back("Case_Female_Freq_Genotype22");
    casecontrolgenoinsert += ",Case_Female_Freq_Genotype22";
    sql += "Case_Female_Count_Genotype11 integer,";
    headers[mytablename].push_back("Case_Female_Count_Genotype11");
    casecontrolgenoinsert += ",Case_Female_Count_Genotype11";
    sql += "Case_Female_Count_Genotype12 integer,";
    headers[mytablename].push_back("Case_Female_Count_Genotype12");
    casecontrolgenoinsert += ",Case_Female_Count_Genotype12";
    sql += "Case_Female_Count_Genotype22 integer,";
    headers[mytablename].push_back("Case_Female_Count_Genotype22");
    casecontrolgenoinsert += ",Case_Female_Count_Genotype22";

    sql += "Control_Male_Freq_Genotype11 REAL,";
    headers[mytablename].push_back("Control_Male_Freq_Genotype11");
    casecontrolgenoinsert += ",Control_Male_Freq_Genotype11";
    sql += "Control_Male_Freq_Genotype12 REAL,";
    headers[mytablename].push_back("Control_Male_Freq_Genotype12");
    casecontrolgenoinsert += ",Control_Male_Freq_Genotype12";
    sql += "Control_Male_Freq_Genotype22 REAL,";
    headers[mytablename].push_back("Control_Male_Freq_Genotype22");
    casecontrolgenoinsert += ",Control_Male_Freq_Genotype22";
    sql += "Control_Male_Count_Genotype11 integer,";
    headers[mytablename].push_back("Control_Male_Count_Genotype11");
    casecontrolgenoinsert += ",Control_Male_Count_Genotype11";
    sql += "Control_Male_Count_Genotype12 integer,";
    headers[mytablename].push_back("Control_Male_Count_Genotype12");
    casecontrolgenoinsert += ",Control_Male_Count_Genotype12";
    sql += "Control_Male_Count_Genotype22 integer,";
    headers[mytablename].push_back("Control_Male_Count_Genotype22");
    casecontrolgenoinsert += ",Control_Male_Count_Genotype22";
    sql += "Control_Female_Freq_Genotype11 REAL,";
    headers[mytablename].push_back("Control_Female_Freq_Genotype11");
    casecontrolgenoinsert += ",Control_Female_Freq_Genotype11";
    sql += "Control_Female_Freq_Genotype12 REAL,";
    headers[mytablename].push_back("Control_Female_Freq_Genotype12");
    casecontrolgenoinsert += ",Control_Female_Freq_Genotype12";
    sql += "Control_Female_Freq_Genotype22 REAL,";
    headers[mytablename].push_back("Control_Female_Freq_Genotype22");
    casecontrolgenoinsert += ",Control_Female_Freq_Genotype22";
    sql += "Control_Female_Count_Genotype11 integer,";
    headers[mytablename].push_back("Control_Female_Count_Genotype11");
    casecontrolgenoinsert += ",Control_Female_Count_Genotype11";
    sql += "Control_Female_Count_Genotype12 integer,";
    headers[mytablename].push_back("Control_Female_Count_Genotype12");
    casecontrolgenoinsert += ",Control_Female_Count_Genotype12";
    sql += "Control_Female_Count_Genotype22 integer,";
    headers[mytablename].push_back("Control_Female_Count_Genotype22");
    casecontrolgenoinsert += ",Control_Female_Count_Genotype22";
    casecontrolgenoinsert += ") VALUES (NULL";
    sql = sql.replace(sql.size() - 1, 1, ")");

        //Controller::execute_sql(db, sql);
    myQuery.transaction();
    Controller::execute_sql(myQuery, sql);
    myQuery.commit();
    }

    if(options.doGender()){
        mytablename = base + "_gender_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("Gender");
        primary_table[mytablename].push_back("LOCI");
        string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        genderinsert = "INSERT INTO " + mytablename + " (id, fkey";
        for(int i = 0; i < maxalleles; i++){
            sql += "Allele" + getString<int>(i + 1) + " varchar(10),";
            genderinsert += ",Allele" + getString<int>(i + 1);
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Overall_Male_Allele" + getString<int>(i + 1) + "_freq REAL,";
            headers[mytablename].push_back("Overall_Male_Allele" + getString<int>(i + 1) + "_freq");
            genderinsert += ",Overall_Male_Allele" + getString<int>(i + 1) + "_freq";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Overall_Male_Allele" + getString<int>(i + 1) + "_count integer,";
            headers[mytablename].push_back("Overall_Male_Allele" + getString<int>(i + 1) + "_count");
            genderinsert += ",Overall_Male_Allele" + getString<int>(i + 1) + "_count";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Overall_Female_Allele" + getString<int>(i + 1) + "_freq REAL,";
            headers[mytablename].push_back("Overall_Female_Allele" + getString<int>(i + 1) + "_freq");
            genderinsert += ",Overall_Female_Allele" + getString<int>(i + 1) + "_freq";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Overall_Female_Allele" + getString<int>(i + 1) + "_count integer,";
            headers[mytablename].push_back("Overall_Female_Allele" + getString<int>(i + 1) + "_count");
            genderinsert += ",Overall_Female_Allele" + getString<int>(i + 1) + "_count";
        }

        genderinsert += ") VALUES (NULL";
        sql = sql.replace(sql.size() - 1, 1, ")");

        //Controller::execute_sql(db, sql);
        myQuery.transaction();
        Controller::execute_sql(myQuery, sql);
        myQuery.commit();

        mytablename = base + "_genotype_gender_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("Gender Genotype");
        primary_table[mytablename].push_back("LOCI");
        sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        sql += "Genotype11 varchar(20),";
    sql += "Genotype12 varchar(20),";
    sql += "Genotype22 varchar(20),";
        gendergenoinsert = "INSERT INTO " + mytablename + " (id, fkey, Genotype11, Genotype12, Genotype22";
    sql += "Overall_Male_Freq_Genotype11 REAL,";
    headers[mytablename].push_back("Overall_Male_Freq_Genotype11");
    gendergenoinsert += ",Overall_Male_Freq_Genotype11";
    sql += "Overall_Male_Freq_Genotype12 REAL,";
    headers[mytablename].push_back("Overall_Male_Freq_Genotype12");
    gendergenoinsert += ",Overall_Male_Freq_Genotype12";
    sql += "Overall_Male_Freq_Genotype22 REAL,";
    headers[mytablename].push_back("Overall_Male_Freq_Genotype22");
    gendergenoinsert += ",Overall_Male_Freq_Genotype22";
    sql += "Overall_Male_Count_Genotype11 integer,";
    headers[mytablename].push_back("Overall_Male_Count_Genotype11");
    gendergenoinsert += ",Overall_Male_Count_Genotype11";
    sql += "Overall_Male_Count_Genotype12 integer,";
    headers[mytablename].push_back("Overall_Male_Count_Genotype12");
    gendergenoinsert += ",Overall_Male_Count_Genotype12";
    sql += "Overall_Male_Count_Genotype22 integer,";
    headers[mytablename].push_back("Overall_Male_Count_Genotype22");
    gendergenoinsert += ",Overall_Male_Count_Genotype22";
    sql += "Overall_Female_Freq_Genotype11 REAL,";
    headers[mytablename].push_back("Overall_Female_Freq_Genotype11");
    gendergenoinsert += ",Overall_Female_Freq_Genotype11";
    sql += "Overall_Female_Freq_Genotype12 REAL,";
    headers[mytablename].push_back("Overall_Female_Freq_Genotype12");
    gendergenoinsert += ",Overall_Female_Freq_Genotype12";
    sql += "Overall_Female_Freq_Genotype22 REAL,";
    headers[mytablename].push_back("Overall_Female_Freq_Genotype22");
    gendergenoinsert += ",Overall_Female_Freq_Genotype22";
    sql += "Overall_Female_Count_Genotype11 integer,";
    headers[mytablename].push_back("Overall_Female_Count_Genotype11");
    gendergenoinsert += ",Overall_Female_Count_Genotype11";
    sql += "Overall_Female_Count_Genotype12 integer,";
    headers[mytablename].push_back("Overall_Female_Count_Genotype12");
    gendergenoinsert += ",Overall_Female_Count_Genotype12";
    sql += "Overall_Female_Count_Genotype22 integer,";
    headers[mytablename].push_back("Overall_Female_Count_Genotype22");
    gendergenoinsert += ",Overall_Female_Count_Genotype22";
    sql = sql.replace(sql.size() - 1, 1, ")");

    gendergenoinsert += ") VALUES (NULL";
    //Controller::execute_sql(db, sql);
    myQuery.transaction();
    Controller::execute_sql(myQuery, sql);
    myQuery.commit();
    }

    if(options.doParental()){
        mytablename = base + "_parental_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("Parental");
        primary_table[mytablename].push_back("LOCI");
        string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        parentalinsert = "INSERT INTO " + mytablename + " (id, fkey";
        for(int i = 0; i < maxalleles; i++){
            sql += "Allele" + getString<int>(i + 1) + " varchar(10),";
            parentalinsert += ",Allele" + getString<int>(i + 1);
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Parent_Male_Allele" + getString<int>(i + 1) + "_freq REAL,";
            headers[mytablename].push_back("Parent_Male_Allele" + getString<int>(i + 1) + "_freq");
            parentalinsert += ",Parent_Male_Allele" + getString<int>(i + 1) + "_freq";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Parent_Male_Allele" + getString<int>(i + 1) + "_count integer,";
            headers[mytablename].push_back("Parent_Male_Allele" + getString<int>(i + 1) + "_count");
            parentalinsert += ",Parent_Male_Allele" + getString<int>(i + 1) + "_count";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Parent_Female_Allele" + getString<int>(i + 1) + "_freq REAL,";
            headers[mytablename].push_back("Parent_Female_Allele" + getString<int>(i + 1) + "_freq");
            parentalinsert += ",Parent_Female_Allele" + getString<int>(i + 1) + "_freq";
        }
        for(int i = 0; i < maxalleles; i++){
            sql += "Parent_Female_Allele" + getString<int>(i + 1) + "_count integer,";
            headers[mytablename].push_back("Parent_Female_Allele" + getString<int>(i + 1) + "_count");
            parentalinsert += ",Parent_Female_Allele" + getString<int>(i + 1) + "_count";
        }

        sql = sql.replace(sql.size() - 1, 1, ")");

        parentalinsert += ") VALUES (NULL";
        //Controller::execute_sql(db, sql);
        myQuery.transaction();
        Controller::execute_sql(myQuery, sql);
        myQuery.commit();

        mytablename = base + "_genotype_parental_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("Parental Genotype");
        primary_table[mytablename].push_back("LOCI");
        sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        sql += "Genotype11 varchar(20),";
    sql += "Genotype12 varchar(20),";
    sql += "Genotype22 varchar(20),";
    parentalgenoinsert = "INSERT INTO " + mytablename + " (id, fkey, Genotype11, Genotype12, Genotype22";
    sql += "Parent_Male_Freq_Genotype11 REAL,";
    headers[mytablename].push_back("Parent_Male_Freq_Genotype11");
    parentalgenoinsert += ",Parent_Male_Freq_Genotype11";
    sql += "Parent_Male_Freq_Genotype12 REAL,";
    headers[mytablename].push_back("Parent_Male_Freq_Genotype12");
    parentalgenoinsert += ",Parent_Male_Freq_Genotype12";
    sql += "Parent_Male_Freq_Genotype22 REAL,";
    headers[mytablename].push_back("Parent_Male_Freq_Genotype22");
    parentalgenoinsert += ",Parent_Male_Freq_Genotype22";
    sql += "Parent_Male_Count_Genotype11 integer,";
    headers[mytablename].push_back("Parent_Male_Count_Genotype11");
    parentalgenoinsert += ",Parent_Male_Count_Genotype11";
    sql += "Parent_Male_Count_Genotype12 integer,";
    headers[mytablename].push_back("Parent_Male_Count_Genotype12");
    parentalgenoinsert += ",Parent_Male_Count_Genotype12";
    sql += "Parent_Male_Count_Genotype22 integer,";
    headers[mytablename].push_back("Parent_Male_Count_Genotype22");
    parentalgenoinsert += ",Parent_Male_Count_Genotype22";
    sql += "Parent_Female_Freq_Genotype11 REAL,";
    headers[mytablename].push_back("Parent_Female_Freq_Genotype11");
    parentalgenoinsert += ",Parent_Female_Freq_Genotype11";
    sql += "Parent_Female_Freq_Genotype12 REAL,";
    headers[mytablename].push_back("Parent_Female_Freq_Genotype12");
    parentalgenoinsert += ",Parent_Female_Freq_Genotype12";
    sql += "Parent_Female_Freq_Genotype22 REAL,";
    headers[mytablename].push_back("Parent_Female_Freq_Genotype22");
    parentalgenoinsert += ",Parent_Female_Freq_Genotype22";
    sql += "Parent_Female_Count_Genotype11 integer,";
    headers[mytablename].push_back("Parent_Female_Count_Genotype11");
    parentalgenoinsert += ",Parent_Female_Count_Genotype11";
    sql += "Parent_Female_Count_Genotype12 integer,";
    headers[mytablename].push_back("Parent_Female_Count_Genotype12");
    parentalgenoinsert += ",Parent_Female_Count_Genotype12";
    sql += "Parent_Female_Count_Genotype22 integer,";
    headers[mytablename].push_back("Parent_Female_Count_Genotype22");
    parentalgenoinsert += ",Parent_Female_Count_Genotype22";
    sql = sql.replace(sql.size() - 1, 1, ")");

    parentalgenoinsert += ") VALUES (NULL";
	//Controller::execute_sql(db, sql);
    myQuery.transaction();
    Controller::execute_sql(myQuery, sql);
    myQuery.commit();

    }
    if (options.doGroupFile()) {
        mytablename = base + "_group_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("Groups");
        primary_table[mytablename].push_back("LOCI");
        sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        groupinsert = "INSERT INTO " + mytablename + " (id, fkey";
        for(int i = 0; i < maxalleles; i++){
            sql += "Allele" + getString<int>(i + 1) + " varchar(10),";
            groupinsert += ",Allele" + getString<int>(i + 1);
        }
        map<string, vector<Methods::Sample*> >::iterator giter;
        map<string, vector<Methods::Sample*> > groups = options.getGroups();
        for (giter = groups.begin(); giter != groups.end(); giter++) {
            string mygroup = giter->first;
            for (int i = 0; i < maxalleles; i++) {
                sql += mygroup + "_Allele" + getString<int>(i + 1) + "_freq REAL,";
                headers[mytablename].push_back(mygroup + "_Allele" + getString<int>(i + 1) + "_freq");
                groupinsert += "," + mygroup + "_Allele" + getString<int>(i + 1) + "_freq";
            }
            for (int i = 0; i < maxalleles; i++) {
                sql += mygroup + "_Allele" + getString<int>(i + 1) + "_count integer,";
                headers[mytablename].push_back(mygroup + "_Allele" + getString<int>(i + 1) + "_count");
                groupinsert += "," + mygroup + "_Allele" + getString<int>(i + 1) + "_count";
            }
        }
        sql = sql.replace(sql.size() - 1, 1, ")");
        groupinsert += ") VALUES (NULL";
        //Controller::execute_sql(db, sql);
        myQuery.transaction();
        Controller::execute_sql(myQuery, sql);
        myQuery.commit();

        mytablename = base + "_genotype_group_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("Groups Genotype");
        primary_table[mytablename].push_back("LOCI");
        sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        sql += "Genotype11 varchar(20),";
        sql += "Genotype12 varchar(20),";
        sql += "Genotype22 varchar(20),";
        groupgenoinsert = "INSERT INTO " + mytablename + " (id, fkey, Genotype11, Genotype12, Genotype22";
        for (giter = groups.begin(); giter != groups.end(); giter++) {
            string mygroup = giter->first;

            sql += mygroup + "_Freq_Genotype11 REAL,";
            headers[mytablename].push_back(mygroup + "_Freq_Genotype11");
            groupgenoinsert += "," + mygroup + "_Freq_Genotype11";
            sql += mygroup + "_Freq_Genotype12 REAL,";
            headers[mytablename].push_back(mygroup + "_Freq_Genotype12");
            groupgenoinsert += "," + mygroup + "_Freq_Genotype12";
            sql += mygroup + "_Freq_Genotype22 REAL,";
            headers[mytablename].push_back(mygroup + "_Freq_Genotype22");
            groupgenoinsert += "," + mygroup + "_Freq_Genotype22";
            sql += mygroup + "_Count_Genotype11 REAL,";
            headers[mytablename].push_back(mygroup + "_Count_Genotype11");
            groupgenoinsert += "," + mygroup + "_Count_Genotype11";
            sql += mygroup + "_Count_Genotype12 REAL,";
            headers[mytablename].push_back(mygroup + "_Count_Genotype12");
            groupgenoinsert += "," + mygroup + "_Count_Genotype12";
            sql += mygroup + "_Count_Genotype22 REAL,";
            headers[mytablename].push_back(mygroup + "_Count_Genotype22");
            groupgenoinsert += "," + mygroup + "_Count_Genotype22";
        }
        sql = sql.replace(sql.size() - 1, 1, ")");
        groupgenoinsert += ") VALUES (NULL";
        myQuery.transaction();
        Controller::execute_sql(myQuery, sql);
        myQuery.commit();
    }
}

    void ProcessAlleleFrequency::run(DataSetObject* ds)
    {
    	data_set = ds;
    	if (options.doGroupFile())
    	{
    		options.readGroups(data_set->get_samples());
    	}
    	processtest();
    	return;
    }
#endif
#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
