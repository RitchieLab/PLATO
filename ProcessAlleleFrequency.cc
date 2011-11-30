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
#ifndef MAC
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
#include <Helper.h>

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
				if ((dEquals(minfreq, 0) || dEquals(majfreq, 1))
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
				if (dLess(minfreq, options.getThreshMarkersLow())
						&& options.doThreshMarkersLow()) {
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}
				if (dGreater(minfreq, options.getThreshMarkersHigh())
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
				if ((dEquals(minfreq, 0) || dEquals(majfreq, 1))
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
				if (dLess(minfreq, options.getThreshMarkersLow())
						&& options.doThreshMarkersLow()) {
					//cout << (*markers)[i]->getRSID() << "\t" << minfreq << endl;
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}
				if (dGreater(minfreq, options.getThreshMarkersHigh())
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
				if ((dEquals(minfreq, 0) || dEquals(majfreq, 1))
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
				if (fLess(mark->getMAF(), options.getThreshMarkersLow())
						&& options.doThreshMarkersLow()) {
					mark->setEnabled(false);
					if (inc == false) {
						orig_num_markers++;
						inc = true;
					}
				}
				if (fGreater(mark->getMAF(), options.getThreshMarkersHigh())
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
void ProcessAlleleFrequency::processtest() {

	map<string, double> group_avg;
	int total_snps = 0;

	string afname = opts::_OUTPREFIX_ + "allele_freq" + options.getOut()
			+ ".txt";//+ getString<int>(order) + ".txt";
	if (!overwrite) {
		afname += "." + getString<int> (order);
	}
	string gfname = opts::_OUTPREFIX_ + "allele_freq_genotype"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite) {
		gfname += "." + getString<int> (order);
	}
	string gafname;
	string ggfname;
	string gafnameavg;
	if (options.doGroupFile()) {
		gafname = opts::_OUTPREFIX_ + "allele_freq_group" + options.getOut()
				+ ".txt";//+ getString<int>(order) + ".txt";
		if (!overwrite) {
			gafname += "." + getString<int> (order);
		}
		ggfname = opts::_OUTPREFIX_ + "allele_freq_genotype_group"
				+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if (!overwrite) {
			ggfname += "." + getString<int> (order);
		}
		gafnameavg = opts::_OUTPREFIX_ + "allele_freq_group_avg"
				+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if (!overwrite) {
			gafnameavg += "." + getString<int> (order);
		}
	}

	string pfname = opts::_OUTPREFIX_ + "allele_freq_parental"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite) {
		pfname += "." + getString<int> (order);
	}
	string pgfname = opts::_OUTPREFIX_ + "allele_freq_parental_genotype"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite) {
		pgfname += "." + getString<int> (order);
	}
	string gendfname = opts::_OUTPREFIX_ + "allele_freq_gender"
			+ options.getOut() + ".txt"; //getString<int>(order) + ".txt";
	if (!overwrite) {
		gendfname += "." + getString<int> (order);
	}
	string gendgfname = opts::_OUTPREFIX_ + "allele_freq_gender_genotype"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite) {
		gendgfname += "." + getString<int> (order);
	}
	string ccfname = opts::_OUTPREFIX_ + "allele_freq_casecontrol"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite) {
		ccfname += "." + getString<int> (order);
	}
	string ccgfname = opts::_OUTPREFIX_ + "allele_freq_casecontrol_genotype"
			+ options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if (!overwrite) {
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
	if (!myoutput) {
		opts::printLog("Error opening: " + afname + ".  Exiting!\n");
		throw MethodException("");
	}
	if (!mygeno) {
		opts::printLog("Error opening: " + gfname + ".  Exiting!\n");
		throw MethodException("");
	}

	ofstream gmyoutput; //group file output
	ofstream gmygeno; //group file output
	ofstream gmyoutputavg; //avg groups
	if (options.doGroupFile()) {
		gmyoutput.open(gafname.c_str());
		gmygeno.open(ggfname.c_str());
		opts::addFile("Marker", stepname, gafname);
		opts::addFile("Marker", stepname, ggfname);
		if (!gmyoutput) {
			opts::printLog("Error opening: " + gafname + ".  Exiting!\n");
			throw MethodException("");
		}
		if (!gmygeno) {
			opts::printLog("Error opening: " + ggfname + ".  Exiting!\n");
			throw MethodException("");
		}
		gmyoutput.precision(4);
		gmygeno.precision(4);

		gmyoutputavg.open(gafnameavg.c_str());
		opts::addFile("Batch", stepname, gafnameavg);
		if(!gmyoutputavg){
			throw MethodException("Error opening: " + gafnameavg + ".\n");
		}
		gmyoutputavg << "Batch\tMAF\tN\n";
		opts::addHeader(gafnameavg, "MAF");
		opts::addHeader(gafnameavg, "N");
	}

	int msize = data_set->num_loci();

	int maxalleles = 0;
	for (int i = 0; i < msize; i++) {
		if (data_set->get_locus(i)->isEnabled()) {
			if (data_set->get_locus(i)->getNumAlleles() > maxalleles) {
				maxalleles = data_set->get_locus(i)->getNumAlleles();
			}
		}
	}

	myoutput.precision(4);
	mygeno.precision(4);
	if (data_set->get_locus(0)->getDetailHeaders().size() > 0) {
		myoutput << "Chrom\trsID\tProbeID\tbploc\t"
				<< data_set->get_locus(0)->getDetailHeaders();
		mygeno << "Chrom\trsID\tProbeID\tbploc\t"
				<< data_set->get_locus(0)->getDetailHeaders()
				<< "\tGenotype11\tGenotype12\tGenotype22\tOverall_Freq_Genotype11\tOverall_Freq_Genotype12\tOverall_Freq_Genotype22\tOverall_Count_Genotype11\tOverall_Count_Genotype12\tOverall_Count_Genotype22\tCase_Freq_Genotype11\tCase_Freq_Genotype12\tCase_Freq_Genotype22\tCase_Count_Genotype11\tCase_Count_Genotype12\tCase_Count_Genotype22\tControl_Freq_Genotype11\tControl_Freq_Genotype12\tControl_Freq_Genotype22\tControl_Count_Genotype11\tControl_Count_Genotype12\tControl_Count_Genotype22";
		if (options.doGroupFile()) {
			gmyoutput << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders();
			gmygeno << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders()
					<< "\tGenotype11\tGenotype12\tGenotype22";
		}
	} else {
		myoutput << "Chrom\trsID\tProbeID\tbploc";
		mygeno
				<< "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22\tOverall_Freq_Genotype11\tOverall_Freq_Genotype12\tOverall_Freq_Genotype22\tOverall_Count_Genotype11\tOverall_Count_Genotype12\tOverall_Count_Genotype22\tCase_Freq_Genotype11\tCase_Freq_Genotype12\tCase_Freq_Genotype22\tCase_Count_Genotype11\tCase_Count_Genotype12\tCase_Count_Genotype22\tControl_Freq_Genotype11\tControl_Freq_Genotype12\tControl_Freq_Genotype22\tControl_Count_Genotype11\tControl_Count_Genotype12\tControl_Count_Genotype22";
		if (options.doGroupFile()) {
			gmyoutput << "Chrom\trsID\tProbeID\tbploc";
			gmygeno
					<< "Chrom\trsID\tProbeID\tbploc\tGenotype11\tGenotype12\tGenotype22";
		}
	}
	if (options.doGroupFile()) {
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

	for (int i = 0; i < maxalleles; i++) {
		myoutput << "\tAllele" << (i + 1);
		opts::addHeader(afname, "Allele" + getString<int> (i + 1));
		if (options.doGroupFile()) {
			gmyoutput << "\tAllele" << (i + 1);
			opts::addHeader(gafname, "Allele" + getString<int> (i + 1));
		}
	}
	for (int i = 0; i < maxalleles; i++) {
		myoutput << "\t" << "Overall_Allele" << (i + 1) << "_freq";
		opts::addHeader(afname, "Overall_Allele" + getString<int> (i + 1)
				+ "_freq");
	}
	for (int i = 0; i < maxalleles; i++) {
		myoutput << "\t" << "Overall_Allele" << (i + 1) << "_count";
		opts::addHeader(afname, "Overall_Allele" + getString<int> (i + 1)
				+ "_count");
	}
	for (int i = 0; i < maxalleles; i++) {
		myoutput << "\t" << "Case_Allele" << (i + 1) << "_freq";
		opts::addHeader(afname, "Case_Allele" + getString<int> (i + 1)
				+ "_freq");
	}
	for (int i = 0; i < maxalleles; i++) {
		myoutput << "\t" << "Case_Allele" << (i + 1) << "_count";
		opts::addHeader(afname, "Case_Allele" + getString<int> (i + 1)
				+ "_count");
	}
	for (int i = 0; i < maxalleles; i++) {
		myoutput << "\t" << "Control_Allele" << (i + 1) << "_freq";
		opts::addHeader(afname, "Control_Allele" + getString<int> (i + 1)
				+ "_freq");
	}
	for (int i = 0; i < maxalleles; i++) {
		myoutput << "\t" << "Control_Allele" << (i + 1) << "_count";
		opts::addHeader(afname, "Control_Allele" + getString<int> (i + 1)
				+ "_count");
	}

	if (options.doGroupFile()) {
		map<string, vector<Sample*> >::iterator giter;
		map<string, vector<Sample*> > groups = options.getGroups();
		for (giter = groups.begin(); giter != groups.end(); giter++) {
			string mygroup = giter->first;
			for (int i = 0; i < maxalleles; i++) {
				gmyoutput << "\t" << mygroup << "_Allele" << (i + 1) << "_freq";
				opts::addHeader(gafname, mygroup + "_Allele" + getString<int> (
						i + 1) + "_freq");
			}
			for (int i = 0; i < maxalleles; i++) {
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

	if (options.doParental()) {
		paren.open(pfname.c_str(), ios::out);
		pareng.open(pgfname.c_str(), ios::out);
		if (!paren) {
			opts::printLog("Error opening: " + pfname + ". Exiting!\n");
			throw MethodException("");
		}
		if (!pareng) {
			opts::printLog("Error opening: " + pgfname + ". Exiting!\n");
			throw MethodException("");
		}
		opts::addFile("Marker", stepname, pfname);
		opts::addFile("Marker", stepname, pgfname);
		paren.precision(4);
		pareng.precision(4);
		if (data_set->get_locus(0)->getDetailHeaders().size() > 0) {
			paren << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders();
			pareng << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders()
					<< "\tGenotype11\tGenotype12\tGenotype22\tParent_Male_Freq_Genotype11\tParent_Male_Freq_Genotype12\tParent_Male_Freq_Genotype22\tParent_Male_Count_Genotype11\tParent_Male_Count_Genotype12\tParent_Male_Count_Genotype22\tParent_Female_Freq_Genotype11\tParent_Female_Freq_Genotype12\tParent_Female_Freq_Genotype22\tParent_Female_Count_Genotype11\tParent_Female_Count_Genotype12\tParent_Female_Count_Genotype22\n";
		} else {
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

		for (int i = 0; i < maxalleles; i++) {
			paren << "\t" << "Allele" << (i + 1);
		}
		for (int i = 0; i < maxalleles; i++) {
			paren << "\t" << "Parent_Male_Allele" << (i + 1) << "_freq";
			opts::addHeader(pfname, "Parent_Male_Allele" + getString<int> (i
					+ 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++) {
			paren << "\t" << "Parent_Male_Allele" << (i + 1) << "_count";
			opts::addHeader(pfname, "Parent_Male_Allele" + getString<int> (i
					+ 1) + "_count");
		}
		for (int i = 0; i < maxalleles; i++) {
			paren << "\t" << "Parent_Female_Allele" << (i + 1) << "_freq";
			opts::addHeader(pfname, "Parent_Female_Allele" + getString<int> (i
					+ 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++) {
			paren << "\t" << "Parent_Female_Allele" << (i + 1) << "_count";
			opts::addHeader(pfname, "Parent_Female_Allele" + getString<int> (i
					+ 1) + "_count");
		}
		paren << endl;
	}
	if (options.doGender()) {
		gend.open(gendfname.c_str(), ios::out);
		gendg.open(gendgfname.c_str(), ios::out);
		opts::addFile("Marker", stepname, gendfname);
		opts::addFile("Marker", stepname, gendgfname);
		if (!gend) {
			opts::printLog("Error opening: " + gendfname + ". Exiting!\n");
			throw MethodException("");
		}
		if (!gendg) {
			opts::printLog("Error opening: " + gendgfname + ". Exiting!\n");
			throw MethodException("");
		}
		gend.precision(4);
		gendg.precision(4);
		if (data_set->get_locus(0)->getDetailHeaders().size() > 0) {
			gend << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders();
			gendg << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders()
					<< "\tGenotype11\tGenotype12\tGenotype22\tOverall_Male_Freq_Genotype11\tOverall_Male_Freq_Genotype12\tOverall_Male_Freq_Genotype22\tOverall_Male_Count_Genotype11\tOverall_Male_Count_Genotype12\tOverall_Male_Count_Genotype22\tOverall_Female_Freq_Genotype11\tOverall_Female_Freq_Genotype12\tOverall_Female_Freq_Genotype22\tOverall_Female_Count_Genotype11\tOverall_Female_Count_Genotype12\tOverall_Female_Count_Genotype22\n";
		} else {
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

		for (int i = 0; i < maxalleles; i++) {
			gend << "\t" << "Allele" << (i + 1);
		}
		for (int i = 0; i < maxalleles; i++) {
			gend << "\t" << "Overall_Male_Allele" << (i + 1) << "_freq";
			opts::addHeader(gendfname, "Overall_Male_Allele" + getString<int> (
					i + 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++) {
			gend << "\t" << "Overall_Male_Allele" << (i + 1) << "_count";
			opts::addHeader(gendfname, "Overall_Male_Allele" + getString<int> (
					i + 1) + "_count");
		}
		for (int i = 0; i < maxalleles; i++) {
			gend << "\t" << "Overall_Female_Allele" << (i + 1) << "_freq";
			opts::addHeader(gendfname, "Overall_Female_Allele"
					+ getString<int> (i + 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++) {
			gend << "\t" << "Overall_Female_Allele" << (i + 1) << "_count";
			opts::addHeader(gendfname, "Overall_Female_Allele"
					+ getString<int> (i + 1) + "_count");
		}
		gend << endl;
	}
	if (options.doCaseControl()) {
		cc.open(ccfname.c_str(), ios::out);
		ccg.open(ccgfname.c_str(), ios::out);
		opts::addFile("Marker", stepname, ccfname);
		opts::addFile("Marker", stepname, ccgfname);
		if (!cc) {
			opts::printLog("Error opening: " + ccfname + ". Exiting!\n");
			throw MethodException("");
		}
		if (!ccg) {
			opts::printLog("Error opening: " + ccgfname + ". Exiting!\n");
			throw MethodException("");
		}
		cc.precision(4);
		ccg.precision(4);

		if (data_set->get_locus(0)->getDetailHeaders().size() > 0) {
			cc << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders();
			ccg << "Chrom\trsID\tProbeID\tbploc\t"
					<< data_set->get_locus(0)->getDetailHeaders()
					<< "\tGenotype11\tGenotype12\tGenotype22\tCase_Male_Freq_Genotype11\tCase_Male_Freq_Genotype12\tCase_Male_Freq_Genotype22\tCase_Male_Count_Genotype11\tCase_Male_Count_Genotype12\tCase_Male_Count_Genotype22\tCase_Female_Freq_Genotype11\tCase_Female_Freq_Genotype12\tCase_Female_Freq_Genotype22\tCase_Female_Count_Genotype11\tCase_Female_Count_Genotype12\tCase_Female_Count_Genotype22\tControl_Male_Freq_Genotype11\tControl_Male_Freq_Genotype12\tControl_Male_Freq_Genotype22\tControl_Male_Count_Genotype11\tControl_Male_Count_Genotype12\tControl_Male_Count_Genotype22\tControl_Female_Freq_Genotype11\tControl_Female_Freq_Genotype12\tControl_Female_Freq_Genotype22\tControl_Female_Count_Genotype11\tControl_Female_Count_Genotype12\tControl_Female_Count_Genotype22\n";
		} else {
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

		for (int i = 0; i < maxalleles; i++) {
			cc << "\t" << "Allele" << (i + 1);
		}
		for (int i = 0; i < maxalleles; i++) {
			cc << "\t" << "Case_Male_Allele" << (i + 1) << "_freq";
			opts::addHeader(ccfname, "Case_Male_Allele"
					+ getString<int> (i + 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++) {
			cc << "\t" << "Case_Male_Allele" << (i + 1) << "_count";
			opts::addHeader(ccfname, "Case_Male_Allele"
					+ getString<int> (i + 1) + "_count");
		}
		for (int i = 0; i < maxalleles; i++) {
			cc << "\t" << "Case_Female_Allele" << (i + 1) << "_freq";
			opts::addHeader(ccfname, "Case_Female_Allele" + getString<int> (i
					+ 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++) {
			cc << "\t" << "Case_Female_Allele" << (i + 1) << "_count";
			opts::addHeader(ccfname, "Case_Female_Allele" + getString<int> (i
					+ 1) + "_count");
		}
		for (int i = 0; i < maxalleles; i++) {
			cc << "\t" << "Control_Male_Allele" << (i + 1) << "_freq";
			opts::addHeader(ccfname, "Control_Male_Allele" + getString<int> (i
					+ 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++) {
			cc << "\t" << "Control_Male_Allele" << (i + 1) << "_count";
			opts::addHeader(ccfname, "Control_Male_Allele" + getString<int> (i
					+ 1) + "_count");
		}
		for (int i = 0; i < maxalleles; i++) {
			cc << "\t" << "Control_Female_Allele" << (i + 1) << "_freq";
			opts::addHeader(ccfname, "Control_Female_Allele" + getString<int> (
					i + 1) + "_freq");
		}
		for (int i = 0; i < maxalleles; i++) {
			cc << "\t" << "Control_Female_Allele" << (i + 1) << "_count";
			opts::addHeader(ccfname, "Control_Female_Allele" + getString<int> (
					i + 1) + "_count");
		}
		cc << endl;
	}

	//begin processing
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

	int prev_base = 0;
	int prev_chrom = -1;
	for (int k = 0; k < msize; k++) {
		if (data_set->get_locus(k)->isEnabled() && isValidMarker(data_set->get_locus(k),
				&options, prev_base, prev_chrom)) {

			total_snps++;

			//perform calculations
			af.calcOne(data_set->get_locus(k));
			if (options.doGroupFile()) {
				af.calcOneGroups(data_set->get_locus(k));
			}
			doFilter(data_set->get_locus(k), &af);

			if (data_set->get_locus(k)->isMicroSat()) {
				myoutput << data_set->get_locus(k)->toString();
				if (options.doGroupFile()) {
					gmyoutput << data_set->get_locus(k)->toString();
				}
				if (options.doParental()) {
					paren << data_set->get_locus(k)->toString();
				}
				if (options.doGender()) {
					gend << data_set->get_locus(k)->toString();
				}
				if (options.doCaseControl()) {
					cc << data_set->get_locus(k)->toString();
				}
				int total_o = 0;
				int total_ca = 0;
				int total_con = 0;
				int numalleles = data_set->get_locus(k)->getNumAlleles();
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
					myoutput << "\t" << data_set->get_locus(k)->getAllele(a);
					if (options.doGroupFile()) {
						gmyoutput << "\t" << data_set->get_locus(k)->getAllele(a);
					}
					if (options.doParental()) {
						paren << "\t" << data_set->get_locus(k)->getAllele(a);
					}
					if (options.doGender()) {
						gend << "\t" << data_set->get_locus(k)->getAllele(a);
					}
					if (options.doCaseControl()) {
						cc << "\t" << data_set->get_locus(k)->getAllele(a);
					}
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
						myoutput << "\tNA";
						if (options.doGroupFile()) {
							gmyoutput << "\tNA";
						}
						if (options.doParental()) {
							paren << "\tNA";
						}
						if (options.doGender()) {
							gend << "\tNA";
						}
						if (options.doCaseControl()) {
							cc << "\tNA";
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
					myoutput << "\t" << freq;
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
						myoutput << "\tNA";
					}
				}
				for (int a = 0; a < numalleles; a++) {
					if (useoverall) {
						myoutput << "\t" << af.getMicroCount(a);
					} else {
						myoutput << "\t" << af.getMicroCountP(a);
					}
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
						myoutput << "\tNA";
					}
				}
				//case
				for (int a = 0; a < numalleles; a++) {
					float freq = ((float) af.getMicroCountCa(a)
							/ (float) total_ca);
					myoutput << "\t" << freq;
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
						myoutput << "\tNA";
					}
				}
				for (int a = 0; a < numalleles; a++) {
					myoutput << "\t" << af.getMicroCountCa(a);
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
						myoutput << "\tNA";
					}
				}
				//control
				for (int a = 0; a < numalleles; a++) {
					float freq = ((float) af.getMicroCountCon(a)
							/ (float) total_con);
					myoutput << "\t" << freq;
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
						myoutput << "\tNA";
					}
				}
				for (int a = 0; a < numalleles; a++) {
					myoutput << "\t" << af.getMicroCountCon(a);
				}
				if (maxalleles > numalleles) {
					for (int b = 0; b < (maxalleles - numalleles); b++) {
						myoutput << "\tNA";
					}
				}
				//groups?
				mygeno << data_set->get_locus(k)->toString();
				if (options.doGroupFile()) {
					gmygeno << data_set->get_locus(k)->toString();
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
							gmyoutput << "\t" << freq;
						}
						if (maxalleles > numalleles) {
							for (int b = 0; b < (maxalleles - numalleles); b++) {
								gmyoutput << "\tNA";
							}
						}
						for (int a = 0; a < numalleles; a++) {
							gmyoutput << "\t" << af.getGroupMicroCount(mygroup, a);
						}
						if (maxalleles > numalleles) {
							for (int b = 0; b < (maxalleles - numalleles); b++) {
								gmyoutput << "\tNA";
							}
						}
						gmygeno << "\tNA\tNA\tNA\tNA\tNA\tNA";
					}
					gmyoutput << endl;
				}
				myoutput << endl;

				for (int l = 0; l < 21; l++) {
					mygeno << "\tNA";
					if (options.doGroupFile()) {
						gmygeno << "\tNA";
					}
				}

				mygeno << endl;
				if (options.doGroupFile()) {
					gmygeno << endl;
				}

				if (options.doParental()) {
					pareng << data_set->get_locus(k)->toString();
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
						paren << "\t" << freq;
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							paren << "\tNA";
						}
					}
					for (int a = 0; a < numalleles; a++) {
						paren << "\t" << af.getMicroCountPM(a);
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							paren << "\tNA";
						}
					}
					//Female
					for (int a = 0; a < numalleles; a++) {
						float freq = ((float) af.getMicroCountPF(a)
								/ (float) total_pf);
						paren << "\t" << freq;
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							paren << "\tNA";
						}
					}
					for (int a = 0; a < numalleles; a++) {
						paren << "\t" << af.getMicroCountPF(a);
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							paren << "\tNA";
						}
					}
					paren << endl;

					for (int l = 0; l < 15; l++) {
						pareng << "\tNA";
					}
					pareng << endl;

				}
				if (options.doGender()) {
					gendg << data_set->get_locus(k)->toString();
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
						gend << "\t" << freq;
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							gend << "\tNA";
						}
					}
					for (int a = 0; a < numalleles; a++) {
						if (useoverall) {
							gend << "\t" << af.getMicroCountM(a);
						} else {
							gend << "\t" << af.getMicroCountPM(a);
						}
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							gend << "\tNA";
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
						gend << "\t" << freq;
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							gend << "\tNA";
						}
					}
					for (int a = 0; a < numalleles; a++) {
						if (useoverall) {
							gend << "\t" << af.getMicroCountF(a);
						} else {
							gend << "\t" << af.getMicroCountPF(a);
						}
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							gend << "\tNA";
						}
					}
					gend << endl;

					for (int l = 0; l < 15; l++) {
						gendg << "\tNA";
					}
					gendg << endl;

				}
				if (options.doCaseControl()) {
					ccg << data_set->get_locus(k)->toString();
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
						cc << "\t" << freq;
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							cc << "\tNA";
						}
					}
					for (int a = 0; a < numalleles; a++) {
						cc << "\t" << af.getMicroCountCaM(a);
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							cc << "\tNA";
						}
					}
					//Case Female
					for (int a = 0; a < numalleles; a++) {
						float freq = ((float) af.getMicroCountCaF(a)
								/ (float) total_caf);
						cc << "\t" << freq;
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							cc << "\tNA";
						}
					}
					for (int a = 0; a < numalleles; a++) {
						cc << "\t" << af.getMicroCountCaF(a);
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							cc << "\tNA";
						}
					}
					//Control Male
					for (int a = 0; a < numalleles; a++) {
						float freq = ((float) af.getMicroCountConM(a)
								/ (float) total_conm);
						cc << "\t" << freq;
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							cc << "\tNA";
						}
					}
					for (int a = 0; a < numalleles; a++) {
						cc << "\t" << af.getMicroCountConM(a);
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							cc << "\tNA";
						}
					}
					//Control Female
					for (int a = 0; a < numalleles; a++) {
						float freq = ((float) af.getMicroCountConF(a)
								/ (float) total_conf);
						cc << "\t" << freq;
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							cc << "\tNA";
						}
					}
					for (int a = 0; a < numalleles; a++) {
						cc << "\t" << af.getMicroCountConF(a);
					}
					if (maxalleles > numalleles) {
						for (int b = 0; b < (maxalleles - numalleles); b++) {
							cc << "\tNA";
						}
					}
					cc << endl;

					for (int l = 0; l < 27; l++) {
						ccg << "\tNA";
					}
					ccg << endl;

				}
			} else { //not microsats
				myoutput << data_set->get_locus(k)->toString() << "\t"
						<< data_set->get_locus(k)->getAllele1() << "\t"
						<< data_set->get_locus(k)->getAllele2();
				if (options.doGroupFile()) {
					gmyoutput << data_set->get_locus(k)->toString() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "\t"
							<< data_set->get_locus(k)->getAllele2();
				}
				for (int l = 2; l < maxalleles; l++) {
					myoutput << "\tNA";
					if (options.doGroupFile()) {
						gmyoutput << "\tNA";
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
				myoutput << "\t" << majfreq << "\t" << minfreq;
				for (int l = 2; l < maxalleles; l++) {
					myoutput << "\tNA";
				}
				if (useoverall) {
					myoutput << "\t" << af.getAone_count() << "\t" << af.getAtwo_count();
				} else {
					myoutput << "\t" << af.getAoneP_count() << "\t" << af.getAtwoP_count();
				}
				for (int t = 2; t < maxalleles; t++) {
					myoutput << "\tNA";
				}
				//case
				majfreq = af.getAoneCa_freq();
				minfreq = 1.0f - majfreq;
				myoutput << "\t" << majfreq << "\t" << minfreq;
				for (int l = 2; l < maxalleles; l++) {
					myoutput << "\tNA";
				}
				myoutput << "\t" << af.getAoneCa_count() << "\t" << af.getAtwoCa_count();
				for (int t = 2; t < maxalleles; t++) {
					myoutput << "\tNA";
				}
				//control
				majfreq = af.getAoneCon_freq();
				minfreq = 1.0f - majfreq;
				myoutput << "\t" << majfreq << "\t" << minfreq;
				for (int t = 2; t < maxalleles; t++) {
					myoutput << "\tNA";
				}
				myoutput << "\t" << af.getAoneCon_count() << "\t" << af.getAtwoCon_count();
				for (int t = 2; t < maxalleles; t++) {
					myoutput << "\tNA";
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
						gmyoutput << "\t" << freq << "\t" << freq2;
						for (int b = 2; b < maxalleles; b++) {
							gmyoutput << "\tNA";
						}
						gmyoutput << "\t" << af.getGroupAone_count(mygroup) << "\t"
								<< af.getGroupAtwo_count(mygroup);
						for (int b = 2; b < maxalleles; b++) {
							gmyoutput << "\tNA";
						}

						if(freq < freq2){
							group_avg[mygroup] += freq;
						}
						else{
							group_avg[mygroup] += freq2;
						}
					}
					gmyoutput << endl;
				}

				myoutput << endl;

				mygeno << data_set->get_locus(k)->toString() << "\t"
						<< data_set->get_locus(k)->getAllele1() << "_"
						<< data_set->get_locus(k)->getAllele1() << "\t"
						<< data_set->get_locus(k)->getAllele1() << "_"
						<< data_set->get_locus(k)->getAllele2() << "\t"
						<< data_set->get_locus(k)->getAllele2() << "_"
						<< data_set->get_locus(k)->getAllele2();
				if (options.doGroupFile()) {
					gmygeno << data_set->get_locus(k)->toString() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "_"
							<< data_set->get_locus(k)->getAllele1() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "_"
							<< data_set->get_locus(k)->getAllele2() << "\t"
							<< data_set->get_locus(k)->getAllele2() << "_"
							<< data_set->get_locus(k)->getAllele2();
				}
				//overall
				float freq1 = 0.0f;
				float freq2 = 0.0f;
				float freq3 = 0.0f;
				if (useoverall) {
					freq1 = ((float) af.getAonehomo()) / (af.getPop());
					freq2 = ((float) af.getHet()) / (af.getPop());
					freq3 = ((float) af.getAtwohomo()) / (af.getPop());
					mygeno << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomo() << "\t" << af.getHet()
							<< "\t" << af.getAtwohomo() << "\t";
				} else {
					freq1 = ((float) af.getAonehomoP()) / (af.getPopP());
					freq2 = ((float) af.getHetP()) / (af.getPopP());
					freq3 = ((float) af.getAtwohomoP()) / (af.getPopP());
					mygeno << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoP() << "\t" << af.getHetP()
							<< "\t" << af.getAtwohomoP() << "\t";
				}

				//			mygeno << freq1 << "\t" << freq2 << "\t" << freq3 << "\t" << a1_homo_countP << "\t"
				//				<< a12_countP << "\t"
				//				<< a2_homo_countP << "\t";
				//case overall
				freq1 = ((float) af.getAonehomoCa()) / (af.getPopCa());
				freq2 = ((float) af.getHetCa()) / (af.getPopCa());
				freq3 = ((float) af.getAtwohomoCa()) / (af.getPopCa());
				mygeno << freq1 << "\t" << freq2 << "\t" << freq3 << "\t"
						<< af.getAonehomoCa() << "\t" << af.getHetCa() << "\t"
						<< af.getAtwohomoCa() << "\t";
				//control overall
				freq1 = ((float) af.getAonehomoCon()) / (af.getPopCon());
				freq2 = ((float) af.getHetCon()) / (af.getPopCon());
				freq3 = ((float) af.getAtwohomoCon()) / (af.getPopCon());
				mygeno << freq1 << "\t" << freq2 << "\t" << freq3 << "\t"
						<< af.getAonehomoCon() << "\t" << af.getHetCon() << "\t"
						<< af.getAtwohomoCon();

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
						gmygeno << "\t" << freq1 << "\t" << freq2 << "\t"
								<< freq3 << "\t" << af.getGroupAonehomo(mygroup)
								<< "\t" << af.getGroupHet(mygroup) << "\t"
								<< af.getGroupAtwohomo(mygroup);
					}
					gmygeno << endl;
				}
				mygeno << endl;

				if (options.doParental()) {
					paren << data_set->get_locus(k)->toString() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "\t"
							<< data_set->get_locus(k)->getAllele2();
					for (int l = 2; l < maxalleles; l++) {
						paren << "\tNA";
					}
					//parent male
					float majfreq = af.getAonePM_freq();
					float minfreq = 1.0f - majfreq;
					if(af.getAonePM_count() == 0 && af.getAtwoPM_count() == 0){
						majfreq = 0;
						minfreq = 0;
					}
					paren << "\t" << majfreq << "\t" << minfreq;
					for (int t = 2; t < maxalleles; t++) {
						paren << "\tNA";
					}
					paren << "\t" << af.getAonePM_count() << "\t" << af.getAtwoPM_count();
					for (int t = 2; t < maxalleles; t++) {
						paren << "\tNA";
					}
					//parent female
					majfreq = af.getAonePF_freq();
					minfreq = 1.0f - majfreq;
					paren << "\t" << majfreq << "\t" << minfreq;
					for (int t = 2; t < maxalleles; t++) {
						paren << "\tNA";
					}
					paren << "\t" << af.getAonePF_count() << "\t" << af.getAtwoPF_count();
					for (int t = 2; t < maxalleles; t++) {
						paren << "\tNA";
					}
					paren << endl;

					pareng << data_set->get_locus(k)->toString() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "_"
							<< data_set->get_locus(k)->getAllele1() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "_"
							<< data_set->get_locus(k)->getAllele2() << "\t"
							<< data_set->get_locus(k)->getAllele2() << "_"
							<< data_set->get_locus(k)->getAllele2();
					float freq1 = ((float) af.getAonehomoPM()) / (af.getPopPM());
					float freq2 = ((float) af.getHetPM()) / (af.getPopPM());
					float freq3 = ((float) af.getAtwohomoPM()) / (af.getPopPM());
					pareng << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoPM() << "\t" << af.getHetPM()
							<< "\t" << af.getAtwohomoPM();
					freq1 = ((float) af.getAonehomoPF()) / (af.getPopPF());
					freq2 = ((float) af.getHetPF()) / (af.getPopPF());
					freq3 = ((float) af.getAtwohomoPF()) / (af.getPopPF());
					pareng << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoPF() << "\t" << af.getHetPF()
							<< "\t" << af.getAtwohomoPF();
					pareng << endl;
				}
				if (options.doGender()) {
					//overall male
					gend << data_set->get_locus(k)->toString() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "\t"
							<< data_set->get_locus(k)->getAllele2();
					for (int l = 2; l < maxalleles; l++) {
						gend << "\tNA";
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
					gend << "\t" << majfreq << "\t" << minfreq;
					for (int t = 2; t < maxalleles; t++) {
						gend << "\tNA";
					}
					if (useoverall) {
						gend << "\t" << af.getAoneM_count() << "\t" << af.getAtwoM_count();
					} else {
						gend << "\t" << af.getAonePM_count() << "\t" << af.getAtwoPM_count();
					}
					for (int t = 2; t < maxalleles; t++) {
						gend << "\tNA";
					}
					//overall female
					if (useoverall) {
						majfreq = af.getAoneF_freq();
					} else {
						majfreq = af.getAonePF_freq();
					}
					minfreq = 1.0f - majfreq;
					gend << "\t" << majfreq << "\t" << minfreq;
					for (int t = 2; t < maxalleles; t++) {
						gend << "\tNA";
					}
					if (useoverall) {
						gend << "\t" << af.getAoneF_count() << "\t" << af.getAtwoF_count();
					} else {
						gend << "\t" << af.getAonePF_count() << "\t" << af.getAtwoPF_count();
					}
					for (int t = 2; t < maxalleles; t++) {
						gend << "\tNA";
					}
					gend << endl;

					gendg << data_set->get_locus(k)->toString() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "_"
							<< data_set->get_locus(k)->getAllele1() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "_"
							<< data_set->get_locus(k)->getAllele2() << "\t"
							<< data_set->get_locus(k)->getAllele2() << "_"
							<< data_set->get_locus(k)->getAllele2();
					if (useoverall) {
						float freq1 =
								((float) af.getAonehomoM()) / (af.getPopM());
						float freq2 = ((float) af.getHetM()) / (af.getPopM());
						float freq3 =
								((float) af.getAtwohomoM()) / (af.getPopM());
						gendg << "\t" << freq1 << "\t" << freq2 << "\t"
								<< freq3 << "\t" << af.getAonehomoM() << "\t"
								<< af.getHetM() << "\t" << af.getAtwohomoM();
						freq1 = ((float) af.getAonehomoF()) / (af.getPopF());
						freq2 = ((float) af.getHetF()) / (af.getPopF());
						freq3 = ((float) af.getAtwohomoF()) / (af.getPopF());
						gendg << "\t" << freq1 << "\t" << freq2 << "\t"
								<< freq3 << "\t" << af.getAonehomoF() << "\t"
								<< af.getHetF() << "\t" << af.getAtwohomoF();
					} else {
						float freq1 = ((float) af.getAonehomoPM())
								/ (af.getPopPM());
						float freq2 = ((float) af.getHetPM()) / (af.getPopPM());
						float freq3 = ((float) af.getAtwohomoPM())
								/ (af.getPopPM());
						gendg << "\t" << freq1 << "\t" << freq2 << "\t"
								<< freq3 << "\t" << af.getAonehomoPM() << "\t"
								<< af.getHetPM() << "\t" << af.getAtwohomoPM();
						freq1 = ((float) af.getAonehomoPF()) / (af.getPopPF());
						freq2 = ((float) af.getHetPF()) / (af.getPopPF());
						freq3 = ((float) af.getAtwohomoPF()) / (af.getPopPF());
						gendg << "\t" << freq1 << "\t" << freq2 << "\t"
								<< freq3 << "\t" << af.getAonehomoPF() << "\t"
								<< af.getHetPF() << "\t" << af.getAtwohomoPF();
					}
					gendg << endl;

				}
				if (options.doCaseControl()) {
					//case male
					cc << data_set->get_locus(k)->toString() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "\t"
							<< data_set->get_locus(k)->getAllele2();
					for (int l = 2; l < maxalleles; l++) {
						cc << "\tNA";
					}
					float majfreq = af.getAoneCaM_freq();
					float minfreq = 1.0f - majfreq;
					if(af.getAoneCaM_count() == 0 && af.getAtwoCaM_count() == 0){
						majfreq = 0;
						minfreq = 0;
					}
					cc << "\t" << majfreq << "\t" << minfreq;
					for (int t = 2; t < maxalleles; t++) {
						cc << "\tNA";
					}
					cc << "\t" << af.getAoneCaM_count() << "\t" << af.getAtwoCaM_count();
					for (int t = 2; t < maxalleles; t++) {
						cc << "\tNA";
					}
					//case female
					majfreq = af.getAoneCaF_freq();
					minfreq = 1.0f - majfreq;
					cc << "\t" << majfreq << "\t" << minfreq;
					for (int t = 2; t < maxalleles; t++) {
						cc << "\tNA";
					}
					cc << "\t" << af.getAoneCaF_count() << "\t" << af.getAtwoCaF_count();
					for (int t = 2; t < maxalleles; t++) {
						cc << "\tNA";
					}
					//control male
					majfreq = af.getAoneConM_freq();
					minfreq = 1.0f - majfreq;
					if(af.getAoneConM_count() == 0 && af.getAtwoConM_count() == 0){
						majfreq = 0;
						minfreq = 0;
					}
					cc << "\t" << majfreq << "\t" << minfreq;
					for (int t = 2; t < maxalleles; t++) {
						cc << "\tNA";
					}
					cc << "\t" << af.getAoneConM_count() << "\t" << af.getAtwoConM_count();
					for (int t = 2; t < maxalleles; t++) {
						cc << "\tNA";
					}
					//control female
					majfreq = af.getAoneConF_freq();
					minfreq = 1.0f - majfreq;
					cc << "\t" << majfreq << "\t" << minfreq;
					for (int t = 2; t < maxalleles; t++) {
						cc << "\tNA";
					}
					cc << "\t" << af.getAoneConF_count() << "\t" << af.getAtwoConF_count();
					for (int t = 2; t < maxalleles; t++) {
						cc << "\tNA";
					}
					cc << endl;

					ccg << data_set->get_locus(k)->toString() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "_"
							<< data_set->get_locus(k)->getAllele1() << "\t"
							<< data_set->get_locus(k)->getAllele1() << "_"
							<< data_set->get_locus(k)->getAllele2() << "\t"
							<< data_set->get_locus(k)->getAllele2() << "_"
							<< data_set->get_locus(k)->getAllele2();
					float freq1 = ((float) af.getAonehomoCaM())
							/ (af.getPopCaM());
					float freq2 = ((float) af.getHetCaM()) / (af.getPopCaM());
					float freq3 = ((float) af.getAtwohomoCaM())
							/ (af.getPopCaM());
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoCaM() << "\t" << af.getHetCaM()
							<< "\t" << af.getAtwohomoCaM();
					freq1 = ((float) af.getAonehomoCaF()) / (af.getPopCaF());
					freq2 = ((float) af.getHetCaF()) / (af.getPopCaF());
					freq3 = ((float) af.getAtwohomoCaF()) / (af.getPopCaF());
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoCaF() << "\t" << af.getHetCaF()
							<< "\t" << af.getAtwohomoCaF();
					freq1 = ((float) af.getAonehomoConM()) / (af.getPopConM());
					freq2 = ((float) af.getHetConM()) / (af.getPopConM());
					freq3 = ((float) af.getAtwohomoConM()) / (af.getPopConM());
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoConM() << "\t"
							<< af.getHetConM() << "\t" << af.getAtwohomoConM();
					freq1 = ((float) af.getAonehomoConF()) / (af.getPopConF());
					freq2 = ((float) af.getHetConF()) / (af.getPopConF());
					freq3 = ((float) af.getAtwohomoConF()) / (af.getPopConF());
					ccg << "\t" << freq1 << "\t" << freq2 << "\t" << freq3
							<< "\t" << af.getAonehomoConF() << "\t"
							<< af.getHetConF() << "\t" << af.getAtwohomoConF();
					ccg << endl;

				}
			}

			//filter Markers
			//doFilter(data_set->get_locus(k), &af);
		}
	}

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

}

/*
 *Function: process
 *Description:
 *Main method to begin the whole process.  Flags samples then diverts work to processtest
 *
 *
 */
void ProcessAlleleFrequency::process(DataSet* ds) {
	data_set = ds;
	if (options.doGroupFile()) {
		options.readGroups(data_set->get_samples());
	}

	processtest();
	return;

}

