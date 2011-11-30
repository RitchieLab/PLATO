/**********************************************************************************
*                       Case-Control ChiSquare Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Performs a chi-square based Hardy-Weinberg equilibrium test on all markers.
* The Allele Frequency module is called first in order to get up-to-date
* frequency calculations before calculating the Hardy-Weinberg
*
* Files generated:
*
*File: ProcessCaConChisq.cc
**********************************************************************************/


#include <iostream>
#include <fstream>
//#include "Chrom.h"
//#include <ChiSquareAllelic.h>
//#include <ChiSquareArmitage.h>
//#include <FisherExact.h>
//#include "helper.h"
#include <cdflib.h>
#include <General.h>
#include <Options.h>
#include <Helper.h>
#include <MethodException.h>
#include "ProcessCaConChisq.h"

string ProcessCaConChisq::stepname = "chisquare";

/*
 * Function: PrintSummary
 * Description:
 * Outputs results of chisquare step
 * Resets Marker flags
 */
void ProcessCaConChisq::PrintSummary(){
	string pfname = opts::_OUTPREFIX_ + "chisquare" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		pfname += "." + getString<int>(order);
	}
	ofstream pvals (pfname.c_str());
	if(!pvals){
		opts::printLog("Error opening " + pfname + " for writing results! Exiting!\n");
		throw MethodException("");
	}
	opts::addFile("Marker", stepname, pfname);
	pvals.precision(4);



    pvals << "Chrom"
		  << "\trsID"
		  << "\tProbeID"
		  << "\tbploc";
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		pvals << "\t" << data_set->get_locus(0)->getDetailHeaders();
	}
	pvals  << "\tArmitage_Chi"
		  << "\tArmitage_Pval"
		  << "\tAllelic_Chi"
		  << "\tAllelic_Pval"
		  << "\tAllelic_Exact_Pval"
		  << "\tGenotypic_Chi"
		  << "\tGenotypic_Pval"
		  << "\tGenotypic_Exact_Pval"
		  << "\tA1:A2"
		  << "\tOdds_Ratio"
		  << "\tOR_Allele"
		  << "\tL" << (options.getCI()*100)
		  << "\tU" << (options.getCI()*100)
		  << endl;
	opts::addHeader(pfname, "Armitage_Chi");
	opts::addHeader(pfname, "Armitage_Pval");
	opts::addHeader(pfname, "Allelic_Chi");
	opts::addHeader(pfname, "Allelic_Pval");
	opts::addHeader(pfname, "Allelic_Exact_Pval");
	opts::addHeader(pfname, "Genotypic_Chi");
	opts::addHeader(pfname, "Genotypic_Pval");
	opts::addHeader(pfname, "Genotypic_Exact_Pval");
	opts::addHeader(pfname, "A1:A2");
	opts::addHeader(pfname, "Odds_Ratio");
	opts::addHeader(pfname, "OR_Allele");
	opts::addHeader(pfname, "L" + getString<double>(options.getCI()*100));
	opts::addHeader(pfname, "U" + getString<double>(options.getCI()*100));

	int msize = data_set->num_loci();


	for(int i = 0; i < msize; i++){
		if(data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged()){
			if(data_set->get_locus(i)->isMicroSat()){
				pvals << data_set->get_locus(i)->toString()
					<< "\tNA"
					<< "\tNA"
					<< "\tNA"
					<< "\tNA"
					<< "\tNA"
					<< "\tNA"
					<< "\tNA"
					<< "\tNA"
					<< "\tNA"
					<< "\tNA"
					<< "\tNA" << endl;
				continue;
			}

			pvals << data_set->get_locus(i)->toString() << "\t"
				<< chi_arm[i] << "\t"
				<< pval_arm[i] << "\t"
				<< chi_allele[i] << "\t"
				<< pval_allele[i] << "\t"
				<< pval_allele_exact[i] << "\t"
				<< chi_geno[i] << "\t"
				<< pval_geno[i] << "\t"
				<< pval_geno_exact[i] << "\t"
				<< data_set->get_locus(i)->getAllele1() << ":" << data_set->get_locus(i)->getAllele2() << "\t"
				<< odds_ratio[i] << "\t"
				<< data_set->get_locus(i)->getAllele1() << "\t"
				<< ci_l[i] << "\t"
				<< ci_u[i]
				<< endl;

			if(options.doGroupFile()){
			}
			data_set->get_locus(i)->setFlag(false);
		}
	}

	if(pvals.is_open()){
		pvals.close();
	}
}

/*
 * Function: FilterSummary
 * Description:
 * Outputs remaining marker count
 */
void ProcessCaConChisq::FilterSummary(){
    //string pfname = opts::_OUTPREFIX_ + "post_casecontrol_chisquare_filter_summary_" + getStringInt(order) + ".txt";
	//ofstream myoutput (pfname.c_str());
	//if(!myoutput){
	//	cerr << "Error opening post_casecontrol_chisquare_filter_summary.txt!  Exiting!" << endl;
	//	exit(1);
	//}
    //myoutput.precision(4);


	opts::printLog("Options:\t" + options.toString() + "\n");
    //myoutput << "Threshold:\t" << threshold << endl;
	opts::printLog("Markers Passed:\t" + getString<int>((opts::_MARKERS_WORKING_ - orig_num_markers)) + " (" +
	        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
	        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
    //myoutput << "Markers Passed:\t" << (opts::_MARKERS_WORKING_ - orig_num_markers) << " (" <<
	//        ((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0 <<
	//        "%) of " << opts::_MARKERS_WORKING_ << endl;
	opts::_MARKERS_WORKING_ -= orig_num_markers;
//    myoutput << "Families Passed:\t" << families->getSize() << " (" <<
//	        ((float)families->getSize() / (float)orig_num_families) * 100.0 <<
//	        "%) of " << orig_num_families << endl;
//    myoutput << "Individuals Passed:\t" << families->getNumInds() << " (" <<
//	        ((float) families->getNumInds() / (float) orig_num_individuals) * 100.0 <<
//	        "%) of " << orig_num_individuals << endl;
    //myoutput.close();

	resize(0);
}

/*
 * Function: resize
 * Description:
 * Resizes result vectors to specified count
 */
void ProcessCaConChisq::resize(int i){
	chi_geno.resize(i);
	chi_allele.resize(i);
	pval_geno.resize(i);
	pval_allele.resize(i);
	ci_l.resize(i);
	ci_u.resize(i);
	chi_arm.resize(i);
	pval_arm.resize(i);
	pval_allele_exact.resize(i);
	pval_geno_exact.resize(i);
	odds_ratio.resize(i);
}

/*
 * Function: filter
 * Description:
 * Filters markers based on the armitage pvalue
 */
void ProcessCaConChisq::filter(){

	//false if out in P or C
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int msize = data_set->num_loci();
		for(int i = 0; i < msize; i++){
			if(data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged()){
				bool inc = false;
				if(options.doThreshMarkersLow() && dLess(pval_arm[i], options.getThreshMarkersLow())){
					data_set->get_locus(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && dGreater(pval_arm[i], options.getThreshMarkersHigh())){
					data_set->get_locus(i)->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}
		}
	}
}

/*
 * Function: process
 * Description:
 * Main step processing function
 *
 */
void ProcessCaConChisq::process(DataSet* ds){
	data_set = ds;
	resize(ds->num_loci());
	CaConChisq chisq(data_set);
	chisq.setOptions(options);

	if(options.doGroupFile()){
		options.readGroups(ds->get_samples());
	}

	if(options.doGroupFile()){
		string gpfname = opts::_OUTPREFIX_ + "chisquare_groups" + options.getOut() + ".txt";
		ofstream gpvals;
		gpvals.open(gpfname.c_str());
		if(!gpvals){
			opts::printLog("Error opening " + gpfname + " for writing results! Exiting!\n");
			exit(1);
		}
		opts::addFile("Marker", stepname, gpfname);
		gpvals.precision(4);
		gpvals << "Chrom\trsID\tProbeID\tbploc";
		if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
			gpvals << "\t" << data_set->get_locus(0)->getDetailHeaders();
		}
		map<string, vector<Sample*> > groups = options.getGroups();
		map<string, vector<Sample*> >::iterator giter;
		for(giter = groups.begin(); giter != groups.end(); giter++){
			string group = giter->first;
			gpvals  << "\t" << group << "_Armitage_Chi"
		  		<< "\t" << group << "_Armitage_Pval"
		  		<< "\t" << group << "_Allelic_Chi"
		  		<< "\t" << group << "_Allelic_Pval"
		  		<< "\t" << group << "_Allelic_Exact_Pval"
		  		<< "\t" << group << "_Genotypic_Chi"
		  		<< "\t" << group << "_Genotypic_Pval"
		  		<< "\t" << group << "_Genotypic_Exact_Pval"
		  		<< "\t" << group << "_A1:A2"
		  		<< "\t" << group << "_Odds_Ratio"
		  		<< "\t" << group << "_OR_Allele"
				<< "\t" << group << "_L" << (options.getCI() * 100)
				<< "\t" << group << "_U" << (options.getCI() * 100);
			opts::addHeader(gpfname, group + "_Armitage_Chi");
			opts::addHeader(gpfname, group + "_Armitage_Pval");
			opts::addHeader(gpfname, group + "_Allelic_Chi");
			opts::addHeader(gpfname, group + "_Allelic_Pval");
			opts::addHeader(gpfname, group + "_Allelic_Exact_Pval");
			opts::addHeader(gpfname, group + "_Genotypic_Chi");
			opts::addHeader(gpfname, group + "_Genotypic_Pval");
			opts::addHeader(gpfname, group + "_Genotypic_Exact_Pval");
			opts::addHeader(gpfname, group + "_A1:A2");
			opts::addHeader(gpfname, group + "_Odds_Ratio");
			opts::addHeader(gpfname, group + "_OR_Allele");
			opts::addHeader(gpfname, group + "_L" + getString<double>(options.getCI() * 100));
			opts::addHeader(gpfname, group + "_U" + getString<double>(options.getCI() * 100));
		}
		gpvals << endl;


		int prev_base = 0;
		int prev_chrom = -1;
		for(int i = 0; i < data_set->num_loci(); i++){
			if(data_set->get_locus(i)->isEnabled() && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
				gpvals << data_set->get_locus(i)->toString();
				if(data_set->get_locus(i)->isMicroSat()){
					for(giter = groups.begin(); giter != groups.end(); giter++){
						gpvals << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
					}
					gpvals << endl;
					continue;
				}
				chisq.calculate(i);
				for(giter = groups.begin(); giter != groups.end(); giter++){

					gpvals << "\t" << chisq.getGroupArmitageChi(giter->first) << "\t"
						<< chisq.getGroupArmitagePval(giter->first) << "\t"
						<< chisq.getGroupAllelicChi(giter->first) << "\t"
						<< chisq.getGroupAllelicPval(giter->first) << "\t"
						<< chisq.getGroupAllelicExactPval(giter->first) << "\t"
						<< chisq.getGroupGenotypicChi(giter->first) << "\t"
						<< chisq.getGroupGenotypicPval(giter->first) << "\t"
						<< chisq.getGroupGenotypicExactPval(giter->first) << "\t"
						<< data_set->get_locus(i)->getAllele1() << ":" << data_set->get_locus(i)->getAllele2() << "\t"
						<< chisq.getGroupOddsRatio(giter->first) << "\t"
						<< data_set->get_locus(i)->getAllele1() << "\t"
						<< chisq.getGroupConfIntervalLower(giter->first) << "\t"
						<< chisq.getGroupConfIntervalUpper(giter->first);
				}
				gpvals << endl;
			}

		}
	}

	int msize = data_set->num_loci();
	int prev_base = 0;
	int prev_chrom = -1;
	for(int i = 0; i < msize; i++){
		if(data_set->get_locus(i)->isEnabled() && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
			if(data_set->get_locus(i)->isMicroSat()){
				continue;
			}
			chisq.calculate(i);
			chi_geno[i] = chisq.getGenotypicChi();
			pval_geno[i] = chisq.getGenotypicPval();
			pval_geno_exact[i] = chisq.getGenotypicExactPval();
			chi_allele[i] = chisq.getAllelicChi();
			pval_allele[i] = chisq.getAllelicPval();
			pval_allele_exact[i] = chisq.getAllelicExactPval();
			odds_ratio[i] = chisq.getOddsRatio();
			ci_l[i] = chisq.getConfIntervalLower();
			ci_u[i] = chisq.getConfIntervalUpper();
			chi_arm[i] = chisq.getArmitageChi();
			pval_arm[i] = chisq.getArmitagePval();
		}
	}

}


