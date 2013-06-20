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

#include "ProcessCaConChisq.h"

#include <iostream>

#include <Options.h>
#include <CaConChisq.h>
#include <Helpers.h>
#include <MethodException.h>
#include <MultComparison.h>

/*
#include <StepOptions.h>
#include <DataSet.h>


#include <fstream>
#include <cdflib.h>
#include <General.h>
*/



using Methods::opts;
using Methods::DataSet;
using Methods::CaConChisq;
using Methods::MethodException;
using Methods::MultComparison;
using Methods::Helpers;
using Methods::Sample;

using std::ofstream;
using std::string;
using std::vector;

string ProcessCaConChisq::stepname = ProcessCaConChisq::doRegister("chisquare");

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

	string mcomp = opts::_OUTPREFIX_ + "chisquare_comparisons" + options.getOut() + ".txt";
	if(!overwrite){
		mcomp += "." + getString<int>(order);
	}
	ofstream mult;


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

	int msize = good_markers.size();//data_set->num_loci();

	if(options.doMultCompare()){
		mult.open(mcomp.c_str());
		if(!mult){
			opts::printLog("Error opening " + mcomp + " for writing results! Exiting!\n");
			throw MethodException("");
		}
		   mult << "Chrom"
				  << "\trsID"
				  << "\tProbeID"
				  << "\tbploc";
			if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
				mult << "\t" << data_set->get_locus(0)->getDetailHeaders();
			}

			mult  << "\tCALC"
				  << "\tOriginal_Pval"
				  << "\tGC"
				  << "\tBONF"
				  << "\tHOLM"
				  << "\tSIDAK_SS"
				  << "\tSIDAK_SD"
				  << "\tFDR_BH"
				  << "\tFDR_BY"
				  << endl;
			opts::addHeader(mcomp, "CALC");
			opts::addHeader(mcomp, "Original_Pval");
			opts::addHeader(mcomp, "GC");
			opts::addHeader(mcomp, "BONF");
			opts::addHeader(mcomp, "HOLM");
			opts::addHeader(mcomp, "SIDAK_SS");
			opts::addHeader(mcomp, "SIDAK_SD");
			opts::addHeader(mcomp, "FDR_BH");
			opts::addHeader(mcomp, "FDR_BY");
	}

	for(int i = 0; i < msize; i++){
		if(good_markers[i]->isEnabled()){//data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged()){
			if(good_markers[i]->isMicroSat()){//data_set->get_locus(i)->isMicroSat()){
				pvals << good_markers[i]->toString()//data_set->get_locus(i)->toString()
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
				if(options.doMultCompare()){
					mult << good_markers[i]->toString()//data_set->get_locus(i)->toString()
						<< "\tNA"
						<< "\tNA"
						<< "\tNA"
						<< "\tNA"
						<< "\tNA"
						<< "\tNA"
						<< "\tNA"
						<< "\tNA"
						<< "\tNA"
						<< endl;
				}
				continue;
			}

			pvals << good_markers[i]->toString() << "\t"//data_set->get_locus(i)->toString() << "\t"
				<< chi_arm[i] << "\t"
				<< pval_arm[i] << "\t"
				<< chi_allele[i] << "\t"
				<< pval_allele[i] << "\t"
				<< pval_allele_exact[i] << "\t"
				<< chi_geno[i] << "\t"
				<< pval_geno[i] << "\t"
				<< pval_geno_exact[i] << "\t"
				<< good_markers[i]->getAllele1() << ":" << good_markers[i]->getAllele2() << "\t"
				<< odds_ratio[i] << "\t"
				<< good_markers[i]->getAllele1() << "\t"
				<< ci_l[i] << "\t"
				<< ci_u[i]
				<< endl;

			if(options.doMultCompare()){
				MultComparison mc(options);
				vector<double> chivals;
				vector<int> tcnt;
//				chivals.push_back(chi_arm[i]);
//				chivals.push_back(chi_allele[i]);
//				chivals.push_back(chi_geno[i]);
				//tcnt.push_back(arm_df[i]);
				//tcnt.push_back(allele_df[i]);
				//tcnt.push_back(geno_df[i]);
				mc.calculate(chi_arm, tcnt);
				mult << good_markers[i]->toString() << "\t"
					<< "ARM\t"
					<< pval_arm[i] << "\t"
					<< mc.get_genomic_control(i) << "\t"
					<< mc.get_bonferroni(i) << "\t"
					<< mc.get_holm(i) << "\t"
					<< mc.get_sidak_single_step(i) << "\t"
					<< mc.get_sidak_step_down(i) << "\t"
					<< mc.get_fdr_bh(i) << "\t"
					<< mc.get_fdr_by(i)
					<< endl;
				mc.calculate(chi_allele, tcnt);
				mult << good_markers[i]->toString() << "\t"
					<< "ALLELIC\t"
					<< pval_allele[i] << "\t"
					<< mc.get_genomic_control(i) << "\t"
					<< mc.get_bonferroni(i) << "\t"
					<< mc.get_holm(i) << "\t"
					<< mc.get_sidak_single_step(i) << "\t"
					<< mc.get_sidak_step_down(i) << "\t"
					<< mc.get_fdr_bh(i) << "\t"
					<< mc.get_fdr_by(i)
					<< endl;
				mc.calculate(chi_geno, tcnt);
				mult << good_markers[i]->toString() << "\t"
					<< "GENO\t"
					<< pval_geno[i] << "\t"
					<< mc.get_genomic_control(i) << "\t"
					<< mc.get_bonferroni(i) << "\t"
					<< mc.get_holm(i) << "\t"
					<< mc.get_sidak_single_step(i) << "\t"
					<< mc.get_sidak_step_down(i) << "\t"
					<< mc.get_fdr_bh(i) << "\t"
					<< mc.get_fdr_by(i)
					<< endl;

			}

			if(options.doGroupFile()){
			}
			good_markers[i]->setFlag(false);
		}
	}

	if(pvals.is_open()){
		pvals.close();
	}
}

/*
 * Function: filter
 * Description:
 * Filters markers based on the armitage pvalue
 */
void ProcessCaConChisq::filter(){

	//false if out in P or C
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int msize = good_markers.size();//data_set->num_loci();
		for(int i = 0; i < msize; i++){
			if(good_markers[i]->isEnabled()){//data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged()){
				bool inc = false;
				if(options.doThreshMarkersLow() && Helpers::dLess(pval_arm[i], options.getThreshMarkersLow())){
					good_markers[i]->setEnabled(false);//data_set->get_locus(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && Helpers::dGreater(pval_arm[i], options.getThreshMarkersHigh())){
					good_markers[i]->setEnabled(false);//data_set->get_locus(i)->setEnabled(false);
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
	good_markers = Helpers::findValidMarkers(data_set->get_markers(), &options);

	//resize(good_markers.size());//ds->num_loci());
	CaConChisq chisq(data_set);
	chisq.setOptions(options);

	if(options.doGroupFile()){
		options.readGroups(ds->get_samples());
	}

	string fnamesv = opts::_OUTPREFIX_ + "chisquare_synthview" + options.getOut() + ".txt";
	if(!overwrite){
		fnamesv += "." + getString<int>(order);
	}
	ofstream svout;
	if(options.doOutputSynthView()){
		svout.open(fnamesv.c_str());
		if(!svout){
			opts::printLog("Error opening " + fnamesv + " for writing results! Exiting!\n");
			throw MethodException("Error opening " + fnamesv + " for writing results! Exiting!\n");
		}
		svout.precision(4);
    	svout << "SNP\tChromosome\tLocation";
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

			if(options.doOutputSynthView()){
				svout << "\t" << group << ":ARM:pval" << "\t" << group << ":AL:pval" << "\t" << group << ":ALE:pval";
			}
		}
		if(options.doOutputSynthView()){
			svout << endl;
		}
		gpvals << endl;


//		int prev_base = 0;
//		int prev_chrom = -1;
		for(int i = 0; i < (int)good_markers.size(); i++){//data_set->num_loci(); i++){
			if(good_markers[i]->isEnabled()){//data_set->get_locus(i)->isEnabled() && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
				if(options.doOutputSynthView()){
					svout << good_markers[i]->getRSID() << "\t" << good_markers[i]->getChrom() << "\t" << good_markers[i]->getBPLOC();
				}

				gpvals << good_markers[i]->toString();//data_set->get_locus(i)->toString();
				if(good_markers[i]->isMicroSat()){//data_set->get_locus(i)->isMicroSat()){
					for(giter = groups.begin(); giter != groups.end(); giter++){
						gpvals << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
					}
					gpvals << endl;
					continue;
				}


				chisq.calculate(good_markers[i]);//i);
				for(giter = groups.begin(); giter != groups.end(); giter++){

					gpvals << "\t" << chisq.getGroupArmitageChi(giter->first) << "\t"
						<< chisq.getGroupArmitagePval(giter->first) << "\t"
						<< chisq.getGroupAllelicChi(giter->first) << "\t"
						<< chisq.getGroupAllelicPval(giter->first) << "\t"
						<< chisq.getGroupAllelicExactPval(giter->first) << "\t"
						<< chisq.getGroupGenotypicChi(giter->first) << "\t"
						<< chisq.getGroupGenotypicPval(giter->first) << "\t"
						<< chisq.getGroupGenotypicExactPval(giter->first) << "\t"
						<< good_markers[i]->getAllele1() << ":" << good_markers[i]->getAllele2() << "\t"
						<< chisq.getGroupOddsRatio(giter->first) << "\t"
						<< good_markers[i]->getAllele1() << "\t"
						<< chisq.getGroupConfIntervalLower(giter->first) << "\t"
						<< chisq.getGroupConfIntervalUpper(giter->first);

					if(options.doOutputSynthView()){
						svout << "\t" << chisq.getGroupArmitagePval(giter->first) << "\t" << chisq.getGroupAllelicPval(giter->first) << "\t" << chisq.getGroupAllelicExactPval(giter->first);
					}
				}
				gpvals << endl;
				if(options.doOutputSynthView()){
					svout << endl;
				}
			}

		}
	}
//#else

	int msize = good_markers.size();//data_set->num_loci();
//	int prev_base = 0;
//	int prev_chrom = -1;
	for(int i = 0; i < msize; i++){
		if(good_markers[i]->isEnabled()){//data_set->get_locus(i)->isEnabled() && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
			if(good_markers[i]->isMicroSat()){//data_set->get_locus(i)->isMicroSat()){
				continue;
			}
			chisq.calculate(good_markers[i]);//i);
			chi_geno[i] = chisq.getGenotypicChi();
			geno_df[i] = chisq.getGenotypicDF();
			pval_geno[i] = chisq.getGenotypicPval();
			pval_geno_exact[i] = chisq.getGenotypicExactPval();
			chi_allele[i] = chisq.getAllelicChi();
			pval_allele[i] = chisq.getAllelicPval();
			allele_df[i] = chisq.getAllelicDF();
			pval_allele_exact[i] = chisq.getAllelicExactPval();
			odds_ratio[i] = chisq.getOddsRatio();
			ci_l[i] = chisq.getConfIntervalLower();
			ci_u[i] = chisq.getConfIntervalUpper();
			chi_arm[i] = chisq.getArmitageChi();
			pval_arm[i] = chisq.getArmitagePval();
			arm_df[i] = chisq.getArmitageDF();
		}
	}
}



