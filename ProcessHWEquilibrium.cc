/**********************************************************************************
*                       Hardy-Weinberg Equilibrium Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Performs a chi-square based Hardy-Weinberg equilibrium test on all markers.
* A hybrid approach is used.  If a genotype count is < 5, then it uses the exact calculation.
* Otherwise it uses the classical calculation.
* The Allele Frequency module is called first in order to get up-to-date
* frequency calculations before calculating the Hardy-Weinberg
*
*
*File: HWEquilibrium.cc
**********************************************************************************/


#include <iostream>
#include <fstream>
#include <algorithm>
#include "ProcessHWEquilibrium.h"
#include "Chrom.h"
//#include "ChiSquare.h"
//#include "helper.h"
#include <cdflib.h>
#include <General.h>
#include <Helper.h>

string ProcessHWEquilibrium::stepname = "hwe";

/*
 * Function: PrintSummary
 * Description:
 * Resets marker flags
 */
void ProcessHWEquilibrium::PrintSummary(){
	for(int i = 0; i < data_set->num_loci(); i++){
		data_set->get_locus(i)->setFlag(false);
	}
	return;
}

/*
 * Function: FilterSummary
 * Description:
 * Outputs number of remaining markers
 */
void ProcessHWEquilibrium::FilterSummary(){
	int fsize = data_set->get_families()->size();
	for(int i = 0; i < fsize; i++){
		data_set->get_pedigree(i)->setFlagAFCM(false);
		data_set->get_pedigree(i)->setFlagAFCF(false);
	}


	int msize = data_set->num_inds();
	opts::printLog("Options:\t" + options.toString() + "\n");
    opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
	        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
	        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessHWEquilibrium::filter(){
	return;
}

/*
 * Function: doFilter
 * Description:
 * Filters markers based on Overall pvalue
 */
void ProcessHWEquilibrium::doFilter(Marker* mark, HWEquilibrium* hwe){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		bool inc = false;
		if(options.doThreshMarkersLow() && ((dLess(hwe->getOverall(), options.getThreshMarkersLow()) && hwe->getOverall() != -1))){// || (dLess(hw_C,options.getThreshMarkersLow()) && hw_C != -1))){
			mark->setEnabled(false);
			inc = true;
		}
		if(options.doThreshMarkersHigh() && ((dGreater(hwe->getOverall(), options.getThreshMarkersHigh()) && hwe->getOverall() != -1))){// || (dGreater(hw_C, options.getThreshMarkersHigh()) && hw_C != -1))){
			mark->setEnabled(false);
			inc = true;
		}
		if(inc){
			orig_num_markers++;
		}
	}

}



/*
 * Function: process
 * Description:
 * Performs the HWE process.
 */
void ProcessHWEquilibrium::process(DataSet* ds){
	data_set = ds;

	useoverall = false;
	if(options.doRandomChild() || options.doAll() || options.doAllChildren()){
		useoverall = true;
	}

	HWEquilibrium hwe(data_set);
	hwe.setOptions(options);
	AlleleFrequency af(data_set);
	af.setOptions(options);

	//	af->process_hw(samples, families, markers, marker_map);

	int msize = data_set->num_loci();

	if(options.doHWEPT()){
	    string fname = opts::_OUTPREFIX_ + "hw" + options.getOut() + ".txt";
	    if(!overwrite){
	        fname += "." + getString<int>(order);
	    }
	    ofstream pvals (fname.c_str());
	    if(!pvals){
	        opts::printLog("Error opening " + fname + "!  Exiting!\n");
	        throw MethodException("");
	    }
	    pvals.precision(4);

	    opts::addFile("Marker", stepname, fname);
	    pvals << "Chrom\trsID\tProbeID\tbploc";
	    if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
	        pvals << "\t" << data_set->get_locus(0)->getDetailHeaders();
	    }
	    pvals << "\tPT_pval_case\tPT_pval_control\tP_pval_genotypic\tP_pval_allelic\tP_pval_case\tP_pval_control\n";
	    opts::addHeader(fname, "PT_pval_case");
	    opts::addHeader(fname, "PT_pval_control");
	    opts::addHeader(fname, "P_pval_genotypic");
	    opts::addHeader(fname, "P_pval_allelic");
	    opts::addHeader(fname, "P_pval_case");
	    opts::addHeader(fname, "P_pval_control");

		int prev_base = 0;
		int prev_chrom = -1;
		for(int i = 0; i < msize; i++){
			if(data_set->get_locus(i)->isEnabled() && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
				hwe.calculateHWEPT(data_set->get_locus(i));
				pvals << data_set->get_locus(i)->toString() << "\t" << hwe.getCasePvalHWEPT() << "\t" << hwe.getControlPvalHWEPT() << "\t" << hwe.getGenotypicPval() << "\t" << hwe.getAllelicPval() << "\t" << hwe.getCaseGeneralHWEPT() << "\t" << hwe.getControlGeneralHWEPT() << endl;
			}
		}
		pvals.close();
		return;
	}

	string fname1 = opts::_OUTPREFIX_ + "hw" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	string fname2 = opts::_OUTPREFIX_ + "hw_parental" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	string fname3 = opts::_OUTPREFIX_ + "hw_gender" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname3 += "." + getString<int>(order);
	}
	string fname4 = opts::_OUTPREFIX_ + "hw_casecontrol" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname4 += "." + getString<int>(order);
	}
	ofstream paren;
	ofstream gend;
	ofstream cc;
	ofstream pvals (fname1.c_str());
	if(!pvals){
		opts::printLog("Error opening " + fname1 + "!  Exiting!\n");
		throw MethodException("");
	}
	pvals.precision(4);
	opts::addFile("Marker",stepname,fname1);
	if(options.doParental()){
		paren.open(fname2.c_str());
		if(!paren){
			opts::printLog("Error opening " + fname2 + "! Exiting!\n");
			throw MethodException("");
		}
		paren.precision(4);
		opts::addFile("Marker",stepname,fname2);
	}
	if(options.doGender()){
		gend.open(fname3.c_str());
		if(!gend){
			opts::printLog("Error opening " + fname3 + "! Exiting!\n");
			throw MethodException("");
		}
		gend.precision(4);
		opts::addFile("Marker",stepname,fname3);
	}
	if(options.doCaseControl()){
		cc.open(fname4.c_str());
		if(!cc){
			opts::printLog("Error opening " + fname4 + "! Exiting!\n");
			throw MethodException("");
		}
		cc.precision(4);
		opts::addFile("Marker",stepname,fname4);
	}

    pvals << "Chrom"
		  << "\trsID"
		  << "\tProbeID"
		  << "\tbploc";
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		pvals << "\t" << data_set->get_locus(0)->getDetailHeaders();
	}
	pvals  << "\tGenotype11"
		  << "\tGenotype12"
		  << "\tGenotype22"
		  << "\tOverall_Pvalue"
		  << "\tOverall_Total_Count"
		  << "\tOverall_Obs_Genotype11"
		  << "\tOverall_Exp_Genotype11"
		  << "\tOverall_Obs_Genotype12"
		  << "\tOverall_Exp_Genotype12"
		  << "\tOverall_Obs_Genotype22"
		  << "\tOverall_Exp_Genotype22"
		  << "\tCase_Pvalue"
		  << "\tCase_Total_Count"
		  << "\tCase_Obs_Genotype11"
		  << "\tCase_Exp_Genotype11"
		  << "\tCase_Obs_Genotype12"
		  << "\tCase_Exp_Genotype12"
		  << "\tCase_Obs_Genotype22"
		  << "\tCase_Exp_Genotype22"
		  << "\tControl_Pvalue"
		  << "\tControl_Total_Count"
		  << "\tControl_Obs_Genotype11"
		  << "\tControl_Exp_Genotype11"
		  << "\tControl_Obs_Genotype12"
		  << "\tControl_Exp_Genotype12"
		  << "\tControl_Obs_Genotype22"
		  << "\tControl_Exp_Genotype22"
		  << endl;

	opts::addHeader(fname1, "Overall_Pvalue");
	opts::addHeader(fname1, "Overall_Total_Count");
	opts::addHeader(fname1, "Overall_Obs_Genotype11");
	opts::addHeader(fname1, "Overall_Exp_Genotype11");
	opts::addHeader(fname1, "Overall_Obs_Genotype12");
	opts::addHeader(fname1, "Overall_Exp_Genotype12");
	opts::addHeader(fname1, "Overall_Obs_Genotype22");
	opts::addHeader(fname1, "Overall_Exp_Genotype22");
	opts::addHeader(fname1, "Case_Pvalue");
	opts::addHeader(fname1, "Case_Total_Count");
	opts::addHeader(fname1, "Case_Obs_Genotype11");
	opts::addHeader(fname1, "Case_Exp_Genotype11");
	opts::addHeader(fname1, "Case_Obs_Genotype12");
	opts::addHeader(fname1, "Case_Exp_Genotype12");
	opts::addHeader(fname1, "Case_Obs_Genotype22");
	opts::addHeader(fname1, "Case_Exp_Genotype22");
	opts::addHeader(fname1, "Control_Pvalue");
	opts::addHeader(fname1, "Control_Total_Count");
	opts::addHeader(fname1, "Control_Obs_Genotype11");
	opts::addHeader(fname1, "Control_Exp_Genotype11");
	opts::addHeader(fname1, "Control_Obs_Genotype12");
	opts::addHeader(fname1, "Control_Exp_Genotype12");
	opts::addHeader(fname1, "Control_Obs_Genotype22");
	opts::addHeader(fname1, "Control_Exp_Genotype22");

	if(options.doParental()){
    paren << "Chrom"
		  << "\trsID"
		  << "\tProbeID"
		  << "\tbploc";
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		paren << "\t" << data_set->get_locus(0)->getDetailHeaders();
	}
	paren  << "\tGenotype11"
		  << "\tGenotype12"
		  << "\tGenotype22"
		  << "\tParent_Male_Pvalue"
		  << "\tParent_Male_Total_Count"
		  << "\tParent_Male_Obs_Genotype11"
		  << "\tParent_Male_Exp_Genotype11"
		  << "\tParent_Male_Obs_Genotype12"
		  << "\tParent_Male_Exp_Genotype12"
		  << "\tParent_Male_Obs_Genotype22"
		  << "\tParent_Male_Exp_Genotype22"
		  << "\tParent_Female_Pvalue"
		  << "\tParent_Female_Total_Count"
		  << "\tParent_Female_Obs_Genotype11"
		  << "\tParent_Female_Exp_Genotype11"
		  << "\tParent_Female_Obs_Genotype12"
		  << "\tParent_Female_Exp_Genotype12"
		  << "\tParent_Female_Obs_Genotype22"
		  << "\tParent_Female_Exp_Genotype22"
		  << endl;

	opts::addHeader(fname2, "Parent_Male_Pvalue");
	opts::addHeader(fname2, "Parent_Male_Total_Count");
	opts::addHeader(fname2, "Parent_Male_Obs_Genotype11");
	opts::addHeader(fname2, "Parent_Male_Exp_Genotype11");
	opts::addHeader(fname2, "Parent_Male_Obs_Genotype12");
	opts::addHeader(fname2, "Parent_Male_Exp_Genotype12");
	opts::addHeader(fname2, "Parent_Male_Obs_Genotype22");
	opts::addHeader(fname2, "Parent_Male_Exp_Genotype22");
	opts::addHeader(fname2, "Parent_Female_Pvalue");
	opts::addHeader(fname2, "Parent_Female_Total_Count");
	opts::addHeader(fname2, "Parent_Female_Obs_Genotype11");
	opts::addHeader(fname2, "Parent_Female_Exp_Genotype11");
	opts::addHeader(fname2, "Parent_Female_Obs_Genotype12");
	opts::addHeader(fname2, "Parent_Female_Exp_Genotype12");
	opts::addHeader(fname2, "Parent_Female_Obs_Genotype22");
	opts::addHeader(fname2, "Parent_Female_Exp_Genotype22");
	}
	if(options.doGender()){
    gend << "Chrom"
		  << "\trsID"
		  << "\tProbeID"
		  << "\tbploc";
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		gend << "\t" << data_set->get_locus(0)->getDetailHeaders();
	}
	gend  << "\tGenotype11"
		  << "\tGenotype12"
		  << "\tGenotype22"
		  << "\tOverall_Male_Pvalue"
		  << "\tOverall_Male_Total_Count"
		  << "\tOverall_Male_Obs_Genotype11"
		  << "\tOverall_Male_Exp_Genotype11"
		  << "\tOverall_Male_Obs_Genotype12"
		  << "\tOverall_Male_Exp_Genotype12"
		  << "\tOverall_Male_Obs_Genotype22"
		  << "\tOverall_Male_Exp_Genotype22"
		  << "\tOverall_Female_Pvalue"
		  << "\tOverall_Female_Total_Count"
		  << "\tOverall_Female_Obs_Genotype11"
		  << "\tOverall_Female_Exp_Genotype11"
		  << "\tOverall_Female_Obs_Genotype12"
		  << "\tOverall_Female_Exp_Genotype12"
		  << "\tOverall_Female_Obs_Genotype22"
		  << "\tOverall_Female_Exp_Genotype22"
		  << endl;
	opts::addHeader(fname3, "Overall_Male_Pvalue");
	opts::addHeader(fname3, "Overall_Male_Total_Count");
	opts::addHeader(fname3, "Overall_Male_Obs_Genotype11");
	opts::addHeader(fname3, "Overall_Male_Exp_Genotype11");
	opts::addHeader(fname3, "Overall_Male_Obs_Genotype12");
	opts::addHeader(fname3, "Overall_Male_Exp_Genotype12");
	opts::addHeader(fname3, "Overall_Male_Obs_Genotype22");
	opts::addHeader(fname3, "Overall_Male_Exp_Genotype22");
	opts::addHeader(fname3, "Overall_Female_Pvalue");
	opts::addHeader(fname3, "Overall_Female_Total_Count");
	opts::addHeader(fname3, "Overall_Female_Obs_Genotype11");
	opts::addHeader(fname3, "Overall_Female_Exp_Genotype11");
	opts::addHeader(fname3, "Overall_Female_Obs_Genotype12");
	opts::addHeader(fname3, "Overall_Female_Exp_Genotype12");
	opts::addHeader(fname3, "Overall_Female_Obs_Genotype22");
	opts::addHeader(fname3, "Overall_Female_Exp_Genotype22");
	}
	if(options.doCaseControl()){
    cc << "Chrom"
		  << "\trsID"
		  << "\tProbeID"
		  << "\tbploc";
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		cc << "\t" << data_set->get_locus(0)->getDetailHeaders();
	}
	cc << "\tGenotype11"
		  << "\tGenotype12"
		  << "\tGenotype22"
		  << "\tCase_Male_Pvalue"
		  << "\tCase_Male_Total_Count"
		  << "\tCase_Male_Obs_Genotype11"
		  << "\tCase_Male_Exp_Genotype11"
		  << "\tCase_Male_Obs_Genotype12"
		  << "\tCase_Male_Exp_Genotype12"
		  << "\tCase_Male_Obs_Genotype22"
		  << "\tCase_Male_Exp_Genotype22"
		  << "\tCase_Female_Pvalue"
		  << "\tCase_Female_Total_Count"
		  << "\tCase_Female_Obs_Genotype11"
		  << "\tCase_Female_Exp_Genotype11"
		  << "\tCase_Female_Obs_Genotype12"
		  << "\tCase_Female_Exp_Genotype12"
		  << "\tCase_Female_Obs_Genotype22"
		  << "\tCase_Female_Exp_Genotype22"
		  << "\tControl_Male_Pvalue"
		  << "\tControl_Male_Total_Count"
		  << "\tControl_Male_Obs_Genotype11"
		  << "\tControl_Male_Exp_Genotype11"
		  << "\tControl_Male_Obs_Genotype12"
		  << "\tControl_Male_Exp_Genotype12"
		  << "\tControl_Male_Obs_Genotype22"
		  << "\tControl_Male_Exp_Genotype22"
		  << "\tControl_Female_Pvalue"
		  << "\tControl_Female_Total_Count"
		  << "\tControl_Female_Obs_Genotype11"
		  << "\tControl_Female_Exp_Genotype11"
		  << "\tControl_Female_Obs_Genotype12"
		  << "\tControl_Female_Exp_Genotype12"
		  << "\tControl_Female_Obs_Genotype22"
		  << "\tControl_Female_Exp_Genotype22"
		  << endl;

	opts::addHeader(fname4, "Case_Male_Pvalue");
	opts::addHeader(fname4, "Case_Male_Total_Count");
	opts::addHeader(fname4, "Case_Male_Obs_Genotype11");
	opts::addHeader(fname4, "Case_Male_Exp_Genotype11");
	opts::addHeader(fname4, "Case_Male_Obs_Genotype12");
	opts::addHeader(fname4, "Case_Male_Exp_Genotype12");
	opts::addHeader(fname4, "Case_Male_Obs_Genotype22");
	opts::addHeader(fname4, "Case_Male_Exp_Genotype22");
	opts::addHeader(fname4, "Case_Female_Pvalue");
	opts::addHeader(fname4, "Case_Female_Total_Count");
	opts::addHeader(fname4, "Case_Female_Obs_Genotype11");
	opts::addHeader(fname4, "Case_Female_Exp_Genotype11");
	opts::addHeader(fname4, "Case_Female_Obs_Genotype12");
	opts::addHeader(fname4, "Case_Female_Exp_Genotype12");
	opts::addHeader(fname4, "Case_Female_Obs_Genotype22");
	opts::addHeader(fname4, "Case_Female_Exp_Genotype22");
	opts::addHeader(fname4, "Control_Male_Pvalue");
	opts::addHeader(fname4, "Control_Male_Total_Count");
	opts::addHeader(fname4, "Control_Male_Obs_Genotype11");
	opts::addHeader(fname4, "Control_Male_Exp_Genotype11");
	opts::addHeader(fname4, "Control_Male_Obs_Genotype12");
	opts::addHeader(fname4, "Control_Male_Exp_Genotype12");
	opts::addHeader(fname4, "Control_Male_Obs_Genotype22");
	opts::addHeader(fname4, "Control_Male_Exp_Genotype22");
	opts::addHeader(fname4, "Control_Female_Pvalue");
	opts::addHeader(fname4, "Control_Female_Total_Count");
	opts::addHeader(fname4, "Control_Female_Obs_Genotype11");
	opts::addHeader(fname4, "Control_Female_Exp_Genotype11");
	opts::addHeader(fname4, "Control_Female_Obs_Genotype12");
	opts::addHeader(fname4, "Control_Female_Exp_Genotype12");
	opts::addHeader(fname4, "Control_Female_Obs_Genotype22");
	opts::addHeader(fname4, "Control_Female_Exp_Genotype22");
	}
	int prev_base = 0;
	int prev_chrom = -1;

	for(int i = 0; i < msize; i++){
		if(data_set->get_locus(i)->isEnabled() && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
			hwe.calculate(i);
			af.calculate(i);

			pvals << data_set->get_locus(i)->toString();
			if(data_set->get_locus(i)->isMicroSat()){
				for(int l = 0; l < 27; l++){
					pvals << "\tNA";
				}
			}
			else{
				//overall default = founders
				pvals << "\t" << data_set->get_locus(i)->getAllele1() << "_" << data_set->get_locus(i)->getAllele1() << "\t"
				<< data_set->get_locus(i)->getAllele1() << "_" << data_set->get_locus(i)->getAllele2() << "\t"
				<< data_set->get_locus(i)->getAllele2() << "_" << data_set->get_locus(i)->getAllele2() << "\t"
				<< hwe.getOverall() << "\t";
				if(useoverall){
					pvals << af.getPop() << "\t"
					<< af.getAonehomo() << "\t"
					<< af.getAonehomo_exp() << "\t"
					<< af.getHet() << "\t"
					<< af.getHet_exp() << "\t"
					<< af.getAtwohomo() << "\t"
					<< af.getAtwohomo_exp() << "\t";
				}
				else{
					pvals << af.getPopP() << "\t"
					<< af.getAonehomoP() << "\t"
					<< af.getAonehomoP_exp() << "\t"
					<< af.getHetP() << "\t"
					<< af.getHetP_exp() << "\t"
					<< af.getAtwohomoP() << "\t"
					<< af.getAtwohomoP_exp() << "\t";
				}
				//cases
				pvals << hwe.getCase() << "\t"
				<< af.getPopCa() << "\t"
				<< af.getAonehomoCa() << "\t"
				<< af.getAonehomoCa_exp() << "\t"
				<< af.getHetCa() << "\t"
				<< af.getHetCa_exp() << "\t"
				<< af.getAtwohomoCa() << "\t"
				<< af.getAtwohomoCa_exp() << "\t"
				//controls
				<< hwe.getControl() << "\t"
				<< af.getPopCon() << "\t"
				<< af.getAonehomoCon() << "\t"
				<< af.getAonehomoCon_exp() << "\t"
				<< af.getHetCon() << "\t"
				<< af.getHetCon_exp() << "\t"
				<< af.getAtwohomoCon() << "\t"
				<< af.getAtwohomoCon_exp();
			}
			pvals << endl;
			if(options.doParental()){
				paren << data_set->get_locus(i)->toString();
				if(data_set->get_locus(i)->isMicroSat()){
					for(int l = 0; l < 19; l++){
						paren << "\tNA";
					}
				}
				else{
					paren << "\t" << data_set->get_locus(i)->getAllele1() << "_" << data_set->get_locus(i)->getAllele1() << "\t"
					<< data_set->get_locus(i)->getAllele1() << "_" << data_set->get_locus(i)->getAllele2() << "\t"
					<< data_set->get_locus(i)->getAllele2() << "_" << data_set->get_locus(i)->getAllele2() << "\t"
					//parent male
					<< hwe.getParentalMale() << "\t"
					<< af.getPopPM() << "\t"
					<< af.getAonehomoPM() << "\t"
					<< af.getAonehomoPM_exp() << "\t"
					<< af.getAtwohomoPM() << "\t"
					<< af.getAtwohomoPM_exp() << "\t"
					<< af.getHetPM() << "\t"
					<< af.getHetPM_exp() << "\t"
					//parent female
					<< hwe.getParentalFemale() << "\t"
					<< af.getPopPF() << "\t"
					<< af.getAonehomoPF() << "\t"
					<< af.getAonehomoPF_exp() << "\t"
					<< af.getAtwohomoPF() << "\t"
					<< af.getAtwohomoPF_exp() << "\t"
					<< af.getHetPF() << "\t"
					<< af.getHetPF_exp();
				}
				paren << endl;
			}
			if(options.doGender()){
				gend << data_set->get_locus(i)->toString();
				if(data_set->get_locus(i)->isMicroSat()){
					for(int l = 0; l < 19; l++){
						gend << "\tNA";
					}
				}
				else{
					gend << "\t" << data_set->get_locus(i)->getAllele1() << "_" << data_set->get_locus(i)->getAllele1() << "\t"
					<< data_set->get_locus(i)->getAllele1() << "_" << data_set->get_locus(i)->getAllele2() << "\t"
					<< data_set->get_locus(i)->getAllele2() << "_" << data_set->get_locus(i)->getAllele2() << "\t";
					if(useoverall){
						//parent male
						gend << hwe.getOverallMale() << "\t"
						<< af.getPopM() << "\t"
						<< af.getAonehomoM() << "\t"
						<< af.getAonehomoM_exp() << "\t"
						<< af.getAtwohomoM() << "\t"
						<< af.getAtwohomoM_exp() << "\t"
						<< af.getHetM() << "\t"
						<< af.getHetM_exp() << "\t"
						//parent female
						<< hwe.getOverallFemale() << "\t"
						<< af.getPopF() << "\t"
						<< af.getAonehomoF() << "\t"
						<< af.getAonehomoF_exp() << "\t"
						<< af.getAtwohomoF() << "\t"
						<< af.getAtwohomoF_exp() << "\t"
						<< af.getHetF() << "\t"
						<< af.getHetF_exp();
					}
					else{
						//parent male
						gend << hwe.getParentalMale() << "\t"
						<< af.getPopPM() << "\t"
						<< af.getAonehomoPM() << "\t"
						<< af.getAonehomoPM_exp() << "\t"
						<< af.getAtwohomoPM() << "\t"
						<< af.getAtwohomoPM_exp() << "\t"
						<< af.getHetPM() << "\t"
						<< af.getHetPM_exp() << "\t"
						//parent female
						<< hwe.getParentalFemale() << "\t"
						<< af.getPopPF() << "\t"
						<< af.getAonehomoPF() << "\t"
						<< af.getAonehomoPF_exp() << "\t"
						<< af.getAtwohomoPF() << "\t"
						<< af.getAtwohomoPF_exp() << "\t"
						<< af.getHetPF() << "\t"
						<< af.getHetPF_exp();
					}
				}
				gend << endl;
			}
			if(options.doCaseControl()){
				cc << data_set->get_locus(i)->toString();
				if(data_set->get_locus(i)->isMicroSat()){
					for(int l = 0; l < 35; l++){
						cc << "\tNA";
					}
				}
				else{
					cc << "\t" << data_set->get_locus(i)->getAllele1() << "_" << data_set->get_locus(i)->getAllele1() << "\t"
					<< data_set->get_locus(i)->getAllele1() << "_" << data_set->get_locus(i)->getAllele2() << "\t"
					<< data_set->get_locus(i)->getAllele2() << "_" << data_set->get_locus(i)->getAllele2() << "\t"
					<< hwe.getCaseMale() << "\t"
					<< af.getPopCaM() << "\t"
					<< af.getAonehomoCaM() << "\t"
					<< af.getAonehomoCaM_exp() << "\t"
					<< af.getAtwohomoCaM() << "\t"
					<< af.getAtwohomoCaM_exp() << "\t"
					<< af.getHetCaM() << "\t"
					<< af.getHetCaM_exp() << "\t"

					<< hwe.getCaseFemale() << "\t"
					<< af.getPopCaF() << "\t"
					<< af.getAonehomoCaF() << "\t"
					<< af.getAonehomoCaF_exp() << "\t"
					<< af.getAtwohomoCaF() << "\t"
					<< af.getAtwohomoCaF_exp() << "\t"
					<< af.getHetCaF() << "\t"
					<< af.getHetCaF_exp() << "\t"

					<< hwe.getControlMale() << "\t"
					<< af.getPopConM() << "\t"
					<< af.getAonehomoConM() << "\t"
					<< af.getAonehomoConM_exp() << "\t"
					<< af.getAtwohomoConM() << "\t"
					<< af.getAtwohomoConM_exp() << "\t"
					<< af.getHetConM() << "\t"
					<< af.getHetConM_exp() << "\t"

					<< hwe.getControlFemale() << "\t"
					<< af.getPopConF() << "\t"
					<< af.getAonehomoConF() << "\t"
					<< af.getAonehomoConF_exp() << "\t"
					<< af.getAtwohomoConF() << "\t"
					<< af.getAtwohomoConF_exp() << "\t"
					<< af.getHetConF() << "\t"
					<< af.getHetConF_exp();
				}
				cc << endl;
			}

			doFilter(data_set->get_locus(i), &hwe);
		}
	}

	if(pvals.is_open()){
		pvals.close();
	}

}

