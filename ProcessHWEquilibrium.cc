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
#include <cdflib.h>
#include <General.h>
#include <Helpers.h>
#ifdef PLATOLIB
#include "Controller.h"
#endif

using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif
string ProcessHWEquilibrium::stepname = "hw";

#ifdef PLATOLIB
ProcessHWEquilibrium::ProcessHWEquilibrium(string bn, int pos, Database* pdb)
{
	name = "HWE";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
}
#endif

/*
 * Function: PrintSummary
 * Description:
 * Resets marker flags
 */
void ProcessHWEquilibrium::PrintSummary(){
	for(int i = 0; i < (int)data_set->num_loci(); i++){
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


////	int msize = data_set->num_inds();
	opts::printLog("Options:\t" + options.toString() + "\n");
    opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
	        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
	        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessHWEquilibrium::filter()
{
	#ifdef PLATOLIB
		for(int i = 0; i < (int)data_set->num_loci(); i++)
		{
			data_set->get_locus(i)->setFlag(false);
		}

		for(int i = 0; i < (int)data_set->get_families()->size(); i++)
		{
			data_set->get_pedigree(i)->setFlagAFCM(false);
			data_set->get_pedigree(i)->setFlagAFCF(false);
		}
	#endif
	return;
}//end method filter

/*
 * Function: doFilter
 * Description:
 * Filters markers based on Overall pvalue
 */
void ProcessHWEquilibrium::doFilter(Marker* mark, HWEquilibrium* hwe){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		bool inc = false;
		if(options.doThreshMarkersLow() && ((Helpers::dLess(hwe->getOverall(), options.getThreshMarkersLow()) && hwe->getOverall() != -1))){// || (dLess(hw_C,options.getThreshMarkersLow()) && hw_C != -1))){
			mark->setEnabled(false);
			inc = true;
		}
		if(options.doThreshMarkersHigh() && ((Helpers::dGreater(hwe->getOverall(), options.getThreshMarkersHigh()) && hwe->getOverall() != -1))){// || (dGreater(hw_C, options.getThreshMarkersHigh()) && hw_C != -1))){
			mark->setEnabled(false);
			inc = true;
		}
		if(inc){
			orig_num_markers++;
		}
	}
}//end method doFilter



/*
 * Function: process
 * Description:
 * Performs the HWE process.
 */
void ProcessHWEquilibrium::process(DataSet* ds){
	data_set = ds;
	#ifdef PLATOLIB
		create_tables();
		Query myQuery(*db);
	#endif

	useoverall = false;
	if(options.doRandomChild() || options.doAll() || options.doAllChildren()){
		useoverall = true;
	}

	good_markers = Helpers::findValidMarkers(data_set->get_markers(), &options);

	HWEquilibrium hwe(data_set);
	hwe.setOptions(options);
	AlleleFrequency af(data_set);
	af.setOptions(options);

	//	af->process_hw(samples, families, markers, marker_map);

	int msize = good_markers.size();//data_set->num_loci();

	if(options.doHWEPT()){
	#ifndef PLATOLIB
		//this section deals with setting up the output files and creating headers
		//not necessary for the Viewer...
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
	#endif
//		int prev_base = 0;
//		int prev_chrom = -1;
		#ifdef PLATOLIB
			myQuery.transaction();
		#endif
		for(int i = 0; i < msize; i++){
			if(good_markers[i]->isEnabled()){// && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
				hwe.calculateHWEPT(data_set->get_locus(i));
				#ifdef PLATOLIB
				string sql = hweptinsert;
				sql += "," + getString<int>(data_set->get_locus(i)->getSysprobe());
				sql += "," + getString<float>(hwe.getCasePvalHWEPT());
				sql += "," + getString<float>(hwe.getControlPvalHWEPT());
				sql += "," + getString<float>(hwe.getGenotypicPval());
				sql += "," + getString<float>(hwe.getAllelicPval());
				sql += "," + getString<float>(hwe.getCaseGeneralHWEPT());
				sql += "," + getString<float>(hwe.getControlGeneralHWEPT());
				sql += ")";
				Controller::execute_sql(myQuery, sql);
				#else
				pvals.precision(8);

				pvals << fixed << data_set->get_locus(i)->toString() << "\t" << hwe.getCasePvalHWEPT() << "\t"
					<< hwe.getControlPvalHWEPT() << "\t" << hwe.getGenotypicPval() << "\t"
					<< hwe.getAllelicPval() << "\t" << hwe.getCaseGeneralHWEPT() << "\t"
					<< hwe.getControlGeneralHWEPT() << endl;
				#endif
			}
		}
		#ifdef PLATOLIB
		myQuery.commit();
		#else
		pvals.close();
		#endif
		return;
	}

#ifndef PLATOLIB
	//more stuff dealing with file creation and headers. Not in Viewer
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
#endif
//	int prev_base = 0;
//	int prev_chrom = -1;
	#ifdef PLATOLIB
		myQuery.transaction();
	#endif
	for(int i = 0; i < msize; i++){
		if(good_markers[i]->isEnabled()){// && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
			hwe.calculate(good_markers[i]);
			af.calculate(good_markers[i]);

			#ifdef PLATOLIB
			string insert = defaultinsert;
			insert += "," + getString<int>(good_markers[i]->getSysprobe());
			#else
			pvals << good_markers[i]->toString();
			#endif
			if(good_markers[i]->isMicroSat()){
				for(int l = 0; l < 27; l++){
					#ifdef PLATOLIB
					insert += ",NULL";
					#else
					pvals << "\tNA";
					#endif
				}
			}
			else{
				//overall default = founders
				#ifdef PLATOLIB
				insert += ",'" + good_markers[i]->getAllele1() + "_" + good_markers[i]->getAllele1() + "'";
				insert += ",'" + good_markers[i]->getAllele1() + "_" + good_markers[i]->getAllele2() + "'";
				insert += ",'" + good_markers[i]->getAllele2() + "_" + good_markers[i]->getAllele2() + "'";
				insert += "," + getString<float>(hwe.getOverall());
				#else
				pvals << "\t" << good_markers[i]->getAllele1() << "_" << good_markers[i]->getAllele1() << "\t"
				<< good_markers[i]->getAllele1() << "_" << good_markers[i]->getAllele2() << "\t"
				<< good_markers[i]->getAllele2() << "_" << good_markers[i]->getAllele2() << "\t"
				<< hwe.getOverall() << "\t";
				#endif
				if(useoverall){
					#ifdef PLATOLIB
					insert += "," + getString<int>(af.getPop());
					insert += "," + getString<int>(af.getAonehomo());
					insert += "," + getString<float>(af.getAonehomo_exp());
					insert += "," + getString<int>(af.getHet());
					insert += "," + getString<float>(af.getHet_exp());
					insert += "," + getString<int>(af.getAtwohomo());
					insert += "," + getString<float>(af.getAtwohomo_exp());
					#else
					pvals << af.getPop() << "\t"
					<< af.getAonehomo() << "\t"
					<< af.getAonehomo_exp() << "\t"
					<< af.getHet() << "\t"
					<< af.getHet_exp() << "\t"
					<< af.getAtwohomo() << "\t"
					<< af.getAtwohomo_exp() << "\t";
					#endif
				}
				else{
					#ifdef PLATOLIB
					insert += "," + getString<int>(af.getPopP());
					insert += "," + getString<int>(af.getAonehomoP());
					insert += "," + getString<float>(af.getAonehomoP_exp());
					insert += "," + getString<int>(af.getHetP());
					insert += "," + getString<float>(af.getHetP_exp());
					insert += "," + getString<int>(af.getAtwohomoP());
					insert += "," + getString<float>(af.getAtwohomoP_exp());
					#else
					pvals << af.getPopP() << "\t"
					<< af.getAonehomoP() << "\t"
					<< af.getAonehomoP_exp() << "\t"
					<< af.getHetP() << "\t"
					<< af.getHetP_exp() << "\t"
					<< af.getAtwohomoP() << "\t"
					<< af.getAtwohomoP_exp() << "\t";
					#endif
				}
				#ifdef PLATOLIB
				//cases
				insert += "," + getString<float>(hwe.getCase());
				insert += "," + getString<int>(af.getPopCa());
				insert += "," + getString<int>(af.getAonehomoCa());
				insert += "," + getString<float>(af.getAonehomoCa_exp());
				insert += "," + getString<int>(af.getHetCa());
				insert += "," + getString<float>(af.getHetCa_exp());
				insert += "," + getString<int>(af.getAtwohomoCa());
				insert += "," + getString<float>(af.getAtwohomoCa_exp());
				//controls
				insert += "," + getString<float>(hwe.getControl());
				insert += "," + getString<int>(af.getPopCon());
				insert += "," + getString<int>(af.getAonehomoCon());
				insert += "," + getString<float>(af.getAonehomoCon_exp());
				insert += "," + getString<int>(af.getHetCon());
				insert += "," + getString<float>(af.getHetCon_exp());
				insert += "," + getString<int>(af.getAtwohomoCon());
				insert += "," + getString<float>(af.getAtwohomoCon_exp());
				#else
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
				#endif
			}
			#ifdef PLATOLIB
			insert += ")";
			Controller::execute_sql(myQuery, insert);
			#else
			pvals << endl;
			#endif
			if(options.doParental()){
				#ifdef PLATOLIB
				insert = parentalinsert;
				insert += "," + getString<int>(good_markers[i]->getSysprobe());
				#else
				paren << good_markers[i]->toString();
				#endif
				if(good_markers[i]->isMicroSat()){
					for(int l = 0; l < 19; l++){
						#ifdef PLATOLIB
						insert += ",NULL";
						#else
						paren << "\tNA";
						#endif
					}
				}
				else{
#ifdef PLATOLIB
					insert += ",'" + good_markers[i]->getAllele1() + "_" + good_markers[i]->getAllele1() + "'";
					insert += ",'" + good_markers[i]->getAllele1() + "_" + good_markers[i]->getAllele2() + "'";
					insert += ",'" + good_markers[i]->getAllele2() + "_" + good_markers[i]->getAllele2() + "'";
					//parent male
					insert += "," + getString<float>(hwe.getParentalMale());
					insert += "," + getString<int>(af.getPopPM());
					insert += "," + getString<int>(af.getAonehomoPM());
					insert += "," + getString<float>(af.getAonehomoPM_exp());
					insert += "," + getString<int>(af.getHetPM());
					insert += "," + getString<float>(af.getHetPM_exp());
					insert += "," + getString<int>(af.getAtwohomoPM());
					insert += "," + getString<float>(af.getAtwohomoPM_exp());
					//parent female
					insert += "," + getString<float>(hwe.getParentalFemale());
					insert += "," + getString<int>(af.getPopPF());
					insert += "," + getString<int>(af.getAonehomoPF());
					insert += "," + getString<float>(af.getAonehomoPF_exp());
					insert += "," + getString<int>(af.getHetPF());
					insert += "," + getString<float>(af.getHetPF_exp());
					insert += "," + getString<int>(af.getAtwohomoPF());
					insert += "," + getString<float>(af.getAtwohomoPF_exp());
				#else
					paren << "\t" << good_markers[i]->getAllele1() << "_" << good_markers[i]->getAllele1() << "\t"
					<< good_markers[i]->getAllele1() << "_" << good_markers[i]->getAllele2() << "\t"
					<< good_markers[i]->getAllele2() << "_" << good_markers[i]->getAllele2() << "\t"
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
					#endif
				}
				#ifdef PLATOLIB
				insert += ")";
				Controller::execute_sql(myQuery, insert);
				#else
				paren << endl;
				#endif
			}
			if(options.doGender()){
				#ifdef PLATOLIB
				insert = genderinsert;
				insert += "," + getString<int>(good_markers[i]->getSysprobe());
				#else
				gend << good_markers[i]->toString();
				#endif
				if(good_markers[i]->isMicroSat()){
					for(int l = 0; l < 19; l++){
						#ifdef PLATOLIB
						insert += ",NULL";
						#else
						gend << "\tNA";
						#endif
					}
				}
				else{
					#ifdef PLATOLIB
					insert += ",'" + good_markers[i]->getAllele1() + "_" + good_markers[i]->getAllele1() + "'";
					insert += ",'" + good_markers[i]->getAllele1() + "_" + good_markers[i]->getAllele2() + "'";
					insert += ",'" + good_markers[i]->getAllele2() + "_" + good_markers[i]->getAllele2() + "'";
					#else
					gend << "\t" << good_markers[i]->getAllele1() << "_" << good_markers[i]->getAllele1() << "\t"
					<< good_markers[i]->getAllele1() << "_" << good_markers[i]->getAllele2() << "\t"
					<< good_markers[i]->getAllele2() << "_" << good_markers[i]->getAllele2() << "\t";
					#endif
					if(useoverall){
						#ifdef PLATOLIB
						//parent male
						insert += "," + getString<float>(hwe.getOverallMale());
						insert += "," + getString<int>(af.getPopM());
						insert += "," + getString<int>(af.getAonehomoM());
						insert += "," + getString<float>(af.getAonehomoM_exp());
						insert += "," + getString<int>(af.getHetM());
						insert += "," + getString<float>(af.getHetM_exp());
						insert += "," + getString<int>(af.getAtwohomoM());
						insert += "," + getString<float>(af.getAtwohomoM_exp());
						//parent female
						insert += "," + getString<float>(hwe.getOverallFemale());
						insert += "," + getString<int>(af.getPopF());
						insert += "," + getString<int>(af.getAonehomoF());
						insert += "," + getString<float>(af.getAonehomoF_exp());
						insert += "," + getString<int>(af.getHetF());
						insert += "," + getString<float>(af.getHetF_exp());
						insert += "," + getString<int>(af.getAtwohomoF());
						insert += "," + getString<float>(af.getAtwohomoF_exp());
						#else
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
						#endif
					}
					else{
						#ifdef PLATOLIB
						//parent male
						insert += "," + getString<float>(hwe.getParentalMale());
						insert += "," + getString<int>(af.getPopPM());
						insert += "," + getString<int>(af.getAonehomoPM());
						insert += "," + getString<float>(af.getAonehomoPM_exp());
						insert += "," + getString<int>(af.getHetPM());
						insert += "," + getString<float>(af.getHetPM_exp());
						insert += "," + getString<int>(af.getAtwohomoPM());
						insert += "," + getString<float>(af.getAtwohomoPM_exp());
						//parent female
						insert += "," + getString<float>(hwe.getParentalFemale());
						insert += "," + getString<int>(af.getPopPF());
						insert += "," + getString<int>(af.getAonehomoPF());
						insert += "," + getString<float>(af.getAonehomoPF_exp());
						insert += "," + getString<int>(af.getHetPF());
						insert += "," + getString<float>(af.getHetPF_exp());
						insert += "," + getString<int>(af.getAtwohomoPF());
						insert += "," + getString<float>(af.getAtwohomoPF_exp());
						#else
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
						#endif
					}
				}
				#ifdef PLATOLIB
				insert += ")";
				Controller::execute_sql(myQuery, insert);
				#else
				gend << endl;
				#endif
			}
			if(options.doCaseControl()){
				#ifdef PLATOLIB
				insert == ccinsert;
				insert += "," + getString<int>(good_markers[i]->getSysprobe());
				#else
				cc << good_markers[i]->toString();
				#endif
				if(good_markers[i]->isMicroSat()){
					for(int l = 0; l < 35; l++){
						#ifdef PLATOLIB
						ccinsert += ",NULL";
						#else
						cc << "\tNA";
						#endif
					}
				}
				else{
					#ifdef PLATOLIB
					insert += ",'" + good_markers[i]->getAllele1() + "_" + good_markers[i]->getAllele1() + "'";
					insert += ",'" + good_markers[i]->getAllele1() + "_" + good_markers[i]->getAllele2() + "'";
					insert += ",'" + good_markers[i]->getAllele2() + "_" + good_markers[i]->getAllele2() + "'";

					insert += "," + getString<float>(hwe.getCaseMale());
					insert += "," + getString<int>(af.getPopCaM());
					insert += "," + getString<int>(af.getAonehomoCaM());
					insert += "," + getString<float>(af.getAonehomoCaM_exp());
					insert += "," + getString<int>(af.getHetCaM());
					insert += "," + getString<float>(af.getHetCaM_exp());
					insert += "," + getString<int>(af.getAtwohomoCaM());
					insert += "," + getString<float>(af.getAtwohomoCaM_exp());

					insert += "," + getString<float>(hwe.getCaseFemale());
					insert += "," + getString<int>(af.getPopCaF());
					insert += "," + getString<int>(af.getAonehomoCaF());
					insert += "," + getString<float>(af.getAonehomoCaF_exp());
					insert += "," + getString<int>(af.getHetCaF());
					insert += "," + getString<float>(af.getHetCaF_exp());
					insert += "," + getString<int>(af.getAtwohomoCaF());
					insert += "," + getString<float>(af.getAtwohomoCaF_exp());

					insert += "," + getString<float>(hwe.getControlMale());
					insert += "," + getString<int>(af.getPopConM());
					insert += "," + getString<int>(af.getAonehomoConM());
					insert += "," + getString<float>(af.getAonehomoConM_exp());
					insert += "," + getString<int>(af.getHetConM());
					insert += "," + getString<float>(af.getHetConM_exp());
					insert += "," + getString<int>(af.getAtwohomoConM());
					insert += "," + getString<float>(af.getAtwohomoConM_exp());

					insert += "," + getString<float>(hwe.getControlFemale());
					insert += "," + getString<int>(af.getPopConF());
					insert += "," + getString<int>(af.getAonehomoConF());
					insert += "," + getString<float>(af.getAonehomoConF_exp());
					insert += "," + getString<int>(af.getHetConF());
					insert += "," + getString<float>(af.getHetConF_exp());
					insert += "," + getString<int>(af.getAtwohomoConF());
					insert += "," + getString<float>(af.getAtwohomoConF_exp());
					#else
					cc << "\t" << good_markers[i]->getAllele1() << "_" << good_markers[i]->getAllele1() << "\t"
					<< good_markers[i]->getAllele1() << "_" << good_markers[i]->getAllele2() << "\t"
					<< good_markers[i]->getAllele2() << "_" << good_markers[i]->getAllele2() << "\t"
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
					#endif
				}
#ifdef PLATOLIB
				insert += ")";
				Controller::execute_sql(myQuery, insert);
#else
				cc << endl;
#endif
			}

			doFilter(good_markers[i], &hwe);
		}
	}
#ifdef PLATOLIB
	myQuery.commit();
#else
	if(pvals.is_open())
	{
		pvals.close();
	}
#endif

}//end method process()

#ifdef PLATOLIB
void ProcessHWEquilibrium::dump2db(){}

void ProcessHWEquilibrium::create_tables()
{
	Query myQuery(*db);
	myQuery.transaction();
    for(int i = 0; i < (int)tablename.size(); i++)
    {
        Controller::drop_table(db, tablename[i]);
    }
    headers.clear();
    tablename.clear();
    primary_table.clear();

    string tempbatch = batchname;
    for(int i = 0; i < (int)tempbatch.size(); i++)
    {
        if(tempbatch[i] == ' ')
        {
            tempbatch[i] = '_';
        }
    }
    string mytablename = tempbatch + "_";
    tempbatch = name;
    for(int i = 0; i < (int)tempbatch.size(); i++)
    {
        if(tempbatch[i] == ' ')
        {
            tempbatch[i] = '_';
        }
    }
    string base = mytablename + tempbatch;

    if(options.doHWEPT())
    {
        mytablename = base + "_hwept_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("HWE PT");
        primary_table[mytablename].push_back(Vars::LOCUS_TABLE);


        string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        sql += "PT_pval_case REAL,";
        sql += "PT_pval_control REAL,";
        sql += "P_pval_genotypic REAL,";
        sql += "P_pval_allelic REAL,";
        sql += "P_pval_case REAL,";
        sql += "P_pval_control REAL";
        sql += ")";
        hweptinsert = "INSERT INTO " + mytablename + "(id, fkey, PT_pval_case, PT_pval_control, P_pval_genotypic, P_pval_allelic, P_pval_case, P_pval_control) VALUES (NULL";
        Controller::execute_sql(myQuery, sql);

        return;
    }

    mytablename = base + "_" + getString<int>(position);
    tablename.push_back(mytablename);
    tablenicknames.push_back("");
    primary_table[mytablename].push_back(Vars::LOCUS_TABLE);


    string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
    sql += "fkey integer not null,";
    sql += "Genotype11 varchar(20),";
    sql += "Genotype12 varchar(20),";
    sql += "Genotype22 varchar(20),";
    sql += "Overall_Pvalue REAL,";
    sql += "Overall_Total_Count integer,";
    sql += "Overall_Obs_Genotype11 integer,";
    sql += "Overall_Exp_Genotype11 REAL,";
    sql += "Overall_Obs_Genotype12 integer,";
    sql += "Overall_Exp_Genotype12 REAL,";
    sql += "Overall_Obs_Genotype22 integer,";
    sql += "Overall_Exp_Genotype22 REAL,";
    sql += "Case_Pvalue REAL,";
    sql += "Case_Total_Count integer,";
    sql += "Case_Obs_Genotype11 integer,";
    sql += "Case_Exp_Genotype11 REAL,";
    sql += "Case_Obs_Genotype12 integer,";
    sql += "Case_Exp_Genotype12 REAL,";
    sql += "Case_Obs_Genotype22 integer,";
    sql += "Case_Exp_Genotype22 REAL,";
    sql += "Control_Pvalue REAL,";
    sql += "Control_Total_Count integer,";
    sql += "Control_Obs_Genotype11 integer,";
    sql += "Control_Exp_Genotype11 REAL,";
    sql += "Control_Obs_Genotype12 integer,";
    sql += "Control_Exp_Genotype12 REAL,";
    sql += "Control_Obs_Genotype22 integer,";
    sql += "Control_Exp_Genotype22 REAL";
    sql += ")";
    Controller::execute_sql(myQuery, sql);
    defaultinsert = "INSERT INTO " + mytablename + " (id, fkey, Genotype11, Genotype12, Genotype22";
    defaultinsert += ",Overall_Pvalue, Overall_Total_Count, Overall_Obs_Genotype11, Overall_Exp_Genotype11,";
    defaultinsert += "Overall_Obs_Genotype12, Overall_Exp_Genotype12, Overall_Obs_Genotype22, Overall_Exp_Genotype22";
    defaultinsert += ",Case_Pvalue, Case_Total_Count, Case_Obs_Genotype11, Case_Exp_Genotype11,";
    defaultinsert += "Case_Obs_Genotype12, Case_Exp_Genotype12, Case_Obs_Genotype22, Case_Exp_Genotype22";
    defaultinsert += ",Control_Pvalue, Control_Total_Count, Control_Obs_Genotype11, Control_Exp_Genotype11,";
    defaultinsert += "Control_Obs_Genotype12, Control_Exp_Genotype12, Control_Obs_Genotype22, Control_Exp_Genotype22";
    defaultinsert += ") VALUES (NULL";

    if(options.doParental())
    {
        mytablename = base + "_parental_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("Parental");
        primary_table[mytablename].push_back(Vars::LOCUS_TABLE);

        string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        sql += "Genotype11 varchar(20),";
        sql += "Genotype12 varchar(20),";
        sql += "Genotype22 varchar(20),";
        sql += "Parent_Male_Pvalue REAL,";
        sql += "Parent_Male_Total_Count integer,";
        sql += "Parent_Male_Obs_Genotype11 integer,";
        sql += "Parent_Male_Exp_Genotype11 REAL,";
        sql += "Parent_Male_Obs_Genotype12 integer,";
        sql += "Parent_Male_Exp_Genotype12 REAL,";
        sql += "Parent_Male_Obs_Genotype22 integer,";
        sql += "Parent_Male_Exp_Genotype22 REAL,";
        sql += "Parent_Female_Pvalue REAL,";
        sql += "Parent_Female_Total_Count integer,";
        sql += "Parent_Female_Obs_Genotype11 integer,";
        sql += "Parent_Female_Exp_Genotype11 REAL,";
        sql += "Parent_Female_Obs_Genotype12 integer,";
        sql += "Parent_Female_Exp_Genotype12 REAL,";
        sql += "Parent_Female_Obs_Genotype22 integer,";
        sql += "Parent_Female_Exp_Genotype22 REAL";
        sql += ")";
        Controller::execute_sql(myQuery, sql);
        parentalinsert = "INSERT INTO " + mytablename + " (id, fkey, Genotype11, Genotype12, Genotype22";
        parentalinsert += ",Parent_Male_Pvalue, Parent_Male_Total_Count, Parent_Male_Obs_Genotype11, Parent_Male_Exp_Genotype11,";
        parentalinsert += "Parent_Male_Obs_Genotype12, Parent_Male_Exp_Genotype12, Parent_Male_Obs_Genotype22, Parent_Male_Exp_Genotype22";
        parentalinsert += ",Parent_Female_Pvalue, Parent_Female_Total_Count, Parent_Female_Obs_Genotype11, Parent_Female_Exp_Genotype11,";
        parentalinsert += "Parent_Female_Obs_Genotype12, Parent_Female_Exp_Genotype12, Parent_Female_Obs_Genotype22, Parent_Female_Exp_Genotype22";
        parentalinsert += ") VALUES (NULL";
    }

    if(options.doGender())
    {
        mytablename = base + "_gender_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("Gender");
        primary_table[mytablename].push_back(Vars::LOCUS_TABLE);

        string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        sql += "Genotype11 varchar(20),";
        sql += "Genotype12 varchar(20),";
        sql += "Genotype22 varchar(20),";
        sql += "Overall_Male_Pvalue REAL,";
        sql += "Overall_Male_Total_Count integer,";
        sql += "Overall_Male_Obs_Genotype11 integer,";
        sql += "Overall_Male_Exp_Genotype11 REAL,";
        sql += "Overall_Male_Obs_Genotype12 integer,";
        sql += "Overall_Male_Exp_Genotype12 REAL,";
        sql += "Overall_Male_Obs_Genotype22 integer,";
        sql += "Overall_Male_Exp_Genotype22 REAL,";
        sql += "Overall_Female_Pvalue REAL,";
        sql += "Overall_Female_Total_Count integer,";
        sql += "Overall_Female_Obs_Genotype11 integer,";
        sql += "Overall_Female_Exp_Genotype11 REAL,";
        sql += "Overall_Female_Obs_Genotype12 integer,";
        sql += "Overall_Female_Exp_Genotype12 REAL,";
        sql += "Overall_Female_Obs_Genotype22 integer,";
        sql += "Overall_Female_Exp_Genotype22 REAL";
        sql += ")";
        Controller::execute_sql(myQuery, sql);
        genderinsert = "INSERT INTO " + mytablename + " (id, fkey, Genotype11, Genotype12, Genotype22";
        genderinsert += ",Overall_Male_Pvalue, Overall_Male_Total_Count, Overall_Male_Obs_Genotype11, Overall_Male_Exp_Genotype11,";
        genderinsert += "Overall_Male_Obs_Genotype12, Overall_Male_Exp_Genotype12, Overall_Male_Obs_Genotype22, Overall_Male_Exp_Genotype22";
        genderinsert += ",Overall_Female_Pvalue, Overall_Female_Total_Count, Overall_Female_Obs_Genotype11, Overall_Female_Exp_Genotype11,";
        genderinsert += "Overall_Female_Obs_Genotype12, Overall_Female_Exp_Genotype12, Overall_Female_Obs_Genotype22, Overall_Female_Exp_Genotype22";
        genderinsert += ") VALUES (NULL";

    }

    if(options.doCaseControl())
    {
        mytablename = base + "_casecontrol_" + getString<int>(position);
        tablename.push_back(mytablename);
        tablenicknames.push_back("Case/Control");
        primary_table[mytablename].push_back(Vars::LOCUS_TABLE);

        string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
        sql += "fkey integer not null,";
        sql += "Genotype11 varchar(20),";
        sql += "Genotype12 varchar(20),";
        sql += "Genotype22 varchar(20),";
        sql += "Case_Male_Pvalue REAL,";
        sql += "Case_Male_Total_Count integer,";
        sql += "Case_Male_Obs_Genotype11 integer,";
        sql += "Case_Male_Exp_Genotype11 REAL,";
        sql += "Case_Male_Obs_Genotype12 integer,";
        sql += "Case_Male_Exp_Genotype12 REAL,";
        sql += "Case_Male_Obs_Genotype22 integer,";
        sql += "Case_Male_Exp_Genotype22 REAL,";
        sql += "Case_Female_Pvalue REAL,";
        sql += "Case_Female_Total_Count integer,";
        sql += "Case_Female_Obs_Genotype11 integer,";
        sql += "Case_Female_Exp_Genotype11 REAL,";
        sql += "Case_Female_Obs_Genotype12 integer,";
        sql += "Case_Female_Exp_Genotype12 REAL,";
        sql += "Case_Female_Obs_Genotype22 integer,";
        sql += "Case_Female_Exp_Genotype22 REAL,";
        sql += "Control_Male_Pvalue REAL,";
        sql += "Control_Male_Total_Count integer,";
        sql += "Control_Male_Obs_Genotype11 integer,";
        sql += "Control_Male_Exp_Genotype11 REAL,";
        sql += "Control_Male_Obs_Genotype12 integer,";
        sql += "Control_Male_Exp_Genotype12 REAL,";
        sql += "Control_Male_Obs_Genotype22 integer,";
        sql += "Control_Male_Exp_Genotype22 REAL,";
        sql += "Control_Female_Pvalue REAL,";
        sql += "Control_Female_Total_Count integer,";
        sql += "Control_Female_Obs_Genotype11 integer,";
        sql += "Control_Female_Exp_Genotype11 REAL,";
        sql += "Control_Female_Obs_Genotype12 integer,";
        sql += "Control_Female_Exp_Genotype12 REAL,";
        sql += "Control_Female_Obs_Genotype22 integer,";
        sql += "Control_Female_Exp_Genotype22 REAL";
        sql += ")";
        Controller::execute_sql(myQuery, sql);
        ccinsert = "INSERT INTO " + mytablename + " (id, fkey, Genotype11, Genotype12, Genotype22";
        ccinsert += ",Case_Male_Pvalue, Case_Male_Total_Count, Case_Male_Obs_Genotype11, Case_Male_Exp_Genotype11,";
        ccinsert += "Case_Male_Obs_Genotype12, Case_Male_Exp_Genotype12, Case_Male_Obs_Genotype22, Case_Male_Exp_Genotype22";
        ccinsert += ",Case_Female_Pvalue, Case_Female_Total_Count, Case_Female_Obs_Genotype11, Case_Female_Exp_Genotype11,";
        ccinsert += "Case_Female_Obs_Genotype12, Case_Female_Exp_Genotype12, Case_Female_Obs_Genotype22, Case_Female_Exp_Genotype22";
        ccinsert += ",Control_Male_Pvalue, Control_Male_Total_Count, Control_Male_Obs_Genotype11, Control_Male_Exp_Genotype11,";
        ccinsert += "Control_Male_Obs_Genotype12, Control_Male_Exp_Genotype12, Control_Male_Obs_Genotype22, Control_Male_Exp_Genotype22";
        ccinsert += ",Control_Female_Pvalue, Control_Female_Total_Count, Control_Female_Obs_Genotype11, Control_Female_Exp_Genotype11,";
        ccinsert += "Control_Female_Obs_Genotype12, Control_Female_Exp_Genotype12, Control_Female_Obs_Genotype22, Control_Female_Exp_Genotype22";
        ccinsert += ") VALUES (NULL";

    }
    myQuery.commit();
}

void ProcessHWEquilibrium::run(DataSetObject* ds)
{
	process(ds);
}

#endif

#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
