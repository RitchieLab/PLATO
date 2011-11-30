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
#include "HWEquilibrium.h"
#include "cdflib.h"
#include "General.h"
#include "Helpers.h"
#include "CaConChisq.h"
#include "ChiSquare.h"
#include "ChiSquareAllelic.h"
#include "cdflib.h"
namespace Methods{
string HWEquilibrium::stepname = "hwe";

/*
 * Function: PrintSummary
 * Description:
 * Resets marker flags
 */
void HWEquilibrium::PrintSummary(){
	for(int i = 0; i < (int)markers->size(); i++){
		(*markers)[i]->setFlag(false);
	}
	return;
}

/*
 * Function: FilterSummary
 * Description:
 * Outputs number of remaining markers
 */
void HWEquilibrium::FilterSummary(){
	int fsize = families->size();
	for(int i = 0; i < fsize; i++){
		(*families)[i]->setFlagAFCM(false);
		(*families)[i]->setFlagAFCF(false);
	}


	opts::printLog("Options:\t" + options.toString() + "\n");
    opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
	        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
	        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

	if(af != NULL){
		delete(af);
	}
}

void HWEquilibrium::filter(){
	return;
}

/*
 * Function: doFilter
 * Description:
 * Filters markers based on Overall pvalue
 */
void HWEquilibrium::doFilter(Marker* mark){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		bool inc = false;
		if(options.doThreshMarkersLow() && ((Helpers::dLess(hw_O, options.getThreshMarkersLow()) && hw_O != -1))){
			mark->setEnabled(false);
			inc = true;
		}
		if(options.doThreshMarkersHigh() && ((Helpers::dGreater(hw_O, options.getThreshMarkersHigh()) && hw_O != -1))){
			mark->setEnabled(false);
			inc = true;
		}
		if(inc){
			orig_num_markers++;
		}
	}

}

void HWEquilibrium::calculateHWEPT(Marker* mark){
	hw_Ca = 0;
	hw_Con = 0;
	geno_chi = 0;
	ChiSquareAllelic chi;
	vector<vector<int> > chis;
	chis.push_back(vector<int>(0));
	chis.push_back(vector<int>(0));
	vector<bool> sample_status(samples->size(), false);
	for(int i = 0; i < (int)samples->size(); i++){
		if((*samples)[i]->getPheno() == 2){
			sample_status[i] = true;
		}
		else if((*samples)[i]->getPheno() == 1){
			sample_status[i] = false;
		}
	}

	af->calcOne(mark);
	CaConChisq CCC;
	geno_chi = CCC.calcChiGeno(af->getAonehomoCa(), af->getHetCa(), af->getAtwohomoCa(), af->getAonehomoCon(), af->getHetCon(), af->getAtwohomoCon());
	geno_pval = -1;
	if(geno_chi >= 0){
				double pvalue, df = 2;
				pvalue = -1;
				if(geno_chi > -1){
					pvalue = Helpers::p_from_chi(geno_chi, df);
				}
				geno_pval = pvalue;
			}
			chis[0].push_back(af->getAoneCa_count());
			chis[0].push_back(af->getAtwoCa_count());
			chis[1].push_back(af->getAoneCon_count());
			chis[1].push_back(af->getAtwoCon_count());
			allele_chi = chi.chisquare(chis);
			allele_p = -1;
			if(allele_chi > -1){
				double pvalue, df = 1;
				pvalue = -1;
				if(allele_chi > -1){
					pvalue = Helpers::p_from_chi(allele_chi, df);
				}
				allele_p = pvalue;
			}

			hw_Ca = calcHW(af->getAoneCa_freq(), af->getAtwoCa_freq(), af->getAonehomoCa(), af->getHetCa(), af->getAtwohomoCa(), af->getPopCa());
			hw_Con = calcHW(af->getAoneCon_freq(), af->getAtwoCon_freq(), af->getAonehomoCon(), af->getHetCon(), af->getAtwohomoCon(), af->getPopCon());

			vector<double> results = hwePT(sample_status, mark);
			int ca_pt = 0;
			int con_pt = 0;
			for(int p = 0; p < options.getPermutations(); p++){
				vector<bool> temp_status = sample_status;
				random_shuffle(temp_status.begin(), temp_status.end());
				vector<double> temp = hwePT(temp_status, mark);
				if(temp[0] < results[0]){
					ca_pt += 0;
				}
				else{
					ca_pt += 1;
				}
				if(temp[1] < results[1]){
					con_pt += 0;
				}
				else{
					con_pt += 1;
				}
			}

			ca_pval = (double)ca_pt/(double)options.getPermutations();
			con_pval = (double)con_pt/(double)options.getPermutations();
}

/*
 * Function: processHWEPT
 * Description:
 * Perform HW calc based on Kelli alg.
 */
void HWEquilibrium::processHWEPT(){
	string fname = opts::_OUTPREFIX_ + "hw" + options.getOut() + ".txt";
	if(!overwrite){
		fname += "." + getString<int>(order);
	}
	ofstream pvals (fname.c_str());
	if(!pvals){
		opts::printLog("Error opening " + fname + "!  Exiting!\n");
		throw MethodException("Error opening " + fname + "!  Exiting!\n");
	}
	pvals.precision(4);

	opts::addFile("Marker", stepname, fname);
	pvals << "Chrom\trsID\tProbeID\tbploc";
	if((*markers)[0]->getDetailHeaders().size() > 0){
		pvals << "\t" << (*markers)[0]->getDetailHeaders();
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
	int msize = markers->size();

	ChiSquareAllelic chi;
	vector<vector<int> > chis;
	chis.push_back(vector<int>(0));
	chis.push_back(vector<int>(0));
	vector<bool> sample_status(samples->size(), false);
	for(int i = 0; i < (int)samples->size(); i++){
		if((*samples)[i]->getPheno() == 2){
			sample_status[i] = true;
		}
		else if((*samples)[i]->getPheno() == 1){
			sample_status[i] = false;
		}
	}

	for(int i = 0; i < msize; i++){
		Marker* mark = (*markers)[i];

		if(mark->isEnabled() && Helpers::isValidMarker(mark, &options, prev_base, prev_chrom)){
			af->calcOne(mark);
			CaConChisq CCC;
			double geno_chi = CCC.calcChiGeno(af->getAonehomoCa(), af->getHetCa(), af->getAtwohomoCa(), af->getAonehomoCon(), af->getHetCon(), af->getAtwohomoCon());
			double geno_pval = -1;
			if(geno_chi >= 0){
				double pvalue, df = 2;
				pvalue = -1;
				if(geno_chi > -1){
					pvalue = Helpers::p_from_chi(geno_chi, df);
				}
				geno_pval = pvalue;
			}
			chis[0].push_back(af->getAoneCa_count());
			chis[0].push_back(af->getAtwoCa_count());
			chis[1].push_back(af->getAoneCon_count());
			chis[1].push_back(af->getAtwoCon_count());
			double allele_chi = chi.chisquare(chis);
			double allele_p = -1;
			if(allele_chi > -1){
				double pvalue, df = 1;
				pvalue = -1;
				if(allele_chi > -1){
					pvalue = Helpers::p_from_chi(allele_chi, df);
				}
				allele_p = pvalue;
			}

			hw_Ca = calcHW(af->getAoneCa_freq(), af->getAtwoCa_freq(), af->getAonehomoCa(), af->getHetCa(), af->getAtwohomoCa(), af->getPopCa());
			hw_Con = calcHW(af->getAoneCon_freq(), af->getAtwoCon_freq(), af->getAonehomoCon(), af->getHetCon(), af->getAtwohomoCon(), af->getPopCon());

			vector<double> results = hwePT(sample_status, mark);
			double ca_pt = 0;
			double con_pt = 0;
			for(int p = 0; p < options.getPermutations(); p++){
				vector<bool> temp_status = sample_status;
				random_shuffle(temp_status.begin(), temp_status.end());
				vector<double> temp = hwePT(temp_status, mark);
				if(temp[0] < results[0]){
					ca_pt += temp[0];
				}
				else{
					ca_pt += results[0];
				}
				if(temp[1] < results[1]){
					con_pt += temp[1];
				}
				else{
					con_pt += results[1];
				}
			}

			ca_pt = ca_pt/options.getPermutations();
			con_pt = con_pt/options.getPermutations();
			double ca_pval = 0;
			double con_pval = 0;
			if(ca_pt >= 0){
				double pvalue, df = 2;
				pvalue = -1;
				if(ca_pt > -1){
					pvalue = Helpers::p_from_chi(ca_pt, df);
				}
				ca_pval = pvalue;
			}
			if(con_pt >= 0){
				double pvalue, df = 2;
				pvalue = -1;
				if(con_pt > -1){
					pvalue = Helpers::p_from_chi(con_pt, df);
				}
				con_pval = pvalue;
			}

			pvals.precision(8);
			pvals << fixed << mark->toString() << "\t" << ca_pval << "\t" << con_pval << "\t" << geno_pval << "\t" << allele_p << "\t" << hw_Ca << "\t" << hw_Con << endl;
		}
	}
	pvals.close();
}

/*
 * Function: hwePT
 * Description:
 * Performs the hwePT calculation
 */
vector<double> HWEquilibrium::hwePT(vector<bool> sample_status, Marker* mark){
	//calculate counts
	int ca1homo = 0;
	int ca2homo = 0;
	int cahet = 0;
	int ca1count = 0;
	int ca2count = 0;
	int con1homo = 0;
	int con2homo = 0;
	int conhet = 0;
	int con1count = 0;
	int con2count = 0;
	int conzero = 0;
	int cazero = 0;

	int mloc = mark->getLoc();
	for(int i = 0; i < (int)samples->size(); i++){
		Sample* samp = (*samples)[i];
		if(samp->isEnabled() && (samp->getPheno() == 2 || samp->getPheno() == 1)){
			if(sample_status[i]){
				if(samp->getAone(mloc)){
					if(samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
						ca2count+=2;
						ca2homo++;
					}
					else if(samp->getAtwo(mloc) && samp->getAmissing(mloc)){
						cazero++;
					}
				}
				else{
					if(samp->getAtwo(mloc)){
						ca1count++;
						ca2count++;
						cahet++;
					}
					else{
						ca1count+=2;
						ca1homo++;
					}
				}
			}
			else{
				if(samp->getAone(mloc)){
					if(samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
						con2count+=2;
						con2homo++;
					}
					else if(samp->getAtwo(mloc) && samp->getAmissing(mloc)){
						conzero++;
					}
				}
				else{
					if(samp->getAtwo(mloc)){
						con1count++;
						con2count++;
						conhet++;
					}
					else{
						con1count+=2;
						con1homo++;
					}
				}
			}
		}
	}
	double ca_tot = ca1homo + cahet + ca2homo;
	double con_tot = con1homo + conhet + con2homo;
	double p_case = (ca2homo + (0.5*cahet))/ca_tot;
	double p_cont = (con2homo + (0.5*conhet))/con_tot;
	double ptot = options.getPrevalance() * p_case + (1-options.getPrevalance())*p_cont;
	double p1 = ptot * ptot;
	double p2 = 2*ptot*(1-ptot);
	double p3 = (1-ptot) * (1-ptot);
	double ca_chi3 = (((ca1homo / ca_tot) * (ca1homo / ca_tot)) / p3);
	double ca_chi2 = (((cahet / ca_tot) * (cahet / ca_tot)) / p2);
	double ca_chi1 = (((ca2homo / ca_tot) * (ca2homo / ca_tot)) / p1);
	double ca_chi = ca_chi3 + ca_chi2 + ca_chi1;
	ca_chi = (ca_chi - 1) * ca_tot;
	double con_chi3 = (((con1homo / con_tot) * (con1homo / con_tot)) / p3);
	double con_chi2 = (((conhet / con_tot) * (conhet / con_tot)) / p2);
	double con_chi1 = (((con2homo / con_tot) * (con2homo / con_tot)) / p1);
	double con_chi = con_chi3 + con_chi2 + con_chi1;
	con_chi = (con_chi - 1) * con_tot;
	vector<double> results;
	results.push_back(ca_chi);
	results.push_back(con_chi);
	return results;
}

/*
 * Function: doFilterHWEPT
 * Description:
 * Filters snps based on HWEPT calculations and thresholds
 */
void HWEquilibrium::doFilterHWEPT(Marker* mark){
}


void HWEquilibrium::calculate(Marker* mark){
	if(options.doHWEPT()){
		calculateHWEPT(mark);
	}
	else{
		af->calcOne(mark);
		if(useoverall){
			if(mark->getChrom() != opts::_CHRX_){
				hw_O = calcHW(af->getAone_freq(), af->getAtwo_freq(), af->getAonehomo(), af->getHet(), af->getAtwohomo(), af->getPop());
				hw_OM = calcHW(af->getAoneM_freq(), af->getAtwoM_freq(), af->getAonehomoM(), af->getHetM(), af->getAtwohomoM(), af->getPopM());
				hw_OF = calcHW(af->getAoneF_freq(), af->getAtwoF_freq(), af->getAonehomoF(), af->getHetF(), af->getAtwohomoF(), af->getPopF());
			}
			else{
				hw_O = calcHW(af->getAoneF_freq(), af->getAtwoF_freq(), af->getAonehomoF(), af->getHetF(), af->getAtwohomoF(), af->getPopF());
				hw_OM = calcHW(af->getAoneF_freq(), af->getAtwoF_freq(), af->getAonehomoF(), af->getHetF(), af->getAtwohomoF(), af->getPopF());
				hw_OF = calcHW(af->getAoneF_freq(), af->getAtwoF_freq(), af->getAonehomoF(), af->getHetF(), af->getAtwohomoF(), af->getPopF());
			}

		}
		else{
			if(mark->getChrom() != opts::_CHRX_){
				hw_O = calcHW(af->getAoneP_freq(), af->getAtwoP_freq(), af->getAonehomoP(), af->getHetP(), af->getAtwohomoP(), af->getPopP());
				hw_OM = calcHW(af->getAonePM_freq(), af->getAtwoPM_freq(), af->getAonehomoPM(), af->getHetPM(), af->getAtwohomoPM(), af->getPopPM());
				hw_OF = calcHW(af->getAonePF_freq(), af->getAtwoPF_freq(), af->getAonehomoPF(), af->getHetPF(), af->getAtwohomoPF(), af->getPopPF());
			}
			else{
				hw_O = calcHW(af->getAonePF_freq(), af->getAtwoPF_freq(), af->getAonehomoPF(), af->getHetPF(), af->getAtwohomoPF(), af->getPopPF());
				hw_OM = calcHW(af->getAonePF_freq(), af->getAtwoPF_freq(), af->getAonehomoPF(), af->getHetPF(), af->getAtwohomoPF(), af->getPopPF());
				hw_OF = calcHW(af->getAonePF_freq(), af->getAtwoPF_freq(), af->getAonehomoPF(), af->getHetPF(), af->getAtwohomoPF(), af->getPopPF());
			}

		}
		if(mark->getChrom() != opts::_CHRX_){
			hw_P = calcHW(af->getAoneP_freq(), af->getAtwoP_freq(), af->getAonehomoP(), af->getHetP(), af->getAtwohomoP(), af->getPopP());
			hw_PM = calcHW(af->getAonePM_freq(), af->getAtwoPM_freq(), af->getAonehomoPM(), af->getHetPM(), af->getAtwohomoPM(), af->getPopPM());
			hw_PF = calcHW(af->getAonePF_freq(), af->getAtwoPF_freq(), af->getAonehomoPF(), af->getHetPF(), af->getAtwohomoPF(), af->getPopPF());
			hw_C = calcHW(af->getAoneC_freq(), af->getAtwoC_freq(), af->getAonehomoC(), af->getHetC(), af->getAtwohomoC(), af->getPopC());
			hw_CM = calcHW(af->getAoneCM_freq(), af->getAtwoCM_freq(), af->getAonehomoCM(), af->getHetCM(), af->getAtwohomoCM(), af->getPopCM());
			hw_CF = calcHW(af->getAoneCF_freq(), af->getAtwoCF_freq(), af->getAonehomoCF(), af->getHetCF(), af->getAtwohomoCF(), af->getPopCF());
			hw_Ca = calcHW(af->getAoneCa_freq(), af->getAtwoCa_freq(), af->getAonehomoCa(), af->getHetCa(), af->getAtwohomoCa(), af->getPopCa());
			hw_CaM = calcHW(af->getAoneCaM_freq(), af->getAtwoCaM_freq(), af->getAonehomoCaM(), af->getHetCaM(), af->getAtwohomoCaM(), af->getPopCaM());
			hw_CaF = calcHW(af->getAoneCaF_freq(), af->getAtwoCaF_freq(), af->getAonehomoCaF(), af->getHetCaF(), af->getAtwohomoCaF(), af->getPopCaF());
			hw_Con = calcHW(af->getAoneCon_freq(), af->getAtwoCon_freq(), af->getAonehomoCon(), af->getHetCon(), af->getAtwohomoCon(), af->getPopCon());
			hw_ConM = calcHW(af->getAoneConM_freq(), af->getAtwoConM_freq(), af->getAonehomoConM(), af->getHetConM(), af->getAtwohomoConM(), af->getPopConM());
			hw_ConF = calcHW(af->getAoneConF_freq(), af->getAtwoConF_freq(), af->getAonehomoConF(), af->getHetConF(), af->getAtwohomoConF(), af->getPopConF());
		}
		else{
			hw_P = calcHW(af->getAonePF_freq(), af->getAtwoPF_freq(), af->getAonehomoPF(), af->getHetPF(), af->getAtwohomoPF(), af->getPopPF());
			hw_PM = calcHW(af->getAonePF_freq(), af->getAtwoPF_freq(), af->getAonehomoPF(), af->getHetPF(), af->getAtwohomoPF(), af->getPopPF());
			hw_PF = calcHW(af->getAonePF_freq(), af->getAtwoPF_freq(), af->getAonehomoPF(), af->getHetPF(), af->getAtwohomoPF(), af->getPopPF());
			hw_C = calcHW(af->getAoneCF_freq(), af->getAtwoCF_freq(), af->getAonehomoCF(), af->getHetCF(), af->getAtwohomoCF(), af->getPopCF());
			hw_CM = calcHW(af->getAoneCF_freq(), af->getAtwoCF_freq(), af->getAonehomoCF(), af->getHetCF(), af->getAtwohomoCF(), af->getPopCF());
			hw_CF = calcHW(af->getAoneCF_freq(), af->getAtwoCF_freq(), af->getAonehomoCF(), af->getHetCF(), af->getAtwohomoCF(), af->getPopCF());
			hw_Ca = calcHW(af->getAoneCaF_freq(), af->getAtwoCaF_freq(), af->getAonehomoCaF(), af->getHetCaF(), af->getAtwohomoCaF(), af->getPopCaF());
			hw_CaM = calcHW(af->getAoneCaF_freq(), af->getAtwoCaF_freq(), af->getAonehomoCaF(), af->getHetCaF(), af->getAtwohomoCaF(), af->getPopCaF());
			hw_CaF = calcHW(af->getAoneCaF_freq(), af->getAtwoCaF_freq(), af->getAonehomoCaF(), af->getHetCaF(), af->getAtwohomoCaF(), af->getPopCaF());
			hw_Con = calcHW(af->getAoneConF_freq(), af->getAtwoConF_freq(), af->getAonehomoConF(), af->getHetConF(), af->getAtwohomoConF(), af->getPopConF());
			hw_ConM = calcHW(af->getAoneConF_freq(), af->getAtwoConF_freq(), af->getAonehomoConF(), af->getHetConF(), af->getAtwohomoConF(), af->getPopConF());
			hw_ConF = calcHW(af->getAoneConF_freq(), af->getAtwoConF_freq(), af->getAonehomoConF(), af->getHetConF(), af->getAtwohomoConF(), af->getPopConF());
		}

	}
}


/*
 * Function: process
 * Description:
 * Performs the HWE process.
 */
void HWEquilibrium::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
	useoverall = false;
	af = new AlleleFrequency(samples, families);
	af->setOptions(options);
	af->setRank(rank);
	if(options.doRandomChild() || options.doAll() || options.doAllChildren()){
		af->flagSamples();
		useoverall = true;
	}
	if(options.doHWEPT()){
		processHWEPT();
		return;
	}

	int msize = markers->size();

	string fname1 = opts::_OUTPREFIX_ + "hw" + options.getOut() + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	string fname2 = opts::_OUTPREFIX_ + "hw_parental" + options.getOut() + ".txt";
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	string fname3 = opts::_OUTPREFIX_ + "hw_gender" + options.getOut() + ".txt";
	if(!overwrite){
		fname3 += "." + getString<int>(order);
	}
	string fname4 = opts::_OUTPREFIX_ + "hw_casecontrol" + options.getOut() + ".txt";
	if(!overwrite){
		fname4 += "." + getString<int>(order);
	}
	ofstream paren;
	ofstream gend;
	ofstream cc;
	ofstream pvals (fname1.c_str());
	if(!pvals){
		opts::printLog("Error opening " + fname1 + "!  Exiting!\n");
		throw MethodException("Error opening " + fname1 + "!  Exiting!\n");
	}
	pvals.precision(4);
	opts::addFile("Marker",stepname,fname1);
	if(options.doParental()){
		paren.open(fname2.c_str());
		if(!paren){
			opts::printLog("Error opening " + fname2 + "! Exiting!\n");
			throw MethodException("Error opening " + fname2 + "! Exiting!\n");
		}
		paren.precision(4);
		opts::addFile("Marker",stepname,fname2);
	}
	if(options.doGender()){
		gend.open(fname3.c_str());
		if(!gend){
			opts::printLog("Error opening " + fname3 + "! Exiting!\n");
			throw MethodException("Error opening " + fname3 + "! Exiting!\n");
		}
		gend.precision(4);
		opts::addFile("Marker",stepname,fname3);
	}
	if(options.doCaseControl()){
		cc.open(fname4.c_str());
		if(!cc){
			opts::printLog("Error opening " + fname4 + "! Exiting!\n");
			throw MethodException("Error opening " + fname4 + "! Exiting!\n");
		}
		cc.precision(4);
		opts::addFile("Marker",stepname,fname4);
	}

    pvals << "Chrom"
		  << "\trsID"
		  << "\tProbeID"
		  << "\tbploc";
	if((*markers)[0]->getDetailHeaders().size() > 0){
		pvals << "\t" << (*markers)[0]->getDetailHeaders();
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
	if((*markers)[0]->getDetailHeaders().size() > 0){
		paren << "\t" << (*markers)[0]->getDetailHeaders();
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
	if((*markers)[0]->getDetailHeaders().size() > 0){
		gend << "\t" << (*markers)[0]->getDetailHeaders();
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
	if((*markers)[0]->getDetailHeaders().size() > 0){
		cc << "\t" << (*markers)[0]->getDetailHeaders();
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
		if((*markers)[i]->isEnabled() && Helpers::isValidMarker((*markers)[i], &options, prev_base, prev_chrom)){
			af->calcOne((*markers)[i]);
			if(useoverall){
				hw_O = calcHW(af->getAone_freq(), af->getAtwo_freq(), af->getAonehomo(), af->getHet(), af->getAtwohomo(), af->getPop());
				hw_OM = calcHW(af->getAoneM_freq(), af->getAtwoM_freq(), af->getAonehomoM(), af->getHetM(), af->getAtwohomoM(), af->getPopM());
				hw_OF = calcHW(af->getAoneF_freq(), af->getAtwoF_freq(), af->getAonehomoF(), af->getHetF(), af->getAtwohomoF(), af->getPopF());
			}
			else{
				hw_O = calcHW(af->getAoneP_freq(), af->getAtwoP_freq(), af->getAonehomoP(), af->getHetP(), af->getAtwohomoP(), af->getPopP());
				hw_OM = calcHW(af->getAonePM_freq(), af->getAtwoPM_freq(), af->getAonehomoPM(), af->getHetPM(), af->getAtwohomoPM(), af->getPopPM());
				hw_OF = calcHW(af->getAonePF_freq(), af->getAtwoPF_freq(), af->getAonehomoPF(), af->getHetPF(), af->getAtwohomoPF(), af->getPopPF());
			}
			hw_P = calcHW(af->getAoneP_freq(), af->getAtwoP_freq(), af->getAonehomoP(), af->getHetP(), af->getAtwohomoP(), af->getPopP());
			hw_PM = calcHW(af->getAonePM_freq(), af->getAtwoPM_freq(), af->getAonehomoPM(), af->getHetPM(), af->getAtwohomoPM(), af->getPopPM());
			hw_PF = calcHW(af->getAonePF_freq(), af->getAtwoPF_freq(), af->getAonehomoPF(), af->getHetPF(), af->getAtwohomoPF(), af->getPopPF());
			hw_C = calcHW(af->getAoneC_freq(), af->getAtwoC_freq(), af->getAonehomoC(), af->getHetC(), af->getAtwohomoC(), af->getPopC());
			hw_CM = calcHW(af->getAoneCM_freq(), af->getAtwoCM_freq(), af->getAonehomoCM(), af->getHetCM(), af->getAtwohomoCM(), af->getPopCM());
			hw_CF = calcHW(af->getAoneCF_freq(), af->getAtwoCF_freq(), af->getAonehomoCF(), af->getHetCF(), af->getAtwohomoCF(), af->getPopCF());
			hw_Ca = calcHW(af->getAoneCa_freq(), af->getAtwoCa_freq(), af->getAonehomoCa(), af->getHetCa(), af->getAtwohomoCa(), af->getPopCa());
			hw_CaM = calcHW(af->getAoneCaM_freq(), af->getAtwoCaM_freq(), af->getAonehomoCaM(), af->getHetCaM(), af->getAtwohomoCaM(), af->getPopCaM());
			hw_CaF = calcHW(af->getAoneCaF_freq(), af->getAtwoCaF_freq(), af->getAonehomoCaF(), af->getHetCaF(), af->getAtwohomoCaF(), af->getPopCaF());
			hw_Con = calcHW(af->getAoneCon_freq(), af->getAtwoCon_freq(), af->getAonehomoCon(), af->getHetCon(), af->getAtwohomoCon(), af->getPopCon());
			hw_ConM = calcHW(af->getAoneConM_freq(), af->getAtwoConM_freq(), af->getAonehomoConM(), af->getHetConM(), af->getAtwohomoConM(), af->getPopConM());
			hw_ConF = calcHW(af->getAoneConF_freq(), af->getAtwoConF_freq(), af->getAonehomoConF(), af->getHetConF(), af->getAtwohomoConF(), af->getPopConF());

			pvals << (*markers)[i]->toString();
			if((*markers)[i]->isMicroSat()){
				for(int l = 0; l < 27; l++){
					pvals << "\tNA";
				}
			}
			else{
				//overall default = founders
				pvals << "\t" << (*markers)[i]->getAllele1() << "_" << (*markers)[i]->getAllele1() << "\t"
				<< (*markers)[i]->getAllele1() << "_" << (*markers)[i]->getAllele2() << "\t"
				<< (*markers)[i]->getAllele2() << "_" << (*markers)[i]->getAllele2() << "\t"
				<< hw_O << "\t";
				if(useoverall){
					pvals << af->getPop() << "\t"
					<< af->getAonehomo() << "\t"
					<< af->getAonehomo_exp() << "\t"
					<< af->getHet() << "\t"
					<< af->getHet_exp() << "\t"
					<< af->getAtwohomo() << "\t"
					<< af->getAtwohomo_exp() << "\t";
				}
				else{
					pvals << af->getPopP() << "\t"
					<< af->getAonehomoP() << "\t"
					<< af->getAonehomoP_exp() << "\t"
					<< af->getHetP() << "\t"
					<< af->getHetP_exp() << "\t"
					<< af->getAtwohomoP() << "\t"
					<< af->getAtwohomoP_exp() << "\t";
				}
				//cases
				pvals << hw_Ca << "\t"
				<< af->getPopCa() << "\t"
				<< af->getAonehomoCa() << "\t"
				<< af->getAonehomoCa_exp() << "\t"
				<< af->getHetCa() << "\t"
				<< af->getHetCa_exp() << "\t"
				<< af->getAtwohomoCa() << "\t"
				<< af->getAtwohomoCa_exp() << "\t"
				//controls
				<< hw_Con << "\t"
				<< af->getPopCon() << "\t"
				<< af->getAonehomoCon() << "\t"
				<< af->getAonehomoCon_exp() << "\t"
				<< af->getHetCon() << "\t"
				<< af->getHetCon_exp() << "\t"
				<< af->getAtwohomoCon() << "\t"
				<< af->getAtwohomoCon_exp();
			}
			pvals << endl;
			if(options.doParental()){
				paren << (*markers)[i]->toString();
				if((*markers)[i]->isMicroSat()){
					for(int l = 0; l < 19; l++){
						paren << "\tNA";
					}
				}
				else{
					paren << "\t" << (*markers)[i]->getAllele1() << "_" << (*markers)[i]->getAllele1() << "\t"
					<< (*markers)[i]->getAllele1() << "_" << (*markers)[i]->getAllele2() << "\t"
					<< (*markers)[i]->getAllele2() << "_" << (*markers)[i]->getAllele2() << "\t"
					//parent male
					<< hw_PM << "\t"
					<< af->getPopPM() << "\t"
					<< af->getAonehomoPM() << "\t"
					<< af->getAonehomoPM_exp() << "\t"
					<< af->getAtwohomoPM() << "\t"
					<< af->getAtwohomoPM_exp() << "\t"
					<< af->getHetPM() << "\t"
					<< af->getHetPM_exp() << "\t"
					//parent female
					<< hw_PF << "\t"
					<< af->getPopPF() << "\t"
					<< af->getAonehomoPF() << "\t"
					<< af->getAonehomoPF_exp() << "\t"
					<< af->getAtwohomoPF() << "\t"
					<< af->getAtwohomoPF_exp() << "\t"
					<< af->getHetPF() << "\t"
					<< af->getHetPF_exp();
				}
				paren << endl;
			}
			if(options.doGender()){
				gend << (*markers)[i]->toString();
				if((*markers)[i]->isMicroSat()){
					for(int l = 0; l < 19; l++){
						gend << "\tNA";
					}
				}
				else{
					gend << "\t" << (*markers)[i]->getAllele1() << "_" << (*markers)[i]->getAllele1() << "\t"
					<< (*markers)[i]->getAllele1() << "_" << (*markers)[i]->getAllele2() << "\t"
					<< (*markers)[i]->getAllele2() << "_" << (*markers)[i]->getAllele2() << "\t";
					if(useoverall){
						//parent male
						gend << hw_OM << "\t"
						<< af->getPopM() << "\t"
						<< af->getAonehomoM() << "\t"
						<< af->getAonehomoM_exp() << "\t"
						<< af->getAtwohomoM() << "\t"
						<< af->getAtwohomoM_exp() << "\t"
						<< af->getHetM() << "\t"
						<< af->getHetM_exp() << "\t"
						//parent female
						<< hw_OF << "\t"
						<< af->getPopF() << "\t"
						<< af->getAonehomoF() << "\t"
						<< af->getAonehomoF_exp() << "\t"
						<< af->getAtwohomoF() << "\t"
						<< af->getAtwohomoF_exp() << "\t"
						<< af->getHetF() << "\t"
						<< af->getHetF_exp();
					}
					else{
						//parent male
						gend << hw_OM << "\t"
						<< af->getPopPM() << "\t"
						<< af->getAonehomoPM() << "\t"
						<< af->getAonehomoPM_exp() << "\t"
						<< af->getAtwohomoPM() << "\t"
						<< af->getAtwohomoPM_exp() << "\t"
						<< af->getHetPM() << "\t"
						<< af->getHetPM_exp() << "\t"
						//parent female
						<< hw_OF << "\t"
						<< af->getPopPF() << "\t"
						<< af->getAonehomoPF() << "\t"
						<< af->getAonehomoPF_exp() << "\t"
						<< af->getAtwohomoPF() << "\t"
						<< af->getAtwohomoPF_exp() << "\t"
						<< af->getHetPF() << "\t"
						<< af->getHetPF_exp();
					}
				}
				gend << endl;
			}
			if(options.doCaseControl()){
				cc << (*markers)[i]->toString();
				if((*markers)[i]->isMicroSat()){
					for(int l = 0; l < 35; l++){
						cc << "\tNA";
					}
				}
				else{
					cc << "\t" << (*markers)[i]->getAllele1() << "_" << (*markers)[i]->getAllele1() << "\t"
					<< (*markers)[i]->getAllele1() << "_" << (*markers)[i]->getAllele2() << "\t"
					<< (*markers)[i]->getAllele2() << "_" << (*markers)[i]->getAllele2() << "\t"
					<< hw_CaM << "\t"
					<< af->getPopCaM() << "\t"
					<< af->getAonehomoCaM() << "\t"
					<< af->getAonehomoCaM_exp() << "\t"
					<< af->getAtwohomoCaM() << "\t"
					<< af->getAtwohomoCaM_exp() << "\t"
					<< af->getHetCaM() << "\t"
					<< af->getHetCaM_exp() << "\t"

					<< hw_CaF << "\t"
					<< af->getPopCaF() << "\t"
					<< af->getAonehomoCaF() << "\t"
					<< af->getAonehomoCaF_exp() << "\t"
					<< af->getAtwohomoCaF() << "\t"
					<< af->getAtwohomoCaF_exp() << "\t"
					<< af->getHetCaF() << "\t"
					<< af->getHetCaF_exp() << "\t"

					<< hw_ConM << "\t"
					<< af->getPopConM() << "\t"
					<< af->getAonehomoConM() << "\t"
					<< af->getAonehomoConM_exp() << "\t"
					<< af->getAtwohomoConM() << "\t"
					<< af->getAtwohomoConM_exp() << "\t"
					<< af->getHetConM() << "\t"
					<< af->getHetConM_exp() << "\t"

					<< hw_ConF << "\t"
					<< af->getPopConF() << "\t"
					<< af->getAonehomoConF() << "\t"
					<< af->getAonehomoConF_exp() << "\t"
					<< af->getAtwohomoConF() << "\t"
					<< af->getAtwohomoConF_exp() << "\t"
					<< af->getHetConF() << "\t"
					<< af->getHetConF_exp();
				}
				cc << endl;
			}

			doFilter((*markers)[i]);
		}
	}

	if(pvals.is_open()){
		pvals.close();
	}

}

/*
 * Function: calcHW_exact
 * Description:
 * Performs the exact HW calculation
 */
float HWEquilibrium::calcHW_exact(int obs_hets, int obs_hom1, int obs_hom2){
/*
 *   This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
 *        Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
 *             Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
 *
 *                 Written by Jan Wigginton
 *                     */
    if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
    {
		opts::printLog("FATAL ERROR - SNP-HWE: Current genotype configuration (" + getString<int>(obs_hets) + " " + getString<int>(obs_hom1) + " " + getString<int>(obs_hom2) + " ) includes a negative count\n");
        //exit(1);
		throw MethodException("FATAL ERROR - SNP-HWE: Current genotype configuration (" + getString<int>(obs_hets) + " " + getString<int>(obs_hom1) + " " + getString<int>(obs_hom2) + " ) includes a negative count\n");
    }

    int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
    int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

    int rare_copies = 2 * obs_homr + obs_hets;
    int genotypes   = obs_hets + obs_homc + obs_homr;
    if(genotypes <= 0){
        return -1;
    }

    double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
    if (het_probs == NULL)
    {
		opts::printLog("FATAL ERROR - SNP-HWE: Current genotype configuration (" + getString<int>(obs_hets) + " " + getString<int>(obs_hom1) + " " + getString<int>(obs_hom2) + " ) includes a negative count\n");
        //exit(1);
		throw MethodException("FATAL ERROR - SNP-HWE: Current genotype configuration (" + getString<int>(obs_hets) + " " + getString<int>(obs_hom1) + " " + getString<int>(obs_hom2) + " ) includes a negative count\n");
    }

    int i;
    for (i = 0; i <= rare_copies; i++)
        het_probs[i] = 0.0;
       /* start at midpoint */
        int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
        /* check to ensure that midpoint and rare alleles have same parity */
        if ((rare_copies & 1) ^ (mid & 1))
            mid++;

        int curr_hets = mid;
        int curr_homr = (rare_copies - mid) / 2;
        int curr_homc = genotypes - curr_hets - curr_homr;

        het_probs[mid] = 1.0;
        double sum = het_probs[mid];
        for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
        {
            het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
            sum += het_probs[curr_hets - 2];

           /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
            curr_homr++;
            curr_homc++;
        }

        curr_hets = mid;
        curr_homr = (rare_copies - mid) / 2;
        curr_homc = genotypes - curr_hets - curr_homr;
        for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
        {
            het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc /((curr_hets + 2.0) * (curr_hets + 1.0));
            sum += het_probs[curr_hets + 2];

           /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
            curr_homr--;
            curr_homc--;
        }

        for (i = 0; i <= rare_copies; i++)
            het_probs[i] /= sum;

        double p_hwe = 0.0;
        /*  p-value calculation for p_hwe  */
        for (i = 0; i <= rare_copies; i++)
        {
            if (het_probs[i] > het_probs[obs_hets])
               continue;
            p_hwe += het_probs[i];
        }
        p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;
        free(het_probs);
     return (float) p_hwe;
}

/*
 * Function: calcHW
 * Description:
 * Performs classic HW calculation
 */
float HWEquilibrium::calcHW(float aonefreq, float atwofreq, int aoneobs, int hetobs, int atwoobs, int population){
	bool DEBUG = false;
	if(DEBUG){
        cout << "Initial:\t" << aonefreq << "\t" << atwofreq << "\t" << aoneobs << "\t" << atwoobs << "\t" << hetobs << "\t" << population << endl;
    }
	DEBUG = false;
    float p2 = pow(aonefreq, 2);
    float q2 = pow(atwofreq,2);
    float pq = 2 * aonefreq * atwofreq;
    if(DEBUG){
        cout << "P2,Q2,PQ:\t" << p2 << "\t" << q2 << "\t" << pq << endl;
    }
    float expp2 = p2 * (float)population;
    float expq2 = q2 * (float)population;
    float exppq = pq * (float)population;
    if(expq2 < 5 || expp2 < 5 || exppq < 5){
        return calcHW_exact(hetobs, aoneobs, atwoobs);
    }
    if(DEBUG){
        cout << "Exp:\t" << expp2 << "\t" << expq2 << "\t" << exppq << endl;
    }
    float chi1 = 0.0;
    float chi2 = 0.0;
    float chi3 = 0.0;

        chi1 = (pow(((float)aoneobs - expp2), 2)) / expp2;
        chi2 = (pow(((float)atwoobs - expq2),2)) / expq2;
        chi3 = (pow(((float)hetobs - exppq),2)) / exppq;


    double chi = chi1 + chi2 + chi3;

    if(DEBUG){
        cout << "ChiTot:\t" << chi1 << "\t" << chi2 << "\t" << chi3 << "\t=" << chi << endl;
    }
	double results, df = 1;
    results = Helpers::p_from_chi(chi, df);
	if(DEBUG){
      cout << "Results: " << results << endl;
    }
    return ((float)results);
}

}
