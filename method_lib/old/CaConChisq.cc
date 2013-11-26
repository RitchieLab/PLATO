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
*File: CaConChisq.cc
**********************************************************************************/


#include <iostream>
#include <fstream>
#include "CaConChisq.h"
#include "ChiSquareAllelic.h"
#include "ChiSquareArmitage.h"
#include "FisherExact.h"
#include "cdflib.h"
#include "General.h"
#include "Options.h"
#include "Helpers.h"
namespace Methods{
string CaConChisq::stepname = "chisquare";

/*DEPRECATED
 *
 * Function: PrintSummary
 * Description:
 * Outputs results of chisquare step
 * Resets Marker flags
 */
void CaConChisq::PrintSummary(){
	string pfname = opts::_OUTPREFIX_ + "chisquare" + options.getOut() + ".txt";
	if(!overwrite){
		pfname += "." + getString<int>(order);
	}
	ofstream pvals (pfname.c_str());
	if(!pvals){
		opts::printLog("Error opening " + pfname + " for writing results! Exiting!\n");
		throw MethodException("Error opening " + pfname + " for writing results! Exiting!\n");
	}
	opts::addFile("Marker", stepname, pfname);
	pvals.precision(4);



    pvals << "Chrom"
		  << "\trsID"
		  << "\tProbeID"
		  << "\tbploc";
	if((*markers).at(0)->getDetailHeaders().size() > 0){
		pvals << "\t" << (*markers).at(0)->getDetailHeaders();
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

	int msize = markers->size();


	for(int i = 0; i < msize; i++){
		if((*markers).at(i)->isEnabled() && !(*markers).at(i)->isFlagged()){
			if((*markers).at(i)->isMicroSat()){
				pvals << (*markers).at(i)->toString()
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

			pvals << (*markers).at(i)->toString() << "\t"
				<< chi_arm.at(i) << "\t"
				<< pval_arm.at(i) << "\t"
				<< chi_allele.at(i) << "\t"
				<< pval_allele.at(i) << "\t"
				<< pval_allele_exact.at(i) << "\t"
				<< chi_geno.at(i) << "\t"
				<< pval_geno.at(i) << "\t"
				<< pval_geno_exact.at(i) << "\t"
				<< (*markers).at(i)->getAllele1() << ":" << (*markers).at(i)->getAllele2() << "\t"
				<< odds_ratio.at(i) << "\t"
				<< (*markers).at(i)->getAllele1() << "\t"
				<< ci_l.at(i) << "\t"
				<< ci_u.at(i)
				<< endl;

			if(options.doGroupFile()){
			}
			(*markers).at(i)->setFlag(false);
		}
	}

	if(pvals.is_open()){
		pvals.close();
	}
}

/*DEPRECATED
 *
 * Function: FilterSummary
 * Description:
 * Outputs remaining marker count
 */
void CaConChisq::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>((opts::_MARKERS_WORKING_ - orig_num_markers)) + " (" +
	        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
	        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

	resize(0);
	if(af != NULL){
		delete(af);
	}
}

/*
 * Function: resize
 * Description:
 * Resizes result vectors to specified count
 */
void CaConChisq::resize(int i){
	chi_geno.resize(i);
	chi_allele.resize(i);
	pval_geno.resize(i);
	pval_allele.resize(i);
	ci_l.resize(i);
	ci_u.resize(i);
}

/*DEPRECATED
 *
 * Function: filter
 * Description:
 * Filters markers based on the armitage pvalue
 */
void CaConChisq::filter(){

	//false if out in P or C
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int msize = markers->size();
		for(int i = 0; i < msize; i++){
			if((*markers).at(i)->isEnabled() && !(*markers).at(i)->isFlagged()){
				bool inc = false;
				if(options.doThreshMarkersLow() && Helpers::dLess(pval_arm.at(i), options.getThreshMarkersLow())){
					(*markers).at(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && Helpers::dGreater(pval_arm.at(i), options.getThreshMarkersHigh())){
					(*markers).at(i)->setEnabled(false);
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
 * performs chisquare calculation on one marker
 * calculates Armitage, Allelic, Allelic Exact, and Genotypic
 * Calculates by group as well if option is specified
 */
void CaConChisq::calculate(Marker* mark){
	int ssize = samples->size();

	unsigned int cases = 0;
	unsigned int controls = 0;
	vector<vector<int> > chitotals;
	chitotals.push_back(vector<int>(0));
	chitotals.push_back(vector<int>(0));
	FisherExact fish;
	ChiSquareAllelic chi;
	ChiSquareArmitage arm;
	double zt = Helpers::ltqnorm(1 - (1 - options.getCI()) / 2);

	gchi_arm_one.clear();
	gpval_arm_one.clear();
	gchi_allele_one.clear();
	gpval_allele_one.clear();
	gpval_allele_exact_one.clear();
	gchi_geno_one.clear();
	gpval_geno_one.clear();
	gpval_geno_exact_one.clear();
	godds_ratio_one.clear();
	gci_l_one.clear();
	gci_u_one.clear();

	AlleleFrequency gaf(samples, families);
	gaf.setRank(rank);
	gaf.setOptions(options);

	if(options.doGroupFile()){
		map<string, vector<Sample*> > groups = options.getGroups();
		map<string, vector<Sample*> >::iterator giter;
		for(giter = groups.begin(); giter != groups.end(); giter++){
			string mygroup = giter->first;
			vector<Sample*> mysamps = giter->second;
			map<string, vector<Family*> > gfams = options.getGroupFamilies();
			vector<Family*> myfams = gfams[mygroup];
			gaf.resetDataSet(&mysamps, &myfams);

			gaf.calcOne(mark);
			cases = 0;
			controls = 0;
			for(int s = 0; s < (int)mysamps.size(); s++){
				Sample* samp = mysamps.at(s);
				if(samp->getPheno() == 2){
					cases++;
				}
				else if(samp->getPheno() == 1){
					controls++;
				}
			}
			gchi_geno_one[mygroup] = calcChiGeno(gaf.getAonehomoCa(), gaf.getHetCa(), gaf.getAtwohomoCa(), gaf.getAonehomoCon(), gaf.getHetCon(), gaf.getAtwohomoCon());
			gpval_geno_one[mygroup] = -1;
			if(gchi_geno_one[mygroup] >= 0){
	             double pvalue, df = 2;
				 if((gaf.getAonehomoCa() == 0 && gaf.getAonehomoCon() == 0) ||
						 (gaf.getHetCa() == 0 && gaf.getHetCon() == 0) ||
						 (gaf.getAtwohomoCa() == 0 && gaf.getAtwohomoCon() == 0)){
					 df = 1;
				 }
	             pvalue = -1;
	             if(gchi_geno_one[mygroup] > -1 && !isnan(gchi_geno_one[mygroup]) && !isinf(gchi_geno_one[mygroup])){
	            	 pvalue = Helpers::p_from_chi(gchi_geno_one[mygroup], df);

	             }

				gpval_geno_one[mygroup] = pvalue;
			}

			gpval_geno_exact_one[mygroup] = -1;
			chitotals.at(0).push_back(gaf.getAoneCa_count());
			chitotals.at(0).push_back(gaf.getAtwoCa_count());
			chitotals.at(1).push_back(gaf.getAoneCon_count());
			chitotals.at(1).push_back(gaf.getAtwoCon_count());
			gchi_allele_one[mygroup] = chi.chisquare(chitotals);
			gpval_allele_one[mygroup] = -1;
			if(gchi_allele_one[mygroup] > -1){
            	double pvalue, df = 1;
             	pvalue = -1;
             	if(gchi_allele_one[mygroup] > -1 && !isnan(gchi_allele_one[mygroup]) && !isinf(gchi_allele_one[mygroup])){
             		pvalue = Helpers::p_from_chi(gchi_allele_one[mygroup], df);
             	}
				gpval_allele_one[mygroup] = pvalue;
			}
			gpval_allele_exact_one[mygroup] = fish.fisher_2_2(gaf.getAoneCa_count(), gaf.getAtwoCa_count(), gaf.getAoneCon_count(), gaf.getAtwoCon_count());
			if(gaf.getAtwoCa_count() == 0 ||  gaf.getAoneCon_count() == 0){
				godds_ratio_one[mygroup] = -1;
			}
			else{
				godds_ratio_one[mygroup] = (double)((double)gaf.getAoneCa_count() * (double)gaf.getAtwoCon_count()) / (double)((double)gaf.getAtwoCa_count() * (double)gaf.getAoneCon_count());
			}
			chitotals.at(0).clear();
			chitotals.at(1).clear();

			chitotals.at(0).push_back(gaf.getAonehomoCon());
			chitotals.at(0).push_back(gaf.getHetCon());
			chitotals.at(0).push_back(gaf.getAtwohomoCon());
			chitotals.at(1).push_back(gaf.getAonehomoCa());
			chitotals.at(1).push_back(gaf.getHetCa());
			chitotals.at(1).push_back(gaf.getAtwohomoCa());
			gchi_arm_one[mygroup] = arm.armitage(chitotals);
			gpval_arm_one[mygroup] = -1;
			if(gchi_arm_one[mygroup] >= 0){
	            double pvalue, df = 1;
           		pvalue = -1;
           		if(gchi_arm_one[mygroup] > -1 && !isnan(gchi_arm_one[mygroup]) && !isinf(gchi_arm_one[mygroup])){
           			pvalue = Helpers::p_from_chi(gchi_arm_one[mygroup], df);
           		}

				gpval_arm_one[mygroup] = pvalue;
			}
			chitotals.at(0).clear();
			chitotals.at(1).clear();
			double lOR = log(godds_ratio_one[mygroup]);
			double SE = sqrt(1/(double)gaf.getAoneCa_count() + 1/(double)gaf.getAtwoCa_count() + 1/(double)gaf.getAoneCon_count() + 1/(double)gaf.getAtwoCon_count());
			gci_l_one[mygroup] = exp(lOR - options.getCI() * SE);
			gci_u_one[mygroup] = exp(lOR + options.getCI() * SE);
		}
	}


	cases = 0;
	controls = 0;
	for(int i = 0; i < ssize; i++){
		if((*samples).at(i)->isEnabled()){
			if((*samples).at(i)->getPheno() == 2){
				cases++;
			}
			else if((*samples).at(i)->getPheno() == 1){
				controls++;
			}
		}
	}
	chi_geno_one = -1;
	chi_allele_one = -1;
	chi_arm_one = -1;
	pval_arm_one = -1;
	pval_allele_one = -1;
	pval_geno_exact_one = -1;
	pval_geno_one = -1;
	pval_allele_exact_one = -1;
	odds_ratio_one = -1;
	ci_l_one = -1;
	ci_u_one = -1;

	if(mark->isMicroSat()){
		return;
	}
	af->calcOne(mark);

	chi_geno_one = calcChiGeno(af->getAonehomoCa(), af->getHetCa(), af->getAtwohomoCa(), af->getAonehomoCon(), af->getHetCon(), af->getAtwohomoCon());
	pval_geno_one = -1;
	if(chi_geno_one >= 0){
    	double pvalue, df = 2;
		if((af->getAonehomoCa() == 0 && af->getAonehomoCon() == 0) ||
				(af->getHetCa() == 0 && af->getHetCon() == 0) ||
				(af->getAtwohomoCa() == 0 && af->getAtwohomoCon() == 0)){
			df = 1;
		}
	    pvalue = -1;
	    if(chi_geno_one > -1 && !isnan(chi_geno_one) && !isinf(chi_geno_one)){
	    	pvalue = Helpers::p_from_chi(chi_geno_one, df);
	    }
	    geno_df = df;
		pval_geno_one = pvalue;
	}

	pval_geno_exact_one = -1;
	chitotals.at(0).push_back(af->getAoneCa_count());
	chitotals.at(0).push_back(af->getAtwoCa_count());
	chitotals.at(1).push_back(af->getAoneCon_count());
	chitotals.at(1).push_back(af->getAtwoCon_count());
	chi_allele_one = chi.chisquare(chitotals);
	pval_allele_one = -1;
	if(chi_allele_one > -1){
	    double pvalue, df = 1;
	    pvalue = -1;
		if(chi_allele_one > -1 && !isnan(chi_allele_one) && !isinf(chi_allele_one)){
			pvalue = Helpers::p_from_chi(chi_allele_one, df);
		}
		allele_df = df;
		pval_allele_one = pvalue;
	}
	pval_allele_exact_one = fish.fisher_2_2(af->getAoneCa_count(), af->getAtwoCa_count(), af->getAoneCon_count(), af->getAtwoCon_count());
	allele_exact_df = 1;
	if(af->getAtwoCa_count() == 0 ||  af->getAoneCon_count() == 0){
		odds_ratio_one = -1;
		ci_l_one = -999;
		ci_u_one = -999;
	}
	else{
		odds_ratio_one = (double)((double)af->getAoneCa_count() * (double)af->getAtwoCon_count()) / (double)((double)af->getAtwoCa_count() * (double)af->getAoneCon_count());
		double lOR = log(odds_ratio_one);
		double SE = sqrt(1/(double)af->getAoneCa_count() + 1/(double)af->getAtwoCa_count() + 1/(double)af->getAoneCon_count() + 1/(double)af->getAtwoCon_count());
		ci_l_one = exp(lOR - zt * SE);
		ci_u_one = exp(lOR + zt * SE);
	}
	//		}
	chitotals.at(0).clear();
	chitotals.at(0).clear();

	chitotals.at(0).push_back(af->getAonehomoCon());
	chitotals.at(0).push_back(af->getHetCon());
	chitotals.at(0).push_back(af->getAtwohomoCon());
	chitotals.at(1).push_back(af->getAonehomoCa());
	chitotals.at(1).push_back(af->getHetCa());
	chitotals.at(1).push_back(af->getAtwohomoCa());
	chi_arm_one = arm.armitage(chitotals);
	pval_arm_one = -1;
	if(chi_arm_one >= 0){
		double pvalue, df = 1;
	    pvalue = -1;
	    if(chi_arm_one > -1 && !isnan(chi_arm_one) && !isinf(chi_arm_one)){
	    	pvalue = Helpers::p_from_chi(chi_arm_one, df);
	    }
	    arm_df = df;
		pval_arm_one = pvalue;
	}
	chitotals.at(0).clear();
	chitotals.at(1).clear();

}

/*DEPRECATED
 * see calculate
 *
 * Function: process
 * Description:
 * Main step processing function
 *
 */
void CaConChisq::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

	af = new AlleleFrequency(samples, families);

	af->setRank(rank);
	if(options.doGroupFile()){
		options.readGroups(samples);
	}
	af->setOptions(options);

	int msize = markers->size();
	int ssize = samples->size();

	unsigned int cases = 0;
	unsigned int controls = 0;
	vector<vector<int> > chitotals;
	chitotals.push_back(vector<int>(0));
	chitotals.push_back(vector<int>(0));
	FisherExact fish;
	ChiSquareAllelic chi;
	ChiSquareArmitage arm;
	double zt = Helpers::ltqnorm(1 - (1 - options.getCI()) / 2);

	if(options.doGroupFile()){
		string gpfname = opts::_OUTPREFIX_ + "chisquare_groups" + options.getOut() + ".txt";
		ofstream gpvals;
		gpvals.open(gpfname.c_str());
		if(!gpvals){
			opts::printLog("Error opening " + gpfname + " for writing results! Exiting!\n");
			throw MethodException("Error opening " + gpfname + " for writing results! Exiting!\n");
		}
		opts::addFile("Marker", stepname, gpfname);
		gpvals.precision(4);
		gpvals << "Chrom\trsID\tProbeID\tbploc";
		if((*markers).at(0)->getDetailHeaders().size() > 0){
			gpvals << "\t" << (*markers).at(0)->getDetailHeaders();
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
		for(int i = 0; i < msize; i++){
			if((*markers).at(i)->isEnabled() && Helpers::isValidMarker((*markers).at(i), &options, prev_base, prev_chrom)){
				gpvals << (*markers).at(i)->toString();
				if((*markers).at(i)->isMicroSat()){
					for(giter = groups.begin(); giter != groups.end(); giter++){
						gpvals << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
					}
					gpvals << endl;
					continue;
				}
				for(giter = groups.begin(); giter != groups.end(); giter++){
					string mygroup = giter->first;
					vector<Sample*> mysamps = giter->second;
					map<string, vector<Family*> > gfams = options.getGroupFamilies();
					vector<Family*> myfams = gfams[mygroup];
					AlleleFrequency* gaf = new AlleleFrequency(&mysamps, &myfams);
					gaf->setRank(rank);
					gaf->setOptions(options);

					gaf->calcOne((*markers).at(i));
					cases = 0;
					controls = 0;
					for(int s = 0; s < (int)mysamps.size(); s++){
						Sample* samp = mysamps.at(s);
						if(samp->getPheno() == 2){
							cases++;
						}
						else if(samp->getPheno() == 1){
							controls++;
						}
					}
					gchi_geno = calcChiGeno(gaf->getAonehomoCa(), gaf->getHetCa(), gaf->getAtwohomoCa(), gaf->getAonehomoCon(), gaf->getHetCon(), gaf->getAtwohomoCon());
					gpval_geno = -1;
					if(gchi_geno >= 0){
			             double pvalue, df = 2;
			             pvalue = -1;
			             if(gchi_geno > -1){
			            	 pvalue = Helpers::p_from_chi(gchi_geno, df);
			             }

						gpval_geno = pvalue;
					}

					gpval_geno_exact = -1;
					chitotals.at(0).push_back(gaf->getAoneCa_count());
					chitotals.at(0).push_back(gaf->getAtwoCa_count());
					chitotals.at(1).push_back(gaf->getAoneCon_count());
					chitotals.at(1).push_back(gaf->getAtwoCon_count());
					gchi_allele = chi.chisquare(chitotals);
					gpval_allele = -1;
					if(gchi_allele > -1){
		            	double pvalue, df = 1;
		             	pvalue = -1;
		             	if(gchi_allele > -1){
		             		pvalue = Helpers::p_from_chi(gchi_allele, df);
		             	}
						gpval_allele = pvalue;
					}
					gpval_allele_exact = fish.fisher_2_2(gaf->getAoneCa_count(), gaf->getAtwoCa_count(), gaf->getAoneCon_count(), gaf->getAtwoCon_count());
					if(gaf->getAtwoCa_count() == 0 ||  gaf->getAoneCon_count() == 0){
						godds_ratio = -1;
					}
					else{
						godds_ratio = (double)((double)gaf->getAoneCa_count() * (double)gaf->getAtwoCon_count()) / (double)((double)gaf->getAtwoCa_count() * (double)gaf->getAoneCon_count());
					}
					chitotals.at(0).clear();
					chitotals.at(1).clear();

					chitotals.at(0).push_back(gaf->getAonehomoCon());
					chitotals.at(0).push_back(gaf->getHetCon());
					chitotals.at(0).push_back(gaf->getAtwohomoCon());
					chitotals.at(1).push_back(gaf->getAonehomoCa());
					chitotals.at(1).push_back(gaf->getHetCa());
					chitotals.at(1).push_back(gaf->getAtwohomoCa());
					gchi_arm = arm.armitage(chitotals);
					gpval_arm = -1;
					if(gchi_arm >= 0){
			            double pvalue, df = 1;
	             		pvalue = -1;
	             		if(gchi_arm > -1){
	             			pvalue = Helpers::p_from_chi(gchi_arm, df);
	             		}

						gpval_arm = pvalue;
					}
					chitotals.at(0).clear();
					chitotals.at(1).clear();
					double lOR = log(godds_ratio);
					double SE = sqrt(1/(double)gaf->getAoneCa_count() + 1/(double)gaf->getAtwoCa_count() + 1/(double)gaf->getAoneCon_count() + 1/(double)gaf->getAtwoCon_count());
					double OR_lower = exp( lOR - options.getCI() * SE);
					double OR_upper = exp(lOR + options.getCI() * SE);

					gpvals << "\t" << gchi_arm << "\t"
						<< gpval_arm << "\t"
						<< gchi_allele << "\t"
						<< gpval_allele << "\t"
						<< gpval_allele_exact << "\t"
						<< gchi_geno << "\t"
						<< gpval_geno << "\t"
						<< gpval_geno_exact << "\t"
						<< (*markers).at(i)->getAllele1() << ":" << (*markers).at(i)->getAllele2() << "\t"
						<< godds_ratio << "\t"
						<< (*markers).at(i)->getAllele1() << "\t"
						<< OR_lower << "\t"
						<< OR_upper;
					gchi_arm = -1;
					gpval_arm = -1;
					gchi_allele = -1;
					gpval_allele = -1;
					gpval_allele_exact = -1;
					gchi_geno = -1;
					gpval_geno = -1;
					gpval_geno_exact = -1;
					godds_ratio = -1;
					delete gaf;
				}
				gpvals << endl;
			}

		}
	}

	cases = 0;
	controls = 0;
	for(int i = 0; i < ssize; i++){
		if((*samples).at(i)->isEnabled()){
			if((*samples).at(i)->getPheno() == 2){
				cases++;
			}
			else if((*samples).at(i)->getPheno() == 1){
				controls++;
			}
		}
	}
	chi_geno.resize(msize);
	chi_allele.resize(msize);
	chi_arm.resize(msize);
	pval_arm.resize(msize);
	pval_allele.resize(msize);
	pval_geno_exact.resize(msize);
	pval_geno.resize(msize);
	pval_allele_exact.resize(msize);
	odds_ratio.resize(msize);
	ci_l.resize(msize);
	ci_u.resize(msize);

	int prev_base = 0;
	int prev_chrom = -1;
	for(int i = 0; i < msize; i++){
		if((*markers).at(i)->isEnabled() && Helpers::isValidMarker((*markers).at(i), &options, prev_base, prev_chrom)){
			if((*markers).at(i)->isMicroSat()){
				continue;
			}
			af->calcOne((*markers).at(i));
				chi_geno.at(i) = calcChiGeno(af->getAonehomoCa(), af->getHetCa(), af->getAtwohomoCa(), af->getAonehomoCon(), af->getHetCon(), af->getAtwohomoCon());
				pval_geno.at(i) = -1;
				if(chi_geno.at(i) >= 0){
		             double pvalue, df = 2;
		             pvalue = -1;
		             if(chi_geno.at(i) > -1){
		            	 pvalue = Helpers::p_from_chi(chi_geno.at(i), df);
		             }

					pval_geno.at(i) = pvalue;
				}

				pval_geno_exact.at(i) = -1;
				chitotals.at(0).push_back(af->getAoneCa_count());
				chitotals.at(0).push_back(af->getAtwoCa_count());
				chitotals.at(1).push_back(af->getAoneCon_count());
				chitotals.at(1).push_back(af->getAtwoCon_count());
				chi_allele.at(i) = chi.chisquare(chitotals);
				pval_allele.at(i) = -1;
				if(chi_allele.at(i) > -1){
		             double pvalue,df = 1;
		             pvalue = -1;
		             if(chi_allele.at(i) > -1){
		            	 pvalue = Helpers::p_from_chi(chi_allele.at(i), df);
		             }
					pval_allele.at(i) = pvalue;
				}
				pval_allele_exact.at(i) = fish.fisher_2_2(af->getAoneCa_count(), af->getAtwoCa_count(), af->getAoneCon_count(), af->getAtwoCon_count());
				if(af->getAtwoCa_count() == 0 ||  af->getAoneCon_count() == 0){
					odds_ratio.at(i) = -1;
					ci_l.at(i) = -999;
					ci_u.at(i) = -999;
				}
				else{
					odds_ratio.at(i) = (double)((double)af->getAoneCa_count() * (double)af->getAtwoCon_count()) / (double)((double)af->getAtwoCa_count() * (double)af->getAoneCon_count());
					double lOR = log(odds_ratio.at(i));
					double SE = sqrt(1/(double)af->getAoneCa_count() + 1/(double)af->getAtwoCa_count() + 1/(double)af->getAoneCon_count() + 1/(double)af->getAtwoCon_count());
					ci_l.at(i) = exp( lOR - zt * SE);
					ci_u.at(i) = exp(lOR + zt * SE);
				}
			chitotals.at(0).clear();
			chitotals.at(1).clear();

			chitotals.at(0).push_back(af->getAonehomoCon());
			chitotals.at(0).push_back(af->getHetCon());
			chitotals.at(0).push_back(af->getAtwohomoCon());
			chitotals.at(1).push_back(af->getAonehomoCa());
			chitotals.at(1).push_back(af->getHetCa());
			chitotals.at(1).push_back(af->getAtwohomoCa());
			chi_arm.at(i) = arm.armitage(chitotals);
			pval_arm.at(i) = -1;
			if(chi_arm.at(i) >= 0){
	             double pvalue, df = 1;
	             pvalue = -1;
	             if(chi_arm.at(i) > -1){
	            	 pvalue = Helpers::p_from_chi(chi_arm.at(i), df);
	             }

				pval_arm.at(i) = pvalue;
			}
			chitotals.at(0).clear();
			chitotals.at(1).clear();
		}
	}

	if(af){
		delete(af);
		af = NULL;
	}
}

/*
 * calculate chi value from genotypic information
 */
double CaConChisq::calcChiGeno(int cahomo1, int cahet, int cahomo2, int conhomo1, int conhet, int conhomo2){
	double chi = 0;

	float cahomo1exp;
	float cahetexp;
	float cahomo2exp;
	float conhomo1exp;
	float conhetexp;
	float conhomo2exp;

	int homo1 = cahomo1 + conhomo1;
	int het = cahet + conhet;
	int homo2 = cahomo2 + conhomo2;
	int ca = cahomo1 + cahet + cahomo2;
	int con = conhomo1 + conhet + conhomo2;

	cahomo1exp = (float) (((float) homo1 * (float) ca) / (float)(ca + con));
	conhomo1exp = (float) (((float) homo1 * (float) con) / (float)(ca + con));
	cahetexp = (float) (((float) het * (float) ca) / (float)(ca + con));
	conhetexp = (float) (((float) het * (float) con) / (float)(ca + con));
	cahomo2exp = (float) (((float) homo2 * (float) ca) / (float)(ca + con));
	conhomo2exp = (float) (((float) homo2 * (float) con) / (float)(ca + con));

	if(cahomo1exp > 0){
		chi += (pow(((float)cahomo1 - cahomo1exp),2)) / cahomo1exp;
	}
	if(conhomo1exp > 0){
		chi += (pow(((float)conhomo1 - conhomo1exp),2)) / conhomo1exp;
	}
	if(cahetexp > 0){
		chi += (pow(((float)cahet - cahetexp),2)) / cahetexp;
	}
	if(conhetexp > 0){
		chi += (pow(((float)conhet - conhetexp),2)) / conhetexp;
	}
	if(cahomo2exp > 0){
		chi += (pow(((float)cahomo2 - cahomo2exp),2)) / cahomo2exp;
	}
	if(conhomo2exp > 0){
		chi += (pow(((float)conhomo2 - conhomo2exp),2)) / conhomo2exp;
	}

	return chi;
}

}
