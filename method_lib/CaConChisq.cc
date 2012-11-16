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
//#include "Chrom.h"
#include "ChiSquareAllelic.h"
#include "ChiSquareArmitage.h"
#include "FisherExact.h"
//#include "helper.h"
#include "cdflib.h"
#include "General.h"
#include "Options.h"
#include "Helper.h"
namespace Methods{
string CaConChisq::stepname = "chisquare";

/*
 * Function: PrintSummary
 * Description:
 * Outputs results of chisquare step
 * Resets Marker flags
 */
void CaConChisq::PrintSummary(){
	string pfname = opts::_OUTPREFIX_ + "chisquare" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		pfname += "." + getString<int>(order);
	}
	ofstream pvals (pfname.c_str());
	if(!pvals){
		opts::printLog("Error opening " + pfname + " for writing results! Exiting!\n");
		exit(1);
	}
	opts::addFile("Marker", stepname, pfname);
	pvals.precision(4);



    pvals << "Chrom"
		  << "\trsID"
		  << "\tProbeID"
		  << "\tbploc";
	if((*markers)[0]->getDetailHeaders().size() > 0){
		pvals << "\t" << (*markers)[0]->getDetailHeaders();
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
		if((*markers)[i]->isEnabled() && !(*markers)[i]->isFlagged()){
			if((*markers)[i]->isMicroSat()){
				pvals << (*markers)[i]->toString()
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

			pvals << (*markers)[i]->toString() << "\t"
				<< chi_arm[i] << "\t"
				<< pval_arm[i] << "\t"
				<< chi_allele[i] << "\t"
				<< pval_allele[i] << "\t"
				<< pval_allele_exact[i] << "\t"
				<< chi_geno[i] << "\t"
				<< pval_geno[i] << "\t"
				<< pval_geno_exact[i] << "\t"
				<< (*markers)[i]->getAllele1() << ":" << (*markers)[i]->getAllele2() << "\t"
				<< odds_ratio[i] << "\t"
				<< (*markers)[i]->getAllele1() << "\t"
				<< ci_l[i] << "\t"
				<< ci_u[i]
				<< endl;

			if(options.doGroupFile()){
			}
			(*markers)[i]->setFlag(false);
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
void CaConChisq::FilterSummary(){
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

/*
 * Function: filter
 * Description:
 * Filters markers based on the armitage pvalue
 */
void CaConChisq::filter(){

	//false if out in P or C
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int msize = markers->size();
		for(int i = 0; i < msize; i++){
			if((*markers)[i]->isEnabled() && !(*markers)[i]->isFlagged()){
				bool inc = false;
				if(options.doThreshMarkersLow() && dLess(pval_arm[i], options.getThreshMarkersLow())){
					(*markers)[i]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && dGreater(pval_arm[i], options.getThreshMarkersHigh())){
					(*markers)[i]->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}
		}
	}
}

void CaConChisq::calculate(Marker* mark){
	//int msize = markers->size();
	int ssize = samples->size();

	unsigned int cases = 0;
	unsigned int controls = 0;
	vector<vector<int> > chitotals;
	chitotals.push_back(vector<int>(0));
	chitotals.push_back(vector<int>(0));
	FisherExact fish;
	ChiSquareAllelic chi;
	ChiSquareArmitage arm;
	double zt = ltqnorm(1 - (1 - options.getCI()) / 2);

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

	if(options.doGroupFile()){
		map<string, vector<Sample*> > groups = options.getGroups();
		map<string, vector<Sample*> >::iterator giter;

		for(giter = groups.begin(); giter != groups.end(); giter++){
			string mygroup = giter->first;
			vector<Sample*> mysamps = giter->second;
			map<string, vector<Family*> > gfams = options.getGroupFamilies();
			vector<Family*> myfams = gfams[mygroup];//generateFamilySet(&mysamps);
			AlleleFrequency* gaf = new AlleleFrequency(&mysamps, &myfams);
			gaf->setRank(rank);
			gaf->setOptions(options);

			gaf->calcOne(mark);
			cases = 0;
			controls = 0;
			for(int s = 0; s < (int)mysamps.size(); s++){
				Sample* samp = mysamps[s];
				if(samp->getPheno() == 2){
					cases++;
				}
				else if(samp->getPheno() == 1){
					controls++;
				}
			}
			gchi_geno_one[mygroup] = calcChiGeno(gaf->getAonehomoCa(), gaf->getHetCa(), gaf->getAtwohomoCa(), gaf->getAonehomoCon(), gaf->getHetCon(), gaf->getAtwohomoCon());
			gpval_geno_one[mygroup] = -1;
			if(gchi_geno_one[mygroup] >= 0){
	             double pvalue, p, bound, df = 2;
				 if((gaf->getAonehomoCa() == 0 && gaf->getAonehomoCon() == 0) ||
						 (gaf->getHetCa() == 0 && gaf->getHetCon() == 0) ||
						 (gaf->getAtwohomoCa() == 0 && gaf->getAtwohomoCon() == 0)){
					 df = 1;
				 }
	             pvalue = -1;
	             int code = 1, status;
	             if(gchi_geno_one[mygroup] > -1){
	                 cdfchi(&code, &p, &pvalue, &(gchi_geno_one[mygroup]), &df, &status, &bound);
	             }

				gpval_geno_one[mygroup] = pvalue;
			}

			gpval_geno_exact_one[mygroup] = -1;
			chitotals[0].push_back(gaf->getAoneCa_count());
			chitotals[0].push_back(gaf->getAtwoCa_count());
			chitotals[1].push_back(gaf->getAoneCon_count());
			chitotals[1].push_back(gaf->getAtwoCon_count());
			gchi_allele_one[mygroup] = chi.chisquare(chitotals);
			gpval_allele_one[mygroup] = -1;
			if(gchi_allele_one[mygroup] > -1){
            	double pvalue, p, bound, df = 1;
             	pvalue = -1;
             	int code = 1, status;
             	if(gchi_allele_one[mygroup] > -1){
                	cdfchi(&code, &p, &pvalue, &(gchi_allele_one[mygroup]), &df, &status, &bound);
             	}
				gpval_allele_one[mygroup] = pvalue;
			}
			gpval_allele_exact_one[mygroup] = fish.fisher_2_2(gaf->getAoneCa_count(), gaf->getAtwoCa_count(), gaf->getAoneCon_count(), gaf->getAtwoCon_count());
			if(gaf->getAtwoCa_count() == 0 ||  gaf->getAoneCon_count() == 0){
				godds_ratio_one[mygroup] = -1;
			}
			else{
				godds_ratio_one[mygroup] = (double)((double)gaf->getAoneCa_count() * (double)gaf->getAtwoCon_count()) / (double)((double)gaf->getAtwoCa_count() * (double)gaf->getAoneCon_count());
			}
			chitotals[0].clear();
			chitotals[1].clear();

			chitotals[0].push_back(gaf->getAonehomoCon());
			chitotals[0].push_back(gaf->getHetCon());
			chitotals[0].push_back(gaf->getAtwohomoCon());
			chitotals[1].push_back(gaf->getAonehomoCa());
			chitotals[1].push_back(gaf->getHetCa());
			chitotals[1].push_back(gaf->getAtwohomoCa());
			gchi_arm_one[mygroup] = arm.armitage(chitotals);
			gpval_arm_one[mygroup] = -1;
			if(gchi_arm_one[mygroup] >= 0){
	            double pvalue, p, bound, df = 1;
           		pvalue = -1;
           		int code = 1, status;
           		if(gchi_arm_one[mygroup] > -1){
               		cdfchi(&code, &p, &pvalue, &(gchi_arm_one[mygroup]), &df, &status, &bound);
           		}

				gpval_arm_one[mygroup] = pvalue;
			}
			chitotals[0].clear();
			chitotals[1].clear();
			double lOR = log(godds_ratio_one[mygroup]);
			double SE = sqrt(1/(double)gaf->getAoneCa_count() + 1/(double)gaf->getAtwoCa_count() + 1/(double)gaf->getAoneCon_count() + 1/(double)gaf->getAtwoCon_count());
			gci_l_one[mygroup] = exp(lOR - options.getCI() * SE);
			gci_u_one[mygroup] = exp(lOR + options.getCI() * SE);

			delete gaf;
		}
	}


	cases = 0;
	controls = 0;
	for(int i = 0; i < ssize; i++){
		if((*samples)[i]->isEnabled()){
			if((*samples)[i]->getPheno() == 2){
				cases++;
			}
			else if((*samples)[i]->getPheno() == 1){
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
    	double pvalue, p, bound, df = 2;
		if((af->getAonehomoCa() == 0 && af->getAonehomoCon() == 0) ||
				(af->getHetCa() == 0 && af->getHetCon() == 0) ||
				(af->getAtwohomoCa() == 0 && af->getAtwohomoCon() == 0)){
			df = 1;
		}
	    pvalue = -1;
	    int code = 1, status;
	    if(chi_geno_one > -1){
	    	cdfchi(&code, &p, &pvalue, &chi_geno_one, &df, &status, &bound);
		}
	    geno_df = df;
		pval_geno_one = pvalue;
	}

	//		}
	//		else{
				//chi_geno[i] = -1;
	pval_geno_exact_one = -1;
	//		}
	//		if(af->getAoneCon_count(i) > 4 && af->getAtwoCon_count(i) > 4 && af->getAoneCa_count(i) > 4 && af->getAtwoCa_count(i) > 4){
				//cout << (*markers)[i]->getProbeID() << "\t" << af->getAoneCa_count(i) << "\t" << af->getAtwoCa_count(i) << "\t" << af->getAoneCon_count(i) << "\t" << af->getAtwoCon_count(i) << endl;
	chitotals[0].push_back(af->getAoneCa_count());
	chitotals[0].push_back(af->getAtwoCa_count());
	chitotals[1].push_back(af->getAoneCon_count());
	chitotals[1].push_back(af->getAtwoCon_count());
	chi_allele_one = chi.chisquare(chitotals);
	pval_allele_one = -1;
	if(chi_allele_one > -1){
	    double pvalue, p, bound, df = 1;
	    pvalue = -1;
		int code = 1, status;
		if(chi_allele_one > -1){
			cdfchi(&code, &p, &pvalue, &chi_allele_one, &df, &status, &bound);
		}
		allele_df = df;
		pval_allele_one = pvalue;
	}
	//		}
	//		else{
				//chi_allele[i] = -1;
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
	chitotals[0].clear();
	chitotals[1].clear();

	chitotals[0].push_back(af->getAonehomoCon());
	chitotals[0].push_back(af->getHetCon());
	chitotals[0].push_back(af->getAtwohomoCon());
	chitotals[1].push_back(af->getAonehomoCa());
	chitotals[1].push_back(af->getHetCa());
	chitotals[1].push_back(af->getAtwohomoCa());
	chi_arm_one = arm.armitage(chitotals);
	pval_arm_one = -1;
	if(chi_arm_one >= 0){
		double pvalue, p, bound, df = 1;
	    pvalue = -1;
	    int code = 1, status;
	    if(chi_arm_one > -1){
	    	cdfchi(&code, &p, &pvalue, &chi_arm_one, &df, &status, &bound);
	    }
	    arm_df = df;
		pval_arm_one = pvalue;
	}
	chitotals[0].clear();
	chitotals[1].clear();

}

/*
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
	//af->process(samples, families, markers, marker_map);

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
	double zt = ltqnorm(1 - (1 - options.getCI()) / 2);

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
		if((*markers)[0]->getDetailHeaders().size() > 0){
			gpvals << "\t" << (*markers)[0]->getDetailHeaders();
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
			if((*markers)[i]->isEnabled() && isValidMarker((*markers)[i], &options, prev_base, prev_chrom)){
				gpvals << (*markers)[i]->toString();
				if((*markers)[i]->isMicroSat()){
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
					vector<Family*> myfams = gfams[mygroup];//generateFamilySet(&mysamps);
					AlleleFrequency* gaf = new AlleleFrequency(&mysamps, &myfams);
					gaf->setRank(rank);
					gaf->setOptions(options);

					gaf->calcOne((*markers)[i]);
					cases = 0;
					controls = 0;
					for(int s = 0; s < (int)mysamps.size(); s++){
						Sample* samp = mysamps[s];
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
			             double pvalue, p, bound, df = 2;
			             pvalue = -1;
			             int code = 1, status;
			             if(gchi_geno > -1){
			                 cdfchi(&code, &p, &pvalue, &gchi_geno, &df, &status, &bound);
			             }

						gpval_geno = pvalue;
					}

					gpval_geno_exact = -1;
					chitotals[0].push_back(gaf->getAoneCa_count());
					chitotals[0].push_back(gaf->getAtwoCa_count());
					chitotals[1].push_back(gaf->getAoneCon_count());
					chitotals[1].push_back(gaf->getAtwoCon_count());
					gchi_allele = chi.chisquare(chitotals);
					gpval_allele = -1;
					if(gchi_allele > -1){
		            	double pvalue, p, bound, df = 1;
		             	pvalue = -1;
		             	int code = 1, status;
		             	if(gchi_allele > -1){
		                	cdfchi(&code, &p, &pvalue, &gchi_allele, &df, &status, &bound);
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
					chitotals[0].clear();
					chitotals[1].clear();

					chitotals[0].push_back(gaf->getAonehomoCon());
					chitotals[0].push_back(gaf->getHetCon());
					chitotals[0].push_back(gaf->getAtwohomoCon());
					chitotals[1].push_back(gaf->getAonehomoCa());
					chitotals[1].push_back(gaf->getHetCa());
					chitotals[1].push_back(gaf->getAtwohomoCa());
					gchi_arm = arm.armitage(chitotals);
					gpval_arm = -1;
					if(gchi_arm >= 0){
			            double pvalue, p, bound, df = 1;
	             		pvalue = -1;
	             		int code = 1, status;
	             		if(gchi_arm > -1){
	                 		cdfchi(&code, &p, &pvalue, &gchi_arm, &df, &status, &bound);
	             		}

						gpval_arm = pvalue;
					}
					chitotals[0].clear();
					chitotals[1].clear();
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
						<< (*markers)[i]->getAllele1() << ":" << (*markers)[i]->getAllele2() << "\t"
						<< godds_ratio << "\t"
						<< (*markers)[i]->getAllele1() << "\t"
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
		if((*samples)[i]->isEnabled()){
			if((*samples)[i]->getPheno() == 2){
				cases++;
			}
			else if((*samples)[i]->getPheno() == 1){
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
		if((*markers)[i]->isEnabled() && isValidMarker((*markers)[i], &options, prev_base, prev_chrom)){
			if((*markers)[i]->isMicroSat()){
				continue;
			}
			af->calcOne((*markers)[i]);
			//orig_num_markers++;
	//		if(af->getAonehomoCa(i) > 4 && af->getAtwohomoCa(i) > 4 && af->getHetCa(i) > 4 && af->getAonehomoCon(i) > 4 && af->getAtwohomoCon(i) > 4 && af->getHetCon(i) > 4){
				//cout << i << ": " << "Case: H1=" << af->getAonehomoCa(i) << " HET=" << af->getHetCa(i) << " H2=" << af->getAtwohomoCa(i) << ":: Control: H1=" << af->getAonehomoCon(i) << " HET=" << af->getHetCon(i) << " H2=" << af->getAtwohomoCon(i) << endl;
				chi_geno[i] = calcChiGeno(af->getAonehomoCa(), af->getHetCa(), af->getAtwohomoCa(), af->getAonehomoCon(), af->getHetCon(), af->getAtwohomoCon());
				pval_geno[i] = -1;
				if(chi_geno[i] >= 0){
		             double pvalue, p, bound, df = 2;
		             pvalue = -1;
		             int code = 1, status;
		             if(chi_geno[i] > -1){
		                 cdfchi(&code, &p, &pvalue, &chi_geno[i], &df, &status, &bound);
		             }

					pval_geno[i] = pvalue;
				}

	//		}
	//		else{
				//chi_geno[i] = -1;
				pval_geno_exact[i] = -1;
	//		}
	//		if(af->getAoneCon_count(i) > 4 && af->getAtwoCon_count(i) > 4 && af->getAoneCa_count(i) > 4 && af->getAtwoCa_count(i) > 4){
				//cout << (*markers)[i]->getProbeID() << "\t" << af->getAoneCa_count(i) << "\t" << af->getAtwoCa_count(i) << "\t" << af->getAoneCon_count(i) << "\t" << af->getAtwoCon_count(i) << endl;
				chitotals[0].push_back(af->getAoneCa_count());
				chitotals[0].push_back(af->getAtwoCa_count());
				chitotals[1].push_back(af->getAoneCon_count());
				chitotals[1].push_back(af->getAtwoCon_count());
				chi_allele[i] = chi.chisquare(chitotals);
				pval_allele[i] = -1;
				if(chi_allele[i] > -1){
		             double pvalue, p, bound, df = 1;
		             pvalue = -1;
		             int code = 1, status;
		             if(chi_allele[i] > -1){
		                 cdfchi(&code, &p, &pvalue, &chi_allele[i], &df, &status, &bound);
		             }
					pval_allele[i] = pvalue;
				}
	//		}
	//		else{
				//chi_allele[i] = -1;
				pval_allele_exact[i] = fish.fisher_2_2(af->getAoneCa_count(), af->getAtwoCa_count(), af->getAoneCon_count(), af->getAtwoCon_count());
				if(af->getAtwoCa_count() == 0 ||  af->getAoneCon_count() == 0){
					odds_ratio[i] = -1;
					ci_l[i] = -999;
					ci_u[i] = -999;
				}
				else{
					odds_ratio[i] = (double)((double)af->getAoneCa_count() * (double)af->getAtwoCon_count()) / (double)((double)af->getAtwoCa_count() * (double)af->getAoneCon_count());
					double lOR = log(odds_ratio[i]);
					double SE = sqrt(1/(double)af->getAoneCa_count() + 1/(double)af->getAtwoCa_count() + 1/(double)af->getAoneCon_count() + 1/(double)af->getAtwoCon_count());
					ci_l[i] = exp( lOR - zt * SE);
					ci_u[i] = exp(lOR + zt * SE);
				}
	//		}
			chitotals[0].clear();
			chitotals[1].clear();

			chitotals[0].push_back(af->getAonehomoCon());
			chitotals[0].push_back(af->getHetCon());
			chitotals[0].push_back(af->getAtwohomoCon());
			chitotals[1].push_back(af->getAonehomoCa());
			chitotals[1].push_back(af->getHetCa());
			chitotals[1].push_back(af->getAtwohomoCa());
			chi_arm[i] = arm.armitage(chitotals);
			pval_arm[i] = -1;
			if(chi_arm[i] >= 0){
	             double pvalue, p, bound, df = 1;
	             pvalue = -1;
	             int code = 1, status;
	             if(chi_arm[i] > -1){
	                 cdfchi(&code, &p, &pvalue, &chi_arm[i], &df, &status, &bound);
	             }

				pval_arm[i] = pvalue;
			}
			chitotals[0].clear();
			chitotals[1].clear();
		}
	}

/*	markers->resetAlleleInfo();

    FAM::iterator fam_iter;
    FAM* myfams = families->getList();
    MKR* mymarkers = markers->getList();

	AlleleFrequency* af = new AlleleFrequency();
	af->setRank(rank);
	af->process(samples, families, markers, marker_map);

	MKR::iterator m_iter;
	for(m_iter = mymarkers->begin(); m_iter != mymarkers->end(); m_iter++){
		m_iter->second.calcHW();
	}
*/
	if(af){
		delete(af);
		af = NULL;
	}
}

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
	//cout << "  EXP: Case: H1=" << cahomo1exp << " HET=" << cahetexp << " H2=" << cahomo2exp << " :: Controls: H1=" << conhomo1exp << " HET=" << conhetexp << " H2=" << conhomo2exp << "\n";

	if(cahomo1exp > 0){
	//	cout << "CaH1=" << ((pow(((float)cahomo1 - cahomo1exp),2)) / cahomo1exp) << ", ";
		chi += (pow(((float)cahomo1 - cahomo1exp),2)) / cahomo1exp;
	}
	if(conhomo1exp > 0){
	//	cout << "ConH1=" << ((pow(((float)conhomo1 - conhomo1exp),2)) / conhomo1exp) << ", ";
		chi += (pow(((float)conhomo1 - conhomo1exp),2)) / conhomo1exp;
	}
	if(cahetexp > 0){
	//	cout << "CaHet=" << ((pow(((float)cahet - cahetexp),2)) / cahetexp) << ", ";
		chi += (pow(((float)cahet - cahetexp),2)) / cahetexp;
	}
	if(conhetexp > 0){
	//	cout << "ConHet=" << ((pow(((float)conhet - conhetexp),2)) / conhetexp) << ", ";
		chi += (pow(((float)conhet - conhetexp),2)) / conhetexp;
	}
	if(cahomo2exp > 0){
	//	cout << "CaH2=" << ((pow(((float)cahomo2 - cahomo2exp),2)) / cahomo2exp) << ", ";
		chi += (pow(((float)cahomo2 - cahomo2exp),2)) / cahomo2exp;
	}
	if(conhomo2exp > 0){
	//	cout << "ConH2=" << ((pow(((float)conhomo2 - conhomo2exp),2)) / conhomo2exp) << endl;
		chi += (pow(((float)conhomo2 - conhomo2exp),2)) / conhomo2exp;
	}

	return chi;
}

/*float CaConChisq::calcHW_exact(int obs_hets, int obs_hom1, int obs_hom2){
*
 *   This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
 *        Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
 *             Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
 *
 *                 Written by Jan Wigginton
 *
    if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
    {
        printf("FATAL ERROR - SNP-HWE: Current genotype configuration (%d  %d %d ) includes a"
        " negative count", obs_hets, obs_hom1, obs_hom2);
        exit(EXIT_FAILURE);
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
        printf("FATAL ERROR - SNP-HWE: Unable to allocate array for heterozygote probabilities" );
        exit(EXIT_FAILURE);
    }

    int i;
    for (i = 0; i <= rare_copies; i++)
        het_probs[i] = 0.0;
       // start at midpoint
        int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
        // check to ensure that midpoint and rare alleles have same parity
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

           // 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
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

           // add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
            curr_homr--;
            curr_homc--;
        }

        for (i = 0; i <= rare_copies; i++)
            het_probs[i] /= sum;

        * alternate p-value calculation for p_hi/p_lo
 *         double p_hi = het_probs[obs_hets];
 *                 for (i = obs_hets + 1; i <= rare_copies; i++)
 *                             p_hi += het_probs[i];
 *
 *                                     double p_lo = het_probs[obs_hets];
 *                                             for (i = obs_hets - 1; i >= 0; i--)
 *                                                         p_lo += het_probs[i];
 *
 *                                                                 double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
 *                                                                         *

        double p_hwe = 0.0;
        //  p-value calculation for p_hwe
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


float CaConChisq::calcHW(float aonefreq, float atwofreq, int aoneobs, int hetobs, int atwoobs, int population){
	bool DEBUG = false;
	if(DEBUG){
        cout << "Initial:\t" << aonefreq << "\t" << atwofreq << "\t" << aoneobs << "\t" << atwoobs << "\t" << hetobs << "\t" << population << endl;
    }
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

    if(expp2 > 0){
        chi1 = (pow(((float)aoneobs - expp2), 2)) / expp2;
    }
    if(expq2 > 0){
        chi2 = (pow(((float)atwoobs - expq2),2)) / expq2;
    }
    if(exppq > 0){
        chi3 = (pow(((float)hetobs - exppq),2)) / exppq;
    }


    double chi = chi1 + chi2 + chi3;

    if(DEBUG){
        cout << "ChiTot:\t" << chi1 << "\t" << chi2 << "\t" << chi3 << "\t=" << chi << endl;
    }
    //float results = _subchisqrprob(1,chi);
	double p, results, bound, df = 1;
	int code = 1, status;
    cdfchi(&code, &p, &results, &chi, &df, &status, &bound);//chiprobP(chi,1);//ChiSquare::pfromchi(chi, 1);
    if(DEBUG){
      cout << "Results: " << results << endl;
    }
    return ((float)results);
}
*/
/*
void CaConChisq::process(Connection* con, Families* f, Markers* m){
	markers = m;
	families = f;

    if(markers == NULL){
        cout << "Markers are empty...Obtaining Markers from DB..." << endl;
        markers = new Markers();
        markers->fillFromDB(con);
        cout << "Markers obtained: " << markers->getSize() << endl;
        *m = *markers;
    }
    if(families == NULL){
        cout << "Families are empty...Obtaining Families from DB..." << endl;
        families = new Families();
        families->fillFromDB(con);
        cout << "Families obtained: " << families->getSize() << endl;
        *f = *families;
    }


	markers->resetAlleleInfo();

    FAM::iterator fam_iter;
    FAM* myfams = families->getList();
    MKR* mymarkers = markers->getList();

	AlleleFrequency* af = new AlleleFrequency();
	af->setRank(rank);
	af->hw_process(con, families, markers);

	MKR::iterator m_iter;
//	ofstream myoutput ("hwinit.txt", ios::out | ios::app);
//	myoutput.precision(4);
	for(m_iter = mymarkers->begin(); m_iter != mymarkers->end(); m_iter++){
		m_iter->second.calcHW();
		AlleleInfo* a = m_iter->second.getAlleleInfo();
		//myoutput << m_iter->second.getRSID() << "\t" << a->getParentHW() << "\t" << a->getChildHW() << endl;
	}
//	myoutput.close();
}
*/
}
