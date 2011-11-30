#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <unistd.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#include <fenv.h>
#include <algorithm>
#include <MultComparison.h>
#include "ProcessRunTDT.h"
#include "Chrom.h"
#include <General.h>
//#include "ChiSquare.h"
#include <Helper.h>
#include <cdflib.h>

using namespace std;
using namespace Methods;
string ProcessRunTDT::stepname = "tdt";

void ProcessRunTDT::PrintSummary(){
	string fname1 = opts::_OUTPREFIX_ + "tdt" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	ofstream output (fname1.c_str(), ios::out);
	if(!output){
		opts::printLog("Error opening " + fname1 + ". Exiting!\n");
		throw MethodException("");
	}
	opts::addFile("Marker", stepname, fname1);
	output.precision(4);
	output << "Chrom\trsID\tProbeID\tbploc";
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		output << "\t" << data_set->get_locus(0)->getDetailHeaders();
	}
	output << "\tChi_square\tTDT_pvalue\tTDT_neglog(pvalue)\tNum_Fams\tT:U\tA1:A2\tOR\tL" + getString<double>(options.getCI()*100) + "\t" + "U" + getString<double>(options.getCI()*100) << endl;
	opts::addHeader(fname1, "Chi_square");
	opts::addHeader(fname1, "TDT_pvalue");
	opts::addHeader(fname1, "TDT_neglog(pvalue)");
	opts::addHeader(fname1, "Num_Fams");
	opts::addHeader(fname1, "T:U");
	opts::addHeader(fname1, "A1:A2");
	opts::addHeader(fname1, "OR");
	opts::addHeader(fname1, "L" + getString<double>(options.getCI()*100));
	opts::addHeader(fname1, "U" + getString<double>(options.getCI()*100));

	ofstream compout;
	if(options.doMultCompare()){
		string fname2 = opts::_OUTPREFIX_ + "tdt_comparisons" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if(!overwrite){
			fname2 += "." + getString<int>(order);
		}
		compout.open(fname2.c_str());
		if(!compout){
			opts::printLog("Error opening " + fname2 + ". Exiting!\n");
			throw MethodException("");
		}
		compout.precision(4);
		compout << "Chrom\trsID\tProbeID\tbploc";
		compout  << "\tCALC"
			  << "\tOriginal_Pval"
			  << "\tGC"
			  << "\tBONF"
			  << "\tHOLM"
			  << "\tSIDAK_SS"
			  << "\tSIDAK_SD"
			  << "\tFDR_BH"
			  << "\tFDR_BY"
			  << endl;
		opts::addHeader(fname2, "CALC");
		opts::addHeader(fname2, "Original_Pval");
		opts::addHeader(fname2, "GC");
		opts::addHeader(fname2, "BONF");
		opts::addHeader(fname2, "HOLM");
		opts::addHeader(fname2, "SIDAK_SS");
		opts::addHeader(fname2, "SIDAK_SD");
		opts::addHeader(fname2, "FDR_BH");
		opts::addHeader(fname2, "FDR_BY");

	}

	int msize = good_markers.size();//data_set->num_loci();

	double zt = ltqnorm(1 - (1 - options.getCI()) / 2);

	int prev_base = 0;
	int prev_chrom = -1;
	for(int i = 0; i < msize; i++){
		if(good_markers[i]->isEnabled()){//data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged() && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
			output << good_markers[i]->toString() << "\t";
				//<< ((float)(*markers)[i]->getBPLOC()/1000000.0f) << "\t"
			output << chi[i] << "\t"
				<< pval[i] << "\t"
				<< (double)abs(log10(pval[i])) << "\t";
			//output.precision(4);
			output << fams_used[i] << "\t"
				<< trans[i] << ":" << untrans[i] << "\t";
			if(!good_markers[i]->isMicroSat()){
				output << good_markers[i]->getAllele1() << ":" << good_markers[i]->getAllele2();
			}
			else{
				output << "NA";
			}
			double OR = (double)trans[i] / (double)untrans[i];
			output << "\t" << OR;
			double OR_lower = exp(log(OR) - zt * sqrt(1/trans[i] + 1/untrans[i]));
			double OR_upper = exp(log(OR) + zt * sqrt(1/trans[i] + 1/untrans[i]));
			output << "\t" << OR_lower << "\t" << OR_upper;

				output << endl;
/*
			if(options.doMultCompare()){
				MultComparison mc(options);
				vector<double> chivals;
				chivals.push_back(chi[i]);
				vector<int> tcnt;
				mc.calculate(chivals, tcnt);
				compout << data_set->get_locus(i)->toString() << "\t"
					<< "TDT\t"
					<< pval[i] << "\t"
					<< mc.get_genomic_control(0) << "\t"
					<< mc.get_bonferroni(0) << "\t"
					<< mc.get_holm(0) << "\t"
					<< mc.get_sidak_single_step(0) << "\t"
					<< mc.get_sidak_step_down(0) << "\t"
					<< mc.get_fdr_bh(0) << "\t"
					<< mc.get_fdr_by(0)
					<< endl;
			}
*/
				//output.precision(100);
		}
		good_markers[i]->setFlag(false);
	}
	if(options.doMultCompare()){
		MultComparison mc(options);
		vector<int> tcnt;
		mc.calculate(chi, tcnt);
		for(int i = 0; i < good_markers.size(); i++){
		compout << good_markers[i]->toString() << "\t"
			<< "TDT\t"
			<< pval[i] << "\t"
			<< mc.get_genomic_control(i) << "\t"
			<< mc.get_bonferroni(i) << "\t"
			<< mc.get_holm(i) << "\t"
			<< mc.get_sidak_single_step(i) << "\t"
			<< mc.get_sidak_step_down(i) << "\t"
			<< mc.get_fdr_bh(i) << "\t"
			<< mc.get_fdr_by(i)
			<< endl;
		}
	}

	if(output.is_open()){
		output.close();
	}
	if(compout.is_open()){
		compout.close();
	}
}

void ProcessRunTDT::filter(){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int size = good_markers.size();//data_set->num_loci();
		int prev_base = 0;
		int prev_chrom = -1;
		for(int i = 0; i < size; i++){
			if(good_markers[i]->isEnabled()){//data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged() && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
				bool inc = false;
				if(options.doThreshMarkersHigh() && dGreater(pval[i], options.getThreshMarkersHigh())){
					good_markers[i]->setEnabled(false);//data_set->get_locus(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersLow() && dLess(pval[i],options.getThreshMarkersLow())){
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

void ProcessRunTDT::setThreshold(string thresh){
	options.setUp(thresh);
}

void ProcessRunTDT::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
	    getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
    opts::_MARKERS_WORKING_ -= orig_num_markers;
}

void ProcessRunTDT::process(DataSet* ds){
	data_set = ds;

	good_markers = findValidMarkers(data_set->get_markers(), &options);

	RunTDT tdt(data_set);
	tdt.setOptions(options);

	int msize = good_markers.size();//data_set->num_loci();
	chi.resize(msize);
	pval.resize(msize);
	fams_used.resize(msize);
	maf.resize(msize);
	trans.resize(msize);
	untrans.resize(msize);
////	int ssize = data_set->num_inds();

	int prev_base = 0;
	int prev_chrom = -1;
	for(int m = 0; m < msize; m++){
		if(good_markers[m]->isEnabled()){//data_set->get_locus(m)->isEnabled() && isValidMarker(data_set->get_locus(m), &options, prev_base, prev_chrom)){
			tdt.calculate(good_markers[m]);

			pval[m] = tdt.getPval();
			chi[m] = tdt.getChi();
			fams_used[m] = tdt.getFamsUsed();
			maf[m] = -1;
			trans[m] = tdt.getTransmitted();
			untrans[m] = tdt.getUntransmitted();
		}
	}
}

