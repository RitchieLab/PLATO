/**********************************************************************************
*                       Marker Genotype Efficiency Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all genotypes and generates a genotype efficiency for all markers.
*
*
* Files generated:
*	percent_breakdown_by_marker.txt
*	percent_breakdown_by_chrom.txt
*       post_marker_geno_eff_filter_summary.txt
*
*File: LD.cc
**********************************************************************************/

#include <iostream>

#include "ProcessCMH.h"
#include <MultComparison.h>
#include <Options.h>
#include <Helpers.h>
#include <CMH.h>
#include <vector>

using Methods::CMH;
using Methods::DataSet;
using Methods::opts;
using Methods::Marker;
using Methods::MultComparison;

using std::string;
using std::getString;
using std::ofstream;
using std::vector;

const string ProcessCMH::stepname = ProcessCMH::doRegister("cmh");

void ProcessCMH::PrintSummary(){
	int msize = data_set->num_loci();
	for(int m = 0; m < msize; m++){
		data_set->get_locus(m)->setFlag(false);
	}

	int ssize = data_set->num_inds();
	for(int m = 0; m < ssize; m++){
		data_set->get_sample(m)->setFlag(false);
	}
}

void ProcessCMH::doFilter(Methods::Marker* mark, double value){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		if(mark->isEnabled()){// && !mark->isFlagged()){
			bool inc = false;
			if(options.doThreshMarkersLow() && Helpers::dLess(value, options.getThreshMarkersLow())){
				mark->setEnabled(false);
				inc = true;
			}
			if(options.doThreshMarkersHigh() && Helpers::dGreater(value, options.getThreshMarkersHigh())){
				mark->setEnabled(false);
				inc = true;
			}
			if(inc){
				orig_num_markers++;
			}
		}
	}
}


void ProcessCMH::process(DataSet* ds){
	data_set = ds;
	int msize = data_set->num_loci();

	CMH cmh;

	if(options.doGroupFile()){
		options.readGroups(ds->get_samples());
	}

	if(options.doGroupFile()){
	}

	options.readClusters(data_set->get_samples());
	cmh.setOverwrite(overwrite);


	cmh.setOptions(options);
	cmh.resetDataSet(data_set);

	string f = opts::_OUTPREFIX_ + "cmh" + options.getOut() + ".txt";
	if (!overwrite) {
		f += "." + getString<int>(order);
	}
	ofstream MHOUT;
	MHOUT.open(f.c_str(), ios::out);
	if (!MHOUT) {
		throw MethodException("Could not open " + f + " for output.\n");
	}


	if(options.doMultCompare()){
	}

	if(options.getCMH2x2xK()){
		MHOUT << "CHR\tSNP\tA1\tA2\tBP\tCHISQ\tP\tOR\tL";
		MHOUT << getString<double>(options.getCI() * 100.0);
		MHOUT << "\t" << "U";
		MHOUT << getString<double>(options.getCI() * 100.0);
		if (options.doBreslowDay())
			MHOUT << "\tCHISQ_BD" << "\t" << "P_BD";
		MHOUT << "\n";
	}
	else{
		MHOUT << "CHR" << "\t" << "SNP" << "\t" << "CHISQ" << "\t" << "P"
				<< "\n";
	}


	MHOUT.precision(4);
	vector<double> chis;
	vector<double> pvals;
	chis.resize(msize, 0);
	pvals.resize(msize, 0);
	int prev_base = 0;
	int prev_chrom = -1;
	vector<Marker*> good_markers = Helpers::findValidMarkers(ds->get_markers(), &options);
	msize = good_markers.size();

	for(int m = 0; m < (int)msize; m++)
	{
		Marker* mark = good_markers[m];//data_set->get_locus(m);

		if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
			cmh.calculate(mark);
			MHOUT << mark->getChrom() << "\t" << mark->getRSID();
			if(options.getCMH2x2xK()){
				MHOUT << "\t" << mark->getAllele1()
						<< "\t" << mark->getAllele2() << "\t"
						<< mark->getBPLOC() << "\t";

				chis[m] = cmh.get_chisq();
				pvals[m] = cmh.get_pval();

				if (Helpers::realnum(cmh.get_chisq()))
				{
					MHOUT << cmh.get_chisq() << "\t" << cmh.get_pval()
							<< "\t";
				}
				else{
					MHOUT << "NA" << "\t" << "NA" << "\t";

				}
				if (Helpers::realnum(cmh.getOR()))
				{
					MHOUT << cmh.getOR() << "\t";

				}
				else{
					MHOUT << "NA" << "\t";

				}
				if (Helpers::realnum(cmh.getOR_lower()))
				{
					MHOUT << cmh.getOR_lower() << "\t";
				}
				else
				{
					MHOUT << "NA" << "\t ";
				}
				if (Helpers::realnum(cmh.getOR_upper()))
				{
					MHOUT << cmh.getOR_upper();
				}
				else
				{
					MHOUT << "NA";
				}

				if(options.doBreslowDay()){
					if (cmh.get_pvalbd() > -1)
					{
						MHOUT << "\t" << cmh.get_chisqbd() << "\t" << cmh.get_pvalbd();
					}
					else
					{
						MHOUT << "\t" << "NA" << "\t" << "NA";
					}
				}
			}
			else
			{
				MHOUT << "\t" << cmh.get_chisq() << "\t" << cmh.get_pval();
				chis[m] = cmh.get_chisq();
				pvals[m] = cmh.get_pval();
			}
			MHOUT << endl;

//			doFilter(mark, cmh.get_pval());
		}
	}

	if(options.doMultCompare()){
		string fcomp = opts::_OUTPREFIX_ + "cmh_comparisons" + options.getOut() + ".txt";
		if (!overwrite) {
			fcomp += "." + getString<int>(order);
		}
		ofstream COMP;

		COMP.open(fcomp.c_str(), ios::out);
		if(!COMP){
			throw MethodException("Could not open " + fcomp + " for output.\n");
		}

		   COMP << "Chrom"
				  << "\trsID"
				  << "\tProbeID"
				  << "\tbploc";
			if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
				COMP << "\t" << data_set->get_locus(0)->getDetailHeaders();
			}

			COMP  << "\tCALC"
				  << "\tOriginal_Pval"
				  << "\tGC"
				  << "\tBONF"
				  << "\tHOLM"
				  << "\tSIDAK_SS"
				  << "\tSIDAK_SD"
				  << "\tFDR_BH"
				  << "\tFDR_BY"
				  << endl;
			opts::addFile("Marker", stepname, fcomp);

			opts::addHeader(fcomp, "CALC");
			opts::addHeader(fcomp, "Original_Pval");
			opts::addHeader(fcomp, "GC");
			opts::addHeader(fcomp, "BONF");
			opts::addHeader(fcomp, "HOLM");
			opts::addHeader(fcomp, "SIDAK_SS");
			opts::addHeader(fcomp, "SIDAK_SD");
			opts::addHeader(fcomp, "FDR_BH");
			opts::addHeader(fcomp, "FDR_BY");


		MultComparison mc(options);
		vector<int> tcnt;
		mc.calculate(chis, tcnt);

		prev_base = 0;
		prev_chrom = -1;
		for(int m = 0; m < (int)msize; m++){//data_set->num_loci(); m++){
			Marker* mark = good_markers[m];//data_set->get_locus(m);

			if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
				COMP << mark->toString() << "\tCMH\t";
				COMP << pvals[m] << "\t"
				<< mc.get_genomic_control(m) << "\t"
				<< mc.get_bonferroni(m) << "\t"
				<< mc.get_holm(m) << "\t"
				<< mc.get_sidak_single_step(m) << "\t"
				<< mc.get_sidak_step_down(m) << "\t"
				<< mc.get_fdr_bh(m) << "\t"
				<< mc.get_fdr_by(m)
				<< endl;

			}
		}
		if(COMP.is_open()){
			COMP.close();
		}

	}

	prev_base = 0;
	prev_chrom = -1;
	for(int m = 0; m < (int)msize; m++){//data_set->num_loci(); m++){
		Marker* mark = good_markers[m];//data_set->get_locus(m);

		if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
			doFilter(mark, pvals[m]);
		}
	}

	if(MHOUT.is_open()){
		MHOUT.close();
	}

}//end method process()


