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


#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include "ProcessCMH.h"
#include <Options.h>
#include <General.h>
#include <Helper.h>
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"


void ProcessCMH::FilterSummary(){

	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

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

void ProcessCMH::filter(){
}


void ProcessCMH::process(DataSet* ds){
	data_set = ds;

	options.readClusters(data_set->get_samples());

	CMH cmh;
	cmh.setOptions(options);
	cmh.setRank(rank);
	cmh.setOverwrite(overwrite);

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


	int prev_base = 0;
	int prev_chrom = -1;

	for(unsigned int m = 0; m < data_set->num_loci(); m++){
		Marker* mark = data_set->get_locus(m);

		if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
			cmh.calculate(mark);
			MHOUT << mark->getChrom() << "\t" << mark->getRSID();
			if(options.getCMH2x2xK()){
				MHOUT << "\t" << mark->getAllele1()
						<< "\t" << mark->getAllele2() << "\t"
						<< mark->getBPLOC() << "\t";

				if (realnum(cmh.get_chisq()))
					MHOUT << cmh.get_chisq() << "\t" << cmh.get_pval()
							<< "\t";
				else
					MHOUT << "NA" << "\t" << "NA" << "\t";

				if (realnum(cmh.getOR()))
					MHOUT << cmh.getOR() << "\t";
				else
					MHOUT << "NA" << "\t";

				if (realnum(cmh.getOR_lower()))
					MHOUT << cmh.getOR_lower() << "\t";
				else
					MHOUT << "NA" << "\t ";

				if (realnum(cmh.getOR_upper()))
					MHOUT << cmh.getOR_upper();
				else
					MHOUT << "NA";

				if(options.doBreslowDay()){
					if (cmh.get_pvalbd() > -1)
						MHOUT << "\t" << cmh.get_chisqbd() << "\t" << cmh.get_pvalbd();
					else
						MHOUT << "\t" << "NA" << "\t" << "NA";
				}
			}
			else{
				MHOUT << "\t" << cmh.get_chisq() << "\t" << cmh.get_pval();
			}
			MHOUT << endl;
		}
	}
}
