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
#include <MultComparison.h>
#include "ProcessLogReg.h"
#include <Options.h>
#include <General.h>
#include <Helper.h>
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
using namespace Methods;

void ProcessLogReg::FilterSummary(){

	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessLogReg::PrintSummary(){
	int msize = data_set->num_loci();
	for(int m = 0; m < msize; m++){
		data_set->get_locus(m)->setFlag(false);
	}

}

void ProcessLogReg::filter(){
}

void ProcessLogReg::doFilter(Methods::Marker* mark, double value){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		if(mark->isEnabled() && !mark->isFlagged()){
			bool inc = false;
			if(options.doThreshMarkersLow() && dLess(value, options.getThreshMarkersLow())){
				mark->setEnabled(false);
				inc = true;
			}
			if(options.doThreshMarkersHigh() && dGreater(value, options.getThreshMarkersHigh())){
				mark->setEnabled(false);
				inc = true;
			}
			if(inc){
				orig_num_markers++;
			}
		}
	}
}


void ProcessLogReg::process(DataSet* ds){
	data_set = ds;

	//check if new covariate file is listed...or covariate name.
	//create vector of covariate indexes to use if specified.

    string fname = opts::_OUTPREFIX_ + "logreg" + options.getOut() + ".txt";
    if(!overwrite){
        fname += "." + getString<int>(order);
    }
    ofstream lrout (fname.c_str());
    if(!lrout){
        opts::printLog("Error opening " + fname + "!  Exiting!\n");
        throw MethodException("");
    }
    lrout.precision(4);

    lrout << "Chrom\trsID\tProbeID\tBPLOC\tReference_Allele\tTest\tNMISS\tOR\tSE\tL" << getString<double>(options.getCI()*100) << "\t" << "U" << getString<double>(options.getCI()*100) << "\tSTAT\tPvalue\n";

//    DataSet* trimmed_data = new DataSet();
//    trimmed_data->set_markers(ds->get_markers());
//    trimmed_data->set_covariates(ds->get_covariates());
//    trimmed_data->set_traits(ds->get_traits());
//	for(int s = 0; s < data_set->num_inds(); s++){
//		Sample* samp = data_set->get_sample(s);
//		if(samp->isEnabled() && (samp->getPheno() == 1 || samp->getPheno() == 2)){
//			trimmed_data->add_ind(samp);
//		}
//	}
//	trimmed_data->recreate_family_vector();
//	trimmed_data->set_affection_vectors();


	LogisticRegression lr;
	lr.set_parameters(&options);
	ds->set_missing_covalues(-99999);
	lr.resetDataSet(ds);
	if(ds->num_covariates() == 0 && ds->num_traits() == 0)
		lr.setFullInteraction(true);

	lr.setModelType(options.getLRModelType());

	double zt = ltqnorm(1.0 - (1.0 - options.getCI()) / 2.0);
	int prev_base = 0;
	int prev_chrom = -1;
	InputFilter ct_filter;
	vector<string> use_covs = options.getCovars();
	vector<string> use_traits = options.getTraits();
	if(options.doCovarsName()){
		ct_filter.add_covariate_list(&use_covs);
		ct_filter.add_covariate_filter(InputFilter::IncludeCovariateFilter);
	}
//	if(options.doTraitsName()){
//		ct_filter.add_trait_list(&use_traits);
//		ct_filter.add_trait_filter(InputFilter::IncludeTraitFilter);
//	}
	bool cov_use = true;
//	bool trait_use = true;
//	if(options.doCovarsName() && !options.doTraitsName()){
//		trait_use = false;
//	}
//	if(!options.doCovarsName() && options.doTraitsName()){
//		cov_use = false;
//	}

	vector<double> chis(ds->num_loci(), 0);
	vector<double> pvals(ds->num_loci(), 0);

	for(int m = 0; m < (int)ds->num_loci(); m++){
		Marker* mark = ds->get_locus(m);
		if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
			int nmiss = 0;
			for(int s = 0; s < ds->num_inds(); s++){
				Sample* samp = ds->get_sample(s);
				if(samp->isEnabled() && !samp->getAmissing(mark->getLoc()) && (samp->getPheno() == 1 || samp->getPheno() == 2)){
					nmiss++;
				}
			}
			vector<unsigned int> model;
			vector<unsigned int> covs;
			vector<unsigned int> traits;
			model.push_back(m);
			if(cov_use){
				for(int c = 0; c < ds->num_covariates(); c++){
					bool use = true;
					for(int f = 0; f < ct_filter.num_covariate_filters(); f++){
						use = ct_filter.run_covariate_filter(f, ds->get_covariate_name(c));
					}
					if(use){
						covs.push_back(c);
					}
				}
			}
//			if(trait_use){
//				for(int c = 0; c < ds->num_traits(); c++){
//					bool use = true;
//					for(int f = 0; f < ct_filter.num_trait_filters(); f++){
//						use = ct_filter.run_trait_filter(f, ds->get_trait_name(c));
//					}
//					if(use){
//						traits.push_back(c);
//					}
//				}
//			}
			if(covs.size() == 0){// && traits.size() == 0){
				lr.calculate(model);
			}
			else{
				lr.calculate(model, covs, traits);
			}
			vector<double> coefs = lr.getCoefficients();
			vector<double> ses = lr.getCoeffStandardErr();

			for(unsigned int c = 0; c < model.size(); c++){
				lrout << mark->toString() << "\t" << mark->getReferent() << "\t" << options.getLRModelType();
				lrout << "\t" << nmiss;
				lrout << "\t" << exp(coefs[c]);
				double se = ses[c];
				double Z = coefs[c] / se;
				lrout << "\t" << se
					<< "\t" << exp(coefs[c] - zt * se)
					<< "\t" << exp(coefs[c] + zt * se)
					<< "\t" << Z;
				double zz = Z*Z;
				double pvalue, p, bound, df = 1;
				int code = 1, status;
				//cdfchi(&code, &p, &pvalue, &zz, &df, &status, &bound);
				pvalue = p_from_chi(zz, df);
				lrout << "\t" << pvalue;
				chis[m] = zz;
				pvals[m] = pvalue;
//				if(covs.size() == 0 && traits.size() == 0){
//					lrout << "\t" << lr.getFullInteractionP();
//				}
//				else{
//					lrout << "\t" << lr.getOverallP();
//				}

				lrout << endl;
			}
			int buffer = model.size();
			for(unsigned int c = 0; c < covs.size(); c++){
				lrout << mark->toString() << "\t" << mark->getReferent() << "\t" << ds->get_covariate_name(covs[c]);
				lrout << "\t" << exp(coefs[buffer + c]);
				double se = ses[buffer + c];
				double Z = coefs[buffer + c] / se;
				lrout << "\t" << se
					<< "\t" << exp(coefs[buffer + c] - zt * se)
					<< "\t" << exp(coefs[buffer + c] + zt * se)
					<< "\t" << Z;
				double zz = Z*Z;
				double pvalue, p, bound, df = 1;
				int code = 1, status;
				//cdfchi(&code, &p, &pvalue, &zz, &df, &status, &bound);
				pvalue = p_from_chi(zz, df);
				lrout << "\t" << pvalue;
//				if(covs.size() == 0 && traits.size() == 0){
//					lrout << "\t" << lr.getFullInteractionP();
//				}
//				else{
//					lrout << "\t" << lr.getOverallP();
//				}
				lrout << endl;
			}
//			buffer += covs.size();
//			for(unsigned int c = 0; c < traits.size(); c++){
//				lrout << mark->toString() << "\t" << mark->getReferent() << "\t" << ds->get_trait_name(traits[c]);
//				lrout << "\t" << exp(coefs[buffer + c]);
//				double se = ses[buffer + c];
//				double Z = coefs[buffer + c] / se;
//				lrout << "\t" << se
//					<< "\t" << exp(coefs[buffer + c] - zt * se)
//					<< "\t" << exp(coefs[buffer + c] + zt * se)
//					<< "\t" << Z;
//				double zz = Z*Z;
//				double pvalue, p, bound, df = 1;
//				int code = 1, status;
//				cdfchi(&code, &p, &pvalue, &zz, &df, &status, &bound);
//				lrout << "\t" << pvalue;
///				if(covs.size() == 0 && traits.size() == 0){
///					lrout << "\t" << lr.getFullInteractionP();
///				}
///				else{
///					lrout << "\t" << lr.getOverallP();
///				}
//				lrout << endl;
//			}
		}
	}

	if(options.doMultCompare()){
		string fcomp = opts::_OUTPREFIX_ + "logreg_comparisons" + options.getOut() + ".txt";
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
			if(ds->get_locus(0)->getDetailHeaders().size() > 0){
				COMP << "\t" << ds->get_locus(0)->getDetailHeaders();
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
		for(unsigned int m = 0; m < ds->num_loci(); m++){
			Marker* mark = ds->get_locus(m);

			if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
				COMP << mark->toString() << "\tLOGREG\t";
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
	for(unsigned int m = 0; m < ds->num_loci(); m++){
		Marker* mark = ds->get_locus(m);

		if(mark->isEnabled() && isValidMarker(mark, &options, prev_base, prev_chrom)){
			doFilter(mark, pvals[m]);
		}
	}

	lrout.close();
}
