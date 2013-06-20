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


#include "ProcessLogReg.h"
#include <LogisticRegression.h>
#include <iostream>
#include <vector>
#include <map>

#include <Options.h>
#include <MultComparison.h>
#include <Helpers.h>
#include <MethodException.h>
#include <InputFilter.h>

using std::string;
using std::ofstream;
using std::vector;
using std::map;
using Methods::DataSet;
using Methods::Marker;
using Methods::LogisticRegression;
using Methods::MultComparison;
using Methods::opts;
using Methods::Helpers;
using Methods::Sample;
using Methods::MethodException;
using Methods::InputFilter;

string const ProcessLogReg::stepname = ProcessLogReg::doRegister("logreg");

void ProcessLogReg::PrintSummary(){
	int msize = data_set->num_loci();
	for(int m = 0; m < msize; m++){
		data_set->get_locus(m)->setFlag(false);
	}

}

void ProcessLogReg::doFilter(Methods::Marker* mark, double value){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		if(mark->isEnabled() && !mark->isFlagged()){
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


void ProcessLogReg::process(DataSet* ds)
{
	data_set = ds;
	vector<int> good_markers = Helpers::findValidMarkersIndexes(data_set->get_markers(), &options);

	if(options.doGroupFile()){
		opts::printLog("Reading group information [" + options.getGroupFile() + "]\n");
		options.readGroups(ds->get_samples());
	}

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

    lrout << "Chrom\trsID\tProbeID\tBPLOC\t";
    if(options.doGroupFile()){
    	lrout << "GRP\t";
    }
    lrout << "Reference_Allele\tTest\tNMISS\tBETA\tOR\tSE\tL" << getString<double>(options.getCI()*100) << "\t" << "U" << getString<double>(options.getCI()*100) << "\tSTAT\tPvalue\n";

    string fnamesv = opts::_OUTPREFIX_ + "logreg_synthview" + options.getOut() + ".txt";
    if(!overwrite){
    	fnamesv += "." + getString<int>(order);
    }
    ofstream lrsvout;
    if(options.doOutputSynthView()){
    	lrsvout.open(fnamesv.c_str());
    	if(!lrsvout){
    		opts::printLog("Error opening " + fnamesv + "! Exiting!\n");
    		throw MethodException("Error opening " + fnamesv + "! Exiting!\n");
    	}
    	lrsvout.precision(4);
    	lrsvout << "SNP\tChromosome\tLocation";
    }
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
	lr.resetDataSet(ds);
	lr.set_parameters(&options);
	ds->set_missing_covalues(-99999);
	map<string, vector<Sample*> > groups = options.getGroups();
	map<string, vector<Sample*> >::iterator group_iter;
	if(groups.size() == 0){
		groups["GROUP_1"] = *(ds->get_samples());
	}
	if(options.doOutputSynthView()){
		for(group_iter = groups.begin(); group_iter != groups.end(); group_iter++){
			lrsvout << "\t" << group_iter->first << ":pval" << "\t" << group_iter->first << ":beta" << "\t" << group_iter->first << ":N";
		}
		lrsvout << endl;
	}


	vector<vector<double> > gchis(ds->num_loci());
	vector<vector<double> > gpvals(ds->num_loci());
	vector<vector<int> > gnmiss(ds->num_loci());
	vector<vector<double> > gcoefs(ds->num_loci());
	vector<vector<double> > gexpcoefs(ds->num_loci());
	vector<vector<double> > gupper(ds->num_loci());
	vector<vector<double> > glower(ds->num_loci());
	vector<vector<double> > gstat(ds->num_loci());

	vector<double> chis(ds->num_loci(), 0);
	vector<double> pvals(ds->num_loci(), 0);
	int prev_base = 0;
	int prev_chrom = -1;

	if(ds->num_covariates() == 0 && ds->num_traits() == 0)
		lr.setFullInteraction(true);

	lr.setModelType(options.getLRModelType());

	double zt = Helpers::ltqnorm(1.0 - (1.0 - options.getCI()) / 2.0);
	InputFilter ct_filter;

	vector<string> use_covs = options.getCovars();
	vector<string> use_traits = options.getTraits();
	if(options.doCovarsName()){
		ct_filter.add_covariate_list(&use_covs);
		ct_filter.add_covariate_filter(InputFilter::IncludeCovariateFilter);
	}

	bool cov_use = true;



	for(int m = 0; m < (int)good_markers.size(); m++){//(int)ds->num_loci(); m++){
		Marker* mark = ds->get_locus(good_markers[m]);//ds->get_locus(m);
		if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
			if(options.doOutputSynthView())
			{
				lrsvout << mark->getRSID() << "\t" << mark->getChrom() << "\t" << mark->getBPLOC();
			}

	for(group_iter = groups.begin(); group_iter != groups.end(); group_iter++)
	{
		DataSet* tempds = new DataSet;

		tempds->set_markers(ds->get_markers());
		tempds->set_samples(&group_iter->second);
		tempds->recreate_family_vector();
		tempds->set_affection_vectors();
		tempds->set_marker_map(ds->get_marker_map());
		tempds->set_covariates(ds->get_covariates());
		tempds->set_covariate_map(ds->get_covariate_map());
		tempds->set_traits(ds->get_traits());
		tempds->set_trait_map(ds->get_trait_map());

		lr.resetDataSet(tempds);

			int nmiss = 0;
			for(int s = 0; s < tempds->num_inds(); s++)
			{
				Sample* samp = tempds->get_sample(s);
				if(samp->isEnabled() && !samp->getAmissing(mark->getLoc()) && (samp->getPheno() == 1 || samp->getPheno() == 2))
				{
					nmiss++;
				}
			}
			vector<unsigned int> model;
			vector<unsigned int> covs;
			vector<unsigned int> traits;
			model.push_back(good_markers[m]);
			if(cov_use)
			{
				for(int c = 0; c < ds->num_covariates(); c++)
				{
					bool use = true;
					for(int f = 0; f < ct_filter.num_covariate_filters(); f++)
					{
						use = ct_filter.run_covariate_filter(f, tempds->get_covariate_name(c));
					}
					if(use){
						covs.push_back(c);
					}
				}
			}


			if(covs.size() == 0){// && traits.size() == 0){
			//	cout << "on " << m << endl;
				lr.calculate(model);
			}
			else{
				lr.calculate(model, covs, traits);
			}
			vector<double> coefs = lr.getCoefficients();
			//cout << "coefs\n";
			vector<double> ses = lr.getCoeffStandardErr();
			//cout << "stderr\n";
			for(unsigned int c = 0; c < model.size(); c++){

				lrout << mark->toString() << "\t";
				if(options.doGroupFile()){
					lrout << group_iter->first << "\t";
				}
				lrout << mark->getReferent() << "\t" << options.getLRModelType();
				lrout << "\t" << nmiss;
				lrout << "\t" << coefs[c];
				lrout << "\t" << exp(coefs[c]);
				double se = ses[c];

				double Z = coefs[c] / se;

				lrout << "\t" << se
					<< "\t" << exp(coefs[c] - zt * se)
					<< "\t" << exp(coefs[c] + zt * se)
					<< "\t" << Z;
				double zz = Z*Z;
				double pvalue = 1 , df = 1;
//				int code = 1, status;
				//cdfchi(&code, &p, &pvalue, &zz, &df, &status, &bound);
				//cout << "pre-p_from_chi: " << se << " : " << Z << " : " << zz << " : " << df << endl;
				if(se > 0)
				{
					if(!isinf(zz)){
						pvalue = Helpers::p_from_chi(zz, df);
						lrout << "\t" << pvalue;
					}
					else{
						lrout << "\tinf";
					}
				}
				else{
					lrout << "\tnan";
				}
				chis[m] = zz;
				pvals[m] = pvalue;

				gchis[m].push_back(zz);
				gpvals[m].push_back(pvalue);

//				if(covs.size() == 0 && traits.size() == 0){
//					lrout << "\t" << lr.getFullInteractionP();
//				}
//				else{
//					lrout << "\t" << lr.getOverallP();
//				}
				lrout << endl;

				if(options.doOutputSynthView()){
					lrsvout << "\t" << pvalue << "\t" << coefs[c] << "\t" << nmiss;
				}

				//cout << "calc done " << c << endl;
			}
			int buffer = model.size();
			for(int c = 0; c < (int)covs.size(); c++)
			{
				lrout << mark->toString() << "\t";
				if(options.doGroupFile()){
					lrout << group_iter->first << "\t";
				}
				lrout << mark->getReferent() << "\t" << ds->get_covariate_name(covs[c]);
				lrout << "\t" << nmiss;
				lrout << "\t" << coefs[buffer + c];
				lrout << "\t" << exp(coefs[buffer + c]);

				double se = ses[buffer + c];
				double Z = coefs[buffer + c] / se;

				lrout << "\t" << se
					<< "\t" << exp(coefs[buffer + c] - zt * se)
					<< "\t" << exp(coefs[buffer + c] + zt * se)
					<< "\t" << Z;

				double zz = Z*Z;
				double pvalue, df = 1;
//				int code = 1, status;
				//cdfchi(&code, &p, &pvalue, &zz, &df, &status, &bound);
				if(se > 0){
					if(!isinf(zz)){
						pvalue = Helpers::p_from_chi(zz, df);
						lrout << "\t" << pvalue;
					}
					else{
						lrout << "\tinf";
					}
				}
				else{
					lrout << "\tnan";
				}
//				if(covs.size() == 0 && traits.size() == 0){
//					lrout << "\t" << lr.getFullInteractionP();
//				}
//				else{
//					lrout << "\t" << lr.getOverallP();
//				}
				lrout << endl;
			}

			delete tempds;
		}//end group_iter
	if(options.doOutputSynthView()){
		lrsvout << endl;
	}
	}
	}

	if(options.doMultCompare())
	{
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
		for(unsigned int m = 0; m < good_markers.size(); m++){//ds->num_loci(); m++){
			Marker* mark = ds->get_locus(good_markers[m]);//ds->get_locus(m);

			if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
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
	for(unsigned int m = 0; m < good_markers.size(); m++){//ds->num_loci(); m++){
		Marker* mark = ds->get_locus(good_markers[m]);//ds->get_locus(m);

		if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
			doFilter(mark, pvals[m]);
		}
	}

	lrout.close();
}


