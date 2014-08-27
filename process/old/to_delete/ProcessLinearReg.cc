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


#include "ProcessLinearReg.h"
#include <LinearRegression.h>

#include <iostream>
#include <vector>
#include <map>

#include <Options.h>
#include <MultComparison.h>
#include <Helpers.h>
#include <MethodException.h>
//#include <General.h>

//using namespace Methods;
using std::string;
using std::ofstream;
using std::vector;
using std::map;
using Methods::DataSet;
using Methods::Marker;
using Methods::LinearRegression;
using Methods::MultComparison;
using Methods::opts;
using Methods::Helpers;
using Methods::Sample;
using Methods::MethodException;

const string ProcessLinearReg::stepname = ProcessLinearReg::doRegister("linear-reg");

void ProcessLinearReg::PrintSummary(){
	int msize = data_set->num_loci();
	for(int m = 0; m < msize; m++){
		data_set->get_locus(m)->setFlag(false);
	}

}

void ProcessLinearReg::doFilter(Methods::Marker* mark, double value){
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


void ProcessLinearReg::process(DataSet* ds){
	data_set = ds;


		//MEMORY_LEAK
		if(options.doGroupFile()){
			options.readGroups(ds->get_samples());
		}
		//END_MEMORY_LEAK

	vector<Marker*> good_markers = Helpers::findValidMarkers(data_set->get_markers(), &options);
	//MEMORY_LEAK
	map<string, vector<Sample*> > groups = options.getGroups();
	map<string, vector<Sample*> >::iterator group_iter;
	//END_MEMORY_LEAK

	//check if new covariate file is listed...or covariate name.
	//create vector of covariate indexes to use if specified.
	if(options.doCovars())
	{
		opts::printLog("Number of covariates used in Linear Model: " + getString<int>(options.getCovars().size()) + "\n");
	}

    string fname = opts::_OUTPREFIX_ + "linearreg" + options.getOut() + ".txt";
    if(!overwrite){
        fname += "." + getString<int>(order);
    }
    ofstream lrout (fname.c_str());
    if(!lrout){
        opts::printLog("Error opening " + fname + "!  Exiting!\n");
        throw MethodException("");
    }
    lrout.precision(4);

    //MEMORY_LEAK
    lrout << "Chrom\trsID\tProbeID\tBPLOC\t";
    if(options.doGroupFile()){
    	lrout << "GRP\t";
    }
    lrout << "Reference_Allele\tTest\tNMISS\tBeta\tOR\tSE\tSTAT\tPvalue\n";

    string fnamesv = opts::_OUTPREFIX_ + "linearreg_synthview" + options.getOut() + ".txt";
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
    //END_MEMORY_LEAK

	LinearRegression lr;
	lr.setOptions(options);

	//MEMORY_LEAK
	if(groups.size() == 0){
		groups["GROUP_1"] = *(ds->get_samples());
	}
	if(options.doOutputSynthView()){
		for(group_iter = groups.begin(); group_iter != groups.end(); group_iter++){
			lrsvout << "\t" << group_iter->first << ":pval" << "\t" << group_iter->first << ":beta" << "\t" << group_iter->first << ":N";
		}
		lrsvout << endl;
	}
	//lr.resetDataSet(ds);

	//END_MEMORY_LEAK

	vector<double> chis(ds->num_loci(), 0);
	vector<double> main_pvals(ds->num_loci(), 0);

	int prev_base = 0;
	int prev_chrom = -1;
	int msize = good_markers.size();

	for(int i = 0; i < msize; i++){
		Marker* mark = good_markers[i];//ds->get_locus(i);
		if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
			//MEMORY_LEAK
			if(options.doOutputSynthView()){
				lrsvout << mark->getRSID() << "\t" << mark->getChrom() << "\t" << mark->getBPLOC();
			}
			for(group_iter = groups.begin(); group_iter != groups.end(); group_iter++){
				DataSet* tempds = new DataSet();
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

			//END_MEMORY_LEAK	

				lr.calculate(mark);//i);
				vector<double>pvals = lr.getPvalues();
				vector<double>coefs = lr.getCoefs();
				vector<string>labels = lr.getLabels();
				vector<double> vars = lr.getVar();
				vector<double> zs = lr.getZs();


				for(int l = 1; l < (int)labels.size(); l++)
				{
					bool okay = vars[l] < 1e-20 || !Helpers::realnum(vars[l]) ? false : true;
					double se = 0;
					if(okay){
						se = sqrt(vars[l]);
					}
					lrout << mark->getChrom() << "\t" << mark->getRSID() << "\t" << mark->getProbeID() << "\t";
					lrout << mark->getBPLOC();
					if(options.doGroupFile()){
						lrout << "\t" << group_iter->first;
					}
					lrout << "\t" << mark->getReferent() << "\t";
					lrout << labels[l] << "\t" << lr.getCalcMissing() << "\t" << coefs[l] << "\t" << exp(coefs[l]) << "\t" << se << "\t" << zs[l] << "\t" << pvals[l] << endl;

					if(options.doOutputSynthView()){
						//do not display this information if the current label belongs to a covariate...
						//this is an inefficient way to do this, but quick to implement
						vector<string>* covariatesList = ds->get_covariates();
						std::vector<string>::iterator covIterator = std::find_if((*covariatesList).begin(), (*covariatesList).end(), Methods::FindString(labels[l]));
						if(covIterator == (*covariatesList).end())
						{
							lrsvout << "\t" << pvals[l] << "\t" << coefs[l] << "\t" << lr.getCalcMissing();
						}
					}
				}//end loop through labels
				chis[i] = lr.getStatistic();
				main_pvals[i] = pvals[1];
			}//end loop through Groups
			if(options.doOutputSynthView()){
				lrsvout << endl;
			}
		}
//		lrout << options.getLinRModelType() << "\t";
//		lrout << lr.getCalcMissing() << "\t" << lr.getTestCoef() << "\t" << lr.getZ() << "\t" << lr.getPValue() << endl;
	}//end loop through Markers

	if(options.doMultCompare()){
		string fcomp = opts::_OUTPREFIX_ + "linearreg_comparisons" + options.getOut() + ".txt";
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
		for(unsigned int m = 0; m < good_markers.size(); m++){//data_set->num_loci(); m++){
			Marker* mark = good_markers[m];//data_set->get_locus(m);

			if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
				COMP << mark->toString() << "\tLINREG\t";
				COMP << main_pvals[m] << "\t"
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
	for(unsigned int m = 0; m < good_markers.size(); m++){//data_set->num_loci(); m++){
		Marker* mark = good_markers[m];//data_set->get_locus(m);

		if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
			doFilter(mark, main_pvals[m]);
		}
	}

	if(lrout.is_open()){
		lrout.close();
	}
}


