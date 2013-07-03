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
#include "config.h"
#ifdef HAVE_MALLOC_H
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
#include <Helpers.h>
#ifdef PLATOLIB
#include "Controller.h"
#endif
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif
#ifdef PLATOLIB
ProcessLogReg::ProcessLogReg(string bn, int pos, Database* pdb)
{
	name = "Logistic Regression";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
}
#endif

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

void ProcessLogReg::filter()
{
#ifdef PLATOLIB
	for(int m = 0; m < (int)data_set->num_loci(); m++)
	{
		data_set->get_locus(m)->setFlag(false);
	}
#endif
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



void ProcessLogReg::outputResult(ostream& lrout, ostream& lrsvout, LogisticRegression& lr,
	vector<unsigned int>& model, vector<unsigned int>& covs, string groupName, int modelnum,
	DataSet* ds, Marker* mark){
			vector<double> coefs = lr.getCoefficients();
			vector<double> ses = lr.getCoeffStandardErr();

			for(unsigned int c = 0; c < model.size(); c++){

				lrout << mark->toString() << "\t";
				if(options.doGroupFile()){
					lrout << groupName << "\t";
				}
				lrout << mark->getReferent() << "\t" << options.getLRModelType();
				lrout << "\t" << _nmiss;
				lrout << "\t" << coefs[c];
				lrout << "\t" << exp(coefs[c]);
				
				double se = ses[c];
				double Z = coefs[c] / se;

				lrout << "\t" << se
					<< "\t" << exp(coefs[c] - _zt * se)
					<< "\t" << exp(coefs[c] + _zt * se)
					<< "\t" << Z;
				double zz = Z*Z;
				double pvalue = 1 , df = 1;

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
				_chis[modelnum] = zz;
				_pvals[modelnum] = pvalue;

				lrout << endl;

				if(options.doOutputSynthView()){
					lrsvout << "\t" << pvalue << "\t" << coefs[c] << "\t" << _nmiss;
				}

			}	
			
			int buffer = model.size();
			for(int c = 0; c < (int)covs.size(); c++)
			{
				if(mark)
					lrout << mark->toString() << "\t";
				else
					lrout << ds->get_covariate_name(covs[0]) << "\t";
					
				if(options.doGroupFile()){
					lrout << groupName << "\t";
				}
				if(mark)
					lrout << mark->getReferent() << "\t";
					
				lrout << ds->get_covariate_name(covs[c]);

				lrout << "\t";
				lrout << "\t" << coefs[buffer + c];
				lrout << "\t" << exp(coefs[buffer + c]);

				double se = ses[buffer + c];
				double Z = coefs[buffer + c] / se;

				lrout << "\t" << se
					<< "\t" << exp(coefs[buffer + c] - _zt * se)
					<< "\t" << exp(coefs[buffer + c] + _zt * se)
					<< "\t" << Z;

				double zz = Z*Z;
				double pvalue, df = 1;
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

				lrout << endl;
			}
			
			if(mark)
				lrout << mark->toString() << "\t";
			else
				lrout << ds->get_covariate_name(covs[0]) << "\t";
				
			if(options.doGroupFile()){
					lrout << groupName << "\t";
			}
			lrout << "overall";
			lrout << "\t" << _nmiss;
			lrout << "\t-----";
			lrout << "\t-----";
			lrout << "\t-----\t-----\t-----";
			lrout << "\t" << lr.getOverallScore();
			lrout << "\t" << lr.getOverallP() << endl;
			
			
			
}

DataSet* ProcessLogReg::getTempDataSet(DataSet* ds, map<string, vector<Sample*> >::iterator group_iter){
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
		return tempds;	
}



void ProcessLogReg::addCovsTraits(vector<unsigned int>& covs, vector<unsigned int>& traits,
	DataSet* ds, bool cov_use, InputFilter& ct_filter, DataSet* tempds){
	
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

	#ifdef PLATOLIB
	if(trait_use){
		for(int c = 0; c < tempds->num_traits(); c++){
			bool use = true;
			for(int f = 0; f < ct_filter.num_trait_filters(); f++){
				use = ct_filter.run_trait_filter(f, tempds->get_trait_name(c));
			}
			if(use){
				traits.push_back(c);
			}
		}
	}
	#endif
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

		if(options.runCovarsOnly()){
 			lrout << "Variable\t";
		}
		else{
	    lrout << "Chrom\trsID\tProbeID\tBPLOC\t";
	  }
	  if(options.doGroupFile()){
  	  lrout << "GRP\t";
    }
    
    if(!options.runCovarsOnly()){
    	lrout << "Reference_Allele\t";
  	}
    lrout << "Test\tNMISS\tBETA\tOR\tSE\tL" << getString<double>(options.getCI()*100) << "\t" << "U" << getString<double>(options.getCI()*100) << "\tSTAT\tPvalue\n";

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

	int vectorSize;
	if(options.runCovarsOnly()){
		vectorSize=ds->num_covariates();
	}
	else{
		vectorSize=ds->num_loci();
	}
	
	_chis.assign(vectorSize, 0);
	_pvals.assign(vectorSize, 0);
	
	int prev_base = 0;
	int prev_chrom = -1;

	if(ds->num_covariates() == 0 && ds->num_traits() == 0)
		lr.setFullInteraction(true);

	lr.setModelType(options.getLRModelType());

	_zt = Helpers::ltqnorm(1.0 - (1.0 - options.getCI()) / 2.0);
	InputFilter ct_filter;

	vector<string> use_covs = options.getCovars();
	vector<string> use_traits = options.getTraits();
	if(options.doCovarsName()){
		ct_filter.add_covariate_list(&use_covs);
		ct_filter.add_covariate_filter(InputFilter::IncludeCovariateFilter);
	}
	bool cov_use = false;

	if(options.runCovarsOnly()){
		int numCovars = ds->num_covariates();
		vector<unsigned int> addCovars;
		setCovariates(addCovars);
		set<unsigned int> covarMap;
		vector<unsigned int>::iterator iter;
		for(iter = addCovars.begin(); iter != addCovars.end(); iter++){
			covarMap.insert(*iter);
		}
		
		vector<string>* get_covariates();
		double missingCoValue = ds->get_missing_covalue();
		for(int c=0; c < numCovars; c++){
			for(group_iter = groups.begin(); group_iter != groups.end(); group_iter++)
			{
				DataSet* tempds = getTempDataSet(ds, group_iter);
				lr.resetDataSet(tempds);				
				_nmiss = 0;
				for(int s = 0; s < tempds->num_inds(); s++)
				{
					Sample* samp = tempds->get_sample(s);
					if(samp->isEnabled() && !samp->getCovariate(c) != missingCoValue && (samp->getPheno() == 1 || samp->getPheno() == 2))
					{
						_nmiss++;
					}
				}
				
				vector<unsigned int> model;
				vector<unsigned int> covs;
				vector<unsigned int> traits;
				covs.push_back(c);
				if(!addCovars.empty()){
					if(covarMap.find(c) == covarMap.end()){
						covs.insert(covs.end(), addCovars.begin(), addCovars.end());
					}
					else{
						for(vector<unsigned int>::iterator iter=addCovars.begin(); iter!=addCovars.end();
							iter++){
							if(*iter != c)
							{
								covs.push_back(*iter);
							}
						}
					}
				}	
				
				addCovsTraits(covs,traits,ds,cov_use,ct_filter,tempds);
				lr.calculate(model, covs, traits);			
				outputResult(lrout, lrsvout, lr, model, covs, group_iter->first, c, ds, NULL);
				delete tempds;			
			}
		}
	}
	else{
		for(int m = 0; m < (int)good_markers.size(); m++){//(int)ds->num_loci(); m++){
			Marker* mark = ds->get_locus(good_markers[m]);//ds->get_locus(m);
			if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
				if(options.doOutputSynthView())
				{
					lrsvout << mark->getRSID() << "\t" << mark->getChrom() << "\t" << mark->getBPLOC();
				}

		for(group_iter = groups.begin(); group_iter != groups.end(); group_iter++)
		{
			DataSet* tempds = getTempDataSet(ds, group_iter);

			lr.resetDataSet(tempds);

				#ifndef PLATOLIB
				_nmiss = 0;
				for(int s = 0; s < tempds->num_inds(); s++)
				{
					Sample* samp = tempds->get_sample(s);
					if(samp->isEnabled() && !samp->getAmissing(mark->getLoc()) && (samp->getPheno() == 1 || samp->getPheno() == 2))
					{
						_nmiss++;
					}
				}
				#endif
				vector<unsigned int> model;
				vector<unsigned int> covs;
				vector<unsigned int> traits;
				model.push_back(good_markers[m]);
				addCovsTraits(covs, traits, ds,cov_use,ct_filter,tempds);

				if(covs.size() == 0){
					lr.calculate(model);
				}
				else{
					lr.calculate(model, covs, traits);
				}

				outputResult(lrout, lrsvout, lr, model, covs, group_iter->first, m, ds, mark);

				delete tempds;
			}//end group_iter
			if(options.doOutputSynthView()){
				lrsvout << endl;
			}
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
		mc.calculate(_chis, tcnt);

		prev_base = 0;
		prev_chrom = -1;
		for(unsigned int m = 0; m < good_markers.size(); m++){//ds->num_loci(); m++){
			Marker* mark = ds->get_locus(good_markers[m]);//ds->get_locus(m);

			if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
				COMP << mark->toString() << "\tLOGREG\t";
				COMP << _pvals[m] << "\t"
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
			doFilter(mark, _pvals[m]);
		}
	}

	lrout.close();
}

///
/// Set covariates for use in regression models
///
void ProcessLogReg::setCovariates(vector<unsigned int>& covars){
  vector<string> use_covs = options.getCovars();

  if(options.doCovarsName()){
    // convert to indexes and add to covars
    for(unsigned int i =0; i<use_covs.size(); i++){
			int index = data_set->get_covariate_index(use_covs[i]);
			if(index != -1){
	      covars.push_back(data_set->get_covariate_index(use_covs[i]));
	    }
	    else{
	    	throw MethodException("\nERROR: " + use_covs[i] + " not found in covariate file\n\n");
	    }
    }
  }
  else{
    // convert string numbers to numbers
    for(unsigned int i =0; i<use_covs.size(); i++){
			if(use_covs.at(i).find('-') != string::npos){
				vector<string> range;
				General::Tokenize(use_covs.at(i), range, "-");
				int start = atoi(range[0].c_str())-1;
				int last = atoi(range[1].c_str())-1;
				for(int j=start; j<=last; j++){
					covars.push_back(j);
				}
			}
			else{
	      covars.push_back(atoi(use_covs[i].c_str())-1);
	    }
    }
  }
}



#ifdef PLATOLIB
void ProcessLogReg::dump2db(){}
void ProcessLogReg::create_tables()
{
	Query myQuery(*db);
	myQuery.transaction();
    for(unsigned int i = 0; i < tablename.size(); i++){
        Controller::drop_table(db, tablename[i]);
    }
    headers.clear();
    tablename.clear();
    primary_table.clear();

    string tempbatch = batchname;
    for(unsigned int i = 0; i < tempbatch.size(); i++){
        if(tempbatch[i] == ' '){
            tempbatch[i] = '_';
        }
    }
    string mytablename = tempbatch + "_";
    tempbatch = name;
    for(unsigned int i = 0; i < tempbatch.size(); i++){
        if(tempbatch[i] == ' '){
            tempbatch[i] = '_';
        }
    }
    string base = mytablename + tempbatch;


    mytablename = base + "_" + getString<int>(position);
    tablename.push_back(mytablename);
    tablenicknames.push_back("");
    primary_table[mytablename].push_back(Vars::LOCUS_TABLE);

    string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
    sql += "fkey integer not null,";
    sql += "Reference_Allele varchar(10),";
    sql += "Test varchar(20),";
    sql += "odds_ratio REAL,";
    headers[mytablename].push_back("odds_ratio");
    sql += "SE REAL,";
    headers[mytablename].push_back("SE");
    sql += "L" + getString<double>(options.getCI()*100) + " REAL,";
    headers[mytablename].push_back("L" + getString<double>(options.getCI()*100));
    sql += "U" + getString<double>(options.getCI()*100) + " REAL,";
    headers[mytablename].push_back("U" + getString<double>(options.getCI()*100));
    sql += "STAT REAL,";
    headers[mytablename].push_back("STAT");
    sql += "Pvalue REAL";
    headers[mytablename].push_back("Pvalue");
    sql += ")";
    defaultinsert = "INSERT INTO " + mytablename + " (id, fkey, Reference_Allele, Test, odds_ratio, SE,";
    defaultinsert += "L" + getString<double>(options.getCI()*100) + ",";
    defaultinsert += "U" + getString<double>(options.getCI()*100) + ",";
    defaultinsert += "STAT, Pvalue) VALUES (NULL";
    Controller::execute_sql(myQuery, sql);
    myQuery.commit();
}
void ProcessLogReg::run(DataSetObject* ds)
{
	process(ds);
}
#endif

#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
