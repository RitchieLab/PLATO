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
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include <LinearRegression.h>
#include <LinRegression.h>
#include <MultComparison.h>
#include "ProcessLinearReg.h"
#include <Options.h>
#include <General.h>
#include <Helpers.h>
#ifdef PLATOLIB
#include "Controller.h"
#endif
using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif

#ifdef PLATOLIB
ProcessLinearReg::ProcessLinearReg(string bn, int pos, Database* pdb)
{
	name = "Linear Regression";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
}
#endif

void ProcessLinearReg::FilterSummary(){

	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessLinearReg::PrintSummary(){
	int msize = data_set->num_loci();
	for(int m = 0; m < msize; m++){
		data_set->get_locus(m)->setFlag(false);
	}

}

void ProcessLinearReg::filter()
{
#ifdef PLATOLIB
	for (int m = 0; m < (int)data_set->num_loci(); m++)
	{
		data_set->get_locus(m)->setFlag(false);
	}
#endif
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


DataSet* ProcessLinearReg::getTempDataSet(DataSet* ds, map<string, vector<Sample*> >::iterator group_iter){
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


void ProcessLinearReg::addCovsTraits(vector<unsigned int>& covs, vector<unsigned int>& traits,
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


///
/// Set covariates for use in regression models
///
void ProcessLinearReg::setCovariates(vector<unsigned int>& covars){
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


void ProcessLinearReg::processEnvironmental(ostream& lrout)
{

	if(options.doGroupFile()){
		opts::printLog("Reading group information [" + options.getGroupFile() + "]\n");
		options.readGroups(data_set->get_samples());
	}

    lrout.precision(4);

 		lrout << "Variable\t";
	  if(options.doGroupFile()){
  	  lrout << "GRP\t";
    }
    lrout << "Test\tNMISS\tBETA\tOR\tSE\tL" << getString<double>(options.getCI()*100) 
    	<< "\t" << "U" << getString<double>(options.getCI()*100) << "\tSTAT\tPvalue\n";


	LinRegression lr;
	lr.resetDataSet(data_set);
	lr.set_parameters(&options);
	data_set->set_missing_covalues(-99999);
	map<string, vector<Sample*> > groups = options.getGroups();
	map<string, vector<Sample*> >::iterator group_iter;
	if(groups.size() == 0){
		groups["GROUP_1"] = *(data_set->get_samples());
	}
	
	int prev_base = 0;
	int prev_chrom = -1;

	lr.setModelType(options.getLRModelType());

	double zt = Helpers::ltqnorm(1.0 - (1.0 - options.getCI()) / 2.0);
	InputFilter ct_filter;

	vector<string> use_covs = options.getCovars();
	vector<string> use_traits = options.getTraits();
	if(options.doCovarsName()){
		ct_filter.add_covariate_list(&use_covs);
		ct_filter.add_covariate_filter(InputFilter::IncludeCovariateFilter);
	}
	bool cov_use = false;

	int numCovars = data_set->num_covariates();
	vector<unsigned int> addCovars;
	setCovariates(addCovars);
	set<unsigned int> covarMap;
	vector<unsigned int>::iterator iter;
	for(iter = addCovars.begin(); iter != addCovars.end(); iter++){
		covarMap.insert(*iter);
	}
		
	vector<string>* get_covariates();
	int nmiss;
	double missingCoValue = data_set->get_missing_covalue();
	for(int c=0; c < numCovars; c++){
		for(group_iter = groups.begin(); group_iter != groups.end(); group_iter++)
		{
			DataSet* tempds = getTempDataSet(data_set, group_iter);
			lr.resetDataSet(tempds);				
			nmiss = 0;
			for(int s = 0; s < tempds->num_inds(); s++)
			{
				Sample* samp = tempds->get_sample(s);
				if(samp->isEnabled() && !samp->getCovariate(c) != missingCoValue && (samp->getPheno() == 1 || samp->getPheno() == 2))
				{
					nmiss++;
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
				
			addCovsTraits(covs,traits,data_set,cov_use,ct_filter,tempds);
			lr.calculate(model, covs);			
			
			vector<double> coefs = lr.getCoefficients();
			vector<double> ses = lr.getCoeffStandardErr();
			
			for(int c=0; c < (int)covs.size(); c++)
			{
				lrout << data_set->get_covariate_name(covs[0]) << "\t";
				if(options.doGroupFile()){
					lrout <<  group_iter->first << "\t";
				}
				lrout << data_set->get_covariate_name(covs[c]);
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
				double pvalue, df = 1;
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
				lrout << "\n";				
			}
			
			lrout << data_set->get_covariate_name(covs[0]) << "\t";	
			if(options.doGroupFile()){
					lrout << group_iter->first << "\t";
			}
			lrout << "overall";
			lrout << "\t" << nmiss;
			lrout << "\t-----";
			lrout << "\t-----";
			lrout << "\t-----\t-----\t-----\t";
			lrout << lr.getOverallScore();
			lrout << "\t" << lr.getOverallP() << endl;			
			
			delete tempds;			
		}
	}

}


void ProcessLinearReg::process(DataSet* ds){
	data_set = ds;

	#ifdef PLATOLIB
		create_tables();
	#endif

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

#ifndef PLATOLIB
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
#endif

	// run only environmental variables (in covariate file)
	if(options.runCovarsOnly()){
		processEnvironmental(lrout);
		lrout.close();
		return;
	}
	
#ifndef PLATOLIB
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
#endif

	LinearRegression lr;
	lr.setOptions(options);
#ifdef PLATOLIB
	ds->set_missing_covalues(-99999);
#endif

	//MEMORY_LEAK
	if(groups.size() == 0){
		groups["GROUP_1"] = *(ds->get_samples());
	}
#ifndef PLATOLIB
	if(options.doOutputSynthView()){
		for(group_iter = groups.begin(); group_iter != groups.end(); group_iter++){
			lrsvout << "\t" << group_iter->first << ":pval" << "\t" << group_iter->first << ":beta" << "\t" << group_iter->first << ":N";
		}
		lrsvout << endl;
	}
#endif
	//lr.resetDataSet(ds);

	//END_MEMORY_LEAK

	vector<double> chis(ds->num_loci(), 0);
	vector<double> main_pvals(ds->num_loci(), 0);

	int prev_base = 0;
	int prev_chrom = -1;
	int msize = good_markers.size();
	#ifdef PLATOLIB
		Query myQuery(*db);
		myQuery.transaction();
	#endif

	for(int i = 0; i < msize; i++){
		Marker* mark = good_markers[i];//ds->get_locus(i);
		if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
#ifndef PLATOLIB
			//MEMORY_LEAK
			if(options.doOutputSynthView()){
				lrsvout << mark->getRSID() << "\t" << mark->getChrom() << "\t" << mark->getBPLOC();
			}
#endif
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
					#ifdef PLATOLIB
					string sql = defaultinsert;
					sql += "," + getString<int>(mark->getSysprobe());
					sql += ",'" + mark->getReferent() + "'";
					sql += ",'" + labels[l] + "'";
					sql += "," + getString<int>(lr.getCalcMissing());
					sql += ",";
					sql += (isnan(coefs[l]) || isinf(coefs[l])) ? "NULL," : (getString<double>(coefs[l]) + ",");
					sql += (isnan(zs[l]) || isinf(zs[l])) ? "NULL," : (getString<double>(zs[l]) + ",");
					sql += (isnan(pvals[l]) || isinf(pvals[l])) ? "NULL," : (getString<double>(pvals[l]));
					sql += ")";
					cout << "SQL: " << sql << endl;
					Controller::execute_sql(myQuery, sql);
					#else
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
						std::vector<string>::iterator covIterator = std::find_if((*covariatesList).begin(), (*covariatesList).end(), FindString(labels[l]));
						if(covIterator == (*covariatesList).end())
						{
							lrsvout << "\t" << pvals[l] << "\t" << coefs[l] << "\t" << lr.getCalcMissing();
						}
					}
					#endif
				}//end loop through labels
				#ifndef PLATOLIB
				chis[i] = lr.getStatistic();
				main_pvals[i] = pvals[1];
				#endif
				delete tempds;
			}//end loop through Groups
#ifndef PLATOLIB
			if(options.doOutputSynthView()){
				lrsvout << endl;
			}
#endif
		}
		#ifdef PLATOLIB
				hasresults = true;
		#endif
//		lrout << options.getLinRModelType() << "\t";
//		lrout << lr.getCalcMissing() << "\t" << lr.getTestCoef() << "\t" << lr.getZ() << "\t" << lr.getPValue() << endl;
	}//end loop through Markers
	#ifdef PLATOLIB
				myQuery.commit();
	#endif

	#ifndef PLATOLIB
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
	#endif
}

#ifdef PLATOLIB
void ProcessLinearReg::dump2db(){}
void ProcessLinearReg::create_tables()
{
	Query myQuery(*db);
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
    sql += "Nmiss integer,";
    headers[mytablename].push_back("Nmiss");
    sql += "Beta REAL,";
    headers[mytablename].push_back("Beta");
    sql += "STAT REAL,";
    headers[mytablename].push_back("STAT");
    sql += "Pvalue REAL";
    headers[mytablename].push_back("Pvalue");
    sql += ")";
    defaultinsert = "INSERT INTO " + mytablename + " (id, fkey, Reference_Allele, Test, Nmiss, Beta, STAT, Pvalue) VALUES (NULL";
    myQuery.transaction();
    Controller::execute_sql(myQuery, sql);
    myQuery.commit();
}//end method create_tables()

void ProcessLinearReg::run(DataSetObject* ds)
{
	process(ds);
}

#endif

#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
