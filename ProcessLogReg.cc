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
#include <Helpers.h>
#include "Controller.h"
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif

ProcessLogReg::ProcessLogReg(string bn, int pos, Database* pdb)
{
	name = "Logistic Regression";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
}

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


void ProcessLogReg::process(DataSet* ds)
{
	data_set = ds;
	vector<int> good_markers = Helpers::findValidMarkersIndexes(data_set->get_markers(), &options);
	cout << "good_markers contains " << getString<int>(good_markers.size()) << " items \n";

	#ifdef PLATOLIB
		Query myQuery(*db);
		create_tables();
		cout << "finished creating tables" << "\n";
	#else
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

    lrout << "Chrom\trsID\tProbeID\tBPLOC\tReference_Allele\tTest\tNMISS\tBETA\tOR\tSE\tL" << getString<double>(options.getCI()*100) << "\t" << "U" << getString<double>(options.getCI()*100) << "\tSTAT\tPvalue\n";

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
	#endif

	LogisticRegression lr;
	lr.set_parameters(&options);
	ds->set_missing_covalues(-99999);
	lr.resetDataSet(ds);
	if(ds->num_covariates() == 0 && ds->num_traits() == 0)
		lr.setFullInteraction(true);

	lr.setModelType(options.getLRModelType());

	double zt = Helpers::ltqnorm(1.0 - (1.0 - options.getCI()) / 2.0);
//	int prev_base = 0;
//	int prev_chrom = -1;
	InputFilter ct_filter;
	vector<string> use_covs = options.getCovars();
	vector<string> use_traits = options.getTraits();
	if(options.doCovarsName()){
		ct_filter.add_covariate_list(&use_covs);
		ct_filter.add_covariate_filter(InputFilter::IncludeCovariateFilter);
	}
#ifdef PLATOLIB
	if(options.doTraitsName()){
		ct_filter.add_trait_list(&use_traits);
		ct_filter.add_trait_filter(InputFilter::IncludeTraitFilter);
	}
	myQuery.transaction();
#endif
	bool cov_use = true;
#ifdef PLATOLIB
	bool trait_use = true;
	if(options.doCovarsName() && !options.doTraitsName())
	{
	trait_use = false;
	}
	if(!options.doCovarsName() && options.doTraitsName())
	{
		cov_use = false;
	}
#else
	vector<double> chis(ds->num_loci(), 0);
	vector<double> pvals(ds->num_loci(), 0);
#endif

	for(int m = 0; m < (int)good_markers.size(); m++){//(int)ds->num_loci(); m++){
		Marker* mark = ds->get_locus(good_markers[m]);//ds->get_locus(m);
		if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
			#ifndef PLATOLIB
			int nmiss = 0;
			for(int s = 0; s < ds->num_inds(); s++){
				Sample* samp = ds->get_sample(s);
				if(samp->isEnabled() && !samp->getAmissing(mark->getLoc()) && (samp->getPheno() == 1 || samp->getPheno() == 2)){
					nmiss++;
				}
			}
			#endif
			vector<unsigned int> model;
			vector<unsigned int> covs;
			vector<unsigned int> traits;
			model.push_back(good_markers[m]);
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
			#ifdef PLATOLIB
			if(trait_use){
				for(int c = 0; c < ds->num_traits(); c++){
					bool use = true;
					for(int f = 0; f < ct_filter.num_trait_filters(); f++){
						use = ct_filter.run_trait_filter(f, ds->get_trait_name(c));
					}
					if(use){
						traits.push_back(c);
					}
				}
			}
			#endif
			if(covs.size() == 0){// && traits.size() == 0){
			//	cout << "on " << m << endl;
				lr.calculate(model);
			}
			else{
				cout << "Right here 8\n";
				lr.calculate(model, covs, traits);
			}
			vector<double> coefs = lr.getCoefficients();
			//cout << "coefs\n";
			vector<double> ses = lr.getCoeffStandardErr();
			//cout << "stderr\n";
			for(unsigned int c = 0; c < model.size(); c++){
				cout << "Loop iteration # " << getString<int>(c + 1) << "\n";
				#ifdef PLATOLIB
				string sql = defaultinsert;
				cout << "A\n";
				sql += "," + getString<int>(mark->getSysprobe());
				cout << "B\n";
				sql += ",'" + mark->getReferent() + "'";
				cout << "C\n";
				sql += ",'" + options.getLRModelType() + "',";
				cout << "D\n";
				sql += (isnan(exp(coefs[c])) || isinf(exp(coefs[c]))) ? "NULL," : (getString<double>(exp(coefs[c])) + ",");
				cout << "E\n";
				#else
				lrout << mark->toString() << "\t" << mark->getReferent() << "\t" << options.getLRModelType();
				lrout << "\t" << nmiss;
				lrout << "\t" << coefs[c];
				lrout << "\t" << exp(coefs[c]);
				#endif
				double se = ses[c];
				cout << "F\n";
				cout << "coefs[c]: " << getString<double>(coefs[c]) << "\n";
				cout << "se: " << getString<double>(se) << "\n";
				double Z = coefs[c] / se;
				cout << "Z: " << getString<double>(Z) << "\n";
				cout << "G\n";
				#ifdef PLATOLIB
				sql += (isnan(se) || isinf(se)) ? "NULL," : (getString<double>(se) + ",");
				cout << "H\n";
				sql += (isnan(exp(coefs[c] - zt * se)) || isinf(exp(coefs[c] - zt * se))) ? "NULL," : (getString<double>(exp(coefs[c] - zt * se)) + ",");
				cout << "I\n";
				sql += (isnan(exp(coefs[c] + zt * se)) || isinf(exp(coefs[c] + zt * se))) ? "NULL," : (getString<double>(exp(coefs[c] + zt * se)) + ",");
				cout << "J\n";
				sql += (isnan(Z) || isinf(Z)) ? "NULL," : (getString<double>(Z) + ",");
				cout << "K\n";
				#else
				lrout << "\t" << se
					<< "\t" << exp(coefs[c] - zt * se)
					<< "\t" << exp(coefs[c] + zt * se)
					<< "\t" << Z;
				#endif
				double zz = Z*Z;
				double pvalue = 1 , df = 1;
//				int code = 1, status;
				//cdfchi(&code, &p, &pvalue, &zz, &df, &status, &bound);
				//cout << "pre-p_from_chi: " << se << " : " << Z << " : " << zz << " : " << df << endl;
				#ifdef PLATOLIB
				cout << "About to call Helpers::p_from_chi(): \n";
				cout << "zz: " << getString<double>(zz) << "\n";
				cout << "df: " << getString<double>(df) << "\n";
				//having problems with zz being nan (must be > 0)
				if ((!isnan(zz)) && (zz > 0))
				{
					//not going to have problem running p_from_chi
					pvalue = Helpers::p_from_chi(zz, df);
				}
				else
				{
					//TODO:  pvalue is null or < 0...
					pvalue = 0;
				}
                cout << "L\n";
                sql += (isnan(pvalue) || isinf(pvalue)) ? "NULL" : getString<double>(pvalue);
                cout << "M\n";
				#else
				if(se > 0)
				{
					pvalue = Helpers::p_from_chi(zz, df);
					lrout << "\t" << pvalue;
				}
				else{
					lrout << "\tnan";
				}
				chis[m] = zz;
				pvals[m] = pvalue;
				#endif
//				if(covs.size() == 0 && traits.size() == 0){
//					lrout << "\t" << lr.getFullInteractionP();
//				}
//				else{
//					lrout << "\t" << lr.getOverallP();
//				}
				#ifdef PLATOLIB
				sql += ")";
				cout << "The first SQL Statement is: " << sql << "\n";
				Controller::execute_sql(myQuery, sql);
				#else
				lrout << endl;
				#endif
				//cout << "calc done " << c << endl;
			}
			int buffer = model.size();
			cout << "Right here 1\n";
			for(int c = 0; c < (int)covs.size(); c++)
			{
				cout<< "In this for loop" << "\n";
				#ifdef PLATOLIB
				string sql = defaultinsert;
				sql += "," + getString<int>(mark->getSysprobe());
				sql += ",'" + mark->getReferent() + "'";
				sql += ",'" + ds->get_covariate_name(covs[c]) + "',";
				sql += (isnan(exp(coefs[buffer + c])) || isinf(exp(coefs[buffer + c]))) ? "NULL," : (getString<double>(exp(coefs[buffer + c])) + ",");
				#else
				lrout << mark->toString() << "\t" << mark->getReferent() << "\t" << ds->get_covariate_name(covs[c]);
				lrout << "\t" << coefs[buffer + c];
				lrout << "\t" << exp(coefs[buffer + c]);
				#endif
				double se = ses[buffer + c];
				double Z = coefs[buffer + c] / se;
				#ifndef PLATOLIB
				lrout << "\t" << se
					<< "\t" << exp(coefs[buffer + c] - zt * se)
					<< "\t" << exp(coefs[buffer + c] + zt * se)
					<< "\t" << Z;
				#endif
				double zz = Z*Z;
				double pvalue, df = 1;
//				int code = 1, status;
				//cdfchi(&code, &p, &pvalue, &zz, &df, &status, &bound);
				#ifdef PLATOLIB
				cout << "About to call Helpers::p_from_chi\n";
				cout << "zz is: " << getString<double>(zz) << "\n";
				cout << "df is: " << getString<double>(df) << "\n";
				pvalue = Helpers::p_from_chi(zz, df);
				sql += (isnan(pvalue) || isinf(pvalue)) ? "NULL" : getString<double>(pvalue);
				#else
				if(se > 0){
					pvalue = Helpers::p_from_chi(zz, df);
					lrout << "\t" << pvalue;
				}
				else{
					lrout << "\tnan";
				}
				#endif
//				if(covs.size() == 0 && traits.size() == 0){
//					lrout << "\t" << lr.getFullInteractionP();
//				}
//				else{
//					lrout << "\t" << lr.getOverallP();
//				}
				#ifdef PLATOLIB
				sql += ")";
				cout << "The Second SQL Satement is: " << sql;
				Controller::execute_sql(myQuery, sql);
				#else
				lrout << endl;
				#endif
			}
			#ifdef PLATOLIB
			cout << "Right here 2\n";
			buffer += covs.size();
			for(int c = 0; c < (int)traits.size(); c++)
			{
				string sql = defaultinsert;
				sql += "," + getString<int>(mark->getSysprobe());
				sql += ",'" + mark->getReferent() + "'";
				sql += ",'" + ds->get_trait_name(traits[c]) + "',";
				sql += (isnan(exp(coefs[buffer + c])) || isinf(exp(coefs[buffer + c]))) ? "NULL," : (getString<double>(exp(coefs[buffer + c])) + ",");
				double se = ses[buffer + c];
				double Z = coefs[buffer + c] / se;
				sql += (isnan(se) || isinf(se)) ? "NULL," : (getString<double>(se) + ",");
				sql += (isnan(exp(coefs[buffer + c] - zt * se)) || isinf(exp(coefs[buffer + c] - zt * se))) ? "NULL," : (getString<double>(exp(coefs[buffer + c] - zt * se)) + ",");
				sql += (isnan(exp(coefs[buffer + c] + zt * se)) || isinf(exp(coefs[buffer + c] + zt * se))) ? "NULL," : (getString<double>(exp(coefs[buffer + c] + zt * se)) + ",");
				sql += (isnan(Z) || isinf(Z)) ? "NULL," : (getString<double>(Z) + ",");
				//double zz = Z*Z;
				//double pvalue, df = 1;
				//int code = 1, status;
				//pvalue = Helpers::p_from_chi(zz, df);
				//lrout << "\t" << pvalue;
				//if(covs.size() == 0 && traits.size() == 0){
				//	lrout << "\t" << lr.getFullInteractionP();
				//}
				//else{
				//	lrout << "\t" << lr.getOverallP();
				//}
				//lrout << endl;
				sql += ")";
				cout << "The Third SQL Statement is: " << sql << "\n";
				Controller::execute_sql(myQuery, sql);
			}
			#endif
			cout << "Right here 3\n";
		}
		cout << "Right here 4\n";
	}
	cout << "Right here 5\n";
	#ifndef PLATOLIB
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
#endif
cout << "Right here 6\n";
#ifdef PLATOLIB
	myQuery.commit();
	cout << "Right here 7\n";
	hasresults = true;
#else
	lrout.close();
#endif
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
