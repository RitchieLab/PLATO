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
#include <LinearRegression.h>
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


void ProcessLinearReg::process(DataSet* ds){
	data_set = ds;

	#ifdef PLATOLIB
		create_tables();
	#endif

	vector<Marker*> good_markers = Helpers::findValidMarkers(data_set->get_markers(), &options);

	//check if new covariate file is listed...or covariate name.
	//create vector of covariate indexes to use if specified.

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

    lrout << "Chrom\trsID\tProbeID\tBPLOC\tReference_Allele\tTest\tNMISS\tBeta\tOR\tSE\tSTAT\tPvalue\n";
#endif

	LinearRegression lr;
	lr.setOptions(options);
	lr.resetDataSet(ds);
#ifdef PLATOLIB
	ds->set_missing_covalues(-99999);
#endif

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
				lrout << mark->getBPLOC() << "\t" << mark->getReferent() << "\t";
				lrout << labels[l] << "\t" << lr.getCalcMissing() << "\t" << coefs[l] << "\t" << exp(coefs[l]) << "\t" << se << "\t" << zs[l] << "\t" << pvals[l] << endl;
				#endif
			}
			#ifndef PLATOLIB
			chis[i] = lr.getStatistic();
			main_pvals[i] = pvals[1];
			#endif
		}
		#ifdef PLATOLIB
				hasresults = true;
		#endif
//		lrout << options.getLinRModelType() << "\t";
//		lrout << lr.getCalcMissing() << "\t" << lr.getTestCoef() << "\t" << lr.getZ() << "\t" << lr.getPValue() << endl;
	}
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
