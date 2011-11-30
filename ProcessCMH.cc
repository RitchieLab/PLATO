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
#include "ProcessCMH.h"
#include <Options.h>
#include <General.h>
#include <Helpers.h>
using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

string ProcessCMH::stepname = "cmh";
#ifdef PLATOLIB
ProcessCMH::ProcessCMH(string bn, int pos, Database* pdb)
{
	name = "CMH";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
}
#endif

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

void ProcessCMH::filter(){}

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

#ifdef PLATOLIB
	options.readClustersFromString(data_set->get_samples());
	cmh.setOverwrite(true);
	#else
	options.readClusters(data_set->get_samples());
	cmh.setOverwrite(overwrite);
	#endif


	cmh.setOptions(options);
	cmh.setRank(rank);
	cmh.resetDataSet(data_set);

	#ifndef PLATOLIB
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
	#else
	create_tables();
	Query myQuery(*db);
	myQuery.transaction();

	#endif
	int prev_base = 0;
	int prev_chrom = -1;
	vector<Marker*> good_markers = Helpers::findValidMarkers(ds->get_markers(), &options);
	msize = good_markers.size();

	for(int m = 0; m < (int)msize; m++)
	{
		Marker* mark = good_markers[m];//data_set->get_locus(m);

		if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_base, prev_chrom)){
			cmh.calculate(mark);
			#ifdef PLATOLIB
			string sql = defaultinsert;
			sql += "," + getString<int>(mark->getSysprobe());
			#else
			MHOUT << mark->getChrom() << "\t" << mark->getRSID();
			#endif
			if(options.getCMH2x2xK()){
				#ifdef PLATOLIB

				#else
				MHOUT << "\t" << mark->getAllele1()
						<< "\t" << mark->getAllele2() << "\t"
						<< mark->getBPLOC() << "\t";

				chis[m] = cmh.get_chisq();
				pvals[m] = cmh.get_pval();
				#endif

				if (Helpers::realnum(cmh.get_chisq()))
				{
				#ifdef PLATOLIB
					sql += "," + getString<double>(cmh.get_chisq());
					sql += "," + getString<double>(cmh.get_pval());
					#else
					MHOUT << cmh.get_chisq() << "\t" << cmh.get_pval()
							<< "\t";
					#endif
				}
				else{
					#ifdef PLATOLIB
					sql += ", NULL, NULL";

					#else
					MHOUT << "NA" << "\t" << "NA" << "\t";

					#endif
				}
				if (Helpers::realnum(cmh.getOR()))
				{
					#ifdef PLATOLIB
					sql += "," + getString<double>(cmh.getOR());

					#else
					MHOUT << cmh.getOR() << "\t";

					#endif
				}
				else{
					#ifdef PLATOLIB
					sql += ",NULL";

					#else
					MHOUT << "NA" << "\t";

					#endif
				}
				if (Helpers::realnum(cmh.getOR_lower()))
				{
					#ifdef PLATOLIB
					sql += "," + getString<double>(cmh.getOR_lower());
					#else
					MHOUT << cmh.getOR_lower() << "\t";
					#endif
				}
				else
				{
					#ifdef PLATOLIB
					sql += ",NULL";
					#else
					MHOUT << "NA" << "\t ";
					#endif
				}
				if (Helpers::realnum(cmh.getOR_upper()))
				{
					#ifdef PLATOLIB
					sql += "," + getString<double>(cmh.getOR_upper());
					#else
					MHOUT << cmh.getOR_upper();
					#endif
				}
				else
				{
					#ifdef PLATOLIB
					sql += ",NULL";
					#else
					MHOUT << "NA";
					#endif
				}

				if(options.doBreslowDay()){
					if (cmh.get_pvalbd() > -1)
					{
						#ifdef PLATOLIB
						sql += "," + getString<double>(cmh.get_chisqbd());
						sql += "," + getString<double>(cmh.get_pvalbd());
						#else
						MHOUT << "\t" << cmh.get_chisqbd() << "\t" << cmh.get_pvalbd();
						#endif
					}
					else
					{
						#ifdef PLATOLIB
						sql += ",NULL,NULL";
						#else
						MHOUT << "\t" << "NA" << "\t" << "NA";
						#endif
					}
				}
			}
			else
			{
				#ifdef PLATOLIB
				sql += "," + getString<double>(cmh.get_chisq());
				sql += "," + getString<double>(cmh.get_pval());
				#else
				MHOUT << "\t" << cmh.get_chisq() << "\t" << cmh.get_pval();
				chis[m] = cmh.get_chisq();
				pvals[m] = cmh.get_pval();
				#endif
			}
			#ifdef PLATOLIB
			sql += ")";
			Controller::execute_sql(myQuery, sql);
			#else
			MHOUT << endl;
			#endif

//			doFilter(mark, cmh.get_pval());
		}
	}
	#ifdef PLATOLIB
	myQuery.commit();
	#endif

	#ifndef PLATOLIB
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
	#endif

	#ifdef PLATOLIB
	hasresults = true;
	#endif
}//end method process()

#ifdef PLATOLIB
void ProcessCMH::dump2db(){}

void ProcessCMH::create_tables()
{
	Query myQuery(*db);

    for(int i = 0; i < (int)tablename.size(); i++){
        Controller::drop_table(db, tablename[i]);
    }
    headers.clear();
    tablename.clear();
    primary_table.clear();

    string tempbatch = batchname;
    for(int i = 0; i < (int)tempbatch.size(); i++){
        if(tempbatch[i] == ' '){
            tempbatch[i] = '_';
        }
    }
    string mytablename = tempbatch + "_";
    tempbatch = name;
    for(int i = 0; i < (int)tempbatch.size(); i++){
        if(tempbatch[i] == ' '){
            tempbatch[i] = '_';
        }
    }

    mytablename += tempbatch + "_" + getString<int>(position);
    tablename.push_back(mytablename);
    tablenicknames.push_back("");
    primary_table[mytablename].push_back(Vars::LOCUS_TABLE);


    string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
    sql += "fkey integer not null,";
    sql += "CHISQ REAL,";
    headers[mytablename].push_back("CHISQ");
    sql += "Pvalue REAL";
    headers[mytablename].push_back("Pvalue");
    if(options.getCMH2x2xK()){
        sql += ", odds_ratio REAL,";
        headers[mytablename].push_back("odds_ratio");
        sql += "L" + getString<double>(options.getCI() * 100.0) + " REAL,";
        headers[mytablename].push_back("L" + getString<double>(options.getCI() * 100.0));
        sql += "U" + getString<double>(options.getCI() * 100.0) + " REAL";
        headers[mytablename].push_back("U" + getString<double>(options.getCI() * 100.0));
        if(options.doBreslowDay()){
            sql += ",CHISQ_BD REAL,";
            headers[mytablename].push_back("CHISQ_BD");
            sql += "Pvalue_BD REAL";
            headers[mytablename].push_back("Pvalue_BD");
        }
    }
    sql += ")";
    defaultinsert = "INSERT INTO " + mytablename + " (id, fkey, CHISQ, Pvalue";
    if(options.getCMH2x2xK()){
        defaultinsert += ", odds_ratio, ";
        defaultinsert += "L" + getString<double>(options.getCI() * 100.0) + ", ";
        defaultinsert += "U" + getString<double>(options.getCI() * 100.0);
        if(options.doBreslowDay()){
            defaultinsert += ", CHISQ_BD, Pvalue_BD";
        }
    }
    defaultinsert += ") VALUES (NULL";

    myQuery.transaction();
    Controller::execute_sql(myQuery, sql);
    myQuery.commit();
}//end method create_tables()

void ProcessCMH::run(DataSetObject* ds)
{
	//merged the code from the viewer's run method into the process method above
	cout << "In Run method, about to call Process() \n";
	process(ds);
	cout << "Finished process() method \n";
}//end method run()
#endif

#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
