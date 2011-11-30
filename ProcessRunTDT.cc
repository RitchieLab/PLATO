#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <unistd.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#include <fenv.h>
#include <algorithm>
#include <MultComparison.h>
#include "ProcessRunTDT.h"
#include "Chrom.h"
#include <General.h>
//#include "ChiSquare.h"
#include <Helpers.h>
#include <cdflib.h>
#ifdef PLATOLIB
#include "Controller.h"
#endif

using namespace std;
using namespace Methods;
#ifdef PLATOLIB
namespace PlatoLib
{
#endif

string ProcessRunTDT::stepname = "tdt";
#ifdef PLATOLIB
ProcessRunTDT::ProcessRunTDT(string bn, int pos, Database* pdb)
{
    name = "TDT";
    batchname = bn;
    position = pos;
    hasresults = false;
    db = pdb;
}
#endif

void ProcessRunTDT::PrintSummary(){
	string fname1 = opts::_OUTPREFIX_ + "tdt" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	ofstream output (fname1.c_str(), ios::out);
	if(!output){
		opts::printLog("Error opening " + fname1 + ". Exiting!\n");
		throw MethodException("");
	}
	opts::addFile("Marker", stepname, fname1);
	output.precision(4);
	output << "Chrom\trsID\tProbeID\tbploc";
	if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
		output << "\t" << data_set->get_locus(0)->getDetailHeaders();
	}
	output << "\tChi_square\tTDT_pvalue\tTDT_neglog(pvalue)\tNum_Fams\tT:U\tA1:A2\tOR\tL" + getString<double>(options.getCI()*100) + "\t" + "U" + getString<double>(options.getCI()*100) << endl;
	opts::addHeader(fname1, "Chi_square");
	opts::addHeader(fname1, "TDT_pvalue");
	opts::addHeader(fname1, "TDT_neglog(pvalue)");
	opts::addHeader(fname1, "Num_Fams");
	opts::addHeader(fname1, "T:U");
	opts::addHeader(fname1, "A1:A2");
	opts::addHeader(fname1, "OR");
	opts::addHeader(fname1, "L" + getString<double>(options.getCI()*100));
	opts::addHeader(fname1, "U" + getString<double>(options.getCI()*100));

	ofstream svout;
	ofstream groupout;
	map<string, vector<Sample*> > groups = options.getGroups();
	map<string, vector<Sample*> >::iterator group_iter;

	if(options.doGroupFile()){
		string fnameg = opts::_OUTPREFIX_ + "tdt_groups" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if(!overwrite){
			fnameg += "." + getString<int>(order);
		}
		groupout.open(fnameg.c_str());
		if(!groupout){
			opts::printLog("Error opening " + fnameg + ". Exiting!\n");
			throw MethodException("");
		}
//		opts::addFile("Marker", stepname, fname1);
		groupout.precision(4);
		groupout << "Chrom\trsID\tProbeID\tbploc";
		if(data_set->get_locus(0)->getDetailHeaders().size() > 0){
			groupout << "\t" << data_set->get_locus(0)->getDetailHeaders();
		}
		groupout << "\tGRP\tChi_square\tTDT_pvalue\tTDT_neglog(pvalue)\tNum_Fams\tT:U\tA1:A2\tOR\tL" + getString<double>(options.getCI()*100) + "\t" + "U" + getString<double>(options.getCI()*100) << endl;

		if(options.doOutputSynthView()){
			string fnamesv = opts::_OUTPREFIX_ + "tdt_synthview" + options.getOut() + ".txt";
			if(!overwrite){
				fnamesv += "." + getString<int>(order);
			}
			svout.open(fnamesv.c_str());
			if(!svout){
				opts::printLog("Error opening " + fnamesv + ". Exiting!\n");
				throw MethodException("");
			}
			svout.precision(4);
	    	svout << "SNP\tChromosome\tLocation";

	    	for(group_iter = groups.begin(); group_iter != groups.end(); group_iter++){
	    		svout << "\t" << group_iter->first << ":pval";
	    	}
	    	svout << endl;
		}
	}
	ofstream compout;
	if(options.doMultCompare()){
		string fname2 = opts::_OUTPREFIX_ + "tdt_comparisons" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if(!overwrite){
			fname2 += "." + getString<int>(order);
		}
		compout.open(fname2.c_str());
		if(!compout){
			opts::printLog("Error opening " + fname2 + ". Exiting!\n");
			throw MethodException("");
		}
		compout.precision(4);
		compout << "Chrom\trsID\tProbeID\tbploc";
		compout  << "\tCALC"
			  << "\tOriginal_Pval"
			  << "\tGC"
			  << "\tBONF"
			  << "\tHOLM"
			  << "\tSIDAK_SS"
			  << "\tSIDAK_SD"
			  << "\tFDR_BH"
			  << "\tFDR_BY"
			  << endl;
		opts::addHeader(fname2, "CALC");
		opts::addHeader(fname2, "Original_Pval");
		opts::addHeader(fname2, "GC");
		opts::addHeader(fname2, "BONF");
		opts::addHeader(fname2, "HOLM");
		opts::addHeader(fname2, "SIDAK_SS");
		opts::addHeader(fname2, "SIDAK_SD");
		opts::addHeader(fname2, "FDR_BH");
		opts::addHeader(fname2, "FDR_BY");

	}

	int msize = good_markers.size();//data_set->num_loci();

	double zt = Helpers::ltqnorm(1 - (1 - options.getCI()) / 2);

//	int prev_base = 0;
//	int prev_chrom = -1;
	for(int i = 0; i < msize; i++){
		if(good_markers[i]->isEnabled()){//data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged() && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
			output << good_markers[i]->toString() << "\t";
				//<< ((float)(*markers)[i]->getBPLOC()/1000000.0f) << "\t"
			output << chi[i] << "\t"
				<< pval[i] << "\t"
				<< (double)abs(log10(pval[i])) << "\t";
			//output.precision(4);
			output << fams_used[i] << "\t"
				<< trans[i] << ":" << untrans[i] << "\t";
			if(!good_markers[i]->isMicroSat()){
				output << good_markers[i]->getAllele1() << ":" << good_markers[i]->getAllele2();
			}
			else{
				output << "NA";
			}
			double OR = (double)trans[i] / (double)untrans[i];
			output << "\t" << OR;
			double OR_lower = exp(log(OR) - zt * sqrt(1/trans[i] + 1/untrans[i]));
			double OR_upper = exp(log(OR) + zt * sqrt(1/trans[i] + 1/untrans[i]));
			output << "\t" << OR_lower << "\t" << OR_upper;

				output << endl;
				if(options.doOutputSynthView()){
					svout << good_markers[i]->getRSID() << "\t" << good_markers[i]->getChrom() << "\t" << good_markers[i]->getBPLOC();
				}
				if(options.doGroupFile()){
					int gcount = 0;
					for(group_iter = groups.begin(); group_iter != groups.end(); group_iter++){
						groupout << good_markers[i]->toString() << "\t";
						groupout << group_iter->first << "\t";
						groupout << gchi[i][gcount] << "\t"
						<< gpval[i][gcount] << "\t"
						<< (double)abs(log10(gpval[i][gcount])) << "\t"
						<< gfams_used[i][gcount] << "\t"
						<< gtrans[i][gcount] << ":"
						<< guntrans[i][gcount] << "\t";
						if(!good_markers[i]->isMicroSat()){
							groupout << good_markers[i]->getAllele1() << ":" << good_markers[i]->getAllele2();
						}
						else{
							groupout << "NA";
						}
						OR = (double)gtrans[i][gcount] / (double)guntrans[i][gcount];
						groupout << "\t" << OR;
						OR_lower = exp(log(OR) - zt * sqrt(1/gtrans[i][gcount] + 1/guntrans[i][gcount]));
						OR_upper = exp(log(OR) + zt * sqrt(1/gtrans[i][gcount] + 1/guntrans[i][gcount]));
						groupout << "\t" << OR_lower << "\t" << OR_upper;
						groupout << endl;
						if(options.doOutputSynthView()){
							svout << "\t" << gpval[i][gcount];
						}
						gcount++;
					}
				}
				if(options.doOutputSynthView()){
					svout << endl;
				}
/*
			if(options.doMultCompare()){
				MultComparison mc(options);
				vector<double> chivals;
				chivals.push_back(chi[i]);
				vector<int> tcnt;
				mc.calculate(chivals, tcnt);
				compout << data_set->get_locus(i)->toString() << "\t"
					<< "TDT\t"
					<< pval[i] << "\t"
					<< mc.get_genomic_control(0) << "\t"
					<< mc.get_bonferroni(0) << "\t"
					<< mc.get_holm(0) << "\t"
					<< mc.get_sidak_single_step(0) << "\t"
					<< mc.get_sidak_step_down(0) << "\t"
					<< mc.get_fdr_bh(0) << "\t"
					<< mc.get_fdr_by(0)
					<< endl;
			}
*/
				//output.precision(100);
		}
		good_markers[i]->setFlag(false);
	}
	if(options.doMultCompare()){
		MultComparison mc(options);
		vector<int> tcnt;
		mc.calculate(chi, tcnt);
		for(int i = 0; i < (int)good_markers.size(); i++){
		compout << good_markers[i]->toString() << "\t"
			<< "TDT\t"
			<< pval[i] << "\t"
			<< mc.get_genomic_control(i) << "\t"
			<< mc.get_bonferroni(i) << "\t"
			<< mc.get_holm(i) << "\t"
			<< mc.get_sidak_single_step(i) << "\t"
			<< mc.get_sidak_step_down(i) << "\t"
			<< mc.get_fdr_bh(i) << "\t"
			<< mc.get_fdr_by(i)
			<< endl;
		}
	}

	if(output.is_open()){
		output.close();
	}
	if(compout.is_open()){
		compout.close();
	}
	if(groupout.is_open()){
		groupout.close();
	}
	if(svout.is_open()){
		svout.close();
	}
}

void ProcessRunTDT::filter()
{
	#ifndef PLATOLIB
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int size = good_markers.size();//data_set->num_loci();
//		int prev_base = 0;
//		int prev_chrom = -1;
		for(int i = 0; i < size; i++){
			if(good_markers[i]->isEnabled()){//data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged() && isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
				bool inc = false;
				if(options.doThreshMarkersHigh() && Helpers::dGreater(pval[i], options.getThreshMarkersHigh())){
					good_markers[i]->setEnabled(false);//data_set->get_locus(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersLow() && Helpers::dLess(pval[i],options.getThreshMarkersLow())){
					good_markers[i]->setEnabled(false);//data_set->get_locus(i)->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}
		}
	}
	#endif
}//end method filter()

void ProcessRunTDT::setThreshold(string thresh){
	options.setUp(thresh);
}

void ProcessRunTDT::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
	    getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
    opts::_MARKERS_WORKING_ -= orig_num_markers;
}

void ProcessRunTDT::process(DataSet* ds)
{
	data_set = ds;
	good_markers = Helpers::findValidMarkers(data_set->get_markers(), &options);

	RunTDT tdt(data_set);
	tdt.setOptions(options);
	if(options.doGroupFile()){
		options.readGroups(data_set->get_samples());
	}

	int msize = good_markers.size();//data_set->num_loci();
	chi.resize(msize);
	pval.resize(msize);
	fams_used.resize(msize);
	maf.resize(msize);
	trans.resize(msize);
	untrans.resize(msize);

	gpval.resize(msize);
	gchi.resize(msize);
	gfams_used.resize(msize);
	gmaf.resize(msize);
	gtrans.resize(msize);
	guntrans.resize(msize);
////	int ssize = data_set->num_inds();

//	int prev_base = 0;
//	int prev_chrom = -1;
	for(int m = 0; m < msize; m++)
	{
		if(good_markers[m]->isEnabled())
		{//data_set->get_locus(m)->isEnabled() && isValidMarker(data_set->get_locus(m), &options, prev_base, prev_chrom)){
				tdt.calculate(good_markers[m]);

				pval[m] = tdt.getPval();
				chi[m] = tdt.getChi();
				fams_used[m] = tdt.getFamsUsed();
				maf[m] = -1;
				trans[m] = tdt.getTransmitted();
				untrans[m] = tdt.getUntransmitted();
			if(options.doGroupFile()){
				map<string, vector<Sample*> > groups = options.getGroups();
				map<string, vector<Sample*> >::iterator group_iter;
				for(group_iter = groups.begin(); group_iter != groups.end(); group_iter++){
					tdt.calculate(good_markers[m], group_iter->second);
					gpval[m].push_back(tdt.getPval());
					gchi[m].push_back(tdt.getChi());
					gfams_used[m].push_back(tdt.getFamsUsed());
					gmaf[m].push_back(-1);
					gtrans[m].push_back(tdt.getTransmitted());
					guntrans[m].push_back(tdt.getUntransmitted());
				//	cout << group_iter->first << ":" << gpval[m][0] << ":" << gfams_used[m][0] << endl;
				}
			}
		}
	}
}//end method process(DataSet* ds)

#ifdef PLATOLIB
void ProcessRunTDT::create_tables()
{
	Query myQuery(*db);
	myQuery.transaction();
    for(int i = 0; i < (int)tablename.size(); i++){
        Controller::drop_table(db, tablename[i]);
    }
    headers.clear();
    primary_table.clear();
    tablename.clear();

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

    string sql = "CREATE TABLE " + tablename.at(0) + " (id integer primary key,";
    sql += "fkey integer not null,";
    sql += "chi_square REAL,";
    headers[tablename.at(0)].push_back("chi_square");
    sql += "pvalue REAL,";
    headers[tablename.at(0)].push_back("pvalue");
    sql += "neglog_pvalue REAL,";
    headers[tablename.at(0)].push_back("neglog_pvalue");
    sql += "num_fams integer,";
    headers[tablename.at(0)].push_back("num_fams");
    sql += "transmitted integer,";
    headers[tablename.at(0)].push_back("transmitted");
    sql += "untransmitted integer,";
    headers[tablename.at(0)].push_back("untransmitted");
    sql += "a1 varchar(10),";
    sql += "a2 varchar(10),";
    sql += "odds_ratio REAL,";
    headers[tablename.at(0)].push_back("odds_ratio");
    sql += "L" + getString<double>(options.getCI() * 100) + " REAL,";
    headers[tablename.at(0)].push_back("L" + getString<double>(options.getCI() * 100));
    sql += "U" + getString<double>(options.getCI() * 100) + " REAL)";
    headers[tablename.at(0)].push_back("U" + getString<double>(options.getCI() * 100));
    cout << sql << endl;

    Controller::execute_sql(myQuery, sql);
    myQuery.commit();
}//end method create_tables()

void ProcessRunTDT::dump2db()
{
    create_tables();

    string sql = "";

    int msize = data_set->num_loci();

    double zt = Helpers::ltqnorm(1 - (1 - options.getCI()) / 2);//Methods::ltqnorm(1 - (1 - options.getCI()) / 2);

    int prev_base = 0;
    int prev_chrom = -1;
    cout << "START TDT INSERT" << endl;
//    sql = "BEGIN;\n";

    Query myQuery(*db);
    myQuery.transaction();

    for(int i = 0; i < msize; i++){
        if(data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged() && Helpers::isValidMarker(data_set->get_locus(i), &options, prev_base, prev_chrom)){
            sql = "INSERT INTO " + tablename.at(0) + " (id, fkey, chi_square, pvalue, neglog_pvalue, num_fams, transmitted, untransmitted, a1, a2, odds_ratio, ";
            sql += "L" + getString<double>(options.getCI() * 100) + ",";
            sql += "U" + getString<double>(options.getCI() * 100) + ") VALUES (NULL, ";
            sql += getString<int>(data_set->get_locus(i)->getSysprobe()) + ",";
            //sql += "(SELECT id from LOCI WHERE name = '" + data_set->get_locus(i)->getRSID() + "' ";
            //sql += "AND chrom = '" + getString<int>(data_set->get_locus(i)->getChrom()) + "' ";
            //sql += "AND bploc = '" + getString<int>(data_set->get_locus(i)->getBPLOC()) + "'), ";
            sql += (isnan(chi[i]) || isinf(chi[i])) ? "NULL," : getString<double>(chi[i]) + ", ";
            sql += (isnan(pval[i]) || isinf(pval[i])) ? "NULL," : getString<double>(pval[i]) + ", ";
            sql += (isnan((double)std::abs(log10(pval[i]))) || isinf((double)std::abs(log10(pval[i])))) ? "NULL," : getString<double>((double)std::abs(log10(pval[i]))) + ", ";
            sql += getString<int>(fams_used[i]) + ", ";
            sql += getString<int>(trans[i]) + ", ";
            sql += getString<int>(untrans[i]) + ", ";
            sql += "'" + data_set->get_locus(i)->getAllele1() + "', ";
            sql += "'" + data_set->get_locus(i)->getAllele2() + "', ";
            double OR = (double)trans[i] / (double)untrans[i];
            double OR_lower = exp(log(OR) - zt * sqrt(1/trans[i] + 1/untrans[i]));
            double OR_upper = exp(log(OR) + zt * sqrt(1/trans[i] + 1/untrans[i]));
            sql += (isnan(OR) || isinf(OR)) ? "NULL," : getString<double>(OR) + ", ";
            sql += (isnan(OR_lower) || isinf(OR_lower))? "NULL," : getString<double>(OR_lower) + ", ";
            sql += (isnan(OR_upper) || isinf(OR_upper)) ? "NULL)" : getString<double>(OR_upper) + ")";
            Controller::execute_sql(myQuery, sql);
        }
        data_set->get_locus(i)->setFlag(false);
    }

    myQuery.commit();

    cout << "END TDT INSERT" << endl;
    hasresults = true;
}//end method dump2db()

void ProcessRunTDT::run(DataSetObject* ds)
{
	process(ds);
}//end method run(DataSetObject* ds)

#endif

#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
