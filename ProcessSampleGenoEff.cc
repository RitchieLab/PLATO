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
#include <vector>
#include <list>
#include <map>
#include "ProcessSampleGenoEff.h"
#include <Sample.h>
#include <Family.h>
//#include "Chrom.h"
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

string ProcessSampleGenoEff::stepname = ProcessSampleGenoEff::doRegister("sample-geno-eff");

#ifdef PLATOLIB
ProcessSampleGenoEff::ProcessSampleGenoEff(string bn, int pos, Database* pdb)
{
	// TODO Auto-generated constructor stub
	name = "Sample Efficiency";
        batchname = bn;
        position = pos;
        hasresults = false;
        db = pdb;
}
#endif

void ProcessSampleGenoEff::process(DataSet* ds){
	data_set = ds;

	SampleGenoEff sge(data_set);
	sge.setOptions(options);

	int ssize = data_set->num_inds();
	zeros.resize(ssize, 0);
	total.resize(ssize, 0);
	if(opts::_ENZYMES_){
//		enzyme_zeros.resize(ssize, 0);
//		enzyme_total.resize(ssize, 0);
	}
	orig_num_samples = 0;
	//for(int j = 0; j < ssize; j++){
	//	if((*samples)[j]->isEnabled()){
	//		orig_num_samples++;
	//	}
	//}

	for(int i = 0; i < ssize; i++){
		if(data_set->get_sample(i)->isEnabled()){
			sge.calculate(i);
			zeros[i] += sge.getZeros();
			total[i] += sge.getTotal();
		}
	}
}//end method process(DataSet* ds)

void ProcessSampleGenoEff::PrintSummary(){

	string fname1 = opts::_OUTPREFIX_ + "sample_geno_eff" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	ofstream indeff (fname1.c_str());
	if(!indeff.is_open()){
		opts::printLog("Unable to open " + fname1 + "\n");
		throw MethodException("");
	}
	opts::addFile("Sample",stepname, fname1);
	string sdetails = "";
	if(opts::_SAMPDESC_.length() > 0){
		sdetails = data_set->get_sample(0)->getDetailHeaders();
	}
	opts::addHeader(fname1, "%GenoEff_All");
	indeff << "FamID\t"
		   << "IndID\tCenter\tSex\t"
		   << "Affection_Status\t"
		   << "Plate\tWell\t%GenoEff_All";
/*	vector<string> enzymes;
	if(opts::_ENZYMES_){
		int prev_base = 0;
		int prev_chrom = -1;
		int msize = markers->size();
		for(int i = 0; i < msize; i++){
			if(data_set->get_locus(i)->isEnabled() && !data_set->get_locus(i)->isFlagged() && isValidMarker(data_set->get_locus(i), &options, &prev_base, &prev_chrom)){
				vector<string>::iterator e_iter = find(enzymes.begin(), enzymes.end(), data_set->get_locus(i)->getEnzyme());
				if(e_iter == enzymes.end()){
					enzymes.push_back(data_set->get_locus(i)->getEnzyme());
				}
			}
		}
		if(enzymes.size() > 1){
			sort(enzymes.begin(), enzymes.end());
		}
		for(int i = 0; i < enzymes.size(); i++){
			indeff << "\tGenoEff_" << enzymes[i];
			opts::addHeader(fname1, "GenoEff_" + enzymes[i]);
		}
	}
*/
	if(opts::_SAMPDESC_.length() > 0){
		indeff << "\t" << sdetails;
	}
	indeff << endl;

	indeff.precision(4);

	int ssize = data_set->num_inds();

	for(int i = 0; i < ssize; i++){
		if(data_set->get_sample(i)->isEnabled()){
			float percent = 0.0f;

			if(total[i] > 0){
				percent = (1.0f - ((float)zeros[i]/(float)total[i]));// * 100.0f;
			}

			indeff << data_set->get_sample(i)->getFamID() << "\t"
				   << data_set->get_sample(i)->getInd() << "\t"
				   << data_set->get_sample(i)->getFamily()->getCenter() << "\t";
			if(data_set->get_sample(i)->getSex()){
				indeff << "M\t";
			}
			else{
				indeff << "F\t";
			}
			indeff << data_set->get_sample(i)->getPheno() << "\t";
			//if((*samples)[i]->getAffected()){
			//	indeff << "Y\t";
			//}
			//else{
			//	indeff << "N\t";
			//}
			indeff << data_set->get_sample(i)->getPlate() << "\t"
				   << data_set->get_sample(i)->getWell() << "\t"
				   << percent;
/*			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < enzymes.size(); e++){
					float epercent = 0.0f;
					if(enzyme_total[i][enzymes[e]] > 0){
						epercent = (1.0f - ((float)enzyme_zeros[i][enzymes[e]] / (float)enzyme_total[i][enzymes[e]])) * 100.0f;
					}
					indeff << "\t" << epercent;
				}
			}
*/
			if(opts::_SAMPDESC_.length() > 0){
				indeff << "\t" << data_set->get_sample(i)->getDetails();
			}
			indeff << endl;
		}
	}

	int msize = data_set->num_loci();
	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}

	if(indeff.is_open()){
		indeff.close();
	}

}

void ProcessSampleGenoEff::filter(){
	if(options.doThreshSamplesLow() || options.doThreshSamplesHigh()){
		int ssize = data_set->num_inds();

		for(int i = 0; i < ssize; i++){
			if(data_set->get_sample(i)->isEnabled()){
				float percent = 0.0f;
				bool inc = false;
				if(total[i] > 0){
					percent = (1.0f - ((float)zeros[i]/(float)total[i]));// * 100.0f;
				}

				if(options.doThreshSamplesLow() && Helpers::dLess(percent, options.getThreshSamplesLow())){
					data_set->get_sample(i)->setEnabled(false);
					inc = true;
				}
				if(options.doThreshSamplesHigh() && Helpers::dGreater(percent, options.getThreshSamplesHigh())){
					data_set->get_sample(i)->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_samples++;
				}
			}
		}
	}
}//end method filter()

void ProcessSampleGenoEff::FilterSummary(){
////	int cutsamps = 0;

	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Samples Passed:\t" + getString<int>(opts::_SAMPLES_WORKING_ - orig_num_samples) + " (" +
		getString<float>(((float) (opts::_SAMPLES_WORKING_ - orig_num_samples) / (float) opts::_SAMPLES_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_SAMPLES_WORKING_) + "\n");
	opts::_SAMPLES_WORKING_ -= orig_num_samples;

}//end method FilterSummary()

#ifdef PLATOLIB
void ProcessSampleGenoEff::create_tables()
{
	Query myQuery(*db);
	myQuery.transaction();

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
    primary_table[mytablename].push_back(Vars::SAMPLE_TABLE);

    string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
    sql += "fkey integer not null,";
    sql += "percent_GenoEff_All REAL";
    headers[mytablename].push_back("percent_GenoEff_All");
    sql += ")";
    Controller::execute_sql(myQuery, sql);
}//end method create_tables()

void ProcessSampleGenoEff::dump2db()
{
    create_tables();

    Query myQuery(*db);
    myQuery.transaction();

    for(int i = 0; i < data_set->num_inds(); i++){
        if(data_set->get_sample(i)->isEnabled()){
            float percent = 0.0f;
            if(total[i] > 0){
                percent = (1.0f - ((float)zeros[i]/(float)total[i]));
            }
            string sql = "INSERT INTO " + tablename.at(0) + " (id, fkey, percent_GenoEff_All) VALUES (NULL";
            sql += "," + getString<int>(data_set->get_sample(i)->getSysid());
            sql += "," + getString<float>(percent);
            sql += ")";
            Controller::execute_sql(myQuery, sql);
        }
    }
    myQuery.commit();

    hasresults = true;
    for(int i = 0; i < (int)data_set->num_loci(); i++){
        data_set->get_locus(i)->setFlag(false);
    }
}//end method dump2db()

void ProcessSampleGenoEff::run(DataSetObject* ds)
{

	process(ds);
}
#endif

#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
