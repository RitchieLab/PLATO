/**********************************************************************************
*                       Run of Homozygosity Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Performs a run of homozygosity search based on Lencz et al (PNAS 0710021104)
* Performs a ROH based on straight forward common homozygous spans.
*
*
*File: Homozygous.cc
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
#include "ProcessHomozygous.h"
#include <Options.h>
#include <General.h>
#include <Helpers.h>
#include <cdflib.h>
#include "Controller.h"

using namespace Methods;

#ifdef PLATOLIB
namespace PlatoLib
{
#endif

string ProcessHomozygous::stepname = "homozygous";

ProcessHomozygous::ProcessHomozygous(string bn, int pos, Database* pdb, string projPath)
{
	name = "Homozygous";
	batchname = bn;
	position = pos;
	hasresults = false;
	db = pdb;
	projectPath = projPath;
}

void ProcessHomozygous::FilterSummary(){

	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;
}

void ProcessHomozygous::PrintSummary(){
	if(options.doHomozygPermute()){
		return;
	}
	if(!options.doHomozygWGHA()){
		string filename = opts::_OUTPREFIX_ + "homozygous_marker" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if(!overwrite){
			filename += "." + getString<int>(order);
		}
		ofstream homo (filename.c_str());
		if(!homo.is_open()){
			opts::printLog("Unable to open " + filename + "\n");
			throw MethodException("");
		}
		opts::addFile("Marker", stepname, filename);
		homo << "Chrom\trsID\tProbeID\tbploc\tUnAff\tAff\tTotal_Maj\tTotal_Min\tTotal\tTotal%\n";
		int msize = good_markers.size();//data_set->num_loci();
		opts::addHeader(filename, "UnAff");
		opts::addHeader(filename, "Aff");
		opts::addHeader(filename, "Total_Maj");
		opts::addHeader(filename, "Total_Min");
		opts::addHeader(filename, "Total");
		opts::addHeader(filename, "Total%");

		for(int i = 0; i < msize; i++){
			Marker* m = good_markers[i];//data_set->get_locus(i);
			float per = (((float)homoallcount[i] / (float)opts::_SAMPLES_WORKING_));// * 100.0f);
			homo.precision(4);
			homo << m->getChrom() << "\t" << m->getRSID() << "\t" << m->getProbeID() << "\t"
				<< m->getBPLOC() << "\t" << homounaffcount[i] << "\t" << homoaffcount[i] << "\t"
				<< homomajallcount[i] << "\t" << homominallcount[i] << "\t"
				<< (homounaffcount[i] + homoaffcount[i]) << "\t" << per << endl;
		}
		if(homo.is_open()){
			homo.close();
		}
	}

}

void ProcessHomozygous::filter(){}

void ProcessHomozygous::process(DataSet* ds){
	data_set = ds;
	good_markers = Helpers::findValidMarkers(data_set->get_markers(), &options);

	Homozygous hom(data_set);
	hom.setOptions(options);
	hom.setOverwrite(this->overwrite);
	hom.setOrder(this->order);
	hom.calculate();

	if(!options.doHomozygPermute() && !options.doHomozygWGHA()){
		homounaffcount = hom.getHomoUnaffCount();
		homoaffcount = hom.getHomoAffCount();
		homomajallcount = hom.getHomoMajAllCount();
		homominallcount = hom.getHomoMinAllCount();
		homoallcount = hom.getHomoAllCount();
	}
#ifdef PLATOLIB
	for(int i = 0; i < (int)(hom.get_filenames().size()); i++)
	{
		filenames.push_back(hom.get_filenames().at(i));
	}
#endif

}
#ifdef PLATOLIB
void ProcessHomozygous::create_tables()
{
	Query myQuery(*db);
	myQuery.transaction();
    for(unsigned int i = 0; i < tablename.size(); i++)
    {
        Controller::drop_table(db, tablename[i]);
    }
    headers.clear();
    tablename.clear();
    primary_table.clear();

    string tempbatch = batchname;
    for(unsigned int i = 0; i < tempbatch.size(); i++)
    {
        if(tempbatch[i] == ' '){
            tempbatch[i] = '_';
        }
    }
    string mytablename = tempbatch + "_";
    tempbatch = name;
    for(unsigned int i = 0; i < tempbatch.size(); i++)
    {
        if(tempbatch[i] == ' '){
            tempbatch[i] = '_';
        }
    }
    string base = mytablename + tempbatch;

    mytablename = base + "_marker_" + getString<int>(position);
    tablename.push_back(mytablename);
    tablenicknames.push_back("");
    primary_table[mytablename].push_back(Vars::LOCUS_TABLE);


    string sql = "CREATE TABLE " + mytablename + " (id integer primary key,";
    sql += "fkey integer not null,";
    sql += "UnAff integer,";
    headers[mytablename].push_back("UnAff");
    sql += "Aff integer,";
    headers[mytablename].push_back("Aff");
    sql += "Total_Maj integer,";
    headers[mytablename].push_back("Total_Maj");
    sql += "Total_Min integer,";
    headers[mytablename].push_back("Total_Min");
    sql += "Total integer,";
    headers[mytablename].push_back("Total");
    sql += "TotalPercent real";
    headers[mytablename].push_back("TotalPercent");
    sql += ")";
    defaultinsert = "INSERT INTO " + mytablename + " (id, fkey, UnAff, Aff, Total_Maj, Total_Min, Total, TotalPercent) VALUES (NULL";

    Controller::execute_sql(myQuery, sql);
    myQuery.commit();
}//end method create_tables

void ProcessHomozygous::dump2db()
{
    if(options.doHomozygPermute())
    {
        return;
    }
    if(!options.doHomozygWGHA())
    {
        create_tables();

        Query myQuery(*db);
        myQuery.transaction();

        int msize = data_set->num_loci();

        for(int i = 0; i < msize; i++)
        {
            //Marker* m = data_set->get_locus(i);
            float per = (((float)homoallcount[i] / (float)opts::_SAMPLES_WORKING_));// * 100.0f);

            string sql = defaultinsert;
            sql += "," + getString<int>(data_set->get_locus(i)->getSysprobe());
            sql += "," + getString<int>(homounaffcount[i]);
            sql += "," + getString<int>(homoaffcount[i]);
            sql += "," + getString<int>(homomajallcount[i]);
            sql += "," + getString<int>(homominallcount[i]);
            sql += "," + getString<int>(homounaffcount[i] + homoaffcount[i]);
            sql += "," + ((isnan(per) || isinf(per)) ? "NULL" : getString<float>(per));
            sql += ")";
            Controller::execute_sql(myQuery, sql);
            data_set->get_locus(i)->setFlag(false);
        }
        myQuery.commit();
    }
}//end method dump2db

void ProcessHomozygous::resize(int i){}

void ProcessHomozygous::run(DataSetObject* ds)
{
	FixOutputName();
	process(ds);
}

void ProcessHomozygous::FixOutputName()
{
#ifdef WIN
	options.setOverrideOut(projectPath);
#else
	options.setOverrideOut(projectPath);
#endif
}

#endif
#ifdef PLATOLIB
}//end namespace PlatoLib
#endif
