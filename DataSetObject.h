/*
 * DataSetObject.h
 *
 *  Created on: Jun 25, 2010
 *      Author: cozartc
 */

#ifndef DATASETOBJECT_H_
#define DATASETOBJECT_H_
#include <DataSet.h>
#include <string>
#include <StepOptions.h>
#include "Vars.h"

using namespace std;
#ifdef PLATOLIB
using namespace PlatoLib;
#endif

class DataSetObject : public Methods::DataSet{
public:
	DataSetObject(){DataSet(); pedfile = mapfile = tpedfile = tfamfile = bin_bim = bin_bed = bin_fam = ""; type = PEDMAP;}
	virtual ~DataSetObject();

	string get_pedfile(){return pedfile;}
	void set_pedfile(string s){pedfile = s; options.setPedFile(pedfile);}
	string get_mapfile(){return mapfile;}
	void set_mapfile(string s){mapfile = s; options.setMapFile(mapfile);}
	string get_tpedfile(){return tpedfile;}
	void set_tpedfile(string s){tpedfile = s; options.setTPedFile(tpedfile);}
	string get_tfamfile(){return tfamfile;}
	void set_tfamfile(string s){tfamfile = s; options.setTFamFile(tfamfile);}
	string get_bin_bim(){return bin_bim;}
	void set_bin_bim(string s){
		bin_bim = s;
		string temp = bin_bim;

	}
	string get_bin_bed(){return bin_bed;}
	void set_bin_bed(string s){bin_bed = s;}
	string get_bin_fam(){return bin_fam;}
	void set_bin_fam(string s){bin_fam = s;}
	void set_bin_prefix(string s){options.setBinInput(s);}
	void set_type(DataTypes dt){type = dt;}
	DataTypes get_type(){return type;}

	Methods::StepOptions get_options(){return options;}
	vector<string> get_files();

	string summarize();
	int num_loci_enabled();
	int num_inds_enabled();

	void load_covariates(vector<string>, map<string, vector<double> >);

protected:
	string pedfile;
	string mapfile;
	string tpedfile;
	string tfamfile;
	string bin_bim;
	string bin_bed;
	string bin_fam;

	vector<string> mycovlist;
	map<string, vector<double> > mycovmap;

	Methods::StepOptions options;
	DataTypes type;
};

#endif /* DATASETOBJECT_H_ */
