/*
 * DataSetObject.cpp
 *
 *  Created on: Jan 20, 2009
 *      Author: gilesjt
 */

#include <Helper.h>
#include "DataSetObject.h"


DataSetObject::~DataSetObject() {
	// TODO Auto-generated destructor stub
}

vector<string> DataSetObject::get_files(){
	vector<string> files;

	if(type == PEDMAP){
		files.push_back(pedfile);
		files.push_back(mapfile);
	}
	else if(type == TPEDFAM){
		files.push_back(tpedfile);
		files.push_back(tfamfile);
	}
	else if(type == BINARY){
		files.push_back(bin_bim);
		files.push_back(bin_bed);
		files.push_back(bin_fam);
	}
	return files;
}

string DataSetObject::summarize(){
	string summary = "";

	if(type == PEDMAP){
		summary += "PED File: ";
		summary += "\t" + pedfile;
		summary += "\n";
		summary += "MAP File: ";
		summary += "\t" + mapfile + "\n";
	}
	else if(type == TPEDFAM){
	}
	else if(type == BINARY){
	}

	summary += "Loci: " + getString<int>(this->num_loci()) + "\n";
	summary += "Samples: " + getString<int>(this->num_inds()) + "\n";
	summary += "Pedigrees: " + getString<int>(this->num_pedigrees()) + "\n";
	summary += "Affected: " + getString<int>(this->num_affected()) + "\n";
	summary += "Unaffected: " + getString<int>(this->num_unaffected()) + "\n";
	summary += "Males: " + getString<int>(this->num_males()) + "\n";
	summary += "Females: " + getString<int>(this->num_females()) + "\n";

	return summary;
}
