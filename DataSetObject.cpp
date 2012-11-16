/*
 * DataSetObject.cpp
 *
 *  Created on: Jun 25, 2010
 *      Author: cozartc
 */

//#include <QString>
//#include <QStringList>
#include <Sample.h>
#include <Helpers.h>
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

int DataSetObject::num_loci_enabled(){
    int count = 0;
    for(int i = 0; i < (int)this->num_loci(); i++){
        if(this->get_locus(i)->isEnabled()){
            count++;
        }
    }
    return count;
}

int DataSetObject::num_inds_enabled(){
    int count = 0;
    for(int i = 0; i < this->num_inds(); i++){
        if(this->get_sample(i)->isEnabled()){
            count++;
        }
    }
    return count;
}

//iterate over samples in the dataset and pull out the covariates that are in the mapping for those that are only existing in the dataset
//if it doesn't have a covariate, fill in with missing value.
void DataSetObject::load_covariates(vector<string> list, map<string, vector<double> > mapping){
//    mycovlist = list;
//    mycovmap = mapping;
//    this->set_covariates(&mycovlist);
    for(int i = 0; i < (int)list.size(); i++){
        this->add_covariate(list[i]);
    }

    for(int s = 0; s < this->num_inds(); s++){
        Methods::Sample* samp = this->get_sample(s);
        string key = samp->getFamID() + "#" + samp->getInd();
        map<string, vector<double> >::iterator iter = mapping.find(key);
        if(iter != mapping.end()){
            vector<double> values = iter->second;
            for(int l = 0; l < (int)values.size(); l++){
                samp->addCovariate(values[l]);
            }
        }
        else{
            for(int l = 0; l < (int)list.size(); l++){
                samp->addCovariate(-99999);
            }
        }
    }
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
            summary += "TPED File: ";
            summary += "\t" + tpedfile;
            summary += "\n";
            summary += "TFAM File: ";
            summary += "\t" + tfamfile + "\n";
	}
	else if(type == BINARY){
            summary += "BIM File: ";
            summary += "\t" + bin_bim;
            summary += "\n";
            summary += "BED File: ";
            summary += "\t" + bin_bed;
            summary += "\n";
            summary += "FAM File: ";
            summary += "\t" + bin_fam;
            summary += "\n";
        }

	summary += "Loci: " + getString<int>(this->num_loci()) + "\n";
        summary += "Loci enabled: " + getString<int>(this->num_loci_enabled()) + "\n";
	summary += "Samples: " + getString<int>(this->num_inds()) + "\n";
        summary += "Samples enabled: " + getString<int>(this->num_inds_enabled()) + "\n";
	summary += "Pedigrees: " + getString<int>(this->num_pedigrees()) + "\n";
//        summary += "Pedigrees enabled: " + getString<int>(this->num_pedigrees_enabled()) + "\n";
	summary += "Affected: " + getString<int>(this->num_affected()) + "\n";
//        summary += "Affected enabled: " + getString<int>(this->num_affected_enabled()) + "\n";
	summary += "Unaffected: " + getString<int>(this->num_unaffected()) + "\n";
//        summary += "Unaffected enabled: " + getString<int>(this->num_unaffected_enabled()) + "\n";
	summary += "Males: " + getString<int>(this->num_males()) + "\n";
	summary += "Females: " + getString<int>(this->num_females()) + "\n";
        summary += "Covariates: " + getString<int>(this->num_covariates()) + "\n";

	return summary;
}
