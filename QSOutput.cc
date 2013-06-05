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
*File: QSOutput.cc
**********************************************************************************/

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include "QSOutput.h"
#include "Chrom.h"
#include "General.h"


void QSOutput::FilterSummary(){
}

void QSOutput::PrintSummary(){

}

void QSOutput::filter(){
}

Sample* QSOutput::find_sample(string i, string f){
	int ssize = samples->size();
	for(int s = 0; s < ssize; s++){
		if((*samples)[s]->getInd() == i && (*samples)[s]->getFamID() == f){
			return ((*samples)[s]);
		}
	}
	return NULL;
}

bool QSOutput::find_marker(int p){
	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		if((*markers)[i]->getSysprobe() == p){
			return true;
		}
	}
	return false;
}

int QSOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

string QSOutput::map_allele(string a){
	if(a == "A"){
		return "1";
	}
	else if(a == "C"){
		return "2";
	}
	else if(a == "G"){
		return "3";
	}
	else if(a == "T"){
		return "4";
	}
	else{
		return "0";
	}
}

