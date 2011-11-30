/**********************************************************************************
*                       Marker Genotype Efficiency Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
*Creates Plink based binary genotype files. (.bed, .fam, .bim)
*See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
*
*File: BINOutput.cc
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
#include "BINOutput.h"
#include "General.h"
#include "Helpers.h"

namespace Methods{

//DEPRECATED
void BINOutput::FilterSummary(){
}

//DEPRECATED
void BINOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}

}

//DEPRECATED
void BINOutput::filter(){
}

//DEPRECATED
Sample* BINOutput::find_sample(string i, string f){
	int ssize = samples->size();
	for(int s = 0; s < ssize; s++){
		if((*samples)[s]->getInd() == i && (*samples)[s]->getFamID() == f){
			return ((*samples)[s]);
		}
	}
	return NULL;
}

//DEPRECATED
bool BINOutput::find_marker(string p){
	int msize = markers->size();
	for(int i = 0; i < msize; i++){
		if((*markers)[i]->getProbeID() == p){
			return true;
		}
	}
	return false;
}

//DEPRECATED
int BINOutput::get_marker_loc(int i){
	int msize = markers->size();
	for(int m = 0; m < msize; m++){
		if((*markers)[m]->getLoc() == i){
			return (*markers)[m]->getLoc();
		}
	}
	return -1;
}

void BINOutput::process(vector<Sample*>* ss, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = ss;
	marker_map = mm;

	writeBit(ss, f, m, mm, &options, overwrite, order);

}

/*
 *Function: writeBit
 *Parameters: sample vector, family vector, marker vector, marker map vector
 *Description:
 *Writes genotype (bed) file as binary genotype file based on Plink output.
 *Family output file (fam) is first 6 columns of original ped file
 *Map file (bim) is extended map file including alleles.
 *
 *
 */
void BINOutput::writeBit(vector<Sample*>* samples, vector<Family*>* families, vector<Marker*>* markers, vector<int>* marker_map, StepOptions* options, bool overwrite, int order){
	string fname = opts::_OUTPREFIX_ + "binary" + options->getOut();
	if(!overwrite){
		fname += "." + getString<int>(order);
	}
	fname += ".fam";

	//output .fam file first
	opts::printLog("Writing family information to " + fname + "\n");
	ofstream BIT(fname.c_str(), ios::out);
	if(!BIT){
		opts::printLog("Error opening " + fname + " for family file creation.  Exiting!\n");
		throw MethodException("Error opening " + fname + " for family file creation.  Exiting!\n");
	}
	filenames.push_back(fname);

	int ssize = samples->size();
	map<string, string> dummy;

	string male = "\t0\t0\t1\t0";
	string female = "\t0\t0\t2\t0";

	for(int s = 0; s < ssize; s++){
		Sample* samp = (*samples)[s];
		if(!samp->isEnabled() && !samp->isExcluded() && !options->doIncDisabledSamples()){
			continue;
		}
		int valid_parents = 0;
		if(samp->getDad() != NULL){
			valid_parents++;
		}
		if(samp->getMom() != NULL){
			valid_parents++;
		}

		BIT << samp->getFamID() << " "
			<< samp->getInd() << " ";

		if((samp->getDad() == NULL || !samp->getDad()->isEnabled()) && samp->getDadID() != "0" && options->doRemMissingParents()){
			BIT << "0 ";
		}
		else if(samp->getDad() == NULL && samp->getDadID() != "0" && options->doDummyMissingParents()){
			BIT << samp->getDadID() << " ";
			dummy[samp->getFamID() + " " + samp->getDadID()] =  male;
		}
		else if(samp->getDad() == NULL && samp->getMom() != NULL && options->doDummyIncompleteParentIds()){
			string did = "DUM999999";
			if(samp->getDadID() != "0"){
				did = samp->getDadID();
			}
			BIT << did << " ";
			dummy[samp->getFamID() + " " + did] = male;
		}
		else if(options->doZeroIncompleteTrioIds() && valid_parents != 2){
			BIT << "0 ";
		}
		else{
			BIT << samp->getDadID() << " ";
		}
		if((samp->getMom() == NULL || !samp->getMom()->isEnabled()) && samp->getMomID() != "0" && options->doRemMissingParents()){
			BIT << "0 ";
		}
		else if(samp->getMom() == NULL && samp->getMomID() != "0" && options->doDummyMissingParents()){
			BIT << samp->getMomID() << " ";
			dummy[samp->getFamID() + " " + samp->getMomID()] = female;
		}
		else if(samp->getMom() == NULL && samp->getDad() != NULL && options->doDummyIncompleteParentIds()){
			string mid = "DUM999998";
			if(samp->getMomID() != "0"){
				mid = samp->getMomID();
			}
			BIT << mid << " ";
			dummy[samp->getFamID() + " " + mid] = female;
		}
		else if(options->doZeroIncompleteTrioIds() && valid_parents != 2){
			BIT << "0 ";
		}
		else{
		    BIT << samp->getMomID() << " ";
		}
		if(samp->getSex()){
			BIT << "1 ";
		}
		else{
			BIT << "2 ";
		}
		if(options->getUsePheno()){
			BIT << samp->getPheno(options->getPhenoLoc());
		}
		else{
			BIT << samp->getPheno();
		}



		BIT << endl;

	}

	if(options->doDummyMissingParents() || options->doDummyIncompleteParentIds()){

		map<string, string>::iterator iter;
		for(iter = dummy.begin(); iter != dummy.end(); iter++){
			BIT << iter->first << iter->second << endl;
		}
	}

	BIT.clear();
	BIT.close();

	//now output snp file
	fname = opts::_OUTPREFIX_ + "binary" + options->getOut();
	if(!overwrite){
		fname += "." + getString<int>(order);
	}
	fname += ".bim";

	opts::printLog("Writing map information to " + fname + "\n");
	BIT.open(fname.c_str(), ios::out);
	if(!BIT){
		opts::printLog("Error opening " + fname + " for map file creation.  Exiting!\n");
		throw MethodException("Error opening " + fname + " for map file creation.  Exiting!\n");
	}
	filenames.push_back(fname);

	int msize = markers->size();
	int prev_base = 0;
	int prev_chrom = -1;
	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, options);
	int gsize = good_markers.size();

	for(int m = 0; m < gsize; m++){
		Marker* mark = good_markers[m];
		if(mark->isEnabled()){

			BIT << mark->getChrom() << "\t";
			BIT << mark->getProbeID() << "\t"
				<< "0\t" //centimorgan column
				<< mark->getBPLOC() << "\t";
			if(mark->getAllele1().length() == 0){
				BIT << "0\t";
			}
			else{
				BIT << mark->getAllele1() << "\t";
			}
			if(mark->getAllele2().length() == 0){
				BIT << "0\t";
			}
			else{
				BIT << mark->getAllele2() << "\t";
			}
			BIT << endl;
		}
	}
	BIT.clear();
	BIT.close();

	//here's the fun part...the binary file.  Please see the url mentioned at the top for references.
	fname = opts::_OUTPREFIX_ + "binary" + options->getOut();
	if(!overwrite){
		fname += "." + getString<int>(order);
	}
	fname += ".bed";
	opts::printLog("Writing genotype bitfile to " + fname + "\n");

	BIT.open(fname.c_str(), ios::out | ios::binary);
	if(!BIT){
		opts::printLog("Error opening " + fname + " for genotype file creation.  Exiting!\n");
		throw MethodException("Error opening " + fname + " for genotype file creation.  Exiting!\n");
	}
	filenames.push_back(fname);

	bitset<8> b;
	char ch[1];

  	b.reset();
    b.set(2);  b.set(3);  b.set(5);  b.set(6);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);

	b.reset();
	b.set(0);  b.set(1);  b.set(3);  b.set(4);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);

	b.reset();
    ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);

	for(int s = 0; s < ssize; s++){
		Sample* samp = (*samples)[s];
		if(!samp->isEnabled() && !samp->isExcluded() && !options->doIncDisabledSamples()){
			continue;
		}
		prev_base = 0;
		prev_chrom = -1;
		for(int m = 0; m < gsize;){
			bitset<8> b;
			b.reset();
			int c = 0;
			while(c < 8 && m < gsize){
				Marker* mark = good_markers[m];
				if(mark->isEnabled() && !mark->isMicroSat()){
					int loc = mark->getLoc();

					if(samp->getAone(loc) && samp->getAtwo(loc) && samp->getAmissing(loc)){
						b.set(c);
						c++;
						c++;
					}else{

						if(samp->getAone(loc))
							b.set(c);

						c++;
						if(samp->getAtwo(loc))
							b.set(c);
						c++;
					}
				}
				m++;
			}

			char ch[1];
			ch[0] = (char)b.to_ulong();
			BIT.write(ch, 1);
		}
	}
	if(options->doDummyMissingParents() || options->doDummyIncompleteParentIds()){

		map<string, string>::iterator iter;
		for(iter = dummy.begin(); iter != dummy.end(); iter++){
			prev_base = 0;
			prev_chrom = -1;
			for(int m = 0; m < msize;){
				bitset<8> b;
				b.reset();
				int c = 0;
				while(c < 8 && m < gsize){
					Marker* mark = good_markers[m];
					if(mark->isEnabled() && !mark->isMicroSat()){

						b.set(c);

						c++;
						c++;
					}
					m++;
				}

				char ch[1];
				ch[0] = (char)b.to_ulong();
				BIT.write(ch, 1);
			}
		}
	}

	BIT.close();
}


int BINOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}


}
