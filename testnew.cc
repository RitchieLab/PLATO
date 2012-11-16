#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include "Helper.h"
#include "Marker.h"
#include "Family.h"
#include "Sample.h"
using namespace std;

void readPed(string, vector<Sample*>*, vector<Family*>*, vector<Marker*>*);
void readMap(string, vector<Marker*>*);
bool readString(FILE*, string*);
void assignLinks(vector<Family*>*);

int main (int argc, char* argv[]){
	if(argc < 3){
		exit(1);
	}
	string ped = argv[1];
	string map = argv[2];

	vector<Marker*> markers;
	markers.resize(0);
	readMap(map, &markers);
	vector<Sample*> samples;
	vector<Family*> families;
	samples.resize(0);
	families.resize(0);

	readPed(ped, &samples, &families, &markers);
	assignLinks(&families);

	cout << "Sample size: " << samples.size() << endl;
	cout << "Marker size: " << markers.size() << endl;
	cout << "First sample one size: " << samples[1]->getAoneSize() << endl;
	cout << "Sample " << samples[1]->getFamID() << " " << samples[1]->getInd();
	if(samples[1]->getDad() != NULL){
		cout << " dad = " << samples[1]->getDad()->getInd();
	}
	if(samples[1]->getMom() != NULL){
		cout << " mom = " << samples[1]->getMom()->getInd();
	}
	cout << "Family size: " << families.size() << endl;
	cout << "First family " << families[1]->getFamID() << " total inds: " << families[1]->getTotalInds() << endl;
	
	vector<Family*>::iterator f_iter;
	for(f_iter = families.begin(); f_iter != families.end(); f_iter++){
		cout << "Family: " << (*f_iter)->getFamID() << "\t" << (*f_iter)->getSamples()->size() << endl;
		vector<Sample*>::iterator s_iter;
		for(s_iter = (*f_iter)->getSamples()->begin(); s_iter != (*f_iter)->getSamples()->end(); s_iter++){
			cout << "\tSample: " << (*s_iter)->getInd() << "\tMom: " << (*s_iter)->getMomID() << "\tDad: " << (*s_iter)->getDadID() << "\tChildren: " << (*s_iter)->getChildren()->size() << endl;
		}
	}
	//string temp;
	//cin >> temp;
}

void assignLinks(vector<Family*>* families){
	vector<Family*>::iterator f_iter;

	for(f_iter = families->begin(); f_iter != families->end(); f_iter++){
		vector<Sample*>::iterator s_iter;
		vector<Sample*>* samps = (*f_iter)->getSamples();
		vector<Sample*> child;
		for(s_iter = samps->begin(); s_iter != samps->end(); s_iter++){
			if((*s_iter)->getDadID() != "0" && (*s_iter)->getMomID() != "0"){
				vector<Sample*>::iterator temp_iter = find_if(samps->begin(), samps->end(), FindSampleByID((*s_iter)->getDadID()));
				if(temp_iter != samps->end()){
					(*s_iter)->setDad((*temp_iter));
				}
				temp_iter = find_if(samps->begin(), samps->end(), FindSampleByID((*s_iter)->getMomID()));
				if(temp_iter != samps->end()){
					(*s_iter)->setMom((*temp_iter));
				}
				child.push_back((*s_iter));
			}
		}
		for(s_iter = samps->begin(); s_iter != samps->end(); s_iter++){
			if((*s_iter)->getDadID() == "0" && (*s_iter)->getMomID() == "0"){
				vector<Sample*>::iterator temp_iter;
				for(temp_iter = child.begin(); temp_iter != child.end(); temp_iter++){
					vector<Sample*>::iterator has_child = find_if((*s_iter)->getChildren()->begin(), (*s_iter)->getChildren()->end(), FindSampleByID((*temp_iter)->getInd()));
					if(has_child == (*s_iter)->getChildren()->end()){
						(*s_iter)->addChild((*temp_iter));
					}
				}
			}
		}
		child.clear();
	}
}

void readMap(string file, vector<Marker*>* markers){
	ifstream input;
	input.open(file.c_str(), ios::in);

	if(!input){
		cerr << "Error opening map file: " << file << endl;
		exit(1);
	}

	while(!input.eof()){
		char buf[256];
		input.getline(buf, 256, '\n');

		string line = buf;

		if(line == ""){
			continue;
		}

		string temp;
		stringstream s(line);
		vector<string> elems;
		while(s >> temp){
			elems.push_back(temp);
		}

		if(elems.size() == 0){
			continue;
		}
		else if(elems.size() > 3){
			cerr << "Map file line has more than 3 elements: " + line << endl;
			exit(1);
		}

		string chr = elems[0];
		string probe_id = elems[1];
		int bploc = atoi(elems[2].c_str());

		Marker* m = new Marker(chr, probe_id, bploc);
		
		markers->push_back(m);
	}

	input.clear();
	input.close();
}

bool readString(FILE* fp, string* s){
	bool done = false;
	*s="";
	while (1)
	{
		char ch = fgetc(fp);
		if ( ch==' ' || ch == '\t' )
		{
			if (done)
				return true;
		}
		else if ( ch=='\n' || ch=='\r' || feof(fp) ) 
			return false;
		else
		{
			*s += ch;
			done = true;
		}
	}
}

void readPed(string file, vector<Sample*>* samples, vector<Family*>* families, vector<Marker*>* markers){
	FILE* input;
	input = fopen(file.c_str(), "r");

	while(!feof(input)){
		Sample* samp = new Sample();
		int f = 0;
		string temp = "";
		if(readString(input, &temp)){
		   samp->setFamID(temp);
	   	   f++;
		   temp = "";
		}
		
		if(samp->getFamID() == ""){
			continue;
		}

		/*check for comments?*/
		string sex = "";
		string pheno = "";
		if(readString(input, &temp)){
			samp->setInd(temp);
		   	f++;
			temp = "";
		}
		if(readString(input, &temp)){
			samp->setDadID(temp);
		   	f++;
			temp = "";
		}
		if(readString(input, &temp)){ 
			samp->setMomID(temp);
			f++;
			temp = "";
		}
		if(readString(input, &sex)) f++;
		if(readString(input, &pheno)) f++;

		if(sex == "1"){
			samp->setSex(true);
		}
		else if(sex == "2"){
			samp->setSex(false);
		}

		samp->resizeAlleles(markers->size());
		
		int gn = 0;
		int i = 0;
		bool linedone = false;
		bool fatal = false;

		string fmsg;
		while(!linedone){
			string one = "";
			string two = "";

			while(1){
				char ch = fgetc(input);

				if(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input)){
					if(ch == '\n' || ch == '\r' || feof(input)){
						linedone = true;
					}

					if(one.length() > 0){
						gn++;
						break;
					}
					if(ch == '\n' || ch == '\r' || feof(input)){
						break;
					}
				}
				else{
					one += ch;
				}
			}
			if(!linedone){
				while(1){
					char ch = fgetc(input);
					if(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input)){
						if(ch == '\n' || ch == '\r' || feof(input)){
							linedone = true;
						}
						if(two.length() > 0){
							gn++;
							break;
						}
						if(ch == '\n' || ch == '\r' || feof(input)){
							break;
						}
					}
					else{
						two += ch;
					}
				}
				if(linedone && one.length() == 0 && two.length() == 0){
					break;
				}

				Marker* m = (*markers)[i];
				if(one != "0"){
					if(one != m->getAllele1() && one != m->getAllele2()){
						if(m->getAllele1() == ""){
							m->setAllele1(one);
						}
						else if(m->getAllele2() == ""){
							m->setAllele2(one);
						}
					}
				}
				
				if(two != one){
					if(two != "0"){
						if(two != m->getAllele1() && two != m->getAllele2()){
							if(m->getAllele1() == ""){
							   	m->setAllele1(two);
							}
							else if(m->getAllele2() == ""){
								m->setAllele2(two);
							}
						}
					}		
				}

				if(one == m->getAllele1() && two == m->getAllele1()){
					samp->addAone(i, false);
					samp->addAtwo(i, false);
				}
				else if(one != "0" && two != "0" && one != two){
					samp->addAone(i, false);
					samp->addAtwo(i, true);
				}
				else if(one == m->getAllele2() && two == m->getAllele2()){
					samp->addAone(i, true);
					samp->addAtwo(i, true);
				}
				else if(one == "0" || two == "0"){
					samp->addAone(i, true);
					samp->addAtwo(i, false);
				}

				i++;
			}/*end !linedone*/
		}/*end while(1)*/

		vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

		if(f_iter != (*families).end()){
			(*f_iter)->AddInd(samp);
			samp->setFamily((*f_iter));
		}
		else{
			Family* fam = new Family();
/*			fam->Setcenter()*/
			fam->setFamID(samp->getFamID());
			fam->AddInd(samp);
			samp->setFamily(fam);
			families->push_back(fam);
		}
		
		samples->push_back(samp);
	}/*end while(eof)*/

	fclose(input);
}

