#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include "Globals.h"
#include "Options.h"
#include "Finalize.h"
#include "Helper.h"

using namespace std;

void Finalize::finish(vector<Marker*>* markers, vector<Sample*>* samples, vector<Family*>* families){
	string mname = opts::_OUTPREFIX_ + "Remaining_makers.txt";
	string sname = opts::_OUTPREFIX_ + "Remaining_samples.txt";
	string fname = opts::_OUTPREFIX_ + "Remaining_families.txt";

	ofstream mout (mname.c_str());
	ofstream sout (sname.c_str());
	ofstream fout (fname.c_str());

	if(!mout){
		opts::printLog("Error opening: " + mname + ".  Exiting!");
		exit(1);
	}
	if(!sout){
		opts::printLog("Error opening: " + sname + ".  Exiting!");
		exit(1);
	}
	if(!fout){
		opts::printLog("Error opening: " + fname + ".  Exiting!");
		exit(1);
	}

	for(int m = 0; m < (int)markers->size(); m++){
		if((*markers)[m]->isEnabled()){
			mout << (*markers)[m]->toString() << endl;
		}
	}

	for(int s = 0; s < (int)samples->size(); s++){
		if((*samples)[s]->isEnabled()){
			sout << (*samples)[s]->toString() << endl;
		}
	}

	for(int f = 0; f < (int)families->size(); f++){
		if((*families)[f]->isEnabled()){
			fout << (*families)[f]->toString() << endl;
		}
	}

	if(mout.is_open()){
		mout.close();
	}
	if(sout.is_open()){
		sout.close();
	}
	if(fout.is_open()){
		fout.close();
	}
}
