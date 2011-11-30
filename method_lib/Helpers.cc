#include <iostream>
#include <string>
#include <vector>
#include "Helpers.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <string>
#include <sstream>
#include <list>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/algorithm/string.hpp>
#include <map>
#include "cdflib.h"
#include "Family.h"
#include "Marker.h"
#include "Sample.h"
#include "StepOptions.h"
#include "DataSet.h"
#include "InputFilter.h"
#include "MethodException.h"
using boost::math::complement;


namespace Methods{

const double Helpers::a[] =
  {
	      -3.969683028665376e+01,
		      2.209460984245205e+02,
			      -2.759285104469687e+02,
				      1.383577518672690e+02,
					      -3.066479806614716e+01,
						       2.506628277459239e+00
								     };

const double Helpers::b[] =
  {
	      -5.447609879822406e+01,
		      1.615858368580409e+02,
			      -1.556989798598866e+02,
				      6.680131188771972e+01,
					      -1.328068155288572e+01
							    };

const double Helpers::c[] =
  {
	      -7.784894002430293e-03,
		      -3.223964580411365e-01,
			      -2.400758277161838e+00,
				      -2.549732539343734e+00,
					      4.374664141464968e+00,
						       2.938163982698783e+00
								     };
const double Helpers::d[] =
  {
	      7.784695709041462e-03,
		      3.224671290700398e-01,
			      2.445134137142996e+00,
				      3.754408661907416e+00
						    };


double Helpers::DoubDif(double a, double b)
{
	double c = Abs(a);
	double d = Abs(b);

	d = Max(c, d);

	return d == 0.0 ? 0.0 : Abs(a - b) / d;
}
float Helpers::FloatDif(float a, float b)
{
	float c = Abs(a);
	float d = Abs(b);

	d = Max(c, d);

	return d == 0.0 ? 0.0 : Abs(a - b) / d;
}
bool Helpers::float_comp(float a, float b){
   if(FloatDif(a, b) <= FTOLERANCE){
   return true;
   }
	return false;
}
bool Helpers::double_comp(double a, double b){
   if(DoubDif(a, b) <= TOLERANCE){
   return true;
   }
	return false;
}
bool Helpers::fEquals(float a, float b){
	if(float_comp(a,b)){
		return true;
	}
	return false;
}
bool Helpers::dEquals(double a, double b){
	if(double_comp(a,b)){
		return true;
	}
	return false;
}
bool Helpers::fLess(float a, float b){
	if(float_comp(a, b)){
		return false;
	}
	if(a < b){
		return true;
	}
	return false;
}
bool Helpers::dLess(double a, double b){
	if(double_comp(a, b)){
		return false;
	}
	if(a < b){
		return true;
	}
	return false;
}
bool Helpers::fLessOrEqual(float a, float b){
	if(float_comp(a, b)){
		return true;
	}
	if(a < b){
		return true;
	}
	return false;
}
bool Helpers::dLessOrEqual(double a, double b){
	if(double_comp(a, b)){
		return true;
	}
	if(a < b){
		return true;
	}
	return false;
}
bool Helpers::fGreater(float a, float b){
	if(float_comp(a, b)){
		return false;
	}
	if(a > b){
		return true;
	}
	return false;
}
bool Helpers::dGreater(double a, double b){
	if(double_comp(a, b)){
		return false;
	}
	if(a > b){
		return true;
	}
	return false;
}
bool Helpers::fGreaterOrEqual(float a, float b){
	if(float_comp(a, b)){
		return true;
	}
	if(a > b){
		return true;
	}
	return false;
}
bool Helpers::dGreaterOrEqual(double a, double b){
	if(double_comp(a, b)){
		return true;
	}
	if(a > b){
		return true;
	}
	return false;
}

bool Helpers::fileExists(const string& fileName){
	fstream fin;
	fin.open(fileName.c_str(), ios::in);
	if(fin.is_open()){
		fin.close();
		return true;
	}
	fin.close();
	return false;
}

bool Helpers::isAlphaNum( string s ){
	for(unsigned int i = 0; i < s.size(); i++){
		if(!isdigit(s[i])){
			return true;
		}
	}
	return false;
}

double Helpers::p_from_chi(double chi, double df){
	if(chi <= 0){
		return -1;
	}
	return cdf(complement(boost::math::chi_squared(df), chi));
}

void Helpers::printFamsToDigit(vector<Methods::Family*>* families, string name, StepOptions options){
	string fname = opts::_OUTPREFIX_ + name + "_family_remap" + options.getOut() + ".txt";
	opts::printLog("Family Map file being written. [" + fname + "]\n");
	ofstream out (fname.c_str());
	if(!out.is_open()){
		opts::printLog("Unable to open " + fname + " for output!\n");
		//exit(1);
		throw MethodException("Unable to open " + fname + " for output!\n");
	}
	out << "FamID_Digit\tSampID_Digit\tDadID_Digit\tMomID_Digit\tFamID_Orig\tSampID_Orig\tDadID_Orig\tMomID_Orig\n";
	for(unsigned int i = 0; i < families->size(); i++){
		Methods::Family* fam = (*families)[i];
		if(fam->isEnabled()){
			vector<Methods::Sample*>* samps = fam->getSamples();
			for(unsigned int s = 0; s < samps->size(); s++){
				Methods::Sample* samp = (*samps)[s];
				if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
					out << fam->getFamID_digit() << "\t" << samp->getInd_digit() << "\t" << samp->getDadID_digit() << "\t" << samp->getMomID_digit() << "\t" << samp->getFamIDOrig() << "\t" << samp->getIndOrig() << "\t" << samp->getDadIDOrig() << "\t" << samp->getMomIDOrig() << endl;
				}
			}
		}
	}

	out.close();
}

void Helpers::remapFamsToDigit(vector<Methods::Family*>* families){
	for(unsigned int i = 0; i < families->size(); i++){
		Methods::Family* fam = (*families)[i];
		fam->setFamDigit(i + 1);
		bool badsamples = true;
		//vector<Methods::Sample*>* samps = fam->getSamples();

		if(badsamples){
			vector<Methods::Sample*>* founders = fam->getFounders();
			vector<Methods::Sample*>* nonfounders = fam->getNonFounders();
			for(unsigned int j = 0; j < founders->size(); j++){
				Methods::Sample* samp = (*founders)[j];
				samp->setIndDigit(j + 1);
			}
			int nonfstart = founders->size() + 1;
			for(unsigned int j = 0; j < nonfounders->size(); j++){
				Methods::Sample* samp = (*nonfounders)[j];
				samp->setIndDigit(j + nonfstart);
			}
			for(unsigned int j = 0; j < nonfounders->size(); j++){
				Methods::Sample* samp = (*nonfounders)[j];
				samp->setParentDigits();
			}
		}
	}
}

bool Helpers::isValidMarker(Methods::Marker* mark, StepOptions* options, int& prev_base, int& prev_chrom){
	if(mark->isEnabled()){
		if(options->doChrom()){
			if(!options->checkChrom(mark->getChrom())){
				return false;
			}
			if(!options->checkBp(mark->getBPLOC())){
				return false;
			}
		}
		if(options->doBpSpace()){
			if(prev_base == 0){
				prev_base = mark->getBPLOC();
				prev_chrom = mark->getChrom();
			}
			else{
				if(mark->getChrom() == prev_chrom && ((mark->getBPLOC() - prev_base) < options->getBpSpace())){
					return false;
				}
				prev_base = mark->getBPLOC();
				prev_chrom = mark->getChrom();
			}
		}
		return true;
	}
	return false;
}

vector<Methods::Marker*> Helpers::findValidMarkers(vector<Methods::Marker*>* marks, StepOptions options){
	int prev_base = 0;
	int prev_chrom = -1;
	vector<Methods::Marker*> good_markers;
	int msize = marks->size();
	for(int m = 0; m < msize; m++){
		Methods::Marker* mark = (*marks)[m];
		if(mark->isEnabled()){
			if(options.getAutosomeOnly() && ((mark->getChrom() >= opts::_CHRX_) || (mark->getChrom() < 1))){
				continue;

			}
			if(options.doChrom()){
				if(!options.checkChrom(mark->getChrom())){
					continue;
				}
				if(!options.checkBp(mark->getBPLOC())){
					continue;
				}
			}
			if(options.doBpSpace()){
				if(prev_base == 0){
					prev_base = mark->getBPLOC();
					prev_chrom = mark->getChrom();
				}
				else{
					if(mark->getChrom() == prev_chrom && ((mark->getBPLOC() - prev_base) < options.getBpSpace())){
						continue;
					}
					prev_base = mark->getBPLOC();
					prev_chrom = mark->getChrom();
				}
			}
			good_markers.push_back(mark);
		}
	}

	return good_markers;
}

vector<Methods::Marker*> Helpers::findValidMarkers(vector<Methods::Marker*>* marks, StepOptions* options){
	int prev_base = 0;
	int prev_chrom = -1;
	vector<Methods::Marker*> good_markers;
	int msize = marks->size();
	for(int m = 0; m < msize; m++){
		Methods::Marker* mark = (*marks)[m];
		if(mark->isEnabled()){
			if(options->getAutosomeOnly() && ((mark->getChrom() >= opts::_CHRX_) || (mark->getChrom() < 1))){
				continue;
			}

			if(options->doChrom()){
				if(!options->checkChrom(mark->getChrom())){
					continue;
				}
				if(!options->checkBp(mark->getBPLOC())){
					continue;
				}
			}
			if(options->doBpSpace()){
				if(prev_base == 0){
					prev_base = mark->getBPLOC();
					prev_chrom = mark->getChrom();
				}
				else{
					if(mark->getChrom() == prev_chrom && ((mark->getBPLOC() - prev_base) < options->getBpSpace())){
						continue;
					}
					prev_base = mark->getBPLOC();
					prev_chrom = mark->getChrom();
				}
			}
			good_markers.push_back(mark);
		}
	}

	return good_markers;
}

vector<int> Helpers::findValidMarkersIndexes(vector<Methods::Marker*>* marks, StepOptions options){
	int prev_base = 0;
	int prev_chrom = -1;
	vector<int> good_markers;
	int msize = marks->size();
	for(int m = 0; m < msize; m++){
		Methods::Marker* mark = (*marks)[m];
		if(mark->isEnabled()){
			if(options.getAutosomeOnly() && ((mark->getChrom() >= opts::_CHRX_) || (mark->getChrom() < 1))){
				continue;
			}

			if(options.doChrom()){
				if(!options.checkChrom(mark->getChrom())){
					continue;
				}
				if(!options.checkBp(mark->getBPLOC())){
					continue;
				}
			}
			if(options.doBpSpace()){
				if(prev_base == 0){
					prev_base = mark->getBPLOC();
					prev_chrom = mark->getChrom();
				}
				else{
					if(mark->getChrom() == prev_chrom && ((mark->getBPLOC() - prev_base) < options.getBpSpace())){
						continue;
					}
					prev_base = mark->getBPLOC();
					prev_chrom = mark->getChrom();
				}
			}
			good_markers.push_back(m);
		}
	}

	return good_markers;
}

vector<int> Helpers::findValidMarkersIndexes(vector<Methods::Marker*>* marks, StepOptions* options){
	int prev_base = 0;
	int prev_chrom = -1;
	vector<int> good_markers;
	int msize = marks->size();
	for(int m = 0; m < msize; m++){
		Methods::Marker* mark = (*marks)[m];
		if(mark->isEnabled()){
			if(options->getAutosomeOnly() && ((mark->getChrom() >= opts::_CHRX_) || (mark->getChrom() < 1))){
				continue;

			}

			if(options->doChrom()){
				if(!options->checkChrom(mark->getChrom())){
					continue;
				}
				if(!options->checkBp(mark->getBPLOC())){
					continue;
				}
			}
			if(options->doBpSpace()){
				if(prev_base == 0){
					prev_base = mark->getBPLOC();
					prev_chrom = mark->getChrom();
				}
				else{
					if(mark->getChrom() == prev_chrom && ((mark->getBPLOC() - prev_base) < options->getBpSpace())){
						continue;
					}
					prev_base = mark->getBPLOC();
					prev_chrom = mark->getChrom();
				}
			}
			good_markers.push_back(m);
		}
	}

	return good_markers;
}

vector<vector<Sample*> > Helpers::generateSampleSets(DataSet* ds, StepOptions* options){
	vector<vector<Sample*> > sample_lists;

	int numsets = options->getSetsSamps();
	int numsamps = options->getRandSamps();
	bool repeats = options->getRandSampsRepeat();
	vector<float> percents = options->getPercentSamps();
	float per_cases = options->getPercentCases();
	float per_controls = options->getPercentControls();

	vector<Sample*> cases = ds->get_affected_vector();
	vector<Sample*> controls = ds->get_unaffected_vector();
	vector<Sample*>* samples = ds->get_samples();

	vector<bool> case_used(cases.size(), false);
	vector<bool> controls_used(controls.size(), false);
	vector<bool> samples_used(samples->size(), false);

	int numused = 0;
	int ssize = ds->num_inds();
	int enabled_size = 0;
	for(int i =0; i < (int)samples->size(); i++){
		if((*samples)[i]->isEnabled()){
			enabled_size++;
		}
		else{
			samples_used[i] = true;
			numused++;
		}
	}

	int enabled_cases = 0;
	for(int i = 0; i < (int)cases.size(); i++){
		if(cases[i]->isEnabled()){
			enabled_cases++;
		}
		else{
			case_used[i] = true;
		}
	}

	int enabled_controls = 0;
	for(int i = 0; i < (int)controls.size(); i++){
		if(controls[i]->isEnabled()){
			enabled_controls++;
		}
		else{
			controls_used[i] = true;
		}
	}

	bool cases_controls = false;
	if(per_cases > 0 || per_controls > 0){
		cases_controls = true;
	}

	vector<int> samps_per_set;
	int num_cases = 0;
	int num_controls = 0;
	if(percents.size() > 0){
		for(int i = 0; i < (int)percents.size(); i++){
			samps_per_set.push_back((int)(enabled_size * (percents[i] / 100.0f)));
		}
	}
	else{
		for(int i = 0; i < numsets; i++){
			samps_per_set.push_back(numsamps);
		}
	}

	if(per_cases > 0 && numsamps > 0){
		num_cases = (int)(numsamps * (per_cases / 100.0f));
		if(num_cases > enabled_cases){
			num_cases = enabled_cases;
		}
		if(num_cases > numsamps){
			num_cases = numsamps;
		}
	}
	if(per_controls > 0 && numsamps > 0){
		num_controls = (int)(numsamps * (per_controls / 100.0f));
		if(num_controls > enabled_controls){
			num_controls = enabled_controls;
		}
		if(num_controls > numsamps){
			num_controls = numsamps;
		}
	}
	if((num_cases + num_controls) > numsamps){
		throw MethodException("Number of cases + controls exceeds sample count.\n");
	}

	int random;
	sample_lists.resize(numsets);
	for(int i = 0; i < numsets; i++){
		vector<bool> set_used(ssize, false);
		if(num_cases > 0){
			int used_cases = 0;
			while(used_cases < num_cases && (int) sample_lists[i].size() < numsamps){
				random = int (rand() % ssize);
				if(!samples_used[random] && !set_used[random] && (*samples)[random]->getAffected()){
					sample_lists[i].push_back((*samples)[random]);
					set_used[random] = true;
					used_cases++;
					if(!repeats){
						samples_used[random] = true;
					}
				}
			}
		}
		if(num_controls > 0){
			int used_controls = 0;
			while(used_controls < num_controls && (int) sample_lists[i].size() < numsamps){
				random = int (rand() % ssize);
				if(!samples_used[random] && !set_used[random] && !(*samples)[random]->getAffected()){
					sample_lists[i].push_back((*samples)[random]);
					set_used[random] = true;
					used_controls++;
					if(!repeats){
						samples_used[random] = true;
					}
				}
			}
		}
		if(num_cases > 0 && num_controls > 0){
			continue;
		}
		while((int) sample_lists[i].size() < numsamps && numused < ssize){
			random = int(rand() % ssize);
			if(!samples_used[random] && !set_used[random]){
				if(num_cases > 0 && (*samples)[random]->getAffected()){
					continue;
				}
				if(num_controls > 0 && !(*samples)[random]->getAffected()){
					continue;
				}
				sample_lists[i].push_back((*samples)[random]);
				set_used[random] = true;
				if(!repeats){
			    	samples_used[random] = true;
			    }
			}
		}
	}
	return sample_lists;
}


vector<Methods::Marker*> Helpers::findRandomMarkers(vector<Methods::Marker*> good_markers, vector<Methods::Marker*>* used_markers, StepOptions* options){
	vector<Methods::Marker*> marks;
	int num = options->getRandomMarkers();
	int random;

	while((int) marks.size() < num && used_markers->size() < good_markers.size()){
		random = int(rand() % good_markers.size());
		Methods::Marker* mark = (good_markers)[random];
		vector<Methods::Marker*>::iterator found = find(used_markers->begin(), used_markers->end(), mark);
		if(found == used_markers->end()){
			used_markers->push_back(mark);
			marks.push_back(mark);
		}
	}
	stable_sort(marks.begin(), marks.end(), less<Methods::Marker*>());
	if(options->doRandomRepeat()){
		used_markers->clear();
	}
	return marks;

}

void Helpers::readCovariateFile(string file, DataSet* ds, StepOptions options, InputFilter* filters){
	map<string, Methods::Sample*> smap;
	map<string, Methods::Sample*>::iterator siter;

	for(int i = 0; i < ds->num_inds(); i++){
		Methods::Sample* samp = ds->get_sample(i);
		smap[samp->getFamID() + "#" + samp->getInd()]= samp;
	}

	vector<bool> use_map;

	ifstream in(file.c_str(), ios::in);
	if(!in){
		opts::printLog("Error opening covariate file: " + file + "\n");
		throw MethodException("Error opening covariate file: " + file + "\n");
	}

	int count = 0;
	while(!in.eof()){
		count++;
		string line;
		getline(in, line);
		if(line == ""){
			continue;
		}
		if(line.at(0) == '#'){
			continue;
		}

		vector<string> elems;
		stringstream ss(line);
		string tok = "";
		while(ss >> tok){
			elems.push_back(tok);
		}

		if(elems.size() <= 2){
			in.close();
			opts::printLog("Covariate file requires 3 or more elements (FamID IndID Cov1 Cov2 etc.)\n");
			throw MethodException("Covariate file requires 3 or more elements(FamID IndID Cov1 Cov2 etc.)\n");
		}
		if(count == 1){//header found
			for(unsigned int i = 2; i < elems.size(); i++){
				bool use = true;
				if(filters != NULL){
					for(int f = 0; f < filters->num_covariate_filters(); f++){
						use = filters->run_covariate_filter(f, elems[i]);
					}
				}
				use_map.push_back(use);
				if(use){
					ds->add_covariate(elems[i]);
				}
				opts::_COVS_FOUND_++;
			}

			continue;
		}

		string famid = elems[0];
		string indid = elems[1];

		Methods::Sample* samp = smap[famid + "#" + indid];
		if(!samp){
			continue;
		}

		for(unsigned int i = 2; i < elems.size(); i++){
			if(use_map[i - 2]){
				double value = -1;
			   	try{
					if(elems[i] == options.getCovarMissing()){
						try{
							for(unsigned int c = 0; c < elems[i].size(); c++){
								if(!isdigit(elems[i][c]) && elems[i][c] != '.' && elems[i][c] != '-'){
									throw "oops!";
								}
							}
							value = (double) atof(options.getCovarMissing().c_str());
						}
						catch(...){
							value = options.getDefaultCovarMissing();
						}
					}
					else{
						value = (double) atof(elems[i].c_str());
					}
				}catch(...){
					throw MethodException("Cannot convert " + elems[i] + " to number on line " + getString<int>(count) + " in file: " + file + "\n");
				}
				samp->addCovariate(value);

			}
		}

	}

	for(int i = 0; i < ds->num_inds(); i++){
		Methods::Sample* samp = ds->get_sample(i);
		if(samp->getCovariateVector().size() == 0){
			for(int c = 0; c < ds->num_covariates(); c++){
				samp->addCovariate(ds->get_missing_covalue());
			}
		}
	}
	in.close();
}


void Helpers::readTraitFile(string file, DataSet* ds, StepOptions options, InputFilter* filters){
	map<string, Methods::Sample*> smap;
	map<string, Methods::Sample*>::iterator siter;

	for(int i = 0; i < ds->num_inds(); i++){
		Methods::Sample* samp = ds->get_sample(i);
		smap[samp->getFamID() + "#" + samp->getInd()]= samp;
	}

	vector<bool> use_map;

	ifstream in(file.c_str(), ios::in);
	if(!in){
		opts::printLog("Error opening trait file: " + file + "\n");
		throw MethodException("Error opening trait file: " + file + "\n");
	}

	int count = 0;
	while(!in.eof()){
		count++;
		string line;
		getline(in, line);
		if(line == ""){
			continue;
		}
		if(line.at(0) == '#'){
			continue;
		}

		vector<string> elems;
		stringstream ss(line);
		string tok = "";
		while(ss >> tok){
			elems.push_back(tok);
		}

		if(elems.size() <= 2){
			in.close();
			opts::printLog("Trait file requires 3 or more elements (FamID IndID Trait1 Trait2 etc.)\n");
			throw MethodException("Trait file requires 3 or more elements(FamID IndID Trait1 Trait2 etc.)\n");
		}

		//2-22-2011 added this to check for header before reading in trait file...
		//If no header line present, process first row as a header, and still process as a normal row
		string firstToken = elems[0];	//if this is the header line, this should be equal to some form of 'FamID'
		if(count == 1)
		{
			if(!(boost::iequals(firstToken, "famid") || boost::iequals(firstToken, "fid")))
			{
				opts::printLog("Trait file must have a header line (FamID IndID Trait1 Trait2, etc.)");
				throw MethodException("Trait file must have a header line (FamID IndID Trait1 Trait2, etc.)");
			}

			for(unsigned int i = 2; i < elems.size(); i++)
			{
				bool use = true;
				if(filters != NULL)
				{
					for(int f = 0; f < filters->num_trait_filters(); f++)
					{
						use = filters->run_trait_filter(f, elems[i]);
					}
				}
				use_map.push_back(use);
				if(use)
				{
					ds->add_trait(elems[i]);
				}
				opts::_TRAITS_FOUND_++;
			}
			continue;
		}

		string famid = elems[0];
		string indid = elems[1];

		Methods::Sample* samp = smap[famid + "#" + indid];
		if(!samp){
			continue;
		}

		for(unsigned int i = 2; i < elems.size(); i++)
		{
			if(use_map[i - 2])
			{
				double value = -1;
			   	try
			   	{
					if(elems[i] == options.getTraitMissing())
					{
						try
						{
							for(unsigned int c = 0; c < elems[i].size(); c++)
							{
								if(!isdigit(elems[i][c]) && elems[i][c] != '.' && elems[i][c] != '-')
								{
									throw "oops!";
								}
							}
							value = (double) atof(options.getTraitMissing().c_str());
						}
						catch(...)
						{
							value = options.getDefaultTraitMissing();
						}
					}
					else
					{
						value = (double) atof(elems[i].c_str());
					}
				}
			   	catch(...)
			   	{
					throw MethodException("Cannot convert " + elems[i] + " to number on line " + getString<int>(count) + " in file: " + file + "\n");
				}
				samp->addTrait(value);

			}
		}

	}

	for(int i = 0; i < ds->num_inds(); i++){
		Methods::Sample* samp = ds->get_sample(i);
		if(samp->getTraitVector().size() == 0){
			for(int c = 0; c < ds->num_traits(); c++){
				samp->addTrait(ds->get_missing_covalue());
			}
		}
	}
	in.close();
}

//Added 02-23-2011 to support new -update-ids function
void Helpers::readIDFile(string file, DataSet* ds)
{
	map <string, string> sampleMap;
	map<string, string>::iterator siter;

	ifstream in(file.c_str(), ios::in);
	if(!in)
	{
		opts::printLog("Error opening ID file: " + file + "\n");
		throw MethodException("Error opening ID file: " + file + "\n");
	}

	//loop through the ID input file, building up a map of OldFamID#OldIndID->NewFamID#NewIndID
	int count = 0;
	while(!in.eof())
	{
		count++;
		string line;
		getline(in, line);
		if(line == "")
		{
			continue;
		}
		if(line.at(0) == '#')
		{
			continue;
		}

		vector<string> elems;
		stringstream ss(line);
		string tok = "";
		while(ss >> tok)
		{
			elems.push_back(tok);
		}

		if(elems.size() != 4)
		{
			in.close();
			opts::printLog("ID file requires 4 elements (OldFamID OldIndID NewFamID NewIndID)\n");
			throw MethodException("ID file requires 4 elements (OldFamID OldIndID NewFamID NewIndID)\n");
		}
		if(count == 1)
		{
			string firstToken = elems[0];
			boost::to_lower(firstToken);
			if(boost::contains(firstToken, "famid"))
			{
				//this is a header line, skip it during processing
				continue;
			}
		}

		//build up the strings to place in the map...
		string oldIDString;
		string newIDString;
		oldIDString = elems[0] + " " + elems[1];
		newIDString = elems[2] + " " + elems[3];
		sampleMap[oldIDString] = newIDString;
	}

	//get pointer to samples vector from dataset
	vector<Sample*>* samples = ds->get_samples();

	//loop through samples vector, updating IDs as we go...
	string currentFamID;
	string currentIndID;
	string searchTerm;
	string newIDString;
	string newFamID;
	string newIndID;
	string currToken;
	for(unsigned int i = 0; i < samples->size(); i++)
	{
		Sample* currSample = samples->at(i);
		currentFamID = currSample->getFamID();
		currentIndID = currSample->getInd();

		searchTerm = currentFamID + " " + currentIndID;
		siter = sampleMap.find(searchTerm);
		vector<string> ids;

		string tok = "";
		if(siter != sampleMap.end())
		{
			//this sample was in the ID list file, update the id's...
			newIDString = sampleMap[searchTerm];
			stringstream idStream(newIDString);
			while (idStream >> tok)
			{
				ids.push_back(tok);
			}
			if(ids.size() < 2)
			{
				opts::printLog("Error processing ID file, Missing NewFamID or NewIndID\n");
				throw MethodException("Error processing ID file, Missing NewFamID or NewIndID\n");
			}
			newFamID = ids.at(0);
			newIndID = ids.at(1);

			if(newFamID != "")
			{
				currSample->setInd(newIndID);
			}
			if(newIndID != "")
			{
				currSample->setFamID(newFamID);
			}
		}
	}
}

string Helpers::stringToLowerCase(string in)
{
	string str = "";
	cout << "in: " << in << endl;
	for (int i=0;i<(int)in.length();i++)
	{
		cout << "in[i]: " << endl;
		try
		{
			str[i] = tolower(in[i]);
		}
		catch (...)
		{
			//there was a problem with tolower() above, just copy the character over...
			str[i] = in[i];
		}
		cout << "str[i]: " << str[i] << endl;
	}
	cout << "str: " << str << endl;
	return str;
}


void Helpers::readCovariates(string file, vector<Methods::Sample*>* samples, vector<string>* covheaders){
	map<string, Methods::Sample*> smap;
	map<string, Methods::Sample*>::iterator siter;

	for(unsigned int i = 0; i < samples->size(); i++){
		Methods::Sample* samp = (*samples)[i];
		smap[samp->getFamID() + "#" + samp->getInd()]= samp;
	}

	ifstream in(file.c_str(), ios::in);
	if(!in){
		opts::printLog("Error opening covariate file: " + file + "\n");
		//exit(1);
		throw MethodException("Error opening covariate file: " + file + "\n");
	}
	int count = 0;
	int numcovs = 0;
	while(!in.eof()){
		count++;
		string line;
		getline(in, line);
		if(line == ""){
			continue;
		}
		if(line.at(0) == '#'){
			continue;
		}

		vector<string> elems;
		stringstream ss(line);
		string tok;
		while(ss >> tok){
		   elems.push_back(tok);
		}
		if(elems.size() <= 2){
			opts::printLog("Line: " + getString<int>(count) + " has incorrect number of columns!");
			in.close();
			throw MethodException("Line: " + getString<int>(count) + " has incorrect number of columns!");
		}
		if(numcovs == 0){
			numcovs = elems.size() - 2;
			if(covheaders == NULL){
				opts::cov_loc.resize(numcovs);
				for(int c = 0; c < numcovs; c++){
					opts::cov_loc[c] = "COV"+getString<int>(c+1);
				}
			}
			else{
				covheaders->resize(numcovs);
				for(int c = 0; c < numcovs; c++){
					(*covheaders)[c] = "COV"+getString<int>(c+1);
				}
			}
		}
		else if((int) elems.size() != numcovs + 2){
			opts::printLog("Line: " + getString<int>(count) + " has incorrect number of columns!");
			in.close();
			throw MethodException("Line: " + getString<int>(count) + " has incorrect number of columns!");
		}
		string fam = elems[0];
		string ind = elems[1];

		siter = smap.find(fam +"#"+ind);
		if(siter != smap.end()){
			Methods::Sample* samp = siter->second;

			samp->resizeCovariates(numcovs);
			for(unsigned int c = 2; c < elems.size(); c++){
				try{
					samp->setCovariate(atof(elems[c].c_str()), (c - 2));
				}catch(...){
					opts::printLog("Column " + getString<int>(c) +" on line: " + getString<int>(count) + " is not a number!?\n");
					throw MethodException("Column " + getString<int>(c) +" on line: " + getString<int>(count) + " is not a number!?\n");
				}
			}
		}
		else if(fam == "FamID" && ind == "IndID"){
			for(unsigned int c = 2; c < elems.size(); c++){
				if(covheaders == NULL){
					opts::cov_loc[c-2] = elems[c];
				}
				else{
					(*covheaders)[c-2] = elems[c];
				}
			}
		}
	}
	in.close();
}

/*
 *Function: readString
 *Description:
 *Reads the next string from a file, space/tab delimited
 *
 *
 */
bool Helpers::readString(FILE* fp, string* s){
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

void Helpers::readZeroGenoFile(string file){
	ifstream PED;
	PED.open(file.c_str());
	if(!PED){
		throw MethodException("Error opening zero genotype file: " + file + ".  Exiting!\n");
	}
	PED.clear();
		string fam = "";
		string ind = "";
		string snp = "";
		string line = "";
	while(getline(PED, line)){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 3){
			throw MethodException("Zero genotype file column size != 3: " + line + " Exitting!!\n");
		}
		fam = tokens[0];
		ind = tokens[1];
		snp = tokens[2];

		opts::zerogenoinfo[(fam + "#" + ind)].push_back(snp);
	}

	PED.clear();
	PED.close();

}

/*
 * Function: zeroSingleGenos
 * Description:
 * Zeros single genotypes based on sample and snp per the file specified by -zero-genos file.txt
 *
 */
void Helpers::zeroSingleGenos(vector<Methods::Marker*>* markers, vector<Methods::Sample*>* samples){
	opts::printLog("Zeroing out requested genotypes from file: " + opts::_ZEROGENOFILE_ + "\n");
	//create sample map
	map<string, Methods::Sample*> smap;
	for(unsigned int i = 0; i < samples->size(); i++){
		Methods::Sample* samp = (*samples)[i];
		smap[samp->getFamID() + "#" + samp->getInd()] = samp;
	}

	map<string, vector<string> >::iterator iter;
	map<string, Methods::Sample* >::iterator sampiter;
	vector<Methods::Marker*>::iterator markiter;

	int samples_done = 0;
	int genos_done = 0;

	for(iter = opts::zerogenoinfo.begin(); iter != opts::zerogenoinfo.end(); iter++){
		string key = iter->first;
		vector<string> value = iter->second;
		sampiter = smap.find(key);
		if(sampiter != smap.end()){
			Methods::Sample* mysamp = sampiter->second;
			bool done = false;
			for(unsigned int m = 0; m < value.size(); m++){
				string snp = value[m];
				markiter = find_if(markers->begin(), markers->end(), FindMarker(snp));
				if(markiter != markers->end()){
					Methods::Marker* mark = *markiter;
					int mloc = mark->getLoc();
					if(!mark->isMicroSat()){
						mysamp->addAone(mloc,true);
						mysamp->addAtwo(mloc,false);
					}
					else{
						mysamp->addAbone(mloc, -1);
						mysamp->addAbtwo(mloc, -1);
					}
					if(!done){
						samples_done++;
						done = true;
					}
					genos_done++;
				}
			}
		}
	}

	opts::printLog("Total samples that had at least 1 genotype zeroed: " + getString<int>(samples_done) + "\n");
	opts::printLog("Total genotypes that were zeroed: " + getString<int>(genos_done) + "\n");
}

/*
 *Function: assignLinks
 *Description:
 *Creates the family structure by assigning each sample their respective
 *parents, siblings and children
 *
 */
void Helpers::assignLinks(vector<Methods::Family*>* families){
    vector<Methods::Family*>::iterator f_iter;
	int fsize = families->size();

	for(int f = 0; f < fsize; f++){
		Methods::Family* fam = (*families)[f];
		vector<Methods::Sample*>* samps = fam->getSamples();
		bool good = false;
		bool excluded = false;
		int ssize = samps->size();

		for(int s = 0; s < ssize; s++){
			Methods::Sample* samp = (*samps)[s];
			if(samp->isEnabled()){
				good = true;
			}
			else if(samp->isExcluded()){
				excluded = true;
			}
			if(samp->getDadID() == "0" && samp->getMomID() == "0"){
				fam->addFounder(samp);
				samp->setFounder(true);
			}
			else{
				fam->addNonFounder(samp);
				if(samp->getDadID() != "0"){
                	vector<Methods::Sample*>::iterator temp_iter = find_if(samps->begin(), samps->end(), FindSampleByID(samp->getDadID()));
                	if(temp_iter != samps->end()){
						Methods::Sample* dad = (*temp_iter);
	 					samp->setDad(dad);
						Methods::Sample* last_child = dad->getLastChild();
						if(last_child != NULL && last_child->getSib() == NULL){
							if(last_child->getSib() == NULL){
								last_child->setSib(samp);
							}
							last_child->setPatSib(samp);
						}
						dad->addChild(samp);
                	}
				}
				if(samp->getMomID() != "0"){
                	vector<Methods::Sample*>::iterator temp_iter = find_if(samps->begin(), samps->end(), FindSampleByID(samp->getMomID()));
                	if(temp_iter != samps->end()){
						Methods::Sample* mom = (*temp_iter);
	 					samp->setMom(mom);
						Methods::Sample* last_child = mom->getLastChild();
						if(last_child != NULL && last_child->getSib() == NULL){
							if(last_child->getSib() == NULL){
								last_child->setSib(samp);
							}
							last_child->setMatSib(samp);
						}
						mom->addChild(samp);
                	}
				}
			}
		}
		fam->setEnabled(good);
		if(!good && excluded){
			fam->setExcluded(excluded);
		}
	}

}

/*
 *Function: reorderAlleles
 *Return: none
 *Parameters: vector of samples, vector of markers
 *Description:
 *Calculates the raw allele counts and remaps the sample bool vector 1 to be the
 *minor allele
 *
 */
void Helpers::reorderAlleles(vector<Methods::Sample*>* samples, vector<Methods::Marker*>* markers){
	int msize = markers->size();
	int ssize = samples->size();

	for(int m = 0; m < msize; m++){
		int mloc = (*markers)[m]->getLoc();
		if((*markers)[m]->getNumAlleles() == 0){
			(*markers)[m]->addAllele("0");
			(*markers)[m]->addAllele("0");
		}
		if((*markers)[m]->getNumAlleles() == 1){
			(*markers)[m]->addAllele("0");
		}
		if((*markers)[m]->getNumAlleles() <= 2){
			int a1 = 0;
			int a2 = 0;
			for(int s = 0; s < ssize; s++){
				if((*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc) && (*samples)[s]->getAmissing(mloc)){
					continue;
				}
				if((*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc)){
					a2+=2;
				}
				else if(!(*samples)[s]->getAone(mloc) && !(*samples)[s]->getAtwo(mloc)){
					a1+=2;
				}
				else if(!(*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc)){
					a1++;
					a2++;
				}
			}
			if(a2 < a1 && (*markers)[m]->getAllele2() != "0"){
				string temp = (*markers)[m]->getAllele1();
				(*markers)[m]->resetAllele1((*markers)[m]->getAllele2());
				(*markers)[m]->resetAllele2(temp);
				for(int s = 0; s < ssize; s++){
					if((*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc) && !(*samples)[s]->getAmissing(mloc)){
						(*samples)[s]->addAone(mloc, false);
						(*samples)[s]->addAtwo(mloc, false);
					}
					else if(!(*samples)[s]->getAone(mloc) && !(*samples)[s]->getAtwo(mloc)){
						(*samples)[s]->addAone(mloc, true);
						(*samples)[s]->addAtwo(mloc, true);
					}
				}
			}
			// set referent allele to be minor allele unless it has already been set by map file
			if((*markers)[m]->getReferent() == ""){
			  (*markers)[m]->setReferent((*markers)[m]->getAllele1());
			}
			// set referent index in Marker for quicker lookup when running analysis
			if((*markers)[m]->getReferent() == (*markers)[m]->getAllele1()){
			  (*markers)[m]->setReferentIndex(0);
			}

			else if((*markers)[m]->getReferent() == (*markers)[m]->getAllele2()){
			  (*markers)[m]->setReferentIndex(1);
			}
			else{
			  throw MethodException("No match for referent allele specified in map file for marker " +
			    (*markers)[m]->getRSID());
			}

		}
	}

}

/*
 * Function: map_allele
 * Description:
 * Maps a string allele to/from ACGT to/from 1234 or to 1/2 or performs custom mapping per user specifications
 */
string Helpers::map_allele(Methods::Marker* mark, string a, StepOptions *options){
	if(options->doAllele1234()){
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
	}
	if(options->doAlleleACGT()){
		if(a == "1"){
			return "A";
		}
		else if(a == "2"){
			return "C";
		}
		else if(a == "3"){
			return "G";
		}
		else if(a == "4"){
			return "T";
		}
	}
	if(options->doAlleleCustom()){
		return options->getAlleleMapping(a);
	}
	if(options->doAllele12()){
		if(mark->getAllele1() == a){
			return "1";
		}
		if(mark->getAllele2() == a){
			return "2";
		}
	}

	return a;
}

/*
 * Function: readCustomAlleles
 * Description:
 * Reads file containing custom allele mappings for conversion
 * Column1 = From column
 * Column2 = To column
 */
map<string, string> Helpers::readCustomAlleles(string f){
	map<string, string> mapping;
	ifstream PED;
	PED.open(f.c_str());
	if(!PED){
		opts::printLog("Error opening custom allele mapping file: " + f + ".  Exiting!\n");
		throw MethodException("Error opening custom allele mapping file: " + f + ".  Exiting!\n");
	}
	PED.clear();
		string from = "";
		string to = "";
		string line = "";
	while(getline(PED, line)){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 2){
			opts::printLog("Custom allele mapping file column size != 2: " + line + " Exiting!!\n");
			throw MethodException("Custom allele mapping file column size != 2: " + line + " Exiting!!\n");
		}
		from = tokens[0];
		to = tokens[1];

		mapping[from] = to;
	}

	PED.clear();
	PED.close();
	return mapping;
}

/*
 *Function: remapSamples
 *Description:
 *If -micro-sats is enabled, when a 3rd allele is found for a marker, this function
 *moves the Sample genotype storage for the specific marker from boolean storage to integer storage.
 *
 */
void Helpers::remapSamples(vector<Methods::Sample*>* samples, vector<Methods::Marker*>* markers, vector<int>* marker_map, int loc){

	for(unsigned int i = 0; i < samples->size(); i++){
		Methods::Sample* samp = (*samples)[i];
		if(!samp->haveMicroSat(loc)){
			samp->addMicroSat(loc);
			if(!samp->getAone(loc) && !samp->getAtwo(loc)){
				samp->addAbone(loc, 0);
				samp->addAbtwo(loc, 0);
			}
			else if(samp->getAone(loc) && samp->getAtwo(loc) && !samp->getAmissing(loc)){
				samp->addAbone(loc, 1);
				samp->addAbtwo(loc, 1);
			}
			else if(!samp->getAone(loc) && samp->getAtwo(loc)){
				samp->addAbone(loc, 0);
				samp->addAbtwo(loc, 1);
			}
			else{
				samp->addAbone(loc, -1);
				samp->addAbtwo(loc, -1);
			}
			samp->addAone(loc, true);
			samp->addAtwo(loc, true);
			samp->addAmissing(loc, true);
		}
	}
}

/*
 *Function: readPed
 *Description:
 *Reads ped file and stores data into sample, family, and marker vectors
 *@throws MethodException
 */
void Helpers::readPedM(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, vector<Methods::Marker*>* markers, vector<int>* marker_map, StepOptions options){
	map<string, vector<string> > descinfo;
	vector<string> exclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> fexclude;
	vector<string> finclude;
	vector<string> descheaders;

    FILE* input;
    input = fopen(options.getPedFile().c_str(), "r");
	if(!input){
		throw MethodException("Error opening pedfile: " + options.getPedFile() + ".");
		return;
	}
	int onind = -1;
    while(!feof(input)){
		onind++;
        Methods::Sample* samp = new Sample();
        int f = 0;
        string temp = "";
        if(readString(input, &temp)){
           samp->setFamID(temp);
              f++;
             temp = "";
         }

		string ftemp = samp->getFamID();
		if(samp->getFamID() == ""){
			delete(samp);
			continue;
		}
        if(ftemp.at(0) == '#'){
			delete(samp);
			while(fgetc(input) != '\n' && !feof(input)){}
            continue;
        }
        string sex = "";
        string pheno = "";
        if(readString(input, &temp)){
            samp->setInd(temp);
            f++;
            temp = "";
        }
		samp->setEnabled(true);

		if(exclude.size() > 0){
			vector<string>::iterator found = find(exclude.begin(), exclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found != exclude.end()){
				if(!opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					continue;
				}
				else{
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
		}
		if(sinclude.size() > 0){
			vector<string>::iterator found = find(sinclude.begin(), sinclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found == sinclude.end()){
				if(!opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					continue;
				}
				else{
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
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

		if(pheno == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		float mypheno = atof(pheno.c_str());
		samp->setPheno(mypheno);
        if(opts::pedinfo.size() > 0){
            map<string, Methods::Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
            if(sfind != opts::pedinfo.end()){
				Methods::Sample* sfound = sfind->second;
                samp->setDadID(sfound->getDadID());
                samp->setMomID(sfound->getMomID());
                samp->setSex(sfound->getSex());
                samp->setPheno(sfound->getPheno());
            }
        }
		if(samp->getPheno() != 0.0f && samp->getPheno() != 1.0f && samp->getPheno() != 2.0f){
			opts::_BINTRAIT_ = false;
		}

		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(unsigned int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					samp->assignDetail(descheaders[i], "NA");
				}
			}

		}
        samp->resizeAlleles(markers->size());

        int gn = 0;
        int i = 0;
        bool linedone = false;
        //bool fatal = false;

        string fmsg;
        while(!linedone){
            string one = "";
            string two = "";

            while(1)
            {
                char ch = fgetc(input);

                if(ch == '/' || ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input))
                {
                    if(ch == '\n' || ch == '\r' || feof(input))
                    {
                        linedone = true;
                    }

                    if(one.length() > 0)
                    {
                        gn++;
                        break;
                    }
                    if(ch == '\n' || ch == '\r' || feof(input))
                    {
                        break;
                    }
                }
                else
                {
                    one += ch;
                }
            }
            if(!linedone)
            {
                while(1)
                {
                    char ch = fgetc(input);
                    if(ch == '/' || ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input))
                    {
                        if(ch == '\n' || ch == '\r' || feof(input))
                        {
                            linedone = true;
                        }
                        if(two.length() > 0)
                        {
                            gn++;
                            break;
                        }
                        if(ch == '\n' || ch == '\r' || feof(input))
                        {
                            break;
                        }
                    }
                    else
                    {
                        two += ch;
                    }
                }
                if(linedone && one.length() == 0 && two.length() == 0)
                {
                    break;
                }

				if(i > (int)markers->size()){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>((2 * markers->size()) + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
				Methods::Marker* m = (*markers)[(*marker_map)[i]];
				if(m->isEnabled()){
					int oldallelecount = m->getNumAlleles();
	                if(one != "0"){
						//new
						if(m->getAlleleLoc(one) < 0){
							m->addAllele(one);
						}

        	        }

            	    if(two != one){
                	    if(two != "0"){
							//new
							if(m->getAlleleLoc(two) < 0){
								m->addAllele(two);
							}

                    	}
	                }


					if(m->getNumAlleles() <= 2){
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
					}
					else if(opts::_MICROSATS_){
						samp->addMicroSat(i);
						int loc1 = m->getAlleleLoc(one);
						int loc2 = m->getAlleleLoc(two);

						samp->addAbone(i, loc1);
						samp->addAbtwo(i, loc2);
						if(oldallelecount <= 2){
							remapSamples(samples, markers, marker_map, i);
						}
					}
					else if(m->getNumAlleles() > 2 && !opts::_MICROSATS_){
						opts::printLog("More than 2 unique alleles found for map location: " + getString<int>(i) + ", line: " + getString<int>(onind + 1) + ".  Microsatellites not specified.\n");
						throw MethodException("More than 2 unique alleles found for map location: " + getString<int>(i) + ", line: " + getString<int>(onind + 1) + ".  Microsatellites not specified.\n");
					}
				}
				else{
					samp->addAone(i,true);
					samp->addAtwo(i, false);
				}
    	    	i++;
				if(i > (int)markers->size()){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>((2 * markers->size()) + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
            }/*end !linedone*/
        }/*end while(1)*/
		if(gn != (int)(2* markers->size())){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>(((2 * markers->size()) + 6));
					text += " columns but found ";
					text += getString<int>((f + gn));
					text += "\n";
					throw MethodException(text);
		}

        vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Methods::Family* fam = new Family();
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
        }
		if(fexclude.size() > 0){
			vector<string>::iterator found = find(fexclude.begin(), fexclude.end(), samp->getFamID());
			if(found != fexclude.end()){
				samp->setEnabled(false);
				vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		if(finclude.size() > 0){
			vector<string>::iterator found = find(finclude.begin(), finclude.end(), samp->getFamID());
			if(found == finclude.end()){
				samp->setEnabled(false);
				vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
    }/*end while(eof)*/

    fclose(input);
}

/*
 * Function: removeBeginWhiteSpace
 * Removes all white space at beginning of a string
 *
 * return: cleaned string
 */
string Helpers::removeBeginWhiteSpace(string l){
	if(l.size() == 0){
		return l;
	}
	while((l.at(0) == ' ' || l.at(0) == '\t') && l.size() > 0){
		l.erase(0);
	}
	return l;
}

/*
 * Function: readSampleFile
 * Reads a file containing a list of Samples (ie: FamID & IndID only) and places each one as a Sample class into the specified vector
 * return: void
 */
void Helpers::readSampleFile(string file, vector<Methods::Sample*>* mlist){
	if(mlist == NULL){
		throw MethodException("NULL vector passed to readSampleFile...\n");
	}
	ifstream input;
	input.open(file.c_str());
	if(!input){
		opts::printLog("Error opening sample list file: " + file + ".  Exiting!\n");
		throw MethodException("Error opening sample list file: " + file + ".\n");
	}
	string line = "";
	while(getline(input, line)){
		line = removeBeginWhiteSpace(line);
		if(line.size() > 0){
			if(line[0] == '#'){
				continue;
			}
		}
		else{
			continue;
		}
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 2){
			opts::printLog("Sample list file column size != 2: " + line + "\n");
			throw MethodException("Sample list file column size != 2: " + line + "\n");
		}
		Methods::Sample* m = new Sample();
		m->setFamID(tokens[0]);
		m->setInd(tokens[1]);
		mlist->push_back(m);
	}
	input.close();
}

/*
 * Function: readFamilyFile
 * Reads a file containing a list of Family IDs and places each one as a Family class into the specified vector
 * return: void
 */
void Helpers::readFamilyFile(string file, vector<Methods::Family*>* mlist){
	if(mlist == NULL){
		throw MethodException("NULL vector passed to readFamilyFile...\n");
	}
	ifstream input;
	input.open(file.c_str());
	if(!input){
		opts::printLog("Error opening family list file: " + file + ".  Exiting!\n");
		throw MethodException("Error opening family list file: " + file + ".\n");
	}
	string line = "";
	while(getline(input, line)){
		line = removeBeginWhiteSpace(line);
		if(line.size() > 0){
			if(line[0] == '#'){
				continue;
			}
		}
		else{
			continue;
		}
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 1){
			opts::printLog("Family list file column size != 1: " + line + "\n");
			throw MethodException("Family list file column size != 1: " + line + "\n");
		}
		Methods::Family* m = new Family();
		m->setFamID(tokens[0]);
		mlist->push_back(m);
	}
	input.close();
}

void Helpers::readCovTraitFile(string file, vector<string>* list){
	if(list == NULL){
		throw MethodException("NULL vector passed to readCovTraitFile...\n");
	}

	ifstream input;
	input.open(file.c_str());
	if(!input){
		opts::printLog("Error opening cov/trait list file: " + file + ". Exiting!\n");
		throw MethodException("Error opening cov/trait list file: " + file + ".\n");
	}

	string line = "";
	while(getline(input, line)){
		line = removeBeginWhiteSpace(line);
		if(line.size() > 0){
			if(line[0] == '#'){
				continue;
			}
		}
		else{
			continue;
		}
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 1){
			opts::printLog("Cov/Trait list file column size != 1: " + line + "\n");
			throw MethodException("Cov/Trait list file column size != 1: " + line + "\n");
		}
		list->push_back(tokens[0]);
	}
	input.close();
}

/*
 * Function: readLocusFile
 * Reads a file containing a list of loci (ie: rsid only) and places each one as a Marker class into the specified vector
 * return: void
 */
void Helpers::readLocusFile(string file, vector<Methods::Marker*>* mlist){
	if(mlist == NULL){
		throw MethodException("NULL vector passed to readLocusFile...\n");
	}
	ifstream input;
	input.open(file.c_str());
	if(!input){
		opts::printLog("Error opening locus list file: " + file + ".  Exiting!\n");
		throw MethodException("Error opening locus list file: " + file + ".\n");
	}
	string line = "";
	while(getline(input, line)){
		line = removeBeginWhiteSpace(line);
		if(line.size() > 0){
			if(line[0] == '#'){
				continue;
			}
		}
		else{
			continue;
		}
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 1){
			opts::printLog("Locus list file column size != 1: " + line + "\n");
			throw MethodException("Locus list file column size != 1: " + line + "\n");
		}
		Methods::Marker* m = new Marker();
		m->setRSID(tokens[0]);
		mlist->push_back(m);
	}
	input.close();
}

/*
 * Function: readFreqFile
 * Reads a frequency file in order to set pre-defined minor allele frequencies.  Format of the file is
 * snpid followed by the frequency (2 columns).  Returns a map<string, float>.
 */

map<string, float> Helpers::readFreqFile(string file){
	map<string, float> freqs;

	ifstream finput;
	finput.open(file.c_str(), ios::in);

	if(!finput){
		throw MethodException("Error opening marker frequency file: " + file + "\n");
	}

	string line = "";

	int freqline = 1;
	while(getline(finput, line)){
		if(line.size() == 0){
			continue;
		}
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 2){
			throw MethodException("Marker frequency file column size != 2 on line: " + line + "\n");
		}
		if(tokens[0].at(0) == '#'){
			continue;
		}

		float freq = -1.0f;
		istringstream b(tokens[1]);
		if(!(b >> freq)){
			throw MethodException(tokens[1] + " is not a valid number on line: " + getString<int>(freqline));
		}
		freqs[tokens[0]] = freq;
		freqline++;
	}

	if(finput){
		finput.close();
	}

	return freqs;
}

/*
 * Function: readMapM
 * Reads a map file and filters based on vector of filter functions passed.  No filtering if vector is null or empty
 * Filter function needs to have as parameters for Markers: Marker and vector of markers to compare against.
 * Filter for covs & traits takes string name and vector of traits/covs to compare against.
 *
 * Return: bool
 */

void Helpers::readMapM(DataSet* ds, StepOptions options,  InputFilter* filters){
	vector<Methods::Marker*>* markers = ds->get_markers();
	vector<int>* marker_map = ds->get_marker_map();
	vector<string>* covariates = ds->get_covariates();
	vector<string>* traits = ds->get_traits();
	map<int, string>* master_map = ds->get_master_map();
	int marker_count = 0;
	//int cov_count = 0;
	//int trait_count = 0;

	map<string, vector<string> > descinfo;
	vector<string> descheaders;
	map<string,int> exclude;
	map<string,int> include;
	map<string,float> frequencies = options.getFrequencies();

	ifstream input;
    input.open(options.getMapFile().c_str(), ios::in);

    if(!input){
		opts::printLog("Error opening map file: " + options.getMapFile() + "\n");
        throw MethodException("Error opening map file: " + options.getMapFile() + "\n");
    }
	int count = 0;
	int num_cols = 3 + options.getMapContainsReferent();

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
        else if((int)elems.size() > num_cols){
			opts::printLog("Map file line has more than " + getString<int>(num_cols) + " elements on line: " + line + "\n");
			throw MethodException("Map file has more than " + getString<int>(num_cols) + " elements on line: " + line + "\n");
        }
		else if(elems.size() < 3 && elems.size() > 0){
			opts::printLog("Map file line has fewer than 3 elements on line: " + line + "\n");
			throw MethodException("Map file has fewer than 3 elements on line: " + line + "\n");
		}

        string chr = elems[0];
        string probe_id = elems[1];
        int bploc = atoi(elems[2].c_str());
        string ref_allele = "";
        if(elems.size() > 3){
          ref_allele = elems[3];
        }

		if(chr == "C"){
			bool cuse = true;
			if(filters != NULL){
				for(int f = 0; f < filters->num_covariate_filters(); f++){
					cuse = filters->run_covariate_filter(f, probe_id);
				}
			}
			if(cuse){
				covariates->push_back(probe_id);
				(*master_map)[count] = chr;
			}
			else{
				(*master_map)[count] = "E";
			}
			count++;
			opts::_COVS_FOUND_++;
			continue;
		}
		else if(chr == "T"){
			bool cuse = true;
			if(filters != NULL){
				for(int f = 0; f < filters->num_trait_filters(); f++){
					cuse = filters->run_trait_filter(f, probe_id);
				}
			}
			if(cuse){
				traits->push_back(probe_id);
				(*master_map)[count] = chr;
			}
			else{
				(*master_map)[count] = "E";
			}
			count++;
			opts::_TRAITS_FOUND_++;
			continue;
		}

		bool use = true;
		opts::_MARKERS_FOUND_++;

        Methods::Marker* m = new Marker(chr, probe_id, bploc);
		if(opts::_AUTOONLY_ && ((m->getChrom() >= opts::_CHRX_) || (m->getChrom() < 1))){
			use = false;
		}
		m->setEnabled(use);
		m->setLoc(marker_count);
		m->setRSID(probe_id);
		m->setReferent(ref_allele);

		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe_id);
			if(found != frequencies.end()){
				m->setMAF(frequencies[probe_id]);
				m->setFreqFlag(true);
			}
		}
		if(opts::_MAPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[probe_id];
			for(unsigned int i = 1; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					m->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					m->assignDetail(descheaders[i], "NA");
				}
			}
		}
        markers->push_back(m);
		(*master_map)[count] = "M";
		marker_count++;
		count++;
    }

    input.clear();
    input.close();

	marker_map->resize(markers->size());

	//put markers in chrom/bploc order
	stable_sort(markers->begin(), markers->end(), less<Methods::Marker*>());

	for(unsigned int i =0; i < markers->size(); i++){
		(*marker_map)[(*markers)[i]->getLoc()] = i;
	}

	//filter marker, set enabled to T/F based on filters.
	if(filters != NULL){
		for(int f = 0; f < filters->num_locus_filters(); f++){
			filters->run_locus_filter(f, markers);
		}
	}


}

/*
 *Function: readMap
 *Description: Takes DataSet as parameter
 *Reads the map file and performs the appropriate inclusion/exclusion of markers
 *
 */
void Helpers::readMapM(DataSet* ds, StepOptions options){
	readMapM(ds, options, NULL);
}

/*
 *Function: readMap
 *Description:
 *Reads the map file and performs the appropriate inclusion/exclusion of markers
 *
 */
void Helpers::readMapM(vector<Methods::Marker*>* markers, vector<int>* marker_map, StepOptions options){
	map<string, vector<string> > descinfo;
	vector<string> descheaders;
	//vector<string> exclude;
	//vector<string> include;
	map<string,int> exclude;
	map<string,int> include;
	map<string,float> frequencies = options.getFrequencies();

	ifstream input;
    input.open(options.getMapFile().c_str(), ios::in);

    if(!input){
		throw MethodException("Error opening map file: " + options.getMapFile() + ".\n");
    }
	int count = 0;
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
			throw MethodException("Map file line has more than 3 elements: " + line + "\n");
        }
		else if(elems.size() < 3 && elems.size() > 0){
			throw MethodException("Map file line has fewer than 3 elements: " + line + "\n");
		}


        string chr = elems[0];
        string probe_id = elems[1];
        int bploc = atoi(elems[2].c_str());

		bool use = true;

		if(exclude.size() > 0){
			map<string, int>::iterator found = exclude.find(probe_id);
			if(found != exclude.end()){
				use = false;
			}
		}
		if(include.size() > 0){
			map<string, int>::iterator found = include.find(probe_id);
			if(found == include.end()){
				use = false;
			}
		}

        Methods::Marker* m = new Marker(chr, probe_id, bploc);
		if(opts::_AUTOONLY_ && ((m->getChrom() >= opts::_CHRX_) || (m->getChrom() < 1))){
			use = false;
		}
		m->setEnabled(use);
		m->setLoc(count);
		m->setRSID(probe_id);


		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe_id);
			if(found != frequencies.end()){
				m->setMAF(frequencies[probe_id]);
				m->setFreqFlag(true);
			}
		}
		if(opts::_MAPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[probe_id];
			for(unsigned int i = 1; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					m->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					m->assignDetail(descheaders[i], "NA");
				}
			}
		}
        markers->push_back(m);
		count++;
    }

    input.clear();
    input.close();

	marker_map->resize(markers->size());

	//put markers in chrom/bploc order
	stable_sort(markers->begin(), markers->end(), less<Methods::Marker*>());

	for(unsigned int i =0; i < markers->size(); i++){
		(*marker_map)[(*markers)[i]->getLoc()] = i;
	}

}


void Helpers::readTPedM(DataSet* ds, StepOptions options, InputFilter* filters){
	vector<Methods::Marker*>* markers = ds->get_markers();
	vector<Methods::Sample*>* samples = ds->get_samples();
	vector<int>* marker_map = ds->get_marker_map();
	vector<string>* covs = ds->get_covariates();
	vector<string>* traits = ds->get_traits();

	map<string, vector<string> > descinfo;
	vector<string> descheaders;
	map<string,int> exclude;
	map<string,int> include;
	map<string,float> frequencies = options.getFrequencies();
	vector<bool> includeme;

	ifstream input;
    input.open(options.getTPedFile().c_str(), ios::in);

    if(!input){
		opts::printLog("Error opening ped file: " + options.getTPedFile() + "\n");
		throw MethodException("Error opening ped file: " + options.getTPedFile() + "\n");
    }
	int count = 0;
    while(!input.eof()){
		string chr;
		string probe_id;
		string cm;
		string bp;
		string referent = "";
		input >> chr;
		if(chr==""){
			continue;
		}
		input >> probe_id >> cm >> bp;
		if(options.getMapContainsReferent()){
			input >> referent;
		}
		input.ignore(1000000,'\n');


		bool use = true;
		if(chr == "C"){
			covs->push_back(probe_id);
			opts::_COVS_FOUND_++;
			includeme.push_back(use);
			count++;
			continue;
		}
		if(chr == "T"){
			traits->push_back(probe_id);
			opts::_TRAITS_FOUND_++;
			includeme.push_back(use);
			count++;
			continue;
		}
		opts::_MARKERS_FOUND_++;

		int bploc = atoi(bp.c_str());
        Methods::Marker* m = new Marker(chr, probe_id, bploc);
		if(opts::_AUTOONLY_ && ((m->getChrom() >= opts::_CHRX_) || (m->getChrom() < 1))){
			use = false;
		}

		includeme.push_back(use);

		m->setEnabled(use);
		m->setLoc(count);
		m->setRSID(probe_id);
		m->setReferent(referent);


		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe_id);
			if(found != frequencies.end()){
				m->setMAF(frequencies[probe_id]);
				m->setFreqFlag(true);
			}
		}
		if(opts::_MAPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[probe_id];
			for(unsigned int i = 1; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					m->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					m->assignDetail(descheaders[i], "NA");
				}
			}
		}
        markers->push_back(m);
		count++;
    }

    input.clear();
    input.close();

	if(filters != NULL){
		for(int f = 0; f < filters->num_locus_filters(); f++){
			filters->run_locus_filter(f, markers);
		}
		for(unsigned int i = 0; i < markers->size(); i++){
			if(!(*markers)[i]->isEnabled()){
				includeme[i] = false;
			}
		}
	}

	marker_map->resize(markers->size());

	//put markers in chrom/bploc order
	stable_sort(markers->begin(), markers->end(), less<Methods::Marker*>());

	for(unsigned int i =0; i < markers->size(); i++){
		(*marker_map)[(*markers)[i]->getLoc()] = i;
	}

	count = 0;
	for(unsigned int i = 0; i < samples->size(); i++){
		(*samples)[i]->resizeAlleles(markers->size());
		(*samples)[i]->resizeCovariates(covs->size());
		(*samples)[i]->resizeTraits(traits->size());
	}

	FILE* PED;
	PED = fopen(options.getTPedFile().c_str(), "r");
	int i =0;

	while(!feof(PED)){
		string dummy;
		string chr = "";
		int f =0;
		if(readString(PED, &chr)) f++;
		dummy = chr;
		if(dummy==""){
			continue;
		}

		if(dummy.substr(0,1) == "#"){
			while(fgetc(PED) != '\n' && !feof(PED)){}
			continue;
		}

		if(!includeme[i]){
			while(fgetc(PED) != '\n' && !feof(PED)){}
			i++;
			continue;
		}

		if(readString(PED,&dummy)) f++; //probe
		if(readString(PED,&dummy)) f++; //cm
		if(readString(PED,&dummy)) f++; //bploc
		if(options.getMapContainsReferent()){
			if(readString(PED, &dummy)) f++; //referent
		}

		int gn = 0;
		int c=0; //ind count
		int cloc = 0;
		int tloc = 0;
		bool linedone = false;
		string fmsg;
		while(!linedone){
			Methods::Sample* samp = (*samples)[c];
			cout << "Working on Sample: " << samp->toString() << endl;
			string one = "";
			string two = "";
			while(1){
				char ch = fgetc(PED);
				if(ch == '/' || ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(PED)){
					if(ch == '\n' || ch == '\r' || feof(PED)){
						linedone = true;
					}
					if(one.length() > 0){
						gn++;
						break;
					}
					if(ch == '\n' || ch == '\r' || feof(PED)){
						break;
					}
				}
				else{
					one += ch;
				}
			}
			if(chr == "C"){
				double value = -1;
			   	try{
					if(one == options.getCovarMissing()){
						value = (double) atof(options.getCovarMissing().c_str());
					}
					else{
						value = (double) atof(one.c_str());
					}
				}catch(...){
					throw MethodException("Cannot convert " + one + " to number on line " + getString<int>(i + 1) + " for covariate <" + (*covs)[cloc] + "> in file: " + options.getTPedFile() + "\n");
				}
				samp->setCovariate(value, cloc);
				cloc++;
				continue;
			}
			else if(chr == "T"){
				double value = -1;
			   	try{
					if(one == options.getTraitMissing()){
						value = (double) atof(options.getTraitMissing().c_str());
					}
					else{
						value = (double) atof(one.c_str());
					}
				}catch(...){
					throw MethodException("Cannot convert " + one + " to number on line " + getString<int>(i + 1) + " for trait <" + (*traits)[tloc] + "> in file: " + options.getTPedFile() + "\n");
				}
				samp->setTrait(value, tloc);
				tloc++;
				continue;
			}

			if(!linedone){
				while(1){
					char ch = fgetc(PED);
					if(ch == '/' || ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(PED)){
						if(ch == '\n' || ch == '\r' || feof(PED)){
							linedone = true;
						}
						if(two.length() > 0){
							gn++;
							break;
						}
						if(ch == '\n' || ch == '\r' || feof(PED)){
							break;
						}
					}
					else{
						two += ch;
					}
				}
			}

			if(linedone && one.length() == 0 && two.length() == 0){
				break;
			}

			if(includeme[i]){
				int k = (*marker_map)[i];
				Methods::Marker* mark = (*markers)[k];

				int oldallelecount = mark->getNumAlleles();
                if(one != opts::_NOCALL_){
					//new
					if(mark->getAlleleLoc(one) < 0){
						mark->addAllele(one);
					}
        	    }

            	if(two != one){
                    if(two != opts::_NOCALL_){
						//new
						if(mark->getAlleleLoc(two) < 0){
							mark->addAllele(two);
						}
	                }
				}

				if(mark->getNumAlleles() <= 2){
					cout << "Adding normal Marker, i= " << getString<int>(i) << endl;
   		            if(one == mark->getAllele1() && two == mark->getAllele1()){
   	         	        samp->addAone(i, false);
   	    	            samp->addAtwo(i, false);
						samp->addAmissing(i, false);
		            }
	                else if(one != opts::_NOCALL_ && two != opts::_NOCALL_ && one != two){
	       	            samp->addAone(i, false);
 		      	        samp->addAtwo(i, true);
	               		samp->addAmissing(i, false);
					}
		            else if(one == mark->getAllele2() && two == mark->getAllele2()){
	                    samp->addAone(i, true);
	       	            samp->addAtwo(i, true);
	           	    	samp->addAmissing(i, false);
					}
	               	else if(one == opts::_NOCALL_ || two == opts::_NOCALL_){
		                samp->addAone(i, true);
	                    samp->addAtwo(i, true);
						samp->addAmissing(i, true);
	       	        }
				}
				else if(opts::_MICROSATS_){
					cout << "adding microsatellite, i= " << getString<int>(i) << endl;
					cout << "ONE: " << one << " TWO: "<< two << endl;
					samp->addMicroSat(i);
					int loc1 = mark->getAlleleLoc(one);
					int loc2 = mark->getAlleleLoc(two);
					samp->addAbone(i, loc1);
					samp->addAbtwo(i, loc2);
					if(oldallelecount <= 2){
						remapSamples(samples, markers, marker_map, i);
					}
				}
				else if(mark->getNumAlleles() > 2 && !opts::_MICROSATS_){
					opts::printLog("More than 2 unique alleles found for map location: " + getString<int>(k) + ", line: " + getString<int>(c + 1) + ".  Microsatellites not specified.\n");
					throw MethodException("More than 2 unique alleles found for map location: " + getString<int>(k) + ", line: " + getString<int>(c + 1) + ".  Microsatellites not specified.\n");
				}
			}
			c++;
			if(c > (int)samples->size()){
				opts::printLog("Problem with line " + getString<int>(i+1) + " in " + options.getTPedFile() + "\n");
				opts::printLog("Expecting 4 + 2 * " + getString<int>(samples->size()) + " = " +
						getString<int>(4+2*samples->size()) + " columns, but found more\n");
				throw MethodException("Problem with line " + getString<int>(i+1) + " in " + options.getTPedFile() + "\nExpecting 4 + 2 * " + getString<int>(samples->size()) + " = " + getString<int>(4+2*samples->size()) + " columns, but found more\n");
			}
		}//line done? next sample

		if(gn != (int)(2 * samples->size())){
			opts::printLog("Problem with line " + getString<int>(i+1) + " in " + options.getTPedFile() + "\n");
			opts::printLog("Expecting 4 + 2 * " + getString<int>(samples->size()) + " = " +
					getString<int>(4+2*samples->size()) + " columns, but found more\n");
			throw MethodException("Problem with line " + getString<int>(i+1) + " in " + options.getTPedFile() + "\nExpecting 4 + 2 * " + getString<int>(samples->size()) + " = " + getString<int>(4+2*samples->size()) + " columns, but found more\n");
		}
		i++;
	}
	fclose(PED);

	reorderAlleles(samples, markers);
}


void Helpers::readTPedM(DataSet* ds, StepOptions options){
	readTPedM(ds, options, NULL);
}


void Helpers::readTFamM(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, StepOptions options, InputFilter* filters){
	map<string, vector<string> > descinfo;
	vector<string> sexclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> finclude;
	vector<string> fexclude;
	vector<string> descheaders;

	opts::printLog("Reading family information from " + options.getTFamFile() + "\n");

	ifstream PED;
	PED.open(options.getTFamFile().c_str());
	if(!PED){
		opts::printLog("Error opening family information file: " + options.getTFamFile() + ".  Exiting!\n");
		throw MethodException("Error opening family information file: " + options.getTFamFile() + ".  Exiting!\n");
	}
	PED.clear();
		string fam = "";
		string ind = "";
		string dad = "";
		string mom = "";
		string sex = "";
		string aff = "";
		string line = "";
	while(getline(PED, line)){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 6){
			opts::printLog("Family information file column size != 6: " + line + "\n");
			throw MethodException("Family information file column size != 6: " + line + "\n");
		}
		fam = tokens[0];
		ind = tokens[1];
		dad = tokens[2];
		mom = tokens[3];
		sex = tokens[4];
		aff = tokens[5];
		Methods::Sample* samp = new Sample();
		samp->setFamID(fam);
		samp->setInd(ind);
		samp->setEnabled(true);
		opts::_SAMPLES_FOUND_++;

		if(filters != NULL){
			for(int f = 0; f < filters->num_sample_filters(); f++){
				bool res = filters->run_sample_filter(f, samp);
				if(!res && !opts::_KEEP_EXC_SAMPLES_){
					samp->setEnabled(false);
				}
				else if(!res && opts::_KEEP_EXC_SAMPLES_){
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
		}
		samp->setDadID(dad);
		samp->setMomID(mom);
		if(sex == "1"){
			samp->setSex(true);
		}
		else if(sex == "2"){
			samp->setSex(false);
		}

		if(aff == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atof(aff.c_str()));
        if(opts::pedinfo.size() > 0){
            map<string, Methods::Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
            if(sfind != opts::pedinfo.end()){
				Methods::Sample* sfound = sfind->second;
                samp->setDadID(sfound->getDadID());
                samp->setMomID(sfound->getMomID());
                samp->setSex(sfound->getSex());
                samp->setPheno(sfound->getPheno());
            }
        }
		if(samp->getPheno() != 0.0f && samp->getPheno() != 1.0f && samp->getPheno() != 2.0f){
			opts::_BINTRAIT_ = false;
		}

		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(unsigned int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
				    samp->assignDetail(descheaders[i], "NA");
				}
			}

		}
        vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Methods::Family* fam = new Family();
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
			opts::_FAMILIES_FOUND_++;
        }

		if(filters != NULL){
			for(int f = 0; f < filters->num_family_filters(); f++){
				bool res = filters->run_family_filter(f, samp->getFamily());
				if(!res){
					samp->setEnabled(false);
					samp->getFamily()->setEnabled(false);
				}
				else{
					samp->setEnabled(true);
					samp->getFamily()->setEnabled(true);
				}
			}
		}
        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
	}

	PED.clear();
	PED.close();

	assignLinks(families);
}

void Helpers::readTFamM(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, StepOptions options){
	readTFamM(samples, families, options, NULL);
}

void Helpers::readTFamM(DataSet* ds, StepOptions options, InputFilter* filters){
	readTFamM(ds->get_samples(), ds->get_families(), options, filters);
}

void Helpers::readTFamM(DataSet* ds, StepOptions options){
	readTFamM(ds->get_samples(), ds->get_families(), options, NULL);
}




/*
 *Function: readBin
 *Paramters: sample vector, family vector, marker vector, marker_map vector
 *Description:
 *Reads set of binary input files (bim, bed, fam)
 *Based on Plink
 *
 */
void Helpers::readBinM(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, vector<Methods::Marker*>* markers, vector<int>* marker_map, StepOptions options, InputFilter* filters){
	map<string, vector<string> > mdescinfo;
	map<string,int> exclude;
	map<string,int> include;
	vector<string> mdescheaders;
	map<string,float> frequencies = options.getFrequencies();

	opts::printLog("Reading map from " + options.getBinInput() + ".bim\n");

	ifstream MAP((options.getBinInput() + ".bim").c_str(), ios::in);
	if(!MAP){
		opts::printLog("Error opening map file: " + options.getBinInput() + ".bim.  Exiting!\n");
		throw MethodException("Error opening map file: " + options.getBinInput() + ".bim.\n");
	}

	int count = 0;

		string chrom = "";
		string probe = "";
		int bploc = 0;
		string a1 = "";
		string a2 = "";
		string centi = "";
		string rsid = "";
		string enzyme = "";
		string line = "";
		string referent = "";
	int num_cols = 6 + options.getMapContainsReferent();
	while(getline(MAP, line)){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if((int)tokens.size() != num_cols){
			opts::printLog(".bim file column size != " + getString<int>(num_cols) + ": " + line + " stopping!!\n");
			throw MethodException(".bim file column size != " + getString<int>(num_cols) + ": " + line + " stopping!!\n");
		}
		chrom = tokens[0];
		probe = tokens[1];
		centi = tokens[2];
		bploc = atoi(tokens[3].c_str());
		a1 = tokens[4];
		a2 = tokens[5];
		if(options.getMapContainsReferent()){
			referent = tokens[6];
		}
		if(rsid == "."){
			rsid = "";
		}
		if(enzyme == "."){
			enzyme = "";
			opts::_ENZYMES_ = false;
		}
		bool use = true;

		Methods::Marker* m = new Marker(chrom, probe, bploc);
		if(opts::_AUTOONLY_ && ((m->getChrom() >= opts::_CHRX_) || (m->getChrom() < 1))){
			use = false;
		}

		m->setEnabled(use);
		m->setLoc(count);
		m->setAllele1(a1);
		m->setAllele2(a2);
		m->setReferent(referent);
		m->setEnzyme(enzyme);
		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe);
			if(found == frequencies.end()){
				found = frequencies.find(m->getRSID());
				if(found != frequencies.end()){
					m->setMAF(frequencies[m->getRSID()]);
					m->setFreqFlag(true);
				}
			}
			else{
				m->setMAF(frequencies[m->getRSID()]);
				m->setFreqFlag(true);
			}
		}
        if(opts::_MAPDESC_.length() > 0 && mdescinfo.size() > 0){
			vector<string> tokens = mdescinfo[probe];
            for(unsigned int i = 1; i < mdescheaders.size(); i++){
	            if(tokens.size() == mdescheaders.size()){
		            m->assignDetail(mdescheaders[i], tokens[i]);
                }
                else{
	                m->assignDetail(mdescheaders[i], "NA");
                }
            }
        }

		markers->push_back(m);

		count++;
	}
	MAP.close();
	opts::_MARKERS_FOUND_ = markers->size();

	if(filters != NULL){
		for(int f = 0; f < filters->num_locus_filters(); f++){
			filters->run_locus_filter(f, markers);
		}
	}


	marker_map->resize(markers->size());
	stable_sort(markers->begin(), markers->end(), less<Methods::Marker*>());

	for(unsigned int i = 0; i < markers->size(); i++){
		(*marker_map)[(*markers)[i]->getLoc()] = i;
	}


	map<string, vector<string> > descinfo;
	vector<string> sexclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> finclude;
	vector<string> fexclude;
	vector<string> descheaders;

	opts::printLog("Reading family information from " + options.getBinInput() + ".fam\n");

	ifstream PED;
	PED.open((options.getBinInput() + ".fam").c_str());
	if(!PED){
		opts::printLog("Error opening family information file: " + options.getBinInput() + ".fam.  Exiting!\n");
		throw MethodException("Error opening family information file: " + options.getBinInput() + ".fam.\n");
	}
	PED.clear();
		string fam = "";
		string ind = "";
		string dad = "";
		string mom = "";
		string sex = "";
		string aff = "";
		line = "";
	while(getline(PED, line)){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 6){
			opts::printLog("Family information file column size != 6: " + line + " Exitting!!\n");
			throw MethodException("Family information file column size != 6: " + line + " Exitting!!\n");
		}
		fam = tokens[0];
		ind = tokens[1];
		dad = tokens[2];
		mom = tokens[3];
		sex = tokens[4];
		aff = tokens[5];
		Methods::Sample* samp = new Sample();
		samp->setFamID(fam);
		samp->setInd(ind);
		samp->setEnabled(true);
		if(filters != NULL){
			for(int f = 0; f < filters->num_sample_filters(); f++){
				bool res = filters->run_sample_filter(f, samp);
				if(!res){
					samp->setEnabled(false);
					if(opts::_KEEP_EXC_SAMPLES_){
						samp->setExcluded(true);
					}
				}
			}
		}
		samp->setDadID(dad);
		samp->setMomID(mom);
		if(sex == "1"){
			samp->setSex(true);
		}
		else if(sex == "2"){
			samp->setSex(false);
		}

		if(aff == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atof(aff.c_str()));
        if(opts::pedinfo.size() > 0){
            map<string, Methods::Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
            if(sfind != opts::pedinfo.end()){
				Methods::Sample* sfound = sfind->second;
                samp->setDadID(sfound->getDadID());
                samp->setMomID(sfound->getMomID());
                samp->setSex(sfound->getSex());
                samp->setPheno(sfound->getPheno());
            }
        }
		if(samp->getPheno() != 0.0f && samp->getPheno() != 1.0f && samp->getPheno() != 2.0f){
			opts::_BINTRAIT_ = false;
		}

		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(unsigned int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
				    samp->assignDetail(descheaders[i], "NA");
				}
			}

		}
        vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Methods::Family* fam = new Family();
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
        }
		if(filters != NULL){
			for(int f = 0; f < filters->num_family_filters(); f++){
				bool res = filters->run_family_filter(f, samp->getFamily());
				if(!res){
					samp->setEnabled(false);
					samp->getFamily()->setEnabled(false);
					if(opts::_KEEP_EXC_SAMPLES_){
						samp->setExcluded(true);
					}
				}
			}
		}

		samp->resizeAlleles(markers->size());
        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
	}

	PED.clear();
	PED.close();

	opts::_FAMILIES_FOUND_ = families->size();
	opts::_SAMPLES_FOUND_ = samples->size();

	bool ind_major = false;
	bool snp_major = false;
	opts::printLog("Reading genotype bitfile from " + options.getBinInput() + ".bed\n");

	ifstream BIT;
	BIT.open((options.getBinInput()+".bed").c_str(), ios::in | ios::binary);
	if(!BIT){
		opts::printLog("Error opening genotype bitfile: " + options.getBinInput() + ".bed.  Exiting!\n");
		throw MethodException("Error opening genotype bitfile: " + options.getBinInput() + ".bed.\n");
	}
	char temp[1];
	//header
	BIT.read(temp, 1);
	bitset<8> tb;
	tb = temp[0];
	if((tb[2] && tb[3] && tb[5] && tb[6]) && !(tb[0] || tb[1] || tb[4] || tb[7])){
		BIT.read(temp, 1);
		tb = temp[0];
		if((tb[0] && tb[1] && tb[3] && tb[4]) && !(tb[2] || tb[5] || tb[6] || tb[7])){
			BIT.read(temp, 1);
			tb = temp[0];
			if(!tb[0]){
			opts::printLog("IND major mode\n");
				ind_major = true;
			}
			else{
			opts::printLog("SNP major mode\n");
				snp_major = true;
			}
		}
		else{
			opts::printLog("Incorrect bit file version (2nd code)!\n");
			throw MethodException("Incorrect bit file version (2nd code)!\n");
		}
	}
	else{
		opts::printLog("Incorrect bit file version (1st code)!\n");
		throw MethodException("Incorrect bit file version (1st code)!\n");
	}

	int ssize = samples->size();
	int msize = markers->size();

	if(ind_major){
		for(int s = 0; s < ssize; s++){
			Methods::Sample* samp = (*samples)[s];
			for(int m = 0; m < msize;){
				char ch[1];
				BIT.read(ch, 1);
				if(!BIT)
				{
					opts::printLog("Problem with the bed file.\n");
					throw MethodException("Problem with the bed file.\n");
				}
				bitset<8> b;
				b = ch[0];
				int c = 0;
				while(c < 7 && m < msize)
				{
					Methods::Marker* mark = (*markers)[m];
					if(mark->isEnabled())
					{
						int mloc = (*markers)[m]->getLoc();
						samp->addAone(mloc, b[c++]);
						samp->addAtwo(mloc, b[c++]);
						if(samp->getAone(mloc) && !samp->getAtwo(mloc))
						{
							samp->addAtwo(mloc, true);
							samp->addAmissing(mloc, true);
						}
					}
					else
					{
						c+=2;
					}
					m++;
				}
			}
		}
	}
	else{//SNP_MAJOR
		for(int m = 0; m < msize; m++){
			Methods::Marker* mark = (*markers)[m];
			int mloc = mark->getLoc();
			for(int s = 0; s < ssize;){
				char ch[1];
				BIT.read(ch, 1);
				if(!BIT){
					opts::printLog("Problem with the bed file.\n");
					throw MethodException("Problem with the bed file.\n");
				}
				bitset<8> b;
				b = ch[0];
				int c = 0;
				while(c < 7 && s < ssize){
					Methods::Sample* samp = (*samples)[s];
					if(samp->isEnabled() || (samp->isExcluded() && opts::_KEEP_EXC_SAMPLES_)){
						samp->addAone(mloc, b[c++]);
						samp->addAtwo(mloc, b[c++]);
						if(samp->getAone(mloc) && !samp->getAtwo(mloc)){
							samp->addAtwo(mloc, true);
							samp->addAmissing(mloc, true);
						}
					}
					else{
						c+=2;
					}
					s++;
				}
			}
		}
	}
	BIT.clear();
	BIT.close();

	assignLinks(families);
	reorderAlleles(samples, markers);
}

void Helpers::readBinM(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, vector<Methods::Marker*>* markers, vector<int>* marker_map, StepOptions options){
	readBinM(samples, families, markers, marker_map, options, NULL);
}

void Helpers::readBinM(DataSet* ds, StepOptions options, InputFilter* filters){
	readBinM(ds->get_samples(), ds->get_families(), ds->get_markers(), ds->get_marker_map(), options, filters);
}

void Helpers::readBinM(DataSet* ds, StepOptions options){
	readBinM(ds, options, NULL);
}

void Helpers::readPedInfo(){

	opts::printLog("Reading Pedigree information from " + opts::_PEDINFO_ + "\n");

	ifstream PED;
	PED.open(opts::_PEDINFO_.c_str());
	if(!PED){
		opts::printLog("Error opening pedigree information file: " + opts::_PEDINFO_ + ".  Exiting!\n");
		throw MethodException("Error opening pedigree information file: " + opts::_PEDINFO_ + ".  Exiting!\n");
	}
	PED.clear();
		string fam = "";
		string ind = "";
		string dad = "";
		string mom = "";
		string sex = "";
		string aff = "";
		string line = "";
	while(getline(PED, line)){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 6){
			opts::printLog("Pedigree information file column size != 6: " + line + " Exiting!!\n");
			throw MethodException("Pedigree information file column size != 6: " + line + " Exiting!!\n");
		}
		fam = tokens[0];
		ind = tokens[1];
		dad = tokens[2];
		mom = tokens[3];
		sex = tokens[4];
		aff = tokens[5];
		Methods::Sample* samp = new Sample();
		samp->setFamID(fam);
		samp->setInd(ind);
		samp->setEnabled(true);
		samp->setDadID(dad);
		samp->setMomID(mom);
		if(sex == "1"){
			samp->setSex(true);
		}
		else if(sex == "2"){
			samp->setSex(false);
		}

		if(aff == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atof(aff.c_str()));

		opts::pedinfo[(samp->getFamID() + "#" + samp->getInd())] = samp;
	}

	PED.clear();
	PED.close();

}



 /*
 *Function: readPed
 *Description:
 *Reads ped file and stores data into sample, family, and marker vectors
 *
 */
void Helpers::readPedM_3vec(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, vector<Methods::Marker*>* markers, vector<int>* marker_map, StepOptions options){
	map<string, vector<string> > descinfo;
	vector<string> exclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> fexclude;
	vector<string> finclude;
	vector<string> descheaders;
    FILE* input;
    input = fopen(options.getPedFile().c_str(), "r");
	if(!input){
		throw MethodException("Error opening pedfile: " + options.getPedFile() + ".");
	}
	int onind = -1;
    while(!feof(input)){
		onind++;
        Methods::Sample* samp = new Sample();
        int f = 0;
        string temp = "";
        if(readString(input, &temp)){
           samp->setFamID(temp);
              f++;
             temp = "";
         }

		string ftemp = samp->getFamID();
		if(samp->getFamID() == ""){
			delete(samp);
			continue;
		}
        if(ftemp.at(0) == '#'){
			delete(samp);
			while(fgetc(input) != '\n' && !feof(input)){}
            continue;
        }
        string sex = "";
        string pheno = "";
        if(readString(input, &temp)){
            samp->setInd(temp);
            f++;
            temp = "";
        }
		samp->setEnabled(true);

		if(exclude.size() > 0){
			vector<string>::iterator found = find(exclude.begin(), exclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found != exclude.end()){
				if(!opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					continue;
				}
				else{
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
		}
		if(sinclude.size() > 0){
			vector<string>::iterator found = find(sinclude.begin(), sinclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found == sinclude.end()){
				if(!opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					continue;
				}
				else{
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
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

		if(pheno == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atof(pheno.c_str()));
        if(opts::pedinfo.size() > 0){
            map<string, Methods::Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
            if(sfind != opts::pedinfo.end()){
				Methods::Sample* sfound = sfind->second;
                samp->setDadID(sfound->getDadID());
                samp->setMomID(sfound->getMomID());
                samp->setSex(sfound->getSex());
                samp->setPheno(sfound->getPheno());
            }
        }
		if(samp->getPheno() != 0.0f && samp->getPheno() != 1.0f && samp->getPheno() != 2.0f){
			opts::_BINTRAIT_ = false;
		}

		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(unsigned int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					samp->assignDetail(descheaders[i], "NA");
				}
			}

		}
        samp->resizeAlleles(markers->size());

        int gn = 0;
        int i = 0;
        bool linedone = false;

        string fmsg;
        while(!linedone){
            string one = "";
            string two = "";

            while(1){
                char ch = fgetc(input);

                if(ch == '/' || ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input)){
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
                    if(ch == '/' || ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input)){
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

				if(i > (int)markers->size()){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>((2 * markers->size()) + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
				Methods::Marker* m = (*markers)[(*marker_map)[i]];
				if(m->isEnabled()){
					int oldallelecount = m->getNumAlleles();
	                if(one != "0"){
						//new
						if(m->getAlleleLoc(one) < 0){
							m->addAllele(one);
						}

        	        }

            	    if(two != one){
                	    if(two != "0"){
							//new
							if(m->getAlleleLoc(two) < 0){
								m->addAllele(two);
							}

                    	}
	                }

          // smd -- changes here for reading in and keeping genotype
          // information as 3 bit vectors -- doesn't keep phased data yet
					if(m->getNumAlleles() <= 2){
   		 	            if(one == m->getAllele1() && two == m->getAllele1()){
   	    	     	        samp->addAone(i, false);
   	     		            samp->addAtwo(i, false);
   	     		            samp->addAmissing(i, false);
		                }
	    	            else if(one != "0" && two != "0" && one != two){
	        	            samp->addAone(i, false);
 		           	        samp->addAtwo(i, true);
 		           	        samp->addAmissing(i, false);
	                	}
		                else if(one == m->getAllele2() && two == m->getAllele2()){
	    	                samp->addAone(i, true);
	        	            samp->addAtwo(i, true);
	        	            samp->addAmissing(i, false);
	            	    }
	                	else if(one == "0" || two == "0"){
		                    samp->addAone(i, true);
	    	                samp->addAtwo(i, true);
	    	                samp->addAmissing(i, true);
	        	        }
					}
					else if(opts::_MICROSATS_){
						samp->addMicroSat(i);
						int loc1 = m->getAlleleLoc(one);
						int loc2 = m->getAlleleLoc(two);

						samp->addAbone(i, loc1);
						samp->addAbtwo(i, loc2);
						if(oldallelecount <= 2){
							remapSamples(samples, markers, marker_map, i);
						}
					}
					else if(m->getNumAlleles() > 2 && !opts::_MICROSATS_){
						throw MethodException("More than 2 unique alleles found for map location: " + getString<int>(i) + ", line: " + getString<int>(onind + 1) + ".  Microsatellites not specified.\n");
					}
				}
				else{
					samp->addAone(i,true);
					samp->addAtwo(i, true);
					samp->addAmissing(i, true);
				}
    	    	i++;
				if(i > (int)markers->size()){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>((2 * markers->size()) + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
            }/*end !linedone*/
        }/*end while(1)*/
		if(gn != (int)(2* markers->size())){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>(((2 * markers->size()) + 6));
					text += " columns but found ";
					text += getString<int>((f + gn));
					text += "\n";
					throw MethodException(text);
		}

        vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Methods::Family* fam = new Family();
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
        }
		if(fexclude.size() > 0){
			vector<string>::iterator found = find(fexclude.begin(), fexclude.end(), samp->getFamID());
			if(found != fexclude.end()){
				samp->setEnabled(false);
				vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		if(finclude.size() > 0){
			vector<string>::iterator found = find(finclude.begin(), finclude.end(), samp->getFamID());
			if(found == finclude.end()){
				samp->setEnabled(false);
				vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
    }/*end while(eof)*/

    fclose(input);
}


/*
 * Function: readMapMdr
 * Reads a map file intended for MDR input and filters based on vector of filter functions passed.  No filtering
 * if vector is null or empty
 * Filter function needs to have as parameters for Markers: Marker and vector of markers to compare against.
 * Filter for covs & traits takes string name and vector of traits/covs to compare against.
 *
 * Return: bool
 */

void Helpers::readMapMdr(DataSet* ds, StepOptions options,  InputFilter* filters){
	vector<Methods::Marker*>* markers = ds->get_markers();
	vector<int>* marker_map = ds->get_marker_map();
	vector<string>* covariates = ds->get_covariates();
	vector<string>* traits = ds->get_traits();
	map<int, string>* master_map = ds->get_master_map();
	int marker_count = 0;

	map<string, vector<string> > descinfo;
	vector<string> descheaders;
	map<string,int> exclude;
	map<string,int> include;
	map<string,float> frequencies = options.getFrequencies();

	ifstream input;
    input.open(options.getMdrMapFile().c_str(), ios::in);

    if(!input){
		opts::printLog("Error opening map file: " + options.getMdrMapFile() + "\n");
        throw MethodException("Error opening map file: " + options.getMdrMapFile() + "\n");
    }
	int count = 0;
	int num_cols = 3 + options.getMapContainsReferent() + 2;  //+2 for alleles minor, then major

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

        bool no_alleles = false;
        bool referent = options.getMapContainsReferent();
        if(referent && (int)elems.size() == 6){
        	no_alleles = false;
        }
        else if(referent && (int)elems.size() == (num_cols - 2)){
        	no_alleles = true;
        }
        else if(!referent && (int)elems.size() == 3){
        	no_alleles = true;
        }
        else if(!referent && (int)elems.size() == 5){
        	no_alleles = false;
        }


        if(elems.size() == 0){
            continue;
        }
        if((int)elems.size() > num_cols){
			opts::printLog("Map file line has more than " + getString<int>(num_cols) + " elements on line: " + line + "\n");
			throw MethodException("Map file has more than " + getString<int>(num_cols) + " elements on line: " + line + "\n");
        }
		else if(elems.size() < 3 && elems.size() > 0){
			opts::printLog("Map file line has fewer than 3 elements on line: " + line + "\n");
			throw MethodException("Map file has fewer than 3 elements on line: " + line + "\n");
		}

        string chr = elems[0];
        string probe_id = elems[1];
        int bploc = atoi(elems[2].c_str());
        string ref_allele = "";
        string allele1 = "A";
        string allele2 = "B";
        if(elems.size() > 3 && elems.size() < 5){
          ref_allele = elems[3];
        }
        if(!no_alleles && referent && elems.size() == 6){
        	allele1 = elems[4];
        	allele2 = elems[5];
        }
        else if(!no_alleles && !referent && elems.size() == 5){
        	allele1 = elems[3];
        	allele2 = elems[4];
        }

		if(chr == "C"){
			bool cuse = true;
			if(filters != NULL){
				for(int f = 0; f < filters->num_covariate_filters(); f++){
					cuse = filters->run_covariate_filter(f, probe_id);
				}
			}
			if(cuse){
				covariates->push_back(probe_id);
				(*master_map)[count] = chr;
			}
			else{
				(*master_map)[count] = "E";
			}
			count++;
			opts::_COVS_FOUND_++;
			continue;
		}
		else if(chr == "T"){
			bool cuse = true;
			if(filters != NULL){
				for(int f = 0; f < filters->num_trait_filters(); f++){
					cuse = filters->run_trait_filter(f, probe_id);
				}
			}
			if(cuse){
				traits->push_back(probe_id);
				(*master_map)[count] = chr;
			}
			else{
				(*master_map)[count] = "E";
			}
			count++;
			opts::_TRAITS_FOUND_++;
			continue;
		}

		bool use = true;
		opts::_MARKERS_FOUND_++;

        Methods::Marker* m = new Marker(chr, probe_id, bploc);
		if(opts::_AUTOONLY_ && ((m->getChrom() >= opts::_CHRX_) || (m->getChrom() < 1))){
			use = false;
		}
		m->setEnabled(use);
		m->setLoc(marker_count);
		m->setRSID(probe_id);
		m->setReferent(ref_allele);
		if(!no_alleles){
			m->setAllele1(allele1);
			m->setAllele2(allele2);
		}

		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe_id);
			if(found != frequencies.end()){
				m->setMAF(frequencies[probe_id]);
				m->setFreqFlag(true);
			}
		}
		if(opts::_MAPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[probe_id];
			for(unsigned int i = 1; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					m->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					m->assignDetail(descheaders[i], "NA");
				}
			}
		}
        markers->push_back(m);
		(*master_map)[count] = "M";
		marker_count++;
		count++;
    }

    input.clear();
    input.close();

	marker_map->resize(markers->size());

	//put markers in chrom/bploc order
	stable_sort(markers->begin(), markers->end(), less<Methods::Marker*>());

	for(unsigned int i =0; i < markers->size(); i++){
		(*marker_map)[(*markers)[i]->getLoc()] = i;
	}

	//filter marker, set enabled to T/F based on filters.
	if(filters != NULL){
		for(int f = 0; f < filters->num_locus_filters(); f++){
			filters->run_locus_filter(f, markers);
		}
	}


}



/*
 *Function: readMdr (MDR Style input)
 *Description:
 *Reads ped file and stores data into sample, family, and marker vectors
 *
 */
void Helpers::readMdr(DataSet* set, StepOptions options, InputFilter* filters){
	bool reset_todigit = false;
	if(opts::_TODIGIT_){
		reset_todigit = true;
		opts::_TODIGIT_ = false;
	}
	if(set == NULL){
		throw MethodException("NULL DataSet structure passed to readPedM_3vec_set...\n");
	}
	map<string, vector<string> > descinfo;
	vector<string> exclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> fexclude;
	vector<string> finclude;
	vector<string> descheaders;


  vector<Methods::Sample*>* samples = set->get_samples();
  vector<Methods::Marker*>* markers = set->get_markers();
  vector<Methods::Family*>* families = set->get_families();
  vector<string>* covariates = set->get_covariates();
  vector<string>* traits = set->get_traits();
  vector<int>* marker_map = set->get_marker_map();
	map<int, string>* master_map = set->get_master_map();

	int columns = markers->size() * 1 + (master_map->size() - markers->size());

  bool missingData = false;

    FILE* input;
    input = fopen(options.getMdrPedFile().c_str(), "r");
	if(!input){
		throw MethodException("Error opening pedfile: " + options.getMdrPedFile() + ".");
	}
	int onind = -1;
    while(!feof(input)){
		onind++;
        Methods::Sample* samp = new Sample();
        int f = 0;
        string temp = "";
        if(readString(input, &temp)){
           samp->setFamID(temp);
              f++;
             temp = "";
         }

		string ftemp = samp->getFamID();
		if(samp->getFamID() == ""){
			delete(samp);
			continue;
		}
        if(ftemp.at(0) == '#'){
			delete(samp);
			while(fgetc(input) != '\n' && !feof(input)){}
            continue;
        }
		opts::_SAMPLES_FOUND_++;
        /*check for comments?*/
        string sex = "";
        string pheno = "";
        if(readString(input, &temp)){
            samp->setInd(temp);
            f++;
            temp = "";
        }
		samp->setEnabled(true);

		if(filters != NULL){
			bool deleted = false;
			for(int f = 0; f < filters->num_sample_filters(); f++){
				bool res = filters->run_sample_filter(f, samp);
				if(!res && !opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					deleted = true;
					break;
				}
				else if(!res && opts::_KEEP_EXC_SAMPLES_){
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
			if(deleted){
				continue;
			}
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

		if(pheno == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atof(pheno.c_str()));
        if(opts::pedinfo.size() > 0){
            map<string, Methods::Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
            if(sfind != opts::pedinfo.end()){
				Methods::Sample* sfound = sfind->second;
                samp->setDadID(sfound->getDadID());
                samp->setMomID(sfound->getMomID());
                samp->setSex(sfound->getSex());
                samp->setPheno(sfound->getPheno());
            }
        }
		if(samp->getPheno() != 0.0f && samp->getPheno() != 1.0f && samp->getPheno() != 2.0f){
			opts::_BINTRAIT_ = false;
		}

		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(unsigned int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					samp->assignDetail(descheaders[i], "NA");
				}
			}

		}
        samp->resizeAlleles(markers->size());
		samp->resizeCovariates(covariates->size());
		samp->resizeTraits(traits->size());

        int gn = 0;
        int i = 0;
		int overall = 0;
		int cloc = 0;
		int tloc = 0;
        bool linedone = false;

        string fmsg;
        while(!linedone){
            string one = "";

            while(1){
                char ch = fgetc(input);

                if(ch == '/' || ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input)){
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
			if((*master_map)[overall] == "C"){
				double value = -1;
			   	try{
					if(one == options.getCovarMissing()){
						try{
							for(unsigned int c = 0; c < one.size(); c++){
								if(!isdigit(one[c]) && one[c] != '.' && one[c] != '-'){
									throw "oops!";
								}
							}
							value = (double) atof(options.getCovarMissing().c_str());
						}
						catch(...){
							value = options.getDefaultCovarMissing();
						}
					}
					else{
						value = (double) atof(one.c_str());
					}
				}catch(...){
					throw MethodException("Cannot convert " + one + " to number on line " + getString<int>(onind + 1) + " for covariate <" + (*covariates)[cloc] + "> in file: " + options.getPedFile() + "\n");
				}
				samp->setCovariate(value, cloc);
				overall++;
				cloc++;
				continue;
			}
			else if((*master_map)[overall] == "T"){
				double value = -1;
			   	try{
					if(one == options.getTraitMissing()){
						try{
							for(unsigned int c = 0; c < one.size(); c++){
								if(!isdigit(one[c]) && one[c] != '.' && one[c] != '-'){
									throw "oops!";
								}
							}
							value = (double) atof(options.getTraitMissing().c_str());
						}catch(...){
							value = options.getDefaultTraitMissing();
						}
					}
					else{
						value = (double) atof(one.c_str());
					}
				}catch(...){
					throw MethodException("Cannot convert " + one + " to number on line " + getString<int>(onind + 1) + " for trait <" + (*traits)[tloc] + "> in file: " + options.getPedFile() + "\n");
				}
				samp->setTrait(value, tloc);
				overall++;
				tloc++;
				continue;
			}
			//covariate or trait exclude
			else if((*master_map)[overall] == "E"){
				overall++;
				continue;
			}
				if(i > columns){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>(columns + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
				Methods::Marker* m = (*markers)[(*marker_map)[i]];
				if(m->isEnabled()){
	                if(one != opts::_NOCALL_){
	                	if(m->getAllele1() == "" && m->getAllele2() == ""){
	                		m->setAllele1("A");
	                		m->setAllele2("B");
	                	}
	                	//new
        	        }
          // smd -- changes here for reading in and keeping genotype
          // information as 3 bit vectors -- doesn't keep phased data yet
					if(m->getNumAlleles() <= 2){
						if(one == "0"){
   	    	     	        samp->addAone(i, false);
   	     		            samp->addAtwo(i, false);
   	     		            samp->addAmissing(i, false);
		                }
						else if(one == "1"){
	        	            samp->addAone(i, false);
 		           	        samp->addAtwo(i, true);
 		           	        samp->addAmissing(i, false);
	                	}
						else if(one == "2"){
	    	                samp->addAone(i, true);
	        	            samp->addAtwo(i, true);
	        	            samp->addAmissing(i, false);
	            	    }
	                	else if(one == opts::_NOCALLMDR_){
		                    samp->addAone(i, true);
	    	                samp->addAtwo(i, true);
	    	                samp->addAmissing(i, true);
	    	                missingData = true;
	        	        }
					}
				}
				else{
					samp->addAone(i,true);
					samp->addAtwo(i, true);
					samp->addAmissing(i, true);
				}
				overall++;
    	    	i++;
				if(i > columns){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>(columns + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
				if(overall > columns){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>(columns + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
        }/*end while(1)*/
		if(gn != columns){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>(columns + 6);
					text += " columns but found ";
					text += getString<int>((f + gn));
					text += "\n";
					throw MethodException(text);
		}

        vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Methods::Family* fam = new Family();
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
			opts::_FAMILIES_FOUND_++;
        }

		if(filters != NULL){
			bool deleted = false;
			for(int f = 0; f < filters->num_family_filters(); f++){
				bool res = filters->run_family_filter(f, samp->getFamily());
				if(!res && !opts::_KEEP_EXC_SAMPLES_){
					families->erase(families->end() - 1);
					delete(samp->getFamily());
					delete(samp);
					deleted = true;
					break;
				}
				else if(!res && opts::_KEEP_EXC_SAMPLES_){
					samp->getFamily()->setEnabled(false);
					samp->getFamily()->setExcluded(true);
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
			if(deleted){
				continue;
			}
		}

        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
    }/*end while(eof)*/

    fclose(input);

    assignLinks(families);
	reorderAlleles(samples, markers);
    // set DataSet conditions
    set->missing_data_present(missingData);
    set->set_missing_value(3);
    set->set_max_locus(2);
    set->set_max_allele(1);
    set->set_affection_vectors();
	if(reset_todigit){
		opts::_TODIGIT_ = true;
	}

}


double Helpers::ltqnorm(double p){
  double q, r;


  if (p < 0 || p > 1){
	return 0.0;
  }
  else if (p == 0){
	return -HUGE_VAL /* minus "infinity" */;
  }
  else if (p == 1){
	return HUGE_VAL /* "infinity" */;
  }
  else if (p < LOW){
	/* Rational approximation for lower region */
	q = sqrt(-2*log(p));
	return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
    ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
  else if (p > HIGH){
	/* Rational approximation for upper region */
	q  = sqrt(-2*log(1-p));
	return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
		    ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
  else{
    /* Rational approximation for central region */
    q = p - 0.5;
    r = q*q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
	    (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  }
}

//new 12-06-2010
//This method is part of the process of taking input in the form of long format filesets
//These filesets are comprised of a Map file, a Fam/Ped file, and an Lgen file
//
void Helpers::readLgenFile(DataSet* set, StepOptions options, InputFilter* filters)
{
	//new 12-09-2010
	//Need to create a hash of Markers
	map<string, Methods::Marker*> markers_map;
	vector<Methods::Marker*>* markers = set->get_markers();
	for(unsigned int i = 0; i < markers->size(); i++)
	{
		markers_map[(*markers)[i]->getRSID()] = (*markers)[i];
	}

	//Create a hash of Samples
	vector<Methods::Sample*>* samples = set->get_samples();
	map<string, Sample*> samples_map;
	for (unsigned int i = 0; i < samples->size(); i++)
	{
		samples_map[(*samples)[i]->getFamID() + "#" + (*samples)[i]->getIndOrig()] = (*samples)[i];
	}

	//process each line of the file
	FILE* input;
	input = fopen(opts::_LGENFILE_.c_str(), "r");
	if(!input){
		throw MethodException("Error opening lgen file: " + opts::_LGENFILE_ + ".");
		return;
	}

	bool reset_todigit = false;
	if(opts::_TODIGIT_){
		reset_todigit = true;
		opts::_TODIGIT_ = false;
	}

	vector<int>* marker_map = set->get_marker_map();

	bool missingData = false;

	int onind = -1;
	while(!feof(input))
	{
		onind++;

	//find the current line's Sample
	string temp = "";
	readString(input, &temp);

	string familyID = temp;
	if(familyID == "")
	{
		continue;
	}
	if (familyID.at(0) == '#')
	{
		while(fgetc(input) != '\n' && !feof(input)){}
		continue;
	}

	temp = "";
	readString(input, &temp);

	string sampleID = temp;
	if(sampleID == "")
	{
		continue;
	}

	vector<Methods::Sample* >::iterator sampiter;
	vector<Methods::Marker*>::iterator markiter;
	Methods::Sample* mysamp = new Methods::Sample();
	Methods::Marker* mymark = new Methods::Marker();

	mysamp = samples_map[familyID + "#" + sampleID];

	//find the current line's Marker
	temp = "";
	readString(input, &temp);
	string markerID = temp;
	if(markerID =="")
	{
		continue;
	}

	mymark = markers_map[markerID];
	int mloc = mymark->getLoc();
	string one = "";
	string two = "";

	//get the genotype strings (one and two) from the current line
	//check to see if the user specified -compound-genotypes (no space between the alleles)...
	if(opts::_COMPOUND_GENOTYPES_)
	{
		//the alleles do not have a space between them, need to seperate them...
		temp = "";
		readString(input, &temp);
		string genotype = temp;
		if (genotype =="")
			continue;
		if(genotype.length() != 2)
			continue;

		one = genotype[0];
		two = genotype[1];
	}
	else //the alleles are seperated by a space
	{
		temp = "";
		readString(input, &temp);
		one = temp;
		if(one == "")
		{
			continue;
		}

		temp = "";
		readString(input, &temp);
		two = temp;
		if(two == "")
		{
			continue;
		}
	}
	//add to the Sample's genotype vectors
	if(mymark->isEnabled())
	{
		int oldallelecount = mymark->getNumAlleles();
		if(one != opts::_NOCALL_)
		{
			//new
			if(mymark->getAlleleLoc(one) < 0)
			{
				mymark->addAllele(one);
			}

		}

		if(two != one)
		{
			if(two != opts::_NOCALL_)
			{
				//new
				if(mymark->getAlleleLoc(two) < 0)
				{
					mymark->addAllele(two);
				}

			}
		}
		// smd -- changes here for reading in and keeping genotype
		// information as 3 bit vectors -- doesn't keep phased data yet
		if(mymark->getNumAlleles() <= 2)
		{
			if(one == mymark->getAllele1() && two == mymark->getAllele1())
			{
				mysamp->addAone(mloc, false);
				mysamp->addAtwo(mloc, false);
				mysamp->addAmissing(mloc, false);
			}
			else if(one != opts::_NOCALL_ && two != opts::_NOCALL_ && one != two)
			{
				mysamp->addAone(mloc, false);
				mysamp->addAtwo(mloc, true);
				mysamp->addAmissing(mloc, false);
			}
			else if(one == mymark->getAllele2() && two == mymark->getAllele2())
			{
				mysamp->addAone(mloc, true);
				mysamp->addAtwo(mloc, true);
				mysamp->addAmissing(mloc, false);
			}
			else if(one == opts::_NOCALL_ || two == opts::_NOCALL_)
			{
				mysamp->addAone(mloc, true);
				mysamp->addAtwo(mloc, true);
				mysamp->addAmissing(mloc, true);
				missingData = true;
			}
		}
		else if(opts::_MICROSATS_)
		{
			mysamp->addMicroSat(mloc);
			int loc1 = mymark->getAlleleLoc(one);
			int loc2 = mymark->getAlleleLoc(two);

			mysamp->addAbone(mloc, loc1);
			mysamp->addAbtwo(mloc, loc2);
			if(oldallelecount <= 2)
			{
				remapSamples(samples, markers, marker_map, mloc);
			}
		}
		else if(mymark->getNumAlleles() > 2 && !opts::_MICROSATS_)
		{
			throw MethodException("More than 2 unique alleles found for map location: " + getString<int>(mloc) + ", line: " + getString<int>(onind + 1) + ".  Microsatellites not specified.\n");
		}
	}
	else
	{
		mysamp->addAone(mloc,true);
		mysamp->addAtwo(mloc, true);
		mysamp->addAmissing(mloc, true);
	}

	}//end while(!feof(input))
	reorderAlleles(samples, markers);
	// set DataSet conditions
    set->missing_data_present(missingData);
    set->set_missing_value(3);
    set->set_max_locus(2);
    set->set_max_allele(1);
    set->set_affection_vectors();
	if(reset_todigit){
		opts::_TODIGIT_ = true;
	}
}// end Helpers::readLgenFile

//new 12-10-2010
//This method reads in a file containing reference alleles for use with the -lgen-file option
//This will override the referent allele, if any, included as part of a map file
void Helpers::readReferenceFile(DataSet* set, StepOptions options, InputFilter* filters)
{
	//process each line of the file
	FILE* input;
	input = fopen(opts::_REFERENCE_FILE_.c_str(), "r");
	if(!input){
		throw MethodException("Error opening reference file: " + opts::_REFERENCE_FILE_ + ".");
		return;
	}

	//build a map of the markers keyed on their rsid for quick access when reading the reference allele file
	map<string, Methods::Marker*> markers_map;
	vector<Methods::Marker*>* markers = set->get_markers();
	for(int i = 0; i < (int)markers->size(); i++)
	{
		markers_map[(*markers)[i]->getRSID()] = (*markers)[i];
	}

	string temp = "";
	while(!feof(input))
	{
		temp = "";
		readString(input, &temp);

		string rsid = temp;
		if(rsid == "")
		{
			continue;
		}
		if (rsid.at(0) == '#')
		{
			while(fgetc(input) != '\n' && !feof(input)){}
			continue;
		}

		temp = "";
		readString(input, &temp);

		string refAllele = temp;
		if(refAllele == "")
		{
			continue;
		}

		Methods::Marker* myMark = markers_map[rsid];

		myMark->setReferent(refAllele);
	}
}

void Helpers::readPedM_3vec_set(DataSet* set, StepOptions options, InputFilter* filters){
	bool reset_todigit = false;
	if(opts::_TODIGIT_){
		reset_todigit = true;
		opts::_TODIGIT_ = false;
	}
	if(set == NULL){
		throw MethodException("NULL DataSet structure passed to readPedM_3vec_set...\n");
	}
	map<string, vector<string> > descinfo;
	vector<string> exclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> fexclude;
	vector<string> finclude;
	vector<string> descheaders;
	cout << "MAP: " << options.getMapFile() << endl;
	cout << "PED: " << options.getPedFile() << endl;


  vector<Methods::Sample*>* samples = set->get_samples();
  vector<Methods::Marker*>* markers = set->get_markers();
  vector<Methods::Family*>* families = set->get_families();
  vector<string>* covariates = set->get_covariates();
  vector<string>* traits = set->get_traits();
  vector<int>* marker_map = set->get_marker_map();
	map<int, string>* master_map = set->get_master_map();

	int columns = markers->size() * 2 + (master_map->size() - markers->size());

  bool missingData = false;

    FILE* input;
    input = fopen(options.getPedFile().c_str(), "r");
	if(!input){
		throw MethodException("Error opening pedfile: " + options.getPedFile() + ".");
	}
	int onind = -1;
    while(!feof(input)){
		onind++;
        Methods::Sample* samp = new Sample();
        int f = 0;
        string temp = "";
        if(readString(input, &temp)){
           samp->setFamID(temp);
              f++;
             temp = "";
         }

		string ftemp = samp->getFamID();
		if(samp->getFamID() == ""){
			delete(samp);
			continue;
		}
        if(ftemp.at(0) == '#'){
			delete(samp);
			while(fgetc(input) != '\n' && !feof(input)){}
            continue;
        }
		opts::_SAMPLES_FOUND_++;
        /*check for comments?*/
        string sex = "";
        string pheno = "";
        if(readString(input, &temp)){
            samp->setInd(temp);
            f++;
            temp = "";
        }
		samp->setEnabled(true);

		if(filters != NULL){
			bool deleted = false;
			for(int f = 0; f < filters->num_sample_filters(); f++){
				bool res = filters->run_sample_filter(f, samp);
				if(!res && !opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					deleted = true;
					break;
				}
				else if(!res && opts::_KEEP_EXC_SAMPLES_){
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
			if(deleted){
				continue;
			}
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

		if(pheno == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atof(pheno.c_str()));
        if(opts::pedinfo.size() > 0){
            map<string, Methods::Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
            if(sfind != opts::pedinfo.end()){
				Methods::Sample* sfound = sfind->second;
                samp->setDadID(sfound->getDadID());
                samp->setMomID(sfound->getMomID());
                samp->setSex(sfound->getSex());
                samp->setPheno(sfound->getPheno());
            }
        }
		if(samp->getPheno() != 0.0f && samp->getPheno() != 1.0f && samp->getPheno() != 2.0f){
			opts::_BINTRAIT_ = false;
		}

		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(unsigned int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					samp->assignDetail(descheaders[i], "NA");
				}
			}

		}
        samp->resizeAlleles(markers->size());
		samp->resizeCovariates(covariates->size());
		samp->resizeTraits(traits->size());

        int gn = 0;
        int i = 0;
		int overall = 0;
		int cloc = 0;
		int tloc = 0;
        bool linedone = false;

        //new 12-07-2010
        if(! opts::_LGENFILE_.length() > 0) //if there is an lgen file, there will not be any genotype data in the ped file
        {

        string fmsg;
        while(!linedone)
        {
            string one = "";
            string two = "";

            while(1){
                char ch = fgetc(input);

                if(ch == '/' || ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input)){
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
			if((*master_map)[overall] == "C"){
				double value = -1;
			   	try{
					if(one == options.getCovarMissing()){
						try{
							for(unsigned int c = 0; c < one.size(); c++){
								if(!isdigit(one[c]) && one[c] != '.' && one[c] != '-'){
									throw "oops!";
								}
							}
							value = (double) atof(options.getCovarMissing().c_str());
						}
						catch(...){
							value = options.getDefaultCovarMissing();
						}
					}
					else{
						value = (double) atof(one.c_str());
					}
				}catch(...){
					throw MethodException("Cannot convert " + one + " to number on line " + getString<int>(onind + 1) + " for covariate <" + (*covariates)[cloc] + "> in file: " + options.getPedFile() + "\n");
				}
				samp->setCovariate(value, cloc);
				overall++;
				cloc++;
				continue;
			}
			else if((*master_map)[overall] == "T"){
				double value = -1;
			   	try{
					if(one == options.getTraitMissing()){
						try{
							for(unsigned int c = 0; c < one.size(); c++){
								if(!isdigit(one[c]) && one[c] != '.' && one[c] != '-'){
									throw "oops!";
								}
							}
							value = (double) atof(options.getTraitMissing().c_str());
						}catch(...){
							value = options.getDefaultTraitMissing();
						}
					}
					else{
						value = (double) atof(one.c_str());
					}
				}catch(...){
					throw MethodException("Cannot convert " + one + " to number on line " + getString<int>(onind + 1) + " for trait <" + (*traits)[tloc] + "> in file: " + options.getPedFile() + "\n");
				}
				samp->setTrait(value, tloc);
				overall++;
				tloc++;
				continue;
			}
			//covariate or trait exclude
			else if((*master_map)[overall] == "E"){
				overall++;
				continue;
			}
            if(!linedone){
                while(1){
                    char ch = fgetc(input);
                    if(ch == '/' || ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input)){
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

				if(i > columns){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + options.getPedFile() + "\n";
					text += "Expecting ";
					text += getString<int>(columns + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
				Methods::Marker* m = (*markers)[(*marker_map)[i]];
				if(m->isEnabled()){
					int oldallelecount = m->getNumAlleles();
	                if(one != opts::_NOCALL_){
						//new
						if(m->getAlleleLoc(one) < 0){
							m->addAllele(one);
						}

        	        }

            	    if(two != one){
                	    if(two != opts::_NOCALL_){
							//new
							if(m->getAlleleLoc(two) < 0){
								m->addAllele(two);
							}
                    	}
	                }

				  // smd -- changes here for reading in and keeping genotype
				  // information as 3 bit vectors -- doesn't keep phased data yet
					if(m->getNumAlleles() <= 2)
					{
   		 	            if(one == m->getAllele1() && two == m->getAllele1())
   		 	            {
   	    	     	        samp->addAone(i, false);
   	     		            samp->addAtwo(i, false);
   	     		            samp->addAmissing(i, false);
		                }
	    	            else if(one != opts::_NOCALL_ && two != opts::_NOCALL_ && one != two)
	    	            {
	        	            samp->addAone(i, false);
 		           	        samp->addAtwo(i, true);
 		           	        samp->addAmissing(i, false);
	                	}
		                else if(one == m->getAllele2() && two == m->getAllele2())
		                {
	    	                samp->addAone(i, true);
	        	            samp->addAtwo(i, true);
	        	            samp->addAmissing(i, false);
	            	    }
	                	else if(one == opts::_NOCALL_ || two == opts::_NOCALL_)
	                	{
		                    samp->addAone(i, true);
	    	                samp->addAtwo(i, true);
	    	                samp->addAmissing(i, true);
	    	                missingData = true;
	        	        }
					}
					else if(opts::_MICROSATS_)
					{
						samp->addMicroSat(i);
						int loc1 = m->getAlleleLoc(one);
						int loc2 = m->getAlleleLoc(two);

						samp->addAbone(i, loc1);
						samp->addAbtwo(i, loc2);
						if(oldallelecount <= 2)
						{
							remapSamples(samples, markers, marker_map, i);
						}
					}
					else if(m->getNumAlleles() > 2 && !opts::_MICROSATS_)
					{
						throw MethodException("More than 2 unique alleles found for map location: " + getString<int>(i) + ", line: " + getString<int>(onind + 1) + ".  Microsatellites not specified.\n");
					}
				}
				else
				{
					samp->addAone(i,true);
					samp->addAtwo(i, true);
					samp->addAmissing(i, true);
				}
				overall++;
    	    	i++;
				if(i > columns){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + options.getPedFile() + "\n";
					text += "Expecting ";
					text += getString<int>(columns + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
				if(overall > columns){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + options.getPedFile() + "\n";
					text += "Expecting ";
					text += getString<int>(columns + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					throw MethodException(text);
				}
            }/*end !linedone*/
        }/*end while(1)*/
		if(gn != columns){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + options.getPedFile() + "\n";
					text += "Expecting ";
					text += getString<int>(columns + 6);
					text += " columns but found ";
					text += getString<int>((f + gn));
					text += "\n";
					throw MethodException(text);
		}

		//new 12-07-2010
		//this closes out the new check for lgen input
        }//end if (!opts::_LGENFILE.length() > 0)

        vector<Methods::Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Methods::Family* fam = new Family();
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
			opts::_FAMILIES_FOUND_++;
        }

		if(filters != NULL){
			bool deleted = false;
			for(int f = 0; f < filters->num_family_filters(); f++){
				bool res = filters->run_family_filter(f, samp->getFamily());
				if(!res && !opts::_KEEP_EXC_SAMPLES_){
					families->erase(families->end() - 1);
					delete(samp->getFamily());
					delete(samp);
					deleted = true;
					break;
				}
				else if(!res && opts::_KEEP_EXC_SAMPLES_){
					samp->getFamily()->setEnabled(false);
					samp->getFamily()->setExcluded(true);
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
			if(deleted){
				continue;
			}
		}

        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
    }/*end while(eof)*/

    fclose(input);

    assignLinks(families);
    if (! opts::_LGENFILE_.length() > 0)
    {
    	reorderAlleles(samples, markers);
    }
    // set DataSet conditions
    set->missing_data_present(missingData);
    set->set_missing_value(3);
    set->set_max_locus(2);
    set->set_max_allele(1);
    set->set_affection_vectors();
	if(reset_todigit){
		opts::_TODIGIT_ = true;
	}

}

void Helpers::readPedM_3vec_set(DataSet* set, StepOptions options){
	if(set == NULL){
		throw MethodException("NULL DataSet structure passed to readPedM_3vec_Set...\n");
	}
	readPedM_3vec_set(set, options, NULL);
}

bool Helpers::realnum(double d){
	double zero = 0;
	if(d != d || d == 1/zero || d == -1/zero)
		return false;
	else
		return true;
}

double Helpers::normdist(double z){
	 double sqrt2pi = 2.50662827463;
	 double t0, z1, p0 ;
	 t0 = 1 / (1 + 0.2316419 * fabs(z));
	 z1 = exp(-0.5 * z*z ) / sqrt2pi;
	 p0 = z1 * t0
	    * (0.31938153 +
	    t0 * (-0.356563782 +
	    t0 * (1.781477937 +
	    t0 * (-1.821255978 +
	    1.330274429 * t0))));
	 return z >= 0 ? 1 - p0 : p0 ;

}

//pulled from Plink
double Helpers::pT(double T, double df){
	if(!realnum(T)){
		return -9;
	}

	T = fabs(T);

	double p, q;
	int st = 0;      //error variable
	int w = 1;       //function variable
	double bnd = 1;  //boundary function

	// NCP is set to 0
	cdft(&w, &p, &q, &T, &df, &st, &bnd);

	//Check status
	if(st != 0){ return -9; }

	// Return two-sided p-value
	return 2 * q;
}

double Helpers::chiprobP(double chi, double df){
	if(chi < 0){
		return -1;
	}
	double pvalue = -1;
	double p, bound;
	int code = 1, status;

	cdfchi(&code, &p, &pvalue, &chi, &df, &status, &bound);
	return pvalue;
}

double Helpers::SQR(double a){
	return a*a;
}

double Helpers::pythag(const double a, const double b){
  double absa,absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void Helpers::svdcmp(vector<vector<double> > & a,
		        vector<double> & w,
				        vector<vector<double> > &v){
  bool flag;
  int i,its,j,jj,k,l = 0,nm;
  double anorm,c,f,g,h,s,scale,x,y,z;
  double volatile temp;

  int m=a.size();
  if (m==0){
	  cerr << "Internal problem in SVD function (no observations left?)\n";
	  throw MethodException("Internal problem in SVD function (no observations left?)\n");
  }
  int n=a[0].size();

  vector<double> rv1(n);
  g=scale=anorm=0.0;
  for (i=0;i<n;i++) {
    l=i+2;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a[k][i]);
      if (scale != 0.0) {
        for (k=i;k<m;k++) {
	      a[k][i] /= scale;
	      s += a[k][i]*a[k][i];
    	}
		f=a[i][i];
		g = -SIGN(sqrt(s),f);
		h=f*g-s;
		a[i][i]=f-g;
		for (j=l-1;j<n;j++) {
		  for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
		  f=s/h;
		  for (k=i;k<m;k++) a[k][j] += f*a[k][i];
		}
		for (k=i;k<m;k++) a[k][i] *= scale;
	  }
	}
	w[i]=scale *g;
	g=s=scale=0.0;
	if (i+1 <= m && i+1 != n) {
	  for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
	  if (scale != 0.0) {
	    for (k=l-1;k<n;k++) {
	      a[i][k] /= scale;
	      s += a[i][k]*a[i][k];
	    }
	    f=a[i][l-1];
	    g = -SIGN(sqrt(s),f);////////////!!!!!!!!!!!!!
		h=f*g-s;
		a[i][l-1]=f-g;
		for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
    	for (j=l-1;j<m;j++) {
		  for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
		  for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
		}
	    for (k=l-1;k<n;k++) a[i][k] *= scale;
	  }
	}
	anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));///////!!!!!!!!!!!!!
  }
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g != 0.0) {
        for (j=l;j<n;j++)
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l;j<n;j++) {
          for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
	      for (k=l;k<n;k++) v[k][j] += s*v[k][i];
        }
      }
      for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=MIN(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) a[i][j]=0.0;
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j<n;j++) {
        for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
        f=(s/a[i][i])*g;
	    for (k=i;k<m;k++) a[k][j] += f*a[k][i];
	  }
	  for (j=i;j<m;j++) a[j][i] *= g;
	} else for (j=i;j<m;j++) a[j][i]=0.0;
	++a[i][i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its<30;its++) {
      flag=true;
      for (l=k;l>=0;l--) {
        nm=l-1;
	    temp=fabs(rv1[l])+anorm;
	    if (temp == anorm) {
	      flag=false;
	      break;
	    }
	    temp=fabs(w[nm])+anorm;
	    if (temp == anorm) break;
	  }
	  if (flag) {
    	c=0.0;
	    s=1.0;
		for (i=l;i<k+1;i++) {
		  f=s*rv1[i];
		  rv1[i]=c*rv1[i];
		  temp = fabs(f)+anorm;
		  if (temp == anorm) break;
		  g=w[i];
		  h=pythag(f,g);
		  w[i]=h;
		  h=1.0/h;
		  c=g*h;
		  s = -f*h;
		  for (j=0;j<m;j++) {
		    y=a[j][nm];
		    z=a[j][i];
		    a[j][nm]=y*c+z*s;
		    a[j][i]=z*c-y*s;
		  }
		}
	  }
	  z=w[k];
	  if (l == k) {
	    if (z < 0.0) {
	      w[k] = -z;
	      for (j=0;j<n;j++) v[j][k] = -v[j][k];
	    }
        break;
	  }
	  if (its == 29){
		  cerr << "SVD function cannot converge: multicollinearity issues?\n";
		  throw MethodException("SVD function cannot converge: multicollinearity issues?\n");
	  }
	  x=w[l];
	  nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rv1[i];
	    y=w[i];
	    h=s*g;
	    g=c*g;
	    z=pythag(f,h);/////////!!!!!!!!!!!!!!!!!
        rv1[j]=z;
	    c=f/z;
	    s=h/z;
	    f=x*c+g*s;
	    g=g*c-x*s;
	    h=y*s;
	    y *= c;
	    for (jj=0;jj<n;jj++) {
	      x=v[jj][j];
	      z=v[jj][i];
	      v[jj][j]=x*c+z*s;
	      v[jj][i]=z*c-x*s;
	    }
	    z=pythag(f,h);
	    w[j]=z;
	    if (z) {
	      z=1.0/z;
	      c=f*z;
	      s=h*z;
	    }
	    f=c*g+s*y;
	    x=c*y-s*g;
	    for (jj=0;jj<m;jj++) {
	      y=a[jj][j];
	      z=a[jj][i];
	      a[jj][j]=y*c+z*s;
	      a[jj][i]=z*c-y*s;
	    }
	  }
	  rv1[l]=0.0;
	  rv1[k]=f;
	  w[k]=x;
	}
  }
}


vector< vector<double> > Helpers::svd_inverse(vector< vector<double> > & u){
  const double eps = 1e-12;
  std::string msga;
  std::string msgb;
  std::stringstream ss;

  if (u.size() == 0){
    cerr << "Internal problem: matrix with no rows (inverse function)\n";
    throw MethodException("Internal problem: matrix with no rows (inverse function)\n");
  }
  if (u.size() != u[0].size() ){
    cerr << "Internal problem: Cannot invert non-square matrix\n" << u.size() << " : " << u[0].size() << "\n";
   ss << u.size();
   msga = ss.str();
   ss << u[0].size();
   msgb = ss.str();

    throw MethodException("Internal problem: Cannot invert non-square matrix\n" + msga + " : " + msgb + "\n");
  }
  int n = u.size();

  vector<double> w(n,0);

  vector<vector<double> > v(n);
  for (int i=0; i<n; i++)
    v[i].resize(n,0);

  svdcmp(u,w,v);//////!!!!!!!!!!!!!!!!!!!!!!!

  // Look for singular values
  double wmax = 0;
  for (int i=0; i<n; i++)
    wmax = w[i] > wmax ? w[i] : wmax;
  double wmin = wmax * eps;
  for (int i=0; i<n; i++)
  {
    w[i] = w[i] < wmin ? 0 : 1/w[i];
  }

  vector<vector<double> > r(n);
  for (int i=0; i<n; i++)
  {
	r[i].resize(n,0);
    for (int j=0; j<n; j++)
      u[i][j] = u[i][j] * w[j];
  }


  for (int i=0; i<n; i++)
	  for (int j=0; j<n; j++)
		  for (int k=0; k<n; k++)
			  r[i][j] += u[i][k] * v[j][k];

  return r;
}

void Helpers::sizeMatrix(vector<vector<double> > &m, int r, int c){
	m.clear();
	m.resize(r);
	for(int i = 0; i < r; i++){
		m[i].resize(c,0);
	}
}


void Helpers::svbksb(vector<vector<double> > &u, vector<double> &w, vector<vector<double> > &v,
		vector<double> &b, vector<double> &x){
	  int jj,j,i;
	  double s;

	  int m=u.size();
	  int n=u[0].size();
	  vector<double> tmp(n);
	  for (j=0;j<n;j++) {
	    s=0.0;
	    if (w[j] != 0.0) {
	      for (i=0;i<m;i++) s += u[i][j]*b[i];
	      s /= w[j];
	    }
	    tmp[j]=s;
	  }

	  for (j=0;j<n;j++) {
	    s=0.0;
	    for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
	    x[j]=s;
	  }

}

void Helpers::multMatrix(vector<vector<double> > & a,
        vector<vector<double> > & b,
        vector<vector<double> > & c)
{

  int ar = a.size();
  int br = b.size();
  if (ar == 0 || br == 0)
    throw MethodException("Internal error: multiplying 0-sized matrices");

  int ac = a[0].size();
  int bc = b[0].size();
  if ( ac != br )
    throw MethodException("Internal error: non-conformable matrices in multMatrix()");

  int cr = ar;
  int cc = bc;

  c.clear();
  sizeMatrix(c,cr,cc);

  for (int i=0; i<ar; i++)
    for (int j=0; j<bc; j++)
      for (int k=0; k<ac; k++)
    c[i][j] += a[i][k] * b[k][j];

}

};
