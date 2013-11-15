/*
 * DataLoader.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "DataLoader.h"

#include "DataSet.h"
#include "Family.h"
#include "Sample.h"
#include "Marker.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>

using std::multimap;
using std::map;
using std::string;
using std::ifstream;
using std::stringstream;

namespace Methods{

void DataLoader::readMap(const string& fn){

	DataSet& dataset = *ds_ptr;

	ifstream input(fn);

    if(!input.is_open()){
    	throw std::invalid_argument("Error opening map file: " + fn);
    }

	string line;

    while(getline(input, line)){

    	stringstream s(line);

    	Marker* newMarker=parseMap(s);

    	// We return -1 if it's a line to ignore
    	if(newMarker != reinterpret_cast<Marker*>(-1)){
    		_marker_incl.push_back(newMarker != 0);
    	}

    }

    input.close();

}

void DataLoader::readPed(const string& fn){
	// a listing of ID -> mother's ID (and father's ID, accordingly)
	map<string, string> mom_map;
	map<string, string> dad_map;

	// a listing of ID -> samples
	map<string, Sample*> id_map;

	// a multiset of family ID -> individual ID
	multimap<string, string> family_map;

    ifstream input(fn);

	if(!input.is_open()){
		throw std::invalid_argument("Error opening pedfile: " + fn);
	}

	string line;

	int lineno = -1;
	while(getline(input, line)){
		++lineno;

		stringstream s(line);

		string fid;
		s >> fid;

		// skip empty lines and those beginning with #
		if(fid.length() == 0 || fid[0] == '#'){
			continue;
		}

		string iid;
		string dad_id;
		string mom_id;
		string sex;
		string pheno;

		if(_ped_fid){
			s >> iid;
		}else{
			iid = fid;
			fid = "";
		}

		if(_ped_parents){
			s >> dad_id;
			s >> mom_id;
		}

		if(_ped_gender){
			s >> sex;
		}

		if(_ped_pheno){
			s >> pheno;
		}

		bool use = true;

		// check for filters at this point, before we add the sample
		if(use){
			Sample* samp;
			if(_ped_fid){
				samp = ds_ptr->addSample(fid, iid, _marker_incl.count());
				family_map.insert(fid, iid);
			}else{
				samp = ds_ptr->addSample(iid, _marker_incl.count());
			}

			if(_ped_gender){
				if(sex == "1"){
					samp->setGender(true);
				}else if(sex == "2"){
					samp->setGender(false);
				}// anything other than 1/2 is unknown gender
			}

			if(_ped_pheno){
				if(pheno == (_ped_control1 ? "1" : "0")){
					samp->setAffected(false);
				}else if(pheno == (_ped_control1 ? "2" : "1")){
					samp -> setAffected(true);
				}
			}

			if(_ped_parents){
				mom_map[iid] = mom_id;
				dad_map[iid] = dad_id;
			}

			if (_ped_genotype) {
				// Now, load up the genotypes
				DataSet::marker_iterator mi = ds_ptr->beginMarker();
				DataSet::marker_iterator mi_end = ds_ptr->endMarker();
				for (int i = 0; i < _marker_incl.size(); i++) {
					string g1, g2;
					s >> g1;
					s >> g2;
					if (_marker_incl[i]) {
						if (mi == mi_end) {
							//This is a problem!!
							throw std::logic_error("Error: marker not loaded!");
						}

						parseSample(*mi, samp, g1, g2);

						// Advance the marker
						++mi;
					}
				}
				if (mi != mi_end) {
					throw std::logic_error("Error: not enough markers!");
				}
			}

		}

		_sample_incl.push_back(use);

	}

	input.close();

	// If we had family IDs, add them now
	multimap<string, string>::const_iterator fam_itr = family_map.begin();
	Family* curr_fam_ptr = 0;
	while(fam_itr != family_map.end()){
		// Make sure the sample was created first
		map<string, Sample*>::const_iterator samp_itr = id_map.find((*fam_itr).second);

		if (samp_itr != id_map.end()){
			if(!curr_fam_ptr || curr_fam_ptr->getFamID() != (*fam_itr).first) {
				curr_fam_ptr = ds_ptr->addFamily((*fam_itr).first);
			}

			// Check for founder status of the current ID.
			map<string, string>::const_iterator mom_itr = mom_map.find((*fam_itr).second);
			map<string, string>::const_iterator dad_itr = dad_map.find((*fam_itr).second);

			//If both are not found, then we have a founder!
			if((mom_itr == mom_map.end() || id_map.find((*mom_itr).second) == id_map.end) &&
			   (dad_itr == dad_map.end() || id_map.find((*dad_itr).second) == id_map.end)){
				curr_fam_ptr->addFounder((*samp_itr).second);
			}else{
				curr_fam_ptr->addNonFounder((*samp_itr).second);
			}
		}
		++fam_itr;
	}

	// If we had parental information, add it here
	map<string, string>::const_iterator mom_itr = mom_map.begin();
	while(mom_itr != mom_map.end()){
		// Get the mother and child samples from the ID map
		map<string, Sample*>::const_iterator mom_id_itr = id_map.find((*mom_itr).second);
		map<string, Sample*>::const_iterator kid_id_itr = id_map.find((*mom_itr).first);

		if(mom_id_itr != id_map.end() && kid_id_itr != id_map.end()){
			(*mom_id_itr).second->addChild((*kid_id_itr).second);
			(*kid_id_itr).second->addMother((*mom_id_itr).second);
		}
		++mom_itr;
	}

	map<string, string>::const_iterator dad_itr = dad_map.begin();
	while(dad_itr != dad_map.end()){
		// Get the mother and child samples from the ID map
		map<string, Sample*>::const_iterator dad_id_itr = id_map.find((*dad_itr).second);
		map<string, Sample*>::const_iterator kid_id_itr = id_map.find((*dad_itr).first);

		if(dad_id_itr != id_map.end() && kid_id_itr != id_map.end()){
			(*dad_id_itr).second->addChild((*kid_id_itr).second);
			(*kid_id_itr).second->addFather((*dad_id_itr).second);
		}
		++dad_itr;
	}
}

void DataLoader::readBinPed(const string& fn){

	ifstream BIT(fn, std::ios::in | std::ios::binary);

	if(!BIT.good()){
		throw std::invalid_argument("Error opening genotype bitfile: " + fn);
	}

	unsigned char byte_val;
	unsigned short magic;
	//header
	BIT.read(&magic, 2);

	bool snp_major = false;
	if(magic == 27675){
		BIT.read(&byte_val, 1);
		snp_major = static_cast<bool>(byte_val);
	} else {
		//throw std::invalid_argument("Incorrect magic number in Binary PED file: " + fn);
		// Maybe just warn here - assume v0.99 BED
		BIT.seekg(0);
	}

	const unsigned char GENO_1 = 1;
	const unsigned char GENO_2 = 2;

	if(snp_major){
		// OK, pad the sample size to a multiple of 4 (just "ignore" those samples)
		while(_sample_incl.size() % 4 != 0){
			_sample_incl.push_back(false);
		}
		for(int m = 0; m < _marker_incl.size(); m++){
			if(!_marker_incl[m]){
				// If I'm not including this marker, seek the appropriate number
				// of bytes.
				BIT.seekg(_sample_incl.size() / 4, BIT.cur);
			} else {
				DataSet::sample_iterator si = ds_ptr->sampleBegin();
				// I KNOW that _sample_incl.size() is a multiple of 4!
				for (int s = 0; s * 4 < _sample_incl.size(); s++) {
					char ch;
					BIT.read(&ch, 1);
					if (!BIT.good()) {
						throw std::logic_error(
								"Problem with the bed file - unexpected EOF.");
					}

					for (int i = 0; i < 4; i++) {
						if (_sample_incl[s * 4 + i]) {
							if (si == ds_ptr->sampleEnd()) {
								throw std::logic_error(
										"Problem with the bed file - not enough samples.");
							}
							if ((ch & GENO_2) && !(ch & GENO_1)) {
								(*si)->appendMissingGenotype;
							} else {
								(*si)->appendGenotype((ch & GENO_2) >> 1, ch
										& GENO_1);
							}

							++si;
						}
						// shift bits two places and go again!
						ch >>= 2;
					}
				}
			}
		}
	} else {
		while(_marker_incl.size() % 4 != 0){
			_marker_incl.push_back(false);
		}
		DataSet::sample_iterator si = ds_ptr->sampleBegin();
		// OK, must be individual-major then!
		for (int s = 0; s < _sample_incl.size(); s++) {
			if(!_sample_incl[s]){
				BIT.seekg(_marker_incl.size() / 4, BIT.cur);
			} else {
				if (si == ds_ptr->sampleEnd()) {
					throw std::logic_error("Problem with the bed file - not enough samples.");
				}

				// I KNOW that _marker_incl.size() is a multiple of 4!
				for (int m = 0; m * 4 < _marker_incl.size(); s++) {
					char ch;
					BIT.read(&ch, 1);
					if (!BIT.good()) {
						throw std::logic_error("Problem with the bed file - unexpected EOF.");
					}

					for (int i = 0; i < 4; i++) {
						if (_marker_incl[m * 4 + i]) {
							if ((ch & GENO_2) && !(ch & GENO_1)) {
								(*si)->appendMissingGenotype;
							} else {
								(*si)->appendGenotype((ch & GENO_2) >> 1, ch & GENO_1);
							}
						}
						// shift bits two places and go again!
						ch >>= 2;
					}
				}
				++si;
			}
		}

	}

	BIT.close();

}

void DataLoader::readTPed(const string& fn){
	// We'll assume that we've already read the "tfam", which is really just a
	// "ped" with no genotype information.
    ifstream input(fn);

	if(!input.is_open()){
		throw std::invalid_argument("Error opening pedfile: " + fn);
	}

	string line;

	int lineno = -1;
	while(getline(input, line)){
		++lineno;

		stringstream s(line);
		Marker* newMarker = parseMap(s);

		if(newMarker != reinterpret_cast<Marker*>(-1)){
			_marker_incl.push_back(newMarker != 0);

			if(newMarker){
				DataSet::sample_iterator si = ds_ptr->beginSample();
				for(int i=0; i<_sample_incl.size(); i++){
					// now, go through and parse each sample.
					string g1, g2;
					s >> g1;
					s >> g2;
					if(!s.good()){
						throw std::logic_error("Error: not enough samples in TPed");
					}

					if(_sample_incl[i]){
						if(si == ds_ptr->endSample()){
							throw std::logic_error("Error: too many samples in TPed");
						}

						parseSample(newMarker, *si, g1, g2);
						++si;
					}
				}
			}
		}


	}
}

Marker* DataLoader::parseMap(stringstream& ss) {
	string chr;
	ss >> chr;

	// skip empty lines or those that begin with "#"
	if (chr.size() == 0 || chr[0] == '#') {
		return reinterpret_cast<Marker*>(-1);
	}
	string id;
	int bploc;
	int distance;
	string ref;
	string alt;

	ss >> id;
	if (_map_distance) {
		ss >> distance;
	}
	ss >> bploc;
	if (_map_ref) {
		ss >> ref;
	}
	if (_map_alt) {
		ss >> alt;
	}

	// If location is negative, exclude the marker
	bool use = bploc > 0;

	// check for other conditions of use here!

	Marker* m = 0;
	if (use) {
		m = dataset.addMarker(chr, static_cast<unsigned int> (bploc),
				id);

		if (_map_ref) {
			m->addAllele(ref);
			m->setRefAllele(ref);
		}
		if (_map_alt) {
			m->addAllele(alt);
			m->setAltAllele(alt);
		}
	}
	return m;
}

void DataLoader::parseSample(Marker* m, Sample* s, const string& g1, const string& g2){

	if (g1 == _ped_missing_geno || g2 == _ped_missing_geno) {
		samp->appendMissingGenotype();
	} else {
		unsigned char g_idx1 =
				static_cast<unsigned char> (m->addAllele(g1));
		unsigned char g_idx2 =
				static_cast<unsigned char> (m->addAllele(g2));

		samp->appendGenotype(g1, g2);
	}

}

}

