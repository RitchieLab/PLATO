/*
 * DataLoader.cpp
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#include "DataLoader.h"

#include "data/DataSet.h"
#include "data/Family.h"
#include "data/Sample.h"
#include "data/Marker.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>

using std::multimap;
using std::map;
using std::string;
using std::ifstream;
using std::stringstream;

using PLATO::Data::Marker;
using PLATO::Data::Sample;
using PLATO::Data::DataSet;
using PLATO::Data::Family;

namespace po=boost::program_options;
using po::value;
using po::bool_switch;

namespace PLATO{
namespace Utility{

DataLoader::DataLoader() : _ped_genotype(true), _ped_missing_geno("0"), input(UNKNOWN) {}

po::options_description& DataLoader::addOptions(po::options_description& opts){
	po::options_description data_opts("Data Input Options");

	data_opts.add_options()
		("file", value<string>(&file_base), "Filename base for reading PED/MAP files")
		("ped", value<string>(&ped_fn), "Filename of a PED file")
		("map", value<string>(&map_fn), "Filename of a MAP file")
		("bfile",value<string>(&bfile_base), "Filename base for reading BED/BIM/FAM files")
		("bed", value<string>(&bed_fn), "Filename of a BED file")
		("bim", value<string>(&bim_fn), "Filename of a BIM file")
		("fam", value<string>(&fam_fn), "Filename of a FAM file")
		("tfile",value<string>(&tfile_base), "Filename base for reading TPED/TFAM files")
		("tped", value<string>(&tped_fn), "Filename of a TPED file")
		("tfam", value<string>(&tfam_fn), "Filename of a TFAM file")
		("no-sex", bool_switch(&_ped_no_gender), "PED file does not contain gender")
		("no-parents", bool_switch(&_ped_no_parents), "PED file does not contain parental information")
		("no-fid", bool_switch(&_ped_no_fid), "PED file does not contain fid")
		("no-pheno", bool_switch(&_ped_no_pheno), "PED file does not contain phenotype information")
		("map3", bool_switch(&_map_no_distance), "Specify 3-column MAP file format (no genetic distance)")
		("map-ref", bool_switch(&_map_ref), "Map file contains referent allele information")
		("map-alt", bool_switch(&_map_alt), "Map file contains alternate allele information")
		("control0", bool_switch(&_ped_control0), "Case/Control status is encoded as 0/1, not 1/2")
		("quant", bool_switch(&_quant), "Phenotype is a quantitative value, not case/control status (unparseable is converted to missing)");

	return opts.add(data_opts);
}

void DataLoader::parseOptions(const po::variables_map& vm){
	if(vm.count("file")){
		if(ped_fn.size() == 0){
			ped_fn = file_base + ".ped";
		}
		if(map_fn.size() == 0){
			map_fn = file_base + ".map";
		}
	}

	if(vm.count("bfile")){
		if(bed_fn.size() == 0){
			bed_fn = bfile_base + ".bed";
		}
		if(bim_fn.size() == 0){
			bim_fn = bfile_base + ".bim";
		}
		if(fam_fn.size() == 0){
			fam_fn = bfile_base + ".fam";
		}
	}

	if(vm.count("tfile")){
		if(tped_fn.size() == 0){
			tped_fn = tfile_base + ".tped";
		}
		if(tfam_fn.size() == 0){
			tfam_fn = tfile_base + ".tfam";
		}
	}

	// decide what the style of input we're reading
	if(ped_fn.size() && map_fn.size()){
		input = PED;
	}else if(bed_fn.size() && bim_fn.size() && fam_fn.size()){
		input = BED;
	}else if(tped_fn.size() && tfam_fn.size()){
		input = TPED;
	}

}

void DataLoader::read(DataSet& ds){

	ds_ptr = &ds;

	// first, reading PED/MAP files
	switch(input){
	case PED:
		readMap(map_fn);
		readPed(ped_fn);
		break;
	case BED:
		// the PED file no longer contains genotype information
		_ped_genotype = false;
		// the BIM file contains ref and alt information
		_map_ref = _map_alt = true;

		// NOTE: all the options dictating the map file format are still in play!
		readMap(bim_fn);
		readPed(fam_fn);
		readBinPed(bed_fn);
		break;
	case TPED:
		// No genotype information in the PED
		_ped_genotype = false;

		readPed(tfam_fn);
		readTPed(tped_fn);
		break;
	default:
		throw std::logic_error("Unknown input type");
	}

}


void DataLoader::readMap(const string& fn){

	ifstream input(fn.c_str());

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

    ifstream input(fn.c_str());

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

		if(!_ped_no_fid){
			s >> iid;
		}else{
			iid = fid;
			fid = "";
		}

		if(!_ped_no_parents){
			s >> dad_id;
			s >> mom_id;
		}

		if(!_ped_no_gender){
			s >> sex;
		}

		if(!_ped_no_pheno){
			s >> pheno;
		}

		bool use = true;

		// check for filters at this point, before we add the sample
		if(use){
			Sample* samp;
			if(!_ped_no_fid){
				samp = ds_ptr->addSample(fid, iid, _marker_incl.count());
				family_map.insert(std::make_pair(fid, iid));
			}else{
				samp = ds_ptr->addSample(iid, _marker_incl.count());
			}

			if(!_ped_no_gender){
				if(sex == "1"){
					samp->setGender(true);
				}else if(sex == "2"){
					samp->setGender(false);
				}// anything other than 1/2 is unknown gender
			}

			if(!_ped_no_pheno){
				if(_quant){
					float f;
					if(stringstream(pheno) >> f){
						samp->setPheno(f);
					}
				} else {
					// kept this way to preserve "unknown" status
					if(pheno == (_ped_control0 ? "0" : "1")){
						samp->setAffected(false);
					}else if(pheno == (_ped_control0 ? "1" : "2")){
						samp->setAffected(true);
					}
				}
			}

			if(!_ped_no_parents){
				mom_map[iid] = mom_id;
				dad_map[iid] = dad_id;
			}

			if (_ped_genotype) {
				// Now, load up the genotypes
				DataSet::marker_iterator mi = ds_ptr->beginMarker();
				DataSet::marker_iterator mi_end = ds_ptr->endMarker();
				for (unsigned int i = 0; i < _marker_incl.size(); i++) {
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
			if((mom_itr == mom_map.end() || id_map.find((*mom_itr).second) == id_map.end()) &&
			   (dad_itr == dad_map.end() || id_map.find((*dad_itr).second) == id_map.end())){
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

	ifstream BIT(fn.c_str(), std::ios::in | std::ios::binary);

	if(!BIT.good()){
		throw std::invalid_argument("Error opening genotype bitfile: " + fn);
	}

	char byte_val;
	char magic[2];
	//header
	BIT.read(magic, 2);

	bool snp_major = false;
	if(magic[0] == 108 && magic[1] == 27){
		BIT.read(&byte_val, 1);
		snp_major = static_cast<bool>(byte_val);
	} else {
		//throw std::invalid_argument("Incorrect magic number in Binary PED file: " + fn);
		// Maybe just warn here - assume v0.99 BED
		BIT.seekg(0);
	}

	const unsigned char GENO_1 = 2;
	const unsigned char GENO_2 = 1;

	if(snp_major){
		// OK, pad the sample size to a multiple of 4 (just "ignore" those samples)
		while(_sample_incl.size() % 4 != 0){
			_sample_incl.push_back(false);
		}
		for(unsigned int m = 0; m < _marker_incl.size(); m++){
			if(!_marker_incl[m]){
				// If I'm not including this marker, seek the appropriate number
				// of bytes.
				BIT.seekg(_sample_incl.size() / 4, BIT.cur);
			} else {
				DataSet::sample_iterator si = ds_ptr->beginSample();
				// I KNOW that _sample_incl.size() is a multiple of 4!
				for (unsigned int s = 0; s * 4 < _sample_incl.size(); s++) {
					char ch;
					BIT.read(&ch, 1);
					if (!BIT.good()) {
						throw std::logic_error(
								"Problem with the bed file - unexpected EOF.");
					}

					for (int i = 0; i < 4; i++) {
						if (_sample_incl[s * 4 + i]) {
							if (si == ds_ptr->endSample()) {
								throw std::logic_error(
										"Problem with the bed file - not enough samples.");
							}
							if ((ch & GENO_2) && !(ch & GENO_1)) {
								(*si)->appendMissingGenotype();
							} else {
								(*si)->appendGenotype((ch & GENO_1) >> 1, ch
										& GENO_2);
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
		DataSet::sample_iterator si = ds_ptr->beginSample();
		// OK, must be individual-major then!
		for (unsigned int s = 0; s < _sample_incl.size(); s++) {
			if(!_sample_incl[s]){
				BIT.seekg(_marker_incl.size() / 4, BIT.cur);
			} else {
				if (si == ds_ptr->endSample()) {
					throw std::logic_error("Problem with the bed file - not enough samples.");
				}

				// I KNOW that _marker_incl.size() is a multiple of 4!
				for (unsigned int m = 0; m * 4 < _marker_incl.size(); s++) {
					char ch;
					BIT.read(&ch, 1);
					if (!BIT.good()) {
						throw std::logic_error("Problem with the bed file - unexpected EOF.");
					}

					for (int i = 0; i < 4; i++) {
						if (_marker_incl[m * 4 + i]) {
							if ((ch & GENO_2) && !(ch & GENO_1)) {
								(*si)->appendMissingGenotype();
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
    ifstream input(fn.c_str());

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
				for(unsigned int i=0; i<_sample_incl.size(); i++){
					// now, go through and parse each sample.
					string g1, g2;
					if(!s.good()){
						throw std::logic_error("Error: not enough samples in TPed");
					}
					s >> g1;
					s >> g2;


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
	if (!_map_no_distance) {
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
		m = ds_ptr->addMarker(chr, static_cast<unsigned int> (bploc),
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
		s->appendMissingGenotype();
	} else {
		unsigned char g_idx1 =
				static_cast<unsigned char> (m->addAllele(g1));
		unsigned char g_idx2 =
				static_cast<unsigned char> (m->addAllele(g2));

		s->appendGenotype(g_idx1, g_idx2);
	}

}

}
}
