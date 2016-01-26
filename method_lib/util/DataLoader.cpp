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

#include "Logger.h"
#include "InputManager.h"
#include "ICompressedFile.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <set>
#include <limits>
#include <deque>
#include <algorithm>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

using std::multimap;
using std::map;
using std::set;
using std::string;
using std::ifstream;
using std::stringstream;
using std::vector;
using std::deque;
using std::find;
using std::pair;

using PLATO::Data::Marker;
using PLATO::Data::Sample;
using PLATO::Data::DataSet;
using PLATO::Data::Family;
using PLATO::Utility::Logger;

namespace po=boost::program_options;
namespace fs=boost::filesystem;
using po::value;
using po::bool_switch;
using boost::regex;
using boost::regex_match;
using boost::algorithm::split;

namespace PLATO{
namespace Utility{

// Define the "Group Separator" string as the separation between FID and IID!
// Do this b/c it's legal ASCII, but darn near impossible to include on the cmd line
const string DataLoader::sampl_field_sep = "\x1d";

DataLoader::DataLoader() : _map_others(false), _ped_genotype(true),
		_ped_missing_geno("0"), input(UNKNOWN) {}

po::options_description& DataLoader::addOptions(po::options_description& opts){
	po::options_description data_opts("Data Input Options");

	po::options_description plink_opts("PLINK file input options");

	plink_opts.add_options()
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
		("lfile", value<string>(&lfile_base), "Filename base for reading LGEN/MAP/FAM files")
		("lgen", value<string>(&lgen_fn), "Filename of an LGEN file")
		("no-sex", bool_switch(&_ped_no_gender), "PED file does not contain gender")
		("no-parents", bool_switch(&_ped_no_parents), "PED file does not contain parental information")
		("no-fid", bool_switch(&_ped_no_fid), "PED file does not contain fid")
		("no-pheno", bool_switch(&_ped_no_pheno), "PED file does not contain phenotype information")
		("map3", bool_switch(&_map_no_distance), "Specify 3-column MAP file format (no genetic distance)")
		("map-ref", bool_switch(&_map_ref), "Map file contains referent allele information")
		("map-alt", bool_switch(&_map_alt), "Map file contains alternate allele information")
		("control0", bool_switch(&_ped_control0), "Case/Control status is encoded as 0/1, not 1/2")
		("quant", bool_switch(&_quant), "Phenotype is a quantitative value, not case/control status (unparseable is converted to missing)")
		;

	data_opts.add(plink_opts);

	po::options_description beagle_opts("BEAGLE file input options");

	beagle_opts.add_options()
		("beagle-prefix", value<string>(&beagle_prefix), "Prefix for all BEAGLE files")
		("beagle-suffix", value<string>(&beagle_suffix)->default_value("bgl"), "Suffix for all BEAGLE genotype files")
		("marker-suffix", value<string>(&marker_suffix)->default_value("markers"), "Suffix for all BEAGLE marker files")
		("beagle-files", value<vector<string> >(&beagle_fns)->composing(), "Comma-separated list of all BEAGLE genotype files")
		("marker-files", value<vector<string> >(&marker_fns)->composing(), "Comma-separated list of all BEAGLE marker files")
		("beagle-chroms", value<vector<string> >(&chroms)->composing(), "Comma-separated list of chromosomes corresponding to genotype/marker files")
		("trio", bool_switch(&bgl_trio), "Genotype data given is trio data")
		("pair", bool_switch(&bgl_pair), "Genotype data given is pair data")
		("bgl-phased", bool_switch(&bgl_phased), "Genotype data is phased")
		("bgl-poly", bool_switch(&bgl_poly), "Genotype data is polyallelic")
		("beagle-missing", value<string>(&bgl_missing_allele)->default_value("0"), "String for missing allele")
		;

	data_opts.add(beagle_opts);

	po::options_description vcf_opts("VCF file input options");

	vcf_opts.add_options()
		("vcf-file", value<string>(&vcf_fn), "VCF file to load")
		("no-filter-marker", bool_switch(&no_filter_marker), "Load markers that fail a filter in the FILTER field")
		("no-filter-geno", bool_switch(&no_filter_geno), "Load genotypes that fail a filter in the FT annotation field")
		("vcf-phased", bool_switch(&vcf_phased), "Flag indicating that the VCF file is phased")
		("vcf-poly", bool_switch(&vcf_poly), "Flag indicating that the VCF file is polyalleleic")
		;

	data_opts.add(vcf_opts);

	po::options_description filter_opts("Pre-filtering options");
	filter_opts.add_options()
		("chrom", value<vector<string> >(&chrom_list)->composing(), "Chromosome(s) to include")
		("bp-window", value<vector<string> >(&bp_window)->composing(), "Base pair window(s) to include")
		("incl-marker", value<vector<string> >(&incl_marker_str)->composing(), "Marker ID(s) to include")
		("excl-marker", value<vector<string> >(&excl_marker_str)->composing(), "Marker ID(s) to exclude")
		("incl-marker-fn", value<vector<string> >(&incl_marker_fns)->composing(), "File containing marker ID(s) to include")
		("excl-marker-fn", value<vector<string> >(&excl_marker_fns)->composing(), "File containing marker ID(s) to exclude")
		("incl-sample", value<vector<string> >(&incl_sample_str)->composing(), "Sample(s) to include")
		("excl-sample", value<vector<string> >(&excl_sample_str)->composing(), "Sample(s) to exclude")
		("incl-sample-fn", value<vector<string> >(&incl_marker_fns)->composing(), "File containing Sample(s) to include")
		("excl-sample-fn", value<vector<string> >(&excl_marker_fns)->composing(), "File containing Sample(s) to exclude")
		;

	data_opts.add(filter_opts);

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

	if(vm.count("lfile")){
		if(lgen_fn.size() == 0){
			lgen_fn = lfile_base + ".lgen";
		}
		if(map_fn.size() == 0){
			map_fn = lfile_base + ".map";
		}
		if(fam_fn.size() == 0){
			fam_fn = lfile_base + ".fam";
		}
	}

	_fns_provided[0] = beagle_fns.size() > 0;
	_fns_provided[1] = marker_fns.size() > 0;
	_fns_provided[2] = chroms.size() > 0;

	if(_fns_provided.any()){
		if(_fns_provided.count() != _fns_provided.size()){
			Logger::log_err("ERROR: If you provide filenames for BEAGLE files, you must provide all genotype, marker, and chromosome files", true);
		}

		vector<string> tmp_fnvec;
		InputManager::parseInput(beagle_fns, tmp_fnvec);
		beagle_fns = tmp_fnvec;
		tmp_fnvec.clear();
		InputManager::parseInput(marker_fns, tmp_fnvec);
		marker_fns = tmp_fnvec;
		tmp_fnvec.clear();
		InputManager::parseInput(chroms, tmp_fnvec);
		chroms = tmp_fnvec;

		if(beagle_fns.size() != marker_fns.size() || marker_fns.size() != chroms.size()){
			Logger::log_err("ERROR: genotypes, markers, and chromosome lists must be the same size", true);
		}
		if(beagle_prefix.size() > 0){
			Logger::log_err("WARNING: beagle files were given along with a beagle prefix, ignoring the prefix directive");
			beagle_prefix = "";
		}

	} else if(beagle_prefix.size() > 0){

		// OK, we need to build the file lists from directory listings
		// also, we need to make sure that the chromosomes match!

		// first, build a file object from the prefix (to include any directories)
		fs::path prefix_dir = fs::absolute(beagle_prefix).parent_path();
		if(!fs::is_directory(prefix_dir)){
			Logger::log_err("ERROR: could not find the directory containing BEAGLE files", true);
		}

		set<string> beagle_chroms;
		set<string> marker_chroms;
		// now, list all files in the directory of the above path
		fs::directory_iterator end_itr;

		regex esc("[\\^\\.\\$\\|\\(\\)\\[\\]\\*\\+\\?\\/\\\\]");
		const string rep("\\\\\\1&");
		string bgl_prefix_regex = boost::regex_replace(fs::path(beagle_prefix).filename().string() + ".", esc, rep,boost::match_default | boost::format_sed);
		string bgl_suffix_regex = boost::regex_replace("." + beagle_suffix, esc, rep, boost::match_default | boost::format_sed);
		string bgl_marker_regex = boost::regex_replace("." + marker_suffix, esc, rep, boost::match_default | boost::format_sed);

		regex bgl_re(bgl_prefix_regex + "(.*)" + bgl_suffix_regex + "\\.?(z|gz|bz)?");
		regex marker_re(bgl_prefix_regex + "(.*)" + bgl_marker_regex + "\\.?(z|gz|bz)?");

		boost::smatch match_obj;

		for(fs::directory_iterator dir_itr(prefix_dir); dir_itr != end_itr; dir_itr++){
			if(fs::is_regular_file(dir_itr->status())){
				// check for regex of beagle file
				if(regex_match(dir_itr->path().filename().string(), match_obj, bgl_re)){
					beagle_chroms.insert(match_obj[1]);
				// check for regex of marker file
				} else if (regex_match(dir_itr->path().filename().string(), match_obj, marker_re)){
					marker_chroms.insert(match_obj[1]);
				}
			}
		}

		if(!(beagle_chroms == marker_chroms)){
			Logger::log_err("ERROR: Chromosomes for marker and genotype files do not match", true);
		} else if(beagle_chroms.size() == 0){
			Logger::log_err("ERROR: could not find BEAGLE genotype or marker files", true);
		}else {

			chroms.clear();
			chroms.reserve(beagle_chroms.size());
			chroms.insert(chroms.begin(), beagle_chroms.begin(), beagle_chroms.end());
			beagle_fns.clear();
			beagle_fns.reserve(beagle_chroms.size());
			marker_fns.clear();
			marker_fns.reserve(beagle_chroms.size());
			for(unsigned int i=0; i<chroms.size(); i++){
				beagle_fns.push_back(beagle_prefix + "." + chroms[i] + "." + beagle_suffix);
				marker_fns.push_back(beagle_prefix + "." + chroms[i] + "." + marker_suffix);
			}
		}
	}

	// decide what the style of input we're reading
	if(ped_fn.size() && map_fn.size()){
		input = PED;
	}else if(bed_fn.size() && bim_fn.size() && fam_fn.size()){
		input = BED;
	}else if(tped_fn.size() && tfam_fn.size()){
		input = TPED;
	}else if(lgen_fn.size() && map_fn.size() && fam_fn.size()){
		input = LGEN;
	}else if(beagle_fns.size()){
		input = BEAGLE;
	}else if(vcf_fn.size()){
		input = VCF;
	}

	// load up some pre-filters, please!
	if(chrom_list.size() > 0){
		vector<string> all_chroms;
		InputManager::parseInput(chrom_list,all_chroms);
		for(unsigned int i=0; i<all_chroms.size(); i++){
			unsigned short chr = InputManager::chrStringToInt(all_chroms[i]);
			if(chr != 0){
				chrom_ids.insert(chr);
			}
		}
		// no sense taking up unneeded space
		chrom_list.clear();
	}

	if(bp_window.size() > 0){
		vector<string> all_window;
		InputManager::parseInput(bp_window, all_window);

		// fix this, please!!!
		vector<string> single_window;
		for(unsigned int i=0; i<all_window.size(); i++){
			boost::algorithm::split(single_window, all_window[i], boost::is_any_of("-"));
			if(single_window.size() != 2){
				Logger::log_err("WARNING: Base pair window '" + all_window[i] +
						"' does not have exactly one dash (-), ignoring.");
			} else {
				unsigned int start = std::numeric_limits<unsigned int>::max();
				unsigned int stop = std::numeric_limits<unsigned int>::max();
				bool good = true;
				if(single_window[0] != ""){
					stringstream(single_window[0].c_str()) >> start;
					if(start == std::numeric_limits<unsigned int>::max()){
						good = false;
						Logger::log_err("WARNING: could not parse '" + single_window[0] + "', ignoring window.");
					}
				} else {
					start = 0;
				}

				if(single_window[1] != ""){
					stringstream(single_window[1].c_str()) >> stop;
					if(stop == std::numeric_limits<unsigned int>::max()){
						good = false;
						Logger::log_err("WARNING: could not parse '" + single_window[0] + "', ignoring window.");
					}
				}

				if(good){
					bp_window_set.insert(std::make_pair(start, stop));
				}
			}
		}

		if(!bp_window_set.empty()){
			// OK, now make sure everything is not overlapping
			set<pair<unsigned int, unsigned int> >::iterator bp_it = bp_window_set.begin();
			set<pair<unsigned int, unsigned int> >::iterator bp_next = ++bp_window_set.begin();

			while(bp_next != bp_window_set.end()){
				if((*bp_it).second >= (*bp_next).first){
					// insert a new element (which I KNOW will go in between
					// bp_it and bp_next!
					bp_window_set.insert(std::make_pair((*bp_it).first, (*bp_next).second));
					// erase and increment bp_it and bp_next
					bp_window_set.erase(bp_it++);
					bp_window_set.erase(bp_next++);
					// NOTE: at this point bp_it will point to the newly inserted
					// element, and bp_next will point to the next element!
				} else{
					++bp_next;
					++bp_it;
				}
			}
		}
	}

	readSampleList(incl_sample_str, incl_sample_set);
	readSampleList(excl_sample_str, excl_sample_set);
	readSampleFile(incl_sample_fns, incl_sample_set);
	readSampleFile(excl_sample_fns, excl_sample_set);

	InputManager::parseInput(incl_marker_str, incl_marker_set);
	readMarkerFile(incl_marker_fns, incl_marker_set);
	InputManager::parseInput(excl_marker_str, excl_marker_set);
	readMarkerFile(excl_marker_fns, excl_marker_set);
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
	case LGEN:
		_ped_genotype = false;
		readMap(map_fn);
		readPed(fam_fn);
		readLGen(lgen_fn);
		break;
	case BEAGLE:
		readBeagle();
		break;
	case VCF:
		readVCF();
		break;
	default:
		Logger::log_err("Unknown input type", true);
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

		bool use = filterSample(iid, fid);

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
							Logger::log_err("Error: marker not loaded!", true);
						}

						parseSample(*mi, samp, g1, g2, false);

						// Advance the marker
						++mi;
					}
				}
				if (mi != mi_end) {
					Logger::log_err("Error: not enough markers!", true);
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
		Logger::log_err<std::invalid_argument>("Error opening genotype bitfile: " + fn, true);
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
		Logger::log_err("WARNING: BED magic number not found, assuming v0.99 BED file");
		BIT.seekg(0);
	}

	const unsigned char GENO_1 = 2;
	const unsigned char GENO_2 = 1;

	if(snp_major){
		// OK, pad the sample size to a multiple of 4 (just "ignore" those samples)
		while(_sample_incl.size() % 4 != 0){
			_sample_incl.push_back(false);
		}
		DataSet::marker_iterator mi = ds_ptr->beginMarker();
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
						Logger::log_err("Problem with the bed file - unexpected EOF.", true);
					}

					for (int i = 0; i < 4; i++) {
						if (_sample_incl[s * 4 + i]) {
							if (si == ds_ptr->endSample()) {
								Logger::log_err("Problem with the bed file - not enough samples.", true);
							}
							if ((ch & GENO_2) && !(ch & GENO_1)) {
								(*si)->setMissingGenotype(**mi);
							} else {
								(*si)->setGenotype(**mi, (ch & GENO_1) >> 1, ch
										& GENO_2);
							}

							++si;
						}
						// shift bits two places and go again!
						ch >>= 2;
					}
				}
				++mi;
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
					Logger::log_err("Problem with the bed file - not enough samples.", true);
				}

				// I KNOW that _marker_incl.size() is a multiple of 4!
				DataSet::marker_iterator mi = ds_ptr->beginMarker();
				for (unsigned int m = 0; m * 4 < _marker_incl.size(); s++) {
					char ch;
					BIT.read(&ch, 1);
					if (!BIT.good()) {
						Logger::log_err("Problem with the bed file - unexpected EOF.", true);
					}

					for (int i = 0; i < 4; i++) {
						if (_marker_incl[m * 4 + i]) {
							if ((ch & GENO_2) && !(ch & GENO_1)) {
								(*si)->setMissingGenotype(**mi);
							} else {
								(*si)->setGenotype(**mi,(ch & GENO_2) >> 1, ch & GENO_1);
							}
						}
						// shift bits two places and go again!
						ch >>= 2;
						if(mi != ds_ptr->endMarker()){
							++mi;
						}
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
		Logger::log_err<std::invalid_argument>("Error opening TPEDfile: " + fn, true);
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

void DataLoader::readLGen(const std::string& fn){
	ifstream input(fn.c_str());

	if(!input.is_open()){
		Logger::log_err<std::invalid_argument>("Error opening LGEN file: " + fn, true);
	}

	string line;

	int lineno = 0;
	while(getline(input, line)){
		++lineno;

		boost::algorithm::trim(line);
		if(line[0] != '#'){
			stringstream s(line);
			string fid, iid, mid, g1, g2;

			if(!_ped_no_fid){
				s >> fid;
			}
			s >> iid >> mid >> g1 >> g2;
			Sample* samp;
			if(!_ped_no_fid){
				samp = ds_ptr->getSample(fid, iid);
			} else {
				samp = ds_ptr->getSample(iid);
			}
			if(!samp){
				Logger::log_err("ERROR: Sample not found on line " + boost::lexical_cast<string>(lineno), true);
			}
			Marker* mark = ds_ptr->getMarker(mid);
			if(!mark){
				Logger::log_err("ERROR: Marker not found on line " + boost::lexical_cast<string>(lineno), true);
			}

			parseSample(mark, samp, g1, g2, false);
		}

	}

	input.close();
}

void DataLoader::readBeagle(){

	string currChrom;
	ICompressedFile genof;
	ICompressedFile markerf;
	_map_others = _map_ref = _map_alt = true;
	_ped_missing_geno = bgl_missing_allele;
	_map_no_distance = true;

	ds_ptr->setBiallelic(!bgl_poly);
	ds_ptr->setPhased(bgl_phased);

	// read all the marker files first
	for(unsigned int i=0; i<chroms.size(); i++){
		currChrom = chroms[i];
		markerf.open(marker_fns[i].c_str());

		string markerline;
		while(getline(markerf, markerline)){
			boost::algorithm::trim(markerline);
			if(markerline.size() > 0 && markerline[0] != '#'){
				stringstream s(currChrom + " " + markerline);
				Marker* m = parseMap(s);
				if(m != reinterpret_cast<Marker*>(-1)){
					_marker_incl.push_back(m);
				}
			}
		}
		markerf.close();

	}

	map<string, Sample*> incl_ids;
	vector<Sample*> curr_samps;

	// now read all the genotype information
	Sample* samp = 0;
	bool affect = false;
	for(unsigned int i=0; i<chroms.size(); i++){
		genof.open(beagle_fns[i].c_str());
		curr_samps.clear();
		string line;
		int lineno = 0;
		while(getline(genof, line)){
			boost::algorithm::trim(line);
			if(line.size() > 0 && line[0] != '#'){

				++lineno;
				stringstream s(line);

				// On the 1st line, we MUST have the ID line
				if(lineno == 1){
					if(line[0] != 'I' && line[0] != 'i'){
						Logger::log_err("ERROR: PLATO requires the ID line in all BEAGLE files", true);
					}

					unsigned int fam_sz = 2 + bgl_trio;
					string id1, id2; // these better match!!
					s >> id1 >> id2; // 1st 2 are "I id", which we don't care about
					while( (s >> id1 >> id2) ){
						if(id1 != id2){
							Logger::log_err("ERROR: consecutive IDs do not match in " + beagle_fns[i], true);
						}

						if(i == 0){
							// Append either a Sample pointer or null pointer
							// to the list of curent samples
							curr_samps.push_back(filterSample(id1) ? ds_ptr->addSample(id1, _marker_incl.count()) : 0);
							// Map the id to the sample pointer (make sure to
							// account for the null pointer, which indicates
							// we excluded this sample from analysis!)
							incl_ids[id1] = curr_samps[curr_samps.size() - 1];

							if(bgl_pair || bgl_trio){
								if(curr_samps.size() % (fam_sz) == 0){
									Family* f = ds_ptr->addFamily(boost::lexical_cast<string>(curr_samps.size() / fam_sz));
									if(curr_samps[curr_samps.size() - fam_sz] && samp){
										curr_samps[curr_samps.size() - fam_sz]->addChild(samp);
										samp->addMother(curr_samps[curr_samps.size() - fam_sz]);
									}
									if(samp){
										f->addNonFounder(samp);
									}
									if(curr_samps[curr_samps.size() - fam_sz]){
										f->addFounder(curr_samps[curr_samps.size() - fam_sz]);
									}

									if(bgl_trio){
										if(curr_samps[curr_samps.size() - 2]){
											if(samp){
												samp->addFather(curr_samps[curr_samps.size() - 2]);
												curr_samps[curr_samps.size() - 2]->addChild(samp);
											}
											f->addFounder(curr_samps[curr_samps.size() - 2]);
										}
									}
								}
							} else if (samp) {
								ds_ptr->addFamily(id1)->addFounder(samp);
							}


						} else{
							map<string, Sample*>::const_iterator id_itr = incl_ids.find(id1);
							if(id_itr == incl_ids.end()){
								Logger::log_err("ERROR: ID '" + id1 + "' in file '"
												+ beagle_fns[i] + "' not found in first genotype file", true);
							} else{
								curr_samps.push_back((*id_itr).second);
							}

						}
					}
				}

				// load up the affection status
				if(i==0 && !affect && (line[0] == 'A' || line[0] == 'a')){
					string a1, a2;
					s >> a1 >> a2; // ignore the 1st 2 strings - we don't need them!
					for(unsigned int j=0; j<curr_samps.size(); j++){
						s >> a1 >> a2;
						if(a1 != a2){
							Logger::log_err("WARNING: discordant affection status, making unknown");
						} else if(curr_samps[j]){
							if(a1 == "1"){
								curr_samps[j]->setAffected(false);
							} else if(a1 == "2"){
								curr_samps[j]->setAffected(true);
							}
						}
					}

				} else if(i==0 && (line[0] == 'T' || line[0] == 't' || line[0] == 'A' || line[0] == 'a')){
					string miss_val = "";
					if(line[0] == 'a' || line[0] == 'A'){
						Logger::log_err("WARNING: multiple affection statuses seen, loading as a trait");
						miss_val = "0";
					}
					string trait_id, t1, t2;
					s >> trait_id;
					s >> trait_id;

					for(unsigned int j=0; j<curr_samps.size(); j++){
						s >> t1 >> t2;
						float tval = std::numeric_limits<float>::quiet_NaN();
						if(t1 != t2){
							Logger::log_err("WARNING: discordant trait values, making unknown");
						} else if (!(stringstream(t1) >> tval)){
							tval = std::numeric_limits<float>::quiet_NaN();
						}
						if(curr_samps[j]){
							ds_ptr->addTrait(trait_id,curr_samps[j],tval);
						}
					}

				} else if (line[0] == 'M' || line[0] == 'm'){
					string mid, g1, g2;
					s >> mid; // ignore the "M"
					s >> mid;

					Marker* curr_mark = ds_ptr->getMarker(mid);
					if(curr_mark != 0){
						for(unsigned int j=0; j<curr_samps.size(); j++){
							s >> g1 >> g2;
							if(curr_samps[j]){
								parseSample(curr_mark, curr_samps[j], g1, g2, false, bgl_phased);
							}
						}

					} else {
						Logger::log_err("WARNING: unknown marker found on line " +
								boost::lexical_cast<string>(lineno) + " of '" +
								beagle_fns[i] + "', ignoring");
					}
				} else if (i == 0 && lineno != 1 && !(line[0] == 'a' || line[0] == 'A' || line[0] == 'T' || line[0] == 't')){
					Logger::log_err("WARNING: unknown line identifier on line " +
							boost::lexical_cast<string>(lineno) + " of '" + beagle_fns[i] +
							"': PLATO only recognizes 'A', 'T', and 'M' lines");
				}
			}
		}
		genof.close();
	}

}

void DataLoader::readVCF(){

	Logger::log_err("INFO: The VCF loader in PLATO only uses some of the "
			"information in a VCF file: details may be lost.  "
			"Please keep the original VCF file for completeness.");

	ICompressedFile vcf_f(vcf_fn.c_str());

	// have we already warned about flags and VCF inconsistencies?
	bool poly_warn = vcf_poly;
	bool phase_warn = false;

	string geno_sep = vcf_phased ? "|" : "/";
	string alt_geno_sep = vcf_phased ? "/" : "|";

	ds_ptr->setBiallelic(!vcf_poly);
	ds_ptr->setPhased(vcf_phased);

	string curr_line;
	deque<Sample*> samp_list;
	vector<bool> incl_sample;

	unsigned int lineno = 0;
	unsigned int n_fields;
	while(getline(vcf_f, curr_line)){
		++lineno;
		vector<string> fields;
		if(curr_line.size() > 2 && curr_line[0] == '#' && curr_line[1] != '#'){
			// In this case, we are looking at the "#CHROM" line (we hope!)
			split(fields, curr_line, boost::is_any_of("\t"));
			n_fields = fields.size();

			// Sanity check here! Note that I never use the "QUAL" or "INFO"
			// fields, so I'm not going to check for them!
			if(n_fields < 9 || fields[0] != "#CHROM" || fields[1] != "POS" ||
					fields[2] != "ID" || fields[3] != "REF" || fields[4] != "ALT" ||
					fields[6] != "FILTER" || fields[8] != "FORMAT"){
				Logger::log_err("ERROR: VCF header line is malformed, please check your VCF input.", true);
			}

			// by default, include all samples
			incl_sample.resize(n_fields - 9, true);

			for(unsigned int i=9; i<n_fields; i++){
				if(filterSample(fields[i])){
					samp_list.push_back(ds_ptr->addSample(fields[i], 0));
				}else{
					incl_sample[i] = false;
				}
			}

		} else if(curr_line.size() > 2 && curr_line[0] != '#'){
			// In this case, we're looking at a marker

			split(fields, curr_line, boost::is_any_of("\t"));
			if(fields.size() != n_fields){
				Logger::log_err("ERROR: Mismatched number of fields on line "
						+ boost::lexical_cast<string>(lineno), true);
			}

			string chr, id, ref, alt, filter, format;
			unsigned int bploc = 0;
			int g1, g2;
			chr = fields[0];
			bploc = boost::lexical_cast<unsigned int>(fields[1]);
			id = fields[2];
			ref = fields[3];
			alt = fields[4];
			filter = fields[6];
			format = fields[8];

			// check for inclusion based on the marker-level data
			bool use = filterMarker(chr, bploc, id);
			use &= (no_filter_marker || filter == "." || filter == "PASS");

			if(use){
				Marker* m = ds_ptr->addMarker(chr,bploc,id);
				m->addAllele(ref);
				vector<string> alt_list;
				split(alt_list, alt, boost::is_any_of(","));
				if(!poly_warn && alt_list.size() > 1){
					Logger::log_err("WARNING: VCF appears polyallelic, but "
							"--vcf-poly was not given.  Alternate alleles will "
							"be collapsed to the first alternate.");
					poly_warn = true;
				}
				for(unsigned int i=0; i<alt_list.size(); i++){
					m->addAllele(alt_list[i]);
				}

				// parse the format string
				unsigned int gt_idx;
				unsigned int ft_idx;

				vector<string> format_list;
				split(format_list, format, boost::is_any_of(":"));
				gt_idx = (find(format_list.begin(), format_list.end(), "GT") - format_list.begin());
				ft_idx = (find(format_list.begin(), format_list.end(), "FT") - format_list.begin());
				ft_idx = (ft_idx == format_list.size()) ? static_cast<unsigned int>(-1) : ft_idx;

				if(gt_idx == format_list.size()){
					Logger::log_err("ERROR: No 'GT' format on line " +
							boost::lexical_cast<string>(lineno) +
							", cannot continue.", true);
				}

				deque<Sample*>::const_iterator si = samp_list.begin();
				for(unsigned int i=0; i<n_fields-9; i++){
					if(incl_sample[i]){
						// parse the individual call
						vector<string> geno_list;
						split(geno_list, fields[i+9], boost::is_any_of(":"));

						if(no_filter_geno || ft_idx == static_cast<unsigned int>(-1) ||
						   geno_list[ft_idx] == "PASS"){
							vector<string> call_list;
							split(call_list, geno_list[gt_idx], boost::is_any_of(geno_sep));
							if(call_list.size() == 1){
								split(call_list, geno_list[gt_idx], boost::is_any_of(alt_geno_sep));
								if(!phase_warn && call_list.size() != 1){
									Logger::log_err("WARNING: VCF phasing is "
											"inconsistent with arguments passed "
											"to the loader, genotype data may "
											"have unexpected results.");
								}
							}


							if(call_list.size() != 2){
								Logger::log_err("WARNING: Non-diploid found on line " +
										boost::lexical_cast<string>(lineno) + ", setting to missing");
								(*si)->appendMissingGenotype();
							} else {
								// If we can't convert to an integer, set to -1 (missing genotype)
								(*si)->appendGenotype(
									static_cast<unsigned char>((stringstream(call_list[0]) >> g1) ? g1 : -1),
									static_cast<unsigned char>((stringstream(call_list[1]) >> g2) ? g2 : -1));
							}
						} else {
							(*si)->appendMissingGenotype();
						}
						++si;
					}

				} // end iterating over genotypes

			} // end if(use)

		} // end if (header) else if (regular line)
	} // end while(getline)

	vcf_f.close();


}

Marker* DataLoader::parseMap(stringstream& ss) const {
	string chr;
	ss >> chr;

	// skip empty lines or those that begin with "#"
	if (chr.size() == 0 || chr[0] == '#') {
		return reinterpret_cast<Marker*>(-1);
	}
	string id;
	int bploc;
	float distance;
	string ref;
	string alt;

	ss >> id;
	if (!_map_no_distance) {
		ss >> distance;
	}

	static bool warn_bploc = false;
	if(! (ss >> bploc) ){
		if(!warn_bploc){
			Logger::log_err("WARNING: Parse error when reading base pair location; bploc forced to 1.  Potential 3-column map file (use --map3)?");
		}
		bploc = 1;
	}

	if (_map_ref) {
		ss >> ref;
	}
	if (_map_alt) {
		ss >> alt;
	}

	// check for conditions of use here!
	bool use = filterMarker(chr, bploc, id);

	Marker* m = 0;
	if (use) {
		m = ds_ptr->addMarker(chr, static_cast<unsigned int> (bploc),
				id);

		if (_map_ref) {
			m->addAllele(ref);
			//m->setRefAllele(ref);
		}
		if (_map_alt) {
			m->addAllele(alt);
			//m->setAltAllele(alt);
		}
		if(_map_others) {
			while(ss >> alt){
				m->addAllele(alt);
			}
		}
	}
	return m;
}

void DataLoader::parseSample(Marker* m, Sample* s, const string& g1, const string& g2, bool append, bool phased) const{

	if ((g1 == _ped_missing_geno) + (g2 == _ped_missing_geno) >= (1 + phased)) {
		if(append){
			s->appendMissingGenotype();
		} else {
			s->setMissingGenotype(*m);
		}
	} else {
		unsigned char g_idx1 = (g1 == _ped_missing_geno) ? Sample::missing_allele :
				static_cast<unsigned char> (m->addAllele(g1));
		unsigned char g_idx2 = (g2 == _ped_missing_geno) ? Sample::missing_allele :
				static_cast<unsigned char> (m->addAllele(g2));

		if(append){
			s->appendGenotype(g_idx1, g_idx2);
		} else {
			s->setGenotype(*m, g_idx1, g_idx2);
		}
	}

}

bool DataLoader::filterMarker(const string& chrom, int bploc, const string& id) const{
	// do not use if bploc <= 0
	bool toret = bploc > 0;

	// do not use if chromosome is not in the given chrom list
	toret = toret && (chrom_ids.empty() || chrom_ids.count(InputManager::chrStringToInt(chrom)));

	// do not use if not found in the given include list
	toret = toret && (incl_marker_set.empty() || incl_marker_set.count(id));

	// do not use if found in the given exclude list
	toret = toret && (excl_marker_set.empty() || !excl_marker_set.count(id));

	// do not use if it falls outside the given base pair range
	if(toret && !bp_window_set.empty()){
		// at this point, assume false - prove me wrong!
		toret = false;

		// get the element strictly larger than the current base pair
		set<pair<unsigned int, unsigned int> >::const_iterator bp_itr =
			bp_window_set.upper_bound(pair<unsigned int, unsigned int>(bploc, bploc));

		// subtract one (now, the 1st element of the pair is <= bploc)
		// and make sure that the second element is >= bploc
		if(bp_itr != bp_window_set.begin()){
			--bp_itr;
			if(static_cast<unsigned int>(bploc) <= (*bp_itr).second){
				toret = true;
			}
		}
	}

	return toret;
}

bool DataLoader::filterSample(const string& id, const string& fid) const{

	// Make sure we are in the include list
	bool toret = (incl_sample_set.empty() || incl_sample_set.count(id + sampl_field_sep + fid));

	// Make sure we aren't in the exclude list!
	toret = toret && (excl_sample_set.empty() || !excl_sample_set.count(id + sampl_field_sep + fid));

	return toret;
}

void DataLoader::readSampleList(const vector<string>& in_list, set<string>& out_set){
	vector<string> all_samps;
	InputManager::parseInput(in_list, all_samps);

	for(unsigned int i=0; i<all_samps.size(); i++){
		addSampleToSet(all_samps[i], out_set);
	}
}

void DataLoader::readSampleFile(const vector<string>& fn_list, set<string>& out_set){
	string samp;
	for(unsigned int i=0; i<fn_list.size(); i++){
		ifstream f(fn_list[i].c_str());
		while(getline(f, samp)){
			boost::algorithm::trim(samp);
			addSampleToSet(samp, out_set);
		}
	}
}

void DataLoader::readMarkerFile(const vector<string>& fn_list, set<string>& out_set){
	string marker_id;
	for(unsigned int i=0; i<fn_list.size(); i++){
		ifstream f(fn_list[i].c_str());
		while(f >> marker_id){
			out_set.insert(marker_id);
		}
	}
}

void DataLoader::addSampleToSet(const string& samp, set<string>& out_set){
	stringstream ss(samp);
	string id, fid;
	if(ss >> id){
		fid = id;
		if(ss >> fid){
			Logger::log_err("WARNING: Sample '" + samp +
					"' has more than an 'FID IID', "
					"ignoring everything after 2nd space!");
		}
	}
	out_set.insert(id + sampl_field_sep + fid);
}

}
}
