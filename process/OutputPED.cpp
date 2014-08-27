/*
 * OutputPED.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: jrw32
 */

#include "OutputPED.h"

#include <fstream>

#include "data/DataSet.h"
#include "data/Sample.h"
#include "data/Marker.h"

#include "util/InputManager.h"
#include "util/Logger.h"

namespace po=boost::program_options;

using PLATO::Data::DataSet;
using PLATO::Data::Sample;
using PLATO::Data::Marker;
using PLATO::Utility::InputManager;
using PLATO::Utility::Logger;

using std::string;
using std::ofstream;

namespace PLATO{
namespace ProcessLib{

const string OutputPED::stepname = OutputPED::doRegister("output-ped");

po::options_description& OutputPED::appendOptions(po::options_description& opts){
	po::options_description ped_opts("PED Outuput Options");

	ped_opts.add_options()
		("file,f", po::value<string>(&base_fn)->default_value("plato"), "Base fileneme for ped and map files")
		("ped", po::value<string>(&ped_fn), "Filename of PED file")
		("map", po::value<string>(&map_fn), "Filename of MAP file");

	return opts.add(ped_opts);
}

void OutputPED::parseOptions(const po::variables_map& vm){
	if(ped_fn.size() == 0){
		ped_fn = base_fn + ".ped";
	}

	if(map_fn.size() == 0){
		map_fn = base_fn + ".map";
	}
}

void OutputPED::process(DataSet& ds){

	ofstream ped_f(ped_fn.c_str());
	ofstream map_f(map_fn.c_str());

	if(!ped_f.good() || !map_f.good()){
		Logger::log_err("Unable to open PED/MAP file(s)", true);
	}

	// Let's print the map file first
	DataSet::const_marker_iterator m_itr = ds.beginMarker();
	DataSet::const_marker_iterator m_end = ds.endMarker();

	while(m_itr != m_end){
		printMAPInfo(map_f, **m_itr, false);
		map_f << std::endl;
		++m_itr;
	}

	map_f.close();

	// OK, now time to print out the PED file
	DataSet::const_sample_iterator s_itr = ds.beginSample();
	DataSet::const_sample_iterator s_end = ds.endSample();
	while(s_itr != s_end){
		printPEDHeader(ped_f, **s_itr);

		DataSet::const_marker_iterator ms_itr = ds.beginMarker();
		std::pair<unsigned char, unsigned char> geno;
		while(ms_itr != m_end){
			if((*s_itr)->isMissing(**ms_itr)){
				ped_f << "\t0\t0";
			} else {
				geno = (*s_itr)->getGeno(**ms_itr);
				ped_f << "\t" << (*ms_itr)->getAllele(geno.first)
					  << "\t" << (*ms_itr)->getAllele(geno.second);
			}
			++ms_itr;
		}

		ped_f << std::endl;
		++s_itr;
	}

	ped_f.close();

}

}
}
