/*
 * OutputPED.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: jrw32
 */

#include "OutputPED.h"

#include <fstream>

#include "util/DataSet.h"
#include "util/Sample.h"
#include "util/Marker.h"

#include "InputManager.h"

namespace po=boost::program_options;

using Methods::DataSet;
using Methods::Sample;
using Methods::Marker;
using std::string;
using std::ofstream;

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
		throw std::logic_error("Unable to open PED/MAP file(s)");
	}

	// Let's print the map file first
	DataSet::const_marker_iterator m_itr = ds.beginMarker();
	DataSet::const_marker_iterator m_end = ds.endMarker();

	while(m_itr != m_end){

		int myChrom = (*m_itr)->getChrom();

		map_f << InputManager::chrIntToString((*m_itr)->getChrom()) << "\t"
			  << (*m_itr)->getID() << "\t" << 0 << "\t"
			  << (*m_itr)->getLoc() << std::endl;

		++m_itr;
	}

	map_f.close();

	// OK, now time to print out the PED file
	DataSet::const_sample_iterator s_itr = ds.beginSample();
	DataSet::const_sample_iterator s_end = ds.endSample();
	while(s_itr != s_end){
		Sample* s = *s_itr;

		// printing information about the sample
		ped_f << s->getFID() << "\t" << s->getID() << "\t";
		ped_f << (s->getFather() == 0 ? "0" : s->getFather()->getID()) << "\t";
		ped_f << (s->getMother() == 0 ? "0" : s->getMother()->getID()) << "\t";
		ped_f << (s->isGenderKnown() ? (s->isFemale() + 1) : 0) << "\t";
		ped_f << (s->isAffectedKnown() ? (s->isAffected() + 1) : -9);

		DataSet::const_marker_iterator ms_itr = ds.beginMarker();
		std::pair<unsigned char, unsigned char> geno;
		string allele;
		while(ms_itr != m_end){
			geno = s->getGeno((*ms_itr)->getIndex());
			allele = (*ms_itr)->getAllele(geno.first);
			ped_f << "\t" << (allele == Marker::getMissingAllele() ? "0" : allele);
			allele = (*ms_itr)->getAllele(geno.second);
			ped_f << "\t" << (allele == Marker::getMissingAllele() ? "0" : allele);
			++ms_itr;
		}

		ped_f << std::endl;
		++s_itr;
	}

	ped_f.close();

}

}
