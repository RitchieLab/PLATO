/*
 * OutputPED.cpp
 *
 *  Created on: Nov 27, 2013
 *      Author: jrw32
 */

#include "OutputTPED.h"

#include <fstream>

#include "data/DataSet.h"
#include "data/Sample.h"
#include "data/Marker.h"

#include "util/InputManager.h"

namespace po=boost::program_options;

using PLATO::Data::DataSet;
using PLATO::Data::Sample;
using PLATO::Data::Marker;
using PLATO::Utility::InputManager;

using std::string;
using std::ofstream;

namespace PLATO{
namespace ProcessLib{

const string OutputTPED::stepname = OutputTPED::doRegister("output-tped");

po::options_description& OutputTPED::appendOptions(po::options_description& opts){
	po::options_description ped_opts("TPED Output Options");

	ped_opts.add_options()
		("file,f", po::value<string>(&base_fn)->default_value("plato"), "Base fileneme for transposed ped and map files")
		("tped", po::value<string>(&tped_fn), "Filename of TPED file")
		("tfam", po::value<string>(&tfam_fn), "Filename of TFAM file");

	return opts.add(ped_opts);
}

void OutputTPED::parseOptions(const po::variables_map& vm){
	if(tped_fn.size() == 0){
		tped_fn = base_fn + ".tped";
	}

	if(tfam_fn.size() == 0){
		tfam_fn = base_fn + ".tfam";
	}
}

void OutputTPED::process(DataSet& ds){

	ofstream tped_f(tped_fn.c_str());
	ofstream tfam_f(tfam_fn.c_str());

	if(!tped_f.good() || !tfam_f.good()){
		throw std::logic_error("Unable to open TPED/TFAM file(s)");
	}

	// Print the TFAM file first
	DataSet::const_sample_iterator s_itr = ds.beginSample();
	DataSet::const_sample_iterator s_end = ds.endSample();
	while(s_itr != s_end){
		printPEDHeader(tfam_f, **s_itr);
		tfam_f << std::endl;
		++s_itr;
	}

	tfam_f.close();

	// OK, now for the TPED file
	DataSet::const_marker_iterator m_itr = ds.beginMarker();
	DataSet::const_marker_iterator m_end = ds.endMarker();

	while(m_itr != m_end){

		DataSet::const_sample_iterator sm_itr = ds.beginSample();
		std::pair<unsigned char, unsigned char> geno;

		printMAPInfo(tped_f, **m_itr, false);

		while(sm_itr != s_end){
			if((*sm_itr)->isMissing(**m_itr)){
				tped_f << "\t0\t0";
			} else {
				geno = (*sm_itr)->getGeno(**m_itr);
				tped_f << "\t" << (*m_itr)->getAllele(geno.first)
					   << "\t" << (*m_itr)->getAllele(geno.second);
			}
			++sm_itr;
		}

		tped_f << std::endl;

		++m_itr;
	}

	tped_f.close();

}

}
}
