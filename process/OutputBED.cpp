/*
 * OutputBED.cpp
 *
 *  Created on: Dec 3, 2013
 *      Author: jrw32
 */

#include "OutputBED.h"

#include <fstream>

#include "util/DataSet.h"
#include "util/Sample.h"
#include "util/Marker.h"

namespace po=boost::program_options;

using std::string;
using std::ofstream;

using Methods::DataSet;
using Methods::Sample;
using Methods::Marker;

namespace ProcessLib {

const string OutputBED::stepname = OutputBED::doRegister("output-bed");

po::options_description& OutputBED::appendOptions(po::options_description& opts){
	po::options_description bed_opts("Binary PED output options");

	bed_opts.add_options()
		("file,f", po::value<string>(&base_fn)->default_value("plato"), "Base filename for BED/BIM/FAM files")
		("bed", po::value<string>(&bed_fn), "Filename for BED file")
		("bim", po::value<string>(&bim_fn), "Filename for BIM file")
		("fam", po::value<string>(&fam_fn), "Filename for FAM file")
		("individual-major", po::bool_switch(&_ind_major), "Output BED in individual-major mode");

	return opts.add(bed_opts);
}

void OutputBED::parseOptions(const po::variables_map& vm){
	if(bed_fn.size() == 0){
		bed_fn = base_fn + ".bed";
	}

	if(bim_fn.size() == 0){
		bim_fn = base_fn + ".bim";
	}

	if(fam_fn.size() == 0){
		fam_fn = base_fn + ".fam";
	}
}

void OutputBED::process(DataSet& ds){

	ofstream bed_f(bed_fn.c_str(), std::ios::out | std::ios::binary);
	ofstream bim_f(bim_fn.c_str());
	ofstream fam_f(fam_fn.c_str());

	if(!(bed_f.good() && bim_f.good() && fam_f.good())){
		throw std::logic_error("Could not open BED/BIM/FAM files");
	}

	// Let's print the bim file first
	DataSet::const_marker_iterator m_itr = ds.beginMarker();
	DataSet::const_marker_iterator m_end = ds.endMarker();

	while(m_itr != m_end){
		printMAPInfo(bim_f, **m_itr, true);
		bim_f << std::endl;
		++m_itr;
	}

	bim_f.close();

	// Now, for the fam file
	DataSet::const_sample_iterator s_itr = ds.beginSample();
	DataSet::const_sample_iterator s_end = ds.endSample();
	while(s_itr != s_end){
		printPEDHeader(fam_f, **s_itr);
		fam_f << std::endl;
		++s_itr;
	}

	fam_f.close();

	// yech! now for the BED file

	// We'll constantly use this unsigned char and write a byte at a time
	char ch = 0;

	//First, print the magic number
	ch = 108;
	bed_f.write(&ch, 1);
	ch = 27;
	bed_f.write(&ch, 1);

	// Now, are we individual-major or snp-major?
	ch = !_ind_major;
	bed_f.write(&ch, 1);

	ch = 0;

	if(_ind_major){
		s_itr = ds.beginSample();
		while(s_itr != s_end){

			unsigned int n_mark=4;
			m_itr = ds.beginMarker();

			while(m_itr != m_end){

				ch |= (getBinaryGeno(**s_itr, **m_itr) << ((4 - n_mark)*2)) ;

				if(--n_mark == 0){
					bed_f.write(&ch, 1);
					ch = 0;
					n_mark = 4;
				}

				++m_itr;
			}

			if(n_mark != 4){
				bed_f.write(&ch, 1);
			}

			++s_itr;
		}
	} else {

		m_itr = ds.beginMarker();
		while (m_itr != m_end) {

			unsigned int n_samp = 4;
			s_itr = ds.beginSample();
			while (s_itr != s_end) {

				ch |= (getBinaryGeno(**s_itr, **m_itr) << ((4 - n_samp)*2)) ;

				if (--n_samp == 0) {
					bed_f.write(&ch, 1);
					// probably unnecessary, but I like it
					ch = 0;
					n_samp = 4;
				}

				++s_itr;
			}

			if (n_samp != 4) {
				bed_f.write(&ch, 1);
			}

			++m_itr;
		}
	}

	bed_f.close();

}

unsigned char OutputBED::getBinaryGeno(const Sample& s, const Marker& m) const {

	unsigned char g_char = s.getAdditiveGeno(m);

	switch(g_char){
	case 0:
		// homozygous referent (00)
		g_char = 0;
		break;
	case 1:
		// heterozygous (10)
		g_char = 2;
		break;
	case 2:
		// homozygous alternate (11)
		g_char = 3;
		break;
	default:
		// missing (01)
		g_char = 1;
	}

	return g_char;

}
}
