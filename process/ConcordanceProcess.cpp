/*
 * ConcordanceProcess.cpp
 *
 *  Created on: Nov 26, 2013
 *      Author: jrw32
 */

#include "ConcordanceProcess.h"

#include <map>
#include <fstream>
#include <iostream>

#include "data/DataSet.h"
#include "data/Marker.h"
#include "data/Sample.h"

#include "util/Logger.h"

using std::string;
using std::map;
using std::ofstream;
using std::endl;
using std::pair;

using PLATO::Utility::DataLoader;
using PLATO::Utility::Logger;
using PLATO::Data::DataSet;
using PLATO::Data::Sample;
using PLATO::Data::Marker;


namespace po=boost::program_options;

namespace PLATO{
namespace ProcessLib{

const std::string ConcordanceProcess::stepname = ConcordanceProcess::doRegister("concordance");

po::options_description& ConcordanceProcess::appendOptions(po::options_description& opts){
	po::options_description concord_opts("Concordance Options");

	concord_opts.add_options()
		("prefix", po::value<string>(&_prefix)->default_value("concordance"), "Prefix for concordance files")
		("error", po::value<string>(&_error_fn)->default_value("error"), "Suffix for concordance error file")
		("sample", po::value<string>(&_sample_fn)->default_value("sample"), "Suffix for sample concordance file")
		("marker", po::value<string>(&_marker_fn)->default_value("marker"), "Suffix for marker concordance file")
		("marker-mismatch", po::value<string>(&_mmismatch_fn)->default_value("marker-mismatch"), "Suffix for marker mismatch file")
		("sample-mismatch", po::value<string>(&_smismatch_fn)->default_value("sample-mismatch"), "Suffix for sample mismatch file")
		("extension", po::value<string>(&_ext)->default_value("txt"), "Extension for concordance files")
		("inc-missing", po::bool_switch(&_incl_missing), "Consider missing genotypes in the concordance checks")
		("sep", po::value<string>(&_sep)->default_value("\t"), "Separator for output files")
	;


	return DataLoader::addOptions(opts);
}

void ConcordanceProcess::parseOptions(const po::variables_map& vm){
	DataLoader::parseOptions(vm);

	if(_mmismatch_fn != ""){
		_mmismatch_fn = _prefix + "." + _mmismatch_fn + "." + _ext;
	}
	if(_smismatch_fn != ""){
		_smismatch_fn = _prefix + "." + _smismatch_fn + "." + _ext;
	}
	if(_error_fn != ""){
		_error_fn = _prefix + "." + _error_fn + "." + _ext;
	}
	if(_sample_fn != ""){
		_sample_fn = _prefix + "." + _sample_fn + "." + _ext;
	}
	if(_marker_fn != ""){
		_marker_fn = _prefix + "." + _marker_fn + "." + _ext;
	}
}

void ConcordanceProcess::process(DataSet& ds){
	DataSet alt_ds;
	read(alt_ds);

	// OK, now I can check the concordance between ds and alt_ds!

	// Check to make sure that the phasing agrees!
	if(ds.isPhased() != alt_ds.isPhased()){
		Logger::log_err("ERROR: Cannot compare datasets with mismatched phasing", true);
	}

	// these give the correspondence between markers and samples in either dataset
	map<const Sample*, const Sample*> sample_map;
	map<const Marker*, const Marker*> marker_map;

	// gives a map of samples -> (# errors, # checked)
	// note that we don't need a similar marker file b/c we'll print it on the fly
	map<const Sample*, pair<unsigned int, unsigned int> > sample_count;

	// First, find all of the matching markers and samples
	// let's sort all the markers and samples to make it easier on us
	ds.sortMarkers();
	alt_ds.sortMarkers();
	ds.sortSamples();
	alt_ds.sortSamples();

	ofstream* m_mismatch_f = 0;
	if(_mmismatch_fn != ""){
		m_mismatch_f = new ofstream(_mmismatch_fn.c_str());
		(*m_mismatch_f) << "Chr" << _sep << "ID" << _sep << "bploc" << _sep
				<< "in_base" << endl;
	}
	DataSet::const_marker_iterator mi = ds.beginMarker();
	DataSet::const_marker_iterator ami = alt_ds.beginMarker();
	DataSet::const_marker_iterator me = ds.endMarker();
	DataSet::const_marker_iterator ame = alt_ds.endMarker();

	while(mi != me && ami != ame){
		if(ami == ame || *mi < *ami){
			outputMarkerMismatch(m_mismatch_f, *mi, true);
			++mi;
		} else if (mi == me || *ami < *mi){
			outputMarkerMismatch(m_mismatch_f, *ami, false);
			++ami;
		} else {
			marker_map[*mi] = *ami;
			++mi;
			++ami;
		}
	}
	if(m_mismatch_f){
		m_mismatch_f->close();
		delete m_mismatch_f;
	}

	ofstream* s_mismatch_f = 0;
	if(_smismatch_fn != ""){
		s_mismatch_f = new ofstream(_smismatch_fn.c_str());
		(*s_mismatch_f) << "FID" << _sep << "IID" << _sep << "Sex" << _sep
				<< "Status" << _sep << "in_base" << endl;
	}
	DataSet::const_sample_iterator si = ds.beginSample();
	DataSet::const_sample_iterator asi = alt_ds.beginSample();
	DataSet::const_sample_iterator se = ds.endSample();
	DataSet::const_sample_iterator ase = alt_ds.endSample();

	while(si != se && asi != ase){
		if(asi == ase || *si < *asi){
			outputSampleMismatch(s_mismatch_f, *si, true);
			++si;
		} else if (mi == me || *asi < *si){
			outputSampleMismatch(s_mismatch_f, *asi, false);
			++asi;
		} else {
			sample_map[*si] = *asi;
			sample_count[*si] = std::make_pair(0,0);
			++si;
			++asi;
		}
	}
	if(s_mismatch_f){
		s_mismatch_f->close();
		delete s_mismatch_f;
	}

	// At this point, we have the markers and samples ready to go.  Let's start
	// a concordance check!
	ofstream* error_f = 0;
	ofstream* sample_f = 0;
	ofstream* marker_f = 0;

	if(_error_fn != ""){
		error_f = new ofstream(_error_fn.c_str());
		(*error_f) << "FID" << _sep << "IID" << _sep << "Chr" << _sep
				<< "rsID" << _sep << "bploc" << _sep << "Original_Genotype"
				<< _sep << "New_Genotype" << endl;
	}

	if(_sample_fn != ""){
		sample_f = new ofstream(_sample_fn.c_str());
		(*sample_f) << "FID" << _sep << "IID" << _sep << "Sex" << _sep
				<< "Status" << _sep << "Errors" << _sep << "Total_Compared"
				<< _sep << "%Concordance" << endl;
	}

	if(_marker_fn != ""){
		marker_f = new ofstream(_marker_fn.c_str());
		(*marker_f) << "Chr" << _sep << "rsID" << "bploc" << _sep << "Errors"
				<< _sep << "Total_Compared" << _sep << "%Concordance" << endl;
	}

	map<const Marker*, const Marker*>::const_iterator mmi = marker_map.begin();
	bool err = false;
	bool cmp = false;
	while(mmi != marker_map.end()){
		map<const Sample*, const Sample*>::const_iterator ssi = sample_map.begin();
		map<const Sample*, pair<unsigned int, unsigned int> >::iterator sci = sample_count.begin();
		unsigned int n_err = 0;
		unsigned int n_comp = 0;
		while(ssi != sample_map.end()){
			unsigned char res = match(error_f, (*mmi).first, (*ssi).first,
					(*mmi).second, (*ssi).second, ds.isPhased());

			cmp = _incl_missing || (res & 2);
			err = cmp && (res & 1);

			n_err += err;
			n_comp += cmp;
			(*sci).second.first += err;
			(*sci).second.second += cmp;

			++sci;
			++ssi;
		}

		if(marker_f){
			printMarker(*marker_f, *((*mmi).first));
			(*marker_f) << _sep << n_err << _sep << n_comp << _sep
					    << 100 * n_err / static_cast<float>(n_comp) << endl;
		}

		++mmi;
	}


	if(error_f){
		error_f->close();
		delete error_f;
	}

	if(marker_f){
		marker_f->close();
		delete marker_f;
	}

	if(sample_f){
		// Before we close, let's go ahead and output the sample concordance file
		map<const Sample*, pair<unsigned int, unsigned int> >::const_iterator sc_itr = sample_count.begin();
		while(sc_itr != sample_count.end()){

			printSample(*sample_f, *((*sc_itr).first));
			(*sample_f) << _sep
					    << (*sc_itr).second.first << _sep
					    << (*sc_itr).second.second << _sep
					    << 100 * (*sc_itr).second.first / static_cast<float>((*sc_itr).second.second) << endl;
		}

		sample_f->close();
		delete sample_f;
	}
}

unsigned char ConcordanceProcess::match(ofstream* error_f, const Marker* m, const Sample* s,
		   const Marker* a_m, const Sample* a_s, bool phased) const{

	pair<unsigned char, unsigned char> g1 = s->getGeno(*m);
	pair<unsigned char, unsigned char> g2 = a_s->getGeno(*a_m);

	bool err = false;
	bool missing = s->isMissing(*m) || a_s->isMissing(*a_m);

	const string& a11 = m->getAllele(g1.first);
	const string& a12 = m->getAllele(g1.second);
	const string& a21 = a_m->getAllele(g2.first);
	const string& a22 = a_m->getAllele(g2.second);

	if(phased){
		err = ((a11 != a21) || (a12 != a22));
	} else {
		if(a11 == a21){
			err = (a12 != a22);
		} else if (a11 == a22){
			err = (a12 != a21);
		} else {
			err = true;
		}
	}

	if(err && error_f){
		printSample(*error_f, *s, false);
		(*error_f) << _sep;
		printMarker(*error_f, *m);
		(*error_f) << _sep << a11 << (phased ? "|" : "/") << a12 << _sep
				<< a21 << (phased ? "|" : "/") << a22 << endl;
	}

	return (missing * 2) | (err * 1);
}

void ConcordanceProcess::outputMarkerMismatch(ofstream* f, const Marker* m, bool in_base) const {
	if(f){
		printMarker(*f, *m);
		(*f) << _sep << in_base << endl;
	}
}
void ConcordanceProcess::outputSampleMismatch(std::ofstream* f, const Data::Sample* s, bool in_base) const {
	if(f){
		printSample(*f, *s);
		(*f) << _sep << in_base << endl;
	}
}

void ConcordanceProcess::printMarker(ofstream& f, const Marker& m) const {
	f << m.getChromStr() << _sep << m.getID() << _sep << m.getLoc();
}

void ConcordanceProcess::printSample(ofstream& f, const Sample& s, bool extra_data) const{
	 f << s.getFID() << _sep << s.getID();

	 if(extra_data){
		f << _sep << (s.isGenderKnown() ? (s.isMale() ? "M" : "F") : "U") << _sep
		  << (s.isAffectedKnown() ? (s.isAffected() ? "Y" : "N") : "U");
	 }
}

}
}
