#include "Missing.h"

#include <cmath>

#include "data/DataSet.h"
#include "data/Marker.h"
#include "data/Sample.h"

using std::string;

using PLATO::Data::DataSet;
using PLATO::Data::Marker;
using PLATO::Data::Sample;

namespace PLATO{
namespace Utility{

double Missing::sampleMissing(const DataSet& ds, const Sample& s){
	unsigned int n_marker = 0;
	unsigned int n_missing = 0;

	DataSet::const_marker_iterator mi = ds.beginMarker();
	while(mi != ds.endMarker()){
		if(s.getAdditiveGeno(**mi) == Sample::missing_allele){
			++n_missing;
		}
		++mi;
		++n_marker;
	}

	return n_missing / static_cast<double>(n_marker);

}

double Missing::markerMissing(const DataSet& ds, const Marker& m){
	unsigned int n_sample = 0;
	unsigned int n_missing = 0;

	DataSet::const_sample_iterator si = ds.beginSample();
	while(si != ds.endSample()){
		if((*si)->getAdditiveGeno(m) == Sample::missing_allele){
			++n_missing;
		}
		++n_sample;
		++si;
	}

	return n_missing / static_cast<double>(n_sample);

}

double Missing::traitMissing(const DataSet& ds, const string& t){
	unsigned int n_sample = 0;
	unsigned int n_missing = 0;

	DataSet::const_sample_iterator si = ds.beginSample();
	while(si != ds.endSample()){
		if(std::isnan(ds.getTrait(t, *si))){
			++n_missing;
		}
		++n_sample;
		++si;
	}

	return n_missing / static_cast<double>(n_sample);

}

}
}
