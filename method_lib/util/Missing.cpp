#include "Missing.h"

#include "data/DataSet.h"
#include "data/Marker.h"
#include "data/Sample.h"

using PLATO::Data::DataSet;
using PLATO::Data::Marker;
using PLATO::Data::Sample;

namespace PLATO{
namespace Utility{

double Missing::sampleMissing(const DataSet& ds, const Sample& s){
	int n_marker;
	int n_missing;

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
	int n_sample;
	int n_missing;

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

}
}
