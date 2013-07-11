#include "GenotypeReader.h"

#include<vector>

using std::string;
using std::map;
using std::vector;

namespace po = boost::program_options;
using po::value;

GenotypeReader::GenotypeReader(const string& prefix, const po::variables_map& vm) : _prefix(prefix){

}

GenotypeReader::~GenotypeReader(){

}

po::options_description& GenotypeReader::addOptions(po::options_description& opts){

	po::options_description input("General Input Options");

	input.add_options()
			("map", value<string>(),
					"Specify .map file")
			("fam", value<string>(),
					"Specify .fam file")
			("bim", value<string>(),
					"Spefify .bim file");

	po::options_description filters("Filtering Options");

	filters.add_options()
			("chr", value<unsigned short>(&_chrom)->default_value(0),
					"Select a particluar chromosome"),
			("from-bp", value<unsigned int>(&_start)->default_value(0),
					"Select markers starting at this base pair position")
			("to-bp", value<unsigned int>(&_stop)->default_value(static_cast<unsigned int>(-1)),
					"Select markers ending at this base pair position")
			("snps", value<string>(),
					"Select comma-delimited list of SNPs, allowing for ranges")
			("window", value<unsigned int>()->default_value(0),
					"Select a window around given SNPs")
			("snp-exc", value<vector<string> >()->composing(),
					"Comma-separated list of markers to exclude")
			("snp-exc-file", value<vector<string> >()->composing(),
					"Files containing markers to exclude")
			("ind-exc", value<vector<string> >()->composing(),
					"Comma-separated list of individuals to exlude")
			("ind-exc-file", value<vector<string> >()->composing(),
					"Files containing individuals to exclude")
			("fam-exc", value<vector<string> >()->composing(),
					"Comma-separated list of families to exclude")
			("fam-exc-file", value<vector<string> >()->composing(),
					"Files containing families to exclude")
			("snp-inc", value<vector<string> >()->composing(),
					"Comma-separated list of markers to include")
			("snp-inc-file", value<vector<string> >()->composing(),
					"Files containing markers to include")
			("ind-inc", value<vector<string> >()->composing(),
					"Comma-separated list of individuals to include")
			("ind-inc-file", value<vector<string> >()->composing(),
					"Files containing individuals to include")
			("fam-inc", value<vector<string> >()->composing(),
					"Comma-separated list of families to include")
			("fam-inc-file", value<vector<string> >()->composing(),
					"Files containing families to include")
			;

	return desc;
}
