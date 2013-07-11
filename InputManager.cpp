#include "InputManager.h"

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

namespace po=boost::program_options;

using std::string;
using std::map;

static InputManager::IChromMap InputManager::s_chr_map;
static map<unsigned short, string> InputManager::s_chrint_map;
static const string InputManager::s_missing_chr_str("");

po::options_description& InputManager::addOptions(po::options_description& opts){
	po::options_description global("Chromosome Mapping Options");
	global.add_options()
			("chroms", value<unsigned short>()->default_value(22),
					"Number of autosomal chromosomes")
			("extra-chroms", value<string>()->default_value("X,Y,XY,M-MT"),
					"Comma-separated list of extra chromosomes to include.  For each chromosome, a dash separates aliases.");

	opts.add(global);

	return opts;
}

void InputManager::parseGlobalOptions(const po::variables_map& vm){

	unsigned short curr_chr = 0;
	unsigned short max_chrom = vm["chroms"].as<unsigned short>();

	s_chr_map.clear();
	// Add chromosome numbers to the map
	while(curr_chr < max_chrom){
		++curr_chr;
		s_chr_map[boost::lexical_cast<string>(curr_chr)] = curr_chr;
	}

	// Add the extra chromosomes
	string chr_sep = ",";
	string alias_sep = "-";
	boost::tokenizer < boost::char_separator<char> > tok(
			vm["extra-chroms"].as<string> (), sep);
	for (boost::tokenizer<boost::char_separator<char> >::iterator beg =
			tok.begin(); beg != tok.end(); ++beg) {

		++curr_chr;

		boost::tokenizer < boost::char_separator<char> > alias_tok(
				*beg, alias_sep);
		for (boost::tokenizer<boost::char_separator<char> >::iterator a_beg =
				aliase_tok.begin(); a_beg != alias_tok.end(); ++alias_beg) {
			s_chr_map[*a_beg] = curr_chr;
		}

	}
}
