#include "GenotypeReader.h"

#include<vector>
#include<boost/algorithm/string.hpp>

using std::string;
using std::map;
using std::vector;

using boost::algorithm::split;
using boost::algorithm::is_any_of;

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

DataSet& GenotypeReader::readBim(const string& fn, DataSet& ds) const{


	map<string, vector<string> > mdescinfo;
	map<string,int> exclude;
	map<string,int> include;
	vector<string> mdescheaders;
	map<string,float> frequencies = options.getFrequencies();


	Logger& log = Logger::getLogger();
	log.log("Reading bim from " + fn);

	ifstream MAP(fn.c_str(), ios::in);
	if(!MAP){
		log.log_err("Error opening bim file: " + fn + ", exiting!");
		throw exception("Error opening bim file: " + fn);
	}


	string line;
	vector<string> split_line;

	while(getline(MAP, line)){
		split_line.clear();
		split(split_line, line, is_any_of(" \t"), boost::token_compress_on);

		if()

	}


	int count = 0;

		string chrom = "";
		string probe = "";
		int bploc = 0;
		string a1 = "";
		string a2 = "";
		string centi = "";
		string rsid = "";
		string enzyme = "";
		string line = "";
		string referent = "";
	int num_cols = 6 + options.getMapContainsReferent();
	while(getline(MAP, line)){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if((int)tokens.size() != num_cols){
			opts::printLog(".bim file column size != " + getString<int>(num_cols) + ": " + line + " stopping!!\n");
			throw MethodException(".bim file column size != " + getString<int>(num_cols) + ": " + line + " stopping!!\n");
		}
		chrom = tokens.at(0);
		probe = tokens.at(1);
		centi = tokens.at(2);
		bploc = atoi(tokens.at(3).c_str());
		a1 = tokens.at(4);
		a2 = tokens.at(5);
		if(options.getMapContainsReferent()){
			referent = tokens.at(6);
		}
		if(rsid == "."){
			rsid = "";
		}
		if(enzyme == "."){
			enzyme = "";
			opts::_ENZYMES_ = false;
		}
		bool use = true;

		Methods::Marker* m = new Marker(chrom, probe, bploc);
		if(opts::_AUTOONLY_ && ((m->getChrom() >= opts::_CHRX_) || (m->getChrom() < 1))){
			use = false;
		}

		m->setEnabled(use);
		m->setLoc(count);
		m->setAllele1(a1);
		m->setAllele2(a2);
		m->setReferent(referent);
		m->setEnzyme(enzyme);
		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe);
			if(found == frequencies.end()){
				found = frequencies.find(m->getRSID());
				if(found != frequencies.end()){
					m->setMAF(frequencies[m->getRSID()]);
					m->setFreqFlag(true);
				}
			}
			else{
				m->setMAF(frequencies[m->getRSID()]);
				m->setFreqFlag(true);
			}
		}
        if(opts::_MAPDESC_.length() > 0 && mdescinfo.size() > 0){
			vector<string> tokens = mdescinfo[probe];
            for(unsigned int i = 1; i < mdescheaders.size(); i++){
	            if(tokens.size() == mdescheaders.size()){
		            m->assignDetail(mdescheaders[i], tokens[i]);
                }
                else{
	                m->assignDetail(mdescheaders[i], "NA");
                }
            }
        }

		markers->push_back(m);

		count++;
	}
	MAP.close();
	opts::_MARKERS_FOUND_ = markers->size();

	if(filters != NULL){
		for(int f = 0; f < filters->num_locus_filters(); f++){
			filters->run_locus_filter(f, markers);
		}
	}


	marker_map->resize(markers->size());
	stable_sort(markers->begin(), markers->end(), less<Methods::Marker*>());

	for(unsigned int i = 0; i < markers->size(); i++){
		(*marker_map)[(*markers).at(i)->getLoc()] = i;
	}

	return ds;
}
}
