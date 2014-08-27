/*
 * BatchProcess.cpp
 *
 *  Created on: Nov 26, 2013
 *      Author: jrw32
 */

#include "BatchProcess.h"

#include <fstream>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "data/DataSet.h"

namespace po=boost::program_options;

using std::string;
using std::ifstream;
using std::vector;

using PLATO::Data::DataSet;

namespace PLATO{
namespace ProcessLib{

const std::string BatchProcess::stepname = BatchProcess::doRegister("batch");

po::options_description& BatchProcess::appendOptions(po::options_description& opts){
	po::options_description batch_opts("Batch Processing Options");

	batch_opts.add_options()
		("file,f", po::value<string>(&batch_fn)->required(),"Batch filename");

	return opts.add(batch_opts);
}

void BatchProcess::process(DataSet& ds){
	// open up that batch file...
	ifstream input(batch_fn.c_str());

    if(!input.is_open()){
    	throw std::invalid_argument("Error opening batch file: " + batch_fn);
    }

	string line;

	ProcessFactory& fact = ProcessFactory::getFactory();
	vector<Process*> process_list;

    while(getline(input, line)){
    	vector<string> cmd_args;
    	boost::split(cmd_args,line,boost::is_any_of(" \t"));

    	// ignore empty lines and lines that begin with "#"
    	if(cmd_args.size() > 0 && cmd_args[0].size() > 0 && cmd_args[0][0] != '#'){
    		string cmd = cmd_args[0];
    		cmd_args.erase(cmd_args.begin());

    		Process* p = fact.Create(cmd);
    		if(p == 0){
    			throw std::invalid_argument("Error: unrecognized command: " + cmd);
    		}
    		po::options_description subopts;
			p->addOptions(subopts);

			// now, parse the remaining options, stopping at the name of the next command
			po::variables_map subvm;
			po::parsed_options subparsed = po::command_line_parser(cmd_args).options(subopts).run();
			po::store(subparsed, subvm);
			po::notify(subvm);

			p->parseOptions(subvm);

    		process_list.push_back(p);
    	}
    }

    // Now that we have all of our processes ready to go, let's run them!
    for(unsigned int i=0; i<process_list.size(); i++){
    	process_list[i]->run(ds);
    }

    // OK, now delete all of our processes
    for(unsigned int i=0; i<process_list.size(); i++){
    	delete process_list[i];
    }
    process_list.clear();


}

}
}
