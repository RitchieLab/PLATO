/***************************************************************
 *		wasp - Genome Association Study Pipeline
 * Written by: Justin Giles
 *	      Vanderbilt University
 *	      Center for Human Genetics Research
 *
 *wasp.cc - Main class for wasp application
 *
 * See README file for complete overview.
 ***************************************************************/

#include "main.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "ProcessFactory.h"
#include "Process.h"
#include "InputManager.h"

#include "util/CommandLineParser.h"
#include "util/DataSet.h"

namespace po = boost::program_options;

using std::cout;
using std::string;
using std::vector;

using Methods::DataSet;

/*
 *Function: Initialize()
 *Return: void
 *Description:
 *Initializes the batch file options enumeration.
 */
void Initialize() {
	srand((unsigned) time(0));
}

int main(int argc, char** argv) {

	string logfn;

	po::options_description cmd_opts("Command Line Options");
	cmd_opts.add_options()
		("help,h","Display this help message")
		("list-command,L", "List the available commands")
		("help-command,C", po::value<string>(), "Get help for a particular command")
		("logfile,f", po::value<string>(&logfn)->default_value("plato.log"),"Name of the log file of PLATO output")
		("version,v", "Print the version of PLATO and exit");

	InputManager::addOptions(cmd_opts);

	po::variables_map vm;
	po::parsed_options parsed = CommandLineParser(argc, argv).options(cmd_opts).run();
	po::store(parsed, vm);
	po::notify(vm);

	InputManager::parseGlobalOptions(vm);

	// OK, check for a few things now:
	if(vm.count("help")){
		cout << cmd_opts;
		return 0;
	}

	if(vm.count("list-command")){
		print_steps();
		return 0;
	}

	if(vm.count("help-command")){
		ProcessFactory& f = ProcessFactory::getFactory();
		string cmd_str = vm["help-command"].as<string>();
		Process* p = f.Create(cmd_str);
		if(!p){
			cout << "ERROR: Unrecognized command '" << cmd_str << "'\n\n";
			cout << "Please use the 'list-command' option to see a list of commands available\n";
			return 1;
		}else{
			po::options_description opts;
			p->addOptions(opts);
			p->printHelp(cout);
			delete p;
			return 0;
		}
	}

	if(vm.count("version")){
		std::cout << PACKAGE_STRING << "\n";
		std::cout << "(c) Ritchie Lab, 2013\n";
		std::cout << "To report bugs, please email " << PACKAGE_BUGREPORT << "\n";
		return 0;
	}

	// TODO: set up the logger here

	vector<string> unrec_opt = po::collect_unrecognized(parsed.options, po::include_positional);
	vector<Process*> process_list;

	// OK, now we create the list of processes that we want to create
	while(unrec_opt.size() > 0){
		ProcessFactory& f = ProcessFactory::getFactory();

		string cmd_str = unrec_opt[0];

		Process* p = f.Create(cmd_str);
		if(!p){
			cout << "ERROR: Unrecognized command '" << cmd_str << "'\n\n";
			cout << "Please use the 'list-command' option to see a list of commands available\n";
			return 1;
		}else{
			// set up the options specific to the process
			po::options_description subopts;
			p->addOptions(subopts);

			// pop off the first string, which is the name of the command
			unrec_opt.erase(unrec_opt.begin());

			// now, parse the remaining options, stopping at the name of the next command
			po::variables_map subvm;
			po::parsed_options subparsed = CommandLineParser(unrec_opt).options(subopts).run();
			po::store(subparsed, subvm);
			po::notify(subvm);

			p->parseOptions(subvm);

			unrec_opt = po::collect_unrecognized(subparsed.options, po::include_positional);
		}

		process_list.push_back(p);

	}

	// Create a global DataSet for our use
	DataSet global_ds;

	// Go through each step and run on the global dataset
	for(unsigned int i=0; i<process_list.size(); i++){
		process_list[i]->run(global_ds);
	}

	// OK, now delete all of our processes
	for(unsigned int i=0; i<process_list.size(); i++){
		delete process_list[i];
	}
	process_list.clear();

	// And we're done! (dataset will clean up on its own)
	return 0;
}

/*
 *Function: print_steps
 *Return: void
 *Parameters: STEPS object
 *Description:
 *Outputs valid batch file steps and their descriptions when the -S command line argument is used
 */
void print_steps() {
	unsigned int field = 0;

	ProcessFactory& f = ProcessFactory::getFactory();

	for (ProcessFactory::const_iterator s_iter = f.begin(); s_iter != f.end(); s_iter++) {
		if (s_iter->first.size() > field) {
			field = s_iter->first.size();
		}
	}
	field += 5;
	cout << std::left << std::setw(field) << "Step:" << "Description:" << std::endl;
	cout << std::left << std::setw(field) << "----------" << "---------------" << std::endl;

	for (ProcessFactory::const_iterator s_iter = f.begin(); s_iter != f.end(); s_iter++) {

		Process* p = f.Create((*s_iter).first);
		cout << std::left << std::setw(field) << p->getName() << p->getDesc() << std::endl;
		delete p;
		//cout << "\tThresh: " << mystep.getThreshold() << endl;
	}
}

/*
 *Function: flipStrand
 *Description:
 *Flips the strand of the specified markers
 *
 */

/*
void flipStrand(vector<Marker*>* markers) {
	if (opts::_FLIPFILE_.length() > 0) {
		map<string, int> exclude;
		opts::printLog("Flipping markers found in: " + opts::_FLIPFILE_ + "\n");
		ifstream einput;
		einput.open(opts::_FLIPFILE_.c_str(), ios::in);
		if (!einput) {
			opts::printLog("Error opening marker strand flip file: "
					+ opts::_FLIPFILE_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";

		while (getline(einput, line)) {
			vector<string> tokens = General::ParseDelimitedLine(line);
			if (tokens.size() != 1) {
				opts::printLog("Strand flip markers column size != 1: " + line
						+ " Skipping line!!\n");
				continue;
			}
			exclude[tokens[0]] = 1;
			//exclude.push_back(probe);
		}
		einput.close();

		int msize = markers->size();

		for (int m = 0; m < msize; m++) {
			Marker* mark = (*markers)[m];
			if (!mark->isMicroSat()) {
				map<string, int>::iterator found = exclude.find(
						mark->getProbeID());
				if (found != exclude.end()) {
					string a1 = mark->getAllele1();
					string a2 = mark->getAllele2();

					if (a1 == "A") {
						mark->resetAllele1("T");
					} else if (a1 == "C") {
						mark->resetAllele1("G");
					} else if (a1 == "G") {
						mark->resetAllele1("C");
					} else if (a1 == "T") {
						mark->resetAllele1("A");
					} else if (a1 == "1") {
						mark->resetAllele1("4");
					} else if (a1 == "2") {

static void Initialize();
						mark->resetAllele1("3");
					} else if (a1 == "3") {
						mark->resetAllele1("2");
					} else if (a1 == "4") {
						mark->resetAllele1("1");
					}

					if (a2 == "A") {
						mark->resetAllele2("T");
					} else if (a2 == "C") {
						mark->resetAllele2("G");
					} else if (a2 == "G") {
						mark->resetAllele2("C");
					} else if (a2 == "T") {
						mark->resetAllele2("A");
					} else if (a2 == "1") {
						mark->resetAllele2("4");
					} else if (a2 == "2") {
						mark->resetAllele2("3");
					} else if (a2 == "3") {
						mark->resetAllele2("2");
					} else if (a2 == "4") {
						mark->resetAllele2("1");
					}
				}
			}
		}
	}

}

*/

/*
 *Function: printFamilies
 *Description:
 *outputs the structure of all families
 *
 */
/*
void printFamilies(vector<Family*>* families) {
	int fsize = families->size();

	string fname = opts::_OUTPREFIX_ + "family_structure.txt";
	ofstream fout(fname.c_str());
	if (!fout.is_open()) {
		opts::printLog("Unable to open " + fname + "\n");
		exit(1);
	}

	fname = opts::_OUTPREFIX_ + "family_data.txt";
	ofstream fdout(fname.c_str());
	if (!fdout.is_open()) {
		opts::printLog("Unable to open " + fname + "\n");
		exit(1);
	}

	for (int f = 0; f < fsize; f++) {
		Family* fam = (*families)[f];
		if (fam->isEnabled()) {
			fout << "Family: " << fam->getFamID() << endl;
			fout
					<< "----------------------------------------------------------"
					<< endl;
			vector<Sample*>* founders = fam->getFounders();
			int size = founders->size();
			if (size > 0) {
				for (int fs = 0; fs < size; fs++) {
					Sample* founder = (*founders)[fs];
					if (founder->isEnabled() || (founder->isExcluded()
							&& opts::_KEEP_EXC_SAMPLES_)) {
						fout << "Founder: " << founder->getInd() << endl;
						fout << descendTree(founder, 0) << endl;
					}
				}
			} else {
				fout << "No Founders...\n";
				vector<Sample*>* samps = fam->getSamples();
				int ssize = samps->size();
				for (int ss = 0; ss < ssize; ss++) {
					Sample* samp = (*samps)[ss];
					if (samp->isEnabled() || (samp->isExcluded()
							&& opts::_KEEP_EXC_SAMPLES_)) {
						fout << descendTree(samp, 0) << endl;
					}
				}
			}
			fout << endl << endl;
		}
	}

	fdout
			<< "FamID\tMultigenerational?\tMutigen_with_affecteds\tTotal_Generations\tTotal_Num_Affected\tGeneration1-N_Affected Counts\n";
	//output summary levels & affecteds
	for (int f = 0; f < fsize; f++) {
		Family* fam = (*families)[f];
		if (fam->isEnabled()) {
			map<int, vector<Sample*> > levels;
			fdout << fam->getFamID();
			vector<Sample*>* founders = fam->getFounders();
			int size = founders->size();
			if (size > 0) {
				//vector<int> levels(size);
				for (int fs = 0; fs < size; fs++) {
					Sample* founder = (*founders)[fs];
					if (founder->isEnabled() || (founder->isExcluded()
							&& opts::_KEEP_EXC_SAMPLES_)) {
						//levels[fs] = descendTree3(founder, 0);
						map<int, vector<Sample*> > levelstemp = descendTree3(
								founder, 1);
						//if(founder->getPheno() == 2){// && !founder->isFlagged()){
						map<int, vector<Sample*> >::iterator titer;
						for (titer = levelstemp.begin(); titer
								!= levelstemp.end(); titer++) {
							vector<Sample*> mysamps = titer->second;
							for (unsigned int ms = 0; ms < mysamps.size(); ms++) {
								levels[titer->first].push_back(mysamps[ms]);
							}
						}
						levels[1].push_back(founder);
						//							founder->setFlag(true);
						//}
					}
				}
				//sort(levels.begin(), levels.end());
				//fdout << "\t" << levels[levels.size() - 1];
			} else {
				vector<Sample*>* samps = fam->getSamples();
				int ssize = samps->size();
				//vector<int> levels(ssize);
				//map<int, int> levels;
				for (int ss = 0; ss < ssize; ss++) {
					Sample* samp = (*samps)[ss];
					if (samp->isEnabled() || (samp->isExcluded()
							&& opts::_KEEP_EXC_SAMPLES_)) {
						//levels[ss] = descendTree3(samp, 0);
						map<int, vector<Sample*> > levelstemp = descendTree3(
								samp, 1);
						//if(samp->getPheno() == 2){// && !samp->isFlagged()){
						map<int, vector<Sample*> >::iterator titer;
						for (titer = levelstemp.begin(); titer
								!= levelstemp.end(); titer++) {
							vector<Sample*> mysamps = titer->second;
							for (unsigned int ms = 0; ms < mysamps.size(); ms++) {
								levels[titer->first].push_back(mysamps[ms]);
							}
						}
						levels[1].push_back(samp);
						//							samp->setFlag(true);
						//}
					}
				}

				//sort(levels.begin(), levels.end());
				//fdout << "\t" << levels[levels.size() - 1];
			}
			//			vector<Sample*>* fsamps = fam->getSamples();
			int affected = 0;
			//			for(int s = 0; s < fsamps->size(); s++){
			//				if((*fsamps)[s]->getPheno() == 2 && (*fsamps)[s]->isFlagged()){
			//affected++;
			//					(*fsamps)[s]->setFlag(false);
			//				}
			//			}
			map<int, vector<Sample*> >::iterator iter = levels.end();
			if (levels.size() > 1) {
				fdout << "\tY";
			} else {
				fdout << "\tN";
			}
			for (iter = levels.begin(); iter != levels.end(); iter++) {
				vector<Sample*> mysamps = iter->second;
				map<Sample*, bool> newsamps;
				for (unsigned int ms = 0; ms < mysamps.size(); ms++) {
					newsamps[mysamps[ms]] = mysamps[ms]->getAffected();
				}
				map<Sample*, bool>::iterator niter;
				for (niter = newsamps.begin(); niter != newsamps.end(); niter++) {

					if (niter->second) {
						bool useme = true;
						map<int, vector<Sample*> >::iterator tempiter;
						for (tempiter = levels.begin(); tempiter
								!= levels.end(); tempiter++) {
							if (tempiter->first > iter->first) {
								vector<Sample*> findsamples = tempiter->second;
								vector<Sample*>::iterator findme = find(
										findsamples.begin(), findsamples.end(),
										niter->first);
								if (findme != findsamples.end()) {
									useme = false;
									break;
								}
							}
						}
						if (useme)
							affected++;
					}
				}
			}
			string affcounts = "";
			int numlevelsofaff = 0;
			for (iter = levels.begin(); iter != levels.end(); iter++) {
				int laffected = 0;
				vector<Sample*> mysamps = iter->second;
				map<Sample*, bool> newsamps;
				for (unsigned int ms = 0; ms < mysamps.size(); ms++) {
					newsamps[mysamps[ms]] = mysamps[ms]->getAffected();
				}
				bool incme = false;
				map<Sample*, bool>::iterator niter;
				for (niter = newsamps.begin(); niter != newsamps.end(); niter++) {
					if (niter->second) {
						bool useme = true;
						map<int, vector<Sample*> >::iterator tempiter;
						for (tempiter = levels.begin(); tempiter
								!= levels.end(); tempiter++) {
							if (tempiter->first > iter->first) {
								vector<Sample*> findsamples = tempiter->second;
								vector<Sample*>::iterator findme = find(
										findsamples.begin(), findsamples.end(),
										niter->first);
								if (findme != findsamples.end()) {
									useme = false;
									break;
								}
							}
						}
						if (useme) {
							laffected++;
							if (!incme) {
								numlevelsofaff++;
								incme = true;
							}
						}
					}
				}
				affcounts += "\t" + getString<int> (laffected);
			}
			if (numlevelsofaff > 1) {
				fdout << "\tY";
			} else {
				fdout << "\tN";
			}
			fdout << "\t" << levels.size();
			fdout << "\t" << affected;
			fdout << affcounts;
			fdout << endl;
		}
	}

	fout.close();
	fdout.close();

}
*/

/*
 *Function: descendTree
 *Description:
 *Recursively moves through pedigree until the child leaves are found to find structure
 */
/*
string descendTree(Sample* sample, int level) {
	string infospace = "     ";
	for (int i = 0; i < (level * 2); i++) {
		infospace += "     ";
	}
	vector<Sample*>* children = sample->getChildren();
	if (children->size() == 0) {
		string v = infospace + "->Ind: " + sample->getInd() + "\n" + infospace
				+ "->Mom: " + sample->getMomID() + "\n" + infospace + "->Dad: "
				+ sample->getDadID() + "\n";
		if (sample->getSib() != NULL) {
			v += infospace + "->Sibling: " + sample->getSib()->getInd() + "\n";
		}
		v += infospace + "->Children: 0\n";
		return v;
	}
	int csize = children->size();

	string val = "";
	string v = infospace + "->Ind: " + sample->getInd() + "\n" + infospace
			+ "->Mom: " + sample->getMomID() + "\n" + infospace + "->Dad: "
			+ sample->getDadID() + "\n";
	if (sample->getSib() != NULL) {
		v += infospace + "->Sibling: " + sample->getSib()->getInd() + "\n";
	}
	v += infospace + "->Children: " + getString<int> (
			sample->getChildren()->size()) + "\n";
	val += v;
	for (int c = 0; c < csize; c++) {
		string temp = descendTree((*children)[c], (level + 1));

		val += temp;
	}
	return val;
}
*/

/*
 *Function: descendTree3
 *Description: returns # of levels
 *Recursively moves through pedigree until the child leaves are found to find structure
 */
/*
map<int, vector<Sample*> > descendTree3(Sample* sample, int level) {
	map<int, vector<Sample*> > values;

	vector<Sample*>* children = sample->getChildren();
	if (children->size() == 0) {
		return values;
		//return level;
	}
	int csize = children->size();
	vector<int> levels(csize);
	for (int c = 0; c < csize; c++) {
		map<int, vector<Sample*> > tempvalues = descendTree3((*children)[c],
				(level + 1));
		map<int, vector<Sample*> >::iterator iter;
		for (iter = tempvalues.begin(); iter != tempvalues.end(); iter++) {
			vector<Sample*> samps = iter->second;
			for (unsigned int is = 0; is < samps.size(); is++) {
				values[iter->first].push_back(samps[is]);
			}
		}
		//		if((*children)[c]->getPheno() == 2){// && !(*children)[c]->isFlagged()){
		values[(level + 1)].push_back((*children)[c]);
		//(*children)[c]->setFlag(true);
		//		}
	}
	//	sort(levels.begin(), levels.end());
	//	level = levels[levels.size() - 1];
	return values;
}
*/

/*
 *Function: compileOutputs
 *Description:
 *Compiles all QC & analysis outputs into single marker, sample, and family files
 *Performed when all steps are complete.
 */
/*
void compileOutputs(vector<Marker*>* markers, vector<Family*>* families,
		vector<Sample*>* samples) {
	opts::printLog("Compiling QC and Analysis output files...\n");

	map<string, vector<string> > filenames = opts::getFilenames();
	if (filenames["Marker"].size() > 0) {
		map<Marker*, vector<string> > marker_output;
		vector<string> all_columns;
		map<string, Marker*> snpmap;
		for (unsigned int i = 0; i < markers->size(); i++) {
			Marker* mark = (*markers)[i];
			string chr = getString<int> (mark->getChrom());
			string bp = getString<int> (mark->getBPLOC());
			snpmap[chr + "#" + bp] = mark;
		}

		cout << "Markers initialized...\n";

		for (unsigned int i = 0; i < filenames["Marker"].size(); i++) {
			string file = filenames["Marker"][i];
			string step = opts::filesteps[file];

			for (unsigned int j = 0; j < opts::fileheaders[file].size(); j++) {
				cout << "Pushing " << file << " : " << j << endl;
				all_columns.push_back("(" + file + ")"
						+ opts::fileheaders[file][j]);
			}
		}
		string outfile = opts::_OUTPREFIX_ + "marker_summary.txt";
		opts::printLog("Working on compiling SNP information...[" + outfile
				+ "]\n");
		for (unsigned int i = 0; i < filenames["Marker"].size(); i++) {
			string filename = filenames["Marker"][i];
			string step = opts::filesteps[filename];
			vector<string> filecols = opts::fileheaders[filename];

			opts::printLog("\tParsing " + filename + "\n");

			ifstream input;
			input.open(filename.c_str(), ios::in);

			if (!input) {
				cerr << "Error opening file: " << filename << endl;
				exit(1);
			}

			string header = "";
			getline(input, header);

			vector<string> columns;
			General::Tokenize(header, columns, "\t");
			int chrloc = -1;
			int bploc = -1;
			for (unsigned int l = 0; l < columns.size(); l++) {
				if (columns[l] == "Chrom") {
					chrloc = l;
				}
				if (columns[l] == "bploc") {
					bploc = l;
				}
				if (chrloc > -1 && bploc > -1) {
					break;
				}
			}
			string line;
			while (getline(input, line)) {
				int count = 0;
				vector<string> elems;
				string token;
				istringstream isstream(line);
				while (getline(isstream, token, '\t')) {
					elems.push_back(token);
				}

				//General::Tokenize(line, elems, "\t");
				if (elems.size() == 0) {
					continue;
				}
				if (elems.size() != columns.size()) {
					cout << "Error on line: " << count << endl;
					exit(1);
				}
				int chrom = atoi(elems[chrloc].c_str());
				int bp = atoi(elems[bploc].c_str());
				string key = elems[chrloc] + "#" + elems[bploc];

				Marker* itermark = snpmap[key];//(*markers)[count];
				///				while(itermark->getChrom() != chrom && itermark->getBPLOC() != bploc){
				///					if(count + 1 >= markers->size()){
				///						itermark = NULL;
				///						break;
				///					}
				///					itermark = (*markers)[++count];
				///				}
				//vector<Marker*>::iterator itermark = find_if(markers->begin(), markers->end(), FindMarkerByChromBploc(chrom, bp));
				if (itermark == NULL) {
					cout << "Cannot find marker with chrom = " << chrom
							<< " and bploc = " << bp << endl;
					exit(1);
				}
				Marker* mark = itermark;
				map<Marker*, vector<string> >::iterator data =
						marker_output.find(mark);
				if (data == marker_output.end()) {
					vector<string> mydata;
					mydata.resize(all_columns.size(), "N/A");
					marker_output[mark] = mydata;
				}
				for (unsigned int j = 0; j < filecols.size(); j++) {
					vector<string>::iterator realcol = find(columns.begin(),
							columns.end(), filecols[j]);
					vector<string>::iterator allloc = find(all_columns.begin(),
							all_columns.end(), "(" + filename + ")"
									+ filecols[j]);
					if (realcol != columns.end()) {
						int myloc = realcol - columns.begin();
						int myallloc = allloc - all_columns.begin();
						marker_output[mark][myallloc] = elems[myloc];
					}
				}
				count++;
			}
			input.close();

		}
		ofstream out(outfile.c_str());
		out << "Chrom\trsID\tProbeID\tbploc";
		for (unsigned int i = 0; i < all_columns.size(); i++) {
			out << "\t" << all_columns[i];
		}
		out << "\n";
		map<Marker*, vector<string> >::iterator iter;
		for (iter = marker_output.begin(); iter != marker_output.end(); iter++) {
			Marker* mark = iter->first;
			vector<string> data = iter->second;
			out << mark->getChrom() << "\t" << mark->getRSID() << "\t"
					<< mark->getProbeID() << "\t" << mark->getBPLOC();
			for (unsigned int i = 0; i < data.size(); i++) {
				out << "\t" << data[i];
			}
			out << "\n";
		}
		out.close();
	}
	if (filenames["Family"].size() > 0) {
		map<Family*, vector<string> > family_output;
		vector<string> all_columns;
		for (unsigned int i = 0; i < filenames["Family"].size(); i++) {
			string file = filenames["Family"][i];
			string step = opts::filesteps[file];

			for (unsigned int j = 0; j < opts::fileheaders[file].size(); j++) {
				all_columns.push_back("(" + file + ")"
						+ opts::fileheaders[file][j]);
			}
		}
		string outfile = opts::_OUTPREFIX_ + "family_summary.txt";
		opts::printLog("Working on compiling Family information...[" + outfile
				+ "]\n");
		for (unsigned int i = 0; i < filenames["Family"].size(); i++) {
			string filename = filenames["Family"][i];
			string step = opts::filesteps[filename];
			vector<string> filecols = opts::fileheaders[filename];

			opts::printLog("\tParsing " + filename + "\n");

			ifstream input;
			input.open(filename.c_str(), ios::in);

			if (!input) {
				cerr << "Error opening file: " << filename << endl;
				exit(1);
			}

			string header = "";
			getline(input, header);

			vector<string> columns;
			General::Tokenize(header, columns, "\t");
			int famloc = -1;
			for (unsigned int l = 0; l < columns.size(); l++) {
				if (columns[l] == "FamID") {
					famloc = l;
				}
				if (famloc > -1) {
					break;
				}
			}
			if (famloc < 0) {
				cout << "FamID column not found!\n";
				exit(1);
			}
			string line;
			int count = 1;
			while (getline(input, line)) {
				count++;
				vector<string> elems;
				string token;
				istringstream isstream(line);
				while (getline(isstream, token, '\t')) {
					elems.push_back(token);
				}
				//General::Tokenize(line, elems, "\t");
				if (elems.size() == 0) {
					continue;
				}
				string famid = elems[famloc];
				vector<Family*>::iterator iterfam = find_if(families->begin(),
						families->end(), FindFamily(famid));
				if (iterfam == families->end()) {
					cout << "Cannot find family with famid = " << famid << endl;
					exit(1);
				}
				Family* fam = *iterfam;
				map<Family*, vector<string> >::iterator data =
						family_output.find(fam);
				if (data == family_output.end()) {
					vector<string> mydata;
					mydata.resize(all_columns.size(), "N/A");
					family_output[fam] = mydata;
				}
				for (unsigned int j = 0; j < filecols.size(); j++) {
					vector<string>::iterator realcol = find(columns.begin(),
							columns.end(), filecols[j]);
					vector<string>::iterator allloc = find(all_columns.begin(),
							all_columns.end(), "(" + filename + ")"
									+ filecols[j]);
					if (realcol != columns.end()) {
						int myloc = realcol - columns.begin();
						int myallloc = allloc - all_columns.begin();
						family_output[fam][myallloc] = elems[myloc];
					}
				}
			}
			input.close();

		}
		ofstream out(outfile.c_str());
		out << "FamID\tNumInds";
		for (unsigned int i = 0; i < all_columns.size(); i++) {
			out << "\t" << all_columns[i];
		}
		out << "\n";
		map<Family*, vector<string> >::iterator iter;
		for (iter = family_output.begin(); iter != family_output.end(); iter++) {
			Family* fam = iter->first;
			vector<string> data = iter->second;
			out << fam->getFamID() << "\t" << fam->getSamples()->size();
			for (unsigned int i = 0; i < data.size(); i++) {
				out << "\t" << data[i];
			}
			out << "\n";
		}
		out.close();
	}
	if (filenames["Sample"].size() > 0) {
		map<Sample*, vector<string> > sample_output;
		vector<string> all_columns;
		for (unsigned int i = 0; i < filenames["Sample"].size(); i++) {
			string file = filenames["Sample"][i];
			string step = opts::filesteps[file];

			for (unsigned int j = 0; j < opts::fileheaders[file].size(); j++) {
				all_columns.push_back("(" + file + ")"
						+ opts::fileheaders[file][j]);
			}
		}
		string outfile = opts::_OUTPREFIX_ + "sample_summary.txt";
		opts::printLog("Working on compiling Sample information...[" + outfile
				+ "]\n");
		for (unsigned int i = 0; i < filenames["Sample"].size(); i++) {
			string filename = filenames["Sample"][i];
			string step = opts::filesteps[filename];
			vector<string> filecols = opts::fileheaders[filename];

			opts::printLog("\tParsing " + filename + "\n");

			ifstream input;
			input.open(filename.c_str(), ios::in);

			if (!input) {
				cerr << "Error opening file: " << filename << endl;
				exit(1);
			}

			string header = "";
			getline(input, header);

			vector<string> columns;
			General::Tokenize(header, columns, "\t");
			int famloc = -1;
			int indloc = -1;
			for (unsigned int l = 0; l < columns.size(); l++) {
				if (columns[l] == "FamID") {
					famloc = l;
				}
				if (columns[l] == "IndID") {
					indloc = l;
				}
				if (famloc > -1 && indloc > -1) {
					break;
				}
			}
			string line;
			int count = 1;
			while (getline(input, line)) {
				count++;
				vector<string> elems;
				string token;
				istringstream isstream(line);
				while (getline(isstream, token, '\t')) {
					elems.push_back(token);
				}
				//General::Tokenize(line, elems, "\t");
				if (elems.size() == 0) {
					continue;
				}
				string famid = elems[famloc];
				string indid = elems[indloc];
				vector<Sample*>::iterator itersamp = find_if(samples->begin(),
						samples->end(), FindSampleByFamAndID(famid, indid));
				if (itersamp == samples->end()) {
					cout << "Cannot find sample with famid = " << famid
							<< " and indid = " << indid << endl;
					exit(1);
				}
				Sample* samp = *itersamp;
				map<Sample*, vector<string> >::iterator data =
						sample_output.find(samp);
				if (data == sample_output.end()) {
					vector<string> mydata;
					mydata.resize(all_columns.size(), "N/A");
					sample_output[samp] = mydata;
				}
				for (unsigned int j = 0; j < filecols.size(); j++) {
					vector<string>::iterator realcol = find(columns.begin(),
							columns.end(), filecols[j]);
					vector<string>::iterator allloc = find(all_columns.begin(),
							all_columns.end(), "(" + filename + ")"
									+ filecols[j]);
					if (realcol != columns.end()) {
						int myloc = realcol - columns.begin();
						int myallloc = allloc - all_columns.begin();
						sample_output[samp][myallloc] = elems[myloc];
					}
				}
			}
			input.close();

		}
		ofstream out(outfile.c_str());
		out << "FamID\tIndID\tSex\tAffection_Status";
		for (unsigned int i = 0; i < all_columns.size(); i++) {
			out << "\t" << all_columns[i];
		}
		out << "\n";
		map<Sample*, vector<string> >::iterator iter;
		for (iter = sample_output.begin(); iter != sample_output.end(); iter++) {
			Sample* samp = iter->first;
			vector<string> data = iter->second;
			out << samp->getFamID() << "\t" << samp->getInd() << "\t";
			if (samp->getSex()) {
				out << "M\t";
			} else {
				out << "F\t";
			}
			out << samp->getPheno();

			for (unsigned int i = 0; i < data.size(); i++) {
				out << "\t" << data[i];
			}
			out << "\n";
		}
		out.close();
	}
	if (filenames["Batch"].size() > 0) {
		map<string, vector<string> > batch_output;
		vector<string> all_columns;
		for (unsigned int i = 0; i < filenames["Batch"].size(); i++) {
			string file = filenames["Batch"][i];
			string step = opts::filesteps[file];

			for (unsigned int j = 0; j < opts::fileheaders[file].size(); j++) {
				all_columns.push_back("(" + file + ")"
						+ opts::fileheaders[file][j]);
			}
		}
		string outfile = opts::_OUTPREFIX_ + "batch_summary.txt";
		opts::printLog("Working on compiling Batch information...[" + outfile
				+ "]\n");
		for (unsigned int i = 0; i < filenames["Batch"].size(); i++) {
			string filename = filenames["Batch"][i];
			string step = opts::filesteps[filename];
			vector<string> filecols = opts::fileheaders[filename];

			opts::printLog("\tParsing " + filename + "\n");

			ifstream input;
			input.open(filename.c_str(), ios::in);

			if (!input) {
				cerr << "Error opening file: " << filename << endl;
				exit(1);
			}

			string header = "";
			getline(input, header);

			vector<string> columns;
			General::Tokenize(header, columns, "\t");
			int famloc = -1;
			for (unsigned int l = 0; l < columns.size(); l++) {
				if (columns[l] == "Batch") {
					famloc = l;
				}
				if (famloc > -1) {
					break;
				}
			}
			if (famloc < 0) {
				cout << "Batch column not found!\n";
				exit(1);
			}
			string line;
			int count = 1;
			while (getline(input, line)) {
				count++;
				vector<string> elems;
				string token;
				istringstream isstream(line);
				while (getline(isstream, token, '\t')) {
					elems.push_back(token);
				}
				//General::Tokenize(line, elems, "\t");
				if (elems.size() == 0) {
					continue;
				}
				string batchid = elems[famloc];
				//				vector<>::iterator iterfam = find_if(families->begin(), families->end(), FindFamily(famid));
				//				if(iterfam == families->end()){
				//					cout << "Cannot find family with famid = " << famid << endl;
				//					exit(1);
				//				}
				//				Family* fam = *iterfam;
				map<string, vector<string> >::iterator data =
						batch_output.find(batchid);
				if (data == batch_output.end()) {
					vector<string> mydata;
					mydata.resize(all_columns.size(), "N/A");
					batch_output[batchid] = mydata;
				}
				for (unsigned int j = 0; j < filecols.size(); j++) {
					vector<string>::iterator realcol = find(columns.begin(),
							columns.end(), filecols[j]);
					vector<string>::iterator allloc = find(all_columns.begin(),
							all_columns.end(), "(" + filename + ")"
									+ filecols[j]);
					if (realcol != columns.end()) {
						int myloc = realcol - columns.begin();
						int myallloc = allloc - all_columns.begin();
						batch_output[batchid][myallloc] = elems[myloc];
					}
				}
			}
			input.close();

		}
		ofstream out(outfile.c_str());
		out << "Batch";
		for (unsigned int i = 0; i < all_columns.size(); i++) {
			out << "\t" << all_columns[i];
		}
		out << "\n";
		map<string, vector<string> >::iterator iter;
		for (iter = batch_output.begin(); iter != batch_output.end(); iter++) {
			string batchid = iter->first;
			vector<string> data = iter->second;
			out << batchid;
			for (unsigned int i = 0; i < data.size(); i++) {
				out << "\t" << data[i];
			}
			out << "\n";
		}
		out.close();
	}

}

*/
