#ifndef PLATO_MAIN_H
#define PLATO_MAIN_H

#include "config.h"
#include <deque>
#include <map>
#include <utility>

#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>

// this is the entry point for the "master" (this is essentially the "old" main)
int master_main(int, char**);

// handles both MPI and regular invocations
int main(int, char**);

void MPISendResponses(std::map<int, std::deque<std::pair<unsigned int, const char*> > >& resp_queue_map, boost::mutex& resp_mutex);

void MPIProbeInput(std::map<int, std::deque<std::pair<unsigned int, const char*> > >& resp_queue_map, boost::mutex& resp_mutex, bool& running);

void print_steps();

/*
void startProcess(vector<Process*>&, int, InputFilter*);
void parseInput(const string&, vector<Process*>& order_out);
void usage();
void error_check();
void parseParameters();
void print_help();
void runStep(Process*, DataSet*);
void flipStrand(vector<Marker*>*);
void printFamilies(vector<Family*>*);
void printOptions();
string descendTree(Sample*, int);
map<int, vector<Sample*> > descendTree3(Sample*, int);
void compileOutputs(vector<Marker*>*, vector<Family*>*, vector<Sample*>*);
void webcheck(vector<string>, map<string, vector<string> >);
map<string, vector<string> > getBatchArgs(string);
vector<vector<Process*> > optimize(vector<Process*>&);
*/

#endif
