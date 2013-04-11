#ifndef wasp_H
#define wasp_H

#include <DataSet.h>
#include <Options.h>
#include <Helpers.h>
#include <InputFilter.h>
#include <StepOptions.h>
#include <MethodException.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <map>

#include <vector>

#include "config.h"

#include "Finalize.h"

#include "ProcessFactory.h"
#include "Process.h"

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>

#define TOTAL_STEPS 15
using namespace Methods;

int main(int, char**);

static void Initialize();
void print_steps();
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


#endif
