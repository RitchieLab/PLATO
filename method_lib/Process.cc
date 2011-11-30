#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Process.h"
using namespace std;

void Process::~Process(){}
void Process::PrintSummary(){}
void Process::filter(){};
void Process::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){};
void Process::setThreshold(string s){};
void Process::FilterSummary(){};
void Process::setRank(int r){};
void Process::setDBOUT(){};
void Process::setMarkerList(){};
void Process::setStratify(){};
void Process::setOrder(int o){};
void Process::setOverwrite(bool v){};
bool Process::hasIncExc(){};
