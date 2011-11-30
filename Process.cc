#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Process.h"
using namespace std;

Process::~Process(){}
void Process::PrintSummary(){}
void Process::filter(){};
//void Process::process(Connection* con, Families* families, Markers* markers){};
//void Process::process(Families* families, Markers* markers){};
void Process::process(DataSet* ds){};
void Process::setThreshold(string s){};
void Process::FilterSummary(){};
void Process::setRank(int r){};
//void Process::updateFamsMarks(Families* f, Markers* m){};
void Process::setDBOUT(){};
void Process::setMarkerList(){};
void Process::setStratify(){};
void Process::setOrder(int o){};
void Process::setOverwrite(bool v){};
bool Process::hasIncExc(){};
