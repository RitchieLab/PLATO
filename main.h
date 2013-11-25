#ifndef PLATO_MAIN_H
#define PLATO_MAIN_H

int main(int, char**);

static void Initialize();
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
