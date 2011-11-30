#ifndef wasp_H
#define wasp_H

#include <DataSet.h>
#include <Options.h>
#include <Helper.h>
#include <InputFilter.h>
#include <StepOptions.h>
#include <MethodException.h>
#include <stdio.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <map>
#include "sockets.h"
#include "Step.h"
#include "ExampleModule.h"
#include "ProcessKinship.h"
#include "ProcessEarth.h"
#include "ProcessMendelianErrors.h"
#include "ProcessMarkerGenoEff.h"
#include "ProcessAlleleFrequency.h"
#include "ProcessSampleGenoEff.h"
#include "ProcessHWEquilibrium.h"
#include "ProcessGenderCheck.h"
//#include "QualityScore.h"
#include "PercentByFamily.h"
#include "ProcessPEDOutput.h"
#include "ProcessBINOutput.h"
#include "ProcessTPEDOutput.h"
#include "ProcessGRROutput.h"
#include "ProcessSuperlinkOutput.h"
#include <Globals.h>
//#include "PEDOutput.h"
//#include "QSOutput.h"
#include "ProcessRunTDT.h"
#include "ProcessCaConChisq.h"
#include "ProcessSTRUCTOutput.h"
#include "ProcessPHASEOutput.h"
#include "ProcessBEAGLEOutput.h"
#include "ProcessLAPISOutput.h"
#include "ProcessMDROutput.h"
#include "ProcessHomozygous.h"
#include "ProcessLD.h"
#include "ProcessPowerMarkerOutput.h"
#include "ProcessDeletions.h"
#include "ProcessMitoCheck.h"
#include "ProcessQTDTOutput.h"
#include "ProcessFBATOutput.h"
#include "ProcessPDT2Output.h"
#include "ProcessConcordance.h"
#include "ProcessLogReg.h"
#include "ProcessLinearReg.h"
#include "ProcessEigenstratOutput.h"
#include "ProcessIBS.h"
#include "ProcessCMH.h"
#include "ProcessMDR.h"
#include "ProcessMDRPDT.h"
#include "ProcessClusterMissing.h"
#include "ProcessFilterProcess.h"
#include "ProcessFst.h"
#include "ProcessEpistasis.h"
#include "Finalize.h"
#include "Process.h"

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>

//#include "Finalize.h"
//#include "Setup.h"
//#include "Options.h"
//#include "Helper.h"

//string USER = "";
//string PASS = "";
//string DBNAME = "";
#define TOTAL_STEPS 15
using namespace Methods;

typedef map<string,Step> STEPS;
typedef vector<Step> ORDER;

int main(int, char**);
static void Initialize();
STEPS initializeSteps();
Step initializeSteps(string);
void print_steps(STEPS);
//void startProcess(ORDER*, STEPS*, void *, int, int, vector<string>, vector<string>);
void startProcess(ORDER*, void *, int, InputFilter*);
//ORDER parseInput(string, STEPS*);
ORDER parseInput(string);
//vector<string> readMarkers(string);
//vector<string> readCenters(string);
void usage();
void error_check();
void parseParameters();
//void Tokenize(const string&, vector<string>&, const string&);
//vector<string> ParseDelimitedLine(string);
void print_help();
void runStep(Step, DataSet*);
//void fill_temp_marker_table(Connection*, Markers*, int);
//void writeBit(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
//void readBin(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
//void readPed(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
//void readMap(vector<Marker*>*, vector<int>*);
//void readTPed(vector<Marker*>*, vector<Sample*>*, vector<int>*);
//void readTFam(vector<Sample*>*, vector<Family*>*);
//bool readString(FILE*, string*);
//void assignLinks(vector<Family*>*);
//void reorderAlleles(vector<Sample*>*, vector<Marker*>*);
//void remapSamples(vector<Sample*>*, vector<Marker*>*, vector<int>*, int);
void flipStrand(vector<Marker*>*);
void printFamilies(vector<Family*>*);
string descendTree(Sample*, int);
map<int, vector<Sample*> > descendTree3(Sample*, int);
void compileOutputs(vector<Marker*>*, vector<Family*>*, vector<Sample*>*);
void webcheck(vector<string>, map<string, vector<string> >);
map<string, vector<string> > getBatchArgs(string);


//typedef list<int> ORDER;
/*
int SR_NUMPERNODE = 10;
int SR_MARKERDATA = 20;
int SR_ALLELEDATA = 30;
//#define SR_FAMILYDATA 40
int SR_FAMILYDATA = 40;
int SR_INDDATA = 50;
int SR_CHREFFDATA = 60;
int SR_ENZYMEDATA = 70;
int DONE = 999;
int RECEIVE = 1;
int _RESTRICT_ = 0;
int _RESTRICT_SAMPS_ = 0;
int _CASECONTROL_ = 0;
int _DBOUTPUT_ = 0;
int _MARKERLIST_ = 0;
int _STRATIFY_ = 0;
int _DBINPUT_ = 0;
string _PEDFILE_ = "";
string _MAPFILE_ = "";
string _MAPDESC_ = "";
string _SAMPDESC_ = "";
string _MARKEXCL_ = "";
string _SAMPEXCL_ = "";
*/
#endif
