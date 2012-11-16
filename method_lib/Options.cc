#include <iostream>
#include <string>
#include "Options.h"

using namespace std;
namespace Methods{
		 string opts::_WASPVER_ = "0.83";
		 int opts::_MARKERS_FOUND_ = 0;
		 int opts::_MARKERS_WORKING_ = 0;
		 int opts::_COVS_FOUND_ = 0;
		 int opts::_COVS_WORKING_ = 0;
		 int opts::_TRAITS_FOUND_ = 0;
		 int opts::_TRAITS_WORKING_ = 0;
		 int opts::_SAMPLES_FOUND_ = 0;
		 int opts::_SAMPLES_WORKING_ = 0;
		 int opts::_FAMILIES_FOUND_ = 0;
		 int opts::_FAMILIES_WORKING_ = 0;
		 int opts::_STRUCTSPACE_ = 0;
		 int opts::_GRRSPACE_ = 0;
		 int opts::_BP_HIGH_ = 0;
		 int opts::_BP_LOW_ = 0;
		 int opts::_CHROM_ = 0;
		 int opts::_BP_SPACE_ = 0;
		 int opts::_CHRX_ = 23;
		 int opts::_CHRY_ = 24;
		 int opts::_CHRXY_ = 25;
		 int opts::_MITO_ = 26;
		 bool opts::_MAP_INCLUDES_REF_ = false;
		 bool opts::_TODIGIT_ = false;
		 bool opts::_COMPILE_OUTPUTS_ = true;
		 bool opts::_KEEP_EXC_SAMPLES_ = false;
		 bool opts::_DOG_ = false;
		 bool opts::_OUTPUTREMAINS_ = false;
         bool opts::_RESTRICT_ = false;
         bool opts::_RESTRICT_SAMPS_ = false;
         bool opts::_CASECONTROL_ = false;
         bool opts::_DBOUTPUT_ = false;
         bool opts::_MARKERLIST_ = false;
         bool opts::_STRATIFY_ = false;
         bool opts::_DBINPUT_ = false;
		 bool opts::_ENZYMES_ = false;
		 bool opts::_QSCUT_ = false;
		 bool opts::_MAKEBIN_ = false;
		 bool opts::_STRUCTALL_ = false;
		 bool opts::_MICROSATS_ = false;
		 bool opts::_FLIPSTRAND_ = false;
		 bool opts::_BP_LOW_LIMIT_ = false;
		 bool opts::_BP_HIGH_LIMIT_ = false;
		 bool opts::_CHROM_LIMIT_ = false;
		 bool opts::_BP_SPACE_LIMIT_ = false;
		 bool opts::_PRINTFAMS_ = false;
		 bool opts::_FREQ_FILE_EXISTS_ = false;
		 bool opts::_COVFILE_EXISTS_ = false;
		 bool opts::_TRAITFILE_EXISTS_ = false;
		 string opts::_SAMPLEBPRANGEFILTER_ = "";
		 string opts::_COVAR_MISSING_ = "-99999";
		 string opts::_TRAIT_MISSING_ = "-99999";
		 string opts::_FREQ_FILE_ = "";
		 string opts::_PEDINFO_ = "";
         string opts::_PEDFILE_ = "";
         string opts::_MAPFILE_ = "";
		 string opts::_TPEDFILE_ = "";
		 string opts::_FAMFILE_ = "";
		 string opts::_MDRFILE_ = "";
		 string opts::_MDRPEDFILE_ = "";
		 string opts::_MDRMAPFILE_ = "";
         string opts::_MAPDESC_ = "";
         string opts::_SAMPDESC_ = "";
         string opts::_MARKEXCL_ = "";
         string opts::_SAMPEXCL_ = "";
		 string opts::_ZEROGENOFILE_ = "";
		 string opts::_NOCALL_ = "0";
		 string opts::_NOCALLMDR_ = "-1";
         string opts::USER = "";
         string opts::PASS = "";
         string opts::DBNAME = "";
		 string opts::_QSCUTOFF_ = "1";
		 string opts::_QSFILE_ = "";
		 string opts::_EXCSAMPS_ = "";
		 string opts::_EXCCOVS_ = "";
		 string opts::_EXCTRAITS_ = "";
		 string opts::_INCCOVS_ = "";
		 string opts::_INCTRAITS_ = "";
		 string opts::_EXCMARKERS_ = "";
		 string opts::_EXCFAMILIES_ = "";
		 string opts::_INCMARKERS_ = "";
		 string opts::_INCCENTERS_ = "";
		 string opts::_INCSAMPLES_ = "";
		 string opts::_INCFAMILIES_ = "";
		 string opts::_BINPREFIX_ = "";
		 string opts::_STRUCT_STRAT_FILE_ = "";
		 string opts::_OUTPREFIX_ = "";
		 string opts::_FLIPFILE_ = "";
		 string opts::_COVFILE_ = "";
		 string opts::_TRAITFILE_="";
		 ofstream opts::_LOG_;
		 vector<string> opts::cov_loc;
		 vector<string> opts::trait_loc;
		 vector<int> opts::trait_covs;
		 map<string, vector<string> > opts::filenames;
		 map<string, vector<string> > opts::fileheaders;
		 map<string, string> opts::filesteps;
		 map<string, Sample*> opts::pedinfo;
		 map<string, vector<string> > opts::zerogenoinfo;

		 void opts::printLog(string s){
			 opts::_LOG_ << s;
			 opts::_LOG_.flush();

			 cout << s;
			 cout.flush();
		 }

		 void opts::addFile(string type, string step, string s){
			 filenames[type].push_back(s);
			 filesteps[s] = step;
		 }

		 void opts::addHeader(string file, string head){
			 fileheaders[file].push_back(head);
		 }

}
