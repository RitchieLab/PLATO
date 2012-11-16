#ifndef OPTIONS_H
#define OPTIONS_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "Sample.h"



using namespace std;

namespace Methods{
class Sample;
class opts {
	public:
		static const int SR_NUMPERNODE = 10;
		static const int SR_MARKERDATA = 20;
		static const int SR_ALLELEDATA = 30;
		static const int SR_FAMILYDATA = 40;
		static const int SR_INDDATA = 50;
		static const int SR_CHREFFDATA = 60;
		static const int SR_ENZYMEDATA = 70;
		static const int DONE = 999;
		static const int RECEIVE = 1;
		static const int MAXROWS = 150000;
		static int _STRUCTSPACE_;
		static int _GRRSPACE_;
		static int _MARKERS_FOUND_;
		static int _MARKERS_WORKING_;
		static int _COVS_FOUND_;
		static int _COVS_WORKING_;
		static int _TRAITS_FOUND_;
		static int _TRAITS_WORKING_;
		static int _SAMPLES_FOUND_;
		static int _SAMPLES_WORKING_;
		static int _FAMILIES_FOUND_;
		static int _FAMILIES_WORKING_;
		static int _BP_HIGH_;
		static int _BP_LOW_;
		static int _CHROM_;
		static int _BP_SPACE_;
		static int _CHRX_;
		static int _CHRY_;
		static int _CHRXY_;
		static int _MITO_;
		static int _NUMTHREADS_;
		static int _PHENO_MISS_;
		static bool _MAP_INCLUDES_REF_;
		static bool _TODIGIT_;
		static bool _COMPILE_OUTPUTS_;
		static bool _KEEP_EXC_SAMPLES_;
		static bool _DOG_;
		static bool _OUTPUTREMAINS_;
		static bool _RESTRICT_;
		static bool _RESTRICT_SAMPS_;
		static bool _CASECONTROL_;
		static bool _DBOUTPUT_;
		static bool _MARKERLIST_;
		static bool _STRATIFY_;
		static bool _DBINPUT_;
		static bool _ENZYMES_;
		static bool _QSCUT_;
		static bool _MAKEBIN_;
		static bool _STRUCTALL_;
		static bool _MICROSATS_;
		static bool _FLIPSTRAND_;
		static bool _BP_LOW_LIMIT_;
		static bool _BP_HIGH_LIMIT_;
		static bool _CHROM_LIMIT_;
		static bool _BP_SPACE_LIMIT_;
		static bool _PRINTFAMS_;
		static bool _FREQ_FILE_EXISTS_;
		static bool _COVFILE_EXISTS_;
		static bool _TRAITFILE_EXISTS_;
		static bool _GROUP_EXISTS_;
		static bool _THREADS_;
		static bool _BINTRAIT_;
		static bool _AUTOONLY_;
		//lgen input option
		static bool _COMPOUND_GENOTYPES_;
		static string _SAMPLEBPRANGEFILTER_;
		static string _COVAR_MISSING_;
		static string _TRAIT_MISSING_;
		static string _WASPVER_;
		static string _FREQ_FILE_;
		static string _PEDINFO_;
		static string _PEDFILE_;
		//lgen
		static string _LGENFILE_;
		static string _REFERENCE_FILE_;

		static string _MAPFILE_;
		static string _FAMFILE_;
		static string _TPEDFILE_;
		static string _MDRFILE_;
		static string _MDRPEDFILE_;
		static string _MDRMAPFILE_;
		static string _MAPDESC_;
		static string _SAMPDESC_;
		static string _MARKEXCL_;
		static string _SAMPEXCL_;
		static string _ZEROGENOFILE_;
		static string USER;
		static string PASS;
		static string DBNAME;
		static string _QSCUTOFF_;
		static string _QSFILE_;
		static string _EXCSAMPS_;
		static string _EXCCOVS_;
		static string _EXCTRAITS_;
		static string _INCCOVS_;
		static string _INCTRAITS_;
		static string _EXCMARKERS_;
		static string _EXCFAMILIES_;
		static string _INCMARKERS_;
		static string _INCCENTERS_;
		static string _INCSAMPLES_;
		static string _INCFAMILIES_;
		static string _BINPREFIX_;
		static string _STRUCT_STRAT_FILE_;
		static string _OUTPREFIX_;
		static string _FLIPFILE_;
		static string _COVFILE_;
		static string _TRAITFILE_;
		static string _NOCALL_;
		static string _NOCALLMDR_;
		static ofstream _LOG_;
		static vector<string> cov_loc;
		static vector<string> trait_loc;
		static vector<int> trait_covs;
		static map<string, vector<string> > filenames;
		static map<string, vector<string> > fileheaders;
		static map<string, string> filesteps;
		static void printLog(string s);
		static void addFile(string type, string step, string s);
		static void addHeader(string file, string head);
		static map<string, vector<string> > getFilenames(){return filenames;};
		static map<string, Sample*> pedinfo;
		static map<string, vector<string> > zerogenoinfo;


};
};
#endif
