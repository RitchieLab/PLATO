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
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <string.h>
#include <string>
#include <iomanip>
#include <fstream>
#include <bitset>
#include "Globals.h"

#ifdef _USEMPI_
#include <mpi.h>
#endif
#include <vector>
#include "General.h"
#include "wasp.h"


using namespace std;

//enum for batch file step switch
enum StepValue{ 
					   e_examplemodule,    //examplemodule
					   e_marker_geno_eff,  //marker-geno-eff
					   e_sample_geno_eff,  //sample-geno-eff
					   e_family_geno_eff,  //family-geno-eff
					   e_allele_freq,      //allele-freq
					   e_mendelian_error,  //mendelian-error
					   e_hw,               //hw
					   e_gender_error,     //gender-error
					   e_tdt,              //tdt
					   e_grr_output,       //grr-output
					   e_ped_output,       //ped-output
					   e_chisquare, //ca-con-chisquare
					   e_structure_output, //structure-output
					   e_phase_output,	   //phase-output
					   e_beagle_output,	   //beagle-output
					   e_lapis_output,	   //lapis-output
					   e_superlink_output, //superlink-output
					   e_homozygous, 	   //homozygous
					   e_ld, 	           //ld
					   e_powermarker_output,//powermarker-output
					   e_deletions,	       //deletions
					   e_fbat_output,	   //fbat-output
					   e_qtdt_output,	   //qtdt-output
					   e_mdr_output, 	   //mdr-output
					   e_pdt2_output,	   //output-pdt2
					   e_concordance,	   //concordance
					   e_mitocheck         //mito-check
					 };
enum cmdArgs{
	a_ped,
	a_map,
	a_bin_input,
	a_incmarkers,
	a_excmarkers,
	a_incsamples,
	a_excsamples,
	a_incfamilies,
	a_excfamilies,
	a_chrom,
	a_bp_min,
	a_bp_max,
	a_bp_span,
	a_pedinfo,
	a_freq_file
	
};
static map<string, StepValue> s_mapStepValues;

int main(int, char**);
static void Initialize();
STEPS initializeSteps();
Step initializeSteps(string);
void print_steps(STEPS);
//void startProcess(ORDER*, STEPS*, void *, int, int, vector<string>, vector<string>);
void startProcess(ORDER*, void *, int, int, vector<string>, vector<string>);
//ORDER parseInput(string, STEPS*);
ORDER parseInput(string);
vector<string> readMarkers(string);
vector<string> readCenters(string);
void usage();
void error_check();
void parseParameters();
//void Tokenize(const string&, vector<string>&, const string&);
//vector<string> ParseDelimitedLine(string);
void print_help();
//void fill_temp_marker_table(Connection*, Markers*, int);
void writeBit(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
void readBin(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
void readPed(vector<Sample*>*, vector<Family*>*, vector<Marker*>*, vector<int>*);
void readMap(vector<Marker*>*, vector<int>*);
void readTPed(vector<Marker*>*, vector<Sample*>*, vector<int>*);
void readTFam(vector<Sample*>*, vector<Family*>*);
//bool readString(FILE*, string*);
//void assignLinks(vector<Family*>*);
//void reorderAlleles(vector<Sample*>*, vector<Marker*>*);
//void remapSamples(vector<Sample*>*, vector<Marker*>*, vector<int>*, int);
void flipStrand(vector<Marker*>*);
void printFamilies(vector<Family*>*);
string descendTree(Sample*, int);
int descendTree2(Sample*, int);
map<int, vector<Sample*> > descendTree3(Sample*, int);
void compileOutputs(vector<Marker*>*, vector<Family*>*, vector<Sample*>*);
void webcheck(vector<string>, map<string, vector<string> >);
map<string, vector<string> > getBatchArgs(string);


/*
 *Function: Initialize()
 *Return: void
 *Description:
 *Initializes the batch file options enumeration.
 */
void Initialize(){
	s_mapStepValues["examplemodule"] = e_examplemodule;
	s_mapStepValues["marker-geno-eff"] = e_marker_geno_eff;
	s_mapStepValues["sample-geno-eff"] = e_sample_geno_eff;
	s_mapStepValues["family-geno-eff"] = e_family_geno_eff;
	s_mapStepValues["allele-freq"] = e_allele_freq;
	s_mapStepValues["mendelian-error"] = e_mendelian_error;
	s_mapStepValues["hw"] = e_hw;
	s_mapStepValues["gender-error"] = e_gender_error;
	s_mapStepValues["tdt"] = e_tdt;
	s_mapStepValues["output-grr"] = e_grr_output;
	s_mapStepValues["output-ped"] = e_ped_output;
	s_mapStepValues["chisquare"] = e_chisquare;
	s_mapStepValues["output-structure"] = e_structure_output;
	s_mapStepValues["output-phase"] = e_phase_output;
	s_mapStepValues["output-beagle"] = e_beagle_output;
	s_mapStepValues["output-superlink"] = e_superlink_output;
	s_mapStepValues["homozygous"] = e_homozygous;
	s_mapStepValues["ld"] = e_ld;
	s_mapStepValues["output-powermarker"] = e_powermarker_output;
	s_mapStepValues["deletions"] = e_deletions;
	s_mapStepValues["concordance"] = e_concordance;
	s_mapStepValues["mito-check"] = e_mitocheck;
	s_mapStepValues["output-lapis"] = e_lapis_output;
	s_mapStepValues["output-fbat"] = e_fbat_output;
	s_mapStepValues["output-mdr"] = e_mdr_output;
	s_mapStepValues["output-qtdt"] = e_qtdt_output;
	s_mapStepValues["output-pdt2"] = e_pdt2_output;
}

int
main (int argc, char* argv[])
{

	Initialize();
	int myrank, nprocs, resultlen;
	char* hostname;
	myrank = 0;
	nprocs = 0;
#ifdef _USEMPI_
	//set up MPI stuff
	MPI::Init(argc, argv);
	hostname = new char[MPI_MAX_PROCESSOR_NAME];
	myrank = MPI::COMM_WORLD.Get_rank();
	nprocs = MPI::COMM_WORLD.Get_size();

	MPI::Get_processor_name(hostname, resultlen);
#endif
	STEPS steps;
	ORDER proc_order;
	vector<string> usemarkers;
	vector<string> usecenters;
	vector<string> usechroms;
	vector<long> usembrange;
	

	//Parse command line arguments
	//assume first argument is -h, -S or the batch file
	if(argc >= 2){
	   	string arg = argv[1];
	   	if(arg == "-h"){
			print_help();
			exit(1);
		}
		if(arg == "-S"){
			steps = initializeSteps();
			print_steps(steps);
			exit(1);
		}
		vector<string> arguments;
		
		for(int i = 2; i < argc; i++){
			string list = "";
			list.assign(argv[i]);
			arguments.push_back(list);
		}
		vector<string>::iterator niter = find(arguments.begin(), arguments.end(), "-noweb");
		map<string, vector<string> > batchargs = getBatchArgs(arg);
		if(niter == arguments.end()){
			webcheck(arguments, batchargs);
		}

		for(int j = 0; j < arguments.size(); j++){
			//pedfile
			if(arguments[j] == "-ped"){
				if(j + 1 >= arguments.size()){
					cerr << "-ped option requires specifying a filename. (Ex: -ped pedinput.ped).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-ped option requires specifying a filename. (Ex: -ped pedinput.ped).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-ped option missing file name?  Halting." << endl;
					exit(1);
				}
				
				opts::_PEDFILE_ = arguments[++j];
			}
			else if(arguments[j] == "-noweb"){
			}
			else if(arguments[j] == "-missing-geno"){
				if(j + 1 >= arguments.size()){
					cerr << "-missing-geno option requires specifying a missing code. (Ex: -missing-geno N). Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-missing-geno option requires specifying a missing code. (Ex: -missing-geno N). Halting!" << endl;
					exit(1);
				}
				//if(test[0] == '-'){
				//	cerr << "-missing-geno option missing code?  Halting." << endl;
				//	exit(1);
				//}
				opts::_NOCALL_ = arguments[++j];
			}
			else if(arguments[j] == "-pedinfo"){
				if(j + 1 >= arguments.size()){
					cerr << "-pedinfo option requires specifying a filename. (Ex: -pedinfo pedinfo.txt). Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-pedinfo option requires specifying a filename. (Ex: -pedinfo pedinfo.txt). Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-pedinfo option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_PEDINFO_ = arguments[++j];
			}
			else if(arguments[j] == "-tped"){
				if(j + 1 >= arguments.size()){
					cerr << "-tped option requires specifying a filename. (Ex: -tped pedinput.ped).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-tped option requires specifying a filename. (Ex: -tped pedinput.ped).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-tped option missing file name?  Halting." << endl;
					exit(1);
				}
				
				opts::_TPEDFILE_ = arguments[++j];
			}
			else if(arguments[j] == "-tfam"){
				if(j + 1 >= arguments.size()){
					cerr << "-tfam option requires specifying a filename. (Ex: -tfam family.map).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-tfam option requires specifying a filename. (Ex: -tfam family.map).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-tfam option missing file name?  Halting." << endl;
					exit(1);
				}
				
				opts::_FAMFILE_ = arguments[++j];
			}
			//dog flag
			else if(arguments[j] == "-dog"){
				opts::_DOG_ = true;
				opts::_CHRX_ = 39;
				opts::_CHRY_ = 40;
			}
			//map file
			else if(arguments[j] == "-map"){
				if(j + 1 >= arguments.size()){
					cerr << "-map option requires specifying a filename. (Ex: -map pedmap.map).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-map option requires specifying a filename. (Ex: -map pedmap.map).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-map option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_MAPFILE_ = arguments[++j];
			}
			//zero genotype file
			else if(arguments[j] == "-zero-geno-file"){
				if(j + 1 >= arguments.size()){
					cerr << "-zero-geno-file option requires specifying a filename. (Ex: -zero-geno-file myzeros.txt).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-zero-geno-file option requires specifying a filename. (Ex: -zero-geno-file myzeros.txt).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-zero-geno-file option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_ZEROGENOFILE_ = arguments[++j];
			}
			//map descriptive file
			else if(arguments[j] == "-mapdesc"){
				if(j + 1 >= arguments.size()){
					cerr << "-mapdesc option requires specifying a filename. (Ex: -mapdesc pedmap.desc).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-mapdesc option requires specifying a filename. (Ex: -mapdesc pedmap.desc).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-mapdesc option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_MAPDESC_ = arguments[++j];
			}
			//sample descriptive file
			else if(arguments[j] == "-sampdesc"){
				if(j + 1 >= arguments.size()){
					cerr << "-sampdesc option requires specifying a filename. (Ex: -sampdesc sampleinfo.txt).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-sampdesc option requires specifying a filename. (Ex: -sampdesc sampleinfo.txt).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-sampdesc option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_SAMPDESC_ = arguments[++j];
			}
			//exclude samples
			else if(arguments[j] == "-excsamples"){
				if(j + 1 >= arguments.size()){
					cerr << "-excsamples option requires specifying a filename. (Ex: -excsamples excludelist.txt).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-excsamples option requires specifying a filename. (Ex: -excsamples excludelist.txt).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-excsamples option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_EXCSAMPS_ = arguments[++j];
			}
			//exclude markers
			else if(arguments[j] == "-excmarkers"){
				if(j + 1 >= arguments.size()){
					cerr << "-excmarkers option requires specifying a filename. (Ex: -excmarkers excludelist.txt).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-excmarkers option requires specifying a filename. (Ex: -excmarkers excludelist.txt).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-excmarkers option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_EXCMARKERS_ = arguments[++j];
			}
			//exclude families
			else if(arguments[j] == "-excfamilies"){
				if(j + 1 >= arguments.size()){
					cerr << "-excfamilies option requires specifying a filename. (Ex: -excfamilies excludelist.txt).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-excfamilies option requires specifying a filename. (Ex: -excfamilies excludelist.txt).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-excfamilies option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_EXCFAMILIES_ = arguments[++j];
			}
			//include markers
			else if(arguments[j] == "-incmarkers"){
				if(j + 1 >= arguments.size()){
					cerr << "-incmarkers option requires specifying a filename. (Ex: -incmarkers includelist.txt).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-incmarkers option requires specifying a filename. (Ex: -incmarkers includelist.txt).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-incmarkers option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_INCMARKERS_ = arguments[++j];
			}
			//include centers (not working right?)
			else if(arguments[j] == "-inccenters"){
				if(j + 1 >= arguments.size()){
					cerr << "-map option requires specifying a filename. (Ex: -map pedmap.map).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-map option requires specifying a filename. (Ex: -map pedmap.map).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-map option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_INCCENTERS_ = arguments[++j];
			}
			//include samples
			else if(arguments[j] == "-incsamples"){
				if(j + 1 >= arguments.size()){
					cerr << "-incsamples option requires specifying a filename. (Ex: -incsamples includelist.txt).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-incsamples option requires specifying a filename. (Ex: -incsamples includelist.txt).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-incsamples option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_INCSAMPLES_ = arguments[++j];
			}
			//include families
			else if(arguments[j] == "-incfamilies"){
				if(j + 1 >= arguments.size()){
					cerr << "-incfamilies option requires specifying a filename. (Ex: -incfamilies includelist.txt).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-incfamilies option requires specifying a filename. (Ex: -incfamilies includelist.txt).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-incfamilies option missing file name?  Halting." << endl;
					exit(1);
				}
				opts::_INCFAMILIES_ = arguments[++j];
			}
			//binary input files based on plink input
			else if(arguments[j] == "-bin-input"){
				if(j + 1 >= arguments.size()){
					cerr << "-bin-input option requires specifying a filename prefix. (Ex: -bin-input data).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-bin-input option requires specifying a filename prefix. (Ex: -bin-input data).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-bin-input option missing filename prefix?  Halting." << endl;
					exit(1);
				}
				opts::_BINPREFIX_ = arguments[++j];
			}
			//output binary genotype file set
			else if(arguments[j] == "-make-bin"){
				if(j + 1 >= arguments.size()){
					cerr << "-make-bin option requires specifying a filename prefix. (Ex: -make-bin data).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-make-bin option requires specifying a filename prefix. (Ex: -make-bin data).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-make-bin option missing filename prefix?  Halting." << endl;
					exit(1);
				}
				opts::_BINPREFIX_ = arguments[++j];
				opts::_MAKEBIN_ = true;
			}
			//out option...sets global prefix for output files
			else if(arguments[j] == "-out"){
				if(j + 1 >= arguments.size()){
					cerr << "-out option requires specifying a filename prefix. (Ex: -out run1).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-out option requires specifying a filename prefix. (Ex: -out run1).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-out option missing filename prefix?  Halting." << endl;
					exit(1);
				}
				opts::_OUTPREFIX_ = arguments[++j];
				opts::_OUTPREFIX_ += ".";
			}
			//enable micro satellites
			else if(arguments[j] == "-micro-sats"){
				opts::_MICROSATS_ = true;
			}
			//flip strands
			else if(arguments[j] == "-flip"){
				opts::_FLIPSTRAND_ = true;

				if(j + 1 >= arguments.size()){
					cerr << "-flip option requires specifying a filename. (Ex: -flip strand_info.txt).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-flip option requires specifying a filename. (Ex: -flip strand_info.txt).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-flip option missing filename?  Halting." << endl;
					exit(1);
				}
				opts::_FLIPFILE_ = arguments[++j];
			}
			//specify specific chromosome
			else if(arguments[j] == "-chrom"){
				if(j + 1 >= arguments.size()){
					cerr << "-chrom option requires specifying a chromosome. (Ex: -chrom 12).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-chrom option requires specifying a chromosome. (Ex: -chrom 12).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-chrom option missing value?  Halting." << endl;
					exit(1);
				}
				opts::_CHROM_LIMIT_ = true;
				string val = arguments[++j];
				if(val != "X" && val != "x" && val != "Y" && val != "y"){
					opts::_CHROM_ = atoi(val.c_str());
				}
				else if(val == "X" || val == "x"){
					opts::_CHROM_ = opts::_CHRX_;
					vector<string>::iterator df = find(arguments.begin(), arguments.end(), "-dog");
					if(df != arguments.end()){
						opts::_CHROM_ = 39;
					}
				}
				else if(val == "Y" || val == "y"){
					opts::_CHROM_ = opts::_CHRY_;
					vector<string>::iterator df = find(arguments.begin(), arguments.end(), "-dog");
					if(df != arguments.end()){
						opts::_CHROM_ = 40;
					}
				}
				else if(val == "Xy" || val == "XY" || val == "xY" || val == "xy"){
					opts::_CHROM_ = opts::_CHRXY_;
					vector<string>::iterator df = find(arguments.begin(), arguments.end(), "-dog");
					if(df != arguments.end()){
						opts::_CHROM_ = 41;
					}
				}
			}
			//specify bploc min
			else if(arguments[j] == "-bp-min"){

				if(j + 1 >= arguments.size()){
					cerr << "-bp-min option requires specifying a positive number. (Ex: -bp-min 1234567).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-bp-min option requires specifying a positive number. (Ex: -bp-min 1234567).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-bp-min option missing value?  Halting." << endl;
					exit(1);
				}
				opts::_BP_LOW_LIMIT_ = true;
				opts::_BP_LOW_ = atoi(arguments[++j].c_str());
			}
			//specify bploc max
			else if(arguments[j] == "-bp-max"){
				if(j + 1 >= arguments.size()){
					cerr << "-bp-max option requires specifying a value. (Ex: -bp-max 987654321).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-bp-max option requires specifying a value. (Ex: -bp-max 987654321).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-bp-max option missing value?  Halting." << endl;
					exit(1);
				}
				opts::_BP_HIGH_LIMIT_ = true;
				opts::_BP_HIGH_ = atoi(arguments[++j].c_str());
			}
			//specify bploc spacing
			else if(arguments[j] == "-bp-space"){
				if(j + 1 >= arguments.size()){
					cerr << "-bp-space option requires specifying a value. (Ex: -bp-space 1000000).  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-bp-space option requires specifying a value. (Ex: -bp-space 1000000).  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-bp-space option missing value?  Halting." << endl;
					exit(1);
				}
				opts::_BP_SPACE_LIMIT_ = true;
				opts::_BP_SPACE_ = atoi(arguments[++j].c_str());
			}
			//enable predefined minor allele frequencies
			else if(arguments[j] == "-freq-file"){
				if(j + 1 >= arguments.size()){
					cerr << "-freq-file option requires specifying a filename. (Ex: -freq-file mymaffile.txt.  Halting!" << endl;
					exit(1);
				}
				string test = arguments[j+1];
				if(test.length() == 0){
					cerr << "-freq-file option requires specifying a filename. (Ex: -freq-file mymaffile.txt.  Halting!" << endl;
					exit(1);
				}
				if(test[0] == '-'){
					cerr << "-freq-file option requires specifying a filename. (Ex: -freq-file mymaffile.txt.  Halting!" << endl;
					exit(1);
				}
				opts::_FREQ_FILE_EXISTS_ = true;
				opts::_FREQ_FILE_ = arguments[++j];
			}
			//enable output of remaining samples, markers, families
			else if(arguments[j] == "-remaining"){
				opts::_OUTPUTREMAINS_ = true;
			}
			//print family structures and exit
			else if(arguments[j] == "-print-families"){
				opts::_PRINTFAMS_ = true;
			}
			//disable summary file production
			else if(arguments[j] == "-no-summary"){
				opts::_COMPILE_OUTPUTS_ = false;
			}
			else{
				cout << "Argument: " << arguments[j] << " unknown.  Exitting.\n";
				exit(1);
			}
		}
		string logoutput = opts::_OUTPREFIX_ + "wasp.log";
//		ofstream log(logoutput.c_str());
		opts::_LOG_.open(logoutput.c_str());
		if(!opts::_LOG_){
			cerr << "Error opening wasp.log.  Exiting!\n";
			exit(1);
		}
		string args = "";
		for(int j = 0; j < argc; j++){
			string t = argv[j];
			args += t + " ";
//			log << argv[j] << " ";
		}
		args += "\n";
		opts::printLog(args);
//		log << endl;
//		log.close();	

		if(opts::_PEDFILE_.length() == 0 && opts::_BINPREFIX_.length() == 0 && opts::_TPEDFILE_.length() == 0){
			opts::printLog("No pedfile specified and no binary input files specified.  Exiting!\n");
			exit(0);
		}
		error_check();
		//steps = initializeSteps();
		if(!opts::_MAKEBIN_){
			proc_order = parseInput(arg);//, &steps);
		}
		
	}
	else{
		usage();
		exit(1);
	}
		//Begin batch file steps and reading in data
		startProcess(&proc_order, NULL, myrank, nprocs, usemarkers, usecenters);
#ifdef _USEMPI_
	MPI::Finalize();
#endif
    return 1;
}

/*
 *Function: error_check
 *Return: void
 *Description: 
 *Performs basic error checking on completeness of command line arguments.
 *
 */
void error_check(){
	if((opts::_BP_LOW_LIMIT_ || opts::_BP_HIGH_LIMIT_) && !opts::_CHROM_LIMIT_){
		opts::printLog("-chrom is required when specifying -bp-min and/or -bp-max.  Exiting!\n");
		exit(0);
	}
}

/*depricated*/
vector<string> readMarkers(string f){
	ifstream input;
	input.open(f.c_str(), ios::in);

	if(!input){
		cerr << "Error opening marker file: " << f << endl;
		exit(1);
	}

	string rsid = "";
	vector<string> mrks;
	while(input >> rsid){
		mrks.push_back(rsid);
	}
	input.close();
	
	return mrks;
}

/*depricated*/
vector<string> readCenters(string f){
	ifstream input;
	input.open(f.c_str(), ios::in);
	if(!input){
		cerr << "Error opening center file: " << f << endl;
		exit(1);
	}

	string center = "";
	vector<string> cntrs;
	while(input >> center){
		cntrs.push_back(center);
	}

	input.close();
	
	return cntrs;
}

/*
 *Function: print_help
 *Return: void
 *Description:
 *Outputs general help statement when -h option is used.
 *
 */
void print_help(){
	int w = 5;
	cout << "wasp - Whole-Genome Associstion Study Pipeline (Version: " << opts::_WASPVER_ << ")" << endl << endl;
	cout << "Center for Human Genetics Research" << endl;
	cout << "Vanderbilt University Medical Center" << endl << endl;
	cout << "!!!!  For most up-to-date options and details  !!!!" << endl;
	cout << "!!!!  visit http://chgr.mc.vanderbilt.edu/wasp !!!!" << endl << endl;
	cout << "usage: wasp <batchfile> [<option1> <option2>...]" << endl << endl;
	cout << "   Special options: " << endl;
	cout << "     -h                 Prints this help" << endl;
	cout << "     -S                 Displays list of valid processing steps" << endl << endl;
	cout << "   Options: " << endl;
	cout << "     -ped {pedfile}     Specify pedfile <fam> <ind> <dad> <mom> <sex> <aff> <genotypes...>" << endl;
	cout << "     -map {mapfile}     Specify mapfile <chrom> <rsid> <bploc>" << endl << endl;
	cout << "     -tped {tpedfile}   Specify tpedfile <chr> <snp> <cm> <bploc> <genotypes...>" << endl;
	cout << "     -tfam {famfile}    Specify family map file for transposed set" << endl;
	cout << "              Format: <fam> <ind> <dad> <mom> <gend> <aff>" << endl;
	cout << "     -bin-input {file_prefix} Read inputs in binary format (based on -make-bin output) (optional)" << endl;
	cout << "     -micro-sats        Use if maker(s) have more than 2 alleles." << endl;
	cout << "     -mapdesc {mapdescriptive_file} Specify more details about map file (optional)" << endl;
	cout << "                                    Mapfile descriptor format: Chrom Marker Basepair_location Other_marker_id Enzyme" << endl;
	cout << "     -sampdesc {sample_Descriptive_file} Specify more details about samples (optional)" << endl;
	cout << "                                         Sample descriptor format: Family Ind Center Plate Well" << endl << endl;
	cout << "     -excsamples {samplefile} File containing samples to exclude (optional)" << endl;
	cout << "     -excmarkers {markerfile} File containing markers to exclude (optional)" << endl << endl;
	cout << "     -incsamples {samplefile} File containing samples to include (optional)" << endl;
	cout << "     -incmarkers {markerfile} File containing markers to include (optional)" << endl << endl;
	cout << "     -make-bin {name}         Convert Ped/Map files into binary input file format with {name} as the file prefix (optional)" << endl;
	cout << endl;
	cout << "     -out {name}         Prepends all output files with {name} (optional)" << endl;
	cout << endl;
	cout << "     -freq-file {frequencyfile} File containing predefined minor allele frequencies (optional)" << endl;
	cout << "                                Format: Marker_id Minor_allele_frequency" << endl;
	cout << "     -flip {filename}  Specifies a list of markers to flip to the opposite strand (optional)" << endl;
	cout << endl;
	cout << "     -chrom {chrom}  Limit markers read to specified chromosome (optional)" << endl;
	cout << "     -bp-min {number}  Specifies minimum chromosome specific base-pair location to load (optional)" << endl;
	cout << "     -bp-max {number} Specifies maximum chromosome specific base-pair location to load (optional)" << endl;
	cout << "     -bp-space {number} Specifies to load only markers spaced by {number} base-pair across all chromosomes. (optional)" << endl;
	cout << endl;
	cout << endl;
	cout << "   Step specific options (to be used in the batch file):" << endl;
	cout << "     -chrom {chrom}  Valid on all steps." << endl;
	cout << "     -bp-min {number}  Valid on all steps." << endl;
	cout << "     -bp-max {number}  Valid on all steps." << endl;
	cout << "     -bp-space {number}  Valid on all steps." << endl;
	cout << "     -thresh-markers-max {number}  Valid on marker-geno-eff, family-geno-eff, allele-freq, mendelian-error, hw, gender-error, tdt, chisquare.." << endl;
	cout << "     -thresh-markers-min {number}  Valid on marker-geno-eff, family-geno-eff, allele-freq, mendelian-error, hw, gender-error, tdt, chisquare." << endl;
	cout << "     -thresh-samples-max {chrom}  Valid on sample-geno-eff, mendelian-error." << endl;
	cout << "     -thresh-samples-min {chrom}  Valid on sample-geno-eff, mendelian-error." << endl;
	cout << "     -thresh-families-max {chrom}  Valid on mendelian-error." << endl;
	cout << "     -thresh-families-min {chrom}  Valid on mendelian-error." << endl;
	cout << "     -zero  Valid on mendelian-error." << endl;
	cout << "     -zero-l2  Valid on mendelian-error." << endl;
	cout << "     -zero-l2-fams Valid on mendelian-error." << endl;
	cout << "     -penetrance-file {filename} Valid on output-superlink." << endl;
	cout << "     -rm-mono  Valid on allele-freq." << endl;
	cout << "     -rm-het-only Valid on allele-freq." << endl;
	cout << "     -center-file {filename}  Valid on output-lapis." << endl;
	cout << "     -out {fileprefix}  Valid on all steps." << endl;
	cout << "     -covar-file {filename}" << endl;
	cout << "     -covars-name {list of headers}" << endl;
	cout << "     -covars-number {list of header numbers}" << endl;
	cout << "     -trait-file {filename}" << endl;
	cout << "     -traits-name {list of headers}" << endl;
	cout << "     -traits-number {list of header numbers}" << endl;
	cout << "     -covar-missing {missing value}" << endl;
	cout << "     -trait-missing {missing value}" << endl;
	cout << "     -group-file {filename}" << endl;
	cout << "     -ped {filename}  Valid on concordance" << endl;
	cout << "     -map {filename}  Valid on concordance" << endl;
	cout << "     -tped {filename} Valid on concordance" << endl;
	cout << "     -tfam {filename} Valid on concordance" << endl;
	cout << "     -bin-input {fileprefix} Valid on concordance" << endl;
	cout << "     -strat-file {filename}  Valid on output-structure." << endl;
	cout << "     -parents-only  Valid on output-structure." << endl;
	cout << "     -trio  Valid on output-phase." << endl;
	cout << "     -disease {name}  Valid on output-beagle." << endl;
	cout << "     -homozyg-zeros {number}  Valid on homozygous." << endl;
	cout << "     -homozyg-span {number}  Valid on homozygous." << endl;
	cout << "     -homozyg-seq-prob {number}  Valid on homozygous." << endl;
	cout << "     -homozyg-wgha  Alternate homozygous span algorithm.  Valid on homozygous." << endl;
	cout << "     -homozyg-min-samp {number}  Valid on homozygous." << endl;
	cout << "     -ld-pairwise {number} {number} {number} Valid on ld." << endl;
	cout << "     -ld-vif {number} {number} {number}  Valid on ld." << endl;
	cout << "     -ld-chop  Valid on ld." << endl;
	cout << "     -deletion-span {number}   Valid on deletions." << endl;
	cout << "     -deletion {number}  Valid on deletions." << endl;
	cout << "     -filter-overall  Valid on allele-freq." << endl;
	cout << "     -filter-file  Valid on allele-freq." << endl;
	cout << "     -parental    Provides further breakdown. Valid on allele-freq and hw." << endl;
	cout << "     -gender      Provides further breakdown. Valid on allele-freq and hw." << endl;
	cout << "     -casecontrol Provides further breakdown. Valid on allele-freq and hw." << endl;
	cout << "     -all  Override overall calc to include all samples.  Valid on allele-freq and hw." << endl;
	cout << "     -all-children  Override to include only all children.  Valid on allele-freq and hw." << endl;
	cout << "     -random-child  Override to include one random child from each family." << endl;
	cout << "                    Valid on allele-freq and hw." << endl;
	cout << "     -unk-spouses   Override to use only spouses with unknown phenotype." << endl;
	cout << "                    Valid on allele-freq and hw." << endl;
	cout << "     -unaff-spouses-only  Override to use only spouses with unaffected phenotype." << endl;
	cout << "                    Valid on allele-freq and hw." << endl;
	cout << "     -no-summary    Disables production of marker, sample, family summary files." << endl;
	cout << endl;
	cout << endl << endl;
	cout << "BATCH FILE FORMATTING:" << endl << endl;
	cout << "The batch file consists of a list of steps and options for the process" << endl;
	cout << "to follow in order of appearance in the file." << endl << endl;
	cout << "Example batch file:" << endl;
	cout << "marker-geno-eff -thresh-markers-min 90" << endl;
	cout << "allele-freq -thresh-markers-min 0.1" << endl;
	cout << endl;
	
}

/*
 *Function: print_steps
 *Return: void
 *Parameters: STEPS object
 *Description:
 *Outputs valid batch file steps and their descriptions when the -S command line argument is used
 */
void print_steps(STEPS s){
	STEPS::iterator s_iter;
	int field = 0;
	for(s_iter = s.begin(); s_iter != s.end(); s_iter++){
		Step mystep =  s_iter->second;
		if(s_iter->first.size() > field){
			field = s_iter->first.size();
		}
	}
	field+=5;
	cout << left << setw(field) << "Step:" << "Description:" << endl;
	cout << left << setw(field) << "----------" << "---------------" << endl;
	for(s_iter = s.begin(); s_iter != s.end(); s_iter++){
		Step mystep =  s_iter->second;
		cout << left << setw(field) << s_iter->first <<  mystep.getName() << endl;
		//cout << "\tThresh: " << mystep.getThreshold() << endl;
	}
}

/*
 *Function: usage
 *Description:
 *outputs general usage statement
 */
void usage(){
	cout << "wasp - Whole-Genome Associstion Study Pipeline (Version: " << opts::_WASPVER_ << ")" << endl << endl;
	cout << "usage: wasp <batchfile> [<option1> <option2>....]" << endl
		 << endl;
	cout << "For a list of valid steps to be inserted into the batch file:" << endl;
	cout << "\t\t> wasp -S" << endl << endl;
	cout << "For help: " << endl;
	cout << "\t\t> wasp -h" << endl << endl;
}

/*
 *Function: startProcess
 *Parameters:
 *ORDER object
 *database connection (null)
 *rank
 *
 *Description:
 *Starts wasp process
 */
void startProcess(ORDER* order, void* con, int myrank, int nprocs, vector<string> usemarkers, vector<string> usecenters){
	ORDER::iterator o_iter;
	//Families* families = new Families();
	//Markers*  markers = new Markers();

/*
	if(usemarkers.size() > 0 && myrank == 0){
		markers->setUse(usemarkers);
	}
	if(usecenters.size() > 0 && myrank == 0){
		families->setCenters(usecenters);
	}
	if(_RESTRICT_ && myrank == 0){
		markers->setRestrict();
	}
	if(_RESTRICT_SAMPS_ && myrank == 0){
		families->setRestrict();
	}
	if(_CASECONTROL_ && myrank == 0){
		families->setCaseControl();
		markers->setCaseControl();
	}
*/


	vector<Family*> families;
	vector<Sample*> samples;
	vector<Marker*> markers;	
	vector<int> marker_map;
#ifndef _USEMPI_
		//read zerogenofile
		if(opts::_ZEROGENOFILE_.length() > 0){
			readZeroGenoFile();
		}
		//read Binary input
		if(opts::_BINPREFIX_.length() > 0 && !opts::_MAKEBIN_){
			if(opts::_PEDINFO_.length() > 0){
				opts::printLog("Reading Pedigree informatino file: " + opts::_PEDINFO_ + "\n");
				readPedInfo();
			}
			opts::printLog("Reading data using Binary inputs: prefix: " + opts::_BINPREFIX_ + "\n");
			if(opts::_MICROSATS_){
				opts::printLog("You specified microsatellite markers exist.  This option is not compatible with binary input files.  Please use the standard PED file format for your input when using microsatellite markers.  Exiting...\n");
				exit(1);
			}
			readBin(&samples, &families, &markers, &marker_map);
			assignLinks(&families);
			reorderAlleles(&samples, &markers);
			if(opts::_FLIPSTRAND_){
				flipStrand(&markers);
			}
		}
		//read ped & map file
		else if(opts::_PEDFILE_.length() > 0 && opts::_MAPFILE_.length() > 0 && opts::_TPEDFILE_.length() == 0){
			if(opts::_MAPFILE_.length() > 0){
				opts::printLog("Reading MAP file: " + opts::_MAPFILE_ + "\n");
				readMap(&markers, &marker_map);
			}
			if(opts::_PEDFILE_.length() > 0){
				if(opts::_PEDINFO_.length() > 0){
					opts::printLog("Reading Pedigree informatino file: " + opts::_PEDINFO_ + "\n");
					readPedInfo();
				}
				opts::printLog("Reading PED file: " + opts::_PEDFILE_ + "\n");
				readPed(&samples, &families, &markers, &marker_map);
				assignLinks(&families);
				reorderAlleles(&samples, &markers);
				if(opts::_FLIPSTRAND_){
					flipStrand(&markers);
				}
			}
		}
		else if(opts::_PEDFILE_.length() == 0 && opts::_TPEDFILE_.length() > 0 && opts::_FAMFILE_.length() > 0){
			if(opts::_FAMFILE_.length() > 0){
				if(opts::_PEDINFO_.length() > 0){
					opts::printLog("Reading Pedigree informatino file: " + opts::_PEDINFO_ + "\n");
					readPedInfo();
				}
				opts::printLog("Reading Pedigree file: " + opts::_FAMFILE_ + "\n");
				readTFam(&samples, &families);
				assignLinks(&families);
			}
			if(opts::_TPEDFILE_.length() > 0){
				opts::printLog("Reading TPED file: " + opts::_TPEDFILE_ + "\n");
				readTPed(&markers, &samples, &marker_map);
				reorderAlleles(&samples, &markers);
				if(opts::_FLIPSTRAND_){
					flipStrand(&markers);
				}
			}
		}
		else if(!opts::_DBINPUT_){
			cerr << "No input method specified!" << endl;
			exit(1);
		}
		//int goodsamps = 0;
		//int goodfams = 0;
		//int goodmarkers = 0;
		for(int i = 0; i < samples.size(); i++){
			Sample* samp = samples[i];
			if(samp->isEnabled()){
				//goodsamps++;
				opts::_SAMPLES_WORKING_++;
			}
		}
		for(int i = 0; i < families.size(); i++){
			Family* fam = families[i];
			if(fam->isEnabled()){
				//goodfams++;
				opts::_FAMILIES_WORKING_++;
			}
		}
		for(int i = 0; i < markers.size(); i++){
			Marker* mark = markers[i];
			if(mark->isEnabled()){
				//goodmarkers++;
				opts::_MARKERS_WORKING_++;
			}
		}
		string text = "Markers found: ";
		text += getString<int>(opts::_MARKERS_FOUND_); 
		text += ", "; 
		text += getString<int>(opts::_MARKERS_WORKING_);
	   	text += " are enabled!\n";
		opts::printLog(text);
		text = "Samples found: ";
	   	text += getString<int>(opts::_SAMPLES_FOUND_);
	   	text += ", ";
	   	text += getString<int>(opts::_SAMPLES_WORKING_);
	   	text += " are enabled!\n";
		opts::printLog(text);
		text = "Families found: ";
	    text += getString<int>(opts::_FAMILIES_FOUND_);
	   	text += ", ";
	   	text += getString<int>(opts::_FAMILIES_WORKING_);
	    text += " are enabled!\n";
		int founders = 0;
		int nonfounders = 0;
		for(int f = 0; f < families.size(); f++){
			Family* fam = families[f];
			if(fam->isEnabled()){
				vector<Sample*>* founder = fam->getFounders();
				vector<Sample*>* nonfound = fam->getNonFounders();
				for(int s = 0; s < founder->size(); s++){
					Sample* samp = (*founder)[s];
					if(samp->isEnabled()){
						founders++;
					}
				}
				for(int s = 0; s < nonfound->size(); s++){
					Sample* samp = (*nonfound)[s];
					if(samp->isEnabled()){
						nonfounders++;
					}
				}
			}
		}
		text += "Founders found: ";
		text += getString<int>(founders);
		text += "\n";
		opts::printLog(text);
		opts::printLog("Non-founders found: " + getString<int>(nonfounders) + "\n");

		if(opts::zerogenoinfo.size() > 0){
			zeroSingleGenos(&markers, &samples);
		}

		if(opts::_MAKEBIN_){
			writeBit(&samples, &families, &markers, &marker_map);
			exit(1);
		}
		if(opts::_PRINTFAMS_){
			printFamilies(&families);
			exit(1);
		}
	//}		
	//set_me_up->summary();
#endif
		
#ifdef _USEMPI_
	//set up MPI structure datatypes
	MPI_Datatype famtype, enzymetype, indtype, markertype, alleletype, chrefftype;
	MPI_Datatype ftype[2] = {MPI_INT, MPI_INT};
	int fblock[2] = {1, 1};
	MPI_Aint fdisp[2] = {0, sizeof(int)};
	MPI_Type_struct(2, fblock, fdisp, ftype, &famtype);
	MPI_Type_commit(&famtype);

	MPI_Datatype etype[5] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_CHAR};
	int eblock[5] = {1, 1, 1, 1, 10};
	MPI_Aint edisp[5] = {
		0,
		sizeof(int),
		2*sizeof(int),
		3*sizeof(int),
		4*sizeof(int),
	};
	MPI_Type_struct(5, eblock, edisp, etype, &enzymetype);
	MPI_Type_commit(&enzymetype);
	
	MPI_Datatype itype[21] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_INT, MPI_CHAR, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_CHAR, MPI_CHAR, MPI_CHAR};
	int iblock[21] = {1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 21, 6};
	MPI_Aint idisp[21] = {
		0, 
		sizeof(int), 
		2*sizeof(int), 
		3*sizeof(int), 
		4*sizeof(int), 
		5*sizeof(int), 
		(5*sizeof(int))+sizeof(float), 
		(6*sizeof(int))+sizeof(float), 
		(6*sizeof(int))+sizeof(float)+(2*sizeof(char)), 
		(7*sizeof(int))+sizeof(float)+(2*sizeof(char)), 
		(8*sizeof(int))+sizeof(float)+(2*sizeof(char)), 
		(9*sizeof(int))+sizeof(float)+(2*sizeof(char)), 
		(10*sizeof(int))+sizeof(float)+(2*sizeof(char)), 
		(11*sizeof(int))+sizeof(float)+(2*sizeof(char)), 
		(12*sizeof(int))+sizeof(float)+(2*sizeof(char)),
		(13*sizeof(int))+sizeof(float)+(2*sizeof(char)),
		(14*sizeof(int))+sizeof(float)+(2*sizeof(char)),
		(15*sizeof(int))+sizeof(float)+(2*sizeof(char)),
		(16*sizeof(int))+sizeof(float)+(2*sizeof(char)),
		(16*sizeof(int))+sizeof(float)+(6*sizeof(char)),
		(16*sizeof(int))+sizeof(float)+(27*sizeof(char))
	};
	MPI_Type_struct(21, iblock, idisp, itype, &indtype);
	MPI_Type_commit(&indtype);
	
	MPI_Datatype mtype[30] = {
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_DOUBLE, 
		MPI_DOUBLE, 
		MPI_INT, 
		MPI_FLOAT,
		MPI_FLOAT,
		MPI_INT,
		MPI_INT,
		MPI_INT, //array 50
		MPI_INT, //array 50
		MPI_INT, //array 50
		MPI_INT, //array 50
		MPI_INT, //array 50
		MPI_INT, //array 50
		//MPI_INT, //array 50
		MPI_CHAR,
		MPI_CHAR,
		MPI_CHAR,
		MPI_CHAR,
		MPI_INT};
	int mblock[30] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 50, 50, 50, 50, 50, 50, 4, 30, 3, 30, 1};
	//int mblock[14] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 30, 3};
	//MPI_Aint mdisp[18] = {
	MPI_Aint mdisp[30] = {
		0, //MPI_INT
		sizeof(int), //MPI_INT
		(2*sizeof(int)), //MPI_INT
		(3*sizeof(int)), //MPI_INT
		(4*sizeof(int)), //MPI_INT
		(5*sizeof(int)), //MPI_INT
		(6*sizeof(int)), //MPI_INT
		(7*sizeof(int)), //MPI_INT
		(8*sizeof(int)), //MPI_INT
		(9*sizeof(int)), //MPI_FLOAT
		(9*sizeof(int))+sizeof(float), //MPI_FLOAT
		(9*sizeof(int))+(2*sizeof(float)), //MPI_FLOAT
		(9*sizeof(int))+(3*sizeof(float)), //MPI_DOUBLE
		(9*sizeof(int))+(3*sizeof(float))+sizeof(double), //MPI_DOUBLE
		(9*sizeof(int))+(3*sizeof(float))+(2*sizeof(double)), //MPI_INT
		(10*sizeof(int))+(3*sizeof(float))+(2*sizeof(double)), //MPI_FLOAT
		(10*sizeof(int))+(4*sizeof(float))+(2*sizeof(double)), //MPI_FLOAT
		(10*sizeof(int))+(5*sizeof(float))+(2*sizeof(double)), //MPI_INT
		(11*sizeof(int))+(5*sizeof(float))+(2*sizeof(double)), //MPI_INT
		(12*sizeof(int))+(5*sizeof(float))+(2*sizeof(double)), //MPI_INT ARRAY
		(62*sizeof(int))+(5*sizeof(float))+(2*sizeof(double)), //MPI_INT ARRAY
		(112*sizeof(int))+(5*sizeof(float))+(2*sizeof(double)), //MPI_INT ARRAY
		(162*sizeof(int))+(5*sizeof(float))+(2*sizeof(double)), //MPI_INT ARRAY
		(212*sizeof(int))+(5*sizeof(float))+(2*sizeof(double)), //MPI_INT ARRAY
		(262*sizeof(int))+(5*sizeof(float))+(2*sizeof(double)), //MPI_INT ARRAY
		//(157*sizeof(int))+(3*sizeof(float))+(2*sizeof(double)), //MPI_INT ARRAY
		(312*sizeof(int))+(5*sizeof(float))+(2*sizeof(double)),//MPI_CHAR
		(312*sizeof(int))+(5*sizeof(float))+(2*sizeof(double))+(4*sizeof(char)),//MPI_CHAR
		(312*sizeof(int))+(5*sizeof(float))+(2*sizeof(double))+(34*sizeof(char)), //MPI_CHAR
		(312*sizeof(int))+(5*sizeof(float))+(2*sizeof(double))+(37*sizeof(char)), //MPI_CHAR
		(312*sizeof(int))+(5*sizeof(float))+(2*sizeof(double))+(67*sizeof(char))
		//(7*sizeof(int))+(3*sizeof(float))+(2*sizeof(double)),//MPI_CHAR
		//(7*sizeof(int))+(3*sizeof(float))+(2*sizeof(double))+(30*sizeof(char)) //MPI_CHAR
		};
	MPI_Type_struct(30, mblock, mdisp, mtype, &markertype);
	//MPI_Type_struct(14, mblock, mdisp, mtype, &markertype);
	MPI_Type_commit(&markertype);
	
	MPI_Datatype atype[93] = {
		MPI_INT, 
		MPI_CHAR, 
		MPI_CHAR, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_FLOAT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT, 
		MPI_INT};
	int ablock[93] = {1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Aint adisp[93] = {
		0, //mpi_int
		sizeof(int), //mpi_char
		sizeof(int)+(2*sizeof(char)), //mpi_char
		sizeof(int)+(4*sizeof(char)), //mpi_float
		sizeof(int)+(4*sizeof(char))+sizeof(float), //mpi_float
		sizeof(int)+(4*sizeof(char))+(2*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(3*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(4*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(5*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(6*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(7*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(8*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(9*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(10*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(11*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(12*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(13*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(14*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(15*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(16*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(17*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(18*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(19*sizeof(float)), //mpi_float
		sizeof(int)+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(2*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(3*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(4*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(5*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(6*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(7*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(8*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(9*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(10*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(11*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(12*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(13*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(14*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(15*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(16*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(17*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(18*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(19*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(20*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(21*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(22*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(23*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(24*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(25*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(26*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(27*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(28*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(29*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(30*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(31*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(32*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(33*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(34*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(35*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(36*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(37*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(38*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(39*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(40*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(41*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(42*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(43*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(44*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(45*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(46*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(47*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(48*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(49*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(50*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_int
		(51*sizeof(int))+(4*sizeof(char))+(20*sizeof(float)), //mpi_float
		(51*sizeof(int))+(4*sizeof(char))+(21*sizeof(float)), //mpi_float
		(51*sizeof(int))+(4*sizeof(char))+(22*sizeof(float)), //mpi_float
		(51*sizeof(int))+(4*sizeof(char))+(23*sizeof(float)), //mpi_float
		(51*sizeof(int))+(4*sizeof(char))+(24*sizeof(float)), //mpi_float
		(51*sizeof(int))+(4*sizeof(char))+(25*sizeof(float)), //mpi_float
		(51*sizeof(int))+(4*sizeof(char))+(26*sizeof(float)), //mpi_float
		(51*sizeof(int))+(4*sizeof(char))+(27*sizeof(float)), //mpi_float
		(51*sizeof(int))+(4*sizeof(char))+(28*sizeof(float)), //mpi_float
		(51*sizeof(int))+(4*sizeof(char))+(29*sizeof(float)), //mpi_float
		(51*sizeof(int))+(4*sizeof(char))+(30*sizeof(float)), //mpi_int
		(52*sizeof(int))+(4*sizeof(char))+(30*sizeof(float)), //mpi_int
		(53*sizeof(int))+(4*sizeof(char))+(30*sizeof(float)), //mpi_int
		(54*sizeof(int))+(4*sizeof(char))+(30*sizeof(float)), //mpi_int
		(55*sizeof(int))+(4*sizeof(char))+(30*sizeof(float)), //mpi_int
		(56*sizeof(int))+(4*sizeof(char))+(30*sizeof(float)), //mpi_int
		(57*sizeof(int))+(4*sizeof(char))+(30*sizeof(float)), //mpi_int
		(58*sizeof(int))+(4*sizeof(char))+(30*sizeof(float)), //mpi_int
		(59*sizeof(int))+(4*sizeof(char))+(30*sizeof(float)), //mpi_int
		(60*sizeof(int))+(4*sizeof(char))+(30*sizeof(float)) //mpi_int
	};
	MPI_Type_struct(93, ablock, adisp, atype, &alleletype);
	MPI_Type_commit(&alleletype);
		
	MPI_Datatype ctype[5] = {
		MPI_INT, 
		MPI_INT,  
		MPI_INT, 
		MPI_INT, 
		MPI_CHAR};
	int cblock[5] = {1, 1, 1, 1, 3};
	MPI_Aint cdisp[5] = {
		0,
		sizeof(int),
		2*sizeof(int), 
		3*sizeof(int),
		4*sizeof(int)};
	MPI_Type_struct(5, cblock, cdisp, ctype, &chrefftype);
	MPI_Type_commit(&chrefftype);
	
	
	//obtain initial offering of Families and Markers by MASTER
	if(myrank == 0){
		cout << "My rank is: " << myrank << endl;
		cout << "Obtaining families from DB..." << endl;
		families->fillFromDB(con);
		cout << "Families obtained: " << families->getSize() << endl;
		cout << "Obtaining markers from DB..." << endl;
		markers->fillFromDB(con);
		cout << "Markers obtained: " << markers->getSize() << endl;
		
		fam_data* fams = new fam_data[families->getSize()];
		ind_data* inds = new ind_data[families->getNumInds()];
		//don't need this yet.
		//chr_eff_data* chreff = new chr_eff_data[families->getNumInds() * ];	

		FAM::iterator fam_iter;
		FAM* myfams = families->getList();
		int famcount = 0;
		int indcount = 0;
		for(fam_iter = myfams->begin(); fam_iter != myfams->end(); fam_iter++){
			fams[famcount] = fam_iter->getFamStruct();

			INDS::iterator ind_iter;
			INDS members = fam_iter->getIndList();
			for(ind_iter = members.begin(); ind_iter != members.end(); ind_iter++){

				inds[indcount] = ind_iter->getIndStruct();
//cout << ind_iter->getFamID() << "\t" << ind_iter->getInd() << "\t" << ind_iter->getPlate() << "\t" << ind_iter->getWell() << endl;
//if(indcount < 3){
//cout << "MASTER: " << inds[indcount].famid << "= famid\t" << inds[indcount].ind << "= ind\t" << inds[indcount].sex << "= sex\t" << inds[indcount].mymother << "= mymom\t" << inds[indcount].myfather << "= mydad\t" << inds[indcount].mother << "= mom\t" << inds[indcount].father << "= dad\t" << inds[indcount].child << "=child" << endl;  
//}
				indcount++;
			}
			famcount++;
		}
		
		int num_markers = markers->getSize();
		int num_per_node = (int)(num_markers / (nprocs - 1));
		
		int remain = num_markers;
		int curr_count = 0;	
		MKR* mymarkers = markers->getList();
		MKR::iterator m_iter;
		m_iter = mymarkers->begin();
		int familysize = families->getSize();
		int familyinds = families->getNumInds();

		//MASTER sends SLAVES all family and marker information (INITIAL)
		for(int i = 1; i < nprocs; i++){
			int receive = -1;
			MPI_Status status;
cout << "Master sending FAMILY CUE to: " << i << endl;
			MPI_Send(&SR_FAMILYDATA, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
cout << "Master sending FAMILY SIZE to: " << i << endl;
			MPI_Send(&familysize, 1, MPI_INT, i, SR_FAMILYDATA, MPI_COMM_WORLD);
cout << "Master sending FAMILIES to: " << i << endl;
			MPI_Send(fams, families->getSize(), famtype, i, SR_FAMILYDATA, MPI_COMM_WORLD);
			MPI_Recv(&receive, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD, &status);
cout << "Master sending IND CUE to: " << i << endl;
			MPI_Send(&SR_INDDATA, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
cout << "Master sending IND SIZE to: " << i << endl;
			MPI_Send(&familyinds, 1, MPI_INT, i, SR_INDDATA, MPI_COMM_WORLD);
cout << "Master sending INDS to: " << i << endl;
			MPI_Send(inds, familyinds, indtype, i, SR_INDDATA, MPI_COMM_WORLD);
			MPI_Recv(&receive, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD, &status);
			int remainder = remain - num_per_node;
			if(remainder < num_per_node){
				num_per_node += remainder;
			}
			
			marker_data* mkrs = new marker_data[num_per_node];
			allele_struct_data* alleles = new allele_struct_data[num_per_node];
			//cout << "Master sending NUM_PER_NODE " << num_per_node << " to " << i << endl;
			//MPI_Send(&num_per_node, 1, MPI_INT, i, 10, MPI_COMM_WORLD);	
			
			for(m_iter; m_iter != mymarkers->end(); m_iter++){
				if(curr_count == num_per_node){
					break;
				}
//		if(first){
//			cout << "MASTER CREATING " << m_iter->first << " (" << curr_count << ")" << endl;
//		}		
				mkrs[curr_count] = m_iter->second.getMarkerStruct();
				alleles[curr_count] = m_iter->second.getAlleleStruct();		
				remain--;
				curr_count++;
			}

//	first = false;
			curr_count = 0;
//			int num_per_remain = num_per_node;
//			int send_count = 0;
//			int asend_count = 0;
//			while(num_per_remain > 0){
//				int size = 1000;
//				if((num_per_remain - 1000) < 0){
//					size = num_per_remain;
//				}
cout << "Master sending MARKER CUE to: " << i << endl;
				MPI_Send(&SR_MARKERDATA, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
cout << "Master sending NUM_PER_NODE to: " << i << endl;
				MPI_Send(&num_per_node, 1, MPI_INT, i, SR_MARKERDATA, MPI_COMM_WORLD);
//				MPI_Send(&size, 1, MPI_INT, i, SR_MARKERDATA, MPI_COMM_WORLD);
cout << "Master sending MARKERS to: " << i << endl;
//				marker_data* tempmkrs = new marker_data[size];
//				for(int t = 0; t < size; t++){
//					tempmkrs[t] = mkrs[send_count];
//					send_count++;
//				}
//if(i == 1){
//	cout << "SENDING TO 1, MARKERS -> " << tempmkrs[0].sysprobe << " to " << tempmkrs[size - 1].sysprobe << endl;
//}
				MPI_Send(mkrs, num_per_node, markertype, i, SR_MARKERDATA, MPI_COMM_WORLD);
				MPI_Recv(&receive, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD, &status);
//				MPI_Send(tempmkrs, size, markertype, i, SR_MARKERDATA, MPI_COMM_WORLD);
//				delete(tempmkrs);
cout << "Master sending ALLELE CUE to: " << i << endl;
				MPI_Send(&SR_ALLELEDATA, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
cout << "Master sending ALLELE NUM to: " << i << endl;
				MPI_Send(&num_per_node, 1, MPI_INT, i, SR_ALLELEDATA, MPI_COMM_WORLD);
//				MPI_Send(&size, 1, MPI_INT, i, SR_ALLELEDATA, MPI_COMM_WORLD);
cout << "Master sending ALLELES to: " << i << endl;
//				allele_struct_data* tempallele = new allele_struct_data[size];
//				for(int a = 0; a < size; a++){
//					tempallele[a] = alleles[asend_count];
//					asend_count++;
//				}
				MPI_Send(alleles, num_per_node, alleletype, i, SR_ALLELEDATA, MPI_COMM_WORLD);
				MPI_Recv(&receive, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD, &status);
//				MPI_Send(tempallele, size, alleletype, i, SR_ALLELEDATA, MPI_COMM_WORLD);
//				delete(tempallele);
//			}
			delete [] mkrs;
			delete [] alleles;
cout << "Master sending DONE to: " << i << endl;
			MPI_Send(&DONE, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
		}
		delete [] fams;
		delete [] inds;
	}
#endif
	
	
	//Start processing steps by slaves.  MASTER gathers all data from slaves and does summary work
	for(o_iter = order->begin(); o_iter != order->end(); o_iter++){
		Step current_step = (Step)(*o_iter);//(*steps)[(*o_iter)];
#ifndef _USEMPI_
		opts::printLog("Working on " + current_step.getName() + "\n");
		//current_step.process(con, families, markers);
				current_step.process(&samples, &families, &markers, &marker_map);
			current_step.PrintSummary();
			current_step.filter();
			current_step.FilterSummary();
			current_step.close();
#endif
#ifdef _USEMPI_		
		//if slave
		if(myrank != 0){
			int receive = -1;

			if(families != NULL){
				delete(families);
				families = new Families();
				if(_CASECONTROL_){
					families->setCaseControl();
				}
			}
			if(markers != NULL){
				delete(markers);
				markers = new Markers();
				markers->setCaseControl();
			}
			if(families == NULL){
				families = new Families();
				if(_CASECONTROL_){
					families->setCaseControl();
				}
			}
			if(markers == NULL){
				markers = new Markers();
			}
			
			MPI_Status status;
			//receive markers
			marker_data* mdata;
			allele_struct_data* adata;
			fam_data* fdata;
			ind_data* idata;
			chr_eff_data* cdata;
			while(receive != DONE){
				MPI_Recv(&receive, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received CUE (" << receive << ")" << endl;

				if(receive == SR_MARKERDATA){
					int size;
					MPI_Recv(&size, 1, MPI_INT, 0, SR_MARKERDATA, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received MARKER SIZE (" << size << ")" << endl;
					mdata = new marker_data[size];
					MPI_Recv(mdata, size, markertype, 0, SR_MARKERDATA, MPI_COMM_WORLD, &status);
					int count_recv;
					MPI_Get_count(&status, markertype, &count_recv);
cout << "SLAVE(" << myrank << ") received MARKERS -> " << count_recv << endl;
//cout << "SLAVE (" << myrank << ") FIRST MARKER = " << mdata[0].sysprobe << " LAST MARKER = " << mdata[size - 1].sysprobe << endl;

					MPI_Send(&DONE, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD);
				}
				else if(receive == SR_ALLELEDATA){
					int size;
					MPI_Recv(&size, 1, MPI_INT, 0, SR_ALLELEDATA, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received ALLELE SIZE (" << size << ")" << endl;
					adata = new allele_struct_data[size];
					MPI_Recv(adata, size, alleletype, 0, SR_ALLELEDATA, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received ALLELES" << endl;
					if(sizeof(mdata) > 0){
						for(int i = 0; i < size; i++){
							Marker* newmarker = new Marker(mdata[i], adata[i]);
							markers->addMarker(newmarker);
//if(myrank == 1){
//	cout << "MARKER: adding (" << mdata[i].sysprobe << ") " << newmarker->getSysprobe() << " (" << i << ") zerocount_case = " << mdata[i].zero_calls_case << " zerocount_control = " << mdata[i].zero_calls_control;
//	cout << "  TOTAL: " << markers->getSize() << endl;
//}
							delete newmarker;	
						}
					}
//					if(o_iter == order->begin()){
//cout << "SLAVE(" << myrank << ") generating temporary table TODO_MARKERS_" << myrank << endl;
//						generate_temp_marker_table(con, myrank);
//						fill_temp_marker_table(con, markers, myrank);
//					}
					delete [] adata;
					delete [] mdata;
					MPI_Send(&DONE, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD);
				}
/*				else if(receive == SR_ENZYMEDATA){
					int size;
					MPI_Recv(&size, 1, MPI_INT, 0, SR_ENZYMEDATA, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received ENZYME SIZE (" << size << ")" << endl;
					edata = new enzyme_count[size];
					MPI_Recv(edata, size, enzymetype, 0, SR_ENZYMEDATA, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received ENZYMES" << endl;
					for(int i = 0; i < size; i++){
						Family*
				}
*/
				else if(receive == SR_FAMILYDATA){
					int size;
					MPI_Recv(&size, 1, MPI_INT, 0, SR_FAMILYDATA, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received FAM SIZE (" << size << ")" << endl;
					fdata = new fam_data[size];
					MPI_Recv(fdata, size, famtype, 0, SR_FAMILYDATA, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received FAMS" << endl;
					for(int i = 0; i < size; i++){
						Family* newfamily = new Family(fdata[i]);
						families->addFamily(newfamily);
						delete newfamily;
					}
					delete [] fdata;
					MPI_Send(&DONE, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD);
				}
				else if(receive == SR_INDDATA){
					int size;
					MPI_Recv(&size, 1, MPI_INT, 0, SR_INDDATA, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received IND SIZE (" << size << ")" << endl;
					idata = new ind_data[size];
					MPI_Recv(idata, size, indtype, 0, SR_INDDATA, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received INDS" << endl;
					for(int i = 0; i < size; i++){
						families->addInd(idata[i]);
//						if(myrank == 1){
//							cout << "IND entering: " << i << " - " << idata[i].famid << " - " << idata[i].ind;
//							cout << " TOTAL: " << families->getNumInds() << endl;
//						}
					}
					delete [] idata;
					MPI_Send(&DONE, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD);
				}
				else if(receive == SR_CHREFFDATA){
					int size;
					MPI_Recv(&size, 1, MPI_INT, 0, SR_CHREFFDATA, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received CHR SIZE (" << size << ")" << endl;
					if(size > 0){
						cdata = new chr_eff_data[size];
						MPI_Recv(cdata, size, chrefftype, 0, SR_CHREFFDATA, MPI_COMM_WORLD, &status);
cout << "SLAVE(" << myrank << ") received CHREFF" << endl;

						for(int i = 0; i < size; i += NUMCHROMS){
							for(int j = 0; j < NUMCHROMS; j++){
								families->addChrEff(cdata[i+j]);
							}
						}
						delete [] cdata;
					}
					MPI_Send(&DONE, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD);
				}
				else if(receive == SR_NUMPERNODE){
				}
			}
		}
		
		if(myrank == 0){
			cout << "Working on " << current_step.getName() << endl;
		}
   	//Environment* env = Environment::createEnvironment(
   	//	Environment::DEFAULT);
	//Connection* con = getDBConnection(env);
	//Markers submarkers[];
	//Families* node_fams;
	//*node_fams = *families;
	


		//SLAVE does work
		fam_data* fams;
		ind_data* inds;
		chr_eff_data* chreffs;
		chr_eff_data* chreffs_real;
		marker_data* mkrs;
		allele_struct_data* alldat;
		enzyme_count* enzdata;
		if(myrank != 0){
			cout << myrank << " Families: " << families->getSize() << " Inds: " << families->getNumInds() << " Markers: " << markers->getSize() << endl;
			current_step.setRank(myrank);			

			current_step.process(con, families, markers);

			//SLAVE prepares to send results back to MASTER			
			int indsize = families->getNumInds();
			int famsize = families->getSize();
			int markersize = markers->getSize();
			int chrsize = indsize * NUMCHROMS;
			vector<string> enzymes = markers->getEnzymes();
			int enzymesize = enzymes.size();
			int enzdatapoints = 8;
			
			FAM::iterator fam_iter;
			FAM* myfams = families->getList();

			int famcount = 0;
			int indcount = 0;
			int chrcount = 0;
			int enzymecount = 0;

			fams = new fam_data[famsize];
			inds = new ind_data[indsize];
			chreffs = new chr_eff_data[chrsize];
			enzdata = new enzyme_count[indsize * enzymesize * enzdatapoints];
			for(fam_iter = myfams->begin(); fam_iter != myfams->end(); fam_iter++){
				fams[famcount] = fam_iter->getFamStruct();
				for(int i = 0; i < enzymes.size(); i++){
					enzdata[enzymecount++] = fam_iter->getEnzMEStruct(enzymes[i]);
					enzdata[enzymecount++] = fam_iter->getEnzZeroEffStruct(enzymes[i]);
					enzdata[enzymecount++] = fam_iter->getEnzTotEffStruct(enzymes[i]);
/*if(myrank == 1){
cout << fam_iter->getFamID() << endl;
cout << "--->" << enzdata[enzymecount - 1].famid << "\t" << enzdata[enzymecount - 1].ind << "\t" << enzdata[enzymecount - 1].type << "\t" << enzdata[enzymecount - 1].enzyme << "\t" << enzdata[enzymecount - 1].count << endl;
}*/
				}
				
				INDS::iterator ind_iter;
				INDS members = fam_iter->getIndList();
				for(ind_iter = members.begin(); ind_iter != members.end(); ind_iter++){
					inds[indcount] = ind_iter->getIndStruct();
					for(int i = 0; i < enzymes.size(); i++){
						enzdata[enzymecount++] = ind_iter->getEnzZeroEffStruct(enzymes[i]);
						enzdata[enzymecount++] = ind_iter->getEnzTotalEffStruct(enzymes[i]);
						enzdata[enzymecount++] = ind_iter->getEnzHetStruct(enzymes[i]);
						enzdata[enzymecount++] = ind_iter->getEnzTotGendStruct(enzymes[i]);
						enzdata[enzymecount++] = ind_iter->getEnzMEStruct(enzymes[i]);
/*if(myrank == 1){
cout << ind_iter->getFamID() << "\t" << ind_iter->getInd() << endl;
cout << "--->" << enzdata[enzymecount - 1].famid << "\t" << enzdata[enzymecount - 1].ind << "\t" << enzdata[enzymecount - 1].type << "\t" << enzdata[enzymecount - 1].enzyme << "\t" << enzdata[enzymecount - 1].count << endl;
}*/
					}
					CHEFFS::iterator ch_iter;
					CHEFFS* cheffs = ind_iter->getChEffs();
					for(ch_iter = cheffs->begin(); ch_iter != cheffs->end(); ch_iter++){
						chreffs[chrcount] = ind_iter->getChEffStruct(ch_iter->first);
						chrcount++;
					}
					indcount++;
				}
				famcount++;
			}
			
			chreffs_real = new chr_eff_data[chrcount];
			for(int t = 0; t < chrcount; t++){
				chreffs_real[t] = chreffs[t];
			}
			
			int receive = -1;
			MPI_Status status;
//cout << "SLAVE (" << myrank << ") sending MASTER FAM CUE" << endl;
//			MPI_Send(&SR_FAMILYDATA, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD);
cout << "SLAVE (" << myrank << ") sending MASTER FAM size: " << famsize << endl;
			MPI_Send(&famsize, 1, MPI_INT, 0, SR_FAMILYDATA, MPI_COMM_WORLD);
cout << "SLAVE (" << myrank << ") sending MASTER FAM DATA" << endl;
			MPI_Send(fams, famsize, famtype, 0, SR_FAMILYDATA, MPI_COMM_WORLD);
			MPI_Recv(&receive, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD, &status);

			
//cout << "SLAVE (" << myrank << ") sending MASTER IND CUE" << endl;
//			MPI_Send(&SR_INDDATA, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD);
cout << "SLAVE (" << myrank << ") sending MASTER IND SIZE: " << indsize << endl;
			MPI_Send(&indsize, 1, MPI_INT, 0, SR_INDDATA, MPI_COMM_WORLD);
cout << "SLAVE (" << myrank << ") sending MASTER INDS" << endl;
			MPI_Send(inds, indsize, indtype, 0, SR_INDDATA, MPI_COMM_WORLD);
			MPI_Recv(&receive, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD, &status);
			

//cout << "SLAVE (" << myrank << ") sending MASTER CHR CUE" << endl;
//			MPI_Send(&SR_CHREFFDATA, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD);
cout << "SLAVE (" << myrank << ") sending MASTER CHR SIZE: " << chrcount << endl;
			MPI_Send(&chrcount, 1, MPI_INT, 0, SR_CHREFFDATA, MPI_COMM_WORLD);
			if(chrcount > 0){
cout << "SLAVE (" << myrank << ") sending MASTER CHRS" << endl;
//cout << myrank << " -> " << chreffs[0].chrom << "\t" << chreffs[0].famid << "\t" << chreffs[0].ind << "\t" << chreffs[0].zero << "\t" << chreffs[0].total << endl;
				MPI_Send(chreffs, chrcount, chrefftype, 0, SR_CHREFFDATA, MPI_COMM_WORLD);
			}
			MPI_Recv(&receive, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD, &status);

cout << "SLAVE (" << myrank << ") sending MASTER ENZYME SIZE: " << enzymecount << endl;
			MPI_Send(&enzymecount, 1, MPI_INT, 0, SR_ENZYMEDATA, MPI_COMM_WORLD);
			if(enzymecount > 0){
cout << "SLAVE (" << myrank << ") sending MASTER ENZYME DATA" << endl;
//for(int z = 0; z < 15; z++){
//cout << "SLAVE " << enzdata[z].famid << "\t" << enzdata[z].ind << "\t" << enzdata[z].type << "\t" << enzdata[z].enzyme << "\t" << enzdata[z].count << endl;
//}
				MPI_Send(enzdata, enzymecount, enzymetype, 0, SR_ENZYMEDATA, MPI_COMM_WORLD);
			}
			MPI_Recv(&receive, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD, &status);

			
			
			MKR::iterator my_m_iter;
			MKR* mymarkers = markers->getList();

			mkrs = new marker_data[markersize];
			int markercount = 0;
			for(my_m_iter = mymarkers->begin(); my_m_iter != mymarkers->end(); my_m_iter++){
				mkrs[markercount] = my_m_iter->second.getMarkerStruct();
if(mkrs[markercount].zero_calls_case > 0 || mkrs[markercount].zero_calls_control > 0){

cout << "SENDTOMASTER (" << myrank << ") " << mkrs[markercount].dbsnp_rsid << "\t" << mkrs[markercount].chrom << "\t" << mkrs[markercount].zero_calls_case << "\t" << mkrs[markercount].zero_calls_control << "\t" << mkrs[markercount].total_calls_case << "\t" << mkrs[markercount].total_calls_control << endl;
}

				markercount++;
			}

//cout << "SLAVE (" << myrank << ") sending MASTER MARKER CUE" << endl;
//			MPI_Send(&SR_MARKERDATA, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD);
cout << "SLAVE (" << myrank << ") sending MASTER MARKER SIZE: " << markersize << endl;
			MPI_Send(&markersize, 1, MPI_INT, 0, SR_MARKERDATA, MPI_COMM_WORLD);
cout << "SLAVE (" << myrank << ") sending MASTER MARKER DATA" << endl;
			MPI_Send(mkrs, markersize, markertype, 0, SR_MARKERDATA, MPI_COMM_WORLD);
			MPI_Recv(&receive, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD, &status);
		
			int allelecount = 0;
			alldat = new allele_struct_data[markersize];
			for(my_m_iter = mymarkers->begin(); my_m_iter != mymarkers->end(); my_m_iter++){
				alldat[allelecount] = my_m_iter->second.getAlleleStruct();
				allelecount++;
			}
		
//cout << "SLAVE (" << myrank << ") sending MASTER ALLELE CUE" << endl;
//			MPI_Send(&SR_ALLELEDATA, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD);
cout << "SLAVE (" << myrank << ") sending MASTER ALLELE SIZE: " << markersize << endl;
			MPI_Send(&markersize, 1, MPI_INT, 0, SR_ALLELEDATA, MPI_COMM_WORLD);
cout << "SLAVE (" << myrank << ") sending MASTER ALLELE DATA" << endl;
//cout << "..vals -- First: " << alldat[0].sysprobe << " Last: " << alldat[markersize - 1].sysprobe << endl;
			MPI_Send(alldat, markersize, alleletype, 0, SR_ALLELEDATA, MPI_COMM_WORLD);
			MPI_Recv(&receive, 1, MPI_INT, 0, RECEIVE, MPI_COMM_WORLD, &status);
			
			//delete(fams);
			//delete(inds);
			//delete(mkrs);
			
			//delete(chreffs);
			//delete(alldat);
			
			//delete(families);
			//delete(markers);
//cout << "SLAVE all data deleted!!" << endl;			
		}//end slave send
		else{ //begin MASTER join
			MPI_Status status;
			if(markers != NULL){
				delete(markers);
				markers = NULL;
				markers = new Markers();
			}
			if(markers == NULL){
				markers = new Markers();
			}
			Families* tempfamilies = new Families();
			FAM* orig_fams = families->getList();
			FAM::iterator fiter;
			for(fiter = orig_fams->begin(); fiter != orig_fams->end(); fiter++){
				tempfamilies->addFamily(&(*fiter));
			}
			for(int i = 1; i < nprocs; i++){
				int size;
				MPI_Recv(&size, 1, MPI_INT, i, SR_FAMILYDATA, MPI_COMM_WORLD, &status);
cout << "MASTER received FAM SIZE (" << size << ") from SLAVE (" << i << ")" << endl;
				fam_data* fams = new fam_data[size];
				MPI_Recv(fams, size, famtype, i, SR_FAMILYDATA, MPI_COMM_WORLD, &status);
cout << "MASTER received FAM DATA (" << size << ") from SLAVE (" << i << ")" << endl;
cout << "MASTER iterating families" << endl;
				for(int f = 0; f < size; f++){
					families->merge(fams[f], tempfamilies);
				}
				delete [] fams;
cout << "MASTER done with families" << endl;
				MPI_Send(&DONE, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
		
				MPI_Recv(&size, 1, MPI_INT, i, SR_INDDATA, MPI_COMM_WORLD, &status);
cout << "MASTER received IND SIZE (" << size << ") from SLAVE (" << i << ")" << endl;
				ind_data* inds = new ind_data[size];
				MPI_Recv(inds, size, indtype, i, SR_INDDATA, MPI_COMM_WORLD, &status);
cout << "MASTER received IND DATA (" << size << ") from SLAVE (" << i << ")" << endl;
				for(int f = 0; f < size; f++){
/*cout << endl << endl;
cout << inds[f].famid << endl;
cout << inds[f].ind << endl;
cout << inds[f].gender_errors << endl;
cout << inds[f].zero_calls << endl;
cout << inds[f].total_calls << endl;
cout << inds[f].geno_eff << endl;
cout << inds[f].enabled << endl;
cout << inds[f].sex << endl;
cout << inds[f].mother << endl;
cout << inds[f].father << endl;
cout << inds[f].child << endl;
cout << inds[f].mymother << endl;
cout << inds[f].myfather << endl;
cout << inds[f].mendelian_errors << endl;
cout << inds[f].locale << endl;
cout << inds[f].plate << endl;
cout << inds[f].well << endl;
*/
					families->merge(inds[f], tempfamilies);
				}
				delete [] inds;
				MPI_Send(&DONE, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);

				MPI_Recv(&size, 1, MPI_INT, i, SR_CHREFFDATA, MPI_COMM_WORLD, &status);
cout << "MASTER received CHR SIZE (" << size << ") from SLAVE (" << i << ")" << endl;
				if(size > 0){
					chr_eff_data* chref = new chr_eff_data[size];
					try{
						MPI_Recv(chref, size, chrefftype, i , SR_CHREFFDATA, MPI_COMM_WORLD, &status);
					}catch(MPI::Exception e){
cout << "my error!!" << endl;
						cout << e.Get_error_code() << " -- " << e.Get_error_string() << endl;
						exit(0);
					}
				//	cout << status->Get_error() << endl;
cout << "MASTER received CHR DATA (" << size << ") from SLAVE (" << i << ")" << endl;
					for(int c = 0; c < size; c++){
						if(chref[c].ind > 0){
							families->merge(chref[c], tempfamilies);
//cout << "MASTERRECEIVE (" << i << ") " << chref[c].famid << "\t" << chref[c].ind << "\t" << chref[c].chrom << "\t" << chref[c].zero << "\t" << chref[c].total << endl;
						}
					}
					delete [] chref;
				}
				MPI_Send(&DONE, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);


				MPI_Recv(&size, 1, MPI_INT, i, SR_ENZYMEDATA, MPI_COMM_WORLD, &status);
cout << "MASTER received ENZYME SIZE (" << size << ") from SLAVE (" << i << ")" << endl;
				if(size > 0){
					enzyme_count* enz = new enzyme_count[size];
					MPI_Recv(enz, size, enzymetype, i , SR_ENZYMEDATA, MPI_COMM_WORLD, &status);
					for(int e = 0; e < size; e++){
//if(enz[e].type == ENZ_FME || enz[e].type == ENZ_IME){
//cout << endl << endl;
//cout << enz[e].famid << endl << enz[e].ind << endl << enz[e].type << endl << enz[e].enzyme << endl << enz[e].count << endl;
//}
						families->merge(enz[e], tempfamilies);
					}
					delete [] enz;
				}
				MPI_Send(&DONE, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);

				
				MPI_Recv(&size, 1, MPI_INT, i, SR_MARKERDATA, MPI_COMM_WORLD, &status);
cout << "MASTER received MARKER SIZE (" << size << ") from SLAVE (" << i << ")" << endl;
				marker_data* mkrs = new marker_data[size];
				MPI_Recv(mkrs, size, markertype, i, SR_MARKERDATA, MPI_COMM_WORLD, &status);
cout << "MASTER received MARKER DATA (" << size << ") from SLAVE (" << i << ")" << endl;
				for(int m = 0; m < size; m++){
					markers->merge(mkrs[m]);
//cout << "MASTERRECEIVE (" << i << ") " << mkrs[m].dbsnp_rsid << "\t" << mkrs[m].enzyme << "\t" << mkrs[m].probe << "\t" << mkrs[m].chrom << "\t" << mkrs[m].pval << "\t" << mkrs[m].minor_allele_freq_tdt << endl;
				}
				delete [] mkrs;
				MPI_Send(&DONE, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);

				MPI_Recv(&size, 1, MPI_INT, i, SR_ALLELEDATA, MPI_COMM_WORLD, &status);
cout << "MASTER received ALLELE SIZE (" << size << ") from SLAVE (" << i << ")" << endl;
				allele_struct_data* asd = new allele_struct_data[size];
				MPI_Recv(asd, size, alleletype, i, SR_ALLELEDATA, MPI_COMM_WORLD, &status);
cout << "MASTER received ALLELE DATA (" << size << ") from SLAVE (" << i << ")" << endl;
				for(int a = 0; a < size; a++){
					markers->merge(asd[a]);
				}
				delete [] asd;
				MPI_Send(&DONE, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);

			}
sleep(5);
			delete(tempfamilies);
			cout << "Markers size: " << markers->getSize() << endl;
			cout << "Family size: " << families->getSize() << endl;
			cout << "Ind size: " << families->getNumInds() << endl;
			current_step.updateFamsMarks(families, markers);
			current_step.PrintSummary(con);
			current_step.filter(con);
			current_step.FilterSummary();
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(myrank > 0){
			delete [] fams;
			delete [] inds;
			delete [] chreffs;
			delete [] chreffs_real;
			delete [] enzdata;
			delete [] mkrs;
			delete [] alldat;
			if(families != NULL){
				delete(families);
				families = NULL;
			}
			if(markers != NULL){
				delete(markers);
				markers = NULL;
			}
cout << "SLAVE (" << myrank << ") families & markers cleaned up!" << endl;			
		}
		//MASTER check for next step and resend filtered information back to slaves for next step
		if(myrank == 0){
			o_iter++;
			if(o_iter != order->end()){
				//send updated info
				fam_data* fams = new fam_data[families->getSize()];
				ind_data* inds = new ind_data[families->getNumInds()];
				int chrsize = families->getNumInds() * NUMCHROMS;
				chr_eff_data* chreffs = new chr_eff_data[chrsize];
				chr_eff_data* chreffs_real;

				FAM::iterator fam_iter;
				FAM* myfams = families->getList();
				int famcount = 0;
				int indcount = 0;
				int chrcount = 0;
				for(fam_iter = myfams->begin(); fam_iter != myfams->end(); fam_iter++){
					fams[famcount] = fam_iter->getFamStruct();

					INDS::iterator ind_iter;
					INDS members = fam_iter->getIndList();
					for(ind_iter = members.begin(); ind_iter != members.end(); ind_iter++){
						inds[indcount] = ind_iter->getIndStruct();
						CHEFFS::iterator ch_iter;
						CHEFFS* cheffs = ind_iter->getChEffs();
						for(ch_iter = cheffs->begin(); ch_iter != cheffs->end(); ch_iter++){
							chreffs[chrcount] = ind_iter->getChEffStruct(ch_iter->first);
							chrcount++;
						}
						indcount++;
					}
					famcount++;
				}
				chreffs_real = new chr_eff_data[chrcount];
				for(int t = 0; t < chrcount; t++){
					chreffs_real[t] = chreffs[t];
				}
				
				int num_markers = markers->getSize();
				int num_per_node = (int)(num_markers / (nprocs - 1));
				
				int remain = num_markers;
				int curr_count = 0;	
				MKR* mymarkers = markers->getList();
				MKR::iterator m_iter;
				m_iter = mymarkers->begin();
				int familysize = families->getSize();
				int familyinds = families->getNumInds();

				for(int i = 1; i < nprocs; i++){
					int receive = -1;
					MPI_Status status;

cout << "Master sending FAMILY CUE to: " << i << endl;
					MPI_Send(&SR_FAMILYDATA, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
cout << "Master sending FAMILY SIZE (" << familysize << ") to: " << i << endl;
					MPI_Send(&familysize, 1, MPI_INT, i, SR_FAMILYDATA, MPI_COMM_WORLD);
cout << "Master sending FAMILIES (" << families->getSize() << ") to: " << i << endl;
					MPI_Send(fams, families->getSize(), famtype, i, SR_FAMILYDATA, MPI_COMM_WORLD);
					MPI_Recv(&receive, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD, &status);

cout << "Master sending IND CUE to: " << i << endl;
					MPI_Send(&SR_INDDATA, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
cout << "Master sending IND SIZE (" << familyinds << ") to: " << i << endl;
					MPI_Send(&familyinds, 1, MPI_INT, i, SR_INDDATA, MPI_COMM_WORLD);
cout << "Master sending INDS to: " << i << endl;
					MPI_Send(inds, familyinds, indtype, i, SR_INDDATA, MPI_COMM_WORLD);
					MPI_Recv(&receive, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD, &status);

cout << "Master sending CHREFF CUE to: " << i << endl;
					MPI_Send(&SR_CHREFFDATA, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
cout << "Master sending CHRSIZE (" << chrcount << ") to: " << i << endl;
					MPI_Send(&chrcount, 1, MPI_INT, i, SR_CHREFFDATA, MPI_COMM_WORLD);
					if(chrcount > 0){
cout << "Master sending CHRS to: " << i << endl;
						MPI_Send(chreffs, chrcount, chrefftype, i, SR_CHREFFDATA, MPI_COMM_WORLD);
					}
					MPI_Recv(&receive, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD, &status);
					
					int remainder = remain - num_per_node;
					if(remainder < num_per_node){
						num_per_node += remainder;
					}
			
					marker_data* mkrs = new marker_data[num_per_node];
					allele_struct_data* alleles = new allele_struct_data[num_per_node];
			
					for(m_iter; m_iter != mymarkers->end(); m_iter++){
						if(curr_count == num_per_node){
							break;
						}
						mkrs[curr_count] = m_iter->second.getMarkerStruct();
						alleles[curr_count] = m_iter->second.getAlleleStruct();		
						remain--;
						curr_count++;
					}

					curr_count = 0;
cout << "Master sending MARKER CUE to: " << i << endl;
					MPI_Send(&SR_MARKERDATA, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
cout << "Master sending NUM_PER_NODE to: " << i << endl;
					MPI_Send(&num_per_node, 1, MPI_INT, i, SR_MARKERDATA, MPI_COMM_WORLD);
cout << "Master sending MARKERS to: " << i << endl;
					MPI_Send(mkrs, num_per_node, markertype, i, SR_MARKERDATA, MPI_COMM_WORLD);
					MPI_Recv(&receive, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD, &status);

cout << "Master sending ALLELE CUE to: " << i << endl;
					MPI_Send(&SR_ALLELEDATA, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
cout << "Master sending ALLELE NUM to: " << i << endl;
					MPI_Send(&num_per_node, 1, MPI_INT, i, SR_ALLELEDATA, MPI_COMM_WORLD);
cout << "Master sending ALLELES to: " << i << endl;
					MPI_Send(alleles, num_per_node, alleletype, i, SR_ALLELEDATA, MPI_COMM_WORLD);
					MPI_Recv(&receive, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD, &status);
					delete [] mkrs;
					delete [] alleles;
cout << "Master sending DONE to: " << i << endl;
					MPI_Send(&DONE, 1, MPI_INT, i, RECEIVE, MPI_COMM_WORLD);
				}
				delete [] fams;
				delete [] inds;
				delete [] chreffs;
				delete [] chreffs_real;
			}
			o_iter--;

		}
#endif
	}

	if(myrank == 0){
		if(opts::_OUTPUTREMAINS_){
			Finalize* f = new Finalize();
			f->finish(&markers, &samples, &families);
			delete(f);
		}
		if(opts::_COMPILE_OUTPUTS_ && opts::filenames.size() > 0){
			compileOutputs(&markers, &families, &samples);
		}
	}
	int fsize = families.size();
	int ssize = samples.size();
	int msize = markers.size();
	for(int i = 0; i < msize; i++){
		delete(markers[i]);
	}
	for(int i = 0; i < ssize; i++){
		delete(samples[i]);
	}
	for(int i = 0; i < fsize; i++){
		delete(families[i]);
	}
	
	//cout << myrank << " leaving process function" << endl;
}

/*
 * Function: getBatchArgs
 * Return: map<string, vector<string> > (string = step, vector = all step options, even repeating)
 * Parameters: Batch file
 * Description:
 * Places all steps and step options in a map to check against the web error checker
 */
map<string, vector<string> > getBatchArgs(string f){
	map<string, vector<string> > batchargs;

	ifstream infile;
	infile.open(f.c_str(), ios::in);
	if(!infile){
		cerr << "Error opening batch file: " << f << endl;
		exit(1);
	}
	string line = "";
	while(getline(infile, line)){
		if(line == ""){
			continue;
		}
		while(line.length() > 0 && (line.at(line.length() - 1) == ' ' || line.at(line.length() - 1) == '\t')){
			line.erase(line.length() - 1);
		}
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() == 0){
			continue;
		}
		else if(tokens[0] == "#"){
			continue;
		}
		else{
			string step = tokens[0];
			batchargs[step].push_back("NA");
			for(int i = 1; i < tokens.size(); i++){
				batchargs[step].push_back(tokens[i]);
			}
		}
	}
	infile.close();

	return batchargs;	

}

/*
 *Function: parseInput
 *Return: ORDER object (batch file contents)
 *Parameters: batch file
 *Description:
 *Parses batch file and sets up steps for processing.
 */
ORDER parseInput(string file){
	ifstream inFile;
	ORDER process;

	const char* temp = file.c_str();
	inFile.open(temp, ios::in);
	if(!inFile){
		cerr << "Error opening batch file: " << file << endl;
		exit(1);
	}

//	string step;
//	string thresh;
	int count = 0;	
	string line;
	map<string, int> allsteps;
	while(getline(inFile, line)){
		if(line == ""){
			continue;
		}
		count++;
		bool overwrite = true;
		while(line.length() > 0 && (line.at(line.length() - 1) == ' ' || line.at(line.length() - 1) == '\t')){
			line.erase(line.length() - 1);
		}
		map<string, int>::iterator found = allsteps.find(line);
		if(found != allsteps.end()){
			overwrite = false;
		}
		allsteps[line]++;
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() == 0){
			continue;
		}
		else if(tokens[0][0] == '#'){
			continue;
		}
		else{
			string step = tokens[0];
			string thresh = "";
			for(int i = 1; i < tokens.size(); i++){
				thresh += tokens[i] + " ";
			}
			Step s = initializeSteps(step);
			s.setThreshold(thresh);
			s.setOrder(count);
			s.setOverwrite(overwrite);
			if(s.hasIncExc()){
				opts::_KEEP_EXC_SAMPLES_ = true;
			}
			process.push_back(s);
		}
	}
	inFile.close();

	
	return process;
}

/*
 *Function: initializeSteps
 *Parameter: batch step string
 *Return: Step object
 *Description:
 *Creates a Step object based on each step listed in the batch file
 *
 */
Step initializeSteps(string i){
		Step *newstep = NULL;
		Process* tempproc = NULL;
		map<string, StepValue>::iterator found = s_mapStepValues.find(i);
		if(found == s_mapStepValues.end()){
			opts::printLog("Step: " + i + " not recognized.  Exiting.\n");
			exit(1);
		}
		switch (s_mapStepValues[i]){
			case e_examplemodule:
				newstep = new Step("Example Module", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ExampleModule();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_marker_geno_eff:
				newstep = new Step("Marker Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new MarkerGenoEff();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_sample_geno_eff:
				newstep = new Step("Individual Sample Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new SampleGenoEff();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;

			case e_family_geno_eff:
				newstep = new Step("Family Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PercentByFamily();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_concordance:
				newstep = new Step("Concordance Check", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new Concordance();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
	
			case e_allele_freq:
				newstep = new Step("Allele Frequencies", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new AlleleFrequency();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;

			case e_mendelian_error:
				newstep = new Step("Mendelian Errors", "", true);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new MendelianErrors();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
	
			case e_hw:
				newstep = new Step("Hardy-Weinberg Calculations", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new HWEquilibrium();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
	
			case e_gender_error:
				newstep = new Step("Gender Correctness (using X-chromosome markers)", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new GenderCheck();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			
//			case 8:
//				newstep = new Step("Genotype Quality Score Analysis (Platform Specific)", "", false);
//				if(tempproc != NULL){
//					delete(tempproc);
//				}
//				tempproc = new QualityScore();
//				if(opts::_DBOUTPUT_){
//					tempproc->setDBOUT();
//				}
//				if(opts::_MARKERLIST_){
//					tempproc->setMarkerList();
//				}
//				if(opts::_STRATIFY_){
//					tempproc->setStratify();
//				}
//				//if(opts::_QSFILE_.length() <= 0){
//				//	cerr << "Quality Score file not specified!  Exiting...\n";
//				//	exit(1);
//				//}
//				newstep->setProcess(tempproc);
//				break;

			case e_tdt:
				newstep = new Step("TDT", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new RunTDT();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;

			case e_superlink_output:
				newstep = new Step("Superlink Output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new SuperlinkOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_grr_output:
				newstep = new Step("GRR Output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new GRROutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
//			case 11:
//				newstep = new Step("PED Output", "", false);
//				if(tempproc != NULL){
//					delete(tempproc);
//				}
//				tempproc = new PEDOutput();
//				if(opts::_DBOUTPUT_){
//					tempproc->setDBOUT();
//				}
//				if(opts::_MARKERLIST_){
//					tempproc->setMarkerList();
//				}
//				if(opts::_STRATIFY_){
//					tempproc->setStratify();
//				}
//				newstep->setProcess(tempproc);
//				break;
//			case 12:
//				newstep = new Step("Quality Score Output", "", false);
//				if(tempproc != NULL){
//					delete(tempproc);
//				}
//				tempproc = new QSOutput();
//				if(opts::_DBOUTPUT_){
//					tempproc->setDBOUT();
//				}
//				if(opts::_MARKERLIST_){
//					tempproc->setMarkerList();
//				}
//				if(opts::_STRATIFY_){
//					tempproc->setStratify();
//				}
//				newstep->setProcess(tempproc);
//				break;
			case e_ped_output:
				newstep = new Step("PED file output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PartialOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_chisquare:
				newstep = new Step("Chisquare test", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new CaConChisq();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;
			case e_structure_output:
				newstep = new Step("Create Structure input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new STRUCTOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_phase_output:
				newstep = new Step("Create Phase input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PHASEOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_beagle_output:
				newstep = new Step("Create BEAGLE input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new BEAGLEOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_lapis_output:
				newstep = new Step("Create LAPIS input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new LAPISOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_qtdt_output:
				newstep = new Step("Create QTDT input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new QTDTOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_mdr_output:
				newstep = new Step("Create MDR input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new MDROutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_fbat_output:
				newstep = new Step("Create FBAT input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new FBATOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_homozygous:
				newstep = new Step("Homozygous spans", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new Homozygous();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_ld:
				newstep = new Step("LD calculations", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new LD();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_powermarker_output:
				newstep = new Step("Create PowerMarker input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PowerMarkerOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_deletions:
				newstep = new Step("Deletion Detection", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new Deletions();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_pdt2_output:
				newstep = new Step("Create PDT2 input file.", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PDT2Output();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			case e_mitocheck:
				newstep = new Step("Mitochondrial Error Checking (Chrom 25)", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new MitoCheck();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				break;	
			default:
				opts::printLog("Step case " + i + " not found!! Exiting...\n");
				exit(1);
				break;
		}
	return *newstep;
}

/*Function: initializeSteps
 *Return: STEPS object
 *
 *Description:
 *Returns a list of valid steps.  Used when -S command line argument is specified
 *
 */
STEPS initializeSteps(){
	STEPS steps;
	map<string, StepValue>::iterator s_iter;
    for(s_iter = s_mapStepValues.begin(); s_iter != s_mapStepValues.end(); s_iter++){
		Step *newstep = NULL;
		Process* tempproc = NULL;
		switch (s_mapStepValues[s_iter->first]){
			case e_examplemodule:
				newstep = new Step("Example Module", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new ExampleModule();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_marker_geno_eff:
				newstep = new Step("Marker Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new MarkerGenoEff();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_sample_geno_eff:
				newstep = new Step("Individual Sample Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new SampleGenoEff();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;

			case e_family_geno_eff:
				newstep = new Step("Family Genotyping Efficiency", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PercentByFamily();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_concordance:
				newstep = new Step("Concordance Check", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new Concordance();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
	
			case e_allele_freq:
				newstep = new Step("Allele Frequencies", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new AlleleFrequency();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;

			case e_mendelian_error:
				newstep = new Step("Mendelian Errors", "", true);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new MendelianErrors();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
	
			case e_hw:
				newstep = new Step("Hardy-Weinberg Calculations", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new HWEquilibrium();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
	
			case e_gender_error:
				newstep = new Step("Gender Correctness (using X-chromosome markers)", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new GenderCheck();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
/*			
			case 8:
				newstep = new Step("Genotype Quality Score Analysis (Platform Specific)", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new QualityScore();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				//if(opts::_QSFILE_.length() <= 0){
				//	cerr << "Quality Score file not specified!  Exiting...\n";
				//	exit(1);
				//}
				newstep->setProcess(tempproc);
				steps[i] = *newstep;
				delete newstep;
				break;
*/
			case e_tdt:
				newstep = new Step("TDT", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new RunTDT();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;

			case e_superlink_output:
				newstep = new Step("Create Superlink input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new SuperlinkOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;

				break;
			case e_grr_output:
				newstep = new Step("Create GRR input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new GRROutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;

				break;
/*			case 11:
				newstep = new Step("PED Output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PEDOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[i] = *newstep;
				delete newstep;
				break;
			case 12:
				newstep = new Step("Quality Score Output", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new QSOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[i] = *newstep;
				delete newstep;
				break;
*/
			case e_ped_output:
				newstep = new Step("Create PED file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PartialOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_chisquare:
				newstep = new Step("Chisquare test", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new CaConChisq();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;
			case e_structure_output:
				newstep = new Step("Create Structure input files", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new STRUCTOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				delete newstep;
				break;	
            case e_phase_output:
                newstep = new Step("Create Phase input file", "", false);
                if(tempproc != NULL){
                    delete(tempproc);
                }
                tempproc = new PHASEOutput();
                if(opts::_DBOUTPUT_){
                    tempproc->setDBOUT();
                }
                if(opts::_MARKERLIST_){
                    tempproc->setMarkerList();
                }
                if(opts::_STRATIFY_){
                   tempproc->setStratify();
                }
                newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
               break;
			case e_beagle_output:
				newstep = new Step("Create BEAGLE input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new BEAGLEOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;	
			case e_lapis_output:
				newstep = new Step("Create LAPIS input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new LAPISOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;	
			case e_qtdt_output:
				newstep = new Step("Create QTDT input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new QTDTOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;	
			case e_mdr_output:
				newstep = new Step("Create MDR input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new MDROutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;	
			case e_fbat_output:
				newstep = new Step("Create FBAT input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new FBATOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;	
			case e_homozygous:
				newstep = new Step("Homozygous spans", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new Homozygous();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;	
			case e_ld:
				newstep = new Step("LD calculations", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new LD();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;	
			case e_powermarker_output:
				newstep = new Step("Create PowerMarker input file", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PowerMarkerOutput();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;	
			case e_deletions:
				newstep = new Step("Deletion Detection", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new Deletions();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;	
			case e_pdt2_output:
				newstep = new Step("Create PDT2 input file.", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new PDT2Output();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;	
			case e_mitocheck:
				newstep = new Step("Mitochondrial Error Checking (Chrom 25)", "", false);
				if(tempproc != NULL){
					delete(tempproc);
				}
				tempproc = new MitoCheck();
				if(opts::_DBOUTPUT_){
					tempproc->setDBOUT();
				}
				if(opts::_MARKERLIST_){
					tempproc->setMarkerList();
				}
				if(opts::_STRATIFY_){
					tempproc->setStratify();
				}
				newstep->setProcess(tempproc);
				steps[s_iter->first] = *newstep;
				break;	
			default:
				cerr << "Step case " << s_iter->first<< " not found!! Exiting..." << endl;
				exit(1);
		}
	}
	return steps;
}


/*
 *Function: reorderAlleles
 *Return: none
 *Parameters: vector of samples, vector of markers
 *Description:
 *Calculates the raw allele counts and remaps the sample bool vector 1 to be the 
 *minor allele
 *
 */
/*void reorderAlleles(vector<Sample*>* samples, vector<Marker*>* markers){
	int msize = markers->size();
	int ssize = samples->size();

	for(int m = 0; m < msize; m++){
		int mloc = (*markers)[m]->getLoc();
		if((*markers)[m]->getNumAlleles() == 0){
			(*markers)[m]->addAllele("0");
			(*markers)[m]->addAllele("0");
		}
		if((*markers)[m]->getNumAlleles() == 1){
			(*markers)[m]->addAllele("0");
		}
		if((*markers)[m]->getNumAlleles() <= 2){
			int a1 = 0;
			int a2 = 0;
			for(int s = 0; s < ssize; s++){
				if((*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc)){
					a2+=2;
				}
				else if(!(*samples)[s]->getAone(mloc) && !(*samples)[s]->getAtwo(mloc)){
					a1+=2;
				}
				else if(!(*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc)){
					a1++;
					a2++;
				}
			}
			if(a2 < a1 && (*markers)[m]->getAllele2() != "0"){
				string temp = (*markers)[m]->getAllele1();
				//(*markers)[m]->setAllele1((*markers)[m]->getAllele2());
				//(*markers)[m]->setAllele2(temp);
				(*markers)[m]->resetAllele1((*markers)[m]->getAllele2());
				(*markers)[m]->resetAllele2(temp);
				for(int s = 0; s < ssize; s++){
					if((*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc)){
						(*samples)[s]->addAone(mloc, false);
						(*samples)[s]->addAtwo(mloc, false);
					}
					else if(!(*samples)[s]->getAone(mloc) && !(*samples)[s]->getAtwo(mloc)){
						(*samples)[s]->addAone(mloc, true);
						(*samples)[s]->addAtwo(mloc, true);
					}	
				}
			}
		}
	}
}*/

/*
 *Function: assignLinks
 *Description:
 *Creates the family structure by assigning each sample their respective
 *parents, siblings and children
 *
 */
/*void assignLinks(vector<Family*>* families){
    vector<Family*>::iterator f_iter;

	int fsize = families->size();

	for(int f = 0; f < fsize; f++){
		Family* fam = (*families)[f];
		vector<Sample*>* samps = fam->getSamples();
		bool good = false;
		bool excluded = false;
		int ssize = samps->size();

		for(int s = 0; s < ssize; s++){
			Sample* samp = (*samps)[s];
			if(samp->isEnabled()){
				good = true;
			}
			else if(samp->isExcluded()){
				excluded = true;
			}
			if(samp->getDadID() == "0" && samp->getMomID() == "0"){
				fam->addFounder(samp);
				samp->setFounder(true);
			}
			else{
				fam->addNonFounder(samp);
				if(samp->getDadID() != "0"){
                	vector<Sample*>::iterator temp_iter = find_if(samps->begin(), samps->end(), FindSampleByID(samp->getDadID()));
                	if(temp_iter != samps->end()){
						Sample* dad = (*temp_iter);
	 					samp->setDad(dad);
						Sample* last_child = dad->getLastChild();
						if(last_child != NULL && last_child->getSib() == NULL){
							last_child->setSib(samp);
						}
						dad->addChild(samp);
                	}
				}
				if(samp->getMomID() != "0"){
                	vector<Sample*>::iterator temp_iter = find_if(samps->begin(), samps->end(), FindSampleByID(samp->getMomID()));
                	if(temp_iter != samps->end()){
						Sample* mom = (*temp_iter);
	 					samp->setMom(mom);
						Sample* last_child = mom->getLastChild();
						if(last_child != NULL && last_child->getSib() == NULL){
							last_child->setSib(samp);
						}
						mom->addChild(samp);
                	}
				}
			}
		}
		fam->setEnabled(good);
		if(!good && excluded){
			fam->setExcluded(excluded);
		}
	}
}*/


/*
 *Function: readMap
 *Description:
 *Reads the map file and performs the appropriate inclusion/exclusion of markers
 *
 */
void readMap(vector<Marker*>* markers, vector<int>* marker_map){
	map<string, vector<string> > descinfo;
	vector<string> descheaders;
	//vector<string> exclude;
	//vector<string> include;
	map<string,int> exclude;
	map<string,int> include;
	map<string,float> frequencies;

	//read in frequency file
	if(opts::_FREQ_FILE_EXISTS_){
		opts::printLog("Marker Minor Allele Frequencies being read from: " + opts::_FREQ_FILE_ + "\n");
		ifstream finput;
		finput.open(opts::_FREQ_FILE_.c_str(), ios::in);
		if(!finput){
			opts::printLog("Error opening marker frequency file: " + opts::_FREQ_FILE_ + "\n");
			exit(1);
		}
		string line = "";
		int freqline = 1;
		while(getline(finput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2){
				opts::printLog("Marker frequency file column size != 2 on line: " + line + " Skipping line!!\n");
				continue;
			}
			try{
				float freq = -1.0f;
				for(int c = 0; c < tokens[1].size(); c++){
					if(!isdigit(tokens[1][c]) && tokens[1][c] != '.'){
						throw "Number format exception";
					}
				}
				freq = std::atof(tokens[1].c_str());
				frequencies[tokens[0]] = freq;	
			}catch(...){
				opts::printLog(tokens[1] + " is not a valid value for line: " + getString<int>(freqline) + "\n");
				exit(1);
			}
			freqline++;
		}
	}
	//read in map descriptive file
	if(opts::_MAPDESC_.length() > 0){
		opts::printLog("Marker description being used: " + opts::_MAPDESC_ + "\n");
		ifstream dinput;
		dinput.open(opts::_MAPDESC_.c_str(), ios::in);
		if(!dinput){
			opts::printLog("Error opening marker description file: " + opts::_MAPDESC_ + "\n");
			exit(1);
		}
		string head = "";
		getline(dinput, head, '\n');
		descheaders = General::ParseDelimitedLine(head);
		while(!dinput.eof()){
			string id = "";
			string line = "";
			getline(dinput, line, '\n');
			if(line == ""){
				continue;
			}
			if(line[0] == '#'){
				continue;
			}

			vector<string> tokens;
			tokens = General::ParseDelimitedLine(line);
			if(tokens.size() < 2){
				opts::printLog("Marker description column size is < 2: " + line + "\n");
				exit(1);
			}
			if(tokens.size() != descheaders.size()){
				opts::printLog("Line is not the same size as header line: " + line + "\n");
				exit(1);
			}

			id = tokens[0];
			descinfo[id] = tokens;
		}
		dinput.close();
	}
	//read in marker exclusion file
	if(opts::_EXCMARKERS_.length() > 0){
		opts::printLog("Excluding markers from file: " + opts::_EXCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker exclusion file: " + opts::_EXCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";
		
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Excluding markers column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			exclude[tokens[0]] = 1;
			//exclude.push_back(probe);
		}
		einput.close();
	}
	//read in marker inclusion file
	if(opts::_INCMARKERS_.length() > 0){
		opts::printLog("Including markers from file: " + opts::_INCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker inclusion file: " + opts::_INCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";
		while(getline(einput,line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Including markers column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			include[tokens[0]] = 1;
			//include.push_back(probe);
		}
		einput.close();
	}
		
	
	ifstream input;
    input.open(opts::_MAPFILE_.c_str(), ios::in);

    if(!input){
		opts::printLog("Error opening map file: " + opts::_MAPFILE_ + "\n");
        exit(1);
    }
	int count = 0;
	int prev_chrom = -1;
	int prev_bploc = -1;
    while(!input.eof()){
        char buf[256];
        input.getline(buf, 256, '\n');

        string line = buf;

        if(line == ""){
            continue;
        }

        string temp;
        stringstream s(line);
        vector<string> elems;
        while(s >> temp){
            elems.push_back(temp);
        }

        if(elems.size() == 0){
            continue;
        }
        else if(elems.size() > 3){
			opts::printLog("Map file line has more than 3 elements: " + line + "\n");
            exit(1);
        }
		else if(elems.size() < 3 && elems.size() > 0){
			opts::printLog("Map file line has fewer than 3 elements: " + line + "\n");
			exit(1);
		}
			

        string chr = elems[0];
        string probe_id = elems[1];
        int bploc = atoi(elems[2].c_str());

		bool use = true;
		opts::_MARKERS_FOUND_++;
		if(exclude.size() > 0){
			//vector<string>::iterator found = find(exclude.begin(), exclude.end(), probe_id);
			map<string, int>::iterator found = exclude.find(probe_id);
			if(found != exclude.end()){
				use = false;
			}
		}
		if(include.size() > 0){
			//vector<string>::iterator found = find(include.begin(), include.end(), probe_id);
			map<string, int>::iterator found = include.find(probe_id);
			if(found == include.end()){
				use = false;
			}
		}

        Marker* m = new Marker(chr, probe_id, bploc);
		if(opts::_CHROM_LIMIT_){
			if(m->getChrom() != opts::_CHROM_){
				use = false;
			}
			else{
				if(opts::_BP_LOW_LIMIT_){
					if(bploc < opts::_BP_LOW_){
						use = false;
					}
				}
				if(opts::_BP_HIGH_LIMIT_){
					if(bploc > opts::_BP_HIGH_){
						use = false;
					}
				}
			}
		}
		if(opts::_BP_SPACE_LIMIT_){
			if(prev_bploc == -1){
				prev_bploc = bploc;
				prev_chrom = m->getChrom();
			}
			else{
				if(m->getChrom() == prev_chrom && ((bploc - prev_bploc) < opts::_BP_SPACE_)){
					use = false;
				}
				else{
					prev_bploc = bploc;
					prev_chrom = m->getChrom();
				}
			}
		}
		m->setEnabled(use);
		m->setLoc(count);
		m->setRSID(probe_id);
		
	
		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe_id);
			if(found != frequencies.end()){
				m->setMAF(frequencies[probe_id]);
				m->setFreqFlag(true);
			}
		}	
		if(opts::_MAPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[probe_id];
			for(int i = 1; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					m->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					m->assignDetail(descheaders[i], "NA");
				}
			}
		}
        markers->push_back(m);
		count++;
    }

    input.clear();
    input.close();

	marker_map->resize(markers->size());

	//put markers in chrom/bploc order
	stable_sort(markers->begin(), markers->end(), less<Marker*>());

	for(int i =0; i < markers->size(); i++){
		(*marker_map)[(*markers)[i]->getLoc()] = i;//(*markers)[i]->getLoc();
	}
	
}

/*
 *Function: readString
 *Description:
 *Reads the next string from a file, space/tab delimited
 *
 *
 */
/*bool readString(FILE* fp, string* s){
    bool done = false;
    *s="";
    while (1)
    {
        char ch = fgetc(fp);
        if ( ch==' ' || ch == '\t' )
        {
            if (done)
                return true;
        }
        else if ( ch=='\n' || ch=='\r' || feof(fp) )
            return false;
        else
        {
            *s += ch;
            done = true;
        }
    }
}*/


/*
 *Function: writeBit
 *Parameters: sample vector, family vector, marker vector, marker map vector
 *Description:
 *Writes genotype (bed) file as binary genotype file based on Plink output.
 *Family output file (fam) is first 6 columns of original ped file
 *Map file (bim) is extended map file including alleles.
 *
 *
 */
void writeBit(vector<Sample*>* samples, vector<Family*>* families, vector<Marker*>* markers, vector<int>* marker_map){
	opts::printLog("Writing family information to " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".fam\n");
	ofstream BIT((opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".fam").c_str(), ios::out);
	if(!BIT){
		opts::printLog("Error opening " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".fam for family file creation.  Exiting!\n");
		exit(1);
	}

	int ssize = samples->size();

	for(int s = 0; s < ssize; s++){
		Sample* samp = (*samples)[s];

		BIT << samp->getFamID() << " "
			<< samp->getInd() << " "
			<< samp->getDadID() << " "
			<< samp->getMomID() << " ";
		if(samp->getSex()){
			BIT << "1 ";
		}
		else{
			BIT << "2 ";
		}

		BIT << samp->getPheno();
		//if(samp->getAffected()){
		//	BIT << "2";
		//}
		//else{
		//	BIT << "0";
		//}
		BIT << endl;
	}
	BIT.clear();
	BIT.close();


	opts::printLog("Writing map information to " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bim\n");
	BIT.open((opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bim").c_str(), ios::out);
	if(!BIT){
		opts::printLog("Error opening " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bim for map file creation.  Exiting!\n");
		exit(1);
	}

	int msize = markers->size();

	for(int m = 0; m < msize; m++){
		Marker* mark = (*markers)[m];

		BIT << mark->getChrom() << "\t";
		BIT << mark->getProbeID() << "\t"
			<< "0\t" //centimorgan column
			<< mark->getBPLOC() << "\t";
		if(mark->getAllele1().length() == 0){
			BIT << "0\t";
		}
		else{
			BIT << mark->getAllele1() << "\t";
		}
		if(mark->getAllele2().length() == 0){
			BIT << "0\t";
		}
		else{
			BIT << mark->getAllele2() << "\t";
		}
//		if(mark->getRSID() == ""){
//			BIT << ".\t";
//		}
//		else{
//			BIT << mark->getRSID() << "\t";
//		}
//		if(mark->getEnzyme() == ""){
//			BIT << ".";
//		}
//		else{
//			BIT << mark->getEnzyme();
//		}
		BIT << endl;
	}
	BIT.clear();
	BIT.close();

	opts::printLog("Writing genotype bitfile to " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bed\n");
	BIT.open((opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bed").c_str(), ios::out | ios::binary);
	if(!BIT){
		opts::printLog("Error opening " + opts::_OUTPREFIX_ + opts::_BINPREFIX_ + ".bed for genotype file creation.  Exiting!\n");
		exit(1);
	}

	bitset<8> b;
	char ch[1];

  	b.reset();
    b.set(2);  b.set(3);  b.set(5);  b.set(6);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);

	b.reset();
	b.set(0);  b.set(1);  b.set(3);  b.set(4);
	ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);
				
	b.reset();
    ch[0] = (char)b.to_ulong();
	BIT.write(ch,1);

	for(int s = 0; s < ssize; s++){
		Sample* samp = (*samples)[s];
		
		for(int m = 0; m < msize;){
			bitset<8> b;
			b.reset();
			int c = 0;
			while(c < 8 && m < msize){
				Marker* mark = (*markers)[m];
				int loc = mark->getLoc();

				if(samp->getAone(loc)) 
					b.set(c);

				c++;
				if(samp->getAtwo(loc))
					b.set(c);
				c++;
				m++;
			}

			char ch[1];
			ch[0] = (char)b.to_ulong();
			BIT.write(ch, 1);
		}
	}
	BIT.close();
}

/*
 *Function: readBin
 *Paramters: sample vector, family vector, marker vector, marker_map vector
 *Description:
 *Reads set of binary input files (bim, bed, fam)
 *Based on Plink
 *
 */

void readBin(vector<Sample*>* samples, vector<Family*>* families, vector<Marker*>* markers, vector<int>* marker_map){
	//vector<string> exclude;
	map<string, vector<string> > mdescinfo;
	map<string,int> exclude;
	map<string,int> include;
	vector<string> mdescheaders;
	//vector<string> include;
	map<string,float> frequencies;
	if(opts::_FREQ_FILE_EXISTS_){
		opts::printLog("Marker Minor Allele Frequencies being read from: " + opts::_FREQ_FILE_ + "\n");
		ifstream finput;
		finput.open(opts::_FREQ_FILE_.c_str(), ios::in);
		if(!finput){
			opts::printLog("Error opening marker frequency file: " + opts::_FREQ_FILE_ + "\n");
			exit(1);
		}
		string line = "";
		int freqline = 1;
		while(getline(finput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2){
				opts::printLog("Marker frequency file column size != 2 on line: " + line + " Skipping line!!\n");
				continue;
			}
			try{
				float freq = -1.0f;
				for(int c = 0; c < tokens[1].size(); c++){
					if(!isdigit(tokens[1][c]) && tokens[1][c] != '.'){
						throw "Number format exception";
					}
				}
				freq = std::atof(tokens[1].c_str());
				frequencies[tokens[0]] = freq;	
			}catch(...){
				opts::printLog(tokens[1] + " is not a valid value for line: " + getString<int>(freqline) + "\n");
				exit(1);
			}
			freqline++;
		}
	}
	if(opts::_MAPDESC_.length() > 0){
        opts::printLog("Marker description being used: " + opts::_MAPDESC_ + "\n");
        ifstream dinput;
        dinput.open(opts::_MAPDESC_.c_str(), ios::in);
        if(!dinput){
            opts::printLog("Error opening marker description file: " + opts::_MAPDESC_ + "\n");
            exit(1);
        }
        string head = "";
        getline(dinput, head, '\n');
        mdescheaders = General::ParseDelimitedLine(head);
        while(!dinput.eof()){
            string id = "";
            string line = "";
            getline(dinput, line, '\n');
            if(line == ""){
                continue;
            }
            if(line[0] == '#'){
                continue;
            }
            vector<string> tokens;
            tokens = General::ParseDelimitedLine(line);
            if(tokens.size() < 2){
                opts::printLog("Marker description column size is < 2: " + line + "\n");
                exit(1);
            }
            if(tokens.size() != mdescheaders.size()){
                opts::printLog("Line is not the same size as header line: " + line + "\n");
                exit(1);
			}
            id = tokens[0];
	        mdescinfo[id] = tokens;
	    }
	    dinput.close();
	}
	if(opts::_EXCMARKERS_.length() > 0){
		opts::printLog("Excluding markers from file: " + opts::_EXCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker exclusion file: " + opts::_EXCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Marker exclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			exclude[tokens[0]] = 1;
			//exclude.push_back(probe);
		}
		einput.close();
	}
	if(opts::_INCMARKERS_.length() > 0){
		opts::printLog("Including markers from file: " + opts::_INCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker exclusion file: " + opts::_INCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Marker inclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			include[tokens[0]] = 1;
			//include.push_back(probe);
		}
		einput.close();
	}
	
	opts::printLog("Reading map from " + opts::_BINPREFIX_ + ".bim\n");

	ifstream MAP((opts::_BINPREFIX_ + ".bim").c_str(), ios::in);
	if(!MAP){
		opts::printLog("Error opening map file: " + opts::_BINPREFIX_ + ".bim.  Exiting!\n");
		exit(1);
	}
	//MAP.clear();

	int count = 0;

		string chrom = "";
		string probe = "";
		int bploc = 0;
		string a1 = "";
		string a2 = "";
		string centi = "";
		string rsid = "";
		string enzyme = "";
		string line = "";
	int prev_bploc = -1;
	int prev_chrom = -1;
	while(getline(MAP, line)){ 
			//>> chrom 
			//>> probe
			//>> bploc
			//>> a1
			//>> a2
			//>> rsid
			//>> enzyme){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 6){
			opts::printLog(".bim file column size != 6: " + line + " stopping!!\n");
			exit(0);
		}
		chrom = tokens[0];
		probe = tokens[1];
		centi = tokens[2];
		bploc = atoi(tokens[3].c_str());
		a1 = tokens[4];
		a2 = tokens[5];
		//rsid = tokens[5];
		//enzyme = tokens[6];
		//opts::_ENZYMES_ = true;
		if(rsid == "."){
			rsid = "";
		}
		if(enzyme == "."){
			enzyme = "";
			opts::_ENZYMES_ = false;
		}
		bool use = true;
		opts::_MARKERS_FOUND_++;

        if(exclude.size() > 0){
            //vector<string>::iterator found = find(exclude.begin(), exclude.end(), probe);
            map<string,int>::iterator found = exclude.find(probe);
            if(found != exclude.end()){
                use = false;
            }
			else{
				found = exclude.find(rsid);
				if(found != exclude.end()){
					use = false;
				}
			}
        }
        if(include.size() > 0){
            //vector<string>::iterator found = find(include.begin(), include.end(), probe);
            map<string,int>::iterator found = include.find(probe);
            if(found == include.end()){
                use = false;
				found = include.find(rsid);
				if(found == include.end()){
					use = false;
				}
				else{
					use = true;
				}
            }
        }

		Marker* m = new Marker(chrom, probe, bploc);
		if(opts::_CHROM_LIMIT_){
			if(m->getChrom() != opts::_CHROM_){
				use = false;
			}
			else{
				if(opts::_BP_LOW_LIMIT_){
					if(bploc < opts::_BP_LOW_){
						use = false;
					}
				}
				if(opts::_BP_HIGH_LIMIT_){
					if(bploc > opts::_BP_HIGH_){
						use = false;
					}
				}
			}
		}
        if(opts::_BP_SPACE_LIMIT_){
            if(prev_bploc == -1){
                prev_bploc = bploc;
                prev_chrom = m->getChrom();
            }
            else{
                if(m->getChrom() == prev_chrom && ((bploc - prev_bploc) < opts::_BP_SPACE_)){
                    use = false;
                }
                else{
                    prev_bploc = bploc;
                    prev_chrom = m->getChrom();
                }
            }
        }

		m->setEnabled(use);
		m->setLoc(count);
		m->setAllele1(a1);
		m->setAllele2(a2);
		m->setRSID(rsid);
		if(rsid == "." || rsid == ""){
			m->setRSID(probe);
		}
		m->setEnzyme(enzyme);
		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe);
			if(found == frequencies.end()){
				found = frequencies.find(m->getRSID());
				if(found != frequencies.end()){
					m->setMAF(frequencies[m->getRSID()]);
					m->setFreqFlag(true);
				}
			}
			else{
				m->setMAF(frequencies[m->getRSID()]);
				m->setFreqFlag(true);
			}
		}
        if(opts::_MAPDESC_.length() > 0 && mdescinfo.size() > 0){
			vector<string> tokens = mdescinfo[probe];
            for(int i = 1; i < mdescheaders.size(); i++){
	            if(tokens.size() == mdescheaders.size()){
		            m->assignDetail(mdescheaders[i], tokens[i]);
                }
                else{
	                m->assignDetail(mdescheaders[i], "NA");
                }
            }
        }

		markers->push_back(m);

		count++;
	}
	MAP.close();

	marker_map->resize(markers->size());
	stable_sort(markers->begin(), markers->end(), less<Marker*>());

	for(int i = 0; i < markers->size(); i++){
		(*marker_map)[(*markers)[i]->getLoc()] = i;
	}


	map<string, vector<string> > descinfo;
	vector<string> sexclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> finclude;
	vector<string> fexclude;
	vector<string> descheaders;
	if(opts::_SAMPDESC_.length() > 0){
		opts::printLog("Sample description being used: " + opts::_SAMPDESC_ + "\n");
		ifstream dinput;
		dinput.open(opts::_SAMPDESC_.c_str(), ios::in);
		if(!dinput){
			opts::printLog("Error opening sample description file: " + opts::_SAMPDESC_ + "\n");
			exit(1);
		}
        string head = "";
        getline(dinput, head, '\n');
        descheaders = General::ParseDelimitedLine(head);
		while(!dinput.eof()){
			string fam = "";
			string id = "";
			string center = "";
			string plate = "";
			string well = "";
//old        	char buf[256];
//old	        dinput.getline(buf, 256, '\n');
            string line = "";
			getline(dinput, line, '\n');

    	    //old string line = buf;

        	if(line == ""){
            	continue;
        	}
			vector<string> tokens = General::ParseDelimitedLine(line);
            if(tokens.size() < 2){
				opts::printLog("Sample description column size is < 2: " + line + " Skipping line!!\n");
	            continue;
	        }
	        if(tokens.size() != descheaders.size()){
				opts::printLog("Line is not the same size as header line: " + line + " Skipping line!!\n");
	            continue;
	        }
	
			fam = tokens[0];
			id = tokens[1];
			descinfo[fam + " " + id] = tokens;
//			center = tokens[2];
//			plate = tokens[3];
//			well = tokens[4];
//			descinfo[fam + " " + id] = center + " " + plate + " " + well;

//			descinfo.push_back(line);	
		}
		dinput.close();
	}
	if(opts::_EXCSAMPS_.length() > 0){
		opts::printLog("Excluding samples from file: " + opts::_EXCSAMPS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCSAMPS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening sample exclusion file: " + opts::_EXCSAMPS_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2 && line != ""){
				opts::printLog("Sample exclusion column size != 2: " + line + " Skipping line!!\n");
				continue;
			}
			sexclude.push_back(tokens[0] + " " + tokens[1]);
		}
		einput.close();
	}
	if(opts::_EXCFAMILIES_.length() > 0){
		opts::printLog("Excluding families from file: " + opts::_EXCFAMILIES_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCFAMILIES_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening family exclusion file: " + opts::_EXCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Family exclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			fexclude.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCCENTERS_.length() > 0){
		opts::printLog("Including centers from file: " + opts::_INCCENTERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCCENTERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening center inclusion file: " + opts::_INCCENTERS_ + "\n");
			exit(1);
		}
		string center = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Center inclusion columns != 1: " + line + " Skipping line!!\n");
				continue;
			}
			inccenters.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCSAMPLES_.length() > 0){
		opts::printLog("Including samples from file: " + opts::_INCSAMPLES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCSAMPLES_.c_str(), ios::in);
		if(!iinput){
			opts::printLog("Error opening sample inclusion file: " + opts::_INCSAMPLES_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		while(getline(iinput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2 && line != ""){
				opts::printLog("Sample inclusion columns != 2: " + line + " Skipping line!!\n");
				continue;
			}
			sinclude.push_back(tokens[0] + " " + tokens[1]);
		}
		iinput.close();
	}
	if(opts::_INCFAMILIES_.length() > 0){
		opts::printLog("Including families from file: " + opts::_INCFAMILIES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCFAMILIES_.c_str(), ios::in);
		if(!iinput){
			opts::printLog("Error opening family inclusion file: " + opts::_INCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(iinput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Family inclusion columns != 1: " + line + " Skipping line!!\n");
				continue;
			}
			finclude.push_back(tokens[0]);
		}
		iinput.close();
	}
	
	opts::printLog("Reading family information from " + opts::_BINPREFIX_ + ".fam\n");

	ifstream PED;
	PED.open((opts::_BINPREFIX_ + ".fam").c_str());
	if(!PED){
		opts::printLog("Error opening family information file: " + opts::_BINPREFIX_ + ".fam.  Exiting!\n");
		exit(1);
	}
	PED.clear();
		string fam = "";
		string ind = "";
		string dad = "";
		string mom = "";
		string sex = "";
		string aff = "";
		line = "";
	while(getline(PED, line)){
	   	//	>> fam
		//	>> ind
		//	>> dad
		//	>> mom
		//	>> sex
		//	>> aff){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 6){
			opts::printLog("Family information file column size != 6: " + line + " Exitting!!\n");
			exit(0);
		}
		fam = tokens[0];
		ind = tokens[1];
		dad = tokens[2];
		mom = tokens[3];
		sex = tokens[4];
		aff = tokens[5];
		Sample* samp = new Sample();
		samp->setFamID(fam);
		samp->setInd(ind);
		samp->setEnabled(true);
		opts::_SAMPLES_FOUND_++;
		if(sexclude.size() > 0){
			vector<string>::iterator found = find(sexclude.begin(), sexclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found != sexclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				if(opts::_KEEP_EXC_SAMPLES_){
					samp->setExcluded(true);
				}
			}
		}
		if(sinclude.size() > 0){
			vector<string>::iterator found = find(sinclude.begin(), sinclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found == sinclude.end()){
				samp->setEnabled(false);
				if(opts::_KEEP_EXC_SAMPLES_){
					samp->setExcluded(true);
				}
			}
		}

		samp->setDadID(dad);
		samp->setMomID(mom);
		if(sex == "1"){
			samp->setSex(true);
		}
		else if(sex == "2"){
			samp->setSex(false);
		}

		if(aff == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atoi(aff.c_str()));
		if(opts::pedinfo.size() > 0){
			map<string, Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
			if(sfind != opts::pedinfo.end()){
				Sample* sfound = sfind->second;
				samp->setDadID(sfound->getDadID());
				samp->setMomID(sfound->getMomID());
				samp->setSex(sfound->getSex());
				samp->setPheno(sfound->getPheno());
			}
		}
		
		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
				    samp->assignDetail(descheaders[i], "NA");
				}
			}

			//if(tokens.size() > 0){
			//	center = tokens[0];
			//}
			//if(opts::_INCCENTERS_.length() > 0 && inccenters.size() > 0){
			//	bool found = false;
			//	for(int cent = 0; cent < inccenters.size(); cent++){
			//		if(inccenters[cent] == center){
			//			found = true;
			//		}
			//	}
			//	if(!found){
			//		samp->setEnabled(false);
			//	}
//
			//}
			//samp->getFamily()->setCenter(tokens[2]);
			//if(tokens.size() > 0){
			//	samp->setPlate(tokens[1]);
			//	samp->setWell(tokens[2]);
			//}
		}
        vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Family* fam = new Family();
/*          fam->Setcenter()*/
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
			opts::_FAMILIES_FOUND_++;
        }
		if(fexclude.size() > 0){
			vector<string>::iterator found = find(fexclude.begin(), fexclude.end(), samp->getFamID());
			if(found != fexclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		if(finclude.size() > 0){
			vector<string>::iterator found = find(finclude.begin(), finclude.end(), samp->getFamID());
			if(found == finclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		samp->resizeAlleles(markers->size());
        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
	}	

	PED.clear();
	PED.close();

	bool ind_major = false;
	bool snp_major = false;
	opts::printLog("Reading genotype bitfile from " + opts::_BINPREFIX_ + ".bed\n");

	ifstream BIT;
	BIT.open((opts::_BINPREFIX_+".bed").c_str(), ios::in | ios::binary);
	if(!BIT){
	opts::printLog("Error opening genotype bitfile: " + opts::_BINPREFIX_ + ".bed.  Exiting!\n");
		exit(1);
	}
	char temp[1];
	//header
	BIT.read(temp, 1);
	bitset<8> tb;
	tb = temp[0];
	if((tb[2] && tb[3] && tb[5] && tb[6]) && !(tb[0] || tb[1] || tb[4] || tb[7])){
		BIT.read(temp, 1);
		tb = temp[0];
		if((tb[0] && tb[1] && tb[3] && tb[4]) && !(tb[2] || tb[5] || tb[6] || tb[7])){
			BIT.read(temp, 1);
			tb = temp[0];
			if(!tb[0]){
			opts::printLog("IND major mode\n");
				ind_major = true;
			}
			else{
			opts::printLog("SNP major mode\n");
				snp_major = true;
			}
		}
		else{
			opts::printLog("Incorrect bit file version (2nd code)!\n");
			exit(1);
		}
	}
	else{
	opts::printLog("Incorrect bit file version (1st code)!\n");
		exit(1);
	}
	
	int ssize = samples->size();
	int msize = markers->size();

	if(ind_major){
		for(int s = 0; s < ssize; s++){
			Sample* samp = (*samples)[s];
			for(int m = 0; m < msize;){
				char ch[1];
				//cout << samples->size() << "\t" << samp->getFamID() << "\t" << samp->getInd() << "\t" << s << "\t" << (*markers)[m]->getProbeID() << "\t" << m << endl;
				BIT.read(ch, 1);
				if(!BIT){
				opts::printLog("Problem with the bed file.\n");
					exit(1);
				}
				bitset<8> b;
				b = ch[0];
				int c = 0;
				while(c < 7 && m < msize){
					Marker* mark = (*markers)[m];
					if(mark->isEnabled()){
						int mloc = (*markers)[m]->getLoc();
						samp->addAone(mloc, b[c++]);
						samp->addAtwo(mloc, b[c++]);
					}
					else{
//					samp->addAone(m, true);
//					samp->addAtwo(m, false);
						c+=2;
					}
					m++;
				}
			}
		}	
	}
	else{//SNP_MAJOR
		for(int m = 0; m < msize; m++){
			Marker* mark = (*markers)[m];
			int mloc = mark->getLoc();
			for(int s = 0; s < ssize;){
				char ch[1];
				//cout << samples->size() << "\t" << samp->getFamID() << "\t" << samp->getInd() << "\t" << s << "\t" << (*markers)[m]->getProbeID() << "\t" << m << endl;
				BIT.read(ch, 1);
				if(!BIT){
					opts::printLog("Problem with the bed file.\n");
					exit(1);
				}
				bitset<8> b;
				b = ch[0];
				int c = 0;
				while(c < 7 && s < ssize){
					Sample* samp = (*samples)[s];
					if(samp->isEnabled() || (samp->isExcluded() && opts::_KEEP_EXC_SAMPLES_)){
						samp->addAone(mloc, b[c++]);
						samp->addAtwo(mloc, b[c++]);
					}
					else{
//					samp->addAone(m, true);
//					samp->addAtwo(m, false);
						c+=2;
					}
					s++;
				}
			}
		}	
	}
	BIT.clear();
	BIT.close();
}

/*
 *Function: readPed
 *Description:
 *Reads ped file and stores data into sample, family, and marker vectors
 *
 */
void readPed(vector<Sample*>* samples, vector<Family*>* families, vector<Marker*>* markers, vector<int>* marker_map){
	map<string, vector<string> > descinfo;
	vector<string> exclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> fexclude;
	vector<string> finclude;
	vector<string> descheaders;
	if(opts::_SAMPDESC_.length() > 0){
		opts::printLog("Sample description being used: " + opts::_SAMPDESC_ + "\n");
		ifstream dinput;
		dinput.open(opts::_SAMPDESC_.c_str(), ios::in);
		if(!dinput){
		opts::printLog("Error opening sample description file: " + opts::_SAMPDESC_ + ".  Exiting!\n");
			exit(1);
		}
		string head = "";
		getline(dinput, head, '\n');
		descheaders = General::ParseDelimitedLine(head);
		while(!dinput.eof()){
			string fam = "";
			string id = "";
			string center = "";
			string plate = "";
			string well = "";
     //old   	char buf[256];
	 //old       dinput.getline(buf, 256, '\n');
	 		string line = "";
			getline(dinput, line, '\n');

//old    	    string line = buf;

        	if(line == ""){
            	continue;
        	}
//			vector<string> temp_tokens;
//			Tokenize(line, temp_tokens, " ");
			vector<string> tokens;
			tokens = General::ParseDelimitedLine(line);
			if(tokens.size() < 2){
				opts::printLog("Sample description column size is < 2: " + line + " Skipping line!!\n");
				continue;
			}
			if(tokens.size() != descheaders.size()){
			opts::printLog("Line is not the same size as header line: " + line + " Skipping line!!\n");
				continue;
			}

			fam = tokens[0];
			id = tokens[1];
			descinfo[fam + " " + id] = tokens;
//			for(int i = 0; i < temp_tokens.size(); i++){
//				Tokenize(temp_tokens[i], tokens, "\t");
//			}

/*old			if(tokens.size() != 5){
				cerr << "Sample description column size != 5: " << line << " Skipping line!!" << endl;
				continue;
			}
			fam = tokens[0];
			id = tokens[1];
			center = tokens[2];
			plate = tokens[3];
			well = tokens[4];
			descinfo[fam + " " + id] = center + " " + plate + " " + well;
*/

//			descinfo.push_back(line);	
		}
		dinput.close();
	}
	if(opts::_EXCSAMPS_.length() > 0){
		opts::printLog("Excluding samples from file: " + opts::_EXCSAMPS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCSAMPS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening sample exclusion file: " + opts::_EXCSAMPS_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		//while(einput >> fam >> ind){
		while(getline(einput,line)){
//			vector<string> temp_tokens;
			vector<string> tokens;
			tokens = General::ParseDelimitedLine(line);
//			Tokenize(line, temp_tokens, " ");
//			for(int i = 0; i < temp_tokens.size(); i++){
//				Tokenize(temp_tokens[i], tokens, "\t");
//			}
			if(tokens.size() != 2 && line != ""){
			opts::printLog("Sample Exclusion column size != 2: " + line + " Skipping line!!\n");
				continue;
			}
			exclude.push_back(tokens[0] + " " + tokens[1]);
		}
		einput.close();
	}
	if(opts::_EXCFAMILIES_.length() > 0){
	opts::printLog("Excluding families from file: " + opts::_EXCFAMILIES_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCFAMILIES_.c_str(), ios::in);
		if(!einput){
		opts::printLog("Error opening family exclusion file: " + opts::_EXCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(einput, line)){
//			vector<string> temp_tokens;
			vector<string> tokens;
			tokens = General::ParseDelimitedLine(line);
//			Tokenize(line, temp_tokens, " ");
//			for(int i = 0; i < temp_tokens.size(); i++){
//				Tokenize(temp_tokens[i], tokens, "\t");
//			}
			if(tokens.size() != 1 && line != ""){
			opts::printLog("Family Exclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			fexclude.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCCENTERS_.length() > 0){
	opts::printLog("Including centers from file: " + opts::_INCCENTERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCCENTERS_.c_str(), ios::in);
		if(!einput){
		opts::printLog("Error opening center inclusion file: " + opts::_INCCENTERS_ + "\n");
			exit(1);
		}
		string center = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
			opts::printLog("Center Inclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			inccenters.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCSAMPLES_.length() > 0){
	opts::printLog("Including samples from file: " + opts::_INCSAMPLES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCSAMPLES_.c_str(), ios::in);
		if(!iinput){
		opts::printLog("Error opening sample inclusion file: " + opts::_INCSAMPLES_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		while(getline(iinput,line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2 && line != ""){
			opts::printLog("Sample inclusion column size != 2:" + line + " Skipping line!!\n");
				continue;
			}
			sinclude.push_back(tokens[0] + " " + tokens[1]);
		}
		iinput.close();
	}
	if(opts::_INCFAMILIES_.length() > 0){
	opts::printLog("Including families from file: " + opts::_INCFAMILIES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCFAMILIES_.c_str(), ios::in);
		if(!iinput){
		opts::printLog("Error opening family inclusion file: " + opts::_INCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(iinput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
			opts::printLog("Family inclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			finclude.push_back(tokens[0]);
		}
		iinput.close();
	}
	
    FILE* input;
    input = fopen(opts::_PEDFILE_.c_str(), "r");
	if(!input){
	opts::printLog("Error opening pedfile: " + opts::_PEDFILE_ + ".  Exiting!\n");
		exit(1);
	}
	int onind = -1;
    while(!feof(input)){
		onind++;
        Sample* samp = new Sample();
        int f = 0;
        string temp = "";
        if(readString(input, &temp)){
           samp->setFamID(temp);
              f++;
             temp = "";
         }

		string ftemp = samp->getFamID();
		if(samp->getFamID() == ""){
			delete(samp);
			continue;
		}
        if(ftemp.at(0) == '#'){
			delete(samp);
			while(fgetc(input) != '\n' && !feof(input)){}
            continue;
        }
		opts::_SAMPLES_FOUND_++;
        /*check for comments?*/
        string sex = "";
        string pheno = "";
        if(readString(input, &temp)){
            samp->setInd(temp);
            f++;
            temp = "";
        }
		samp->setEnabled(true);

		if(exclude.size() > 0){
			vector<string>::iterator found = find(exclude.begin(), exclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found != exclude.end()){
				if(!opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					continue;
				}
				else{
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
		}
		if(sinclude.size() > 0){
			vector<string>::iterator found = find(sinclude.begin(), sinclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found == sinclude.end()){
				if(!opts::_KEEP_EXC_SAMPLES_){
					delete(samp);
					while(fgetc(input) != '\n' && !feof(input)){}
					continue;
				}
				else{
					samp->setEnabled(false);
					samp->setExcluded(true);
				}
			}
		}
        if(readString(input, &temp)){
            samp->setDadID(temp);
            f++;
            temp = "";
        }
        if(readString(input, &temp)){
            samp->setMomID(temp);
            f++;
            temp = "";
        }
        if(readString(input, &sex)) f++;
        if(readString(input, &pheno)) f++;

        if(sex == "1"){
            samp->setSex(true);
        }
        else if(sex == "2"){
            samp->setSex(false);
        }

		if(pheno == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atoi(pheno.c_str()));

		if(opts::pedinfo.size() > 0){
			map<string, Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
			if(sfind != opts::pedinfo.end()){
				Sample* sfound = sfind->second;
				samp->setDadID(sfound->getDadID());
				samp->setMomID(sfound->getMomID());
				samp->setSex(sfound->getSex());
				samp->setPheno(sfound->getPheno());
			}
		}

		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					samp->assignDetail(descheaders[i], "NA");
				}
			}
			//if(tokens.size() > 0){
			//	center = tokens[0];
			//}
			//if(opts::_INCCENTERS_.length() > 0 && inccenters.size() > 0){
			//	bool found = false;
			//	for(int cent = 0; cent < inccenters.size(); cent++){
			//		if(inccenters[cent] == center){
			//			found = true;
			//		}
			//	}
			//	if(!found){
			//		delete(samp);
			//		while(fgetc(input) != '\n' && !feof(input)){}
			//		continue;
			//	}

			//}
			//samp->getFamily()->setCenter(tokens[2]);
			//if(tokens.size() > 0){
			//	samp->setPlate(tokens[1]);
			//	samp->setWell(tokens[2]);
			//}
		}
        samp->resizeAlleles(markers->size());

        int gn = 0;
        int i = 0;
        bool linedone = false;
        bool fatal = false;

        string fmsg;
        while(!linedone){
            string one = "";
            string two = "";

            while(1){
                char ch = fgetc(input);

                if(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input)){
                    if(ch == '\n' || ch == '\r' || feof(input)){
                        linedone = true;
                    }

                    if(one.length() > 0){
                        gn++;
                        break;
                    }
                    if(ch == '\n' || ch == '\r' || feof(input)){
                        break;
                    }
                }
                else{
                    one += ch;
                }
            }
            if(!linedone){
                while(1){
                    char ch = fgetc(input);
                    if(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(input)){
                        if(ch == '\n' || ch == '\r' || feof(input)){
                            linedone = true;
                        }
                        if(two.length() > 0){
                            gn++;
                            break;
                        }
                        if(ch == '\n' || ch == '\r' || feof(input)){
                            break;
                        }
                    }
                    else{
                        two += ch;
                    }
                }
                if(linedone && one.length() == 0 && two.length() == 0){
                    break;
                }
                
				//Marker* m = (*markers)[i];
				if(i > markers->size()){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>((2 * markers->size()) + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					opts::printLog(text);
					exit(0);
				}
				Marker* m = (*markers)[(*marker_map)[i]];
				if(m->isEnabled()){
					int oldallelecount = m->getNumAlleles();
	                if(one != opts::_NOCALL_){
						//new
						if(m->getAlleleLoc(one) < 0){
							m->addAllele(one);
						}
/*old
    	                if(one != m->getAllele1() && one != m->getAllele2()){
        	                if(m->getAllele1() == ""){
            	                m->setAllele1(one);
                	        }
                    	    else if(m->getAllele2() == ""){
                        	    m->setAllele2(one);
	                        }
    	                }
old*/
        	        }

            	    if(two != one){
                	    if(two != opts::_NOCALL_){
							//new
							if(m->getAlleleLoc(two) < 0){
								m->addAllele(two);
							}
/*old                    	    if(two != m->getAllele1() && two != m->getAllele2()){
                        	    if(m->getAllele1() == ""){
                            	    m->setAllele1(two);
	                            }
    	                        else if(m->getAllele2() == ""){
        	                        m->setAllele2(two);
            	                }
                	        }
old*/
                    	}
	                }


					if(m->getNumAlleles() <= 2){
   		 	            if(one == m->getAllele1() && two == m->getAllele1()){
   	    	     	        samp->addAone(i, false);
   	     		            samp->addAtwo(i, false);
		                }
	    	            else if(one != opts::_NOCALL_ && two != opts::_NOCALL_ && one != two){
	        	            samp->addAone(i, false);
 		           	        samp->addAtwo(i, true);
	                	}
		                else if(one == m->getAllele2() && two == m->getAllele2()){
	    	                samp->addAone(i, true);
	        	            samp->addAtwo(i, true);
	            	    }
	                	else if(one == opts::_NOCALL_ || two == opts::_NOCALL_){
		                    samp->addAone(i, true);
	    	                samp->addAtwo(i, false);
	        	        }
					}
					else if(opts::_MICROSATS_){
						samp->addMicroSat(i);
						int loc1 = m->getAlleleLoc(one);
						int loc2 = m->getAlleleLoc(two);
						//samp->addAbone(loc1);
						//samp->addAbtwo(loc2);
						samp->addAbone(i, loc1);
						samp->addAbtwo(i, loc2);
						if(oldallelecount <= 2){
							remapSamples(samples, markers, marker_map, i);
						}
					}
					else if(m->getNumAlleles() > 2 && !opts::_MICROSATS_){
						string alleles = "";
						for(int aa = 0; aa < m->getAlleles().size(); aa++){
							alleles += " " + m->getAllele(aa);
						}
						opts::printLog("More than 2 unique alleles found for SNP: " + m->getProbeID() + ", PED File line: " + getString<int>(onind + 1) + ".  Microsatellites not specified. (Offending alleles: " + alleles + "\n");
						exit(1);
					}
				}
				else{
					samp->addAone(i,true);
					samp->addAtwo(i, false);
				}
    	    	i++;
				if(i > markers->size()){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>((2 * markers->size()) + 6);
					text += " columns but found ";
					text += getString<int>(f + gn);
					text += "\n";
					opts::printLog(text);
					
					exit(0);
				}
            }/*end !linedone*/
        }/*end while(1)*/
		if(gn != (2* markers->size())){
					string text = "Problem with line: ";
					text += getString<int>(onind + 1);
					text += " in file: " + opts::_PEDFILE_ + "\n";
					text += "Expecting ";
					text += getString<int>(((2 * markers->size()) + 6));
					text += " columns but found ";
					text += getString<int>((f + gn));
					text += "\n";
					opts::printLog(text);
			exit(0);
		}

        vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Family* fam = new Family();
/*          fam->Setcenter()*/
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
			opts::_FAMILIES_FOUND_++;
        }
		if(fexclude.size() > 0){
			vector<string>::iterator found = find(fexclude.begin(), fexclude.end(), samp->getFamID());
			if(found != fexclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				families->erase(families->end() - 1);
				delete(samp->getFamily());
				delete(samp);
				continue;
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		if(finclude.size() > 0){
			vector<string>::iterator found = find(finclude.begin(), finclude.end(), samp->getFamID());
			if(found == finclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
    }/*end while(eof)*/

    fclose(input);
}

/*
 *Function: remapSamples
 *Description:
 *If -micro-sats is enabled, when a 3rd allele is found for a marker, this function
 *moves the Sample genotype storage for the specific marker from boolean storage to integer storage.
 *
 */
/*void remapSamples(vector<Sample*>* samples, vector<Marker*>* markers, vector<int>* marker_map, int loc){

	Marker* m = (*markers)[(*marker_map)[loc]];
	for(int i = 0; i < samples->size(); i++){
		Sample* samp = (*samples)[i];
		if(!samp->haveMicroSat(loc)){
			samp->addMicroSat(loc);
			if(!samp->getAone(loc) && !samp->getAtwo(loc)){
				samp->addAbone(loc, 0);
				samp->addAbtwo(loc, 0);
			}
			else if(samp->getAone(loc) && samp->getAtwo(loc)){
				samp->addAbone(loc, 1);
				samp->addAbtwo(loc, 1);
			}
			else if(!samp->getAone(loc) && samp->getAtwo(loc)){
				samp->addAbone(loc, 0);
				samp->addAbtwo(loc, 1);
			}
			else{
				cout << "Remapping sample: " << samp->toString() << "\t-1/-1\n";
				samp->addAbone(loc, -1);
				samp->addAbtwo(loc, -1);
			}
			samp->addAone(loc, true);
			samp->addAtwo(loc, false);
		}
	}
}*/

/*
 *Function: flipStrand
 *Description:
 *Flips the strand of the specified markers
 *
 */
void flipStrand(vector<Marker*>* markers){
	if(opts::_FLIPFILE_.length() > 0){
		map<string, int> exclude;
		opts::printLog("Flipping markers found in: " + opts::_FLIPFILE_ + "\n");
		ifstream einput;
		einput.open(opts::_FLIPFILE_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker strand flip file: " + opts::_FLIPFILE_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";
		
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Strand flip markers column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			exclude[tokens[0]] = 1;
			//exclude.push_back(probe);
		}
		einput.close();

		int msize = markers->size();

		for(int m = 0; m < msize; m++){
			Marker* mark = (*markers)[m];
			if(!mark->isMicroSat()){
				map<string, int>::iterator found = exclude.find(mark->getProbeID());
    	        if(found != exclude.end()){
					string a1 = mark->getAllele1();
					string a2 = mark->getAllele2();

					if(a1 == "A"){
						mark->resetAllele1("T");
					}
					else if(a1 == "C"){
						mark->resetAllele1("G");
					}
					else if(a1 == "G"){
						mark->resetAllele1("C");
					}
					else if(a1 == "T"){
						mark->resetAllele1("A");
					}
					else if(a1 == "1"){
						mark->resetAllele1("4");
					}
					else if(a1 == "2"){
						mark->resetAllele1("3");
					}
					else if(a1 == "3"){
						mark->resetAllele1("2");
					}
					else if(a1 == "4"){
						mark->resetAllele1("1");
					}

					if(a2 == "A"){
						mark->resetAllele2("T");
					}
					else if(a2 == "C"){
						mark->resetAllele2("G");
					}
					else if(a2 == "G"){
						mark->resetAllele2("C");
					}
					else if(a2 == "T"){
						mark->resetAllele2("A");
					}
					else if(a2 == "1"){
						mark->resetAllele2("4");
					}
					else if(a2 == "2"){
						mark->resetAllele2("3");
					}
					else if(a2 == "3"){
						mark->resetAllele2("2");
					}
					else if(a2 == "4"){
						mark->resetAllele2("1");
					}
				}
			}
		}
	}
	
}


/*
 *Function: printFamilies
 *Description:
 *outputs the structure of all families
 *
 */
void printFamilies(vector<Family*>* families){
	int fsize = families->size();

	string fname = opts::_OUTPREFIX_ + "family_structure.txt";
	ofstream fout (fname.c_str());
	if(!fout.is_open()){
		opts::printLog("Unable to open " + fname + "\n");
		exit(1);
	}
	fname = opts::_OUTPREFIX_ + "family_data.txt";
	ofstream fdout (fname.c_str());
	if(!fdout.is_open()){
		opts::printLog("Unable to open " + fname + "\n");
		exit(1);
	}
	//output structure
	for(int f = 0; f < fsize; f++){
		Family* fam = (*families)[f];
		if(fam->isEnabled()){
			fout << "Family: " << fam->getFamID() << endl;
			fout << "----------------------------------------------------------" << endl;
			vector<Sample*>* founders = fam->getFounders();
			int size = founders->size();
			if(size > 0){
				for(int fs = 0; fs < size; fs++){
					Sample* founder = (*founders)[fs];
					if(founder->isEnabled() || (founder->isExcluded() && opts::_KEEP_EXC_SAMPLES_)){
						fout << "Founder: " << founder->getInd() << endl;
						fout << descendTree(founder, 0) << endl;
					}
				}
			}
			else{
				fout << "No Founders...\n";
				vector<Sample*>* samps = fam->getSamples();
				int ssize = samps->size();
				for(int ss = 0; ss < ssize; ss++){
					Sample* samp = (*samps)[ss];
					if(samp->isEnabled() || (samp->isExcluded() && opts::_KEEP_EXC_SAMPLES_)){
						fout << descendTree(samp, 0) << endl;
					}
				}
			}
			fout << endl << endl;
		}
	}

	fdout << "FamID\tMultigenerational?\tMutigen_with_affecteds\tTotal_Generations\tTotal_Num_Affected\tGeneration1-N_Affected Counts\n";
	//output summary levels & affecteds
	for(int f = 0; f < fsize; f++){
		Family* fam = (*families)[f];
		if(fam->isEnabled()){
			map<int, vector<Sample*> > levels;
			fdout << fam->getFamID();
			vector<Sample*>* founders = fam->getFounders();
			int size = founders->size();
			if(size > 0){
				//vector<int> levels(size);
				for(int fs = 0; fs < size; fs++){
					Sample* founder = (*founders)[fs];
					if(founder->isEnabled() || (founder->isExcluded() && opts::_KEEP_EXC_SAMPLES_)){
						//levels[fs] = descendTree3(founder, 0);
						map<int, vector<Sample*> > levelstemp = descendTree3(founder, 1);
						//if(founder->getPheno() == 2){// && !founder->isFlagged()){
						map<int, vector<Sample*> >::iterator titer;
						for(titer = levelstemp.begin(); titer != levelstemp.end(); titer++){
							vector<Sample*> mysamps = titer->second;
							for(int ms = 0; ms < mysamps.size(); ms++){
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
			}
			else{
				vector<Sample*>* samps = fam->getSamples();
				int ssize = samps->size();
				//vector<int> levels(ssize);
				//map<int, int> levels;
				for(int ss = 0; ss < ssize; ss++){
					Sample* samp = (*samps)[ss];
					if(samp->isEnabled() || (samp->isExcluded() && opts::_KEEP_EXC_SAMPLES_)){
						//levels[ss] = descendTree3(samp, 0);
						map<int, vector<Sample*> > levelstemp = descendTree3(samp, 1);
						//if(samp->getPheno() == 2){// && !samp->isFlagged()){
						map<int, vector<Sample*> >::iterator titer;
						for(titer = levelstemp.begin(); titer != levelstemp.end(); titer++){
							vector<Sample*> mysamps = titer->second;
							for(int ms = 0; ms < mysamps.size(); ms++){
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
			vector<Sample*>* fsamps = fam->getSamples();
			int affected = 0;
//			for(int s = 0; s < fsamps->size(); s++){
//				if((*fsamps)[s]->getPheno() == 2 && (*fsamps)[s]->isFlagged()){
					//affected++;
//					(*fsamps)[s]->setFlag(false);
//				}
//			}
			map<int, vector<Sample*> >::iterator iter = levels.end();
			if(levels.size() > 1){
				fdout << "\tY";
			}
			else{
				fdout << "\tN";
			}
			for(iter = levels.begin(); iter != levels.end(); iter++){
				vector<Sample*> mysamps = iter->second;
				map<Sample*, bool> newsamps;
				for(int ms = 0; ms < mysamps.size(); ms++){
					newsamps[mysamps[ms]] = mysamps[ms]->getAffected();
				}
				map<Sample*, bool>::iterator niter;
				for(niter = newsamps.begin(); niter != newsamps.end(); niter++){

					if(niter->second){
						bool useme = true;
						map<int, vector<Sample*> >::iterator tempiter;
						for(tempiter = levels.begin(); tempiter != levels.end(); tempiter++){
							if(tempiter->first > iter->first){
								vector<Sample*> findsamples = tempiter->second;
								vector<Sample*>::iterator findme = find(findsamples.begin(), findsamples.end(), niter->first);
								if(findme != findsamples.end()){
									useme = false;
									break;
								}
							}
						}
						if(useme)
							affected++;
					}
				}
			}
			string affcounts = "";
			int numlevelsofaff = 0;
			for(iter = levels.begin(); iter != levels.end(); iter++){
				int laffected = 0;
				vector<Sample*> mysamps = iter->second;
				map<Sample*, bool> newsamps;
				for(int ms = 0; ms < mysamps.size(); ms++){
					newsamps[mysamps[ms]] = mysamps[ms]->getAffected();
				}
				bool incme = false;
				map<Sample*, bool>::iterator niter;
				for(niter = newsamps.begin(); niter != newsamps.end(); niter++){
					if(niter->second){
						bool useme = true;
						map<int, vector<Sample*> >::iterator tempiter;
						for(tempiter = levels.begin(); tempiter != levels.end(); tempiter++){
							if(tempiter->first > iter->first){
								vector<Sample*> findsamples = tempiter->second;
								vector<Sample*>::iterator findme = find(findsamples.begin(), findsamples.end(), niter->first);
								if(findme != findsamples.end()){
									useme = false;
									break;
								}
							}
						}
						if(useme){
							laffected++;
							if(!incme){
								numlevelsofaff++;
								incme = true;
							}
						}
					}
				}
				affcounts += "\t" + getString<int>(laffected);
			}
			if(numlevelsofaff > 1){
				fdout << "\tY";
			}
			else{
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

/*
 *Function: descendTree
 *Description:
 *Recursively moves through pedigree until the child leaves are found to find structure
 */
string descendTree(Sample* sample, int level){
	string infospace = "     ";
	for(int i = 0; i < (level * 2); i++){
		infospace += "     ";
	}
	vector<Sample*>* children = sample->getChildren();
	if(children->size() == 0){
		string v = infospace + "->Ind: " + sample->getInd() + "\n" + infospace + "->Mom: " + sample->getMomID() + "\n" + infospace + "->Dad: " + sample->getDadID() + "\n";
		if(sample->getSib() != NULL){
			v += infospace + "->Sibling: " + sample->getSib()->getInd() + "\n";
		}
	   	v += infospace + "->Children: 0\n";
		return v;
	}
	int csize = children->size();

	string val = "";
	string v = infospace + "->Ind: " + sample->getInd() + "\n" + infospace + "->Mom: " + sample->getMomID() + "\n" + infospace + "->Dad: " + sample->getDadID() + "\n";
	if(sample->getSib() != NULL){
		v += infospace + "->Sibling: " + sample->getSib()->getInd() + "\n";
	}
	v += infospace + "->Children: " + getString<int>(sample->getChildren()->size()) + "\n";
	val += v;
	for(int c = 0; c < csize; c++){
		string temp = descendTree((*children)[c], (level + 1));
		 
		val += temp;
	}
	return val;
}

/*
 *Function: descendTree2
 *Description: returns # of levels
 *Recursively moves through pedigree until the child leaves are found to find structure
 */
int descendTree2(Sample* sample, int level){
	vector<Sample*>* children = sample->getChildren();
	if(children->size() == 0){
		return level;
	}
	int csize = children->size();
	vector<int> levels(csize);
	for(int c = 0; c < csize; c++){
		levels[c] = descendTree2((*children)[c], (level + 1));
	}
	sort(levels.begin(), levels.end());
	level = levels[levels.size() - 1];
	return level;
}

/*
 *Function: descendTree3
 *Description: returns # of levels
 *Recursively moves through pedigree until the child leaves are found to find structure
 */
map<int, vector<Sample*> > descendTree3(Sample* sample, int level){
	map<int, vector<Sample*> > values;

	vector<Sample*>* children = sample->getChildren();
	if(children->size() == 0){
		return values;
		//return level;
	}
	int csize = children->size();
	vector<int> levels(csize);
	for(int c = 0; c < csize; c++){
		map<int, vector<Sample*> > tempvalues = descendTree3((*children)[c], (level + 1));
		map<int, vector<Sample*> >::iterator iter;
		for(iter = tempvalues.begin(); iter != tempvalues.end(); iter++){
			vector<Sample*> samps = iter->second;
			for(int is = 0; is < samps.size(); is++){
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

/*
 *Function: compileOutputs
 *Description:
 *Compiles all QC & analysis outputs into single marker, sample, and family files
 *Performed when all steps are complete.
 */
void compileOutputs(vector<Marker*>* markers, vector<Family*>* families, vector<Sample*>* samples){
	opts::printLog("Compiling QC and Analysis output files...\n");


	map<string, vector<string> > filenames = opts::getFilenames();
	if(filenames["Marker"].size() > 0){
		map<Marker*, vector<string> > marker_output;
		vector<string> all_columns;
		map<string, Marker*> snpmap;
		for(int i = 0; i < markers->size(); i++){
			Marker* mark = (*markers)[i];
			string chr = getString<int>(mark->getChrom());
			string bp = getString<int>(mark->getBPLOC());
			snpmap[chr + "#" + bp] = mark;
		}

		for(int i = 0; i < filenames["Marker"].size(); i++){
			string file = filenames["Marker"][i];
			string step = opts::filesteps[file];

			for(int j = 0; j < opts::fileheaders[file].size(); j++){
				all_columns.push_back("(" + file + ")" + opts::fileheaders[file][j]);
			}
		}
		string outfile = opts::_OUTPREFIX_ + "marker_summary.txt";
		opts::printLog("Working on compiling SNP information...[" + outfile + "]\n");
		for(int i = 0; i < filenames["Marker"].size(); i++){
			string filename = filenames["Marker"][i];
			string step = opts::filesteps[filename];
			vector<string> filecols = opts::fileheaders[filename];

			opts::printLog("\tParsing " + filename + "\n");

	
			ifstream input;
			input.open(filename.c_str(), ios::in);
	
			if(!input){
				cerr << "Error opening file: " << filename << endl;
				exit(1);
			}
	
			string header = "";
		   	getline(input, header);

			vector<string> columns;
			General::Tokenize(header, columns, "\t");
			int chrloc = -1;
			int bploc = -1;
			for(int l = 0; l < columns.size(); l++){
				if(columns[l] == "Chrom"){
					chrloc = l;
				}
				if(columns[l] == "bploc"){
					bploc = l;
				}
				if(chrloc > -1 && bploc > -1){
					break;
				}
			}
			string line;
			while(getline(input, line)){
				int count = 0;
				vector<string> elems;
				string token;
				istringstream isstream (line);
				while(getline(isstream, token, '\t')){
					elems.push_back(token);
				}
				//General::Tokenize(line, elems, "\t");
				if(elems.size() == 0){
					continue;
				}
				if(elems.size() != columns.size()){
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
				if(itermark == NULL){
					cout << "Cannot find marker with chrom = " << chrom << " and bploc = " << bp << endl;
					exit(1);
				}
				Marker* mark = itermark;
				map<Marker*, vector<string> >::iterator data = marker_output.find(mark);
				if(data == marker_output.end()){
					vector<string> mydata;
					mydata.resize(all_columns.size(), "N/A");
					marker_output[mark] = mydata;
				}
				for(int j = 0; j < filecols.size(); j++){
					vector<string>::iterator realcol = find(columns.begin(), columns.end(), filecols[j]);
					vector<string>::iterator allloc = find(all_columns.begin(), all_columns.end(), "(" +filename+")"+filecols[j]);
					if(realcol != columns.end()){
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
		for(int i = 0; i < all_columns.size(); i++){
			out << "\t" << all_columns[i];
		}
		out << "\n";
		map<Marker*, vector<string> >::iterator iter;
		for(iter = marker_output.begin(); iter != marker_output.end(); iter++){
			Marker* mark = iter->first;
			vector<string> data = iter->second;
			out << mark->getChrom() << "\t"
				<< mark->getRSID() << "\t"
				<< mark->getProbeID() << "\t"
				<< mark->getBPLOC();
			for(int i = 0; i < data.size(); i++){
				out << "\t" << data[i];
			}
			out << "\n";
		}
		out.close();
	}
	if(filenames["Family"].size() > 0){
	map<Family*, vector<string> > family_output;
		vector<string> all_columns;
		for(int i = 0; i < filenames["Family"].size(); i++){
			string file = filenames["Family"][i];
			string step = opts::filesteps[file];

			for(int j = 0; j < opts::fileheaders[file].size(); j++){
				all_columns.push_back("(" + file + ")" + opts::fileheaders[file][j]);
			}
		}
		string outfile = opts::_OUTPREFIX_ + "family_summary.txt";
		opts::printLog("Working on compiling Family information...[" + outfile + "]\n");
		for(int i = 0; i < filenames["Family"].size(); i++){
			string filename = filenames["Family"][i];
			string step = opts::filesteps[filename];
			vector<string> filecols = opts::fileheaders[filename];

			opts::printLog("\tParsing " + filename + "\n");

	
			ifstream input;
			input.open(filename.c_str(), ios::in);
	
			if(!input){
				cerr << "Error opening file: " << filename << endl;
				exit(1);
			}
	
			string header = "";
		   	getline(input, header);

			vector<string> columns;
			General::Tokenize(header, columns, "\t");
			int famloc = -1;
			for(int l = 0; l < columns.size(); l++){
				if(columns[l] == "FamID"){
					famloc = l;
				}
				if(famloc > -1){
					break;
				}
			}
			if(famloc < 0){
				cout << "FamID column not found!\n";
				exit(1);
			}
			string line;
			int count = 1;
			while(getline(input, line)){
				count++;
				vector<string> elems;
				string token;
				istringstream isstream (line);
				while(getline(isstream, token, '\t')){
					elems.push_back(token);
				}
				//General::Tokenize(line, elems, "\t");
				if(elems.size() == 0){
					continue;
				}
				string famid = elems[famloc];
				vector<Family*>::iterator iterfam = find_if(families->begin(), families->end(), FindFamily(famid));
				if(iterfam == families->end()){
					cout << "Cannot find family with famid = " << famid << endl;
					exit(1);
				}
				Family* fam = *iterfam;
				map<Family*, vector<string> >::iterator data = family_output.find(fam);
				if(data == family_output.end()){
					vector<string> mydata;
					mydata.resize(all_columns.size(), "N/A");
					family_output[fam] = mydata;
				}
				for(int j = 0; j < filecols.size(); j++){
					vector<string>::iterator realcol = find(columns.begin(), columns.end(), filecols[j]);
					vector<string>::iterator allloc = find(all_columns.begin(), all_columns.end(), "(" +filename+")"+filecols[j]);
					if(realcol != columns.end()){
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
		for(int i = 0; i < all_columns.size(); i++){
			out << "\t" << all_columns[i];
		}
		out << "\n";
		map<Family*, vector<string> >::iterator iter;
		for(iter = family_output.begin(); iter != family_output.end(); iter++){
			Family* fam = iter->first;
			vector<string> data = iter->second;
			out << fam->getFamID() << "\t"
				<< fam->getSamples()->size();
			for(int i = 0; i < data.size(); i++){
				out << "\t" << data[i];
			}
			out << "\n";
		}
		out.close();
	}
	if(filenames["Sample"].size() > 0){
	map<Sample*, vector<string> > sample_output;
		vector<string> all_columns;
		for(int i = 0; i < filenames["Sample"].size(); i++){
			string file = filenames["Sample"][i];
			string step = opts::filesteps[file];

			for(int j = 0; j < opts::fileheaders[file].size(); j++){
				all_columns.push_back("(" + file + ")" + opts::fileheaders[file][j]);
			}
		}
		string outfile = opts::_OUTPREFIX_ + "sample_summary.txt";
		opts::printLog("Working on compiling Sample information...[" + outfile + "]\n");
		for(int i = 0; i < filenames["Sample"].size(); i++){
			string filename = filenames["Sample"][i];
			string step = opts::filesteps[filename];
			vector<string> filecols = opts::fileheaders[filename];

			opts::printLog("\tParsing " + filename + "\n");

	
			ifstream input;
			input.open(filename.c_str(), ios::in);
	
			if(!input){
				cerr << "Error opening file: " << filename << endl;
				exit(1);
			}
	
			string header = "";
		   	getline(input, header);

			vector<string> columns;
			General::Tokenize(header, columns, "\t");
			int famloc = -1;
			int indloc = -1;
			for(int l = 0; l < columns.size(); l++){
				if(columns[l] == "FamID"){
					famloc = l;
				}
				if(columns[l] == "IndID"){
					indloc = l;
				}
				if(famloc > -1 && indloc > -1){
					break;
				}
			}
			string line;
			int count = 1;
			while(getline(input, line)){
				count++;
				vector<string> elems;
				string token;
				istringstream isstream (line);
				while(getline(isstream, token, '\t')){
					elems.push_back(token);
				}
				//General::Tokenize(line, elems, "\t");
				if(elems.size() == 0){
					continue;
				}
				string famid = elems[famloc];
				string indid = elems[indloc];
				vector<Sample*>::iterator itersamp = find_if(samples->begin(), samples->end(), FindSampleByFamAndID(famid, indid));
				if(itersamp == samples->end()){
					cout << "Cannot find sample with famid = " << famid << " and indid = " << indid << endl;
					exit(1);
				}
				Sample* samp = *itersamp;
				map<Sample*, vector<string> >::iterator data = sample_output.find(samp);
				if(data == sample_output.end()){
					vector<string> mydata;
					mydata.resize(all_columns.size(), "N/A");
					sample_output[samp] = mydata;
				}
				for(int j = 0; j < filecols.size(); j++){
					vector<string>::iterator realcol = find(columns.begin(), columns.end(), filecols[j]);
					vector<string>::iterator allloc = find(all_columns.begin(), all_columns.end(), "(" +filename+")"+filecols[j]);
					if(realcol != columns.end()){
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
		for(int i = 0; i < all_columns.size(); i++){
			out << "\t" << all_columns[i];
		}
		out << "\n";
		map<Sample*, vector<string> >::iterator iter;
		for(iter = sample_output.begin(); iter != sample_output.end(); iter++){
			Sample* samp = iter->first;
			vector<string> data = iter->second;
			out << samp->getFamID() << "\t"
				<< samp->getInd() << "\t";
			if(samp->getSex()){
				out << "M\t";
			}
			else{
				out << "F\t";
			}
			out << samp->getPheno(); 
			
			for(int i = 0; i < data.size(); i++){
				out << "\t" << data[i];
			}
			out << "\n";
		}
		out.close();
	}
}

void readTPed(vector<Marker*>* markers, vector<Sample*>* samples, vector<int>* marker_map){
	map<string, vector<string> > descinfo;
	vector<string> descheaders;
	//vector<string> exclude;
	//vector<string> include;
	map<string,int> exclude;
	map<string,int> include;
	map<string,float> frequencies;
	vector<bool> includeme;

	//read in frequency file
	if(opts::_FREQ_FILE_EXISTS_){
		opts::printLog("Marker Minor Allele Frequencies being read from: " + opts::_FREQ_FILE_ + "\n");
		ifstream finput;
		finput.open(opts::_FREQ_FILE_.c_str(), ios::in);
		if(!finput){
			opts::printLog("Error opening marker frequency file: " + opts::_FREQ_FILE_ + "\n");
			exit(1);
		}
		string line = "";
		int freqline = 1;
		while(getline(finput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2){
				opts::printLog("Marker frequency file column size != 2 on line: " + line + " Skipping line!!\n");
				continue;
			}
			try{
				float freq = -1.0f;
				for(int c = 0; c < tokens[1].size(); c++){
					if(!isdigit(tokens[1][c]) && tokens[1][c] != '.'){
						throw "Number format exception";
					}
				}
				freq = std::atof(tokens[1].c_str());
				frequencies[tokens[0]] = freq;	
			}catch(...){
				opts::printLog(tokens[1] + " is not a valid value for line: " + getString<int>(freqline) + "\n");
				exit(1);
			}
			freqline++;
		}
	}
	//read in map descriptive file
	if(opts::_MAPDESC_.length() > 0){
		opts::printLog("Marker description being used: " + opts::_MAPDESC_ + "\n");
		ifstream dinput;
		dinput.open(opts::_MAPDESC_.c_str(), ios::in);
		if(!dinput){
			opts::printLog("Error opening marker description file: " + opts::_MAPDESC_ + "\n");
			exit(1);
		}
		string head = "";
		getline(dinput, head, '\n');
		descheaders = General::ParseDelimitedLine(head);
		while(!dinput.eof()){
			string id = "";
			string line = "";
			getline(dinput, line, '\n');
			if(line == ""){
				continue;
			}
			if(line[0] == '#'){
				continue;
			}

			vector<string> tokens;
			tokens = General::ParseDelimitedLine(line);
			if(tokens.size() < 2){
				opts::printLog("Marker description column size is < 2: " + line + "\n");
				exit(1);
			}
			if(tokens.size() != descheaders.size()){
				opts::printLog("Line is not the same size as header line: " + line + "\n");
				exit(1);
			}

			id = tokens[0];
			descinfo[id] = tokens;
		}
		dinput.close();
	}
	//read in marker exclusion file
	if(opts::_EXCMARKERS_.length() > 0){
		opts::printLog("Excluding markers from file: " + opts::_EXCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker exclusion file: " + opts::_EXCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";
		
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Excluding markers column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			exclude[tokens[0]] = 1;
			//exclude.push_back(probe);
		}
		einput.close();
	}
	//read in marker inclusion file
	if(opts::_INCMARKERS_.length() > 0){
		opts::printLog("Including markers from file: " + opts::_INCMARKERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCMARKERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening marker inclusion file: " + opts::_INCMARKERS_ + "\n");
			exit(1);
		}
		string probe = "";
		string line = "";
		while(getline(einput,line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1){
				opts::printLog("Including markers column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			include[tokens[0]] = 1;
			//include.push_back(probe);
		}
		einput.close();
	}
		
	
	FILE* input;
    input = fopen(opts::_TPEDFILE_.c_str(), "r");

    if(!input){
		opts::printLog("Error opening ped file: " + opts::_TPEDFILE_ + "\n");
        exit(1);
    }
	int count = 0;
	int prev_chrom = -1;
	int prev_bploc = -1;
    while(!feof(input)){
		int f = 0;
		string chr;
		string probe_id;
		string cm;
		string bp;
		if(readString(input, &chr)) f++;

		if(chr==""){
			continue;
		}

		if(readString(input, &probe_id)) f++;
		if(readString(input, &cm)) f++;
		if(readString(input, &bp)) f++;
		
		while(fgetc(input) != '\n' && !feof(input)){}

		int bploc = atoi(bp.c_str());

		bool use = true;
		opts::_MARKERS_FOUND_++;
		if(exclude.size() > 0){
			//vector<string>::iterator found = find(exclude.begin(), exclude.end(), probe_id);
			map<string, int>::iterator found = exclude.find(probe_id);
			if(found != exclude.end()){
				use = false;
			}
		}
		if(include.size() > 0){
			//vector<string>::iterator found = find(include.begin(), include.end(), probe_id);
			map<string, int>::iterator found = include.find(probe_id);
			if(found == include.end()){
				use = false;
			}
		}

        Marker* m = new Marker(chr, probe_id, bploc);
		if(opts::_CHROM_LIMIT_){
			if(m->getChrom() != opts::_CHROM_){
				use = false;
			}
			else{
				if(opts::_BP_LOW_LIMIT_){
					if(bploc < opts::_BP_LOW_){
						use = false;
					}
				}
				if(opts::_BP_HIGH_LIMIT_){
					if(bploc > opts::_BP_HIGH_){
						use = false;
					}
				}
			}
		}
		if(opts::_BP_SPACE_LIMIT_){
			if(prev_bploc == -1){
				prev_bploc = bploc;
				prev_chrom = m->getChrom();
			}
			else{
				if(m->getChrom() == prev_chrom && ((bploc - prev_bploc) < opts::_BP_SPACE_)){
					use = false;
				}
				else{
					prev_bploc = bploc;
					prev_chrom = m->getChrom();
				}
			}
		}
		includeme.push_back(use);
		if(!use){
	//		delete(m);
	//		continue;
		}
		m->setEnabled(use);
		m->setLoc(count);
		m->setRSID(probe_id);
		
	
		if(frequencies.size() > 0){
			map<string,float>::iterator found = frequencies.find(probe_id);
			if(found != frequencies.end()){
				m->setMAF(frequencies[probe_id]);
				m->setFreqFlag(true);
			}
		}	
		if(opts::_MAPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[probe_id];
			for(int i = 1; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					m->assignDetail(descheaders[i], tokens[i]);
				}
				else{
					m->assignDetail(descheaders[i], "NA");
				}
			}
		}
        markers->push_back(m);
		count++;
    }

    //input.clear();
    fclose(input);

	marker_map->resize(markers->size());

	//put markers in chrom/bploc order
	stable_sort(markers->begin(), markers->end(), less<Marker*>());

	for(int i =0; i < markers->size(); i++){
		(*marker_map)[(*markers)[i]->getLoc()] = i;//(*markers)[i]->getLoc();
	}

	count = 0;
	for(int i = 0; i < samples->size(); i++){
		(*samples)[i]->resizeAlleles(markers->size());
	}
	
	FILE* PED;
	PED = fopen(opts::_TPEDFILE_.c_str(), "r");
	int i =0;
	
	while(!feof(PED)){
		string dummy;
		string tempdummy;
		int f =0;
		if(readString(PED, &dummy)) f++;

		if(dummy==""){
			continue;
		}

		if(dummy.substr(0,1) == "#"){
			while(fgetc(PED) != '\n' && !feof(PED)){}
			continue;
		}

		if(!includeme[i]){
			while(fgetc(PED) != '\n' && !feof(PED)){}
			i++;
			continue;
		}

		if(readString(PED,&tempdummy)) f++; //probe
		if(readString(PED,&dummy)) f++; //cm
		if(readString(PED,&dummy)) f++; //bploc

		int gn = 0;
		int c=0; //ind count
		bool linedone = false;
		bool fatal = false;
		string fmsg;
		while(!linedone){
			Sample* samp = (*samples)[c];
			string one = "";
			string two = "";
			while(1){
				char ch = fgetc(PED);
				if(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(PED)){
					if(ch == '\n' || ch == '\r' || feof(PED)){
						linedone = true;
					}
					if(one.length() > 0){
						gn++;
						break;
					}
					if(ch == '\n' || ch == '\r' || feof(PED)){
						break;
					}
				}
				else{
					one += ch;
				}
			}

			if(!linedone){
				while(1){
					char ch = fgetc(PED);
					if(ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || feof(PED)){
						if(ch == '\n' || ch == '\r' || feof(PED)){
							linedone = true;
						}
						if(two.length() > 0){
							gn++;
							break;
						}
						if(ch == '\n' || ch == '\r' || feof(PED)){
							break;
						}
					}
					else{
						two += ch;
					}
				}
			}

			if(linedone && one.length() == 0 && two.length() == 0){
				break;
			}

			if(includeme[i]){
				int k = (*marker_map)[i];
				Marker* mark = (*markers)[k];

				int oldallelecount = mark->getNumAlleles();
                if(one != opts::_NOCALL_){
					//new
					if(mark->getAlleleLoc(one) < 0){
						mark->addAllele(one);
					}
        	    }

            	if(two != one){
                    if(two != opts::_NOCALL_){
						//new
						if(mark->getAlleleLoc(two) < 0){
							mark->addAllele(two);
						}
	                }
				}
//cout << mark->toString() << "\t" << mark->getAllele1() << mark->getAllele2() << "\t" << samp->toString() << "\t" << one << two;
				if(mark->getNumAlleles() <= 2){
   		            if(one == mark->getAllele1() && two == mark->getAllele1()){
   	         	        samp->addAone(i, false);
   	    	            samp->addAtwo(i, false);
		            }
	                else if(one != opts::_NOCALL_ && two != opts::_NOCALL_ && one != two){
	       	            samp->addAone(i, false);
 		      	        samp->addAtwo(i, true);
	               	}
		            else if(one == mark->getAllele2() && two == mark->getAllele2()){
	                    samp->addAone(i, true);
	       	            samp->addAtwo(i, true);
	           	    }
	               	else if(one == opts::_NOCALL_ || two == opts::_NOCALL_){
		                samp->addAone(i, true);
	                    samp->addAtwo(i, false);
	       	        }
//cout << "\t" << samp->getAone(i) << samp->getAtwo(i) << endl;
				}
				else if(opts::_MICROSATS_){
					samp->addMicroSat(i);
					int loc1 = mark->getAlleleLoc(one);
					int loc2 = mark->getAlleleLoc(two);
					//samp->addAbone(loc1);
					//samp->addAbtwo(loc2);
					samp->addAbone(i, loc1);
					samp->addAbtwo(i, loc2);
					if(oldallelecount <= 2){
						remapSamples(samples, markers, marker_map, i);
					}
				}
				else if(mark->getNumAlleles() > 2 && !opts::_MICROSATS_){
					opts::printLog("More than 2 unique alleles found for map location: " + getString<int>(k) + ", line: " + getString<int>(c + 1) + " [" + samp->toString() + ":::" + mark->toString() + " " + one + "/" + two + ".  Microsatellites not specified.\n");
					exit(1);
				}
			}
			c++;
			if(c > samples->size()){
				opts::printLog("Problem with line " + getString<int>(i+1) + " in " + opts::_TPEDFILE_ + "\n");
				opts::printLog("Expecting 4 + 2 * " + getString<int>(samples->size()) + " = " +
						getString<int>(4+2*samples->size()) + " columns, but found more\n");
				exit(1);
			}
		}//line done? next snp

		if(gn != 2 * samples->size()){
			opts::printLog("Problem with line " + getString<int>(i+1) + " in " + opts::_TPEDFILE_ + "\n");
			opts::printLog("Expecting 4 + 2 * " + getString<int>(samples->size()) + " = " +
					getString<int>(4+2*samples->size()) + " columns, but found more\n");
			exit(1);
		}
		i++;
	}
	fclose(PED);
}

void readTFam(vector<Sample*>* samples, vector<Family*>* families){
	map<string, vector<string> > descinfo;
	vector<string> sexclude;
	vector<string> inccenters;
	vector<string> sinclude;
	vector<string> finclude;
	vector<string> fexclude;
	vector<string> descheaders;
	if(opts::_SAMPDESC_.length() > 0){
		opts::printLog("Sample description being used: " + opts::_SAMPDESC_ + "\n");
		ifstream dinput;
		dinput.open(opts::_SAMPDESC_.c_str(), ios::in);
		if(!dinput){
			opts::printLog("Error opening sample description file: " + opts::_SAMPDESC_ + "\n");
			exit(1);
		}
        string head = "";
        getline(dinput, head, '\n');
        descheaders = General::ParseDelimitedLine(head);
		while(!dinput.eof()){
			string fam = "";
			string id = "";
			string center = "";
			string plate = "";
			string well = "";
//old        	char buf[256];
//old	        dinput.getline(buf, 256, '\n');
            string line = "";
			getline(dinput, line, '\n');

    	    //old string line = buf;

        	if(line == ""){
            	continue;
        	}
			vector<string> tokens = General::ParseDelimitedLine(line);
            if(tokens.size() < 2){
				opts::printLog("Sample description column size is < 2: " + line + " Skipping line!!\n");
	            continue;
	        }
	        if(tokens.size() != descheaders.size()){
				opts::printLog("Line is not the same size as header line: " + line + " Skipping line!!\n");
	            continue;
	        }
	
			fam = tokens[0];
			id = tokens[1];
			descinfo[fam + " " + id] = tokens;
//			center = tokens[2];
//			plate = tokens[3];
//			well = tokens[4];
//			descinfo[fam + " " + id] = center + " " + plate + " " + well;

//			descinfo.push_back(line);	
		}
		dinput.close();
	}
	if(opts::_EXCSAMPS_.length() > 0){
		opts::printLog("Excluding samples from file: " + opts::_EXCSAMPS_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCSAMPS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening sample exclusion file: " + opts::_EXCSAMPS_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2 && line != ""){
				opts::printLog("Sample exclusion column size != 2: " + line + " Skipping line!!\n");
				continue;
			}
			sexclude.push_back(tokens[0] + " " + tokens[1]);
		}
		einput.close();
	}
	if(opts::_EXCFAMILIES_.length() > 0){
		opts::printLog("Excluding families from file: " + opts::_EXCFAMILIES_ + "\n");
		ifstream einput;
		einput.open(opts::_EXCFAMILIES_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening family exclusion file: " + opts::_EXCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Family exclusion column size != 1: " + line + " Skipping line!!\n");
				continue;
			}
			fexclude.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCCENTERS_.length() > 0){
		opts::printLog("Including centers from file: " + opts::_INCCENTERS_ + "\n");
		ifstream einput;
		einput.open(opts::_INCCENTERS_.c_str(), ios::in);
		if(!einput){
			opts::printLog("Error opening center inclusion file: " + opts::_INCCENTERS_ + "\n");
			exit(1);
		}
		string center = "";
		string line = "";
		while(getline(einput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Center inclusion columns != 1: " + line + " Skipping line!!\n");
				continue;
			}
			inccenters.push_back(tokens[0]);
		}
		einput.close();
	}
	if(opts::_INCSAMPLES_.length() > 0){
		opts::printLog("Including samples from file: " + opts::_INCSAMPLES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCSAMPLES_.c_str(), ios::in);
		if(!iinput){
			opts::printLog("Error opening sample inclusion file: " + opts::_INCSAMPLES_ + "\n");
			exit(1);
		}
		string fam = "";
		string ind = "";
		string line = "";
		while(getline(iinput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 2 && line != ""){
				opts::printLog("Sample inclusion columns != 2: " + line + " Skipping line!!\n");
				continue;
			}
			sinclude.push_back(tokens[0] + " " + tokens[1]);
		}
		iinput.close();
	}
	if(opts::_INCFAMILIES_.length() > 0){
		opts::printLog("Including families from file: " + opts::_INCFAMILIES_ + "\n");
		ifstream iinput;
		iinput.open(opts::_INCFAMILIES_.c_str(), ios::in);
		if(!iinput){
			opts::printLog("Error opening family inclusion file: " + opts::_INCFAMILIES_ + "\n");
			exit(1);
		}
		string fam = "";
		string line = "";
		while(getline(iinput, line)){
			vector<string> tokens = General::ParseDelimitedLine(line);
			if(tokens.size() != 1 && line != ""){
				opts::printLog("Family inclusion columns != 1: " + line + " Skipping line!!\n");
				continue;
			}
			finclude.push_back(tokens[0]);
		}
		iinput.close();
	}
	
	opts::printLog("Reading family information from " + opts::_FAMFILE_ + "\n");

	ifstream PED;
	PED.open((opts::_FAMFILE_).c_str());
	if(!PED){
		opts::printLog("Error opening family information file: " + opts::_FAMFILE_ + ".  Exiting!\n");
		exit(1);
	}
	PED.clear();
		string fam = "";
		string ind = "";
		string dad = "";
		string mom = "";
		string sex = "";
		string aff = "";
		string line = "";
	while(getline(PED, line)){
	   	//	>> fam
		//	>> ind
		//	>> dad
		//	>> mom
		//	>> sex
		//	>> aff){
		vector<string> tokens = General::ParseDelimitedLine(line);
		if(tokens.size() != 6){
			opts::printLog("Family information file column size != 6: " + line + " Exitting!!\n");
			exit(0);
		}
		fam = tokens[0];
		ind = tokens[1];
		dad = tokens[2];
		mom = tokens[3];
		sex = tokens[4];
		aff = tokens[5];
		Sample* samp = new Sample();
		samp->setFamID(fam);
		samp->setInd(ind);
		samp->setEnabled(true);
		opts::_SAMPLES_FOUND_++;
		if(sexclude.size() > 0){
			vector<string>::iterator found = find(sexclude.begin(), sexclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found != sexclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				if(opts::_KEEP_EXC_SAMPLES_){
					samp->setExcluded(true);
				}
			}
		}
		if(sinclude.size() > 0){
			vector<string>::iterator found = find(sinclude.begin(), sinclude.end(), samp->getFamID() + " " + samp->getInd());
			if(found == sinclude.end()){
				samp->setEnabled(false);
				if(opts::_KEEP_EXC_SAMPLES_){
					samp->setExcluded(true);
				}
			}
		}

		samp->setDadID(dad);
		samp->setMomID(mom);
		if(sex == "1"){
			samp->setSex(true);
		}
		else if(sex == "2"){
			samp->setSex(false);
		}

		if(aff == "2"){
			samp->setAffected(true);
		}
		else{
			samp->setAffected(false);
		}
		samp->setPheno(atoi(aff.c_str()));
		if(opts::pedinfo.size() > 0){
			map<string, Sample*>::iterator sfind = opts::pedinfo.find(samp->getFamID() + "#" + samp->getInd());
			if(sfind != opts::pedinfo.end()){
				Sample* sfound = sfind->second;
				samp->setDadID(sfound->getDadID());
				samp->setMomID(sfound->getMomID());
				samp->setSex(sfound->getSex());
				samp->setPheno(sfound->getPheno());
			}
		}
		
		string center = "";
		if(opts::_SAMPDESC_.length() > 0 && descinfo.size() > 0){
			vector<string> tokens = descinfo[samp->getFamID() + " " + samp->getInd()];
			for(int i = 2; i < descheaders.size(); i++){
				if(tokens.size() == descheaders.size()){
					samp->assignDetail(descheaders[i], tokens[i]);
				}
				else{
				    samp->assignDetail(descheaders[i], "NA");
				}
			}

			//if(tokens.size() > 0){
			//	center = tokens[0];
			//}
			//if(opts::_INCCENTERS_.length() > 0 && inccenters.size() > 0){
			//	bool found = false;
			//	for(int cent = 0; cent < inccenters.size(); cent++){
			//		if(inccenters[cent] == center){
			//			found = true;
			//		}
			//	}
			//	if(!found){
			//		samp->setEnabled(false);
			//	}
//
			//}
			//samp->getFamily()->setCenter(tokens[2]);
			//if(tokens.size() > 0){
			//	samp->setPlate(tokens[1]);
			//	samp->setWell(tokens[2]);
			//}
		}
        vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(),FindFamily(samp->getFamID()));

        if(f_iter != (*families).end()){
            (*f_iter)->AddInd(samp);
            samp->setFamily((*f_iter));
        }
        else{
            Family* fam = new Family();
/*          fam->Setcenter()*/
            fam->setFamID(samp->getFamID());
            fam->AddInd(samp);
			fam->setCenter(center);
			fam->setEnabled(true);
            samp->setFamily(fam);
            families->push_back(fam);
			fam->setLoc((families->size() - 1));
			opts::_FAMILIES_FOUND_++;
        }
		if(fexclude.size() > 0){
			vector<string>::iterator found = find(fexclude.begin(), fexclude.end(), samp->getFamID());
			if(found != fexclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		if(finclude.size() > 0){
			vector<string>::iterator found = find(finclude.begin(), finclude.end(), samp->getFamID());
			if(found == finclude.end()){
				//cout << "Disabling sample: " << samp->getFamID() << "\t" << samp->getInd() << endl;
				samp->setEnabled(false);
				vector<Family*>::iterator f_iter = find_if(families->begin(), families->end(), FindFamily(samp->getFamID()));
				if(f_iter != (*families).end()){
					(*f_iter)->setEnabled(false);
				}
			}
		}
		//samp->resizeAlleles(markers->size());
        samples->push_back(samp);
		samp->setLoc((samples->size() - 1));
	}	

	PED.clear();
	PED.close();
}




//////////////////////////////////////////////////////////////////
// Borrowed from Plink v0.99s                                   //
//                                                              //
//////////////////////////////////////////////////////////////////


#define  PORT_NUM                80     
#define  IP_ADDR    "160.129.37.40"  
#define  GET_STRING "GET /wasp/files/version.txt HTTP/1.1\nHost: chgr.mc.vanderbilt.edu\nConnection: close\n\n"


void webcheck(vector<string> a, map<string, vector<string> > b)
{

#ifdef SKIP
	opts::printLog("Web-check not implemented on this system...\n");
  return;
#else
  
  opts::printLog("Web-based version check ( -noweb to skip )\n");
  
  vector<string> tokens = socketConnection( 
					    IP_ADDR,
					    PORT_NUM,
					    GET_STRING);
					    
  bool print = false;
  bool print2 = false;
  bool version_okay = true;
  vector<string>::iterator fiter;
  map<string, vector<string> >::iterator siter;

  for (int i=0; i<tokens.size(); i++)
    {

      if (tokens[i]=="END") break;

      if (tokens[i]=="END-MESSAGE")
	{
	  print2=false;
	  continue;
	}

      if (tokens[i]=="WARN")
	{
	  if ( i < tokens.size()-1 ) 
	    {
	      i++;
		  fiter = find(a.begin(), a.end(), tokens[i]);
	      if ( fiter != a.end())
		{
		  opts::printLog("\n*** ALERT ***\n*** A warning flag has been set for: "+tokens[i]+
			   "\n*** See http://chgr.mc.vanderbilt.edu/wasp/\n");
		}
	    }
	  continue;
	}

	  if(tokens[i] == "WARNSTEP"){
	  	if(i < tokens.size() - 1){
			i++;
			siter = b.find(tokens[i]);
			if(siter != b.end()){
				opts::printLog("\n*** ALERT ***\n*** A warning flag has been set for STEP: " + tokens[i] + 
						"\n*** See http://chgr.mc.vanderbilt.edu/wasp/\n");
			}
		}
		continue;
	  }

	  if(tokens[i] == "WARNSTEPOPTION"){
		  if(i < tokens.size() - 1){
			  i++;
			  for(siter = b.begin(); siter != b.end(); siter++){
				  vector<string> second = siter->second;
				  fiter = find(second.begin(), second.end(), tokens[i]);
				  if(fiter != second.end()){
						opts::printLog("\n*** ALERT ***\n*** A warning flag has been set for STEP OPTION: " + tokens[i] + 
							"\n*** See http://chgr.mc.vanderbilt.edu/wasp/\n");
				  	break;
				  }
			  }
		  }
		  continue;
	  }
      
      if (tokens[i]=="FATAL")
	{
	  if ( i < tokens.size()-1 ) 
	    {
	      i++;
		  fiter = find(a.begin(), a.end(), tokens[i]);
	      if ( fiter != a.end()){
		opts::printLog("A serious warning flag has been set for: "+tokens[i]+
		    "\nWasp has been instructed to stop"+
 	            "\nPlease see http://chgr.mc.vanderbilt.edu/wasp/\n");
		  	exit(1);
		  }
	    }
	  continue;
	}

	  if(tokens[i] == "FATALSTEP"){
	  	if(i < tokens.size() - 1){
			i++;
			siter = b.find(tokens[i]);
			if(siter != b.end()){
		opts::printLog("A serious warning flag has been set for STEP: "+tokens[i]+
		    "\nWasp has been instructed to stop"+
 	            "\nPlease see http://chgr.mc.vanderbilt.edu/wasp/\n");
		  	exit(1);
			}
		}
		continue;
	  }

	  if(tokens[i] == "FATALSTEPOPTION"){
		  if(i < tokens.size() - 1){
			  i++;
			  for(siter = b.begin(); siter != b.end(); siter++){
				  vector<string> second = siter->second;
				  fiter = find(second.begin(), second.end(), tokens[i]);
				  if(fiter != second.end()){
		opts::printLog("A serious warning flag has been set for STEP OPTION: "+tokens[i]+
		    "\nWasp has been instructed to stop"+
 	            "\nPlease see http://chgr.mc.vanderbilt.edu/wasp/\n");
		  	exit(1);
				  }
			  }
		  }
		  continue;
	  }

      if (tokens[i]=="MESSAGE-ALL")
	{
	  print2=true;
	  continue;
	}

      // Display any other messages
      // Either conditional on old version (print)
      // or a broadcast to all users (print2)

      if ( ( print && !version_okay) || print2 ) 
	{
	  if (tokens[i]=="\\n")
	    opts::printLog("\n");
	  else
	    opts::printLog(tokens[i]+" ");
	}

      // Check version code
      if (tokens[i]=="WASPVER") 
	{
	  print=true;
	  if ( i < tokens.size() - 1) 
	    {
	      if (tokens[i+1] == opts::_WASPVER_)
		opts::printLog(" OK, v"+opts::_WASPVER_+" is current\n");
	      else
		{
		  opts::printLog("\n\n          *** UPDATE REQUIRED ***\n\n");
		  opts::printLog("\tThis version        : "+opts::_WASPVER_+"\n");
		  opts::printLog("\tMost recent version : "+tokens[i+1]+"\n\n");
		  opts::printLog("Please upgrade your version of Wasp as soon as possible!\n"); 
		  opts::printLog("  (visit the above website for free download)\n\n");
		  version_okay=false;
		}

	      // Skip the version number
	      i++;
	    }
	}

    }

  // did we get the information we needed?
  if (!print) opts::printLog(" problem connecting to web\n");

  opts::printLog("\n");
#endif  
}

