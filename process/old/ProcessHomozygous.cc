/**********************************************************************************
*                       Run of Homozygosity Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Performs a run of homozygosity search based on Lencz et al (PNAS 0710021104)
* Performs a ROH based on straight forward common homozygous spans.
*
*
*File: Homozygous.cc
**********************************************************************************/


#include "ProcessHomozygous.h"
#include <Options.h>
#include <Homozygous.h>
#include <MethodException.h>
#include <Helpers.h>

using std::string;
using std::ofstream;
using std::getString;

using Methods::opts;
using Methods::Marker;
using Methods::Homozygous;
using Methods::MethodException;
using Methods::DataSet;
using Methods::Helpers;

const string ProcessHomozygous::stepname = ProcessHomozygous::doRegister("homozygous");

void ProcessHomozygous::PrintSummary(){
	if(options.doHomozygPermute()){
		return;
	}
	if(!options.doHomozygWGHA()){
		string filename = opts::_OUTPREFIX_ + "homozygous_marker" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if(!overwrite){
			filename += "." + getString<int>(order);
		}
		ofstream homo (filename.c_str());
		if(!homo.is_open()){
			opts::printLog("Unable to open " + filename + "\n");
			throw MethodException("");
		}
		opts::addFile("Marker", stepname, filename);
		homo << "Chrom\trsID\tProbeID\tbploc\tUnAff\tAff\tTotal_Maj\tTotal_Min\tTotal\tTotal%\n";
		int msize = good_markers.size();//data_set->num_loci();
		opts::addHeader(filename, "UnAff");
		opts::addHeader(filename, "Aff");
		opts::addHeader(filename, "Total_Maj");
		opts::addHeader(filename, "Total_Min");
		opts::addHeader(filename, "Total");
		opts::addHeader(filename, "Total%");

		for(int i = 0; i < msize; i++){
			Marker* m = good_markers[i];//data_set->get_locus(i);
			float per = (((float)homoallcount[i] / (float)opts::_SAMPLES_WORKING_));// * 100.0f);
			homo.precision(4);
			homo << m->getChrom() << "\t" << m->getRSID() << "\t" << m->getProbeID() << "\t"
				<< m->getBPLOC() << "\t" << homounaffcount[i] << "\t" << homoaffcount[i] << "\t"
				<< homomajallcount[i] << "\t" << homominallcount[i] << "\t"
				<< (homounaffcount[i] + homoaffcount[i]) << "\t" << per << endl;
		}
		if(homo.is_open()){
			homo.close();
		}
	}

}

void ProcessHomozygous::process(DataSet* ds){
	data_set = ds;
	good_markers = Helpers::findValidMarkers(data_set->get_markers(), &options);

	Homozygous hom(data_set);
	hom.setOptions(options);
	hom.setOverwrite(this->overwrite);
	hom.setOrder(this->order);
	hom.calculate();

	if(!options.doHomozygPermute() && !options.doHomozygWGHA()){
		homounaffcount = hom.getHomoUnaffCount();
		homoaffcount = hom.getHomoAffCount();
		homomajallcount = hom.getHomoMajAllCount();
		homominallcount = hom.getHomoMinAllCount();
		homoallcount = hom.getHomoAllCount();
	}

}
