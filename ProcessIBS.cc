/**********************************************************************************
*                       Marker Genotype Efficiency Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all genotypes and generates a genotype efficiency for all markers.
*
*
* Files generated:
*	percent_breakdown_by_marker.txt
*	percent_breakdown_by_chrom.txt
*       post_marker_geno_eff_filter_summary.txt
*
*File: PartialOutput.cc
**********************************************************************************/


#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include "ProcessIBS.h"
#include "Chrom.h"
#include <General.h>
#include <Helpers.h>
using namespace Methods;
string ProcessIBS::stepname = "ibs";

void ProcessIBS::FilterSummary(){
}

void ProcessIBS::PrintSummary(){
	int msize = data_set->num_loci();

	for(int i = 0; i < msize; i++){
		data_set->get_locus(i)->setFlag(false);
	}

}

void ProcessIBS::filter(){
}


void ProcessIBS::process(DataSet* ds){
	data_set = ds;
	vector<Marker*> good_markers = Helpers::findValidMarkers(data_set->get_markers(), &options);

	IBS ibs;
	ibs.resetDataSet(data_set);
	ibs.set_parameters(&options);

	int ssize = data_set->num_inds();
////	int msize = data_set->num_loci();

	if (options.getDoIBSAllPairs()) {
		string filename = opts::_OUTPREFIX_ + "ibs" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if (!overwrite) {
			filename += "." + getString<int> (order);
		}

		ofstream out(filename.c_str());
		if (!out.is_open()) {
			opts::printLog("Unable to open " + filename + "\n");
			throw MethodException("");
		}

		out.precision(4);
		out << "FamID1\tInd1\tFamID2\tInd2\tNum_Comps\tIBS_AVG\n";
		for (int s1 = 0; s1 < ssize - 1; s1++) {
			for (int s2 = s1 + 1; s2 < ssize; s2++) {
				Sample* samp1 = data_set->get_sample(s1);
				Sample* samp2 = data_set->get_sample(s2);
				double avg = ibs.calcPairAverage(s1, s2);
				int comps = ibs.getComparisons();
				out << samp1->getFamID() << "\t" << samp1->getInd() << "\t";
				out << samp2->getFamID() << "\t" << samp2->getInd() << "\t";
				out << comps << "\t" << avg << endl;
			}
		}
		if (out.is_open()) {
			out.close();
		}
	}

	if(options.getDoIBSPairs()){
		string filename2 = opts::_OUTPREFIX_ + "ibs_raw" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if(!overwrite){
			filename2 += "." + getString<int>(order);
		}

		ofstream rawout (filename2.c_str());
		if(!rawout.is_open()){
			opts::printLog("Unable to open " + filename2 + "\n");
			throw MethodException("");
		}

		rawout.precision(4);
		rawout << "SNP\tFamID1\tInd1\tFamID2\tInd2\tIBS\n";
		vector<Sample*>* samps = data_set->get_samples();
		map<string, vector<string> > pairs = options.get_ibs_pairs();
		map<string, vector<string> >::iterator iter;
		int msize = good_markers.size();//data_set->num_loci();
		for(iter = pairs.begin(); iter != pairs.end(); iter++){
			string s1 = (string)iter->first;
			vector<string> tokens = General::ParseDelimitedLine(s1);

			vector<string> comps = (vector<string>)iter->second;
			vector<Sample*>::iterator found = find_if(samps->begin(), samps->end(), FindSampleByFamAndID(tokens[0], tokens[1]));

			int samp1loc = -1;
			int samp2loc = -1;
			if(found != samps->end()){
				samp1loc = found - samps->begin();
			}
			if(samp1loc == -1){
				continue;
			}
			for(int s = 0; s < (int)comps.size(); s++){
				string s2 = comps[s];
				vector<string> toks = General::ParseDelimitedLine(s2);
				vector<Sample*>::iterator found2 = find_if(samps->begin(), samps->end(), FindSampleByFamAndID(toks[0], toks[1]));

				if(found2 != samps->end()){
					samp2loc = found2 - samps->begin();
				}
				if(samp2loc == -1){
					continue;
				}
				for(int m = 0; m < msize; m++){
					int value = ibs.calcPairLocus(samp1loc, samp2loc, good_markers[m]);
					rawout << good_markers[m]->getProbeID() << "\t" << s1 << "\t" << s2 << "\t" << value << endl;
				}
			}
		}
		if(rawout.is_open()){
			rawout.close();
		}
	}

	if(options.getDoIBSTrioPairs() || options.getDoIBSAllTrioPairs()){
		string filename2 = opts::_OUTPREFIX_ + "ibs_trios" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
		if(!overwrite){
			filename2 += "." + getString<int>(order);
		}

		ofstream rawout (filename2.c_str());
		if(!rawout.is_open()){
			opts::printLog("Unable to open " + filename2 + "\n");
			throw MethodException("");
		}

		rawout.precision(4);
		rawout << "FamID1\tFamID2\tPaternal_avg\tMaternal_avg\tITBS_avg\tNsnps\n";

		ofstream fullraw;
		if(options.getIbsTriosRaw()){
			string filename3 = opts::_OUTPREFIX_ + "ibs_trios_raw" + options.getOut() + ".txt";
			if(!overwrite){
				filename3 += "." + getString<int>(order);
			}
			fullraw.open(filename3.c_str());
			if(!fullraw.is_open()){
				opts::printLog("Unable to open " + filename3 + "\n");
				throw MethodException("");
			}
			fullraw << "SNP\tFamID1\tFamID2\tPaternal\tMaternal\tIBS\n";
		}

		vector<Family*>* fams = data_set->get_families();
		map<string, vector<string> > pairs = options.get_ibs_trio_pairs();
		map<string, vector<string> >::iterator iter;

		int msize = good_markers.size();//data_set->num_loci();
		if(pairs.size() > 0){
			for(iter = pairs.begin(); iter != pairs.end(); iter++){
				string s1 = (string)iter->first;
				vector<string> tokens = General::ParseDelimitedLine(s1);

				vector<string> comps = (vector<string>)iter->second;
				vector<Family*>::iterator found = find_if(fams->begin(), fams->end(), FindFamily(tokens[0]));

				int fam1loc = -1;
				int fam2loc = -1;
				if(found != fams->end()){
					fam1loc = found - fams->begin();
				}
				if(fam1loc == -1){
					continue;
				}
				for(int s = 0; s < (int)comps.size(); s++){
					string s2 = comps[s];
					vector<string> toks = General::ParseDelimitedLine(s2);
					vector<Family*>::iterator found2 = find_if(fams->begin(), fams->end(), FindFamily(toks[0]));

					if(found2 != fams->end()){
						fam2loc = found2 - fams->begin();
					}
					if(fam2loc == -1){
						continue;
					}
//					int prev_bploc = 0;
//					int prev_chrom = -1;
					double paternal = 0;
					double maternal = 0;
					int total_snps = 0;
					for(int m = 0; m < msize; m++){
						Marker* mark = good_markers[m];//data_set->get_locus(m);
						if(mark->isEnabled() && mark->getChrom() != opts::_CHRX_ &&
								mark->getChrom() != opts::_CHRY_ && mark->getChrom() != opts::_CHRXY_){// &&
//								isValidMarker(mark, &options, prev_bploc, prev_chrom)){
							vector<double> value = ibs.calcTriosLocus(fam1loc, fam2loc, mark);
							paternal += value[0];
							maternal += value[1];
							total_snps++;

							if(options.getIbsTriosRaw()){
								fullraw << mark->getRSID() << "\t"
									    << data_set->get_pedigree(fam1loc)->getFamID() << "\t"
									    << data_set->get_pedigree(fam2loc)->getFamID() << "\t"
									    << value[0] << "\t"
									    << value[1] << "\t"
									    << (value[0] + value[1])
									    << endl;
							}
						}
					}

					rawout << data_set->get_pedigree(fam1loc)->getFamID() << "\t" << data_set->get_pedigree(fam2loc)->getFamID() << "\t";
					rawout << (paternal / total_snps) << "\t" << (maternal / total_snps) << "\t" << ((paternal + maternal) / total_snps)
						<< "\t" << total_snps << endl;
				}
			}
		}
		else{
			int fsize = fams->size();

			for(int f1 = 0; f1 < fsize - 1; f1++){
				for(int f2 = f1 + 1; f2 < fsize; f2++){
					Family* fam1 = data_set->get_pedigree(f1);
					Family* fam2 = data_set->get_pedigree(f2);
//					int prev_bploc = 0;
//					int prev_chrom = -1;
					double paternal = 0;
					double maternal = 0;
					int total_snps = 0;
					for(int m = 0; m < msize; m++){
						Marker* mark = good_markers[m];//data_set->get_locus(m);
						if(mark->isEnabled()){// && isValidMarker(mark, &options, prev_bploc, prev_chrom)){
							vector<double> value = ibs.calcTriosLocus(f1, f2, mark);
							paternal += value[0];
							maternal += value[1];
							total_snps++;
							if(options.getIbsTriosRaw()){
								fullraw << mark->getRSID() << "\t"
									    << fam1->getFamID() << "\t"
									    << fam2->getFamID() << "\t"
									    << value[0] << "\t"
									    << value[1] << "\t"
									    << (value[0] + value[1])
									    << endl;
							}

						}
					}

					rawout << fam1->getFamID() << "\t" << fam2->getFamID() << "\t";
					rawout << (paternal / total_snps) << "\t" << (maternal / total_snps) << "\t" << ((paternal + maternal) / total_snps)
						<< "\t" << total_snps << endl;
				}
			}
		}
		if(rawout.is_open()){
			rawout.close();
		}
		if(fullraw.is_open()){
			fullraw.close();
		}
	}

}

