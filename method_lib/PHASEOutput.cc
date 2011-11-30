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
*File: PHASEOutput.cc
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
#include <iomanip>
#include <string>
#include <list>
#include <algorithm>
#include <map>
#include "PHASEOutput.h"
#include "General.h"
#include "Helpers.h"

namespace Methods{
void PHASEOutput::FilterSummary(){
}

void PHASEOutput::PrintSummary(){
	int msize = markers->size();
	int ssize = samples->size();

	for(int i = 0; i < ssize; i++){
		(*samples).at(i)->setFlag(false);
	}

	for(int i = 0; i < msize; i++){
		(*markers).at(i)->setFlag(false);
	}

}

void PHASEOutput::filter(){
}

void PHASEOutput::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

   	int ssize = samples->size();
	int msize = markers->size();

	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	msize = good_markers.size();

	vector<int> chrom_counts;
	chrom_counts.resize(26,0);
	for(int i = 0; i < msize; i++){
		Marker* mark = good_markers.at(i);
		if(mark->isEnabled()){
			chrom_counts.at(mark->getChrom() - 1)++;
		}
	}
	int parents = 0;
	int stotal = 0;
	vector<bool> samp_flags(ssize, false);
	for(int i = 0; i < ssize; i++){
		Sample* samp = (*samples).at(i);
		int sloc = samp->getLoc();
		if((samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) ||
				(!samp->isEnabled() && options.doIncDisabledSamples())) && !samp_flags.at(sloc)){
			if(samp->getDad() != NULL && samp->getMom() != NULL && (samp->getDad()->isEnabled() ||
					(samp->getDad()->isExcluded() && options.doIncExcludedSamples()) ||
					(!samp->getDad()->isEnabled() && options.doIncDisabledSamples())) &&
					(samp->getMom()->isEnabled() || (samp->getMom()->isExcluded() && options.doIncExcludedSamples()) ||
					(!samp->getMom()->isEnabled() && options.doIncDisabledSamples())) &&
					!samp_flags.at(samp->getDad()->getLoc()) && !samp_flags.at(samp->getMom()->getLoc()) &&
					samp->getFamily()->getTotalInds() == 3){
				samp_flags.at(samp->getDad()->getLoc()) = true;
				samp_flags.at(samp->getMom()->getLoc()) = true;
				parents += 2;
			}
		}
		if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
			stotal++;
		}
	}

	string marker_types = "";

	for(int c = 0; c < (int)chrom_counts.size(); c++){
		if(chrom_counts.at(c) > 0){
			if(options.getChrom() == -1 || (options.getChrom() == (c + 1))){
				string fname1 = opts::_OUTPREFIX_ + "input_phase_chr" + getString<int>(c + 1) + options.getOut() + ".txt";
				if(options.getOverrideOut().size() > 0){
					fname1 = options.getOverrideOut() + "_chr" + getString<int>(c + 1) + ".txt";
				}
				if(!overwrite){
					fname1 += "." + getString<int>(order);
				}
				filenames.push_back(fname1);
				ofstream str (fname1.c_str());
				if(!str.is_open()){
					opts::printLog("Unable to open " + fname1 + " for output!\n");
					throw MethodException("Unable to open " + fname1 + " for output!\n");
				}
				fname1 = opts::_OUTPREFIX_ + "input_phase_chr" + getString<int>(c+ 1) + "_map" + options.getOut() + ".txt";
				if(options.getOverrideOut().size() > 0){
					fname1 = options.getOverrideOut() + "_chr" + getString<int>(c + 1) + "_map.txt";
				}
				if(!overwrite){
					fname1 += "." + getString<int>(order);
				}
				filenames.push_back(fname1);
				ofstream str_map (fname1.c_str());
				if(!str_map.is_open()){
					opts::printLog("Unable to open " + fname1 + " for output!\n");
					throw MethodException("Unable to open " + fname1 + " for output!\n");
				}
				if(options.doTriosOnly()){
					str << parents << "\n";
				}
				else{
					str << stotal << "\n";
				}
				str << chrom_counts.at(c) << "\n";
				str << "P";
				marker_types = "";
				for(int m = 0; m < msize; m++){
					Marker* mark = good_markers.at(m);
					if(mark->getChrom() == (c + 1) && mark->isEnabled()){
						if(!mark->isMicroSat()){
							marker_types += "S";
						}
						else{
							marker_types += "M";
						}
						str << " " << mark->getBPLOC();
						str_map << mark->getChrom() << "\t" << mark->getRSID() << "\t" << mark->getBPLOC() << "\n";
					}
				}
				str_map.close();
				str << "\n";
				str << marker_types << "\n";


				if(!options.doTriosOnly()){
					for(int s = 0; s < ssize; s++){
						Sample* samp = (*samples).at(s);
						string a1 = "";
						string a2 = "";
						if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
							str << samp->getFamID() << "_" << samp->getInd() << "\n";
							for(int m = 0; m < msize; m++){
								Marker* mark = good_markers.at(m);
								if(mark->getChrom() == (c + 1) && mark->isEnabled()){
									int mloc = mark->getLoc();
									if(!mark->isMicroSat()){
										if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
											a1 += "? ";
											a2 += "? ";
										}
										else if(!samp->getAone(mloc) && !samp->getAtwo(mloc)){
											a1 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
											a2 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
										}
										else if(!samp->getAone(mloc) && samp->getAtwo(mloc)){
											a1 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
											a2 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
										}
										else if(samp->getAone(mloc) && samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
											a1 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
											a2 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
										}
										else{
											a1 += "? ";
											a2 += "? ";
										}
									}//end !microsat
									else{
										if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
											a1 += "-1 ";
											a2 += "-1 ";
										}
										else if(samp->getAbone(mloc) == -1){
											a1 += "-1 ";
											a2 += "-1 ";
										}
										else{
											a1 += Helpers::map_allele(mark, mark->getAllele(samp->getAbone(mloc)), &options) + " ";
											a2 += Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(mloc)), &options) + " ";
										}
									}//end microsat
								}
							}
							a1 = a1.erase((a1.size() - 1), 1);
							a2 = a2.erase((a2.size() - 1), 1);

							str << a1 << "\n" << a2 << "\n";
						}
					}
				}
				else{
					//reset flags
					for(int s = 0; s < ssize; s++){
						samp_flags.at(s) = false;
					}
					//parents first
					for(int s = 0; s < ssize; s++){
						Sample* samp = (*samples).at(s);
						if((samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) ||
									(!samp->isEnabled() && options.doIncDisabledSamples())) && samp->getDad() != NULL && samp->getMom() != NULL && (samp->getDad()->isEnabled() ||
									(samp->getDad()->isExcluded() && options.doIncExcludedSamples()) || (!samp->getDad()->isEnabled() && options.doIncDisabledSamples())) &&
									(samp->getMom()->isEnabled() || (samp->getMom()->isExcluded() && options.doIncExcludedSamples()) ||
									(!samp->getMom()->isEnabled() && options.doIncDisabledSamples())) && samp->getFamily()->getTotalInds() == 3)
						{
							Sample* dad = samp->getDad();
							Sample* mom = samp->getMom();

							string da1 = "";
							string da2 = "";
							string ma1 = "";
							string ma2 = "";

							for(int m = 0; m < msize; m++){
								Marker* mark = good_markers.at(m);
								if(mark->getChrom() == (c + 1) && mark->isEnabled()){
									int mloc = mark->getLoc();
									if(!mark->isMicroSat()){
										if((dad->isExcluded() && options.doZeroExcluded()) || (!dad->isEnabled() && options.doIncDisabledSamples())){
											da1 += "? ";
											da2 += "? ";
										}
										else if(!dad->getAone(mloc) && !dad->getAtwo(mloc)){
											da1 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
											da2 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
										}
										else if(!dad->getAone(mloc) && dad->getAtwo(mloc)){
											da1 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
											da2 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
										}
										else if(dad->getAone(mloc) && dad->getAtwo(mloc) && !dad->getAmissing(mloc)){
											da1 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
											da2 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
										}
										else{
											da1 += "? ";
											da2 += "? ";
										}

										if((mom->isExcluded() && options.doZeroExcluded()) || (!mom->isEnabled() && options.doIncDisabledSamples())){
											ma1 += "? ";
											ma2 += "? ";
										}
										else if(!mom->getAone(mloc) && !mom->getAtwo(mloc)){
											ma1 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
											ma2 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
										}
										else if(!mom->getAone(mloc) && mom->getAtwo(mloc)){
											ma1 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
											ma2 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
										}
										else if(mom->getAone(mloc) && mom->getAtwo(mloc) && !mom->getAmissing(mloc)){
											ma1 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
											ma2 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
										}
										else{
											ma1 += "? ";
											ma2 += "? ";
										}
									}//end !microsat
									else{
										if((dad->isExcluded() && options.doZeroExcluded()) || (!dad->isEnabled() && options.doIncDisabledSamples())){
											da1 += "-1 ";
											da2 += "-1 ";
										}
										else if(dad->getAbone(mloc) == -1){
											da1 += "-1 ";
											da2 += "-1 ";
										}
										else{
											da1 += Helpers::map_allele(mark, mark->getAllele(dad->getAbone(mloc)), &options) + " ";
											da2 += Helpers::map_allele(mark, mark->getAllele(dad->getAbtwo(mloc)), &options) + " ";
										}
										if((mom->isExcluded() && options.doZeroExcluded()) || (!mom->isEnabled() && options.doIncDisabledSamples())){
											ma1 += "-1 ";
											ma2 += "-1 ";
										}
										else if(mom->getAbone(mloc) == -1){
											ma1 += "-1 ";
											ma2 += "-1 ";
										}
										else{
											ma1 += Helpers::map_allele(mark, mark->getAllele(mom->getAbone(mloc)), &options) + " ";
											ma2 += Helpers::map_allele(mark, mark->getAllele(mom->getAbtwo(mloc)), &options) + " ";
										}
									}//end microsat
								}
							}
							da1 = da1.erase((da1.size() - 1), 1);
							da2 = da2.erase((da2.size() - 1), 1);
							ma1 = ma1.erase((ma1.size() - 1), 1);
							ma2 = ma2.erase((ma2.size() - 1), 1);

							str << mom->getFamID() << "_" << mom->getInd() << "\n";
							str << ma1 << "\n" << ma2 << "\n";
							str << dad->getFamID() << "_" << dad->getInd() << "\n";
							str << da1 << "\n" << da2 << "\n";
							samp_flags.at(samp->getLoc()) = true;
						}
					}

					//do children
					for(int s = 0; s < ssize; s++){
						Sample* samp = (*samples).at(s);
						if(samp_flags.at(samp->getLoc())){
							string a1 = "";
							string a2 = "";
							str << samp->getFamID() << "_" << samp->getInd() << "\n";
							for(int m = 0; m < msize; m++){
								Marker* mark = good_markers.at(m);
								if(mark->getChrom() == (c + 1) && mark->isEnabled()){
									int mloc = mark->getLoc();
									if(!mark->isMicroSat()){
										if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
											a1 += "? ";
											a2 += "? ";
										}
										else if(!samp->getAone(mloc) && !samp->getAtwo(mloc)){
											a1 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
											a2 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
										}
										else if(!samp->getAone(mloc) && samp->getAtwo(mloc)){
											a1 += Helpers::map_allele(mark, mark->getAllele1(), &options) + " ";
											a2 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
										}
										else if(samp->getAone(mloc) && samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
											a1 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
											a2 += Helpers::map_allele(mark, mark->getAllele2(), &options) + " ";
										}
										else{
											a1 += "? ";
											a2 += "? ";
										}
									}//end !microsat
									else{
										if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
											a1 += "-1 ";
											a2 += "-1 ";
										}
										else if(samp->getAbone(mloc) == -1){
											a1 += "-1 ";
											a2 += "-1 ";
										}
										else{
											a1 += Helpers::map_allele(mark, mark->getAllele(samp->getAbone(mloc)), &options) + " ";
											a2 += Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(mloc)), &options) + " ";
										}
									}//end microsat
								}
							}
							a1 = a1.erase((a1.size() - 1), 1);
							a2 = a2.erase((a2.size() - 1), 1);

							str << a1 << "\n" << a2 << "\n";
						}
					}
				}
				str.close();
			}
		}
	}

}


int PHASEOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}

}
