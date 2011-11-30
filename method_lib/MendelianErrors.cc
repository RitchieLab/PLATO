/**********************************************************************************
*                       Mendelian Error Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all family genotypes and calculates mendelian errors.  Generates
* counts by family and marker.  Also generates a list of deletion candidates.
* Errors are checked from Father to Child and Mother to child for a total of up
* to 3 mendelian errors per trio (or more if larger family structure).
* Performs initial scan and removes bad markers based on threshold.  A secondary
* scan is then performed.
*
*File: MendelianErrors.cc
**********************************************************************************/


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include "Globals.h"
#include <string>
#include <list>
#include <map>
#include <algorithm>
#include <bitset>
#include "MendelianErrors.h"
//#include "Chrom.h"
#include "General.h"
#include "Helpers.h"

using namespace std;
namespace Methods{
string MendelianErrors::stepname = "mendelian-error";

void MendelianErrors::Tokenize(const string& str, vector<string>& tokens, const string& delimiter){
	string::size_type lastPos = str.find_first_not_of(delimiter, 0);
	string::size_type pos = str.find_first_of(delimiter, lastPos);

	while(string::npos != pos || string::npos != lastPos){
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		lastPos = str.find_first_not_of(delimiter, pos);
		pos = str.find_first_of(delimiter, lastPos);
	}
}

void MendelianErrors::setThreshold(string thresh){
	options.setUp(thresh);
}

void MendelianErrors::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::printLog("Families Passed:\t" + getString<int>(opts::_FAMILIES_WORKING_ - orig_num_families) + " (" +
        getString<float>(((float)(opts::_FAMILIES_WORKING_ - orig_num_families) / (float)opts::_FAMILIES_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_FAMILIES_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;
	opts::_FAMILIES_WORKING_ -= orig_num_families;

}

void MendelianErrors::calcThreshold(){
	/*int nummarkers = markers->getSize();
	int numfams = families->getSize();

	if(marker_thresh == -1){
		marker_thresh = (int)((float)nummarkers * error_rate * 0.33) + 1;
	}
	if(fam_thresh == -1){
		fam_thresh = (int)((float) numfams * error_rate * 0.33) + 1;
	}
	*/
}

void MendelianErrors::PrintSummary(){
	string fname1 = opts::_OUTPREFIX_ + "mendelian_error_family" + options.getOut() + ".txt";//+ getString<int>(order) + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	string fname2 = opts::_OUTPREFIX_ + "mendelian_error_individual" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	string fname3 = opts::_OUTPREFIX_ + "mendelian_error_marker" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname3 += "." + getString<int>(order);
	}
	ofstream myoutputf (fname1.c_str());
	ofstream myoutputi (fname2.c_str());
	ofstream myoutputm (fname3.c_str());
	if(!myoutputf){
		opts::printLog("Error opening " + fname1 + ".  Exiting!\n");
		//exit(1);
		throw MethodException("Error opening " + fname1 + ".  Exiting!\n");
	}
	opts::addFile("Family",stepname, fname1);
	if(!myoutputi){
		opts::printLog("Error opening " + fname2 + ".  Exiting!\n");
		//exit(1);
		throw MethodException("Error opening " + fname2 + ".  Exiting!\n");
	}
	opts::addFile("Sample",stepname,fname2);
	if(!myoutputm){
		opts::printLog("Error opening " + fname3 + ".  Exiting!\n");
		//exit(1);
		throw MethodException("Error opening " + fname3 + ".  Exiting!\n");
	}
	opts::addFile("Marker",stepname,fname3);
	myoutputi << "FamID\tIndID\tCenter\tSex\tAffection Satus\tPlate\tWell\tME_count_All";
	opts::addHeader(fname2, "ME_count_All");

	myoutputf << "FamID\tNumInds\tCenter\tME_count_All";
	opts::addHeader(fname1, "ME_count_All");

	myoutputm << "Chrom\trsID\tProbeID\tbploc\t";
	if((*markers)[0]->getDetailHeaders().size() > 0){
		myoutputm << (*markers)[0]->getDetailHeaders() << "\t";
	}
	myoutputm << "ME_count_All" << endl;
	opts::addHeader(fname3, "ME_count_All");

	int msize = markers->size();
	int fsize = families->size();
	int ssize = samples->size();

	vector<string> enzymes;
    if(opts::_ENZYMES_){
        int msize = markers->size();
        for(int i = 0; i < msize; i++){
            if((*markers)[i]->isEnabled() && !(*markers)[i]->isFlagged()){
				if(options.doChrom()){
					if(!options.checkChrom((*markers)[i]->getChrom())){
					    continue;
				    }
				    if(!options.checkBp((*markers)[i]->getBPLOC())){
					    continue;
				    }
				}

                vector<string>::iterator e_iter = find(enzymes.begin(), enzymes.end(), (*markers)[i]->getEnzyme());
                if(e_iter == enzymes.end()){
                    enzymes.push_back((*markers)[i]->getEnzyme());
                }
            }
        }
        if(enzymes.size() > 1){
            sort(enzymes.begin(), enzymes.end());
        }
        for(int i = 0; i < (int)enzymes.size(); i++){
            myoutputi << "\tME_count_" << enzymes[i];
			myoutputf << "\tME_count_" << enzymes[i];
			opts::addHeader(fname1, "ME_count_" + enzymes[i]);
			opts::addHeader(fname2, "ME_count_" + enzymes[i]);
        }
    }
	string sdetails = "";
	if(opts::_SAMPDESC_.length() > 0){
		sdetails = (*samples)[0]->getDetailHeaders();
	}
	myoutputi << "\t" << sdetails;
	myoutputi << endl;
	myoutputf << endl;

	for(int i = 0; i < msize; i++){
		if((*markers)[i]->isEnabled() && !(*markers)[i]->isFlagged()){
			if(options.doChrom()){
				if(!options.checkChrom((*markers)[i]->getChrom())){
				    continue;
			    }
			    if(!options.checkBp((*markers)[i]->getBPLOC())){
				    continue;
			    }
			}
			myoutputm << (*markers)[i]->toString() << "\t"
					  << merrors[i] << endl;
		}
	}

	for(int i = 0; i < fsize; i++){
		if((*families)[i]->isEnabled()){
			myoutputf << (*families)[i]->getFamID() << "\t"
				<< (*families)[i]->getSamples()->size() << "\t"
				<< (*families)[i]->getCenter() << "\t"
				<< ferrors[i];
			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < (int)enzymes.size(); e++){
					myoutputf << "\t" << fenzyme[i][enzymes[e]];
				}
			}
			myoutputf << endl;
		}
	}

	for(int i = 0; i < ssize; i++){
		if((*samples)[i]->isEnabled()){
			myoutputi << (*samples)[i]->getFamID() << "\t"
				<< (*samples)[i]->getInd() << "\t"
				<< (*samples)[i]->getFamily()->getCenter() << "\t";
			if((*samples)[i]->getSex()){
				myoutputi << "M\t";
			}
			else{
				myoutputi << "F\t";
			}
			if((*samples)[i]->getPheno() == 2){//getAffected()){
				myoutputi << "Y\t";
			}
			else if((*samples)[i]->getPheno() == 1){
				myoutputi << "N\t";
			}
			else{
				myoutputi << "U\t";
			}
			myoutputi << (*samples)[i]->getPlate() << "\t"
				<< (*samples)[i]->getWell() << "\t"
				<< serrors[i];
			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < (int)enzymes.size(); e++){
					myoutputi << "\t" << senzyme[i][enzymes[e]];
				}
			}
			if(opts::_SAMPDESC_.length() > 0){
				myoutputi << "\t" << (*samples)[i]->getDetails();
			}
			myoutputi << endl;
		}
	}

	if(myoutputm.is_open()){
		myoutputm.close();
	}
	if(myoutputi.is_open()){
		myoutputi.close();
	}
	if(myoutputf.is_open()){
		myoutputf.close();
	}

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}
}

void MendelianErrors::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

	calcThreshold();

	good_markers = Helpers::findValidMarkers(markers, &options);

	perform_evaluation(false);
	mid_process_output();
	if(options.zeroGenos()){
		zeroErrors();
	}
	filter_markers();
	//resetCounts();
	//perform_evaluation(false);
}

void MendelianErrors::zeroErrors(){

	int esize = error_map.size();
	for(int s = 0; s < esize; s++){
		if(error_map[s].size() > 0){
			for(int m = 0; m < (int)error_map[s].size(); m++){
				Marker* aloc = error_map[s][m];
				int mloc = aloc->getLoc();//(*markers)[aloc]->getLoc();
				(*samples)[s]->addAone(mloc, true);
				(*samples)[s]->addAtwo(mloc, false);
				if(aloc->isMicroSat()){//(*markers)[aloc]->isMicroSat()){
					(*samples)[s]->addAbone(mloc, -1);
					(*samples)[s]->addAbtwo(mloc, -1);
				}
			}
		}
	}
}

void MendelianErrors::perform_evaluation(bool output){
	ofstream myoutput;
	//ofstream qsoutput;
	ofstream erroroutput;
	ofstream level2output;
	if(output){
		ostringstream rstr;
		rstr << rank;
		string filename = opts::_OUTPREFIX_ + "ME_mid_possible_deletion_" + getString<int>(order) + ".txt";
		myoutput.open(filename.c_str(), ios::out | ios::app);
		if(!myoutput){
			opts::printLog("Error opening ME_mid_possible_deletion.txt.  Exiting!\n");
			//exit(1);
			throw MethodException("Error opening ME_mid_possible_deletion.txt.  Exiting!\n");
		}
		myoutput << "Chrom\trsID\tProbeID\tbploc\tFamID\tFather_Geno\tMother_Geno\tChild_Geno\tPossible_Deletion?" << endl;
		//qsoutput.open("ME_quality_scores.txt", ios::out | ios::app);
		//qsoutput << "Chrom\trsID\tProbeID\tbploc\tEnzyme\tFamID\tIndID\tSex\tPlate\tWell\tGenotype\tQuality_Score" << endl;
	}
	string fname = opts::_OUTPREFIX_ + "mendelian_error_errors" + options.getOut() + ".txt";// + getString<int>(order) + ".txt";
	while(!overwrite && Helpers::fileExists(fname)){
		fname += "." + getString<int>(order);
	}
	errors_file_name = fname;
	erroroutput.open(fname.c_str(), ios::out);
	if(!erroroutput){
		opts::printLog("Error opening " + fname + ".  Exiting!\n");
		//exit(1);
		throw MethodException("Error opening " + fname + ".  Exiting!\n");

	}
	string fname2 = opts::_OUTPREFIX_ + "mendelian_error_level2" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	while(!overwrite && Helpers::fileExists(fname2)){
		fname2 += "." + getString<int>(order);
	}
	level2output.open(fname2.c_str(), ios::out);
	level2_file_name = fname2;
	if(!level2output){
		opts::printLog("Error opening " + fname2 + ".  Exiting!\n");
		//exit(1);
		throw MethodException("Error opening " + fname2 + ".  Exiting!\n");
	}
	level2output << "(*** denotes probable location of error)\n";
	erroroutput << "Chrom\trsID\tProbeID\tbploc\tFamID\tIndID\tSex\tGenotype\tP_ID\t\tP_Sex\tP_Genotype" << endl;

	int fsize = families->size();
	int msize = markers->size();
	int ssize = samples->size();
	ferrors.resize(fsize);
	merrors.resize(msize);
	serrors.resize(ssize);
	error_map.resize(ssize);
	if(opts::_ENZYMES_){
		senzyme.resize(ssize);
		fenzyme.resize(fsize);
	}
	//for(int s = 0; s < ssize; s++){
	//	Sample* samp = (*samples)[s];
	//	if(samp->isEnabled()){
	//		samp->resizeDeletionCat(msize, '0');
	//	}
	//}

//	int prev_base = 0;
//	int prev_chrom = -1;

	if(good_markers.size() == 0){
		good_markers = Helpers::findValidMarkers(markers, &options);
	}

	msize = good_markers.size();

	for(int m = 0; m < msize; m++){
		vector<Family*> family_inc;

		Marker* mark = good_markers[m];//(*markers)[m];
		if(mark->isEnabled() && mark->getChrom() <= opts::_CHRX_){
/*			if(options.doChrom()){
				if(!options.checkChrom((*markers)[m]->getChrom())){
				    continue;
			    }
			    if(!options.checkBp((*markers)[m]->getBPLOC())){
				    continue;
			    }
			}
            if(options.doBpSpace()){
	            if(prev_base == 0){
		            prev_base = (*markers)[m]->getBPLOC();
		            prev_chrom = (*markers)[m]->getChrom();
		        }
            	else{
            		if((*markers)[m]->getChrom() == prev_chrom && (((*markers)[m]->getBPLOC() - prev_base) < options.getBpSpace())){
						(*markers)[m]->setFlag(true);
            			continue;
            		}
            		prev_base = (*markers)[m]->getBPLOC();
            		prev_chrom = (*markers)[m]->getChrom();
            	}
            }
*/
			int mloc = mark->getLoc();

			for(int s = 0; s < ssize; s++){
				if((*samples)[s]->isEnabled()){
					bool inc = false;
					bool dadinc = false;
					bool childinc = false;
					bool mominc = false;
					if((*samples)[s]->getDadID() != "0" && (*samples)[s]->getDad() != NULL){
						if(((*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc) && (*samples)[s]->getAmissing(mloc) &&
							!mark->isMicroSat()) ||
							((*samples)[s]->getDad()->getAone(mloc) && (*samples)[s]->getDad()->getAtwo(mloc) && (*samples)[s]->getDad()->getAmissing(mloc) &&
							 !mark->isMicroSat()) ||
							(mark->isMicroSat() && (*samples)[s]->getAbone(mloc) == -1) ||
							(mark->isMicroSat() && (*samples)[s]->getDad()->getAbone(mloc) == -1)){
							//continue;
						}
						else{
							if(mark->getChrom() != opts::_CHRX_ || (!(*samples)[s]->getSex() && mark->getChrom() == opts::_CHRX_)){
								if((!mark->isMicroSat() &&
							   		(*samples)[s]->getAone(mloc) != (*samples)[s]->getDad()->getAone(mloc) &&
							   		(*samples)[s]->getAtwo(mloc) != (*samples)[s]->getDad()->getAone(mloc) &&
							   		(*samples)[s]->getAone(mloc) != (*samples)[s]->getDad()->getAtwo(mloc) &&
							   		(*samples)[s]->getAtwo(mloc) != (*samples)[s]->getDad()->getAtwo(mloc)) ||
							   		(mark->isMicroSat() &&
							   		(*samples)[s]->getAbone(mloc) != (*samples)[s]->getDad()->getAbone(mloc) &&
							   		(*samples)[s]->getAbtwo(mloc) != (*samples)[s]->getDad()->getAbone(mloc) &&
							   		(*samples)[s]->getAbone(mloc) != (*samples)[s]->getDad()->getAbtwo(mloc) &&
							   		(*samples)[s]->getAbtwo(mloc) != (*samples)[s]->getDad()->getAbtwo(mloc))){
									vector<Family*>::iterator found = find(family_inc.begin(), family_inc.end(), (*samples)[s]->getFamily());
									if(!inc && found == family_inc.end()){
										ferrors[(*samples)[s]->getFamily()->getLoc()]++;
										if(opts::_ENZYMES_){
											fenzyme[(*samples)[s]->getFamily()->getLoc()][mark->getEnzyme()]++;
										}
										merrors[m]++;
										family_inc.push_back((*samples)[s]->getFamily());
										inc = true;
									}
									if(!childinc){
										serrors[s]++;
										if(opts::_ENZYMES_){
											senzyme[s][mark->getEnzyme()]++;
										}
										childinc = true;
									}
									if(!dadinc){
										serrors[(*samples)[s]->getDad()->getLoc()]++;
										if(opts::_ENZYMES_){
											senzyme[(*samples)[s]->getDad()->getLoc()][mark->getEnzyme()]++;
										}
										dadinc = true;
									}
									//cout << "we1\n";
									write_error(erroroutput, mark, (*samples)[s], (*samples)[s]->getDad(), mloc);
								}
							}
						}//end else
					}

					if((*samples)[s]->getMomID() != "0" && (*samples)[s]->getMom() != NULL){
						if(((!mark->isMicroSat()) &&
							(((*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc) && (*samples)[s]->getAmissing(mloc)) ||
							((*samples)[s]->getMom()->getAone(mloc) && (*samples)[s]->getMom()->getAtwo(mloc) && (*samples)[s]->getMom()->getAmissing(mloc)))) ||
							((mark->isMicroSat()) && (((*samples)[s]->getAbone(mloc) == -1 || (*samples)[s]->getAbone(mloc) == -1))) || ((mark->isMicroSat()) && (((*samples)[s]->getMom()->getAbone(mloc) == -1)))){
							//continue;
						}
						else{
							if((!mark->isMicroSat() &&
						   		(*samples)[s]->getAone(mloc) != (*samples)[s]->getMom()->getAone(mloc) &&
						   		(*samples)[s]->getAtwo(mloc) != (*samples)[s]->getMom()->getAone(mloc) &&
						   		(*samples)[s]->getAone(mloc) != (*samples)[s]->getMom()->getAtwo(mloc) &&
						   		(*samples)[s]->getAtwo(mloc) != (*samples)[s]->getMom()->getAtwo(mloc)) ||
						   		(mark->isMicroSat() &&
						   		(*samples)[s]->getAbone(mloc) != (*samples)[s]->getMom()->getAbone(mloc) &&
						   		(*samples)[s]->getAbtwo(mloc) != (*samples)[s]->getMom()->getAbone(mloc) &&
						   		(*samples)[s]->getAbone(mloc) != (*samples)[s]->getMom()->getAbtwo(mloc) &&
						   		(*samples)[s]->getAbtwo(mloc) != (*samples)[s]->getMom()->getAbtwo(mloc))){
								vector<Family*>::iterator found = find(family_inc.begin(), family_inc.end(), (*samples)[s]->getFamily());
								if(!inc && found == family_inc.end()){
									ferrors[(*samples)[s]->getFamily()->getLoc()]++;
									if(opts::_ENZYMES_){
										fenzyme[(*samples)[s]->getFamily()->getLoc()][mark->getEnzyme()]++;
									}
									merrors[m]++;
									family_inc.push_back((*samples)[s]->getFamily());
									inc = true;
								}
								if(!childinc){
									serrors[s]++;
									if(opts::_ENZYMES_){
										senzyme[s][mark->getEnzyme()]++;
									}
									childinc = true;
								}
								if(!mominc){
									serrors[(*samples)[s]->getMom()->getLoc()]++;
									if(opts::_ENZYMES_){
										senzyme[(*samples)[s]->getMom()->getLoc()][mark->getEnzyme()]++;
									}
									mominc = true;
								}
								//cout << "we2:" << (*markers)[m]->toString() << "\t" << (*samples)[s]->toString() << "\t" << (*markers)[m]->getLoc() << ":" << mloc << "\n";
								write_error(erroroutput, mark, (*samples)[s], (*samples)[s]->getMom(), mloc);
							}
						}//end else
					}
					if((*samples)[s]->getDadID() != "0" && (*samples)[s]->getDad() != NULL && (*samples)[s]->getMomID() != "0" && (*samples)[s]->getMom() != NULL){
						if(((!mark->isMicroSat()) &&
							(((*samples)[s]->getAone(mloc) && (*samples)[s]->getAtwo(mloc) && (*samples)[s]->getAmissing(mloc)) ||
							((*samples)[s]->getDad()->getAone(mloc) && (*samples)[s]->getDad()->getAtwo(mloc) && (*samples)[s]->getDad()->getAmissing(mloc)) ||
							((*samples)[s]->getMom()->getAone(mloc) && (*samples)[s]->getMom()->getAtwo(mloc) && (*samples)[s]->getMom()->getAmissing(mloc)))) ||
						   ((mark->isMicroSat()) &&
							(((*samples)[s]->getAbone(mloc) == -1 || (*samples)[s]->getDad()->getAbone(mloc) == -1 || (*samples)[s]->getMom()->getAbone(mloc) == -1)))){
							//continue;
						}
						else{
							if(mark->getChrom() != opts::_CHRX_ || (!(*samples)[s]->getSex() && mark->getChrom() == opts::_CHRX_)){
								if((!mark->isMicroSat() &&
									(*samples)[s]->getDad()->getAone(mloc) == (*samples)[s]->getDad()->getAtwo(mloc) &&
									(*samples)[s]->getMom()->getAone(mloc) == (*samples)[s]->getMom()->getAtwo(mloc) &&
									(*samples)[s]->getMom()->getAone(mloc) == (*samples)[s]->getDad()->getAone(mloc) &&
									((*samples)[s]->getAone(mloc) != (*samples)[s]->getMom()->getAone(mloc) ||
									 (*samples)[s]->getAtwo(mloc) != (*samples)[s]->getMom()->getAone(mloc))) ||
									(mark->isMicroSat() && ((
									(*samples)[s]->getDad()->getAbone(mloc) == (*samples)[s]->getDad()->getAbtwo(mloc) &&
									(*samples)[s]->getMom()->getAbone(mloc) == (*samples)[s]->getMom()->getAbtwo(mloc) &&
									(*samples)[s]->getMom()->getAbone(mloc) == (*samples)[s]->getDad()->getAbone(mloc) &&
									((*samples)[s]->getAbone(mloc) != (*samples)[s]->getMom()->getAbone(mloc) ||
									 (*samples)[s]->getAbtwo(mloc) != (*samples)[s]->getMom()->getAbone(mloc))) ||
										(((*samples)[s]->getAbone(mloc) != (*samples)[s]->getMom()->getAbone(mloc) &&
										  (*samples)[s]->getAbone(mloc) != (*samples)[s]->getMom()->getAbtwo(mloc) &&
										  (*samples)[s]->getAbone(mloc) != (*samples)[s]->getDad()->getAbone(mloc) &&
										  (*samples)[s]->getAbone(mloc) != (*samples)[s]->getDad()->getAbtwo(mloc)
										 ) ||
										 ((*samples)[s]->getAbtwo(mloc) != (*samples)[s]->getMom()->getAbone(mloc) &&
										  (*samples)[s]->getAbtwo(mloc) != (*samples)[s]->getMom()->getAbtwo(mloc) &&
										  (*samples)[s]->getAbtwo(mloc) != (*samples)[s]->getDad()->getAbone(mloc) &&
										  (*samples)[s]->getAbtwo(mloc) != (*samples)[s]->getDad()->getAbtwo(mloc)
										 )
										  )))){

									vector<Family*>::iterator found = find(family_inc.begin(), family_inc.end(), (*samples)[s]->getFamily());
							//	cout << "Possible increase at " << (*markers)[m]->getChrom() << "\t" << (*markers)[m]->getProbeID() << endl;
									if(!inc && found == family_inc.end()){
										ferrors[(*samples)[s]->getFamily()->getLoc()]++;
										if(opts::_ENZYMES_){
											fenzyme[(*samples)[s]->getFamily()->getLoc()][mark->getEnzyme()]++;
										}
										merrors[m]++;
										family_inc.push_back((*samples)[s]->getFamily());
										inc = true;
									}
									if(!childinc){
										serrors[s]++;
										if(opts::_ENZYMES_){
											senzyme[s][mark->getEnzyme()]++;
										}
										childinc = true;
									}
									if(!mominc){
										serrors[(*samples)[s]->getMom()->getLoc()]++;
										if(opts::_ENZYMES_){
											senzyme[(*samples)[s]->getMom()->getLoc()][mark->getEnzyme()]++;
										}
										mominc = true;
									}
									if(!dadinc){
										serrors[(*samples)[s]->getDad()->getLoc()]++;
										if(opts::_ENZYMES_){
											senzyme[(*samples)[s]->getDad()->getLoc()][mark->getEnzyme()]++;
										}
										dadinc = true;
									}
									//cout << "we3\n";
									write_error(erroroutput, mark, (*samples)[s], (*samples)[s]->getDad(), mloc);
									//cout << "we4\n";
									write_error(erroroutput, mark, (*samples)[s], (*samples)[s]->getMom(), mloc);

								}
							}
						}//end else
					}
					if(childinc){
						(*samples)[s]->getFamily()->setMeError(true);
						error_map[s].push_back(mark);
					}
					if(dadinc){
						(*samples)[s]->getFamily()->setMeError(true);
						error_map[(*samples)[s]->getDad()->getLoc()].push_back(mark);
					}
					if(mominc){
						(*samples)[s]->getFamily()->setMeError(true);
						error_map[(*samples)[s]->getMom()->getLoc()].push_back(mark);
					}
				}
			}// end level 1

			//check sibships
		//	for(int f = 0; f < fsize; f++){
		//		Family* fam = (*families)[f];
		//		if(fam->isEnabled()){
		//			vector<Sample*>* founders = fam->getFounders();
		//		}
		//	}
			//end sibships

			//level 2
			for(int f = 0; f < fsize; f++){
				Family* fam = (*families)[f];
				if(fam->isEnabled() && !fam->hasMeError() && fam->getSamples()->size() > 1){
					bool realgeno = getPossibleGenos(fam, mark);
					if(!realgeno){
						continue;
					}
//					vector<Family*>* nuclears = generateNuclearFams(fam);
//					for(int n = 0; n < nuclears->size(); n++){
//					}
					int removed = 0;
					do{
						removed = 0;
						removed += checkTypes(fam, mark);
					//	level2output << "After checkTypes: " << removed << endl; //
					//	printError(level2output, fam, mark); //
						genotypeElimination2(fam, mark, removed);
					//	level2output << "removed: " << removed << endl; //
					//	printError(level2output, fam, mark); //
					//	level2output << removed << "\t" << fam->toString() << "\t" << mark->toString() << endl; //
						//printError(level2output, fam, mark); //
					}while(removed > 0);
					if(checkErrors(fam, mark)){
						printError(level2output, fam, mark);
						if(options.zeroL2Genos()){
							vector<Sample*>* fsamps = fam->getSamples();
							for(int fs = 0; fs < (int)fsamps->size(); fs++){
								Sample* mysamp = (*fsamps)[fs];
								if(options.zeroL2FamGenos()){
									if(!mark->isMicroSat()){
										mysamp->addAone(mark->getLoc(), true);
										mysamp->addAtwo(mark->getLoc(), true);
										mysamp->addAmissing(mark->getLoc(), true);
									}
									else{
										mysamp->addAbone(mark->getLoc(), -1);
										mysamp->addAbtwo(mark->getLoc(), -1);
									}
								}
								else if(mysamp->IsMEerror()){
									if(!mark->isMicroSat()){
										mysamp->addAone(mark->getLoc(), true);
										mysamp->addAtwo(mark->getLoc(), true);
										mysamp->addAmissing(mark->getLoc(), true);
									}
									else{
										mysamp->addAbone(mark->getLoc(), -1);
										mysamp->addAbtwo(mark->getLoc(), -1);
									}
								}
							}
						}
					}
				}
				fam->setMeError(false);
			}
			//end level 2
		}
	}

	if(myoutput && myoutput.is_open()){
		myoutput.close();
	}
}

void MendelianErrors::genotypeElimination2(Family* fam, Marker* mark, int &removed){
	vector<Sample*>* samps = fam->getSamples();
	for(int s = 0; s < (int)samps->size(); s++){
		Sample* samp = (*samps)[s];
		if(samp->isEnabled() && !samp->getSex() && samp->getChildren()->size() > 0){
			vector<Sample*>* children = samp->getChildren();
			int enableddads = 0;
			for(int ch = 0; ch < (int)children->size(); ch++){
				Sample* child = (*children)[ch];
				Sample* dad = child->getDad();
				Sample* mom = samp;
				if(dad && dad->isEnabled() && mom->isEnabled()){
					enableddads++;
					if(!mark->isMicroSat()){
						vector<bool> da1 = dad->getAonePossible();
						vector<bool> da2 = dad->getAtwoPossible();
						vector<bool> da3 = dad->getAmissingPossible();
						vector<bool> ma1 = mom->getAonePossible();
						vector<bool> ma2 = mom->getAtwoPossible();
						vector<bool> ma3 = mom->getAmissingPossible();
						vector<bool> ca1 = child->getAonePossible();
						vector<bool> ca2 = child->getAtwoPossible();
						vector<bool> ca3 = child->getAmissingPossible();

						vector<string> alleles = mark->getAlleles();
						for(int mg = 0; mg < (int)ma1.size(); mg++){
							for(int dg = 0; dg < (int)da1.size(); dg++){
								vector<vector<bool> > zygotes;
								vector<bool> geno1(3,false);
								vector<bool> geno2(3,false);
								vector<bool> geno3(3,false);
								vector<bool> geno4(3,false);
								geno1[0] = da1[dg];
								geno1[1] = ma1[mg];
								geno2[0] = da1[dg];
								geno2[1] = ma2[mg];
								geno3[0] = da2[dg];
								geno3[1] = ma1[mg];
								geno4[0] = da2[dg];
								geno4[1] = ma2[mg];
								vector<vector<bool> >::iterator found = find(zygotes.begin(), zygotes.end(), geno1);
								if(found == zygotes.end()){
									zygotes.push_back(geno1);
								}
								found = find(zygotes.begin(), zygotes.end(), geno2);
								if(found == zygotes.end()){
									zygotes.push_back(geno2);
								}
								found = find(zygotes.begin(), zygotes.end(), geno3);
								if(found == zygotes.end()){
									zygotes.push_back(geno3);
								}
								found = find(zygotes.begin(), zygotes.end(), geno4);
								if(found == zygotes.end()){
									zygotes.push_back(geno4);
								}
							//	cout << dad->toString() << "(" << alleles[da1[dg]] << "/" << alleles[da2[dg]] << ")\t" << mom->toString() << "(" << alleles[ma1[mg]] << "/" << alleles[ma2[mg]] << ")" << endl;
							//	cout << "ZYGOTES:\n";
							//	for(int z = 0; z < zygotes.size(); z++){
							//		cout << alleles[zygotes[z][0]] << "/" << alleles[zygotes[z][1]] << " ";
							//	}
							//	cout << endl;
							//	cout << child->toString() << " genos found: ";

								for(int c = 0; c < (int)ca1.size(); c++){
									vector<bool> cgeno(3,false);
									cgeno[0] = ca1[c];
									cgeno[1] = ca2[c];
									found = find(zygotes.begin(), zygotes.end(), cgeno);
									if(found != zygotes.end()){
										child->addMEsaved(c);
										mom->addMEsaved(mg);
										dad->addMEsaved(dg);
							//			if(!cgeno[0]){
							//				cout << alleles[0];
							//			}
							//			else{
							//				cout << alleles[1];
							//			}
							//			if(!cgeno[1]){
							//				cout << "/" << alleles[0] << " ";
							//			}
							//			else{
							//				cout << "/" << alleles[1] << " ";
							//			}
									}
								}//end foreach childgeno
							//	cout << endl << endl;
							}//end foreach dadgeno
						}//end foreach momgeno
					}
					else{
						vector<int> da1 = dad->getAbonePossible();
						vector<int> da2 = dad->getAbtwoPossible();
						vector<int> ma1 = mom->getAbonePossible();
						vector<int> ma2 = mom->getAbtwoPossible();
						vector<int> ca1 = child->getAbonePossible();
						vector<int> ca2 = child->getAbtwoPossible();

						vector<string> alleles = mark->getAlleles();
						for(int mg = 0; mg < (int)ma1.size(); mg++){
							for(int dg = 0; dg < (int)da1.size(); dg++){
								vector<vector<int> > zygotes;
								vector<int> geno1(2,-1);
								vector<int> geno2(2,-1);
								vector<int> geno3(2,-1);
								vector<int> geno4(2,-1);
								geno1[0] = da1[dg];
								geno1[1] = ma1[mg];
								geno2[0] = da1[dg];
								geno2[1] = ma2[mg];
								geno3[0] = da2[dg];
								geno3[1] = ma1[mg];
								geno4[0] = da2[dg];
								geno4[1] = ma2[mg];
								vector<vector<int> >::iterator found = find(zygotes.begin(), zygotes.end(), geno1);
								if(found == zygotes.end()){
									zygotes.push_back(geno1);
								}
								found = find(zygotes.begin(), zygotes.end(), geno2);
								if(found == zygotes.end()){
									zygotes.push_back(geno2);
								}
								found = find(zygotes.begin(), zygotes.end(), geno3);
								if(found == zygotes.end()){
									zygotes.push_back(geno3);
								}
								found = find(zygotes.begin(), zygotes.end(), geno4);
								if(found == zygotes.end()){
									zygotes.push_back(geno4);
								}
						//		cout << dad->toString() << "(" << alleles[da1[dg]] << "/" << alleles[da2[dg]] << ")\t" << mom->toString() << "(" << alleles[ma1[mg]] << "/" << alleles[ma2[mg]] << ")" << endl;
						//		cout << "ZYGOTES:\n";
						//		for(int z = 0; z < zygotes.size(); z++){
						//			cout << alleles[zygotes[z][0]] << "/" << alleles[zygotes[z][1]] << " ";
						//		}
						//		cout << endl;
						//		cout << child->toString() << " genos found: ";
								for(int c = 0; c < (int)ca1.size(); c++){
									vector<int> cgeno(2,-1);
									cgeno[0] = ca1[c];
									cgeno[1] = ca2[c];
									//if(find_match(zygotes, cgeno)){
									found = find(zygotes.begin(), zygotes.end(), cgeno);
									if(found != zygotes.end()){
										child->addMEsaved(c);
										mom->addMEsaved(mg);
										dad->addMEsaved(dg);
										//cout << alleles[cgeno[0]] << "/" << alleles[cgeno[1]] << " ";
									}
									else{
									}
								}//end foreach childgeno
								//cout << endl << endl;
							}//end foreach dadgeno
						}//end foreach momgeno
					}//end ms else
				}//end mom & dad enabled
			}//end foreach children
			//foreach child, remove child, dad
			vector<Sample*> dads;
			for(int ch = 0; ch < (int)children->size(); ch++){
				Sample* child = (*children)[ch];
				Sample* dad = child->getDad();
				if(child->isEnabled() && dad && dad->isEnabled()){
					vector<Sample*>::iterator found = find(dads.begin(), dads.end(), dad);
					if(found == dads.end()){
						dads.push_back(dad);
					}
					if(!mark->isMicroSat()){
						vector<bool> temp1;
						vector<bool> temp2;
						vector<bool> temp3;
						removed += child->getAonePossible().size() - child->getMEsavedCount();
						vector<bool> ca1 = child->getAonePossible();
						vector<bool> ca2 = child->getAtwoPossible();
						vector<bool> ca3 = child->getAmissingPossible();
						for(int i = 0; i < child->getMEsavedCount(); i++){
							temp1.push_back(ca1[child->getMEsaved(i)]);
							temp2.push_back(ca2[child->getMEsaved(i)]);
							temp3.push_back(ca3[child->getMEsaved(i)]);
						}
						child->setAonePossible(temp1);
						child->setAtwoPossible(temp2);
						child->setAmissingPossible(temp3);

						child->clearMEsaved();
					}
					else{
						vector<int> temp1;
						vector<int> temp2;
						removed += child->getAbonePossible().size() - child->getMEsavedCount();
						vector<int> ca1 = child->getAbonePossible();
						vector<int> ca2 = child->getAbtwoPossible();
						for(int i = 0; i < child->getMEsavedCount(); i++){
							temp1.push_back(ca1[child->getMEsaved(i)]);
							temp2.push_back(ca2[child->getMEsaved(i)]);
						}
						child->setAbonePossible(temp1);
						child->setAbtwoPossible(temp2);

						child->clearMEsaved();
					}
				}
			}
			for(int dd = 0; dd < (int)dads.size(); dd++){
				Sample* dad = dads[dd];
				if(!mark->isMicroSat()){
					vector<bool> temp1;
					vector<bool> temp2;
					vector<bool> temp3;
					vector<bool> ca1 = dad->getAonePossible();
					vector<bool> ca2 = dad->getAtwoPossible();
					vector<bool> ca3 = dad->getAmissingPossible();
					removed += dad->getAonePossible().size() - dad->getMEsavedCount();
					temp1.clear();
					temp2.clear();
					temp3.clear();
					//cout << "Dad ME saved count: " << dad->getMEsavedCount() << endl;
					for(int i = 0; i < dad->getMEsavedCount(); i++){
						temp1.push_back(ca1[dad->getMEsaved(i)]);
						temp2.push_back(ca2[dad->getMEsaved(i)]);
						temp3.push_back(ca3[dad->getMEsaved(i)]);
					}
					dad->setAonePossible(temp1);
					dad->setAtwoPossible(temp2);
					dad->setAmissingPossible(temp3);
				}
				else{
					vector<int> temp1;
					vector<int> temp2;
					removed += dad->getAbonePossible().size() - dad->getMEsavedCount();
					vector<int> ca1 = dad->getAbonePossible();
					vector<int> ca2 = dad->getAbtwoPossible();
					temp1.clear();
					temp2.clear();
					for(int i = 0; i < dad->getMEsavedCount(); i++){
						temp1.push_back(ca1[dad->getMEsaved(i)]);
						temp2.push_back(ca2[dad->getMEsaved(i)]);
					}
					dad->setAbonePossible(temp1);
					dad->setAbtwoPossible(temp2);
				}
				dad->clearMEsaved();
			}
			//remove mom stuff
			if(enableddads > 0){
				if(!mark->isMicroSat()){
					vector<bool> temp1;
					vector<bool> temp2;
					vector<bool> temp3;
					removed += samp->getAonePossible().size() - samp->getMEsavedCount();
					vector<bool> ma1 = samp->getAonePossible();
					vector<bool> ma2 = samp->getAtwoPossible();
					vector<bool> ma3 = samp->getAmissingPossible();
					for(int i = 0; i < samp->getMEsavedCount(); i++){
						temp1.push_back(ma1[samp->getMEsaved(i)]);
						temp2.push_back(ma2[samp->getMEsaved(i)]);
						temp3.push_back(ma3[samp->getMEsaved(i)]);
					}
					samp->setAonePossible(temp1);
					samp->setAtwoPossible(temp2);
					samp->setAmissingPossible(temp3);
					samp->clearMEsaved();
				}
				else{
					vector<int> temp1;
					vector<int> temp2;
					removed += samp->getAbonePossible().size() - samp->getMEsavedCount();
					vector<int> ma1 = samp->getAbonePossible();
					vector<int> ma2 = samp->getAbtwoPossible();
					for(int i = 0; i < samp->getMEsavedCount(); i++){
						temp1.push_back(ma1[samp->getMEsaved(i)]);
						temp2.push_back(ma2[samp->getMEsaved(i)]);
					}
					samp->setAbonePossible(temp1);
					samp->setAbtwoPossible(temp2);
					samp->clearMEsaved();
				}
			}
		}//end samp enabled & have mom and dad
	}//end foreach parent
}

bool MendelianErrors::find_match(vector<vector<int> > zygotes, vector<int> geno){
	bool found = false;
	for(int i = 0; i < (int)zygotes.size(); i++){
		if(zygotes[i][0] == geno[0] && zygotes[i][1] == geno[1]){
			found = true;
			break;
		}
	}
	return found;
}

void MendelianErrors::printError(ofstream &level2output, Family* fam, Marker* mark){
	level2output << "Marker: " << mark->toString() << endl;
	level2output << "---------------------------\n";
	vector<Sample*>* samps = fam->getSamples();
	for(int i = 0; i < (int)samps->size(); i++){
		Sample* samp = (*samps)[i];
		if(samp->isEnabled()){
			if(samp->IsMEerror()){
				level2output << "***";
			}
			level2output << samp->toString() << "\tPossible Genotypes: " << samp->possibleToString(mark->getAlleles()) << "\tActual Genotype: " << samp->genoToString(mark->isMicroSat(), mark->getLoc(), mark->getAlleles()) << endl;
		}
	}
	level2output << "\n\n";
}

bool MendelianErrors::checkErrors(Family* fam, Marker* mark){
	bool bad = true;
	int enabled = 0;
	int zero = 0;
	vector<Sample*>* samps = fam->getSamples();
	for(int i = 0; i < (int)samps->size(); i++){
		Sample* samp = (*samps)[i];
		samp->setME2error(false);
		for(int m = 0; m < (int)error_map[samp->getLoc()].size(); m++){
			if(error_map[samp->getLoc()][m] == mark){//(*markers)[error_map[samp->getLoc()][m]] == mark){
				return false;
			}
		}
	}
	for(int i = 0; i < (int)samps->size(); i++){
		Sample* samp = (*samps)[i];
		if(samp->isEnabled()){
			samp->setME2error(false);

			enabled++;
			if(!mark->isMicroSat()){
				if(samp->getAonePossible().size() == 0){
					bad = true;
					samp->setME2error(true);
					break;
				}
				bool a1 = samp->getAone(mark->getLoc());
				bool a2 = samp->getAtwo(mark->getLoc());
				bool a3 = samp->getAmissing(mark->getLoc());
				if(a1 && a2 && a3){
					zero++;
					continue;
				}
				bool het = false;
				if(!a1 && a2){
					het = true;
				}
				vector<bool> sa1 = samp->getAonePossible();
				vector<bool> sa2 = samp->getAtwoPossible();
				vector<bool> sa3 = samp->getAmissingPossible();
				for(int a = 0; a < (int)sa1.size(); a++){
					if(het && sa1[a] && sa2[a] && sa3[a]){
						bad = false;
						break;
					}
					if(a1 && sa1[a] && a2 && sa2[a] && !a3 && !sa3[a]){
						bad = false;
						break;
					}
					if(!a1 && !sa1[a] && a2 && sa2[a]){
						bad = false;
						break;
					}
					if(a1 && sa1[a] && !a2 && !sa2[a]){
						bad = false;
						break;
					}
					if(!a1 && !sa1[a] && !a2 && !sa2[a]){
						bad = false;
						break;
					}
					if(het && sa1[a] && !sa2[a] && !a1 && a2){
						bad = false;
						break;
					}
					bad = true;
				}
				if(bad){
					samp->setME2error(true);
					break;
				}
			}
			else{
				if(samp->getAbonePossible().size() == 0){
					bad = true;
					samp->setME2error(true);
					break;
				}
				int a1 = samp->getAbone(mark->getLoc());
				int a2 = samp->getAbtwo(mark->getLoc());
				if(a1 == -1 && a2 == -1){
					zero++;
					continue;
				}
				bool het = false;
				if(a1 != a2){
					het = true;
				}
				vector<int> sa1 = samp->getAbonePossible();
				vector<int> sa2 = samp->getAbtwoPossible();
				for(int a = 0; a < (int)sa1.size(); a++){
					if(sa1[a] == a1 && sa2[a] == a2){
						bad = false;
						break;
					}
					if(het && sa1[a] == a2 && sa2[a] == a1){
						bad = false;
						break;
					}
				}
				if(bad){
					samp->setME2error(true);
					break;
				}
			}
		}
	}
	if(zero == enabled){
		bad = false;
	}
	return bad;
}

//performs checks if any parent or child is untyped for the original genotype
int MendelianErrors::checkTypes(Family* fam, Marker* mark){
	//cout << "Starting checkTypes\n";
	int removed = 0;
	int mloc = mark->getLoc();
	if(mark->isMicroSat()){
		vector<Sample*>* samps = fam->getSamples();
		for(int s = 0; s < (int)samps->size(); s++){
			Sample* samp = (*samps)[s];
			if(samp->isEnabled()){
				samp->clearMEsaved();
				bool child = false;
				Sample* dad = samp->getDad();
				Sample* mom = samp->getMom();
				//untyped child vs typed parents
				if((dad || mom) && (samp->getAbone(mloc) == -1 && samp->getAbtwo(mloc) == -1)){
					if(dad != NULL){
						child = true;
						if(dad->getAbone(mloc) != -1 && dad->getAbtwo(mloc) != -1){
							vector<int> da1 = dad->getAbonePossible();
							vector<int> da2 = dad->getAbtwoPossible();
							for(int ca = 0; ca < (int)samp->getAbonePossible().size(); ca++){
								vector<int>::iterator found = find(da1.begin(), da1.end(), samp->getAbonePossibleVal(ca));
								if(found == da1.end()){
									found = find(da2.begin(), da2.end(), samp->getAbonePossibleVal(ca));
									if(found == da2.end()){
										samp->remMSPossible(ca);
										removed++;
										ca--;
									}
								}

							}
						}
					}
					if(mom != NULL){
						child = true;
						if(mom->getAbone(mloc) != -1 && mom->getAbtwo(mloc) != -1){
							vector<int> ma1 = mom->getAbonePossible();
							vector<int> ma2 = mom->getAbtwoPossible();
							for(int ca = 0; ca < (int)samp->getAbtwoPossible().size(); ca++){
								vector<int>::iterator found = find(ma1.begin(), ma1.end(), samp->getAbtwoPossibleVal(ca));
								if(found == ma1.end()){
									found = find(ma2.begin(), ma2.end(), samp->getAbtwoPossibleVal(ca));
									if(found == ma2.end()){
										samp->remMSPossible(ca);
										removed++;
										ca--;
									}
								}
							}
						}
					}
				}
				//untyped parent vs typed children
				if(!dad && !mom && !child){
					if(samp->getAbone(mloc) == -1 && samp->getAbtwo(mloc) == -1){
						vector<Sample*>* children = samp->getChildren();
						for(int c = 0; c < (int)children->size(); c++){
							Sample* child = (*children)[c];
							if(child->getAbone(mloc) != -1 && child->getAbtwo(mloc) != -1){
								vector<int> ca1 = child->getAbonePossible();
								vector<int> ca2 = child->getAbtwoPossible();
								for(int sa = 0; sa < (int)samp->getAbonePossible().size(); sa++){
									if(samp->getSex()){
										vector<int>::iterator found = find(ca1.begin(), ca1.end(), samp->getAbonePossibleVal(sa));
										if(found == ca1.end()){
											found = find(ca1.begin(), ca1.end(), samp->getAbtwoPossibleVal(sa));
											if(found == ca1.end()){
												samp->remMSPossible(sa);
												sa--;
												removed++;
											}
										}
									}
									else{
										vector<int>::iterator found = find(ca2.begin(), ca2.end(), samp->getAbonePossibleVal(sa));
										if(found == ca2.end()){
											found = find(ca2.begin(), ca2.end(), samp->getAbtwoPossibleVal(sa));
											if(found == ca2.end()){
												samp->remMSPossible(sa);
												sa--;
												removed++;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else{
		vector<Sample*>* samps = fam->getSamples();
		for(int s = 0; s < (int)samps->size(); s++){
			Sample* samp = (*samps)[s];
			if(samp->isEnabled()){
				samp->clearMEsaved();
				bool child = false;
				Sample* dad = samp->getDad();
				Sample* mom = samp->getMom();
				//untyped child vs typed parents
				if((dad || mom) && (samp->getAone(mloc) && samp->getAtwo(mloc) && samp->getAmissing(mloc))){
					if(dad != NULL){
						child = true;
						if(!(dad->getAone(mloc) && dad->getAtwo(mloc) && dad->getAmissing(mloc))){
							vector<bool> da1 = dad->getAonePossible();
							vector<bool> da2 = dad->getAtwoPossible();
							vector<bool> da3 = dad->getAmissingPossible();
							for(int ca = 0; ca < (int)samp->getAonePossible().size(); ca++){
								vector<bool>::iterator found = find(da1.begin(), da1.end(), samp->getAonePossibleVal(ca));
								if(found == da1.end()){
									found = find(da2.begin(), da2.end(), samp->getAonePossibleVal(ca));
									if(found == da2.end()){
										samp->remPossible(ca);
										removed++;
										ca--;
									}
								}

							}
						}
					}
					if(mom != NULL){
						child = true;
						if(!(mom->getAone(mloc) && mom->getAtwo(mloc) && mom->getAmissing(mloc))){
							vector<bool> ma1 = mom->getAonePossible();
							vector<bool> ma2 = mom->getAtwoPossible();
							vector<bool> ma3 = mom->getAmissingPossible();
							for(int ca = 0; ca < (int)samp->getAtwoPossible().size(); ca++){
								vector<bool>::iterator found = find(ma1.begin(), ma1.end(), samp->getAtwoPossibleVal(ca));
								if(found == ma1.end()){
									found = find(ma2.begin(), ma2.end(), samp->getAtwoPossibleVal(ca));
									if(found == ma2.end()){
										samp->remPossible(ca);
										removed++;
										ca--;
									}
								}
							}
						}
					}
				}
				//untyped parent vs typed children
				if(!dad && !mom && !child){
					if(samp->getAone(mloc) && samp->getAtwo(mloc) && samp->getAmissing(mloc)){
						vector<Sample*>* children = samp->getChildren();
						for(int c = 0; c < (int)children->size(); c++){
							Sample* child = (*children)[c];
							if(!(child->getAone(mloc) && child->getAtwo(mloc) && child->getAmissing(mloc))){
								vector<bool> ca1 = child->getAonePossible();
								vector<bool> ca2 = child->getAtwoPossible();
								vector<bool> ca3 = child->getAmissingPossible();
								for(int sa = 0; sa < (int)samp->getAonePossible().size(); sa++){
									if(samp->getSex()){
										vector<bool>::iterator found = find(ca1.begin(), ca1.end(), samp->getAonePossibleVal(sa));
										if(found == ca1.end()){
											found = find(ca1.begin(), ca1.end(), samp->getAtwoPossibleVal(sa));
											if(found == ca1.end()){
												samp->remPossible(sa);
												sa--;
												removed++;
											}
										}
									}
									else{
										vector<bool>::iterator found = find(ca2.begin(), ca2.end(), samp->getAonePossibleVal(sa));
										if(found == ca2.end()){
											found = find(ca2.begin(), ca2.end(), samp->getAtwoPossibleVal(sa));
											if(found == ca2.end()){
												samp->remPossible(sa);
												sa--;
												removed++;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return removed;
}


void MendelianErrors::determineZygotes(Family* fam, Marker* mark){
////	int mloc = mark->getLoc();

	vector<Sample*>* samps = fam->getSamples();

	for(int s = 0; s < (int)samps->size(); s++){
		Sample* samp = (*samps)[s];
		if(samp->isEnabled()){
			Sample* dad = samp->getDad();
			Sample* mom = samp->getMom();
			if(!mark->isMicroSat()){
				vector<bool> da1;
				vector<bool> da2;
				vector<bool> ma1;
				vector<bool> ma2;
				if(dad != NULL){
					da1 = dad->getAonePossible();
					da2 = dad->getAtwoPossible();
				}
				if(mom != NULL){
					ma1 = mom->getAonePossible();
					ma2 = mom->getAtwoPossible();
				}

				vector<bool> dadalleles;
				vector<bool> momalleles;
				for(int p = 0; p < (int)da1.size(); p++){
					if(!(da1[p] && !da2[p])){
						vector<bool>::iterator found = find(dadalleles.begin(), dadalleles.end(), da1[p]);
						if(found == dadalleles.end()){
							dadalleles.push_back(da1[p]);
						}
						found = find(dadalleles.begin(), dadalleles.end(), da2[p]);
						if(found == dadalleles.end()){
							dadalleles.push_back(da2[p]);
						}
					}
				}
				for(int p = 0; p < (int)ma1.size(); p++){
					if(!(ma1[p] && !ma2[p])){
						vector<bool>::iterator found = find(momalleles.begin(), momalleles.end(), ma1[p]);
						if(found == momalleles.end()){
							momalleles.push_back(ma1[p]);
						}
						found = find(momalleles.begin(), momalleles.end(), ma2[p]);
						if(found == momalleles.end()){
							momalleles.push_back(ma2[p]);
						}
					}
				}
				vector<bitset<2> > zygotes;
				for(int d = 0; d < (int)dadalleles.size(); d++){
					for(int m = 0; m < (int)momalleles.size(); m++){
						if(!(dadalleles[d] && !momalleles[m])){
							bitset<2> b;
							b.reset();
							if(momalleles[m]){
								b.set(1);
							}
							if(dadalleles[d]){
								b.set(0);
							}
							vector<bitset<2> >::iterator found = find(zygotes.begin(), zygotes.end(), b);
							if(found == zygotes.end()){
								zygotes.push_back(b);
							}
						}
						if(!(momalleles[m] && !dadalleles[d])){
							bitset<2> b;
							b.reset();
							if(momalleles[m]){
								b.set(0);
							}
							if(dadalleles[d]){
								b.set(1);
							}
							vector<bitset<2> >::iterator found = find(zygotes.begin(), zygotes.end(), b);
							if(found == zygotes.end()){
								zygotes.push_back(b);
							}
						}
					}
				}
				for(int m = 0; m < (int)momalleles.size(); m++){
					for(int d = 0; d < (int)dadalleles.size(); d++){
						if(!(dadalleles[d] && !momalleles[m])){
							bitset<2> b;
							b.reset();
							if(momalleles[m]){
								b.set(1);
							}
							if(dadalleles[d]){
								b.set(0);
							}
							vector<bitset<2> >::iterator found = find(zygotes.begin(), zygotes.end(), b);
							if(found == zygotes.end()){
								zygotes.push_back(b);
							}
						}
						if(!(momalleles[m] && !dadalleles[d])){
							bitset<2> b;
							b.reset();
							if(momalleles[m]){
								b.set(0);
							}
							if(dadalleles[d]){
								b.set(1);
							}
							vector<bitset<2> >::iterator found = find(zygotes.begin(), zygotes.end(), b);
							if(found == zygotes.end()){
								zygotes.push_back(b);
							}
						}
					}
				}

			}
			else{//microsatellites
				vector<int> da1;
				vector<int> da2;
				vector<int> ma1;
				vector<int> ma2;
				if(dad != NULL){
					da1 = dad->getAbonePossible();
					da2 = dad->getAbtwoPossible();
				}
				if(mom != NULL){
					ma1 = mom->getAbonePossible();
					ma2 = mom->getAbtwoPossible();
				}
				vector<int> dadalleles;
				vector<int> momalleles;
				for(int p = 0; p < (int)da1.size(); p++){
					if(!(da1[p] && !da2[p])){
						vector<int>::iterator found = find(dadalleles.begin(), dadalleles.end(), da1[p]);
						if(found == dadalleles.end()){
							dadalleles.push_back(da1[p]);
						}
						found = find(dadalleles.begin(), dadalleles.end(), da2[p]);
						if(found == dadalleles.end()){
							dadalleles.push_back(da2[p]);
						}
					}
				}
				for(int p = 0; p < (int)ma1.size(); p++){
					if(!(ma1[p] && !ma2[p])){
						vector<int>::iterator found = find(momalleles.begin(), momalleles.end(), ma1[p]);
						if(found == momalleles.end()){
							momalleles.push_back(ma1[p]);
						}
						found = find(momalleles.begin(), momalleles.end(), ma2[p]);
						if(found == momalleles.end()){
							momalleles.push_back(ma2[p]);
						}
					}
				}

				//generate zygotes
				vector<vector<int> > zygotes;
				for(int d = 0; d < (int)dadalleles.size(); d++){
					for(int m = 0; m < (int)momalleles.size(); m++){
						if(dadalleles[d] != -1 && momalleles[m] != -1){
							vector<int> geno;
							geno.resize(2,-1);
							geno[0] = dadalleles[d];
							geno[1] = momalleles[m];
							vector<vector<int> >::iterator found = find(zygotes.begin(), zygotes.end(), geno);
							if(found == zygotes.end()){
								zygotes.push_back(geno);
							}
						}
					}
				}
				for(int m = 0; m < (int)momalleles.size(); m++){
					for(int d = 0; d < (int)dadalleles.size(); d++){
						if(dadalleles[d] != -1 && momalleles[m] != -1){
							vector<int> geno;
							geno.resize(2,-1);
							geno[0] = dadalleles[d];
							geno[1] = momalleles[m];
							vector<vector<int> >::iterator found = find(zygotes.begin(), zygotes.end(), geno);
							if(found == zygotes.end()){
								zygotes.push_back(geno);
							}
						}
					}
				}
				if(dad){
					dad->addZygotes(zygotes);
				}
				if(mom){
					mom->addZygotes(zygotes);
				}
			}
		}
	}
}

bool MendelianErrors::getPossibleGenos(Family* fam, Marker* mark){
	bool realgeno = false;
	int mloc = mark->getLoc();
	vector<Sample*>* samps = fam->getSamples();
	for(int s = 0; s < (int)samps->size(); s++){
		Sample* samp = (*samps)[s];
		if(samp->isEnabled()){
			samp->clearPossible();
			samp->clearMEsaved();
			if(!mark->isMicroSat()){
				if(samp->getAone(mloc)){
					if(samp->getAtwo(mloc) && !samp->getAmissing(mloc)){
						realgeno = true;
						samp->addPossible(true, true, false);
					}
					else if(!samp->getAtwo(mloc)){
						realgeno = true;
						samp->addPossible(true, true, false);
						samp->addPossible(false,false, false);
						samp->addPossible(false, true, false);
						samp->addPossible(true, false, false);
					//	samp->addPossible(true, true, true);
					}
					else{
						samp->addPossible(true, true, false);
						samp->addPossible(false,false, false);
						samp->addPossible(false, true, false);
						samp->addPossible(true, false, false);
						samp->addPossible(true, true, true);
					}
				}
				else{
					if(!samp->getAtwo(mloc)){
						realgeno = true;
						samp->addPossible(false, false, false);
					}
					else{
						realgeno = true;
						samp->addPossible(true, true, false);
						samp->addPossible(false,false, false);
						samp->addPossible(false, true, false);
						samp->addPossible(true, false, false);
					//	samp->addPossible(true, true, true);
					}
				}
			}
			else{
				if(samp->getAbone(mloc) != -1 && samp->getAbtwo(mloc) != -1){
					realgeno = true;
					if(samp->getAbone(mloc) == samp->getAbtwo(mloc)){
						samp->addMSPossible(samp->getAbone(mloc), samp->getAbone(mloc));
					}
					else{
						samp->addMSPossible(samp->getAbone(mloc), samp->getAbone(mloc));
						samp->addMSPossible(samp->getAbone(mloc), samp->getAbtwo(mloc));
						samp->addMSPossible(samp->getAbtwo(mloc), samp->getAbtwo(mloc));
						samp->addMSPossible(samp->getAbtwo(mloc), samp->getAbone(mloc));
					}
				}
				else{
					int num = mark->getNumAlleles();
					for(int i = 0; i < num; i++){
						for(int j = 0; j < num; j++){
							samp->addMSPossible(i, j);
						}
					}
				}
			}
		}
	}
	return realgeno;
}

void MendelianErrors::filter_markers(){
	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int msize = good_markers.size();//markers->size();

		for(int i = 0; i < msize; i++){
			if(good_markers[i]->isEnabled()){// && !(*markers)[i]->isFlagged()){
/*				if(options.doChrom()){
					if(!options.checkChrom((*markers)[i]->getChrom())){
				    	continue;
			    	}
			    	if(!options.checkBp((*markers)[i]->getBPLOC())){
				    	continue;
			    	}
				}
*/
				bool inc = false;
				if(options.doThreshMarkersHigh() && merrors[i] > options.getThreshMarkersHigh()){
					(*markers)[i]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersLow() && merrors[i] < options.getThreshMarkersLow()){
					(*markers)[i]->setEnabled(false);
					inc = true;
				}

				if(inc){
					orig_num_markers++;
				}
			}
		}
	}
}

void MendelianErrors::mid_process_output(){
}

void MendelianErrors::resetCounts(){
}

void MendelianErrors::filter(){
	if(options.doThreshFamiliesLow() || options.doThreshFamiliesHigh()){
		int fsize = families->size();

		for(int i = 0; i < fsize; i++){
			if((*families)[i]->isEnabled()){
				bool inc = false;
				if(options.doThreshFamiliesLow() && ferrors[i] < options.getThreshFamiliesLow()){
					(*families)[i]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshFamiliesHigh() && ferrors[i] > options.getThreshFamiliesHigh()){
					(*families)[i]->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_families++;
				}
			}
		}
	}
}

void MendelianErrors::write_error(ofstream &erroroutput, Marker* m, Sample* s, Sample* p, int mloc){
	if(!erroroutput){
		opts::printLog("File handle for Mendelian error output is corrupted.  Exiting!\n");
		//exit(1);
		throw MethodException("File handle for Mendelian error output is corrupted.  Exiting!\n");
	}
	//cout << m->getChrom() << "\t"
	//	<< m->getRSID() << "\t"
	//	<< m->getProbeID() << "\t"
	//	<< m->getEnzyme() << "\t"
	//	<< s->getFamID() << "\t"
	//	<< s->getInd() << "\t";
	erroroutput << m->toString() << "\t"
		<< s->getFamID() << "\t"
		<< s->getInd() << "\t";
	if(s->getSex()){
	//	cout << "M\t";
		erroroutput << "M\t";
	}
	else{
	//	cout << "F\t";
		erroroutput << "F\t";
	}
	if(!m->isMicroSat()){
		if(s->getAone(mloc) && s->getAtwo(mloc)){
			erroroutput << m->getAllele2() << "/" << m->getAllele2();
		}
		else if(!s->getAone(mloc) && !s->getAtwo(mloc)){
			erroroutput << m->getAllele1() << "/" << m->getAllele1();
		}
		else if(!s->getAone(mloc) && s->getAtwo(mloc)){
			erroroutput << m->getAllele1() << "/" << m->getAllele2();
		}
		else if(s->getAone(mloc) && !s->getAtwo(mloc)){
			erroroutput << "0/0";
		}
//		if(s->getAone(mloc)){
//			erroroutput << m->getAllele2() << "/";
//		}
//		else{
//			erroroutput << m->getAllele1() << "/";
//		}
//		if(s->getAtwo(mloc)){
//			erroroutput << m->getAllele2() << "\t";
//		}
//		else{
//			erroroutput << m->getAllele1() << "\t";
//		}
	}
	else{
	//	cout << m->getAllele(s->getAbone(mloc)) << "/";
	//	cout << m->getAllele(s->getAbtwo(mloc)) << "\t";
		erroroutput << m->getAllele(s->getAbone(mloc)) << "/";
		erroroutput << m->getAllele(s->getAbtwo(mloc)) << "\t";
	}
	//cout << p->getInd() << "\t";
	erroroutput << "\t" << p->getInd() << "\t";
	if(p->getSex()){
	//	cout << "M\t";
		erroroutput << "M\t";
	}
	else{
	//	cout << "F\t";
		erroroutput << "F\t";
	}
	if(!m->isMicroSat()){
		if(p->getAone(mloc) && p->getAtwo(mloc)){
			erroroutput << m->getAllele2() << "/" << m->getAllele2();
		}
		else if(!p->getAone(mloc) && !p->getAtwo(mloc)){
			erroroutput << m->getAllele1() << "/" << m->getAllele1();
		}
		else if(!p->getAone(mloc) && p->getAtwo(mloc)){
			erroroutput << m->getAllele1() << "/" << m->getAllele2();
		}
		else if(p->getAone(mloc) && !p->getAtwo(mloc)){
			erroroutput << "0/0";
		}
		//if(p->getAtwo(mloc)){
		//	erroroutput << m->getAllele2() << "";
		//}
		//else{
		//	erroroutput << m->getAllele1() << "";
		//}
	}
	else{
	//	cout << m->getAllele(p->getAbone(mloc)) << "/";
	//	cout << m->getAllele(p->getAbtwo(mloc)) << "";
		erroroutput << m->getAllele(p->getAbone(mloc)) << "/";
		erroroutput << m->getAllele(p->getAbtwo(mloc)) << "";
	}
	erroroutput << endl;
}

}
