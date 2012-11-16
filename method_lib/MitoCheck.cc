/**********************************************************************************
*                       Mito Check Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Iterates over all Male genotypes on the X Chromosome and generates a count of
* all genotyping errors (not homozygous).  Outputs include counts by individual
* and counts by marker.
*
* Performs an initial scan.  Filters out any Markers with counts exceeding
* threshold.  Performs a secondary scan with bad markers removed.
*
* Files generated:
*	gender_errors_ind.txt
*	gender_errors_markers.txt
*       post_gender_error_filter_summary.txt
*       post_gender_error_filter_summary_chrom.txt
*
*File: MitoCheck.cc
**********************************************************************************/


#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include "MitoCheck.h"
//#include "Chrom.h"
#include "General.h"
#include "Helper.h"
namespace Methods{
string MitoCheck::stepname = "mito-check";

void MitoCheck::calcThreshold(){
/*	int nummarkers = markers->getNumOfType("X");
	int numfams = families->getSize();
	int nummales = families->getNumMales();

	if(marker_thresh == -1.0){
		marker_thresh = (int)((float)nummales * error_rate) + 1;
	}
	if(ind_thresh == -1.0){
		ind_thresh = (int)((float)nummarkers * error_rate) + 1;
	}
*/
}



void MitoCheck::Tokenize(const string& str, vector<string>& tokens, const string& delimiter){
    string::size_type lastPos = str.find_first_not_of(delimiter, 0);
    string::size_type pos = str.find_first_of(delimiter, lastPos);

    while(string::npos != pos || string::npos != lastPos){
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiter, pos);
        pos = str.find_first_of(delimiter, lastPos);
    }
}

void MitoCheck::setThreshold(string thresh){
	options.setUp(thresh);
}

void MitoCheck::PrintSummary(){
	int msize = markers->size();
	int ssize = samples->size();

	string fname1 = opts::_OUTPREFIX_ + "mito_check_marker" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	ofstream myoutput (fname1.c_str());
	if(!myoutput){
		opts::printLog("Error opening " + fname1 + "!  Exiting!\n");
		exit(1);
	}
	opts::addFile("Marker",stepname,fname1);
	myoutput.precision(4);
	myoutput << "Chrom\trsID\tProbeID\tbploc";
	if((*markers)[0]->getDetailHeaders().size() > 0){
		myoutput << "\t" << (*markers)[0]->getDetailHeaders();
	}
	myoutput << "\tNum_errors" << endl;
	opts::addHeader(fname1, "Num_errors");

	for(int i = 0; i < msize; i++){
		if((*markers)[i]->isEnabled() && (*markers)[i]->getChrom() == opts::_MITO_ && !(*markers)[i]->isFlagged()){
			if(options.doChrom()){
				if(!options.checkChrom((*markers)[i]->getChrom())){
					continue;
				}
				if(!options.checkBp((*markers)[i]->getBPLOC())){
					continue;
				}
			}

			myoutput << (*markers)[i]->toString() << "\t"
				<< merrors[i] << endl;
		}
	}
	string fname2 = opts::_OUTPREFIX_ + "mito_check_ind" + options.getOut() + ".txt";//getString<int>(order) + ".txt";
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	ofstream hetout (fname2.c_str());
	if(!hetout){
		opts::printLog("Error opening " + fname2 + "!  Exiting!!\n");
		exit(1);
	}
	opts::addFile("Sample",stepname,fname2);
	hetout.precision(4);
	hetout << "FamID\t";
	hetout << "IndID\tCenter\tSex\t";
	hetout << "Affection_Status\t";
	hetout << "Plate\tWell\tNum_errors";
	opts::addHeader(fname2, "Num_errors");

    string sdetails = "";
    if(opts::_SAMPDESC_.length() > 0){
        sdetails = (*samples)[0]->getDetailHeaders();
    	hetout << "\t" << sdetails;
	}

    hetout << endl;

	for(int i = 0; i < ssize; i++){
		if((*samples)[i]->isEnabled()){
			hetout << (*samples)[i]->toString() << "\t"
				<< (*samples)[i]->getFamily()->getCenter() << "\t";
			if((*samples)[i]->getSex()){
				hetout << "M\t";
			}
			else{
				hetout << "F\t";
			}
			if((*samples)[i]->getPheno() == 2){
				hetout << "Y\t";
			}
			else if((*samples)[i]->getPheno() == 1){
				hetout << "N\t";
			}
			else{
				hetout << "U\t";
			}

			hetout << (*samples)[i]->getPlate() << "\t"
				<< (*samples)[i]->getWell() << "\t"
				<< serrors[i];
			if(opts::_SAMPDESC_.length() > 0){
				hetout << "\t" << (*samples)[i]->getDetails();
			}
			hetout << endl;
		}
	}
	string fname3 = opts::_OUTPREFIX_ + "mito_check_errors_" + getString<int>(order) + ".txt";
	ofstream errorout (fname3.c_str());
	if(!errorout){
		opts::printLog("Error opening " + fname3 + "!  Exiting!!\n");
		exit(1);
	}
	errorout << "Chrom\trsID\tProbeID\tbploc\tFamID\tMom\tMom_Genotype\tChild\tChild_Genotype\n";

	for(int i = 0; i < (int)error_map.size(); i++){
		Sample* child = (*samples)[i];
		Sample* mom = child->getMom();
		for(int j = 0; j < (int)error_map[i].size(); j++){
			Marker* mark = (*markers)[error_map[i][j]];
			int loc = mark->getLoc();
			string ma1;
			string ma2;
			string ca1;
			string ca2;
			if(!mark->isMicroSat()){
				if(child->getAone(loc)){
					if(child->getAtwo(loc)){
						ca1 = mark->getAllele2();
						ca2 = mark->getAllele2();
					}
					else{
						ca1 = "0";
						ca2 = "0";
					}
				}
				else{
					if(!child->getAtwo(loc)){
						ca1 = mark->getAllele1();
						ca2 = mark->getAllele1();
					}
					else{
						ca1 = mark->getAllele1();
						ca2 = mark->getAllele2();
					}
				}
				if(mom->getAone(loc)){
					if(mom->getAtwo(loc)){
						ma1 = mark->getAllele2();
						ma2 = mark->getAllele2();
					}
					else{
						ca1 = "0";
						ca2 = "0";
					}
				}
				else{
					if(!mom->getAtwo(loc)){
						ma1 = mark->getAllele1();
						ma2 = mark->getAllele1();
					}
					else{
						ma1 = mark->getAllele1();
						ma2 = mark->getAllele2();
					}
				}
			}
			else{
				ca1 = mark->getAllele(child->getAbone(loc));
				ca2 = mark->getAllele(child->getAbtwo(loc));
				ma1 = mark->getAllele(child->getAbone(loc));
				ma2 = mark->getAllele(child->getAbtwo(loc));
			}
			errorout << mark->toString() << "\t" << child->getFamID() << "\t" << mom->getInd() << "\t" << ma1 << "_" << ma2 << "\t" << child->getInd() << "\t" << ca1 << "_" << ca2 << endl;
		}

	}

	if(hetout.is_open()){
		hetout.close();
	}
	if(myoutput.is_open()){
		myoutput.close();
	}

	if(errorout.is_open()){
		errorout.close();
	}
	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}
}

void MitoCheck::FilterSummary(){
	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::printLog("Individuals Passed:\t" + getString<int>(opts::_SAMPLES_WORKING_ - orig_num_samples) + " (" +
	    getString<float>(((float) (opts::_SAMPLES_WORKING_ - orig_num_samples) / (float) opts::_SAMPLES_WORKING_) * 100.0) +
	    "%) of " + getString<int>(opts::_SAMPLES_WORKING_) + "\n");

	opts::_MARKERS_WORKING_ -= orig_num_markers;
	opts::_SAMPLES_WORKING_ -= orig_num_samples;

}

void MitoCheck::filter(){

	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int msize = markers->size();
		for(int m = 0; m < msize; m++){
			if((*markers)[m]->isEnabled() && (*markers)[m]->getChrom() == opts::_MITO_ && !(*markers)[m]->isFlagged()){
				if(options.doChrom()){
					if(!options.checkChrom((*markers)[m]->getChrom())){
					    continue;
				    }
				    if(!options.checkBp((*markers)[m]->getBPLOC())){
					    continue;
				    }
				}

				bool inc = false;
				if(options.doThreshMarkersLow() && merrors[m] < options.getThreshMarkersLow()){
					(*markers)[m]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshMarkersHigh() && merrors[m] > options.getThreshMarkersHigh()){
					(*markers)[m]->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_markers++;
				}
			}
		}
	}

	if(options.doThreshSamplesLow() || options.doThreshSamplesHigh()){
		int ssize = samples->size();
		for(int i = 0; i < ssize; i++){
			if((*samples)[i]->isEnabled()){
				bool inc = false;
				if(options.doThreshSamplesLow() && serrors[i] < options.getThreshSamplesLow()){
					(*samples)[i]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshSamplesHigh() && serrors[i] > options.getThreshSamplesHigh()){
					(*samples)[i]->setEnabled(false);
					inc = true;
				}
				if(inc){
					orig_num_samples++;
				}
			}
		}
	}
}

void MitoCheck::filter_markers(){


/*	MKR* mymarkers = markers->getList();
	MKR::iterator m_iter;

ofstream myoutput ("gender_errors_markers.txt", ios::out | ios::app);
	for(m_iter = mymarkers->begin(); m_iter != mymarkers->end();){
		myoutput << m_iter->second.getChrom() << "\t" <<  m_iter->second.getRSID() << "\t" << m_iter->second.getBPLOC() << "\t" << m_iter->second.getMitoErrors() << endl;
		if(m_iter->second.getMitoErrors() > marker_thresh){
			mymarkers->erase(m_iter++);
		}
		else{
			++m_iter;
		}
	}
	myoutput.close();
*/
}


void MitoCheck::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

	calcThreshold();

	perform_evaluation(false);
}
void MitoCheck::perform_evaluation(bool dofams){
	int msize = markers->size();
//	int fsize = families->size();
	int ssize = samples->size();

	merrors.resize(msize);
	serrors.resize(ssize);
	error_map.resize(ssize);

	int prev_base = 0;
	int prev_chrom = -1;
	for(int m = 0; m < msize; m++){
		if((*markers)[m]->isEnabled()){
			if((*markers)[m]->getChrom() != opts::_MITO_){
				continue;
			}
			if(options.doChrom()){
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


			int mloc = (*markers)[m]->getLoc();
			bool micro = (*markers)[m]->isMicroSat();
			for(int i = 0; i < ssize; i++){
				if((*samples)[i]->isEnabled()){
					if(!micro){
						if((*samples)[i]->getAone(mloc) && (*samples)[i]->getAtwo(mloc) && (*samples)[i]->getAmissing(mloc)){
							continue;
						}
						if(!(*samples)[i]->getAone(mloc) && (*samples)[i]->getAtwo(mloc)){
							serrors[i]++;
							merrors[m]++;
						}
						if((*samples)[i]->getMom()){
							Sample* mom = (*samples)[i]->getMom();
							if(mom->getAone(mloc) == mom->getAtwo(mloc) && !mom->getAmissing(mloc) && (*samples)[i]->getAone(mloc) == (*samples)[i]->getAtwo(mloc) && mom->getAone(mloc) != (*samples)[i]->getAone(mloc)){
								serrors[i]++;
								merrors[m]++;
								error_map[i].push_back(m);
							}
						}
					}
					else{
						if((*samples)[i]->getAbone(mloc) == -1){
							continue;
						}
						if((*samples)[i]->getAbone(mloc) != (*samples)[i]->getAbtwo(mloc)){
							serrors[i]++;
							merrors[m]++;
						}
						if((*samples)[i]->getMom()){
							Sample* mom = (*samples)[i]->getMom();
							if(mom->getAbone(mloc) == mom->getAbtwo(mloc) && (*samples)[i]->getAbone(mloc) == (*samples)[i]->getAbtwo(mloc) && mom->getAbone(mloc) == (*samples)[i]->getAbone(mloc)){
								serrors[i]++;
								merrors[m]++;
								error_map[i].push_back(m);
							}
						}
					}
				}
			}
		}
	}

}

}
