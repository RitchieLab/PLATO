/**********************************************************************************
*                       Gender Check Module
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
*File: Heterozygosity.cc
**********************************************************************************/


#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include "Heterozygosity.h"
#include "Chrom.h"
#include <General.h>
#include <Helper.h>


void Heterozygosity::calcThreshold(){
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



void Heterozygosity::Tokenize(const string& str, vector<string>& tokens, const string& delimiter){
    string::size_type lastPos = str.find_first_not_of(delimiter, 0);
    string::size_type pos = str.find_first_of(delimiter, lastPos);

    while(string::npos != pos || string::npos != lastPos){
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = str.find_first_not_of(delimiter, pos);
        pos = str.find_first_of(delimiter, lastPos);
    }
}

void Heterozygosity::setThreshold(string thresh){
	options.setUp(thresh);
	//vector<string> tokens;
//    Tokenize(thresh, tokens, ":");
//    if(tokens.size() == 2){
//        ind_thresh = std::atof(tokens[0].c_str());
//        marker_thresh = std::atof(tokens[1].c_str());
//    }
//    else{
//		opts::printLog("Incorrect number of Gender Error thresholds.  Value should be #:#\n");
//        exit(1);
//    }
}

void Heterozygosity::PrintSummary(){
	int msize = markers->size();
	int ssize = samples->size();

	string fname1 = opts::_OUTPREFIX_ + "GE_Male_ChrX_error_marker_" + getString<int>(order) + ".txt";
	ofstream myoutput (fname1.c_str());
	if(!myoutput){
		opts::printLog("Error opening GE_Male_ChrX_error_marker.txt!  Exiting!\n");
		exit(1);
	}
	opts::addFile(fname1);
	myoutput.precision(4);
	myoutput << "Chrom\trsID\tProbeID\tbploc\tEnzyme\tCount_Geno_all_Males\tCount_Geno_HET_Males\t%HET_Males" << endl;

	for(int i = 0; i < msize; i++){
		if((*markers)[i]->isEnabled() && (*markers)[i]->getChrom() == otps::_CHRX_ && !(*markers)[i]->isFlagged()){
			if(options.doChrom()){
				if(!options.checkChrom((*markers)[i]->getChrom())){
					continue;
				}
				if(!options.checkBp((*markers)[i]->getBPLOC())){
					continue;
				}
			}

			myoutput << (*markers)[i]->toString() << "\t"
				<< mtotal[i] << "\t"
				<< merrors[i] << "\t";
			if(mtotal[i] > 0){
				myoutput << (((float)merrors[i]/(float)mtotal[i]));// * 100.0f);
			}
			else{
				myoutput << "0";
			}
			myoutput << endl;
		}
	}
	string fname2 = opts::_OUTPREFIX_ + "GE_Individual_ChrX_HET_" + getString<int>(order) + ".txt";
	ofstream hetout (fname2.c_str());
	if(!hetout){
		opts::printLog("Error opening GE_Individual_ChrX_HET.txt.  Exiting!!\n");
		exit(1);
	}
	opts::addFile(fname2);
	hetout.precision(4);
	hetout << "FamID\t";
	hetout << "IndID\tCenter\tSex\t";
	hetout << "Affection_Status\t";
	hetout << "Plate\tWell\tCount_Geno_total_all\tCount_Geno_HET_all\t%HET_all";

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
        for(int i = 0; i < enzymes.size(); i++){
			hetout << "\tCount_Geno_total_" << enzymes[i] << "\tCount_Geno_HET_" << enzymes[i] << "\t%HET_" << enzymes[i];
        }
    }
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
				<< stotal[i] << "\t"
				<< shets[i] << "\t";
			if(stotal[i] > 0){
				hetout << (((float)shets[i]/(float)stotal[i]));// * 100.0f);
			}
			else{
				hetout << "0";
			}

			if(opts::_ENZYMES_ && enzymes.size() > 0){
				for(int e = 0; e < enzymes.size(); e++){
					hetout << "\t" << senzyme_tot[i][enzymes[e]]
						<< "\t" << senzyme_hets[i][enzymes[e]];
					if(senzyme_tot[i][enzymes[e]] > 0){
						hetout << "\t" << ((1.0f - ((float)senzyme_hets[i][enzymes[e]]/(float)senzyme_tot[i][enzymes[e]])));// * 100.0f);
					}
					else{
						hetout << "\t" << "0";
					}
				}
			}
			if(opts::_SAMPDESC_.length() > 0){
				hetout << "\t" << (*samples)[i]->getDetails();
			}
			hetout << endl;
		}
	}

	if(hetout.is_open()){
		hetout.close();
	}
	if(myoutput.is_open()){
		myoutput.close();
	}

	for(int i = 0; i < msize; i++){
		(*markers)[i]->setFlag(false);
	}
}

void Heterozygosity::FilterSummary(){
//    string fname = opts::_OUTPREFIX_ + "post_gender_error_filter_summary_" + getStringInt(order) + ".txt";
//	ofstream myoutput (fname.c_str());
  //  myoutput.precision(4);
//	ofstream chromoutput ("post_gender_error_filter_summary_chrom.txt");
//	chromoutput.precision(4);
//	chromoutput << "Chrom\t#Markers Passed" << endl;

	opts::printLog("Options:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
        getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
        "%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
    //myoutput << "Threshold:\t" << "N/A" << endl;
    //myoutput << "Families Passed:\t" << (opts::_FAMILIES_WORKING_ - orig_num_families) << " (" <<
    //    ((float)(opts::_FAMILIES_WORKING_ - orig_num_families) / (float)opts::_FAMILIES_WORKING_) * 100.0 <<
    //    "%) of " << opts::_FAMILIES_WORKING_ << endl;
	opts::printLog("Individuals Passed:\t" + getString<int>(opts::_SAMPLES_WORKING_ - orig_num_samples) + " (" +
	    getString<float>(((float) (opts::_SAMPLES_WORKING_ - orig_num_samples) / (float) opts::_SAMPLES_WORKING_) * 100.0) +
	    "%) of " + getString<int>(opts::_SAMPLES_WORKING_) + "\n");

//	myoutput << "Threshold:\t" << marker_thresh << endl;
//    myoutput << "Markers Passed:\t" << (opts::_MARKERS_WORKING_ - orig_num_markers) << " (" <<
//        ((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0 <<
//        "%) of " << opts::_MARKERS_WORKING_ << endl;
    //myoutput << "Threshold:\t" << "N/A" << endl;
    //myoutput << "Families Passed:\t" << (opts::_FAMILIES_WORKING_ - orig_num_families) << " (" <<
    //    ((float)(opts::_FAMILIES_WORKING_ - orig_num_families) / (float)opts::_FAMILIES_WORKING_) * 100.0 <<
    //    "%) of " << opts::_FAMILIES_WORKING_ << endl;
//	myoutput << "Threshold:\t" << ind_thresh << endl;
//    myoutput << "Individuals Passed:\t" << (opts::_SAMPLES_WORKING_ - orig_num_samples) << " (" <<
//	    ((float) (opts::_SAMPLES_WORKING_ - orig_num_samples) / (float) opts::_SAMPLES_WORKING_) * 100.0 <<
//	    "%) of " << opts::_SAMPLES_WORKING_ << endl;

 //   myoutput.close();
	opts::_MARKERS_WORKING_ -= orig_num_markers;
	opts::_SAMPLES_WORKING_ -= orig_num_samples;

}

void Heterozygosity::filter(){

	if(options.doThreshMarkersLow() || options.doThreshMarkersHigh()){
		int msize = markers->size();
		for(int m = 0; m < msize; m++){
			if((*markers)[m]->isEnabled() && !(*markers)[m]->isFlagged()){
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
			if((*samples)[i]->isEnabled() && (*samples)[i]->getSex()){
				bool inc = false;
				if(options.doThreshSamplesLow() && shets[i] < options.getThreshSamplesLow()){
					(*samples)[i]->setEnabled(false);
					inc = true;
				}
				if(options.doThreshSamplesHigh() && shets[i] > options.getThreshSamplesHigh()){
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

void Heterozygosity::filter_markers(){


/*	MKR* mymarkers = markers->getList();
	MKR::iterator m_iter;

ofstream myoutput ("gender_errors_markers.txt", ios::out | ios::app);
	for(m_iter = mymarkers->begin(); m_iter != mymarkers->end();){
		myoutput << m_iter->second.getChrom() << "\t" <<  m_iter->second.getRSID() << "\t" << m_iter->second.getBPLOC() << "\t" << m_iter->second.getGenderErrors() << endl;
		if(m_iter->second.getGenderErrors() > marker_thresh){
			mymarkers->erase(m_iter++);
		}
		else{
			++m_iter;
		}
	}
	myoutput.close();
*/
}

void Heterozygosity::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm){
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;

	calcThreshold();

	perform_evaluation(false);
}
void Heterozygosity::perform_evaluation(bool dofams){
	int msize = markers->size();
	int fsize = families->size();
	int ssize = samples->size();

	merrors.resize(msize);
	shets.resize(ssize);
	mtotal.resize(msize);
	stotal.resize(ssize);
	senzyme_hets.resize(ssize);
	senzyme_tot.resize(ssize);

	int prev_base = 0;
	int prev_chrom = -1;
	for(int m = 0; m < msize; m++){
		if((*markers)[m]->isEnabled()){
			if((*markers)[m]->getChrom() != opts::_CHRX_){
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
						if((*samples)[i]->getAone(mloc) && !(*samples)[i]->getAtwo(mloc)){
						continue;
						}
						if(!(*samples)[i]->getAone(mloc) && (*samples)[i]->getAtwo(mloc)){
							shets[i]++;
							if(opts::_ENZYMES_){
								senzyme_hets[i][(*markers)[m]->getEnzyme()]++;
							}
							if((*samples)[i]->getSex()){
								merrors[m]++;
							}
						}
					}
					else{
						if((*samples)[i]->getAbone(mloc) == -1){
							continue;
						}
						if((*samples)[i]->getAbone(mloc) != (*samples)[i]->getAbtwo(mloc)){
							shets[i]++;
							if(opts::_ENZYMES_){
								senzyme_hets[i][(*markers)[m]->getEnzyme()]++;
							}
							if((*samples)[i]->getSex()){
								merrors[m]++;
							}
						}
					}
					if((*samples)[i]->getSex()){
						mtotal[m]++;
					}
					stotal[i]++;
					if(opts::_ENZYMES_){
						senzyme_tot[i][(*markers)[m]->getEnzyme()]++;
					}
				}
			}
		}
	}

}
