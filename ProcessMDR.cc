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
*File: LD.cc
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
#include "ProcessMDR.h"
#include <Options.h>
#include <General.h>
#include <Helpers.h>
//#include "Markers.h"
//#include "Chrom.h"
//#include "Families.h"
using namespace Methods;

void ProcessMDR::FilterSummary(){

	opts::printLog("Threshold:\t" + options.toString() + "\n");
	opts::printLog("Markers Passed:\t" + getString<int>(opts::_MARKERS_WORKING_ - orig_num_markers) + " (" +
		getString<float>(((float)(opts::_MARKERS_WORKING_ - orig_num_markers) / (float)opts::_MARKERS_WORKING_) * 100.0) +
		"%) of " + getString<int>(opts::_MARKERS_WORKING_) + "\n");
	opts::_MARKERS_WORKING_ -= orig_num_markers;

}

void ProcessMDR::PrintSummary(){
	int msize = data_set->num_loci();
	for(int m = 0; m < msize; m++){
		data_set->get_locus(m)->setFlag(false);
	}

}

void ProcessMDR::filter(){
}


void ProcessMDR::process(DataSet* ds){
	data_set = ds;

	//check if new covariate file is listed...or covariate name.
	//create vector of covariate indexes to use if specified.

    string fname = opts::_OUTPREFIX_ + "mdr" + options.getOut() + ".txt";
    if(!overwrite){
        fname += "." + getString<int>(order);
    }
    ofstream mdrout (fname.c_str());
    if(!mdrout){
        opts::printLog("Error opening " + fname + "!  Exiting!\n");
        throw MethodException("");
    }
    mdrout.precision(4);

    mdrout << "Chrom\trsID\tProbeID\tBPLOC\tTrain_Balanced_Accuracy\tTest_Balanced_Accuracy\tTrain_False_Neg\tTest_False_Neg\tTrain_False_Pos\tTest_False_Pos\tTrain_True_Neg\tTest_True_Neg\tTrain_True_Pos\tTest_True_Pos\n";

    //Set up covariate indexes
    vector<unsigned int> traits; //not used.
    vector<unsigned int> covariates;
    vector<string> cov_use = options.getCovars();

    for(unsigned int i = 0; i < cov_use.size(); i++){
    	int index = data_set->get_covariate_index(cov_use[i]);
    	if(index > -1){
    		covariates.push_back((unsigned int)index);
    	}
    }

    //initialize training & testing groups by iterating over first 2 sets of groups in Group File
    //if no group file, then randomly choose which samples to use, trying to get the same ratio of
    //case to control in each set.
    DataSet* training_data_set = new DataSet();
    DataSet* testing_data_set = new DataSet();
    training_data_set->set_markers(data_set->get_markers());
    testing_data_set->set_markers(data_set->get_markers());
    training_data_set->set_covariates(data_set->get_covariates());
    testing_data_set->set_covariates(data_set->get_covariates());
    if(options.doGroupFile()){
    	options.readGroups(data_set->get_samples());
    	map<string, vector<Sample*> > groups = options.getGroups();
    	int count = 0;
    	map<string, vector<Sample*> >::iterator iter;
    	for(iter = groups.begin(); iter != groups.end(); iter++){
    		if(count == 2){
    			break;
    		}
    		if(count == 0){
    			training_data_set->set_samples(&(iter->second));
    		}
    		else{
    			testing_data_set->set_samples(&(iter->second));
    		}
    		count++;
    	}
    	training_data_set->recreate_family_vector();
    	training_data_set->set_affection_vectors();
    	testing_data_set->recreate_family_vector();
    	testing_data_set->set_affection_vectors();
    }
    else{
    	opts::printLog("Randomly generating two case/control datasets.\n");
    	vector<DataSet*> new_sets = data_set->generate_case_control_subsets(2);
    	if(new_sets.size() == 2){
    		training_data_set = new_sets[0];
    		testing_data_set = new_sets[1];
    	}
    	else{
    		throw MethodException("Could not create training & testing data sets.");
    	}

    }
	opts::printLog("Training Dataset Summary:");
	string summary = "Samples: " + getString<int>(training_data_set->num_inds()) + "\n";
	summary += "Cases: " + getString<int>(training_data_set->num_affected()) + "\n";
	summary += "Controls: " + getString<int>(training_data_set->num_unaffected()) + "\n";
	summary += "Males: " + getString<int>(training_data_set->num_males()) + "\n";
	summary += "Females: " + getString<int>(training_data_set->num_females()) + "\n";
	opts::printLog(summary);

	opts::printLog("Testing Dataset Summary:");
	summary = "Samples: " + getString<int>(testing_data_set->num_inds()) + "\n";
	summary += "Cases: " + getString<int>(testing_data_set->num_affected()) + "\n";
	summary += "Controls: " + getString<int>(testing_data_set->num_unaffected()) + "\n";
	summary += "Males: " + getString<int>(testing_data_set->num_males()) + "\n";
	summary += "Females: " + getString<int>(testing_data_set->num_females()) + "\n";
	opts::printLog(summary);


    MDR mdr_train;
    mdr_train.resetDataSet(training_data_set);
    mdr_train.set_parameters(&options);

    MDR mdr_test;
    mdr_test.resetDataSet(testing_data_set);
    mdr_test.set_parameters(&options);


    opts::printLog("Performing MDR calculation on training and testing data sets.\n");

    int prev_base = 0;
    int prev_chrom = -1;
    int msize = data_set->num_loci();
    for(int m = 0; m < msize; m++){
    	Marker* mark = data_set->get_locus(m);
    	if(mark->isEnabled() && Helpers::isValidMarker(mark, &options, prev_chrom, prev_base)){
    		vector<unsigned int> loci;
    		loci.push_back((unsigned int)m);
    		if(covariates.size() == 0){
    			mdr_train.calculate(loci);
    			mdr_test.calculate(loci);
    		}
    		else{
    			mdr_train.calculate(loci,covariates,traits);
    			mdr_test.calculate(loci, covariates,traits);
    		}

    		vector<unsigned int> trainaffcells = mdr_train.getAffectedCellTotals();
    		vector<unsigned int> trainunaffcells = mdr_train.getUnaffectedCellTotals();
    		vector<unsigned int> testaffcells = mdr_test.getAffectedCellTotals();
    		vector<unsigned int> testunaffcells = mdr_test.getUnaffectedCellTotals();

    		map<vector<string>, unsigned short > mapping;
    		int cellsize = testaffcells.size();
////    		int traincellsize = trainaffcells.size();

    		vector<string> rules;

    	    map<vector<string>, int> testing_aff_map;
    	    map<vector<string>, int> testing_unaff_map;
    	    map<vector<string>, int> training_aff_map;
    	    map<vector<string>, int> training_unaff_map;

    	    for(int i = 0; i < cellsize; i++){
    			vector<unsigned int> genos = mdr_test.genotypes_at_index(i, (loci.size() + covariates.size() + traits.size()));
    			vector<string> conv = testing_data_set->convert_geno_tostring(genos, loci, covariates, traits);
    			testing_aff_map[conv] = testaffcells[i];
    			testing_unaff_map[conv] = testunaffcells[i];

    			genos = mdr_train.genotypes_at_index(i, (loci.size() + covariates.size() + traits.size()));
    			conv = training_data_set->convert_geno_tostring(genos, loci, covariates, traits);
    			training_aff_map[conv] = trainaffcells[i];
    			training_unaff_map[conv] = trainunaffcells[i];

    	    }

    	    map<vector<string>, int>::iterator giter;
    	    map<vector<string>, int>::iterator findme;
    	    for(giter = training_aff_map.begin(); giter != training_aff_map.end(); giter++){
    	    	vector<string> conv = giter->first;
    	    	int trainaff = training_aff_map[giter->first];
    	    	int trainunaff = training_unaff_map[giter->first];

    	    	if(trainaff == 0 && trainunaff == 0){
    	    		continue;
    	    	}

    	    	double unaffratio = 0;
    	    	double affratio = 0;

    	    	if(trainunaff > 0){
    	    		unaffratio = (double)((double) trainunaff / training_data_set->num_unaffected());
    	    	}
    	    	if(trainaff > 0){
    	    		affratio = (double)((double) trainaff / training_data_set->num_affected());
    	    	}

    	    	string rule = "";
    	    	for(int g = 0; g < (int)conv.size(); g++){
    	    		rule += conv[g] + "  ";
    	    	}
    	    	if(Helpers::dGreater(affratio, unaffratio)){
    	    		rule += "CASE";
    	    		mapping[conv] = 2;
    	    	}
    	    	else if(Helpers::dLess(affratio, unaffratio)){
    	    		rule += "CONTROL";
    	    		mapping[conv] = 1;
    	    	}
    	    	else{
    	    		rule += "CNC";
    	    		mapping[conv] = 0;
    	    	}
    	    	rules.push_back(rule);
    	    }

////    	    for(int s = 0; s < testing_data_set->num_inds(); s++){
////    	    	Sample* samp = testing_data_set->get_sample(s);

////    	    	vector<unsigned int> genos;

////    	    }

    		mdrout << mark->toString();

    		mdrout << "\t" << mdr_train.getBalancedAcc();
    		mdrout << "\t" << mdr_test.getBalancedAcc();
    		mdrout << "\t" << mdr_train.getFalseNegatives();
    		mdrout << "\t" << mdr_test.getFalseNegatives();
    		mdrout << "\t" << mdr_train.getFalsePositives();
    		mdrout << "\t" << mdr_test.getFalsePositives();
    		mdrout << "\t" << mdr_train.getTrueNegatives();
    		mdrout << "\t" << mdr_test.getTrueNegatives();
    		mdrout << "\t" << mdr_train.getTruePositives();
    		mdrout << "\t" << mdr_test.getTruePositives();
    		mdrout << endl;

    	}
    }

    if(training_data_set != NULL){
    	delete training_data_set;
    }
    if(testing_data_set != NULL){
    	delete testing_data_set;
    }

    if(mdrout.is_open()){
    	mdrout.close();
    }
}
