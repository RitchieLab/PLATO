/**********************************************************************************
*                       Marker Genotype Efficiency Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Example modlue framework
*
* Files generated:
*
*File: ExampleModule.cc
**********************************************************************************/


#include <iostream>
#include <Options.h>

#include "ExampleModule.h" //////CHANGE TO REAL MODULE NAME

using std::cout;
using Methods::opts;

const string ExampleModule::stepname = ExampleModule::doRegiser("example");

//PrintSummary()
//used to output results after processing the data
void ExampleModule::PrintSummary(){
	//output for process goes here
	//
	cout << "Hi, I'm in the ExampleModule PrintSummary method\n";
}

//filter()
//used to filter markers, samples, families after processing the data
void ExampleModule::filter(){
	//if markers, samples, or families need to be filtered, do it here (ie: use "setEnabled(flase);" method in classes)
	//
	cout << "Hi, I'm in the ExampleModule filter method\n";
}

//process()
//main method to get the process going and to the work
void ExampleModule::process(DataSet* ds){
	data_set = ds;
	//main processing of step goes here.
	cout << "Hi, I'm in the ExampleModule process method\n";
	cout << "I'm going to process " << data_set->num_loci() << " markers, " << data_set->num_pedigrees() << " families, and " << data_set->num_inds() << " samples!\n";

}

