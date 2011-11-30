/**********************************************************************************
*                       Superlink Input Module
*
*Written by: Justin Giles
*            Vanderbilt University
*            Center for Human Genetics Research
*
* Generates Superlink input files
*
*
*
*File: SuperlinkOutput.cc
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
#include <time.h>
#include "SuperlinkOutput.h"
#include "General.h"
#include "Helpers.h"
namespace Methods{
string SuperlinkOutput::stepname = "output-superlink";

void SuperlinkOutput::FilterSummary(){
}

void SuperlinkOutput::PrintSummary(){
	int msize = markers->size();

	for(int i = 0; i < msize; i++){
		(*markers).at(i)->setFlag(false);
	}

}

void SuperlinkOutput::filter(){
}

void SuperlinkOutput::process(vector<Sample*>* s, vector<Family*>* f, vector<Marker*>* m, vector<int>* mm)
{
	markers = m;
	families = f;
	samples = s;
	marker_map = mm;
	int msize = markers->size();
	int fsize = families->size();
	int numgoodmarkers = 0;
	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	numgoodmarkers = good_markers.size();
	msize = good_markers.size();
	int numfams = 0;
	bool alldigit = true;
	for(int f = 0; f < fsize; f++){
		Family* fam = (*families).at(f);
		if(Helpers::isAlphaNum(fam->getFamID())){
			alldigit = false;
		}
		else{
			fam->setAlphanumeric(false);
		}
		vector<Sample*>* fsamps = fam->getSamples();
		int fssize = fsamps->size();
		for(int s = 0; s < fssize; s++){
			Sample* samp = (*fsamps).at(s);
			if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())){
				numfams++;
				break;
			}
		}
	}

		opts::printLog("Remapping families to digit format.\n");
		Helpers::remapFamsToDigit(families);
		Helpers::printFamsToDigit(families, "input_superlink", options);
	string fname1 = opts::_OUTPREFIX_ + "input_superlink_locus" + options.getOut() + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname1 = options.getOverrideOut() + ".locus";
	}
	if(!overwrite){
		fname1 += "." + getString<int>(order);
	}
	filenames.push_back(fname1);
	ofstream locus (fname1.c_str());
	if(!locus.is_open()){
		opts::printLog("Unable to open " + fname1 + " for output!\n");
		throw MethodException("Unable to open " + fname1 + " for output!\n");
	}
	string fname2 = opts::_OUTPREFIX_ + "input_superlink_ped" + options.getOut() + ".txt";
	if(options.getOverrideOut().size() > 0){
		fname2 = options.getOverrideOut() + ".ped";
	}
	if(!overwrite){
		fname2 += "." + getString<int>(order);
	}
	filenames.push_back(fname2);
	ofstream ped (fname2.c_str());

    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
	string date = asctime(timeinfo);
	date = date.erase((date.size() - 1), 1);
	locus << (numgoodmarkers + 1) << " 0 0 5 0" << endl;
	//filler
	locus << "0 0.0 0.0 0" << endl;

	//print order of markers
	for(int i = 0; i < (numgoodmarkers + 1); i++){
		locus << (i+1) << " ";
	}
	locus << endl;

	locus << "###########\n#Insert affection status locus information here" << endl << "#Including penetrance information\n############" << endl;

	if (!opts::_FREQ_FILE_EXISTS_)
	{
		//the user did not specify a list of frequencies on the command line, so run AlleleFrequency to calculate them...
		AlleleFrequency* af = new AlleleFrequency(samples, families);
		af->setOptions(options);
		af->flagSamples();
		for(int m = 0; m < msize; m++){
			Marker* mark = good_markers.at(m);
			if(mark->isEnabled()){
				locus << "3 " << mark->getNumAlleles() << " #" << mark->getProbeID() << "#" << endl;
				af->calcOne(mark);
				if(mark->getNumAlleles() < 3){
					locus << af->getAone_freq() << " " << af->getAtwo_freq() << endl;
				}
				else{
					locus << af->getMicroFreq(0);
					for(int a = 1; a < mark->getNumAlleles(); a++){
						locus << " " << af->getMicroFreq(a);
					}
					locus << endl;
				}
			}
		}
		delete(af);
		locus.close();
	}
	else
	{
		//the user provided a list of frequencies, use the ones from the specified file...
		//the frequencies map contains snpid, frequency pairs
		//loop through all markers (snps), and print out frequency associated with each...
		for(int m = 0; m < msize; m++)
		{
			Marker* mark = good_markers.at(m);
			if(mark->isEnabled())
			{
				//this is a valid marker, needs printing...
				locus << "3 " << mark->getNumAlleles() << " #" << mark->getProbeID() << "#" << endl;
				locus << getString<float>(mark->getMAF()) << " #" << getString<float>(1- mark->getMAF());
			}
		}
	}

	for(int f = 0; f < fsize; f++)
	{
		Family* fam = (*families).at(f);
		if(fam->isEnabled() || (fam->isExcluded() && options.doIncExcludedSamples()) || (!fam->isEnabled() && options.doIncDisabledSamples()))
		{
			vector<Sample*>* fSamples = fam->getSamples();
			int fssize = fSamples->size();
			if(fSamples->size() == 0)
			{
				continue;
			}

			//need to instantiate the samplesUsed map to all false...
			for(int i = 0; i < fssize; i++)
			{
				samplesUsed[(*fSamples).at(i)->getInd_digit()] = false;
			}
			//find the TRUE founders so that we know that we are starting at the top of the pedigree
			vector<Sample*> trueFounders = FindTrueFounders(fSamples);

			//the PrintGenerationsToPed method recursively prints each generation and finds the next
			PrintGenerationsToPed(trueFounders, &ped);

			//now check to be sure that all samples have been used
			for(int s = 0; s < fssize; s++)
			{
				Sample* samp = (*fSamples).at(s);
				if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples()))
				{
					if (!samplesUsed.at(samp->getInd_digit()))
					{
						//this sample not used yet, print to PED file...
						PrintPedigreeInfo(samp, &ped);
					}
				}
			}
		}
	}
	ped.close();
}

//This method will loop through the sample list passed in, printing the results in the form of a Superlink Pedigree file
//It will then find the next generation of samples, and recurse.
void SuperlinkOutput::PrintGenerationsToPed(vector<Sample*> samps, ofstream* ped)
{
	if (samps.size() > 0)
	{
		for(unsigned int i = 0; i < samps.size(); i++)
		{
			//print the information for the current True Founder
			PrintPedigreeInfo((samps).at(i), ped);
			samplesUsed.at(samps.at(i)->getInd_digit()) = true;
			//mark the current sample as used
		}
			PrintGenerationsToPed(FindNextGeneration(samps), ped);
	}
}

//FindNextGeneration will find the next generation of samples given a current sample vector
//the next generation will include all children and partners of those children
vector<Sample*> SuperlinkOutput::FindNextGeneration(vector<Sample*> samps)
{
	//the vector being passed in represents a "generation" or level of a pedigree
	//this method will find the entire level below that one...
	int ssize = samps.size();
	vector<Sample*> children;
	vector<Sample*> grandchildren;
	vector<int> childrenIDs;

	for(int i = 0; i < ssize; i++)
	{
		vector<Sample*>* currChildren = (samps).at(i)->getChildren();
		//add all children to the children vector...

		for(unsigned int j = 0; j < (*currChildren).size(); j++)
		{
			Sample* currChild = (*currChildren).at(j);
			if(!samplesUsed.at(currChild->getInd_digit()))
			{
				children.push_back(currChild);
				(samplesUsed).at(currChild->getInd_digit()) = true;
				childrenIDs.push_back(currChild->getInd_digit());
			}
			vector<Sample*>* currGrandchildren = currChild->getChildren();
			Sample* currGrandchild;

			for(unsigned int g = 0; g < currGrandchildren->size(); g++)
			{
				grandchildren.push_back((*currGrandchildren).at(g));
				currGrandchild = (*currGrandchildren).at(g);
				if ((!samplesUsed.at(currGrandchild->getMom()->getInd_digit())) && currGrandchild->getMom() != NULL)
				{
					//mom needs to be added to the children vector...
					childrenIDs.push_back(currGrandchild->getMom()->getInd_digit());
					children.push_back(currGrandchild->getMom());
					(samplesUsed).at(currGrandchild->getMom()->getInd_digit()) = true;
				}
				if ((!samplesUsed.at(currGrandchild->getDad()->getInd_digit())) && currGrandchild->getDad() != NULL)
				{
					//mom needs to be added to the children vector...
					childrenIDs.push_back(currGrandchild->getDad()->getInd_digit());
					children.push_back(currGrandchild->getDad());
					(samplesUsed).at(currGrandchild->getDad()->getInd_digit()) = true;
				}
			}//end loop through grandchildren
		}//end loop through children
	}//end loop through parameter sample list
	return children;
}//end FindNextGeneration

//FindTrueFounders will find only those founders whose partners also have no parents
vector<Sample*> SuperlinkOutput::FindTrueFounders(vector<Sample*>* samps)
{
	//loop through all samples, pick out the TRUE founders (samples with no parents, and whose partners also have no parents)...
	int ssize = samps->size();
	vector<Sample*> trueFounders;

	for(int s = 0; s < ssize; s++)
	{
		Sample* samp = (*samps).at(s);
		if((samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples())) && samp->isFounder())
		{
			vector<Sample*>* children = samp->getChildren();
			bool isTrueFounder = true;
			for(unsigned int c = 0; c < children->size(); c++)
			{
				Sample* dad = (*children).at(c)->getDad();
				Sample* mom = (*children).at(c)->getMom();

				if ( (dad != NULL && !dad->isFounder()) || (mom != NULL && !mom->isFounder()))
				{
					//at least one parent has a sample and is NOT a founder, so the current founder is not a TRUE founder
					isTrueFounder = false;
					break;
				}
			}//end for children
			if(isTrueFounder)
			{
				//add to list of True Founders...
				(trueFounders).push_back(samp);
			}
		}
	}
	return trueFounders;
}

//PrintPedigreeInfo
void SuperlinkOutput::PrintPedigreeInfo(Sample* samp, ofstream* ped)
{
	vector<Marker*> good_markers = Helpers::findValidMarkers(markers, &options);
	int msize = good_markers.size();
	if(samp->isEnabled() || (samp->isExcluded() && options.doIncExcludedSamples()) || (!samp->isEnabled() && options.doIncDisabledSamples()))
	{
		(*ped) << samp->getFamily()->getFamID_digit() << "\t" << samp->getInd_digit() << " " << samp->getDadID_digit() << " " << samp->getMomID_digit() << " ";
		vector<Sample*>* children = samp->getChildren();
		if(children->size() == 0)
		{
			(*ped) << "0 ";
		}
		else
		{
			(*ped) << (*children).at(0)->getInd_digit() << " ";
		}
		//paternal sib
		Sample* sib = samp->getPatSib();
		if(sib == NULL)
		{
			(*ped) << "0 ";
		}
		else
		{
			(*ped) << sib->getInd_digit() << " ";
		}
		//maternal sib
		sib = samp->getMatSib();
		if(sib == NULL)
		{
			(*ped) << "0 ";
		}
		else
		{
			(*ped) << sib->getInd_digit() << " ";
		}

		//gender
		if(samp->getSex())
		{
			(*ped) << "1";
		}
		else
		{
			(*ped) << "2";
		}
		//filler
		(*ped) << " 0";
		//Disease status
		if(samp->getPheno() == 2)
		{
			(*ped) << " " << "2";
		}
		else if(samp->getPheno() == 1)
		{
			(*ped) << " " << "1";
		}
		else
		{
			(*ped) << " " << "0";
		}

		//penetrance info
		(*ped) << " " << options.findPenetranceCode(samp->getFamID() + " " + samp->getInd());
		for(int m = 0; m < msize; m++)
		{
			Marker* mark = good_markers.at(m);
			if(mark->isEnabled())
			{
				int mloc = mark->getLoc();
				if((samp->isExcluded() && options.doZeroExcluded()) || (!samp->isEnabled() && options.doZeroDisabled()))
				{
					(*ped) << " 0 0";
					continue;
				}
				if(!mark->isMicroSat())
				{
					if(!samp->getAone(mloc) && !samp->getAtwo(mloc))
					{
						(*ped) << " " << Helpers::map_allele(mark, mark->getAllele1(), &options);
						(*ped) << " " << Helpers::map_allele(mark, mark->getAllele1(), &options);
					}
					else if(!samp->getAone(mloc) && samp->getAtwo(mloc))
					{
						(*ped) << " " << Helpers::map_allele(mark, mark->getAllele1(), &options);
						(*ped) << " " << Helpers::map_allele(mark, mark->getAllele2(), &options);
					}
					else if(samp->getAone(mloc) && samp->getAtwo(mloc) && !samp->getAmissing(mloc))
					{
						(*ped) << " " << Helpers::map_allele(mark, mark->getAllele2(), &options);
						(*ped) << " " << Helpers::map_allele(mark, mark->getAllele2(), &options);
					}
					else
					{
						(*ped) << " 0 0";
					}
				}//end !microsat
				else
				{
					if(samp->getAbone(mloc) == -1)
					{
						(*ped) << " 0 0";
					}
					else
					{
						(*ped) << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbone(mloc)), &options);
						(*ped) << " " << Helpers::map_allele(mark, mark->getAllele(samp->getAbtwo(mloc)), &options);
					}
				}//end microsat
			}
		}
		(*ped) << endl;
	}
}//end PrintPedigreeInfo(Sample* samp)

int SuperlinkOutput::map_sex(char c){
	if(c == 'M'){
		return 1;
	}
	return 2;
}
}
