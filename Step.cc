#include <stdio.h>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include "Step.h"
//#include "Process.h"
using namespace Methods;
#ifdef PLATOLIB
using namespace PlatoLib;
#endif
string Step::getName(){
	return name;
}

void Step::setProcess(Process* p){
	myprocess = p;
}

void Step::setProcess(void* p){
	myprocess = (Process*) p;
}

string Step::getThreshold(){
	return threshold;
}

bool Step::getMultiThresh(){
	return multi_thresh;
}

void Step::setThreshold(string val){
	threshold = val;
	myprocess->setThreshold(threshold);
}

void Step::setOverwrite(bool val){
	myprocess->setOverwrite(val);
}

void Step::setOrder(int val){
	order = val;
	myprocess->setOrder(order);
}

bool Step::hasIncExc(){
	return myprocess->hasIncExc();
}
