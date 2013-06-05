#include <stdio.h>
#include <math.h>
#include "config.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <map>
#include <vector>
#include <list>
#include "Options.h"
#include "Helpers.h"
#include "Globals.h"
#include "Sample.h"
using namespace std;
namespace Methods{
string Sample::getFamID(){
	if(opts::_TODIGIT_){
		return pfam->getFamID();
	}
	return famid;
}

string Sample::getInd(){
    if(opts::_TODIGIT_ && id_digit != -1){
	    return(getString<int>(id_digit));
    }
	return id;
}
string Sample::getMomID(){
	if(opts::_TODIGIT_ && mom_digit != -1){
		return getString<int>(mom_digit);
	}
	return mom;
}
string Sample::getDadID(){
	if(opts::_TODIGIT_ && dad_digit != -1){
		return getString<int>(dad_digit);
	}
	return dad;
}

string Sample::toString(){
	if(opts::_TODIGIT_ && id_digit != -1){
		return (getString<int>(pfam->getFamID_digit()) + "\t" + getString<int>(id_digit));
	}
	else{
		return (famid + "\t" + id);
	}
}

Sample* Sample::getLastChild(){
	int ssize = children.size();
	if(ssize == 0){
		return NULL;
	}
	return children.at(ssize - 1);
}

}
