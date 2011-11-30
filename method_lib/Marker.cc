#include <stdio.h>
#include <math.h>
#ifndef MAC
#include <malloc.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <string>
#include <list>
#include "Globals.h"
#include "Marker.h"
#include "Helpers.h"

using namespace std;
namespace Methods{
string Marker::toString(){
	 return (getString<int>(chrom) + "\t" + rsid + "\t" + probe_id + "\t" + getString<int>(bploc) + getDetails());
}
}
