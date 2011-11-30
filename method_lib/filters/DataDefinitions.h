//DataDefinitions.h

#ifndef datahandling__DATADEFINITIONS_H__
#define datahandling__DATADEFINITIONS_H__

#include <vector>

namespace Filters{

typedef std::vector<unsigned int> LocusVec;

// Size of rs number character arrays in binary files
#define RSNUMSIZE 20;

// Maximum size of combination of loci in a result
#define MAXIMUMLOCSIZE = 20;

}

#endif


