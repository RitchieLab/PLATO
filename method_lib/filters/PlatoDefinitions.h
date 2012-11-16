//PlatoDefinitions.h

#ifndef __PLATODEFINITIONS_H__
#define __PLATODEFINITIONS_H__

#include <map>
#include <string>
#include <vector>
using namespace std;

namespace Filters{
/// Enumeration specifying the
enum ListOperation{
  // pass list on to next filter
  PassList,
  // removes designated combination sizes from list
  RemoveFromList,
  // add designated combianation sizes from list
  AddToList
};

/// Passes on information on list operations for the current filter
struct ListOptions{
  ListOperation op;
  vector<unsigned int> lociCombinations;
}; 
}

#endif
