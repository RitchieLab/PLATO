//Filter.cpp

#include "Filter.h"
namespace Filters{
///
/// Constructor
/// @param filterName name of this filter
///
Filter::Filter(string filterName){
  name = filterName;
  filterID = "Not Set";
  sort_priority = lowest;
}


// destructor clears any subfilters
Filter::~Filter(){
  for(unsigned int currSub=0; currSub < subfilters.size(); currSub++){
    delete subfilters[currSub];
  }
}


///
/// Copy Constructor
/// @param origFilter original Filter object
///
Filter::Filter(const Filter & origFilter){
  copy(origFilter);
}


///
/// Copies 
/// @param origFilter Filter to copy
/// @return none
/// 
void Filter::copy(const Filter & origFilter){
  for(unsigned int subIndex=0; subIndex < origFilter.subfilters.size();
    subIndex++){
      Filter * newSubFilter = FilterFactory::create_filter(origFilter.subfilters[subIndex]->getName());
      *newSubFilter = *origFilter.subfilters[subIndex];
      subfilters.push_back(newSubFilter);
  }
  filterID = origFilter.filterID;
  name = origFilter.name;
}


///
/// operator= overloaded for Filter object
/// @param origFilter filter to copy
/// @return this object
///
Filter & Filter::operator=(const Filter & origFilter){
  
  unsigned int subIndex;
  // clear any existing Filters
  for(subIndex=0; subIndex < subfilters.size(); subIndex++){
    delete subfilters[subIndex];
  }
  copy(origFilter);
  
  return *this;
}

///
/// Equivalent to NORMDIST function in excel
/// @param xx
/// @param mean
/// @param standard_dev
/// @return double fraction that will be < xx
///
double Filter::normdist(double xx, double mean, double standard_dev)
{
    double res;
    double x=(xx - mean) / standard_dev;
    if (x == 0)
    {
        res=0.5;
    }
    else
    {
        double oor2pi = 1/(sqrt(double(2) * 3.14159265358979323846));
        double t = 1 / (double(1) + 0.2316419 * fabs(x));
        t *= oor2pi * exp(-0.5 * x * x) 
             * (0.31938153   + t 
             * (-0.356563782 + t
             * (1.781477937  + t 
             * (-1.821255978 + t * 1.330274429))));
        if (x >= 0)
        {
            res = double(1) - t;
        }
        else
        {
            res = t;
        }
    }
    return res;
}


///
/// Returns number of seconds difference between 
/// 2 timespec time structures
/// @param x
/// @param y
/// @return time difference in seconds
///
// double Filter::timespecdiff(timespec* x, timespec* y){
// 
//   timespec result;
//   /* Perform the carry for the later subtraction by updating y. */
//   if (x->tv_nsec < y->tv_nsec) {
//     int nsec = (y->tv_nsec - x->tv_nsec) / 1000000000 + 1;
//     y->tv_nsec -= 1000000000 * nsec;
//     y->tv_sec += nsec;
//   }
//   if (x->tv_nsec - y->tv_nsec > 1000000000) {
//     int nsec = (x->tv_nsec - y->tv_nsec) / 1000000000;
//     y->tv_nsec += 1000000000 * nsec;
//     y->tv_sec -= nsec;
//   }
//      
//   /* Compute the time remaining to wait.
//       tv_nsec is certainly positive. */
//   result.tv_sec = x->tv_sec - y->tv_sec;
//   result.tv_nsec = x->tv_nsec - y->tv_nsec;
//   
//   return result.tv_sec + (double(result.tv_nsec) / 1000000000);
//   
//   /* Return 1 if result is negative. */
// //   return x->tv_sec < y->tv_sec;
// 
// }

///
/// Adds subfilter to current filter that will be used as needed
/// @param subfilter additional sub filter for this filter
/// @return none
///
void Filter::add_subfilter(Filter * subfilter){
  subfilters.push_back(subfilter);
}
}

