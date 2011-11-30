// ComboGenerator.cpp
#include<iostream>
using namespace std;
#include "ComboGenerator.h"

namespace Methods{
///
/// Constructor initializes parameters
///
ComboGenerator::ComboGenerator(){
  initialize();
}

///
/// Initializes variables for generator
/// @return none
///
void ComboGenerator::initialize(){
  ComboInterval = 1000;
  AlreadyStarted = false;
  NumLoci = 0;
}

///
/// Uses modified code by Donald Knuth and Glenn C. Rhoads to generate
/// all possible combinations from COMBOSTART to COMBOEND in intervals specified in COMBOINTERVAL
/// That is, N = COMBOEND, C = COMBOSTART. Summation of N choose C  -->  N choose N.
/// Code is confusing but fast
/// @return true if all combinations completed
///
bool ComboGenerator::GenerateCombinations(){
  // Clear the ComboList!!! Very important, or Memory will overflow
  ComboList.clear();

  counter = 0;  // counter for the number of combinations created

  if(AlreadyStarted)   // If combination generator has already been started
  {
    counter--;
    goto resume;    // Ugh.. a GOTO... but it works!
  }
  else{
    kdec = ComboEnd;
    
    c = new int[ComboEnd+3];

  }

  AlreadyStarted = true;  // If it wasn't already started, it is now

  init:

      if(kdec == (ComboStart-1))
      {
        delete [] c;
        AlreadyStarted = false;
        return true;
      }

      for (int i=1; i <= kdec; i++)
      {  c[i] = i;  }
      c[kdec+1] = NumLoci+1;
      c[kdec+2] = 0;
      j = kdec;
      
      if(NumLoci == ComboEnd){
        j=0;
      }
      
      // Add a new combination to the ComboList
      ComboList.push_back(std::vector <unsigned int>());

  visit:
      for (int i=kdec; i >= 1; i--)
      {
    // Add an element to the new combination
    ComboList.at(counter).push_back(c[i]-1);
      }


      if(counter >= ComboInterval)  // If you exceed the interval limit
      {
          return false;
      }

  resume:
            // Add a new combination to the ComboList
      ComboList.push_back(std::vector <unsigned int>());

      counter++;
      if (j > 0) {x = j+1; goto incr;}

      if (c[1] + 1 < c[2])
         {
         c[1] += 1;
         goto visit;
         }

      j = 2;

   do_more:
      c[j-1] = j-1;
      x = c[j] + 1;
      if (x == c[j+1]) {j++; goto do_more;}

      if (j > kdec)
      {
          kdec--;

    // Remove the last empty combination before quitting
          ComboList.pop_back();
          goto init;
      }

   incr:
      c[j] = x;
      j--;
      goto visit;
}


///
/// Advances parameters without adding any combinations
/// to the combination list
///
bool ComboGenerator::AdvanceParameters(){
  counter = 0;  // counter for the number of combinations created

  if(AlreadyStarted)   // If combination generator has already been started
  {
    counter--;
    goto resume;    // Ugh.. a GOTO... but it works!
  }
  else{
    kdec = ComboEnd;
    x = 0;
    j = kdec;
    
    if(NumLoci == ComboEnd){
      j=0;
    }
    
    c = new int[ComboEnd+3];

    for(int i=1; i <= kdec; i++)
    {
      c[i] = i;
    }
    c[kdec+1] = NumLoci+1;
    c[kdec+2] = 0;
  }

  AlreadyStarted = true;  // If it wasn't already started, it is now

  init:

      if(kdec == (ComboStart-1))
      {
        delete [] c;
        AlreadyStarted = false;
        --counter;
        return true;
      }

      for (int i=1; i <= kdec; i++)
      {  c[i] = i;  }
      c[kdec+1] = NumLoci+1;
      c[kdec+2] = 0;
      j = kdec;
      
      if(NumLoci == ComboEnd){
        j=0;
      }
      

  visit:


      if(counter >= ComboInterval)  // If you exceed the interval limit
      {
          return false;
      }

  resume:
            // Add a new combination to the ComboList

      counter++;
      if (j > 0) {x = j+1; goto incr;}

      if (c[1] + 1 < c[2])
         {
         c[1] += 1;
         goto visit;
         }

      j = 2;

   do_more:
      c[j-1] = j-1;
      x = c[j] + 1;
      if (x == c[j+1]) {j++; goto do_more;}

      if (j > kdec)
      {
          kdec--;

    // Remove the last empty combination before quitting
          goto init;
      }

   incr:
      c[j] = x;
      j--;
      goto visit;  
}


///
/// Creates initial state 
///
void ComboGenerator::initialize_state(){
  kdec = ComboEnd;
  x = 0;
  j = kdec;
  
  c = new int[ComboEnd+3];

  for(int i=1; i <= kdec; i++){
    c[i] = i;
  }
  c[kdec+1] = NumLoci+1;
  c[kdec+2] = 0;
 
}


///
/// Generates combinations in amount passed
/// @param new_interval Combinations to create
/// 
bool ComboGenerator::GenerateCombinations(int new_interval){
  int old_interval = GetComboInterval();
  SetComboInterval(new_interval);
  bool done = GenerateCombinations();
  SetComboInterval(old_interval);
  return done;
}


///
///Sets the interval for generating combinations.  Used so that
///don't excessively use memory 
//@param cmbInterval interval to use in producing combinations
//@return none
///
void ComboGenerator::SetComboInterval(int cmbInterval){
  ComboInterval = cmbInterval;
}


///
///Sets the start combination number and ending combination
///@param  combStart minimum size of combination
///@param  combEnd maximum size of combination
///@return  none
///
void ComboGenerator::ComboEnds(int combStart, int combEnd){
  ComboStart = combStart;
  ComboEnd = combEnd;
}

///
///Sets total number of loci for generator
///@param nLoci number of loci
///@return none
///
void ComboGenerator::SetLoci(int nLoci){
  NumLoci = nLoci;
}


///
/// Calculates number of combinations 
/// @param n Number of items
/// @param k Combination size
///
double ComboGenerator::calc_combinations(double n, unsigned int k){
  if(k > n)
    return 0.0;
    
  if(k > n/2)
    k = n-k;
    
  long double accum = 1;
  for(unsigned i=1; i<=k; i++)
    accum = accum * (n-k+i)/i;
  
  return accum;
    
}

}
