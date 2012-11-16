//ChiSquareArmitage.cpp

#include "ChiSquareArmitage.h"
namespace Methods{
/// Calculate chi square <br>
/// Assumes binary status with each individual being either a 0 or 1
/// @param chiTotals two dimensional vector with each row being a status and
///        each column a genotype.  Values represent totals for
///        that status and genotype
/// Ret: chi square
///
double ChiSquareArmitage::armitage(vector<vector<int> > & chiTotals){
  int affected=0, unaffected=0; // total row counts
  //double total_pop = cases + controls;
  vector<double> totalInCell; // total column counts
  vector<vector<double> > freqsInCell;
  freqsInCell.push_back(vector<double>(0));
  freqsInCell.push_back(vector<double>(0));
  for(unsigned int cell=0; cell<chiTotals[0].size(); cell++){
    unaffected+=chiTotals[0][cell];
    affected+=chiTotals[1][cell];
    totalInCell.push_back(chiTotals[0][cell] + chiTotals[1][cell]);
  }

  double totalinds = (double)(unaffected + affected);

  for(unsigned int cell = 0; cell < chiTotals[0].size(); cell++){
	  freqsInCell[0].push_back((chiTotals[0][cell]/(double)unaffected));
	  freqsInCell[1].push_back((chiTotals[1][cell]/(double)affected));
  }
  
  double chiscore=0.0;

  //slager scheid
  double numer = (((double)unaffected / (double)totalinds * (double)chiTotals[1][1]) - ((double)affected / (double)totalinds * (double)chiTotals[0][1])) + 2*(((double)unaffected / (double)totalinds * (double)chiTotals[1][2]) - ((double)affected / (double)totalinds * (double)chiTotals[0][2]));
  double denom = (double)affected * (double)unaffected * ( ((double)totalinds * ((double)totalInCell[1] + 4*(double)totalInCell[2]) - (((double)totalInCell[1] + 2*(double)totalInCell[2]) * ((double)totalInCell[1] + 2*(double)totalInCell[2]))) / ((double)totalinds * (double)totalinds * (double)totalinds));

  chiscore = (numer * numer) / denom;
  
  /*  English way
  double numer = (freqsInCell[1][1] * 0.5 + freqsInCell[1][2]) - (freqsInCell[0][1]*0.5 + freqsInCell[0][2]);
  numer = numer * numer;
  double denom = (1/(double)cases + 1/(double)controls) * (1/((double) total_pop * (double) total_pop)) * ((total_pop * (0.25 * totalInCell[1] + totalInCell[2])) - ((0.5 * totalInCell[1]+ totalInCell[2]) * (0.5 * totalInCell[1] + totalInCell[2])));

  cout << "numer = " << numer << " denom = " << denom << " un = " << unaffected << " aff = " << affected << " caf0 = " << freqsInCell[1][0] << " caf1 = " << freqsInCell[1][1] << " caf2 = " << freqsInCell[1][2] << " conf0 = " << freqsInCell[0][0] << " conf1 = " << freqsInCell[0][1] << " conf2 = " << freqsInCell[0][2] << endl;
  
  chiscore = numer / denom;
  */
  
  // need to take square root -- adjusted on account of Ben's review of formulas 3/13/09 -- smd
  // returned to original so that score can be used to get p value against standard distribution 4/30/09 -- smd
  return chiscore;
}
}
