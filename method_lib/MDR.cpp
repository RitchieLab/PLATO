//MDR.cpp

#include "MDR.h"
namespace Methods{
///
/// Constructor
///
MDR::MDR(){
  dataset = NULL;
  initialize();
}


///
/// Alternative constructor
///
MDR::MDR(DataSet* ds){
  resetDataSet(ds);
  initialize();
}

///
/// Initialize starting variables
///
void MDR::initialize(){

  maxLocusValue = 2;
  missingValue = maxLocusValue+1;
  LociComboLimit = 5;
  LociComboMin = 1;
  tieCellValue = 1;
  set_threshold = 1.0;
  missingCoValue = -9999;

}


///
/// Reset data and set markers vector pointer
/// @param ds DataSet pointer
///
void MDR::resetDataSet(DataSet* ds){
  dataset = ds;
  calculate_set_threshold();
  missingValue = dataset->get_missing_value();
  markers = dataset->get_markers();
  missingCoValue = int(dataset->get_missing_covalue());
  setIndexConverter(dataset->get_max_locus(), dataset->missing_data_present(),
    missingValue);
}



///
/// Sets FlatIndex object to match current analysis
///
// void MDR::setIndexConverter(){
//   int numGenosPerLocus = dataset->get_max_locus()+1;
//   if(dataset->missing_data_present())
//     numGenosPerLocus++;
//   indexConverter.set_genos_per_locus(numGenosPerLocus);
//
//   indexConverter.set_included_indexes(LociComboMin, LociComboLimit,
//     !dataset->missing_data_present(), missingValue);
//   includedIndexes = indexConverter.get_included_indexes();
// }

///
/// Sets FlatIndex object to match current analysis
/// @param max_locus Maximum locus value for data
/// @param any_missing_data true when there are missing data in set
/// @param missing_val Value of missing data in set
///
void MDR::setIndexConverter(int max_locus, bool any_missing_data, int missing_val){

  int numGenosPerLocus = max_locus+1;
  if(any_missing_data)
    numGenosPerLocus++;

  indexConverter.set_genos_per_locus(numGenosPerLocus);

  indexConverter.set_included_indexes(LociComboMin, LociComboLimit,
    !any_missing_data, missing_val);

  includedIndexes = indexConverter.get_included_indexes();

}



///
/// Calculates overall threshold for determining high
/// and low risk cells in a model
///
void MDR::calculate_set_threshold(){

  vector<int> statusTotals(2,0);
  int numInds = dataset->num_inds();

  for(int currInd=0; currInd < numInds; currInd++)
//     statusTotals[dataset[currInd].status]++;
    statusTotals[dataset->get_sample(currInd)->getAffected()]++;

  set_threshold = 0.0;

  if(statusTotals[0] > 0)
    set_threshold = float(statusTotals[1])/statusTotals[0];
}


///
/// Calculates threshold for the indicated loci.  Useful when
/// missing data occurs in the set.
/// @param vector<unsigned int>
///
float MDR::calculate_model_threshold(vector<unsigned int> & lociComb){
  vector<int> statusTotals(2,0);
  int numInds = dataset->num_inds();
  unsigned int currLoc, numLoci = lociComb.size();
  bool missingPresent;

  for(int currInd=0; currInd < numInds; currInd++){
    missingPresent = false;
    for(currLoc=0; currLoc < numLoci; currLoc++)
      if(dataset->get_sample(currInd)->get_genotype(lociComb[currLoc])==missingValue){
        missingPresent=true;
        break;
      }

    if(!missingPresent)
      statusTotals[dataset->get_sample(currInd)->getAffected()]++;
  }

  float threshold=0;
  if(statusTotals[0] > 0)
    threshold = float(statusTotals[1])/statusTotals[0];

  return threshold;
}



///
/// calculate ratio of affected to unaffected for this model
/// @param loci vector pass empty if no markers included in analysis
/// @param covars vector pass empty if no covariates included in analysis
/// @param traits pass empty if no traits in analysis
/// @return threshold
///
float MDR::calculate_model_threshold(vector<unsigned int> loci,
  vector<unsigned int>& covars, vector<unsigned int> & traits){

  // check each individual at each snp, covariate or trait and
  // exclude from threshold calculation if individual is missing
  // at any of these
  vector<int> statusTotals(2,0);
  int numInds = dataset->num_inds();
  unsigned int curr, numLoci = loci.size(), numCovars = covars.size(), numTraits = traits.size();
  bool missingPresent;

  for(int currInd=0; currInd < numInds; currInd++){
    missingPresent = false;
    for(curr=0; curr < numLoci; curr++)
//       if(dataset[currInd][lociComb[currLoc]]==missingValue){
      if(dataset->get_sample(currInd)->get_genotype(loci[curr])==missingValue){
        missingPresent=true;
        break;
      }

    if(!missingPresent)
      for(curr=0; curr < numCovars; curr++)
        if(dataset->get_sample(currInd)->getCovariate(covars[curr])== missingCoValue){
          missingPresent = true;
          break;
        }

    if(!missingPresent)
      for(curr=0; curr < numTraits; curr++)
        if(dataset->get_sample(currInd)->getTrait(traits[curr])==missingCoValue){
          missingPresent = true;
          break;
        }

    if(!missingPresent)
      statusTotals[dataset->get_sample(currInd)->getAffected()]++;
  }

  float threshold=0;
  if(statusTotals[0] > 0)
    threshold = float(statusTotals[1])/statusTotals[0];

  return threshold;


}


///
/// Calculates MDR values
/// @param loci vector of loci that are in order of
///
void MDR::calculate(vector<unsigned int> loci){

  // assume loci are in marker_map order so need to alter to order contained
  // in samples
  convert_loci_indexes(loci);

  // make sure
  if(loci.size() > LociComboLimit){
    LociComboLimit = loci.size();
    setIndexConverter(dataset->get_max_locus(), dataset->missing_data_present(),
      missingValue);
  }

  distributeInds(loci);

  if(calc_each_threshold){
    mod.threshold = calculate_model_threshold(loci);
  }
  else{
    mod.threshold = set_threshold;
  }

  runMDR(loci.size());

}



///
/// Converts loci from marker_map order to order stored in
/// the samples
/// @param loci
///
void MDR::convert_loci_indexes(vector<unsigned int> & loci){
  vector<unsigned int>::iterator iter;
  for(iter=loci.begin(); iter!=loci.end(); iter++)
    *iter = (*markers)[*iter]->getLoc();
}


///
/// Checks traits and covariates passed and
/// sets up FlatIndex object to work with appropriately
/// sized vectors.
/// @param covars vector<unsigned int> Covariate indexes
/// @param traits vector<unsigned int> Trait indexes
/// @return highest value found in a covariate or trait
///
int MDR::set_converter_covars(vector<unsigned int> & covars,
  vector<unsigned int> & traits){

  int highest_value = dataset->get_max_locus();
  int num_inds = dataset->num_inds();

  bool missing_present = dataset->missing_data_present();

  for(unsigned int i=0; i<covars.size(); i++){
    for(int j=0; j<num_inds; j++){
      if((*dataset)[j]->getCovariate(covars[i]) > highest_value){
        highest_value = int((*dataset)[j]->getCovariate(covars[i]));
      }
      if((*dataset)[j]->getCovariate(covars[i]) == missingCoValue){
        missing_present = true;
      }
    }
  }

  for(unsigned int i=0; i<traits.size(); i++){
    for(int j=0; j<num_inds; j++){
      if((*dataset)[j]->getTrait(traits[i]) > highest_value){
        highest_value = int((*dataset)[j]->getTrait(traits[i]));
      }
      if((*dataset)[j]->getTrait(traits[i]) == missingCoValue){
        missing_present = true;
      }
    }
  }

  if(highest_value == missingCoValue && highest_value > (int)dataset->get_max_locus()){
    highest_value--;
  }

  setIndexConverter(highest_value, missing_present,
    highest_value + 1);

  return highest_value;

}


///
/// Calculates MDR on loci, covariates and traits.  When
/// only running analysis on loci, it is faster to use the
/// calculate method that only accepts a list of loci.
/// @param loci vector pass empty if no markers included in analysis
/// @param covars vector pass empty if no covariates included in analysis
/// @param traits pass empty if no traits in analysis
///
void MDR::calculate(vector<unsigned int> loci,
  vector<unsigned int>& covars, vector<unsigned int> & traits){

  // need to establish the size of the grid to use in calculating
  // MDR on this combination
  // Only necessary if any covars or traits in the analysis
  int high_value = dataset->get_max_locus();

  if(loci.size() > LociComboLimit)
    LociComboLimit = loci.size() + covars.size() + traits.size();

  if(covars.size() > 0 || traits.size() > 0){
    high_value = set_converter_covars(covars, traits);
  }


  // assume loci are in marker_map order so need to alter to order contained
  // in samples
  convert_loci_indexes(loci);

  // distribute individuals to create totals
  distributeInds(loci, covars, traits, high_value+1);

  // set the threshold based on whether threshold should be set based
  // or based on this individual model (when missing data in set)
  if(calc_each_threshold){
    mod.threshold = calculate_model_threshold(loci, covars, traits);
  }
  else{
    mod.threshold = set_threshold;
  }

  runMDR(loci.size() + covars.size() + traits.size());

}


///
/// Calculates MDR on model
/// @param model_size Size of model being analyzed
///
void MDR::runMDR(unsigned int model_size){

  calculateStats(model_size);

  mod.balanced_accuracy = .5 * (float(mod.classhigh)/(mod.classhigh + mod.misclasslow)
    + float(mod.classlow)/(mod.classlow + mod.misclasshigh));

}


///
/// Calculate corrrect and incorrectly classified individuals
///       for further use
/// @param combSize Number of loci being analyzed by model
///
void MDR::calculateStats(unsigned int combSize){

  vector<int> & includedCells = includedIndexes[combSize];
  unsigned int numCells = includedCells.size();

  float calculatedRisk;

  mod.classhigh = 0;
  mod.misclasshigh = 0;
  mod.classlow = 0;
  mod.misclasslow = 0;
  mod.totaltiecells = 0;
//cout << "numCells = " << numCells << " combSize = " << combSize << " includedindexes = " << includedIndexes.size() << endl;
  for(unsigned int currCell=0; currCell<numCells; currCell++){
    // check that unaffected total in cell is greater than zero
    if(mod.unaffected[includedCells[currCell]]>0)
      calculatedRisk = float(mod.affected[includedCells[currCell]])
        /mod.unaffected[includedCells[currCell]];
    else if(mod.affected[includedCells[currCell]] > 0)
      calculatedRisk = mod.threshold + 1;
    else
      continue; //skip to next cell
    int cellstatus;
    if(calculatedRisk > mod.threshold){
      cellstatus = 1;
    }
    else if(calculatedRisk < mod.threshold){
      cellstatus = 0;
    }
    else{
      cellstatus = tieCellValue;
    }

    switch(cellstatus){
      case 1:
        mod.classhigh += mod.affected[includedCells[currCell]];
        mod.misclasshigh += mod.unaffected[includedCells[currCell]];
        break;
      case 0:
        mod.classlow += mod.unaffected[includedCells[currCell]];
        mod.misclasslow += mod.affected[includedCells[currCell]];
        break;
      case -1:
        mod.totaltiecells+= mod.affected[includedCells[currCell]]
          + mod.unaffected[includedCells[currCell]];
      break;
    };
  }

}


///
/// Distributes individuals int totals in model
/// @param loci vector pass empty if no markers included in analysis
/// @param covars vector pass empty if no covariates included in analysis
/// @param traits pass empty if no traits in analysis
/// @param miss Value to use for any missing data in the set
///
void MDR::distributeInds(vector<unsigned int>& loci,
  vector<unsigned int>& covars, vector<unsigned int> & traits, int miss){

  unsigned int combSize = loci.size() + covars.size() + traits.size();
  unsigned int locSize = loci.size();
  unsigned int coSize = covars.size();
  unsigned int traitSize = traits.size();

  // establish a vector of correct size that can be used to distribute individuals
  vector<int> scores(combSize,0);

  int linearSize = indexConverter.get_size_array(combSize);

  mod.affected.assign(linearSize, 0);
  mod.unaffected.assign(linearSize, 0);

  unsigned int currVal, numInds = dataset->num_inds(), curr_score=0;

  for(unsigned currInd=0; currInd < numInds; currInd++){
    curr_score = 0;
    for(currVal=0; currVal < locSize; currVal++){
      if((*dataset)[currInd]->get_genotype(loci[currVal]) == missingValue){
        scores[curr_score++] = miss;
      }
      else
        scores[curr_score++] = (*dataset)[currInd]->get_genotype(loci[currVal]);
    }

    for(currVal=0; currVal < coSize; currVal++){
      if((*dataset)[currInd]->getCovariate(covars[currVal]) == missingCoValue){
        scores[curr_score++] = miss;
      }
      else
        scores[curr_score++] = int((*dataset)[currInd]->getCovariate(covars[currVal]));
    }

    for(currVal=0; currVal < traitSize; currVal++){
      if((*dataset)[currInd]->getTrait(traits[currVal]) == missingCoValue){
        scores[curr_score++] = miss;
      }
      else
        scores[curr_score++] = int((*dataset)[currInd]->getTrait(traits[currVal]));
    }

    // deterimine index and add to appropriate status totals
    if(dataset->get_sample(currInd)->getAffected())
      mod.affected[indexConverter.flatten_indexes(scores)]++;
    else
      mod.unaffected[indexConverter.flatten_indexes(scores)]++;

  }

}



///
/// Distributes individuals into totals in model
/// @param lociComb SNPs to analyze
///
void MDR::distributeInds(vector<unsigned int>& lociComb){

  unsigned int combSize = lociComb.size();
  // establish a vector of correct size that can be used to distribute individuals
  vector<int> genotype(combSize,0);

  int linearSize=indexConverter.get_size_array(lociComb.size());

  mod.affected.assign(linearSize,0);
  mod.unaffected.assign(linearSize,0);

  unsigned int currLoc, numInds = dataset->num_inds();
  // place each individual in proper location in model linear arrays
  // based on location and genoype
  for(unsigned int currInd=0; currInd < numInds; currInd++){
    // create genotype vector
    for(currLoc=0; currLoc < combSize;
      currLoc++){
      genotype[currLoc] = dataset->get_sample(currInd)->get_genotype(lociComb[currLoc]);
    }
    // deterimine index and add to appropriate status totals
    if(dataset->get_sample(currInd)->getAffected())
      mod.affected[indexConverter.flatten_indexes(genotype)]++;
    else
      mod.unaffected[indexConverter.flatten_indexes(genotype)]++;
  }

}

///
/// Sets the parameters using StepOptions class
/// @param options StepOptions containing options
///
void MDR::set_parameters(StepOptions* options){
  calc_each_threshold = !options->getOnlySetThreshold();
}
}
