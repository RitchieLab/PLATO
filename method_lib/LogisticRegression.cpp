//LogisticRegression.cpp

#include "LogisticRegression.h"
#include "ComboGenerator.h"
#include "Helpers.h"
namespace Methods{
///
/// Constructor
///
LogisticRegression::LogisticRegression(){
  set = NULL;
  initialize();
}


///
/// Alternative constructor
///
LogisticRegression::LogisticRegression(DataSet* ds){
  set = ds;
  initialize();
}

///
/// Initialize starting variables
///
void LogisticRegression::initialize(){
  fullInteraction = true;
  modelSize = 0;
  maxLocusValue = 2;
  missingValue = maxLocusValue+1;
  LociComboLimit = 10;
  LociComboMin = 1;
  maxIterations = 20;
  modType = Additive;
  defaultComboInterval = 10000;
  includeInteractions = true;
  coeff_intercept = 0.0;
  missingCoValue = -99999;

  ModelTypeMap["DOMINANT"] = Dominant;
  ModelTypeMap["RECESSIVE"] = Recessive;
  ModelTypeMap["ADDITIVE"] = Additive;

  PiD2 = 3.141592653589793 / 2;
  set_model();

  initialize_interactions();
}


///
/// Generates and stores interaction lists for producing interaction
/// terms when doing multi-locus models
///
void LogisticRegression::initialize_interactions(){
  vector<vector<unsigned int> > blank_vector;

  interaction_lists.clear();

  // indexes 0 and 1 contain empty vectors as placeholders
  // so that index of the interaction list will match
  // number of loci in the model
  interaction_lists.push_back(blank_vector);
  interaction_lists.push_back(blank_vector);

  unsigned int comboIndex;
  for(unsigned int curr_num=2; curr_num <= LociComboLimit; curr_num++){
    interaction_lists.push_back(blank_vector);

    ComboGenerator generator;
    generator.ComboEnds(2, curr_num);
    generator.SetLoci(curr_num);
    generator.SetComboInterval(defaultComboInterval);
    bool done=true;

    do{
      done = generator.GenerateCombinations();
      for(comboIndex=0; comboIndex < generator.ComboList.size(); comboIndex++){
        interaction_lists[curr_num].push_back(generator.ComboList[comboIndex]);
      }
    }while(!done);

  }
}


///
/// Runs logistic regression on dataset referenced in this method <br>
/// @param loci
/// @return none
///
void LogisticRegression::calculate(vector<unsigned int>& loci){
  // determine included indexes to use in running this dataset
//   int numGenosPerLocus = set->get_max_locus()+1;
//   missingValue = numGenosPerLocus;
//   if(set->missing_data_present())
//     numGenosPerLocus++;
//
//   indexConverter.set_genos_per_locus(numGenosPerLocus);
//   indexConverter.set_included_indexes(LociComboMin, LociComboLimit,
//     !set->missing_data_present(), missingValue);
//   includedIndexes = indexConverter.get_included_indexes();

  // when model is too large for current interaction list
  // update interactions
  if(loci.size() > LociComboLimit){
    LociComboLimit = loci.size();
    initialize_interactions();
  }

//  cout << "Point 1\n";
  vector<int> includedCells = includedIndexes[loci.size()];
//  cout << "Point 2\n";
  vector<unsigned int> converted_loci = convert_loc_order(loci);
//  cout << "Point 3\n";

  summarize_data(loci);//converted_loci);
//  cout << "Point 4\n";

  // assume loci are in marker_map order so need to alter to order contained
  // in samples
//   vector<unsigned int>::iterator iter;
//
//   for(iter=loci.begin(); iter!=loci.end(); iter++){
//     *iter = (*markers)[*iter]->getLoc();
//   }

  // run routine that calculates LR values
  calculateLR(summary_data, true, includedCells);
 // cout << "Point 5\n";

}


///
/// Reset data and set markers vector pointer
/// @param ds DataSet pointer
///
void LogisticRegression::resetDataSet(DataSet* ds){
  set = ds;
  markers = set->get_markers();
  missingValue = set->get_missing_value();
  missingCoValue = set->get_missing_covalue();

  int numGenosPerLocus = set->get_max_locus()+1;
  missingValue = numGenosPerLocus;
  if(set->missing_data_present())
    numGenosPerLocus++;

  indexConverter.set_genos_per_locus(numGenosPerLocus);
  indexConverter.set_included_indexes(LociComboMin, LociComboLimit,
    !set->missing_data_present(), missingValue);
  includedIndexes = indexConverter.get_included_indexes();

}


///
/// Runs logistic regression on dataset passed <br>
/// @param resultList ResultSet that may contain results from other filters
/// and will contain results of this sequential replication analysis
/// @param dataset DataSet to be analyzed
/// @return none
///
void LogisticRegression::calculate(vector<unsigned int>& loci, DataSet* dataset){
  resetDataSet(dataset);
  calculate(loci);
}


void LogisticRegression::calculateLR(vector<vector<double> >& data, bool summary_data, vector<int>& includedCells){


   unsigned int nRows = includedCells.size();
  // number of coefficients is equal to the main effects + the interactions
  unsigned int nColumns = data[0].size()-2;
  if(!summary_data)
    nColumns++;

  unsigned int nP = nColumns + 1;
  unsigned int nP1 = nP + 1;
  unsigned int sY0 = 0;
  unsigned int sY1 = 0;
  unsigned int sC = 0;
  double x;
  double v, xij, s,q;

  unsigned int unaffIndex = data[0].size()-2;
  unsigned int affIndex = unaffIndex+1;

  vector<double> X(nRows * ( nColumns + 1 ),0);
  vector<unsigned int> Y0(nRows,0);
  vector<unsigned int> Y1(nRows,0);
  vector<double> xM(nColumns+1,0.0);
  vector<double> xSD(nColumns+1, 0.0);
  //vector<double> coefficients(nP);
  vector<double> Par(nP);
  coefficients.clear();
  standard_errors.clear();
  coefficients.resize(nP-1);
  standard_errors.resize(nP-1);

  vector<double> SEP(nP);
  vector<double> Arr(nP * nP1);

  unsigned int i,j,k;

  // This loop stores the original values in X and adds values into
  // xM and xSD for mean and standard deviation calculations
  for(i=0; i<nRows; i++){
    X[ix(i,0,nColumns+1)] = 1;
    // store predictor values
    for(j=1; j<=nColumns; j++){
      X[ix(i,j,nColumns+1)] = data[includedCells[i]][j-1];
    }

    // status handled differently if data is summary
    if(summary_data)
    {
      Y0[i] = (unsigned int)(data[includedCells[i]][unaffIndex]);
      sY0 += Y0[i];
      Y1[i] = (unsigned int)(data[includedCells[i]][affIndex]);
      sY1 += Y1[i];
    }
    else
    { // not summary
      if(data[includedCells[i]].back() == 0)
      {
        Y0[i] = 1; sY0++;
      }
      else
      {
        Y1[i] = 1; sY1++;
      }
    }

    sC += Y0[i] + Y1[i];

    for(j=1; j<= nColumns; j++)
    {
      x = X[ix(i,j,nColumns+1)];
      xM[j] += (Y0[i] + Y1[i])*x;
		  xSD[j] += (Y0[i] + Y1[i])*x*x;
    }
  }

  // calculate mean and standard deviation
  for (j = 1; j<=nColumns; j++) {
    xM[j]  = xM[j]  / sC;
    xSD[j] = xSD[j] / sC;
    xSD[j] = sqrt( fabs( xSD[j] - xM[j] * xM[j] ) );
  }

  xM[0] = 0;
  xSD[0] = 1;

  // adjusts X values using the mean and standard deviation values
  for (i = 0; i<nRows; i++) {
    for (j = 1; j<=nColumns; j++) {
      X[ix(i,j,nColumns+1)] = ( X[ix(i,j,nColumns+1)] - xM[j] ) / xSD[j];
    }
  }

  Par[0] = log( double(sY1) / sY0 ); // use natural log of the ratio
  for (j = 1; j<=nColumns; j++) { // zero out all the others
    Par[j] = 0;
  }

  double LnV=0,Ln1mV=0, LLn=0;
  double LLp = 2e+10; // stores previous value of LL to check for convergence
  double LL = 1e+10;
  unsigned int numIterations = 0;

  while( fabs(LLp-LL)>0.0000001 ) {
    if(++numIterations > maxIterations){
      break;
    }

    LLp = LL;
    LL = 0;

    // zero out Arr for this iteration
    for (j = 0; j<=nColumns; j++) {
      for (k = j; k<=nColumns+1; k++) {
        Arr[ix(j,k,nColumns+2)] = 0;
      }
    }

	  // add to LL for each row
    for (i = 0; i<nRows; i++) {
      v = Par[0]; // for this ind start with Par[0] value
      for (j = 1; j<=nColumns; j++) {
        v = v + Par[j] * X[ix(i,j,nColumns+1)]; // add the value of Par for column and the value at that column
      }
      if( v>15 ) { LnV = -exp(-v); Ln1mV = -v; q = exp(-v); v=exp(LnV); }
      else {
        if( v<-15 ) { LnV = v; Ln1mV = -exp(v); q = exp(v); v=exp(LnV); }
        else { v = 1 / ( 1 + exp(-v) ); LnV = log(v); Ln1mV = log(1-v); q = v*(1-v); }
      }
      // calculate LL for this ind and add to running total
      LL = LL - 2*Y1[i]*LnV - 2*Y0[i]*Ln1mV;
      for (j = 0; j<=nColumns; j++) {
        xij = X[ix(i,j,nColumns+1)];
        Arr[ix(j,nColumns+1,nColumns+2)] = Arr[ix(j,nColumns+1,nColumns+2)] + xij * ( Y1[i] * (1 - v) + Y0[i] * (-v) );
        for (k=j; k<=nColumns; k++) {
          Arr[ix(j,k,nColumns+2)] = Arr[ix(j,k,nColumns+2)] + xij * X[ix(i,k,nColumns+1)] * q * (Y0[i] + Y1[i]);
        }
      }
    }

	  // when this is the first iteration, set LLn (null model) to be the current value of LL
    if( LLp==1e+10 ) { LLn = LL;}

    for (j = 1; j<=nColumns; j++) {
      for (k=0; k<j; k++) {
        Arr[ix(j,k,nColumns+2)] = Arr[ix(k,j,nColumns+2)];
      }
    }

    for (i=0; i<=nColumns; i++) {
      s = Arr[ix(i,i,nColumns+2)];
      Arr[ix(i,i,nColumns+2)] = 1;
      for (k=0; k<=nColumns+1; k++) {
        Arr[ix(i,k,nColumns+2)] = Arr[ix(i,k,nColumns+2)] / s;
      }
      for (j=0; j<=nColumns; j++) {
        if (i!=j) {
          s = Arr[ix(j,i,nColumns+2)]; Arr[ix(j,i,nColumns+2)] = 0;
          for (k=0; k<=nColumns+1; k++) {
            Arr[ix(j,k,nColumns+2)] = Arr[ix(j,k,nColumns+2)] - s * Arr[ix(i,k,nColumns+2)];
          }
        }
      }
    }
    for( j=0; j<=nColumns; j++) {
      Par[j] = Par[j] + Arr[ix(j,nColumns+1,nColumns+2)];
    }
  } // complete iteration

  // calculate p values for the coefficients
  // interaction coefficient for all loci is the last one
  for(j=1; j<=nColumns; j++) {
    Par[j] = Par[j] / xSD[j];
    Par[0] = Par[0] - Par[j] * xM[j];
  }

  if(isnan(LL)){
    overallPvalue = coeffPvalue = 1.0;
    LLR = 0.0;
  }
  else{
    // calculate coefficient p value
    for(j=1; j<=nColumns; j++){
      SEP[j] = sqrt( Arr[ix(j,j,nP+1)] ) / xSD[j];
    }
    j=nColumns;
    coeffPvalue = norm(fabs(Par[j]/SEP[j]));
    // calculate overall p value
    overallPvalue = ChiSq(fabs(LLn-LL), nColumns);

    coeff_intercept = Par[0];
    // adjust coefficients so that the zero index is now first coefficient
    for(j=1; j<=nColumns; j++) {
      coefficients[j-1] = Par[j];
      standard_errors[j-1] = SEP[j];
    }

    // calculate
    LLR = LLn-LL;
  }

}



///
/// Returns p value based on chi square scored passed and
/// degrees of freedom
/// @param x
/// @param n
/// @returns p value
///
double LogisticRegression::ChiSq(double x, unsigned int n) {

  if(x > 1000 || n>1000){
    double q = norm((pow(x/n,1/3)+2/(9*n)-1)/sqrt(2/(9*n)))/2;
    if(isnan(q)){
      return 0.0;
    }
    if(x>n){
      return q;
    }
    else{
      return 1-q;
    }
  }

  double p = exp(-0.5*x);
  if((n%2)==1) {
    p=p*sqrt(2*x/3.141592653589793);
  }
  unsigned int k=n;
  while(k>=2){
     p=p*x/k;
     k=k-2;
  }
  double t = p;
  unsigned int a=n;
  while(t>1e-15*p){
    a=a+2;
    t=t*x/a;
    p=p+t;
  }
  return 1-p;
}

///
/// Returns encoding for a genotype based on the model set.
/// @param geno Genotype value being checked
/// @param referent_allele Index of the allele that is the referent for the marker
/// @return recoded value for the genotype
///
int LogisticRegression::get_geno_conversion(int geno, int referent_allele){
	return geno_convert[referent_allele][geno];
}


///
/// Sets values to use for each genotype based on model selected <br>
/// Only works for SNPs currently.
///
void LogisticRegression::set_model(){
  // set up 2-D array
  geno_convert.assign(2, vector<unsigned int>(4,0));

//   geno_convert_zero.assign(4,0);
//   geno_convert_one.assign(4,0);

  switch(modType){
    case Dominant:
      geno_convert[0][0] = 1;
      geno_convert[0][1] = 1;
      geno_convert[0][3] = 2;
      geno_convert[1][1] = 1;
      geno_convert[1][2] = 1;
      geno_convert[1][3] = 2;
      maxLocusValue = 1;
      break;
    case Recessive:
//geno_convert[0] = 1;
      geno_convert[0][0] = 1;
      geno_convert[0][3] = 2;
      geno_convert[1][2] = 1;
      geno_convert[1][3] = 2;
      maxLocusValue = 1;
      break;
    case Additive:
      geno_convert[0][1] = 1;
      geno_convert[0][0] = 2;
      geno_convert[0][3] = 3;
      geno_convert[1][1] = 1;
      geno_convert[1][2] = 2;
      geno_convert[1][3] = 3;
      maxLocusValue = 2;
      break;
  }
}


///
/// Calculates p value
/// @param z
/// @return p value
///
double LogisticRegression::norm(double z){
  double q = z * z;
  if(fabs(z)>7){
    return (1-1/q+3/(q*q))*exp(-q/2)/(fabs(z)*sqrt(PiD2));
  }
  else{
    return ChiSq(q, 1);
  }
}


///
/// Fills summary vector with totals for use
/// in logistic regression routine
/// @param genos vector<unsigned int> containing loci to analyze
/// @return
///
void LogisticRegression::summarize_data(vector<unsigned int> & genos)
{
  unsigned int combSize = genos.size();

  // clear and initialize arrays
  // set model
  initialize_summary(combSize);

  unsigned int currInd, currLoc;

  // establish a vector of correct size that can be used to distribute individuals
  vector<int> genotype(combSize,0);

  // determine indexes for unaffected and affected totals
  unsigned int unaffIndex = summary_data[0].size()-2;
  unsigned int numInds = set->num_inds();

  vector<int> ref_alleles(genos.size(), 0);
  // determine which geno_conversion to use for each locus
  for(unsigned int curr_mark=0; curr_mark < genos.size(); curr_mark++)
  {
    ref_alleles[curr_mark] = set->get_locus(genos[curr_mark])->getReferentIndex();
  }

  // add to summary totals for each genotype
  // use current model to convert genotypes
  for(currInd=0; currInd < numInds; currInd++)
  {
	if(!(*set)[currInd]->isEnabled())
		continue;
    for(currLoc=0; currLoc < combSize; currLoc++)
    {
      // genotype[currLoc] = geno_convert[dataset[currInd][genos[currLoc]]];
      // genotype[currLoc] = geno_convert[(*set)[currInd]->get_genotype(genos[currLoc])];
      genotype[currLoc] = geno_convert[ref_alleles[currLoc]][(*set)[currInd]->get_genotype(set->get_locus(genos[currLoc])->getLoc())];
    }

    // increment count based on status of individual
    // summary_data[indexConverter.flatten_indexes(genotype)][unaffIndex+dataset[currInd].status]++;
    summary_data[indexConverter.flatten_indexes(genotype)][unaffIndex+(*set)[currInd]->getAffected()]++;
  }

}

///
/// Sizes and clears arrays for holding data
/// @param currModelSize size of model
/// @return
///
void LogisticRegression::initialize_summary(unsigned int currModelSize){

  unsigned int array_size = indexConverter.get_size_array(currModelSize);

  // check to see if need to resize all parameters
  if(currModelSize != modelSize){
    modelSize = currModelSize;

    // establish vector for holding individual totals
    summary_data.assign(array_size, vector<double>(currModelSize,0));
    // set the values for the predictor variables
      zero_summary(array_size, currModelSize);
  }
  else{
    // zero out totals for summary data
      zero_summary(array_size, currModelSize);
  }

}

///
/// Returns index into array
/// @param j
/// @param k
/// @param nCols
/// @returns index
///
unsigned int LogisticRegression::ix(int j,int k,int nCols){
  return j * nCols + k;
}


///
///  Sets genotype for single array, adds interaction terms if
///  any needed and requested by user, zeroes out totals
///  @param array_size unsigned int number of rows in summary
///  @param currModelSize unsigned int current model size to work with
///  @return
///
void LogisticRegression::zero_summary(unsigned int array_size, unsigned int model_size){

  double product=1;
  for(unsigned int sub_array_index=0; sub_array_index < array_size; sub_array_index++){
    summary_data[sub_array_index] = indexConverter.decode_index_double(sub_array_index, model_size);

    // this version can handle any size model
    if(summary_data[sub_array_index].size() > 1 && includeInteractions){
      for(int i=int(interaction_lists[model_size].size())-1; i >= 0; i--){
        product=1;
        for(int j=int(interaction_lists[model_size][i].size())-1; j >= 0; j--){
          product*= summary_data[sub_array_index][interaction_lists[model_size][i][j]];
        }
        summary_data[sub_array_index].push_back(product);
      }

      // when not full interaction remove the full interaction term from the list
      if(!fullInteraction){
        summary_data[sub_array_index].pop_back();
      }

    }
    summary_data[sub_array_index].push_back(0);
    summary_data[sub_array_index].push_back(0);
  }

}


///
/// Performs logistic regression on the snps and covariates
/// specified.  No interaction terms are included
/// @param loci vector pass empty if no markers included in analysis
/// @param covars vector pass empty if no covariates included in analysis
/// @param traits pass empty if no traits in analysis
///
void LogisticRegression::calculate(vector<unsigned int>& loci, vector<unsigned int>& covars, vector<unsigned int> & traits)
{

  unsigned int numLoci = loci.size();
  unsigned int numCovars = covars.size();
  unsigned int numTraits = traits.size();

  // determine size of row for each sample in dataset
  unsigned int row_size = numLoci + numCovars + numTraits + 1;

  vector<double> row(row_size, 0);

  // convert loci indexes to marker map
  vector<unsigned int> converted_loci = convert_loc_order(loci);

  // establish summary_data vector
  summary_data.clear();
  unsigned int numInds = set->num_inds();
  vector<int> includedCells;

  vector<int> ref_alleles(converted_loci.size(), 0);
  // determine which geno_conversion to use for each locus
  for(unsigned int curr_mark=0; curr_mark < converted_loci.size(); curr_mark++)
  {
    ref_alleles[curr_mark] = set->get_locus(loci[curr_mark])->getReferentIndex();
  }

  unsigned int currInd, currValue, i;
  for(currInd=0; currInd < numInds; currInd++)
  {
	  if(!(*set)[currInd]->isEnabled())
		  continue;
    currValue = 0;
    bool any_missing = false;
    for(i=0; i < numLoci; i++)
    {
      if((*set)[currInd]->get_genotype(converted_loci[i]) != missingValue)
      {
        row[currValue++] = geno_convert[ref_alleles[i]][(*set)[currInd]->get_genotype(converted_loci[i])];
      }
      else
      {
        any_missing = true;
        break;
      }
    }
    for(i=0; i < numCovars; i++)
    {
      //cout << "Adding " << getString<int>(numCovars) << " Covariates for Logistic Regression";
      if((*set)[currInd]->getCovariate(covars[i]) != missingCoValue)
        row[currValue++] = (*set)[currInd]->getCovariate(covars[i]);
      else
      {
        any_missing = true;
        break;
      }
    }
    for(i=0; i < numTraits; i++)
    {
      if((*set)[currInd]->getTrait(traits[i]) != missingCoValue)
        row[currValue++] = (*set)[currInd]->getTrait(traits[i]);
      else
      {
        any_missing = true;
        break;
      }
    }
    if(!any_missing)
    {
      row[currValue] = ((*set)[currInd]->getAffected());
      summary_data.push_back(row);
      includedCells.push_back(summary_data.size()-1);
    }
  }

  calculateLR(summary_data, false, includedCells);

}

///
/// Indexes refer to the order of the markers in map location.
/// To use in sample need to convert them to the indexes in the Samples.
/// @param loci vector of unsigned int that will be converted
/// @return returns vector with updated indexes
///
vector<unsigned int> LogisticRegression::convert_loc_order(vector<unsigned int>& loci){
  // assume loci are in marker_map order so need to alter to order contained
  // in samples
  vector<unsigned int>::iterator iter;
  vector<unsigned int> converted_indexes;

  for(iter=loci.begin(); iter!=loci.end(); iter++){
    converted_indexes.push_back((*markers)[*iter]->getLoc());
  }

  return converted_indexes;
}

vector<unsigned int> LogisticRegression::convert_loc_order(vector<Marker*> loci){
	vector<unsigned int> converted_indexes;
	for(unsigned int i = 0; i < loci.size(); i++){
		converted_indexes.push_back(loci[i]->getLoc());
	}
	return converted_indexes;
}


///
/// Sets the parameters using StepOptions class
/// @param options StepOptions containing options
///
void LogisticRegression::set_parameters(StepOptions* options){

    setFullInteraction(options->getLRFullInteraction());
    setIncludeInteractions(options->getLRIncludeInteractions());
    setMaximumIterations(options->getLRMaximumIterations());
    setModelType(options->getLRModelType());
}


///
/// Sets the model type to use in the calculation
/// @param modelType string containing model type to run
/// @throws Exception on error
///
void LogisticRegression::setModelType(string modelType){
  if(ModelTypeMap.find(modelType) != ModelTypeMap.end()){
    modType = ModelTypeMap[modelType];
  }
  else{
    throw MethodException(modelType + " is not a valid LR model type");
  }
  set_model();
}
}
