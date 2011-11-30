//LogisticRegression.cpp

#include "LogisticRegression.h"
#include "ComboGenerator.h"
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

  vector<int> includedCells = includedIndexes[loci.size()];
  vector<unsigned int> converted_loci = convert_loc_order(loci);

  summarize_data(converted_loci);

  // assume loci are in marker_map order so need to alter to order contained
  // in samples
//   vector<unsigned int>::iterator iter;
//
//   for(iter=loci.begin(); iter!=loci.end(); iter++){
//     *iter = (*markers)[*iter]->getLoc();
//   }

  // run routine that calculates LR values
  calculateLR(summary_data, true, includedCells);

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


//output data
/*
for(size_t i=0; i < data.size();i++){
  cout << data[i][0];
  for(size_t j=1; j < data[i].size(); j++){
   cout << "," << data[i][j];
  }
cout << endl;
}
*/
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
//  vector<double> coefficients(nP);
  vector<double> Par(nP);
  coefficients.resize(nP-1);
  standard_errors.resize(nP-1);

  vector<double> SEP(nP);
  vector<double> Arr(nP * nP1);

  unsigned int i,j,k;

  // This loop stores the original values in X and adds values into
  // xM and xSD for mean and standard deviation calculations
  for(i=0; i<nRows; i++){
    X[ix(i,0,nColumns+1)] = 1;
// cout << "X[" << ix(i,0,nColumns+1) << "]=1" << endl;
    // store predictor values
    for(j=1; j<=nColumns; j++){
      X[ix(i,j,nColumns+1)] = data[includedCells[i]][j-1];
    }

    // status handled differently if data is summary
    if(summary_data){
      Y0[i] = (unsigned int)(data[includedCells[i]][unaffIndex]);
      sY0 += Y0[i];
      Y1[i] = (unsigned int)(data[includedCells[i]][affIndex]);
      sY1 += Y1[i];
    }
    else{ // not summary
      if(data[includedCells[i]].back() == 0){
        Y0[i] = 1; sY0++;
      }
      else{
        Y1[i] = 1; sY1++;
      }
    }

    sC += Y0[i] + Y1[i];

    for(j=1; j<= nColumns; j++){
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
// cout << "X from " << X[ix(i,j,nColumns+1)] << " ";
      X[ix(i,j,nColumns+1)] = ( X[ix(i,j,nColumns+1)] - xM[j] ) / xSD[j];
// cout << X[ix(i,j,nColumns+1)] << endl;
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
// cout << "v=" << v << endl;
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
// cout << "df=" << nColumns ;
    overallPvalue = ChiSq(fabs(LLn-LL), nColumns);
// cout << " pvalue = " << overallPvalue << endl;
// exit(1);

    coeff_intercept = Par[0];
    // adjust coefficients so that the zero index is now first coefficient
//cout << "-----------" << endl;
    for(j=1; j<=nColumns; j++) {
      coefficients[j-1] = Par[j];
standard_errors[j-1] = SEP[j];
//cout << "coef=" << coefficients[j-1] << " stderror=" << standard_errors[j-1] << endl;
    }
//cout << "-------------" << endl;

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
//cout << "in set_model" << endl;
//cout << "modType=" << modType << endl;
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
void LogisticRegression::summarize_data(vector<unsigned int> & genos){
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
  for(unsigned int curr_mark=0; curr_mark < genos.size(); curr_mark++){
    ref_alleles[curr_mark] = set->get_locus(genos[curr_mark])->getReferentIndex();
  }

  // add to summary totals for each genotype
  // use current model to convert genotypes
  for(currInd=0; currInd < numInds; currInd++){
    for(currLoc=0; currLoc < combSize;
      currLoc++){
      // genotype[currLoc] = geno_convert[dataset[currInd][genos[currLoc]]];
      // genotype[currLoc] = geno_convert[(*set)[currInd]->get_genotype(genos[currLoc])];
      genotype[currLoc] = geno_convert[ref_alleles[currLoc]][(*set)[currInd]->get_genotype(genos[currLoc])];
    }

    // increment count based on status of individual
    // summary_data[indexConverter.flatten_indexes(genotype)][unaffIndex+dataset[currInd].status]++;
    summary_data[indexConverter.flatten_indexes(genotype)][unaffIndex+(*set)[currInd]->getAffected()]++;
  }

//output summary data
// cout << "Summary data in summarize_data:" << endl;
// for(uint y=0; y<summary_data.size(); y++){
//   for(uint z=0; z<summary_data[y].size(); z++){
//      cout << summary_data[y][z] << " ";
//   }
//   cout << endl;
// }
// cout << "------- end summary data ---------" << endl;


}

///
/// Sizes and clears arrays for holding data
/// @param currModelSize size of model
/// @return
///
void LogisticRegression::initialize_summary(unsigned int currModelSize){

//   unsigned int counter;

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
void LogisticRegression::calculate(vector<unsigned int>& loci,
  vector<unsigned int>& covars, vector<unsigned int> & traits){

//for(uint i=0; i<geno_convert.size(); i++){
//  cout << "genoconvert i=" << geno_convert[i] << endl;
//}

// cout << "negative calculate" << endl;

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
  for(unsigned int curr_mark=0; curr_mark < converted_loci.size(); curr_mark++){
    ref_alleles[curr_mark] = set->get_locus(converted_loci[curr_mark])->getReferentIndex();
  }

  unsigned int currInd, currValue, i;
  for(currInd=0; currInd < numInds; currInd++){
    currValue = 0;
    bool any_missing = false;
    for(i=0; i < numLoci; i++){
//  cout << "geno=" << (*set)[currInd]->get_genotype(converted_loci[i]) <<  " missingValue = " << missingValue << endl;

//       if((*set)[currInd]->get_genotype(converted_loci[i]) != missingValue)
//         row[currValue++] = geno_convert[(*set)[currInd]->get_genotype(converted_loci[i])];
      if((*set)[currInd]->get_genotype(converted_loci[i]) != missingValue){
        row[currValue++] = geno_convert[ref_alleles[i]][(*set)[currInd]->get_genotype(converted_loci[i])];
      }
      else{
        any_missing = true;
        break;
      }
    }
//  cout << "any_missing=" << any_missing << endl;
//  exit(1);
    for(i=0; i < numCovars; i++){
      if((*set)[currInd]->getCovariate(covars[i]) != missingCoValue)
        row[currValue++] = (*set)[currInd]->getCovariate(covars[i]);
      else{
        any_missing = true;
        break;
      }
    }

    for(i=0; i < numTraits; i++){
      if((*set)[currInd]->getTrait(traits[i]) != missingCoValue)
        row[currValue++] = (*set)[currInd]->getTrait(traits[i]);
      else{
        any_missing = true;
        break;
      }
    }
// cout << "any_missing=" << any_missing << endl;
    if(!any_missing){
      row[currValue] = ((*set)[currInd]->getAffected());
      summary_data.push_back(row);
      includedCells.push_back(summary_data.size()-1);
    }

  }

// output data in summary_data
/*
 for(uint i=0; i<summary_data.size(); i++){
   if(summary_data[i].size() > 0){
     cout << summary_data[i][0];
   }
   for(uint j=1; j<summary_data[i].size(); j++){
     cout << "," << summary_data[i][j];
   }
   cout << endl;
 }
 */
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

void LogicRegression::set_anneal_parameters(int start, int end, int iter, int earlyout, int update){
	anneal_start = start;
	anneal_end = end;
	anneal_iter = iter;
	anneal_earlyout = earlyout;
	anneal_update = update;
	if(iter < 100 && (start != 0 || end != 0)){
		throw MethodException("Logic Regression: Anneal Parameters: Not enough repetitions.");
	}
	if(start < end){
		throw MethodException("Logic Regression: Anneal Parameters: Starting temperature below ending temperature.");
	}

}

void LogicRegression::logreg(vector<int> bin, vector<double> wgt, vector<double> cens, int type, int select,
		vector<int> ntrees, vector<int> nleaves)

void LogicRegression::logreg(vector<int> resp, vector<int> bin, vector<int> sep, vector<double> wgt, vector<double> cens,
		int type, int select, vector<int> ntrees, vector<int> nleaves, vector<int> penalty, vector<int> seed,
		vector<int> kfold, vector<int> nrep, Fitness* oldfit){//, tree.control, mc.control){

	TreeAnnealParameters tree_control;
	//build matrices?  Phenotype + genotype columns, rows are samples

	//tree.control = unknown
	//mc.control = unknown
	//resp = phenotype column
	//bin = binary columns (genotypes, covars, etc)
	//anneal.control = anneal_* vars
	//seed = seed value?
	//penalty = penalty value?
	//type=mdl can be 1,2,3,4,5,0 or C,R,L,S,E,O (<- letter O)
	//select=choice cat be 1,2,3,4,5,6,7 or S,M,C,N,R,G,B
	//sep = separate predictors?

	if(oldfit != NULL){
		if(oldfit->getResp().size() > 0){
			resp = oldfit->getResp();
		}
		if(oldfit->getBinNames().size() > 0){
			binnames = oldfit->getBinNames();
		}
		if(oldfit->getBinary().size() > 0){
			bin = oldfit->getBinary();
		}
		if(oldfit->getNsep() > 0){
			sep = oldfit->getSeparate();
		}
		if(oldfit->getCensor().size() > 0){
			cens = oldfit->getCensor();
		}
		if(oldfit->getWeight().size() > 0){
			wgt = oldfit->getWeight();
		}
		if(oldfit->getNrep().size() > 0){
			nrep = oldfit->getNrep();
		}
		if(oldfit->getKfold().size() > 0){
			kfold = oldfit->getKfold();
		}
		if(oldfit->getPenalty().size() > 0){
			penalty = oldfit->getPenalty();
		}
		if(oldfit->getAnnealControl().size() > 0){
			anneal_control = oldfit->getAnnealControl();
		}
		//mc.control
		//tree.control
		if(oldfit->getNTrees().size() > 0){
			ntrees = oldfit->getNTrees();
			if(ntrees.size() == 1){
				//ntrees <- c(ntrees, ntrees)
				ntrees.clear();
				ntrees.push_back(ntrees[0]);
				ntrees.push_back(ntrees[0]);
			}
		}
		if(oldfit->getNLeaves().size() > 0){
			nleaves = oldfit->getNLeaves();
			if(nleaves.size() == 1){
				//nleavs <- c(nleaves, nleaves)
				nleaves.clear();
				nleaves.push_back(nleaves[0]);
				nleaves.push_back(nleaves[0]);
			}
		}
		if(oldfit->getSeed().size() > 0){
			seed = oldfit->getSeed();
		}
	}
		int mdl = type;
		if(mdl < 0 || mdl > 5){
			throw MethodException("model not implemented");
		}

		int n1 = resp.size();

		int choice = select;
		if(choice < 0 || choice > 7){
			throw MethodException("not a valid selection approach");
		}

		if((mdl == 1 || mdl == 3) && (logreg.binary(resp) == false)){
			throw MethodException("some non binary response data, and a binary fitting function");
		}
		if((mdl == 4 || mdl == 5) && *(min_element(resp.begin(), resp.end())) < 0){
			throw MethodException("survival data needs to be positive");
		}

		//check data integrity?? responses x bindata availalbe?

		int nsep = sep.size();
		if(sep.size() > 0){
			//check integrity?
			if(sep.size() != n1){
				throw MethodException("length of separate predictors not a multiple of response");
			}
		}

		if(nsep > 0 && (mdl == 2 || mdl == 5 || mdl == 4)){
			vector<int> tep = sep;
			int k = 0;
			for(unsigned int i = 0; i < nsep; i++){
				int min = *(min_element(tep.begin(), tep.begin() + i));
				int max = *(max_element(tep.begin(), tep.begin() + i));

				if((min != 0 || max != 1) && k != 2){
					k = 1;
					if(max != min){
						tep[,i] = (tep[,i] - min) / (max - min);
						if(!logreg.binary(tep[,i])){
							k = 2;
						}
					}
					else{
						tep[,i] = 0;
					}
				}
			}

			if(k == 1 && mdl 1= 4){
				sep = tep;
				opts::printLog("separate covariates that were binary recoded to 0/1");
			}
			if(k == 2 && mdl != 4){
				if(mdl == 5){
					throw MethodException("Exponential surviaval models are only implemented if all separate preditors are binary");
				}
				opts::printLog("Logistic regression runs much faster if all separate predictors are binary.");
			}
			if(mdl == 4 && k != 2){
				opts::printLog("The code will be run much faster if you use exponential regression instead"
						+ " of proportional hazards models after a cumulative hazards transformation, "
						+ " as all your separate predictors are binary.");
			}
		}

		if(nesp == 0 && mdl == 4){
			opts::printLog("The code will be run much faster if you use exponential regression instead"
					+ " of proportional hazards models after a cumulative hazards transformation, "
					+ " as all your separate predictors are binary.");
		}

		if(tree == NULL){
			tree_control = logreg.tree.control(n1);
		}
		else{
			tree_control = logreg.tree.control(tree_control, n1);
		}

		if(wgt.size() == 0){
			wgt.resize(n1, 1);
		}
		else{
			if(wgt.size() != n1){
				throw MethodException("weight does not have the right length");
			}
			if(*(min_element(wgt.begin(), wgt.end())) < 0 && mdl != 0){
				throw MethodException("weight has to be larger than 0");
			}
		}

		if(cens.size() == 0){
			cens.resize(n1, 1);
		}
		else{
			if(cens.size() != n1){
				throw MethodException("cens doesn't have right length");
			}
			if((mdl == 4 || mdl == 5) && (logreg.binary(cens) == false)){
				throw MethodException("censoring data has to be 0/1");
			}
		}

		vector<int> ntr;
		vector<int> msz;
		if(choice == 1 || choice == 4 || choice == 7){
			if(ntrees.size() == 0){
				ntr.resize(2,1);
			}
			else{
				ntr.resize(2);
				ntr.push_back(ntrees[0]);
				ntr.push_back(ntrees[1]);
			}
			if(nleaves.size() == 0){
				msz.resize(2,-1);
			}
			else{
				msz.resize(2);
				msz[0] = nleaves[0];
				msz[1] = nleaves[1];

			}
		}
		else{
			if(ntrees.size() == 0){
				ntr.resize(2,1);
			}
			else{
				ntr.resize(2);
				ntr[0] = ntrees[0];
				ntr[1] = ntrees[1];
			}
			if(ntr[1] == -1){
				ntr[1] = ntr[0];
				if(choice == 6 && nleaves.size() == 0){
					nleaves.resize(tree.treesize * ntr[1], 2);
				}
			}
			msz.resize(2);
			msz[0] = nleaves[0];
			msz[1] = nleaves[1];
			if(msz[1] == -1){
				msz[1] = msz[0];
			}
			if(ntr[0] > ntr[1]){
				throw MethodException("upper limit of number of trees smaller than lower limit");
			}
			if(msz[0] > msz[1]){
				throw MethodException("upper limit of number of leaves smaller than lower limit");
			}
			if(msz[1] < 0){
				throw MethodException("number of leaves needs to be at least 0");
			}


		}

		if(ntr[0] < 1){
			throw MethodException("number of trees needs to be at least 1");
		}
		if(mdl == 1 && *(max_element(ntr.begin(), ntr.end())) > 1){
			throw MethodException("for classification only 1 tree is possible");
		}
		if(mdl == 1 && nsep > 0){
			throw MethodException("for classification no separate predictors are possible");
		}
		if(anneal_control.size() == 0){
			anneal_contorl = logreg.anneal.control();
		}
		else{
			anneal_control = logreg.anneal_control(anneal_control);
		}

		if(mc_control.size() == 0){
			mc_control = logreg.mc.control();
		}
		else{
			mc_control = logreg.mc.control(mc_control);
		}

		if(mc_control.size() > 0 && select == 7){
			anneal_control.update = mc_control.update;
		}

		if(penalty < 0 && (choice == 1 || choice == 6)){
			throw MethodException("penalty should be at leaset 0");
		}

		if(choice == 3){
			nrep = kfold;
			if(nrep < 2){
				throw MethodException("kfold needs to be at least 2");
			}
		}

		if(choice > 3 && nrep < 1 && choice < 6){
			throw MethodException("nrep needs to be at least 1");
		}

		int xseed = seed;
		if(xseed == 0){
			//get random seed???
		}

		//ipars??
		//rpars??

		int nkn = tree_control.treesize * 2 - 1;
		int nxx = 2;

		int na = 0;
		int nb = 0;
		int nc = 0;
		int nd = 0;

		if(choice == 1 || choice == 7){
			na = ntr[0] * (nkn * 4 + 3);
			nb = nsep + ntr[0] + 1;
			nc = 2;
		}

		if(choice == 2){
			nd = (ntr[1] - ntr[0] + 1) & (msz[1] - msz[0] + 1);
			na = nd * ntr[1] * (nkn * 4 + 3);
			nb = nd * (nsep + ntr[1] + 1);
			nc = nd;
		}

		if(choice == 6){
			nd = msz[1] + 2;
			na = nd * ntr[1] * (nkn * 4 + 3);
			nb = nd * (nsep + ntr[1] + 1);
			nc = nd;
		}

		if(choice == 3){
			na = nb = 2;
			nd = (ntr[1] - ntr[0] + 1) * (msz[1] - msz[0] + 1);
			nc = nd * 8 * nrep;
		}

		if(choice == 4){
			na = nb = 2;
			nc = nrep + 2;
		}

		if(choice == 7){
			na = 256;
			nb = n2;
			nc = n2 * n2;
			if(abs(mc_control.output) < 2){
				nc = 1;
			}
			nxx = nc * n2;
			if(abs(mc_control.output) < 3){
				nxx = 1;
			}
		}
		vector<int> xtree;
		if(choice != 5){
			xtree.resize(na, -100);
		}

		if(choice == 5){
			nb = 2;
			nd = (ntr[1] - ntr[0] + 1) * (msz[1] - msz[1] + 1);
			nc = (nrep + 2) * nd + 2;
			xtree.clear(); //= NULL

			if(tree_control.treesize != oldfit.getTreeControl().treesize){
				throw MethodException("treesize should match treesize in oldfit");
			}
			if(oldfit->getChoice() != 2){
				throw MethodException("oldfit not an object from a multiple model fit");
			}
			if(oldfit->getChoice() == 2){
				for(int i = 0; i < oldfit->getNModels(); i++){
					tmptree = oldfit->getTree(i);
					if(i == 1){

					}
				}
			}
		}

		Fitness fit = slogreg_(n1, n2, nsep, ipars, rpars, sep, cens, orders, resp, wgt, bin, xtree, rep(-100, nb),
				rep(-100, nc), rep(0,nxx));

		ipars = choice;
		rpars = penalty;

		if(mdl == 0){
			type = "own.scoring";
		}
		else if(mdl == 1){
			type = "classification";
		}
		else if(mdl == 2){
			type = "regression";
		}
		else if(mdl == 3){
			type = "logistic";
		}
		else if(mdl == 4){
			type = "proportional.hazards";
		}
		else if(mdl == 5){
			type = "exponential.survival";
		}

		string chs = "";
		if(choice == 1){
			chs = "single.model";
		}
		else if(choice == 2){
			chs = "multiple.models";
		}
		else if(choice == 3){
			chs = "cross.validation";
		}
		else if(choice == 4){
			chs = "null.model.test";
		}
		else if(choice == 5){
			chs = "randomization.test";
		}
		else if(choice == 7){
			chs = "bayesian";
		}
		else if(choice == 6){
			chs = "greedy";
		}

		if(tree_control.getOpers() != 2 && tree_control.getOpers() != 3){
			tree_control.setOperators("both");
		}
		else if(tree_control.getOpers() == 2){
			tree_control.setOperators("and");
		}
		else if(tree_control.getOpers() == 3){
			tree_control.setOperators("or");
		}


		if(choice == 7){

		}

}




void LogicRegression::eval_logreg(string type){
	if(type == "logregtree"){
	}
	else if(type == "logregmodel"){

	}
}
}
