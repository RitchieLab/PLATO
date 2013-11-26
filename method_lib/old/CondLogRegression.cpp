//CondLogRegression.cpp

#include "CondLogRegression.h"
#include "ComboGenerator.h"

namespace Methods{  
///
/// constructor -- initialize variables
///
CondLogRegression::CondLogRegression(){
  set=NULL;
  initialize();
}


///
/// Alternative constructor
///
CondLogRegression::CondLogRegression(DataSet* ds){
  resetDataSet(ds);
  initialize();
}

///
/// Sets DataSet pointer for running analysis
/// @param ds DataSet 
///
void CondLogRegression::resetDataSet(DataSet* ds){
  set = ds;
  missingValue = set->get_missing_value();
//  param.nvmax = 0;
  select_individuals();
}


///
/// initializes maps used in switches for setting parameters
/// @return none
///
void CondLogRegression::initialize(){
  modelSize = 0;
  param.nvmax=0;
  modType = Additive;
//   pValType = Overall;
  maxLocusValue = 2;
  missingValue = maxLocusValue+1;
  LociComboLimit = 10;
  LociComboMin = 1;
  maxIterations = 20;
  defaultComboInterval = 10000;
  original_strat_size = 0;

  interaction_included = true;
  
  ModelTypeMap["DOMINANT"] = Dominant;
  ModelTypeMap["RECESSIVE"] = Recessive;
  ModelTypeMap["ADDITIVE"] = Additive;
  
  PiD2 = 3.141592653589793 / 2;
  
  vector<Sample*> empty;
  strata_inds.push_back(empty);
  strata_inds.push_back(empty);
  
  set_model();
  initialize_interactions();
}


///
/// Sets values to use for each genotype
/// @return
///
void CondLogRegression::set_model(){
  
  geno_convert.assign(2, vector<unsigned int>(4,0));
  
  switch(modType){
    case Dominant:
      geno_convert.at(0).at(0) = 1;
      geno_convert.at(0).at(1) = 1;
      geno_convert.at(0).at(3) = 2;
      geno_convert.at(1).at(1) = 1;
      geno_convert.at(1).at(2) = 1;
      geno_convert.at(1).at(3) = 2;
      maxLocusValue = 1;
      break;
    case Recessive:
      geno_convert.at(0).at(0) = 1;
      geno_convert.at(0).at(3) = 2;
      geno_convert.at(1).at(2) = 1;
      geno_convert.at(1).at(3) = 2;
      maxLocusValue = 1;
      break;
    case Additive:
      geno_convert.at(0).at(1) = 1;
      geno_convert.at(0).at(0) = 2;
      geno_convert.at(0).at(3) = 3;
      geno_convert.at(1).at(1) = 1;
      geno_convert.at(1).at(2) = 2;
      geno_convert.at(1).at(3) = 3;
      maxLocusValue = 2;
      break;
  }

}


///
/// Generates and stores interaction lists for producing interaction
/// terms when doing multi-locus models
/// @return
///
void CondLogRegression::initialize_interactions(){
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
        interaction_lists[curr_num].push_back(generator.ComboList.at(comboIndex));
      }
    }while(!done);

  }
}

///
/// Sets the parameters using StepOptions class
/// @param options StepOptions containing options
///
void CondLogRegression::set_parameters(StepOptions* options){
    setFullInteraction(options->getCondLRFullInteraction());
    setIncludeInteractions(options->getCondLRIncludeInteractions());
    setMaximumIterations(options->getCondLRMaximumIterations());
    setModelType(options->getCondLRModelType());
}

///
/// Sets the model type to use in the calculation
/// @param modelType string containing model type to run
/// @throws MethodException on error
///
void CondLogRegression::setModelType(string modelType){
  if(ModelTypeMap.find(modelType) != ModelTypeMap.end()){
    modType = ModelTypeMap[modelType];
    set_model();
  }
  else{
    throw MethodException(modelType + " is not a valid Conditional LR model type");
  }
}


///
/// Selects an affected and unaffected from each pedigree <br>
/// Each pair will make up the strata for the conditional LR
/// @return
///
void CondLogRegression::select_individuals(){
  // at each index the two inds make a stratum for use
  // in the conditional LR
  strata_inds.at(0).clear();
  strata_inds.at(1).clear();
  
  include_strata.clear();
  
  unsigned int num_pedigrees = set->num_pedigrees();
  
  bool affected_found = false;
  bool unaffected_found = false;
  unsigned int num_inds, curr_ind;
  Sample* ind;
  vector<Sample*>* samples;
  
  for(unsigned int currped=0; currped < num_pedigrees; currped++){
    affected_found = false;
    unaffected_found = false;
    num_inds = set->get_pedigree(currped)->getTotalInds();
    samples = set->get_pedigree(currped)->getSamples();

    curr_ind = 0;
    while((affected_found != true || unaffected_found != true) && curr_ind < num_inds){
//       ind = (*(set->get_pedigree(currped).get_samples()))[curr_ind];
      ind = (*samples)[curr_ind];
      // to be used must match the status
//        if(ind->getDad() != NULL){
        if(ind->getAffected() && !affected_found){
          strata_inds.at(1).push_back(ind);
          affected_found=true;
        }
        else if(!ind->getAffected() && !unaffected_found){
          strata_inds.at(0).push_back(ind);
          unaffected_found=true;
        }
//        }
      curr_ind++;
    } 
    include_strata.push_back(currped);
  }
  original_strat_size = include_strata.size();
}



///
/// Runs conditional logistic regression on dataset passed 
/// @param loci
/// @param dataset DataSet to be analyzed
/// @return none
///
void CondLogRegression::calculate(vector<unsigned int>& loci, DataSet* dataset){
  resetDataSet(dataset);
  calculate(loci);
}


///
/// Runs conditional logistic regression on loci
/// @param loci
/// @return none
///
void CondLogRegression::calculate(vector<unsigned int>& oloci){

  vector<unsigned int> loci;
  for(vector<unsigned int>::iterator iter=oloci.begin(); iter!=oloci.end(); iter++)
    loci.push_back(set->get_markers()->at(*iter)->getLoc());

  if((unsigned int)(param.nvmax) != (loci.size() + interaction_included)){
    // free any previously allocated memory
    if(param.nvmax > 0)
      free_params(param);
        
    // allocate memory for the arrays
    allocate_arrays(param, loci);
  }

  // when missing data exists in set
  // calculate number of pedigrees to include
  if(set->missing_data_present()){
    adjust_arrays(param, loci);
  }
  
  fill_z_array(param, loci);
    
  likelihood_ratio = run_conditional_lr(param);
  // alter so only does likelihood ratio
//  p_value = ChiSquare::pfromchi(likelihood_ratio, 2);
  p_value = ChiSq(likelihood_ratio, param.nv);

    ////////////////////////////////////
    // covariate matrix is as follows  
    //      snp1   snp2    Int
    // B1   cov1
    // B2   cov2   cov3
    // B3   cov4   cov5    cov6
    /////////////////////////////////////
  
  // push back coefficients and covariates
  coefficients.clear();
  covariates.clear();
  unsigned int end_z = (unsigned int)(param.nv * (float(param.nv+1))/2);
  for(int z=0; z<param.nv; z++){
    coefficients.push_back(param.b[z]);
  }
  for(unsigned int z=0; z< end_z; z++){
    covariates.push_back(param.cov[z]);
  }
  
//   if(num_results > 0)  
//     free_params(param);
}


///
/// Runs conditional LR and returns  p value
/// @param
/// @return p value
///
float CondLogRegression::run_conditional_lr(conditional_lr_parameters& param){

  // need to construct the data structure for the 
  
  // parameters into the 
  // int ns --> number of strata
  // int[] nca  -> size(ns) number  of cases in each stratum
  // int[] nct --> size (ns) number of controls in each stratum
  // int nimax --> maximum number of cases plus controls permitted in any analysis
  // int nmax --> maximum  number of cases plus controls permitted in each stratum
  // int nmax1 --> nmax+1
  // int nvmax = maximum number of variables permitted
  // int nvmax1 --> nvmax*(nvmax+1)/2
  // int nv --> number of parameters to be estimated
  // int [] z --> size (nvmax, nimax) matrix containing covariate values where z[c,ind] is the value of ind at that covariate
  //        individuals follow as first cases, then first control, then second case & second control 
  //        for each stratum
  // int [] ivar --> size(nvmax) first nv positions specify the indices of the nv parameters to be estimated
  // float [] b --> size(nvmax) first nv positions of b give the initial estimates of the parameters
  //               as specified by ivar.
//output_params(param);

  // run conditional LR
	logcch_(&param.ns, param.nca, param.nct, &param.nimax, &param.nmax, &param.nmax1,
	  &param.nvmax, &param.nvmax1, &param.nv, param.z, param.ivar, param.covi, param.cntr,
	  param.w, param.wb, param.wdb, param.wd2b, param.u, param.ins, param.db,
	  param.d2b, param.dl, param.b, param.cov, &param.chi2, &param.st, &param.ifault);

float output = param.chi2;
// int logcch_(integer *ns, integer *nca, integer *nct, integer *nimax, integer *nmax, 
// integer *nmax1, integer *nvmax, integer *nvmax1, integer *nv, float **z__, 
// integer *ivar, float *covi, float **cntr, float *w, float *wb, 
// float **wdb, float **wd2b, float *u, integer *ins, float *db, 
// float *d2b, float *dl, float *b, float *cov, float *chi2, float *st, integer *ifault);  

//   free_params(param);
//exit(1);
  
  return output;  

}


///
/// Allocates memory for parameters that are not dependent on 
/// the number of individuals in the set
/// @param params condtional_lr_parameters
/// 
void CondLogRegression::allocate_arrays(conditional_lr_parameters& params, vector<unsigned int>& loci){
  unsigned int num_loci = loci.size();
    // input parameters that aren't dependent on number of pedigrees
  params.nmax = 2;
  params.nmax1 = params.nmax+1;
  if(!interaction_included)
    params.nvmax = num_loci;    
  else
    params.nvmax = num_loci + interaction_lists[num_loci].size();
  params.nvmax1 = params.nvmax * (params.nvmax+1) / 2; 
  params.nv = params.nvmax;
  params.ivar = new integer[params.nvmax];
  params.b = new doublereal[params.nvmax];
  
  // workspace parameters that aren't dependent on number of pedigrees
  params.covi = new doublereal[params.nvmax1];
  params.w = new doublereal[params.nvmax];
  params.wb = new doublereal[params.nmax1];
  params.wdb = new doublereal[params.nvmax * params.nmax1];
  params.wd2b = new doublereal[params.nvmax1 * params.nmax1];
  params.u = new doublereal[params.nmax];
  params.db = new doublereal[params.nvmax];
  params.d2b = new doublereal[params.nvmax1];
  
  // output params that aren't dependent on number of pedigrees
  params.dl = new doublereal[params.nvmax];
  params.cov = new doublereal[params.nvmax1];
  
  // set sizes initially to be same as size with no missing data
  params.ns = include_strata.size();
  params.nca = new integer[params.ns*params.nmax];
  params.nct = new integer[params.ns*params.nmax];
   
  params.nimax = params.ns*params.nmax;
  
  params.z = new doublereal[params.nvmax * params.nimax];
  
  params.cntr = new doublereal[params.nvmax * params.ns];
  params.ins = new integer[params.ns];  
  
  fill_arrays(params);
}


///
/// Determines which pedigrees will be included for this set
/// of loci and resizes arrays that are dependent on it to match
/// @param params conditional_lr_parameters
/// @param loci SNPs to analyze
/// @return
///
void CondLogRegression::adjust_arrays(conditional_lr_parameters& params, 
  vector<unsigned int>& loci){
  unsigned int num_loci = loci.size();
  unsigned int currloc;
  unsigned int previous_strat_size = include_strata.size();
  include_strata.clear();
  bool include_this_stratum;

  for(unsigned int currstrata=0; currstrata< original_strat_size; currstrata++){
    include_this_stratum = true;
    for(currloc=0; currloc < num_loci; currloc++){
      if(strata_inds.at(0).at(currstrata)->get_genotype(loci.at(currloc)) == (unsigned int)(missingValue) ||
        strata_inds.at(1).at(currstrata)->get_genotype(loci.at(currloc)) == (unsigned int)(missingValue)){
        include_this_stratum = false;
        break;
      }
    }
    if(include_this_stratum){
      include_strata.push_back(currstrata);
    }
  }

  
  if(previous_strat_size != include_strata.size()){
    delete [] params.nca;
    delete [] params.nct;
    delete [] params.z;
    delete [] params.cntr;
    delete [] params.ins;
  
    params.ns = include_strata.size();
    params.nca = new integer[params.ns*params.nmax];
    params.nct = new integer[params.ns*params.nmax];
    
    params.nimax = params.ns*params.nmax;
  
    params.z = new doublereal[params.nvmax * params.nimax];
 
    params.cntr = new doublereal[params.nvmax * params.ns];
    params.ins = new integer[params.ns];      
    
    fill_arrays(params);
  }
}


///
/// Fills arrays that are dependent on number of strata
/// @param params conditional_lr_parameter
///
void CondLogRegression::fill_arrays(conditional_lr_parameters& params){
  // set number of cases and controls in each stratum (pedigree)
  for(int stratum=0; stratum < params.ns; stratum++){
    params.nca[stratum] = 1;
    params.nct[stratum] = 1;
  }  

    
}


///
/// Only includes the pedigree (stratum) when no missing data present
/// @param dataset DataSet with analysis data
/// @return
///
void CondLogRegression::fill_params(conditional_lr_parameters& params, 
  vector<unsigned int>& loci){
    
  unsigned int num_loci = loci.size();
  
  // input parameters that aren't dependent on number of pedigrees
  params.nmax = 2;
  params.nmax1 = params.nmax+1;
  if(!interaction_included)
    params.nvmax = num_loci;    
  else
    params.nvmax = num_loci + interaction_lists[num_loci].size();
  params.nvmax1 = params.nvmax * (params.nvmax+1) / 2; 
  params.nv = params.nvmax;
  params.ivar = new integer[params.nvmax];
  params.b = new doublereal[params.nvmax];
  
  // workspace parameters that aren't dependent on number of pedigrees
  params.covi = new doublereal[params.nvmax1];
  params.w = new doublereal[params.nvmax];
  params.wb = new doublereal[params.nmax1];
  params.wdb = new doublereal[params.nvmax * params.nmax1];
  params.wd2b = new doublereal[params.nvmax1 * params.nmax1];
  params.u = new doublereal[params.nmax];
  params.db = new doublereal[params.nvmax];
  params.d2b = new doublereal[params.nvmax1];
  
  // output params that aren't dependent on number of pedigrees
  params.dl = new doublereal[params.nvmax];
  params.cov = new doublereal[params.nvmax1];
  
  // when missing data, check which individuals will be included
  if(set->missing_data_present()){
    unsigned int currloc;
    include_strata.clear();
    bool include_this_stratum;
    for(unsigned int currstrata=0; currstrata< strata_inds.size(); currstrata++){
      include_this_stratum = true;
      for(currloc=0; currloc < num_loci; currloc++){
        if(strata_inds.at(0).at(currstrata)->get_genotype(loci.at(currloc)) == (unsigned int)(missingValue) ||
          strata_inds.at(1).at(currstrata)->get_genotype(loci.at(currloc)) == (unsigned int)(missingValue)){
          include_this_stratum = false;
          break;
        }
      }
      if(include_this_stratum){
        include_strata.push_back(currstrata);
      }
    }
  }

  params.ns = include_strata.size();
  params.nca = new integer[params.ns*params.nmax];
  params.nct = new integer[params.ns*params.nmax];
  // set number of cases and controls in each stratum (pedigree)
  for(int stratum=0; stratum < params.ns; stratum++){
    params.nca[stratum] = 1;
    params.nct[stratum] = 1;
  }
    
  params.nimax = params.ns*params.nmax;
  
  params.z = new doublereal[params.nvmax * params.nimax];
  
  params.cntr = new doublereal[params.nvmax * params.ns];
  params.ins = new integer[params.ns];  
  
  // fill ivar -- holds the positions of the variables to be included
  for(int var=0; var<params.nv; var++){
    params.ivar[var] = var+1;
    params.b[var] = 0.0;
  }
  
  fill_z_array(params, loci);
  
}

/*    struct conditional_lr_parameters{
      float * z;
      integer * nca;
      integer * nct;
      integer * ivar;
      integer * ins;
      float * b;
      integer ns, nimax, nmax, nvmax, nvmax1, nv, nmax1;
      float *cntr, *w, *wb, *wdb, *wd2b, *u,  *db, *d2b, *dl, *cov, *covi;
      float chi2, st;
      integer ifault;
    }; */
    
void CondLogRegression::output_params(conditional_lr_parameters& params){
  cout << "ns=" << params.ns << endl;
  cout << "nimax=" << params.nimax << endl;
  
  for(int i=0; i<params.ns; i++){
    cout << "nca[" << i << "]=" << params.nca[i] << endl;
  }
  for(int i=0; i<params.ns; i++){
    cout << "nct[" << i << "]=" << params.nct[i] << endl;
  }
  
  unsigned int z_total = params.nvmax * params.nimax;
  for(unsigned int i=0; i<z_total; i++){
    cout << "z[" << i << "]=" << params.z[i] << endl;
  }
  
  for(int i=0; i<params.nvmax; i++){
    cout << "b[" << i << "]=" << params.b[i] << endl;
  }
  
  for(int i=0; i<params.nvmax; i++){
    cout << "ivar[" << i << "]=" << params.ivar[i] << endl;
  }  
  
  cout << "nmax=" << params.nmax << endl;
  cout << "nvmax=" << params.nvmax << endl;
  cout << "nvmax1=" << params.nvmax1 << endl;
  cout << "nv=" << params.nv << endl;
  cout << "nmax1=" << params.nmax1 << endl;
  cout << "ifault=" << params.ifault << endl;
  cout << "st=" << params.st << endl;
  cout << "chi2=" << params.chi2 << endl;
}


///
/// Fills z array using strata to include
/// @param strata Lists indices of strata to include in z array for analysis
/// Strata with any missing data are excluded
/// @return
///
void CondLogRegression::fill_z_array(conditional_lr_parameters& params, 
  vector<unsigned int>& loci){
      
  for(int var=0; var<params.nv; var++){
    params.ivar[var] = var+1;
    params.b[var] = 0.0;
  } 
  
  unsigned int num_loci = loci.size();
  unsigned int curr_z_index = 0;

  unsigned int currloc, interact;
 
  vector<int> ref_alleles(loci.size(), 0);
  // determine which geno_conversion to use for each locus
  for(unsigned int curr_mark=0; curr_mark < loci.size(); curr_mark++){
    ref_alleles[curr_mark] = set->get_locus(loci[curr_mark])->getReferentIndex();
  }
  
  for(unsigned int stratum=0; stratum < include_strata.size(); ++stratum){
    for(currloc=0; currloc < num_loci; currloc++){
//     params.z[curr_z_index++] = geno_convert[strata_inds[1][include_strata[stratum]]->get_genotype(loci[currloc])];
    params.z[curr_z_index++] = geno_convert.at(ref_alleles.at(currloc)).at(strata_inds.at(1).at(include_strata.at(stratum))->get_genotype(loci.at(currloc)));
    }
    if(interaction_included){
      for(interact=0; interact < interaction_lists.at(num_loci).size(); ++interact){
//         params.z[curr_z_index] = geno_convert[strata_inds[1][include_strata[stratum]]->get_genotype(loci[interaction_lists[num_loci][interact][0]])];
        params.z[curr_z_index] = geno_convert.at(ref_alleles.at(interaction_lists.at(num_loci).at(interact).at(0))).at(strata_inds.at(1).at(include_strata.at(stratum))->get_genotype(loci.at(interaction_lists.at(num_loci).at(interact).at(0))));
        for(currloc=1; currloc < interaction_lists.at(num_loci).at(interact).size(); currloc++){
          params.z[curr_z_index] *= 
//             geno_convert[strata_inds[1][include_strata[stratum]]->get_genotype(loci[interaction_lists[num_loci][interact][currloc]])];
            geno_convert.at(ref_alleles.at(interaction_lists.at(num_loci).at(interact).at(0))).at(strata_inds.at(1).at(include_strata.at(stratum))->get_genotype(loci.at(interaction_lists.at(num_loci).at(interact).at(currloc))));
        }
        curr_z_index++;
      }
    }
    for(currloc=0; currloc < num_loci; currloc++){     
//       params.z[curr_z_index++] = geno_convert[strata_inds[0][include_strata[stratum]]->get_genotype(loci[currloc])];
      params.z[curr_z_index++] = geno_convert.at(ref_alleles.at(currloc)).at(strata_inds.at(0).at(include_strata.at(stratum))->get_genotype(loci.at(currloc)));
    }
    if(interaction_included){
      for(interact=0; interact < interaction_lists.at(num_loci).size(); ++interact){
         params.z[curr_z_index] = 
//            geno_convert[strata_inds[0][include_strata[stratum]]->get_genotype(loci[interaction_lists[num_loci][interact][0]])];
           geno_convert.at(ref_alleles.at(interaction_lists.at(num_loci).at(interact).at(0))).at(strata_inds.at(0).at(include_strata.at(stratum))->get_genotype(loci.at(interaction_lists.at(num_loci).at(interact).at(0))));
        for(currloc=1; currloc < interaction_lists.at(num_loci).at(interact).size(); currloc++){
          params.z[curr_z_index] *= 
//             geno_convert[strata_inds[0][include_strata[stratum]]->get_genotype(loci[interaction_lists[num_loci][interact][currloc]])];
            geno_convert.at(ref_alleles.at(interaction_lists.at(num_loci).at(interact).at(currloc))).at(strata_inds.at(0).at(include_strata.at(stratum))->get_genotype(loci.at(interaction_lists.at(num_loci).at(interact).at(currloc))));
        }
        curr_z_index++;
      }
    }
  }
}

//     struct conditional_lr_parameters{
//       int * z;
//       int * nca;
//       int * nct;
//       int * ivar;
//       float * b;
//       int ns, nimax, nmax, nvmax, nvmax1, nv;
//       float* cov1, *cntr, *w, *wb, *wdb, *wd2b, *u, *ins, *db, *d2b, *dl, *cov;
//       int * 
//       float chi2, st;
//       int ifault;
//       
//     };

///
/// Frees memory used for dynamic parameters in conditional LR
/// @param params conditional_lr_parameters&
/// @return
///
void CondLogRegression::free_params(conditional_lr_parameters& params){

//  delete [] params.cov1;
  delete [] params.wb;
  delete [] params.u;
  delete [] params.ins;
  delete [] params.db;
  delete [] params.d2b;
  delete [] params.dl;
  delete [] params.nca;
  delete [] params.nct;
  delete [] params.b;
  delete [] params.cov;
  delete [] params.ivar;
  
  delete [] params.z;
  delete [] params.covi;
  delete [] params.wdb;
  delete [] params.wd2b;
  
  delete [] params.w;
  delete [] params.cntr;
  
}


///
/// Returns p value based on chi square scored passed and
/// degrees of freedom
/// @param x
/// @param n
/// @returns p value
///
double CondLogRegression::ChiSq(double x, unsigned int n) {
//     if(x>1000 | n>1000) { var q=Norm((Power(x/n,1/3)+2/(9*n)-1)/Sqrt(2/(9*n)))/2; if (x>n) {return q} else {return 1-q} }
//     var p=Math.exp(-0.5*x); if((n%2)==1) { p=p*Math.sqrt(2*x/Pi) }
//     var k=n; while(k>=2) { p=p*x/k; k=k-2 }
//     var t=p; var a=n; while(t>1e-15*p) { a=a+2; t=t*x/a; p=p+t }
//     return 1-p
    
// cout << "x=" << x << " n=" << n << endl;    
  if(x > 1000 || n>1000){
// cout << "argument passed to norm from ChiSq=" << (pow(x/n,1/3)+2/(9*n)-1)/sqrt(2/(9*n)) << endl;
// cout << "pow(x/n,1/3)+2/(9*n)-1 = " << pow(x/n,1/3)+2/(9*n)-1 << endl;
// cout << "sqrt(2/(9*n)=" << sqrt(2/(9*n)) << endl;
    double q = norm((pow(x/n,1/3)+2/(9*n)-1)/sqrt(2/(9*n)))/2;
// cout << "q=" << q << endl;
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
/// Calculates p value
/// @param z
/// @return p value
///
double CondLogRegression::norm(double z){
  double q = z * z;
  if(fabs(z)>7){
    return (1-1/q+3/(q*q))*exp(-q/2)/(fabs(z)*sqrt(PiD2));
  }
  else{
    return ChiSq(q, 1);
  }
}

}
