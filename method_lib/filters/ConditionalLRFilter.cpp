//ConditionalLRFilter.cpp

#include "ConditionalLRFilter.h"
#include <math.h>
#include <ComboGenerator.h>
#include <ChiSquare.h>

using namespace Methods;

namespace Filters{

///
/// constructor -- initialize variables
///
ConditionalLRFilter::ConditionalLRFilter():Filter("ConditionalLR"){
  initialize();
}


///
/// Alternative constructor
/// @param filterName name of filter
///
ConditionalLRFilter::ConditionalLRFilter(string filterName):Filter(filterName){
  initialize();
}


///
/// initializes maps used in switches for setting parameters
/// @return none
///
void ConditionalLRFilter::initialize(){
  pThreshold = 1.0;

  pValType = Overall;
  coeff_output = false;

//   ConfigMap["LOGREGRESSTYPE"] = CalcType;
  ConfigMap["THRESHOLD"] = Threshold;
  ConfigMap["MODELTYPE"] = ModelType;
//   ConfigMap["COMBOLIMIT"] = ComboLimit;
//   ConfigMap["COMBOMIN"] = ComboMin;
  ConfigMap["MAXITERATIONS"] = MaximumIterations;
  ConfigMap["INCLUDEINTERACTION"] = InteractionIncluded;
  ConfigMap["PTYPE"] = PType;
  ConfigMap["OUTPUTCOEFF"] = IncludeCoeff;
  ConfigMap["MAXMODELSIZE"] = MaxModelSize;

  PValueTypeMap["OVERALL"] = Overall;
//   PValueTypeMap["COEFFICIENT"] = Coefficient;
  PValueTypeMap["NONE"] = None;

  
  TotalTypeDF["DOMINANT"] = 1;
  TotalTypeDF["RECESSIVE"] = 1;
  TotalTypeDF["ADDITIVE"] = 2;
  
  max_model_size = 2;

}


///
/// Returns estimated run time in seconds 
/// @param num_models Number of models that will be processed total
/// @param dataset DataSet
/// @param resultList ResultList contains list of models to test as time trial
/// @return ProcessEstimate that lists estimated run time and number of models
/// that will pass threshold
///
Filter::ProcessEstimate ConditionalLRFilter::estimate_run_time(double num_models, DataSet & dataset, 
  ResultSet& resultList, int model_size){
  Filter::ProcessEstimate time_estimate;
  clock_t start, end;

  start = clock();
  lr_calculator.resetDataSet(&dataset);
  end = clock();
  
  clock_t reset_time = end - start;
 
  size_t num_run = resultList.size();

  start = clock();
  analyze(resultList, dataset);

  end = clock();

  time_estimate.seconds = double(end - start - reset_time) / (CLOCKS_PER_SEC * num_run);

  time_estimate.seconds *= num_models;
  
  // calculate total that will be passed through
  // depends on threshold (converted to p value if not already)
  float pthresh = pThreshold;
  // need to determine the p value threshold when threshold is not
  // the 
  if(pValType == None){
    int df;
    if(model_size == 1){
     df = 1;
    }
    else{
    // degrees of freedom should be determined by the number of snps and then all
    // the combinations included in that 
      ComboGenerator gen;
      gen.SetLoci(model_size);
      gen.ComboEnds(1,model_size);
      gen.GenerateCombinations();
      df = int(gen.ComboList.size());
    }
    pthresh = ChiSquare::pfromchi(pThreshold, df);
  }
  time_estimate.num_models = num_models * pthresh;
  
  return time_estimate;
}




///
/// Runs logistic regression on dataset passed <br>
/// @param resultList ResultSet that may contain results from other filters
/// and will contain results of this sequential replication analysis
/// @param dataset DataSet to be analyzed
/// @return none
/// @throws FilterExcept when no sequential replication type has been set
///
void ConditionalLRFilter::analyze(ResultSet & resultList, DataSet & dataset){

  // calculate position for p-value or likelihood ratio
  // will always be next value
  unsigned int threshold_index = 0;
  if(resultList.size() > 0){
    threshold_index = resultList.begin()->analysisScores.size();
  }
  bool erase_result = false;
  float result_value;

  lr_calculator.resetDataSet(&dataset);

  for(ResultIter currResult = resultList.begin(); currResult !=resultList.end();){
    erase_result = false;

    lr_calculator.calculate(currResult->genoCombination);
    if(pValType == Overall){
      result_value = lr_calculator.getOverallP();
    }
//     else if(pValType == Coefficient){
//       result_value = lr_calculator.getFullInteractionP();
//     }
    else{
      result_value = lr_calculator.getLogLikelihoodRatio();
    }

    if(pValType != None){
      if(result_value > pThreshold){
        erase_result = true;
      }
    }
    else
      if(result_value < pThreshold)
        erase_result = true;

    if(!erase_result){
      currResult->analysisScores.push_back(result_value);

     // add covariates and coefficients if requested
     if(coeff_output){

       vector<double>::iterator iter;
       vector<double> values = lr_calculator.getCoefficients();
       for(iter = values.begin(); iter != values.end(); iter++)
         currResult->analysisScores.push_back(*iter);

//        values = lr_calculator.getCovariates();
//        for(iter = values.begin(); iter != values.end(); iter++)
//          currResult->analysisScores.push_back(*iter);


     }

      ++currResult;
    }
    else
      resultList.erase(currResult++);
  }

}



///
/// Sets parameters for the filter and throws a FilterExcept when
/// a parameter is unacceptable
/// @param  params map with key being parameter identifier and value being the value for that param
/// @return  none
/// @throws FilterExcept when required parameter is not set or parameter is out of bounds
///
void ConditionalLRFilter::set_params(PARAMS params, DataSet* set){

std::map<string, string>::iterator configIter;
  string value;
  unsigned int maxIterations;

  StepOptions options;

  // check every entry and throw exception when unknown entry encountered
  for(configIter = params.begin(); configIter != params.end(); configIter++){
    switch(ConfigMap[configIter->first]){
      case Threshold:
        pThreshold = Stringmanip::stodouble(configIter->second);
        if(pThreshold < 0 || pThreshold > 1){
          throw FilterExcept(configIter->first +
            " must be greater than zero and less than or equal to one for " +
            name + " filter\n\n");
        }
        break;
      case ModelType:
        modType = Stringmanip::to_upper(configIter->second);
        options.setCondLRModelType(modType);
//         if(ModelTypeMap.find(value) != ModelTypeMap.end()){
//           lr_calculator.setModelType(ModelTypeMap[value]);
//         }
//         else{
//           throw FilterExcept(configIter->second + " is not a valid option for paramater " +
//             configIter->first + " for " + name + " filter\n\n");
//         }
        break;
      case PType:
        value = Stringmanip::to_upper(configIter->second);
        if(PValueTypeMap.find(value) != PValueTypeMap.end()){
          pValType = PValueTypeMap[value];
        }
        else{
          throw FilterExcept(configIter->second + " is not a valid option for paramater " +
            configIter->first + " for " + name + " filter\n\n");
        }
        break;
      case MaximumIterations:
        maxIterations = Stringmanip::stoi(configIter->second);
        if(maxIterations < 1){
          throw FilterExcept(configIter->second + " must be greater than zero for " +
            name + " filter\n\n");
        }
//         lr_calculator.setMaximumIterations(maxIterations);
        options.setCondLRMaximumIterations(maxIterations);
        break;
      case IncludeCoeff:
        value = Stringmanip::to_upper(configIter->second);
        if(value.compare("TRUE") == 0 || value.compare("ON")==0)
          coeff_output = true;
        else
          coeff_output = false;
        break;
      case InteractionIncluded:
        value = Stringmanip::to_upper(configIter->second);
        if(value.compare("TRUE") == 0)
//           lr_calculator.setIncludeInteractions(true);
          options.setCondLRIncludeInteractions(true);
        else
//           lr_calculator.setIncludeInteractions(false);
          options.setCondLRIncludeInteractions(false);
        break;
      case MaxModelSize:
        max_model_size = Stringmanip::stoi(configIter->second);
        break;
      case NoMatch:
        throw FilterExcept(configIter->first + " is not defined as a valid parameter for " +
          name + " filter\n\n");
    };
  }
  lr_calculator.set_parameters(&options);

    // adjust name of this filter to reflect output of coefficients
  if(coeff_output){
 ////unused   int num_coeff;
    if(!options.getCondLRIncludeInteractions()){
      max_coeff = max_model_size;
    }
    else{
      Methods::ComboGenerator generator;
      generator.ComboEnds(1, max_model_size);
      generator.SetLoci(max_model_size);
      generator.SetComboInterval(1000);
 ////unused     bool done=true;
      generator.GenerateCombinations();
      max_coeff = int(generator.ComboList.size());
    }

    for(int i=0; i<max_coeff; i++){
      name += "\tcoeff ";
    }
    name += "\n";

  }

}

}




