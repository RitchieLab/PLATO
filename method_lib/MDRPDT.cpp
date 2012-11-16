//MDRPDT.cpp

#include "MDRPDT.h"

namespace Methods{

///
/// Constructor
///
MDRPDT::MDRPDT(){
  initialize();
}


///
/// Alternative constructer that establishes DataSet
/// @param ds DataSet
///
MDRPDT::MDRPDT(DataSet* ds){
  initialize();
  resetDataSet(ds);
}

///
/// Destructor
///
MDRPDT::~MDRPDT(){
  if(eval != NULL)
    delete eval;
  if(dist)
    delete dist;
}


///
/// Initializes
///
void MDRPDT::initialize(){
  num_ptests = 0;
  random_seed = 7;
  xvCount = 1;
  eval = NULL;
  dist = NULL;
  dataset = NULL;
  avg_mor = 0.0;
  ptests_run = false;
}


///
/// set parameters for method using StepOptions class
///
void MDRPDT::set_parameters(StepOptions* options){
  num_ptests = options->getMDRPDTNumPTests();
  random_seed = options->getMDRPDTRandSeed();
  xvCount = options->getMDRDPTNumCrossvals();
  Utility::Random rnd(random_seed);
  MdrPDT::Evaluation::EvaluationMethod::maxModelSize = options->getMDRPDTMaxCombo();
  MdrPDT::Evaluation::EvaluationMethod::minModelSize = options->getMDRPDTMinCombo();
  if(num_ptests > 0 && !ptests_run && dataset != NULL){
    if(!params_set){
      Load(dataset);
    }    
    RunPTests();
    ptests_run = true;
  }
  params_set = true;
}


///
/// Converts original marker map index into index of the 
/// converted data in the repository
/// @param locus in marker map position
/// @return index in repository
///
int MDRPDT::index_conversion(int loc){
  return repo.getNewMarkerIndex(loc);
}

///
/// runs MDR-PDT on the locus index passed 
/// @param locus index in marker_map position to run PDT on
///
void MDRPDT::calculate(unsigned int locus){
  MdrPDT::PdtModel currModel;
  int repo_index = index_conversion(locus);
  currModel.SetSnpLabels(data.GetSnpLabels());
  currModel.SetDSPCounts(data.GetDSPCounts());
  // original has GetSnp(snpID-1) but as the repo_index is zero-based i don't think it is needed
  currModel.AddSnp(repo_index, data.GetSNP(repo_index)); 
  runPDT(currModel);
}


///
/// runs MDR-PDT on the loci passed in the vector
/// @param loci vector of indexes in marker_map order to run PDT on
///
void  MDRPDT::calculate(vector<unsigned int>& loci){

  MdrPDT::PdtModel currModel;
  vector<unsigned int>::iterator iter;
  int repo_index;
  
  currModel.SetSnpLabels(data.GetSnpLabels());
  currModel.SetDSPCounts(data.GetDSPCounts());
  
  for(iter = loci.begin(); iter != loci.end(); ++iter){
    repo_index = index_conversion(*iter);
    currModel.AddSnp(repo_index+1, data.GetSNP(repo_index)); 
  }
  runPDT(currModel);
}


///
/// Runs PDT on the model.
/// @param model PdtModel to run
///
void MDRPDT::runPDT(MdrPDT::PdtModel& model){
  eval->EvaluateModel(&data, model);
  
//   t_stat = model.GetTrainingT(0);
  t_stats_training.clear();
  mor_values.clear();
  for(int cv=0; cv < xvCount; cv++){
    t_stats_training.push_back(model.GetTrainingT(cv));
    MdrPDT::Evaluation::MatchedOddsRatio mor = model.EvaluateMOR(cv);
    mor_values.push_back(mor.GetRatio());
  }
  if(xvCount > 1){
    t_stats_testing.clear();
    for(int cv=0; cv < xvCount; cv++){
      t_stats_testing.push_back(model.GetTestingT(cv));
    }
  }
  
}


///
/// Returns average Matched Odds Ratio for the model
/// @return average Matched Odds Ratio
///
float MDRPDT::getAvgMatchedOddsRatio(){
  if(mor_values.empty())
    return 0.0;

  vector<float>::iterator iter;
  float total=0.0;
  for(iter = mor_values.begin(); iter != mor_values.end(); iter++){
    total += *iter;
  }

  return total/mor_values.size();

}



///
/// Loads current dataset into repositories for 
/// analysis.  The Markers need to be labelled for inclusion.
/// @param ds DataSet containing current families and data
/// 
void MDRPDT::Load(Methods::DataSet* ds){
  
  Utility::Random rnd(random_seed);
  
  // call load on familiy repository
  repo.Load(ds);
  
  repo.PostLoad();
  
  // initialize the data -- don't randomize
  repo.InitializeData(data, rnd, xvCount, false);

  if(eval != NULL)
    delete eval;
  eval = new MdrPDT::Evaluation::TStatistic(xvCount, repo.GetPedigreeCount(), data.GetIndividualCount());
  
}

///
/// Runs p-tests and creates distribution.
///
void MDRPDT::RunPTests(){
  if(dist)
    delete dist;
  
  dist = new MdrPDT::Distribution::OmnibusDistribution(num_ptests);
  
//   dist = new MdrPDT::Distribution::NTestDistribution(num_ptests);
  
  MdrPDT::Evaluation::TStatistic eval(xvCount, repo.GetPedigreeCount(), data.GetIndividualCount());
  
  int testNumber = 1;
  
  MdrPDT::Evaluation::EvaluationMethod::maxModelSize=2;
  vector<string> exclusions;
  
  
  while(testNumber <= num_ptests){
    Utility::Random rnd(random_seed+testNumber);
    MdrPDT::GenotypeRepository test;
    
    // create the test data partition
    repo.InitializeData(test, rnd, xvCount, true);

    test.InitExclusionList(exclusions);
    
		MdrPDT::Evaluation::TStatistic testEval(xvCount, repo.GetPedigreeCount(), test.GetIndividualCount());
		MdrPDT::Evaluation::TFinalReport ptestWinner(MdrPDT::Evaluation::EvaluationMethod::maxModelSize, xvCount, 1);
		
		testEval.BasicEval(&test, ptestWinner);
		
		string bestID = "";
		float bestScore = 0.0;
		for (int o=MdrPDT::Evaluation::EvaluationMethod::minModelSize-1; o<MdrPDT::Evaluation::EvaluationMethod::maxModelSize;o++) {
			string modelID;
			float curScore = ptestWinner.GetBestMOR(o, modelID);
			if (curScore>bestScore) {
				bestScore=curScore;
				bestID = modelID;
			}
		}
		
    dist->AppendTest(1, testNumber-1, bestID.c_str(), bestScore);
    testNumber++;
  }
  
}

};
