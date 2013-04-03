//Interactions.cpp

// March 7, 2013

#include "Interactions.h"
#include "LinRegression.h"
#include "LogisticRegression.h"
#include <gsl/gsl_cdf.h>
#include "Helpers.h"
#include "AlleleFrequency.h"

namespace Methods{
///
/// basic method
/// @param ds DataSet
///
void Interactions::calculate(DataSet* ds){
  
  data_set = ds;
  
  ofstream inter_out;
  
  lrt_threshold = options.getLRTPval();
  
  // create appropriate type of regression calculator
  
  if (opts::_BINTRAIT_)
	{
	  regressor = new LogisticRegression;  
	  openOutput(inter_out, false);
	}
	else
	{
	  regressor = new LinRegression;
	  openOutput(inter_out, true);
	}
	
	// set covariates if used in model
	SetCovariates();
	
	regressor->resetDataSet(data_set);
	regressor->set_parameters(&options);
	
  // check if lists are from bio file or exhaustive and run calculation
  if(options.getBioSnpFile() != "")
	{
	
	  CalculateBioFile(inter_out, options.getBioSnpFile());
	}
	else
	{
	  CalculateExhaustive(inter_out);
	}
  
  // free memory
  delete regressor;
}


///
/// Set covariates for use in regression models
///
void Interactions::SetCovariates(){
  vector<string> use_covs = options.getCovars();

  if(options.doCovarsName()){
    // convert to indexes and add to covars
    for(unsigned int i =0; i<use_covs.size(); i++){
      covars.push_back(data_set->get_covariate_index(use_covs[i]));
    }
  }
  else{
    // convert string numbers to numbers
    for(unsigned int i =0; i<use_covs.size(); i++)
      covars.push_back(atoi(use_covs[i].c_str())-1);
  }
  
}



///
/// Take SNP pairings from bio file
/// @param inter_out ostream to write to 
/// @param biofiltername name of the biofilter file
///
void Interactions::CalculateBioFile(ostream& inter_out, string biofiltername){
  
  ofstream epi_log;
  openLog(epi_log);
  
  // retrieve the SNP pairs
  options.readBioTextFile(biofiltername);
  map<string, vector<string> > pairs = options.getBioPairs();
  
  MarkerInfo mark1, mark2;
  
  uni_results.clear(); // clear any previous single SNP results
  
  for(map<string, vector<string> >::iterator bio_iter = pairs.begin(); bio_iter != pairs.end(); 
    bio_iter++){
    
    string snp1 = bio_iter->first;

    vector<string> snp2_list = bio_iter->second;
    
    if(!getMarker(snp1, mark1, epi_log))
      continue;
    
    
    for(vector<string>::iterator snp2_iter=snp2_list.begin(); snp2_iter != snp2_list.end();
      snp2_iter++){
      
      if(!getMarker(*snp2_iter, mark2, epi_log))
        continue;
      
      if(mark1.loc_index == mark2.loc_index){
			  epi_log << snp1 << " " << *snp2_iter << " ---> " << "Snps have same index!" << endl;
				continue;        
      }
      
      CalculatePair(mark1, mark2, inter_out);
    }
  }
  
  epi_log.close();
}


///
/// Calculate all pairwise SNP interactions
/// @param inter_out ostream to write to
///
void Interactions::CalculateExhaustive(ostream& inter_out){
  
  vector<int> good_indexes = Helpers::findValidMarkersIndexes(data_set->get_markers(), &options);
  
  ofstream epi_log;
  openLog(epi_log);
  
  MarkerInfo mark1, mark2;
  uni_results.clear();
  
  for(vector<int>::iterator snp1_iter = good_indexes.begin(); snp1_iter != good_indexes.end(); ++snp1_iter){
    if(!getMarker(*snp1_iter, mark1, epi_log))
      continue;   
    for(vector<int>::iterator snp2_iter = snp1_iter+1; snp2_iter != good_indexes.end(); ++snp2_iter){
      if(!getMarker(*snp2_iter, mark2, epi_log))
        continue;
      CalculatePair(mark1, mark2, inter_out);
      
    }
  }
  
  epi_log.close();
  
}


///
/// Calculate the single SNP models, the reduced model and the full model for the pair
/// of SNPs passed into this function.
/// @param snp1 MarkerInfo first marker
/// @param snp2 MarkerInfo second marker
/// 
void Interactions::CalculatePair(MarkerInfo& snp1, MarkerInfo& snp2, ostream& inter_out){

  // output results -- can add checks throughout for thresholds

  // get snp1 results
  UniRegression snp1_results, snp2_results;
  snp1_results = GetSingleRegression(snp1.loc_index);
  snp2_results = GetSingleRegression(snp2.loc_index);
  
  vector<unsigned int> snps;
  bool snp1_invalid = isnan(snp1_results.p_value);
  bool snp2_invalid = isnan(snp2_results.p_value);
  if(!snp1_invalid)
    snps.push_back(snp1.loc_index);
  if(!snp2_invalid)
    snps.push_back(snp2.loc_index);

  double lrt_p_value = 1.0;
  double red_p_value, red_rsq, red_llr, full_p_value,
    full_rsq, full_llr, likelihood_ratio;
  vector<double> red_coeff_p, red_beta, red_se, full_coeff_p, full_beta, full_se;
  
  if(!snps.empty()){
  // calculate reduced model
  regressor->setIncludeInteractions(false);
  regressor->calculate(snps, covars);
  red_coeff_p = regressor->getCoeffPValues();
  red_beta = regressor->getCoefficients();
  red_se = regressor->getCoeffStandardErr();
  red_p_value = regressor->getOverallP();
  red_rsq = regressor->adjusted_rsquared();
  red_llr = regressor->getLLR();
 
  // calculate full model
  if(snp1_invalid || snp2_invalid)
    regressor->setIncludeInteractions(false);
  else
    regressor->setIncludeInteractions(true);
  regressor->calculate(snps, covars);
  full_coeff_p = regressor->getCoeffPValues();
  full_beta = regressor->getCoefficients();
  full_se = regressor->getCoeffStandardErr();
  full_p_value = regressor->getOverallP();
  full_rsq = regressor->adjusted_rsquared();
  full_llr = regressor->getLLR();

  likelihood_ratio = -2 * (red_llr - full_llr);

  lrt_p_value = GetLLRPValue(likelihood_ratio);
  }
  
  if(lrt_p_value > lrt_threshold){
    return;
  }

  string sep=" ";
  string invalid = "NA";
  inter_out << snp1.marker->getRSID() << "_" << snp2.marker->getRSID() << sep <<
    snp1.marker->getRSID() << sep << snp2.marker->getRSID() << sep <<
    snp1.marker->getChrom() << ":" << snp1.marker->getBPLOC() << sep <<
    snp1_results.maf << sep << snp1_results.ngenotypes << sep <<
    snp2.marker->getChrom() << ":" << snp2.marker->getBPLOC() << sep <<
    snp2_results.maf<< sep << snp2_results.ngenotypes << sep;
  
   if(!snp1_invalid && !snp2_invalid)
     inter_out << snp1_results.p_value << sep << snp1_results.beta << sep << snp1_results.se << sep <<
      snp2_results.p_value << sep << snp2_results.beta << sep << snp2_results.se << sep <<
      red_coeff_p[0] << sep << red_beta[0] << sep << red_beta[0] << sep <<
      red_coeff_p[1] << sep << red_beta[1] << sep << red_beta[1] << sep <<
      full_coeff_p[0] << sep << full_beta[0] << sep << full_se[0] << sep <<
      full_coeff_p[1] << sep << full_beta[1] << sep << full_se[1] << sep <<
      full_coeff_p[2] << sep << full_beta[2] << sep << full_se[2] << sep <<
      red_p_value << sep << red_rsq << sep << 
      full_p_value << sep << full_rsq << sep <<
      full_rsq - red_rsq << sep <<
      lrt_p_value << endl;    
    else if(snp1_invalid && snp2_invalid)
      inter_out << invalid << sep << invalid<< sep << invalid<< sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      invalid << sep << invalid<< sep << invalid << sep <<
      invalid << sep << invalid << sep << invalid<< sep <<
      invalid << sep << invalid<< sep << invalid << sep <<
      invalid << sep << invalid<< sep << invalid << sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      invalid << sep << invalid << sep << 
      invalid<< sep << invalid<< sep <<
      invalid << sep <<
      invalid << endl;
    else if(snp1_invalid)
     inter_out << invalid << sep << invalid << sep << invalid << sep <<
      snp2_results.p_value << sep << snp2_results.beta << sep << snp2_results.se << sep <<
      invalid << sep << invalid<< sep << invalid << sep <<
      red_coeff_p[0] << sep << red_beta[0] << sep << red_beta[0] << sep <<
      invalid<< sep << invalid << sep << invalid << sep <<
      full_coeff_p[0] << sep << full_beta[0] << sep << full_se[0] << sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      red_p_value << sep << red_rsq << sep << 
      full_p_value << sep << full_rsq << sep <<
      full_rsq - red_rsq << sep <<
      lrt_p_value << endl;    
    else //snp2 invalid
      inter_out << snp1_results.p_value << sep << snp1_results.beta << sep << snp1_results.se << sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      red_coeff_p[0] << sep << red_beta[0] << sep << red_beta[0] << sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      full_coeff_p[0] << sep << full_beta[0] << sep << full_se[0] << sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      invalid << sep << invalid<< sep << invalid << sep <<
      red_p_value << sep << red_rsq << sep << 
      full_p_value << sep << full_rsq << sep <<
      full_rsq - red_rsq << sep <<
      lrt_p_value << endl;
}


///
/// Calculate p value for likelihood ratio test
///
double Interactions::GetLLRPValue(double llr){
//   return gsl_cdf_chisq_P(llr,1);
  return  gsl_cdf_chisq_Q(llr,1);
}



///
/// Calculates single SNP when necessary
///  @param snp1 locus index for the snp
///
Interactions::UniRegression Interactions::GetSingleRegression(int snp_index){
  if(uni_results.find(snp_index) == uni_results.end())
  {
    regressor->setIncludeInteractions(false);
    vector<unsigned int> snps;
    snps.push_back(snp_index);
    // run this SNP
    UniRegression results;
    regressor->calculate(snps, covars);
    
    results.p_value = (regressor->getCoeffPValues())[0];
    results.beta = (regressor->getCoefficients())[0];
    results.se = ( regressor->getCoeffStandardErr())[0];
    results.ngenotypes = regressor->getNumGenotypes();
    results.maf = calcMAF(snp_index);

    uni_results[snp_index] = results;
  }

  return uni_results[snp_index];
}




///
///  Get marker information based on name
/// @param name of SNP to find
/// @param m MarkerInfo to fill 
/// @param epi_log log to write information
/// @return true if marker is found and good to use
///
bool Interactions::getMarker(string name, MarkerInfo & m, ostream& epi_log){
    
  try{
    m.loc_index = data_set->get_locus_index(name);
	  m.marker = data_set->get_locus(m.loc_index);
	}
	catch(MethodException & ex){
	  epi_log << name << " not found!" << endl;
	  return false;
	}
	
	if(!m.marker->isEnabled()){
	  epi_log << name << " is disabled!" << endl;
	  return false;
	}  
	
	if(opts::_CHRX_ == m.marker->getChrom()){
	  epi_log << name << " is on Chr X and is skipped!" << endl;
	  return false;
	}

  return true;
}


///
///  Get marker information based on name
/// @param name of SNP to find
/// @param m MarkerInfo to fill 
/// @param epi_log log to write information
/// @return true if marker is found and good to use
///
bool Interactions::getMarker(int index, MarkerInfo & m, ostream& epi_log){
    

    m.loc_index = index;
	  m.marker = data_set->get_locus(m.loc_index);
	
	if(!m.marker->isEnabled()){
	  epi_log << m.marker->getRSID() << " is disabled!" << endl;
	  return false;
	}  
	
	if(opts::_CHRX_ == m.marker->getChrom()){
	  epi_log <<  m.marker->getRSID() << " is on Chr X and is skipped!" << endl;
	  return false;
	}

  return true;
}

///
/// Opens output and applies headers to output
///
void Interactions::openOutput(ofstream & out, bool isLinearReg){

  // check for overwriting existing file
  string fname = opts::_OUTPREFIX_ + "interaction" + options.getOut() + ".txt";
  if(!overwrite)
  {
      fname += "." + getString<int>(order);
  }  
  out.open(fname.c_str());
  if(!out)
  {
    throw MethodException("Unable to open " + fname + " for writing!\n");
  }
  
  string rsq;
  if(isLinearReg)
    rsq = "rsq";
  else
    rsq = "pseudo_rsq";
  
  // add headers to output file
  string sep = " ";
  
  out << "Model_ID" << sep << "Var1_ID" << sep << "Var2_ID" << sep << "Var1_Pos" << 
    sep << "Var1_MAF" << sep << "Var1_#" << sep << "Var2_Pos" << sep << "Var2_MAF" <<
    sep << "Var2_#" << sep << "Uni_Var1_P" << sep << "Uni_Var1_beta" << sep <<
    "Uni_Var1_SE" << sep << "Uni_Var2_P" << sep << "Uni_Var2_beta" << sep << 
    "Uni_Var2_SE" << sep << "Red_Var1_P" << sep << "Red_Var1_beta" << sep <<
    "Red_Var1_SE" << sep << "Red_Var2_P" << sep << "Red_Var2_beta" << sep <<
    "Red_Var2_SE" << sep << "Full_Var1_P" << sep << "Full_Var1_beta" << sep <<
    "Full_Var1_SE" << sep << "Full_Var2_P" << sep <<  "Full_Var2_beta" << sep <<
    "Full_Var2_SE" << sep << "Full_INT_P" << sep << "Full_INT_beta" << sep <<
    "Full_INT_SE" << sep << "Red_Model_P" << sep << "Red_Model_" << rsq << sep <<
    "Full_Model_P" << sep << "Full_Model_" << rsq << sep << rsq << "_diff" << sep << 
    "LRT_P" << endl;
  
}

/// 
/// Calculate MAF
/// 
///
double Interactions::calcMAF(int marker_index){

  AlleleFrequency af;
	af.resetDataSet(data_set);
	af.initializeCounts(0);
	
	bool useoverall=false;
	
	if (options.doRandomChild() || options.doAll() || options.doAllChildren()
			|| options.doUnaffSpousesOnly() || options.doUnknownSpouses()) {
		af.setOptions(options);
	} else {
		options.setFoundersOnly();
		af.setOptions(options);
	}
	if(options.doAll() || options.doFilterOverall()){
		useoverall = true;
	}
	
	Marker* mark = data_set->get_locus(marker_index);
	af.calcOne(mark);
	float minfreq;
	
	if (!mark->isMicroSat()){
    float majfreq = 0.0f;
		if (useoverall) {
		  majfreq = af.getAone_freq();
		} else {
		  majfreq = af.getAoneP_freq();
		}
		if(majfreq > 0.5) 
  	  minfreq = 1.0f - majfreq;
  	else
  	  minfreq = majfreq;
		
	}
	
	return minfreq;
}


///
/// Return stream to log file 
/// @return ostream
///
void Interactions::openLog(ofstream& epi_log){
  string fname = opts::_OUTPREFIX_ + "interaction" + options.getOut() + ".txt";
  if(!overwrite)
  {
        fname += "." + getString<int>(order);
  }

  string logname = fname + ".log";
  epi_log.open(logname.c_str());
  
  epi_log << "LRT threshold is " << lrt_threshold << endl;
  
  if(!epi_log)
	{
		throw MethodException("Cannot open " + logname + " for writing!");
	}

}

} // end namespace Methods

