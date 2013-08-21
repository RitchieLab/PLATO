//Interactions.cpp

// March 7, 2013

#include "Interactions.h"
#include "LinRegression.h"
#include "LogisticRegression.h"
#include <gsl/gsl_cdf.h>
#include "Helpers.h"
#include "AlleleFrequency.h"
#include<set>
#include<sstream> 

namespace Methods{
///
/// basic method
/// @param ds DataSet
///
void Interactions::calculate(DataSet* ds){
  
  data_set = ds;
  
  ofstream inter_out;
  
  lrt_threshold = options.getLRTPval();
  string gxeFile = options.getGXEFile();
  bool isGXE = !options.getGXEcovars().empty() || gxeFile != "";
  
  // create appropriate type of regression calculator
  if(PhenoBinary())
	{
	  regressor = new LogisticRegression;
	  openOutput(inter_out, false, isGXE);
	}
	else
	{
	  regressor = new LinRegression;
	  openOutput(inter_out, true, isGXE);
	}
	
	// set covariates if used in model
	SetCovariates();
	
	regressor->resetDataSet(data_set);
	regressor->set_parameters(&options);

  // check if this is a GXE analysis
  if(options.getGXEcovars().empty() and gxeFile == ""){
  	// check if lists are from bio file or exhaustive and run calculation
	  if(options.getBioSnpFile() != "")
		{
	  	CalculateBioFile(inter_out, options.getBioSnpFile());
		}
		else
		{
		  CalculateExhaustive(inter_out);
		}
	}
	else{
		if(gxeFile==""){
			CalculateGXE(inter_out);
		}
		else{
			CalculateGXEFile(inter_out, gxeFile);
		}
	}
  
  // free memory
  delete regressor;
}



///
/// Return true if phenotype is binary
/// @return true for binary phenotype
///
bool Interactions::PhenoBinary(){
  set<string> phenos;
  phenos.insert("0");
  phenos.insert("1");
  phenos.insert("2");
  phenos.insert("-9");
  phenos.insert("-1");
  phenos.insert("-99999");

  if(options.getUsePheno()){
    int index = options.getPhenoLoc();
    if(options.getPhenoName() != ""){
      index = data_set->get_trait_index(options.getPhenoName());
      if(index==-1){
      	throw MethodException("The phenotype " + options.getPhenoName() + " is not found in traits file\n");
      }
    }
    for(int i=0; i<data_set->num_inds(); i++){
      stringstream ss;
      ss << data_set->get_sample(i)->getPheno(index);
      if(phenos.find(ss.str()) == phenos.end())
        return false;
    }
  }
  else{
    for(int i=0; i<data_set->num_inds(); i++){
      stringstream ss;
      ss << data_set->get_sample(i)->getPheno();
      if(phenos.find(ss.str()) == phenos.end())
        return false;
    }
  }
  return true;
}


///
/// Set covariates for use in regression models
///
void Interactions::SetCovariates(){
  vector<string> use_covs = options.getCovars();
	
  if(options.doCovarsName()){
    // convert to indexes and add to covars
    for(unsigned int i =0; i<use_covs.size(); i++){
			int index = data_set->get_covariate_index(use_covs[i]);
			if(index != -1){
	      covars.push_back(data_set->get_covariate_index(use_covs[i]));
	    }
	    else{
	    	throw MethodException("\nERROR: " + use_covs[i] + " not found in covariate file\n\n");
	    }
    }
  }
  else{
    // convert string numbers to numbers
    for(unsigned int i =0; i<use_covs.size(); i++){
			if(use_covs.at(i).find('-') != string::npos){
				vector<string> range;
				General::Tokenize(use_covs.at(i), range, "-");
				int start = atoi(range[0].c_str())-1;
				int last = atoi(range[1].c_str())-1;
				for(int j=start; j<=last; j++){
					covars.push_back(j);
				}
			}
			else{
	      covars.push_back(atoi(use_covs[i].c_str())-1);
	    }
    }
  }
}


///
/// Take gene-environment pairings from GXE file
/// @param inter_out ostream to write to 
/// @param gxefilename name of the GXE file
///
void Interactions::CalculateGXEFile(ostream& inter_out, string gxefilename){
  
  ofstream epi_log;
  openLog(epi_log);
  
  // retrieve the GxE pairs
  // it is gene as the map key
  options.readGXETextFile(gxefilename);
  map<string, vector<string> > pairs = options.getGXEPairs();
  
  MarkerInfo mark;
  // all models will have same covariates 
	modelCovars.push_back(0);
	modelCovars.insert(modelCovars.end(), covars.begin(), covars.end());
  uni_results.clear(); // clear any previous single SNP results
  int envIndex;
  
  map<string, vector<string> >::iterator endIter = pairs.end();
  for(map<string, vector<string> >::iterator gxeIter = pairs.begin(); gxeIter != endIter; 
    gxeIter++){
    
    string snp = gxeIter->first;

    vector<string>& envs = gxeIter->second;
    
    if(!getMarker(snp, mark, epi_log))
      continue;
    
    
    for(vector<string>::iterator envIter=envs.begin(); envIter != envs.end();
      envIter++){
      envIndex = data_set->get_covariate_index(*envIter);
    	modelCovars.front() = envIndex;
			CalculateGXEPair(mark, envIndex, inter_out);
    }
  }
  
  epi_log.close();
}



///
/// Sets covariate indexes from options
///
void Interactions::setGXECovars(){
	vector<string> gxe_covars = options.getGXEcovars();

	if(options.doGXEName()){
    // convert to indexes and add to covars
    for(unsigned int i =0; i<gxe_covars.size(); i++){
      gXecovars.push_back(data_set->get_covariate_index(gxe_covars[i]));
    }
  }
	else{
    for(unsigned int i =0; i<gxe_covars.size(); i++){
			if(gxe_covars.at(i).find('-') != string::npos){
				vector<string> range;
				General::Tokenize(gxe_covars.at(i), range, "-");
				int start = atoi(range[0].c_str())-1;
				int last = atoi(range[1].c_str())-1;
				for(int j=start; j<=last; j++){
					gXecovars.push_back(j);
				}
			}
			else{
	      gXecovars.push_back(atoi(gxe_covars[i].c_str())-1);
	    }
    }		
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
/// Calculate all pairwise combinations of the environmental factors (covariates)
/// and the SNPs
///
void Interactions::CalculateGXE(ofstream& inter_out){
	vector<int> good_indexes = Helpers::findValidMarkersIndexes(data_set->get_markers(), &options);
	ofstream epi_log;
	openLog(epi_log);
	
	MarkerInfo mark;
	uni_results.clear();
	cov_results.clear();
	
	setGXECovars();
	vector<unsigned int>::iterator endIter = gXecovars.end();
	vector<int>::iterator snpsEnd = good_indexes.end();

	// all models will have same covariates 
	modelCovars.push_back(0);
	modelCovars.insert(modelCovars.end(), covars.begin(), covars.end());

	for(vector<int>::iterator snpIter = good_indexes.begin(); snpIter != snpsEnd; ++snpIter){	
		if(!getMarker(*snpIter, mark, epi_log)){
				continue;
		}
		for(vector<unsigned int>::iterator cIter = gXecovars.begin(); cIter != endIter; ++cIter){
			// environmental factor will be first in the covariate list
			modelCovars.front() = *cIter;
			CalculateGXEPair(mark, *cIter, inter_out);
		}
	}
	
	epi_log.close();
}


///
/// Calculate the single SNP model, covariatae model, the reduced model and the full model for the GXE pair
/// @param snp MarkerInfo marker
/// @param environ Index of the covariate for this environmental factor
/// @param inter_out output stream
/// 
void Interactions::CalculateGXEPair(MarkerInfo& snp, int environ, ostream& inter_out){
	// get snp1 and snp2 results
  UniRegression snp_results, covar_results;
  snp_results = GetSingleRegression(snp.loc_index);
  covar_results = GetSingleEnvRegression(environ);
  
  // continue HERE 
  vector<unsigned int> snps, covariates;
  snp_results.valid = !isnan(snp_results.p_value);
  covar_results.valid = !isnan(covar_results.p_value);
  snps.push_back(snp.loc_index);
  
  ComplexResults results;
  if(snp_results.valid and covar_results.valid){
  	CalculateComplexResults(results, snps, modelCovars);
  	if(results.lrt_p_value > lrt_threshold){
    	return;
  	}
  }
  
  string sep=" ";
  string invalid = "NA";
  string environName = data_set->get_covariate_name(environ);
  
  inter_out << snp.marker->getRSID() << "_" << environName << sep <<
    snp.marker->getRSID() << sep << environName << sep <<
    snp.marker->getChrom() << ":" << snp.marker->getBPLOC() << sep <<
    snp_results.maf << sep << snp_results.ngenotypes << sep <<
    covar_results.ngenotypes << sep;  

  OutputPair(results, snp_results, covar_results, inter_out);  
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
/// Calculates full and reduced models and stores in the results parameter
/// @param results ComplexResults to store analysis
/// @param modsnps Indexes for model snps
/// @param modcovars Indexes for model covariates
///
void Interactions::CalculateComplexResults(ComplexResults& results, vector<unsigned int>& modsnps,
	vector<unsigned int>& modcovars){

  // calculate reduced model
  regressor->setIncludeInteractions(false);
  regressor->calculate(modsnps, modcovars);
  results.red_coeff_p = regressor->getCoeffPValues();
  results.red_beta = regressor->getCoefficients();
  results.red_se = regressor->getCoeffStandardErr();
  results.red_p_value = regressor->getOverallP();
  results.red_rsq = regressor->adjusted_rsquared();
  results.red_llr = regressor->getLLR();
 
  regressor->setIncludeInteractions(true);
  regressor->calculate(modsnps, modcovars);
  results.full_coeff_p = regressor->getCoeffPValues();
  results.full_beta = regressor->getCoefficients();
  results.full_se = regressor->getCoeffStandardErr();
  results.full_p_value = regressor->getOverallP();
  results.full_rsq = regressor->adjusted_rsquared();
  results.full_llr = regressor->getLLR();
//   results.likelihood_ratio = -2 * (results.red_llr - results.full_llr);
	results.lrt_p_value = GetLLRPValue(results.likelihood_ratio);
  results.likelihood_ratio = -(results.red_llr - results.full_llr);
  results.lrt_p_value = GetLLRPValue(results.likelihood_ratio);
}



///
/// Outputs results to the stream passed
/// @param complex
/// @param var1
/// @param var2
/// @param inter_out
///
void Interactions::OutputPair(ComplexResults& complex, UniRegression& var1, UniRegression& var2,
	ostream& inter_out){
	
	string sep=" ";
  string invalid = "NA";
   if(var1.valid && var2.valid)
     inter_out << var1.p_value << sep << var1.beta << sep << var1.se << sep <<
      var2.p_value << sep << var2.beta << sep << var2.se << sep <<
      complex.red_coeff_p[0] << sep << complex.red_beta[0] << sep << complex.red_se[0] << sep <<
      complex.red_coeff_p[1] << sep << complex.red_beta[1] << sep << complex.red_se[1] << sep <<
      complex.full_coeff_p[0] << sep << complex.full_beta[0] << sep << complex.full_se[0] << sep <<
      complex.full_coeff_p[1] << sep << complex.full_beta[1] << sep << complex.full_se[1] << sep <<
      complex.full_coeff_p[2] << sep << complex.full_beta[2] << sep << complex.full_se[2] << sep <<
      complex.red_p_value << sep << complex.red_rsq << sep << 
      complex.full_p_value << sep << complex.full_rsq << sep <<
      complex.full_rsq - complex.red_rsq << sep <<
      complex.lrt_p_value << endl;    
    else if(!var1.valid && !var2.valid)
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
    else if(!var1.valid)
     inter_out << invalid << sep << invalid << sep << invalid << sep <<
      var2.p_value << sep << var2.beta << sep << var2.se << sep <<
      invalid << sep << invalid<< sep << invalid << sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      invalid<< sep << invalid << sep << invalid << sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      invalid << sep << invalid << sep << 
      invalid << sep << invalid << sep <<
      invalid << sep <<
      invalid << endl;    
    else //var 2 invalid
      inter_out << var1.p_value << sep << var1.beta << sep << var1.se << sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      invalid<< sep << invalid << sep << invalid << sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      invalid << sep << invalid<< sep << invalid << sep <<
      invalid << sep << invalid << sep << invalid << sep <<
      invalid << sep << invalid<< sep << invalid << sep <<
      invalid << sep << invalid << sep << 
      invalid << sep << invalid << sep <<
      invalid << sep <<
      invalid << endl;
}


///
/// Calculate the single SNP models, the reduced model and the full model for the pair
/// of SNPs passed into this function.
/// @param snp1 MarkerInfo first marker
/// @param snp2 MarkerInfo second marker
/// 
void Interactions::CalculatePair(MarkerInfo& snp1, MarkerInfo& snp2, ostream& inter_out){

  // get snp1 and snp2 results
  UniRegression snp1_results, snp2_results;
  snp1_results = GetSingleRegression(snp1.loc_index);
  snp2_results = GetSingleRegression(snp2.loc_index);

  vector<unsigned int> snps;  
	snp1_results.valid = !isnan(snp1_results.p_value);
  snp2_results.valid = !isnan(snp2_results.p_value);  
  
  if(snp1_results.valid)
    snps.push_back(snp1.loc_index);
  if(snp2_results.valid)
    snps.push_back(snp2.loc_index);
  
  ComplexResults results;
  if(snp1_results.valid and snp2_results.valid){
	  CalculateComplexResults(results, snps, covars);
  	  if(results.lrt_p_value > lrt_threshold){
    		return;
		  }
	}
  
  string sep=" ";
  string invalid = "NA";
  inter_out << snp1.marker->getRSID() << "_" << snp2.marker->getRSID() << sep <<
    snp1.marker->getRSID() << sep << snp2.marker->getRSID() << sep <<
    snp1.marker->getChrom() << ":" << snp1.marker->getBPLOC() << sep <<
    snp1_results.maf << sep << snp1_results.ngenotypes << sep <<
    snp2.marker->getChrom() << ":" << snp2.marker->getBPLOC() << sep <<
    snp2_results.maf<< sep << snp2_results.ngenotypes << sep;
    
  OutputPair(results, snp1_results, snp2_results, inter_out);
}


///
/// Calculate p value for likelihood ratio test
///
double Interactions::GetLLRPValue(double llr){
  return  gsl_cdf_chisq_Q(llr,1);
}


///
///  Calculates single covariate regression when necessary
///  @param env_index locus index for the snp
///
Interactions::UniRegression& Interactions::GetSingleEnvRegression(int env_index){
  if(cov_results.find(env_index) == cov_results.end())
  {
    regressor->setIncludeInteractions(false);
    vector<unsigned int> snps;
    
    UniRegression results;
    regressor->calculate(snps, modelCovars);
    
    results.p_value = (regressor->getCoeffPValues())[0];
    results.beta = (regressor->getCoefficients())[0];
    results.se = ( regressor->getCoeffStandardErr())[0];
    results.ngenotypes = regressor->getNumGenotypes();

    cov_results[env_index] = results;
  }

  return cov_results[env_index];
}



///
/// Calculates single SNP when necessary
/// @param snp_index
///
Interactions::UniRegression& Interactions::GetSingleRegression(int snp_index){
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

  return true;
}

///
/// Opens output and applies headers to output
///
void Interactions::openOutput(ofstream & out, bool isLinearReg, bool isGXE){

	string fname;
	// check for overwriting existing file
  if(opts::_OUTPREFIX_.length() > 0)
  {
  	fname = opts::_OUTPREFIX_ + "interaction.txt";
	}
  else
  {
  	fname = "interaction" + options.getOut() + ".txt";
	}
	
  
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
    sep << "Var1_MAF" << sep << "Var1_#";
  if(!isGXE){
   out << sep<< "Var2_Pos" << sep << "Var2_MAF";
  }
  out << sep << "Var2_#" << sep << "Uni_Var1_P" << sep << "Uni_Var1_beta" << sep <<
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
  string fname;
  
	if(opts::_OUTPREFIX_.length() > 0)
  {
  	fname = opts::_OUTPREFIX_ + "interaction";
	}
  else
  {
  	fname = "interaction" + options.getOut();
	}
	
  
  if(!overwrite)
  {
        fname += "." + getString<int>(order);
  }

  string logname = fname + ".log";
  epi_log.open(logname.c_str());
  
  epi_log << "LRT threshold is " << lrt_threshold << endl;
  for(vector<unsigned int>::iterator iter=covars.begin(); iter != covars.end(); ++iter){
	epi_log << "Covariate included: " << data_set->get_covariate_name(*iter) << endl;
  }  
  
  if(!epi_log)
	{
		throw MethodException("Cannot open " + logname + " for writing!");
	}

}

} // end namespace Methods

