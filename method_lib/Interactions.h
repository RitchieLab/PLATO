#ifndef INTERACTIONS_H
#define INTERACTIONS_H

#include "config.h"
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <vector>
#include <map>
#include "Marker.h"
#include "Family.h"
#include "Sample.h"
#include "Globals.h"
#include "Options.h"
#include "StepOptions.h"
#include "DataSet.h"
#include "Regression.h"


using namespace std;

namespace Methods{
  ///
  /// Runs regression tests on SNPS to evaluate the interaction models.
  /// It runs logistic regression for binary traits (case/control)
  /// and runs linear regression for continuous traits.
  ///
  class Interactions{
  
  public:
    Interactions(){
			data_set = NULL;
			families = NULL;
			markers = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			order = 0;
		};
		
		Interactions(DataSet* ds){
			data_set = ds;
			families = ds->get_families();
			markers = ds->get_markers();
			samples = ds->get_samples();
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			order = 0;
		}
    Interactions(int thresh) : threshold(thresh){
			data_set = NULL;
			families = NULL;
			markers = NULL;
			orig_num_markers = 0;
			orig_num_families = 0;
			orig_num_individuals = 0;
			rank = 0;
			_DBOUTPUT_ = false;
			order = 0;
		};
		void calculate(DataSet* ds);

    struct MarkerInfo{
      Marker * marker;
      int loc_index;
    };
  
    struct UniRegression{
      double p_value, beta, se, ngenotypes, maf;
      bool valid;
    };
    
    struct ComplexResults{
    	double lrt_p_value;
	  	double red_p_value, red_rsq, red_llr, full_p_value,
  	  full_rsq, full_llr, likelihood_ratio;
  		vector<double> red_coeff_p, red_beta, red_se, full_coeff_p, full_beta, full_se;
    };
    
    void setOverwrite(bool v){overwrite = v;}
    void setOrder(int o){order = o;}
    void setOptions(StepOptions o){options = o;}

  private:
  
    void CalculateBioFile(ostream& inter_out, string biofiltername);
    void CalculateGXEFile(ostream& inter_out, string gxefilename);
    void CalculateExhaustive(ostream& out);
    void CalculateGXE(ofstream& inter_out);
    void SetCovariates();
    void setGXECovars();
    
    bool PhenoBinary();
    map<int, UniRegression> uni_results;
    map<int, UniRegression> cov_results;
    
    void CalculatePair(MarkerInfo& snp1, MarkerInfo& snp2, ostream& inter_out, 
    	ostream& epi_log);
    void CalculateGXEPair(MarkerInfo& snp, int environ, ostream& inter_out, 
    	ostream& epi_log);
    void CalculateComplexResults(ComplexResults& results, vector<unsigned int>& modsnps,
    	vector<unsigned int>& modcovars);
    void OutputPair(ComplexResults& complex, UniRegression& var1, UniRegression& var2,
			ostream& inter_out);
    double GetLLRPValue(double llr);
    UniRegression& GetSingleRegression(int snp_index);
    UniRegression& GetSingleEnvRegression(int env_index);
    
    bool getMarker(string name, MarkerInfo & m, ostream& epi_log);
    bool getMarker(int index, MarkerInfo & m, ostream& epi_log);
    
    double calcMAF(int marker_index);
    
    void openOutput(ofstream & out, bool isLinearReg, bool isGXE);
    
    void openLog(ofstream& epi_log);
    
    vector<unsigned int> covars, gXecovars, modelCovars;
  
		DataSet* data_set;
		Regression * regressor;
		vector<Sample*>* samples;
		vector<Family*>* families;
		vector<Marker*>* markers;
		vector<int>* marker_map;
		double lrt_threshold;
		StepOptions options;

		int threshold;
		int orig_num_markers;
		int orig_num_families;
		int orig_num_individuals;
		int rank;
		bool _DBOUTPUT_;
		bool _MARKERLIST_;
		bool _STRATIFY_;
		bool overwrite;
		int order;
		vector<string> filenames;
  
  };
};

#endif
