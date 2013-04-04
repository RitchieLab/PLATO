#ifndef INTERACTIONS_H
#define INTERACTIONS_H

#ifndef MAC
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
    };
    
    void setOverwrite(bool v){overwrite = v;}
    void setOrder(int o){order = o;}
    void setOptions(StepOptions o){options = o;}

  private:
  
    void CalculateBioFile(ostream& inter_out, string biofiltername);
    void CalculateExhaustive(ostream& out);
    void SetCovariates();
    
    bool PhenoBinary();
    map<int, UniRegression> uni_results;
    
    void CalculatePair(MarkerInfo& snp1, MarkerInfo& snp2, ostream& inter_out);
    double GetLLRPValue(double llr);
    UniRegression GetSingleRegression(int snp_index);
    
    bool getMarker(string name, MarkerInfo & m, ostream& epi_log);
    bool getMarker(int index, MarkerInfo & m, ostream& epi_log);
    
    double calcMAF(int marker_index);
    
    void openOutput(ofstream & out, bool isLinearReg);
    
    void openLog(ofstream& epi_log);
    
    vector<unsigned int> covars;
  
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
